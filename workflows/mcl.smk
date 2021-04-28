if 'edge_weight' not in config:
    config['edge_weight'] = 'pearson'

if 'correlation' not in config:
    config['correlation'] = 0.6

if 'inflation' not in config:
    config['inflation'] = [1.4]

if 'ceilnb' not in config:
    config['ceilnb'] = [10000]

if 'regular_size' not in config:
    config['regular_size'] = [15, 200]

if 'small_size' not in config:
    config['small_size'] = 3

if 'top_big' not in config:
    config['top_big'] = 3

if 'n_threads' not in config:
    config['n_threads'] = 1

if 'data_matrix' not in config:
    raise OSError(r"""Please specify the input data matrix for MCL clustering with:

            config['data_matrix'] = 'my-data-matrix.txt'

            or

            snakemake -C data_matrix=my-data-matrix.txt ...
            """)

if 'mcl_output_dir' not in config:
    raise OSError(r"""\nPlease specify the output directory for MCL clustering with:

            config['mcl_output_dir'] = 'my-mcl-output'

            or

            snakemake -C mcl_output_dir=my-mcl-output ...
            """)

rule make_mcl_network:
    input:
        mat= config['data_matrix'],
    output:
        network= temp(os.path.join(config['mcl_output_dir'], 'network.mci')),
        gene_dict= temp(os.path.join(config['mcl_output_dir'], 'gene_dict.txt')),
    params:
        n_threads= config['n_threads'],
        edge_weight= config['edge_weight']
    shell:
        r"""
        mcxarray -data {input.mat} -skipr 1 -skipc 1 -o {output.network} -t {params.n_threads} -co 0.2 \
            --write-binary --{params.edge_weight} -tf 'abs()' -write-tab {output.gene_dict}
        """

rule write_similarity_matrix:
    input:
        network= os.path.join(config['mcl_output_dir'], 'network.mci'),
        gene_dict= os.path.join(config['mcl_output_dir'], 'gene_dict.txt'),
    output:
        matrix= os.path.join(config['mcl_output_dir'], 'similarity_matrix.tsv.gz'),
    shell:
        r"""
        mcxdump -imx {input.network} -tab {input.gene_dict} --dump-upper | gzip > {output.matrix}
        """

rule query_mcl_network:
    input:
        network= os.path.join(config['mcl_output_dir'], 'network.mci'),
    output:
        qry= os.path.join(config['mcl_output_dir'], 'vary_correlation_query.txt'),
    shell:
        r"""
        mcx query -imx {input.network} --vary-correlation -o {output.qry}
        """

rule select_mcl_network:
    input:
        network= os.path.join(config['mcl_output_dir'], 'network.mci'),
    output:
        selected_ntwk= temp(os.path.join(config['mcl_output_dir'], 'correlation_{r}/selected_network.mci')),
    shell:
        r"""
        mcx alter -imx {input.network} -tf 'gq({wildcards.r}), add(-{wildcards.r})' \
            --write-binary -o {output.selected_ntwk}
        """

rule cluster_network:
    input:
        selected_ntwk= os.path.join(config['mcl_output_dir'], 'correlation_{r}/selected_network.mci'),
    output:
        clst= temp(os.path.join(config['mcl_output_dir'], 'correlation_{r}/inflation_{inflation}/ceilnb_{ceilnb}/cluster.mci')),
    params:
        nt= config['n_threads']
    shell:
        r"""
        mcl {input.selected_ntwk} -o {output.clst} -I {wildcards.inflation} -te {params.nt} -tf '#ceilnb({wildcards.ceilnb})'
        """

rule distance_between_clusters:
    input:
        clst= expand(os.path.join(config['mcl_output_dir'], 'correlation_{r}/inflation_{inflation}/ceilnb_{ceilnb}/cluster.mci'), 
                r= config['correlation'], inflation= config['inflation'], ceilnb= config['ceilnb']),
    output:
        dist= os.path.join(config['mcl_output_dir'], 'distance_between_clusters.tsv'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)

n_files <- strsplit('{input.clst}', ' ')
if(length(n_files) == 1) {{
    cat('\nOnly one cluster file found - skipping distance between clusters\n\n')
    fout <- file('{output.dist}', 'w')
    close(fout)
    quit(status= 0)
}}

dists <- fread(cmd= 'clm dist --sort --index {input.clst}', header= FALSE, col.names= c('rand_index', 'jaccard_index', 'adjusted_rand_index', 'n1', 'n2'))
stopifnot(grepl('^rand=', dists$rand_index))
stopifnot(grepl('^jaccard=', dists$jaccard_index))
stopifnot(grepl('^arand=', dists$adjusted_rand_index))

dists[, rand_index := sub('^rand=', '', rand_index)]
dists[, jaccard_index := sub('^jaccard=', '', jaccard_index)]
dists[, adjusted_rand_index := sub('^arand=', '', adjusted_rand_index)]
dists[, n1 := sub('^n1=', '', n1)]
dists[, n2 := sub('^n2=', '', n2)]

dists[, corr_n1 := as.numeric(sub('/.*', '', sub('.*correlation_', '', n1)))]
dists[, inflation_n1 := as.numeric(sub('/.*', '', sub('.*inflation_', '', n1)))]
dists[, ceilnb_n1 := as.numeric(sub('/.*', '', sub('.*ceilnb_', '', n1)))]

dists[, corr_n2 := as.numeric(sub('/.*', '', sub('.*correlation_', '', n2)))]
dists[, inflation_n2 := as.numeric(sub('/.*', '', sub('.*inflation_', '', n2)))]
dists[, ceilnb_n2 := as.numeric(sub('/.*', '', sub('.*ceilnb_', '', n2)))]

setcolorder(dists, c('corr_n1', 'corr_n2', 'inflation_n1', 'inflation_n2', 'ceilnb_n1', 'ceilnb_n2'))

write.table(dists, '{output.dist}', row.names= FALSE, sep= '\t', quote= FALSE)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule dump_clusters:
    input:
        clst= os.path.join(config['mcl_output_dir'], 'correlation_{r}/inflation_{inflation}/ceilnb_{ceilnb}/cluster.mci'),
        gene_dict= os.path.join(config['mcl_output_dir'], 'gene_dict.txt'),
    output:
        clst= os.path.join(config['mcl_output_dir'], 'correlation_{r}/inflation_{inflation}/ceilnb_{ceilnb}/cluster.tsv'),
    shell:
        r"""
        mcxdump -icl {input.clst} --dump-pairs -o {output.clst} -tabr {input.gene_dict}
        """

rule summarise_clusters:
    input:
        clst= expand(os.path.join(config['mcl_output_dir'], 'correlation_{r}/inflation_{inflation}/ceilnb_{ceilnb}/cluster.tsv'), 
                r= config['correlation'], inflation= config['inflation'], ceilnb= config['ceilnb']),
    output:
        smry= os.path.join(config['mcl_output_dir'], 'cluster_summary.tsv'),
    params:
        regular_size= config['regular_size'],
        small_size= config['small_size'],
        top_big= config['top_big'],
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R

library(data.table)

REGULAR_SIZE <- as.numeric(strsplit('{params.regular_size}', ' ')[[1]])
stopifnot(length(REGULAR_SIZE) == 2)
SMALL_SIZE <- as.numeric('{params.small_size}')
TOP_BIG <- as.numeric('{params.top_big}')

clst_files <- strsplit('{input.clst}', ' ')[[1]]

clst <- list()
for(fn in clst_files) {{
    dt <- fread(fn, header= FALSE, col.names= c('cluster_id', 'gene_id'))
    dt[, filename := fn]
    clst[[length(clst) + 1]] <- dt
}}
clst <- rbindlist(clst)

clst_size <- clst[, list(size= .N), by= list(cluster_id, filename)][order(filename, size)]

num_clst <- clst_size[, list(n_clusters= .N, median_size= median(size)), by= filename]
regular <- clst_size[size >= REGULAR_SIZE[1] & size <= REGULAR_SIZE[2]][, list(ngenes_regular= sum(size)), by= filename]
small <- clst_size[size <= SMALL_SIZE][, list(ngenes_small_k= sum(size)), by= filename]
big <- clst_size[cluster_id <= TOP_BIG, list(ngenes_big_k= sum(size)), by= filename]

smry <- merge(num_clst, regular, all.x= TRUE, by= 'filename')
smry <- merge(smry, small, all.x= TRUE, by= 'filename')
smry <- merge(smry, big, all.x= TRUE, by= 'filename')
smry[is.na(smry)] <- 0

n_genes <- length(unique(clst$gene_id))
smry[, pct_genes_big_k := 100 * ngenes_big_k/n_genes]
smry[, pct_genes_small_k := 100 * ngenes_small_k/n_genes]
smry[, pct_genes_regular := 100 * ngenes_regular/n_genes]
smry[, corr := as.numeric(sub('/.*', '', sub('.*correlation_', '', filename)))]
smry[, inflation := as.numeric(sub('/.*', '', sub('.*inflation_', '', filename)))]
smry[, ceilnb := as.numeric(sub('/.*', '', sub('.*ceilnb_', '', filename)))]
smry <- smry[order(-pct_genes_regular)]

setcolorder(smry, c('corr', 'inflation', 'ceilnb', 'n_clusters', 'pct_genes_regular', 'pct_genes_big_k', 'pct_genes_small_k'))

smry[, pct_genes_big_k := sprintf('%.2f', pct_genes_big_k)]
smry[, pct_genes_small_k := sprintf('%.2f', pct_genes_small_k)]
smry[, pct_genes_regular := sprintf('%.2f', pct_genes_regular)]

stopifnot(all(clst_files %in% smry$filename))
stopifnot(length(unique(smry$filename)) == nrow(smry))

write.table(smry, '{output.smry}', row.names= FALSE, quote= FALSE, sep= '\t')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """
