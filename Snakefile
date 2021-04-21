pearson_r_cutoff = config['pearson_r_cutoff']
inflation = config['inflation']
ceilnb = config['ceilnb']

rule all:
    input:
        'mcl/similarity_matrix.tsv.gz',
        'mcl/vary_correlation_query.txt',
        expand('mcl/Pearson_{pearson_r}/I_{inflation}/ceilnb_{ceilnb}/cluster.tsv', 
                pearson_r= pearson_r_cutoff, inflation= inflation, ceilnb= ceilnb),
        'mcl/cluster_summary.tsv',
        'mcl/distance_between_clusters.tsv',
        

rule make_mcl_network:
    input:
        mat= config['data_matrix'],
    output:
        network= temp('mcl/network.mci'),
        gene_dict= temp('mcl/gene_dict.txt'),
    shell:
        r"""
        mcxarray -data {input.mat} -skipr 1 -skipc 1 -o {output.network} \
            --write-binary --pearson -co 0.2 -tf 'abs()' -write-tab {output.gene_dict}
        """

rule write_similarity_matrix:
    input:
        network= 'mcl/network.mci',
        gene_dict= 'mcl/gene_dict.txt',
    output:
        matrix= 'mcl/similarity_matrix.tsv.gz',
    shell:
        r"""
        mcxdump -imx {input.network} -tab {input.gene_dict} --dump-upper | gzip > {output.matrix}
        """

rule query_mcl_network:
    input:
        network= 'mcl/network.mci',
    output:
        qry= 'mcl/vary_correlation_query.txt',
    shell:
        r"""
        mcx query -imx {input.network} --vary-correlation -o {output.qry}
        """

rule select_mcl_network:
    input:
        network= 'mcl/network.mci',
    output:
        selected_ntwk= temp('mcl/Pearson_{pearson_r}/selected_network.mci'),
    shell:
        r"""
        mcx alter -imx {input.network} -tf 'gq({wildcards.pearson_r}), add(-{wildcards.pearson_r})' \
            --write-binary -o {output.selected_ntwk}
        """

rule cluster_network:
    input:
        selected_ntwk= 'mcl/Pearson_{pearson_r}/selected_network.mci',
    output:
        clst= temp('mcl/Pearson_{pearson_r}/I_{inflation}/ceilnb_{ceilnb}/cluster.mci'),
    params:
        nt= config['n_threads']
    shell:
        r"""
        mcl {input.selected_ntwk} -o {output.clst} -I {wildcards.inflation} -te {params.nt} -tf '#ceilnb({wildcards.ceilnb})'
        """

rule distance_between_clusters:
    input:
        clst= expand('mcl/Pearson_{pearson_r}/I_{inflation}/ceilnb_{ceilnb}/cluster.mci', 
                pearson_r= pearson_r_cutoff, inflation= inflation, ceilnb= ceilnb),
    output:
        dist= 'mcl/distance_between_clusters.tsv',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R

library(data.table)

dists <- fread(cmd= 'clm dist --sort --index {input.clst}', header= FALSE, col.names= c('rand_index', 'jaccard_index', 'adjusted_rand_index', 'n1', 'n2'))
stopifnot(grepl('^rand=', dists$rand_index))
stopifnot(grepl('^jaccard=', dists$jaccard_index))
stopifnot(grepl('^arand=', dists$adjusted_rand_index))

dists[, rand_index := sub('^rand=', '', rand_index)]
dists[, jaccard_index := sub('^jaccard=', '', jaccard_index)]
dists[, adjusted_rand_index := sub('^arand=', '', adjusted_rand_index)]
dists[, n1 := sub('^n1=', '', n1)]
dists[, n2 := sub('^n2=', '', n2)]

dists[, corr_n1 := as.numeric(sub('/.*', '', sub('.*Pearson_', '', n1)))]
dists[, inflation_n1 := as.numeric(sub('/.*', '', sub('.*I_', '', n1)))]
dists[, ceilnb_n1 := as.numeric(sub('/.*', '', sub('.*ceilnb_', '', n1)))]

dists[, corr_n2 := as.numeric(sub('/.*', '', sub('.*Pearson_', '', n2)))]
dists[, inflation_n2 := as.numeric(sub('/.*', '', sub('.*I_', '', n2)))]
dists[, ceilnb_n2 := as.numeric(sub('/.*', '', sub('.*ceilnb_', '', n2)))]

setcolorder(dists, c('corr_n1', 'corr_n2', 'inflation_n1', 'inflation_n2', 'ceilnb_n1', 'ceilnb_n2'))

write.table(dists, '{output.dist}', row.names= FALSE, sep= '\t', quote= FALSE)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule dump_clusters:
    input:
        clst= 'mcl/Pearson_{pearson_r}/I_{inflation}/ceilnb_{ceilnb}/cluster.mci',
        gene_dict= 'mcl/gene_dict.txt',
    output:
        clst= 'mcl/Pearson_{pearson_r}/I_{inflation}/ceilnb_{ceilnb}/cluster.tsv',
    shell:
        r"""
        mcxdump -icl {input.clst} --dump-pairs -o {output.clst} -tabr {input.gene_dict}
        """

rule summarise_clusters:
    input:
        clst= expand('mcl/Pearson_{pearson_r}/I_{inflation}/ceilnb_{ceilnb}/cluster.tsv', 
                pearson_r= pearson_r_cutoff, inflation= inflation, ceilnb= ceilnb),
    output:
        smry= 'mcl/cluster_summary.tsv',
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

regular <- clst_size[size >= REGULAR_SIZE[1] & size <= REGULAR_SIZE[2]][, list(ngenes_regular= sum(size)), by= filename]
small <- clst_size[size <= SMALL_SIZE][, list(ngenes_smallk= sum(size)), by= filename]
big <- clst_size[cluster_id <= TOP_BIG, list(ngenes_bigk= sum(size)), by= filename]
num_clst <- clst_size[, list(n_clusters= .N, median_size= median(size)), by= filename]

smry <- merge(regular, merge(
                        num_clst, merge(
                            big, small, by= 'filename'), by= 'filename'), by= 'filename')

n_genes <- length(unique(clst$gene_id))
smry[, pct_genes_bigk := 100 * ngenes_bigk/n_genes]
smry[, pct_genes_smallk := 100 * ngenes_smallk/n_genes]
smry[, pct_genes_regular := 100 * ngenes_regular/n_genes]
smry[, corr := as.numeric(sub('/.*', '', sub('.*Pearson_', '', filename)))]
smry[, inflation := as.numeric(sub('/.*', '', sub('.*I_', '', filename)))]
smry[, ceilnb := as.numeric(sub('/.*', '', sub('.*ceilnb_', '', filename)))]
smry <- smry[order(-pct_genes_regular)]

setcolorder(smry, c('corr', 'inflation', 'ceilnb', 'n_clusters', 'pct_genes_regular', 'pct_genes_bigk', 'pct_genes_smallk'))

write.table(smry, '{output.smry}', row.names= FALSE, quote= FALSE, sep= '\t')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """
