defaults = {'mcl_output_dir': 'mcl', 'pearson_r_cutoff': [0.6, 0.8], 'inflation': [1.4, 2], 'ceilnb': [10000], 'n_threads': 1, 'regular_size': [15, 200], 'small_size': 3, 'top_big': 3}
help_doc = r"""
USAGE
snakemake -C [option=value] [option=value] ...

OPTIONS and DEFAULTS

==== Input/Output ====

data_matrix
    Input data matrix. First row is header and first column is feature
    identifier. For example, a data matrix to cluster genes has one row per
    gene and one column per sample. Required.

mcl_output_dir={mcl_output_dir} 
    Write output in this directory. Created if it does not exists. 

====  Clustering ====

pearson_r_cutoff='{pearson_r_cutoff}'
    List correlation cutoffs. Reset to 0 correlations below this cutoff

inflation='{inflation}'
    List of inflations for MCL clustering

ceilnb='{ceilnb}'
    List of 'ceilnb' value. Nodes will be pruned if they exceed this many
    connections

==== Summary and performance ====

regular_size='{regular_size}'
    List of length 2. Consider as "regular" clusters within this size interval

small_size='{small_size}'
    Consider "small" cluster of this or lower size. Useful to assess whether
    the clustering is over-fragmented

top_big='{top_big}'
    Number of top/biggest clusters. Useful to assess whether few big clusters
    attract most of the features

n_threads=1
    Number of threads for steps allowing multithreading

help=1
    Show this help

EXAMPLES

Cluster testing all combination of 2 cutoffs for correlation, 2 of ceilnb, and
1 inflation:

    snakemake --jobs 5 -C data_matrix=data/logrpkm.tsv \
                          pearson_r_cutoff='[0.6, 0.8]' \
                          ceilnb='[100, 200]' \
                          inflation='[1.4]'

Show this help:

    snakemake --jobs 5 -C help=1
""".format(**defaults)

if 'help' in config:
    print(help_doc)
    sys.exit(0)

if 'data_matrix' not in config:
    raise OSError('''\nPlease specify the input data matrix for MCL clustering. E.g:\n
    snakemake -C data_matrix=my-matrix.txt ...\n''')
if 'mcl_output_dir' not in config:
    config['mcl_output_dir'] = defaults['mcl_output_dir']
if 'pearson_r_cutoff' not in config:
    config['pearson_r_cutoff'] = defaults['pearson_r_cutoff']
if 'inflation' not in config:
    config['inflation'] = defaults['inflation']
if 'ceilnb' not in config:
    config['ceilnb'] = defaults['ceilnb']
if 'regular_size' not in config:
    config['regular_size'] = defaults['regular_size']
if 'small_size' not in config:
    config['small_size'] = defaults['small_size']
if 'top_big' not in config:
    config['top_big'] = defaults['top_big']
if 'n_threads' not in config:
    config['n_threads'] = defaults['n_threads']
 
mcl_output_dir = config['mcl_output_dir']

rule mcl_all:
    input:
        f'{mcl_output_dir}/similarity_matrix.tsv.gz',
        f'{mcl_output_dir}/vary_correlation_query.txt',
        expand(os.path.join(mcl_output_dir, 'Pearson_{pearson_r}/I_{inflation}/ceilnb_{ceilnb}/cluster.tsv'),
                pearson_r= config['pearson_r_cutoff'], inflation= config['inflation'], ceilnb= config['ceilnb']),
        f'{mcl_output_dir}/cluster_summary.tsv',
        f'{mcl_output_dir}/distance_between_clusters.tsv',
        
rule make_mcl_network:
    input:
        mat= config['data_matrix'],
    output:
        network= temp(f'{mcl_output_dir}/network.mci'),
        gene_dict= temp(f'{mcl_output_dir}/gene_dict.txt'),
    shell:
        r"""
        mcxarray -data {input.mat} -skipr 1 -skipc 1 -o {output.network} \
            --write-binary --pearson -co 0.2 -tf 'abs()' -write-tab {output.gene_dict}
        """

rule write_similarity_matrix:
    input:
        network= f'{mcl_output_dir}/network.mci',
        gene_dict= f'{mcl_output_dir}/gene_dict.txt',
    output:
        matrix= f'{mcl_output_dir}/similarity_matrix.tsv.gz',
    shell:
        r"""
        mcxdump -imx {input.network} -tab {input.gene_dict} --dump-upper | gzip > {output.matrix}
        """

rule query_mcl_network:
    input:
        network= f'{mcl_output_dir}/network.mci',
    output:
        qry= f'{mcl_output_dir}/vary_correlation_query.txt',
    shell:
        r"""
        mcx query -imx {input.network} --vary-correlation -o {output.qry}
        """

rule select_mcl_network:
    input:
        network= f'{mcl_output_dir}/network.mci',
    output:
        selected_ntwk= temp(os.path.join(mcl_output_dir, 'Pearson_{pearson_r}/selected_network.mci')),
    shell:
        r"""
        mcx alter -imx {input.network} -tf 'gq({wildcards.pearson_r}), add(-{wildcards.pearson_r})' \
            --write-binary -o {output.selected_ntwk}
        """

rule cluster_network:
    input:
        selected_ntwk= os.path.join(mcl_output_dir, 'Pearson_{pearson_r}/selected_network.mci'),
    output:
        clst= temp(os.path.join(mcl_output_dir, 'Pearson_{pearson_r}/I_{inflation}/ceilnb_{ceilnb}/cluster.mci')),
    params:
        nt= config['n_threads']
    shell:
        r"""
        mcl {input.selected_ntwk} -o {output.clst} -I {wildcards.inflation} -te {params.nt} -tf '#ceilnb({wildcards.ceilnb})'
        """

rule distance_between_clusters:
    input:
        clst= expand(os.path.join(mcl_output_dir, 'Pearson_{pearson_r}/I_{inflation}/ceilnb_{ceilnb}/cluster.mci'), 
                pearson_r= config['pearson_r_cutoff'], inflation= config['inflation'], ceilnb= config['ceilnb']),
    output:
        dist= os.path.join(mcl_output_dir, 'distance_between_clusters.tsv'),
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
        clst= os.path.join(mcl_output_dir, 'Pearson_{pearson_r}/I_{inflation}/ceilnb_{ceilnb}/cluster.mci'),
        gene_dict= os.path.join(mcl_output_dir, 'gene_dict.txt'),
    output:
        clst= os.path.join(mcl_output_dir, 'Pearson_{pearson_r}/I_{inflation}/ceilnb_{ceilnb}/cluster.tsv'),
    shell:
        r"""
        mcxdump -icl {input.clst} --dump-pairs -o {output.clst} -tabr {input.gene_dict}
        """

rule summarise_clusters:
    input:
        clst= expand(os.path.join(mcl_output_dir, 'Pearson_{pearson_r}/I_{inflation}/ceilnb_{ceilnb}/cluster.tsv'), 
                pearson_r= config['pearson_r_cutoff'], inflation= config['inflation'], ceilnb= config['ceilnb']),
    output:
        smry= os.path.join(mcl_output_dir, 'cluster_summary.tsv'),
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
