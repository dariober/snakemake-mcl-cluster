config['data_matrix'] = 'data/data_matrix.tsv'
config['mcl_output_dir'] = 'test_out/mcl'
config['inflation'] = [1.4, 2]

rule all:
    input:
        'test_out/mcl/cluster_summary.tsv',

include: '../workflows/mcl.smk'
