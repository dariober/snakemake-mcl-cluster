config['inflation'] = [1.4, 2.0]

samples = ['sample_1', 'sample_2']

rule all:
    input:
        expand('test_out/{sample}/cluster_summary.tsv', sample= samples),

rule make_data_matrix:
    input:
        mat= 'data/data_matrix.tsv',
    output:
        mat= 'test_out/{sample}/data_matrix.tsv',
    shell:
        r"""
        cp {input.mat} {output.mat}
        """

include: '../workflows/mcl.smk'
