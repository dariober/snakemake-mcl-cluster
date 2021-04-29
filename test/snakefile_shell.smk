rule all:
    input:
        'test_out/mcl/cluster_summary.tsv',

rule mcl:
    input:
        mat= 'data/data_matrix.tsv',
    output:
        'test_out/mcl/cluster_summary.tsv',
    shell:
        r"""
        ../snakemake-mcl-cluster.py -i data/data_matrix.tsv -o test_out -d ../workflows
        """

