rule all:
    input:
        'test_out/mcl_shell/cluster_summary.tsv',

rule mcl:
    input:
        mat= 'data/data_matrix.tsv',
    output:
        'test_out/mcl_shell/cluster_summary.tsv',
    shell:
        r"""
        ../snakemake-mcl-cluster.py -m data/data_matrix.tsv -o test_out/mcl_shell -d ../workflows
        """

