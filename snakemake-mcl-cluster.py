#!/usr/bin/env python3

import subprocess
import argparse
import os
import sys
import urllib.request
import filecmp
import shutil

GITHUB_REPOSITORY = "https://raw.githubusercontent.com/dariober/snakemake-mcl-cluster/51167fe/workflows"

def get_workflowdir(workflow_dir, download_to= None):
    """Return path workflow directory with snakemake files or download them and
    return path to the receiving directory
    """
    expected_files = ['Snakefile', 'mcl.smk'] 

    if workflow_dir.startswith('https://'):
        assert download_to is not None, 'Need a path to a download dir'
        os.makedirs(os.path.join(download_to), exist_ok= True)

        for smk in expected_files:
            url = os.path.join(workflow_dir, smk)
            dst = os.path.join(download_to, smk)
            try:
                tmp = urllib.request.urlretrieve(url, filename= dst)
            except Exception as e:
                sys.stderr.write('Downloading %s to %s\n' % (url, dst))
                sys.stderr.write(str(e) + '\n')
                sys.exit(1)
        workflow_dir = download_to
    
    workflow_dir = os.path.abspath(workflow_dir)

    if not os.path.isdir(workflow_dir):
        sys.stderr.write('\n"%s" is not a directory\n\n' % workflow_dir)
        sys.exit(1)

    for smk in expected_files:
        if not os.path.isfile(os.path.join(workflow_dir, smk)):
            sys.stderr.write('\nWorkflow dir "%s" does not contain %s\n\n' % (workflow_dir, smk))
            sys.exit(1)
    
    return workflow_dir

def copy_data_matrix(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok= True)
    if not os.path.exists(dst) or not filecmp.cmp(src, dst):
        shutil.copyfile(src, dst)

parser = argparse.ArgumentParser(description= """
DESCRIPTION

Run mcl clustering on array of parameters and summarise the clustering
results
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

io_args = parser.add_argument_group('Input/Output', '')
params_args = parser.add_argument_group('Clustering parameters', '')
smry_args = parser.add_argument_group('Only for summary', '')

io_args.add_argument('--data-matrix', '-i',
    required= True,
    help='''Input data matrix with features to cluster as rows
and observations as columns. First column is feature
identifier and first row is header''')

io_args.add_argument('--output-dir', '-o',
    required= True,
    help='Output directory')

io_args.add_argument('--workflow-dir', '-d',
    default= GITHUB_REPOSITORY,
    help= r"""Directory of snakemake workflows. Default to GitHub link
{%(default)s}""")

io_args.add_argument('--jobs', '-j',
    type= int,
    default= 1,
    help='Number of jobs to run in parallel {%(default)s}')

io_args.add_argument('--n-threads', '-t',
    type= int,
    default= 1,
    help='Number of threads for mcl clustering {%(default)s}')

io_args.add_argument('--snakemake', '-s',
    default= '',
    help='String of further options to pass to snakemake. E.g. "-n -p" {%(default)s}')

params_args.add_argument('--edge-weight', '-w',
    default= 'pearson',
    choices= ['pearson', 'spearman', 'cosine'],
    help='Compute edge weights using this metric {%(default)s}')

params_args.add_argument('--correlation', '-r',
    default= [0.6],
    nargs= '+',
    type= float,
    help="""Reset to 0 weights (e.g. Pearson correlations) below this cutoff.
Make sure this cutoff is suitable for the chosen edge weight metric
{%(default)s}""")

params_args.add_argument('--inflation', '-I',
    default= [1.4],
    nargs= '+',
    type= float,
    help='Inflation for mcl algorithm {%(default)s}')

params_args.add_argument('--ceilnb', '-n',
    default= [10000],
    nargs= '+',
    type= int,
    help='Nodes with more than this many connections will be pruned {%(default)s}')

smry_args.add_argument('--regular-size', '-R',
    default= [15, 200],
    type= int,
    nargs= 2,
    help='Clusters of size within this range are "regular" {%(default)s}')

smry_args.add_argument('--small-size', '-S',
    default= 3,
    type= int,
    help='Clusters below this are "small" {%(default)s}')

smry_args.add_argument('--top-big', '-b',
    default= 3,
    type= int,
    help='Number of top, biggest clusters {%(default)s}')

parser.add_argument('--quiet', '-q', action='store_true', help= "Print progress messages only in case of failure")
parser.add_argument('--version', '-v', action='version', version='%(prog)s 0.1.0')

if __name__ == "__main__":

    args = parser.parse_args()

    workflow_dir = get_workflowdir(args.workflow_dir, download_to= os.path.join(args.output_dir, 'workflows'))

    copy_data_matrix(args.data_matrix, os.path.join(args.output_dir, 'mcl/data_matrix.tsv'))

    config = []
    config.append("edge_weight=%s" % args.edge_weight)
    config.append("correlation='[%s]'" % ','.join([str(x) for x in args.correlation]))
    config.append("inflation='[%s]'" % ','.join([str(x) for x in args.inflation]))
    config.append("ceilnb='[%s]'" % ','.join([str(x) for x in args.ceilnb]))
    config.append("regular_size='[%s]'" % ','.join([str(x) for x in args.regular_size]))
    config.append("small_size=%s" % args.small_size)
    config.append("top_big=%s" % args.top_big)
    config.append("n_threads=%s" % args.n_threads)
    config.append("workflow_dir='%s'" % workflow_dir)

    config = ' '.join(config)

    cmd = r"""snakemake --rerun-incomplete --jobs {jobs} \
        --directory {outdir} \
        --config {config} \
        --snakefile '{snakefile}' {opts}
    """.format(jobs= args.jobs, 
            outdir= args.output_dir, 
            config= config, 
            snakefile= os.path.join(workflow_dir, 'Snakefile'), 
            opts= args.snakemake)
    
    log = []
    if not args.quiet:
        print(cmd)
    log.append(cmd)

    p = subprocess.Popen(cmd, shell= True, executable= 'bash', stdout= subprocess.PIPE, stderr= subprocess.STDOUT)
    while p.poll() is None:
        l = p.stdout.readline().decode() # This blocks until it receives a newline.
        log.append(l)
        if not args.quiet:
            sys.stdout.write(l)
    if p.returncode != 0:
        sys.stderr.write(''.join(log) + '\n')
        sys.exit(p.returncode)
