#!/usr/bin/env python3

import unittest
import shutil
import subprocess as sp
import os
import sys
import filecmp

class TestMcl(unittest.TestCase):

    def setUp(self):
        sys.stderr.write('\n' + self.id().split('.')[-1] + ' ') # Print test name
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')

    def tearDown(self):
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')

    def testShowHelp(self):
        p = sp.Popen("../snakemake-mcl-cluster.py --help", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)

    def testDefaults(self):
        p = sp.Popen("../snakemake-mcl-cluster.py -i data/data_matrix.tsv -d ../workflows -o test_out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.isfile('test_out/mcl/cluster_summary.tsv'))
        self.assertTrue(os.path.getsize('test_out/mcl/cluster_summary.tsv') > 100)

    def testPullScriptsFromGitHub(self):
        p = sp.Popen("../snakemake-mcl-cluster.py -i data/data_matrix.tsv -o test_out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.isfile('test_out/mcl/cluster_summary.tsv'))
        self.assertTrue(os.path.getsize('test_out/mcl/cluster_summary.tsv') > 100)

    def testArrayOfParameters(self):
        p = sp.Popen("../snakemake-mcl-cluster.py -j 5 -i data/data_matrix.tsv -d ../workflows -o test_out -I 1.5 2 -r 0.6 0.8 -n 50 100 1000", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.isfile('test_out/mcl/correlation_0.6/inflation_1.5/ceilnb_50/cluster.tsv'))
        self.assertTrue(os.path.isfile('test_out/mcl/correlation_0.6/inflation_1.5/ceilnb_100/cluster.tsv'))
        self.assertTrue(os.path.isfile('test_out/mcl/correlation_0.6/inflation_1.5/ceilnb_1000/cluster.tsv'))
        self.assertTrue(os.path.isfile('test_out/mcl/correlation_0.6/inflation_2.0/ceilnb_50/cluster.tsv'))
        self.assertTrue(os.path.isfile('test_out/mcl/correlation_0.6/inflation_2.0/ceilnb_100/cluster.tsv'))
        self.assertTrue(os.path.isfile('test_out/mcl/correlation_0.6/inflation_2.0/ceilnb_1000/cluster.tsv'))
        self.assertTrue(os.path.isfile('test_out/mcl/correlation_0.8/inflation_1.5/ceilnb_50/cluster.tsv'))
        self.assertTrue(os.path.isfile('test_out/mcl/correlation_0.8/inflation_1.5/ceilnb_100/cluster.tsv'))
        self.assertTrue(os.path.isfile('test_out/mcl/correlation_0.8/inflation_1.5/ceilnb_1000/cluster.tsv'))
        self.assertTrue(os.path.isfile('test_out/mcl/correlation_0.8/inflation_2.0/ceilnb_50/cluster.tsv'))
        self.assertTrue(os.path.isfile('test_out/mcl/correlation_0.8/inflation_2.0/ceilnb_100/cluster.tsv'))
        self.assertTrue(os.path.isfile('test_out/mcl/correlation_0.8/inflation_2.0/ceilnb_1000/cluster.tsv'))
        self.assertTrue(os.path.isfile('test_out/mcl/cluster_summary.tsv'))
        self.assertTrue(os.path.getsize('test_out/mcl/cluster_summary.tsv') > 100)

    def testNothingToBeDone(self):
        p = sp.Popen("../snakemake-mcl-cluster.py -i data/data_matrix.tsv -d ../workflows -o test_out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)

        p = sp.Popen("../snakemake-mcl-cluster.py -i data/data_matrix.tsv -d ../workflows -o test_out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('Nothing to be done' in stdout.decode())

    def testDoNotRecreateExistingFiles(self):
        p = sp.Popen("../snakemake-mcl-cluster.py -i data/data_matrix.tsv -d ../workflows -o test_out -I 0.8", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        mt = os.stat('test_out/mcl/correlation_0.6/inflation_0.8/ceilnb_10000/cluster.tsv').st_mtime

        p = sp.Popen("../snakemake-mcl-cluster.py -i data/data_matrix.tsv -d ../workflows -o test_out -I 0.6 0.8", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.exists('test_out/mcl/correlation_0.6/inflation_0.6/ceilnb_10000/cluster.tsv'))

        mt2 = os.stat('test_out/mcl/correlation_0.6/inflation_0.8/ceilnb_10000/cluster.tsv').st_mtime
        self.assertEqual(mt, mt2)

    def testRerunIfInputChanged(self):
        p = sp.Popen("../snakemake-mcl-cluster.py -i data/data_matrix.tsv -d ../workflows -o test_out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        smry = os.stat('test_out/mcl/cluster_summary.tsv').st_mtime
        mci = os.stat('test_out/mcl/network.mci').st_mtime

        p = sp.Popen("../snakemake-mcl-cluster.py -i data/data_matrix2.tsv -d ../workflows -o test_out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        smry2 = os.stat('test_out/mcl/cluster_summary.tsv').st_mtime
        mci2 = os.stat('test_out/mcl/network.mci').st_mtime

        self.assertNotEqual(smry, smry2)
        self.assertNotEqual(mci, mci2)

    def testDoNothingOnDifferentInputNameWithSameContent(self):
        """This behaviour is kind of undesirable. If two input matrices have
        same content but different names the analysis is not repeated
        """
        assert filecmp.cmp('data/data_matrix.tsv', 'data/data_matrix_copy.tsv')

        p = sp.Popen("../snakemake-mcl-cluster.py -i data/data_matrix.tsv -d ../workflows -o test_out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)

        p = sp.Popen("../snakemake-mcl-cluster.py -i data/data_matrix_copy.tsv -d ../workflows -o test_out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('Nothing to be done' in stdout.decode())

    def testIncludeInExternalSnakemakePipeline(self):
        p = sp.Popen("snakemake -s snakefile_include.smk -j 1", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.isfile('test_out/sample_1/cluster_summary.tsv'))
        self.assertTrue(os.path.getsize('test_out/sample_2/cluster_summary.tsv') > 100)

    def testRunInShellInExternalSnakemakePipeline(self):
        p = sp.Popen("snakemake -s snakefile_shell.smk -j 1", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.isfile('test_out/mcl/cluster_summary.tsv'))
        self.assertTrue(os.path.getsize('test_out/mcl/cluster_summary.tsv') > 100)

if __name__ == '__main__':
    unittest.main()
