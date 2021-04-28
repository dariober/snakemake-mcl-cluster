#!/usr/bin/env python3

import unittest
import shutil
import subprocess as sp
import os

class cnv_facets(unittest.TestCase):

    def tearDown(self):
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')

    def testShowHelp(self):
        p = sp.Popen("../snakemake-mcl-cluster.py --help", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)

    def testDefaults(self):
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')
        p = sp.Popen("../snakemake-mcl-cluster.py -m data/data_matrix.tsv -d ../workflows -o test_out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.isfile('test_out/cluster_summary.tsv'))
        self.assertTrue(os.path.getsize('test_out/cluster_summary.tsv') > 100)

    def testParamArray(self):
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')
        p = sp.Popen("../snakemake-mcl-cluster.py --jobs 5 -m data/data_matrix.tsv -d ../workflows -I 1.5 2 -r 0.6 0.8 -n 50 100 1000 -o test_out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.isfile('test_out/cluster_summary.tsv'))
        self.assertTrue(os.path.getsize('test_out/cluster_summary.tsv') > 100)

if __name__ == '__main__':
    unittest.main()
