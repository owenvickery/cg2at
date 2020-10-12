#!/usr/bin/env python3
import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/..')
import unittest
from unittest.mock import patch
import mock
import numpy as np
import gen, gro, at_mod, at_mod_p, at_mod_np, read_in, g_var




run_dir = os.path.dirname(os.path.realpath(__file__))+'/'
class TestSum(unittest.TestCase):
    def test_add_residue_to_dictionary(self):
        g_var.p_residues, g_var.np_residues = ['PHE'], ['CHOL', 'NA', 'W']
        g_var.cg_water_types = ['W', 'SOL', 'WN', 'WF', 'PW']
        line_test = [{'residue_name':'CHOL'}, {'residue_name':'PHE'}, {'residue_name':'NA'}, {'residue_name':'W'}]
        for l_val, l in enumerate(line_test):
            read_in.add_residue_to_dictionary(l)
        self.assertEqual(g_var.cg_residues, {'CHOL': {}, 'PROTEIN': {}, 'NA': {}, 'SOL': {}})

    def test_check_new_box(self):
        box = [100,100,100] 
        new_box = [90,90,90]
        g_var.args.box = [90,90,90]
        self.assertEqual(read_in.check_new_box([50,50, 50], box, new_box), False)
        self.assertEqual(read_in.check_new_box([95,95,95], box, new_box), True)

    def test_swap(self):
        swap=[{'POPE': {'POPE:POPG': {'NH3': 'GL0', 'resid': 'ALL', 'range': 'ALL'}}},\
              {'POPE': {'POPE:POPG': {'NH3': 'GL0', 'resid': 'ALL', 'range': 'ALL'}}},\
              {'POPE': {'POPE:POPG': {'NH3': 'skip', 'resid': 'ALL', 'range': 'ALL'}}},\
              {'NA+': {'NA+:skip': {'ALL': 'ALL', 'resid': [4000, 4001, 4002], 'range': ['4000-4002']}}},\
              {'NA+': {'NA+:skip': {'ALL': 'ALL', 'resid': [4000, 4001, 4002], 'range': ['4000-4002']}}}]
        test_suite = [['NH3', 'POPE', 50],['C1A', 'POPE', 50],['NH3', 'POPE', 50],['NA+', 'NA+', 4000],['NA+', 'NA+',30]  ]
        test_out = [['GL0', 'POPG'],['C1A', 'POPG'],['SKIP', 'POPG'],['NA+', 'SKIP'],['NA+', 'NA+'] ]
        for test_val, test in enumerate(test_suite):
            g_var.swap_dict = swap[test_val]
            atom, resname = read_in.swap(test[0], test[1], test[2])
            self.assertEqual(atom, test_out[test_val][0])
            self.assertEqual(resname, test_out[test_val][1])



    def test_add_to_cg_database(self):
        pass

    def test_read_initial_cg_pdb(self):
        pass

if __name__ == '__main__':
    unittest.main()
