#!/usr/bin/env python3
import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/..')
import unittest
import numpy as np
import gen, gro, at_mod, at_mod_p, at_mod_np, read_in, g_var

run_dir = os.path.dirname(os.path.realpath(__file__))+'/'
class TestSum(unittest.TestCase):
    def test_time(self):
        t1 = 1602152289.7587886
        t2 = 1602154289.7587886
        result = gen.fix_time(t2, t1)
        self.assertEqual(result, ' 0 hours 33 min 20 sec ' )

    def test_dist(self):
        p1 = [58.274,66.912,14.038]
        p2 = [59.360,66.216,14.859]
        result = gen.calculate_distance(p1, p2)
        self.assertEqual(result, 1.5290039241283893)

    def test_atoms(self):
        atoms = ['H1', '1H', 'HA', 'CA', '1CA', 'CAH' ]   
        result = []
        for atom_name in atoms:
            result.append(gen.is_hydrogen(atom_name))
        self.assertEqual(result, [True, True, True, False, False, False])

    def test_split_swap(self):
        swap = ['POPE,NH3:POPG,GL0', 'POPG:POPE', 'NA+:skip:4000-4002','POPG:skip','GLU,SC2:ASP,skip']
        swap_out = [['ALL', 'ALL'], ['ALL', 'ALL'], [['4000-4002'], [4000, 4001, 4002]], ['ALL', 'ALL'], ['ALL', 'ALL']] 
        for swap_val, swap_type in enumerate(swap):
            result1, result2 = gen.split_swap(swap_type)
            self.assertEqual(result1,swap_out[swap_val][0])
            self.assertEqual(result2,swap_out[swap_val][1])

    def test_sort_swap_group(self):
        out=[{'POPE': {'POPE:POPG': {'NH3': 'GL0', 'resid': 'ALL', 'range': 'ALL'}}},\
             {'POPG': {'POPG:POPE': {'ALL': 'ALL', 'resid': 'ALL', 'range': 'ALL'}}},\
             {'NA+': {'NA+:skip': {'ALL': 'ALL', 'resid': [4000, 4001, 4002], 'range': ['4000-4002']}}}, \
             {'POPG': {'POPG:skip': {'ALL': 'ALL', 'resid': 'ALL', 'range': 'ALL'}}}, \
             {'GLU': {'GLU:ASP': {'SC2': 'skip', 'resid': 'ALL', 'range': 'ALL'}}}]
        for swap_val, swap_types in enumerate([['POPE,NH3:POPG,GL0'], ['POPG:POPE'], ['NA+:skip:4000-4002'], ['POPG:skip'], ['GLU,SC2:ASP,skip']]):
            g_var.args.swap = swap_types
            gen.sort_swap_group()
            self.assertEqual(g_var.swap_dict, out[swap_val])
            g_var.swap_dict={}

    def test_new_box_vec(self):
        box_vec, box = 'CRYST1   100   100   100  90.00  90.00  90.00 P 1           1', [[50, 50, 50],[50, 0, 50]]
        out = [['CRYST1   50.000   50.000   50.000    90.00    90.00    90.00 P 1           1\n', [25., 25., 25.]], ['CRYST1   50.000  100.000   50.000    90.00    90.00    90.00 P 1           1\n', [25.,  0., 25.]]]
        for check_val, check in enumerate(box): 
            result1, result2 = gen.new_box_vec(box_vec, box[check_val])
            self.assertEqual(result1,out[check_val][0])
            self.assertIsNone(np.testing.assert_array_equal(result2, np.array(out[check_val][1])))

    def test_strip_header(self):
        header_types = [ '[BB]', ' [BB]', '[ BB]', '[ BB ]']
        out = ['BB','BB','BB','BB']
        for header_val, header in enumerate(header_types):
            self.assertEqual(gen.strip_header(header), out[header_val])

    def test_sep_fragments_topology(self):
        results = gen.sep_fragments_topology(run_dir+'files_test/PHE/PHE')
        out = {'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, 'GROUPS': {'group_max': 2, 'SC1': 1, 'SC2': 1, 'SC3': 1}, 'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}}
        self.assertEqual(results, out)

    def test_empty_sep_fragments_topology(self):
        empty = {'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': []}, 'GROUPS': {'group_max': 1}, 'CONNECT': {'atoms': {}}}
        results = gen.sep_fragments_topology(run_dir+'files_test/PHE/missing')
        self.assertEqual(results, empty)

    def test_get_fragment_topology(self):
        results = gen.get_fragment_topology('PHE', 'files_test/PHE/PHE.pdb')
        grouped_atoms = {2: {'BB': [1, 2, 10, 11]}, 1: {'SC1': [3, 4, 5], 'SC2': [8, 9], 'SC3': [6, 7]}}
        res_top_out={'PHE': {'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, 'GROUPS': {'BB': 2, 'SC1': 1, 'SC2': 1, 'SC3': 1}, 'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}, 'ATOMS': ['N', 'CA', 'C', 'O']}}
        self.assertEqual(g_var.res_top, res_top_out)
        self.assertEqual(results, grouped_atoms)

    def test_sort_connectivity(self):
        atom_dict, heavy_bond = {1: {'BB': [1, 2, 10, 11]}, 2: {'SC1': [3, 4, 5]}, 3: {'SC2': [6, 7, 8, 9]}}, {3: [2, 4], 2: [3, 1, 10], 4: [3, 5], 5: [4, 6], 6: [5, 7], 7: [6, 9, 8], 9: [7], 1: [2], 10:[2, 11], 11: [10], 8: [7]}  
        out = {1: {2: ['SC1']}, 2: {3: ['BB'], 5: ['SC2']}, 3: {6: ['SC1']}}
        self.assertEqual(gen.sort_connectivity(atom_dict, heavy_bond), out)

    def test_atom_bond_check(self):
        line_sep = [['[', 'atoms', ']'], ['[', 'bonds', ']'], ['[', 'fail', ']']]
        out = [[True, False],[False, True],[False, False]]
        for header_val, header in enumerate(line_sep):
            result1, result2 = gen.atom_bond_check(header)
            self.assertEqual(result1, out[header_val][0])
            self.assertEqual(result2, out[header_val][1])

    def test_fetch_amino_rtp_file_location(self):
        self.assertEqual(gen.fetch_amino_rtp_file_location(run_dir+'files_test/charmm36-mar2019-updated.ff'), [run_dir+'files_test/charmm36-mar2019-updated.ff/merged.rtp'])

    def test_fetch_atom_masses(self):
        self.assertEqual(gen.fetch_atom_masses(run_dir+'files_test'), {'AG': '107.86820', 'AL': '26.98154'})

    def test_fragment_location(self):
        g_var.np_directories=[[run_dir+'files_test/']]
        self.assertEqual(gen.fragment_location('PHE'), run_dir+'files_test/PHE/PHE.pdb')

    # def test_fetch_bond_info(self):
    #     residue = PHE
    #     amino_acid_itp = ['/home/owen/Documents/scripts/cg2at/database/forcefields/charmm36-mar2019-updated.ff/merged.rtp']
    #     location = '/home/owen/Documents/scripts/cg2at/database/fragments/martini_2-2_charmm36/protein/PHE/PHE.pdb'
    #     fetch_bond_info(residue, amino_acid_itp, g_var.at_mass, location)


####### test at_mod_p
    def test_shrink_coordinates(self):
        p1, p2 = np.array([58.274,66.912,12.038]), np.array([59.360,66.216,14.859])
        p1a, p2a = np.array([58.450475 , 66.7989   , 12.4964125]), np.array([59.183525 , 66.3291   , 14.4005875])
        p1s, p2s = at_mod_p.shrink_coordinates(p1, p2)
        self.assertIsNone(np.testing.assert_array_equal(p1s, p1a))
        self.assertIsNone(np.testing.assert_array_equal(p2s, p2a))




if __name__ in ['__main__', 'test_gen']:
    unittest.main()
