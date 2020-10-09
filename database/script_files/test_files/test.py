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
        results = gen.sep_fragments_topology(run_dir+'files_test/database_test/fragments/test_1/protein/PHE/PHE')
        out = {'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, 'GROUPS': {'group_max': 2, 'SC1': 1, 'SC2': 1, 'SC3': 1}, 'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}}
        self.assertEqual(results, out)

    def test_empty_sep_fragments_topology(self):
        empty = {'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': []}, 'GROUPS': {'group_max': 1}, 'CONNECT': {'atoms': {}}}
        results = gen.sep_fragments_topology(run_dir+'files_test/PHE/missing')
        self.assertEqual(results, empty)

    def test_get_fragment_topology(self):
        g_var.res_top = {}
        results = gen.get_fragment_topology('PHE', run_dir+'files_test/database_test/fragments/test_1/protein/PHE/PHE.pdb')
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
        self.assertEqual(gen.fetch_amino_rtp_file_location(run_dir+'files_test/database_test/forcefields/charmm36-mar2019-updated.ff'), [run_dir+'files_test/database_test/forcefields/charmm36-mar2019-updated.ff/merged.rtp'])

    def test_fetch_atom_masses(self):
        self.assertEqual(gen.fetch_atom_masses(run_dir+'files_test'), {'AG': '107.86820', 'AL': '26.98154'})

    def test_fragment_location(self):
        g_var.np_directories=[[run_dir+'files_test/database_test/fragments/test_1/protein/']]
        self.assertEqual(gen.fragment_location('PHE'), run_dir+'files_test/database_test/fragments/test_1/protein/PHE/PHE.pdb')

    def test_path_leaf(self):
        location = [run_dir+'files_test', run_dir+'files_test/']
        out = [run_dir, 'files_test']
        for loc_val, loc in enumerate(location):
            result1, result2 = gen.path_leaf(loc)
            self.assertEqual(result1, out[0])
            self.assertEqual(result2, out[1])

    def test_read_database_directories(self):
        g_var.database_dir = run_dir+'files_test/database_test/'
        gen.read_database_directories()
        self.assertEqual(g_var.forcefield_available, ['charmm36-mar2019-updated.ff', 'empty.ff'])
        self.assertEqual(g_var.fragments_available, ['test_1', 'test_2'])

    @mock.patch('builtins.input', return_value='0 1')
    def test_ask_database_forcefield_list(self, input):
        forcefields = ['charmm36-mar2019-updated.ff', 'empty.ff']
        self.assertEqual(gen.ask_database(forcefields, 'forcefields', True), True)

    @mock.patch('builtins.input', return_value='3')
    def test_ask_database_forcefield_over(self, input):
        forcefields = ['charmm36-mar2019-updated.ff', 'empty.ff']
        self.assertEqual(gen.ask_database(forcefields, 'forcefields', True), True)

    @mock.patch('builtins.input', return_value='3')
    def test_ask_database_forcefield_over(self, input):
        forcefields = ['charmm36-mar2019-updated.ff', 'empty.ff']
        self.assertEqual(gen.ask_database(forcefields, 'forcefields', True), True)
    @mock.patch('builtins.input', return_value='a')
    def test_ask_database_forcefield_str(self, input):
        forcefields = ['charmm36-mar2019-updated.ff', 'empty.ff']
        self.assertEqual(gen.ask_database(forcefields, 'forcefields', True), True)
    @mock.patch('builtins.input', return_value='0')
    def test_ask_database_forcefield_correct(self, input):
        forcefields = ['charmm36-mar2019-updated.ff', 'empty.ff']
        self.assertEqual(gen.ask_database(forcefields, 'forcefields', True), 0)
    @mock.patch('builtins.input', return_value='0 1')
    def test_ask_database_fragments_list(self, input):
        fragments = ['test_1', 'test_2']
        self.assertIsNone(np.testing.assert_array_equal(gen.ask_database(fragments, 'fragments', True), [0, 1]))
    @mock.patch('builtins.input', return_value='3')
    def test_ask_database_fragments_over(self, input):
        fragments = ['test_1', 'test_2']
        self.assertEqual(gen.ask_database(fragments, 'fragments', True), True)
    @mock.patch('builtins.input', return_value='a')
    def test_ask_database_fragments_str(self, input):
        fragments = ['test_1', 'test_2']
        self.assertEqual(gen.ask_database(fragments, 'fragments', True), True)
    @mock.patch('builtins.input', return_value='0 1')
    def test_database_selection(self, input):
        fragments = ['test_1', 'test_2']
        self.assertIsNone(np.testing.assert_array_equal(gen.ask_database(fragments, 'fragments', True), [0, 1]))

    def test_folder_copy_and_check(self):
        loc = run_dir+'files_test/database_test/forcefields/empty.ff'
        test = run_dir+'files_test/database_test/forcefields/_empty.ff'
        self.assertEqual(os.path.exists(test), False)
        self.assertEqual(os.path.exists(loc), True)
        gen.folder_copy_and_check(loc, test)
        self.assertEqual(os.path.exists(test), True)
        os.system('rm -r '+test)

    def test_forcefield_selection_user(self):
        g_var.args.ff = run_dir+'files_test/database_test/forcefields/empty.ff'
        g_var.final_dir = run_dir+'files_test/database_test/'
        gen.forcefield_selection(True)
        self.assertEqual(g_var.forcefield_location,run_dir+'files_test/database_test/forcefields/')
        self.assertEqual(g_var.forcefield, 'empty.ff')
        self.assertEqual(g_var.opt['ff'], 'empty.ff')
        os.system('rm -r '+g_var.final_dir+g_var.forcefield)

    def test_forcefield_selection_specified(self):
        g_var.args.ff = 'empty.ff'
        g_var.final_dir = run_dir+'files_test/database_test/'
        g_var.database_dir = run_dir+'files_test/database_test/'
        g_var.forcefield_available = ['charmm36-mar2019-updated.ff', 'empty.ff']
        gen.forcefield_selection(True)
        self.assertEqual(g_var.forcefield_location,run_dir+'files_test/database_test/forcefields/')
        self.assertEqual(g_var.forcefield, 'empty.ff')
        self.assertEqual(g_var.opt['ff'], 'empty.ff')

    @mock.patch('builtins.input', return_value='1')
    def test_forcefield_selection_number(self, input):
        g_var.forcefield_location, g_var.forcefield = '',''
        g_var.final_dir = run_dir+'files_test/database_test/'
        g_var.database_dir = run_dir+'files_test/database_test/'
        g_var.forcefield_available = ['charmm36-mar2019-updated.ff', 'empty.ff']
        gen.forcefield_selection(True)
        self.assertEqual(g_var.forcefield_location,run_dir+'files_test/database_test/forcefields/')
        self.assertEqual(g_var.forcefield, 'empty.ff')
        self.assertEqual(g_var.opt['ff'], 'empty.ff')

    def test_fetch_residues(self):
        g_var.np_residues, g_var.p_residues, g_var.mod_residues, g_var.o_residues, g_var.np_directories, g_var.p_directories, g_var.mod_directories, g_var.o_directories = [],[],[],[],[],[],[],[]
        frag_location = [run_dir+'files_test/database_test/fragments/', run_dir+'files_test/database_test/fragments/']
        g_var.fragments_available = ['test_1', 'test_2']
        fragment_number = [0]
        gen.fetch_residues(frag_location, g_var.fragments_available, fragment_number, True)
        self.assertIsNone(np.testing.assert_array_equal(g_var.np_residues, ['SOL']))
        self.assertIsNone(np.testing.assert_array_equal(g_var.p_residues, ['PHE']))
        self.assertIsNone(np.testing.assert_array_equal(g_var.mod_residues, []))
        self.assertIsNone(np.testing.assert_array_equal(g_var.o_residues, []))
        self.assertIsNone(np.testing.assert_array_equal(g_var.np_directories, [[run_dir+'files_test/database_test/fragments/test_1/non_protein/', 'SOL']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.p_directories, [[run_dir+'files_test/database_test/fragments/test_1/protein/', 'PHE']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.mod_directories, [[run_dir+'files_test/database_test/fragments/test_1/protein_modified/']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.o_directories, [[run_dir+'files_test/database_test/fragments/test_1/other/']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.opt['fg'], ['test_1']))
        

    def test_fetch_frag_number(self):
        g_var.args.fg = ['test_1']
        fragments_available = ['test_1', 'test_2']
        results = gen.fetch_frag_number(fragments_available)
        self.assertEqual(results, [0])


    def test_fragment_selection(self):
        g_var.np_residues, g_var.p_residues, g_var.mod_residues, g_var.o_residues, g_var.np_directories, g_var.p_directories, g_var.mod_directories, g_var.o_directories = [],[],[],[],[],[],[],[]
        g_var.opt['fg'] = ''
        g_var.args.fg = ['test_1']
        g_var.database_dir = run_dir+'files_test/database_test/'
        g_var.fragments_available = ['test_1', 'test_2']
        gen.fragment_selection(True)
        self.assertIsNone(np.testing.assert_array_equal(g_var.np_residues, ['SOL']))
        self.assertIsNone(np.testing.assert_array_equal(g_var.p_residues, ['PHE']))
        self.assertIsNone(np.testing.assert_array_equal(g_var.mod_residues, []))
        self.assertIsNone(np.testing.assert_array_equal(g_var.o_residues, []))
        self.assertIsNone(np.testing.assert_array_equal(g_var.np_directories, [[run_dir+'files_test/database_test/fragments/test_1/non_protein/', 'SOL']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.p_directories, [[run_dir+'files_test/database_test/fragments/test_1/protein/', 'PHE']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.mod_directories, [[run_dir+'files_test/database_test/fragments/test_1/protein_modified/']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.o_directories, [[run_dir+'files_test/database_test/fragments/test_1/other/']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.opt['fg'], ['test_1']))

    @mock.patch('builtins.input', return_value='0')
    def test_fragment_selection(self, input):
        g_var.np_residues, g_var.p_residues, g_var.mod_residues, g_var.o_residues, g_var.np_directories, g_var.p_directories, g_var.mod_directories, g_var.o_directories = [],[],[],[],[],[],[],[]
        g_var.opt['fg'] = ''
        g_var.args.fg = None
        g_var.database_dir = run_dir+'files_test/database_test/'
        g_var.fragments_available = ['test_1', 'test_2']
        gen.fragment_selection(True)
        self.assertIsNone(np.testing.assert_array_equal(g_var.np_residues, ['SOL']))
        self.assertIsNone(np.testing.assert_array_equal(g_var.p_residues, ['PHE']))
        self.assertIsNone(np.testing.assert_array_equal(g_var.mod_residues, []))
        self.assertIsNone(np.testing.assert_array_equal(g_var.o_residues, []))
        self.assertIsNone(np.testing.assert_array_equal(g_var.np_directories, [[run_dir+'files_test/database_test/fragments/test_1/non_protein/', 'SOL']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.p_directories, [[run_dir+'files_test/database_test/fragments/test_1/protein/', 'PHE']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.mod_directories, [[run_dir+'files_test/database_test/fragments/test_1/protein_modified/']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.o_directories, [[run_dir+'files_test/database_test/fragments/test_1/other/']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.opt['fg'], ['test_1']))

    @mock.patch('builtins.input', return_value='0')
    def test_ask_for_water_model(self, input):
        water = ['tip3p', 'tip4p', 'spc', 'spce']
        directory = [run_dir+'files_test/database_test/fragments/test_1/non_protein/']
        result1, result2 = gen.ask_for_water_model(directory, water)
        self.assertEqual(result1, run_dir+'files_test/database_test/fragments/test_1/non_protein/SOL/')
        self.assertEqual(result2, 'tip3p')

    @mock.patch('builtins.input', return_value='0')
    def test_check_water_molecules(self, input):
        g_var.np_directories = [[run_dir+'files_test/database_test/fragments/test_1/non_protein/', 'SOL']]
        gen.check_water_molecules(True)
        self.assertIsNone(np.testing.assert_array_equal(g_var.water_info, [[run_dir+'files_test/database_test/fragments/test_1/non_protein/SOL', 'tip3p', 'tip4p', 'spc', 'spce']]))
        self.assertEqual(g_var.water_dir, run_dir+'files_test/database_test/fragments/test_1/non_protein/SOL/')
        self.assertEqual(g_var.water, 'tip3p')
        self.assertEqual(g_var.opt['w'], 'tip3p')

    def test_fetch_bond_info_atoms_heavy(self):
        residue, line_sep, residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict = 'PHE', ['N', 'NH1', '-0.470', '0'], [], {}, [], {}, []
        at_mass = {'NH1': 14.007, 'CT1': 12.011, 'CT2': 12.011, 'CA': 12.011, 'C': 12.011, 'O': 15.999}
        residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict = gen.fetch_bond_info_atoms(residue, line_sep, residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict, at_mass)
        self.assertEqual(residue_list, ['PHE'])
        self.assertEqual(atom_conversion, {'N': 1})
        self.assertEqual(H_dict, [])
        self.assertEqual(res_at_mass, {'N': 14.007})
        self.assertEqual(heavy_dict, ['N'])

    def test_fetch_bond_info_atoms_H(self):
        residue, line_sep, residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict = 'PHE', ['HN', 'H', '0.310', '1'], ['PHE'], {'N': 1}, [], {'N': 14.007}, ['N']
        at_mass = {'NH1': 14.007, 'CT1': 12.011, 'CT2': 12.011, 'CA': 12.011, 'C': 12.011, 'O': 15.999}
        residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict = gen.fetch_bond_info_atoms(residue, line_sep, residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict, at_mass)
        self.assertEqual(residue_list, ['PHE'])
        self.assertEqual(atom_conversion, {'N': 1, 'HN': 2})
        self.assertEqual(H_dict, ['HN'])
        self.assertEqual(res_at_mass, {'N': 14.007})
        self.assertEqual(heavy_dict, ['N'])

    def test_fetch_bond_info(self):
        g_var.forcefield_location = run_dir+'files_test/database_test/forcefields/'
        g_var.forcefield = 'charmm36-mar2019-updated.ff'
        g_var.p_residues =  ['PHE']
        g_var.p_directories =[[run_dir+'files_test/database_test/fragments/test_1/protein/', 'PHE']]
        residue = 'PHE'
        amino_acid_itp = [run_dir+'files_test/database_test/forcefields/charmm36-mar2019-updated.ff/merged.rtp']
        location = run_dir+'files_test/database_test/fragments/test_1/protein/PHE/PHE.pdb'
        res_at_mass = {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'CG': 12.011, 'CD1': 12.011, 'CE1': 12.011, 'CZ': 12.011, 'CD2': 12.011, 'CE2': 12.011, 'C': 12.011, 'O': 15.999}
        at_mass = {'NH1': 14.007, 'CT1': 12.011, 'CT2': 12.011, 'CA': 12.011, 'C': 12.011, 'O': 15.999}
        # hydrogen, heavy_bond, residue_list, res_at_mass, amide_hydrogen = gen.fetch_bond_info(residue, amino_acid_itp, at_mass, location)
        # print(hydrogen)
        # print(heavy_bond)
        # print(residue_list)
        # print(res_at_mass)
        # print(amide_hydrogen)
# {'N': ['HN'], 'CA': ['HA'], 'CB': ['HB1', 'HB2'], 'CD1': ['HD1'], 'CD2': ['HD2'], 'CE1': ['HE1'], 'CE2': ['HE2'], 'CZ': ['HZ']}
# {3: [2, 4], 2: [3, 1, 10], 4: [3, 8, 5], 8: [4, 9], 6: [5, 7], 5: [6, 4], 7: [9, 6], 9: [7, 8], 1: [2], 10: [2, 11], 11: [10]}
# ['PHE']
# {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'CG': 12.011, 'CD1': 12.011, 'CE1': 12.011, 'CZ': 12.011, 'CD2': 12.011, 'CE2': 12.011, 'C': 12.011, 'O': 15.999}
# HN




    # def test_fetch_fragment(self):
    #     g_var.forcefield_location = run_dir+'files_test/database_test/forcefields/'
    #     g_var.forcefield = 'charmm36-mar2019-updated.ff'
    #     g_var.p_residues =  ['PHE']
    #     g_var.p_directories =[[run_dir+'files_test/database_test/fragments/test_1/protein/', 'PHE']]
    #     gen.fetch_fragment()







####### test at_mod_p
    def test_shrink_coordinates(self):
        p1, p2 = np.array([58.274,66.912,12.038]), np.array([59.360,66.216,14.859])
        p1a, p2a = np.array([58.450475 , 66.7989   , 12.4964125]), np.array([59.183525 , 66.3291   , 14.4005875])
        p1s, p2s = at_mod_p.shrink_coordinates(p1, p2)
        self.assertIsNone(np.testing.assert_array_equal(p1s, p1a))
        self.assertIsNone(np.testing.assert_array_equal(p2s, p2a))




if __name__ in ['__main__', 'test_gen']:
    unittest.main()
