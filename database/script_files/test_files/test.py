#!/usr/bin/env python3
import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/..')
import unittest
from unittest.mock import patch
import mock
import filecmp
import numpy as np
from scipy.spatial import cKDTree
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
            g_var.swap_dict={}
            gen.sort_swap_group()
            self.assertEqual(g_var.swap_dict, out[swap_val])
            g_var.swap_dict={}

    def test_print_swap_residues(self):
        g_var.args.swap=True
        g_var.swap_dict = {'NA+': {'NA+:skip': {'ALL': 'ALL', 'resid': [4000, 4001, 4002], 'range': ['4000-4002']}}} 
        correct = '\nYou have chosen to swap the following residues\n\n residue  bead              residue     bead       range   \n -------  ----              -------     ----       -----   \n   NA+    ALL      -->       skip       ALL      4000-4002 \n'
        to_print = gen.print_swap_residues()       
        self.assertEqual(to_print, correct)

    def test_new_box_vec(self):
        box_vec, box = 'CRYST1   100   100   100  90.00  90.00  90.00 P 1           1', [[50, 50, 50],[50, 0, 50]]
        out = [['CRYST1   50.000   50.000   50.000    90.00    90.00    90.00 P 1           1\n', [25., 25., 25.]], ['CRYST1   50.000  100.000   50.000    90.00    90.00    90.00 P 1           1\n', [25.,  0., 25.]]]
        for check_val, check in enumerate(box): 
            result1, result2 = gen.new_box_vec(box_vec, box[check_val])
            self.assertEqual(result1,out[check_val][0])
            self.assertIsNone(np.testing.assert_array_equal(result2, np.array(out[check_val][1])))

    def test_strip_header(self):
        header_types = [ '[BB]', ' [BB]', '[ BB]', '[ BB ]', '[ BB AA ]']
        out = ['BB','BB','BB','BB']
        with self.assertRaises(SystemExit) as cm:
            for header_val, header in enumerate(header_types):
                self.assertEqual(gen.strip_header(header), out[header_val])
        self.assertEqual(cm.exception.code, 'There is a issue in one of the fragment headers: \n[ BB AA ]')

    def test_sep_fragments_topology(self):
        results = gen.sep_fragments_topology(run_dir+'database_test/fragments/test_1/protein/PHE/PHE')
        out = {'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, 'GROUPS': {'group_max': 2, 'SC1': 1, 'SC2': 1, 'SC3': 1}, 'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}}
        self.assertEqual(results, out)

    def test_empty_sep_fragments_topology(self):
        empty = {'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': []}, 'GROUPS': {'group_max': 1}, 'CONNECT': {'atoms': {}}}
        results = gen.sep_fragments_topology(run_dir+'PHE/missing')
        self.assertEqual(results, empty)

    def test_get_fragment_topology(self):
        g_var.res_top = {}
        results = gen.get_fragment_topology('PHE', run_dir+'database_test/fragments/test_1/protein/PHE/PHE.pdb')
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
        self.assertEqual(gen.fetch_amino_rtp_file_location(run_dir+'database_test/forcefields/charmm36-mar2019-updated.ff'), [run_dir+'database_test/forcefields/charmm36-mar2019-updated.ff/merged.rtp'])

    def test_fetch_atom_masses(self):
        self.assertEqual(gen.fetch_atom_masses(run_dir+'database_test'), {'AG': '107.86820', 'AL': '26.98154'})

    def test_fragment_location(self):
        g_var.np_directories=[[run_dir+'database_test/fragments/test_1/protein/']]
        self.assertEqual(gen.fragment_location('PHE'), run_dir+'database_test/fragments/test_1/protein/PHE/PHE.pdb')

    def test_path_leaf(self):
        location = [run_dir+'database_test', run_dir+'database_test/']
        out = [run_dir, 'database_test']
        for loc_val, loc in enumerate(location):
            result1, result2 = gen.path_leaf(loc)
            self.assertEqual(result1, out[0])
            self.assertEqual(result2, out[1])

    def test_read_database_directories(self):
        g_var.database_dir = run_dir+'database_test/'
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
        loc = run_dir+'database_test/forcefields/empty.ff'
        test = run_dir+'database_test/forcefields/_empty.ff'
        self.assertEqual(os.path.exists(test), False)
        self.assertEqual(os.path.exists(loc), True)
        gen.folder_copy_and_check(loc, test)
        self.assertEqual(os.path.exists(test), True)
        os.system('rm -r '+test)

    def test_forcefield_selection_user(self):
        g_var.args.ff = run_dir+'database_test/forcefields/empty.ff'
        g_var.final_dir = run_dir+'database_test/'
        gen.forcefield_selection(True)
        self.assertEqual(g_var.forcefield_location,run_dir+'database_test/forcefields/')
        self.assertEqual(g_var.forcefield, 'empty.ff')
        self.assertEqual(g_var.opt['ff'], 'empty.ff')
        os.system('rm -r '+g_var.final_dir+g_var.forcefield)

    def test_forcefield_selection_specified(self):
        g_var.args.ff = 'empty.ff'
        g_var.final_dir = run_dir+'database_test/'
        g_var.database_dir = run_dir+'database_test/'
        g_var.forcefield_available = ['charmm36-mar2019-updated.ff', 'empty.ff']
        gen.forcefield_selection(True)
        self.assertEqual(g_var.forcefield_location,run_dir+'database_test/forcefields/')
        self.assertEqual(g_var.forcefield, 'empty.ff')
        self.assertEqual(g_var.opt['ff'], 'empty.ff')

    @mock.patch('builtins.input', return_value='1')
    def test_forcefield_selection_number(self, input):
        g_var.forcefield_location, g_var.forcefield = '',''
        g_var.final_dir = run_dir+'database_test/'
        g_var.database_dir = run_dir+'database_test/'
        g_var.forcefield_available = ['charmm36-mar2019-updated.ff', 'empty.ff']
        gen.forcefield_selection(True)
        self.assertEqual(g_var.forcefield_location, run_dir+'database_test/forcefields/')
        self.assertEqual(g_var.forcefield, 'empty.ff')
        self.assertEqual(g_var.opt['ff'], 'empty.ff')

    def test_fetch_residues(self):
        g_var.np_residues, g_var.p_residues, g_var.mod_residues, g_var.o_residues, g_var.np_directories, g_var.p_directories, g_var.mod_directories, g_var.o_directories = [],[],[],[],[],[],[],[]
        frag_location = [run_dir+'database_test/fragments/', run_dir+'database_test/fragments/']
        g_var.fragments_available = ['test_1', 'test_2']
        fragment_number = [0]
        gen.fetch_residues(frag_location, g_var.fragments_available, fragment_number, True)
        self.assertIsNone(np.testing.assert_array_equal(g_var.np_residues, ['CHOL', 'ION', 'SOL']))
        self.assertIsNone(np.testing.assert_array_equal(g_var.p_residues, ['PHE']))
        self.assertIsNone(np.testing.assert_array_equal(g_var.mod_residues, []))
        self.assertIsNone(np.testing.assert_array_equal(g_var.o_residues, []))
        self.assertIsNone(np.testing.assert_array_equal(g_var.np_directories, [[run_dir+'database_test/fragments/test_1/non_protein/', 'CHOL', 'ION', 'SOL']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.p_directories, [[run_dir+'database_test/fragments/test_1/protein/', 'PHE']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.mod_directories, [[run_dir+'database_test/fragments/test_1/protein_modified/']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.o_directories, [[run_dir+'database_test/fragments/test_1/other/']]))
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
        g_var.database_dir = run_dir+'database_test/'
        g_var.fragments_available = ['test_1', 'test_2']
        gen.fragment_selection(True)
        self.assertIsNone(np.testing.assert_array_equal(g_var.np_residues, ['CHOL', 'ION', 'SOL']))
        self.assertIsNone(np.testing.assert_array_equal(g_var.p_residues, ['PHE']))
        self.assertIsNone(np.testing.assert_array_equal(g_var.mod_residues, []))
        self.assertIsNone(np.testing.assert_array_equal(g_var.o_residues, []))
        self.assertIsNone(np.testing.assert_array_equal(g_var.np_directories, [[run_dir+'database_test/fragments/test_1/non_protein/', 'CHOL', 'ION', 'SOL']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.p_directories, [[run_dir+'database_test/fragments/test_1/protein/', 'PHE']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.mod_directories, [[run_dir+'database_test/fragments/test_1/protein_modified/']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.o_directories, [[run_dir+'database_test/fragments/test_1/other/']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.opt['fg'], ['test_1']))

    @mock.patch('builtins.input', return_value='0')
    def test_fragment_selection(self, input):
        g_var.np_residues, g_var.p_residues, g_var.mod_residues, g_var.o_residues, g_var.np_directories, g_var.p_directories, g_var.mod_directories, g_var.o_directories = [],[],[],[],[],[],[],[]
        g_var.opt['fg'] = ''
        g_var.args.fg = None
        g_var.database_dir = run_dir+'database_test/'
        g_var.fragments_available = ['test_1', 'test_2']
        gen.fragment_selection(True)
        self.assertIsNone(np.testing.assert_array_equal(g_var.np_residues, ['CHOL', 'ION', 'SOL']))
        self.assertIsNone(np.testing.assert_array_equal(g_var.p_residues, ['PHE']))
        self.assertIsNone(np.testing.assert_array_equal(g_var.mod_residues, []))
        self.assertIsNone(np.testing.assert_array_equal(g_var.o_residues, []))
        self.assertIsNone(np.testing.assert_array_equal(g_var.np_directories, [[run_dir+'database_test/fragments/test_1/non_protein/', 'CHOL', 'ION', 'SOL']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.p_directories, [[run_dir+'database_test/fragments/test_1/protein/', 'PHE']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.mod_directories, [[run_dir+'database_test/fragments/test_1/protein_modified/']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.o_directories, [[run_dir+'database_test/fragments/test_1/other/']]))
        self.assertIsNone(np.testing.assert_array_equal(g_var.opt['fg'], ['test_1']))

    @mock.patch('builtins.input', return_value='0')
    def test_ask_for_water_model(self, input):
        water = ['tip3p', 'tip4p', 'spc', 'spce']
        directory = [run_dir+'database_test/fragments/test_1/non_protein/']
        result1, result2 = gen.ask_for_water_model(directory, water)
        self.assertEqual(result1, run_dir+'database_test/fragments/test_1/non_protein/SOL/')
        self.assertEqual(result2, 'tip3p')

    @mock.patch('builtins.input', return_value='0')
    def test_check_water_molecules(self, input):
        g_var.np_directories = [[run_dir+'database_test/fragments/test_1/non_protein/', 'SOL']]
        gen.check_water_molecules(True)
        self.assertIsNone(np.testing.assert_array_equal(g_var.water_info, [[run_dir+'database_test/fragments/test_1/non_protein/SOL', 'tip3p', 'tip4p', 'spc', 'spce']]))
        self.assertEqual(g_var.water_dir, run_dir+'database_test/fragments/test_1/non_protein/SOL/')
        self.assertEqual(g_var.water, 'tip3p')
        self.assertEqual(g_var.opt['w'], 'tip3p')

    def test_fetch_bond_info_atoms_heavy_linked(self):
        residue, line_sep, residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict = 'PHE', ['N', 'NH1', '-0.470', '0'], [], {}, [], {}, []
        at_mass = {'NH1': 14.007, 'CT1': 12.011, 'CT2': 12.011, 'CA': 12.011, 'C': 12.011, 'O': 15.999}
        residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict = gen.fetch_bond_info_atoms_linked(residue, line_sep, residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict, at_mass)
        self.assertEqual(residue_list, ['PHE'])
        self.assertEqual(atom_conversion, {'N': 1})
        self.assertEqual(H_dict, [])
        self.assertEqual(res_at_mass, {'N': 14.007})
        self.assertEqual(heavy_dict, ['N'])

    def test_fetch_bond_info_atoms_H_linked(self):
        residue, line_sep, residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict = 'PHE', ['HN', 'H', '0.310', '1'], ['PHE'], {'N': 1}, [], {'N': 14.007}, ['N']
        at_mass = {'NH1': 14.007, 'CT1': 12.011, 'CT2': 12.011, 'CA': 12.011, 'C': 12.011, 'O': 15.999}
        residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict = gen.fetch_bond_info_atoms_linked(residue, line_sep, residue_list, atom_conversion, H_dict, res_at_mass, heavy_dict, at_mass)
        self.assertEqual(residue_list, ['PHE'])
        self.assertEqual(atom_conversion, {'N': 1, 'HN': 2})
        self.assertEqual(H_dict, ['HN'])
        self.assertEqual(res_at_mass, {'N': 14.007})
        self.assertEqual(heavy_dict, ['N'])

    def test_fetch_bond_info_atoms_heavy_NP(self):
        residue, line_sep, residue_list, H_dict, res_at_mass, heavy_dict = 'CHOL', ['1', 'CRL1', '1', 'CHOL', 'C3', '1', '0.14', '12.011', ';', 'qtot', '0.14'], [], [], {}, []
        residue_list, H_dict, res_at_mass, heavy_dict = gen.fetch_bond_info_atoms_NP(residue, line_sep, residue_list, H_dict, res_at_mass, heavy_dict)
        self.assertEqual(residue_list, ['CHOL'])
        self.assertEqual(H_dict, [])
        self.assertEqual(res_at_mass, {'C3': 12.011})
        self.assertEqual(heavy_dict, [1])

    def test_fetch_bond_info_atoms_H_NP(self):
        residue, line_sep, residue_list, H_dict, res_at_mass, heavy_dict = 'CHOL', ['2', 'HGA1', '1', 'CHOL', 'H3', '2', '0.09', '1.008', ';', 'qtot', '0.23'], ['CHOL'], [], {'C3': 12.011}, [1]
        residue_list, H_dict, res_at_mass, heavy_dict = gen.fetch_bond_info_atoms_NP(residue, line_sep, residue_list, H_dict, res_at_mass, heavy_dict)
        self.assertEqual(residue_list, ['CHOL'])
        self.assertEqual(H_dict, [2])
        self.assertEqual(res_at_mass, {'C3': 12.011})
        self.assertEqual(heavy_dict, [1])

    def test_get_atomistic_P(self):
        location = run_dir+'database_test/fragments/test_1/protein/PHE/PHE.pdb'
        atom_conversion = gen.get_atomistic(location)
        self.assertEqual(atom_conversion, {'N': 1, 'CA': 2, 'C': 10, 'O': 11, 'CB': 3, 'CG': 4, 'CD1': 5, 'CD2': 8, 'CE2': 9, 'CE1': 6, 'CZ': 7})

    def test_get_atomistic_NP(self):
        correct={'C3': 1, 'H3': 2, 'O3': 3, "H3'": 4, 'C4': 5, 'H4A': 6, 'H4B': 7, 'C5': 8, 'C10': 39, 'C19': 40, 'H19A': 41, 'H19B': 42, 'H19C': 43,\
         'C1': 44, 'H1A': 45, 'H1B': 46, 'C2': 47, 'H2A': 48, 'H2B': 49, 'C6': 9, 'H6': 10, 'C7': 11, 'H7A': 12, 'H7B': 13, 'C8': 14, 'H8': 15, 'C12': 31,\
          'H12A': 32, 'H12B': 33, 'C11': 34, 'H11A': 35, 'H11B': 36, 'C9': 37, 'H9': 38, 'C14': 16, 'H14': 17, 'C15': 18, 'H15A': 19, 'H15B': 20, 'C16': 21,\
           'H16A': 22, 'H16B': 23, 'C17': 24, 'H17': 25, 'C13': 26, 'C18': 27, 'H18A': 28, 'H18B': 29, 'H18C': 30, 'C20': 50, 'H20': 51, 'C21': 52, 'H21A': 53,\
            'H21B': 54, 'H21C': 55, 'C22': 56, 'H22A': 57, 'H22B': 58, 'C23': 59, 'H23A': 60, 'H23B': 61, 'C24': 62, 'H24A': 63, 'H24B': 64, 'C25': 65, 'H25': 66,\
             'C26': 67, 'H26A': 68, 'H26B': 69, 'H26C': 70, 'C27': 71, 'H27A': 72, 'H27B': 73, 'H27C': 74}
        location = run_dir+'database_test/fragments/test_1/non_protein/CHOL/CHOL.pdb'
        atom_conversion = gen.get_atomistic(location)
        self.assertEqual(atom_conversion,correct)
        
    def test_add_to_topology_list_heavy_h_h_bond(self):
        bond = ['N', 'HN']
        heavy_dict = ['N', 'CA', 'CB', 'CG', 'CD1', 'CE1', 'CZ', 'CD2', 'CE2', 'C', 'O']
        H_dict = ['HN', 'HA', 'HB1', 'HB2', 'HD1', 'HE1', 'HZ', 'HD2', 'HE2']
        atom_conversion = {'N': 1, 'CA': 2, 'C': 10, 'O': 11, 'CB': 3, 'CG': 4, 'CD1': 5, 'CD2': 8, 'CE2': 9, 'CE1': 6, 'CZ': 7}
        hydrogen, amide_h = gen.add_to_topology_list(bond[0], bond[1], {}, heavy_dict, H_dict, atom_conversion, 'PHE', ['PHE'])
        self.assertEqual(hydrogen,{'N': ['HN']}) 
        self.assertEqual(amide_h, 'HN')

    def test_add_to_topology_list_heavy_h_heavy_bond(self):
        bond = ['CB', 'CA']
        heavy_dict = ['N', 'CA', 'CB', 'CG', 'CD1', 'CE1', 'CZ', 'CD2', 'CE2', 'C', 'O']
        H_dict = ['HN', 'HA', 'HB1', 'HB2', 'HD1', 'HE1', 'HZ', 'HD2', 'HE2']
        atom_conversion = {'N': 1, 'CA': 2, 'C': 10, 'O': 11, 'CB': 3, 'CG': 4, 'CD1': 5, 'CD2': 8, 'CE2': 9, 'CE1': 6, 'CZ': 7}
        hydrogen, amide_h = gen.add_to_topology_list(bond[0], bond[1], {}, heavy_dict, H_dict, atom_conversion, 'PHE', ['PHE'])
        self.assertEqual(hydrogen,{}) 
        self.assertEqual(amide_h, None)

    def test_add_to_topology_list_heavy_heavy_h_bond(self):
        bond = ['N', 'HN']
        heavy_dict = ['N', 'CA', 'CB', 'CG', 'CD1', 'CE1', 'CZ', 'CD2', 'CE2', 'C', 'O']
        H_dict = ['HN', 'HA', 'HB1', 'HB2', 'HD1', 'HE1', 'HZ', 'HD2', 'HE2']
        atom_conversion = {'N': 1, 'CA': 2, 'C': 10, 'O': 11, 'CB': 3, 'CG': 4, 'CD1': 5, 'CD2': 8, 'CE2': 9, 'CE1': 6, 'CZ': 7}
        hydrogen, amide_h = gen.add_to_topology_list(bond[0], bond[1], {}, heavy_dict, heavy_dict, atom_conversion, 'PHE', ['PHE'])
        self.assertEqual(hydrogen,{}) 
        self.assertEqual(amide_h, None)

    def test_add_to_topology_list_heavy_heavy_heavy_bond(self):
        bond = ['CB', 'CA']
        heavy_dict = ['N', 'CA', 'CB', 'CG', 'CD1', 'CE1', 'CZ', 'CD2', 'CE2', 'C', 'O']
        H_dict = ['HN', 'HA', 'HB1', 'HB2', 'HD1', 'HE1', 'HZ', 'HD2', 'HE2']
        atom_conversion = {'N': 1, 'CA': 2, 'C': 10, 'O': 11, 'CB': 3, 'CG': 4, 'CD1': 5, 'CD2': 8, 'CE2': 9, 'CE1': 6, 'CZ': 7}
        hydrogen, amide_h = gen.add_to_topology_list(bond[0], bond[1], {}, heavy_dict, heavy_dict, atom_conversion, 'PHE', ['PHE'])
        self.assertEqual(hydrogen,{3: [2], 2: [3]}) 
        self.assertEqual(amide_h, None)

    def test_fetch_bond_info_P(self):
        g_var.forcefield_location = run_dir+'database_test/forcefields/'
        g_var.forcefield = 'charmm36-mar2019-updated.ff'
        g_var.p_residues =  ['PHE']
        g_var.p_directories =[[run_dir+'database_test/fragments/test_1/protein/', 'PHE']]
        residue = 'PHE'
        amino_acid_itp = [run_dir+'database_test/forcefields/charmm36-mar2019-updated.ff/merged.rtp']
        location = run_dir+'database_test/fragments/test_1/protein/PHE/PHE.pdb'
        res_at_mass = {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'CG': 12.011, 'CD1': 12.011, 'CE1': 12.011, 'CZ': 12.011, 'CD2': 12.011, 'CE2': 12.011, 'C': 12.011, 'O': 15.999}
        at_mass = {'NH1': 14.007, 'CT1': 12.011, 'CT2': 12.011, 'CA': 12.011, 'C': 12.011, 'O': 15.999}
        hydrogen, heavy_bond, residue_list, res_at_mass, amide_hydrogen = gen.fetch_bond_info(residue, amino_acid_itp, at_mass, location)
        self.assertEqual(hydrogen, {'N': ['HN'], 'CA': ['HA'], 'CB': ['HB1', 'HB2'], 'CD1': ['HD1'], 'CD2': ['HD2'], 'CE1': ['HE1'], 'CE2': ['HE2'], 'CZ': ['HZ']})
        self.assertEqual(heavy_bond, {3: [2, 4], 2: [3, 1, 10], 4: [3, 8, 5], 8: [4, 9], 6: [5, 7], 5: [6, 4], 7: [9, 6], 9: [7, 8], 1: [2], 10: [2, 11], 11: [10]})
        self.assertEqual(residue_list, ['PHE'])
        self.assertEqual(res_at_mass, {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'CG': 12.011, 'CD1': 12.011, 'CE1': 12.011, 'CZ': 12.011, 'CD2': 12.011, 'CE2': 12.011, 'C': 12.011, 'O': 15.999})
        self.assertEqual(amide_hydrogen, 'HN')

    def test_fetch_bond_info_NP(self):
        g_var.forcefield_location = run_dir+'database_test/forcefields/'
        g_var.forcefield = 'charmm36-mar2019-updated.ff'
        g_var.np_residues =  ['CHOL']
        g_var.np_directories = [[run_dir+'database_test/fragments/test_1/non_protein/', 'CHOL',  'SOL']]
        residue = 'CHOL'
        NP_itp = [run_dir+'database_test/fragments/test_1/non_protein/CHOL/CHOL.itp']
        location = run_dir+'database_test/fragments/test_1/non_protein/CHOL/CHOL.pdb'
        at_mass = {'NH1': 14.007, 'CT1': 12.011, 'CT2': 12.011, 'CA': 12.011, 'C': 12.011, 'O': 15.999}

        hydrogen_correct = {1: [2], 3: [4], 5: [6, 7], 9: [10], 11: [12, 13], 14: [15], 16: [17],\
         18: [19, 20], 21: [22, 23], 24: [25], 27: [28, 29, 30], 31: [32, 33], 34: [35, 36], 37: [38],\
          40: [41, 42, 43], 44: [45, 46], 47: [48, 49], 50: [51], 52: [53, 54, 55], 56: [57, 58], 59: [60, 61],\
           62: [63, 64], 65: [66], 67: [68, 69, 70], 71: [72, 73, 74]} 
        heavy_bond_correct = {1: [3, 5, 47], 3: [1], 5: [1, 8], 47: [1, 44], 8: [5, 9, 39], 9: [8, 11], 39: [8, 37, 40, 44],\
         11: [9, 14], 14: [11, 16, 37], 16: [14, 18, 26], 37: [14, 34, 39], 18: [16, 21], 26: [16, 24, 27, 31], 21: [18, 24],\
          24: [21, 26, 50], 50: [24, 52, 56], 27: [26], 31: [26, 34], 34: [31, 37], 40: [39], 44: [39, 47], 52: [50], 56: [50, 59],\
           59: [56, 62], 62: [59, 65], 65: [62, 67, 71], 67: [65], 71: [65]}
        res_at_mass_correct = {'C3': 12.011, 'O3': 15.9994, 'C4': 12.011, 'C5': 12.011, 'C6': 12.011, 'C7': 12.011, 'C8': 12.011,\
         'C14': 12.011, 'C15': 12.011, 'C16': 12.011, 'C17': 12.011, 'C13': 12.011, 'C18': 12.011, 'C12': 12.011, 'C11': 12.011, \
         'C9': 12.011, 'C10': 12.011, 'C19': 12.011, 'C1': 12.011, 'C2': 12.011, 'C20': 12.011, 'C21': 12.011, 'C22': 12.011, \
         'C23': 12.011, 'C24': 12.011, 'C25': 12.011, 'C26': 12.011, 'C27': 12.011}
        hydrogen, heavy_bond, residue_list, res_at_mass, amide_hydrogen = gen.fetch_bond_info(residue, NP_itp, at_mass, location)
        self.assertEqual(hydrogen, hydrogen_correct)
        self.assertEqual(heavy_bond, heavy_bond_correct)
        self.assertEqual(residue_list, ['CHOL'])
        self.assertEqual(res_at_mass, res_at_mass_correct)
        self.assertEqual(amide_hydrogen, None)

    def test_fetch_atoms_water_ion(self):
        at_mass_water = {'OT':15.99940, 'OW':15.99940, 'OWT4':15.99940, 'MWT4':15.99940, 'HT':1, 'HW':1, 'HWT4':1}
        at_mass = gen.fetch_atoms_water_ion(run_dir+'database_test/fragments/test_1/non_protein/SOL/', at_mass_water)
        self.assertEqual(at_mass, {'OW': 15.9994, 'HW1': 1.0, 'HW2': 1.0, 'MW': 15.9994})

    def test_create_ion_list(self):
        gen.create_ion_list(run_dir+'database_test/fragments/test_1/non_protein/ION/ION.pdb')
        self.assertEqual(g_var.ions, ['NA+', 'NA', 'CL-', 'CL', 'K+', 'K'])
        

    def test_fetch_fragment(self):
        g_var.res_top = {}
        res_top_correct = {'C_TERMINAL': 'default', 'N_TERMINAL': 'default', \
        'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, \
        'GROUPS': {'BB': 2, 'SC1': 1, 'SC2': 1, 'SC3': 1}, \
        'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}, \
        'ATOMS': ['N', 'CA', 'C', 'O'], \
        'RESIDUE': ['PHE'], \
        'atom_masses': {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'CG': 12.011, 'CD1': 12.011, 'CE1': 12.011, 'CZ': 12.011, \
                        'CD2': 12.011, 'CE2': 12.011, 'C': 12.011, 'O': 15.999}, 'amide_h': 'HN'}
        g_var.forcefield_location = run_dir+'database_test/forcefields/'
        g_var.forcefield = 'charmm36-mar2019-updated.ff'
        g_var.p_residues =  ['PHE']
        g_var.p_directories =[[run_dir+'database_test/fragments/test_1/protein/', 'PHE']]
        gen.fetch_fragment()
        self.assertEqual(g_var.res_top['PHE'], res_top_correct)
        self.assertEqual(len(g_var.res_top), 3)

    def test_fetch_chain_groups(self):
        test = [['1,2','3,4'], ['chain'], ['all']]
        out = [{1: 0, 2: 0, 3: 1, 4: 1}, 'chain', 'all']
        for scene_val, scene in enumerate(test):
            g_var.group_chains = {}
            g_var.args.group = scene
            gen.fetch_chain_groups()
            self.assertEqual(g_var.group_chains, out[scene_val])
            g_var.group_chains = None

    def test_AnglesToRotMat(self):
        R = gen.AnglesToRotMat(np.array([1.2, 1.2, 1.2]))
        out = np.array([[ 0.13130314, -0.02295255,  0.99107652],[ 0.33773159,  0.94096257, -0.02295255], [-0.93203909,  0.33773159,  0.13130314]])
        for i in range(3):
            np.testing.assert_array_almost_equal(R, out)

    def test_angle_clockwise(self):
        A = np.array([-1.43688546, -0.27613562] ) 
        B = np.array([ 0.90246221, -1.1875225 ])
        angle_result = gen.angle_clockwise(A, B)
        np.testing.assert_array_almost_equal(angle_result, 243.64512301056095)

    def test_print_sequence_info(self):
        correct = 'Summary of coarsegrain PROTEIN chains\n\n chain number  length of chain\n\n ------------  ---------------\n       0            3      \nSequences:\n\nchain: 0\n1        10        20        30        40        50        60        70        \nAAA                                                                             \n\nSummary of atomistic PROTEIN chains\n\n chain number  length of chain\n\n ------------  ---------------\n       0            3      \nSequences:\n\nchain: 0\n1        10        20        30        40        50        60        70        \nAAA                                                                             \n'
        g_var.seq_cg['PROTEIN'], g_var.seq_at['PROTEIN']= {0:'AAA'}, {0:'AAA'}
        to_print = gen.print_sequnce_info('PROTEIN')
        self.assertEqual(to_print,correct)

    def test_write_system_components(self):
        correct = '\n----------------------------------------------------------------------------------------------------\n                               Script has completed, time for a beer                                \n\nmolecules          number          \n---------          ------          \n PROTEIN              1            \n'
        g_var.system = {'PROTEIN':1}
        to_write = gen.write_system_components()
        self.assertEqual(to_write,correct)

    def test_database_information(self):
        output = '\n              The available forcefields within your database are (flag -ff):              \n------------------------------------------------------------------------------------------\n\n                                 charmm36-mar2019-updated                                 \n\n\n          The available fragment libraries within your database are (flag -fg):           \n------------------------------------------------------------------------------------------\n\n                                          test_1                                          \n\n"If all else fails, immortality can always be assured by spectacular error." (John Kenneth Galbraith)\n'
        g_var.forcefield_available =['charmm36-mar2019-updated']
        g_var.fragments_available = ['test_1']
        with self.assertRaises(SystemExit) as cm:
            gen.database_information()
        self.assertEqual(cm.exception.code, output)

    def test_fragments_in_use(self):
        output = '\n\n               The following residues are available in the database: test_1               \n------------------------------------------------------------------------------------------\n\n\n                                   Non protein residues                                   \n                                   --------------------                                   \n                                      CHOL, ION, SOL                                      \n\n                                     Protein residues                                     \n                                     ----------------                                     \n                                           PHE                                            \n\n                                Modified protein residues                                 \n                                -------------------------                                 \n                                                                                          \n\n                                  Other linked residues                                   \n                                  ---------------------                                   \n                                                                                          \n\n                                      Water residues                                      \n                                      --------------                                      \n                                 spc, spce, tip3p, tip4p                                  \n\n------------------------------------------------------------------------------------------\n\n'
        g_var.forcefield_available =['charmm36-mar2019-updated']
        g_var.fragments_available = ['test_1']
        g_var.args.fg=['test_1']
        to_print = gen.fragments_in_use('')
        self.assertEqual(to_print, output)

    def test_create_pdb(self):
        g_var.box_vec = 'test box vec '
        self.assertEqual(os.path.exists(run_dir+'database_test/test.pdb'), False)
        pdb_input = gen.create_pdb(run_dir+'database_test/test.pdb')
        pdb_input.close()
        self.assertEqual(os.path.exists(run_dir+'database_test/test.pdb'), True)
        self.assertTrue(filecmp.cmp(run_dir+'database_test/test.pdb', run_dir+'database_test/test_correct.pdb', shallow=False), msg='test pdb file header incorrect')
        os.remove(run_dir+'database_test/test.pdb')

    def test_mkdir_directory(self):
        self.assertEqual(os.path.exists(run_dir+'database_test/test_dir'), False)
        gen.mkdir_directory(run_dir+'database_test/test_dir')
        self.assertEqual(os.path.exists(run_dir+'database_test/test_dir'), True)
        os.system('rm -r '+run_dir+'database_test/test_dir')

    def test_clean(self):
        g_var.working_dir=run_dir+'database_test/'
        g_var.cg_residues= {'test_clean':1}
        for file_name in ['test_temp.pdb', 'test_tem.pdb', 'test_tem.tpr', 'test_temp.tpr']:
            self.assertEqual(os.path.exists(run_dir+'database_test/test_clean/'+file_name), True)
        gen.clean(True)
        result = [False, True, True, True]
        for file_val, file_name in enumerate(['test_temp.pdb', 'test_tem.pdb', 'test_tem.tpr', 'test_temp.tpr']):
            self.assertEqual(os.path.exists(run_dir+'database_test/test_clean/'+file_name), result[file_val])
            open(run_dir+'database_test/test_clean/'+file_name, 'w').close()

    def test_print_water_selection(self):
        g_var.args.w = 'test'
        water = ['tip3p', 'tip4p', 'spc', 'spce']
        directory = [run_dir+'database_test/fragments/test_1/non_protein/']
        correct  = "\nThe water type test doesn't exist\n\nPlease select a water molecule from below:\n\n     Selection              water_molecule        \n     ---------                ----------          \n         0                      tip3p             \n         1                      tip4p             \n         2                       spc              \n         3                       spce             \n"
        result = gen.print_water_selection(water, directory)
        self.assertEqual(result, correct)
        g_var.args.w = None
        with self.assertRaises(SystemExit) as cm:
            result = gen.print_water_selection([], directory)
        self.assertEqual(cm.exception.code, '\nCannot find any water models in: \n\n'+directory[0]+'SOL/'+'\n')

    def test_trunc_coord(self):
        xyz = [83.97299999999998, 44.467999999999996, 28.028000000000002]
        correct = [83.973, 44.468, 28.028]
        x, y, z = gen.trunc_coord(xyz)
        self.assertEqual([x, y, z], correct)

    def test_pdbatom(self):
        line = 'ATOM     31  N   ALA     3      74.248  59.378  51.435  1.00  0.00           N'
        line_sep = gen.pdbatom(line)
        self.assertEqual(line_sep, {'atom_number': 31, 'atom_name': 'N', 'residue_name': 'ALA', 'chain': ' ', 'residue_id': 3, 'x': 74.248, 'y': 59.378, 'z': 51.435})
        with self.assertRaises(SystemExit) as cm:
            line = 'ATOM     31  N   ALA     3      74.248  59.37833  51.435  1.00  0.00           N'
            line_sep = gen.pdbatom(line)
        self.assertEqual(cm.exception.code, '\npdb line is wrong:\t'+line)

    # def test_find_gromacs(self):
    #     pass

    # def test_correct_number_cpus(self):
    #     pass

    # def test_flags_used(self):
    #     pass

    # def test_print_script_timings(self):
    #     pass

    # def test_cg2at_header(self):
    #     pass

######## read_in file
    def test_add_residue_to_dictionary(self):
        g_var.p_residues, g_var.np_residues, g_var.o_residues = ['PHE'], ['CHOL', 'NA', 'W', 'ION'], ['DNA']
        g_var.cg_water_types = ['W', 'SOL', 'WN', 'WF', 'PW']
        line_test = [{'residue_name':'CHOL'}, {'residue_name':'PHE'}, {'residue_name':'ION'}, {'residue_name':'NA'}, {'residue_name':'W'}, {'residue_name':'DNA'}, {'residue_name':'SKIP'}]
        for l_val, l in enumerate(line_test):
            read_in.add_residue_to_dictionary(l)
        self.assertEqual(g_var.cg_residues, {'CHOL': {}, 'PROTEIN': {}, 'ION': {}, 'SOL': {}, 'NA': {}, 'OTHER': {}})
        with self.assertRaises(SystemExit) as cm:
            read_in.add_residue_to_dictionary({'residue_name':'AAA'})
        self.assertEqual(cm.exception.code, '\nAAA is not in the fragment database!')

    def test_add_to_cg_database(self):
        g_var.cg_residues = {'ION':{}, 'CHOL':{}, 'W':{}, 'SOL':{}, 'DNA':{}, 'PROTEIN':{}, 'OTHER':{}, 'NA':{}}
        g_var.p_residues, g_var.np_residues, g_var.o_residues = ['PHE'], ['CHOL', 'NA', 'W', 'ION'], ['DNA']
        g_var.water = 'tip3p'
        line_test = [{'residue_name':'CHOL', 'atom_name':'test_chol'}, {'residue_name':'PHE', 'atom_name':'test_phe'}, {'residue_name':'ION', 'atom_name':'test_ion'}, {'residue_name':'NA', 'atom_name':'test_na'}, {'residue_name':'W', 'atom_name':'test_w'}, {'residue_name':'DNA', 'atom_name':'test_dna'}]
        output = {'ION': {10: {'test_ion': {}}}, 'CHOL': {10: {'test_ion': {}}}, 'W': {10: {'test_ion': {}}}, 'SOL': {10: {'tip3p': {'residue_name': 'SOL'}}}, 'DNA': {}, 'PROTEIN': {10: {'test_ion': {}}}, 'OTHER': {10: {'test_ion': {}}}, 'NA': {10: {'test_ion': {}}}}
        for l_val, l in enumerate(line_test):
            read_in.add_to_cg_database(l, 10, {'test_ion':{}})   
        self.assertEqual(g_var.cg_residues, output)

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
        test_out = [['GL0', 'POPG'],['C1A', 'POPG'],['SKIP', 'POPG'],['NA+', 'SKIP'],['NA+', 'ION'] ]
        for test_val, test in enumerate(test_suite):
            g_var.swap_dict = swap[test_val]
            atom, resname = read_in.swap(test[0], test[1], test[2])
            self.assertEqual(atom, test_out[test_val][0])
            self.assertEqual(resname, test_out[test_val][1])

    def test_real_box_vectors(self):
        box_vec = 'CRYST1  159.804  124.407  103.403  90.00  90.00  90.00 P 1           1'
        correct1 = np.array([[1.59804000e+02, 7.61773172e-15, 6.33160765e-15], [0.00000000e+00, 1.24407000e+02, 6.33160765e-15], [0.00000000e+00, 0.00000000e+00, 1.03403000e+02]])
        correct2 = np.array([[ 6.25766564e-03,  0.00000000e+00,  0.00000000e+00], [-3.83171510e-19,  8.03813290e-03,  0.00000000e+00], [-3.83171510e-19, -4.92193687e-19,  9.67089930e-03]])
        r_b_vec, r_b_inv = read_in.real_box_vectors(box_vec)
        np.testing.assert_array_almost_equal(r_b_vec, correct1)
        np.testing.assert_array_almost_equal(r_b_inv, correct2)

    def test_brute_mic(self):
        r_b_vec = np.array([[1.59804000e+02, 7.61773172e-15, 6.33160765e-15], [0.00000000e+00, 1.24407000e+02, 6.33160765e-15], [0.00000000e+00, 0.00000000e+00, 1.03403000e+02]])
        p1, p2= np.array([155, 50, 50]), np.array([[5, 50, 50], [160, 50, 50], [50, 50, 50]])
        correct = [[164.804,  50.,     50.   ], [160,  50,  50], [50, 50, 50]]
        for i in range(3):
            result = read_in.brute_mic(p1, p2[i], r_b_vec)
            self.assertIsNone(np.testing.assert_array_equal(result, correct[i]))

    def test_duplicate_chain(self):
        g_var.args.d = ['0:2', '3:45']
        g_var.chain_count = 2
        g_var.atomistic_protein_input_raw = {0:1, 1:2}
        with self.assertRaises(SystemExit) as cm:
            read_in.duplicate_chain(True)
        self.assertEqual(g_var.chain_count, 3)
        self.assertEqual(g_var.atomistic_protein_input_raw, {0: 1, 1: 2, 2: 1})      
        self.assertEqual(cm.exception.code, 'your atomistic chain duplication input is incorrrect')

    def test_filter_input(self):
        input_no_atom = ['CRYST1  159.804  124.407  103.403  90.00  90.00  90.00 P 1           1\n', 'MODEL        1\n', '\n']
        with self.assertRaises(SystemExit) as no_atom:
            read_in.filter_input(input_no_atom)
        self.assertEqual(no_atom.exception.code, 'input coarsegrain structure seems to contain no beads')
        input_no_box =['MODEL        1\n', 'ATOM  17534  W     W   676      54.070   9.145   7.973  1.00  0.00           W\n', '\n']
        g_var.input_directory = run_dir+'/test_inputs/CG/no_box/'
        with self.assertRaises(SystemExit) as no_box:
            read_in.filter_input(input_no_box)
        self.assertEqual(no_box.exception.code, 'The input file is missing the Box vectors')

    def test_read_initial_cg_pdb(self):
        g_var.p_residues, g_var.np_residues = ['PHE'], ['CHOL', 'NA', 'W', 'ION']
        g_var.cg_water_types = ['W', 'SOL', 'WN', 'WF', 'PW']
        g_var.water = 'tip3p'
        g_var.cg_residues = {}
        cg_residues_correct = {'SOL': {0: {'tip3p': {'residue_name': 'SOL', 'coord': np.array([54.07 ,  9.145,  7.973])}}, 1: {'tip3p': {'residue_name': 'SOL', 'coord': np.array([108.503,  91.79 ,  90.375])}}, 2: {'tip3p': {'residue_name': 'SOL', 'coord': np.array([60.503, 61.79 , 60.375])}}, 3: {'tip3p': {'residue_name': 'SOL', 'coord': np.array([70.503, 71.79 , 70.375])}}}, 'ION': {1: {'CL-': {'residue_name': 'ION', 'coord': np.array([108.503,  91.79 ,  90.375])}}, 2: {'CL-': {'residue_name': 'ION', 'coord': np.array([60.503, 61.79 , 60.375])}}, 3: {'NA+': {'residue_name': 'ION', 'coord': np.array([70.503, 71.79 , 70.375])}}}}
        g_var.input_directory = run_dir+'/test_inputs/CG/'
        box_vec = read_in.read_initial_cg_pdb(True)
        self.assertEqual(box_vec, 'CRYST1  159.804  124.407  103.403  90.00  90.00  90.00 P 1           1\n')
        self.assertCountEqual(g_var.cg_residues, cg_residues_correct)

    def test_fix_pbc(self):
        g_var.cg_residues = {'PROTEIN': {0: {'BB': {'residue_name': 'ALA', 'coord': np.array([41.938, 58.822, 52.274])}}, 1: {'BB': {'residue_name': 'ALA', 'coord': np.array([75.016, 60.271, 50.624])}}, 2: {'BB': {'residue_name': 'ALA', 'coord': np.array([77.956, 62.112, 52.127])}}, 3: {'BB': {'residue_name': 'ALA', 'coord': np.array([ 1.118, 63.4  , 50.505])}}, 4: {'BB': {'residue_name': 'ALA', 'coord': np.array([ 3.974, 65.397, 51.969])}}}}
        box_vec, new_box, box_shift = 'CRYST1   80.000  124.000  103.000  90.00  90.00  90.00 P 1           1', 'CRYST1   80.000  124.000  103.000  90.00  90.00  90.00 P 1           1', [0, 0, 0]
        g_var.res_top['ALA']={'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, 'GROUPS': {'BB': 1}, 'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}, 'ATOMS': ['N', 'CA', 'CB', 'C', 'O'], 'RESIDUE': ['ALA'], 'atom_masses': {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'C': 12.011, 'O': 15.999}, 'amide_h': 'HN'}
        correct = {'PROTEIN': {0: {'BB': {'residue_name': 'ALA', 'coord': np.array([41.938, 58.822, 52.274])}}, 1: {'BB': {'residue_name': 'ALA', 'coord': np.array([75.016, 60.271, 50.624])}}, 2: {'BB': {'residue_name': 'ALA', 'coord': np.array([77.956, 62.112, 52.127])}}, 3: {'BB': {'residue_name': 'ALA', 'coord': np.array([81.118, 63.4  , 50.505])}}, 4: {'BB': {'residue_name': 'ALA', 'coord': np.array([83.974, 65.397, 51.969])}}}}
        read_in.fix_pbc(box_vec, new_box, box_shift)
        self.assertCountEqual(g_var.cg_residues, correct)

    def test_read_in_atomistic(self):
        g_var.box_vec = 'CRYST1  159.804  124.407  103.403  90.00  90.00  90.00 P 1           1\n'
        g_var.p_residues, g_var.np_residues = ['ALA'], ['CHOL', 'NA', 'W', 'ION']
        g_var.cg_water_types = ['W', 'SOL', 'WN', 'WF', 'PW']
        g_var.water = 'tip3p'
        g_var.res_top['ALA']={'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, 'GROUPS': {'BB': 1}, 'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}, 'ATOMS': ['N', 'CA', 'CB', 'C', 'O'], 'RESIDUE': ['ALA'], 'atom_masses': {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'C': 12.011, 'O': 15.999}, 'amide_h': 'HN'}
        atomistic_protein_input_correct_long = {0: {0: {1: {'coord': np.array([65.317, 54.293, 52.349]), 'atom': 'N', 'res_type': 'ALA', 'frag_mass': 14.007, 'resid': 0}, 2: {'coord': np.array([65.999, 55.581, 52.443]), 'atom': 'CA', 'res_type': 'ALA', 'frag_mass': 12.011, 'resid': 0}, 3: {'coord': np.array([67.496, 55.408, 52.345]), 'atom': 'C', 'res_type': 'ALA', 'frag_mass': 12.011, 'resid': 0}, 4: {'coord': np.array([68.097, 54.526, 52.969]), 'atom': 'O', 'res_type': 'ALA', 'frag_mass': 15.999, 'resid': 0}, 5: {'coord': np.array([65.561, 56.245, 53.76 ]), 'atom': 'CB', 'res_type': 'ALA', 'frag_mass': 12.011, 'resid': 0}}, 1: {11: {'coord': np.array([68.168, 56.207, 51.587]), 'atom': 'N', 'res_type': 'ALA', 'frag_mass': 14.007, 'resid': 1}, 12: {'coord': np.array([69.616, 56.04 , 51.492]), 'atom': 'CA', 'res_type': 'ALA', 'frag_mass': 12.011, 'resid': 1}, 13: {'coord': np.array([70.315, 57.378, 51.492]), 'atom': 'C', 'res_type': 'ALA', 'frag_mass': 12.011, 'resid': 1}, 14: {'coord': np.array([69.917, 58.326, 50.804]), 'atom': 'O', 'res_type': 'ALA', 'frag_mass': 15.999, 'resid': 1}, 15: {'coord': np.array([69.908, 55.215, 50.227]), 'atom': 'CB', 'res_type': 'ALA', 'frag_mass': 12.011, 'resid': 1}}}, 1: {3: {31: {'coord': np.array([74.248, 59.378, 51.435]), 'atom': 'N', 'res_type': 'ALA', 'frag_mass': 14.007, 'resid': 3}, 32: {'coord': np.array([75.698, 59.21 , 51.435]), 'atom': 'CA', 'res_type': 'ALA', 'frag_mass': 12.011, 'resid': 3}, 33: {'coord': np.array([76.397, 60.545, 51.336]), 'atom': 'C', 'res_type': 'ALA', 'frag_mass': 12.011, 'resid': 3}, 34: {'coord': np.array([76.038, 61.416, 50.536]), 'atom': 'O', 'res_type': 'ALA', 'frag_mass': 15.999, 'resid': 3}, 35: {'coord': np.array([76.058, 58.262, 50.278]), 'atom': 'CB', 'res_type': 'ALA', 'frag_mass': 12.011, 'resid': 3}}}}
        atomistic_protein_input, chain_count = read_in.read_in_atomistic(run_dir+'test_inputs/AT/AT_INPUT_long.pdb')
        self.assertEqual(chain_count, 2)
        self.assertCountEqual(atomistic_protein_input,atomistic_protein_input_correct_long)
        with self.assertRaises(SystemExit) as cm:
            atomistic_protein_input, chain_count = read_in.read_in_atomistic(run_dir+'test_inputs/AT/AT_INPUT_long.pdbhh')
        self.assertEqual(cm.exception.code, 'cannot find atomistic protein : '+run_dir+'test_inputs/AT/AT_INPUT_long.pdbhh')

#### test at_mod_p

    def test_shrink_coordinates(self):
        p1, p2 = np.array([58.274,66.912,12.038]), np.array([59.360,66.216,14.859])
        p1a, p2a = np.array([58.450475 , 66.7989   , 12.4964125]), np.array([59.183525 , 66.3291   , 14.4005875])
        p1s, p2s = at_mod_p.shrink_coordinates(p1, p2)
        self.assertIsNone(np.testing.assert_array_equal(p1s, p1a))
        self.assertIsNone(np.testing.assert_array_equal(p2s, p2a))

#### test at_mod_np
    
    # def test_build_atomistic_system(self):
        # at_mod_np.build_atomistic_system(residue_type, a)

#### test at_mod

    def test_sanity_check_fragments(self):
        g_var.res_top['PHE']={'C_TERMINAL': 'default', 'N_TERMINAL': 'default', \
        'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, \
        'GROUPS': {'BB': 2, 'SC1': 1, 'SC2': 1, 'SC3': 1}, \
        'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}, \
        'ATOMS': ['N', 'CA', 'C', 'O'], \
        'RESIDUE': ['PHE'], \
        'atom_masses': {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'CG': 12.011, 'CD1': 12.011, 'CE1': 12.011, 'CZ': 12.011, \
                        'CD2': 12.011, 'CE2': 12.011, 'C': 12.011, 'O': 15.999}, 'amide_h': 'HN'}
        res = 'PHE'
        cg = {'BB': {'residue_name': 'PHE', 'coord': np.array([84.312, 45.09 , 28.573])}, 'SC1': {'residue_name': 'PHE', 'coord': np.array([82.306, 43.106, 29.565])}, 'SC2': {'residue_name': 'PHE', 'coord': np.array([79.798, 42.108, 29.557])}, 'SC3': {'residue_name': 'PHE', 'coord': np.array([81.894, 40.487, 30.077])}}
        sin_bead = False
        bead_list, atom_list = at_mod.sanity_check_fragments(res, cg, sin_bead)
        self.assertEqual(bead_list, ['BB', 'SC1', 'SC2', 'SC3'])
        self.assertEqual(atom_list, [1, 2, 10, 11, 3, 4, 5, 8, 9, 6, 7])

    def test_get_atomistic(self):
        residue_correct = {2: {'BB': {1: {'coord': np.array([37.827, 15.084,  9.828]), 'atom': 'N', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 14.007}, 2: {'coord': np.array([38.493, 16.128, 10.269]), 'atom': 'CA', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 10: {'coord': np.array([39.816, 15.84 , 10.395]), 'atom': 'C', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 11: {'coord': np.array([40.176, 14.877, 10.872]), 'atom': 'O', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 15.999}}}, 1: {'SC1': {3: {'coord': np.array([37.998, 16.524, 11.52 ]), 'atom': 'CB', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 4: {'coord': np.array([36.657, 16.857, 11.574]), 'atom': 'CG', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 5: {'coord': np.array([35.802, 15.957, 11.799]), 'atom': 'CD1', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}}, 'SC2': {8: {'coord': np.array([36.27 , 18.027, 11.34 ]), 'atom': 'CD2', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 9: {'coord': np.array([35.046, 18.306, 11.331]), 'atom': 'CE2', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}}, 'SC3': {6: {'coord': np.array([34.578, 16.227, 11.781]), 'atom': 'CE1', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 7: {'coord': np.array([34.2  , 17.397, 11.547]), 'atom': 'CZ', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}}}}
        frag_mass_correct = {'BB': [[37.827000000000005, 15.084000000000001, 9.828, 14.007], [38.493, 16.128000000000004, 10.269, 12.011], [39.816, 15.840000000000002, 10.395000000000001, 12.011], [40.176, 14.877, 10.872, 15.999]], 'SC1': [[37.998, 16.524, 11.520000000000001, 12.011], [36.657, 16.857, 11.574, 12.011], [35.802, 15.957, 11.799, 12.011]], 'SC2': [[36.269999999999996, 18.027, 11.34, 12.011], [35.046, 18.306, 11.331, 12.011]], 'SC3': [[34.578, 16.227, 11.781, 12.011], [34.2, 17.397, 11.547, 12.011]]}
        residue, fragment_mass = at_mod.get_atomistic(run_dir+'database_test/fragments/test_1/protein/PHE/PHE.pdb')
        self.assertCountEqual(residue, residue_correct)
        self.assertCountEqual(fragment_mass, frag_mass_correct)


    def test_sanity_check_atoms(self):
        raised = False
        try:
            at_mod.sanity_check_atoms([1, 2, 3, 4, 5], 'ALA')
        except:
            raised = True
        self.assertFalse(raised, 'Exception raised')    

        with self.assertRaises(SystemExit) as cm:
            at_mod.sanity_check_atoms([1, 2, 3, 4, 6], 'ALA')
        self.assertEqual(cm.exception.code, 'atom number '+str(5)+' is missing from fragment library: ALA\n')

    def test_sanity_check_beads(self):
        bead_list = ['BB', 'SC1', 'SC2', 'SC3']
        cg = {'BB': {'residue_name': 'PHE', 'coord': np.array([84.312, 45.09 , 28.573])}, 'SC1': {'residue_name': 'PHE', 'coord': np.array([82.306, 43.106, 29.565])}, 'SC2': {'residue_name': 'PHE', 'coord': np.array([79.798, 42.108, 29.557])}, 'SC3': {'residue_name': 'PHE', 'coord': np.array([81.894, 40.487, 30.077])}}
        res = 'PHE' 
        bead_list_cg = at_mod.sanity_check_beads(bead_list, cg, res)
        self.assertEqual(bead_list_cg, ['BB', 'SC1', 'SC2', 'SC3'])
        cgf = {'FF': {'residue_name': 'PHE', 'coord': np.array([84.312, 45.09 , 28.573])}, 'SC1': {'residue_name': 'PHE', 'coord': np.array([82.306, 43.106, 29.565])}, 'SC2': {'residue_name': 'PHE', 'coord': np.array([79.798, 42.108, 29.557])}, 'SC3': {'residue_name': 'PHE', 'coord': np.array([81.894, 40.487, 30.077])}}
        with self.assertRaises(SystemExit) as cm:
            bead_list_cg = at_mod.sanity_check_beads(bead_list, cgf, res)
        self.assertEqual(cm.exception.code, 'The bead FF is missing from the fragment library: PHE\n')

    def test_sanity_check_protein_other(self):
        g_var.res_top['PHE']={'C_TERMINAL': 'default', 'N_TERMINAL': 'default', \
        'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, \
        'GROUPS': {'BB': 2, 'SC1': 1, 'SC2': 1, 'SC3': 1}, \
        'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}, \
        'ATOMS': ['N', 'CA', 'C', 'O'], \
        'RESIDUE': ['PHE'], \
        'atom_masses': {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'CG': 12.011, 'CD1': 12.011, 'CE1': 12.011, 'CZ': 12.011, \
                        'CD2': 12.011, 'CE2': 12.011, 'C': 12.011, 'O': 15.999}, 'amide_h': 'HN'}
        g_var.cg_residues = {'PROTEIN': {0: {'BB': {'residue_name': 'PHE', 'coord': np.array([84.312, 45.09 , 28.573])}, 'SC1': {'residue_name': 'PHE', 'coord': np.array([82.306, 43.106, 29.565])}, 'SC2': {'residue_name': 'PHE', 'coord': np.array([79.798, 42.108, 29.557])}, 'SC3': {'residue_name': 'PHE', 'coord': np.array([81.894, 40.487, 30.077])}}}}        
        raised = False
        try:
            at_mod.sanity_check_protein_other('PROTEIN')
        except:
            raised = True
        self.assertFalse(raised, 'Exception raised')    

    def test_sanity_check_protein_other_wrong_bead(self):
        g_var.res_top['PHE']={'C_TERMINAL': 'default', 'N_TERMINAL': 'default', \
        'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, \
        'GROUPS': {'BB': 2, 'SC1': 1, 'SC2': 1, 'SC3': 1}, \
        'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}, \
        'ATOMS': ['N', 'CA', 'C', 'O'], \
        'RESIDUE': ['PHE'], \
        'atom_masses': {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'CG': 12.011, 'CD1': 12.011, 'CE1': 12.011, 'CZ': 12.011, \
                        'CD2': 12.011, 'CE2': 12.011, 'C': 12.011, 'O': 15.999}, 'amide_h': 'HN'}
        g_var.cg_residues = {'PROTEIN': {0: { 'SC1': {'residue_name': 'PHE', 'coord': np.array([82.306, 43.106, 29.565])}, 'SC2': {'residue_name': 'PHE', 'coord': np.array([79.798, 42.108, 29.557])}, 'SC3': {'residue_name': 'PHE', 'coord': np.array([81.894, 40.487, 30.077])}}}}        
        with self.assertRaises(SystemExit) as cm:
            at_mod.sanity_check_protein_other('PROTEIN', True)
        self.assertEqual(cm.exception.code, 'number of atomistic fragments: 4 does not equal number of CG beads: 3')    
        
    def test_sanity_check(self):
        g_var.res_top['PHE']={'C_TERMINAL': 'default', 'N_TERMINAL': 'default', \
        'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, \
        'GROUPS': {'BB': 2, 'SC1': 1, 'SC2': 1, 'SC3': 1}, \
        'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}, \
        'ATOMS': ['N', 'CA', 'C', 'O'], \
        'RESIDUE': ['PHE'], \
        'atom_masses': {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'CG': 12.011, 'CD1': 12.011, 'CE1': 12.011, 'CZ': 12.011, \
                        'CD2': 12.011, 'CE2': 12.011, 'C': 12.011, 'O': 15.999}, 'amide_h': 'HN'}
        g_var.cg_residues = {'PROTEIN': {0: {'BB': {'residue_name': 'PHE', 'coord': np.array([84.312, 45.09 , 28.573])}, 'SC1': {'residue_name': 'PHE', 'coord': np.array([82.306, 43.106, 29.565])}, 'SC2': {'residue_name': 'PHE', 'coord': np.array([79.798, 42.108, 29.557])}, 'SC3': {'residue_name': 'PHE', 'coord': np.array([81.894, 40.487, 30.077])}}}}        
        raised = False
        try:
            at_mod.sanity_check()
        except:
            raised = True
        self.assertFalse(raised, 'Exception raised')         

    def test_sanity_check_solvent(self):
        g_var.res_top['SOL'] = {'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': []}, 'GROUPS': {'tip3p': 1, 'tip4p': 2, 'spc': 3, 'spce': 4}, 'CONNECT': {'atoms': {}}, 'RESIDUE': ['SOL'], 'atom_masses': {'OW': 15.9994, 'HW1': 1.008, 'HW2': 1.008, 'MW': 0.0}}
        g_var.cg_residues = {'SOL': {0: {'tip3p': {'residue_name': 'SOL', 'coord': np.array([60.577, 12.72 ,  2.4  ])}}}}
        raised = False
        try:
            at_mod.sanity_check_solvent('SOL')
        except:
            raised = True
        self.assertFalse(raised, 'Exception raised')         

    def test_sanity_check_non_protein(self):
        g_var.res_top['CHOL'] = {'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': []}, 'GROUPS': {'ROH': 1, 'R1': 1, 'R2': 1, 'R3': 2, 'R4': 2, 'R5': 2, 'C1': 3, 'C2': 4}, 'CONNECT': {'atoms': {}}, 'RESIDUE': ['CHOL'], 'atom_masses': {'C3': 12.011, 'O3': 15.9994, 'C4': 12.011, 'C5': 12.011, 'C6': 12.011, 'C7': 12.011, 'C8': 12.011, 'C14': 12.011, 'C15': 12.011, 'C16': 12.011, 'C17': 12.011, 'C13': 12.011, 'C18': 12.011, 'C12': 12.011, 'C11': 12.011, 'C9': 12.011, 'C10': 12.011, 'C19': 12.011, 'C1': 12.011, 'C2': 12.011, 'C20': 12.011, 'C21': 12.011, 'C22': 12.011, 'C23': 12.011, 'C24': 12.011, 'C25': 12.011, 'C26': 12.011, 'C27': 12.011}}
        g_var.cg_residues = {'CHOL': {0: {'ROH': {'residue_name': 'CHOL', 'coord': np.array([60.577, 12.72 ,  2.4  ])}, 'R1': {'residue_name': 'CHOL', 'coord': np.array([60.577, 12.72 ,  2.4  ])}, 'R2': {'residue_name': 'CHOL', 'coord': np.array([60.577, 12.72 ,  2.4  ])}, 'R3': {'residue_name': 'CHOL', 'coord': np.array([60.577, 12.72 ,  2.4  ])}, 'R4': {'residue_name': 'CHOL', 'coord': np.array([60.577, 12.72 ,  2.4  ])}, 'R5': {'residue_name': 'CHOL', 'coord': np.array([60.577, 12.72 ,  2.4  ])}, 'C1': {'residue_name': 'CHOL', 'coord': np.array([60.577, 12.72 ,  2.4  ])}, 'C2': {'residue_name': 'CHOL', 'coord': np.array([60.577, 12.72 ,  2.4  ])}}}}
        raised = False
        try:
            at_mod.sanity_check_non_protein('CHOL')
        except:
            raised = True
        self.assertFalse(raised, 'Exception raised') 

    def test_fix_atom_wrap(self):
        correct = {'PROTEIN': {0: {'BB': {'residue_name': 'PHE', 'coord': np.array([84.312, 45.09 , 28.573])}, 'SC1': {'residue_name': 'PHE', 'coord': np.array([82.306, 43.106, 29.565])}, 'SC2': {'residue_name': 'PHE', 'coord': np.array([79.798, 42.108, 29.557])}, 'SC3': {'residue_name': 'PHE', 'coord': np.array([81.894, 40.487, 30.077])}}}}        

        g_var.cg_residues = {'PROTEIN': {0: {'BB': {'residue_name': 'PHE', 'coord': np.array([84.312, 45.09 , 28.573])}, '1SC': {'residue_name': 'PHE', 'coord': np.array([82.306, 43.106, 29.565])}, 'SC2': {'residue_name': 'PHE', 'coord': np.array([79.798, 42.108, 29.557])}, 'SC3': {'residue_name': 'PHE', 'coord': np.array([81.894, 40.487, 30.077])}}}}        

        bead_list, bead_list_cg, res_type, residue= ['BB', 'SC1', 'SC2', 'SC3'], ['BB', '1SC', 'SC2', 'SC3'], 'PROTEIN', 0
        at_mod.fix_atom_wrap(bead_list, bead_list_cg, res_type, residue)
        self.assertCountEqual(g_var.cg_residues, correct)

    def test_rotate_atom(self):
        coord_test, center_test,xyz_rot_apply_test = np.array([78.09324928, 73.84975034, 66.00837494]), np.array([77.01,  75.365, 64.802]), np.array([[-0.13784372,  0.95043581,  0.27869495], [-0.92459301, -0.02258486, -0.38028633], [-0.35514346, -0.31009948,  0.88187949]])
        correct = [77.83323396, 76.05468437, 66.74400215]
        coord = at_mod.rotate_atom(coord_test, center_test,xyz_rot_apply_test)
        np.testing.assert_array_almost_equal(coord,correct)

    def test_kabsch_rotate(self):
        at, cg = np.array([[-0.82883763,  0.86100236,  0.16477411], [ 0.54816237, -0.86699764, -0.57322589]]), np.array([[ 4.109,   0.132,   2.326 ], [ 0.813,   1.628,  -0.0295]])
        correct = np.array([[-0.40049698, -0.30353545, -0.86456254], [ 0.79299884, -0.58754593, -0.16106715], [-0.45908061, -0.750104,    0.47601363]])
        result = at_mod.kabsch_rotate(at,cg)
        np.testing.assert_array_almost_equal(result, correct)

    def test_find_cross_vector(self):
        ca = [np.array([66.297, 54.97 , 52.774]), np.array([69.286, 56.581, 51.125]), np.array([72.312, 58.261, 52.646])]
        correct = [ 0.48016321, -0.87707549,  0.01348649]
        result = at_mod.find_cross_vector(ca)
        np.testing.assert_array_almost_equal(result, correct)

    def test_noramlised_vector(self):
        c1, c2 = np.array([94.5448748,  70.14307068, 51.5138988 ]), np.array([94.52699527, 69.12609437, 51.06636723])
        correct = np.array([0.01608978, 0.91517598, 0.40273322])
        result = at_mod.noramlised_vector(c1, c2)
        np.testing.assert_array_almost_equal(result, correct)

    def test_align_to_vector(self):
        v1,v2 = np.array([0.01608978, 0.91517598, 0.40273322]), np.array([-0.4341557,   0.81322425, -0.38752439])
        correct = np.array([[ 0.87549255,  0.48315775, -0.00844915], [-0.33766917,  0.5991694,  -0.72593082], [-0.34567663,  0.63840004,  0.68771582]])
        result = at_mod.align_to_vector(v1, v2)
        np.testing.assert_array_almost_equal(result, correct)


    def test_align_at_frag_to_CG_frag(self):
        at, cg, group = np.array([45.0827361,  45.18484289, 45.58242755]), np.array([93.604, 69.2,   50.814]), {'BB': {1: {'coord': np.array([43.8669, 44.4825, 45.6372]), 'atom': 'N', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 14.007}, 2: {'coord': np.array([45.1359, 44.5455, 45.2862]), 'atom': 'CA', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 3: {'coord': np.array([45.1899, 45.0225, 43.9902]), 'atom': 'CB', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 4: {'coord': np.array([45.8379, 45.3195, 46.1592]), 'atom': 'C', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 5: {'coord': np.array([45.4599, 46.3005, 46.5192]), 'atom': 'O', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 15.999}}}
        correct1 = {'BB': {1: {'coord': np.array([92.3881639 , 68.49765711, 50.86877245]), 'atom': 'N', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 14.007}, 2: {'coord': np.array([93.6571639 , 68.56065711, 50.51777245]), 'atom': 'CA', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 3: {'coord': np.array([93.7111639 , 69.03765711, 49.22177245]), 'atom': 'CB', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 4: {'coord': np.array([94.3591639 , 69.33465711, 51.39077245]), 'atom': 'C', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 5: {'coord': np.array([93.9811639 , 70.31565711, 51.75077245]), 'atom': 'O', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 15.999}}} 
        correct2 = np.array([-48.5212639,  -24.01515711,  -5.23157245])
        result1, result2 = at_mod.align_at_frag_to_CG_frag(at, cg, group)
        self.assertCountEqual(result1,correct1)
        np.testing.assert_array_almost_equal(result2, correct2)

    def test_COM(self):
        mass = np.array([[43.8669, 44.4825, 45.6372, 14.007], [45.13590000000001, 44.5455, 45.2862, 12.011], [45.1899, 45.0225, 43.9902, 12.011], [45.8379, 45.3195, 46.1592, 12.011], [45.459900000000005, 46.3005, 46.519200000000005, 15.999]] )
        fragment = {'BB': {1: {'coord': np.array([43.8669, 44.4825, 45.6372]), 'atom': 'N', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 14.007}, 2: {'coord': np.array([45.1359, 44.5455, 45.2862]), 'atom': 'CA', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 3: {'coord': np.array([45.1899, 45.0225, 43.9902]), 'atom': 'CB', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 4: {'coord': np.array([45.8379, 45.3195, 46.1592]), 'atom': 'C', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 5: {'coord': np.array([45.4599, 46.3005, 46.5192]), 'atom': 'O', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 15.999}}}
        com_correct = np.array([45.0827361,  45.18484289, 45.58242755])
        result = at_mod.COM(mass, fragment)
        np.testing.assert_array_almost_equal(result,com_correct)

    def test_COM_error(self):
        mass = np.array([[43.8669, 44.4825, 45.6372, 14.007], [45.13590000000001, 44.5455, 45.2862]])
        fragment = {'BB': {1: {'coord': np.array([43.8669, 44.4825, 45.6372]), 'atom': 'N', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 14.007}, 2: {'coord': np.array([45.1359, 44.5455, 45.2862]), 'atom': 'CA', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}}}
        with self.assertRaises(SystemExit) as cm:
            result = at_mod.COM(mass, fragment)
        self.assertEqual(cm.exception.code,  'missing the mass one of the atoms in ALA')  

    def test_rigid_fit(self):
        group = {'BB': {1: {'coord': np.array([43.8669, 44.4825, 45.6372]), 'atom': 'N', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 14.007}, 2: {'coord': np.array([45.1359, 44.5455, 45.2862]), 'atom': 'CA', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 3: {'coord': np.array([45.1899, 45.0225, 43.9902]), 'atom': 'CB', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 4: {'coord': np.array([45.8379, 45.3195, 46.1592]), 'atom': 'C', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 5: {'coord': np.array([45.4599, 46.3005, 46.5192]), 'atom': 'O', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 15.999}}} 
        frag_mass = {'BB': [[43.8669, 44.4825, 45.6372, 14.007], [45.13590000000001, 44.5455, 45.2862, 12.011], [45.1899, 45.0225, 43.9902, 12.011], [45.8379, 45.3195, 46.1592, 12.011], [45.459900000000005, 46.3005, 46.519200000000005, 15.999]]} 
        resid, cg = 9, {'BB': {'residue_name': 'ALA', 'coord': np.array([93.604, 69.2  , 50.814])}}

        rigid_mass_cg_cor, at_frag_centers_cor, cg_frag_centers_cor= np.array([93.604, 69.2,  50.814]), {'BB': np.array([93.604, 69.2  , 50.814])}, {'BB': np.array([93.604, 69.2  , 50.814])}
        group_cor={'BB': {1: {'coord': np.array([92.3881639 , 68.49765711, 50.86877245]), 'atom': 'N', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 14.007}, 2: {'coord': np.array([93.6571639 , 68.56065711, 50.51777245]), 'atom': 'CA', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 3: {'coord': np.array([93.7111639 , 69.03765711, 49.22177245]), 'atom': 'CB', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 4: {'coord': np.array([94.3591639 , 69.33465711, 51.39077245]), 'atom': 'C', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 5: {'coord': np.array([93.9811639 , 70.31565711, 51.75077245]), 'atom': 'O', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 15.999}}}

        rigid_mass_cg, at_frag_centers, cg_frag_centers, group = at_mod.rigid_fit(group, frag_mass, resid, cg)
        np.testing.assert_array_almost_equal(rigid_mass_cg, rigid_mass_cg_cor)
        self.assertCountEqual(at_frag_centers, at_frag_centers_cor)
        self.assertCountEqual(cg_frag_centers, cg_frag_centers_cor)
        self.assertCountEqual(group, group_cor)

    def test_overlapping_atoms(self):
        coordinates = [[0, 1, 0 ], [1,2,1], [0,1.1,0]]
        correct = [[0, 2]]
        tree = cKDTree(coordinates)
        overlapped = at_mod.overlapping_atoms(tree)
        self.assertEqual(overlapped, correct)

    def test_check_atom_overlap(self):
        coordinates = [[0, 1, 0 ], [1,2,1], [0,1.1,0]]
        result = at_mod.check_atom_overlap(coordinates)
        tree = cKDTree(result)
        overlapped = at_mod.overlapping_atoms(tree)
        self.assertEqual(len(overlapped), 0)

    def test_split_fragment_names(self):
        res_cor, group_cor, bead_cor = {1: {'BB': {}}}, 1, 'BB'
        line, residue, resname = '[ BB ]', {}, 'ALA'
        residue, group, bead = at_mod.split_fragment_names(line, residue, resname)
        self.assertCountEqual(residue,res_cor) 
        self.assertEqual(group, group_cor)
        self.assertEqual(bead, bead_cor )

    def test_get_atomistic(self):
        g_var.res_top['PHE'] = {'C_TERMINAL': 'default', 'N_TERMINAL': 'default', \
        'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, \
        'GROUPS': {'BB': 2, 'SC1': 1, 'SC2': 1, 'SC3': 1}, \
        'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}, \
        'ATOMS': ['N', 'CA', 'C', 'O'], \
        'RESIDUE': ['PHE'], \
        'atom_masses': {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'CG': 12.011, 'CD1': 12.011, 'CE1': 12.011, 'CZ': 12.011, \
                        'CD2': 12.011, 'CE2': 12.011, 'C': 12.011, 'O': 15.999}, 'amide_h': 'HN'}

        residue_cor = {2: {'BB': {1: {'coord': np.array([37.827, 15.084,  9.828]), 'atom': 'N', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 14.007}, 2: {'coord': np.array([38.493, 16.128, 10.269]), 'atom': 'CA', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 10: {'coord': np.array([39.816, 15.84 , 10.395]), 'atom': 'C', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 11: {'coord': np.array([40.176, 14.877, 10.872]), 'atom': 'O', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 15.999}}}, 1: {'SC1': {3: {'coord': np.array([37.998, 16.524, 11.52 ]), 'atom': 'CB', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 4: {'coord': np.array([36.657, 16.857, 11.574]), 'atom': 'CG', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 5: {'coord': np.array([35.802, 15.957, 11.799]), 'atom': 'CD1', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}}, 'SC2': {8: {'coord': np.array([36.27 , 18.027, 11.34 ]), 'atom': 'CD2', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 9: {'coord': np.array([35.046, 18.306, 11.331]), 'atom': 'CE2', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}}, 'SC3': {6: {'coord': np.array([34.578, 16.227, 11.781]), 'atom': 'CE1', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 7: {'coord': np.array([34.2  , 17.397, 11.547]), 'atom': 'CZ', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}}}}
        fragment_cor = {'BB': [[37.827000000000005, 15.084000000000001, 9.828, 14.007], [38.493, 16.128000000000004, 10.269, 12.011], [39.816, 15.840000000000002, 10.395000000000001, 12.011], [40.176, 14.877, 10.872, 15.999]], 'SC1': [[37.998, 16.524, 11.520000000000001, 12.011], [36.657, 16.857, 11.574, 12.011], [35.802, 15.957, 11.799, 12.011]], 'SC2': [[36.269999999999996, 18.027, 11.34, 12.011], [35.046, 18.306, 11.331, 12.011]], 'SC3': [[34.578, 16.227, 11.781, 12.011], [34.2, 17.397, 11.547, 12.011]]}
        residue, fragment_mass = at_mod.get_atomistic(run_dir+'database_test/fragments/test_1/protein/PHE/PHE.pdb')
        self.assertCountEqual(residue, residue_cor)
        self.assertCountEqual(fragment_mass, fragment_cor)

    def test_connectivity(self):
        g_var.sorted_connect['PHE'] = {2: {2: ['SC1']}, 1: {3: ['BB']}}
        cg = {'BB': {'residue_name': 'PHE', 'coord':  np.array([84.312, 45.09 , 28.573])}, 
            'SC1': {'residue_name': 'PHE', 'coord':  np.array([82.306, 43.106, 29.565])}, 
            'SC2': {'residue_name': 'PHE', 'coord':  np.array([79.798, 42.108, 29.557])}, 
            'SC3': {'residue_name': 'PHE', 'coord':  np.array([81.894, 40.487, 30.077])}}
        at_frag_centers = {'SC1':  np.array([82.35866667, 41.30419048, 29.808]), 'SC2':  np.array([81.19766667, 43.02469048, 29.5125]), 'SC3':  np.array([79.92866667, 41.67019048, 29.841])} 
        cg_frag_centers = {'SC1':  np.array([82.306, 43.106, 29.565]), 'SC2':  np.array([79.798, 42.108, 29.557]), 'SC3':  np.array([81.894, 40.487, 30.077])}
        group = {'SC1': {3: {'coord':  np.array([83.53766667, 41.38219048, 29.697     ]), 'atom': 'CB', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011},
                         4: {'coord':  np.array([82.19666667, 41.71519048, 29.751     ]), 'atom': 'CG', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 
                        5: {'coord':  np.array([81.34166667, 40.81519048, 29.976     ]), 'atom': 'CD1', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}}, 
                'SC2': {8: {'coord':  np.array([81.80966667, 42.88519048, 29.517     ]), 'atom': 'CD2', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 
                        9: {'coord':  np.array([80.58566667, 43.16419048, 29.508     ]), 'atom': 'CE2', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}}, 
                'SC3': {6: {'coord':  np.array([80.11766667, 41.08519048, 29.958     ]), 'atom': 'CE1', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}, 
                        7: {'coord':  np.array([79.73966667, 42.25519048, 29.724     ]), 'atom': 'CZ', 'resid': 1, 'res_type': 'PHE', 'frag_mass': 12.011}}}
        group_number = 1
        at_connection_cor =   np.array([  np.array([83.53766667, 41.38219048, 29.697 ]),  np.array([82.35866667, 41.30419048, 29.808 ]),  np.array([81.19766667, 43.02469048, 29.5125]),  np.array([79.92866667, 41.67019048, 29.841     ])])
        cg_connection_cor =   np.array([  np.array([84.312, 45.09 , 28.573]),  np.array([82.306, 43.106, 29.565]),  np.array([79.798, 42.108, 29.557]),  np.array([81.894, 40.487, 30.077])])
        at_connection, cg_connection = at_mod.connectivity(cg, at_frag_centers, cg_frag_centers, group, group_number)
        np.testing.assert_array_almost_equal(at_connection, at_connection_cor)
        np.testing.assert_array_almost_equal(cg_connection, cg_connection_cor)

    def test_BB_connectivity_1st_res(self):
        g_var.res_top['ALA']={'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, 'GROUPS': {'BB': 1}, 'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}, 'ATOMS': ['N', 'CA', 'CB', 'C', 'O'], 'RESIDUE': ['ALA'], 'atom_masses': {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'C': 12.011, 'O': 15.999}, 'amide_h': 'HN'}

        at_connections,cg_connections = [], []
        cg_residues =  {0: {'BB': {'residue_name': 'ALA', 'coord': np.array([66.297, 54.97 , 52.774])}}, 
                        1: {'BB': {'residue_name': 'ALA', 'coord': np.array([69.286, 56.581, 51.125])}}} 
        at_residues = {1: {'coord': np.array([65.0811639 , 54.26765711, 52.82877245]), 'atom': 'N', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 14.007}, 
                       2: {'coord': np.array([66.3501639 , 54.33065711, 52.47777245]), 'atom': 'CA', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 
                       3: {'coord': np.array([66.4041639 , 54.80765711, 51.18177245]), 'atom': 'CB', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 
                       4: {'coord': np.array([67.0521639 , 55.10465711, 53.35077245]), 'atom': 'C', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 
                       5: {'coord': np.array([66.6741639 , 56.08565711, 53.71077245]), 'atom': 'O', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 15.999}} 
        residue_number, BB_bead = 0, 'BB'
        at_connections_cor,cg_connections_cor = np.array([np.array([67.0521639 , 55.10465711, 53.35077245])]), np.array([np.array([69.286, 56.581, 51.125])])
        at_connections,cg_connections, new_chain = at_mod.BB_connectivity(at_connections,cg_connections, cg_residues, at_residues, residue_number, BB_bead)
        np.testing.assert_array_almost_equal(at_connections,at_connections_cor)
        np.testing.assert_array_almost_equal(cg_connections,cg_connections_cor)
        self.assertFalse(new_chain, 'Exception raised')

    def test_BB_connectivity_middle_res(self):
        g_var.res_top['ALA']={'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, 'GROUPS': {'BB': 1}, 'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}, 'ATOMS': ['N', 'CA', 'CB', 'C', 'O'], 'RESIDUE': ['ALA'], 'atom_masses': {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'C': 12.011, 'O': 15.999}, 'amide_h': 'HN'}

        at_connections,cg_connections = [], []
        cg_residues =  {0: {'BB': {'residue_name': 'ALA', 'coord': np.array([72.878, 57.293, 50.09 ])}}, 
                        1: {'BB': {'residue_name': 'ALA', 'coord': np.array([75.867, 58.904, 48.441])}}, 
                        2: {'BB': {'residue_name': 'ALA', 'coord': np.array([78.893, 60.584, 49.962])}}, 
                        3: {'BB': {'residue_name': 'ALA', 'coord': np.array([91.971, 72.032, 58.312])}}}
        at_residues = {1: {'coord': np.array([74.6511639 , 58.20165711, 48.49577245]), 'atom': 'N', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 14.007}, 
                       2: {'coord': np.array([75.9201639 , 58.26465711, 48.14477245]), 'atom': 'CA', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 
                       3: {'coord': np.array([75.9741639 , 58.74165711, 46.84877245]), 'atom': 'CB', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 
                       4: {'coord': np.array([76.6221639 , 59.03865711, 49.01777245]), 'atom': 'C', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 
                       5: {'coord': np.array([76.2441639 , 60.01965711, 49.37777245]), 'atom': 'O', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 15.999}} 
        residue_number, BB_bead = 1, 'BB'
        at_connections_cor,cg_connections_cor = np.array([ np.array([74.6511639 , 58.20165711, 48.49577245]), np.array([76.6221639 , 59.03865711, 49.01777245])]), np.array([ np.array([72.878, 57.293, 50.09 ]), np.array([78.893, 60.584, 49.962])])
        at_connections,cg_connections, new_chain = at_mod.BB_connectivity(at_connections,cg_connections, cg_residues, at_residues, residue_number, BB_bead)
        np.testing.assert_array_almost_equal(at_connections,at_connections_cor)
        np.testing.assert_array_almost_equal(cg_connections,cg_connections_cor)
        self.assertFalse(new_chain, 'Exception raised') 

    def test_BB_connectivity_end_res(self):
        g_var.res_top['ALA']={'C_TERMINAL': 'default', 'N_TERMINAL': 'default', 'CHIRAL': {'atoms': ['CA', 'HA', 'CB', 'N', 'C'], 'CA': {'m': 'HA', 'c1': 'CB', 'c2': 'N', 'c3': 'C'}}, 'GROUPS': {'BB': 1}, 'CONNECT': {'atoms': {'N': -1, 'C': 1}, 'BB': {'atom': ['N', 'C'], 'Con_Bd': ['BB', 'BB'], 'dir': [-1, 1]}}, 'ATOMS': ['N', 'CA', 'CB', 'C', 'O'], 'RESIDUE': ['ALA'], 'atom_masses': {'N': 14.007, 'CA': 12.011, 'CB': 12.011, 'C': 12.011, 'O': 15.999}, 'amide_h': 'HN'}

        at_connections,cg_connections = [], []
        cg_residues =  {0: {'BB': {'residue_name': 'ALA', 'coord': np.array([72.878, 57.293, 50.09 ])}}, 
                        1: {'BB': {'residue_name': 'ALA', 'coord': np.array([75.867, 58.904, 48.441])}}, 
                        2: {'BB': {'residue_name': 'ALA', 'coord': np.array([78.893, 60.584, 49.962])}}, 
                        3: {'BB': {'residue_name': 'ALA', 'coord': np.array([91.971, 72.032, 58.312])}}}
        at_residues = {1: {'coord': np.array([77.6771639 , 59.88165711, 50.01677245]), 'atom': 'N', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 14.007}, 
                       2: {'coord': np.array([78.9461639 , 59.94465711, 49.66577245]), 'atom': 'CA', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 
                       3: {'coord': np.array([79.0001639 , 60.42165711, 48.36977245]), 'atom': 'CB', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 
                       4: {'coord': np.array([79.6481639 , 60.71865711, 50.53877245]), 'atom': 'C', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 12.011}, 
                       5: {'coord': np.array([79.2701639 , 61.69965711, 50.89877245]), 'atom': 'O', 'resid': 1, 'res_type': 'ALA', 'frag_mass': 15.999}}
        residue_number, BB_bead = 2, 'BB'
        at_connections_cor,cg_connections_cor = np.array([np.array([77.6771639 , 59.88165711, 50.01677245])]), np.array([np.array([75.867, 58.904, 48.441])])
        at_connections,cg_connections, new_chain = at_mod.BB_connectivity(at_connections,cg_connections, cg_residues, at_residues, residue_number, BB_bead)
        np.testing.assert_array_almost_equal(at_connections,at_connections_cor)
        np.testing.assert_array_almost_equal(cg_connections,cg_connections_cor)
        self.assertTrue(new_chain, 'Exception raised') 



    # def test_merge_indivdual_chain_pdbs(self):
    #     pass

    def test_index_conversion_generate(self):
        merge = [{'atom_name':'AB'}, {'atom_name':'ABM'},{'atom_name':'MAB'}, {'atom_name':'1AB'},{'atom_name':'MAB'}]
        merge_coords = [[0, 1, 2],[1, 1, 2],[2, 1, 2],[3, 1, 2],[4, 1, 2],[5, 1, 2]]
        coords_cor = [[0, 1, 2], [1, 1, 2], [3, 1, 2]]
        index_cor = {0: 0, 1: 1, 3: 2}
        coords, index_conversion = at_mod.index_conversion_generate(merge, merge_coords)
        self.assertEqual(coords, coords_cor)
        self.assertEqual(index_conversion, index_cor)
    #     pass

    # def test_write_pdb(self):
    #     pass

    # def test_merge_system_pdbs(self):
    #     pass

    # def test_read_in_merged_pdbs(self):
    #     pass

    # def test_check_overlap_chain(self):
    #     pass

    # def test_fetch_chiral_coord(self):
    #     pass

    # def test_fix_chirality(self):
    #     pass

    # def test_check_hydrogens(self):
    #     pass

    # def test_check_ringed_lipids(self):
    #     pass

    # def test_fetch_start_of_residue(self):
    #     pass

    # def test_get_np_resname(self):
    #     pass

    # def test_fix_threaded_lipids(self):
    #     pass

    # def test_(self):
    #     pass

    # def test_(self):
    #     pass

    # def test_(self):
    #     pass

    # def test_(self):
    #     pass

    # def test_(self):
    #     pass

    # def test_(self):
    #     pass



if __name__ == '__main__':
    unittest.main()
