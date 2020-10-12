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

    # @patch('builtins.print')
    # def test_print_swap_residues(self, mock_print):
    #     g_var.args.swap = {'NA+': {'NA+:skip': {'ALL': 'ALL', 'resid': [4000, 4001, 4002], 'range': ['4000-4002']}}} 
    #     gen.print_swap_residues()       
    #     # sys.stdout.write(str( mock_print.call_args_list ) + '\n')
    #     print(str( mock_print.getvalue() ) + '\n') 
    #     print(len(mock_print.call_args_list))
    #     # mock_print.assert_called_with(

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
        self.assertEqual(atom_conversion,{'N': 1, 'CA': 2, 'C': 10, 'O': 11, 'CB': 3, 'CG': 4, 'CD1': 5, 'CD2': 8, 'CE2': 9, 'CE1': 6, 'CZ': 7})

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

######## read_in file
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
        test_out = [['GL0', 'POPG'],['C1A', 'POPG'],['SKIP', 'POPG'],['ION', 'SKIP'],['ION', 'ION'] ]
        for test_val, test in enumerate(test_suite):
            g_var.swap_dict = swap[test_val]
            atom, resname = read_in.swap(test[0], test[1], test[2])
            self.assertEqual(atom, test_out[test_val][0])
            self.assertEqual(resname, test_out[test_val][1])



    def test_add_to_cg_database(self):
        pass

    def test_read_initial_cg_pdb(self):
        pass



####### test at_mod_p
    def test_shrink_coordinates(self):
        p1, p2 = np.array([58.274,66.912,12.038]), np.array([59.360,66.216,14.859])
        p1a, p2a = np.array([58.450475 , 66.7989   , 12.4964125]), np.array([59.183525 , 66.3291   , 14.4005875])
        p1s, p2s = at_mod_p.shrink_coordinates(p1, p2)
        self.assertIsNone(np.testing.assert_array_equal(p1s, p1a))
        self.assertIsNone(np.testing.assert_array_equal(p2s, p2a))




if __name__ == '__main__':
    unittest.main()
