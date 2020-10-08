#!/usr/bin/env python3
import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__))+'/..')
import unittest
# from gen import fix_time
import gen, g_var

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
        swap = 'POPE,NH3:POPG,GL0'# POPG,GL0:POPE,NH3'
        swap='NA+:skip:4000-4100'
        result, result2 = gen.split_swap(swap)
        print(result, result2)
if __name__ in ['__main__', 'test_gen']:
    unittest.main()

