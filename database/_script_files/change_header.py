#!/usr/bin/env python3

import os, sys
import numpy as np
import math
from distutils.dir_util import copy_tree
from shutil import copyfile
import glob
import re
import shlex
import argparse


parser = argparse.ArgumentParser(description='Converts CG representation into an atomistic representation', prog='CG2AT', epilog='Enjoy the program and best of luck!\n')
parser.add_argument('-f', help='choose your forcefield. (Optional)',metavar='charmm36',type=str)
args = parser.parse_args()
options = vars(args)

def sep_fragments_header(line, residue_name):
    line = line.replace('[','')
    line = line.replace(']','')
    line_sep = shlex.split(line)
    residue = {}
    residue['ter']=False
    cont=False
    for top in line_sep:
        if cont == True:
            break
        if top == 'frag':
            cont=True
    if not cont:
        sys.exit('error in header') 
    return top


def read_pdb_in(pdb_file):
    pdb_lines=[]
    with open(pdb_file, 'r') as pdb_input:
        with open(pdb_file+'2', 'w') as pdb_output:
            for line_nr, line in enumerate(pdb_input.readlines()):
                if line.startswith('['):
                    # print(line.strip())
                    header = sep_fragments_header(line.strip(), 'a')
                    pdb_output.write('[ '+header+' ]\n')
                else:
                    pdb_output.write(line)
    os.system('mv '+pdb_file+'2 '+pdb_file)

read_pdb_in(args.f+'/'+args.f+'.pdb')

