#!/usr/bin/env python3
import os, sys
import numpy as np
from subprocess import Popen, PIPE
import subprocess, shlex
from time import gmtime, strftime
import math
import multiprocessing as mp
import argparse
import copy
from shutil import copyfile
from distutils.dir_util import copy_tree
import time
from string import ascii_uppercase
from pathlib import Path
import re
import datetime
import glob
from scipy.spatial import KDTree
import difflib
import gen

parser = argparse.ArgumentParser(description='Converts CG representation into an atomistic representation', epilog='Enjoy the program and best of luck!\n', allow_abbrev=True)
parser.add_argument('-c', help='coarse grain coordinates',metavar='pdb/gro/tpr',type=str, required=True)
parser.add_argument('-a', help='atomistic coordinates',metavar='pdb/gro/tpr',type=str)
parser.add_argument('-v', action="count", default=0, help="increase output verbosity (eg -vv, 3 levels)")
parser.add_argument('-ter', help='interactively choose terminal species', action='store_true')
parser.add_argument('-clean', help='removes all part files from build', action='store_true')
# parser.add_argument('-chains', help='list of chains to rigid body fit together, starts at chain 0',metavar='1 2',type=int, nargs='*')
parser.add_argument('-w', help='choose your solvent, common choices are: tip3p, tip4p, spc and spce. This is optional',metavar='tip3p',type=str)
parser.add_argument('-ff', help='choose your forcefield. This is optional',metavar='charmm36',type=str)
parser.add_argument('-fg', help='choose your fragment library. This is optional',metavar='charmm36',type=str, nargs='*')
args = parser.parse_args()
options = vars(args)
initialisation_time=time.time()

cg, at, verbose, ter, clean, water, force, frag= args.c, args.a, args.v , args.ter, args.clean, args.w, args.ff, args.fg

timestamp =  strftime("%Y-%m-%d_%H-%M-%S", gmtime())

start_dir=os.getcwd()+'/'  ### initial working directory
working_dir=os.getcwd()+'/CG2AT_'+timestamp+'/'   ### working directory 
final_dir=os.getcwd()+'/CG2AT_'+timestamp+'/FINAL/'  ### final directory for run files
input_directory=os.getcwd()+'/CG2AT_'+timestamp+'/INPUT/'  ### contains input run files
database_dir=os.path.dirname(os.path.realpath(__file__))+'/../'


forcefield_available_prov, fragments_available_prov = gen.read_database_directories()

##### select forcefield

try: 
    forcefield_number = forcefield_available_prov.index(args.ff.split('.')[0]+'.ff')
except:
    if args.ff != None: 
        print('Cannot find forcefield: '+args.ff.split('.')[0]+'.ff  please select one from below\n')
    forcefield_number = gen.database_selection(forcefield_available_prov, 'forcefields')
forcefield_location, forcefield=gen.sort_forcefield(forcefield_available_prov, forcefield_number)

### reads in and sorts fragment information

try: 
    fragment_number = []
    for frag in args.fg:
        fragment_number.append(fragments_available_prov.index(frag))
except:
    if args.fg != None: 
        print('Cannot find fragment library: '+frag+' please select library from below\n')
    fragment_number = gen.database_selection(fragments_available_prov, 'fragments')


p_directories_unsorted, mod_directories_unsorted, np_directories_unsorted=gen.fetch_residues(fragments_available_prov, fragment_number)

np_residues, p_residues, mod_residues,np_directories, p_directories, mod_directories = gen.sort_directories(p_directories_unsorted, 
																						mod_directories_unsorted, np_directories_unsorted, args.v)
### reads in water molecules
water_dir, water = gen.check_water_molecules(args.w, np_directories)