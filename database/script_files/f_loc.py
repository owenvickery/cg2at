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
import gen, g_var

forcefield_available_prov, fragments_available_prov = gen.read_database_directories()

##### select forcefield

try: 
    forcefield_number = forcefield_available_prov.index(g_var.ff.split('.')[0]+'.ff')
except:
    if g_var.ff != None: 
        print('Cannot find forcefield: '+g_var.ff.split('.')[0]+'.ff  please select one from below\n')
    forcefield_number = gen.database_selection(forcefield_available_prov, 'forcefields')
forcefield_location, forcefield=gen.sort_forcefield(forcefield_available_prov, forcefield_number)

### reads in and sorts fragment information

try: 
    fragment_number = []
    for frag in g_var.fg:
        fragment_number.append(fragments_available_prov.index(frag))
except:
    if g_var.fg != None: 
        print('Cannot find fragment library: '+frag+' please select library from below\n')
    fragment_number = gen.database_selection(fragments_available_prov, 'fragments')


p_directories_unsorted, mod_directories_unsorted, np_directories_unsorted = gen.fetch_residues(fragments_available_prov, fragment_number)

np_residues, p_residues, mod_residues, np_directories, p_directories, mod_directories = gen.sort_directories(p_directories_unsorted, 
																						mod_directories_unsorted, np_directories_unsorted)
### reads in water molecules
water_dir, water = gen.check_water_molecules(g_var.w, np_directories)

### return backbone information

backbone=gen.fetch_fragment(p_directories)

### finds initial rotation matrices
x_rot, y_rot, z_rot=[],[],[]
for angle in range(0,360, 5):
    angle=np.radians(angle)
    x_rot.append(gen.eulerAnglesToRotationMatrix([angle,0,0]))
    y_rot.append(gen.eulerAnglesToRotationMatrix([0,angle,0]))
    z_rot.append(gen.eulerAnglesToRotationMatrix([0,0,angle]))