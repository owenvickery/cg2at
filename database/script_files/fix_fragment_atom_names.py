#!/usr/bin/env python3

import os, sys
import numpy as np
import argparse
import re
import shlex
from shutil import copyfile


parser = argparse.ArgumentParser(description='fixes atom names in fragment pdb file', epilog='Enjoy the program and best of luck!\n', allow_abbrev=True)
parser.add_argument('-f', help='fragment pdb file',metavar='pdb',type=str, required=True)
parser.add_argument('-i', help='itp of molecule',metavar='itp',type=str, required=True)
args = parser.parse_args()
options = vars(args)


def fetch_itp_atom_names(itp_file):
    atoms=False
    atom_names=[]
    with open(itp_file, 'r') as itp_input:
        for line_nr, line in enumerate(itp_input.readlines()):
            if len(line.split()) >= 3:
                if not line.startswith(';') :
                    if line.split()[0] == '[':
                        if line.split()[1] == 'atoms':
                            atoms=True
                        else:
                            atoms=False               
                    elif atoms:
                        atom_names.append(line.split()[4])
    return atom_names

def read_pdb_in(pdb_file):
    pdb_lines=[]
    with open(pdb_file, 'r') as pdb_input:
        for line_nr, line in enumerate(pdb_input.readlines()):
            pdb_lines.append(line)
    return pdb_lines

def pdbatom(line):
### get information from pdb file
### atom number, atom name, residue name,chain, resid,  x, y, z, backbone (for fragment), connect(for fragment)
    try:
        return dict([('atom_number',int(line[7:11].replace(" ", ""))),('atom_name',str(line[12:16]).replace(" ", "")),('residue_name',str(line[17:21]).replace(" ", "")),\
            ('chain',line[21]),('residue_id',int(line[22:26])), ('x',float(line[30:38])),('y',float(line[38:46])),('z',float(line[46:54])), ('backbone',int(float(line[55:61]))),\
            ('connect',int(float(line[62:67])))])
    except:
        sys.exit('\npdb line is wrong:\t'+line)

def swap_atom_name(atom):
    if atom not in atom_names:
        if len(atom)==4:
            atom_switch=atom[1:]+atom[0]
            if atom_switch in atom_names:
                return atom_switch
    return False

def final_split(seg, switch):
    seg_split = seg.split()
    for seg_split_group in seg_split:
        atom_switch = swap_atom_name(seg_split_group)
        if atom_switch:
            switch[seg_split_group] = atom_switch   
    return switch 

def parse_pdb(pdb_file, pdb_lines):
    with open(pdb_file, 'w') as pdb_output:
        for line in pdb_lines:
            if line.startswith('ATOM'):
                line_sep = pdbatom(line)
                atom_switch = swap_atom_name(line_sep['atom_name'])
                if atom_switch:
                    line= re.sub(line_sep['atom_name'], atom_switch,line)
                pdb_output.write(line)
            else:
                switch = {}
                for seg in shlex.split(line):
                    if seg.count(':') >= 1:
                        seg_split = seg.split(':')
                        for seg_split_group in seg_split:
                            if seg_split_group.count(',') >= 1:
                                seg_split_comma = seg_split_group.split(',')
                                for sg_st_cm in seg_split_comma:
                                    switch = final_split(sg_st_cm, switch)
                            else:
                                switch = final_split(seg_split_group, switch)
                for atom in switch:
                    line = re.sub(atom, switch[atom],line)
                pdb_output.write(line)

atom_names = fetch_itp_atom_names(args.i)
pdb_lines = read_pdb_in(args.f)
copyfile(args.f, '_'+args.f)
parse_pdb(args.f, pdb_lines)
# os.system('gawk -i inplace -F\'hydrogen\' \'{print $1}\' '+args.f)
# os.system('gawk -i inplace -F\'con\' \'{print $1"  "}\' '+args.f)



