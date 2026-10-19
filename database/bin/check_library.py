#!/usr/bin/env python3

import os, sys
import gen, g_var


def write_posre_file(residue, posres):
    with open(residue+'_posre.itp', 'w') as posre_output:
        posre_output.write('; position restraints for '+g_var.args.posre+'\n[ position_restraints ]\n;  i    funct       fcx        fcy        fcz\n')
        for atom in posres[g_var.args.posre]:
            posre_output.write(atom)
        print('Created posres file for '+g_var.args.posre)

def append_ifdef(frag_location):
    with open(frag_location, 'a') as itp_input:
        itp_input.write('\n; Include heavy atom position restraint file\n#ifdef NP\n#include "'+g_var.args.posre+'_posre.itp"\n#endif\n')
    print('Added ifdef flag to '+g_var.args.posre)

def read_itp(filename_itp, footer=False):
    posre = {}
    molecule ={}
    mol_type = False
    atoms=False
    if os.path.exists(filename_itp):
        with open(filename_itp, 'r') as itp_input:
            for line in itp_input.readlines():
                if not (line.isspace() or len(line) == 0) and not line.startswith(';'):
                    line_sep = line.split()
                    if line.startswith('[') and gen.strip_header(line) == 'moleculetype':
                        mol_type = True
                    elif line.startswith('[') and gen.strip_header(line) == 'atoms':
                        atoms, mol_type = True, False
                    elif line.startswith('['):    
                        atoms = False
                        mol_type = False  
                    elif mol_type:
                        molecule[line_sep[0]] = []
                        posre[line_sep[0]] = []
                        mol = line_sep[0]
                    elif atoms:
                        molecule[mol].append(line_sep[4])
                        if len(line.split()) >= 8:
                            if int(float(line.split()[7])) > 1:
                                posre[mol].append('   {0:5}{1:3}{2:11}{3:11}{4:11}\n'.format(line.split()[0],1,50,50,50))
                    if '#ifdef NP' in line:
                        footer=True
        return molecule, posre, footer
    else:
        sys.exit('Cannot find itp file: '+filename_itp)
    

def check_frag_file(directory, molecule):
    exists = []
    wrong = []
    for mol, beads in molecule.items():
        beads.sort()
        bead_list = []
        if os.path.exists(directory+mol+'/'+mol+'.pdb'):
            exists.append(mol)
            with open(directory+mol+'/'+mol+'.pdb', 'r') as pdb_input:
                for line_nr, line in enumerate(pdb_input.readlines()):
                    if line.startswith('['):
                        bead_list.append(gen.strip_header(line))
                bead_list.sort()

            if bead_list == beads:
                print(mol, 'Correct')
            else:
                l1=mol+' Incorrect: '+directory+mol+'/'+mol+'.pdb'
                l2='\nMartini '+mol+': '+', '.join(map(str, beads))
                l3='\nCG2AT2  '+mol+': '+', '.join(map(str, bead_list))
                wrong.append([l1+l2+l3])
    for residue in wrong:
        print(residue)
    return exists

def check_fragments_missing_from_itp(directory, molecule):
    print('\n-------------------')
    print('Fragment files missing from itp\n')
    for filename in os.listdir(directory):
        if filename not in molecule:
            print(filename)

def check_fragments_to_add(directory, molecule, exists):
    print('\n-------------------')
    print('Fragment files to add\n')   
    count = 0 
    for mol, beads in molecule.items():
        if not os.path.exists(directory+mol+'/'+mol+'.pdb'):
            count+=1
            print('Fragment entry required for: ',mol)
    print('Currently missing: ',count)

def add_posres_file():
    print()
    frag_location=gen.fragment_location(g_var.args.posre)[:-4]
    molecule, posre, footer = read_itp(frag_location+'.itp')
    write_posre_file(frag_location, posre)
    if not footer:
        append_ifdef(frag_location+'.itp')
    sys.exit()

def compare_forcefield_to_database():
    molecule, posre, footer = read_itp(g_var.args.compare)
    for directory in g_var.np_directories:
        print('\nChecking: ', directory[0],'\n')
        exists = check_frag_file(directory[0], molecule)
        check_fragments_missing_from_itp(directory[0], molecule)
        check_fragments_to_add(directory[0], molecule, exists)
    sys.exit()

