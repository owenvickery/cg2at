#!/usr/bin/env python3

import os, sys
import numpy as np
from scipy.spatial import cKDTree
import time
import multiprocessing as mp
import gen, g_var, f_loc, at_mod

def get_solvent_names():
    solvent_list = ['ION', 'SOL']
    solvent_dict={}
    for solvent in solvent_list:
        frag_location=at_mod.fragment_location(solvent) ### get fragment location from database
        solvent_type_list, ion_m = at_mod.get_atomistic(frag_location)
        for solvent_type in solvent_type_list:
            for res_type in solvent_type_list[solvent_type]:
                count=0
                for atom in solvent_type_list[solvent_type][res_type]:
                    if solvent_type_list[solvent_type][res_type][atom]['frag_mass'] > 1:
                        count+=1
                if res_type not in solvent_dict:
                    solvent_dict[res_type]=[solvent, count]

    
    return solvent_dict

def convert_AT2CG(at_residues, box_vec):
    solvent_dict = get_solvent_names()
    cg_residues = {}
    for resname in at_residues:
        print('converting residue type: '+resname)
        if resname in solvent_dict or resname == 'SOL':
            cg_residues = convert_solvent(at_residues,resname, solvent_dict, cg_residues)
        else:
            if resname in f_loc.np_residues:
                cg_residues = convert_non_solvent(at_residues, resname, cg_residues)
            elif resname == 'PROTEIN':
                cg_residues = convert_protein(at_residues, resname, cg_residues)
            else:
                print('cannot find: '+resname)
    write_to_pdb(cg_residues, box_vec, solvent_dict)  
    return cg_residues

def print_system_info(cg_residues):
    print('\n{:-<100}'.format(''))
    print('{0:^100}'.format('Script has completed, time for a beer'))
    print('\n{0:^20}{1:^10}'.format('molecules','number'))
    print('{0:^20}{1:^10}'.format('---------','------'))
    for residue_type in cg_residues:
        if residue_type in ['PROTEIN']:
            print('{0:^20}{1:^10}'.format('num of aas', len(cg_residues[residue_type])))
        if residue_type not in ['ION', 'PROTEIN']:
            print('{0:^20}{1:^10}'.format(residue_type, len(cg_residues[residue_type])))
        elif residue_type == 'ION':
            for ion_type in cg_residues[residue_type]:
                print('{0:^20}{1:^10}'.format(ion_type, len(cg_residues[residue_type][ion_type])))

def write_to_pdb(cg_residues, box_vec, solvent_dict):
    pdb_output=gen.create_pdb(g_var.final_dir+'AT2CG_merged.pdb', box_vec)  
    at_count=1
    for residue_type in cg_residues:
        if residue_type not in ['PROTEIN', 'SOL', 'ION']:
            for resid in cg_residues[residue_type]:
                for bead in cg_residues[residue_type][resid]:
                    coord=cg_residues[residue_type][resid][bead]
                    pdb_output.write(g_var.pdbline%((at_count, bead, residue_type,' ',resid,\
                                    coord[0],coord[1],coord[2],1,0))+'\n')
                    at_count+=1
        elif residue_type == 'SOL':
            for resid, coord in enumerate(cg_residues[residue_type]):
                pdb_output.write(g_var.pdbline%((at_count, 'W', residue_type,' ',resid+1,\
                                coord[0],coord[1],coord[2],1,0))+'\n')
                at_count+=1
        elif residue_type == 'ION':
            for ion_type in cg_residues[residue_type]:
                for resid, coord in enumerate(cg_residues[residue_type][ion_type]):
                    pdb_output.write(g_var.pdbline%((at_count, ion_type, residue_type,' ',resid+1,\
                                    coord[0],coord[1],coord[2],1,0))+'\n')
                    at_count+=1
        else:
            for resid in cg_residues[residue_type]:
                for amino_acid in cg_residues[residue_type][resid]:
                    for bead in cg_residues[residue_type][resid][amino_acid]:
                        coord=cg_residues[residue_type][resid][amino_acid][bead]
                        pdb_output.write(g_var.pdbline%((at_count, bead, amino_acid,' ',resid,\
                                        coord[0],coord[1],coord[2],1,0))+'\n')
                        at_count+=1

def convert_protein(at_residues, resname, cg_residues):
    protein_frag = {}
    for residue in f_loc.p_residues:
        frag_location=at_mod.fragment_location(residue) ### get fragment location from database
        residue_type, residue_type_mass = at_mod.get_atomistic(frag_location)
        fragment_num, fragment_name=fetch_fragments(residue_type)
        protein_frag[residue]=fragment_name

    cg_residues[resname]={}
    for residue in at_residues['PROTEIN']:
        cg_residues_temp={}
        cg_residues[resname][residue]={}
        for atom in at_residues[resname][residue]:
            res_type = at_residues[resname][residue][atom]['residue_name']
            if res_type not in cg_residues[resname][residue]:
                cg_residues[resname][residue][res_type]={}
            for frag in protein_frag[res_type]:
                if atom in protein_frag[res_type][frag]:
                    if frag not in cg_residues_temp:
                        cg_residues_temp[frag]={}
                    cg_residues_temp[frag][atom]=dict([('coord',at_residues[resname][residue][atom]['coord'])])
        for fragment in cg_residues_temp:
            frag_com = com_residue(cg_residues_temp[fragment])
            cg_residues[resname][residue][res_type][fragment]=frag_com
    return cg_residues      

def fetch_fragments(residue_type):
    fragment_name={}
    fragment_number={}
    for group in residue_type:
        for frag in residue_type[group]:
            fragment_name[frag]=[]
            fragment_number[frag]=[]
            for atom in residue_type[group][frag]:
                if residue_type[group][frag][atom]['frag_mass'] > 1:
                    fragment_name[frag].append(residue_type[group][frag][atom]['atom'])
                    fragment_number[frag].append(atom)
    return fragment_number, fragment_name

def convert_non_solvent(at_residues, resname, cg_residues):
    frag_location=at_mod.fragment_location(resname) ### get fragment location from database
    residue_type, residue_type_mass = at_mod.get_atomistic(frag_location)
    cg_residues[resname]={}
    fragment_num, fragment_name = fetch_fragments(residue_type)
    for residue in at_residues[resname]:
        cg_residues[resname][residue]={}
        cg_residues_temp={}
        for at_val, atom in enumerate(at_residues[resname][residue]):
            for frag in fragment_num:
                if at_val+1 in fragment_num[frag]:
                    if frag not in cg_residues_temp:
                        cg_residues_temp[frag]={}
                    cg_residues_temp[frag][atom]=dict([('coord',at_residues[resname][residue][atom]['coord'])])
        for fragment in cg_residues_temp:
            frag_com = com_residue(cg_residues_temp[fragment])
            cg_residues[resname][residue][fragment]=frag_com
    return cg_residues      

def convert_solvent(solvent,resname, solvent_dict, cg_residues):
    cg1_solvent=[]
    for residue in solvent[resname]:
        cg1_solvent.append(com_residue(solvent[resname][residue]))
    if resname in solvent_dict:
        atom_name=resname
        residue_size=solvent_dict[resname][1]
        resname = solvent_dict[resname][0]
    if resname == 'SOL':
        if 'residue_size' not in locals():
            residue_size = 4
        cg_residues[resname]=condense_solvent(cg1_solvent, residue_size)
        return cg_residues
    else:
        if resname not in cg_residues:
            cg_residues[resname]={}
        if 'SOL' in cg_residues:
            cg_residues, cg1_solvent = switch_W_to_ION(cg1_solvent, cg_residues)
            cg_residues[resname][atom_name]=cg1_solvent
            return cg_residues
        else:
            cg_residues[resname][atom_name]=cg1_solvent
            return cg_residues
    
def switch_W_to_ION(cg1_solvent, cg_residues):
    tree = cKDTree(cg_residues['SOL'])
    for ion in cg1_solvent:
        radius=0.1
        ndx = tree.query_ball_point(ion, r=radius)
        while len(ndx) <= 1:
            radius+=0.01
            ndx = tree.query_ball_point(ion, r=radius)
        del cg_residues['SOL'][ndx[0]]
        tree = cKDTree(cg_residues['SOL'])
    return cg_residues, cg1_solvent

def com_residue(solvent):
    coord=[]
    for atom in solvent:
        if 'H' not in atom:
            for atom_ch in atom:
                if atom_ch in g_var.mass:
                    coord.append([solvent[atom]['coord'][0],solvent[atom]['coord'][1],solvent[atom]['coord'][2],g_var.mass[atom_ch]])
    return np.average(np.array(coord)[:,:3], axis=0, weights=np.array(coord)[:,3])

def condense_solvent(coordinates, res_size):
#### creates tree of atom coordinates
    tree = cKDTree(coordinates)
    coord_cg = []
    total=len(coordinates)
    while len(coordinates) > 0:
        print('Condensing '+str(total)+' water residues: percentage complete: ',np.round(((total-len(coordinates))/total)*100, 2),'%', end='\r')
        coord_temp = []
        radius=2
        ndx = tree.query_ball_point(coordinates[0], r=radius)
        if len(coordinates) >= res_size: 
            while len(ndx) < res_size:
                radius+=0.01
                ndx = tree.query_ball_point(coordinates[0], r=radius)
        ndx_cut=ndx[:res_size]
        coord_temp = [coordinates[ndx_val] for ndx_val in ndx_cut]
        ndx_cut.sort(reverse = True)
        for ndx_val in ndx_cut:
            del coordinates[ndx_val]
        if len(coordinates) > 0:
            tree = cKDTree(coordinates)
        coord_cg.append(np.mean(coord_temp, axis=0))
    print('                                                                       ', end='\r')
    print('Finished condensing water molecules')       
    return coord_cg

def write_topology(cg_residues):
    os.chdir(g_var.final_dir)
    if not os.path.exists('topol_final.top'):
        with open('topol_final.top', 'w') as topol_write:       
            topol_write.write('; Include forcefield parameters\n#include \"'+g_var.final_dir+f_loc.forcefield+'/'+f_loc.forcefield[:-4]+'.itp\"\n')
            topol_write.write('[ system ]\n; Name\nSomething clever....\n\n[ molecules ]\n; Compound        #mols\n')
            for residue_type in cg_residues:
                if residue_type in ['PROTEIN']:
                    topol_write.write('{0:20}{1:^10}\n'.format('Protein', 'XXX'))
                if residue_type not in ['ION', 'PROTEIN']:
                    topol_write.write('{0:20}{1:^10}\n'.format(residue_type, len(cg_residues[residue_type])))
                elif residue_type == 'ION':
                    for ion_type in cg_residues[residue_type]:
                        topol_write.write('{0:20}{1:^10}\n'.format(ion_type, len(cg_residues[residue_type][ion_type])))