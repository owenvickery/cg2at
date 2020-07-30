#!/usr/bin/env python3

import os, sys
import numpy as np
import copy
from string import ascii_uppercase
import difflib
from scipy.spatial import cKDTree
import multiprocessing as mp
import gen, g_var, f_loc, at_mod, read_in
import math


def build_protein_atomistic_system():   
#### initisation of counters
    chain_count=0
    g_var.backbone_coords[chain_count]=[]
    g_var.o_system[chain_count]=[]
    g_var.coord_atomistic[chain_count]={}
    g_var.seq_cg[chain_count]=[]
    print('Converting other')
    gen.mkdir_directory(g_var.working_dir+'OTHER')  ### make and change to protein directory
#### for each residue in protein
    initial=True
    residue_type={}
    residue_type_mass={}
    for cg_residue_id, residue_number in enumerate(cg_residues):
        resname = cg_residues[residue_number][next(iter(cg_residues[residue_number]))]['residue_name']
        if cg_residue_id == 0:
            g_var.o_system[chain_count].append(terminal_residue(resname)[0]) 
        g_var.coord_atomistic[chain_count][residue_number]={}
        frag_location=gen.fragment_location(resname, f_loc.database_locations) ### get fragment location from database
        residue_type[resname], residue_type_mass[resname] = at_mod.get_atomistic(frag_location)
        for group in residue_type[resname]:
            center, at_frag_centers, cg_frag_centers, group_fit = at_mod.rigid_fit(residue_type[resname][group], residue_type_mass[resname], residue_number, cg_residues[residue_number])
            at_connect, cg_connect = at_mod.connectivity(cg_residues[residue_number], at_frag_centers, cg_frag_centers, group_fit, group)
            if BB_bead in group_fit:
                BB_connect = []
                at_connect, cg_connect, new_chain = BB_connectivity(at_connect,cg_connect, cg_residues, group_fit[BB_bead], residue_number, BB_bead)
                g_var.seq_cg = add_to_sequence(g_var.seq_cg, resname, chain_count)
                g_var.backbone_coords[chain_count].append(np.append(cg_residues[residue_number][BB_bead]['coord'], 1))     
            if len(at_connect) == len(cg_connect) and len(at_connect) != 0:
                xyz_rot_apply=at_mod.kabsch_rotate(np.array(at_connect)-center, np.array(cg_connect)-center)
            elif len(at_connect) == 0:
                xyz_rot_apply = False
                print('Cannot find any connectivity for residue number: '+str(residue_number)+', residue type: '+str(resname)+', group: '+str(group))
            else:
                print('atom connections: '+str(len(at_connect))+' does not equal CG connections: '+str(len(cg_connect)))
                sys.exit('residue number: '+str(residue_number)+', residue type: '+str(resname)+', group: '+group)

            for bead in group_fit:
                for atom in group_fit[bead]:
                    group_fit[bead][atom]['coord'] = at_mod.rotate_atom(group_fit[bead][atom]['coord'], center, xyz_rot_apply)   
                    atom_new = group_fit[bead][atom].copy()
                    g_var.coord_atomistic[chain_count][residue_number][atom] = atom_new
        if new_chain:
            g_var.o_system[chain_count].append(terminal_residue(resname)[1])
            chain_count+=1
            if cg_residue_id+1 != len(cg_residues):
                g_var.backbone_coords[chain_count]=[]
                g_var.coord_atomistic[chain_count]={}
                g_var.o_system[chain_count]=[]
                g_var.o_system[chain_count].append(terminal_residue(cg_residues[residue_number+1][next(iter(cg_residues[residue_number+1]))]['residue_name'])[0])
                g_var.seq_cg[chain_count]=[]




def atomistic_other_residue_types():
    print('Converting Protein')
    gen.mkdir_directory(g_var.working_dir+'OTHER')
    atomistic_fragments={}  #### residue dictionary
#### run through every residue in a particular residue type
    residue_type={}
    residue_type_mass={}
    for cg_resid, cg_res in enumerate(cg_residues['OTHER']):
        atomistic_fragments[cg_resid]={}
        frag_location=gen.fragment_location(cg_residue_type, f_loc.database_locations) ### get fragment location from database
        residue_type[cg_residue_type], residue_type_mass[cg_residue_type] = at_mod.get_atomistic(frag_location)
        for group in residue_type[cg_residue_type]:
            center, at_frag_centers, cg_frag_centers, group_fit = at_mod.rigid_fit(residue_type[cg_residue_type][group], residue_type_mass[cg_residue_type], cg_res, cg_residues[cg_res])
            at_connect, cg_connect = at_mod.connectivity(cg_residues[cg_res], at_frag_centers, cg_frag_centers, group_fit, group)
            if len(at_connect) == len(cg_connect) and len(cg_connect) > 0:
                try:
                    xyz_rot_apply=at_mod.kabsch_rotate(np.array(at_connect)-center, np.array(cg_connect)-center)
                except:
                    sys.exit('There is a issue with residue: '+cg_residue_type+' in group: '+str(group))
            else:
                print('atom connections: '+str(len(at_connect))+' does not match CG connections: '+str(len(cg_connect)))
                sys.exit('residue number: '+str(cg_resid)+', residue type: '+str(cg_residue_type)+', group: '+group)
            for bead in group_fit:
                for atom in group_fit[bead]:
                    group_fit[bead][atom]['coord'] = at_mod.rotate_atom(group_fit[bead][atom]['coord'], center, xyz_rot_apply)   
                    atom_new = group_fit[bead][atom].copy()
                    atomistic_fragments[cg_resid][atom] = atom_new
    return atomistic_fragments   

def non_solvent(system, atomistic_fragments, residue_type):
#### loop through all resids of that residue type 
    if not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_all.pdb') and not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb'):
        pdb_output_all = gen.create_pdb(g_var.working_dir+residue_type+'/'+residue_type+'_all.pdb')     
        pdb_all_coord = []   
        pdb_output_all_temp = []
    for resid in atomistic_fragments[residue_type]:
        # skip = os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_'+str(resid)+'.pdb')
        if not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb'):
            if not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_'+str(resid)+'.pdb'):
                pdb_output_ind = gen.create_pdb(g_var.working_dir+residue_type+'/'+residue_type+'_'+str(resid)+'.pdb')         
                atomistic_fragments[residue_type][resid] = at_mod.check_hydrogens(atomistic_fragments[residue_type][resid])
            ####### check if any atoms in residue overlap #######
                coord=[]
                for atom in atomistic_fragments[residue_type][resid]:
                    coord.append(atomistic_fragments[residue_type][resid][atom]['coord'])
                coord=at_mod.check_atom_overlap(coord)
                for atom_val, atom in enumerate(atomistic_fragments[residue_type][resid]):
                    atomistic_fragments[residue_type][resid][atom]['coord']=coord[atom_val]

                for at_id, atom in enumerate(atomistic_fragments[residue_type][resid]):
                #### write residue out to a pdb file
                    short_line=atomistic_fragments[residue_type][resid][at_id+1]
                    x, y, z = gen.trunc_coord([short_line['coord'][0], short_line['coord'][1],short_line['coord'][2]])
                    pdb_output_ind.write(g_var.pdbline%((at_id+1,short_line['atom'],short_line['res_type'],' ',1,x,y,z,0.00, 0.00))+'\n')
                    if 'pdb_output_all' in locals():
                        x, y, z = gen.trunc_coord([short_line['coord'][0], short_line['coord'][1],short_line['coord'][2]])
                        pdb_all_coord.append([x, y, z])
                        pdb_output_all_temp.append([at_id+1, short_line['atom'],short_line['res_type'],' ',1,x, y, z, 0.00, 0.00])
    if 'pdb_output_all' in locals():
        p_a_c=at_mod.check_atom_overlap(pdb_all_coord)
        for at_val, line in enumerate(pdb_output_all_temp):
            pdb_output_all.write(g_var.pdbline%((line[0],line[0],line[0],' ',1,p_a_c[at_val][0], p_a_c[at_val][1], p_a_c[at_val][2],0.00, 0.00))+'\n')

    system[residue_type]=int(resid)+1
    return system