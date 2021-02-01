#!/usr/bin/env python3

import os, sys
import numpy as np
import math
import gen, g_var, at_mod


def build_atomistic_system(residue_type, a):
    system={}
    atomistic_fragments={}
#### for each residue type covert to atomistic except protein
    # for residue_type in [key for value, key in enumerate(cg_residues) if key not in ['PROTEIN']]:
    if residue_type not in system and residue_type != 'ION':
        system[residue_type] = 0
#### reset counters for each residue type
    print('Converting residue type: ' +residue_type)
#### creates folder for residue type
    gen.mkdir_directory(g_var.working_dir+residue_type)
    if residue_type in ['ION','SOL']:
        if os.path.exists(g_var.working_dir+'SOL'+'/SOL_all.pdb') and residue_type == 'SOL':
            atomistic_fragments[residue_type], system['SOL'] = read_solvent_conversion(residue_type, g_var.cg_residues[residue_type])
        else:
            atomistic_fragments[residue_type], solvent_number = atomistic_non_protein_solvent(residue_type, g_var.cg_residues[residue_type])
        if residue_type in ['ION']:
            system = solvent_ion( system, atomistic_fragments, residue_type)
        elif not os.path.exists(g_var.working_dir+'SOL'+'/SOL_all.pdb'):
            system['SOL'] = solvent_number
            write_solvent( system, atomistic_fragments, residue_type)   
    else:
        if not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_all.pdb'):
            atomistic_fragments[residue_type], system[residue_type] = atomistic_non_protein_non_solvent(residue_type, g_var.cg_residues[residue_type])
            write_solvent(system, atomistic_fragments, residue_type)
        else:
            system[residue_type]=len(cg_residues) 
            
    print('{:<100}'.format(''), end='\r')
    print('Finished converting: '+residue_type)
    return system 

def write_solvent(system, atomistic_fragments, residue_type):
    merge_coord = [line['coord'] for line in atomistic_fragments[residue_type]]
    coord, index_conversion = at_mod.index_conversion_generate(atomistic_fragments[residue_type], merge_coord)
    at_mod.write_pdb(atomistic_fragments[residue_type], coord, index_conversion, g_var.working_dir+residue_type+'/'+residue_type+'_all.pdb')


def solvent_ion(system, atomistic_fragments, residue_type):
    #### creates ion pdb with header
    gen.mkdir_directory(g_var.working_dir+'ION/MIN')
    skip = os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb')
    if not skip:
        pdb_ion = gen.create_pdb(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb')
    # for resid in atomistic_fragments[residue_type]:
    for at_id, atom in enumerate(atomistic_fragments[residue_type]):
    #### write ion coordinate out
        if not skip:
            x, y, z = gen.trunc_coord([atom['coord'][0],atom['coord'][1],atom['coord'][2]])
            pdb_ion.write(g_var.pdbline%((at_id+1,atom['atom_name'],atom['residue_name'],' ',1,x,y,z,0.00, 0.00))+'\n')
        if atom['residue_name'] not in system:
            system[atom['residue_name']]=1
        else:
            system[atom['residue_name']]+=1
    return system

def atomistic_non_protein_non_solvent(cg_residue_type,cg_residues):
    atomistic_fragments={}  #### residue dictionary
#### run through every residue in a particular residue type
    residue_type={}
    residue_type_mass={}
    atomistic_list = []
    for cg_resid, cg_residue in enumerate(cg_residues):
        atomistic_fragments[cg_resid]={}
        frag_location=gen.fragment_location(cg_residue_type) ### get fragment location from database
        residue_type[cg_residue_type], residue_type_mass[cg_residue_type] = at_mod.get_atomistic(frag_location)
        for group in residue_type[cg_residue_type]:
            center, at_frag_centers, cg_frag_centers, group_fit = at_mod.rigid_fit(residue_type[cg_residue_type][group], residue_type_mass[cg_residue_type], cg_residue, cg_residues[cg_residue])
            at_connect, cg_connect = at_mod.connectivity(cg_residues[cg_residue], at_frag_centers, cg_frag_centers, group_fit, group)
            xyz_rot_apply=at_mod.get_rotation(cg_connect, at_connect, center, cg_residue_type, group, cg_resid)
            atomistic_fragments = at_mod.apply_rotations(atomistic_fragments,cg_resid, group_fit, center, xyz_rot_apply)
        atomistic_fragments[cg_resid] = at_mod.check_hydrogens(atomistic_fragments[cg_resid])
        for atom in atomistic_fragments[cg_resid].values():
            atom['atom_name'] = atom.pop('atom')
            atom['residue_name'] = atom.pop('res_type')
            atom['residue_id'] = atom.pop('resid')
            atomistic_list.append(atom)
    return atomistic_list, len(atomistic_fragments)
   

def read_solvent_conversion(cg_residue_type,cg_residues):
    residue_type={}
    residue_type_mass={}
    for cg_resid, cg_residue in enumerate(cg_residues):
        frag_location=gen.fragment_location(cg_residue_type) ### get fragment location from database
        residue_type[cg_residue_type], residue_type_mass[cg_residue_type] = at_mod.get_atomistic(frag_location)
        for res_type in residue_type[cg_residue_type]:
            if g_var.water in residue_type[cg_residue_type][res_type]:
                sol_p_bead = 0
                for atom in residue_type[cg_residue_type][res_type][g_var.water].values():
                    if int(atom['resid_ori']) > sol_p_bead:
                        sol_p_bead = int(atom['resid_ori'])
                return sol_p_bead, sol_p_bead*len(cg_residues)

def atomistic_non_protein_solvent(cg_residue_type,cg_residues): 
    atomistic_fragments={}  #### residue dictionary
#### run through every residue in a particular residue type
    residue_type={}
    residue_type_mass={}
    atomistic_fragments_list = []
    for cg_resid, cg_residue in enumerate(cg_residues):
        frag_location=gen.fragment_location(cg_residue_type) ### get fragment location from database
        residue_type[cg_residue_type], residue_type_mass[cg_residue_type] = at_mod.get_atomistic(frag_location)
        sol_p_bead = 0    
        atomistic_fragments[cg_resid]={}
        for res_type in residue_type[cg_residue_type]:
            if g_var.water in residue_type[cg_residue_type][res_type]:
                center, at_frag_centers, cg_frag_centers, group_fit = at_mod.rigid_fit(residue_type[cg_residue_type][res_type], residue_type_mass[cg_residue_type]
                                                                                       , cg_residue, cg_residues[cg_residue])
                xyz_rot_apply=gen.AnglesToRotMat([np.random.uniform(0, math.pi*2), np.random.uniform(0, math.pi*2), np.random.uniform(0, math.pi*2)])
                # atomistic_fragments = at_mod.apply_rotations(atomistic_fragments,cg_resid, group_fit, center, xyz_rot_apply)
                for bead in at_mod.apply_rotations(atomistic_fragments,cg_resid, group_fit, center, xyz_rot_apply)[cg_resid].values():
                    bead['atom_name'] = bead.pop('atom')
                    bead['residue_name'] = bead.pop('res_type')
                    bead['residue_id'] = bead.pop('resid')
                    atomistic_fragments_list.append(bead)    
                    if bead['resid_ori'] > sol_p_bead:
                        sol_p_bead = bead['resid_ori']

    return atomistic_fragments_list, sol_p_bead*len(cg_residues)