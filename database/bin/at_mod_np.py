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
            atomistic_fragments[residue_type], solvent_number = at_np_solvent(residue_type, g_var.cg_residues[residue_type])
        if residue_type in ['ION']:
            system = solvent_ion( system, atomistic_fragments, residue_type)
        elif not os.path.exists(g_var.working_dir+'SOL'+'/SOL_all.pdb'):
            system['SOL'] = solvent_number
            write_solvent(system, atomistic_fragments, residue_type)   
    else:
        if not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb'):
            atomistic_fragments[residue_type], system[residue_type] = at_np_solvent(residue_type, g_var.cg_residues[residue_type])
            write_solvent(system, atomistic_fragments, residue_type)
        else:
            system[residue_type]=len(g_var.cg_residues[residue_type]) 
            
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

def at_np_solvent(cg_residue_type,cg_residues):   
    atomistic_fragments={}  #### residue dictionary
#### run through every residue in a particular residue type
    residue_type={}
    residue_type_mass={}
    atomistic_fragments_list = []
    for cg_resid, cg_residue in enumerate(cg_residues):
        atomistic_fragments[cg_resid]={}
        frag_location=gen.fragment_location(cg_residue_type) ### get fragment location from database
        residue_type[cg_residue_type], residue_type_mass[cg_residue_type] = at_mod.get_atomistic(frag_location)

        for group in residue_type[cg_residue_type]:
            if convert_question(residue_type[cg_residue_type][group], cg_residue_type, cg_residues[cg_residue]):
                center, at_frag_centers, cg_frag_centers, group_fit = at_mod.rigid_fit(residue_type[cg_residue_type][group], 
                                                                                        residue_type_mass[cg_residue_type], 
                                                                                        cg_residue, cg_residues[cg_residue])
                at_connect, cg_connect = at_mod.connectivity(cg_residues[cg_residue], at_frag_centers, cg_frag_centers, group_fit, group)
                xyz_rot_apply=at_mod.get_rotation(cg_connect, at_connect, center, cg_residue_type, group, cg_resid)

                atomistic_fragments = at_mod.apply_rotations(atomistic_fragments,cg_resid, group_fit, center, xyz_rot_apply)
        if cg_residue_type not in ['SOL','ION']:
            atomistic_fragments[cg_resid] = at_mod.check_hydrogens(atomistic_fragments[cg_resid])
        atomistic_fragments_list, sol_p_bead =  sort_np_dictionary(atomistic_fragments[cg_resid], atomistic_fragments_list)
    if cg_residue_type in ['SOL']:
        return atomistic_fragments_list, sol_p_bead*len(cg_residues)
    else:
        return atomistic_fragments_list, len(atomistic_fragments)

def sort_np_dictionary(atomistic_fragments, atomistic_fragments_list):
    sol_p_bead = 0
    dict_start = np.min(list(atomistic_fragments.keys()))
    for atom_val in range(dict_start, len(atomistic_fragments)+dict_start):
        atom = atomistic_fragments[atom_val]
        atom['atom_name'] = atom.pop('atom')
        atom['residue_name'] = atom.pop('res_type')
        atom['residue_id'] = atom.pop('resid')
        atomistic_fragments_list.append(atom)    
        if atom['resid_ori'] > sol_p_bead:
            sol_p_bead = atom['resid_ori']
    return atomistic_fragments_list, sol_p_bead 


def convert_question(res_type, cg_residue_type, bead):
    if cg_residue_type == 'SOL' and g_var.water in res_type:
        return True
    elif cg_residue_type == 'ION' and list(bead.keys())[0] in res_type:
        return True 
    elif cg_residue_type not in ['SOL', 'ION']:
        return True
    else:
        return False
