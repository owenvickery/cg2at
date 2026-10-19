#!/usr/bin/env python3

import os
import copy
import numpy as np
import gen, g_var, at_mod
from exceptions import CG2ATError



def build_atomistic_system(residue_type):
    system={}
#### for each residue type covert to atomistic except protein
    print('Converting residue type: ' +residue_type)
#### creates folder for residue type
    gen.mkdir_directory(g_var.working_dir+residue_type)
    if not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb'):
        atomistic_fragments, system[residue_type] = at_np_solvent(residue_type, g_var.cg_residues[residue_type])
        write_solvent(atomistic_fragments, residue_type)
    elif residue_type in g_var.sol_residues:
        atomistic_fragments, system[residue_type] = read_solvent_conversion(residue_type, g_var.cg_residues[residue_type])
    else:
        system[residue_type]=len(g_var.cg_residues[residue_type])
            
    print('{:<100}'.format(''), end='\r')
    print('Finished converting: '+residue_type)
    return system 

def write_solvent(atomistic_fragments, residue_type):
    merge_coord = [line['coord'] for line in atomistic_fragments]
    coord, index_conversion = at_mod.index_conversion_generate(atomistic_fragments, merge_coord)
    at_mod.write_pdb(atomistic_fragments, coord, index_conversion, g_var.working_dir+residue_type+'/'+residue_type+'_all.pdb')

def read_solvent_conversion(cg_residue_type,cg_residues):
    residue_type={}
    residue_type_mass={}
    for cg_resid, cg_residue in enumerate(cg_residues):
        frag_location=gen.fragment_location(cg_residue_type) ### get fragment location from database
        residue_type[cg_residue_type], residue_type_mass[cg_residue_type] = at_mod.get_atomistic(frag_location, cg_residue_type)

        for group in residue_type[cg_residue_type]:
            sol_p_bead = 0
            for frag in residue_type[cg_residue_type][group].values():
                for atom in frag.values():
                    if int(atom['resid_ori']) > sol_p_bead:
                        sol_p_bead = int(atom['resid_ori'])
            return sol_p_bead, sol_p_bead*len(cg_residues)
    
    raise CG2ATError('There is an issue with the solvent recalculation')

def at_np_solvent(cg_residue_type,cg_residues):   
    atomistic_fragments={}  #### residue dictionary
#### run through every residue in a particular residue type
    residue_type={}
    residue_type_mass={}
    atomistic_fragments_list = []
    # Fragment location and template are constant for all copies of the same
    # residue type — look them up once before the loop.
    frag_location = gen.fragment_location(cg_residue_type)
    frag_template, frag_mass_template = at_mod.get_atomistic(frag_location, cg_residue_type)

    total = len(cg_residues)
    for cg_resid, cg_residue in enumerate(cg_residues):
        gen.print_progress('Converting ' + cg_residue_type, cg_resid + 1, total)
        atomistic_fragments[cg_resid] = {}
        # Deep-copy the template so each residue gets an independent, mutable copy.
        residue_type[cg_residue_type] = copy.deepcopy(frag_template)
        residue_type_mass[cg_residue_type] = copy.deepcopy(frag_mass_template)
        for group in residue_type[cg_residue_type]:
            center, at_frag_centers, cg_frag_centers, group_fit = at_mod.rigid_fit(residue_type[cg_residue_type][group], 
                                                                                    residue_type_mass[cg_residue_type], 
                                                                                    cg_residue, cg_residues[cg_residue])
            at_connect, cg_connect = at_mod.connectivity(cg_residues[cg_residue], at_frag_centers, cg_frag_centers, group_fit, group)
            xyz_rot_apply=at_mod.get_rotation(cg_connect, at_connect, center, cg_residue_type, group, cg_resid)

            atomistic_fragments = at_mod.apply_rotations(atomistic_fragments,cg_resid, group_fit, center, xyz_rot_apply)
            
        if cg_residue_type in g_var.np_residues:
            atomistic_fragments[cg_resid] = at_mod.check_hydrogens(atomistic_fragments[cg_resid])
        if len(g_var.res_top[cg_residue_type]['CHIRAL']) > 1:
            chiral = at_mod.get_chiral_non_carbonyl(atomistic_fragments[cg_resid])
            atomistic_fragments[cg_resid] = at_mod.correct_chiral_atoms(atomistic_fragments[cg_resid], chiral)            
        atomistic_fragments_list, sol_p_bead =  sort_np_dictionary(atomistic_fragments[cg_resid], atomistic_fragments_list)
    gen.finish_progress('Converting ' + cg_residue_type, total)
    if cg_residue_type in g_var.sol_residues:
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
