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
        atomistic_fragments[residue_type], solvent_number = atomistic_non_protein_solvent(residue_type, g_var.cg_residues[residue_type])
        if residue_type in ['ION']:
            system = solvent_ion( system, atomistic_fragments, residue_type)
        else:
            if solvent_number != 0:
                system['SOL']=solvent_number
            else:
                system = solvent_sol( system, atomistic_fragments, residue_type)   
    else:
        atomistic_fragments[residue_type], solvent_number = atomistic_non_protein_non_solvent(residue_type, g_var.cg_residues[residue_type])
        if solvent_number != 0:
            system[residue_type]=solvent_number
        else:
            system = non_solvent(system, atomistic_fragments, residue_type)
    print('{:<100}'.format(''), end='\r')
    print('Finished converting: '+residue_type)
    return system 

def non_solvent(system, atomistic_fragments, residue_type):
#### loop through all resids of that residue type 
    system[residue_type]=len(atomistic_fragments[residue_type])
    
    if not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_all.pdb') and not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb'):
        NP = {}
        count = 0
        coord = []
        index_conversion = {}
        for resid in atomistic_fragments[residue_type]:
            if not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb'):
                atomistic_fragments[residue_type][resid] = at_mod.check_hydrogens(atomistic_fragments[residue_type][resid])
                for at_id, atom in enumerate(atomistic_fragments[residue_type][resid]):
                    if not atomistic_fragments[residue_type][resid][atom]['atom'].startswith('M'):
                        index_conversion[count] = len(coord)
                        coord.append(atomistic_fragments[residue_type][resid][atom]['coord'])
                    NP[count] = atomistic_fragments[residue_type][resid][atom]
                    count+=1
        
        coord=at_mod.check_atom_overlap(coord)
        pdb_output_all = gen.create_pdb(g_var.working_dir+residue_type+'/'+residue_type+'_all.pdb')
        atom_counter = 1
        for at_val, atom in enumerate(NP):
            if at_val in index_conversion:
                x, y, z = gen.trunc_coord(coord[index_conversion[at_val]])
            else:
                x, y, z = gen.trunc_coord([NP[atom]['coord'][0],NP[atom]['coord'][1],NP[atom]['coord'][2]])
            if atom_counter >= 99_999:
                atom_counter=1
            else:
                atom_counter+=1
            pdb_output_all.write(g_var.pdbline%((atom_counter,NP[atom]['atom'],NP[atom]['res_type'],' ',1, x, y, z,0.00, 0.00))+'\n')
    return system

def solvent_sol(system, atomistic_fragments, residue_type):
    #### creates solvent directory and SOL key in system dictionay otherwise it appends solvent molecules to sol pdb
    sol = {}
    coord=[]
    count = 0 
    index_conversion = {}
    if not os.path.exists(g_var.working_dir+'SOL'+'/SOL_all.pdb'):
        pdb_sol = gen.create_pdb(g_var.working_dir+'SOL'+'/SOL_all.pdb')
        for resid in atomistic_fragments[residue_type]:
        ####### check if any atoms in residue overlap #######
            for atom in atomistic_fragments[residue_type][resid]:
                if not atomistic_fragments[residue_type][resid][atom]['atom'].startswith('M'):
                    index_conversion[count] = len(coord)
                    coord.append(atomistic_fragments[residue_type][resid][atom]['coord'])
                sol[count] = atomistic_fragments[residue_type][resid][atom]
                count+=1
        coord=at_mod.check_atom_overlap(coord)
        atom_counter = 1
        for at_val, atom in enumerate(sol):
            if sol[atom]['frag_mass'] > 1:                  
                system['SOL']+=1
            if at_val in index_conversion:
                x, y, z = gen.trunc_coord(coord[index_conversion[at_val]])
            else:
                x, y, z = gen.trunc_coord([sol[atom]['coord'][0],sol[atom]['coord'][1],sol[atom]['coord'][2]])
            if atom_counter >= 99_999:
                atom_counter=1
            else:
                atom_counter+=1
            pdb_sol.write(g_var.pdbline%((atom_counter,sol[atom]['atom'],sol[atom]['res_type'],' ',1,x, y, z,0.00, 0.00))+'\n')
        return system

def solvent_ion(system, atomistic_fragments, residue_type):
    #### creates ion pdb with header
    gen.mkdir_directory(g_var.working_dir+'ION/MIN')
    skip = os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb')
    if not skip:
        pdb_ion = gen.create_pdb(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb')
    for resid in atomistic_fragments[residue_type]:
        for at_id, atom in enumerate(atomistic_fragments[residue_type][resid]):
        #### write ion coordinate out
            if not skip:
                short_line=atomistic_fragments[residue_type][resid][at_id+1]
                x, y, z = gen.trunc_coord([short_line['coord'][0],short_line['coord'][1],short_line['coord'][2]])
                pdb_ion.write(g_var.pdbline%((at_id+1,short_line['atom'],short_line['res_type'],' ',1,x,y,z,0.00, 0.00))+'\n')
            if atomistic_fragments[residue_type][resid][at_id+1]['res_type'] != 'SOL':
                if atomistic_fragments[residue_type][resid][at_id+1]['res_type'] not in system:
                    system[atomistic_fragments[residue_type][resid][at_id+1]['res_type']]=1
                else:
                    system[atomistic_fragments[residue_type][resid][at_id+1]['res_type']]+=1
    return system

def atomistic_non_protein_non_solvent(cg_residue_type,cg_residues):
    atomistic_fragments={}  #### residue dictionary
#### run through every residue in a particular residue type
    residue_type={}
    residue_type_mass={}
    if not os.path.exists(g_var.working_dir+cg_residue_type+'/'+cg_residue_type+'_all.pdb'):
        for cg_resid, cg_residue in enumerate(cg_residues):
            atomistic_fragments[cg_resid]={}
            frag_location=gen.fragment_location(cg_residue_type) ### get fragment location from database
            residue_type[cg_residue_type], residue_type_mass[cg_residue_type] = at_mod.get_atomistic(frag_location)
            for group in residue_type[cg_residue_type]:
                center, at_frag_centers, cg_frag_centers, group_fit = at_mod.rigid_fit(residue_type[cg_residue_type][group], residue_type_mass[cg_residue_type], cg_residue, cg_residues[cg_residue])
                at_connect, cg_connect = at_mod.connectivity(cg_residues[cg_residue], at_frag_centers, cg_frag_centers, group_fit, group)
                if len(at_connect) == len(cg_connect) and len(cg_connect) > 0:
                    try:
                        xyz_rot_apply=at_mod.kabsch_rotate(np.array(at_connect)-center, np.array(cg_connect)-center)
                    except BaseException:
                        sys.exit('There is a issue with residue: '+cg_residue_type+' in group: '+str(group))
                else:
                    print('atom connections: '+str(len(at_connect))+' does not match CG connections: '+str(len(cg_connect)))
                    sys.exit('residue number: '+str(cg_resid)+', residue type: '+str(cg_residue_type)+', group: '+group)
                for bead in group_fit:
                    for atom in group_fit[bead]:
                        group_fit[bead][atom]['coord'] = at_mod.rotate_atom(group_fit[bead][atom]['coord'], center, xyz_rot_apply)   
                        atom_new = group_fit[bead][atom].copy()
                        atomistic_fragments[cg_resid][atom] = atom_new
        return atomistic_fragments, 0
    else:
        return atomistic_fragments, len(cg_residues)   

def atomistic_non_protein_solvent(cg_residue_type,cg_residues): 
    atomistic_fragments={}  #### residue dictionary
#### run through every residue in a particular residue type
    residue_type={}
    residue_type_mass={}

    for cg_resid, cg_residue in enumerate(cg_residues):
        for bead in cg_residues[cg_residue]:
            fragment = bead
            break
        atomistic_fragments[cg_resid]={}
        frag_location=gen.fragment_location(cg_residue_type) ### get fragment location from database
        residue_type[cg_residue_type], residue_type_mass[cg_residue_type] = at_mod.get_atomistic(frag_location)
        if os.path.exists(g_var.working_dir+'SOL'+'/SOL_all.pdb') and cg_residue_type == 'SOL':
            sol_p_bead = 0
            for atom in residue_type_mass[cg_residue_type][g_var.water]:
                if atom[3] > 1:
                    sol_p_bead+=1
            return sol_p_bead, sol_p_bead*len(cg_residues)
        for res_type in residue_type[cg_residue_type]:
            if fragment in residue_type[cg_residue_type][res_type]:
                center, at_frag_centers, cg_frag_centers, group_fit = at_mod.rigid_fit(residue_type[cg_residue_type][res_type], residue_type_mass[cg_residue_type]
                                                                                       , cg_residue, cg_residues[cg_residue])
                xyz_rot_apply=gen.AnglesToRotMat([np.random.uniform(0, math.pi*2), np.random.uniform(0, math.pi*2), np.random.uniform(0, math.pi*2)])
                for bead in group_fit:
                    for atom in group_fit[bead]:
                        group_fit[bead][atom]['coord'] = at_mod.rotate_atom(group_fit[bead][atom]['coord'], center, xyz_rot_apply)   
                        atom_new = group_fit[bead][atom].copy()
                        atomistic_fragments[cg_resid][atom] = atom_new
    return atomistic_fragments, 0



def merge_minimised(residue_type):
    os.chdir(g_var.working_dir+residue_type+'/MIN')
    print('Merging individual residues : '+residue_type)
#### create merged pdb in min folder
    if not os.path.exists(g_var.working_dir+residue_type+'/MIN/'+residue_type+'_merged.pdb'):
        pdb_output=gen.create_pdb(g_var.working_dir+residue_type+'/MIN/'+residue_type+'_merged.pdb')  
        if residue_type =='SOL':
            resid_range=1
        else:
            resid_range=g_var.system[residue_type]
        merge,merge_coords=[],[]
    #### run through every resid 
        for resid in range(resid_range):
            merge_temp ,dump = at_mod.read_in_merged_pdbs([], [], g_var.working_dir+residue_type+'/MIN/'+residue_type+'_'+str(resid)+'.pdb')
            merge, merge_coords = at_mod.fix_chirality(merge,merge_temp,merge_coords, residue_type)    
        if residue_type !='SOL':
            merge_coords = at_mod.check_atom_overlap(merge_coords)
        for line_val, line in enumerate(merge):
            x, y, z = gen.trunc_coord([merge_coords[line_val][0],merge_coords[line_val][1],merge_coords[line_val][2]])
            pdb_output.write(g_var.pdbline%((int(line['atom_number']), line['atom_name'], line['residue_name'],' ',line['residue_id'],\
                x, y, z,1,0))+'\n')
        pdb_output.write('TER\nENDMDL')
        pdb_output.close()

