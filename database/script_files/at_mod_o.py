#!/usr/bin/env python3

import os, sys
import numpy as np
import copy
import gen, g_var, at_mod, read_in

def build_non_protein_linked_atomistic_system(cg_residues):   
#### initisation of counters
    chain_count=0
    g_var.seq_cg_other[chain_count] = []
    # g_var.o_system[chain_count]=[]
    g_var.o_system[chain_count]=[False,False]
    g_var.other_atomistic[chain_count]={}
    print('Converting other linked residues')
    gen.mkdir_directory(g_var.working_dir+'OTHER')  ### make and change to protein directory
#### for each residue in protein
    residue_type={}
    residue_type_mass={}
    for cg_residue_id, residue_number in enumerate(cg_residues):
        resname = cg_residues[residue_number][next(iter(cg_residues[residue_number]))]['residue_name']
        if cg_residue_id == 0:
            g_var.o_system[chain_count][0]=resname
            # g_var.o_system[chain_count].append(terminal_residue(resname)[0]) 
        g_var.other_atomistic[chain_count][residue_number]={}
        frag_location=gen.fragment_location(resname) ### get fragment location from database
        residue_type[resname], residue_type_mass[resname] = at_mod.get_atomistic(frag_location)
        g_var.seq_cg_other = add_to_sequence(g_var.seq_cg_other, resname, chain_count)
        new_chain = False
        for group in residue_type[resname]:
            skip = False
            if next(iter(residue_type[resname][group])) in g_var.res_top[resname]['CONNECT'] and next(iter(residue_type[resname][group])) not in cg_residues[residue_number]:
                next(iter(residue_type[resname][group]))
                skip = True
            if not skip:
                center, at_frag_centers, cg_frag_centers, group_fit = at_mod.rigid_fit(residue_type[resname][group], residue_type_mass[resname], residue_number, cg_residues[residue_number])
                at_connect, cg_connect = at_mod.connectivity(cg_residues[residue_number], at_frag_centers, cg_frag_centers, group_fit, group)
                for group_bead in group_fit:
                    if group_bead in g_var.res_top[resname]['CONNECT']:
                        at_connect, cg_connect, new_chain = at_mod.BB_connectivity(at_connect,cg_connect, cg_residues, group_fit[group_bead], residue_number, group_bead)
                        # g_var.backbone_coords[chain_count].append(np.append(cg_residues[residue_number][group_bead]['coord'], 1)) 
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
                        g_var.other_atomistic[chain_count][residue_number][atom] = atom_new
        if new_chain:
            # g_var.o_system[chain_count].append(terminal_residue(resname)[1])
            g_var.o_system[chain_count][1]=resname
            chain_count+=1
            if cg_residue_id+1 != len(cg_residues):
                g_var.other_atomistic[chain_count]={}
                # g_var.o_system[chain_count]=[]
                # g_var.o_system[chain_count].append(terminal_residue(cg_residues[residue_number+1][next(iter(cg_residues[residue_number+1]))]['residue_name'])[0])
                g_var.o_system[chain_count]=[resname, False]
                # g_var.o_system[chain_count]=[False, False]
                g_var.seq_cg_other[chain_count]=[]
    if g_var.v >=1:
        print('Sequence of other linked non protein chains')
        print('\n{0:^15}{1:^12}'.format('chain number', 'length of chain')) #   \nchain number\tDelta A\t\tno in pdb\tlength of chain')
        print('\n{0:^15}{1:^12}'.format('------------', '---------------'))
        for chain in g_var.seq_cg_other:
            print('{0:^15}{1:^12}'.format(chain, len(g_var.seq_cg_other[chain])))
        print()
    g_var.system['OTHER']=chain_count
    g_var.system['OTHER']

def add_to_sequence(sequence, residue, chain_count):
    if residue  not in g_var.o_residues:
        sequence[chain_count]+=residue
    else:
        sequence[chain_count]+='X'    
    return sequence    

def finalise_novo_atomistic():
    final_at_residues={}
    final_at = {}
    for chain in g_var.other_atomistic: 
        at_number=0    
        final_at_residues[chain]={}  
        final_at[chain]={}
        coords=[]
        skip = os.path.exists(g_var.working_dir+'OTHER/OTHER_de_novo_'+str(chain)+'.pdb')
        if not skip:
            pdb_output = gen.create_pdb(g_var.working_dir+'OTHER/OTHER_de_novo_'+str(chain)+'.pdb')
        for res_index, residue_id in enumerate(g_var.other_atomistic[chain]):
            g_var.other_atomistic[chain][residue_id] = at_mod.check_hydrogens(g_var.other_atomistic[chain][residue_id])
            offset = min(g_var.other_atomistic[chain][residue_id].keys())
            for at_val, atom in enumerate(g_var.other_atomistic[chain][residue_id]):
                if at_val+offset in g_var.other_atomistic[chain][residue_id]:
                    g_var.other_atomistic[chain][residue_id][at_val+offset]['resid'] = res_index
                    coords.append(g_var.other_atomistic[chain][residue_id][at_val+offset]['coord'])
                    final_at[chain][at_number]=g_var.other_atomistic[chain][residue_id][at_val+offset]
                    at_number+=1
                else:
                    print(residue_id,at_val+1,min(g_var.other_atomistic[chain][residue_id].keys())) #'\n',g_var.other_atomistic[chain][residue_id])
                    sys.exit()
                
            final_at_residues[chain][res_index]=g_var.other_atomistic[chain][residue_id]
        coords = at_mod.check_atom_overlap(coords)
        for atom in final_at[chain]:
            final_at[chain][atom]['coord']=coords[atom]
            final_at_residues[chain][final_at[chain][atom]['resid']][atom]=final_at[chain][atom]
            if not skip:
                x, y, z = gen.trunc_coord(final_at[chain][atom]['coord'])
                pdb_output.write(g_var.pdbline%((atom,final_at[chain][atom]['atom'],final_at[chain][atom]['res_type'],' ',\
                                final_at[chain][atom]['resid'],x,y,z,1,0))+'\n')
    if 'pdb_output' in locals():
        pdb_output.close()
    return final_at_residues

def terminal_residue(resname):
    ter = ['5TER', '3TER']
    if g_var.res_top[resname]['N_TERMINAL'] in ['5TER', 'None']:
        ter[0] = g_var.res_top[resname]['N_TERMINAL']
    if g_var.res_top[resname]['C_TERMINAL'] in ['3TER', 'CT2', 'None']:
        ter[1] = g_var.res_top[resname]['C_TERMINAL']
    return ter