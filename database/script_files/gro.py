#!/usr/bin/env python3

import os, sys
import numpy as np
import subprocess 
import multiprocessing as mp
from shutil import copyfile
from distutils.dir_util import copy_tree
from pathlib import Path
import re
import gen, g_var, f_loc, at_mod



def collect_input(cg, at):
    if not os.path.exists(cg):
        sys.exit('\nPrint cannot find CG input file: '+cg)
    gen.mkdir_directory(g_var.working_dir)
    gen.mkdir_directory(g_var.final_dir)
    gen.mkdir_directory(g_var.input_directory)
#### collates all input files in input directory
    copyfile(cg, g_var.input_directory+cg.split('/')[-1])
    if at != None:
        if not os.path.exists(at):
            sys.exit('Print cannot find AT input file: '+at)
        copyfile(at, g_var.input_directory+at.split('/')[-1])
    os.chdir(g_var.input_directory)
#### converts input files into pdb files 
    gromacs(g_var.gmx+' editconf -f '+cg.split('/')[-1]+' -resnr 0 -o CG_input.pdb')
    if at != None:
        gromacs(g_var.gmx+' editconf -f '+at.split('/')[-1]+' -resnr 0 -o AT_input.pdb')
        return True
    return False

def gromacs(cmd):
#### if the flag gromacs is used every gromacs command will be printed to the terminal 
    if g_var.v >= 3:
        print('\nrunning gromacs: \n '+cmd+'\n')
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    err, out = output.communicate()
    exitcode = output.returncode
    out=out.decode("utf-8")
#### all gromacs outputs will be saved into gromacs_outputs within the folder it is run
    with open('gromacs_outputs', 'a') as checks:
        checks.write(out)
#### standard catch for failed gromacs commands
        if 'File input/output error:' in out:
            sys.exit('\n'+out)
        elif 'Error in user input:' in out:
            sys.exit('\n'+out)
        elif 'did not converge to Fmax ' in out:
            sys.exit('\n'+out)
        elif 'Segmentation fault (core dumped):' in out:
            sys.exit('\n'+out)
        elif 'Fatal error:' in out:
            sys.exit('\n'+out)
        elif 'but did not reach the requested Fmax' in out:
            sys.exit('\n'+out)
        elif 'number of atoms in the topology (' in out:
            sys.exit('\n'+out+'\n\nIf it is only 2 atoms out check cysteine distances, and increase -cys cutoff')


def make_min(residue):#, fragments):
#### makes minimisation folder
    gen.mkdir_directory('min')
#### makes em.mdp file for each residue
    if not os.path.exists('em_'+residue+'.mdp'):
        with open('em_'+residue+'.mdp','w') as em:
            em.write('define = \n integrator = steep\nnsteps = 10000\nemtol = 1000\nemstep = 0.001\ncutoff-scheme = Verlet\n')

def minimise_protein(protein, p_system, user_at_input):
#### makes em.mdp for each chain
    os.chdir(g_var.working_dir+'/PROTEIN')
    gen.mkdir_directory(g_var.working_dir+'FORCEFIELD')
    copy_tree(f_loc.forcefield_location+f_loc.forcefield+'.ff', g_var.working_dir+'PROTEIN/'+f_loc.forcefield+'.ff/.')
    copyfile(f_loc.forcefield_location+'/residuetypes.dat', 'residuetypes.dat')
    make_min('PROTEIN')
    for chain in range(protein):
        pdb2gmx_selections=ask_terminal(chain, p_system)
        minimise_protein_chain(chain, 'novo_', ' << EOF \n1\n'+str(pdb2gmx_selections[0])+'\n'+str(pdb2gmx_selections[1]))
        pdb2gmx_selections = histidine_protonation(chain, 'novo_', pdb2gmx_selections)
        if user_at_input:
            minimise_protein_chain(chain, 'at_rep_user_supplied_', pdb2gmx_selections)
    os.chdir('..')

def histidine_protonation(chain, input, chain_ter):
#### reads protonation state of histidine from itp file
    histidines=[]
    with open('PROTEIN_'+input+str(chain)+'.top', 'r') as top_input:
        for line in top_input.readlines():
            if line.startswith('; residue'):
                if line.split()[5] in ['HSD','HID']:
                    histidines.append(0)
                elif line.split()[5] in ['HSE', 'HIE']:
                    histidines.append(1)
                elif line.split()[5] in ['HSP','HIS1']:
                    histidines.append(2)
    pdb2gmx_selections='-his << EOF \n1'
    for his in histidines:
        pdb2gmx_selections+='\n'+str(his)
    pdb2gmx_selections+='\n'+str(chain_ter[0])+'\n'+str(chain_ter[1])
    return pdb2gmx_selections

def ask_terminal(chain, p_system):
#### default termini is neutral, however if ter flag is supplied you interactively choose termini 
    default_ter=[1,1]
    ter_name=['N terminal','C terminal']
    for ter_val,  ter_residue in enumerate(p_system['terminal_residue'][chain]):
        if not ter_residue:
            if g_var.ter:
                print('\n please select species for '+ter_name[ter_val]+' residue in chain '+str(chain)+' :\n 0: charged\n 1: neutral')
                while True:
                    try:
                        number = int(input('\nplease select terminal species: '))
                        if number in [0,1]:
                            default_ter[ter_val]=number
                            break
                    except KeyboardInterrupt:
                        sys.exit('\nInterrupted')
                    except:
                        print("Oops!  That was a invalid choice")
        else:
            if g_var.ter:
                print('\n The '+ter_name[ter_val]+' residue is non adjustable')
            if ter_val == 0:
                default_ter[ter_val]=3
            else:
                default_ter[ter_val]=4
    return default_ter

def minimise_protein_chain(chain, input, pdb2gmx_selections):
#### pdb2gmx on on protein chain, creates the topologies    
    gromacs(g_var.gmx+' pdb2gmx -f PROTEIN_'+input+str(chain)+'.pdb -o PROTEIN_'+input+str(chain)+'_gmx.pdb -water none \
-p PROTEIN_'+input+str(chain)+'.top  -i PROTEIN_'+input+str(chain)+'_posre.itp -ter '+pdb2gmx_selections+'\nEOF') #### single chains
#### converts the topology file and processes it into a itp file
    convert_topology('PROTEIN_'+input, chain)
#### writes topology overview for each chain 
    write_topol('PROTEIN_'+input, 1, str(chain))
#### writes restraints file for each chain
    write_posres(chain)
#### grompps each protein chain
    gromacs(g_var.gmx+' grompp '+
            '-f em_PROTEIN.mdp '+
            '-p topol_PROTEIN_'+input+str(chain)+'.top '+
            '-c PROTEIN_'+input+str(chain)+'_gmx.pdb '+
            '-o min/PROTEIN_'+input+str(chain)+' '+
            '-maxwarn 1 ')
#### minimises chain
    os.chdir('min')
    gromacs(g_var.gmx+' mdrun -v -deffnm PROTEIN_'+input+str(chain)+' -c PROTEIN_'+input+str(chain)+'.pdb')
    os.chdir('..')  

def write_posres(chain):
#### if not posres file exist create one
    if not os.path.exists(g_var.working_dir+'PROTEIN/PROTEIN_'+str(chain)+'_steered_posre.itp'):
        posres_output = open(g_var.working_dir+'PROTEIN/PROTEIN_'+str(chain)+'_steered_posre.itp', 'w')
        posres_output.write('[ position_restraints ]\n; atom  type      fx      fy      fz\n')
    #### read in each chain from after pdb2gmx 
        with open(g_var.working_dir+'PROTEIN/PROTEIN_novo_'+str(chain)+'_gmx.pdb', 'r') as pdb_input:
            at_counter=0
            for line in pdb_input.readlines():
                if line.startswith('ATOM'):
                    line_sep = gen.pdbatom(line)
                    at_counter+=1
                #### if atom is in the restraint list for that residue add to position restraint file
                    if line_sep['atom_name'] in f_loc.backbone[line_sep['residue_name']]['restraint']:
                        posres_output.write(str(at_counter)+'     1  1000  1000  1000\n')

def steered_md_atomistic_to_cg_coord(chain):
    os.chdir(g_var.working_dir+'PROTEIN')
    gen.mkdir_directory('steered_md')
#### create bog standard mdp file, simulation is only 3 ps in a vaccum so settings should not have any appreciable effect 
    with open('steered_md.mdp', 'w') as steered_md:
        steered_md.write('define = -DPOSRES_STEERED\nintegrator = md\nnsteps = 3000\ndt = 0.001\ncontinuation   = no\n')
        steered_md.write('constraint_algorithm = lincs\nconstraints = h-bonds\nns_type = grid\nnstlist = 25\nrlist = 1\n')
        steered_md.write('rcoulomb = 1\nrvdw = 1\ncoulombtype  = PME\npme_order = 4\nfourierspacing = 0.16\ntcoupl = V-rescale\n')
        steered_md.write('tc-grps = system\ntau_t = 0.1\nref_t = 310\npcoupl = no\npbc = xyz\nDispCorr = no\ngen_vel = yes\ngen_temp = 310\ngen_seed = -1')    
#### run grompp on chain 
    gromacs(g_var.gmx+' grompp '+
            '-f steered_md.mdp '+
            '-p topol_PROTEIN_at_rep_user_supplied_'+str(chain)+'.top '+
            '-c min/PROTEIN_at_rep_user_supplied_'+str(chain)+'.pdb '+
            '-r min/PROTEIN_novo_'+str(chain)+'.pdb '+
            '-o steered_md/PROTEIN_at_rep_user_supplied_'+str(chain)+' -maxwarn 1 ')
#### run mdrun on steered MD
    os.chdir('steered_md')
    gromacs(g_var.gmx+' mdrun -v -deffnm PROTEIN_at_rep_user_supplied_'+str(chain)+' -c PROTEIN_at_rep_user_supplied_'+str(chain)+'.pdb')
#### if no pdb file is created stop script with error message
    if os.path.exists('PROTEIN_at_rep_user_supplied_'+str(chain)+'.pdb'):
        pass
    else:
        sys.exit('Steered MD failed!')

def convert_topology(topol, protein_number):
#### reads in topology 
    if Path(topol+str(protein_number)+'.top').exists():
        read=False
        with open(topol+str(protein_number)+'.itp', 'w') as itp_write:
            for line in open(topol+str(protein_number)+'.top', 'r').readlines():
            #### copies between moleculetype and position restraint section
                if len(line.split()) > 1: 
                    if read == False and line.split()[1] == 'moleculetype':
                        read = True
                    if read == True and line.split()[1] == 'POSRES':
                        read = False
                #### sorts out chain naming
                    if line.split()[0][:-1] == 'Protein_chain_':
                        line= re.sub('Protein_chain_.?', 'protein_'+str(protein_number),line)
            #### writes to itp file copied section          
                if read:
                    itp_write.write(line)
        #### adds position restraint section to end of itp file         
            itp_write.write('#ifdef POSRES\n#include \"PROTEIN_'+str(protein_number)+'_posre.itp\"\n#endif\n') 
            itp_write.write('\n; Include CA Position restraint file\n#ifdef POSRES_STEERED\n#include \"PROTEIN_'+str(protein_number)+'_steered_posre.itp\"\n#endif')
    else:
        sys.exit('cannot find : '+topol+'_'+str(protein_number)+'.top')

def write_topol(residue_type, residue_number, chain):
#### open topology file
    found=False
    with open('topol_'+residue_type+chain+'.top', 'w') as topol_write:
    #### add standard headers may need to be changed dependant on forcefield
        topol_write.write('; Include forcefield parameters\n#include \"'+g_var.working_dir+'FORCEFIELD/'+f_loc.forcefield+'.ff/forcefield.itp\"\n')
        if 'SOL' == residue_type:
            topol_write.write('#include \"'+f_loc.water_dir+f_loc.water+'.itp\"\n\n#include \"'+g_var.working_dir+'/FORCEFIELD/'+f_loc.forcefield+'.ff/ions.itp\"\n\n')
    #### add location of residue topology file absolute file locations
        if residue_type not in ['ION','SOL']:
            for directory in range(len(f_loc.np_directories)):
                if os.path.exists(f_loc.np_directories[directory][0]+residue_type+'/'+residue_type+'.itp'):
                    topol_write.write('#include \"'+f_loc.np_directories[directory][0]+residue_type+'/'+residue_type+'.itp\"\n')
                    found=True
                    break
            if os.path.exists(g_var.working_dir+'/PROTEIN/'+residue_type+chain+'.itp'):
                topol_write.write('#include \"'+residue_type+chain+'.itp\"\n')
                found=True
            if not found:
                sys.exit('cannot find itp : '+residue_type+'/'+residue_type+chain)
    #### topology section headers
        topol_write.write('\n\n[ system ]\n; Name\nSomething clever....\n\n[ molecules ]\n; Compound        #mols\n')
    #### individual number of residues
        if residue_type.split('_')[0] == 'PROTEIN':
             residue_type='PROTEIN_'
        topol_write.write(residue_type+chain+'    '+str(residue_number))


#################################################################   Non protein

def non_protein_minimise(resid, residue_type):
#### in the case of SOL all residues are minimised, whilst in all other cases individual residues are minimised separately
    if residue_type != 'SOL':
        individual = 1
        resid=resid
    else:
        individual=resid
        resid=1
    os.chdir(g_var.working_dir+residue_type)
### write topology and minimisation parts (min folder and em.mdp)
    write_topol(residue_type, individual, '')
    make_min(residue_type)#, fragment_names)
#### spin up multiprocessing for grompp 
    pool = mp.Pool(mp.cpu_count())
    pool_process = pool.map_async(gromacs, [(g_var.gmx+' grompp '+
                                  '-po md_out-'+residue_type+'_temp_'+str(rid)+' '+
                                  '-f em_'+residue_type+'.mdp '+
                                  '-p topol_'+residue_type+'.top '+
                                  '-c '+residue_type+'_'+str(rid)+'.pdb '+
                                  '-o min/'+residue_type+'_temp_'+str(rid)+' -maxwarn 1')
                                  for rid in range(0, resid)]).get()          ## minimisation grompp parallised  
    pool.close()
#### close grompp multiprocessing and change to min directory and spin up mdrun multiprocessing
    os.chdir('min')
    pool = mp.Pool(mp.cpu_count())
    pool.map_async(gromacs, [(g_var.gmx+' mdrun -v -nt 1 -deffnm '+residue_type+'_temp_'+str(rid)+' -c '+residue_type+'_'+str(rid)+'.pdb') \
                            for rid in range(0, resid)]).get()
    pool.close()
    os.chdir(g_var.working_dir)



def minimise_merged(residue_type, np_system):
#### write topology for merged system
    os.chdir(g_var.working_dir+residue_type)
    write_topol(residue_type, np_system[residue_type], '')
#### grompp with merged system
    gromacs(g_var.gmx+' grompp '+
            '-po md_out-'+residue_type+' '+
            '-f em_'+residue_type+'.mdp '+
            '-p topol_'+residue_type+'.top '+
            '-c '+g_var.working_dir+residue_type+'/min/'+residue_type+'_merged.pdb '+
            '-o '+g_var.working_dir+residue_type+'/min/'+residue_type+'_merged_min -maxwarn 1')
#### change to min directory and minimise
    os.chdir('min') 
    gromacs(g_var.gmx+' mdrun -v -deffnm '+residue_type+'_merged_min -c ../'+residue_type+'_merged.pdb')
    os.chdir(g_var.working_dir)



################################################################ Gromacs for merged system

def write_merged_topol(system, protein):
    os.chdir(g_var.working_dir+'MERGED')
    with open('topol_final.top', 'w') as topol_write:
    #### writes topology headers (will probably need updating with other forcefields)
        topol_write.write('; Include forcefield parameters\n#include \"'+g_var.working_dir+'FORCEFIELD/'+f_loc.forcefield+'.ff/forcefield.itp\"\n')
        if 'SOL' in system:
            copyfile(f_loc.water_dir+f_loc.water+'.itp', f_loc.water+'.itp')
            topol_write.write('#include \"'+f_loc.water+'.itp\"')
            topol_write.write('\n#include \"'+g_var.working_dir+'/FORCEFIELD/'+f_loc.forcefield+'.ff/ions.itp\"\n\n')
    #### runs through residue types and copies to MERGED directory and simplifies the names
        for residue_type in system:
            if residue_type not in ['ION','SOL']:
            #### copies 1st itp file it comes across 
                for directory in f_loc.np_directories:
                    if os.path.exists(directory[0]+residue_type+'/'+residue_type+'.itp'):       
                        topol_write.write('#include \"'+residue_type+'.itp\"\n')
                        copyfile(directory[0]+residue_type+'/'+residue_type+'.itp', residue_type+'.itp')
                        break
            #### copies across protein itp files and simplifies the names 
                if residue_type == 'PROTEIN':
                    for protein_unit in range(system[residue_type]): 
                        topol_write.write('#include \"PROTEIN_'+str(protein_unit)+'.itp\"\n')
                        copyfile(g_var.working_dir+'PROTEIN/PROTEIN'+protein+'_'+str(protein_unit)+'.itp', 'PROTEIN_'+str(protein_unit)+'.itp')
                        copyfile(g_var.working_dir+'PROTEIN/PROTEIN_'+str(protein_unit)+'_steered_posre.itp', 'PROTEIN_'+str(protein_unit)+'_steered_posre.itp')
                        copyfile(g_var.working_dir+'PROTEIN/PROTEIN'+protein+'_'+str(protein_unit)+'_posre.itp', 'PROTEIN_'+str(protein_unit)+'_posre.itp')

        topol_write.write('[ system ]\n; Name\nSomething clever....\n\n[ molecules ]\n; Compound        #mols\n')
    #### adds number of residues to the topology
        for residue_type in system:
            if residue_type not in  ['PROTEIN']:
                topol_write.write(residue_type+'    '+str(system[residue_type])+'\n')   
        #### adds monomers separately
            if residue_type == 'PROTEIN':
                for protein_unit in range(system[residue_type]):
                    topol_write.write('PROTEIN_'+str(protein_unit)+'    1\n')    

def minimise_merged_pdbs(system, protein):
    print('Minimising merged atomistic files : '+protein[1:])
    os.chdir(g_var.working_dir+'MERGED')
#### grompps final merged systems
    gromacs(g_var.gmx+' grompp '+
            '-po md_out-merged_cg2at '+
            '-f em_merged_cg2at.mdp '+
            '-p topol_final.top '+
            '-r merged_cg2at'+protein+'.pdb '+
            '-c merged_cg2at'+protein+'.pdb '+
            '-o min/merged_cg2at'+protein+'_minimised '+
            '-maxwarn 1')
    os.chdir('min')
#### runs minimises final systems
    gromacs(g_var.gmx+' mdrun -v -deffnm merged_cg2at'+protein+'_minimised -c merged_cg2at'+protein+'_minimised.pdb')


def alchembed(system):
    
    os.chdir(g_var.working_dir+'MERGED')
    gen.mkdir_directory('alchembed')
#### runs through each chain and run alchembed on each sequentially
    for chain in range(system):
        print('Running alchembed on chain: '+str(chain))
    #### creates a alchembed mdp for each chain 
        with open('alchembed_'+str(chain)+'.mdp', 'w') as alchembed:
            alchembed.write('define = -DPOSRES\nintegrator = sd\nnsteps = 500\ndt = 0.001\ncontinuation = no\nconstraint_algorithm = lincs')
            alchembed.write('\nconstraints = h-bonds\nns_type = grid\nnstlist = 25\nrlist = 1\nrcoulomb = 1\nrvdw = 1\ncoulombtype  = PME')
            alchembed.write('\npme_order = 4\nfourierspacing = 0.16\ntc-grps = system\ntau_t = 0.1\nref_t = 310\npcoupl = no\ncutoff-scheme = Verlet')
            alchembed.write('\npbc = xyz\nDispCorr = no\ngen_vel = yes\ngen_temp = 310\ngen_seed = -1\nfree_energy = yes\ninit_lambda = 0.00')
            alchembed.write('\ndelta_lambda = 1e-3\nsc-alpha = 0.1000\nsc-power = 1\nsc-r-power = 6\ncouple-moltype = protein_'+str(chain))
            alchembed.write('\ncouple-lambda0 = none\ncouple-lambda1 = vdw')
    #### if 1st chain use minimised structure for coordinate input
        if chain == 0:
            gromacs(g_var.gmx+' grompp '+
                    '-po md_out-merged_cg2at_alchembed_'+str(chain)+' '+
                    '-f alchembed_'+str(chain)+'.mdp '+
                    '-p topol_final.top '+
                    '-r min/merged_cg2at_at_rep_user_supplied_minimised.pdb '+
                    '-c min/merged_cg2at_at_rep_user_supplied_minimised.pdb '+
                    '-o alchembed/merged_cg2at_at_rep_user_supplied_alchembed_'+str(chain)+' '+
                    '-maxwarn 1')
    #### if not 1st chain use the previous output of alchembed tfor the input of the next chain 
        else:
            gromacs(g_var.gmx+' grompp '+
                '-po md_out-merged_cg2at_alchembed_'+str(chain)+' '+
                '-f alchembed_'+str(chain)+'.mdp '+
                '-p topol_final.top '+
                '-r min/merged_cg2at_at_rep_user_supplied_minimised.pdb '+
                '-c alchembed/merged_cg2at_at_rep_user_supplied_alchembed_'+str(chain-1)+'.pdb '+
                '-o alchembed/merged_cg2at_at_rep_user_supplied_alchembed_'+str(chain)+' '+
                '-maxwarn 1')          
        os.chdir('alchembed')
    #### run alchembed on the chain of interest
        gromacs(g_var.gmx+' mdrun -v -deffnm merged_cg2at_at_rep_user_supplied_alchembed_'+str(chain)+' -c merged_cg2at_at_rep_user_supplied_alchembed_'+str(chain)+'.pdb')
        os.chdir('..')
#### copy final output to the FINAL folder
    copyfile('alchembed/merged_cg2at_at_rep_user_supplied_alchembed_'+str(chain)+'.pdb', g_var.final_dir+'final_cg2at_at_rep_user_supplied.pdb')
    copyfile('merged_cg2at_no_steered.pdb', g_var.final_dir+'final_cg2at_no_steered.pdb')
