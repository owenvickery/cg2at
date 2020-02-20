#!/usr/bin/env python3

import os, sys
import numpy as np
import subprocess 
import multiprocessing as mp
from shutil import copyfile
from distutils.dir_util import copy_tree
from pathlib import Path
import re
import time
import gen, g_var, f_loc, at_mod


#### collects input structures and creates initial folders
def collect_input(cg, at):
    if not os.path.exists(cg):
        sys.exit('\ncannot find CG input file: '+cg)
    gen.mkdir_directory(g_var.working_dir)
    gen.mkdir_directory(g_var.final_dir)
    gen.mkdir_directory(g_var.input_directory)
    if not g_var.at2cg:
        gen.mkdir_directory(g_var.merged_directory)
#### collates all input files in input directory
    gen.file_copy_and_check(cg, g_var.input_directory+cg.split('/')[-1])
    if at != None:
        if not os.path.exists(at):
            sys.exit('cannot find AT input file: '+at)
        gen.file_copy_and_check(at, g_var.input_directory+at.split('/')[-1])
    os.chdir(g_var.input_directory)
    if cg.split('/')[-1].endswith('.tpr'):
        input_sort(cg, 'conversion')
    else:
        gromacs([g_var.gmx+' editconf -f '+cg.split('/')[-1]+' -resnr 0 -c -o CG_input_temp.pdb', 'CG_input_temp.pdb'])
        gromacs([g_var.gmx+' trjconv -f CG_input_temp.pdb -s '+cg.split('/')[-1]+' -pbc atom -o conversion_input.pdb '+
                        '<< EOF\nSystem\nEOF\n', 'conversion_input.pdb'])

#### converts input files into pdb files 
    if at != None:
        # input_sort(at, 'AT')
        gromacs([g_var.gmx+' editconf -f '+at.split('/')[-1]+' -resnr 0 -o AT_input.pdb', 'AT_input.pdb'])
        return True
    return False

#### fixes pbc complex pbc issue if supplied with a tpr
def input_sort(loc, rep):
    gromacs([g_var.gmx+' make_ndx -f '+loc.split('/')[-1]+' -o '+rep+'_input.ndx << EOF\nq\nEOF\n', rep+'_input.ndx'])
    gromacs([g_var.gmx+' editconf -f '+loc.split('/')[-1]+' -o '+rep+'_input_temp.pdb', rep+'_input_temp.pdb'])
    ndx_f = open(rep+'_input.ndx', "r")
    ndx = ndx_f.read()
    if 'Protein' in ndx:
        gromacs([g_var.gmx+' trjconv -f '+rep+'_input_temp.pdb -s '+loc.split('/')[-1]+' -pbc whole -center -o '+rep+'_input_temp_pbc.pdb -n '+rep+'_input.ndx '+
                        '<< EOF\nProtein\nSystem\nEOF\n', rep+'_input_temp_pbc.pdb'])
    else:
        gromacs([g_var.gmx+' trjconv -f '+rep+'_input_temp.pdb -s '+loc.split('/')[-1]+' -pbc whole -o '+rep+'_input_temp_pbc.pdb -n '+rep+'_input.ndx '+
                        '<< EOF\nSystem\nEOF\n', rep+'_input_temp_pbc.pdb'])        
    gromacs([g_var.gmx+' trjconv -f '+rep+'_input_temp_pbc.pdb -s '+loc.split('/')[-1]+' -pbc atom -o '+rep+'_input_temp_atom.pdb -n '+rep+'_input.ndx '+
                    '<< EOF\nSystem\nEOF\n', rep+'_input_temp_atom.pdb'])
    gromacs([g_var.gmx+' editconf -f '+rep+'_input_temp_atom.pdb -resnr 0 -o '+rep+'_input.pdb -n '+rep+'_input.ndx << EOF\nSystem\nEOF\n', rep+'_input.pdb'])

#### gromacs parser
def gromacs(gro):
    cmd,output = gro[0], gro[1]
    if os.path.exists(output):
        pass
    else:
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
            elif 'Segmentation fault' in out:
                sys.exit('\n'+out)
            elif 'number of atoms in the topology (' in out:
                print('\n'+out+'\n\n')
                print('{0:^90}\n\n{1:^90}\n'.format('***NOTE***','If it is only out by multiples of two, check cysteine distances and increase -cys cutoff'))
                sys.exit('{0:^90}\n\n'.format('A lot of Martini v2-2 disulphide bonds can be up to 10 A (current search cutoff is '+str(g_var.cys)+' A)')) 
            elif 'Fatal error:' in out:
                sys.exit('\n'+out)
    if len(gro) == 4: 
        gro[3].put(gro[2])
        return gro[2]

def gromacs_equilibration(gro):
    cmd,output = gro[0], gro[1]
    if os.path.exists(output):
        return True
    else:
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
            if  'Warning: pressure scaling more than 1%' in out:
                print('pressure coupling failed trying Parrinello-Rahman instead')
                return False
            elif 'Segmentation fault' in out:
                return False
            else:
                return True


def make_min(residue):#, fragments):
#### makes minimisation folder
    gen.mkdir_directory('min')
#### makes em.mdp file for each residue
    if not os.path.exists('em_'+residue+'.mdp'):
        with open('em_'+residue+'.mdp','w') as em:
            em.write('define = \n integrator = steep\nnsteps = 10000\nemtol = 750\nemstep = 0.001\ncutoff-scheme = Verlet\n')

def minimise_protein(protein, p_system, user_at_input, box_vec):
#### makes em.mdp for each chain
    os.chdir(g_var.working_dir+'/PROTEIN')
    gen.folder_copy_and_check(f_loc.forcefield_location+f_loc.forcefield, g_var.working_dir+'PROTEIN/'+f_loc.forcefield+'/.')
    gen.file_copy_and_check(f_loc.forcefield_location+'/residuetypes.dat', 'residuetypes.dat')
    make_min('PROTEIN')
    for chain in range(protein):
        if g_var.v >= 1:
            print('Minimising protein chain: '+str(chain), end='\r')
        pdb2gmx_selections=ask_terminal(chain, p_system)
        pdb2gmx_chain(chain, 'de_novo_', ' << EOF \n1\n'+str(pdb2gmx_selections[0])+'\n'+str(pdb2gmx_selections[1]))
        pdb2gmx_selections = histidine_protonation(chain, 'de_novo_', pdb2gmx_selections)
        at_mod.check_overlap_chain(chain, 'de_novo_', box_vec)
        minimise_protein_chain(chain, 'de_novo_')
        if user_at_input:
            pdb2gmx_chain(chain, 'aligned_', pdb2gmx_selections)
            at_mod.check_overlap_chain(chain, 'aligned_', box_vec)
            minimise_protein_chain(chain, 'aligned_')
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


### interactive terminal residue selection
def ask_ter_question(default_ter, ter_name, ter_val, chain):
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
    return default_ter

def ask_terminal(chain, p_system):
#### default termini is neutral, however if ter flag is supplied you interactively choose termini 
    default_ter=[0,0]
    ter_name=['N terminal','C terminal']
    for ter_val,  ter_residue in enumerate(p_system['terminal_residue'][chain]):
        if not ter_residue:
            if g_var.nt and ter_val==0:
                default_ter[ter_val]=1
            elif g_var.ct and ter_val==1:
                default_ter[ter_val]=1
            elif g_var.ter:
                default_ter = ask_ter_question(default_ter, ter_name, ter_val, chain)          
        else:
            if g_var.ter:
                print('\n The '+ter_name[ter_val]+' residue is non adjustable')
            if ter_val == 0:
                default_ter[ter_val]=3
            else:
                default_ter[ter_val]=4
    return default_ter

def pdb2gmx_chain(chain, input, pdb2gmx_selections):
#### pdb2gmx on on protein chain, creates the topologies    
    gromacs([g_var.gmx+' pdb2gmx -f PROTEIN_'+input+str(chain)+'.pdb -o PROTEIN_'+input+str(chain)+'_gmx.pdb -water none \
    -p PROTEIN_'+input+str(chain)+'.top  -i PROTEIN_'+input+str(chain)+'_posre.itp '+g_var.vs+' -ter '+pdb2gmx_selections+'\nEOF', 'PROTEIN_'+input+str(chain)+'_gmx.pdb']) #### single chains
#### converts the topology file and processes it into a itp file
    convert_topology('PROTEIN_'+input, chain)
#### writes topology overview for each chain 
    write_topol('PROTEIN_'+input, 1, str(chain))
#### writes restraints file for each chain
    write_posres(chain)

def minimise_protein_chain(chain, input):
#### grompps each protein chain
    gromacs([g_var.gmx+' grompp '+
                '-f em_PROTEIN.mdp '+
                '-p topol_PROTEIN_'+input+str(chain)+'.top '+
                '-c PROTEIN_'+input+str(chain)+'_gmx_checked.pdb '+
                '-o min/PROTEIN_'+input+str(chain)+' '+
                '-maxwarn 1 ', 'min/PROTEIN_'+input+str(chain)+'.tpr'])
#### minimises chain
    os.chdir('min')
    gromacs([g_var.gmx+' mdrun -v -nt '+str(g_var.ncpus)+' -deffnm PROTEIN_'+input+str(chain)+' -c PROTEIN_'+input+str(chain)+'.pdb', 'PROTEIN_'+input+str(chain)+'.pdb'])
    os.chdir('..')  


def posres_header(file_write):
        posres_output = open(file_write, 'w')
        posres_output.write('[ position_restraints ]\n; atom  type      fx      fy      fz\n')
        return posres_output

def write_posres(chain):
#### if not posres file exist create one
    steered_posres = posres_header(g_var.working_dir+'PROTEIN/PROTEIN_'+str(chain)+'_steered_posre.itp')
    low_posres = posres_header(g_var.working_dir+'PROTEIN/PROTEIN_'+str(chain)+'_low_posre.itp')
    mid_posres = posres_header(g_var.working_dir+'PROTEIN/PROTEIN_'+str(chain)+'_mid_posre.itp')
    high_posres = posres_header(g_var.working_dir+'PROTEIN/PROTEIN_'+str(chain)+'_high_posre.itp')
    #### read in each chain from after pdb2gmx 
    with open(g_var.working_dir+'PROTEIN/PROTEIN_de_novo_'+str(chain)+'_gmx.pdb', 'r') as pdb_input:
        at_counter=0
        for line in pdb_input.readlines():
            if line.startswith('ATOM'):
                line_sep = gen.pdbatom(line)
                at_counter+=1
            #### if atom is in the restraint list for that residue add to position restraint file
                if line_sep['atom_name'] in f_loc.backbone[line_sep['residue_name']]['posres']:
                    steered_posres.write(str(at_counter)+'     1  1000  1000  1000\n')
                if not gen.is_hydrogen(line_sep['atom_name']):
                    low_posres.write(str(at_counter)+'     1  250  250  250\n')
                    mid_posres.write(str(at_counter)+'     1  1000  1000  1000\n')
                    high_posres.write(str(at_counter)+'     1  6000  6000  6000\n')

def steered_md_atomistic_to_cg_coord(chain):
    os.chdir(g_var.working_dir+'PROTEIN')
    gen.mkdir_directory('steered_md')
#### create bog standard mdp file, simulation is only 3 ps in a vaccum so settings should not have any appreciable effect 
    write_steered_mdp(g_var.working_dir+'PROTEIN/steered_md.mdp', '-DPOSRES_STEERED','Berendsen', 2000, 0.001)
#### run grompp on chain 
    gromacs([g_var.gmx+' grompp '+
                '-f steered_md.mdp '+
                '-p topol_PROTEIN_aligned_'+str(chain)+'.top '+
                '-c min/PROTEIN_aligned_'+str(chain)+'.pdb '+
                '-r min/PROTEIN_de_novo_'+str(chain)+'.pdb '+
                '-o steered_md/PROTEIN_steered_'+str(chain)+' -maxwarn 2 ', 'steered_md/PROTEIN_steered_'+str(chain)+'.tpr'])
#### run mdrun on steered MD
    os.chdir('steered_md')
    gromacs([g_var.gmx+' mdrun -v -nt '+str(g_var.ncpus)+' -pin on -deffnm PROTEIN_steered_'+str(chain)+' -c PROTEIN_steered_'+str(chain)+'.pdb', 
                'PROTEIN_steered_'+str(chain)+'.pdb'])

def convert_topology(topol, protein_number):
#### reads in topology 
    if Path(topol+str(protein_number)+'.top').exists():
        read=False
        if not os.path.exists(topol+str(protein_number)+'.itp'):
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
                        elif line.split()[0] == 'Protein':
                            line= re.sub('Protein', 'protein_'+str(protein_number),line)
                #### writes to itp file copied section          
                    if read:
                        itp_write.write(line)
            #### adds position restraint section to end of itp file         
                itp_write.write('#ifdef POSRES\n#include \"PROTEIN_'+str(protein_number)+'_posre.itp\"\n#endif\n') #_low_posre
                itp_write.write('#ifdef LOWPOSRES\n#include \"PROTEIN_'+str(protein_number)+'_low_posre.itp\"\n#endif\n')
                itp_write.write('#ifdef MIDPOSRES\n#include \"PROTEIN_'+str(protein_number)+'_mid_posre.itp\"\n#endif\n')
                itp_write.write('#ifdef HIGHPOSRES\n#include \"PROTEIN_'+str(protein_number)+'_high_posre.itp\"\n#endif\n')
                itp_write.write('\n; Include CA Position restraint file\n#ifdef POSRES_STEERED\n#include \"PROTEIN_'+str(protein_number)+'_steered_posre.itp\"\n#endif')
    else:
        sys.exit('cannot find : '+topol+'_'+str(protein_number)+'.top')

def write_topol(residue_type, residue_number, chain):
#### open topology file
    found=False
    with open('topol_'+residue_type+chain+'.top', 'w') as topol_write:
    #### add standard headers may need to be changed dependant on forcefield
        topol_write.write('; Include forcefield parameters\n#include \"'+f_loc.forcefield_location+f_loc.forcefield+'/forcefield.itp\"\n')
        if 'SOL' == residue_type:
            topol_write.write('#include \"'+f_loc.water_dir+f_loc.water+'.itp\"\n\n#include \"'+f_loc.forcefield_location+f_loc.forcefield+'/ions.itp\"\n\n')
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
    pool = mp.Pool(g_var.ncpus)
    m = mp.Manager()
    q = m.Queue()
    pool_process = pool.map_async(gromacs, [(g_var.gmx+' grompp '+
                                  '-po md_out-'+residue_type+'_temp_'+str(rid)+' '+
                                  '-f em_'+residue_type+'.mdp '+
                                  '-p topol_'+residue_type+'.top '+
                                  '-c '+residue_type+'_'+str(rid)+'.pdb '+
                                  '-o min/'+residue_type+'_temp_'+str(rid)+' -maxwarn 1', 
                                  'min/'+residue_type+'_temp_'+str(rid)+'.tpr',rid, q) for rid in range(0, resid)])
    while not pool_process.ready():
        report_complete('GROMPP', q.qsize(), resid)
    print('                                                                       ', end='\r')
    print('GROMPP completed on residue type: '+residue_type)       
    pool.close()
    pool.join()
#### close grompp multiprocessing and change to min directory and spin up mdrun multiprocessing
    os.chdir('min')
    m = mp.Manager()
    q = m.Queue()
    pool = mp.Pool(g_var.ncpus)
    pool_process = pool.map_async(gromacs, [(g_var.gmx+' mdrun -v -nt 1 -deffnm '+residue_type+'_temp_'+str(rid)+' -c '+residue_type+'_'+str(rid)+'.pdb', 
                                  residue_type+'_'+str(rid)+'.pdb',rid, q) for rid in range(0, resid)])          ## minimisation grompp parallised  
    while not pool_process.ready():
        report_complete('Minimisation', q.qsize(), resid)
    print('                                                                       ', end='\r')
    print('Minimisation completed on residue type: '+residue_type)
    pool.close()
    os.chdir(g_var.working_dir)

def report_complete(func, size, resid):
    print('                                                                                                       ', end='\r')
    print('Running '+func+' on '+str(resid)+' residues: percentage complete: ',np.round((size/resid)*100,2),'%', end='\r')
    time.sleep(0.1)

def minimise_merged(residue_type, np_system):
#### write topology for merged system
    os.chdir(g_var.working_dir+residue_type)
    write_topol(residue_type, np_system[residue_type], '')
#### grompp with merged system
    gromacs([g_var.gmx+' grompp '+
            '-po md_out-'+residue_type+' '+
            '-f em_'+residue_type+'.mdp '+
            '-p topol_'+residue_type+'.top '+
            '-c '+g_var.working_dir+residue_type+'/min/'+residue_type+'_merged.pdb '+
            '-o '+g_var.working_dir+residue_type+'/min/'+residue_type+'_merged_min -maxwarn 1', g_var.working_dir+residue_type+'/min/'+residue_type+'_merged_min.tpr'])
#### change to min directory and minimise
    os.chdir('min') 
    gromacs([g_var.gmx+' mdrun -v -nt '+str(g_var.ncpus)+' -pin on -deffnm '+residue_type+'_merged_min -c ../'+residue_type+'_merged.pdb', '../'+residue_type+'_merged.pdb'])
    os.chdir(g_var.working_dir)



################################################################ Gromacs for merged system

def write_merged_topol(system, protein):
    os.chdir(g_var.working_dir+'MERGED')
    if not os.path.exists('topol_final.top'):
        with open('topol_final.top', 'w') as topol_write:
        #### writes topology headers (will probably need updating with other forcefields)
            topol_write.write('; Include forcefield parameters\n#include \"'+g_var.final_dir+f_loc.forcefield+'/forcefield.itp\"\n')
            if 'SOL' in system:
                gen.file_copy_and_check(f_loc.water_dir+f_loc.water+'.itp', f_loc.water+'.itp')
                topol_write.write('#include \"'+f_loc.water+'.itp\"')
                topol_write.write('\n#include \"'+g_var.final_dir+f_loc.forcefield+'/ions.itp\"\n\n')
        #### runs through residue types and copies to MERGED directory and simplifies the names
            for residue_type in system:
                if residue_type not in ['ION','SOL']:
                #### copies 1st itp file it comes across 
                    for directory in f_loc.np_directories:
                        if os.path.exists(directory[0]+residue_type+'/'+residue_type+'.itp'):       
                            topol_write.write('#include \"'+residue_type+'.itp\"\n')
                            gen.file_copy_and_check(directory[0]+residue_type+'/'+residue_type+'.itp', residue_type+'.itp')
                            break
                #### copies across protein itp files and simplifies the names 
                    if residue_type == 'PROTEIN':
                        for protein_unit in range(system[residue_type]): 
                            topol_write.write('#include \"PROTEIN_'+str(protein_unit)+'.itp\"\n')
                            gen.file_copy_and_check(g_var.working_dir+'PROTEIN/PROTEIN'+protein+'_'+str(protein_unit)+'.itp', 'PROTEIN_'+str(protein_unit)+'.itp')
                            gen.file_copy_and_check(g_var.working_dir+'PROTEIN/PROTEIN_'+str(protein_unit)+'_steered_posre.itp', 'PROTEIN_'+str(protein_unit)+'_steered_posre.itp')
                            gen.file_copy_and_check(g_var.working_dir+'PROTEIN/PROTEIN_'+str(protein_unit)+'_low_posre.itp', 'PROTEIN_'+str(protein_unit)+'_low_posre.itp')
                            gen.file_copy_and_check(g_var.working_dir+'PROTEIN/PROTEIN_'+str(protein_unit)+'_mid_posre.itp', 'PROTEIN_'+str(protein_unit)+'_mid_posre.itp')
                            gen.file_copy_and_check(g_var.working_dir+'PROTEIN/PROTEIN_'+str(protein_unit)+'_high_posre.itp', 'PROTEIN_'+str(protein_unit)+'_high_posre.itp')
                            gen.file_copy_and_check(g_var.working_dir+'PROTEIN/PROTEIN'+protein+'_'+str(protein_unit)+'_posre.itp', 'PROTEIN_'+str(protein_unit)+'_posre.itp')

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
    gromacs([g_var.gmx+' grompp '+
            '-po md_out-merged_cg2at '+
            '-f em_merged_cg2at.mdp '+
            '-p topol_final.top '+
            '-r merged_cg2at'+protein+'.pdb '+
            '-c merged_cg2at'+protein+'.pdb '+
            '-o min/merged_cg2at'+protein+'_minimised '+
            '-maxwarn 1', 'min/merged_cg2at'+protein+'_minimised.tpr'])
    os.chdir('min')
#### runs minimises final systems
    gromacs([g_var.gmx+' mdrun -nt '+str(g_var.ncpus)+' -v -pin on -deffnm merged_cg2at'+protein+'_minimised -c merged_cg2at'+protein+'_minimised.pdb', 'merged_cg2at'+protein+'_minimised.pdb'])


def alchembed(system, protein_type):
    os.chdir(g_var.working_dir+'MERGED')
    gen.mkdir_directory('alchembed')
#### runs through each chain and run alchembed on each sequentially
    for chain in range(system):
        print('Running alchembed on chain: '+str(chain))
    #### creates a alchembed mdp for each chain 
        if not os.path.exists('alchembed_'+str(chain)+'.mdp'):
            with open('alchembed_'+str(chain)+'.mdp', 'w') as alchembed:
                alchembed.write('define = -DPOSRES\nintegrator = sd\nnsteps = 500\ndt = 0.001\ncontinuation = no\nconstraint_algorithm = lincs')
                alchembed.write('\nconstraints = h-bonds\nns_type = grid\nnstlist = 25\nrlist = 1\nrcoulomb = 1\nrvdw = 1\ncoulombtype  = PME')
                alchembed.write('\npme_order = 4\nfourierspacing = 0.16\ntc-grps = system\ntau_t = 0.1\nref_t = 310\ncutoff-scheme = Verlet')
                alchembed.write('\npcoupl = Berendsen\npcoupltype = semiisotropic\n tau_p = 2.0\nref_p = 1.0 1.0\ncompressibility = 4.5e-5 4.5e-5\n')
                alchembed.write('\npbc = xyz\nDispCorr = no\ngen_vel = yes\ngen_temp = 310\ngen_seed = -1\nfree_energy = yes\ninit_lambda = 0.00')
                alchembed.write('\ndelta_lambda = 1e-3\nsc-alpha = 0.1000\nsc-power = 1\nsc-r-power = 6\ncouple-moltype = protein_'+str(chain))
                alchembed.write('\ncouple-lambda0 = none\ncouple-lambda1 = vdw\nrefcoord_scaling = all')
    #### if 1st chain use minimised structure for coordinate input
        if chain == 0:
            gromacs([g_var.gmx+' grompp '+
                    '-po md_out-merged_cg2at_alchembed_'+str(chain)+' '+
                    '-f alchembed_'+str(chain)+'.mdp '+
                    '-p topol_final.top '+
                    '-r min/merged_cg2at_'+protein_type+'_minimised.pdb '+
                    '-c min/merged_cg2at_'+protein_type+'_minimised.pdb '+
                    '-o alchembed/merged_cg2at_'+protein_type+'_supplied_alchembed_'+str(chain)+' '+
                    '-maxwarn 1', 'alchembed/merged_cg2at_'+protein_type+'_supplied_alchembed_'+str(chain)+'.tpr'])
    #### if not 1st chain use the previous output of alchembed tfor the input of the next chain 
        else:
            gromacs([g_var.gmx+' grompp '+
                '-po md_out-merged_cg2at_alchembed_'+str(chain)+' '+
                '-f alchembed_'+str(chain)+'.mdp '+
                '-p topol_final.top '+
                '-r min/merged_cg2at_'+protein_type+'_minimised.pdb '+
                '-c alchembed/merged_cg2at_'+protein_type+'_supplied_alchembed_'+str(chain-1)+'.pdb '+
                '-o alchembed/merged_cg2at_'+protein_type+'_supplied_alchembed_'+str(chain)+' '+
                '-maxwarn 1', 'alchembed/merged_cg2at_'+protein_type+'_supplied_alchembed_'+str(chain)+'.tpr'])          
        os.chdir('alchembed')
    #### run alchembed on the chain of interest
        gromacs([g_var.gmx+' mdrun -nt '+str(g_var.ncpus)+' -v -pin on -deffnm merged_cg2at_'+protein_type+'_supplied_alchembed_'+str(chain)+
                ' -c merged_cg2at_'+protein_type+'_supplied_alchembed_'+str(chain)+'.pdb', 'merged_cg2at_'+protein_type+'_supplied_alchembed_'+str(chain)+'.pdb'])
        os.chdir('..')
#### copy final output to the FINAL folder
    gen.file_copy_and_check('alchembed/merged_cg2at_'+protein_type+'_supplied_alchembed_'+str(chain)+'.pdb', g_var.merged_directory+'final_cg2at_'+protein_type+'.pdb')
    gen.file_copy_and_check('alchembed/merged_cg2at_'+protein_type+'_supplied_alchembed_'+str(chain)+'.pdb', g_var.final_dir+'final_cg2at_'+protein_type+'.pdb')
    # gen.file_copy_and_check('merged_cg2at_no_steered.pdb', g_var.final_dir+'final_cg2at_no_steered.pdb')

def write_steered_mdp(loc, posres,pc_type, time, timestep):
    if not os.path.exists(loc):
        with open(loc, 'w') as steered_md:
            steered_md.write('define = '+posres+'\nintegrator = md\nnsteps = '+str(time)+'\ndt = '+str(timestep)+'\ncontinuation   = no\nconstraint_algorithm = lincs\n')
            steered_md.write('nstxtcout = 10\nnstenergy = 10\nconstraints = h-bonds\nns_type = grid\nnstlist = 25\nrlist = 1.2\nrcoulomb = 1.2\nrvdw = 1.2\ncoulombtype  = PME\n')
            steered_md.write('pme_order = 4\nfourierspacing = 0.135\ntcoupl = v-rescale\ntc-grps = system\ntau_t = 0.1\nref_t = 310\n')
            steered_md.write('pcoupl = '+pc_type+'\npcoupltype = semiisotropic\ntau_p = 2.0\nref_p = 1.0 1.0\ncompressibility = 4.5e-5 4.5e-5\n')
            steered_md.write('pbc = xyz\nDispCorr = no\ngen_vel = no\nrefcoord_scaling = all\ncutoff-scheme = Verlet')   

def write_nvt_mdp(loc, posres, time, timestep):
    if not os.path.exists(loc):
        with open(loc, 'w') as steered_md:
            steered_md.write('define = '+posres+'\nintegrator = md\nnsteps = '+str(time)+'\ndt = '+str(timestep)+'\ncontinuation   = no\nconstraint_algorithm = lincs\n')
            steered_md.write('constraints = h-bonds\nns_type = grid\nnstlist = 25\nrlist = 1.2\nrcoulomb = 1.2\nrvdw = 1.2\ncoulombtype  = PME\n')
            steered_md.write('pme_order = 4\nfourierspacing = 0.135\ntcoupl = v-rescale\ntc-grps = system\ntau_t = 0.1\nref_t = 310\n')
            steered_md.write('pcoupl = no\npbc = xyz\nDispCorr = no\ngen_vel = yes\ngen_temp = 310\ngen_seed = -1\nrefcoord_scaling = all\ncutoff-scheme = Verlet')   

def reverse_steer(protein_type, fc, input_file ):
    gen.mkdir_directory(g_var.merged_directory+'reverse_steer')
    equil_type = ['Parrinello-Rahman', 'Berendsen']
    for equil_type_val, npt_type in enumerate([fc+'_posres-pr.mdp', fc+'_posres-b.mdp']):
        os.chdir(g_var.merged_directory)
        write_steered_mdp(g_var.merged_directory+npt_type, '-D'+fc.upper()+'POSRES', equil_type[equil_type_val] ,2000, 0.001)  
        gromacs([g_var.gmx+' grompp '+
                ' -po md_out-merged_cg2at_reverse_steer_'+fc+
                ' -t '+input_file+'.cpt '+
                ' -f '+npt_type+' '+
                ' -p topol_final.top '+
                ' -r merged_cg2at_'+protein_type+'.pdb '+
                ' -c '+input_file+'.pdb '+
                ' -o reverse_steer/merged_cg2at_'+protein_type+'_reverse_steer_'+fc+'_'+equil_type[equil_type_val]+' '+
                ' -maxwarn '+str(equil_type_val+2), 'reverse_steer/merged_cg2at_'+protein_type+'_reverse_steer_'+fc+'_'+equil_type[equil_type_val]+'.tpr'])  
        os.chdir('reverse_steer')
        equil = gromacs_equilibration([g_var.gmx+' mdrun -v -nt '+str(g_var.ncpus)+' -pin on -deffnm merged_cg2at_'+protein_type+'_reverse_steer_'+fc+'_'+equil_type[equil_type_val]+
                                     ' -c merged_cg2at_'+protein_type+'_reverse_steer_'+fc+'.pdb -cpo merged_cg2at_'+protein_type+'_reverse_steer_'+fc+'.cpt'
                                     ,'merged_cg2at_'+protein_type+'_reverse_steer_'+fc+'.pdb'])
        if equil:
            break
    if not equil:
        sys.exit('reverse steer failed')

def run_nvt(loc):
    print('Running NVT on de novo system')
    os.chdir(g_var.merged_directory)    
    write_nvt_mdp(g_var.merged_directory+'nvt.mdp', '-DPOSRES', 100, 0.001)
    gen.mkdir_directory(g_var.merged_directory+'NVT')
    gromacs([g_var.gmx+' grompp'+
            ' -po md_out-merged_cg2at_nvt'+
            ' -f nvt.mdp'+
            ' -p topol_final.top'+
            ' -r '+loc+
            ' -c '+loc+
            ' -o NVT/merged_cg2at_de_novo_nvt'+
            ' -maxwarn 2', 'NVT/merged_cg2at_de_novo_nvt.tpr'])   
    os.chdir(g_var.merged_directory+'NVT')
    gromacs([g_var.gmx+' mdrun -v -pin on -nt '+str(g_var.ncpus)+' -deffnm merged_cg2at_de_novo_nvt'+
        ' -c merged_cg2at_de_novo_nvt.pdb', 'merged_cg2at_de_novo_nvt.pdb'])      

def run_npt(input_file):
    print('Running NPT on de novo system')
    os.chdir(g_var.merged_directory)        
    gen.mkdir_directory(g_var.merged_directory+'NPT')
    equil_type = ['Parrinello-Rahman', 'Berendsen']
    for equil_type_val, npt_type in enumerate(['npt-pr.mdp', 'npt-b.mdp']):
        os.chdir(g_var.merged_directory)   
        write_steered_mdp(g_var.merged_directory+npt_type, '-DPOSRES', equil_type[equil_type_val] ,5000, 0.001)
        gromacs([g_var.gmx+' grompp'+
                ' -po md_out-merged_cg2at_npt'+
                ' -f '+npt_type+
                ' -p topol_final.top'+
                ' -t '+input_file+'.cpt '+
                ' -r '+input_file+'.pdb '+
                ' -c '+input_file+'.pdb '+
                ' -o NPT/merged_cg2at_de_novo_npt_'+equil_type[equil_type_val]+
                ' -maxwarn '+str(equil_type_val+2), 'NPT/merged_cg2at_de_novo_npt_'+equil_type[equil_type_val]+'.tpr'])   
        os.chdir(g_var.merged_directory+'NPT')
        equil = gromacs_equilibration([g_var.gmx+' mdrun -v -nt '+str(g_var.ncpus)+' -pin on -deffnm merged_cg2at_de_novo_npt_'+equil_type[equil_type_val]+
                ' -c merged_cg2at_de_novo_npt.pdb -cpo merged_cg2at_de_novo_npt.cpt'
                , 'merged_cg2at_de_novo_npt.pdb'])  
        if equil:
            break
    if not equil:
        sys.exit('NPT run failed')
    gen.file_copy_and_check('merged_cg2at_de_novo_npt.pdb', g_var.final_dir+'final_cg2at_de_novo.pdb') 