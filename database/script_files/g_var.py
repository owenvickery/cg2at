#!/usr/bin/env python3

import os, sys
from time import gmtime, strftime
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description='Converts CG representation into an atomistic representation', prog='CG2AT', epilog='Enjoy the program and best of luck!\n')
group_req = parser.add_mutually_exclusive_group()
group_req.add_argument('-info', help=' provides version, available forcefields and fragments', action='store_true')
group_req.add_argument('-c', help='coarse grain coordinates',metavar='pdb/gro/tpr',type=str)
parser.add_argument('-a', help='atomistic coordinates (Optional)',metavar='pdb/gro/tpr',type=str, nargs='*')
parser.add_argument('-d', help='duplicate atomistic chains. (0:3 1:3 means 3 copies each of chain 0 and 1)',type=str, nargs='*', default=[])
parser.add_argument('-group', help='treat user supplied atomistic chains, as rigid bodies. (0,1 2,3 or all or chain)',type=str, nargs='*')
parser.add_argument('-loc', help='output folder name, (default = CG2AT_timestamp)',metavar='CG2AT',type=str)
parser.add_argument('-o', help='Final output supplied (default = all)', default='all', type=str, choices= ['all', 'align', 'de_novo', 'none'])
parser.add_argument('-w', help='choose your solvent, common choices are: tip3p, tip4p, spc and spce. This is optional',metavar='tip3p',type=str)
parser.add_argument('-ff', help='choose your forcefield. (Optional)',metavar='charmm36',type=str)
parser.add_argument('-fg', help='choose your fragment library. (Optional)',metavar='martini-2-2',type=str, nargs='*')
parser.add_argument('-mod', help='treat fragments individually', action='store_true')
parser.add_argument('-swap', help='creates a swap dictionary supply residues as PIP2,D3A:PVCL2,C3A (Optional)',metavar='PIP2,D3A:PVCL2,C3A',type=str, nargs='*')
parser.add_argument('-v', action="count", default=0, help="increase output verbosity (eg -vv, 3 levels) (Optional)")
group_C = parser.add_mutually_exclusive_group()
group_N = parser.add_mutually_exclusive_group()
parser.add_argument('-ter', help='interactively choose terminal species (Optional)', action='store_true')
group_N.add_argument('-nt', help='choose neutral N terminal state', action='store_true')
group_C.add_argument('-ct', help='choose neutral C terminal state', action='store_true')
parser.add_argument('-messy', help='do not remove part files CG2AT', action='store_true')
parser.add_argument('-gromacs', help='gromacs executable name (Optional)',metavar='gmx_avx',type=str)
parser.add_argument('-cys', help='cutoff for disulphide bonds, sometimes CYS are too far apart (Optional)',metavar='7',type=float, default=7)
parser.add_argument('-silent', help='silent cysteines question', action='store_true')
parser.add_argument('-box', help='box size in Angstrom (0 = use input file) (Optional)',metavar='100',type=float, nargs=3)
parser.add_argument('-vs', help='use virtual sites', action='store_true')
parser.add_argument('-al', help='switches off alchembed (WARNING may cause issues with lipids and rings)', action='store_false')
parser.add_argument('-sf', help='scale factor for fragments, shrinks fragments before fitting to CG',metavar='0.9',type=float, default=0.9)
parser.add_argument('-version', action='version', version='%(prog)s 2.0')
parser.add_argument('-ncpus', help='maximum number of cores to use (default = all)', type=int)
parser.add_argument('-disre', help='switches off the distance restraint matrix for the backbone', action='store_false')
parser.add_argument('-ov', help='amount of overlap allowed between atoms', type=float, default=0.14)

args = parser.parse_args()
opt = vars(args)
opt['input']=''.join([ i+' ' for i in sys.argv])+'\n'

#### if missing structure file print help and quit
if not args.info and args.c == None:
    help_output = parser.print_help(sys.stderr)
    print(parser.print_help(sys.stderr), '\n')
    sys.exit('Error: the following arguments are required: -c\n')

# convert argparser into global variables to be read by the other files

# input/output files  
c = args.c
a = args.a
o = args.o

# forcfield and fragment inputs
w = args.w
ff = args.ff
fg = args.fg
mod = args.mod

# modifiers
ov = args.ov
cys = args.cys
silent = args.silent
swap = args.swap
group = args.group
duplicate = args.d
disre = args.disre
ter = args.ter
nt = args.nt
ct = args.ct
alchembed = args.al
ncpus = args.ncpus
gmx = args.gromacs
box = args.box
v = args.v
messy  = args.messy
info = args.info

#### virtual site information
if args.vs:
    vs = '-vsite h'
    sf = args.sf-0.1
else:
    vs = ''
    sf=args.sf

### hardcoded variables for use elsewhere in the script

topology = {'C_TERMINAL':'', 'N_TERMINAL':'', 'STEER':[], 'CHIRAL':{'atoms':[]}, 'GROUPS':{'group_max':1}, 'CONNECT':{'atoms':{}}}

box_line="CRYST1 %8.3f %8.3f %8.3f %8.2f %8.2f %8.2f P 1           1\n"

pdbline = "ATOM  %5d %4s %4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f"

cg_water_types = ['W', 'SOL', 'WN', 'WF', 'PW']

aas = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 
       'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 
       'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

termini_selections = {'charmm':{'N_TERMINAL':{'PRO':{'NH':1,'NH3+':2, '5TER':3, 'NONE':4},'NORM':{'NH3+':0,'NH2':1,'5TER':2, 'NONE':3},}, 
                                'C_TERMINAL':{'NORM':{'COO-':0,'COOH':1,'CT2':2, '3TER':3, 'NONE':4}}},
                     'opls':{'N_TERMINAL':{'PRO':{'NH':2, 'NH3+':3, 'NONE':5},'NORM':{'NH3+':0,'NH2':2, 'NONE':3}}, 
                                'C_TERMINAL':{'NORM':{'COO-':0,'COOH':2, 'NONE':3}}}
}
### CG2AT folder locations

timestamp =  strftime("%Y-%m-%d_%H-%M-%S", gmtime())

if args.loc != None:
    working_dir_name = args.loc
else:
    working_dir_name =  'CG2AT_'+timestamp




start_dir       = os.getcwd()+'/'  ### initial working directory
working_dir     = os.getcwd()+'/'+working_dir_name+'/'   ### working directory 
final_dir       = os.getcwd()+'/'+working_dir_name+'/FINAL/'  ### final directory for run files
input_directory = os.getcwd()+'/'+working_dir_name+'/INPUT/'  ### contains input run files
merged_directory = os.getcwd()+'/'+working_dir_name+'/MERGED/'  ### contains run files
scripts_dir     = os.path.dirname(os.path.realpath(__file__))+'/' ### contains script files
database_dir    = str(Path(*Path(scripts_dir).parts[:-1]))+'/' ### contains database files
box_vec = ''
user_at_input = False
p_system = {} ## contains the chain termini info e.g. if chain has a non standard temini 0:[True, False]
o_system = {}
system = {}  ## number of system components e.g. PROTEIN:2 POPE:10, POPG:20
backbone_coords = {} ## CG coordinates of the backbone beads
coord_atomistic = {} ## de_novo atomisitic information e.g. coord_atomistic[chain_count][residue_number][atom][info....]
user_cys_bond = {} ## contains resid of disulphide bonds e.g. user_cys_bond[chain][[cys_resid,cys_resid], [cys_resid,cys_resid]])
cg_residues = {} ## dictionary of CG beads eg cg_residues[residue type(POPE)][resid(1)][bead name(BB)][residue_name(PO4)/coordinates(coord)]
seq_cg = {} ## CG sequence e.g. seq_cg[chain][sequence]
seq_at = {} ## user AT sequence e.g. seq_at[chain][sequence]
seq_cg_other = {} ## CG seq for linked non protein residues e.g. DNA
tc = {} ## contains script timings
atomistic_protein_input_raw = {} ## Raw user atomistic info coord_atomistic[chain_count][residue_number][atom][info....]
atomistic_protein_input_aligned = {}
chain_count = 0 ## number of user atomistic chains
other_atomistic ={}
forcefield_available, fragments_available = '',''
forcefield_location, forcefield = '','' ## forcefield info
np_residues, p_residues, mod_residues, o_residues, np_directories, p_directories, mod_directories, o_directories = [],[],[],[],[],[],[],[]  ## fragment info
database_locations = [] ## grouped directories
water_dir, water = [],[] ## water information
swap_dict ={} ## CG residue swap
res_top, sorted_connect, hydrogen, heavy_bond, ions, at_mass = {},{},{},{},[],{}