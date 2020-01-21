#!/usr/bin/env python3

import os, sys
from time import gmtime, strftime
import distutils.spawn
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description='Converts CG representation into an atomistic representation', prog='CG2AT', epilog='Enjoy the program and best of luck!\n', allow_abbrev=True)

parser.add_argument('-a', help='atomistic coordinates (Optional)',metavar='pdb/gro/tpr',type=str)
parser.add_argument('-loc', help='output folder name, (default = CG2AT_timestamp)',metavar='CG2AT',type=str)
parser.add_argument('-v', action="count", default=0, help="increase output verbosity (eg -vv, 3 levels) (Optional)")
group_C = parser.add_mutually_exclusive_group()
group_N = parser.add_mutually_exclusive_group()
parser.add_argument('-ter', help='interactively choose terminal species (Optional)', action='store_true')
group_N.add_argument('-nt', help='choose charged N terminal', action='store_true')
group_C.add_argument('-ct', help='choose charged C terminal state', action='store_true')
group_N.add_argument('-capN', help='cap N terminal with ACE (Optional) not currently working', action='store_true')
group_C.add_argument('-capC', help='cap C terminal with NME (Optional) not currently working', action='store_true')
parser.add_argument('-clean', help='removes all part files from build', action='store_true')
parser.add_argument('-mod', help='treat fragments individually', action='store_true')
parser.add_argument('-w', help='choose your solvent, common choices are: tip3p, tip4p, spc and spce. This is optional',metavar='tip3p',type=str)
parser.add_argument('-ff', help='choose your forcefield. (Optional)',metavar='charmm36',type=str)
parser.add_argument('-fg', help='choose your fragment library. (Optional)',metavar='martini-2-2',type=str, nargs='*')
parser.add_argument('-gromacs', help='gromacs executable name (Optional)',metavar='gmx_avx',type=str)
parser.add_argument('-cys', help='cutoff for disulphide bonds, sometimes CYS are too far apart (Optional)',metavar='7',type=float, default=7)
parser.add_argument('-swap', help='creates a swap dictionary supply residues as PIP2,D3A:PVCL2,C3A (Optional)',metavar='PIP2,D3A:PVCL2,C3A',type=str, nargs='*')
parser.add_argument('-box', help='box size in Angstrom (0 = use input file) (Optional)',metavar='100',type=float, nargs=3)
parser.add_argument('-vs', help='use virtual sites', action='store_true')
parser.add_argument('-al', help='switches off alchembed (WARNING may cause issues with lipids and rings)', action='store_false')
parser.add_argument('-sf', help='scale factor for fragments, shrinks fragments before fitting to CG',metavar='0.9',type=float, default=0.9)
parser.add_argument('-at2cg', help='converts atomistic to coarsegrain ', action='store_true')
parser.add_argument('-version', action='version', version='%(prog)s 2.0')
group_req = parser.add_mutually_exclusive_group()
group_req.add_argument('-c', help='coarse grain coordinates',metavar='pdb/gro/tpr',type=str)
group_req.add_argument('-info', help=' provides version, available forcefields and fragments', action='store_true')
parser.add_argument('-o', help='Final output supplied (default = all)', default='all', type=str, choices= ['all', 'align', 'steer', 'none'])
args = parser.parse_args()
options = vars(args)

#### if missing structure file print help and quit
if not args.info and args.c == None:
    help_output = parser.print_help(sys.stderr)
    print(parser.print_help(sys.stderr), '\n')
    sys.exit('Error: the following arguments are required: -c\n')


# convert argparser into global variables to be read by the other files

#### run atomistic to coarsegrain
at2cg=args.at2cg

# input/output files  
c, a, o = args.c, args.a, args.o
# forcfield and fragment inputs
w, ff, fg, mod = args.w, args.ff, args.fg, args.mod
cys, swap = args.cys, args.swap
ter, nt, ct, capN, capC = args.ter, args.nt, args.ct, args.capN, args.capC
alchembed =  args.al
#### virtual site information
vs=args.vs
if args.vs:
    vst='0.004'
    vs = '-vsite h'
    sf = args.sf-0.1
else:
    vst='0.002'
    vs = ''
    sf=args.sf
box = args.box
v, clean  = args.v, args.clean

#### print information and quit
info = args.info

### hardcoded variables for use elsewhere in the script

variables_to_save={'-c':c,'-a':a, '-w':w, '-ff':ff, '-fg':fg, '-mod':mod, 
                   '-cys':cys, '-swap':swap, '-ter':ter, '-nt':nt, '-ct':ct, 
                   '-vs':args.vs, '-box':box,'-loc':args.loc, '-at2cg':args.at2cg, '-al':args.al}
topology = ['frag', 'group', 'C_ter', 'N_ter', 'posres', 'sul']
box_line="CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00 P 1           1\n"
pdbline = "ATOM  %5d %4s %4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f"
mass = {'H': 1,'C': 12,'N': 14,'O': 16,'P': 31,'M': 0, 'B': 32 ,'S': 32} 
aas = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 
       'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 
       'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

### CG2AT folder locations

timestamp       =  strftime("%Y-%m-%d_%H-%M-%S", gmtime())

if args.loc != None:
    working_dir_name = args.loc
else:
    if args.at2cg:
        working_dir_name =  'AT2CG_'+timestamp
    else:
        working_dir_name =  'CG2AT_'+timestamp

start_dir       = os.getcwd()+'/'  ### initial working directory
working_dir     = os.getcwd()+'/'+working_dir_name+'/'   ### working directory 
final_dir       = os.getcwd()+'/'+working_dir_name+'/FINAL/'  ### final directory for run files
input_directory = os.getcwd()+'/'+working_dir_name+'/INPUT/'  ### contains input run files
merged_directory = os.getcwd()+'/'+working_dir_name+'/MERGED/'  ### contains input run files
scripts_dir     = os.path.dirname(os.path.realpath(__file__))+'/' ### contains script files
database_dir    = str(Path(*Path(scripts_dir).parts[:-1]))+'/' ### contains database files


### finds gromacs installation


if args.gromacs != None:
    gmx=distutils.spawn.find_executable(args.gromacs)
else:
    gmx=distutils.spawn.find_executable('gmx')
if gmx==None or type(gmx) != str:
    for root, dirs, files in os.walk(os.environ.get("GMXBIN")):
        for file_name in files:
            if file_name.startswith('gmx') and file_name.islower() and '.' not in file_name:
                gmx=distutils.spawn.find_executable(file_name)
                if type(gmx) == str and gmx != None :
                    break
                else:
                    gmx=None
        break
    if gmx == None:
        sys.exit('Cannot find gromacs installation')
