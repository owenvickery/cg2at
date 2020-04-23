import os,sys
import argparse


parser = argparse.ArgumentParser(description='Creates posres file for fragments', prog='CG2AT fragments', epilog='Enjoy the program and best of luck!\n')
group_req = parser.add_mutually_exclusive_group()
group_req.add_argument('-f', help='fragment itp file',metavar='POPC.itp',type=str, default='')
group_req.add_argument('-dir', help='apply to all files', action='store_true')

args = parser.parse_args()







def write_posre_file(residue):
    itp_file=[]
    with open(residue+'_posre.itp', 'w') as posre_output:
        posre_output.write('; position restraints for '+residue+'\n[ position_restraints ]\n;  i    funct       fcx        fcy        fcz\n')
        read = False
        with open(residue+'.itp', 'r') as itp_input:
            for line in itp_input.readlines():
                if len(line) > 0:
                    itp_file.append(line)
                    if '[' in line and 'atoms' in line:
                        read = True
                    elif read and '[' in line:
                        read = False
                    elif read and not line.startswith(';') and len(line.split()) >= 7:
                        if int(float(line.split()[7])) > 1:
                            posre_output.write('   {0:5}{1:3}{2:11}{3:11}{4:11}\n'.format(line.split()[0],1,50,50,50))
        print('Created posres file for '+residue)
    return itp_file

def append_ifdef(residue, itp_file):
    if '#ifdef NP\n' not in itp_file:
        with open(residue+'.itp', 'a') as itp_input:
            itp_input.write('\n; Include heavy atom position restraint file\n#ifdef NP\n#include "'+residue+'_posre.itp"\n#endif\n')
        print('Added ifdef flag to '+residue)


if args.dir:
    for folder in os.listdir('.'):
        if not folder.startswith('_') and folder not in ['SOL', 'ION'] and '.' not in folder:
            os.chdir(folder)
            for file in os.listdir('.'):
                if file.endswith('.itp') and 'posre' not in file:
                    print('processing residue '+file)
                    itp_file = write_posre_file(file.split('.')[0])
                    append_ifdef(file.split('.')[0], itp_file)
            os.chdir('..')
else:
    if not os.path.exists(args.f):
        sys.exit('Cannot find '+args.f)
    else:
        itp_file = write_posre_file(args.f.split('.')[0])
        append_ifdef(args.f.split('.')[0], itp_file)    
