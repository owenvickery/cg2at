import sys
import numpy as np
import copy

atoms=[]
with open('modern', 'r') as pdb_input:
    for line in pdb_input.readlines():
        if len(line) > 0 and not line.startswith(';'):
            line_sep=line.split()
            # print(line)
            atoms.append([line_sep[0], line_sep[1], line_sep[2],line_sep[3]])

atoms=np.array(atoms)

with open('to_add', 'r') as pdb_input:
    for line_val, line in enumerate(pdb_input.readlines()):
        if len(line) > 0:
            line_val+=1
            if not line.startswith(';'):
                line_sep=line.split()
                skip = False
                for at_val, at in enumerate(atoms):
                    # if line_sep[0] == at[0] and line_sep[1] == at[1] and line_sep[2] == at[2]:
                    #     print(at_val, at,'\n',line_val,line)
                    #     skip=True
                    # elif line_sep[0] == at[2] and line_sep[1] == at[1] and line_sep[2] == at[0]:
                    #     print(at_val, at,'\n',line_val,line)
                    #     skip=True

                    if line_sep[0] == at[0] and line_sep[1] == at[1] and line_sep[2] == at[2] and line_sep[3] == at[3]:
                        skip=True
                        l=[at_val, at,'\n',line_val,line]
                        # print(at_val, at,'\n',line_val,line)
                    elif line_sep[0] == at[3] and line_sep[1] == at[2] and line_sep[2] == at[1] and line_sep[3] == at[0]:
                        l = [at_val, at,'\n',line_val,line]
                        # print(at_val, at,'\n',line_val,line)
                        skip=True
                if not skip:
                    print(line)

