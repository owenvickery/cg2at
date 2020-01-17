# Procedure of DMPI topology generated with custom PI head group
## Forcefield modification
To paramertise the phosphatidylinositol (PI) head group, the file lipids.rtp in
the forcefield folder (lipid17.ff) has to be modified to include the entry for
PI. The file `residuetypes.dat` also need to modified where a line for the PI
head group has to be added to allow `pdb2gmx` to recognise the head group correctly.

## Lipid topology generation
After the relevant filed are modified. The topology for the new file can be
generated using the `pdb2gmx` function of Gromacs. The `GMXLIB` has to be set to
the current directory first.

```bash
cat lipids.rtp >> ../lipid17.ff/lipids.rtp
cat ../residuetypes.dat >> residuetypes.dat
export GMXLIB=$PWD/../
gmx pdb2gmx -f DMPI.pdb -o DMPI.gro -ff lipid17
```

## Charmm lipid conversion

The script `charmmlipid2amber.py` has been modified to allow PI head group
to be converted from charmm name to amber name.

