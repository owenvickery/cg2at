#!/usr/bin/env python
import argparse
import sys
import networkx as nx
import re
import numpy as np
import itertools

_s = re.compile('\s+')
_p = re.compile('(\d+)\s+(\d+)')
_rings_str = """* place dummy atoms at the center of each rings
*

set ncycles = %s

read sequence POL @NCYCLES
generate DUM warn

%s

scalar wmain set 0
%s

return
"""

def build_topology(psffile):
    g = nx.Graph()
    flag = 0

    for line in open(psffile).readlines():
        if flag == 0 and line.strip().endswith('NATOM'):
            flag = 1
            continue
        if flag == 0 and line.strip().endswith('bonds'):
            flag = 2
            continue
        if flag == 1 and not line.strip(): flag = 0
        if flag == 2 and not line.strip(): break

        if flag == 1:
            num, segid, resid, resname, name = _s.split(line)[1:6]
            if resname.startswith('TIP3'): continue
            if name.startswith('H'): continue
            g.add_node(int(num), **{'segid': segid, 'resname': resname, 'name': name, 'resid': resid})

        if flag == 2:
            for pair in _p.findall(line):
                num1, num2 = map(int, pair)
                if g.has_node(num1) and g.has_node(num2): g.add_edge(num1, num2)
    return g

def build_atomtable(psf, crdfile):
    crds = {}
    flag = 0
    for line in open(crdfile).readlines():
        entries = _s.split(line)
        if flag == 0 and entries[0] != '*':
            flag = 1
            continue
        if flag == 1 and not line.strip(): break

        if flag == 1:
            num, resid, resname, name, x, y, z, segid = _s.split(line.strip())[:8]
            if resname.startswith('TIP3'): continue
            if name.startswith('H'): continue
            if psf.node[int(num)]['name'] != name: raise AtomMismatch("%d %s != %d %s" % (int(num), psf.node[int(num)]['name'], int(num), name))
            crds[int(num)] = np.array((float(x), float(y), float(z)))
    return crds

def build_exclusionlist(penstr):
    exclud_list = []
    for line in open(penstr).readlines():
        entries = _s.split(line)
        if entries[0].startswith('%'): continue
        if entries[0].startswith('-'):
            exclud_list.append(entries[1])
    return exclud_list

class AtomMismatch(Exception):
    pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('psf', metavar='psffile', help='PSF file')
    parser.add_argument('crd', metavar='crdfile', help='CRD file')
    parser.add_argument('penstr', metavar='penstrfile', help='pentest.py output STR file')
    args = parser.parse_args()

    # build connectivity of atoms
    psf = build_topology(args.psf)
    crd = build_atomtable(psf, args.crd)
    exclud_list = build_exclusionlist(args.penstr)

    # sanity check
    if len(psf.nodes()) != len(crd): raise AtomMismatch('Number of atom does not match')

    # print out ring atoms
    molecules = nx.connected_component_subgraphs(psf)
    dummys = []
    fixes = []
    i = 1
    _coor = 'coor set xdir %8.3f ydir %8.3f zdir %8.3f sele segid dum .and. resi %d end'
    _fix = 'scalar wmain set 1.0 sele segid %s .and. resi %s end'
    for m in molecules:
        cycles = nx.cycle_basis(m)
        if not cycles: continue
        for cycle in cycles:
            segid = psf.node[cycle[0]]['segid']
            if segid != 'MEMB' and segid in exclud_list: continue
            atoms = np.array([crd[num] for num in cycle])
            com = np.sum(atoms, axis=0)/len(atoms)
            dummys.append(_coor % (com[0], com[1], com[2], i))
            if _fix % (psf.node[cycle[0]]['segid'], psf.node[cycle[0]]['resid']) not in fixes:
                fixes.append(_fix % (psf.node[cycle[0]]['segid'], psf.node[cycle[0]]['resid']))
            i += 1
    print _rings_str % (len(dummys), "\n".join(dummys), "\n".join(fixes))

if __name__ == '__main__':
    main()
