#!/usr/local/bin/python3
# -*- coding: latin-1 -*-
"""Interconvert PDB internal and external coordinates.

Interconvert PDB Structure data between external (X, Y, Z cartesian)
coordinates and internal (bond length, angle and dihedral angle) measurements.
"""

#
# replicate buildprot with biopython
#

import argparse
import os
import sys
import re

# print(sys.path)

import gzip

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO

from InternalCoordinates import read_PIC, PIC_Chain, report_PIC
from InternalCoordinates import internal_to_atom_coordinates, write_PIC


PDB_repository_base = None

if os.path.isdir('/media/data'):
    PDB_repository_base = '/media/data/pdb/'
elif os.path.isdir('/Volumes/data'):
    PDB_repository_base = '/Volumes/data/pdb/'


arg_parser = argparse.ArgumentParser(
    description='Interconvert .pic (pprotein internal coordinates) and '
                '.pdb (protein data bank) files.')
arg_parser.add_argument(
    'file', nargs='*',
    help='a .pdb or .pic path/filename to read (first model, first chain), '
         'or a PDB idCode with optional chain ID to read from {0} as .ent.gz'
         .format((PDB_repository_base or '[PDB resource not defined - '
                  'please configure before use]')))

arg_parser.add_argument(
    '-f', dest='filelist',
    help='a Dunbrack cullPDB pdb ID list to read from {0} as .ent.gz'
         .format((PDB_repository_base or '[PDB resource not defined - '
                  'please configure before use]')))
arg_parser.add_argument('-wp', help='write pdb file with .pypdb extension',
                        action="store_true")
arg_parser.add_argument('-wi', help='write pic file with .pypdb extension',
                        action="store_true")
arg_parser.add_argument('-ct', help='test conversion supplied pdb/pic to ',
                        action="store_true")

args = arg_parser.parse_args()


print(args)

toProcess = args.file
pdbidre = re.compile(r'(^\d(\w\w)\w)(\w)?$')

if args.filelist:
    flist = open(args.filelist, 'r')
    for aline in flist:
        fields = aline.split()
        pdbidMatch = pdbidre.match(fields[0])
        if pdbidMatch:
            # print(m.group(1) + ' ' + m.group(2))
            # toProcess.append(PDB_repository_base + m.group(2)
            # + '/pdb' + m.group(1) + '.ent.gz' )
            toProcess.append(pdbidMatch.group(0))


if len(toProcess):
    print(len(toProcess), 'entries to process')
else:
    print("no files to process. use '-h' for help")
    sys.exit(0)

PDB_parser = PDBParser(PERMISSIVE=False, QUIET=False)

for target in toProcess:
    pdb_input = False
    pic_input = False
    pdb_structure = None
    # pdb_chain = None
    prot_id = ''

    pdbidMatch = pdbidre.match(target)
    if pdbidMatch is not None:
        assert PDB_repository_base, 'PDB repository base directory missing, '
        'please configure for this host'
        pdb_input = True
        filename = (PDB_repository_base + pdbidMatch.group(2).lower() + '/pdb'
                    + pdbidMatch.group(1).lower() + '.ent.gz')
        prot_id = pdbidMatch.group(1)
    else:
        filename = target
        prot_id = target
        if '.' in prot_id:
            prot_id = prot_id[0:prot_id.find('.')]
        if '/' in prot_id:
            prot_id = prot_id[prot_id.rfind('/')+1:]
        if 'pdb' in prot_id:
            prot_id = prot_id[prot_id.rfind('pdb') + 3:]

    if not pdb_input:
        pdb_structure = read_PIC(
            gzip.open(filename, mode='rt')
            if filename.endswith('.gz') else filename)
        if pdb_structure is not None:
            pic_input = True

    if pdb_structure is None:
        pdb_input = True
        pdb_structure = PDB_parser.get_structure(
            prot_id,
            gzip.open(filename, mode='rt')
            if filename.endswith('.gz') else filename)

    # get specified chain if given, else just pick first for now
    # count atoms to detect pic file if don't know already

    rCount = 0
    aCount = 0

    if pdbidMatch is not None and pdbidMatch.group(3) is not None:
        # have chain specifier for PDBid
        if pdb_structure[0][pdbidMatch.group(3)] is not None:
            pdb_chain = pdb_structure[0][pdbidMatch.group(3)]
            rCount = len(pdb_chain)
        else:
            print('chain ' + pdbidMatch.group(3) + ' not found in ' + filename)
            continue

    if pdb_input:
        print('parsed pdb input', prot_id, filename)
    #    for res in pdb_chain.get_residues():   # pdb_structure.get_residues():
    #        print(res.get_full_id(), res.resname,
    #              'disordered' if res.disordered else '')
        for chn in pdb_structure.get_chains():
            chn.pic = PIC_Chain(chn)
            chn.pic.dihedra_from_atoms()
    else:
        print('parsed pic input ', filename)
        print(report_PIC(pdb_structure))
        internal_to_atom_coordinates(pdb_structure)

    # print(pdb_structure.header['idcode'], pdb_chain.id, ':',
    #      pdb_structure.header['head'])

    if args.wp:
        io = PDBIO()
        io.set_structure(pdb_structure)
        io.save(target + '.PyPDB')
        print('wrote pdb output for', target)

    if args.wi:
        write_PIC(pdb_structure, target + '.PyPIC')
        print('wrote pic output for', target)
print('normal termination')
