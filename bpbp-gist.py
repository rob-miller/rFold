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
from io import StringIO
# import itertools
# print(sys.path)

import gzip
import warnings

from Bio.PDB.PDBParser import PDBParser

from Bio.PDB.pic import read_PIC, report_PIC, structure_rebuild_test
from Bio.PDB.pic import internal_to_atom_coordinates, write_PIC
from Bio.PDB.pic import atom_to_internal_coordinates, PIC_Residue
from Bio.PDB.pic import write_PDB, write_SCAD, PIC_Chain

PDB_repository_base = None

if os.path.isdir('/media/data/pdb'):
    PDB_repository_base = '/media/data/pdb/'
elif os.path.isdir('/Volumes/data/pdb'):
    PDB_repository_base = '/Volumes/data/pdb/'


scale_val = 2

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
arg_parser.add_argument('-skip', dest='skip_count',
                        help='count of pdb ID list entries to skip')
arg_parser.add_argument('-limit', dest='limit_count',
                        help='stop after this many pdb ID list entries')
arg_parser.add_argument('-wp', help='write pdb file with .PyPDB extension',
                        action="store_true")
arg_parser.add_argument('-wi', help='write pic file with .PyPIC extension',
                        action="store_true")
arg_parser.add_argument('-ws', help='write OpenSCAD file with .scad extension',
                        action="store_true")
arg_parser.add_argument('-scale', dest='scale',
                        help='OpenSCAD output: units (usually mm) per '
                        'angstrom, default ' + str(scale_val))
arg_parser.add_argument('-maxp', dest='maxp',
                        help='max N-C peptide bond length for chain breaks,'
                        'default ' + str(PIC_Chain.MaxPeptideBond))
arg_parser.add_argument('-backbone',
                        help='OpenSCAD output: skip sidechains',
                        action="store_true")
arg_parser.add_argument('-t', help='test conversion pdb/pic to pic/pdb',
                        action="store_true")
arg_parser.add_argument('-tv', help='verbose test conversion pdb<>pic',
                        action="store_true")
arg_parser.add_argument('-nh', help='ignore hydrogens on PDB read',
                        action="store_true")
arg_parser.add_argument('-amide',
                        help='only amide proton, skip other Hs on PDB read',
                        action="store_true")
arg_parser.add_argument('-gcb',
                        help='generate GLY C-beta atoms',
                        action="store_true")
arg_parser.add_argument('-rama',
                        help='print psi, phi, omega values',
                        action="store_true")

args = arg_parser.parse_args()

# print(args)

if args.nh:
    PIC_Residue.accept_atoms = PIC_Residue.accept_backbone
if args.amide:
    PIC_Residue.accept_atoms = PIC_Residue.accept_backbone + ('H',)
if args.gcb:
    PIC_Residue.gly_Cbeta = True
if args.maxp:
    PIC_Chain.MaxPeptideBond = float(args.maxp)
if args.scale:
    scale_val = args.scale
if args.skip_count:
    args.skip_count = int(args.skip_count)
if args.limit_count:
    args.limit_count = int(args.limit_count)

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

PDB_parser = PDBParser(PERMISSIVE=True, QUIET=True)

fileNo = 1

for target in toProcess:
    if args.skip_count and fileNo <= args.skip_count:
        fileNo += 1
        continue
    if args.limit_count is not None:
        if args.limit_count <= 0:
            sys.exit(0)
        args.limit_count -= 1
    pdb_input = False
    pic_input = False
    pdb_structure = None
    # pdb_chain = None
    prot_id = ''
    outfile = os.path.basename(target)

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
        try:
            pdb_structure = PDB_parser.get_structure(
                prot_id,
                gzip.open(filename, mode='rt')
                if filename.endswith('.gz') else filename)
        except FileNotFoundError:
            print(filename, "not found")
            continue

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
        if not args.t:
            print(fileNo, 'parsed pdb input', prot_id, filename)
        # print('header:', pdb_structure.header.get('head', 'NONE'))
        # print('idcode:', pdb_structure.header.get('idcode', 'NONE'))
        # print('deposition date:', pdb_structure.header.get(
        #    'deposition_date', 'NONE'))
    #    for res in pdb_chain.get_residues():   # pdb_structure.get_residues():
    #        print(res.get_full_id(), res.resname,
    #              'disordered' if res.disordered else '')
    else:
        print(fileNo, 'parsed pic input ', filename)
        report_PIC(pdb_structure, verbose=True)

    # print(pdb_structure.header['idcode'], pdb_chain.id, ':',
    #      pdb_structure.header['head'])

    if args.wp:
        if pic_input:
            internal_to_atom_coordinates(pdb_structure)
        write_PDB(pdb_structure, outfile + '.PyPDB')
        print('wrote pdb output for', outfile)

    if args.wi:
        if pdb_input:
            # add_PIC(pdb_structure)
            atom_to_internal_coordinates(pdb_structure)
        write_PIC(pdb_structure, outfile + '.PyPIC')
        print('wrote pic output for', outfile)

    if args.t or args.tv:
        sp = StringIO()
        if pdb_input:
            with warnings.catch_warnings(record=True) as w:
                # warnings.simplefilter("error")
                try:
                    r = structure_rebuild_test(pdb_structure, args.tv)
                    warns = (len(w) > 0)
                    if args.tv and warns:
                        for wrn in w:
                            print(wrn.message)
                    print(prot_id, fileNo, r['report'],
                          ('WARNINGS' if warns else ''))
                except Exception as e:
                    print(prot_id, fileNo, 'ERROR FAIL:', type(e), e)

        elif pic_input:
            internal_to_atom_coordinates(pdb_structure)
            write_PDB(pdb_structure, sp)
            sp.seek(0)
            pdb2 = PDB_parser.get_structure(prot_id, sp)
            atom_to_internal_coordinates(pdb2)
            sp2 = StringIO()
            write_PIC(pdb2, sp2)
            sp2.seek(0)
            inf = open(filename, 'r')
            lineCount = 0
            matchCount = 0
            diffCount = 0
            # for line1, line2 in itertools.zip_longest(inf, sp2):
            for line1, line2 in zip(inf, sp2):
                lineCount += 1
                if line1 == line2:
                    matchCount += 1
                else:
                    diffCount += 1
                    if args.tv:
                        print(line1, '!=', line2)
            print(lineCount, matchCount, diffCount)
    if args.rama:
        if pdb_input:
            atom_to_internal_coordinates(pdb_structure)
        for r in pdb_structure.get_residues():
            # print(r.pic.get_dihedral('N:CA:C:O'))
            print(r, r.pic.get_angle('psi'), r.pic.get_angle('phi'),
                  r.pic.get_angle('omg'), r.pic.get_angle('chi2'),
                  r.pic.get_length('0C:1N'))

    if args.ws:
        write_SCAD(pdb_structure, target + '.scad', scale_val, pdbid=prot_id,
                   backboneOnly=args.backbone)

    fileNo += 1

print('normal termination')
