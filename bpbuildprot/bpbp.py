#!/usr/local/bin/python3

#
# replicate buildprot with biopython
#


import argparse
import re
import os
import sys

import gzip

#from Bio.PDB import *

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import *

#class ICChain():
    
#class ICResidue():
#    def __init__(self)

#class ICDihedra():

#class ICHedra():




PDB_repository_base = None

if os.path.isdir('/media/data'):
    PDB_repository_base = '/media/data/pdb/'
elif (os.path.isdir('/Volumes/data')):
    PDB_repository_base = '/Volumes/data/pdb/'
                  

arg_parser = argparse.ArgumentParser(description='Interconvert .pic (pprotein internal coordinates) and '
                                             '.pdb (protein data bank) files.')
arg_parser.add_argument('file', nargs='*',
                        help='a .pdb or .pic path/filename to read, or a PDB idCode with optional chain ID'
                        'to read from {0} as .ent.gz'.format(
                        (PDB_repository_base or '[PDB resource not defined - please configure before use]')))
arg_parser.add_argument('-f', dest='filelist',
                        help='a Dunbrack cullPDB pdb ID list to read from {0} as .ent.gz'.format(
                        (PDB_repository_base or '[PDB resource not defined - please configure before use]')))
arg_parser.add_argument('-wp', help='write pdb file with .pypdb extension', action="store_true")
arg_parser.add_argument('-wi', help='write pic file with .pypdb extension', action="store_true")
arg_parser.add_argument('-ct', help='test conversion supplied pdb/pic to ', action="store_true")

args = arg_parser.parse_args()


print(args)

toProcess=args.file
pdbidre = re.compile(r'(^\d(\w\w)\w)(\w?)$')

if args.filelist:
    flist = open(args.filelist,'r')
    for aline in flist:
        fields = aline.split()
        m = pdbidre.match(fields[0])
        if m:
            # print(m.group(1) + ' ' + m.group(2))
            # toProcess.append(PDB_repository_base + m.group(2) + '/pdb' + m.group(1) + '.ent.gz' )
            toProcess.append(m.group(0))


if len(toProcess):
    print(len(toProcess), 'entries to process')
else:
    print("no files to process. use '-h' for help")
    sys.exit(0)

PDB_parser = PDBParser(PERMISSIVE=False, QUIET=True)

for target in toProcess:
    pdb_input=False
    pic_input=False
    pdb_structure = None
    pdb_chain = None

    m = pdbidre.match(target)
    if m is not None:
        assert PDB_repository_base, 'PDB repository base directory missing, please configure for this host'
        pdb_input=True
        filename = PDB_repository_base + m.group(2).lower() + '/pdb' + m.group(1).lower() + '.ent.gz'
    else:
        filename = target

    if filename.endswith('.gz'):
        pdb_structure = PDB_parser.get_structure('prot', gzip.open(filename, mode='rt'))
    else:
        pdb_structure = PDB_parser.get_structure('prot', filename)

    # for k,v in structure.header.items():
    #     print(k,v)

    if m is not None and m.group(3) is not None:
        if pdb_structure[0][m.group(3)] is not None:
            pdb_input = True
            pdb_chain = pdb_structure[0][m.group(3)]
        else:
            print('chain ' + m.group(3) + ' not found in ' + filename)
    else:
        pdb_chain = pdb_structure[0].get_chains()    # get any chain

    if pdb_chain is None:
        for model in pdb_structure:
            for pdb_chain in model:
                break;
            break;

    if pdb_input:
        rCount = len(pdb_chain)     # does not work if chain not well defined like pic file
    else:
        rCount = 0
        residue = None
        for residue in pdb_chain:
            rCount += 1

        if hasattr(pdb_chain, 'id'):
            pdb_input = true
        else:
            pic_input = True
            aCount = 0
            for atom in residue:
                aCount += 1



    if pdb_input:
        print('parsed pdb input ', filename, ' ok. residues: ', rCount)
        print(pdb_structure.header['idcode'], pdb_chain.id, ':', pdb_structure.header['head'])
    else:
        print('parsed pic input ', filename, ' ok. residues: ', rCount, 'atoms: ', aCount)
        print(pdb_structure.header['idcode'],':', pdb_structure.header['head'])

    if args.wp:
        io = PDBIO()
        io.set_structure(structure)
        io.save(f+'.PyPDB')
