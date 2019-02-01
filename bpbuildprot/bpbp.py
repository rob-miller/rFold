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
                  

parser = argparse.ArgumentParser(description='Parse PDB files on command line or from a file.')
parser.add_argument('file', nargs='*',
                    help='a pdb path/filename to read')
parser.add_argument('-f', dest='filelist',
                    help='a Dunbrack cullPDB pdb ID list to read from {0}'.format(
                        (PDB_repository_base or '[PDB resource not defined - please configure before use]')))
parser.add_argument('-w',help='write pdb file with .pypdb extension',action="store_true")

args = parser.parse_args()


print(args)

toProcess=args.file

if args.filelist:
    assert PDB_repository_base,'PDB repository base directory missing, please configure for this host'
    pdbidre = re.compile(r'(\d(\w\w)\w)(\w?)')
    flist = open(args.filelist,'r')
    for aline in flist:
        fields = aline.split()
        m = pdbidre.match(fields[0].lower())
        if m:
            # print(m.group(1) + ' ' + m.group(2))
            toProcess.append(PDB_repository_base + m.group(2) + '/pdb' + m.group(1) + '.ent.gz' )


parser = PDBParser(PERMISSIVE=False,QUIET=True,get_header=True)

for f in toProcess:
    if f.endswith('.gz'):
        print('trying compressed file ' + f + ' ...')
        structure = parser.get_structure('prot',gzip.open(f,mode='rt'))
    else:
        print('trying ' + f + ' ...')
        structure = parser.get_structure('prot',f)

    for k,v in structure.header.items():
        print(k,v)

    print(structure.header['head'])
    mCount=0
    cCount=0
    rCount=0
    for model in structure:
        mCount += 1
        for chain in model:
            cCount += 1
            for residue in chain:
                rCount += 1
                # print(residue)
                
    print('parsed ' , f , ' ok. models: ' , mCount , '  chains: ' , cCount , ' residues: ' , rCount)
    if args.w:
        io = PDBIO()
        io.set_structure(structure)
        io.save(f+'.PyPDB')
