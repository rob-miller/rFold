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

from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Entity import Entity

from Bio.Data import IUPACData

from Bio.PDB.StructureBuilder import StructureBuilder

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to build proteins from internal coordinates.")
''' 
class PIC_protein:
    def __init__(self,id):
        self.id = id
        self.chains = [{}]
        self.model = 0
        self.chains_maxpos = [{}]

    def load(self,filename):
        di_hedron_re = re.compile(r'^(?P<pdbid>\w+)\s+(?P<chn>\w?)\s*'    # pdbid and chain id
        '(?P<a1>-?[\w_]+):(?P<a2>-?[\w_]+):(?P<a3>-?[\w_]+)'              # 3 atom specifiers for hedron
        '(:(?P<a4>-?[\w_]+))?'                                            # 4th atom speicfier for dihedron
        '\s+'
        '|((?P<len1>\S+)\s+(?P<angle2>\S+)\s+(?P<len3>\S+)\s*$)'          # len-angle-len for hedron
        '((?P<dihedral1>\S+)\s*$)')                                       # dihedral angle for dihedron

        fil = open(filename, 'r')
        for aline in fil:
            m = di_hedron_re.match(aline)
            if m:
                mdl = self.model
                chn = m.group('chn')
                pos = int(re.search(r'\d+',m.group('a1')).group(0))
                #print(mdl,chn,pos)
                try:
                    if self.chains_maxpos[mdl][chn] > pos:
                        # no KeyError means already have this model, chain defined longer than current position
                        # so start new model
                        # print(pos + ' less than maxpos ' + self.chains_maxpos[mdl][chn] +
                        # ' on chain so new model and new chain')
                        mdl += 1
                        self.model = mdl
                        self.chains[mdl][chn] = PicChain(self, self.model, chn)
                except KeyError:  # KeyError means no chain in current model for this chain id
                    self.chains[mdl][chn] = PicChain(self, self.model, chn)

                self.chains[mdl][chn].load(m.groupdict())
                self.chains_maxpos[mdl][chn] = pos

            #else:
                #print('no match on ' + aline)


class PicChain:
    def __init__(self,parent,parent_model,id):
        self.id = id
        self.parent = parent
        self.parent_model = parent_model
        # dict of N,Ca,C atom coordinates for first residue of each ordered chain segment indexed by 1st atom id
        self.initNCaC = {}
        # possibly sparse array (dict) of residues comprising sequence
        self.PIC_residues = {}

    def load(self,kwargs):
        print('bar')

class PIC_residue:
    def __init__(self,**kwargs):
        pass

'''

akre = re.compile(r'(?P<pos>\d+)(?P<res>[A-Za-z])(?P<atm>[A-Za-z]+)')
dhkre = re.compile(r'(?P<a1>-?[\w_]+):(?P<a2>-?[\w_]+):(?P<a3>-?[\w_]+)(:(?P<a4>-?[\w_]+))?')

def split_atom_key(atom_key):
    rslt = akre.match(atom_key)
    return rslt.groups()

def split_dh_key(dh_key):
    rslt = dhkre.match(dh_key)
    return rslt.groups()

def gen_key(di_hedron):
    kl = []
    for k in ['a1','a2','a3']:
        kl.append(di_hedron[k])
    a4 = di_hedron['a4']
    if a4:
        kl.append(a4)
    return ':'.join(kl)


def residue_pic_init(self):
    self.hedra = {}      # array of hedron objects indexed by hedron keys ([resPos res atom] x3)
    self.dihedra = {}    # array of dihedron objects indexed by dihedron keys ([resPos res atom] x4)
    self.h2dhIndex = {}  # map of key(first 3 atom ids) to dihedron for all dihedra in Residue


def residue_pic_load(self,di_hedron):
    dhk = gen_key(di_hedron)
    try:
        if di_hedron['a4'] is not None:   # parse regex defaults this to None instead of KeyError
            self.dihedra[dhk] = Dihedron(di_hedron)
        else:
            self.hedra[dhk] = Hedron(di_hedron)
    except AttributeError:    # start residue created from pic file by PDB_parser not pic initialised
        self.pic_init()
        self.pic_load(di_hedron)


Residue.pic_init = residue_pic_init
Residue.pic_load = residue_pic_load


def chain_init_pic(self):
    # print('chain pic init')
    pass


def chain_pic_load(self,di_hedra):
    sak = split_atom_key(di_hedra['a1'])
    res_id = self._translate_id(int(sak[0]))
    try:
        res = self.__getitem__(res_id)
    except KeyError:
        # print(res_id,IUPACData.protein_letters_1to3[sak[1]].upper())
        # do outside Chain - so Chain() does not know about Residue()
        res = Residue(res_id,IUPACData.protein_letters_1to3[sak[1]].upper(),' ')
        self.add(res)
        res.pic_init()
        #print('created',res, 'for',di_hedra['a1'])
    res.pic_load(di_hedra)


Chain.pic_init = chain_init_pic
Chain.pic_load = chain_pic_load





class Dihedron:
    def __init__(self, dh_dict):  # kwargs because parsed into group.dict
        self.id = dh_dict['a1'] + ':' + dh_dict['a2'] + ':' + dh_dict['a3'] + ':' + dh_dict['a4']
        self.dihedral1 = dh_dict['dihedral1']  # angle formed between 3 atoms

        # no residue or position, just atoms
        self.dclass = ''
        atomre = re.compile(r'^\d+[a-zA-Z](\w+)$')
        for k in ('a1','a2','a3','a4'):
            m = atomre.match(dh_dict[k])
            self.dclass += m.group(1)

        self.id3  = dh_dict['a1'] + ':' + dh_dict['a2'] + ':' + dh_dict['a3']
        self.id32 = dh_dict['a2'] + ':' + dh_dict['a3'] + ':' + dh_dict['a4']

        # 4 matrices specifying hedron space coordinates of constituent atoms, initially atom3 on +Z axis
        self.InitialCoords = []

        # hedra making up this dihedron; set by self:setHedra()
        self.hedron1 = None
        self.hedron2 = None

        self.res = None        # Residue object which includes this dihedron; set by Residue:linkDihedra()
        self.reverse = False   # order of atoms in dihedron is reversed from order of atoms in hedra

        # flag indicting that atom coordinates are up to date (do not need to be recalculated from dihedral1)
        self.updated = True

        #print(self)

    def __str__(self):
        return '4-' + self.id + ' ' + self.dclass + ' ' + str(self.dihedral1)


class Hedron:
    def __init__(self, h_dict):  # kwargs because parsed into group.dict
        self.id = h_dict['a1'] + ':' + h_dict['a2'] + ':' + h_dict['a3']
        self.len1 = h_dict['len1']     # distance between 1st and 2nd atom
        self.angle2 = h_dict['angle2']  # angle formed between 3 atoms
        self.len3 = h_dict['len3']     # distance between 2nd and 3rd atoms

        # no position or residue, just atom chars at end of key
        self.hclass = ''
        atomre = re.compile(r'^\d+[a-zA-Z](\w+)$')
        for k in ('a1','a2','a3'):
            m = atomre.match(h_dict[k])
            self.hclass += m.group(1)

        # 3 matrices specifying hedron space coordinates of constituent atoms, initially atom3 on +Z axis
        self.atoms = []
        # 3 matrices, hedron space coordinates, reversed order - initially atom1 on +Z axis
        self.atomsR = []

        self.updated = True

        #print(self)

    def __str__(self):
        return '3-' + self.id + ' ' + self.hclass






PDB_repository_base = None

if os.path.isdir('/media/data'):
    PDB_repository_base = '/media/data/pdb/'
elif os.path.isdir('/Volumes/data'):
    PDB_repository_base = '/Volumes/data/pdb/'
                  

arg_parser = argparse.ArgumentParser(description='Interconvert .pic (pprotein internal coordinates) and .pdb (protein '
                                                 'data bank) files.')
arg_parser.add_argument('file', nargs='*',
                        help='a .pdb or .pic path/filename to read (first model, first chain),'
                        'or a PDB idCode with optional chain ID to read from {0} as .ent.gz'.format(
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
pdbidre = re.compile(r'(^\d(\w\w)\w)(\w)?$')

if args.filelist:
    flist = open(args.filelist,'r')
    for aline in flist:
        fields = aline.split()
        pdbidMatch = pdbidre.match(fields[0])
        if pdbidMatch:
            # print(m.group(1) + ' ' + m.group(2))
            # toProcess.append(PDB_repository_base + m.group(2) + '/pdb' + m.group(1) + '.ent.gz' )
            toProcess.append(pdbidMatch.group(0))


if len(toProcess):
    print(len(toProcess), 'entries to process')
else:
    print("no files to process. use '-h' for help")
    sys.exit(0)

PDB_parser = PDBParser(PERMISSIVE=False, QUIET=False)

for target in toProcess:
    pdb_input=False
    pic_input=False
    pdb_structure = None
    pdb_chain = None

    pdbidMatch = pdbidre.match(target)
    if pdbidMatch is not None:
        assert PDB_repository_base, 'PDB repository base directory missing, please configure for this host'
        pdb_input=True
        filename = PDB_repository_base + pdbidMatch.group(2).lower() + '/pdb' + pdbidMatch.group(1).lower() + '.ent.gz'
    else:
        filename = target

    pdb_structure = PDB_parser.get_structure('prot',
                                             gzip.open(filename, mode='rt') if filename.endswith('.gz')
                                             else filename)

    # get specified chain if given, else just pick first for now
    # count atoms to detect pic file if don't know already

    rCount = 0
    aCount = 0

    if pdbidMatch is not None and pdbidMatch.group(3) is not None:
        if pdb_structure[0][pdbidMatch.group(3)] is not None:
            pdb_chain = pdb_structure[0][pdbidMatch.group(3)]
            rCount = len(pdb_chain)  # does not work if chain not well defined like pic file
        else:
            print('chain ' + pdbidMatch.group(3) + ' not found in ' + filename)
            continue
    else:
        pdat = {}
        rct = 0
        mc = 0
        for mdl in pdb_structure:
            cc = 0
            for chn in mdl:
                rc = 0
                for res in chn:
                    rc += 1
                rct += rc
                pdat[mc,cc] = rc
                cc += 1
            mc += 1
        if mc>0 or cc>0:
            print(mc , ' models ', cc, ' chains ', rct, ' residues total')
            for mi in range(mc):
                for ci in range(cc):
                    print('model ', mi, ' chain ', ci, ' residues ', pdat[mi,ci])
            print('selecting first chain only')

        pdb_chain = pdb_structure[0].child_list[0]  # first model, first chain
        for residue in pdb_chain:
            rCount += 1
            for atom in residue:
                aCount += 1
                # print(atom)

    if not pdb_input:
        if aCount == 3*rCount:
            pic_input = True
        else:
            pdb_input = True

    if pdb_input:
        print('parsed pdb input ', filename)
        foo = {}
        foo['a1'] = '1KN'
        pdb_chain.pic_load(foo)
    else:
        print('parsed pic input ', filename)
        #print(pdb_chain)
        #print(pdb_chain[1])
        #print(pdb_chain.has_id(1))
        pdb_chain.pic_init()
        di_hedron_re = re.compile(
            r'^(?P<pdbid>\w+)\s+(?P<chn>\w?)\s*'                       # pdbid and chain id
            r'(?P<a1>-?[\w_]+):(?P<a2>-?[\w_]+):(?P<a3>-?[\w_]+)'      # 3 atom specifiers for hedron
            r'(:(?P<a4>-?[\w_]+))?'                                    # 4th atom speicfier for dihedron
            r'\s+'
            r'(((?P<len1>\S+)\s+(?P<angle2>\S+)\s+(?P<len3>\S+)\s*$)|'  # len-angle-len for hedron
            r'((?P<dihedral1>\S+)\s*$))')                               # dihedral angle for dihedron
        picfile = open(filename,'r')
        for aline in picfile:
            m = di_hedron_re.match(aline)
            if m:
                # print(m.groupdict())
                pdb_chain.pic_load(m.groupdict())

    print(pdb_structure.header['idcode'], pdb_chain.id, ':', pdb_structure.header['head'])



    if args.wp:
        io = PDBIO()
        io.set_structure(structure)
        io.save(f+'.PyPDB')
