#!/usr/local/bin/python3

#
# replicate buildprot with biopython
#


import argparse
import re
import os
import sys

import gzip

from Bio import PDB

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import *

from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Entity import Entity

from Bio.Data import IUPACData

from Bio.PDB.StructureBuilder import StructureBuilder

import math

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


# ------------ utility functions

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
    for k in [0,1,2]:
        kl.append(di_hedron[k])
    a4 = di_hedron[3]
    if a4:
        kl.append(a4)
    return ':'.join(kl)


def check_res(res):
    rid = res.id
    if ' ' != rid[0]:
        return False
    if ' ' != rid[2]:
        return False
    return True


def gen_Mrot(angle_rads,axis):
    cosang = math.cos(angle_rads)
    sinang = math.sin(angle_rads)

    if 'z' == axis:
        return numpy.arr([[cosang,-sinang,0,0],
                          [sinang,cosang,0,0],
                          [0,0,1,0],
                          [0,0,0,1]],dtype=float64)
    elif 'y' == axis:
        return numpy.arr([[cosang,0,sinang,0],
                          [0, 1, 0, 0],
                          [-sinang,0,cosang,0],
                          [0,0,0,1]],dtype=float64)
    else:
        return numpy.arr([[1,0,0,0],
                          [0,cosang,-sinang,0],
                          [0,sinang,cosang,0],
                          [0,0,0,1]],dtype=float64)



# ------------ could add to protein/model classes


def internal_to_atom_coordinates(struct):
    for chn in struct.get_chains():
        chn.set_first_last()
    for chn in struct.get_chains():
        chn.link_residues()
    for chn in struct.get_chains():
        chn.render_dihedra()
    for chn in struct.get_chains():
        chn.assemble_residues()


# ------------ add to residue class


def residue_pic_init(self):
    self.hedra = {}      # array of hedron objects indexed by hedron keys ([resPos res atom] x3)
    self.dihedra = {}    # array of dihedron objects indexed by dihedron keys ([resPos res atom] x4)
    self.id3_dh_index = {}  # map of dihedron key (first 3 atom ids) to dihedron for all dihedra in Residue


def residue_pic_load(self,di_hedron):
    dhk = gen_key([di_hedron['a1'],di_hedron['a2'],di_hedron['a3'],di_hedron['a4']])
    try:
        if di_hedron['a4'] is not None:   # parse regex defaults this to None instead of KeyError
            self.dihedra[dhk] = Dihedron(di_hedron)
        else:
            self.hedra[dhk] = Hedron(di_hedron)
    except AttributeError:    # start residue created from pic file by PDB_parser not pic initialised
        self.pic_init()
        self.pic_load(di_hedron)


def residue_link_dehedra(self):
    '''
    k3i = {}
    for k, dihedron in pairs(self['dihedra']) do
        dihedron['res'] = self # each dihedron can find its residue
        k3 = dihedron['key3']
        if not k3i[k3] then
             k3i[k3] = {}
        k3i[k3][  # k3i[k3] +1 ] = dihedron       # map to find each dihedron from atom tokens 1,2,3
    self['key3index'] = k3i
     '''

    id3i = {}
    for dh in self.dihedra.values():
        dh.res = self               # each dihedron can find its residue
        id3 = dh.id3
        if id3 not in id3i:
            id3i[id3] = []
        id3i[id3].append(dh)        # map to find each dihedron from atom tokens 1,2,3
    self.id3_dh_index = id3i



def residue_render_hedra(self):
    for h in self.hedra.values():
        if h.updated:
            h.init_pos()


def residue_render_dihedra(self):
    for d in self.hedra.values():
        if d.updated:
            d.init_pos()


def residue_assemble(self,startpos):
    '''


    '''
    pass


Residue.pic_init = residue_pic_init
Residue.pic_load = residue_pic_load
Residue.link_dihedra = residue_link_dehedra
Residue.render_hedra = residue_render_hedra
Residue.render_dihedra = residue_render_dihedra
Residue.assemble = residue_assemble

# ------------ add to chain class


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


def chain_set_first_last(self):
    # SEQFAIL
    self.firstpos = 999999
    self.lastpos  = -999999

    for res in self.get_residues():
        if check_res(res):
            respos = res.id[1]
            if respos < self.firstpos:
                self.firstpos = respos
            if respos > self.lastpos:
                self.lastpos = respos



def chain_set_prev_post(self, res):
    #print(res.id, res.disordered)
    #if respos > self[firstPos] then res[prev] = self[residues][respos-1]
    #if respos < len(self[residues]) then res[next] = self[residues][respos+1]

    # SEQFAIL
    respos = res.id[1]
    if respos > self.firstpos:
        res.prev = self[respos-1]
    if respos < self.lastpos:
        res.post = self[respos+1]



def chain_link_residues(self):
    # SEQFAIL
    for res in self.get_residues():
        self.set_prev_next(res)
        res.link_dihedra()


def chain_render_dihedra(self):
    # SEQFAIL
    for res in self.get_residues():
        res.render_hedra()
    for res in self.get_residues():
        res.render_dihedra()


def chain_assemble_residues(self, start=False, fin=False):
    '''
    for res in orderec_residues
        if (not start or start<= respos) and (not fin or fin >= respos):
            if res[prev] and res[prev] is ordered:    # if no chain break
                startPos={}
                rp = res[prev]
                akl = res:NCaCKeySplit()  # ordered list - tuple?
                for ak in akl:
                    startPos[ak] = rp[atomCoords][ak]
            else:  # have chain break so use stopred start pos
                startPos = self[initNCaC][respos]
            res[atomCoords] = res:assemble(startPos)
    '''
    pass

Chain.pic_init = chain_init_pic
Chain.pic_load = chain_pic_load
Chain.set_first_last = chain_set_first_last
Chain.set_prev_next = chain_set_prev_post
Chain.link_residues = chain_link_residues
Chain.render_dihedra = chain_render_dihedra
Chain.assemble_residues = chain_assemble_residues


class HedronMatchError(Exception):
    pass


class HedronIncompleteError(Exception):
    pass


class Dihedron:
    def __init__(self, dh_dict):  # kwargs because parsed into group.dict

        self.a14 = [dh_dict['a1'],dh_dict['a2'],dh_dict['a3'],dh_dict['a4']]

        self.id = gen_key(self.a14)
        self.dihedral1 = float(dh_dict['dihedral1'])  # angle formed between two planes

        # no residue or position, just atoms # maybe also grab residue?
        self.dclass = ''
        atomre = re.compile(r'^\d+[a-zA-Z](\w+)$')
        for ak in self.a14:
            m = atomre.match(ak)
            self.dclass += m.group(1)

        self.id3 = self.a14[0] + ':' + self.a14[1] + ':' + self.a14[2]
        self.id32 = self.a14[1] + ':' + self.a14[2] + ':' + self.a14[3]

        # 4 matrices specifying hedron space coordinates of constituent atoms, initially atom3 on +Z axis
        self.InitialCoords = []

        # hedra making up this dihedron; set by self:setHedra()
        self.hedron1 = None
        self.hedron2 = None

        self.h1key = None
        self.h2key = None

        self.initial_coords = None
        self.a4_pre_rotation = None

        self.res = None        # Residue object which includes this dihedron; set by Residue:linkDihedra()
        self.reverse = False   # order of atoms in dihedron is reversed from order of atoms in hedra

        # flag indicting that atom coordinates are up to date (do not need to be recalculated from dihedral1)
        self.updated = True

        #print(self, self.dclass)

    def __str__(self):
        return '4-' + self.id + ' ' + self.dclass + ' ' + str(self.dihedral1)

    # find specified hedron on this residue or its adjacent neighbors
    @staticmethod
    def get_hedron(res,id3):
        hedron = res.hedra.get(id3,None)
        if not hedron and res.hasAttr(prev) and not res.prev.disordered and (' ' == res.prev.id[0]):  # not hetero
            hedron = res.prev.hedra.get(id3,None)
        if not hedron and res.hasAttr(post) and not res.post.disordered and (' ' == res.post.id[0]):
            hedron = res.post.hedra.get(id3,None)
        return hedron

    def set_hedra(self):
        rev = False
        res = self.res
        h1key = gen_key(self.a14[0:2])
        hedron1 = Dihedron.get_hedron(res,h1key)
        if not hedron1:
            rev = True
            h1key = gen_key(self.a14[2:0:-1])
            hedron1 = Dihedron.get_hedron(res,h1key)
            h2key = gen_key(self.a14[3:1:-1])
        else:
            h1key = gen_key(self.a14[1:3])

        if not hedron1:
            raise HedronMatchError(res, "can't find 1st hedron", h1key, self)

        hedron2 = Dihedron.get_hedron(res,h2key)

        if not hedron2:
            raise HedronMatchError(res, "can't find 2nd hedron", h2key, self)

        self.hedron1 = hedron1
        self.h1key = h1key
        self.hedron2 = hedron2
        self.h2key = h2key

        self.reverse = rev

        return rev

    def init_pos(self):

        rev = self.setHedra()
        hedron1 = self.hedron1
        hedron2 = self.hedron2

        '''
# not sure this happens...

   local complete = true
   for i=1,3 do
      if not hedron1['atoms'][i] then complete=false end
   end
   for i=1,3 do
      if not hedron2['atoms'][i] then complete=false end
   end
   if not complete then
      if utils.warn then
         io.stderr:write('dihedron: hedra missing atoms ' .. self:tostring() .. (reverse and ' reverse ' or ' forward ') .. hedron1:tostring() .. ' ' .. hedron2:tostring() .. '\n')
      end
      --io.write('dihedron: hedra missing atoms ' .. self:tostring() .. (reverse and ' reverse ' or ' forward ') .. hedron1:tostring() .. ' ' .. hedron2:tostring() .. '\n')
      return
   end
        '''

        initial = []

        if not rev:
            initial.append(hedron1.atoms[0].copy())
            initial.append(hedron1.atoms[1].copy())
            initial.append(hedron1.atoms[2].copy())

            a4_pre_rotation = hedron2.atomsR[2].copy()
            a4shift = hedron2.len1
        else:
            initial.append(hedron1.atomsR[2].copy())
            initial.append(hedron1.atomsR[1].copy())
            initial.append(hedron1.atomsR[0].copy())

            a4_pre_rotation = hedron2.atoms[0].copy()
            a4shift = hedron2.len3

        a4_pre_rotation[2][0] *= -1                     # a4 to +Z
        a4_pre_rotation[2][0] += a4shift                # hedron2 shift up so a2 at 0,0,0
        mrz = gen_Mrot(math.radians(self.dihedral1),'z')
        initial[3] = mrz * a4_pre_rotation

        self.initial_coords = initial
        self.a4_pre_rotation = a4_pre_rotation

        self.updated = False




class Hedron:
    def __init__(self, h_dict):  # kwargs because parsed into group.dict
        self.id = h_dict['a1'] + ':' + h_dict['a2'] + ':' + h_dict['a3']
        self.len1 = float(h_dict['len1'])     # distance between 1st and 2nd atom
        self.angle2 = float(h_dict['angle2'])  # angle formed between 3 atoms
        self.len3 = float(h_dict['len3'])     # distance between 2nd and 3rd atoms

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

    def init_pos(self):

        # build hedron with a2 on +Z axis, a1 at origin, a0 in -Z at angle n XZ plane
        for i in range(0,3):
            # note this initializes a1 to 0,0,0 origin
            self.atoms.append(numpy.array([[0],[0],[0],[1]], dtype=numpy.float64))  # 4x1 array

        # supplement: angles which add to 180 are supplementary
        sar = math.radians(180.0 - self.angle2)

        # a2 is len3 up from a2 on Z axis, X=Y=0
        self.atoms[2][2][0] = self.len3
        # a0 X is sin( sar ) * len1
        self.atoms[0][0][0] = math.sin(sar) * self.len1
        # a0 Z is -(cos( sar ) * len1) (assume angle2 always obtuse, so a0 is in -Z
        self.atoms[0][2][0] = - (math.cos(sar) * self.len3)

        # same again but 'reversed' : a0 on Z axis, a1 at origin, a2 in -Z
        for i in range(0,3):
            # atom[1] to 0,0,0 origin
            self.atomsR.append(numpy.array([[0],[0],[0],[1]], dtype=numpy.float64))

        # a0r is len1 up from a1 on Z axis, X=Y=0
        self.atomsR[0][2][0] = self.len1
        # a2r X is sin( sar ) * len3
        self.atomsR[2][0][0] = math.sin(sar)* self.len3
        # a2r Z is -(cos( sar ) * len3)
        self.atomsR[2][2][0] = - (math.cos(sar) * self.len3)

        self.updated = False



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
        # foo = {}
        # foo['a1'] = '1KN'
        # pdb_chain.pic_load(foo)
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


        internal_to_atom_coordinates(pdb_structure)

    print(pdb_structure.header['idcode'], pdb_chain.id, ':', pdb_structure.header['head'])



    if args.wp:
        io = PDBIO()
        io.set_structure(structure)
        io.save(f+'.PyPDB')
