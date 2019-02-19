#!/usr/local/bin/python3

#
# replicate buildprot with biopython
#

import argparse
import re
import os
import sys

import gzip

# from Bio import PDB

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.File import as_handle
from Bio.PDB.Polypeptide import three_to_one

# from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
# from Bio.PDB.Atom import Atom
# from Bio.PDB.Entity import Entity

from Bio.Data import IUPACData

# from Bio.PDB.StructureBuilder import StructureBuilder

import math

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy to build proteins from internal coordinates.")

__updated__ = '2019-02-19 21:22:58'
print('ver: ' + __updated__)
print(sys.version)


class PIC_Utils:

    # ------------ utility functions

    akre = re.compile(r'(?P<pos>\d+)(?P<res>[A-Za-z])(?P<atm>[A-Za-z]+)')
    dhkre = re.compile(r'(?P<a1>-?[\w_]+):(?P<a2>-?[\w_]+):(?P<a3>-?[\w_]+)'
                       r'(:(?P<a4>-?[\w_]+))?')
    di_hedron_re = re.compile(
        # pdbid and chain id
        r'^(?P<pdbid>\w+)\s+(?P<chn>\w?)\s*'
        # 3 atom specifiers for hedron
        r'(?P<a1>-?[\w_]+):(?P<a2>-?[\w_]+):(?P<a3>-?[\w_]+)'
        # 4th atom speicfier for dihedron
        r'(:(?P<a4>-?[\w_]+))?'
        r'\s+'
        # len-angle-len for hedron
        r'(((?P<len1>\S+)\s+(?P<angle2>\S+)\s+(?P<len3>\S+)\s*$)|'
        # dihedral angle for dihedron
        r'((?P<dihedral1>\S+)\s*$))')

    @staticmethod
    def split_atom_key(atom_key):
        rslt = PIC_Utils.akre.match(atom_key)
        return rslt.groups()

    @staticmethod
    def split_dh_key(dh_key):
        rslt = PIC_Utils.dhkre.match(dh_key)
        return rslt.groups()

    @staticmethod
    def gen_key(di_hedron):
        kl = di_hedron[0:3]
        if 4 == len(di_hedron):
            kl.append(di_hedron[3])
        return ':'.join(kl)

    @staticmethod
    def check_res(res):
        rid = res.id
        if ' ' != rid[0]:
            return False
        if ' ' != rid[2]:
            return False
        return True

    @staticmethod
    def gen_Mrot(angle_rads, axis):
        cosang = math.cos(angle_rads)
        sinang = math.sin(angle_rads)

        if 'z' == axis:
            return numpy.arr([[cosang, -sinang, 0, 0],
                              [sinang, cosang, 0, 0],
                              [0, 0, 1, 0],
                              [0, 0, 0, 1]], dtype=numpy.float64)
        elif 'y' == axis:
            return numpy.arr([[cosang, 0, sinang, 0],
                              [0, 1, 0, 0],
                              [-sinang, 0, cosang, 0],
                              [0, 0, 0, 1]], dtype=numpy.float64)
        else:
            return numpy.arr([[1, 0, 0, 0],
                              [0, cosang, -sinang, 0],
                              [0, sinang, cosang, 0],
                              [0, 0, 0, 1]], dtype=numpy.float64)

    # ------------ could add to protein/model classes

    @staticmethod
    def internal_to_atom_coordinates(struct):
        for chn in struct.get_chains():
            chn.set_first_last()
        for chn in struct.get_chains():
            chn.link_residues()
        for chn in struct.get_chains():
            chn.render_dihedra()
        for chn in struct.get_chains():
            chn.assemble_residues()

#    @staticmethod
#    def atoms_to_internal_coordinates(chn):
#        chn.pic.dihedra_from_atoms()

    @staticmethod
    def load_pic_file(file):
        with as_handle(file, mode='rU') as handle:
            for aline in handle.readlines():
                m = PIC_Utils.di_hedron_re.match(aline)
                if m:
                    # print(m.groupdict())
                    # pdb_chain.pic_load(m.groupdict())
                    pass


class HedronMatchError(Exception):
    pass


class HedronIncompleteError(Exception):
    pass


class PIC_Dihedron:
    def __init__(self, dh_dict):  # kwargs because parsed into group.dict

        self.a14 = [dh_dict['a1'], dh_dict['a2'], dh_dict['a3'], dh_dict['a4']]
        self.id = PIC_Utils.gen_key(self.a14)
        if 'dihedral1' in dh_dict:
            self.dihedral1 = float(dh_dict['dihedral1'])
        else:
            self.dihedral1 = None

        # no residue or position, just atoms # maybe also grab residue?
        self.dclass = ''
        atomre = re.compile(r'^\d+[a-zA-Z](\w+)$')
        for ak in self.a14:
            m = atomre.match(ak)
            self.dclass += m.group(1)

        self.id3 = self.a14[0] + ':' + self.a14[1] + ':' + self.a14[2]
        self.id32 = self.a14[1] + ':' + self.a14[2] + ':' + self.a14[3]

        # 4 matrices specifying hedron space coordinates of constituent atoms,
        # initially atom3 on +Z axis
        self.InitialCoords = []

        # hedra making up this dihedron; set by self:setHedra()
        self.hedron1 = None
        self.hedron2 = None

        self.h1key = None
        self.h2key = None

        self.initial_coords = None
        self.a4_pre_rotation = None

        # Residue object which includes this dihedron;
        # set by Residue:linkDihedra()
        self.res = None
        # order of atoms in dihedron is reversed from order of atoms in hedra
        self.reverse = False

        # flag indicting that atom coordinates are up to date
        # (do not need to be recalculated from dihedral1)
        self.updated = True

        # print(self, self.dclass)

    def __str__(self):
        return '4-' + self.id + ' ' + self.dclass + ' ' + str(self.dihedral1)

    # find specified hedron on this residue or its adjacent neighbors
    @staticmethod
    def get_hedron(res, id3):
        hedron = res.hedra.get(id3, None)
        if (not hedron and res.rprev and not res.prev.disordered
                and (' ' == res.rprev.id[0])):  # not hetero
            hedron = res.rprev.hedra.get(id3, None)
        if (not hedron and res.rnext and not res.rnext.disordered
                and (' ' == res.rnext.id[0])):
            hedron = res.rnext.hedra.get(id3, None)
        return hedron

    def set_hedra(self):
        rev = False
        res = self.res
        h1key = PIC_Utils.gen_key(self.a14[0:2])
        hedron1 = PIC_Dihedron.get_hedron(res, h1key)
        if not hedron1:
            rev = True
            h1key = PIC_Utils.gen_key(self.a14[2:0:-1])
            hedron1 = PIC_Dihedron.get_hedron(res, h1key)
            h2key = PIC_Utils.gen_key(self.a14[3:1:-1])
        else:
            h1key = PIC_Utils.gen_key(self.a14[1:3])

        if not hedron1:
            raise HedronMatchError(res, "can't find 1st hedron", h1key, self)

        hedron2 = PIC_Dihedron.get_hedron(res, h2key)

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
#   this warns if the dihedron does not have 2 complete 3-atom hedra

   local complete = true
   for i=1,3 do
      if not hedron1['atoms'][i] then complete=false end
   end
   for i=1,3 do
      if not hedron2['atoms'][i] then complete=false end
   end
   if not complete then
      if utils.warn then
         io.stderr:write('dihedron: hedra missing atoms ' .. self:tostring()
         .. (reverse and ' reverse ' or ' forward ') .. hedron1:tostring()
         .. ' ' .. hedron2:tostring() .. '\n')
      end
      --io.write('dihedron: hedra missing atoms ' .. self:tostring()
      .. (reverse and ' reverse ' or ' forward ') .. hedron1:tostring()
      .. ' ' .. hedron2:tostring() .. '\n')
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

        # a4 to +Z
        a4_pre_rotation[2][0] *= -1
        # hedron2 shift up so a2 at 0,0,0
        a4_pre_rotation[2][0] += a4shift
        mrz = PIC_Utils.gen_Mrot(math.radians(self.dihedral1), 'z')
        initial[3] = mrz * a4_pre_rotation

        self.initial_coords = initial
        self.a4_pre_rotation = a4_pre_rotation

        self.updated = False


class Hedron:
    def __init__(self, h_dict):  # kwargs because parsed into group.dict
        self.id = h_dict['a1'] + ':' + h_dict['a2'] + ':' + h_dict['a3']

        if 'len1' in h_dict:
            # distance between 1st and 2nd atom
            self.len1 = float(h_dict['len1'])
            # angle formed between 3 atoms
            self.angle2 = float(h_dict['angle2'])
            # distance between 2nd and 3rd atoms
            self.len3 = float(h_dict['len3'])
        else:
            self.len1 = None
            self.angle2 = None
            self.angle3 = None

        # no position or residue, just atom chars at end of key
        self.hclass = ''
        atomre = re.compile(r'^\d+[a-zA-Z](\w+)$')
        for k in ('a1', 'a2', 'a3'):
            m = atomre.match(h_dict[k])
            self.hclass += m.group(1)

        # 3 matrices specifying hedron space coordinates of constituent atoms,
        # initially atom3 on +Z axis
        self.atoms = []
        # 3 matrices, hedron space coordinates, reversed order
        # initially atom1 on +Z axis
        self.atomsR = []

        self.updated = True

        # print(self)

    def __str__(self):
        return '3-' + self.id + ' ' + self.hclass

    def init_pos(self):

        # build hedron with a2 on +Z axis, a1 at origin,
        # a0 in -Z at angle n XZ plane
        for i in range(0, 3):
            # note this initializes a1 to 0,0,0 origin
            self.atoms.append(numpy.array([[0], [0], [0], [1]],
                                          dtype=numpy.float64))  # 4x1 array

        # supplement: angles which add to 180 are supplementary
        sar = math.radians(180.0 - self.angle2)

        # a2 is len3 up from a2 on Z axis, X=Y=0
        self.atoms[2][2][0] = self.len3
        # a0 X is sin( sar ) * len1
        self.atoms[0][0][0] = math.sin(sar) * self.len1
        # a0 Z is -(cos( sar ) * len1)
        # (assume angle2 always obtuse, so a0 is in -Z)
        self.atoms[0][2][0] = - (math.cos(sar) * self.len3)

        # same again but 'reversed' : a0 on Z axis, a1 at origin, a2 in -Z
        for i in range(0, 3):
            # atom[1] to 0,0,0 origin
            self.atomsR.append(numpy.array([[0], [0], [0], [1]],
                                           dtype=numpy.float64))

        # a0r is len1 up from a1 on Z axis, X=Y=0
        self.atomsR[0][2][0] = self.len1
        # a2r X is sin( sar ) * len3
        self.atomsR[2][0][0] = math.sin(sar) * self.len3
        # a2r Z is -(cos( sar ) * len3)
        self.atomsR[2][2][0] = - (math.cos(sar) * self.len3)

        self.updated = False

# ------------ add to residue class


class PIC_Residue:

    def __init__(self, parent, ndx):
        self.residue = parent
        self.ndx = ndx
        # dict of hedron objects indexed by hedron keys ([resPos res atom] x3)
        self.hedra = {}
        # dict of dihedron objects indexed by dihedron keys
        # ([resPos res atom] x4)
        self.dihedra = {}
        # map of dihedron key (first 3 atom ids) to dihedron
        # for all dihedra in Residue
        # set by Residue.link_dihedra()
        self.id3_dh_index = {}
        # reference to adjacent residues in chain
        self.rprev = None
        self.rnext = None
        # one letter amino acid code
        self.lc = three_to_one(parent.resname).upper()
        # generated from dihedra include i+1 atoms, not residue.atoms
        # or initialised here from parent residue if loaded from coordinates
        self.atom_coords = {}
        rbase = str(ndx) + self.lc
        for atom in parent.get_atoms():
            ak = rbase + atom.name
            # TODO: are 4x1 arrays necessary?
            arr41 = numpy.append(atom.coord, [1])
            self.atom_coords[ak] = numpy.array(arr41, dtype=numpy.float64)

        # print(self.atom_coords)

    def pic_load(self, di_hedron):
        dhk = PIC_Utils.gen_key([di_hedron['a1'], di_hedron['a2'],
                                 di_hedron['a3'], di_hedron['a4']])
        # parse regex defaults this to None instead of KeyError
        if di_hedron['a4'] is not None:
            self.dihedra[dhk] = PIC_Dihedron(di_hedron)
        else:
            self.hedra[dhk] = Hedron(di_hedron)

    def link_dehedra(self):
        id3i = {}
        for dh in self.dihedra.values():
            dh.res = self        # each dihedron can find its residue.pic
            id3 = dh.id3
            if id3 not in id3i:
                id3i[id3] = []
            id3i[id3].append(dh)
        # map to find each dihedron from atom tokens 1-3
        self.id3_dh_index = id3i

    def render_hedra(self):
        for h in self.hedra.values():
            if h.updated:
                h.init_pos()

    def render_dihedra(self):
        for d in self.hedra.values():
            if d.updated:
                d.init_pos()

    def NCaCKeySplit(self):
        rbase = self.ndx + self.lc
        return [rbase + 'N', rbase + 'CA', rbase + 'C']

    def assemble(self, atomCoordsIn, genSCAD=False):
        '''
        --- join dihedrons from N-CA-C and N-CA-CB hedrons, computing protein
                space coordinates for backbone and sidechain atoms
        -- @param atomCoordsIn optional table of atom_token : 4x1 matrix of
                protein space coordinates
        -- @param genSCAD boolean if true, return table of transformation
                matrices for each hedron key (called from Chain:writeSCAD())
        -- @return atomCoords for residue in protein space relative to
                acomCoordsIn OR table of transformation matrices according to
                genSCAD parameter
        '''

        '''
    form queue, start with n-ca-c, o-c-ca, n-ca-cb
            [ o-c-ca not 2nd hedron for any dihedron and thus won't be picked
            up w/o adding here ]
    gen triple keys for current residue
    if no atomCoordsIn, use initial coords from generating dihedral for n-ca-c
            initial positions (dihedron coordinate space)

    while queue not empty
        get triple key
    for each dihedral starting with triple key (1st hedron)
            if have coordinates for all 4 atoms already
                add 2nd hedron key to back of queue
            else if have coordinates for 1st 3 atoms
                compute forward and reverse transform to take 1st 3 atoms
                        to/from dihedron initial coordinate space
                use reverse transform to get position of 4th atom in current
                        coordinates from dihedron initial coordinates
                add 2nd hedron key to back of queue
            else
                ordering failed, put triple key at back of queue and hope next
                    time we have 1st 3 atom positions (should not happen)

    loop terminates (queue drains) as triple keys which do not start any
            dihedra are removed without action
        '''

    pass

    def dihedra_from_atoms(self):
        skbase = str(self.ndx) + self.lc
        sN, sCA, sC, sO, sCB = skbase + 'N', skbase + \
            'CA', skbase + 'C', skbase + 'O', skbase + 'CB'
        pu = PIC_Utils

        if self.rnext:
            # atom_coords, hedra and dihedra for backbone dihedra
            # which reach into next residue
            rn = self.rnext
            nkbase = str(rn.ndx) + rn.lc
            nN, nCA, nC = nkbase + 'N', nkbase + 'CA', nkbase + 'C'

            for ak in (nN, nCA, nC):
                self.atom_coords[ak] = rn.atom_coords[ak]

            self.dihedra[pu.gen_key([sN, sCA, sC, nN])] = PIC_Dihedron(
                {'a1': sN, 'a2': sCA, 'a3': sC, 'a4': nN})   # psi

            self.dihedra[pu.gen_key([sCA, sC, nN, nCA])] = PIC_Dihedron(
                {'a1': sCA, 'a2': sC, 'a3': nN, 'a4': nCA})   # omega i+1

            self.dihedra[pu.gen_key([sC, nN, nCA, nC])] = PIC_Dihedron(
                {'a1': sC, 'a2': nN, 'a3': nCA, 'a4': nC})   # phi i+1

            rn.hedra[pu.gen_key([nN, nCA, nC])] = Hedron(
                {'a1': nN, 'a2': nCA, 'a3': nC})

        # backbone O and C-beta hedra and dihedra within this residue
        self.dihedra[pu.gen_key([sN, sCA, sC, sO])] = PIC_Dihedron(
            {'a1': sN, 'a2': sCA, 'a3': sC, 'a4': sO})
        self.dihedra[pu.gen_key([sO, sC, sCA, sCB])] = PIC_Dihedron(
            {'a1': sO, 'a2': sC, 'a3': sCA, 'a4': sCB})

        sNCaCKey = pu.gen_key([sN, sCA, sC])
        if sNCaCKey not in self.hedra:
            self.hedra[sNCaCKey] = Hedron({'a1': sN, 'a2': sCA, 'a3': sC})

        self.hedra[pu.gen_key([sCA, sC, sO])] = Hedron(
            {'a1': sCA, 'a2': sC, 'a3': sO})

        self.hedra[pu.gen_key([sCB, sCA, sC])] = Hedron(
            {'a1':  sCB, 'a2': sCA, 'a3': sC})

        # if res is not G or A then
        # only needed for sidechain CG residues
        # (not gly or ala or any missing rest of side chain)
        if ((skbase + 'CG') in self.atom_coords or
            (skbase + 'CG1') in self.atom_coords or
            (skbase + 'OG') in self.atom_coords or
            (skbase + 'OG1') in self.atom_coords or
                (skbase + 'SG') in self.atom_coords):
            self.hedra[pu.gen_key([sN, sCA, sCB])] = Hedron(
                {'a1': sN, 'a2': sCA, 'a3': sCB})

        # amide proton N H if present
        sH = skbase + 'H'
        if sH in self.atom_coords:
            self.hedra[pu.gen_key([sH, sN, sCA])] = Hedron(
                {'a1': sH, 'a2': sN, 'a3': sCA})
            self.hedra[pu.gen_key([sC, sCA, sN, sH])] = PIC_Dihedron(
                {'a1': sC, 'a2': sCA, 'a3': sN, 'a4': sH})

        

# Residue.pic_init = residue_pic_init
# Residue.pic_load = residue_pic_load
# Residue.link_dihedra = residue_link_dehedra
# Residue.render_hedra = residue_render_hedra
# Residue.render_dihedra = residue_render_dihedra
# Residue.NCaCKeySplit = residue_NCaCKeySplit
# Residue.assemble = residue_assemble

# ------------ add to chain class


class PIC_Chain:

    def __init__(self, parent):
        self.chain = parent
        self.ordered_aa_pic_list = []
        self.initNCaC = {}
        ndx = 0
        last_res = None
        last_ord_res = None
        for res in parent.get_residues():
            # select only ordered not hetero
            if res.id[0] == ' ' and not res.disordered:
                res.pic = PIC_Residue(res, ndx)
                self.ordered_aa_pic_list.append(res.pic)
                if last_res and last_ord_res == last_res:
                    # no missing or hetatm
                    last_ord_res.pic.rnext = res.pic
                    res.pic.rprev = last_ord_res.pic
                ndx += 1
                last_ord_res = res
            last_res = res

    def pic_load(self, di_hedra):
        # NOT READY
        sak = PIC_Utils.split_atom_key(di_hedra['a1'])
        res_id = self._translate_id(int(sak[0]))
        try:
            res = self.__getitem__(res_id)
        except KeyError:
            # print(res_id,IUPACData.protein_letters_1to3[sak[1]].upper())
            # do outside Chain - so Chain() does not know about Residue()
            res = Residue(res_id,
                          IUPACData.protein_letters_1to3[sak[1]].upper(), ' ')
            self.add(res)
            res.pic_init()
            # print('created',res, 'for',di_hedra['a1'])
        res.pic.pic_load(di_hedra)

    def link_residues(self):
        for rpic in self.ordered_aa_pic_list:
            rpic.link_dihedra()

    def render_dihedra(self):
        for rpic in self.ordered_aa_pic_list:
            rpic.render_hedra()
        for rpic in self.ordered_aa_pic_list:
            rpic.render_dihedra()

    def assemble_residues(self, start=False, fin=False):
        '''
        for res in ordered_residues
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

        for rpic in self.ordered_aa_pic_list:
            respos = rpic.parent.id[1]
            if (not start or start <= respos) and (not fin or fin >= respos):
                if rpic.rprev:
                    startPos = {}
                    rp = rpic.rprev
                    # nb akl for this res n-ca-c in rp dihedra
                    akl = rpic.NCaCKeySplit()
                    for ak in akl:
                        startPos[ak] = rp.atom_coords[ak]
                else:
                    # in theory the atom posns already added by load_structure
                    startPos = self.initNCaC[respos]
                # rtm unlike lua assemble() mods res coords
                rpic.assemble(startPos)

    def dihedra_from_atoms(self):
        for rpic in self.ordered_aa_pic_list:
            rpic.dihedra_from_atoms()


# Chain.pic_init = chain_init_pic
# Chain.pic_load = chain_pic_load
# Chain.set_first_last = chain_set_first_last
# Chain.set_prev_next = chain_set_prev_next
# Chain.link_residues = chain_link_residues
# Chain.render_dihedra = chain_render_dihedra
# Chain.assemble_residues = chain_assemble_residues


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
    pdb_chain = None
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
        pdb_structure = PIC_Utils.load_pic_file(
            gzip.open(filename, mode='rt')
            if filename.endswith('.gz') else filename)
        if pdb_structure:
            pic_input = True

    if not pdb_structure:
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
    else:
        if (len(pdb_structure) > 1) or (len(pdb_structure[0].child_list) > 1):
            print(pdb_structure.header['idcode'],
                  'selecting first model, first chain only')
        pdb_chain = pdb_structure[0].child_list[0]  # first model, first chain

    if pdb_input:
        print('parsed pdb input', prot_id, filename)
    #    for res in pdb_chain.get_residues():   # pdb_structure.get_residues():
    #        print(res.get_full_id(), res.resname,
    #              'disordered' if res.disordered else '')
        pdb_chain.pic = PIC_Chain(pdb_chain)
        pdb_chain.pic.dihedra_from_atoms()
    else:
        print('parsed pic input ', filename)
        PIC_Utils.internal_to_atom_coordinates(pdb_structure)

    print(pdb_structure.header['idcode'], pdb_chain.id, ':',
          pdb_structure.header['head'])

    if args.wp:
        io = PDBIO()
        io.set_structure(pdb_structure)
        io.save(target + '.PyPDB')
