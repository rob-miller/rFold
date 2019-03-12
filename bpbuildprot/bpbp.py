#!/usr/local/bin/python3

#
# replicate buildprot with biopython
#

import argparse
import re
import os
import sys

from collections import deque

# print(sys.path)

import gzip

# from Bio import PDB

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.File import as_handle
from Bio.PDB.Polypeptide import three_to_one

# from Bio.PDB.Chain import Chain
# from Bio.PDB.Residue import Residue
# from Bio.PDB.Atom import Atom
# from Bio.PDB.Entity import Entity

# from Bio.Data import IUPACData

from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.parse_pdb_header import _parse_pdb_header_list
from Bio.PDB.Atom import Atom, DisorderedAtom

# from Bio.PDB.Residue import Residue

import math
from Bio.PDB.vectors import rotaxis2m, calc_dihedral

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy to build proteins from internal coordinates.")

from PIC_Data import pic_data_sidechains

from Bio.PDB.PDBExceptions import PDBException

__updated__ = '2019-03-12 17:24:36'
print('ver: ' + __updated__)
print(sys.version)


# class PIC_Utils:

# ------------ utility functions

# akre = re.compile(
# r'(?P<pos>\d+)(?P<res>[A-Za-z])(?P<atm>[A-Za-z]+)(?P<occ>')
# dhkre = re.compile(r'(?P<a1>-?[\w-]+):(?P<a2>-?[\w-]+):(?P<a3>-?[\w-]+)'
#                    r'(:(?P<a4>-?[\w-]+))?')

#    @staticmethod
#    def split_atom_key(atom_key):
#        rslt = PIC_Utils.akre.match(atom_key)
#        return rslt.groups()

#    @staticmethod
#    def split_dh_key(dh_key):
#        rslt = PIC_Utils.dhkre.match(dh_key)
#        return rslt.groups()

# def check_res(res):
#    rid = res.id
#    if ' ' != rid[0]:
#        return False
#    if ' ' != rid[2]:
#        return False
#    return True


def gen_Mrot(angle_rads, axis):
    cosang = numpy.cos(angle_rads)
    sinang = numpy.sin(angle_rads)

    if 'z' == axis:
        return numpy.array([[cosang, -sinang, 0, 0],
                            [sinang, cosang, 0, 0],
                            [0, 0, 1, 0],
                            [0, 0, 0, 1]], dtype=numpy.float64)
    elif 'y' == axis:
        return numpy.array([[cosang, 0, sinang, 0],
                            [0, 1, 0, 0],
                            [-sinang, 0, cosang, 0],
                            [0, 0, 0, 1]], dtype=numpy.float64)
    else:
        return numpy.array([[1, 0, 0, 0],
                            [0, cosang, -sinang, 0],
                            [0, sinang, cosang, 0],
                            [0, 0, 0, 1]], dtype=numpy.float64)


def gen_Mtrans(xyz):
    return numpy.array([[1, 0, 0, xyz[0]],
                        [0, 1, 0, xyz[1]],
                        [0, 0, 1, xyz[2]],
                        [0, 0, 0, 1]
                        ], dtype=numpy.float64)


# ------------ could add to protein/model classes


def internal_to_atom_coordinates(struct):
    for chn in struct.get_chains():
        if hasattr(chn, 'pic'):
            chnPIC = chn.pic
            # chnPIC.set_first_last()  # TODO: clean up - not used here
            chnPIC.link_residues()
            chnPIC.render_dihedra()
            chnPIC.assemble_residues()
            chnPIC.coords_to_structure()

#    @staticmethod
#    def atoms_to_internal_coordinates(chn):
#        chn.pic.dihedra_from_atoms()


#
# read pic file
#
# file format:
#   (opt) PDB HEADER record
#       - idcode and deposition date recommended but optional
#       - deposition date in PDB format or as changed by Biopython
#   (opt) PDB TITLE record
#   repeat:
#   Biopython Residue Full ID - sets ID of returned structure
#   PDB ATOM records for chain start N, CA, C (optional)
#   PIC Hedra records for residue
#   PIC Dihedra records for residue
#
# returns either
#    Biopython Structure object, Residues with .pic attributes
#      but no coordinates except for chain start N, CA, C atoms
#  or
#    None on parse fail
#

def load_pic_file(file):
    pdb_hdr_re = re.compile(
        r'^HEADER\s{4}(?P<cf>.{1,40})'
        r'(?:\s+(?P<dd>\d\d\d\d-\d\d-\d\d|\d\d-\w\w\w-\d\d))?'
        r'(?:\s+(?P<id>[1-9][0-9A-Z]{3}))?\s*$',
    )
    pdb_ttl_re = re.compile(r'^TITLE\s{5}(?P<ttl>.+)\s*$')
    biop_id_re = re.compile(r'^\(\'(?P<pid>\w*)\',\s(?P<mdl>\d+),\s'
                            r'\'(?P<chn>\w)\',\s\(\'(?P<het>\s|[\w-]+)'
                            r'\',\s(?P<pos>\d+),\s\'(?P<icode>\s|\w)\'\)\)'
                            r'\s(?P<res>[A-Z]{3})'
                            r'\s\[(?P<segid>[a-zA-z\s]{4})\]'
                            r'\s(?P<ndx>\d+)?'
                            r'\s*$')
    pdb_atm_re = re.compile(r'^ATOM\s\s(?:\s*(?P<ser>\d+))\s(?P<atm>[\w\s]{4})'
                            r'(?P<alc>\w|\s)(?P<res>[A-Z]{3})\s(?P<chn>.)'
                            r'(?P<pos>[\s\-\d]{4})(?P<icode>[A-Za-z\s])\s\s\s'
                            r'(?P<x>[\s\-\d\.]{8})(?P<y>[\s\-\d\.]{8})'
                            r'(?P<z>[\s\-\d\.]{8})(?P<occ>[\s\d\.]{6})'
                            r'(?P<tfac>[\s\d\.]{6})\s{6}'
                            r'(?P<segid>[a-zA-z\s]{4})(?P<elm>.{2})'
                            r'(?P<chg>.{2})?\s*$')

    struct_builder = StructureBuilder()

    # init empty header dict
    # - could use to parse HEADERR and TITLE lines except
    # -- deposition_date format changed from original PDB header
    header_dict = _parse_pdb_header_list([])

    curr_SMCS = [None, None, None, None]  # struct model chain seg
    SMCS_init = [
        struct_builder.init_structure,
        struct_builder.init_model,
        struct_builder.init_chain,
        struct_builder.init_seg
    ]

    sb_res = None

    with as_handle(file, mode='rU') as handle:
        for aline in handle.readlines():
            if aline.startswith('HEADER '):
                m = pdb_hdr_re.match(aline)
                if m:
                    header_dict['head'] = m.group('cf')  # classification
                    header_dict['idcode'] = m.group('id')
                    header_dict['deposition_date'] = m.group('dd')

                    # print('HDR: classification: ', m.group('cf').strip(),
                    #      'dep date: ', m.group('dd')
                    #      if m.group('dd') else '',
                    #      'idcode: ', m.group('id')
                    #      if m.group('id') else ''
                    #      )
                else:
                    print('Reading pic file', file, 'HEADER fail: ', aline)
                pass
            elif aline.startswith('TITLE '):
                m = pdb_ttl_re.match(aline)
                if m:
                    header_dict['name'] = m.group('ttl').strip()
                    # print('TTL: ', m.group('ttl').strip())
                else:
                    print('Reading pic file', file, 'TITLE fail:, ', aline)
            elif aline.startswith('('):  # Biopython ID line for Residue
                m = biop_id_re.match(aline)
                if m:
                    # check Structure, Model, Chain, SegID
                    this_SMCS = [m.group(1), int(m.group(2)),
                                 m.group(3), m.group(8)]
                    if curr_SMCS != this_SMCS:
                        for i in range(4):
                            if curr_SMCS[i] != this_SMCS[i]:
                                SMCS_init[i](this_SMCS[i])
                                curr_SMCS[i] = this_SMCS[i]
                                if 0 == i:
                                    # 0 = init structure level so add header
                                    struct_builder.set_header(header_dict)

                    struct_builder.init_residue(
                        m.group('res'), m.group('het'),
                        int(m.group('pos')), m.group('icode'))

                    sb_res = struct_builder.residue
                    sb_res.pic = PIC_Residue(sb_res, m.group('ndx'))
                    # print('res id:', m.groupdict())
                    # print(report_PIC(struct_builder.get_structure()))
                else:
                    print('Reading pic file', file, 'res fail: ', aline)
            elif aline.startswith('ATOM '):
                m = pdb_atm_re.match(aline)
                if m:
                    if sb_res is None:
                        # ATOM without res spec already loaded, not a pic file
                        print('no sb_res - not pic file', aline)
                        return None
                    if (sb_res.resname != m.group('res')
                            or sb_res.id[1] != int(m.group('pos'))):
                        # TODO: better exception here?
                        raise Exception('pic ATOM read confusion:',
                                        sb_res.resname, str(sb_res.id),
                                        aline)

                    coord = numpy.array(
                        (float(m.group('x')), float(m.group('y')),
                         float(m.group('z'))), "f")
                    struct_builder.init_atom(
                        m.group('atm').strip(), coord, float(m.group('tfac')),
                        float(m.group('occ')), m.group('alc'), m.group('atm'),
                        int(m.group('ser')), m.group('elm').strip())

                    # print('atom: ', m.groupdict())
                else:
                    print('Reading pic file', file, 'ATOM fail: ', aline)
            else:
                m = DH_Base.edron_re.match(aline)
                if m:
                    sb_res.pic.load_PIC(m.groupdict())
                elif aline.strip():
                    print('Reading PIC file', file,
                          'parse fail on: .', aline, '.')
                    return None

    struct = struct_builder.get_structure()
    for chn in struct.get_chains():
        chn.pic = PIC_Chain(chn)

    # print(report_PIC(struct_builder.get_structure()))
    return struct


def write_PIC(entity, pdbid=None, chainid=None, s=''):
    try:
        if 'A' == entity.level:
            raise PDBException("No PIC output at Atom level")
        elif 'R' == entity.level:
            if hasattr(entity, 'pic'):
                if not chainid or not pdbid:
                    chain = entity.parent
                    if not chainid:
                        chainid = chain.id
                    if not pdbid:
                        struct = chain.parent.parent
                        pdbid = struct.header.get('idcode', '0PDB')
                s = entity.pic.write_PIC(pdbid, chainid, s)
            else:
                s += PIC_Residue.residue_string(entity)
        elif 'C' == entity.level:
            if not chainid:
                chainid = entity.id
            for res in entity:
                s = write_PIC(res, pdbid, chainid, s)
                pass
        elif 'M' == entity.level:
            for chn in entity:
                s = write_PIC(chn, pdbid, chainid, s)
        elif 'S' == entity.level:
            if not pdbid:
                pdbid = entity.header.get('idcode', None)
            hdr = entity.header.get('head', None)
            dd = entity.header.get('deposition_date', None)
            if hdr:
                s += ('HEADER    {:40}{:8}   {:4}\n'
                      ).format(hdr.upper(), (dd or ''), (pdbid or ''))
            nam = entity.header.get('name', None)
            if nam:
                s += 'TITLE     ' + nam.upper() + '\n'
            for mdl in entity:
                s = write_PIC(mdl, pdbid, chainid, s)
        else:
            raise PDBException("Cannot identify level: " + str(entity.level))
    except KeyError:
        raise Exception("write_PIC: argument is not a Biopython PDB Entity "
                        + str(entity))
    return s


def report_PIC(entity, reportDict=None):
    if reportDict is None:
        reportDict = {'idc': None, 'hdr': 0, 'mdl': 0, 'chn': 0,
                      'res': 0, 'res_e': 0, 'dih': 0, 'hed': 0}
    try:
        if 'A' == entity.level:
            raise PDBException("No PIC output at Atom level")
        elif 'R' == entity.level:
            if hasattr(entity, 'pic'):
                reportDict['res'] += 1
                dlen = len(entity.pic.dihedra)
                hlen = len(entity.pic.hedra)
                if 0 < dlen or 0 < hlen:
                    reportDict['res_e'] += 1
                    reportDict['dih'] += dlen
                    reportDict['hed'] += hlen

        elif 'C' == entity.level:
            reportDict['chn'] += 1
            for res in entity:
                reportDict = report_PIC(res, reportDict)

        elif 'M' == entity.level:
            reportDict['mdl'] += 1
            for chn in entity:
                reportDict = report_PIC(chn, reportDict)

        elif 'S' == entity.level:
            if reportDict['idc'] is None:
                reportDict['idc'] = entity.header.get('idcode', None)

            hdr = entity.header.get('head', None)
            if hdr:
                reportDict['hdr'] += 1
            nam = entity.header.get('name', None)
            if nam:
                reportDict['hdr'] += 1
            for mdl in entity:
                reportDict = report_PIC(mdl, reportDict)
        else:
            raise PDBException("Cannot identify level: " + str(entity.level))
    except KeyError:
        raise Exception("write_PIC: argument is not a Biopython PDB Entity "
                        + str(entity))
    return reportDict


def zdh(lst):
    return dict(zip(['a1', 'a2', 'a3', 'a4'], lst))


def genCBhamelryck(res):
    # https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
    # How do I put a virtual Cβ on a Gly residue?
    n = res['N'].get_vector()
    c = res['C'].get_vector()
    ca = res['CA'].get_vector()
    # center at origin
    n = n - ca
    c = c - ca
    # find rotation matrix that rotates n -120 degrees along the ca-c vector
    rot = rotaxis2m(-math.pi*120.0/180.0, c)
    # apply rotation to ca-n vector
    cb_at_origin = n.left_multiply(rot)
    # put on top of ca atom
    cb = cb_at_origin + ca
    return cb


def genCBjones(res):
    # code by DT Jones, based on
    # Prot. Engr vol. 2, issue 2, p. 121 (1 July, 1988)
    # Model building of disulfide bonds in proteins with known three-
    # dimensional structure
    # Hazes and Dijkstra, pp 119-125
    # https://doi.org/10.1093/protein/2.2.119
    # https://academic.oup.com/peds/article/2/2/119/1484803
    # TODO: update the constants!  Choose for Ala
    n = res['N'].coord
    c = res['C'].coord
    ca = res['CA'].coord

    nca = ca - n
    cca = ca - c
    xx = nca + cca
    yy = numpy.cross(nca, cca)
    kCACBDIST = 1.538
    kTETH_ANG = 0.9128

    sx = kCACBDIST * numpy.cos(kTETH_ANG) / numpy.sqrt(xx.dot(xx))
    sy = kCACBDIST * numpy.sin(kTETH_ANG) / numpy.sqrt(yy.dot(yy))
    cb = ca + (xx * sx + yy * sy)
    return cb


def set_accuracy_95(num):
    return float("{:9.5f}".format(num))


def set_accuracy_83(num):
    return float("{:8.3f}".format(num))


# TODO: rtm changed to remove [0] for assemble()
def get_spherical_coordinates(xyz):
    r = numpy.linalg.norm(xyz)
    if 0 == r:
        return numpy.array([0, 0, 0])
    sign = -1.0 if xyz[1][0] < 0.0 else 1.0
    theta = ((math.pi/2.0 * sign) if 0 == xyz[0][0]
             else numpy.arctan2(xyz[1][0], xyz[0][0]))
    phi = numpy.arccos(xyz[2][0]/r)
    return [r, theta, phi]


# acs[0] on XZ plane
# acs[1] origin
# acs[2] on +Z axis
# return transformation matrix
# if rev, return reverse transformation matrix (to return from coord_space)
def coord_space(acs, rev=False):
    dbg = False
    if dbg:
        for ac in acs:
            print(ac.transpose())

    a0 = acs[0]
    a1 = acs[1]
    a2 = acs[2]

    # tx acs[1] to origin
    tm = gen_Mtrans([-a1[0], -a1[1], -a1[2]])

    # directly translate a2 using a1
    p = a2 - a1
    sc = get_spherical_coordinates(p)

    if dbg:
        print('p', p.transpose())
        print('sc', sc)

    mrz = gen_Mrot(-sc[1], 'z')  # rotate translated a3 -theta about Z
    mry = gen_Mrot(-sc[2], 'y')  # rotate translated a3 -phi about Y

    # mt completes a2-a3 on Z-axis, still need to align a1 with XZ plane
    mt = mry @ mrz @ tm

    if dbg:
        print('mt * a2', (mt @ a2).transpose())

    p = mt @ a0

    # need theta of translated a0
    # sc2 = get_spherical_coordinates(p)
    sign = 1.0 if (p[1][0] < 0.0) else 1.0
    theta2 = ((math.pi/2.0 * sign) if 0 == p[0][0]
              else numpy.arctan2(p[1][0], p[0][0]))
    # rotate a0 -theta2 about Z to align with X
    mrz2 = gen_Mrot(-theta2, 'z')

    mt = mrz2 @ mt

    if not rev:
        return mt

    # generate the reverse transformation

    # rotate a0 theta about Z, reversing alignment with X
    mrz2 = gen_Mrot(theta2, 'z')
    # rotate a2 phi about Y
    mry = gen_Mrot(sc[2], 'y')
    # rotate a2 theta about Z
    mrz = gen_Mrot(sc[1], 'z')
    # translation matrix origin to a1
    tm = gen_Mtrans([a1[0], a1[1], a1[2]])

    # mr = mry @ mrz2
    # mr = mrz @ mr
    # mr = tm @ mr

    mr = tm @ mrz @ mry @ mrz2
    return mt, mr


class HedronMatchError(Exception):
    pass


class HedronIncompleteError(Exception):
    pass


class DihedronIncompleteError(Exception):
    pass


class Atom_Key(object):
    atom_re = re.compile(r'^(?P<respos>-?\d+)(?P<icode>[A-Za-z])?'
                         r'_(?P<resname>[a-zA-Z]+)_(?P<atm>\w+)'
                         r'(?:_(?P<occ>\d\.\d?\d?))?(?:_(?P<altloc>\w))?$')
    # PDB altLoc = Character = [\w ] (any non-ctrl ASCII incl space)
    # PDB iCode = AChar = [A-Za-z]
    fieldNames = ['respos', 'icode', 'resname', 'atm', 'occ', 'altloc']
    fields = {'respos': 0, 'icode': 1, 'resname': 2,
              'atm': 3, 'occ': 4, 'altloc': 5}

    def __init__(self, *args, **kwargs):
        # pass me:
        # (<PIC Residue>, 'CA', ...)
        # ([52, None, 'G', 'CA', ...])
        # (52, None, 'G', 'CA', ...)
        # {respos: 52, icode: None, resname: 'G', ...}
        # 52-G-CA, 52B-G-CA, 52-G-CA_0.33, etc

        akl = []
        self.id = None

        for arg in args:
            if isinstance(arg, PIC_Residue):
                if [] != akl:
                    raise Exception('Atom Key init Residue not first argument')
                akl += arg.rbase
            elif isinstance(arg, Atom):
                if 3 != len(akl):
                    raise Exception('Atom Key init Atom before Residue info')
                akl.append(arg.name)
                occ = arg.occupancy
                akl.append(occ if occ != 1.00 else None)
                altloc = arg.altloc
                akl.append(altloc if altloc != ' ' else None)
            elif isinstance(arg, list):
                akl += arg
            elif '_' in arg:
                # got atom key string, recurse with regex parse
                m = self.atom_re.match(arg)
                if [] != akl:
                    raise Exception(
                        'Atom Key init full key not first argument', arg)
                for fn in Atom_Key.fieldNames:
                    akl.append(m.group(fn))
            else:
                akl.append(arg)

        # did not initialize recursively above so continue here
        akl[0] = str(akl[0])
        for i in range(6):
            if len(akl) <= i:
                akl.append(kwargs.get(self.fieldNames[i], None))

        if akl[4] is not None:
            akl[4] = str(akl[4])

        # now we have len(akl) = 6
        #  start with Residue object if supplied
        #  then args in order
        #  then kwargs or None according to self.fieldNames

        (self.respos, self.icode, self.resname, self.atm, self.occ,
            self.altloc) = akl

        self.akl = tuple(akl)

        self.id = '_'.join(
            [''.join(filter(None, akl[:2])),
                akl[2],
                '_'.join(filter(None, akl[3:]))]
        )

        self._hash = hash(self.akl)
        # pass

    def __repr__(self):
        return self.id

    #def __str__(self):
    #    return self.id

    def __hash__(self):
        return self._hash

    backbone_sort_keys = {'N': 0, 'CA': 1, 'C': 2, 'O': 3}

    def altloc_match(self, other):
        if isinstance(other, type(self)):
            return self.akl[:5] == other.akl[:5]
        else:
            return NotImplemented

    def is_sidechain_gamma(self):
        if 'G' in self.akl[3]:
            return True
        return False

    def _cmp(self, other):
        akl_s = self.akl
        akl_o = other.akl
        atmNdx = self.fields['atm']
        for i in range(6):
            s, o = akl_s[i], akl_o[i]
            if s != o:
                if atmNdx != i:  # only sorting complication is atom level
                    return s, o
                # print('>', s, '<', '>', o, '<')
                if (s in self.backbone_sort_keys
                        and o in self.backbone_sort_keys):
                    return (self.backbone_sort_keys[s],
                            self.backbone_sort_keys[o])
                elif (s in self.backbone_sort_keys
                        and o not in self.backbone_sort_keys):
                    return 0, 1
                elif (s not in self.backbone_sort_keys
                        and o in self.backbone_sort_keys):
                    return 1, 0
                s0, o0 = s[0], o[0]
                if (s0 == 'H' and o0 != 'H'):
                    return 1, 0
                elif (s0 != 'H' and o0 == 'H'):
                    return 0, 1
                else:
                    return s, o
        return 1, 1

    def __eq__(self, other):
        """Test for equality."""
        if isinstance(other, type(self)):
            return self.akl == other.akl
        else:
            return NotImplemented

    def __ne__(self, other):
        """Test for inequality."""
        if isinstance(other, type(self)):
            return self.akl != other.akl
        else:
            return NotImplemented

    def __gt__(self, other):
        """Test greater than."""
        if isinstance(other, type(self)):
            rslt = self._cmp(other)
            return rslt[0] > rslt[1]
        else:
            return NotImplemented

    def __ge__(self, other):
        """Test greater or equal."""
        if isinstance(other, type(self)):
            rslt = self._cmp(other)
            return rslt[0] >= rslt[1]
        else:
            return NotImplemented

    def __lt__(self, other):
        """Test less than."""
        if isinstance(other, type(self)):
            rslt = self._cmp(other)
            return rslt[0] < rslt[1]
        else:
            return NotImplemented

    def __le__(self, other):
        """Test less or equal."""
        if isinstance(other, type(self)):
            rslt = self._cmp(other)
            return rslt[0] <= rslt[1]
        else:
            return NotImplemented


class DH_Base(object):
    edron_re = re.compile(
        # pdbid and chain id
        r'^(?P<pdbid>\w+)\s(?P<chn>[\w|\s])\s'
        # 3 atom specifiers for hedron
        r'(?P<a1>[\w\-\.]+):(?P<a2>[\w\-\.]+):(?P<a3>[\w\-\.]+)'
        # 4th atom speicfier for dihedron
        r'(:(?P<a4>[\w\-\.]+))?'
        r'\s+'
        # len-angle-len for hedron
        r'(((?P<len1>\S+)\s+(?P<angle2>\S+)\s+(?P<len3>\S+)\s*$)|'
        # dihedral angle for dihedron
        r'((?P<dihedral1>\S+)\s*$))')

    def __init__(self, *args, **kwargs):
        # pass me
        #  [ atom key, ... ]
        #  atom key, ...
        #  {'a1': str, 'a2': str, ... }

        aks = []
        for arg in args:
            if isinstance(arg, list):
                aks = arg
            elif isinstance(arg, tuple):
                aks = list(arg)
            else:
                if arg is not None:
                    aks.append(arg)
        if [] == aks:
            # self.aks = [Atom_Key(kwargs['a1']), Atom_Key(kwargs['a2']),
            #            Atom_Key(kwargs['a3'])]
            aks = [kwargs['a1'], kwargs['a2'], kwargs['a3']]

            try:
                if kwargs['a4'] is not None:
                    # self.aks.append(Atom_Key(kwargs['a4']))
                    aks.append(kwargs['a4'])
            except KeyError:
                pass

        # if args are atom key strings instead of Atom_Keys
        for i in range(len(aks)):
            if not isinstance(aks[i], Atom_Key):
                aks[i] = Atom_Key(aks[i])

        self.aks = tuple(aks)
        self.id = DH_Base.gen_key(aks)
        self._hash = hash(self.aks)
        
        # flag indicting that atom coordinates are up to date
        # (do not need to be recalculated from dihedral1)
        self.updated = True

        # no residue or position, just atoms
        self.dh_class = ''
        # same but residue specific
        self.rdh_class = ''

        atmNdx = Atom_Key.fields['atm']
        resNdx = Atom_Key.fields['resname']
        for ak in aks:
            akl = ak.akl
            self.dh_class += akl[atmNdx]
            self.rdh_class += akl[resNdx] + akl[atmNdx]

    def __repr__(self):
        return self.aks

    def __hash__(self):
        return self._hash
        
    def gen_key(lst):
        if isinstance(lst[0], Atom_Key):
            return ':'.join(ak.id for ak in lst)
        else:
            return ':'.join(lst)

# make the keys lists of atom objects and compare should be automatic?

    def _cmp(self, other):
        for ak_s, ak_o in zip(self.aks, other.aks):
            if ak_s != ak_o:
                return ak_s, ak_o
        return 1, 1

    def __eq__(self, other):
        """Test for equality."""
        if isinstance(other, type(self)):
            return self.id == other.id
        else:
            return NotImplemented

    def __ne__(self, other):
        """Test for inequality."""
        if isinstance(other, type(self)):
            return self.id != other.id
        else:
            return NotImplemented

    def __gt__(self, other):
        """Test greater than."""
        if isinstance(other, type(self)):
            rslt = self._cmp(other)
            return rslt[0] > rslt[1]
        else:
            return NotImplemented

    def __ge__(self, other):
        """Test greater or equal."""
        if isinstance(other, type(self)):
            rslt = self._cmp(other)
            return rslt[0] >= rslt[1]
        else:
            return NotImplemented

    def __lt__(self, other):
        """Test less than."""
        if isinstance(other, type(self)):
            rslt = self._cmp(other)
            return rslt[0] < rslt[1]
        else:
            return NotImplemented

    def __le__(self, other):
        """Test less or equal."""
        if isinstance(other, type(self)):
            rslt = self._cmp(other)
            return rslt[0] <= rslt[1]
        else:
            return NotImplemented


'''
    @staticmethod
    def gen_key(edron):
        kl = edron[0:3]
        if 4 == len(edron) and edron[3] is not None:
            kl.append(edron[3])
        return ':'.join(kl)
'''


class PIC_Dihedron(DH_Base):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if 'dihedral1' in kwargs:
            self.dihedral1 = float(kwargs['dihedral1'])
        else:
            self.dihedral1 = None

        # hedra making up this dihedron; set by self:set_hedra()
        self.hedron1 = None
        self.hedron2 = None

        self.h1key = None
        self.h2key = None

        self.id3 = tuple(self.aks[0:3]) #DH_Base.gen_key(self.aks[0:3])
        self.id32 = tuple(self.aks[1:4]) #DH_Base.gen_key(self.aks[1:4])

        # 4 matrices specifying hedron space coordinates of constituent atoms,
        # initially atom3 on +Z axis
        self.initial_coords = None
        self.a4_pre_rotation = None

        # Residue object which includes this dihedron;
        # set by Residue:linkDihedra()
        self.residue = None
        # order of atoms in dihedron is reversed from order of atoms in hedra
        self.reverse = False

        # print(self, self.dclass)

    def __str__(self):
        return ('4-' + str(self.id) + ' ' + self.rdh_class + ' ' +
                str(self.dihedral1) + ' (' + str(self.residue) + ')')

    # find specified hedron on this residue or its adjacent neighbors
    @staticmethod
    def get_hedron(res, id3):
        hedron = res.hedra.get(id3, None)
        if (not hedron and res.rprev
                and (' ' == res.rprev.residue.id[0])):  # not hetero
            hedron = res.rprev.hedra.get(id3, None)
        if (not hedron and res.rnext
                and (' ' == res.rnext.residue.id[0])):
            hedron = res.rnext.hedra.get(id3, None)
        return hedron

    def set_hedra(self):
        rev = False
        res = self.residue
        h1key = self.id3
        hedron1 = PIC_Dihedron.get_hedron(res, h1key)
        if not hedron1:
            rev = True
            h1key = tuple(self.aks[2::-1])
            hedron1 = PIC_Dihedron.get_hedron(res, h1key)
            h2key = tuple(self.aks[3:0:-1])
        else:
            h2key = self.id32

        if not hedron1:
            raise HedronMatchError(res, "can't find 1st hedron", h1key, self)

        hedron2 = PIC_Dihedron.get_hedron(res, h2key)

        if not hedron2:
            # print(res.hedra)
            # print(res.rnext.hedra)
            raise HedronMatchError(res, "can't find 2nd hedron", h2key, self)

        self.hedron1 = hedron1
        self.h1key = h1key
        self.hedron2 = hedron2
        self.h2key = h2key

        self.reverse = rev

        return rev

    def init_pos(self):

        rev = self.set_hedra()
        hedron1 = self.hedron1
        hedron2 = self.hedron2


# not sure this happens...
#   this warns if the dihedron does not have 2 complete 3-atom hedra

        acount = 0
        for a in hedron1.atoms:
            if a is not None:
                acount += 1
        for a in hedron2.atoms:
            if a is not None:
                acount += 1
        if 6 > acount:
            raise Exception('dihedron: hedra missing atoms', self)

        '''
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
        mrz = gen_Mrot(numpy.deg2rad(self.dihedral1), 'z')
        initial.append(mrz @ a4_pre_rotation)

        self.initial_coords = initial
        self.a4_pre_rotation = a4_pre_rotation

        self.updated = False

    def find_bp_atom(self, ak):
        bpa = self.residue.bp_atoms.get(ak, None)
        if bpa is not None:
            return bpa
        if self.residue.rnext:
            bpa = self.residue.rnext.bp_atoms.get(ak, None)
            if bpa is not None:
                return bpa
        if self.residue.rprev:
            bpa = self.residue.rprev.bp_atoms.get(ak, None)
        return bpa

    # get dist angle dist angle dist
    @staticmethod
    def get_dadad(acs):
        a0 = acs[0].squeeze()
        a1 = acs[1].squeeze()
        a2 = acs[2].squeeze()
        a3 = acs[3].squeeze()

        a0a1 = numpy.linalg.norm(a0-a1)
        a1a2 = numpy.linalg.norm(a1-a2)
        a2a3 = numpy.linalg.norm(a2-a3)

        a0a2 = numpy.linalg.norm(a0-a2)
        a1a3 = numpy.linalg.norm(a1-a3)

        sqr_a1a2 = a1a2 * a1a2

        a0a1a2 = numpy.rad2deg(
            numpy.arccos(
                ((a0a1*a0a1) + sqr_a1a2 - (a0a2*a0a2)) / (2*a0a1*a1a2)
            )
        )

        a1a2a3 = numpy.rad2deg(
            numpy.arccos(
                (sqr_a1a2 + (a2a3*a2a3) - (a1a3*a1a3)) / (2*a1a2*a2a3)
            )
        )

        return a0a1, a0a1a2, a1a2, a1a2a3, a2a3

    def dihedron_from_atoms(self):
        # call link_dihedra before this so can find res->atom_coords
        rev = self.set_hedra()
        hed1 = self.hedron1
        hed2 = self.hedron2

        # set_hedra will catch this
        # if (hed1 is None) or (hed2 is None):
        #    raise DihedronIncompleteError(self, 'missing hedra', hed1, hed2)

        atom_coords = self.residue.atom_coords
        aks = self.aks
        acs = []
        estr = ''
        for ak in aks:
            ac = atom_coords[ak]
            if ac is None:
                estr += ak + ' '
            else:
                acs.append(ac)
        if estr != '':
            raise DihedronIncompleteError(
                self, 'missing coordinates for', estr)

        # compare rtm coordspace dihedral
        # vs Vector() dihedral with self.bp_atoms

        # print(self.id)
        mt = coord_space(acs[:3])
        do4 = mt @ acs[3]

        dh1r = numpy.rad2deg(numpy.arctan2(do4[1][0], do4[0][0]))

        a0 = self.find_bp_atom(aks[0])
        a1 = self.find_bp_atom(aks[1])
        a2 = self.find_bp_atom(aks[2])
        a3 = self.find_bp_atom(aks[3])

        NODIHED = -999999
        if a0 and a1 and a2 and a3:
            dh1h = calc_dihedral(
                a0.get_vector(), a1.get_vector(),
                a2.get_vector(), a3.get_vector()
            )
            dh1h = numpy.rad2deg(dh1h)
        else:
            # no gly cb
            dh1h = NODIHED

        dh1h = set_accuracy_95(dh1h)
        dh1r = set_accuracy_95(dh1r)

        if dh1h != dh1r and NODIHED != dh1h:
            print('dihedral disagreement:', self, dh1h, dh1r)

        self.dihedral1 = dh1r

        a0a1, a0a1a2, a1a2, a1a2a3, a2a3 = PIC_Dihedron.get_dadad(acs)

        if not rev:
            hed1.len1 = set_accuracy_95(a0a1)
            hed1.len3 = set_accuracy_95(a1a2)

            hed2.len1 = set_accuracy_95(a1a2)
            hed2.len3 = set_accuracy_95(a2a3)
        else:
            hed1.len3 = set_accuracy_95(a0a1)
            hed1.len1 = set_accuracy_95(a1a2)

            hed2.len3 = set_accuracy_95(a1a2)
            hed2.len1 = set_accuracy_95(a2a3)

        hed1.angle2 = set_accuracy_95(a0a1a2)
        hed2.angle2 = set_accuracy_95(a1a2a3)

        hed1.updated = True
        hed2.updated = True

        self.updated = True

        # print(self)
        # print(hed1)
        # print(hed2)


class Hedron(DH_Base):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # print('initialising', self.id)
        if 'len1' in kwargs:
            # distance between 1st and 2nd atom
            self.len1 = float(kwargs['len1'])
            # angle formed between 3 atoms
            self.angle2 = float(kwargs['angle2'])
            # distance between 2nd and 3rd atoms
            self.len3 = float(kwargs['len3'])
        else:
            self.len1 = None
            self.angle2 = None
            self.len3 = None

        # 3 matrices specifying hedron space coordinates of constituent atoms,
        # initially atom3 on +Z axis
        self.atoms = []
        # 3 matrices, hedron space coordinates, reversed order
        # initially atom1 on +Z axis
        self.atomsR = []

        # print(self)

    def __str__(self):
        return ('3-' + self.id + ' ' + self.rdh_class + ' ' + str(self.len1)
                + ' ' + str(self.angle2) + ' ' + str(self.len3))

    def init_pos(self):

        # build hedron with a2 on +Z axis, a1 at origin,
        # a0 in -Z at angle n XZ plane
        for i in range(0, 3):
            # note this initializes a1 to 0,0,0 origin
            self.atoms.append(numpy.array([[0], [0], [0], [1]],
                                          dtype=numpy.float64))  # 4x1 array

        # supplement: angles which add to 180 are supplementary
        sar = numpy.deg2rad(180.0 - self.angle2)

        # a2 is len3 up from a2 on Z axis, X=Y=0
        self.atoms[2][2][0] = self.len3
        # a0 X is sin( sar ) * len1
        self.atoms[0][0][0] = numpy.sin(sar) * self.len1
        # a0 Z is -(cos( sar ) * len1)
        # (assume angle2 always obtuse, so a0 is in -Z)
        self.atoms[0][2][0] = - (numpy.cos(sar) * self.len1)

        # same again but 'reversed' : a0 on Z axis, a1 at origin, a2 in -Z
        for i in range(0, 3):
            # atom[1] to 0,0,0 origin
            self.atomsR.append(numpy.array([[0], [0], [0], [1]],
                                           dtype=numpy.float64))

        # a0r is len1 up from a1 on Z axis, X=Y=0
        self.atomsR[0][2][0] = self.len1
        # a2r X is sin( sar ) * len3
        self.atomsR[2][0][0] = numpy.sin(sar) * self.len3
        # a2r Z is -(cos( sar ) * len3)
        self.atomsR[2][2][0] = - (numpy.cos(sar) * self.len3)

        self.updated = False

    @staticmethod
    def get_dad(acs):
        a0 = acs[0].squeeze()
        a1 = acs[1].squeeze()
        a2 = acs[2].squeeze()

        a0a1 = numpy.linalg.norm(a0-a1)
        a1a2 = numpy.linalg.norm(a1-a2)

        a0a2 = numpy.linalg.norm(a0-a2)

        a0a1a2 = numpy.rad2deg(
            numpy.arccos(
                ((a0a1*a0a1) + (a1a2*a1a2) - (a0a2*a0a2)) / (2*a0a1*a1a2)
            )
        )
        return a0a1, a0a1a2, a1a2

    def hedron_from_atoms(self, atom_coords):
        aks = self.aks
        acs = []
        estr = ''
        for ak in aks:
            ac = atom_coords[ak]
            if ac is None:
                estr += ak + ' '
            else:
                acs.append(ac)
        if estr != '':
            raise HedronIncompleteError(
                self, 'missing coordinates for', estr)

        len1, angle2, len3 = Hedron.get_dad(acs)
        self.len1 = set_accuracy_95(len1)
        self.angle2 = set_accuracy_95(angle2)
        self.len3 = set_accuracy_83(len3)


# ------------ add to residue class


class PIC_Residue(object):

    def __init__(self, parent, NO_ALTLOC=False):  # , ALT_OK=True):
        # NO_ALTLOC=True will turn off alotloc positions and just use selected
        self.residue = parent
        #self.ndx = ndx
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
        # try:
        #    self.lc = three_to_one(parent.resname).upper()
        # except KeyError:
        #    self.lc = None
        # generated from dihedra include i+1 atoms, not residue.atoms
        # or initialised here from parent residue if loaded from coordinates
        self.atom_coords = {}
        # TODO: remove bp_atoms, just for testing against bp vector()
        self.bp_atoms = {}
        if not NO_ALTLOC:
            self.alt_ids = []

        is20AA = True
        rid = parent.id
        rbase = [rid[1],
                 rid[2] if ' ' != rid[2] else None,
                 parent.resname]
        try:
            rbase[2] = three_to_one(rbase[2]).upper()
        except KeyError:
            is20AA = False

        self.rbase = rbase
        self.lc = rbase[2]

        if is20AA:
            self.NCaCKey = (Atom_Key(self, 'N'),
                            Atom_Key(self, 'CA'),
                            Atom_Key(self, 'C'))

            for atom in parent.get_atoms():
                if hasattr(atom, 'child_dict'):
                    if NO_ALTLOC:
                        self.add_atom(atom.selected_child)
                    else:
                        for atm in atom.child_dict.values():
                            self.add_atom(atm)
                else:
                    self.add_atom(atom)

            # print(self.atom_coords)

    accept_atoms = ('N', 'CA', 'C', 'O', 'H', 'CB', 'CG', 'CG1', 'OG1', 'SG',
                    'CG2', 'CD', 'CD1', 'SD', 'OD1', 'ND1', 'CD2', 'ND2', 'CE',
                    'CE1', 'NE', 'OE1', 'NE1', 'CE2', 'OE2', 'NE2', 'CE3',
                    'CZ', 'NZ', 'CZ2', 'CZ3')

    def add_atom(self, atm):
        if atm.name not in self.accept_atoms:
            # print('skip:', atm.name)
            return
        ak = Atom_Key(self, atm)

        # TODO: are 4x1 arrays necessary?
        arr41 = numpy.append(atm.coord, [1])
        self.atom_coords[ak] = numpy.array(arr41, dtype=numpy.float64)[
            numpy.newaxis].transpose()
        self.bp_atoms[ak] = atm

    def __str__(self):
        return('(' + str(self.residue.id) + ')')

    def load_PIC(self, edron):
        ek = [edron['a1'], edron['a2'], edron['a3']]

        if edron['a4'] is not None:
            ek.append(edron['a4'])
            #self.dihedra[DH_Base.gen_key(ek)] = PIC_Dihedron(ek, **edron)
            self.dihedra[tuple(ek)] = PIC_Dihedron(ek, **edron)
        else:
            #self.hedra[DH_Base.gen_key(ek)] = Hedron(ek, **edron)
            self.hedra[tuple(ek)] = Hedron(ek, **edron)

    def link_dihedra(self):
        id3i = {}
        for dh in self.dihedra.values():
            dh.residue = self        # each dihedron can find its residue.pic
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
        for d in self.dihedra.values():
            if d.updated:
                d.init_pos()

    # def NCaCKeySplit(self):
    #    # rbase = self.ndx + self.lc
    #    return [self.rbase + 'N', self.rbase + 'CA', self.rbase + 'C']

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
        dbg = False

        transformations = {}
        atomCoords = atomCoordsIn
        NCaCKey = self.NCaCKey

        if genSCAD:
            transformations[NCaCKey] = numpy.identity(4, dtype=numpy.float64)
        #gk = DH_Base.gen_key
        q = deque([
            (Atom_Key(self, 'C'), Atom_Key(self, 'CA'), Atom_Key(self, 'N')),
            (Atom_Key(self, 'N'), Atom_Key(self, 'CA'), Atom_Key(self, 'CB')),
            (Atom_Key(self, 'O'), Atom_Key(self, 'C'), Atom_Key(self, 'CA')),
            self.NCaCKey
        ])

        # if list of initial coords not passed as parameter,
        # then use N-CA-C initial coords from creating dihedral
        if atomCoords is None:
            atomCoords = {}
            dlist = self.id3_dh_index[NCaCKey]
            for d in dlist:
                for i, a in enumerate(d.aks):
                    atomCoords[a] = d.initial_coords[i]

# TODO: optimise - need to check both len and not None????

        while q:  # deque is not empty
            if dbg:
                print('assemble loop start q=', q)
            h1k = q.pop()
            # print('  h1k:', h1k)
            dihedra = self.id3_dh_index.get(h1k, None)
            # print('  dihedra:', dihedra)
            if dihedra is not None:
                for d in dihedra:
                    if (4 == len(d.initial_coords)
                            and d.initial_coords[3] is not None):
                        # skip incomplete dihedron if don't have 4th atom due
                        # to missing input data
                        d_h2key = d.hedron2.aks
                        akl = d.aks
                        acount = 0
                        if dbg:
                            print('    process', d, d_h2key, akl)
                        for a in akl:
                            if (a in atomCoords and atomCoords[a] is not None):
                                acount += 1
                        if 4 == acount:
                            # dihedron already done, queue 2nd hedron key
                            q.appendleft(d_h2key)
                            # print('    4- already done, append left')
                        elif 3 == acount:
                            # TODO:[0, 1, x, 3] possible? then bug
                            # print('    3- call coord_space')
                            mt, mtr = coord_space(
                                [atomCoords[a] for a in akl[:3]], True)
                            if genSCAD:
                                transformations[h1k] = mtr
                            if dbg:
                                print(
                                    '        initial_coords[3]=', d.initial_coords[3].transpose())
                            acak3 = mtr @ d.initial_coords[3]
                            if dbg:
                                print('        acak3=', acak3.transpose())
                            for i in range(3):
                                acak3[i][0] = set_accuracy_83(acak3[i][0])
                            atomCoords[akl[3]] = acak3
                            if dbg:
                                print(
                                    '        3- finished, ak:', akl[3], 'coords:', atomCoords[akl[3]].transpose())
                            q.appendleft(d_h2key)
                        else:
                            print('no coords to start', d)
                            pass
                    else:
                        print('no initial coords for', d)
                        pass
        # print('coord_space returning')
        if genSCAD:
            return transformations
        else:
            return atomCoords

    altloc_re = re.compile(r'-([A-Z])\Z')

    def split_akl(self, lst):
        # given a list of Atom_Keys (aks) for a Hedron or Dihedron,
        #  return:
        #       list of matching aks that have id3_dh in this residue
        #            (ak may change if occupancy != 1.00)
        #    or
        #       multiple lists of matching aks expanded for all atom altlocs
        #    or
        #       empty list if any of atom_coord(ak) missing

        altloc_ndx = Atom_Key.fields['altloc']

        # step 1
        # given a list of Atom_Keys (aks)
        #  form a new list of same aks with coords in this residue
        #      plus lists of matching altloc aks in this residue
        edraLst = []
        altlocs = set()
        for ak in lst:
            if ak in self.atom_coords:
                edraLst.append([ak])
            else:
                ak2_lst = []
                for ak2 in self.atom_coords:
                    if ak.altloc_match(ak2):
                        # print(key)
                        ak2_lst.append(ak2)
                        altloc = ak2.akl[altloc_ndx]
                        if altloc is not None:
                            altlocs.add(altloc)
                edraLst.append(ak2_lst)

        # step 2
        # check and finish for
        #   missing atoms
        #   simple case no altlocs
        #
        maxc = 0
        for akl in edraLst:
            lenAKL = len(akl)
            if 0 == lenAKL:
                return []    # atom missing in atom_coords, cannot form object
            elif maxc < lenAKL:
                maxc = lenAKL
        if 1 == maxc:        # simple case no altlocs for any ak in list
            newAKL = []
            for akl in edraLst:
                newAKL.append(akl[0])
            return [newAKL]
        else:
            new_edraLst = []
            for al in altlocs:
                # form complete new list for each altloc
                alhl = []
                for akl in edraLst:
                    if 1 == len(akl):
                        alhl.append(akl[0])  # not all atoms will have altloc
                    else:
                        for ak in akl:
                            if ak.akl[altloc_ndx] == al:
                                alhl.append(ak)
                new_edraLst.append(alhl)
            # print(new_edraLst)
            return new_edraLst

        # can't reach here

    def gen_edra(self, lst):
        # given list of Atom_Keys defining hedron or dihedron
        #  convert to Atom_Keys with coordinates in this residue
        #  add appropriately to self.di/hedra, expand as needed atom altlocs

        if 4 > len(lst):
            dct, obj = self.hedra, Hedron
        else:
            dct, obj = self.dihedra, PIC_Dihedron

        if all(ak in self.atom_coords for ak in lst):
            hl = [lst]
        else:
            hl = self.split_akl(lst)

        for nlst in hl:
            #dct[DH_Base.gen_key(nlst)] = obj(**zdh(nlst))
            dct[tuple(nlst)] = obj(**zdh(nlst))

    def dihedra_from_atoms(self):
        AK = Atom_Key
        S = self

        sN, sCA, sC = AK(S, 'N'), AK(S, 'CA'), AK(S, 'C')
        sO, sCB = AK(S, 'O'), AK(S, 'CB')
        sN = Atom_Key(self, 'N')

        if self.rnext:
            # atom_coords, hedra and dihedra for backbone dihedra
            # which reach into next residue
            rn = self.rnext
            nN, nCA, nC = AK(rn, 'N'), AK(rn, 'CA'), AK(rn, 'C')

            for ak in (nN, nCA, nC):
                if ak in rn.atom_coords:
                    self.atom_coords[ak] = rn.atom_coords[ak]
                else:
                    for rn_ak in rn.atom_coords.keys():
                        if rn_ak.altloc_match(ak):
                            self.atom_coords[rn_ak] = rn.atom_coords[rn_ak]

            self.gen_edra([sN, sCA, sC, nN])  # psi
            self.gen_edra([sCA, sC, nN, nCA])  # omega i+1
            self.gen_edra([sC, nN, nCA, nC])  # phi i+1
            self.gen_edra([sCA, sC, nN])
            self.gen_edra([sC, nN, nCA])
            self.gen_edra([nN, nCA, nC])

        # backbone O and C-beta hedra and dihedra within this residue
        self.gen_edra([sN, sCA, sC, sO])
        self.gen_edra([sO, sC, sCA, sCB])

        if not self.rprev:
            self.gen_edra([sN, sCA, sC])

        self.gen_edra([sCA, sC, sO])
        self.gen_edra([sCB, sCA, sC])

        # if res is not G or A
        # only needed for sidechain CG residues
        # (not gly or ala or any missing rest of side chain)

        if 'G' != self.lc and 'A' != self.lc:
            for ak in self.atom_coords:
                if ak.is_sidechain_gamma():
                    self.gen_edra([sN, sCA, sCB])
                    break

        # amide proton N H if present
        sH = AK(S, 'H')
        self.gen_edra([sH, sN, sCA])
        self.gen_edra([sC, sCA, sN, sH])

        # sidechain hedra and dihedra

        sidechain = pic_data_sidechains.get(self.lc, [])
        for edra in sidechain:
            r_edra = [AK(S, atom) for atom in edra]
            self.gen_edra(r_edra[0:4])  # [4] is label on some table entries

        if sCB not in self.atom_coords:  # add C-beta for Gly
            cb = numpy.append(genCBjones(self.residue), [1])
            self.atom_coords[sCB] = numpy.array(cb, dtype=numpy.float64)[
                numpy.newaxis].transpose()

        # testing C-beta generation against Ala:
        # if 'A' == self.lc:
        #    aj = genCBjones(self.residue)
        #    ah = genCBhamelryck(self.residue)
        #    cac = self.residue['CB'].coord
        #    cav = self.residue['CB'].get_vector()
        #    dj = aj - cac
        #    dh = ah - cav
        #    print('deltaJ=',dj,'deltaH=',dh,aj,cac)

        self.link_dihedra()

        for d in self.dihedra.values():
            # populate values and hedra for dihedron ojects
            d.dihedron_from_atoms()
        for h in self.hedra.values():
            # miss redundant hedra above, needed for some chi1 angles
            if h.len1 is None:
                # print(h)
                h.hedron_from_atoms(self.atom_coords)

    def pdb_atom_string(self, atm):
        res = atm.parent
        chn = res.parent
        s = ('{:6}{:5d} {:4}{:1}{:3} {:1}{:4}{:1}   {:8.3f}{:8.3f}{:8.3f}'
             '{:6.2f}{:6.2f}        {:>4}\n'
             ).format('ATOM', atm.serial_number, atm.fullname, atm.altloc,
                      res.resname, chn.id, res.id[1], res.id[2], atm.coord[0],
                      atm.coord[1], atm.coord[2], atm.occupancy, atm.bfactor,
                      atm.element)
        # print(s)
        return s

    @staticmethod
    def residue_string(res):
        return (str(res.get_full_id())
                + ' ' + res.resname
                + ' [' + res.get_segid() + ']'
                #+ (' ' + str(ndx) if ndx is not None else '')
                + '\n')

    def write_PIC(self, pdbid, chainid, s='', stats=False):
        s += PIC_Residue.residue_string(self.residue)
        if not self.rprev:
            try:
                ts = self.pdb_atom_string(self.residue['N'])
                ts += self.pdb_atom_string(self.residue['CA'])
                ts += self.pdb_atom_string(self.residue['C'])
                s += ts  # only if have all 3 atoms
            except KeyError:
                pass
        # TODO: lose pdbid .. chainid format after lua compare?
        base = pdbid + ' ' + chainid + ' '
        for h in sorted(self.hedra.values()):
            try:
                s += (base + h.id + ' ' + "{:9.5f} {:9.5f} {:9.5f}".format(
                    h.len1, h.angle2, h.len3))
            except KeyError:
                pass
            s += '\n'
        for d in sorted(self.dihedra.values()):
            try:
                s += (base + d.id + ' ' + "{:9.5f}".format(d.dihedral1))
            except KeyError:
                pass
            s += '\n'
        return s

    def coords_to_residue(self):
        seqpos, icode = self.residue.id[1:3]
        seqpos = str(seqpos)
        spNdx, icNdx, atmNdx, occNdx, altlocNdx = [
            Atom_Key.fields.get(k) for k in
            ['respos', 'icode', 'atm', 'occ', 'altloc']
        ]
        Res = self.residue
        ndx = Res.parent.pic.ndx

        for ak in sorted(self.atom_coords):
            # print(ak)
            if (seqpos == ak.akl[spNdx] and
                ((icode == ' ' and ak.akl[icNdx] is None)
                 or icode == ak.akl[icNdx])):

                ac = self.atom_coords[ak]
                atm_coords = ac[:3].transpose()[0]
                akl = ak.akl
                atm, altloc = akl[atmNdx], akl[altlocNdx]

                Atm = None
                newAtom = None

                if Res.has_id(atm):
                    Atm = Res[atm]

                if (Atm is None
                    or (2 == Atm.is_disordered()
                        and not Atm.disordered_has_id(altloc))):
                    #print('new', ak)
                    newAtom = Atom(atm, atm_coords, 0.0,
                                   (1.00 if akl[occNdx] is None else akl[occNdx]),
                                   ' ', atm, ndx, atm[0])
                    ndx += 1
                    if Atm is None:
                        if altloc is None:
                            Res.add(newAtom)
                        else:
                            disordered_atom = DisorderedAtom(atm)
                            Res.add(disordered_atom)
                            disordered_atom.disordered_add(newAtom)
                            Res.flag_disordered()
                    else:
                        Atm.disordered_add(newAtom)
                else:
                    # Atm is not None, might be disordered with altloc
                    #print('update', ak)
                    if (2 == Atm.is_disordered()
                            and Atm.disordered_has_id(altloc)):
                        Atm.disordered_select(altloc)
                    Atm.set_coord(atm_coords)
                    ndx = Atm.get_serial_number()

        Res.parent.pic.ndx = ndx

        # print(ak, ac)


# Residue.pic_init = residue_pic_init
# Residue.pic_load = residue_pic_load
# Residue.link_dihedra = residue_link_dehedra
# Residue.render_hedra = residue_render_hedra
# Residue.render_dihedra = residue_render_dihedra
# Residue.NCaCKeySplit = residue_NCaCKeySplit
# Residue.assemble = residue_assemble

# ------------ add to chain class

def npat(arr):
    arr41 = numpy.append(arr, [1])
    return numpy.array(arr41, dtype=numpy.float64)[
        numpy.newaxis].transpose()


class PIC_Chain:

    def __init__(self, parent):
        self.chain = parent
        self.ordered_aa_pic_list = []
        self.initNCaC = {}
        self.set_residues()  # no effect if no residues loaded

    def set_residues(self):
        #ndx = 0
        last_res = None
        last_ord_res = None
        for res in self.chain.get_residues():
            # select only not hetero
            # if res.id[0] == ' ':  # and not res.disordered:
            if not hasattr(res, 'pic'):
                res.pic = PIC_Residue(res)  # , ndx)
            if res.id[0] == ' ':  # TODO: in set of supported hetatms?
                self.ordered_aa_pic_list.append(res.pic)
                if last_res is not None and last_ord_res == last_res:
                    # no missing or hetatm
                    last_ord_res.pic.rnext = res.pic
                    res.pic.rprev = last_ord_res.pic
                else:
                    # chain break, stash coords to restart
                    respos = str(res.id[1])
                    if ' ' != res.id[2]:
                        respos += res.id[2]  # need icode if set
                    initNCaC = {}
                    rpic = res.pic
                    initNCaC[Atom_Key(rpic, 'N')] = npat(res['N'].coord)
                    initNCaC[Atom_Key(rpic, 'CA')] = npat(res['CA'].coord)
                    initNCaC[Atom_Key(rpic, 'C')] = npat(res['C'].coord)
                    self.initNCaC[respos] = initNCaC
                last_ord_res = res
            else:
                # print('skipping res ' + str(res.id) + ' ' + res.resname
                #       + (' disordered'
                #           if res.disordered else ' not disordered'))
                pass
            last_res = res

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
            respos = rpic.residue.id[1]
            if (not start or start <= respos) and (not fin or fin >= respos):
                if rpic.rprev:
                    startPos = {}
                    rp = rpic.rprev
                    # nb akl for this res n-ca-c in rp dihedra
                    akl = rpic.NCaCKey  # .split(':')  # Split()
                    for ak in akl:
                        startPos[ak] = rp.atom_coords[ak]
                else:
                    # in theory the atom posns already added by load_structure
                    icode = rpic.residue.id[2]
                    NCaCselect = str(respos) + '' if icode == ' ' else icode
                    startPos = self.initNCaC[NCaCselect]
                # rtm unlike lua assemble() mods res coords TODO: correct comment or change back
                rpic.atom_coords = rpic.assemble(startPos)

    def dihedra_from_atoms(self):
        for rpic in self.ordered_aa_pic_list:
            rpic.dihedra_from_atoms()

    def coords_to_structure(self):
        self.ndx = 0
        for res in self.chain.get_residues():
            if hasattr(res, 'pic'):
                res.pic.coords_to_residue()


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
        pdb_structure = load_pic_file(
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
        pic_str = write_PIC(pdb_structure)
        with open(target + '.PyPIC', 'w') as outf:
            print(pic_str, file=outf)
        print('wrote pic output for', target)
print('normal termination')
