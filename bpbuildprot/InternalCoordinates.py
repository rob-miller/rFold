#!/usr/local/bin/python3
"""Interconversion between PDB atom and internal coordinates.

Atom coordinates are cartesian X, Y, Z coordinates as specified in a PDB file.
Internal coordinates are dihedral angles, bond angles and bond lengths.
These routines compute internal coordinates from atom coordinates, read/write
them as files, and regenerate PDB atom coordinate files from internal
coordinates.  Also supported is writing PDB internal coordinate data as
OpenSCAD code to generate STL output for printing protein structures on a
3D printer.
"""

# from Bio import PDB

from collections import deque, namedtuple
import sys
import re

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
from Bio.PDB.vectors import rotaxis2m   # , calc_dihedral

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy to build proteins from internal coordinates.")

from PIC_Data import pic_data_sidechains

from Bio.PDB.PDBExceptions import PDBException

# __updated__ is specifically for the python-coding-tools Visual Studio Code
# extension, which updates the variable on each file save

__updated__ = '2019-03-21 17:38:24'
print('ver: ' + __updated__)
print(sys.version)


def gen_Mrot(angle_rads, axis):
    """Generate a 4x4 single-axis numpy rotation matrix.

    :param float angle_rads: the desired rotation angle in radians
    :param char axis: character specifying the rotation axis
    """
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
    """Generate a 4x4 numpy translation matrix.

    :param list[3] xyz: translation in each axis
    """
    return numpy.array([[1, 0, 0, xyz[0]],
                        [0, 1, 0, xyz[1]],
                        [0, 0, 1, xyz[2]],
                        [0, 0, 0, 1]
                        ], dtype=numpy.float64)


def atom_to_internal_coordinates(struct):
    """Create/update internal coordinates from Atom X,Y,Z coordinates.

    Internal coordinates are bond length, angle and dihedral angles.

    :param Structure struct: Biopython PDB Structure object
    """
    for chn in struct.get_chains():
        if not hasattr(chn, 'pic'):
            chn.pic = PIC_Chain(chn)
        chn.pic.dihedra_from_atoms()


def internal_to_atom_coordinates(struct):
    """Create/update atom coordinates from internal coordinates.

    :param Structure struct: Biopython PDB Structure object with .pic
        attribute containing internal coordinates (bond length, angle,
        and dihedral angles)
    :raises Exception: if any chain does not have .pic attribute
    """
    for chn in struct.get_chains():
        if hasattr(chn, 'pic'):
            chn.pic.internal_to_atom_coordinates()
        else:
            raise Exception(
                'Structure %s Chain %s does not have internal coordinates set'
                % (struct, chn))


def read_PIC(file):
    """Load Protein Internal Coordinate (PIC) data from file.

    PIC file format:
        (optional) PDB HEADER record
            - idcode and deposition date recommended but optional
            - deposition date in PDB format or as changed by Biopython
        (optional) PDB TITLE record
        repeat:
            Biopython Residue Full ID - sets ID of returned structure
            (optional) PDB ATOM records for chain start N, CA, C
            PIC Hedra records for residue
            PIC Dihedra records for residue

    :param Bio.File file: file name or handle
    :returns: Biopython Structure object, Residues with .pic attributes
        but no coordinates except for chain start N, CA, C atoms if supplied,
        or None on parse fail (silent, no exception rasied)
    """
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
    bfac_re = re.compile(r'^BFAC:\s([^\s]+\s+[\-\d\.]+)'
                         r'\s*([^\s]+\s+[\-\d\.]+)?'
                         r'\s*([^\s]+\s+[\-\d\.]+)?'
                         r'\s*([^\s]+\s+[\-\d\.]+)?'
                         r'\s*([^\s]+\s+[\-\d\.]+)?')
    bfac2_re = re.compile(r'([^\s]+)\s+([\-\d\.]+)')
    struct_builder = StructureBuilder()

    # init empty header dict
    # - could use to parse HEADER and TITLE lines except
    #   deposition_date format changed from original PDB header
    header_dict = _parse_pdb_header_list([])

    curr_SMCS = [None, None, None, None]  # struct model chain seg
    SMCS_init = [
        struct_builder.init_structure,
        struct_builder.init_model,
        struct_builder.init_chain,
        struct_builder.init_seg
    ]

    sb_res = None

    with as_handle(file, mode='r') as handle:
        for aline in handle.readlines():
            if aline.startswith('HEADER '):
                m = pdb_hdr_re.match(aline)
                if m:
                    header_dict['head'] = m.group('cf')  # classification
                    header_dict['idcode'] = m.group('id')
                    header_dict['deposition_date'] = m.group('dd')
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
                    # check SMCS = Structure, Model, Chain, SegID
                    this_SMCS = [m.group(1), int(m.group(2)),
                                 m.group(3), m.group(8)]
                    if curr_SMCS != this_SMCS:
                        # init new SMCS level as needed
                        for i in range(4):
                            if curr_SMCS[i] != this_SMCS[i]:
                                SMCS_init[i](this_SMCS[i])
                                curr_SMCS[i] = this_SMCS[i]
                                if 0 == i:
                                    # 0 = init structure so add header
                                    struct_builder.set_header(header_dict)

                    struct_builder.init_residue(
                        m.group('res'), m.group('het'),
                        int(m.group('pos')), m.group('icode'))

                    sb_res = struct_builder.residue
                    sb_res.pic = PIC_Residue(sb_res, m.group('ndx'))
                    # print('res id:', m.groupdict())
                    # print(report_PIC(struct_builder.get_structure()))
                else:
                    print('Reading pic file', file, 'residue fail: ', aline)
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
                        raise Exception('pic ATOM read confusion: %s %s %s'
                                        % (sb_res.resname, str(sb_res.id),
                                            aline))

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
            elif aline.startswith('BFAC: '):
                m = bfac_re.match(aline)
                if m:
                    for bfac_pair in m.groups():
                        if bfac_pair is not None:
                            m2 = bfac2_re.match(bfac_pair)
                            if (m2
                                and sb_res is not None
                                    and hasattr(sb_res, 'pic')):
                                rp = sb_res.pic
                                rp.bfactors[m2.group(1)] = float(m2.group(2))
            else:
                m = Edron_Base.edron_re.match(aline)
                if m:
                    sb_res.pic.load_PIC(m.groupdict())
                elif aline.strip():
                    print('Reading PIC file', file,
                          'parse fail on: .', aline, '.')
                    return None

    struct = struct_builder.get_structure()
    for chn in struct.get_chains():
        chnp = chn.pic = PIC_Chain(chn)
        chnp.link_residues()
        chnp.render_dihedra()

    # print(report_PIC(struct_builder.get_structure()))
    return struct


def write_PIC(entity, file, pdbid=None, chainid=None):
    """Write Protein Internal Coordinates (PIC) to file.

    See read_PIC() for file format.  Recurses to lower entity levels (M, C, R).

    :param Entity entity: Biopython PDB Entity object: S, M, C or R
    :param Bio.File file: file name or handle
    :param str pdbid: PDB idcode, set to 0PDB if cannot be determined
    :param char chainid: PDB Chain ID, set from C level entity.id if needed
    :raises PDBException: if entity level not S, M, C, or R
    :raises Exception: if entity does not have .level attribute
    """
    with as_handle(file, 'w') as fp:
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
                    fp.write(entity.pic.write_PIC(pdbid, chainid))
                else:
                    fp.write(PIC_Residue._residue_string(entity))
            elif 'C' == entity.level:
                if not chainid:
                    chainid = entity.id
                for res in entity:
                    write_PIC(res, fp, pdbid, chainid)
                    pass
            elif 'M' == entity.level:
                for chn in entity:
                    write_PIC(chn, fp, pdbid, chainid)
            elif 'S' == entity.level:
                if not pdbid:
                    pdbid = entity.header.get('idcode', None)
                hdr = entity.header.get('head', None)
                dd = entity.header.get('deposition_date', None)
                if hdr:
                    fp.write(('HEADER    {:40}{:8}   {:4}\n'
                              ).format(hdr.upper(), (dd or ''), (pdbid or '')))
                nam = entity.header.get('name', None)
                if nam:
                    fp.write('TITLE     ' + nam.upper() + '\n')
                for mdl in entity:
                    write_PIC(mdl, fp, pdbid, chainid)
            else:
                raise PDBException("Cannot identify level: "
                                   + str(entity.level))
        except KeyError:
            raise Exception(
                "write_PIC: argument is not a Biopython PDB Entity "
                + str(entity))


def report_PIC(entity, reportDict=None):
    """Generate dict with counts of PIC data elements for each entity level.

    :param Entity entity: Biopython PDB Entity object: S, M, C or R
    :raises PDBException: if entity level not S, M, C, or R
    :raises Exception: if entity does not have .level attribute
    :returns: dict with counts of PIC data elements
    """
    if reportDict is None:
        reportDict = {'idcode': None, 'hdr': 0, 'mdl': 0, 'chn': 0,
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
            if reportDict['idcode'] is None:
                reportDict['idcode'] = entity.header.get('idcode', None)

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


def genCBhamelryck(res):
    """Generate virtual C-beta residue, Hamelryck method.

    Method from
    https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
    'How do I put a virtual C-beta on a Gly residue?'

    :param Residue res: Biopython Residue object
    :returns: numpy 1x3 array, coordinates of virtual C-beta
    """
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
    """Generate virtual C-beta residue, Hazes and Dijkstra / Jones method.

    Based on Prot. Engr vol. 2, issue 2, p. 121 (1 July, 1988)
    'Model building of disulfide bonds in proteins with known three-
    dimensional structure', Hazes and Dijkstra, pp 119-125
    https://doi.org/10.1093/protein/2.2.119
    https://academic.oup.com/peds/article/2/2/119/1484803
    Original code by DT Jones
    TODO: update the constants!  Choose for Ala

    :param Residue res: Biopython Residue object
    :returns: numpy 1x3 array, coordinates of virtual C-beta
    """
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
    """Reduce floating point accuracy to 9.5 (xxxx.xxxxx).

    :param float num: input number
    :returns: float with specified accuracy
    """
    return float("{:9.5f}".format(num))


def set_accuracy_83(num):
    """Reduce floating point accuracy to 8.3 (xxxxx.xxx).

    :param float num: input number
    :returns: float with specified accuracy
    """
    return float("{:8.3f}".format(num))


def get_spherical_coordinates(xyz):
    """Compute spherical coordinates (r, theta, phi) for X,Y,Z point.

    :param array xyz: column vector (3 row x 1 column numpy array)
    :return: tuple of r, theta, phi for input coordinate
    """
    r = numpy.linalg.norm(xyz)
    if 0 == r:
        return numpy.array([0, 0, 0])
    sign = -1.0 if xyz[1][0] < 0.0 else 1.0
    theta = ((math.pi/2.0 * sign) if 0 == xyz[0][0]
             else numpy.arctan2(xyz[1][0], xyz[0][0]))
    phi = numpy.arccos(xyz[2][0]/r)
    return (r, theta, phi)


def coord_space(acs, rev=False):
    """Generate transformation matrix to coordinate space defined by 3 points.

    New coordinate space will have:
        acs[0] on XZ plane
        acs[1] origin
        acs[2] on +Z axis

    :param numpy column array x3 acs: X,Y,Z column input coordinates x3
    :param bool rev: if True, also return reverse transformation matrix
        (to return from coord_space)
    :returns: 4x4 numpy array, x2 if rev=True
    """
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
    # mt = mry @ mrz @ tm  # python 3.5 and later
    mt = mry.dot(mrz.dot(tm))

    if dbg:
        # print('mt * a2', (mt @ a2).transpose())
        print('mt * a2', (mt.dot(a2)).transpose())

    # p = mt @ a0
    p = mt.dot(a0)

    # need theta of translated a0
    # sc2 = get_spherical_coordinates(p)
    sign = -1.0 if (p[1][0] < 0.0) else 1.0
    theta2 = ((math.pi/2.0 * sign) if 0 == p[0][0]
              else numpy.arctan2(p[1][0], p[0][0]))
    # rotate a0 -theta2 about Z to align with X
    mrz2 = gen_Mrot(-theta2, 'z')

    # mt = mrz2 @ mt
    mt = mrz2.dot(mt)

    if not rev:
        return mt

    # rev=True, so generate the reverse transformation

    # rotate a0 theta about Z, reversing alignment with X
    mrz2 = gen_Mrot(theta2, 'z')
    # rotate a2 phi about Y
    mry = gen_Mrot(sc[2], 'y')
    # rotate a2 theta about Z
    mrz = gen_Mrot(sc[1], 'z')
    # translation matrix origin to a1
    tm = gen_Mtrans([a1[0], a1[1], a1[2]])

    # mr = tm @ mrz @ mry @ mrz2
    mr = tm.dot(mrz.dot(mry.dot(mrz2)))

    return mt, mr


class AtomKey(object):
    """Class for dict keys to reference atom coordinates.

    Supports rich comparison and multiple ways to instantiate.
    AtomKeys contain:
     residue position, insertion code, 1 or 3 char residue name,
     atom name, occupancy, and altloc

    Attributes
    ----------
    akl: tuple
        All six fields of AtomKey
    fieldNames : tuple (Class Attribute)
        Mapping of key index positions to names
    fields : namedtuple (Class Attribute)
        Mapping of field names to index positions
    id: str
        '_'-joined AtomKey fields, excluding 'None' fields
    atom_re : compiled regex (Class Attribute)
        A compiled regular expression matching the string form of the key

    Methods
    -------
    altloc_match(other)
        Returns True if this AtomKey matches other AtomKey excluding altloc
        and occupancy fields
    is_sidechain_gamma()
        Returns True if atom name contains 'G', meaning it is gamma element
        of sidechain

    """

    atom_re = re.compile(r'^(?P<respos>-?\d+)(?P<icode>[A-Za-z])?'
                         r'_(?P<resname>[a-zA-Z]+)_(?P<atm>\w+)'
                         r'(?:_(?P<occ>\d\.\d?\d?))?(?:_(?P<altloc>\w))?$')
    # PDB altLoc = Character = [\w ] (any non-ctrl ASCII incl space)
    # PDB iCode = AChar = [A-Za-z]
    fieldNames = ('respos', 'icode', 'resname', 'atm', 'occ', 'altloc')
    fields = namedtuple('fieldsDef', 'respos, icode, resname, '
                        'atm, occ, altloc')(0, 1, 2, 3, 4, 5)

    def __init__(self, *args, **kwargs):
        """Initialize AtomKey with residue and atom data.

        Examples of acceptable input:
            (<PIC_Residue>, 'CA', ...)    : PIC_Residue with atom info
            (<PIC_Residue>, <Atom>)       : PIC_Residue with Biopython Atom
            ([52, None, 'G', 'CA', ...])  : list of ordered data fields
            (52, None, 'G', 'CA', ...)    : multiple ordered arguments
            ({respos: 52, icode: None, atm: 'CA', ...}) : dict with fieldNames
            (respos: 52, icode: None, atm: 'CA', ...) : kwargs with fieldNames
            52_G_CA, 52B_G_CA, 52_G_CA_0.33, 52_G_CA_0.33_B  : id strings
        """
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
                occ = float(arg.occupancy)
                akl.append(occ if occ != 1.00 else None)
                altloc = arg.altloc
                akl.append(altloc if altloc != ' ' else None)
            elif isinstance(arg, list):
                akl += arg
            elif isinstance(arg, dict):
                for k in self.fieldNames:
                    akl.append(arg.get(k, None))
            elif '_' in arg:
                # got atom key string, recurse with regex parse
                m = self.atom_re.match(arg)
                if [] != akl:
                    raise Exception(
                        'Atom Key init full key not first argument: ' + arg)
                for fn in AtomKey.fieldNames:
                    akl.append(m.group(fn))
            else:
                akl.append(arg)

        # process kwargs, initialize occ and altloc to None
        # if not specified above
        for i in range(6):
            if len(akl) <= i:
                akl.append(kwargs.get(self.fieldNames[i], None))

        # tweak local akl to generate id string
        akl[0] = str(akl[0])  # numeric residue position to string

        if akl[4] is not None:
            akl[4] = str(akl[4])  # numeric occupancy to string

        # unused option:
        # (self.respos, self.icode, self.resname, self.atm, self.occ,
        #    self.altloc) = akl

        self.id = '_'.join(
            [''.join(filter(None, akl[:2])),
                akl[2],
                '_'.join(filter(None, akl[3:]))]
        )

        self.akl = tuple(akl)
        self._hash = hash(self.akl)

        # pass

    def __repr__(self):
        """Repr string from id."""
        return self.id

    def __hash__(self):
        """Hash calculated at init from akl tuple."""
        return self._hash

    _backbone_sort_keys = {'N': 0, 'CA': 1, 'C': 2, 'O': 3}
    _sidechain_sort_keys = {'CB': 1,
                            'CG': 2, 'CG1': 2, 'OG': 2, 'OG1': 2, 'SG': 2,
                            'CG2': 3,
                            'CD': 4, 'CD1': 4, 'SD': 4, 'OD1': 4, 'ND1': 4,
                            'CD2': 5, 'ND2': 5, 'OD2': 5,
                            'CE': 6, 'NE': 6, 'CE1': 6, 'OE1': 6, 'NE1': 6,
                            'CE2': 7, 'OE2': 7, 'NE2': 7,
                            'CE3': 8,
                            'CZ': 9, 'CZ2': 9, 'NZ': 9,
                            'NH1': 10, 'OH': 10, 'CZ3': 10,
                            'CH2': 11, 'NH2': 11,
                            'OXT': 12,
                            'H': 13
                            }

    def altloc_match(self, other):
        """Test AtomKey match other discounting occupancy and altloc."""
        if isinstance(other, type(self)):
            return self.akl[:4] == other.akl[:4]
        else:
            return NotImplemented

    def is_sidechain_gamma(self):
        """Test atom name contains 'G', then it is sidechain gamma atom."""
        if 'G' in self.akl[3]:
            return True
        return False

    def _cmp(self, other):
        """Comparison function ranking self vs. other."""
        akl_s = self.akl
        akl_o = other.akl
        atmNdx = self.fields.atm
        occNdx = self.fields.occ
        for i in range(6):
            s, o = akl_s[i], akl_o[i]
            if s != o:
                if atmNdx != i:
                    # only sorting complications at atom level, occ.
                    # otherwise seqpos, insertion code will trigger
                    # before residue name
                    if occNdx == i:
                        tmp = s
                        s = o
                        o = tmp  # swap so higher occupancy comes first
                    if s is None:  # no insert code before named insert code
                        return 0, 1
                    elif o is None:
                        return 1, 0
                    else:
                        return s, o
                # backbone atoms before sidechain atoms
                sb = self._backbone_sort_keys.get(s, None)
                ob = self._backbone_sort_keys.get(o, None)
                if (sb is not None and ob is not None):
                    return sb, ob
                elif (sb is not None and ob is None):
                    return 0, 1
                elif (sb is None and ob is not None):
                    return 1, 0
                # finished backbone and backbone vs. sidechain atoms
                # now hydrogens after sidechain
                s0, o0 = s[0], o[0]
                if (s0 == 'H' and o0 != 'H'):
                    return 1, 0
                elif (s0 != 'H' and o0 == 'H'):
                    return 0, 1
                ss = self._sidechain_sort_keys.get(s, None)
                os = self._sidechain_sort_keys.get(o, None)
                if (ss is not None and os is not None):
                    return ss, os
                elif (ss is not None and os is None):
                    return 0, 1
                elif (ss is None and os is not None):
                    return 1, 0
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


class Edron_Base(object):
    """Base class for Hedron and Dihedron classes.

    Supports rich comparison based on lists of AtomKeys.

    Attributes
    ----------
    aks : tuple
        3 (hedron) or 4 (dihedron) AtomKeys defining this di/hedron
    id : str
        ':'-joined string of AtomKeys for this di/hedron
    atoms_updated : bool
        indicates hedron local atom_coords reflect current di/hedron angle and
        length values in hedron local coordinate space
    dh_class : str
        sequence of atoms (no position or residue) comprising di/hedron
        for statistics
    rdh_class : str
        sequence of residue, atoms comprising di/hedron for statistics
    edron_re : compiled regex (Class Attribute)
        A compiled regular expression matching string IDs for Hedron
        and Dihedron objects

    Methods
    -------
    gen_key([AtomKey, ...] or AtomKey, ...) (Static Method)
        generate a ':'-joined string of AtomKey Ids
    gen_acs(atom_coords)
        generate tuple of atom coords for keys in self.aks

    """

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

    @staticmethod
    def gen_key(lst):
        """Generate string of ':'-joined AtomKey strings from input.

        :param lst: list of AtomKey objects or id strings
        """
        if isinstance(lst[0], AtomKey):
            return ':'.join(ak.id for ak in lst)
        else:
            return ':'.join(lst)

    def __init__(self, *args, **kwargs):
        """Initialize Edron_Base with sequence of AtomKeys.

        Acceptable input:

            [ atom key, ... ]  : list of AtomKeys
            atom key, ...      : sequence of AtomKeys as args
            {'a1': str, 'a2': str, ... }  : dict of AtomKeys as 'a1', 'a2' ...
        """
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
            aks = [kwargs['a1'], kwargs['a2'], kwargs['a3']]

            try:
                if kwargs['a4'] is not None:
                    aks.append(kwargs['a4'])
            except KeyError:
                pass

        # if args are atom key strings instead of AtomKeys
        for i in range(len(aks)):
            if not isinstance(aks[i], AtomKey):
                aks[i] = AtomKey(aks[i])

        self.aks = tuple(aks)
        self.id = Edron_Base.gen_key(aks)
        self._hash = hash(self.aks)

        # flag indicating that atom coordinates are up to date
        # (do not need to be recalculated from dihedral1)
        self.atoms_updated = False

        # no residue or position, just atoms
        self.dh_class = ''
        # same but residue specific
        self.rdh_class = ''

        atmNdx = AtomKey.fields.atm
        resNdx = AtomKey.fields.resname
        for ak in aks:
            akl = ak.akl
            self.dh_class += akl[atmNdx]
            self.rdh_class += akl[resNdx] + akl[atmNdx]

    def gen_acs(self, atom_coords):
        """Generate tuple of atom coord arrays for keys in self.aks.

        :param atom_coords: AtomKey dict of atom coords for residue
        :raises: MissingAtomError any atoms in self.aks missing coordinates
        """
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
            raise MissingAtomError(
                '%s missing coordinates for %s' % (self, estr))
        return tuple(acs)

    def __repr__(self):
        """Tuple of AtomKeys is default repr string."""
        return self.aks

    def __hash__(self):
        """Hash calculated at init from aks tuple."""
        return self._hash

    def _cmp(self, other):
        """Comparison function ranking self vs. other."""
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


class Dihedron(Edron_Base):
    """Class to represent four joined atoms forming a dihedral angle.

    Attributes
    ----------
    dihedral1 : float
        Measurement or specification of dihedral angle
    hedron1, hedron2 : Hedron object references
        The two hedra which form the dihedral angle
    h1key, h2key : tuples of AtomKeys
        Hash keys for hedron1 and hedron2
    id3,id32 : tuples of AtomKeys
        First 3 and second 3 atoms comprising dihedron; hxkey orders may differ
    initial_coords : tuple[4] of numpy arrays [4][1]
        Local atom coords for 4 atoms, [0] on XZ plane, [1] at origin,
        [2] on +Z, [3] rotated by dihedral1
    a4_pre_rotation : numpy array [4][1]
        4th atom of dihedral aligned to XZ plane (dihedral1 not applied)
    pic_residue : PIC_Residue object reference
        PIC_Residue object containing this dihedral
    reverse : bool
        Indicates order of atoms in dihedron is reversed from order of atoms
        in hedra (configured by _set_hedra())

    Methods
    -------
    init_pos()
        Find Hedron objects for self.pic_residue, set initial_coords
        and a4_pre_rotation
    dihedron_from_atoms()
        Compute dihedral and bond lengths, angles from PIC_Residue atom_coords
    set_dihedral()
        Store new dihedral angle and update initial_coords accordingly

    """

    def __init__(self, *args, **kwargs):
        """Initialize Dihedron with sequence of AtomKeys and optional dihedral angle.

        Acceptable input:
            As for Edron_Base, plus optional 'dihedral1' keyworded angle value.
        """
        super().__init__(*args, **kwargs)

        # hedra making up this dihedron; set by self:_set_hedra()
        self.hedron1 = None
        self.hedron2 = None

        self.h1key = None
        self.h2key = None

        self.id3 = tuple(self.aks[0:3])
        self.id32 = tuple(self.aks[1:4])

        # 4 matrices specifying hedron space coordinates of constituent atoms,
        # in this space atom 3 is on on +Z axis
        # see coord_space()
        self.initial_coords = None
        self.a4_pre_rotation = None

        # PIC_Residue object which includes this dihedron;
        # set by Residue:linkDihedra()
        self.pic_residue = None
        # order of atoms in dihedron is reversed from order of atoms in hedra
        self.reverse = False

        if 'dihedral1' in kwargs:
            self.dihedral1 = float(kwargs['dihedral1'])
            # self.init_pos()  # can't do here because need adjacent residues
        else:
            self.dihedral1 = None

        # print(self, self.dclass)

    def __str__(self):
        """Print string for Dihedron object."""
        return ('4-' + str(self.id) + ' ' + self.rdh_class + ' ' +
                str(self.dihedral1) + ' (' + str(self.pic_residue) + ')')

    @staticmethod
    def _get_hedron(pic_res, id3):
        """Find specified hedron on this residue or its adjacent neighbors."""
        hedron = pic_res.hedra.get(id3, None)
        if (not hedron and pic_res.rprev):
            hedron = pic_res.rprev.hedra.get(id3, None)
        if (not hedron and pic_res.rnext):
            hedron = pic_res.rnext.hedra.get(id3, None)
        return hedron

    def _set_hedra(self):
        """Work out hedra keys and set rev flag."""
        rev = False
        res = self.pic_residue
        h1key = self.id3
        hedron1 = Dihedron._get_hedron(res, h1key)
        if not hedron1:
            rev = True
            h1key = tuple(self.aks[2::-1])
            hedron1 = Dihedron._get_hedron(res, h1key)
            h2key = tuple(self.aks[3:0:-1])
        else:
            h2key = self.id32

        if not hedron1:
            raise HedronMatchError(
                "can't find 1st hedron for key %s dihedron %s" % (h1key, self))

        hedron2 = Dihedron._get_hedron(res, h2key)

        if not hedron2:
            # print(res.hedra)
            # print(res.rnext.hedra)
            raise HedronMatchError(
                "can't find 2nd hedron for key %s dihedron %s" % (h2key, self))

        self.hedron1 = hedron1
        self.h1key = h1key
        self.hedron2 = hedron2
        self.h2key = h2key

        self.reverse = rev

        return rev, hedron1, hedron2

    def init_pos(self, updating=False):
        """Set hedron-space atom coords with dihedral1 applied.

        :param updating: bool
            skip _set_hedra if True
        """
        hedron1 = self.hedron1
        if updating and hedron1 is not None:
            rev = self.reverse
            hedron2 = self.hedron2
        else:
            rev, hedron1, hedron2 = self._set_hedra()

        acount = 0
        for a in hedron1.atoms:
            if a is not None:
                acount += 1
        for a in hedron2.atoms:
            if a is not None:
                acount += 1
        if 6 > acount:
            raise MissingAtomError('dihedron: hedra missing atoms: ' + self)

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
        # initial.append(mrz @ a4_pre_rotation)
        initial.append(mrz.dot(a4_pre_rotation))

        self.initial_coords = tuple(initial)
        self.a4_pre_rotation = a4_pre_rotation

        self.atoms_updated = True

    """
    # unused - get Biopython Residue atoms to test dihedral calculations
    def find_bp_atom(self, ak):
        bpa = self.pic_residue.bp_atoms.get(ak, None)
        if bpa is not None:
            return bpa
        if self.pic_residue.rnext:
            bpa = self.pic_residue.rnext.bp_atoms.get(ak, None)
            if bpa is not None:
                return bpa
        if self.pic_residue.rprev:
            bpa = self.pic_residue.rprev.bp_atoms.get(ak, None)
        return bpa
    """

    def set_dihedral(self, dangle_deg):
        """Save new dihedral angle and update initial_coords.

        :param dangle_deg: float
            New dihedral angle in degrees
        """
        self.dihedral1 = dangle_deg
        # if self.atoms_updated:
        if self.hedron1 is not None:  # then we are updating
            initial = list(self.initial_coords)
            mrz = gen_Mrot(numpy.deg2rad(dangle_deg, 'z'))
            # initial[3] = mrz @ self.a4_pre_rotation
            initial[3] = mrz.dot(self.a4_pre_rotation)
            self.initial_coords = tuple(initial)
            # self.atoms_updated = True still valid
        else:
            self.init_pos()

    @staticmethod
    def _get_dadad(acs):
        """Get distance, angle, distance, angle, distance for 4 atoms.

        :param acs: list[4] of numpy [4][1] array
            Atom coordinates
        """
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
        """Compute residue dihedral, bond angles, bond lengths.

        Source data is Biopython Residue.Atom coords.
        Call link_dihedra before this so can find Residue.Atom coords.
        Updates hedron and dihedron values, then all local atom coords
        for both hedra and this dihedron.
        """
        rev, hed1, hed2 = self._set_hedra()

        atom_coords = self.pic_residue.atom_coords
        acs = self.gen_acs(atom_coords)
        mt = coord_space(acs[:3])
        # do4 = mt @ acs[3]
        do4 = mt.dot(acs[3])

        dh1r = numpy.rad2deg(numpy.arctan2(do4[1][0], do4[0][0]))

        """
        # compare rtm coordspace dihedral
        # vs Vector() dihedral with self.bp_atoms

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
        """

        self.dihedral1 = dh1r

        a0a1, a0a1a2, a1a2, a1a2a3, a2a3 = Dihedron._get_dadad(acs)

        if not rev:
            hed1.len1 = set_accuracy_95(a0a1)
            hed1.len3 = hed2.len1 = set_accuracy_95(a1a2)
            hed2.len3 = set_accuracy_95(a2a3)
        else:
            hed1.len3 = set_accuracy_95(a0a1)
            hed1.len1 = hed2.len3 = set_accuracy_95(a1a2)
            hed2.len1 = set_accuracy_95(a2a3)

        hed1.angle2 = set_accuracy_95(a0a1a2)
        hed2.angle2 = set_accuracy_95(a1a2a3)

        # hed1.atoms_updated = False
        # hed2.atoms_updated = False

        hed1.init_pos()
        hed2.init_pos()

        # self.atoms_updated = False

        self.init_pos(True)

        # print(self)
        # print(hed1)
        # print(hed2)


class Hedron(Edron_Base):
    """Class to represent three joined atoms forming a plane.

    Contains atom coordinates in local coordinate space, central atom
    at origin.  Stored in two orientations, with the 3rd (forward) or
    first (reversed) atom on the +Z axis.

    Attributes
    ----------
    len1 : float
        distance between 1st and 2nd atom
    angle2 : float
        angle (degrees) formed by 3 atoms
    len3 : float
        distance between 2nd and 3rd atoms
    atoms : tuple[3] of numpy arrays [4][1]
        3 atoms comprising hedron, 1st on XZ, 2nd at origin, 3rd on +Z
    atomsR : tuple[3] of numpy arrays [4][1]
        atoms reversed, 1st on +Z, 2nd at origin, 3rd on XZ plane

    Methods
    -------
    init_pos()
        Create hedron space atom coordinate numpy arrays.
    hedron_from_atoms()
        Compute length, angle, length for hedron from PIC_Residue atom coords

    """

    def __init__(self, *args, **kwargs):
        """Initialize Hedron with sequence of AtomKeys, kwargs.

        Acceptable input:
            As for Edron_Base, plus optional 'len1', 'angle2', 'len3'
            keyworded values.
        """
        super().__init__(*args, **kwargs)

        # print('initialising', self.id)

        # 3 matrices specifying hedron space coordinates of constituent atoms,
        # initially atom3 on +Z axis
        self.atoms = None
        # 3 matrices, hedron space coordinates, reversed order
        # initially atom1 on +Z axis
        self.atomsR = None

        if 'len1' in kwargs:
            # distance between 1st and 2nd atom
            self.len1 = float(kwargs['len1'])
            # angle formed between 3 atoms
            self.angle2 = float(kwargs['angle2'])
            # distance between 2nd and 3rd atoms
            self.len3 = float(kwargs['len3'])

            self.init_pos()
        else:
            self.len1 = None
            self.angle2 = None
            self.len3 = None

        # print(self)

    def __str__(self):
        """Print string for Hedron object."""
        return ('3-' + self.id + ' ' + self.rdh_class + ' ' + str(self.len1)
                + ' ' + str(self.angle2) + ' ' + str(self.len3))

    def init_pos(self):
        """Initialize Hedron by creating atom coordinate numpy arrays."""
        # build hedron with a2 on +Z axis, a1 at origin,
        # a0 in -Z at angle n XZ plane
        atoms = []
        for _ in range(3):
            # note this initializes a1 to 0,0,0 origin
            atoms.append(numpy.array([[0], [0], [0], [1]],
                                     dtype=numpy.float64))  # 4x1 array

        # supplementary angle radian: angles which add to 180 are supplementary
        sar = numpy.deg2rad(180.0 - self.angle2)

        # a2 is len3 up from a2 on Z axis, X=Y=0
        atoms[2][2][0] = self.len3
        # a0 X is sin( sar ) * len1
        atoms[0][0][0] = numpy.sin(sar) * self.len1
        # a0 Z is -(cos( sar ) * len1)
        # (assume angle2 always obtuse, so a0 is in -Z)
        atoms[0][2][0] = - (numpy.cos(sar) * self.len1)

        self.atoms = tuple(atoms)

        atomsR = []
        # same again but 'reversed' : a0 on Z axis, a1 at origin, a2 in -Z
        for _ in range(3):
            # atom[1] to 0,0,0 origin
            atomsR.append(numpy.array([[0], [0], [0], [1]],
                                      dtype=numpy.float64))

        # a0r is len1 up from a1 on Z axis, X=Y=0
        atomsR[0][2][0] = self.len1
        # a2r X is sin( sar ) * len3
        atomsR[2][0][0] = numpy.sin(sar) * self.len3
        # a2r Z is -(cos( sar ) * len3)
        atomsR[2][2][0] = - (numpy.cos(sar) * self.len3)

        self.atomsR = tuple(atomsR)

        self.atoms_updated = True

    @staticmethod
    def _get_dad(acs):
        """Get distance, angle, distance for 3 atoms.

        :param acs: list[3] of numpy arrays [4][[1]]
        """
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
        """Compute length, angle, length for hedron for residue atom coords."""
        acs = self.gen_acs(atom_coords)

        len1, angle2, len3 = Hedron._get_dad(acs)
        self.len1 = set_accuracy_95(len1)
        self.angle2 = set_accuracy_95(angle2)
        self.len3 = set_accuracy_95(len3)

        # self.atoms_updated = False
        self.init_pos()


class PIC_Residue(object):
    """Class to extend Biopython Residue with internal coordinate data.

    Attributes
    ----------
    residue : Biopython Residue object reference
        The Residue object this extends
    hedra : dict indexed by 3-tuples of AtomKeys
        Hedra forming this residue
    dihedra : dict indexed by 4-tuples of AtomKeys
        Dihedra forming (overlapping) this residue
    rprev, rnext : PIC_Residue objects
        References to adjacent (bonded, not missing) residues in chain
    atom_coords : AtomKey indexed dict of numpy [4][1] arrays
        Local copy of atom homogeneous coordinates [4][1] for work
        distinct from Bopython Residue/Atom
    alt_ids : list of char
        AltLoc IDs from PDB file
    bfactors : dict
        AtomKey indexed B-factors as read from PDB file

    Methods
    -------
    atm241(coord) :
        Convert 1x3 cartesian coords to 4x1 homogeneous coords
    load_PIC(edron) :
        Process parsed (di-/h-)edron data from PIC file
    link_dihedra() :
        Link dihedra to this residue, form id3_dh_index
    render_hedra() :
        Call init_pos for each hedron in hedra
    render_dihedra() :
        Call init_pos for each dihedron in dihedra
    assemble(atomCoordsIn, transforms) :
        Compute atom coordinates for this residue from internal coordinates
    dihedra_from_atoms() :
        Create hedra and dihedra for atom coordinates
    write_PIC(pdbid, chainId, s, stats) :
        Generate PIC format strings for this residue
    coords_to_residue() :
        Convert homogeneous atom_coords to Biopython cartesian Atom coords

    """

    def __init__(self, parent, NO_ALTLOC=False):
        """Initialize PIC_Residue with parent Biopython Residue.

        :param parent: Biopython Residue object
            The Biopython Residue this object extends
        :param NO_ALTLOC: bool default False
            Option to disable processing altloc disordered atoms, use selected.
        """
        # NO_ALTLOC=True will turn off alotloc positions and just use selected
        self.residue = parent
        # self.ndx = ndx
        # dict of hedron objects indexed by hedron keys
        self.hedra = {}
        # dict of dihedron objects indexed by dihedron keys
        self.dihedra = {}
        # map of dihedron key (first 3 atom keys) to dihedron
        # for all dihedra in Residue
        # set by Residue.link_dihedra()
        self.id3_dh_index = {}
        # set of AtomKeys involved in dihedra, used by split_akl
        self.ak_set = set()
        # reference to adjacent residues in chain
        self.rprev = None
        self.rnext = None
        # local copy, homogeneous coordinates for atoms, numpy [4][1]
        # generated from dihedra include some i+1 atoms
        # or initialised here from parent residue if loaded from coordinates
        self.atom_coords = {}
        # TODO: remove bp_atoms, just for testing against bp vector()
        # self.bp_atoms = {}
        # bfactors copied from PDB file
        self.bfactors = {}
        if NO_ALTLOC:
            self.alt_ids = None
        else:
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
            self.NCaCKey = (AtomKey(self, 'N'),
                            AtomKey(self, 'CA'),
                            AtomKey(self, 'C'))

            for atom in parent.get_atoms():
                if hasattr(atom, 'child_dict'):
                    if NO_ALTLOC:
                        self._add_atom(atom.selected_child)
                    else:
                        for atm in atom.child_dict.values():
                            self._add_atom(atm)
                else:
                    self._add_atom(atom)

            # print(self.atom_coords)

    accept_atoms = ('N', 'CA', 'C', 'O', 'H', 'CB', 'CG', 'CG1', 'OG1', 'SG',
                    'CG2', 'CD', 'CD1', 'SD', 'OD1', 'ND1', 'CD2', 'ND2', 'CE',
                    'CE1', 'NE', 'OE1', 'NE1', 'CE2', 'OE2', 'NE2', 'CE3',
                    'CZ', 'NZ', 'CZ2', 'CZ3')

    @staticmethod
    def atm241(coord):
        """Convert 1x3 cartesian coordinates to 4x1 homogeneous coordinates."""
        arr41 = numpy.append(coord, [1])
        return numpy.array(arr41, dtype=numpy.float64)[
            numpy.newaxis].transpose()

    def _add_atom(self, atm):
        """Filter Biopython Atom with accept_atoms; set atom_coords, ak_set."""
        if atm.name not in self.accept_atoms:
            # print('skip:', atm.name)
            return
        ak = AtomKey(self, atm)
        self.atom_coords[ak] = PIC_Residue.atm241(atm.coord)
        self.ak_set.add(ak)
        # self.bp_atoms[ak] = atm

    def __str__(self):
        """Print string is parent Residue ID."""
        return('(' + str(self.residue.id) + ')')

    def load_PIC(self, edron):
        """Process parsed (di-/h-)edron data from PIC file.

        :param edron: parse dictionary from Edron_Base.edron_re
        """
        ek = [edron['a1'], edron['a2'], edron['a3']]

        if edron['a4'] is not None:
            ek.append(edron['a4'])
            self.dihedra[tuple(ek)] = Dihedron(ek, **edron)
        else:
            # self.hedra[Edron_Base.gen_key(ek)] = Hedron(ek, **edron)
            self.hedra[tuple(ek)] = Hedron(ek, **edron)

    def link_dihedra(self):
        """Link dihedra to this residue, form id3_dh_index and ak_set."""
        id3i = {}
        for dh in self.dihedra.values():
            dh.pic_residue = self  # each dihedron can find its pic_residue
            id3 = dh.id3
            if id3 not in id3i:
                id3i[id3] = []
            id3i[id3].append(dh)
            self.ak_set.update(id3)
        # map to find each dihedron from atom tokens 1-3
        self.id3_dh_index = id3i

    def render_dihedra(self):
        # TODO: check!
        """Call init_pos for each dihedron in dihedra."""
        for d in self.dihedra.values():
            # if d.atoms_updated:
            d.init_pos()

    def assemble(self, atomCoordsIn, transforms=False):
        """Compute atom coordinates for this residue from internal coordinates.

        Join dihedrons from N-CA-C and N-CA-CB hedrons, computing protein
        space coordinates for backbone and sidechain atoms

        Algorithm
        ---------
        form double-ended queue, start with n-ca-c, o-c-ca, n-ca-cb
            [ o-c-ca not 2nd hedron for any dihedron and thus won't be picked
            up w/o adding here ]
        generate triple keys for current residue
        if no atomCoordsIn, use initial coords from generating dihedral for
            n-ca-c initial positions (dihedron coordinate space)

        while queue not empty
            get triple key
            for each dihedral starting with triple key (1st hedron)
                if have coordinates for all 4 atoms already
                    add 2nd hedron key to back of queue
                else if have coordinates for 1st 3 atoms
                    compute forward and reverse transform to take 1st 3 atoms
                        to/from dihedron initial coordinate space
                    use reverse transform to get position of 4th atom in
                        current coordinates from dihedron initial coordinates
                    add 2nd hedron key to back of queue
                else
                    ordering failed, put triple key at back of queue and hope
                        next time we have 1st 3 atom positions (should not
                        happen)

        loop terminates (queue drains) as triple keys which do not start any
            dihedra are removed without action

        :param atomCoordsIn: list[3] of numpy arrays (optional)
            N, Ca, C coordinates of previous residue to build this one on to
        :param transforms: bool default False
            Option to return transformation matrices for each hedron instead
            of coordinates.
        :return: atomCoords for residue in protein space relative to
            acomCoordsIn OR table of transformation matrices according to
            transforms parameter
        """
        dbg = False

        transformations = {}
        atomCoords = atomCoordsIn
        NCaCKey = self.NCaCKey

        if transforms:
            transformations[NCaCKey] = numpy.identity(4, dtype=numpy.float64)

        startLst = []
        for lst in [
            (AtomKey(self, 'C'), AtomKey(self, 'CA'), AtomKey(self, 'N')),
            (AtomKey(self, 'N'), AtomKey(self, 'CA'), AtomKey(self, 'CB')),
            (AtomKey(self, 'O'), AtomKey(self, 'C'), AtomKey(self, 'CA')),
            self.NCaCKey
        ]:
            startLst.extend(self._split_akl(lst))

        q = deque(startLst)

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
                            if transforms:
                                transformations[h1k] = mtr
                            if dbg:
                                print(
                                    '        initial_coords[3]=',
                                    d.initial_coords[3].transpose())
                            # acak3 = mtr @ d.initial_coords[3]
                            acak3 = mtr.dot(d.initial_coords[3])
                            if dbg:
                                print('        acak3=', acak3.transpose())
                            for i in range(3):
                                acak3[i][0] = set_accuracy_83(acak3[i][0])
                            atomCoords[akl[3]] = acak3
                            if dbg:
                                print(
                                    '        3- finished, ak:', akl[3],
                                    'coords:', atomCoords[akl[3]].transpose())
                            q.appendleft(d_h2key)
                        else:
                            print('no coords to start', d)
                            pass
                    else:
                        print('no initial coords for', d)
                        pass
        # print('coord_space returning')
        if transforms:
            return transformations
        else:
            return atomCoords

    altloc_re = re.compile(r'-([A-Z])\Z')

    def _split_akl(self, lst):
        """Get AtomKeys for this residue ak_set given generic list of AtomKeys.

        Given a list of AtomKeys (aks) for a Hedron or Dihedron,
          return:
               list of matching aks that have id3_dh in this residue
                    (ak may change if occupancy != 1.00)
            or
               multiple lists of matching aks expanded for all atom altlocs
            or
               empty list if any of atom_coord(ak) missing

        :param lst: list[3] or [4] of AtomKeys
            non-altloc AtomKeys to match to specific AtomKeys for this residue
        """
        altloc_ndx = AtomKey.fields.altloc

        # step 1
        # given a list of AtomKeys (aks)
        #  form a new list of same aks with coords or diheds in this residue
        #      plus lists of matching altloc aks in coords or diheds
        edraLst = []
        altlocs = set()
        for ak in lst:
            if ak in self.ak_set:
                edraLst.append(tuple([ak]))
            else:
                ak2_lst = []
                for ak2 in self.ak_set:
                    if ak.altloc_match(ak2):
                        # print(key)
                        ak2_lst.append(ak2)
                        altloc = ak2.akl[altloc_ndx]
                        if altloc is not None:
                            altlocs.add(altloc)
                edraLst.append(tuple(ak2_lst))

        # step 2
        # check and finish for
        #   missing atoms
        #   simple case no altlocs
        # else form new AtomKey lists covering all altloc permutations
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
            return [tuple(newAKL)]
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
                new_edraLst.append(tuple(alhl))
            # print(new_edraLst)
            return new_edraLst

        # can't reach here

    def _gen_edra(self, lst):
        """Populate hedra/dihedra given edron ID tuple.

        Given list of AtomKeys defining hedron or dihedron
          convert to AtomKeys with coordinates in this residue
          add appropriately to self.di/hedra, expand as needed atom altlocs

        :param lst: tuple of AtomKeys
            Specifies Hedron or Dihedron
        """
        if 4 > len(lst):
            dct, obj = self.hedra, Hedron
        else:
            dct, obj = self.dihedra, Dihedron

        if all(ak in self.atom_coords for ak in lst):
            hl = [lst]
        else:
            hl = self._split_akl(lst)

        for nlst in hl:
            # dct[Edron_Base.gen_key(nlst)] = obj(**zdh(nlst))
            #
            dct[tuple(nlst)] = obj(nlst)

    def dihedra_from_atoms(self):
        """Create hedra and dihedra for atom coordinates."""
        AK = AtomKey
        S = self

        sN, sCA, sC = AK(S, 'N'), AK(S, 'CA'), AK(S, 'C')
        sO, sCB = AK(S, 'O'), AK(S, 'CB')
        # sN = AtomKey(self, 'N')

        if self.rnext:
            # atom_coords, hedra and dihedra for backbone dihedra
            # which reach into next residue
            rn = self.rnext
            nN, nCA, nC = AK(rn, 'N'), AK(rn, 'CA'), AK(rn, 'C')

            for ak in (nN, nCA, nC):
                if ak in rn.atom_coords:
                    self.atom_coords[ak] = rn.atom_coords[ak]
                    self.ak_set.add(ak)
                else:
                    for rn_ak in rn.atom_coords.keys():
                        if rn_ak.altloc_match(ak):
                            self.atom_coords[rn_ak] = rn.atom_coords[rn_ak]
                            self.ak_set.add(rn_ak)

            self._gen_edra([sN, sCA, sC, nN])  # psi
            self._gen_edra([sCA, sC, nN, nCA])  # omega i+1
            self._gen_edra([sC, nN, nCA, nC])  # phi i+1
            self._gen_edra([sCA, sC, nN])
            self._gen_edra([sC, nN, nCA])
            self._gen_edra([nN, nCA, nC])

        # backbone O and C-beta hedra and dihedra within this residue
        self._gen_edra([sN, sCA, sC, sO])
        self._gen_edra([sO, sC, sCA, sCB])

        if not self.rprev:
            self._gen_edra([sN, sCA, sC])

        self._gen_edra([sCA, sC, sO])
        self._gen_edra([sCB, sCA, sC])

        # if res is not G or A
        # only needed for sidechain CG residues
        # (not gly or ala or any missing rest of side chain)

        if 'G' != self.lc and 'A' != self.lc:
            for ak in self.atom_coords:
                if ak.is_sidechain_gamma():
                    self._gen_edra([sN, sCA, sCB])
                    break

        # amide proton N H if present
        sH = AK(S, 'H')
        self._gen_edra([sH, sN, sCA])
        self._gen_edra([sC, sCA, sN, sH])

        # sidechain hedra and dihedra

        sidechain = pic_data_sidechains.get(self.lc, [])
        for edra in sidechain:
            r_edra = [AK(S, atom) for atom in edra]
            self._gen_edra(r_edra[0:4])  # [4] is label on some table entries

        if ('G' == self.lc and sCB not in self.atom_coords):
            # add C-beta for Gly
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

    def _pdb_atom_string(self, atm):
        """Generate PDB ATOM record.

        :param atm: Biopython Atom object reference
        """
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
    def _residue_string(res):
        """Generate PIC Residue string.

        Enough to create Biopython Residue object without actual Atoms.

        :param res: Biopython Residue object reference
        """
        return (str(res.get_full_id())
                + ' ' + res.resname
                + ' [' + res.get_segid() + ']'
                # + (' ' + str(ndx) if ndx is not None else '')
                + '\n')

    def _write_pic_bfac(self, atm, s, col):
        ak = AtomKey(self, atm)
        if (0 == col % 5):
            s += 'BFAC:'
        s += ' ' + ak.id + ' ' + "{:6.2f}".format(atm.get_bfactor())
        col += 1
        if (0 == col % 5):
            s += '\n'
        return s, col

    def write_PIC(self, pdbid, chainid, s='', stats=False):
        """Write PIC format lines for this residue.

        :param str pdbid: PDB idcode string
        :param str chainid: PDB Chain ID character
        :param str s: result string to add to
        :param bool stats: option to output counts TODO: what stats???
        """
        s += PIC_Residue._residue_string(self.residue)
        if not self.rprev:
            try:
                ts = self._pdb_atom_string(self.residue['N'])
                ts += self._pdb_atom_string(self.residue['CA'])
                ts += self._pdb_atom_string(self.residue['C'])
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

        col = 0
        for a in self.residue.get_atoms():
            if hasattr(a, 'child_dict'):
                if self.alt_ids is None:
                    s, col = self._write_pic_bfac(a.selected_child, s, col)
                else:
                    for atm in a.child_dict.values():
                        s, col = self._write_pic_bfac(atm, s, col)
            else:
                s, col = self._write_pic_bfac(a, s, col)
        if (0 != col % 5):
            s += '\n'

        return s

    def coords_to_residue(self):
        """Convert homogeneous PIC atom coords to cartesian Atom coords."""
        seqpos, icode = self.residue.id[1:3]
        seqpos = str(seqpos)
        spNdx, icNdx, resnNdx, atmNdx, occNdx, altlocNdx = AtomKey.fields

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
                    # print('new', ak)
                    occ = akl[occNdx]
                    aloc = akl[altlocNdx]
                    bfac = self.bfactors.get(ak.id, None)
                    newAtom = Atom(atm, atm_coords,
                                   (0.0 if bfac is None else bfac),
                                   (1.00 if occ is None else float(occ)),
                                   (' ' if aloc is None else aloc),
                                   atm, ndx, atm[0])
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
                    # print('update', ak)
                    if (2 == Atm.is_disordered()
                            and Atm.disordered_has_id(altloc)):
                        Atm.disordered_select(altloc)
                    Atm.set_coord(atm_coords)
                    ndx = Atm.get_serial_number()

        Res.parent.pic.ndx = ndx

        # print(ak, ac)


class PIC_Chain:
    """Class to extend Biopython Chain with internal coordinate data.

    Attributes
    ----------
    chain : biopython Chain object reference
        The Chain object this extends
    ordered_aa_pic_list : list of PIC_Residue objects
        PIC_Residue objects PIC algorithms can process (currently no hetatoms)
    initNCaC : AtomKey indexed dictionary of N, Ca, C atom coordinates to start
        chain segments (first residue or after chain break)

    Methods
    -------
    set_residues()
        Add .pic attribute for all Residues, populate ordered_aa_pic_list, set
        PIC_Residue rprev, rnext or initNCaC coordinates
    link_residues()
        Call link_dihedra() on each PIC_Residue (needs rprev, rnext set)
    render_dihedra()
        Call render_hedra() and render_dihedra() on each PIC_Residue
    assemble_residues()
        Generate PIC_Residue atom coords from internal coordinates
    coords_to_structure()
        update Biopython Residue.Atom coords for all Residues with PIC_Residue
        attributes
    internal_to_atom_coordinates()
        Process pic data to Residue/Atom coordinates
    dihedra_from_atoms()
        Calculate dihedrals, angles, bond lengths for Atom data

    """

    def __init__(self, parent):
        """Initialize PIC_Chain object, with or without residue/Atom data.

        :param parent: Biopython Chain object
            Chain object this extends
        """
        self.chain = parent
        self.ordered_aa_pic_list = []
        self.initNCaC = {}
        self.set_residues()  # no effect if no residues loaded

    def set_residues(self):
        """Initialize pic data for loaded Residues.

        Add PIC_Residue as .pic attribute for each Residue in parent Chain;
        populate ordered_aa_pic_list with PIC_Residue references for residues
        which can be built (currently no hetatms); set rprev and rnext on each
        sequential PIC_Residue, populate initNCaC at start and after chain
        breaks.
        """
        # ndx = 0
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
                    for atm in ('N', 'CA', 'C'):
                        initNCaC[AtomKey(rpic, atm)] = PIC_Residue.atm241(
                            res[atm].coord)
                    self.initNCaC[respos] = initNCaC
                last_ord_res = res
            else:
                # print('skipping res ' + str(res.id) + ' ' + res.resname
                #       + (' disordered'
                #           if res.disordered else ' not disordered'))
                pass
            last_res = res

    def link_residues(self):
        """link_dihedra() for each PIC_Residue; needs rprev, rnext set."""
        for rpic in self.ordered_aa_pic_list:
            rpic.link_dihedra()

    def render_dihedra(self):
        """Set dihedron local coords for each PIC_Residue."""
        for rpic in self.ordered_aa_pic_list:
            rpic.render_dihedra()

    def assemble_residues(self, start=False, fin=False):
        """Generate PIC_Residue atom coords from internal coordinates.

        Filter positions between start and fin if set, find appropriate start
        coordinates for each residue and pass to PIC_Residue.assemble()

        :param start, fin: lists
             sequence position, insert code for begin, end of subregion to
             process
        """
        for rpic in self.ordered_aa_pic_list:
            respos, resicode = rpic.residue.id[1:]
            go = True
            if (start and
                (start[0] > respos or
                    (start[0] == respos and start[1] < resicode))):
                go = False
            if (go and fin and
                (fin[0] < respos or
                    (fin[0] == respos and fin[1] > resicode))):
                go = False
            if go:
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
                rpic.atom_coords = rpic.assemble(startPos)
                rpic.ak_set = rpic.atom_coords.keys()

    def coords_to_structure(self):
        """All pic atom_coords to Biopython Residue/Atom coords."""
        self.ndx = 0
        for res in self.chain.get_residues():
            if hasattr(res, 'pic'):
                res.pic.coords_to_residue()

    # TODO: fix to allow only updating parts of structure
    # esp sidechain rotamers

    def internal_to_atom_coordinates(self):
        """Complete process pic data to Residue/Atom coords."""
        # self.link_residues()
        # self.render_dihedra()
        self.assemble_residues()
        self.coords_to_structure()

    def dihedra_from_atoms(self):
        """Calculate dihedrals, angles, bond lengths for Atom data."""
        for rpic in self.ordered_aa_pic_list:
            rpic.dihedra_from_atoms()


# PIC  Exceptions
class HedronMatchError(Exception):
    """Cannot find hedron in residue for given key."""

    pass


class MissingAtomError(Exception):
    """Missing atom coordinates for hedron or dihedron."""

    pass
