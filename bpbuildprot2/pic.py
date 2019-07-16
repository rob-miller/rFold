# Copyright 2019 by Robert T. Miller.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""PIC: Protein Internal Coordinates.

Interconversion between PDB atom and internal coordinates.

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
# import sys
import re
import itertools

from Bio.File import as_handle
from Bio.PDB.Polypeptide import three_to_one
from Bio._py3k import StringIO
from Bio.PDB.PDBIO import PDBIO

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

from ic_data import pic_data_sidechains, pic_data_backbone
from ic_data import residue_atom_bond_state, pic_data_sidechain_extras

from Bio.PDB.PDBExceptions import PDBException

# __updated__ is specifically for the python-coding-tools Visual Studio Code
# extension, which updates the variable on each file save

__updated__ = '2019-05-22 19:27:03'
# print('ver: ' + __updated__)
# print(sys.version)




def write_PDB(entity, file, pdbid=None, chainid=None):
    """Write PDB file with HEADER and TITLE."""
    with as_handle(file, 'w') as fp:
        try:
            if 'S' == entity.level:
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
                io = PDBIO()
                io.set_structure(entity)
                io.save(fp)

            else:
                raise PDBException("level not 'S': "
                                   + str(entity.level))
        except KeyError:
            raise Exception(
                "write_PIC: argument is not a Biopython PDB Entity "
                + str(entity))




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
    rot = rotaxis2m(-math.pi * 120.0 / 180.0, c)
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
    for i in range(3):
        cb[i] = set_accuracy_83(cb[i])
    return cb





