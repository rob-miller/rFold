#!/usr/local/bin/python3
"""Unit tests for pic module (internal coordinates).

Copyright 2019 by Robert T. Miller  All Rights Reserved.
"""

import unittest
import warnings

from Bio.PDB.pic import AtomKey, atom_to_internal_coordinates
from Bio.PDB.pic import write_PIC, read_PIC, compare_residues
from Bio.PDB.pic import internal_to_atom_coordinates, PIC_Residue
from Bio.PDB.PDBParser import PDBParser

from io import StringIO

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB.")

# target = __import__("InternalCoordinates.py")


class test_pic(unittest.TestCase):
    """Test with strings and PDB files."""

    def setUp(self):
        """Load test PDB file here."""
        pass

    def test_AtomKey(self):
        """Test AtomKey inputs and sort."""
        ak = AtomKey(**{'resSeq': 52, 'resname': 'A', 'atm': 'CA'})
        self.assertEqual(ak.id, '52_A_CA', msg='kwargs test')
        ak2 = AtomKey({'resSeq': 52, 'resname': 'A', 'atm': 'CB'})
        self.assertEqual(ak2.id, '52_A_CB', msg='dict test')
        """Test AtomKey sort."""
        self.assertGreater(ak2.id, ak.id, msg='gt sort test')

    def compare_pdb_rebuild(self, pdb):
        sp = StringIO()
        PDB_parser = PDBParser(PERMISSIVE=True, QUIET=True)

        pdb_structure = PDB_parser.get_structure(
            '0PDB', pdb)
        atom_to_internal_coordinates(pdb_structure)
        write_PIC(pdb_structure, sp)
        sp.seek(0)
        pdb2 = read_PIC(sp)
        internal_to_atom_coordinates(pdb2)
        r = compare_residues(pdb_structure, pdb2, verbose=False)
        return r

    def test_noHstruct(self):
        # print(__file__)
        rslt = {'rCount': 835, 'rMatchCount': 835, 'aCount': 6267,
                'disAtmCount': 0, 'aCoordMatchCount': 6086,
                'aFullIdMatchCount': 6086}
        PIC_Residue.accept_atoms = PIC_Residue.accept_backbone
        testRslt = self.compare_pdb_rebuild(
            '/Users/rob/proj/rFold/bpbuildprot/PDB/2XHE.pdb')
        print(testRslt)
        self.assertEqual(rslt, testRslt, msg="no H 2XHE build")

    def test_struct(self):
        rslt = {'rCount': 130, 'rMatchCount': 130, 'aCount': 1855,
                'disAtmCount': 0, 'aCoordMatchCount': 1855,
                'aFullIdMatchCount': 1855}
        PIC_Residue.accept_atoms = (PIC_Residue.accept_backbone
                                    + PIC_Residue.accept_hydrogens)
        testRslt = self.compare_pdb_rebuild(
            '/Users/rob/proj/rFold/bpbuildprot/PDB/2BEG.pdb')
        self.assertEqual(rslt, testRslt, msg="full 2BEG build")


# if __name__ == '__main__':
#    runner = unittest.TextTestRunner(verbosity=2)
#    unittest.main(testRunner=runner)
if __name__ == '__main__':
    unittest.main()
