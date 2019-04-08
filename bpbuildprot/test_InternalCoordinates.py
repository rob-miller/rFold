#!/usr/local/bin/python3
"""Unit tests for InternalCoordinates module.

Copyright 2019 by Robert T. Miller  All Rights Reserved.
"""

import unittest
from InternalCoordinates import AtomKey


# import warnings

"""
try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB.")
"""


# target = __import__("InternalCoordinates.py")


class test_InternalCoordinates(unittest.TestCase):
    """Test with strings for now."""

    def setUp(self):
        """Load test PDB file here."""
        pass

    def test_AtomKey(self):
        #  def test_AtomKey(self):
        """Test AtomKey inputs."""
        ak = AtomKey(**{'respos': 52, 'resname': 'A', 'atm': 'CA'})
        self.assertEqual(ak.id, '52_A_CA', msg='kwargs test')
        ak2 = AtomKey({'respos': 52, 'resname': 'A', 'atm': 'CB'})
        self.assertEqual(ak2.id, '52_A_CB', msg='dict test')
        """Test AtomKey sort."""
        self.assertGreater(ak2.id, ak.id, msg='gt sort test')
        # self.assertEqual(1, 1)

    def test_Hstruct(self):
        

# if __name__ == '__main__':
#    runner = unittest.TextTestRunner(verbosity=2)
#    unittest.main(testRunner=runner)


if __name__ == '__main__':
    unittest.main()
