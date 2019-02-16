//
// OpenSCAD array indices to reference protein data
//

// protein base level
p_pdbid = 0;
p_proteinScale = 1;
p_chainData = 2;

// chain level data
c_chainID = 0;
c_hedra = 1;
c_dihedra = 2;
c_residues = 3;

// hedra definitions
h_len1 = 0;
h_angle2 = 1;
h_len3 = 2;
h_atom1 = 3;
h_atom2 = 4;
h_atom3 = 5;

// dihedra specifications for each residue in sequence, dihedral array
d_dangle1 = 0;
d_h1ndx = 1;
d_h2ndx = 2;
d_reversed = 3;
d_dihedralTransform = 4;

// residueSet: world transform for each residue in sequence array
r_resNdx = 0;
r_resID = 1;
r_resTransform = 2;

