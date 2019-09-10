//protein_scale=10;
topFn = 0;
$fn=topFn;
tubes=true;
support=true;
layerHeight=0.15;
capMode = false;

difference() {
translate([0,0,(capMode ? 0 : bondRadius*protein_scale)]) {
rotate([-90,0,0]) 
chain(protein);
}

/*
//translate([6,-10,0]) cube([30, 38,20]);   // view edge ca 1
translate([-25,-18,0]) cube([85, 38,20]);  // block ca 1
translate([-25,0,0]) cube([28,30,20]);     // block ca 1
rotate([0,0,12]) 
translate([36,0,0]) cube([40,80,20]);

rotate([0,-20,0])
translate([0,30,8]) cube([20,20,20]);
*/
}

//
// peptide.scad
// Copyright(c) 2019 Robert T. Miller.  All rights reserved.
// This file is part of the Biopython distribution and governed by your
// choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
// Please see the LICENSE file that should have been included as part of this
// package.

// This is the support file to build an OpenSCAD (http://www.openscad.org/) model of a protein
// from internal coordinates.
//
// Lines similar to:
//    protein_scale=2;
//    $fn=20;
//    chain(protein);
//
// should be pre-pended to the beginning of this file, and data matrices should be appended
//    below to form a program ready to load into the OpenSCAD application.
//
//  protein_scale is the value supplied when generating the data for build units per PDB angstrom.
//    You may wish to modify it here to adjust the appearance of the model in terms of atom sphere
//    or bond cylinder diamether
//
//  $fn (fragment number) is an OpenSCAD parameter controlling the smoothness of the model surface.
//    Smaller values will render faster, but yield more 'blocky' models.
//
//  chain(protein); is the command to call the chain() function below to render the data in the
//    'protein' array at the end of this file.
//
//  This is intended to be a working example, you are encouraged to modify the OpenSCAD subroutines
//    below to generate a model to your liking.  For more information, start with
//    http://www.openscad.org/cheatsheet/index.html

// output parameters
atomScale=0.8;
defaultAtomRadius = 0.77;

bondRadius = (tubes ? defaultAtomRadius * atomScale : 0.4);

// for $fn=8 which works nice on fdm printer
oRot = 22.5;  // 45/2, rotate to make fn=8 spheres and cylinders flat on build plate
apmFac = cos(180/8);  // apothem factor - multiply by radius for center to octagon side distance
octSide = 2* tan(180/8);  // multiply by radius to get length of octagon side
// for values of topFn:
tfnRot = (topFn ? 90-(180/topFn) : 90-(180/30));

clearance=0.3;  //0.4;   // sliding clearance
pClearance=0.2;  // press-fit clearance (magnets for h-bonds)
shim=0.05;
//extrusionWidth=0.4;
//minThickness=2*extrusionWidth;
bondLenFac = 0.6; //0.7; // 0.66;  // fraction of bond length to extend from atom for each arm of hedron in join
hblen = 1.97;

magR=3/2;    // magnet radius
magL=5;      // magnet length

nozzleDiameter=0.4;
wall = 3*nozzleDiameter;
joinerStep = 1+clearance;

//
// Generate a sphere to represent an atom.
// Colour and size determined for the atom covalent radius specified by the parameter 'a' by lookup
//   in the atomData table below, then scaled by the supplied parameter 'scal'.
module atom(a,scal,clr=0)
{
     ad = atomData[search([a],atomData)[0]];
     color(ad[1]) {
          rotate([0,0,tfnRot]) sphere(r=((ad[2]*atomScale)*scal)+clr);
     }
}

function hFlip(h,rev) = 
                 //   yes reversed                                     :  not reversed
                 //    0    1     2     3     4     5     6      7     :     0     1     2     3    4     5      6      7
                 //  len1  len3  atom1 atom3  a1    a2   a1-a2  a2-a3      len1  len3  atom1 atom3   a1    a3  a1-a2  a2-a3
     (rev ? [ h[2], h[0], h[5], h[3], h[8], h[6], h[10], h[9] ] : [ h[0], h[2], h[3], h[5], h[6], h[8],  h[9], h[10] ]);
     // h[1] = angle2 for both cases


module joinUnit(cOuterLen, cOuterRad, cInnerLen, cInnerRad, male=false) {
/*          rotate([0,0,oRot]) {
               cylinder(h=cInnerLen,r=cInnerRad,center=false,$fn=8);
               cylinder(h=cOuterLen,r=cOuterRad,center=false,$fn=8);
          }
*/
  if (male) {
          //rotate([0,0,tfnRot]) 
          rotate([0,0,oRot]) {
               cylinder(h=cInnerLen, r=cInnerRad, center=false, $fn=8);
               cylinder(h=cOuterLen,r=cOuterRad,center=false, $fn=8);
          }
     } else {
          rotate([0,0,tfnRot]) {
               cylinder(h=cInnerLen,r=cInnerRad,center=false);
               cylinder(h=cOuterLen,r=cOuterRad,center=false);
          }
     }
}

module joiner(bondlen, scal, male=0, ver=0, supportSelect=0) {  // ver = differentiate arms for assembly
     lenfac = bondLenFac;
     jClr = clearance+0.05;
     
     cOuterRad = bondRadius*scal - (2*wall + (male ? jClr/2 : -jClr/2));
     cInnerRad = cOuterRad - joinerStep;  // m/f jClr already in cOuterRad;  - (male ? 0 : -0*jClr/2);

     hArmLen = (bondlen * lenfac);
     lenClr = 0.5*jClr;  // length clearance applied to male and female both, so effective clearance is 2x this value
     cOuterLen = hArmLen * lenfac + (ver ? 0.5 : - 0.5) - (wall+ (male ? lenClr*2 : -lenClr*2  ));

     joinerOffset = (hArmLen * (1 - lenfac)) + (male ? lenClr : -lenClr) - (ver ? 1 : 0);
     //joinerOffset = (hArmLen*0.4) + (male ? 0*jClr : -0*jClr) - (ver ? 1 : 0);

     i=supportSelect-1;
     oside = cOuterRad*octSide;
     wid = oside+2*wall+4*jClr+1;

     if (male) {
	  rotate([0,180,0])
	       translate([0,0,-(bondlen-joinerOffset)]) {
               difference() {
                    joinUnit(cOuterLen, cOuterRad, bondlen, cInnerRad, male=true);
                    if (supportSelect) {
                         rotate([0,0,i*180]) {
                              translate([0,(cOuterRad*apmFac)-0.5*layerHeight,cOuterLen/2]) {
                                   cube([oside+2*shim,layerHeight+shim,cOuterLen+2*shim],center=true);
                              }
                         }
                    }
               }
               if (supportSelect) {
                    rotate([0,0,i*180]) {
                         translate([0,(cOuterRad*apmFac)-0.5*layerHeight,cOuterLen/2]) {
                              for (j=[0:1]) {
                                   rotate([0,(j?60:-60),0])
                                        cube([wid,layerHeight,2*nozzleDiameter],center=true);
                              }
                         }
                    }
               }
          }
     } else {
          translate([0,0,joinerOffset]) {
               joinUnit(cOuterLen, cOuterRad, bondlen, cInnerRad);
               if (supportSelect) {  // extra gap top and bottom because filament sags
                    supHeight = max(5*layerHeight,2*(cOuterRad-cOuterRad*apmFac));  // double because center=true below
                    echo(supHeight, cOuterRad, cOuterRad*apmFac);
                    for(j=[0:1]) {
                         rotate([0,0,j*180]) {
                              translate([0,(cOuterRad*apmFac),cOuterLen/2]) {
                                   cube([oside+2*shim,supHeight+shim,cOuterLen+2*shim],center=true);
                              }
                         }
                    }
               }
          }
          
     }
}

module bond(bl, br, scal, key, atm, ver, supportSel=0) {

     br = br * (key == 4 ? 0.7 : 1);   // handle skinnyBond
     bl = (key == 2 ? bl * bondLenFac : bl);  // make female joiner shorter
     if (key == 3) { // male join
	  joiner(bl, scal, male = true, ver = ver, supportSelect=supportSel);
     } else {  // regular bond / skinny / h-bond / female join
          bhblen = bl +(hblen/2 * protein_scale);
          rotate([0,0,tfnRot]) {
               difference() {
                    union() {
                         cylinder(h=bl,r=br,center=false);
                         if (key == 5) {
                              rotate([0,0,oRot-tfnRot]) cylinder(h=bhblen-1,r=(magR + clearance +wall),center=false, $fn=8);
                         }
                    }
                    atom(atm,scal,-clearance);  // remove overlap with atom to clear area for female join
                    if (key == 5) {
                         translate([0,0,(bhblen-magL)-pClearance])
                              cylinder(h=magL+pClearance+shim, r=magR+pClearance, center=false, $fn=8);
                    }
               }
          }
     }
}

/*
// Generate a 'hedron', one plane of 3 points, consisting of 3 atoms joined by two bonds.
// In some cases the sequence of atoms in the h[] array is reversed, as detailed in the comments.

  bond codes
  1 regular
  2 female join
  3 male join
  4 skinny
  5 h-bond magnet
 */
module hedron(h,rev=0,scal,split=0, supportSel) {
     //   yes reversed                                     :  not reversed
     //    0    1     2     3     4     5     6      7     :     0     1     2     3    4     5      6      7
     //  len1  len3  atom1 atom3  a1    a2   a1-a2  a2-a3      len1  len3  atom1 atom3   a1    a3  a1-a2  a2-a3
     //newh = (rev ? [ h[2], h[0], h[5], h[3], h[8], h[6], h[10], h[9] ] : [ h[0], h[2], h[3], h[5], h[6], h[8],  h[9], h[10] ]);
     // h[1] = angle2 for both cases
     newh = hFlip(h, rev);

     bondRad = bondRadius * scal;
     difference() {
	  union(){
	       if (h[7]) {
		    // central atom at 0,0,0
		    atom(h[4],scal);
	       }

	       if (newh[5] && newh[7] != 2) {  // not female join
		    // comments for non-reversed case
		    // atom 3 is len3 up on +z
		    translate([0,0,newh[1]])
                         difference() {
                         atom(newh[3],scal * (newh[7] == 4 ? 0.7 : 1));
                         if (newh[7] == 5) {  // make room for hbond magnet through atom - this branch not used for backbone N,O
                              translate([0,0,scal*hblen/2-magL-pClearance])
                                   cylinder(h=magL+pClearance,r=magR+pClearance,$fn=8);
                         }
                    }
	       }
	  
	       if (newh[7]) {
		    // atom 2 - atom 3 bond from origin up +z distance len3
		    bond(newh[1], bondRad, scal, newh[7], h[4], ver=1, supportSel=supportSel);
	       }
	       rotate([0, h[1], 0]) {                        // rotate following elements by angle2 about Y
		    if (newh[6]) {
			 bond(newh[0], bondRad, scal, newh[6], h[4], ver=1, supportSel=supportSel);  // h[4] is center atom?
		    }
		    if (newh[4] && newh[6] != 2) {   // if draw atom 2 and atom1-atom2 not joiner
                         translate([0,0,newh[0]]) {
                              difference() {
                                   atom(newh[2],scal * (newh[6] == 4 ? 0.7 : 1));  // put atom1 sphere len1 away on Z
                                   if (newh[6] == 5) {  // make room for hbond magnet through atom
                                        translate([0,0,scal*hblen/2-magL-pClearance])
                                             cylinder(h=magL+pClearance,r=magR+pClearance,$fn=8);
                                   }
                              }
                         }
		    }
	       }
	  }

	  if (split) {
	       // top / bottom half cutter
               // 0 - off, 1 = one side, 2 = other side
               thick = bondRadius * scal;
               Zdim = newh[0];
               Xdim = newh[1];
               echo(Xdim, Zdim);
	       translate([-Xdim,(split == 2 ? -thick : 0),-Zdim]) {
		    cube([2*Xdim,thick,2*Zdim]);
	       }
	  }
          
	  if (newh[7] == 2) {  // female join
	       joiner(newh[1], scal, male=false, ver=1, supportSelect=supportSel);
	  }
          
	  if (newh[6] == 2) {  // female join
	       rotate([0, h[1], 0]) {                        // rotate following elements by angle2 about Y
		    joiner(newh[0], scal, male=false, ver=1, supportSelect=supportSel);
		    translate([0,0,newh[0]])
			 atom(newh[2],scal+0.5,clearance);  // put atom1 sphere len1 away on Z
	       }
	  }
          
	  if (newh[7] == 2 || newh[6] == 2) {  // female join both sides
	       translate([0,0,newh[1]]) atom(newh[3],scal+0.5,clearance);

	  }
     }
}

// hook to call custom routines for specific hedra.  chain is h[h_chain], sequence position is h[h_seqpos1]
module hedronDispatch(h,rev=0,scal) {
     if (capMode) {
          if (h[h_seqpos1] == 1) {
               if (h[h_class] == "NCAC") {
                    hedron(h, rev, scal, 1);
               } else if (h[h_class] == "CBCAC") {
                    color("yellow") {  // ca-cb
                         hedron(h, rev, scal);
                    }
               }
          }
          
     } else {

              if (h[h_seqpos1] == 0) {
              } else if (h[h_seqpos1] == 0) {
              } else if (h[h_seqpos1] == 0) {
     } else if (h[h_class] == "NCAC") {
          
          if (h[h_seqpos1] == 1) {
               hedron(h, rev, scal, 2, (support ? 1 : 0));
          } else if (h[h_seqpos1] == 5) {
               hedron(h, rev, scal, 1, (support ? 2 : 0));
          } else {
               hedron(h, rev, scal, 0, (support ? 2 : 0));
          }
     } else if (h[h_class] == "CBCAC") {
	  color("yellow") {  // ca-cb
               if (h[h_seqpos1] == 1 ) {
               } else if (h[h_seqpos1] == 5 ) {
               } else {
                    hedron(h, rev, scal);
               }
	  }
     } else if (h[h_class] == "HNCA") {
	  color("cyan") {  // h=n
	       hedron(h, rev, scal, 0, (support ? 1 : 0));
	  }
     } else if (h[h_class] == "CACN") {
	  color("darkgray") {  // ca-c-n, ca-c joint
	       hedron(h, rev, scal);
	  }
     } else if (h[h_class] == "CACO") {
           if (h[h_seqpos1] == 5) {
              } else {
	  color("red") {   // c=o
	       hedron(h, rev, scal,0, (support ? 1 : 0));
}

               
	  }
     } else {
          //echo("unrecognised hedra", h[h_class]);
           //color("green")
	  hedron(h, rev, scal, 0, (support ? 1 : 0));
     }
     }
}

// Generate a hedron rotated to specific angle d
module d2(d,hedra,scal)
{
     tz = (d[d_reversed] ? hedra[d[d_h2ndx]][2] : hedra[d[d_h2ndx]][0]);   // get h2 len1 depending on reversed
     rotate(d[d_dangle1]) {                                   // 4. rotate h2 to specified dihedral angle
          translate([0,0,tz]) {                       // 3. translate h2 h2:len1 up +z
               rotate([180, 0, 0]) {                  // 2. rotate h2r about X so h2:a3 in +z and h2:a1 in -z
                    hedronDispatch(hedra[d[d_h2ndx]],(!d[d_reversed]),scal);      // 1. reverse hedron 2 orientation = h2r
               }
          }
     }
}

// Generate two hedra at specified dihedral angle d
module dihedron(d,hedra,scal)
{
     if (d[d_h1new])
          hedronDispatch(hedra[d[d_h1ndx]],d[d_reversed],scal);                // reverse h1 if dihedral reversed
     if (d[d_h2new])
          d2(d,hedra,scal);
}

// Generate a residue consisting of the set of dihedra in the parameter 'r', referring to hedra the
//   table speicified in the parameter 'hedra'.
module residue(r,hedra, scal)
{
     for (d = r) {
          multmatrix(d[d_dihedralTransform]) {
               dihedron(d, hedra, scal);
          }
     }
}

// Generate a chain of residues, each positioned by a supplied rotation/translation matrix.
module chain(protein) 
{
    chnD = protein[p_chainData];
    c = chnD[c_residues];
    dihedra = chnD[c_dihedra];
    hedra = chnD[c_hedra];
     for (r = c) {
          multmatrix(r[r_resTransform]) {
               residue(dihedra[r[r_resNdx]],hedra, protein[p_proteinScale]);
          }
     }
}

//
// OpenSCAD array indices to reference protein data
//

// protein base level
p_pdbid = 0;
p_proteinScale = 1;
p_chainData = 2;

// chain level data
c_chainID = 0;
c_dihedra = 1;
c_hedra = 2;
c_residues = 3;

// hedra definitions
h_len1 = 0;
h_angle2 = 1;
h_len3 = 2;
h_atom1class = 3;
h_atom2class = 4;
h_atom3class = 5;
h_atom1state = 6;
h_atom2state = 7;
h_atom3state = 8;
h_bond1state = 9;
h_bond2state = 10;
h_chain = 11;
h_seqpos1 = 12;  // residue sequence position for first atom in hedra
h_class = 13;

// dihedra specifications for each residue in sequence, dihedral array
d_dangle1 = 0;
d_h1ndx = 1;
d_h2ndx = 2;
d_reversed = 3;
d_h1new = 4;
d_h2new = 5;
d_dihedralTransform = 6;

// residueSet: world transform for each residue in sequence array
r_resNdx = 0;
r_resID = 1;
r_resTransform = 2;


// covalent radii from Heyrovska, Raji : 'Atomic Structures of all the Twenty Essential Amino Acids and a Tripeptide, with Bond Lengths as Sums of Atomic Covalent Radii'
// https://arxiv.org/pdf/0804.2488.pdf

atomData = ( tubes ?
             [ ["Csb","green" , defaultAtomRadius], ["Cres","green" , defaultAtomRadius], ["Cdb","green" , defaultAtomRadius],
               ["Osb","red" , defaultAtomRadius], ["Ores","red" , defaultAtomRadius], ["Odb","red" , defaultAtomRadius],
               ["Nsb","blue" , defaultAtomRadius], ["Nres","blue" , defaultAtomRadius], ["Ndb","blue" , defaultAtomRadius],
               ["Hsb","gray" , defaultAtomRadius],
               ["Ssb","yellow" , defaultAtomRadius] ]
             :
                  
             [ ["Csb","green" , 0.77], ["Cres","green" , 0.72], ["Cdb","green" , 0.67],
               ["Osb","red" , 0.67], ["Ores","red" , 0.635], ["Odb","red" , 0.60],
               ["Nsb","blue" , 0.70], ["Nres","blue" , 0.66], ["Ndb","blue" , 0.62],
               ["Hsb","gray" , 0.37],
               ["Ssb","yellow" , 1.04] ]
     );


//include <../bpbuildprot2/1rtm-2a.PyPIC.scad>;
//include <../bpbuildprot2/1rtm.PyPIC.scad>;
include <1rtm.PyPIC.scad>;

// Protein specific array data below.
/*
   [  //hedra   a1 a2 a3 a1-a2 a2-a3
     [  14.73960, 108.16847,  15.42266, "Nsb", "Csb", "Cdb", 1, 1, 1, 2, 2, "NCAC", 1 ], // (4_A_N, 4_A_CA, 4_A_C)
     [  14.73960, 110.27884,  15.41567, "Nsb", "Csb", "Csb", 1, 0, 0, 1, 0, "NCACB", 2 ], // (4_A_N, 4_A_CA, 4_A_CB)
     [  15.42266, 117.75013,  12.32618, "Csb", "Cdb", "Odb", 0, 0, 1, 0, 5, "CACO", 3 ], // (4_A_CA, 4_A_C, 4_A_O)
     [  15.42266, 119.28867,  12.94834, "Csb", "Cdb", "Nsb", 0, 1, 1, 3, 1, "CACN", 4 ], // (4_A_CA, 4_A_C, 5_A_N)
     [  12.94834, 122.91801,  14.46867, "Cdb", "Nsb", "Csb", 0, 0, 0, 0, 0, "CNCA", 5 ], // (4_A_C, 5_A_N, 5_A_CA)
     [  15.41567, 114.93647,  15.42266, "Csb", "Csb", "Cdb", 1, 0, 0, 4, 0, "CBCAC", 6 ], // (4_A_CB, 4_A_CA, 4_A_C)
     [  10.00209, 116.50485,  14.73960, "Hsb", "Nsb", "Csb", 0, 0, 0, 0, 0, "HNCA", 7 ], // (4_A_H, 4_A_N, 4_A_CA)
     [  14.46867, 110.60090,  15.23981, "Nsb", "Csb", "Cdb", 1, 1, 1, 2, 2, "NCAC", 8 ], // (5_A_N, 5_A_CA, 5_A_C)
     [  14.46867, 114.66653,  15.15796, "Nsb", "Csb", "Csb", 0, 0, 0, 0, 0, "NCACB", 9 ], // (5_A_N, 5_A_CA, 5_A_CB)
     [  15.23981, 121.52916,  12.43648, "Csb", "Cdb", "Odb", 0, 0, 0, 0, 0, "CACO", 10 ], // (5_A_CA, 5_A_C, 5_A_O)
     [  15.15796, 109.99027,  15.23981, "Csb", "Csb", "Cdb", 1, 0, 0, 4, 0, "CBCAC", 11 ], // (5_A_CB, 5_A_CA, 5_A_C)
     [  10.19207, 116.71552,  14.46867, "Hsb", "Nsb", "Csb", 1, 0, 0, 5, 3, "HNCA", 12 ], // (5_A_H, 5_A_N, 5_A_CA)
   ],
*/
