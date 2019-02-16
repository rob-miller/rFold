
include <peptideOffsets.scad>;


// openSCAD control parameters
//$fn=4;

// output parameters
bondRadius=0.4;
clearance=0.25;
extrusionWidth=0.4;
minThickness=2*extrusionWidth;
atomScale=0.8;

// datafile array offsets
     
     
// covalent radii from Heyrovska, Raji : 'Atomic Structures of all the Twenty Essential Amino Acids and a Tripeptide, with Bond Lengths as Sums of Atomic Covalent Radii'
// https://arxiv.org/pdf/0804.2488.pdf

atomData = [ ["Csb","green" , 0.77], ["Cres","green" , 0.72], ["Cdb","green" , 0.67],
                      ["Osb","red" , 0.67], ["Ores","red" , 0.635], ["Odb","red" , 0.60],
                      ["Nsb","blue" , 0.70], ["Nres","blue" , 0.66], ["Ndb","blue" , 0.62],
                      ["Hsb","gray" , 0.37],
                      ["Ssb","yellow" , 1.04] ];


// atom with colour and size
module atom(a,scal)
{
     ad = atomData[search([a],atomData)[0]];
     color(ad[1]) {
          sphere(r=((ad[2]*atomScale)*scal));
     }
}


// /*
///////////// version 0 ... used by chain()

// one plane of 3 points
module hedron(h,rev=0,scal)
{
                 //   yes reversed           :  not reversed
                 //    0    1     2     3    :     0     1     2     3
                 //  len1  len3  atom1 atom3     len1  len3  atom1 atom3
     newh = (rev ? [ h[2], h[0], h[5], h[3] ] : [ h[0], h[2], h[3], h[5] ]);
     scal=scal; 
     
     // central atom at 0,0,0
     atom(h[4],scal);

     // comments for non-reversed case
     // atom 3 is len3 up on +z
     translate([0,0,newh[1]]) atom(newh[3],scal);
     
     // atom 2 - atom 3 bond from origin up +z distance len3
     cylinder(h=newh[1],r=bondRadius*scal,center=false);
     
     // atom 1 - atom 2 bond rotated about Y
     rotate([0, h[1], 0]) {                        // rotate following elements by angle2 about Y
          cylinder(h=newh[0],r=bondRadius*scal,center=false);        // a1-a2 bond is len1
          translate([0,0,newh[0]]) atom(newh[2],scal);  // put atom1 sphere len1 away on Z
     }
}


// hedron rotated to specific angle d
module d2(d,hedra,scal)
{
     tz = (d[d_reversed] ? hedra[d[d_h2ndx]][2] : hedra[d[d_h2ndx]][0]);   // get h2 len1 depending on reversed
     rotate(d[d_dangle1]) {                                   // 4. rotate h2 to specified dihedral angle
          translate([0,0,tz]) {                       // 3. translate h2 h2:len1 up +z
               rotate([180, 0, 0]) {                  // 2. rotate h2r about X so h2:a3 in +z and h2:a1 in -z
                    hedron(hedra[d[d_h2ndx]],(!d[d_reversed]),scal);      // 1. reverse hedron 2 orientation = h2r
               }
          }
     }
}

// two hedra at specified angle d
module dihedron(d,hedra,scal)
{
     //echo(d)
     hedron(hedra[d[d_h1ndx]],d[d_reversed],scal);                // reverse h1 if dihedral reversed
     //hull() {
     d2(d,hedra,scal);
     //}
}





///////////////////////////////////

// x versions used by amideChain()

// one plane of 3 points
// h = hedron data
// rev = h has reversed point order or not
// a1, a2, a3, ax = flags for conditional render at each step - a<n>=bonds, ax=atom sphere; e.g. don't draw atom spheres at joiner points
module hedronx(h,rev=0,a1=1,a2=1,a3=1,ax=1)
{
     
                 //   yes reversed           :  not reversed
                 //    0    1     2     3    :     0     1     2     3
                 //  len1  len3  atom1 atom3     len1  len3  atom1 atom3
     newh = (rev ? [ h[2], h[0], h[5], h[3] ] : [ h[0], h[2], h[3], h[5] ]);
     s = ( (1==ax) ? 1.0 : 0.9 );  // if not drawinng atom then make bond slightly shorter
     scal = 10; //h[6];
          
     if (1==a2)
          // central atom at 0,0,0
          atom(h[4],scal);
     
     // comments for non-reversed case

     if (1==a3) {
          // atom 3 is len3 up on +z
          if (1==ax)
               translate([0,0,newh[1]]) atom(newh[3],scal);
          // atom 2 - atom 3 bond from origin up +z distance len3
          cylinder(h=(newh[1] *s),r=bondRadius*scal,center=false);
     }
     
     if (1==a1)
          // atom 1 - atom 2 bond rotated about Y
          rotate([0, h[1], 0]) {                        // rotate following elements by angle2 about Y
               cylinder(h=(newh[0] *s),r=bondRadius*scal, center=false);        // a1-a2 bond is len1
               if (1==ax)
                    translate([0,0,newh[0]]) atom(newh[2],scal);  // put atom1 sphere len1 away on Z
          }
}

// hedron joiner 

module jhedron(h,rev=0)
{
     //   yes reversed           :  not reversed
     //    0    1     2     3    :     0     1     2     3
     //  len1  len3  atom1 atom3     len1  len3  atom1 atom3
     newh = (rev ? [ h[2], h[0], h[5], h[3] ] : [ h[0], h[2], h[3], h[5] ]);
     scal = 10;  //h[6];
          
     // central atom at 0,0,0
     //atom(h[4],scal);
     
     // comments for non-reversed case

     clen = 0.55;
     // atom 3 is len3 up on +z
     //translate([0,0,newh[1]]) atom(newh[3]);
     // atom 2 - atom 3 bond from origin up +z distance len3
     difference() {
          cylinder(h=newh[1]*clen,r=(bondRadius*scal)+clearance+minThickness, center=false);
          translate([0,0,-(newh[1]*clen*0.5)])
               cylinder(h=newh[1],r=(bondRadius*scal)+clearance,center=false);
               
     }

     // atom 1 - atom 2 bond rotated about Y
     rotate([0, h[1], 0]) {                        // rotate following elements by angle2 about Y
          //cylinder(h=newh[0],r=bondRadius);        // a1-a2 bond is len1
          difference() {
               cylinder(h=newh[1]*clen,r=(bondRadius*scal)+clearance+minThickness, center=false);
               translate([0,0,-(newh[1]*clen*0.5)])
                    cylinder(h=newh[1],r=(bondRadius*scal)+clearance, center=false);
          }
          //translate([0,0,newh[0]]) atom(newh[2]);  // put atom1 sphere len1 away on Z
     }
}

module d2x(d,hedra,a1=1,a2=1,a3=1,ax=1)
{
     tz = (d[3] ? hedra[d[2]][2] : hedra[d[2]][0]);            // get h2 len1 depending on reversed
     rotate(d[0]) {                                            // 4. rotate h2 to specified dihedral angle
          translate([0,0,tz]) {                                // 3. translate h2 h2:len1 up +z
               rotate([180, 0, 0]) {                           // 2. rotate h2r about X so h2:a3 in +z and h2:a1 in -z
                    hedronx(hedra[d[2]],(!d[3]),a1,a2,a3,ax);  // 1. reverse hedron 2 orientation = h2r
               }
          }
     }
}

module dihedronx(d,hedra,a1=1,a2=1,a3=1,a4=1,ax=1)
{
     //echo(d)
     if (1==a1)
          hedronx(hedra[d[1]],d[3],a2,a3,a4,ax);                // reverse h1 if dihedral reversed
     //hull() {
     d2x(d,hedra,a2,a3,a4,ax);
     //}
}




////////////////////////////////////



module amide(a,residues,hedra)
{
     //hedronx(hedra[a[1]],0);

     //translate([0,0,60])
     //jhedron(hedra[a[1]],0);

     //echo(residues[a[0]][a[2]][4]);

//hull() {
          // omega dihedral
          multmatrix(residues[a[0]][a[2]][4]) 
               dihedronx(residues[a[0]][a[2]],hedra,0,0,0,1);
//}
// hull() {
 
          // NCaCO fragment
          multmatrix(residues[a[0]][a[3]][4]) 
               d2x(residues[a[0]][a[3]],hedra,1,1,0);
//}

//hull() {
          // NCaCO fragment
          //multmatrix(residues[a[0]][a[3]][4]) 
          //     d2x(residues[a[0]][a[3]],hedra,1,1,0);

// hull() {
          // CCaNH (r=i+1) fragment
          multmatrix(a[5]) 
               d2x(residues[a[0]+1][a[4]],hedra,1,1,0);
//}     
  //   }
          
//hull() { 
     // omega dihedral
     multmatrix(residues[a[0]][a[2]][4]) 
          dihedronx(residues[a[0]][a[2]],hedra,1,1,0,0,0);
//}
}

module residue(r,hedra, scal)
{
     for (d = r) {
          multmatrix(d[d_dihedralTransform]) {
               dihedron(d, hedra, scal);
          }
     }
}


module chain(protein) 
{
    chnD = protein[p_chainData];
    //a = chnD[cd_amideSet];
    c = chnD[c_residues];
    dihedra = chnD[c_dihedra];
    hedra = chnD[c_hedra];
     for (r = c) {
          multmatrix(r[r_resTransform]) {
               residue(dihedra[r[r_resNdx]],hedra, protein[p_proteinScale]);
          }
     }
}


//module amideChain(a, chains, residues, hedra)
module amideChain(protein)
{
     chnD = protein[p_chainData];
     a = chnD[c_amideSet];
     chains = chnD[c_chain];
     residues = chnD[c_residueSet];
     hedra = chnD[c_hedra];
     
     for (r = a) {
          multmatrix(chains[r[a_resNdx]][c_resTransform] ) {
               amide(r,residues,hedra);
          }
     }
}
