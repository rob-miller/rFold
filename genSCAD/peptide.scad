
bondRadius=0.125;
atomRadius=2*bondRadius;

function atomColor(a) = a == "C" ? "green" : (a == "N" ? "blue" : (a == "O" ? "red" : (a == "S" ? "yellow" : "pink")));

module atom(a)
{
     color(atomColor(a)) {
          sphere(r=atomRadius);
     }
}

module hedron(h,rev=0)
{
     // central atom at 0,0,0
     atom(h[4]);
                 //  len1  len3  atom1 atom3   // len3  len1  atom3 atom1
     newh = (rev ? [ h[0], h[2], h[3], h[5] ] : [ h[2], h[0], h[5], h[3] ]);

     // comments for non-reversed case
     // atom 3 on +z
     translate([0,0,newh[0]]) atom(newh[2]);
     
     // atom 2 - atom 3 bond from origin up +z
     cylinder(h=newh[0],r=bondRadius);
          
     // atom 1 - atom 2 bond rotated about Y
     rotate([0, h[1], 0]) {
          cylinder(h=newh[1],r=bondRadius);  
          translate([0,0,newh[1]]) atom(newh[3]);
     }
}

module dihedron(d,hedra)
{
     hedron(hedra[d[1]],d[3]);                // reverse h1 if dihedral reversed
     rotate(d[0]) {                        // 4. rotate h2 to specified dihedral angle
          translate([0,0,hedra[d[2]][0]]) {       // 3. translate h2 h2:len1 up +z
               rotate([180, 0, 0]) {   // 2. rotate h2r about X so h2:a3 in +z and h2:a1 in -z
                    hedron(hedra[d[2]],!(d[3]));  // 1. reverse hedron 2 orientation = h2r
               }
          }
     }
}


module residue(r,hedra)
{
     for (d = r) {
          multmatrix(d[4]) {
               dihedron(d,hedra);
          }
     }
}

module chain(c, residues, hedra)
{
     for (r = c) {
          multmatrix(r[1]) {
               residue(residues[r[0]],hedra);
          }
     }
}
