#rFold
[rFold](http://rob-miller.github.io/rFold/doc/index.html) is a set of Lua (LuaJIT) classes and submodules to build, manipulate and analyse three-dimensional protein structures.

## Overview

rFold is intended to be a library of functions to assist in protein structure analysis and prediction.  The predictive methodology will be largely based on kinematic and artificial intelligence approaches, but hopefully interested users can pick and choose utilities from the project to meet their own needs.

The intended development roadmap is to create functionality in roughly the following order: Build (generate 3D structure coordinates from internal phi, psi, omega and chi dihedral angles); Collect (catalogue internal coordinate and other data from known protein structures); Move (perturb structures to explore local regions of conformational space); and Evaluate (assess three dimensional residue environments based on data from known structures).

Currently (June, 2016) only the Build module is working,  This software can accurately interconvert between PDB files and fully specified internal coordinates (_phi, psi, omega_ and _chi_ dihedral angles, plus all bond lengths and angles).  Planned functionality will include incorporating average bond lengths and angles and default dihedral angles if not specified.

## External Dependencies

[numlua](https://github.com/carvalho/numlua) -- Numeric Libraries for Lua for matrix support.  I needed to pull the lua_number2int fix for this repository; possibly this should be switched to use [Torch](http://torch.ch/).

[deque](https://github.com/catwell/cw-lua/tree/master/deque) -- Queue implementation used to assemble residues for Build module.

[ldoc](https://github.com/stevedonovan/LDoc) -- Lua documentation generator.

## Build

Accurate creation of 3D coordinate data from protein structure internal coordinates (phi and psi dihedral angles) requires explicit specification of all structure details: bond lengths, angles and dihedral angles.  rFold defines a '.pic' datafile format in which each standalone line specifies 3 atoms and their length, angle, length relationship, or 4 atoms with their dihedral angle relationship.  For example:

1MUD A 1MN:1MCA:1MC   1.45628 114.95350   1.53294

1MUD A 1MN:1MCA:1MC:2QN 144.22629

In addition the coordinates for the first three atoms (usually N, CA and C) of each chain are required from the original PDB file in order to place the chain in the correct coordinate space, if the coordinates are to be compared to the PDB data.

While there is an inevitable buildup of error as residues are assembled along a protein chain, the current implementation with five decimal point precision for the internal coordinates is sufficient to exactly regenerate backbone and sidechain coordinates for PDB chains of over 1200 residues.

Using a May, 2016 `cullpdb_pc20_res2.2_R1.0.curr` file from the [Dunbrack Lab PISCES server](http://dunbrack.fccc.edu/PISCES.php), 97.5% of 5,825 protein chains are regenerated exactly according to a line-by-line comparison test (see buildProtein.lua).  As of June, 2016, the 143 chains generating differences have not yet been investigated.

## Collect

Not yet implemented.  Source tree contains modified DSSP program, changed to output internal coordinate records matching PIC data above and more atom coordinates.

## Move

Not implented.

## Evaluate

Not implemented.  See [Miller, Douthart and Dunker 1994](https://books.google.com/books?id=VmFSNNm7k6cC&lpg=PA22&ots=kq7am9BPse&dq=miller%20douthart%20and%20dunker&pg=PA9) for early work.

## Why Lua / LuaJIT?

In addition to speed and its connection with the [Torch7](http://torch.ch/) Neural Nets / Scientific Computing Framework, I find Lua code remarkably easy to understand and modify long after initial development work.

## License

Copyright 2016 Robert T. Miller

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use these files except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
