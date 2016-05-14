#!/usr/bin/env luajit

parsers = require 'rfold.parsers'

--[[

   dssp output internal coords as lua tables ...
   { triple_id { dist angle dist } }
   { quad_id dihedral }


   triple : gen positions in standard coordinate space

   quad : identify 1st triple
          copy cordinates x3 
          identify 2nd triple but only for last atom
          copy cordinates last atom
          transform last atom copy to <dihedral> rotated position for 1st triple

   -- quad could reference both triples, copy only 4th position
   --- keep copy of a4 initial transformed before mrz -- ready for any dihedral setting
   
   chain <- specified cordinates 1st (or any) triple with triple ID
   --- needs to be first or have to build entire chain and then transform
   ---- else can just set first and build on to that
   ----- no because will render each residue in initialised std coordinate space

   generate sequence of residues from triple IDs

   k1 e2 t3 a4 a5 a6 k7

   residue :
      based on sequence position and residue, gather constituents:
      - reference backbone triples
      - reference backbone dihedrals -- include next res omg, phi if present
      -- only include if 1st dihedral atom in residue
      - reference sidechain triples and dihedrals

      generate all 3d coords in residue coordinate space ... but not the one I want ... at least CA at origin :
      - start with psi
      -- 1st triple in std coord space -- just copy
      -- a4 of 1st dihedral in coord space - just copy  (psi)
      -- dihedral ending in O in same space, just copy O
   
      - while (remaining dihedrals to process)
      -- remove from queue
      -- if 1st 3 atom coordinates known
      --- generate reverse transform for known coordinates
      --- apply to 4th atom and set [ will include positions for residue N+1, save with residue n ]
      ---- store dihedral angle with a4 coordinate to test for changed dihedral angle later
      -- else add to new queue
      -- if processed entire queue with no new coordinates, halt with error
      -- else new queue -> queue and continue

      - create reverse transform for residue n+1 n-ca-c to std space (apply reverse to add n+1 to n)
   
      add residue to chain:
      - get res n's n+1 transform matrix
      - apply to all n+1 atoms
      - concatenate to (n+1)'s n+1 transform matrix

   residue storage:
      - all dihedrals
      -- each dihedral's transform to attach to neighbour?
      - n+1 transform matrix
      - all atom coordinates in space for n-ca-c at origin 

   dihedral storage
      - a4 translated, ready for rotation
      - a4 in position
      - dirty flag when angle changed

   hede storage
      - all 3 atom coordinates in standard coordinate space

   
   n - ca - c - n  ***
                n - ca - c - n
       ca - c - n - ca
            c - n - ca - c
   
   
--]]


local pm = require "rfold.protein"


local args = parsers.parseCmdLine()
local pdbid = parsers.prd(args[1], function (t) pm.load(t) end)
local prot = pm.get(pdbid)
prot:linkResidues()
prot:renderDihedrons()

prot:setStartCoords()
prot:assembleResidues()

print(prot:toPDB())

--prot:setInitialPosition(1)
--prot:assembleChains()



--[[
print(prot:tostring())
local p2 = pm.protein:get(pdbid)
print(p2:tostring())
--]]
