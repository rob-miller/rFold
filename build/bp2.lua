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


local protein = require "rfold.protein"

local hedron_default = {}
local dihedron_default = {}

local function lineByLine(s1,s2)
   local s1t, s2t = {}, {}
   for line in s1:gmatch("[^\r\n]+") do s1t[#s1t+1] = line end
   for line in s2:gmatch("[^\r\n]+") do s2t[#s2t+1] = line end

   local ecnt = 2
   for i,lin in ipairs(s1t) do
      if lin ~= s2t[i] then
         print()
         print(lin)
         print(s2t[i])
         ecnt = ecnt-1
         if 0>ecnt then os.exit() end
      end
   end
   return true
end


local args = parsers.parseCmdLine(
   {
      ['a'] = 'average: generate hedron_default.lua and dihedron_default.lua files of average values from input files',
      ['t'] = 'test: convert (PDB|internal_coordinates) input to (internal_coordinates|PDB) and back again, verify match',
      ['w'] = 'write: converted result to <pdbid>.pdb or <pdbid>.pic in current directory'
   },{
      ['f'] = '<input file> : process files listed in <input file> (1 per line) followed by any on command line'
     },
   'convert PDB files to internal coordinates and vice-versa'
)

local toProcess={}
if args['f'] then
   for i,f in args['f'] do
      io.open(f)
      for line in io.lines do table.insert(toProcess, line) end
   end
end

for i,a in ipairs(args) do table.insert(toProcess, a) end

for i,a in ipairs(toProcess) do

   local pdbid = parsers.parseProteinData(a, function (t) protein.load(t) end)
   local prot = protein.get(pdbid)

   if prot:countDihedra() > 0 then  -- did read internal coordinates, so generate PDB

      prot:setStartCoords()               -- copy residue 1 N, CA, C coordinates from input data to each chain residue 1 initCoords list
      if args['a'] then
         --hedron_default = prot:addHedronData(hedron_default)
         --dihedron_default = prot:addDihedronData(dihedron_default)
      elseif args['t'] then
         io.write('testing ' .. a .. ' ... ')
         io.flush()
         local s0 = prot:writeInternalCoords()  -- get output for internal coordinate data as loaded
         prot:internalToAtomCoords()    -- generate PDB atom coordinates from internal coordinates
         prot:clearInternalCoords()          -- wipe internal coordinates as loaded
         prot:atomsToInternalCoords()    -- generate internal coordinates from PDB atom coordinates
         local s1 = prot:writeInternalCoords()  -- get output for internal coordinate data as generated
         if lineByLine(s0,s1) then
            print('passed.')
         else
            print('failed')
         end
      else
         prot:internalToAtomCoords()
         local s = prot:writePDB()              -- get PDB format text
         if args['w'] then
            io.open(pdbid .. '.gpdb')
            io.write(s)
            io.close()
         else
            print(s)
         end
      end
      
   else                             -- did read PDB, generate internal coordinates
      prot:setStartCoords()            -- copy residue 1 N, CA, C coordinates from input data to each chain residue 1 initCoords list
      prot:atomsToInternalCoords()          -- calculate bond lengths, bond angles and dihedral angles from input coordinates
      if args['a'] then
         --hedron_default = prot:addHedronData(hedron_default)
         --dihedron_default = prot:addDihedronData(dihedron_default)
      elseif args['t'] then
         io.write('testing ' .. a .. ' ... ')
         io.flush()
         local s0 = prot:writePDB(true)     -- PDB format text without REMARK RFOLD record (timestamp may not match)
         prot:clearAtomCoords()
         prot:internalToAtomCoords()
         local s1 = prot:writePDB(true)
         if lineByLine(s0,s1) then
            print('passed.')
         else
            print('failed')
         end
      else
         local s = prot:writeInternalCoords()
         if args['w'] then
            io.open(pdbid .. '.pic')
            io.write(s)
            io.close()
         else
            print(s)
         end
      end         
   end
end


