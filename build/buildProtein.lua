#!/usr/bin/env luajit

--[[
   buildProtein.lua
   
Copyright 2016, 2017 Robert T. Miller

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use these files except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
]]

--- buildProtein.lua [ options ] <PDB file | PIC file> ...
-- 
-- Default action is to read a PDB file and print corresponding 'Protein
-- Internal Coordinates' (PIC) data, or a PIC file and print corresponding
-- PDB format output.
--
-- options:
-- 
-- -f=<list of PDB codes> : read PDB codes and optional chain IDs
--                            from file and use to form filename to read
--                            from PDB_repository (location specified below),
--                            e.g. 7RSA becomes
--                            /media/data/pdb/rs/pdb7rsa.ent.gz 
-- -w : instead of printing, write output to input filename with '.pic'.
--                            or '.gpdb' (generated PDB) extension 
-- -t : test mode - read input file, generate PDB data or internal 
--                            coordinates, attempt to re-generate input file, 
--                            compare to input and report initial differences
--                            if any 
-- 
-- Using the -t and -f option with a cullpdb_pc20_res2.2_R1.0.curr from May, 2016
-- (http://dunbrack.fccc.edu/PISCES.php -- Dunbrack Lab PISCES server),
-- 97.5% of 5,825 protein chains are regenerated correctly (backbone and sidechain
-- atoms) according to the line-by-line comparison test.



local parsers = require 'rfold.parsers'
local protein = require "rfold.protein"
local utils = require 'rfold.utils'

local PDB_repository_base = '/media/data/pdb/'

local args = parsers.parseCmdLine(
   '[ <pdbid[chain id]> | <pdb file> | <pic file> ] ...',
   '\n convert PDB files to internal coordinates and vice-versa. \n PDBid[chain id] file (/xx/<pdbid>.ent.gz) located under ' .. (PDB_repository_base and PDB_repository_base or '[not configured]') .. " \n 'pic file' in 'protein internal coordinates' format as generated by this program.\n",
   {
--      ['a'] = 'average: generate hedron_default.lua and dihedron_default.lua files of average values from input files',
      ['t'] = 'test: convert (PDB|internal_coordinates) input to (internal_coordinates|PDB) and back again, verify match',
      ['w'] = 'write: converted result to <pdbid>.pdb or <pdbid>.pic in current directory',
      ['wa'] = 'writeAll: all intermediate files for -t written in current directory'
   },
   {
      ['f'] = '<input file> : process files listed in <input file> (1 per line) followed by any on command line'
   }
)

local toProcess={}
if args['f'] then
   for i,f in ipairs(args['f']) do
      local fh=io.open(f)
      for line in fh:lines() do
         table.insert(toProcess, line)
      end
   end
end

for i,a in ipairs(args) do table.insert(toProcess, a) end 

for i,a in ipairs(toProcess) do
   if a:match('^IDs%s') then toProcess[i]='' end     -- header line for Dunbrack list : 'IDs         length Exptl.  resolution  R-factor FreeRvalue'
   if a:match('^#') then toProcess[i]='' end         -- comment line starts with '#'
   if a:ematch('^(%d%w%w%w)(%w?)%s*') then           -- looks like pdbcode with optional chain, read as compressed file from PDB_repository_base
      local pdbid, chn = _1:lower(), _2
      local subdir = pdbid:match('^%w(%w%w)%w$')
      toProcess[i] = PDB_repository_base .. '/' .. subdir .. '/pdb' .. pdbid:lower() .. '.ent.gz'
      if (chn) then
         toProcess[i] = toProcess[i] .. ' ' .. chn
      end
      --print(a,toProcess[i])
   end
end

for i,arg in ipairs(toProcess) do
   if '' ~= arg then
      --print(arg)
      if 'quit' == arg then os.exit() end    -- so can insert 'quit' in input file of PDB IDs for testing and debugging
      local file,chain = arg:match('^(%S+)%s?(%w?)%s*$')
      if chain == '' then chain = nil end
      --print(file,chain,arg)
      local pdbid = parsers.parseProteinData(file, function (t) protein.load(t) end, chain)   -- read protein data file either PDB or pic format into object in in-memory table 'proteins' (purged at end of loop)
      local prot = protein.get(pdbid)  -- select protein object from proteins table

      prot:setStartCoords()            -- copy first residue N, CA, C coordinates from input data to each chain residue 1 initCoords list

      -- print('countDihedra: ' .. prot:countDihedra())
      
      if prot:countDihedra() > 0 then  -- did read internal coordinates (pic format data), so generate PDB

         --prot:setStartCoords()               -- copy residue 1 N, CA, C coordinates from input data to each chain residue 1 initCoords list
         if args['a'] then
            -- not yet implemented
            --hedron_default = prot:addHedronData(hedron_default)
            --dihedron_default = prot:addDihedronData(dihedron_default)
         elseif args['t'] then
            io.write('testing ' .. arg .. ' ... ')
            io.flush()
            local s0 = prot:writeInternalCoords()  -- get output for internal coordinate data as loaded
            if args['wa'] then utils.writeFile('in.pic',s0) end
            prot:internalToAtomCoords()    -- generate PDB atom coordinates from internal coordinates
            prot:clearInternalCoords()          -- wipe internal coordinates as loaded
            if args['wa'] then utils.writeFile('gen_intermediate.pdb',prot:writePDB()) end
            prot:atomsToInternalCoords()    -- generate internal coordinates from PDB atom coordinates
            local s1 = prot:writeInternalCoords()  -- get output for internal coordinate data as generated
            if args['wa'] then utils.writeFile('gen_out.pic',s1) end
            local c = utils.lineByLineCompare(s0,s1) 
            if c then
               print('passed. ' .. c .. ' pic lines compared, ' .. prot:report())
            else
               print('failed')
            end
         else
            prot:internalToAtomCoords()
            local s = prot:writePDB()              -- get PDB format text
            if args['w'] then
               utils.writeFile(pdbid .. '.gpdb',s)
            else
               print(s)
            end
         end
         
      else                             -- did read PDB format data, generate internal coordinates
         --prot:setStartCoords()            -- copy first residue N, CA, C coordinates from input data to each chain residue 1 initCoords list
         
         prot:atomsToInternalCoords()          -- calculate bond lengths, bond angles and dihedral angles from input coordinates
         if args['a'] then
            --hedron_default = prot:addHedronData(hedron_default)
            --dihedron_default = prot:addDihedronData(dihedron_default)
         elseif args['t'] then
            io.write('testing ' .. arg .. ' ... ')
            io.flush()
            local s0 = prot:writePDB(true)     -- PDB format text without REMARK RFOLD record (timestamp may not match)
            if args['wa'] then utils.writeFile('in.pdb',s0) end
            prot:clearAtomCoords()
            prot:internalToAtomCoords()
            if args['wa'] then utils.writeFile('gen_intermediate.pic',prot:writeInternalCoords()) end
            local s1 = prot:writePDB(true)
            if args['wa'] then utils.writeFile('gen_out.pdb',s1) end
            local c = utils.lineByLineCompare(s0,s1) 
            if c then
               print('passed: ' .. c .. ' pdb lines compared, ' .. prot:report())
            else
               print('failed')
            end
         else
            local s = prot:writeInternalCoords()
            if args['w'] then
               utils.writeFile(pdbid .. '.pic',s)
            else
               print(s)
            end
         end         
      end
      protein.drop(pdbid)
   end
end

