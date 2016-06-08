#!/usr/bin/env luajit

local parsers = require 'rfold.parsers'
local protein = require "rfold.protein"
local utils = require 'rfold.utils'

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
         if 0>ecnt then return false end
      end
   end
   if #s1t ~= #s2t then
      print(' different linecounts.')
      return false
   end
   --return true
   return #s1t
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

--for k,v in pairs(args) do print(k,v) end
--os.exit()

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
         local c = lineByLine(s0,s1) 
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
      
   else                             -- did read PDB, generate internal coordinates
      prot:setStartCoords()            -- copy first residue N, CA, C coordinates from input data to each chain residue 1 initCoords list
      
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
         local c = lineByLine(s0,s1) 
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
end


