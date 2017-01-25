#!/usr/bin/env luajit

--[[
   rFoldLoadDB.lua
   
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

--- rFoldLoadDB.lua [ options ] [<PDB file | PIC file>] ...
-- 
-- Default action is to read a PDB file and load into configured rFold database.
--
-- options:
-- 
-- -f=<list of PDB codes> : read PDB codes and optional chain IDs
--                            from file and use to form filename to read
--                            from PDB_repository (location specified below),
--                            e.g. 7RSA becomes
--                            /media/data/pdb/rs/pdb7rsa.ent.gz 
--
-- -nd : do not run dssp on PDB files
-- -u  : update if entry in database already (default is to skip on chain by chain basis)
-- 

local parsers = require 'rfold.parsers'
local protein = require "rfold.protein"
local utils = require 'rfold.utils'
local chemdata = require 'rfold.chemdata'

local ps = require 'pipe_simple'

local PDB_repository_base = '/media/data/pdb/'
local mkdssp = 'dssp/rtm-dssp-2.2.1/mkdssp'

local rfpg = require 'rfold.rfPostgresql'  -- autocommit false by default



function verifyDb()
   local c=0
   for _ in pairs(chemdata.atomic_weight) do c = c + 1 end

   local rows = tonumber(rfpg.Q('select count(*) from periodic_table;')[1]);
   --print('periodic_table:', rows)

   if (rows < c) then 
         for k,v in pairs(chemdata.atomic_weight) do
            rfpg.Qcur("insert into periodic_table (atom, weight, electronegativity) values ('" .. k .. "'," .. v .. "," .. chemdata.electronegativity[k] .. ") on conflict (atom) do update set weight=" .. v .. ", electronegativity = " .. chemdata.electronegativity[k] .. ";")
         end
   end

   c=0
   for _ in pairs(chemdata.covalent_radii) do c = c + 1 end
   rows = tonumber(rfpg.Q('select count(*) from atom_class;')[1]);
   if (rows < c) then
         for k,v in pairs(chemdata.covalent_radii) do
            rfpg.Qcur("insert into atom_class (class, r_covalent, v_covalent) values ('" .. k .. "'," .. v .. "," .. chemdata.covalent_volume(k) .. ") on conflict (class) do update set r_covalent=" .. v .. ", v_covalent = " .. chemdata.covalent_volume(k) .. ";")
         end
   end

   for r1,r3 in pairs(chemdata.res3) do
      local pdb_atoms = { 'N', 'CA', 'C', 'O', 'CB', 'OXT' }
      local backbone = chemdata.residue_atom_class['X']
      for i,pa in ipairs(pdb_atoms) do
         if (not (r1 == 'G' and pa == 'CB')) then
            --print('r1',r1)
            --print('pa',pa)
            --print('b[pa]', backbone[pa])
            --print('pa[1]', string.sub(pa,1,1))
            rfpg.Qcur("insert into atoms (name, class, atom) values ('" .. r1 .. pa .. "','" .. backbone[pa] .. "','" .. string.sub(pa,1,1) .. "') on conflict (name) do nothing;")
         end
      end
      local sidechain = (r1 == 'G' or r1 == 'A' or r1 == 'X') and {} or chemdata.residue_atom_class[r1]
      for an,ac in pairs(sidechain) do
         --print('an[1]', string.sub(an,1,1))
         rfpg.Qcur("insert into atoms (name, class, atom) values ('" .. r1 .. an .. "','" .. ac .. "','" .. string.sub(an,1,1) .. "') on conflict (name) do nothing;")
      end
      --print()
   end
   rfpg.Qcur('commit')
end



local args = parsers.parseCmdLine(
   {
      ['nd'] = 'do not run dssp on PDB input file',
      ['u'] = 'update database (default is skip if entry already exists)'
   },{
      ['f'] = '<input file> : process files listed in <input file> (1 per line) followed by any on command line'
     },
   'convert PDB file to internal coordinates and load into database.  repository base: ' .. PDB_repository_base .. ' db: ' .. rfpg.db .. ' on ' .. rfpg.host .. ' ('  ..utils.getHostname() .. ')' 
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
   if a:ematch('^(%d%w%w%w)(%w?)%s+') then           -- looks like pdbcode with optional chain, read as compressed file from PDB_repository_base
      local pdbid, chn = _1:lower(), _2
      local subdir = pdbid:match('^%w(%w%w)%w$')
      toProcess[i] = PDB_repository_base .. '/' .. subdir .. '/pdb' .. pdbid:lower() .. '.ent.gz'
      if (chn) then
         toProcess[i] = toProcess[i] .. ' ' .. chn
      end
   end
end

verifyDb()

for i,arg in ipairs(toProcess) do
   if '' ~= arg then
      --print(arg)
      if 'quit' == arg then os.exit() end    -- so can insert 'quit' in input file of PDB IDs for testing and debugging
      local file,chain = arg:match('^(%S+)%s?(%w?)%s*$')
      if chain == '' then chain = nil end
      
      local pdbid = parsers.parseProteinData(file, function (t) protein.load(t) end, chain)
      local prot = protein.get(pdbid)
      --print(pdbid,prot,prot:countDihedra())
      --print(prot:tostring())


      --- need dihedra and (optionally not) DSSP data
      -- have either PDB coords (pdb file), dihedra data only (pic file), or dihedra and DSSP data (rtm mkdssp output file)
      -- want to test for data consistency each way, essentially repeating buildProtein.lua test mode

      local dsspCount = prot:countDSSPs()  -- useful to know below
      local coordsInternal
      local coords3D
      
      if prot:countDihedra() > 0  then  -- did read internal coordinates, either from .pic file or rtm DSSP output
         
         --- test if we can generate PDB coordinates and match back:
         prot:setStartCoords()                  -- copy residue 1 N, CA, C coordinates from input data to each chain residue 1 initCoords list (not done on loading pic file)
         coordsInternal = prot:writeInternalCoords()  -- get output for internal coordinate data as loaded 
         prot:internalToAtomCoords()            -- generate PDB atom coordinates from internal coordinates (needs dihedron data structures to complete chain)
         prot:clearInternalCoords()             -- wipe internal coordinates as loaded
         prot:atomsToInternalCoords()           -- generate internal coordinates from PDB atom coordinates
         local s1 = prot:writeInternalCoords()  -- get output for internal coordinate data as generated
         
         --print(coordsInternal)
         --print('------------')
         --print(s1)
         local c = utils.lineByLineCompare(coordsInternal,s1) 
         if not c then
            print((0 == dsspCount and 'PIC' or 'DSSP') .. ' file ' .. arg .. ' failed to re-generate matching internal coordinates from 3D coordinates.')
            goto continue
         end
         coords3D = prot:writePDB(true)
      else                             -- did read PDB file
         --- test if we can get internal coordinates and regenerate input pdb file
         prot:atomsToInternalCoords()     -- calculate bond lengths, bond angles and dihedral angles from input coordinates
         coords3D = prot:writePDB(true)   -- loaded PDB format text without REMARK RFOLD record (timestamp may not match)
         prot:setStartCoords()            -- copy first residue N, CA, C coordinates from input data to each chain residue 1 initCoords list
         prot:clearAtomCoords()
         prot:internalToAtomCoords()
         local s1 = prot:writePDB(true)
         --print(coords3D)
         --print('------------')
         --print(s1)
         local c = utils.lineByLineCompare(coords3D,s1)
         if not c then
            print('PDB file ' .. arg .. ' failed to regenerate 3D coordinates from calculated internal coordinates')
            goto continue
         end
         coordsInternal = prot:writeInternalCoords()  
      end

      --- if here then 3D coordinate match internal coordinates and have both
      -- now confirm mkdssp matches calculated internal coordinates

      if 0 == dsspCount and not args['nd'] then -- did read pic or pdb file, check dssp
         -- see if rtm mkdssp gets same internal coordinates
         local dsspstatus,dsspresult,dssperr = ps.pipe_simple(coords3D,mkdssp,unpack({'-i','-'}))
         local pdbid2 = parsers.parseProteinData(dsspresult, function (t) protein.load(t) end, chain, 'dssp_pipe')
         s1 = prot:writeInternalCoords()  -- get output for internal coordinate data as read from dssp
         local c = utils.lineByLineCompare(coordsInternal,s1) 
         if not c then
            print(arg .. ' failed to generate matching DSSP internal coordinates from 3D coordinates.')
            goto continue
         end
      --else    -- otherwise we have DSSP data already, must have read in, and already tested regenerate structure -> regenerate matching internal coordinates
      end
      

      --- if here then we have happy data for prot:
      print('we have happy data for ' .. arg)
      print(prot:tostring())

     -- prot:writeDB(rfpg,args['u'])
      
   end

   ::continue::
end


--[[

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
      protein.drop(pdbid)
   end
end
--]]

--[[
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
--]]
