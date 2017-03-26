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
rfpg.dbReconnect(1)

------- schema initialisation routines --------

function atomPairValueString(a1,a2)
   local t = utils.orderPair(a1,a2)
   local vals = "'{" .. '"' .. t[1] .. '", "' .. t[2] .. '"' .. "}'"
   return vals
end

function bondPairValueString(a1,a3,b1,b2)
   local t = utils.orderPair(a1,a3)  -- order bond pair by end atom sorting
   local vals
   if (a1 == t[1]) then  -- orderis b1=a1-a2, b2=a2-a3
      vals = "'{" .. b1 .. ", " .. b2 .. "}'"
   else                    -- orderis b2=a2-a3, b1=a1-a2
      vals = "'{" .. b2 .. ", " .. b1 .. "}'"
   end
   return vals
end

--- ensure there is an angle_class entry and subsidiary bond_class entries in database for passed angle string specification
-- @param(a) table 3 strings specifying order of residue-atoms comprising angle
-- @return database ID for angle_class
function getOrCreateAngleId(a)
   local angleAtomStr = rfpg.pgArrayOfStrings(a[1],a[2],a[3])
   local qry = rfpg.Q("select id from angle_string where atom_string = '" .. angleAtomStr .. "'")
   local rid
   
   if (not qry) then  -- angle not already configured
      local b = {}
      for k=1,2 do
         local bondAtomStr = rfpg.pgArrayOfStrings(a[k],a[k+1])
         qry = rfpg.Q("select id from bond_string where atom_string = '" .. bondAtomStr .. "'")
         if (qry) then     -- bond already configured for this string
            b[k] = qry[1]
         else              -- need to possibly create bond, or just attach ID to different order for string
            local vals = atomPairValueString(a[k],a[k+1])
            b[k] = rfpg.Q("insert into bond_class (res_atoms) values (" .. vals .. ") on conflict(res_atoms) do update set res_atoms = " .. vals .. " returning id;")[1]  -- id may exist already for different order, update will not change
            rfpg.Qcur("insert into bond_string (atom_string, id) values ('" .. bondAtomStr .. "'," .. b[k] ..')')
         end
      end

      -- now create the angle and set string for ID, as we know it is not already configured if here
      local bpvstr = bondPairValueString(a[1],a[3],b[1],b[2])
      rid = rfpg.Q('insert into angle_class (res_bonds) values (' .. bpvstr .. ') on conflict(res_bonds) do update set res_bonds = ' .. bpvstr .. ' returning id' )[1]
      rfpg.Qcur("insert into angle_string (atom_string,id) values ('" .. angleAtomStr .. "'," .. rid .. ')')
   else
      -- get existing ID for angle string from successful database query
      rid = qry[1]
   end
   return rid
end

--- ensure there is a dihedral_class entry and subsidiary angle_class and bond_class entries in database for passed dihedral angle string specification
-- @param(a) table 4 strings specifying order of residue-atoms comprising dihedral angle
-- @param(res) string indicating residue and adjacent neighbour if involved in bond
-- @param(nameId) int optional database id referencing name of dihedral angle (e.g. psi or chi1)
function rfdb.getOrCreateDihedralId(a,res,nameId)
   local dihedralAtomStr = rfpg.pgArrayOfStrings(a[1],a[2],a[3],a[4])
   local qry = rfpg.Q("select id from dihedral_string where atom_string = '" .. dihedralAtomStr .. "'")
   if (not qry) then
      local angIDs={}
      for j=1,2 do
         angIDs[j] = getOrCreateAngleId({a[j],a[j+1],a[j+2]})
         --rfpg.Q("select id from angle_string where atom_string='" .. rfpg.pgArrayOfStrings(a[j],a[j+1],a[j+2]) .."'")[1]
      end

      -- figure out how bonds and angles got ordered
      local bndNdx={}
      for i=1,3 do
         bndNdx[i] = (a[i] == utils.orderPair(a[i],a[i+1])[1] and '{1,2}' or '{2,1}')
      end
      local angNdx={}
      angNdx[1] = (a[1] == utils.orderPair(a[1],a[3]) and "{" .. bndNdx[1] .. ',' .. bndNdx[2] .. "}" or "{" .. bndNdx[2] .. ',' .. bndNdx[1] .. "}")
      angNdx[2] = (a[2] == utils.orderPair(a[2],a[4]) and "{" .. bndNdx[2] .. ',' .. bndNdx[3] .. "}" or "{" .. bndNdx[3] .. ',' .. bndNdx[2] .. "}")

      -- order angles within dihedral and log
      local dihedNdx
      local dihedIDs
      if (a[1] == utils.orderPair(a[1],a[4])) then
         dihedNdx = '{' .. angNdx[1] .. ',' .. angNdx[2] .. '}'
         dihedIDs = '{' .. angIDs[1] .. ',' .. angIDs[2] .. '}'
      else
         dihedNdx = '{' .. angNdx[2] .. ',' .. angNdx[1] .. '}'
         dihedIDs = '{' .. angIDs[2] .. ',' .. angIDs[1] .. '}'
      end

      local dihedId = rfpg.Q('insert into dihedral_class(res_angles, angle_atom_order, res3' .. (nameId and ',name' or '') .. ") values ('" .. dihedIDs .. "','" .. dihedNdx .. "','" .. res .. (nameId and "'," .. nameId or "'" ) .. ") on conflict(res_angles) do update set res3='" .. res .. "' returning id")[1]
      rfpg.Qcur("insert into dihedral_string (atom_string, id) values ('" .. dihedralAtomStr .. "'," .. dihedId .. ')')
   end
end

--- deterministically populate database tables for bond, angle and dihedral angle classes with serial IDs so IDs do not change with change in protein sets
function verifyDb()
   -- load dihedral_names
   rfpg.Qcur("insert into dihedral_name (name) values ('psi') on conflict (name) do nothing")
   rfpg.Qcur("insert into dihedral_name (name) values ('omega') on conflict (name) do nothing")
   rfpg.Qcur("insert into dihedral_name (name) values ('phi') on conflict (name) do nothing")
   for i=1,5 do
      rfpg.Qcur("insert into dihedral_name (name) values ('chi" .. i .. "') on conflict (name) do nothing")
   end
      
   -- load periodic_table
   local c=0
   for _ in pairs(chemdata.atomic_weight) do c = c + 1 end

   local rows = tonumber(rfpg.Q('select count(*) from periodic_table;')[1]);
   --print('periodic_table:', rows)

   if (rows < c) then
      for k,v in pairs(chemdata.atomic_weight) do
         rfpg.Qcur("insert into periodic_table (atom, weight, electronegativity) values ('" .. k .. "'," .. v .. "," .. chemdata.electronegativity[k] .. ") on conflict (atom) do update set weight=" .. v .. ", electronegativity = " .. chemdata.electronegativity[k] .. ";")
      end
   end

   -- load atom_bond_state (atom classes based on Heyrovska, Raji covalent radii paper https://arxiv.org/pdf/0804.2488.pdf )
   c=0
   for _ in pairs(chemdata.covalent_radii) do c = c + 1 end
   rows = tonumber(rfpg.Q('select count(*) from atom_bond_state;')[1]);
   if (rows < c) then
      for k,v in pairs(chemdata.covalent_radii) do
         rfpg.Qcur("insert into atom_bond_state (state, r_covalent, v_covalent) values ('" .. k .. "'," .. v .. "," .. chemdata.covalent_volume(k) .. ") on conflict (state) do update set r_covalent=" .. v .. ", v_covalent = " .. chemdata.covalent_volume(k) .. ";")
      end
   end

   for r1,r3 in pairs(chemdata.res3) do  -- first round through all residues: sidechain atoms into res_atoms, backbone+Cbeta into res_atoms
      -- load backbone and c-beta (Ala atoms) into res_atoms
      local ala_atoms = { 'N', 'CA', 'C', 'O', 'CB', 'OXT' }
      local backbone = chemdata.residue_atom_bond_state['X']
      for i,pa in ipairs(ala_atoms) do
         --if (not (r1 == 'G' and pa == 'CB')) then
            rfpg.Qcur("insert into res_atoms (name, bond_state, atom) values ('" .. r1 .. pa .. "','" .. backbone[pa] .. "','" .. string.sub(pa,1,1) .. "') on conflict (name) do nothing;")
         --end
      end
      -- load sidechain atoms int res_atoms
      local sidechain = (r1 == 'G' or r1 == 'A' or r1 == 'X') and {} or chemdata.residue_atom_bond_state[r1]
      for an,ac in pairs(sidechain) do
         --print('an[1]', string.sub(an,1,1))
         rfpg.Qcur("insert into res_atoms (name, bond_state, atom) values ('" .. r1 .. an .. "','" .. ac .. "','" .. string.sub(an,1,1) .. "') on conflict (name) do nothing;")
      end
   end
   
   for r1,r3 in pairs(chemdata.res3) do  -- second round through all residues: populate bond_class
      -- load res_bond_class, res_angle_class for backbone and c-beta (ala atoms)
      -- if _ in any atom position in triple, then need to cycle that position for all 20 neighbouring residues
      for r12,r32 in pairs(chemdata.res3) do
         for i,aset in ipairs(chemdata.backbone_angles) do
            local a = {}
            for j=1,3 do
               a[j] = string.gsub(aset[j], '_', r12)              -- substitute _ with each of 20 amino acids (r12)
               if (a[j] == aset[j]) then a[j] = r1 .. a[j] end    -- if no substitution, atom key is for current r1 residue
            end

            getOrCreateAngleId(a)

         end
      end
   end
   
   for r1,r3 in pairs(chemdata.res3) do  -- third round through all residues: populate dihedral_class
      -- same scan for loading dihedral_class
      for r12,r32 in pairs(chemdata.res3) do
         for i,dset in ipairs(chemdata.backbone_dihedrals) do
            local subpos=nil
            local a = {}
            for j=1,4 do
               a[j] = string.gsub(dset[j], '_', r12)              -- substitute _ with each of 20 amino acids (r12)
               if (a[j] == dset[j]) then                          -- if no substitution, atom key is for current r1 residue
                  a[j] = r1 .. a[j]
               else
                  subpos = j
               end
            end

            -- get central residue and neighbour if specified
            local res = (subpos and (2<subpos and "{ NULL," .. '"' .. r1 .. '","' .. r12 .. '"}' or "{" ..'"' .. r12 .. '","' .. r1 .. '", NULL}') or "{ NULL, " .. '"' .. r1 .. '", NULL}')
            local nameId = nil
            if (dset[5]) then nameId = rfpg.Q("select id from dihedral_name where name='" .. dset[5] .."'")[1] end
            getOrCreateDihedralId(a,res,nameId)

         end
      end 
   end

   for r,set in pairs(chemdata.sidechains) do
      for i,v in ipairs(set) do
         local a = {}
         if (v[4]) then -- dihedral angle
            for k=1,4 do
               a[k] = r .. v[k]
            end
            local nameId = nil
            if (v[5]) then nameId = rfpg.Q("select id from dihedral_name where name='" .. v[5] .."'")[1] end
            getOrCreateDihedralId(a,"{ NULL, " .. '"' .. r .. '", NULL}', nameId)
         else  -- regular angle
            for k=1,3 do
               a[k] = r .. v[k]
            end
            getOrCreateAngleId(a)
         end
      end
   end
   
   rfpg.Qcur('commit')
end

--------- end schema initialisation routines ------------------



--------- main ------ 

local args = parsers.parseCmdLine(
   '[ <pdbid[chain id]> | <pdb file> | <pic file> ] ...',
   'convert PDB file to internal coordinates and load into database.  repository base: ' .. PDB_repository_base .. ' db: ' .. rfpg.db .. ' on ' .. rfpg.host .. ' ('  ..utils.getHostname() .. ')',
   {
      ['nd'] = 'do not run dssp on PDB input file',
      ['u'] = 'update database (default is skip if entry already exists)'
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

      prot:setStartCoords()                  -- copy residue 1 N, CA, C coordinates from input data to each chain residue 1 initCoords list (not done on loading pic file)

      if prot:countDihedra() > 0  then  -- did read internal coordinates, either from .pic file or rtm DSSP output
         
         --- test if we can generate PDB coordinates and match back:
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
         --prot:setStartCoords()            -- copy first residue N, CA, C coordinates from input data to each chain residue 1 initCoords list  
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
         --local pdbid2 = parsers.parseProteinData(dsspresult, function (t) protein.load(t) end, chain, 'dssp_pipe')  -- this line re-initialises prot: because reloading same pdbid
         pdbid = protein.stashId(pdbid)  -- so instead we stash protein.proteins[pdbid] entry as new name first
         local pdbid2 = parsers.parseProteinData(dsspresult, function (t) protein.load(t) end, chain, 'dssp_pipe')  -- now we get new structure loaded
         local prot2 = protein.get(pdbid2)
         prot2:setStartCoords()
         s1 = prot2:writeInternalCoords()  -- get output for internal coordinate data as read from dssp
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

      --prot:setStartCoords()
      prot:clearAtomCoords()
      --coordsInternal = prot:writeInternalCoords()
      --print(coordsInternal)
      
      prot:writeDb(rfpg,args['u'])
      
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
