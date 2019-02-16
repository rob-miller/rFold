--[[
   residue.lua
   
Copyright 201,20176 Robert T. Miller

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

--- residue structure classes for rFold 
--
-- **external dependencies** : deque (https://github.com/catwell/cw-lua/tree/master/deque)
--
-- @classmod Residue

local utils = require 'rfold.utils'
local chemdata = require 'rfold.chemdata'
local geom3d = require 'rfold.geom3d'
local hedron = require 'rfold.hedron'
local dihedron = require 'rfold.dihedron'
local deque = require 'deque'

local residue = {}  -- module

local Residue = {}  -- class table

------------------------------------------------------------------------------------------------
-- residue
------------------------------------------------------------------------------------------------
-- class properties

--- array of hedron objects indexed by hedron keys ([resPos res atom] x3)
--@field hedra array of hedron objects

--- array of dihedron objects indexed by dihedron keys ([resPos res atom] x4)
--@field dihedra array of dihedron objects

--- array of { dihedron object, position 1..4 in dihedron specifying atom } for backbone atoms indexed by integer position of print order in PDB format file; set by linkDihedra()
-- @field backbone {dihedron, i}

--- array of { dihedron object, position 1..4 in dihedron specifying atom } for sidechain atoms indexed by integer position of print order in PDB format file; set by linkDihedra()
-- @field sidechain {dihedron, i}

--- array of 4x1 float matrices indexed by atomKey strings (resPos res atom) 
-- @field atomCoords float[4][1]

--- boolean indicating residue is ordered
-- @field ordered boolean

--- residue character GAVLIMFPSTCNQYWDEHKRX
-- @field res char

--- residue sequence position as in PDB or DSSP file
-- @field resn integer

--- DSSP data for this residue as read by parseProteinData(), or subset (ss, ss2, psi, phi, omg, acc) if loaded from database
-- @field dssp array 

--- map of key(first 3 atom ids) to dihedron for all dihedra in Residue
-- @field key3index array


-- @table backboneSort order of backone atoms for pdb files
local backboneSort = { N = 1, CA = 2, C = 3, O = 4 }

-- @table sidechainSort order of sidechain atoms for pdb files
local sidechainSort = { CB = 1,
                        CG = 2, CG1 = 2, OG = 2, OG1 = 2, SG = 2,
                        CG2 = 3,
                        CD = 4, CD1 = 4, SD = 4, OD1 = 4, ND1 = 4,
                        CD2 = 5,  ND2 = 5, OD2 = 5,
                        CE = 6, NE = 6, CE1 = 6, OE1 = 6, NE1 = 6,
                        CE2 = 7, OE2 = 7, NE2 = 7,
                        CE3 = 8,
                        CZ = 9, CZ2 = 9, NZ = 9,
                        NH1 = 10, OH = 10, CZ3 = 10, 
                        CH2 = 11, NH2 = 11,
                        OXT = 12
}


--- Residue class object initialiser (not a class method)
-- @param o table with fields 'res' = 1-letter amino acid code, 'resn' = sequence position of Residue in Chain
-- @return minimally initialised Residue object
function residue.new (o)
   assert ((o and o['res'] and o['resn']), "residue.new() called without 'res' and/or 'resn'")
   setmetatable(o, { __index = Residue })

   if not o['hedra'] then
      o['hedra'] = {}
   end
   if not o['dihedra'] then
      o['dihedra'] = {}
   end

   if not o['backbone'] then
      o['backbone'] = {}
   end
   if not o['sidechain'] then
      o['sidechain'] = {}
   end
   if not o['atomCoords'] then
      o['atomCoords'] = {}
   end

   return o
end

--- comparison function for sorting hedron / dihedron / atom keys
-- @param a 2ECA:2EC:3TN:3TCA or 2ECA 
-- @param b same
-- @return boolean result for '<'
function residue.keysort(a,b)
   if a==b then return false end  -- if = then not <

   local aksan, aksa = a:match('^(-?%d+)%a(%w+):?')
   if (not aksan) or (not aksa) then aksan, aksa = a:match('^(-?%d+)_(%w+):?') end
   local aksbn, aksb = b:match('^(-?%d+)%a(%w+):?')
   if (not aksbn) or (not aksb) then aksbn, aksb = b:match('^(-?%d+)_(%w+):?') end
   --print(a,b,aksan,aksa,aksbn,aksb)

   if aksan ~= aksbn then return aksan < aksbn end   -- seqpos takes precedence
      
   if aksa == aksb then return residue.keysort(a:match('^-?%d+%a%w+:(.+)$'),b:match('^-?%d+%a%w+:(.+)$')) end  -- if first field = then go to next
      
   if backboneSort[aksa] and backboneSort[aksb] then return backboneSort[aksa] < backboneSort[aksb] end
   if sidechainSort[aksa] and sidechainSort[aksb] then return sidechainSort[aksa] < sidechainSort[aksb] end
   if backboneSort[aksa] then return true end  -- backbone aksa, sidechain aksb
   return false -- sidechain aksa, backbone aksb
end





--- callback from file parser, import table data for this Residue according to contents
-- @param t parsed file record data: DSSP record, PDB ATOM record, or hedron / dihedron specification
function Residue:load(t)
   if t['psi'] then -- dssp record
      self['dssp'] = t
   elseif t['tempfact'] then   -- pdb ATOM record
      local ak = t['resn'] .. t['res'] .. t['atom']
      local atom = geom3d.get41mtx()
      atom[1][1] = t['x']
      atom[2][1] = t['y']
      atom[3][1] = t['z']
      self['atomCoords'][ak] = atom
   else
      t['pdbid'] = nil  -- not desired in di,hedron objects
      t['chn'] = nil  
      if t['dihedral1'] then
         --print('dihedron: ' .. self:tostring() .. ' ' .. utils.genKey(t[1],t[2],t[3],t[4]))
         self['dihedra'][ utils.genKey(t[1],t[2],t[3],t[4]) ] = dihedron.new(t)
      elseif t['angle2'] then
         --print('hedron: ' .. self:tostring() .. ' ' ..  utils.genKey(t[1],t[2],t[3]))
         self['hedra'][ utils.genKey(t[1],t[2],t[3]) ] = hedron.new(t)
      else
         for k,v in pairs(t) do print(k,v) end
         assert(nil, self:tostring() .. ' loading ' .. tostring(t) .. ' -- not recognised')
      end
   end
end


--- populate Residue hedra and dihedra tables with keys only for supplied res and resn
-- used to generate database average results, not otherwise in normal use
-- NB: *** neighbouring residues set to wildcard
function Residue:initEmpty()
   for x,tbl in ipairs( { chemdata.backbone_angles, chemdata.backbone_dihedrals, chemdata.sidechains[self['res']] } ) do
      if tbl then
         for i,hd in ipairs(tbl) do   -- hedra or dihedra
            local t = {}
            for k,a in ipairs(hd) do
               if 5 == k then
                  t['name'] = a
               else
                  local resn = self['resn']
                  local res = self['res']
                  if '_' == a:sub(1,1) then         -- handle _CA, _C and _N
                     resn = (k<4 and resn-1 or resn+1)
                     res = '_'   -- *** neighbouring residues set to wildcard
                     a = a:sub(2)  -- trim _ as it is residue
                  end  
                  t[k] = resn .. res .. a
               end
            end
            if t[4] then
               t['dihedral1'] = {}
            else
               t['len1'] = {}
               t['angle2'] = ''
               t['len3'] = {}
            end
            self:load(t)
         end
      end
   end
end

--- populate Residue hedra and dihedra values with average results using rFold db residue selection query
-- add hedron and dihedron tags sd, min, max
-- @param rfpg open database handle
-- @param resSelector rFold database query returning res_id, e.g. "select res_id from dssp where struc='H' and struc2=' X S+   '" limits to residues inside alpha helices
function Residue:getDbStats(rfpg, resSelector)
   for k,d in pairs(self['dihedra']) do
      d:getDbStats(rfpg, resSelector)
   end
   for k,h in pairs(self['hedra']) do
      h:getDbStats(rfpg, resSelector)
   end
end

--- write residue data to rfold database
-- @param rfpg open database handle
-- @param res_id residue ID in residue table
-- @param update optional flag, if false silently skip if [res_id [atom] ] entry exists already in atom_coordinates / dssp / dihedral / angle / bond tables 
function Residue:writeDb(rfpg, res_id, update)
   if self['dssp'] then
      local chk = rfpg.Q("select 1 from dssp where res_id=" .. res_id)
      if (not chk) then
         rfpg.Qcur("insert into dssp (res_id, struc, struc2, psi, phi, omega, acc) values (" .. res_id .. ",'" .. self['dssp']['ss'] .. "','" .. self['dssp']['ss2'] .. "'," .. self['dssp']['psi'] .. "," .. self['dssp']['phi'] .. "," .. self['dssp']['omg'] .. "," .. self['dssp']['acc'] .. ")" )
      elseif update then
         rfpg.Qcur("update dssp set (struc, struc2, psi, phi, omega, acc) = ('" .. self['dssp']['ss'] .. "','" .. self['dssp']['ss2'] .. "'," .. self['dssp']['psi'] .. "," .. self['dssp']['phi'] .. "," .. self['dssp']['omg'] .. "," .. self['dssp']['acc'] .. ") where res_id= " .. res_id  )
      end
   end
   if {} ~= self['atomCoords'] then  -- at chain level wrote atom_oordinates for initNCaC if present; here write all atom coordinates to db if present
      for k,v in pairs(self['atomCoords']) do
         local chk = rfpg.Q("select 1 from atom_coordinates where res_id = " .. res_id .. " and atom = '" .. k .. "' and not init")
         if not chk then
            rfpg.Qcur("insert into atom_coordinates (res_id, atom, x, y, z) values (" .. res_id .. ",'" .. k .. "'," .. v[1][1] .. "," .. v[2][1] .. "," .. v[3][1] .. ")")
         elseif update then
            rfpg.Qcur("update atom_coordinates set (x, y, z) = (" .. v[1][1] .. "," .. v[2][1] .. "," .. v[3][1] .. ") where res_id = " .. res_id .. " and atom = '" .. k .. "'")
         end
      end
   end

   for k,d in pairs(self['dihedra']) do
      d:writeDb(rfpg, res_id, update)
   end

end

--- populate residue object from database using supplied database res_id
-- @param rfpg open database handle
function Residue:dbLoad(rfpg)
   local resid = self['res_id']  -- if res_id not pre-set need to know pdb_no to find it
   --print(resid)
   local cur = rfpg.Qcur('select atom, x, y, z from atom_coordinates where res_id=' .. resid .. ' and not init')
   local atomData = cur:fetch({},'a')
   while atomData do
      local atom = geom3d.get41mtx()
      atom[1][1] = tonumber(atomData['x'])
      atom[2][1] = tonumber(atomData['y'])
      atom[3][1] = tonumber(atomData['z'])
      self['atomCoords'][ atomData['atom'] ] = atom
      --print(atomData['atom'])
      atomData = cur:fetch({},'a')
   end

   cur = rfpg.Qcur('select struc, struc2, psi, phi, omega, acc from dssp where res_id=' .. resid)
   local dsspData = cur:fetch({},'a')
   if dsspData then self['dssp'] = {} end
   while dsspData do
      self['dssp']['ss'] = dsspData['struc']
      self['dssp']['ss2'] = dsspData['struc2']
      self['dssp']['psi'] = dsspData['psi']
      self['dssp']['phi'] = dsspData['phi']
      self['dssp']['omg'] = dsspData['omega']
      self['dssp']['acc'] = dsspData['acc']
      dsspData = cur:fetch({},'a')      
   end

   cur = rfpg.Qcur('select key, angle1, angle2, dangle from dihedral where res_id=' .. resid)
   local dihedralData = cur:fetch({},'a')
   while dihedralData do
      local k = dihedralData['key']
      local d = utils.splitKey(k)  -- generates [1..4] entries for each atom
      d['dihedral1'] = tonumber(dihedralData['dangle'])
      self['dihedra'][k] = dihedron.new(d)
      --self['dihedra'][k]:dbLoad(rfpg)

      local angles = { tonumber(dihedralData['angle1']), tonumber(dihedralData['angle2']) }
      for i,aid in ipairs(angles) do
         local cur2 = rfpg.Qcur('select a.key, a.angle, b1.length as len1, b2.length as len2 from angle a, bond b1, bond b2 where a.angle_id=' .. angles[i] .. ' and a.bond1 = b1.bond_id and a.bond2 = b2.bond_id')
         local hedronData = cur2:fetch({},'a')
         while hedronData do
            --for k,v in pairs(hedronData) do print(k,v) end
            local k = hedronData['key']
            local h = utils.splitKey(k)
            local ak = utils.splitAtomKey(h[1])
            --print(ak[1],self['resn'],k)
            if ak[1] == self['resn'] then
               h['len1'] = hedronData['len1']
               h['angle2'] = hedronData['angle']
               h['len3'] = hedronData['len2']
               self['hedra'][k] = hedron.new(h)
            end
            hedronData = cur2:fetch({},'a')
         end
      end
      
      dihedralData = cur:fetch({},'a')
   end
end

      
--- generate descriptive string for Residue: 1-letter amino acid code, sequence position
-- @return descriptive string
function Residue:tostring()
   return self['res'] .. self['resn']
end

--- count entries in self['hedra'] table
-- @return number of hedra in Residue
function Residue:countHedra()
   local c = 0
   for k,v in pairs(self['hedra']) do
      c = c+1
   end
   return c
end

--- count entries in self['dihedra'] table
-- @return number of dihedra in Residue
function Residue:countDihedra()
   local c = 0
   for k,v in pairs(self['dihedra']) do
      --print('  ' .. k)
      c = c+1
   end
   return c
end

--- delete protein space atom coordinate data from this residue
function Residue:clearAtomCoords()
   self['atomCoords'] = {}
end
   
--- delete internal coordinate data from this chain
function Residue:clearInternalCoords()
   for k,h in pairs(self['hedra']) do -- 23 mar 17 : was ipairs 
      h:clearInternalCoords()
   end
   for k,d in pairs(self['dihedra']) do  -- 23 mar 17 : was ipairs 
      d:clearInternalCoords()
   end
end


--- set hedron space coordinates for each point from supplied length, angle, length
function Residue:renderHedra()
   for k,h in pairs(self['hedra']) do
      if h['updated'] then h:initPos() end
   end
end

--- set dihedron space coordinates for each point from hedra and supplied dihedral angle
function Residue:renderDihedra()
   --for k,d in utils.pairsByKeys(self['dihedra']) do
   for k,d in pairs(self['dihedra']) do
      if d['updated'] then d:initPos() end
   end
end

--- create map of key(first 3 atom ids) to dihedron for all dihedra in Residue; create backbone and sidechain tables of (dihedron, position) for atoms in PDB order
-- @see writePDB
function Residue:linkDihedra()
   local c = 0
   local k3i = {}
   --local k32i = {}
   for k,dihedron in pairs(self['dihedra']) do
      dihedron['res'] = self          -- each dihedron can find its residue
      
      local k3 = dihedron['key3']
      --print(dihedron['key'],k3)
      if not k3i[k3] then
         k3i[k3] = {}
      end
      k3i[k3][ #k3i[k3] +1 ] = dihedron       -- map to find each dihedron from atom tokens 1,2,3 

      --[[
      local k32 = dihedron['key32']
      if not k32i[k32] then
         k32i[k32] = {}
      end
      k32i[k32][ #k32i[k32] +1 ] = dihedron
      --]]
      
      for i=1,4 do                            -- PDB format ordered lists of backbone and sidechain atoms in residue
         --print(dihedron[i])
         local al = utils.splitAtomKey(dihedron[i])
         local n,r,a = tonumber(al[1]),al[2],al[3]

         --print(r,n,a, self['res'],self['resn'], dihedron[i])
         if r == self['res'] and n == self['resn'] then
            local b = backboneSort[ a ]
            local s = sidechainSort[ a ]
            if b and not self['backbone'][b] then
               self['backbone'][b] = { dihedron, i }
               --print('backbone ' .. a .. ' ' .. b .. ' '..  dihedron:tostring() .. ' ' .. i)
            elseif s and not self['sidechain'][s] then
               self['sidechain'][s] = { dihedron, i }
            elseif not (b or s) then
               assert(nil, 'cannot identify atom ' .. r .. n .. ' ' .. a .. ' dihedron ' .. dihedron:tostring() .. ' position ' .. i)
            end
         end
      end
   end

   self['key3index'] = k3i
   --self['key32index'] = k32i
end

-- based on code from DT Jones
local function genGlyCB(n4,ca4,c4)
   if not (n4 and ca4 and c4) then
      if utils.warn then io.stderr:write('genGlyCB without complete NCaC\n') end
      return
   end
         
   local n,ca,c = matrix.new(3),matrix.new(3),matrix.new(3)
   for i=1,3 do
      n[i] = n4[i][1]
      ca[i] = ca4[i][1]
      c[i] = c4[i][1]
   end
   local nca = ca - n
   local cca = ca - c
   local xx = nca + cca
   local yy = matrix.cross(nca,cca)
   local kCACBDIST = 1.538
   local kTETH_ANG = 0.9128
   local sx = kCACBDIST * math.cos(kTETH_ANG) / math.sqrt(matrix.dot(xx,xx));
   local sy = kCACBDIST * math.sin(kTETH_ANG) / math.sqrt(matrix.dot(yy,yy));
   local cb = xx * sx + yy * sy

   local cb4 = geom3d.get41mtx()
   for i=1,3 do
      cb4[i][1] = utils.setAccuracy83(ca[i] + cb[i])
   end

   return cb4
end


--- complete residue, dihedra, hedra data structures from protein space atom coordinates
function Residue:dihedraFromAtoms()

   local skbase = self['resn'] .. self['res']
   local sN, sCA, sC, sO, sCB = skbase .. 'N', skbase .. 'CA', skbase .. 'C', skbase .. 'O', skbase .. 'CB'

   --print(sN, sCA, sC, sO, sCB)
   local nkbase
   local nN, nCA, nCB

   -- atomCoords, hedra and dihedra for backbone dihedra which reach into next residue
   if self['next'] and self['next']['ordered'] then
      local rn = self['next']

      nkbase = rn['resn'] .. rn['res']
      nN, nCA, nC = nkbase .. 'N', nkbase .. 'CA', nkbase .. 'C'
      self['atomCoords'][nN] = rn['atomCoords'][nN]
      self['atomCoords'][nCA] = rn['atomCoords'][nCA]
      self['atomCoords'][nC] = rn['atomCoords'][nC]

      self['dihedra'][utils.genKey(sN,sCA,sC,nN)] = dihedron.new({ sN,sCA,sC,nN })    -- psi 
      self['dihedra'][utils.genKey(sCA,sC,nN,nCA)] = dihedron.new({ sCA,sC,nN,nCA })  -- omega i+1
      self['dihedra'][utils.genKey(sC,nN,nCA,nC)] = dihedron.new({ sC,nN,nCA,nC })    -- phi i+1

      self['hedra'][utils.genKey(sCA,sC,nN)] = hedron.new({ sCA,sC,nN })
      self['hedra'][utils.genKey(sC,nN,nCA)] = hedron.new({ sC,nN,nCA })

      rn['hedra'][utils.genKey(nN,nCA,nC)] = hedron.new({ nN,nCA,nC })
   end

   -- backbone hedra and dihedra within this residue
   self['dihedra'][utils.genKey(sN,sCA,sC,sO)] = dihedron.new({ sN,sCA,sC,sO })
   self['dihedra'][utils.genKey(sO,sC,sCA,sCB)] = dihedron.new({ sO,sC,sCA,sCB })

   local sNCaCkey = utils.genKey(sN,sCA,sC)
   if not self['hedra'][sNCaCkey] then self['hedra'][sNCaCkey] = hedron.new({ sN,sCA,sC }) end

   self['hedra'][utils.genKey(sCA,sC,sO)] = hedron.new({ sCA,sC,sO })
   self['hedra'][utils.genKey(sCB,sCA,sC)] = hedron.new({ sCB,sCA,sC })

   --if ('G' ~= self['res']) and ('A' ~= self['res']) then
   if (self['atomCoords'][skbase .. 'CG']
          or self['atomCoords'][skbase .. 'CG1']
          or self['atomCoords'][skbase .. 'OG']
          or self['atomCoords'][skbase .. 'OG1']
          or self['atomCoords'][skbase .. 'SG'] ) then
      self['hedra'][utils.genKey(sN,sCA,sCB)] = hedron.new({ sN,sCA,sCB })   -- only needed for sidechain CG residues (not gly or ala or any missing rest of side chain)
   end  

   -- terminal OXT if present
   local sOXT = skbase .. 'OXT'
   if self['atomCoords'][sOXT] then
      self['hedra'][utils.genKey(sCA,sC,sOXT)] = hedron.new({ sCA,sC,sOXT })
      self['dihedra'][utils.genKey(sN,sCA,sC,sOXT)] = dihedron.new({ sN,sCA,sC,sOXT })
   end
   
   -- sidechain hedra and dihedra
   local sct = chemdata.sidechains[self['res']]
   if sct then 
      for i,t in ipairs(sct) do
         local nt = utils.addBase(skbase,t)
         if 4 > table.getn(t) then
            self['hedra'][utils.genKey(nt[1],nt[2],nt[3])] = hedron.new(nt)
         else
            if self['atomCoords'][nt[4]] then  -- skip if missing sidechain atom at end
               self['dihedra'][utils.genKey(nt[1],nt[2],nt[3],nt[4])] = dihedron.new(nt)
            end
         end
      end
   elseif 'G' == self['res'] and not self['atomCoords'][sCB] then  -- only G,A do not have entries in chemdata.sidechains
      self['atomCoords'][sCB] = genGlyCB(self['atomCoords'][sN],self['atomCoords'][sCA],self['atomCoords'][sC])
      --print(sCB .. self['atomCoords'][sCB]:transpose():pretty())
   end

   if not self['atomCoords'][sCB] then  -- fabricate CB if not present in coordinate data
      self['atomCoords'][sCB] = genGlyCB(self['atomCoords'][sN],self['atomCoords'][sCA],self['atomCoords'][sC])
   end

   self:linkDihedra()

   --for k,d in utils.pairsByKeys(self['dihedra']) do
   for k,d in pairs(self['dihedra']) do
      d:dihedronFromAtoms(self['atomCoords'])
   end
   for k,h in pairs(self['hedra']) do  -- miss a few redundant hedra needed for some chi1 angles
      if not h['len1'] then h:hedronFromAtoms(self['atomCoords']) end
   end
      
end

--- generate ATOM records for this Residue
--<br>
-- Note: OXT defined in Dihedra like other atoms, if present
-- @todo make use of temperature factor field
-- @param chain chain ID
-- @param ndx ATOM record sequence counter
-- @return string containing sequential ATOM records, ndx
function Residue:writePDB(chain,ndx)
   --print('residue topdb '  .. chain .. ' ' .. ndx)
   local s = ''
   local dihedron
   local position
   local atomName
   local atom
   local ls
   --print('ndx',ndx)
   for i,a in pairs(self['backbone']) do
      if a then
         ndx = ndx + 1
         dihedron = a[1]
         position = a[2]
         --print('res:writePDB chain ' .. chain .. ' ndx ' .. ndx .. ' ' .. dihedron:tostring() .. ' pos ' .. position)
         atomName = utils.splitAtomKey(dihedron[position])[3]   -- dihedron['atomNames'][position]  -- dihedron['atomCoords'][coordsInitial]['names'][position]
         local akl = { dihedron[1], dihedron[2], dihedron[3], dihedron[4] } -- utils.splitKey(dihedron['key'])
         atom = self['atomCoords'][akl[position]] --  dihedron['initialCoords'][position]
         if atom then 
            ls =  utils.atomString(ndx,atomName,self['res'], chain, self['resn'], atom[1][1], atom[2][1],atom[3][1], 1.0, 0.0)
            --print(ls)
            s = s .. ls
         end
      end
   end

      for i = 1,12 do
         local a = self['sidechain'][i]
         if a then
            dihedron = a[1]
            position = a[2]
            atomName = utils.splitAtomKey(dihedron[position])[3] -- dihedron['atomNames'][position]  --  dihedron['atomCoords'][coordsInitial]['names'][position]
            if not ('CB' == atomName and 'G' == self['res']) then  -- don't put gly cbeta in pdb file
               ndx = ndx + 1
               local akl = { dihedron[1], dihedron[2], dihedron[3], dihedron[4] } -- utils.splitKey(dihedron['key'])
               atom = self['atomCoords'][akl[position]] -- dihedron['initialCoords'][position]
               if atom then
                  ls = utils.atomString(ndx,atomName,self['res'], chain, self['resn'], atom[1][1], atom[2][1], atom[3][1], 1.0, 0.0)
                  s = s .. ls
               end
            end
         end
      end
   return s,ndx
end

--- generate string in PIC format for this residue's hedra and dihedra
-- @param pdb pdb id for pic format
-- @param chn pdb chain for pic format
-- @param stats optional boolean report stats (sd, min, max) and name if available
-- @return pic format string of lines for this residue
function Residue:writeInternalCoords( pdb, chn, stats )
   local s = ''
   local base = pdb .. ' ' .. chn .. ' '
   for k,h in utils.pairsByKeys(self['hedra'], residue.keysort) do
      if h['len1'] and h['angle2'] and h['len3'] then
         s = s .. base .. h['key'] .. ' ' .. string.format('%9.5f %9.5f %9.5f', h['len1'], h['angle2'], h['len3'])
         if stats then
            if h['sd'] then
               for i,val in ipairs( { 'len1', 'angle2', 'len3' } ) do
                  --print(s,h['sd'][val], h['min'][val], h['max'][val])
                  --print()
                  s = s .. string.format(' ( ' .. val .. ' sd: %9.5f min: %9.5f max: %9.5f )', h['sd'][val], h['min'][val], h['max'][val])
               end
               s = s .. string.format(' count: %9.0f ', h['count'])
            end
         end
         s = s .. '\n'
      end
   end
   for k,d in utils.pairsByKeys(self['dihedra'], residue.keysort) do
      if d['dihedral1'] then
         s = s .. base .. d['key'] .. ' ' .. string.format('%9.5f', d['dihedral1'])
         if stats then
            if d['sd'] then
               s = s .. string.format(' ( sd: %9.5f min: %9.5f max: %9.5f count: %9.0f )', d['sd'], d['min'], d['max'], d['count'])
            end
            if d['name'] then s = s .. ' ' .. d['name'] end
         end
         s = s .. '\n'
      end
   end
   return s
end

--- create a 4x1 matrix with coordinates for specified atom token as read from DSSP
-- @param atomKey atom token of form (sequence postion)(residue)(atom string)
-- @return 4x1 matrix with protein space coordinates for atom
-- @see Chain:setStartCoords
function Residue:dsspAtom(atomKey)
   local atom = geom3d.get41mtx()
   local atomName = utils.splitAtomKey(atomKey)[3]:lower()
   if not self['dssp'] then assert(nil, 'dsspAtom ' .. atomKey .. ' no DSSP record loaded for residue ' .. self:tostring() ) end
   local datom = self['dssp'][atomName]
   if not datom then assert(nil, 'dsspAtom ' .. atomKey .. ' failed to find ' .. atomName) end

   atom[1][1] = datom['x'] 
   atom[2][1] = datom['y'] 
   atom[3][1] = datom['z']

   return atom
end

--- generate a table of atom tokens for this Residue's N, CA, C backbone atoms
-- @return table of atom tokens { (sequence postion)(residue)(atom string) } for N, CA, C backbone atoms
-- @see Chain:setStartCoords
function Residue:NCaCKeySplit()
   local rbase = self['resn'] .. self['res']
   --local key = utils.genKey(rbase .. 'N', rbase .. 'CA', rbase .. 'C')
   --return utils.splitKey(key)
   return { rbase .. 'N', rbase .. 'CA', rbase .. 'C' }
end

--- join dihedrons from N-CA-C and N-CA-CB hedrons, computing protein space coordinates for backbone and sidechain atoms
-- @param atomCoordsIn optional table of atom_token : 4x1 matrix of protein space coordinates
-- @param genSCAD boolean if true, return table of transformation matrices for each hedron key
-- @return atomCoords for residue in protein space relative to acomCoordsIn OR table of transformation matrices according to genSCAD parameter
function Residue:assemble( atomCoordsIn, genSCAD )
   --[[
   for di,d in pairs(self['dihedra']) do
      print('diheron: ' .. d['key'] .. ' angle: ' .. d['dihedral1'])
   end
   print('genSCAD',genSCAD)
   --]]
--[[
   form queue, start with n-ca-c, o-c-ca, n-ca-cb  [ o-c-ca not 2nd hedron for any dihedron and thus won't be picked up w/o adding here ]
   gen triple keys for current residue
   if no atomCoordsIn, use initial coords from generating dihedral for n-ca-c initial positions (dihedron coordinate space)

   while queue not empty
      get triple key
   for each dihedral starting with triple key (1st hedron)
         if have coordinates for all 4 atoms already
              add 2nd hedron key to back of queue
         else if have coordinates for 1st 3 atoms
              compute forward and reverse transform to take 1st 3 atoms to/from dihedron initial coordinate space
              use reverse transform to get position of 4th atom in current coordinates from dihedron initial coordinates
              add 2nd hedron key to back of queue              
         else
              ordering failed, put triple key at back of queue and hope next time we have 1st 3 atom positions (should not happen)

   loop terminates (queue drains) as triple keys which do not start any dihedra are removed without action
   
--]]

   local atomCoords = atomCoordsIn
   local transformations = {}

   local rbase = self['resn'] .. self['res']
   
   local q = deque:new()
   local NCaCKey = utils.genKey(rbase .. 'N', rbase .. 'CA', rbase .. 'C')
   
   if genSCAD then transformations[NCaCKey] = geom3d.get44mtx() end
   
   q:push_left(NCaCKey)
   q:push_left(utils.genKey(rbase .. 'O', rbase .. 'C', rbase .. 'CA'))
   q:push_left(utils.genKey(rbase .. 'N', rbase .. 'CA', rbase .. 'CB'))

   if not atomCoords then -- if list of initial coords not passed as parameter, use N-CA-C initial coords from creating dihedral
      atomCoords = {}
      local dl = self['key3index'][NCaCKey]
      for di,d in ipairs(dl) do
         local akl = { d[1], d[2], d[3], d[4] }  -- utils.splitKey( d['key'] )
         for ai,a in ipairs(akl) do
            atomCoords[a] = d['initialCoords'][ai]
         end
      end
   end
   
   while not q:is_empty() do
      local h1k = q:pop_right()
      local dihedra = self['key3index'][h1k]
      --print('start h1k ' .. h1k, dihedra)
      if dihedra then
         for di, d in ipairs(dihedra) do
            if d['initialCoords'][4] then  -- (skip incomplete dihedron if don't have 4th atom due to missing input data)
               --print('assemble: ' .. h1k .. ' -> ' .. d['key'])
               local dh2key = d['hedron2']['key']
               local akl = { d[1], d[2], d[3], d[4] }  -- utils.splitKey( d['key'] )
               if atomCoords[akl[1]]
                  and atomCoords[akl[2]]
                  and atomCoords[akl[3]]
                  and atomCoords[akl[4]]
               then
                  -- skip
                  --print('skipping already done ' .. d['key'] .. ' adding hedron ' .. dh2key)
                  q:push_left(dh2key)
               elseif atomCoords[akl[1]]
                  and atomCoords[akl[2]]
                  and atomCoords[akl[3]]
               then
                  local mt, mtr = geom3d.coordSpace( atomCoords[akl[1]], atomCoords[akl[2]], atomCoords[akl[3]], true ) -- get transforms to take 1st hedron to dihedron coordinate space and back
                  atomCoords[akl[4]] = mtr * d['initialCoords'][4]  -- apply back transform to 4th atom's dihedron space coordinates
                  if genSCAD then
                     --print(h1k, mtr:transpose():pretty())
                     transformations[h1k] = mtr
                  end
                  for i=1,3 do atomCoords[akl[4]][i][1] = utils.setAccuracy83(atomCoords[akl[4]][i][1]) end 
                  --print('finished: ' .. d['key'] .. ' adding hedron ' .. dh2key .. ' a4: ' .. akl[4] .. ' -- ' .. atomCoords[akl[4]]:transpose():pretty())
                  q:push_left(dh2key)
               else
                  if utils.warn then 
                     io.stderr:write('no coords to start ' .. d['key'] .. '\n')
                  end
                  --assert(nil,'foo')
                  --q:push_left(h1k)
               end
            end
         end
      end
   end

   --[[
   print(self:tostring())
   for k,a in pairs(atomCoords) do
      print('assemble final:',k, a:transpose():pretty())
   end
   print()
   --]]

   --for x,y in pairs(transformations) do print(x, y:transpose():pretty()) end
   
   if genSCAD then return transformations
   else return atomCoords
   end
      
end

function Residue:printInfo()
   for k,v in pairs(self['dihedra']) do v:printInfo() end
end


function Residue:writeSCADstrings()
   local s = ''
   for k,h in utils.pairsByKeys(self['hedra'], residue.keysort) do
      if h['len1'] and h['angle2'] and h['len3'] then
         s = s .. 'h_' .. h['key']:gsub(':','_') .. ' = [' .. string.format('%9.5f, %9.5f, %9.5f ];\n', h['len1'], h['angle2'], h['len3'])
      end
   end
   for k,d in utils.pairsByKeys(self['dihedra'], residue.keysort) do
      if d['dihedral1'] then
         s = s .. 'd_' .. d['key']:gsub(':','_') .. ' = [ ' .. string.format('%9.5f ];\n', d['dihedral1'])
      end
   end
   return s
end

function Residue:writeSCADhedra()
   local s = ''
   for k,h in utils.pairsByKeys(self['hedra'], residue.keysort) do
      if h['len1'] and h['angle2'] and h['len3'] then
         s = s .. 'h_' .. h['key']:gsub(':','_') .. ' = [' .. string.format('%9.5f, %9.5f, %9.5f ];\n', h['len1'], h['angle2'], h['len3'])
      end
   end
   return s
end

function Residue:writeSCADassembly()
end


return residue
