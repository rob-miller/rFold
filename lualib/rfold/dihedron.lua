--[[
   dihedron.lua
   
Copyright 2016,2017 Robert T. Miller

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

--- dihedron structure class for rFold
--
-- a dihedron consists of two faces, or hedra
-- <br>**external dependencies** : deque (https://github.com/catwell/cw-lua/tree/master/deque)
--
-- @classmod Dihedron

local utils = require 'rfold.utils'
local geom3d = require 'rfold.geom3d'

local dihedron = {} -- module

local Dihedron = {}  -- class table


------------------------------------------------------------------------------------------------
-- dihedron
------------------------------------------------------------------------------------------------
-- class properties:

---  4 atom tokens separated by ':' identifying this dihedron : **6FCA:6FCB:6FCG:6FCD2**
-- @field key string

--- 3 atom tokens separated by ':' identifying first hedron : **6FCA:6FCB:6FCG**
-- @field key3 string

--- 3 atom tokens separated by ':' identifying second hedron : **6FCB:6FCG:6FCD2**
-- @field key32 string

--- dihedral angle for this dihedron
-- @field dihedral1 float

--- flag indicting that atom coordinates are up to date (do not need to be recalculated from dihedral1)
-- @field updated

--- Residue object which includes this dihedron; set by Residue:linkDihedra()
-- @table res

--- first Hedron object; set by self:setHedra()
-- @table hedron1

--- second Hedron object; set by self:setHedra()
-- @table hedron2

--- four 4x1 matrices holding atom coordinates comprising this dihedron
-- @table initialCoords



--- Dihedron class object initialiser (not a class method).
-- <br>
-- The input object 'o' is a table with expected fields: <br>
-- <br> [1..4] = atom tokens for 4 bonded atoms forming dihedral angle<br> 'dihedral1' = dihedral angle formed by 4 atoms
-- <br>
-- @param o table as described above
-- @return minimally initialised Dihedron object
function dihedron.new (o)
   --o = o or {}
   assert(o and o[4],'attempt to instantiate dihedron without atom 4')

   setmetatable(o, { __index = Dihedron })

   o['key'] = utils.genKey(o[1],o[2],o[3],o[4])
   o['key3'] = utils.genKey(o[1],o[2],o[3])
   o['key32'] = utils.genKey(o[2],o[3],o[4])
   
   o['updated'] = true

   if not o['initialCoords'] then
      o['initialCoords'] = {}
   end

   --print('instantiate dihedron ' .. o['key'])
   return o
end


--- generate descriptive string for Dihedron: 4- followed by 'key' for this Dihedron = 4 atom tokens separated by :'s
-- @return descriptive string
function Dihedron:tostring()
   local s = '4-[' .. self['key'] .. ']'
   if self['dihedral1'] then s = s  .. ' ' .. self['dihedral1'] end
   return s
end

local function getHedron(res,key)
   local hedron = res['hedra'][key]
   if not hedron and res['prev'] and res['prev']['ordered'] then hedron = res['prev']['hedra'][key] end
   if not hedron and res['next'] and res['next']['ordered'] then hedron = res['next']['hedra'][key] end
   return hedron
end

--- determine forward or reverse hedra keys based on 4 atom tokens for this dihedron, set object variables for hkey and hedra
-- @return true if hedra keys are reverse ordered from dihedron atom tokens
function Dihedron:setHedra()
   --print('sh: ' .. self:tostring())
   local reverse = false
   local res = self['res']
   local h1key = utils.genKey(self[1], self[2], self[3])
   local hedron1 = getHedron(res,h1key)
   local h2key
   local hedron2
   if not hedron1 then
      reverse=true
      h1key = utils.genKey(self[3], self[2], self[1])
      hedron1 = getHedron(res,h1key)
      h2key = utils.genKey(self[4], self[3], self[2])
   else
      h2key = utils.genKey(self[2], self[3], self[4])
   end

   if not hedron1 then
      --assert(nil,'residue ' .. res['res'] .. res['resn'] .. ' failed to locate 1st hedron ' .. h1key .. ' a4= ' .. self[4])
      io.stderr:write('residue ' .. res['res'] .. res['resn'] .. ' failed to locate 1st hedron ' .. h1key .. ' a4= ' .. self[4] .. '\n')
   end

   hedron2 = getHedron(res,h2key)
   if not hedron2 then 
      --assert(nil,'residue ' .. res['res'] .. res['resn'] .. ' failed to locate 2nd hedron ' .. h2key .. ' a1= ' .. self[1])
      io.stderr:write('residue ' .. res['res'] .. res['resn'] .. ' failed to locate 2nd hedron ' .. h2key .. ' a1= ' .. self[1] .. '\n')
   end

   self['hedron1'] = hedron1
   self['h1key'] = h1key
   self['hedron2'] = hedron2
   self['h2key'] = h2key

   return reverse
end

--- generate dihedron space coordinates for 4 atoms with specified dihedral, first 3 on XZ plane (a1 in -Z), 4th in +Z and rotated
--<br>
-- coordinates stored in 'initialCoords' field, 'updated' set to false
function Dihedron:initPos()

   local reverse = self:setHedra()
   local hedron1 = self['hedron1']
   local hedron2 = self['hedron2']

   if not (hedron1 and hedron2) then return end -- error bad PDB data!
   
   local complete = true
   for i=1,3 do
      if not hedron1['atoms'][i] then complete=false end
   end
   for i=1,3 do
      if not hedron2['atoms'][i] then complete=false end
   end
   if not complete then
      if utils.warn then
         io.stderr:write('dihedron: hedra missing atoms ' .. self:tostring() .. (reverse and ' reverse ' or ' forward ') .. hedron1:tostring() .. ' ' .. hedron2:tostring() .. '\n')
      end
      return
   end
      
   local initial = {}
   local a4preRotation
   local a4shift
   
   if not reverse then
      initial[1] = hedron1['atoms'][1]:copy()
      initial[2] = hedron1['atoms'][2]:copy()
      initial[3] = hedron1['atoms'][3]:copy()
      
      a4preRotation = hedron2['atomsR'][3]:copy()
      a4shift = hedron2['len1']
   else
      initial[1] = hedron1['atomsR'][3]:copy()
      initial[2] = hedron1['atomsR'][2]:copy()
      initial[3] = hedron1['atomsR'][1]:copy()
      
      a4preRotation = hedron2['atoms'][1]:copy()
      a4shift = hedron2['len3']
   end
   
   a4preRotation[3][1] = a4preRotation[3][1] * -1                   -- a4 to +Z
   a4preRotation[3][1] = a4preRotation[3][1] + a4shift      -- hedron2 shift up so a2 at 0,0,0
   local mrz = geom3d.genMrz( math.rad( self['dihedral1'] ) )
   initial[4] = mrz * a4preRotation                              -- dihedral set

   --[[
   initial['names']={}
   for i=1,4 do
      initial['names'][i] = getAtomName(self[i])
   end
   --]]
   
   self['initialCoords'] = initial

   self['a4preRotation'] = a4preRotation

   --print(hedron1:tostring() .. ' ' .. hedron2:tostring() .. ' ' .. self:tostring())
   --print('dip : ' .. ' [' ..  self['initialCoords'][1]:transpose():pretty() .. '][' ..  self['initialCoords'][2]:transpose():pretty() .. '][' ..  self['initialCoords'][3]:transpose():pretty() .. '][' ..  self['initialCoords'][4]:transpose():pretty() .. ']')
   --os.exit()

   self['updated'] = false
end


local function getDistAngleDistAngleDist(a1,a2,a3,a4)
   local a1a2 = geom3d.getDistance3d(a1,a2)
   local a2a3 = geom3d.getDistance3d(a2,a3)
   local a3a4 = geom3d.getDistance3d(a3,a4)

   local a1a3 = geom3d.getDistance3d(a1,a3)
   local a2a4 = geom3d.getDistance3d(a2,a4)

   local a1a2a3 = math.deg( geom3d.getAngleS3(a1a2,a2a3,a1a3) )
   local a2a3a4 = math.deg( geom3d.getAngleS3(a2a3,a3a4,a2a4) )

   return a1a2,a1a2a3,a2a3,a2a3a4,a3a4
end

--- generate self[dihedral1] and self[initialCoords], populate hedra as needed from atomCoords
-- @param atomCoords table of atom token, atom matrix of position in protein coordinate space
function Dihedron:dihedronFromAtoms(atomCoords)
   local reverse = self:setHedra()
   local hedron1 = self['hedron1']
   local hedron2 = self['hedron2']

   if not (hedron1 and hedron2) then return end -- error bad PDB data!
   
   --local initial = {}
   --local a4preRotation
   
   local a1, a2, a3, a4 = atomCoords[self[1]], atomCoords[self[2]], atomCoords[self[3]], atomCoords[self[4]]
   if not (a1 and a2 and a3 and a4) then
      if utils.warn then
         io.stderr:write('dihedron missing coordinates for ')
         for i=1,4 do
            if not atomCoords[self[i]] then io.stderr:write('a'..i..' ') end
         end
         io.stderr:write('  : ' .. self:tostring() .. ' ' .. hedron1:tostring() .. ' ' .. hedron2:tostring() .. '\n')
      end
      return
   end
   local mt = geom3d.coordSpace( a1, a2, a3 )
   local da4 = mt * a4
   self['dihedral1'] = utils.setAccuracy95( math.deg( math.atan2(da4[2][1],da4[1][1]) ) )

   local a1a2, a1a2a3, a2a3, a2a3a4, a3a4 = getDistAngleDistAngleDist(a1,a2,a3,a4)
   
   if not reverse then
      hedron1['len1'] = utils.setAccuracy95(a1a2)
      hedron1['len3'] = utils.setAccuracy95(a2a3)
      
      hedron2['len1'] = utils.setAccuracy95(a2a3)
      hedron2['len3'] = utils.setAccuracy95(a3a4)
   else
      hedron1['len3'] = utils.setAccuracy95(a1a2)
      hedron1['len1'] = utils.setAccuracy95(a2a3)
      
      hedron2['len3'] = utils.setAccuracy95(a2a3)
      hedron2['len1'] = utils.setAccuracy95(a3a4)      
   end

   hedron1['angle2'] = utils.setAccuracy95(a1a2a3)
   hedron2['angle2'] = utils.setAccuracy95(a2a3a4)

   hedron1['updated'] = true
   hedron2['updated'] = true

   self['updated'] = true

   --print(hedron1:tostring() .. ' ' .. hedron2:tostring() .. ' ' .. self:tostring())
   
end

--- delete dihedral angle value for this dihedron, keep key
function Dihedron:clearInternalCoords()
   self['dihedral1'] = nil
end

--- write dihedron data to rfold database
-- @param rfpg open database handle
-- @param res_id residue ID in db residue table
-- @param update optional flag, if false silently skip if entry exists already in dihedral / angle / bond tables 
function Dihedron:writeDb(rfpg, res_id, update)
   --print(self:tostring())
   if not self['dihedral1'] then return end
   local akl = utils.splitKey(self['key'])
   local al = {}
   for i,ak in ipairs(akl) do
      local a = utils.splitAtomKey(ak)
      al[i] = a[2]..a[3]  -- just residue and atom, not residue position
   end
   local as = rfpg.pgArrayOfStrings(al[1],al[2],al[3],al[4])
   local dcid = rfpg.Q("select id from dihedral_string where atom_string='" .. as .. "'")[1]
   local tst = rfpg.Q('select dihedral_id from dihedral where res_id = ' .. res_id .. ' and dihedral_class = ' .. dcid)
   local did
   if not tst then
      did = rfpg.Q('insert into dihedral (res_id, key, dihedral_class, dangle) values (' .. res_id .. ",'" .. self['key'] .. "'," .. dcid .. ',' .. self['dihedral1'] .. ') returning dihedral_id')[1]
   elseif update then
      did = tst[1]
      rfpg.Qcur('update dihedral set dangle = ' .. self['dihedral1'] .. ' where dihedral_id = ' .. did)
   end

   if did then
      if not self['hedron1'] then self:setHedra() end
      local aid1 = self['hedron1']:writeDb(rfpg, res_id, update)
      local aid2 = self['hedron2']:writeDb(rfpg, res_id, update)
      rfpg.Qcur('update dihedral set angle1=' .. aid1 .. ', angle2=' .. aid2 .. ' where dihedral_id=' .. did)
   end
   
end



function Dihedron:printInfo()
   print(self:tostring())
   for k,v in pairs(self) do print(k,v) end
end

return dihedron
