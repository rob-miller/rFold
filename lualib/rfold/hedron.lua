--[[
   hedron.lua
   
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
]]

--- hedron structure class for rFold
-- 
-- a hedron consists of 3 coplanar objects, specified as two bond lengths and the angle between them. 
--
-- @classmod Hedron

--

local utils = require 'rfold.utils'
local geom3d = require 'rfold.geom3d'

local hedron = {} -- module

local Hedron = {}  -- class table



------------------------------------------------------------------------------------------------
-- hedron
------------------------------------------------------------------------------------------------
-- class properties:

--- 3 atom tokens separated by ':' identifying this hedron : **4SCA:4SC:4SO**
-- @field key string

--- distance between first and second atoms
-- @field len1 float

--- angle formed between 3 atoms
-- @field angle2 float

--- distance between second and third atoms
-- @field len3 float

--- three 4x1 matrices specifying coordinates of constituent atoms, initially atom3 on +Z axis
-- @table atoms 

--- three 4x1 matrices specifying coordinates of constituent atoms, reversed order - initially atom1 on +Z axis
-- @table atomsR

--- flag indicting that atom coordinates are up to date (do not need to be recalculated from len1-angle2-len3)
-- @field updated


--- initialiser for hedron class (not a class method).
-- <br><br>
-- The input object 'o' is a table with expected fields:
-- <br> [1..3] = atom string tokens for 3 bonded atoms forming plane<br> 'len1' = atom1 to atom2 distance<br> 'angle2' = angle formed by 3 atoms<br> 'len3' = atom2 to atom3 distance
-- <br>
-- @param o table as defined above
-- @return initialised hedron object
function hedron.new (o)
   --o = o or {}
   assert(o and o[3],'attempt to instantiate hedron without atom 3')
   setmetatable(o, { __index = Hedron })

   o['key'] = utils.genKey(o[1],o[2],o[3])
   o['updated'] = true

   if not o['atoms'] then
      o['atoms'] = {}
   end
   if not o['atomsR'] then
      o['atomsR'] = {}
   end

   --print('instantiate hedron ' .. o['key'])
   
   return o
end


--- generate descriptive string for Hedron. 3- followed by 'key' for this Hedron = 3 atom tokens separated by :'s
-- @return descriptive string
function Hedron:tostring()
   local s = '3-[' .. self['key'] .. ']'
   if (self['len1'] and self['angle2'] and self['len3']) then s = s .. ' ' .. self['len1'] .. ' ' .. self['angle2'] .. ' ' .. self['len3'] end
   return  s
end

--- generate hedron space coordinates for 3 atoms. with specified bond lengths and angle between on XZ plane 
-- (Y=0 for all atoms).
-- <br>initialises [atoms] and [atomsR].
-- <br>Hedron coordinate system: a2 = 0,0,0 ; a3 on Z-axis ; a1 on XZ plane (-Z for angle >90).
-- <br>reverse coordinates swap : a2 = 0,0,0 ; a1 on Z-axis ; a3 on XZ plane.
function Hedron:initPos()

   if not (self['len1'] and self['angle2'] and self['len3']) then
      if utils.warn then
         io.stderr:write('incomplete hedron missing ')
         if not self['len1'] then io.stderr:write('len1 ') end
         if not self['angle2'] then io.stderr:write('angle2 ') end
         if not self['len3'] then io.stderr:write('len3 ') end
         io.stderr:write(self:tostring() .. '\n')
      end
      self['updated'] = false
      return
   end
   self['atoms'][1] = geom3d.get41mtx()   
   self['atoms'][2] = geom3d.get41mtx()    -- a2 to 0,0,0
   self['atoms'][3] = geom3d.get41mtx()

   --local ar = math.rad( self['angle2'] )
   local sar = math.rad(180.0 - self['angle2'])    -- supplement: angles which add to 180 are supplementary

   self['atoms'][3][3][1] = self['len3']   -- a3 is len3 up from a2 on Z axis, X=Y=0
   self['atoms'][1][1][1] = math.sin(sar) * self['len1']      -- a1 X is sin( sar ) * len1
   self['atoms'][1][3][1] = - (math.cos(sar) * self['len1'])  -- a1 Z is -(cos( sar ) * len1) (assume angle2 always obtuse, so a1 is in -Z

   -- same again but a3 at 0,0,0
   self['atomsR'][1] = geom3d.get41mtx()   
   self['atomsR'][2] = geom3d.get41mtx()   -- a2 to 0,0,0
   self['atomsR'][3] = geom3d.get41mtx()    

   self['atomsR'][1][3][1] = self['len1']   -- a1r is len1 up from a2 on Z axis, X=Y=0
   self['atomsR'][3][1][1] = math.sin(sar) * self['len3']    -- a3r X is sin( sar ) * len3
   self['atomsR'][3][3][1] = - (math.cos(sar) * self['len3'])  -- a3r Z is -(cos( sar ) * len3)
   
   self['updated'] = false
end

local function getDistAngleDist(a1,a2,a3)
   local a1a2 = geom3d.getDistance3d(a1,a2)
   local a2a3 = geom3d.getDistance3d(a2,a3)

   local a1a3 = geom3d.getDistance3d(a1,a3)

   local a1a2a3 = math.deg( geom3d.getAngleS3(a1a2,a2a3,a1a3) )

   return a1a2,a1a2a3,a2a3
end

function Hedron:hedronFromAtoms(atomCoords)
   local a1,a2,a3 = atomCoords[self[1]],atomCoords[self[2]],atomCoords[self[3]]
   if not (a1 and a2 and a3) then
      if utils.warn then
         io.stderr:write('hedron missing coordinates for ')
         if not a1 then io.stderr:write('a1 ') end
         if not a2 then io.stderr:write('a2 ') end
         if not a3 then io.stderr:write('a3 ') end
         io.stderr:write(self:tostring() .. '\n')
      end
      return
   end
   local len1,angle2,len3 = getDistAngleDist(a1,a2,a3)
   self['len1'],self['angle2'],self['len3'] = utils.setAccuracy95(len1),utils.setAccuracy95(angle2),utils.setAccuracy95(len3)
end

function Hedron:clearInternalCoords()
   self['len1'] = nil
   self['angle2'] = nil
   self['len3'] = nil
end

--- write hedron data to rfold database
-- @param rfpg open database handle
-- @param res_id residue ID in residue table
-- @param did dihedral id
-- @param update optional flag, if false silently skip if entry exists already in dihedral / angle / bond tables
-- @return angleID identifier for new or updated entry in angle table
function Hedron:writeDb(rfpg, res_id, did, update)
   local akl = utils.splitKey(self['key'])
   local al = {}
   for i,ak in ipairs(akl) do
      local a = utils.splitAtomKey(ak)
      al[i] = a[2]..a[3]
   end
   local as = rfpg.pgArrayOfStrings(al[1],al[2],al[3])
   local acid = rfpg.Q("select id from angle_string where atom_string='" .. as .. "'")[1]
   local tst = rfpg.Q('select angle_id from angle where res_id = ' .. res_id .. ' and angle_class = ' .. acid)
   local aid
   if not tst then
      aid = rfpg.Q('insert into angle (res_id, dihedral_id, angle_class, angle) values (' .. res_id .. ',' .. did .. ',' .. acid .. ',' .. self['angle2'] .. ') returning angle_id')[1]
   elseif update then
      aid = tst[1]
      rfpg.Qcur('update angle set angle = ' .. self['angle2'] .. ' where angle_id = ' .. aid)
   end

   if aid then
      local bids={}
      for i=1,2 do
         local bas = rfpg.pgArrayOfStrings(al[i],al[i+1])
         local bcid = rfpg.Q("select id from bond_string where atom_string='" .. bas .. "'")[1]
         tst = rfpg.Q('select bond_id from bond where res_id=' .. res_id .. ' and angle_id=' .. aid .. ' and bond_class=' .. bcid)
         if not tst then
            bids[i] = rfpg.Q('insert into bond (angle_id, bond_class, res_id, length) values (' .. aid .. ',' .. bcid .. ',' .. res_id .. ',' .. (1==i and self['len1'] or self['len3']) .. ') returning bond_id')[1]
         elseif update then
            bids[i] = tst[1]
            rfpg.Qcur('update bond set length=' .. (1==i and self['len1'] or self['len3']) .. ' where bond_id=' .. bids[i])
         end
      end
      rfpg.Qcur('update angle set bond1=' .. bids[1] .. ', bond2=' .. bids[2] .. ' where angle_id=' .. aid)
   end
   return aid 
end

return hedron
