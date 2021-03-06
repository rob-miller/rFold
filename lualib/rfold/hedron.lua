--[[
   hedron.lua
   
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

--- three 4x1 matrices specifying hedron space coordinates of constituent atoms, initially atom3 on +Z axis
-- @table atoms 

--- three 4x1 matrices specifying hedron space coordinates of constituent atoms, reversed order - initially atom1 on +Z axis
-- @table atomsR

--- flag indicating that atom coordinates are up to date (do not need to be recalculated from len1-angle2-len3)
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

   local classArr={}
   for i=1,3 do classArr[i] = utils.splitAtomKey(o[i])[3] end
   o['class'] = table.concat(classArr)


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

--- get average, std_dev (sample), min, max for len1, angle2, len3 from rfold database (not a class method)
-- calculated over range: total_average of angle for selected subset +/- 1 standard deviation (upper and lower bounds)
-- @param rfpg open database handle
-- @param atom_selector list of atom substrings using postgresql LIKE pattern matching, e.g. { '_N','_CA', '_C' } or { 'GN', 'GCA', 'GC' }
-- @param residue_selector optional rFold database query returning res_id, e.g. "select res_id from dssp where struc='H' and struc2=' X S+   '" limits to residues inside alpha helices
-- @return ( len1 = ( avg, sd, min, max), angle2 =  ( avg, sd, min, max), len3 = ( avg, sd, min, max) )
function hedron.avgDb(rfpg, atom_selector, residue_selector)
   local qry = 'with '

   -- residue subset:
   if residue_selector then
      qry = qry .. 'residue_subset as ( ' .. residue_selector .. ' ), '
   end

   -- atom selector subset:
   qry = qry .. "angle_classes as ( select id from angle_string where atom_string[1] like '" .. atom_selector[1] .. "' and atom_string[2] like '" .. atom_selector[2] .. "' and atom_string[3] like '" .. atom_selector[3] .. "' ), "

   -- initial working subset of all angle, bond1_len, bond2_len for residue and atom string selectors :
   qry = qry .. 'angle_subset as ( select a.angle as angle, b1.length as b1L, b2.length as b2L from angle a, bond b1, bond b2 where a.angle_class in (select id from angle_classes) '
   if residue_selector then
      qry = qry .. 'and a.res_id in (select res_id from residue_subset) '
   end
   qry = qry .. 'and a.bond1=b1.bond_id and a.bond2=b2.bond_id ), '

   -- +/- 1 standard deviation bounds on angle:
   qry = qry .. 'bounds as ( select (avg(angle) - stddev_samp(angle) * 1) as lower_bound, (avg(angle) + stddev_samp(angle) * 1) as upper_bound from angle_subset ) '

   -- columns to return:
   qry = qry .. 'select '
   for i,col in ipairs( { 'b1L', 'angle', 'b2L' } ) do
      for j,fn in ipairs( { 'avg', 'stddev_samp', 'min', 'max' } ) do
         qry = qry .. fn .. '(' .. col .. ') as ' .. fn .. '_' .. col .. ', '
      end
   end
   --- trim last ', '
   --qry = qry:sub(1, -3)
   qry = qry .. 'count(angle) as count_angle '

   -- specify upper/lower bounds subset to calculate over:
   qry = qry .. ' from angle_subset where angle between (select lower_bound from bounds) and (select upper_bound from bounds)'

   --print(qry)
   local rslt = rfpg.Q(qry)

   --[[
   for x,y in ipairs(rslt) do
      print(x,y, y .. ' ')
      if not y then print(x,qry)
         os.exit()
      end
   end
   --]]
   
   return {
      ['len1'] = { rslt[1], rslt[2], rslt[3], rslt[4] },
      ['angle2'] = { rslt[5], rslt[6], rslt[7], rslt[8] },
      ['len3'] = { rslt[9], rslt[10], rslt[11], rslt[12] },
      ['count'] = rslt[13]
   }
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
-- @param res_id residue id
-- @param update optional flag, if false skip with warning if entry exists already in dihedral / angle / bond tables
-- @return angleID identifier for new or updated entry in angle table
function Hedron:writeDb(rfpg, res_id, update)
   --print(self['key'])
   local akl = utils.splitKey(self['key'])
   local al = {}
   for i,ak in ipairs(akl) do
      local a = utils.splitAtomKey(ak)
      al[i] = a[2]..a[3]
   end
   local as = rfpg.pgArrayOfStrings(al[1],al[2],al[3])
   local acid = rfpg.Q("select id from angle_string where atom_string='" .. as .. "'")[1]
   local tst = rfpg.Q('select angle_id from angle where res_id = ' .. res_id .. " and key = '" .. self['key'] .. "' and angle_class = " .. acid .. ' and angle=' .. self['angle2'])
   local aid
   if tst then
      aid = tst[1]   -- already set as we want it, no update
   else
      tst = rfpg.Q('select angle_id from angle where res_id = ' .. res_id .. " and key = '" .. self['key'] .. "' and angle_class = " .. acid )
      if not tst then
         aid = rfpg.Q('insert into angle (res_id, key, angle_class, angle) values (' .. res_id .. ",'" .. self['key'] .. "'," .. acid .. ',' .. self['angle2'] .. ') returning angle_id')[1]
      elseif update then
         aid = tst[1]
         rfpg.Qcur('update angle set angle = ' .. self['angle2'] .. ' where angle_id = ' .. aid)
      else
         aid = tst[1]
         local v = rfpg.Q('select value from angle where angle_id=' .. aid)[1]
         io.stderr('res_id ' .. res_id .. ' hedra key ' .. self['key'] .. ' new angle= ' .. self['angle2'] .. ' already in database with different value= ' .. v .. ' angle_id= ' .. angle_id .. '\n')
      end
   end

   local bids={}
   for i=1,2 do
      local bas = rfpg.pgArrayOfStrings(al[i],al[i+1])
      local bcid = rfpg.Q("select id from bond_string where atom_string='" .. bas .. "'")[1]
      local len = (1==i and self['len1'] or self['len3'])
      
      tst = rfpg.Q('select bond_id from bond where angle_id=' .. aid .. ' and bond_class=' .. bcid .. ' and length =' .. len )
      if tst then
         bids[i] = tst[1]
      else
         tst = rfpg.Q('select bond_id from bond where angle_id=' .. aid .. ' and bond_class=' .. bcid )
         if not tst then
            bids[i] = rfpg.Q('insert into bond (angle_id, bond_class, length) values (' .. aid .. ',' .. bcid .. ',' .. len .. ') returning bond_id')[1]
         elseif update then
            bids[i] = tst[1]
            rfpg.Qcur('update bond set length=' .. len .. ' where bond_id=' .. bids[i])
         else
            bids[i] = tst[1]
            local v = rfpg.Q('select length from bond where bond_id=' .. bids[i])[1]
            io.stderr:write('res_id ' .. res_id .. ' hedra key ' .. self['key'] .. ' bond string= ' .. bas .. ' new len= ' .. len .. ' already in database with different value= ' .. v .. ' bond_id= ' .. bids[i] .. '\n')
            os.exit()
         end
      end
   end
   rfpg.Qcur('update angle set bond1=' .. bids[1] .. ', bond2=' .. bids[2] .. ' where angle_id=' .. aid)

   return aid 
end

--- return val if check is true, else return 0 
local function or0(check,val)
   return val
   --return (check and val or 0)
end

--- populate len1, angle2, len3 values with average results using rFold db residue selection query and atom keys in string (residue will be '_' for postgresql 'like' query if appropriate)
-- add tags sd[len1, angle2, len3], min[...], max[...]
-- @param rfpg open database handle
-- @param resSelector rFold database query returning res_id, e.g. "select res_id from dssp where struc='H' and struc2=' X S+   '" limits to residues inside alpha helices
function Hedron:getDbStats(rfpg, resSelector)
   --print('enter')
   local akl = utils.splitKey(self['key'])
   local al = {}
   for i,ak in ipairs(akl) do
      local a = utils.splitAtomKey(ak)
      al[i] = a[2]..a[3]
   end
   local stats = hedron.avgDb(rfpg,al,resSelector)
   --local c = 

   self['count'] = stats['count']

   -- if no data (c=0) for this angle in database for some reason, return 0's instead of nil's so don't crash elsewhere
   
   self['len1'] = stats['len1'][1]
   self['angle2'] = stats['angle2'][1]
   self['len3'] = stats['len3'][1]

   --[[
   for i,fn in ipairs( {'sd', 'min', 'max'} ) do
      self[fn] = {}
      for j,val in ipairs( {'len1', 'angle2', 'len3'} ) do
         self[fn][val] = stats[val][i+1]
         --print(fn,val,stats[val][i+1])
         --print(fn .. ' ' .. val .. ' ' .. self[fn][val] .. ' ' ..  stats[val][i+1])
      -- heisenbug here : works if the print() below is uncommented, fails repeatably otherwise
         --print()
      end
   end
   --]]

   for i,fn in ipairs( {'sd', 'min', 'max'} ) do
      self[fn] = {}
   end
   for j,val in ipairs( {'len1', 'angle2', 'len3'} ) do
      local stats2 = stats[val]
      for i,fn in ipairs( {'sd', 'min', 'max'} ) do
         self[fn][val] = stats2[i+1]
         --print(fn,val,stats[val][i+1])
         --print(fn .. ' ' .. val .. ' ' .. self[fn][val] .. ' ' ..  stats[val][i+1])
         --print()
      end
   end
   
   
   --[[
   print(collectgarbage("count"))
   print(self['len1'],self['angle2'],self['len3'])

   for i,fn in ipairs( {'sd', 'min', 'max'} ) do
      for i,val in ipairs( { 'len1', 'angle2', 'len3' } ) do
         print(val,self['sd'][val],self['min'][val],self['max'][val])
         print(val .. ' ' .. self['sd'][val] .. ' ' .. self['min'][val] .. ' ' .. self['max'][val])
      end
   end
   --]]
end



return hedron
