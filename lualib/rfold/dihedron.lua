
--- dihedron structure class for rFold
--
-- **external dependencies** : deque (https://github.com/catwell/cw-lua/tree/master/deque)
--
-- @classmod Dihedron

local utils = require 'rfold.utils'
local geom3d = require 'rfold.geom3d'

local dihedron = {} -- module

local Dihedron = {}  -- class table


------------------------------------------------------------------------------------------------
-- dihedron
------------------------------------------------------------------------------------------------

--- Dihedron class object initialiser (not a class method)
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
      assert(nil,'residue ' .. res['res'] .. res['resn'] .. ' failed to locate 1st hedron ' .. h1key .. ' a4= ' .. self[4])
   end

   hedron2 = getHedron(res,h2key)
   if not hedron2 then 
      assert(nil,'residue ' .. res['res'] .. res['resn'] .. ' failed to locate 2nd hedron ' .. h2key .. ' a1= ' .. self[1])
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

   local complete = true
   for i=1,3 do
      if not hedron1['atoms'][i] then complete=false end
   end
   for i=1,3 do
      if not hedron2['atoms'][i] then complete=false end
   end
   if not complete then
      io.stderr:write('dihedron: hedra missing atoms ' .. self:tostring() .. (reverse and ' reverse ' or ' forward ') .. hedron1:tostring() .. ' ' .. hedron2:tostring() .. '\n')
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

   --local initial = {}
   --local a4preRotation
   
   local a1, a2, a3, a4 = atomCoords[self[1]], atomCoords[self[2]], atomCoords[self[3]], atomCoords[self[4]]
   if not a4 then
      io.stderr:write('dihedron missing coordinates for a4: ' .. self:tostring() .. ' ' .. hedron1:tostring() .. ' ' .. hedron2:tostring() .. '\n')
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

function Dihedron:clearInternalCoords()
   self['dihedral1'] = nil
end

return dihedron
