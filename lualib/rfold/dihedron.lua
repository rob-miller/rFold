
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
   return '4-[' .. self['key'] .. ']'
end

--- generate dihedron space coordinates for 4 atoms with specified dihedral, first 3 on XZ plane (a1 in -Z), 4th in +Z and rotated
--<br>
-- coordinates stored in 'initialCoords' field, 'updated' set to false
function Dihedron:initPos()
   local res = self['res']

   local h1key = utils.genKey(self[1], self[2], self[3])
   --print(res:tostring() .. ' ' .. self:tostring() .. ' searching 1: ' .. h1key)
   
   local hedron1 = res['hedra'][h1key]
   if not hedron1 and res['prev'] then hedron1 = res['prev']['hedra'][h1key] end
   if not hedron1 and res['next'] then hedron1 = res['next']['hedra'][h1key] end
   if not hedron1 then return self:initPosR() end
   self['hedron1'] = hedron1
   self['h1key'] = h1key
   --print(res:tostring() .. ' ' .. self:tostring() .. ' found 1: ' .. hedron1:tostring())
   
   local h2key = utils.genKey(self[2], self[3], self[4])
   --print(res:tostring() .. ' ' .. self:tostring() .. ' searching 2: ' .. h2key)

   local hedron2 = res['hedra'][h2key]
   if not hedron2 and res['prev'] then hedron2 = res['prev']['hedra'][h2key] end
   if not hedron2 and res['next'] then hedron2 = res['next']['hedra'][h2key] end
   if not hedron2 then 
      assert(nil,'residue ' .. res['res'] .. res['resn'] .. ' failed to locate 2nd hedron ' .. h2key .. ' a1= ' .. self[1])
   end
   self['hedron2'] = hedron2
   self['h2key'] = h2key
   --print(res:tostring() .. ' ' .. self:tostring() .. ' found 2: ' .. hedron2:tostring())

   local initial = {}
   initial[1] = hedron1['atoms'][1]:copy()
   initial[2] = hedron1['atoms'][2]:copy()
   initial[3] = hedron1['atoms'][3]:copy()

   local a4preRotation = hedron2['atomsR'][3]:copy()
   
   a4preRotation[3][1] = a4preRotation[3][1] * -1                   -- a4 to +Z
   a4preRotation[3][1] = a4preRotation[3][1] + hedron2['len1']      -- hedron2 shift up so a2 at 0,0,0
   
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
   
   -- print('qip : ' .. utils.genKey(self['a1'], self['a2'], self['a3'],self['a4']) .. ' ' .. h1key .. ' ' .. h2key .. ' [' ..  self['a1initial']:transpose():pretty() .. '][' ..  self['a2initial']:transpose():pretty() .. '][' ..  self['a3initial']:transpose():pretty() .. '][' ..  self['a4initial']:transpose():pretty() .. ']')

   self['updated'] = false
end

--- called by initPos if hedron2 is reversed
-- <br>
-- generate coordinates for 4 atoms with specified dihedral, first 3 on XZ plane with Z>=0, 4th in +Z -- hedrons in reverse order from above
-- coordinates stored in 'initialCoords' field, 'updated' set to false
function Dihedron:initPosR()
   
   local res = self['res']

   local h1key = utils.genKey(self[3], self[2], self[1])
   --print(res:tostring() .. ' ' .. self:tostring() .. ' searching reverse 1: ' .. h1key)
   local hedron1 = res['hedra'][h1key]
   if not hedron1 and res['prev'] then hedron1 = res['prev']['hedra'][h1key] end
   if not hedron1 and res['next'] then hedron1 = res['next']['hedra'][h1key] end
   if not hedron1 then
      assert(nil,'residue ' .. res['res'] .. res['resn'] .. ' failed to locate 1st hedron (R) ' .. h1key .. ' a4= ' .. self['a4'])
   end

   self['hedron1'] = hedron1
   self['h1key'] = h1key
   --print(res:tostring() .. ' ' .. self:tostring() .. ' found reverse 1: ' .. hedron1:tostring())
   
   local h2key = utils.genKey(self[4], self[3], self[2])
   --print(res:tostring() .. ' ' .. self:tostring() .. ' reverse searching 2: ' .. h2key)
   local hedron2 = res['hedra'][h2key]
   if not hedron2 and res['prev'] then hedron2 = res['prev']['hedra'][h2key] end
   if not hedron2 and res['next'] then hedron2 = res['next']['hedra'][h2key] end
   if not hedron2 then 
      assert(nil,'residue ' .. res['res'] .. res['resn'] .. ' failed to locate 2nd hedron (R) ' .. h2key .. ' a1= ' .. self[1])
   end
   self['hedron2'] = hedron2
   self['h2key'] = h2key
   --print(res:tostring() .. ' ' .. self:tostring() .. ' found reverse 2: ' .. hedron2:tostring())

   local initial = {}
   initial[1] = hedron1['atomsR'][3]:copy()
   initial[2] = hedron1['atomsR'][2]:copy()
   initial[3] = hedron1['atomsR'][1]:copy()

   local a4preRotation = hedron2['atoms'][1]:copy()
   a4preRotation[3][1] = a4preRotation[3][1] * -1                   -- a4 to +Z
   a4preRotation[3][1] = a4preRotation[3][1] + hedron2['len3']        -- hedron2 shift up so a2 at 0,0,0
   
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
   
   --print('qipR: ' .. utils.genKey(self['a1'], self['a2'], self['a3'],self['a4']) .. ' ' .. h1key .. ' ' .. h2key .. ' [' ..  self['a1initial']:pretty() .. '][' ..  self['a2initial']:pretty() .. '][' ..  self['a3initial']:pretty() .. '][' ..  self['a4initial']:pretty() .. ']') 

   self['updated'] = false
end

return dihedron
