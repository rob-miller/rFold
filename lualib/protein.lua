#!/usr/bin/env luajit

require 'rtmLualib'
require 'rtmGeom3d'

local P = {}

local proteins = {}

P.protein = {}
P.chain = {}
P.residue = {}
P.dihedron = {}
P.hedron = {}


local function genKey(...)
   local key = ''
   for i, s in ipairs{...} do
      if '' ~= key then key = key .. ':' end
      key = key .. s
   end
   return key
end


------------------------------------------------------------------------------------------------
-- hedron
------------------------------------------------------------------------------------------------

function P.hedron:new (o)
   --o = o or {}
   assert(o['a3'],'attempt to instantiate hedron without atom 3')
   setmetatable(o, self)
   self.__index = self

   o['key'] = genKey(o['a1'],o['a2'],o['a3'])
   o['updated'] = true
   
   return o
end

function P.hedron:tostring()
   return '3-[' .. self['key'] .. ']'
end

function P.hedron:initPos()
   -- generate coordinates for 3 atoms with specified bond lengths and angle between on XZ plane (Y=0 for all atoms)
   -- match aacs : a2 = 0,0,0 ; a3 on Z-axis ; a1 on XZ plane
   -- reverse coordinates swap : a2 = 0,0,0 ; a1 on Z-axis ; a3 on XZ plane
   self['a1initial'] = get41mtx()   
   self['a2initial'] = get41mtx()    -- a2 to 0,0,0
   self['a3initial'] = get41mtx()

   self['a3initial'][3][1] = self['len3']   -- a3 is len3 up from a2 on Z axis, X=Y=0

   --local ar = self['angle2'] * d2r
   local sar = (180.0 - self['angle2']) * d2r    -- angles which add to 180 are supplementary

   self['a1initial'][1][1] = math.sin(sar) * self['len1']      -- a1 X is sin( sar ) * len1
   self['a1initial'][3][1] = - (math.cos(sar) * self['len1'])  -- a1 Z is -(cos( sar ) * len1) (assume angle2 always obtuse, so a1 is in -Z

   -- same again but a3 at 0,0,0
   self['a1initialr'] = get41mtx()   
   self['a2initialr'] = get41mtx()   -- a2 to 0,0,0
   self['a3initialr'] = get41mtx()    

   self['a1initialr'][3][1] = self['len1']   -- a1r is len1 up from a2 on Z axis, X=Y=0

   --local ar = self['angle2'] * d2r
   local sar = (180.0 - self['angle2']) * d2r    -- angles which add to 180 are supplementary

   self['a3initialr'][1][1] = math.sin(sar) * self['len3']    -- a3r X is sin( sar ) * len3
   self['a3initialr'][3][1] = - (math.cos(sar) * self['len3'])  -- a3r Z is -(cos( sar ) * len3)

   self['updated'] = false
end


------------------------------------------------------------------------------------------------
-- dihedron
------------------------------------------------------------------------------------------------

function P.dihedron:new (o)
   --o = o or {}
   assert(o['a4'],'attempt to instantiate dihedron without atom 4')

   setmetatable(o, self)
   self.__index = self

   o['key'] = genKey(o['a1'],o['a2'],o['a3'],o['a4'])
   o['updated'] = true

   return o
end

function P.dihedron:tostring()
   return '4-[' .. self['key'] .. ']'
end


function P.dihedron:initPos()
   -- generate coordinates for 4 atoms with specified dihedral, first 3 on XZ plane, 4th in +Z
   local res = self['res']

   local h1key = genKey(self['a1'], self['a2'], self['a3'])
   --print(res:tostring() .. ' ' .. self:tostring() .. ' searching 1: ' .. h1key)
   
   local hedron1 = res['hedra'][h1key]
   if not hedron1 and res['prev'] then hedron1 = res['prev']['hedra'][h1key] end
   if not hedron1 and res['next'] then hedron1 = res['next']['hedra'][h1key] end
   if not hedron1 then return self:initPosR() end
   self['hedron1'] = hedron1
   self['h1key'] = h1key
   --print(res:tostring() .. ' ' .. self:tostring() .. ' found 1: ' .. hedron1:tostring())
   
   local h2key = genKey(self['a2'], self['a3'], self['a4'])
   --print(res:tostring() .. ' ' .. self:tostring() .. ' searching 2: ' .. h2key)

   local hedron2 = res['hedra'][h2key]
   if not hedron2 and res['prev'] then hedron2 = res['prev']['hedra'][h2key] end
   if not hedron2 and res['next'] then hedron2 = res['next']['hedra'][h2key] end
   if not hedron2 then 
      assert(nil,'residue ' .. res['res'] .. res['resn'] .. 'failed to locate 2nd hedronle ' .. h2key .. ' a1= ' .. self['a1'])
   end
   self['hedron2'] = hedron2
   self['h2key'] = h2key
   --print(res:tostring() .. ' ' .. self:tostring() .. ' found 2: ' .. hedron2:tostring())
   
   self['a1initial'] = hedron1['a1initial']:copy()
   self['a2initial'] = hedron1['a2initial']:copy()
   self['a3initial'] = hedron1['a3initial']:copy()

   self['a4initial'] = hedron2['a3initialr']:copy()
   self['a4initial'][3][1] = self['a4initial'][3][1] * -1                   -- a4 to +Z
   self['a4initial'][3][1] = self['a4initial'][3][1] + hedron2['len1']        -- hedron2 shift up so a2 at 0,0,0
   
   local mrz = genMrz( self['dihedral1'] * d2r )
   self['a4rot'] = mrz * self['a4initial']                              -- dihedral set

   -- print('qip : ' .. genKey(self['a1'], self['a2'], self['a3'],self['a4']) .. ' ' .. h1key .. ' ' .. h2key .. ' [' ..  self['a1initial']:transpose():pretty() .. '][' ..  self['a2initial']:transpose():pretty() .. '][' ..  self['a3initial']:transpose():pretty() .. '][' ..  self['a4initial']:transpose():pretty() .. ']')

   self['updated'] = false
end

function P.dihedron:initPosR()
   -- generate coordinates for 4 atoms with specified dihedral, first 3 on XZ plane with Z>=0, 4th in +Z -- hedronles in reverse order from above
   local res = self['res']

   local h1key = genKey(self['a3'], self['a2'], self['a1'])
   --print(res:tostring() .. ' ' .. self:tostring() .. ' searching reverse 1: ' .. h1key)
   local hedron1 = res['hedra'][h1key]
   if not hedron1 and res['prev'] then hedron1 = res['prev']['hedra'][h1key] end
   if not hedron1 and res['next'] then hedron1 = res['next']['hedra'][h1key] end
   if not hedron1 then
      assert(nil,'residue ' .. res['res'] .. res['resn'] .. ' failed to locate 1st hedronle (R) ' .. h1key .. ' a4= ' .. self['a4'])
   end

   self['hedron1'] = hedron1
   self['h1key'] = h1key
   --print(res:tostring() .. ' ' .. self:tostring() .. ' found reverse 1: ' .. hedron1:tostring())
   
   local h2key = genKey(self['a4'], self['a3'], self['a2'])
   --print(res:tostring() .. ' ' .. self:tostring() .. ' reverse searching 2: ' .. h2key)
   local hedron2 = res['hedra'][h2key]
   if not hedron2 and res['prev'] then hedron2 = res['prev']['hedra'][h2key] end
   if not hedron2 and res['next'] then hedron2 = res['next']['hedra'][h2key] end
   if not hedron2 then 
      assert(nil,'residue ' .. res['res'] .. res['resn'] .. ' failed to locate 2nd hedronle (R) ' .. h2key .. ' a1= ' .. self['a1'])
   end
   self['hedron2'] = hedron2
   self['h2key'] = h2key
   --print(res:tostring() .. ' ' .. self:tostring() .. ' found reverse 2: ' .. hedron2:tostring())

   self['a1initial'] = hedron1['a3initialr']:copy()
   self['a2initial'] = hedron1['a2initialr']:copy()
   self['a3initial'] = hedron1['a1initialr']:copy()

   self['a4initial'] = hedron2['a1initial']:copy()
   
   self['a4initial'][3][1] = self['a4initial'][3][1] * -1                   -- a4 to +Z
   self['a4initial'][3][1] = self['a4initial'][3][1] + hedron2['len3']        -- hedron2 shift up so a2 at 0,0,0  -- reverse from above so  use len3
   
   local mrz = genMrz( self['dihedral1'] * d2r )
   self['a4rot'] = mrz * self['a4initial']                               -- dihedral set

   --print('qipR: ' .. genKey(self['a1'], self['a2'], self['a3'],self['a4']) .. ' ' .. h1key .. ' ' .. h2key .. ' [' ..  self['a1initial']:pretty() .. '][' ..  self['a2initial']:pretty() .. '][' ..  self['a3initial']:pretty() .. '][' ..  self['a4initial']:pretty() .. ']') 

   self['updated'] = false
end


------------------------------------------------------------------------------------------------
-- residue
------------------------------------------------------------------------------------------------

function P.residue:new (o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self

   if not o['resn'] then
      o['resn'] = 0
   end
   if not o['res'] then
      o['res'] = ''
   end
   if not o['hedra'] then
      o['hedra'] = {}
   end
   if not o['dihedra'] then
      o['dihedra'] = {}
   end

   return o
end

function P.residue:tostring()
   return self['res'] .. self['resn']
end

function P.residue:countHedra()
   local c = 0
   for k,v in pairs(self['hedra']) do
      c = c+1
   end
   return c
end

function P.residue:countDihedra()
   local c = 0
   for k,v in pairs(self['dihedra']) do
      --print('  ' .. k)
      c = c+1
   end
   return c
end

function P.residue:linkDihedra()
   local c = 0
   for k,v in pairs(self['dihedra']) do
      v['res'] = self
   end
end



------------------------------------------------------------------------------------------------
-- chain
------------------------------------------------------------------------------------------------

function P.chain:new (o)
   o = o or {}
   setmetatable(o, self)
   self.__index = self

   if not o['residues'] then
      o['residues'] = {}
   end
   
   return o
end

function P.chain:initRes( ires, iresn )
   if not self['residues'][iresn] then
      self['residues'][iresn] = P.residue:new{ resn = iresn, res = ires }
      --print('new residue ' .. ires .. ' ' .. iresn)
   else
      --print(' found residue ' .. ires .. ' ' .. iresn)
   end
   if ires ~= self['residues'][iresn]['res'] then
      assert(nil,'stored residue ' .. self['residues'][iresn]['res'] .. iresn .. ' does not match new residue ' .. ires .. iresn)
   end
   return self['residues'][iresn]
end

function P.chain:load(t)
   local res

   local function getRes(s)
      if not s:match('^(%a)(%d+)(%w+)$') then
         assert(nil, 'failed to parse ' .. s .. ' as atom in residue')
      end
      --print('parse: ' .. _1 .. ' ' .. _2 .. ' ' .. _3)
      return _1, tonumber(_2)
   end
   
   if t['psi'] then
      --print('dssp: ' .. t['res'] .. t['resn'])
      res = self:initRes(t['res'], t['resn'])
      res['dssp'] = t
   else
      --print('chain:load not dssp')
      res = self:initRes(getRes(t['a1']))
      if t['dihedral1'] then
         print('dihedron: ' .. res:tostring() .. ' ' .. genKey(t['a1'],t['a2'],t['a3'],t['a4']))
         res['dihedra'][ genKey(t['a1'],t['a2'],t['a3'],t['a4']) ] = P.dihedron:new(t)
         --res:countDihedra()
         --local k = genKey(t['a1'],t['a2'],t['a3'],t['a4'])
         --print(k .. ' : ' .. tostring(res['dihedra'][k]))
      elseif t['angle2'] then
         print('hedron: ' .. res:tostring() .. ' ' ..  genKey(t['a1'],t['a2'],t['a3']))
         res['hedra'][ genKey(t['a1'],t['a2'],t['a3']) ] = P.hedron:new(t)
      else
         assert(nil,'not recognised: ' .. tostring(t))
      end
   end
end

function P.chain:linkResidues()
   for i,v in ipairs(self['residues']) do
      if i > 1 then v['prev'] = self['residues'][i-1] end
      if i < #self['residues'] then v['prev'] = self['residues'][i+1] end
      v:linkDihedra()
   end
end

function P.chain:countHedra()
   local c=0
   for i,v in ipairs(self['residues']) do
      c = c + v:countHedra()
   end
   return c
end

function P.chain:countDihedra()
   local c=0
   for i,v in ipairs(self['residues']) do
      --print(v['res'] .. ' ' .. v['resn'])
      c = c + v:countDihedra()
   end
   return c
end

function P.chain:countDSSPs()
   local c=0
   for i,v in ipairs(self['residues']) do
      if v['dssp'] then c = c + 1 end
   end
   return c
end

function P.chain:seqStr()
   local s=''
   for i,v in ipairs(self['residues']) do
      s = s .. v['res']
   end
   return s
end

function P.chain:update3d()
   for i,r in ipairs(self['residues']) do
      for k,h in pairs(r['hedra']) do
         if h['updated'] then h:initPos() end
      end
   end
   for i,r in ipairs(self['residues']) do
      for k,d in pairs(r['dihedra']) do
         if d['updated'] then d:initPos() end
      end
   end
end   



------------------------------------------------------------------------------------------------
-- protein
------------------------------------------------------------------------------------------------


function P.protein:new (o)
   if o and o['id'] and proteins[o['id']] then
      return proteins[o['id']]
   end
   o = o or {}
   setmetatable(o, self)
   self.__index = self

   if not o['chains'] then
      o['chains'] = {}
   end

   return o
end

function P.protein:get(id)  -- get the protein with this id or init if it does not exist
   if id then
      return P.protein:new{ id = id }
   end
   
   return P.protein:new()
end
   

function P.protein:load (t)
   if not t['pdbid'] then
      assert(nil,'protein:load - no pdbid for ' .. tostring(t))
   end
   if not proteins[t['pdbid']] then
      proteins[t['pdbid']] = self
      self['id'] = t['pdbid']
      --print('load new protein ' .. self['id'])
   else
      --print('load existing protein ' .. self['id'])
   end

   local chn = t['chn'] or ''
   if not self['chains'][chn] then
      self['chains'][chn] = P.chain:new{ id = chn }
      --print('protein:load new chain ' .. chn)
   else
      --print('protein:load found chain ' .. chn)
   end
   self['chains'][chn]:load(t)
end

function P.protein:tostring()
   local s = ''
   for k,v in pairs(self['chains']) do
      s = self.id .. ' ' .. k .. ' ' .. v:seqStr() .. '\n'
      s = s .. '  ' .. #v['residues'] .. ' residues '  .. v:countDSSPs() .. ' dssps ' .. v:countHedra() ..  ' hedra ' .. v:countDihedra() .. ' dihedra'
   end
   return s
end

function P.protein:linkResidues()
   for k,v in pairs(self['chains']) do
      v:linkResidues()
   end
end

function P.protein:update3d()
   for k,v in pairs(self['chains']) do
      v:update3d()
   end
end

   
return P
