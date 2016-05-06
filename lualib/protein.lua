#!/usr/bin/env luajit

require 'rtmLualib'
require 'rtmGeom3d'
local deque = require 'deque'

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
   assert(o[3],'attempt to instantiate hedron without atom 3')
   setmetatable(o, self)
   self.__index = self

   o['key'] = genKey(o[1],o[2],o[3])
   o['updated'] = true

   if not o['atoms'] then
      o['atoms'] = {}
   end
   if not o['atomsR'] then
      o['atomsR'] = {}
   end
   
   return o
end

function P.hedron:tostring()
   return '3-[' .. self['key'] .. ']'
end

function P.hedron:initPos()
   -- generate coordinates for 3 atoms with specified bond lengths and angle between on XZ plane (Y=0 for all atoms)
   -- match aacs : a2 = 0,0,0 ; a3 on Z-axis ; a1 on XZ plane (-Z for angle >90)
   -- reverse coordinates swap : a2 = 0,0,0 ; a1 on Z-axis ; a3 on XZ plane
   self['atoms'][1] = get41mtx()   
   self['atoms'][2] = get41mtx()    -- a2 to 0,0,0
   self['atoms'][3] = get41mtx()

   self['atoms'][3][3][1] = self['len3']   -- a3 is len3 up from a2 on Z axis, X=Y=0

   --local ar = self['angle2'] * d2r
   local sar = (180.0 - self['angle2']) * d2r    -- angles which add to 180 are supplementary

   self['atoms'][1][1][1] = math.sin(sar) * self['len1']      -- a1 X is sin( sar ) * len1
   self['atoms'][1][3][1] = - (math.cos(sar) * self['len1'])  -- a1 Z is -(cos( sar ) * len1) (assume angle2 always obtuse, so a1 is in -Z

   -- same again but a3 at 0,0,0
   self['atomsR'][1] = get41mtx()   
   self['atomsR'][2] = get41mtx()   -- a2 to 0,0,0
   self['atomsR'][3] = get41mtx()    

   self['atomsR'][1][3][1] = self['len1']   -- a1r is len1 up from a2 on Z axis, X=Y=0

   --local ar = self['angle2'] * d2r
   local sar = (180.0 - self['angle2']) * d2r    -- angles which add to 180 are supplementary

   self['atomsR'][3][1][1] = math.sin(sar) * self['len3']    -- a3r X is sin( sar ) * len3
   self['atomsR'][3][3][1] = - (math.cos(sar) * self['len3'])  -- a3r Z is -(cos( sar ) * len3)

   self['updated'] = false
end


------------------------------------------------------------------------------------------------
-- dihedron
------------------------------------------------------------------------------------------------

local coordsInitial = 1
local coordsResidue = 2
local coordsStructure = 3
local coordsAltStructure = 4

local function getAtomName(a)
   return string.match(a,'^%a%d+(%S+)$')
end

function P.dihedron:new (o)
   --o = o or {}
   assert(o[4],'attempt to instantiate dihedron without atom 4')

   setmetatable(o, self)
   self.__index = self

   o['key'] = genKey(o[1],o[2],o[3],o[4])
   o['key3'] = genKey(o[1],o[2],o[3])
   o['key32'] = genKey(o[2],o[3],o[4])
   
   o['updated'] = true

   if not o['atomNames'] then
      o['atomNames'] = {}
   end

   for i=1,4 do
      o['atomNames'][ i ] =  getAtomName(o[i])
   end
   
   if not o['atomCoords'] then
      o['atomCoords'] = {}
   end

   return o
end

function P.dihedron:tostring()
   return '4-[' .. self['key'] .. ']'
end


function P.dihedron:initPos()
   -- generate coordinates for 4 atoms with specified dihedral, first 3 on XZ plane (a1 in -Z), 4th in +Z and rotated
   local res = self['res']

   local h1key = genKey(self[1], self[2], self[3])
   --print(res:tostring() .. ' ' .. self:tostring() .. ' searching 1: ' .. h1key)
   
   local hedron1 = res['hedra'][h1key]
   if not hedron1 and res['prev'] then hedron1 = res['prev']['hedra'][h1key] end
   if not hedron1 and res['next'] then hedron1 = res['next']['hedra'][h1key] end
   if not hedron1 then return self:initPosR() end
   self['hedron1'] = hedron1
   self['h1key'] = h1key
   --print(res:tostring() .. ' ' .. self:tostring() .. ' found 1: ' .. hedron1:tostring())
   
   local h2key = genKey(self[2], self[3], self[4])
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

   --[[
   self['a1initial'] = hedron1['a1initial']:copy()
   self['a2initial'] = hedron1['a2initial']:copy()
   self['a3initial'] = hedron1['a3initial']:copy()

   self['a4initial'] = hedron2['a3initialr']:copy()
   self['a4initial'][3][1] = self['a4initial'][3][1] * -1                   -- a4 to +Z
   self['a4initial'][3][1] = self['a4initial'][3][1] + hedron2['len1']        -- hedron2 shift up so a2 at 0,0,0
   --]]

   local initial = {}
   initial[1] = hedron1['atoms'][1]:copy()
   initial[2] = hedron1['atoms'][2]:copy()
   initial[3] = hedron1['atoms'][3]:copy()

   local a4preRotation = hedron2['atomsR'][3]:copy()
   
   a4preRotation[3][1] = a4preRotation[3][1] * -1                   -- a4 to +Z
   a4preRotation[3][1] = a4preRotation[3][1] + hedron2['len1']      -- hedron2 shift up so a2 at 0,0,0
   
   local mrz = genMrz( self['dihedral1'] * d2r )
   initial[4] = mrz * a4preRotation                              -- dihedral set

   --[[
   initial['names']={}
   for i=1,4 do
      initial['names'][i] = getAtomName(self[i])
   end
   --]]
   
   self['atomCoords'][coordsInitial] = initial                             -- atoms[1] = initial coords;  [2] = residue coordinate system;  [3] = structure coordinate system;  [4]... alternates in structure coordinate system

   self['a4preRotation'] = a4preRotation
   
   -- print('qip : ' .. genKey(self['a1'], self['a2'], self['a3'],self['a4']) .. ' ' .. h1key .. ' ' .. h2key .. ' [' ..  self['a1initial']:transpose():pretty() .. '][' ..  self['a2initial']:transpose():pretty() .. '][' ..  self['a3initial']:transpose():pretty() .. '][' ..  self['a4initial']:transpose():pretty() .. ']')

   self['updated'] = false
end

function P.dihedron:initPosR()
   -- generate coordinates for 4 atoms with specified dihedral, first 3 on XZ plane with Z>=0, 4th in +Z -- hedrons in reverse order from above
   local res = self['res']

   local h1key = genKey(self[3], self[2], self[1])
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
   
   local h2key = genKey(self[4], self[3], self[2])
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

   --[[
   self['a1initial'] = hedron1['a3initialr']:copy()
   self['a2initial'] = hedron1['a2initialr']:copy()
   self['a3initial'] = hedron1['a1initialr']:copy()

   self['a4initial'] = hedron2['a1initial']:copy()   
   self['a4initial'][3][1] = self['a4initial'][3][1] * -1                   -- a4 to +Z
   self['a4initial'][3][1] = self['a4initial'][3][1] + hedron2['len3']        -- hedron2 shift up so a2 at 0,0,0  -- reverse from above so  use len3

   local mrz = genMrz( self['dihedral1'] * d2r )
   self['a4rot'] = mrz * self['a4initial']                               -- dihedral set
   --]]

   local initial = {}
   initial[1] = hedron1['atomsR'][3]:copy()
   initial[2] = hedron1['atomsR'][2]:copy()
   initial[3] = hedron1['atomsR'][1]:copy()

   local a4preRotation = hedron2['atoms'][1]:copy()
   a4preRotation[3][1] = a4preRotation[3][1] * -1                   -- a4 to +Z
   a4preRotation[3][1] = a4preRotation[3][1] + hedron2['len3']        -- hedron2 shift up so a2 at 0,0,0
   
   local mrz = genMrz( self['dihedral1'] * d2r )
   initial[4] = mrz * a4preRotation                              -- dihedral set

   --[[
   initial['names']={}
   for i=1,4 do
      initial['names'][i] = getAtomName(self[i])
   end
   --]]
   
   self['atomCoords'][coordsInitial] = initial                                   -- atoms[1] = initial coords;  [2] = residue coordinate system;  [3] = structure coordinate system;  [4]... alternates in structure coordinate system
   self['a4preRotation'] = a4preRotation
   
   --print('qipR: ' .. genKey(self['a1'], self['a2'], self['a3'],self['a4']) .. ' ' .. h1key .. ' ' .. h2key .. ' [' ..  self['a1initial']:pretty() .. '][' ..  self['a2initial']:pretty() .. '][' ..  self['a3initial']:pretty() .. '][' ..  self['a4initial']:pretty() .. ']') 

   self['updated'] = false
end


------------------------------------------------------------------------------------------------
-- residue
------------------------------------------------------------------------------------------------

local backboneSort = { N = 1, CA = 2, C = 3, O = 4 }
local sidechainSort = { CB = 1,
                        CG = 2, CG1 = 2, OG = 2, OG1 = 2, SG = 2,
                        CG2 = 3,
                        CD = 4, CD1 = 4, SD = 4, OD1 = 4, ND1 = 4,
                        CD2 = 5,  ND2 = 5, OD2 = 5,
                        CE = 6, NE = 6, CE1 = 6, OE1 = 6, NE1 = 6,
                        CE2 = 7, OE2 = 7, NE2 = 7,
                        CE3 = 8,
                        CZ = 9, CZ2 = 9, NZ = 9,
                        CH2 = 10, NH1 = 10, OH = 10,
                        NH2 = 11
}

local res3 = { G = 'GLY', A = 'ALA', V = 'VAL', L = 'LEU', I = 'ILE', M = 'MET', F = 'PHE', P = 'PRO', S = 'SER', T = 'THR', C = 'CYS', N = 'ASN', Q = 'GLN', Y = 'TYR', W = 'TRP',
               D = 'ASP', E = 'GLU', H = 'HIS', K = 'LYS', R = 'ARG' }

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

   if not o['backbone'] then
      o['backbone'] = {}
   end
   if not o['sidechain'] then
      o['sidechain'] = {}
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
   local k3i = {}
   local k32i = {}
   for k,dihedron in pairs(self['dihedra']) do
      dihedron['res'] = self          -- each dihedron can find its residue
      
      local k3 = dihedron['key3']
      local k32 = dihedron['key32']
      --print(dihedron['key'],k3,k32)
      if not k3i[k3] then
         k3i[k3] = {}
      end
      if not k32i[k32] then
         k32i[k32] = {}
      end
      k3i[k3][ #k3i[k3] +1 ] = dihedron       -- hash to find each dihedron from 1,2,3 and 2,3,4 
      k32i[k32][ #k32i[k32] +1 ] = dihedron

      for i=1,4 do                            -- PDB format ordered lists of backbone and sidechain atoms in residue
         local a = dihedron['atomNames'][i]
         local b = backboneSort[ a ]
         local s = sidechainSort[ a ]
         if b and not self['backbone'][b] then
            self['backbone'][b] = { dihedron, i }
            --print('backbone ' .. a .. ' ' .. b .. ' '..  dihedron:tostring() .. ' ' .. i)
         elseif s and not self['sidechain'][s] then
            self['sidechain'][s] = { dihedron, i }
         elseif not (b or s) then
            assert(nil, 'cannot identify atom ' .. a .. ' dihedron ' .. dihedron:tostring() .. ' position ' .. i)
         end
      end
   end

   self['key3index'] = k3i
   self['key32index'] = k32i
end

function P.residue:toPDB(chain,ndx)
   --print('residue topdb '  .. chain .. ' ' .. ndx)
   local s = ''
   local dihedron
   local position
   local atomOffset= coordsResidue --coordsInitial
   local atomName
   local atom
   local ls
   for i,a in ipairs(self['backbone']) do
      if a then
         ndx = ndx + 1
         dihedron = a[1]
         position = a[2]
         --print('res:toPDB chain ' .. chain .. ' ndx ' .. ndx .. ' ' .. dihedron:tostring() .. ' pos ' .. position)
         atomName = dihedron['atomNames'][position]  -- dihedron['atomCoords'][coordsInitial]['names'][position]
         atom = dihedron['atomCoords'][atomOffset][position]
         ls =  string.format('%6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n',
                             'ATOM  ',ndx,atomName,' ',res3[self['res']], chain, self['resn'],
                             ' ', atom[1][1], atom[2][1],atom[3][1],
                             1.0, 0.0, '    ', string.sub(atomName,1,1),'  ')
         --print(ls)
         s = s .. ls
      end
   end
   for i,a in ipairs(self['sidechain']) do
      if a then
         ndx = ndx + 1
         dihedron = a[1]
         position = a[2]
         atomName = dihedron['atomNames'][position]  --  dihedron['atomCoords'][coordsInitial]['names'][position]
         atom = dihedron['atomCoords'][atomOffset][position]
         ls = string.format('%6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n',
                                'ATOM  ',ndx,atomName,' ',res3[self['res']], chain, self['resn'],
                                ' ', atom[1][1], atom[2][1],atom[3][1],
                                1.0, 0.0, '    ', string.sub(atomName,1,1),'  ')
         s = s .. ls
      end
   end

   return s,ndx
end
   
function P.residue:assemble()
   --[[
   for i,d in pairs(self['dihedra']) do
      print(i,d['key'])
   end
   os.exit()
   --]]
   
   local rbase = self['res'] .. self['resn']

   local q = deque:new()
   q:push_left(genKey(rbase .. 'N', rbase .. 'CA', rbase .. 'C'))
   q:push_left(genKey(rbase .. 'O', rbase .. 'C', rbase .. 'CA'))
   q:push_left(genKey(rbase .. 'N', rbase .. 'CA', rbase .. 'CB'))
   
   while not q:is_empty() do
      local h1k = q:pop_right()
      local dihedra1 = self['key3index'][h1k]
      print('start h1k ' .. h1k, dihedra1)
      if dihedra1 then
         for d1ndx, d1 in ipairs(dihedra1) do
            if not d1['atomCoords'][coordsResidue] then
               d1['atomCoords'][coordsResidue] = {}
            end
            --[[
            split dihed key -> atomKeys
            for i=1,4 do
               
            --]]
            local d1key = d1['key']
            local d1h2key = d1['hedron2']['key'] --d1key:gsub('^%w+:','')
            print('adding hedron ' .. d1h2key)
            q:push_left(d1h2key)
            local dihedra2 = self['key3index'][d1h2key]
            print(h1k, d1key, d1h2key, dihedra2 )
            if dihedra2 then
               for d2ndx, d2 in ipairs(dihedra2) do
                  print(h1k, d1key, d1h2key, d2['key'],d2['atomCoords'][coordsResidue] )
                  assert(not d2['atomCoords'][coordsResidue],'d2 already has coordsResidue')
                  local d1h2a = d1['hedron2']['atoms']
                  local mt, mtr = coordSpace(d1h2a[1],d1h2a[2],d1h2a[3],true)
                  local d2CoordsResidue = {}
                  for aindx, ai in ipairs(d2['atomCoords'][coordsInitial]) do
                     d2CoordsResidue[aindx] = mtr * ai
                  end
                  d2['atomCoords'][coordsResidue] = d2CoordsResidue
               end
            end
         end
      end
   end

   print(self:toPDB('A',1))
   os.exit()
   
   --[[

      n-ca-c 1st hedron
      init list with h1

      while list not empty
      get h1
      get dihedrons for h1
      if no residueCoords,
        copy initialCoords for dihedrons
      
        
      for h2s in dihedrons
        get mtr for return from coordsResidue for h2
        for all dihedra2 starting with h2,
          coordsResidue = coordsInitial * m2r
      
        add h2 to list for processing as h1s
      

      start with n ca c = 1st hedron
      find dihedron1 starting with 1st hedron
      1:
      get 2nd hedron
      find dihedron2 starting with 2nd hedron
      get mtr for return from coord space for 2nd hedron
      apply to dihedron2
      go to 1

   --]]
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
      res = self:initRes(getRes(t[1]))
      if t['dihedral1'] then
         --print('dihedron: ' .. res:tostring() .. ' ' .. genKey(t[1],t[2],t[3],t[4]))
         res['dihedra'][ genKey(t[1],t[2],t[3],t[4]) ] = P.dihedron:new(t)
      elseif t['angle2'] then
         --print('hedron: ' .. res:tostring() .. ' ' ..  genKey(t[1],t[2],t[3]))
         res['hedra'][ genKey(t[1],t[2],t[3]) ] = P.hedron:new(t)
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

function P.chain:renderDihedrals()
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

function P.chain:assembleResidues()
   for i,r in ipairs(self['residues']) do
      r:assemble()
   end
end

function P.chain:toPDB(ndx)
   print('chain topdb '  .. ndx)
   local s = ''
   local ls
   for i,v in ipairs(self['residues']) do
      ls,ndx = v:toPDB(self['id'],ndx)
      s = s .. ls
   end
   return s,ndx
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

function P.protein:toPDB()
   print('protein topdb ')   
   local s=''
   local ndx=0
   local ls
   for k,v in pairs(self['chains']) do
      ls, ndx = v:toPDB(ndx)
      s = s .. ls
   end
   return s
end

function P.protein:renderDihedrals()
   for k,v in pairs(self['chains']) do
      v:renderDihedrals()
   end
end

function P.protein:assembleResidues()
   for k,v in pairs(self['chains']) do
      v:assembleResidues()
   end
end

   
return P
