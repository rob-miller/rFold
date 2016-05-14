
--- Residue structure classes for rFold 
--
-- **external dependencies** : deque (https://github.com/catwell/cw-lua/tree/master/deque)
--
-- @classmod Residue

local utils = require 'rfold.utils'
local geom3d = require 'rfold.geom3d'
local hedron = require 'rfold.hedron'
local dihedron = require 'rfold.dihedron'
local deque = require 'deque'

local residue = {}  -- module

local Residue = {}  -- class table

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

   return o
end

--- callback from file parser, import table data for this Residue according to contents
-- @param t parsed file record data: DSSP record, or hedron / dihedron specification
function Residue:load(t)
   if t['psi'] then -- dssp record
      self['dssp'] = t
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
         assert(nil, self:tostring() .. ' loading ' .. tostring(t) .. ' -- not recognised')
      end
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

--- create map of first 3 atoms to dihedron for all dihedra in Residue; create backbone and sidechain tables of (dihedron, position) for atoms in PDB order
-- @see toPDB
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
               assert(nil, 'cannot identify atom ' .. a .. ' dihedron ' .. dihedron:tostring() .. ' position ' .. i)
            end
         end
      end
   end

   self['key3index'] = k3i
   --self['key32index'] = k32i
end

--- generate ATOM records for this Residue
--<br>
-- Note: OXT defined in Dihedra like other atoms, if present
-- @todo make use of temperature factor field
-- @param chain chain ID
-- @param ndx ATOM record sequence counter
-- @return string containing sequential ATOM records, ndx
function Residue:toPDB(chain,ndx)
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
         --print('res:toPDB chain ' .. chain .. ' ndx ' .. ndx .. ' ' .. dihedron:tostring() .. ' pos ' .. position)
         atomName = utils.splitAtomKey(dihedron[position])[3]   -- dihedron['atomNames'][position]  -- dihedron['atomCoords'][coordsInitial]['names'][position]
         local akl = utils.splitKey(dihedron['key'])
         atom = self['atomCoords'][akl[position]] --  dihedron['initialCoords'][position]
         ls =  string.format('%6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n',
                             'ATOM  ',ndx,atomName,' ',utils.res3[self['res']], chain, self['resn'],
                             ' ', atom[1][1], atom[2][1],atom[3][1],
                             1.0, 0.0, '    ', string.sub(atomName,1,1),'  ')
         --print(ls)
         s = s .. ls
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
               local akl = utils.splitKey(dihedron['key'])
               atom = self['atomCoords'][akl[position]] -- dihedron['initialCoords'][position]
               ls = string.format('%6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n',
                                  'ATOM  ',ndx,atomName,' ',utils.res3[self['res']], chain, self['resn'],
                                  ' ', atom[1][1], atom[2][1],atom[3][1],
                                  1.0, 0.0, '    ', string.sub(atomName,1,1),'  ')
               s = s .. ls
            end
         end
      end
   return s,ndx
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
-- @param atomCoordsIn table of atom_token : 4x1 matrix of protein space coordinates
function Residue:assemble( atomCoordsIn )
   --[[
   for di,d in pairs(self['dihedra']) do
      print('diheron: ' .. d['key'] .. ' angle: ' .. d['dihedral1'])
   end
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
   
--]]

   local atomCoords = atomCoordsIn
   
   local rbase = self['resn'] .. self['res']
   
   local q = deque:new()
   local NCaCKey = utils.genKey(rbase .. 'N', rbase .. 'CA', rbase .. 'C')
   q:push_left(NCaCKey)
   q:push_left(utils.genKey(rbase .. 'O', rbase .. 'C', rbase .. 'CA'))
   q:push_left(utils.genKey(rbase .. 'N', rbase .. 'CA', rbase .. 'CB'))

   if not atomCoords then -- if list of initial coords not passed as parameter, use N-CA-C initial coords from creating dihedral
      atomCoords = {}
      local dl = self['key3index'][NCaCKey]
      for di,d in ipairs(dl) do
         local akl = utils.splitKey( d['key'] )
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
            --print('assemble: ' .. h1k .. ' -> ' .. d['key'])
            local dh2key = d['hedron2']['key']
            local akl = utils.splitKey( d['key'] )
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
               --print('finished: ' .. d['key'] .. ' adding hedron ' .. dh2key .. ' a4: ' .. akl[4] .. ' -- ' .. atomCoords[akl[4]]:transpose():pretty())
               q:push_left(dh2key)
            else
               print('no coords to start ' .. d['key'])
               q:push_left(h1k)
            end
         end
      end
   end

   --[[
   for k,a in pairs(atomCoords) do
      print('assemble final:',k, a:transpose():pretty())
   end
   --]]
   
   return atomCoords
end

return residue
