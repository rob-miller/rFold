
--- core protein chain class for rFold
--
-- @classmod Chain

local utils = require 'rfold.utils'
local residue = require 'rfold.residue'

local chain = {}  -- module
local Chain = {}  -- class table

------------------------------------------------------------------------------------------------
-- Chain
------------------------------------------------------------------------------------------------

--- Chain class object initialiser (not a class method)
-- @param o table with field 'id' = chain ID; ' ' or '' is valid
-- @return minimally initialised Chain object
function chain.new (o)
   assert((o and o['id']),'chain.new() called without id')
   assert(o['pdbid'],'chain.new() called without pdbid')
   setmetatable(o, { __index = Chain })

   if not o['residues'] then
      o['residues'] = {}
   end
   
   return o
end

--- find or initialise Residue object for specified amino acid and sequence position
-- @param ires uppercase 1-letter amino acid code
-- @param iresn sequence position of residue in this chain
-- @return Residue object, minimally initialised and referenced from self['residues'][iresn]
function Chain:getResidue( ires, iresn )
   if not self['residues'][iresn] then
      self['residues'][iresn] = residue.new{ resn = iresn, res = ires }
   elseif ires ~= self['residues'][iresn]['res'] then
      assert(nil,'stored residue ' .. self['residues'][iresn]['res'] .. iresn .. ' does not match new residue ' .. ires .. iresn)
   end
   return self['residues'][iresn]
end

--- find or create Residue for passed data table, pass table to it for loading
-- @param t table created by file parser, passed here as callback
function Chain:load(t)
   local res

   if t['resn'] then -- dssp record or pdb atom record
      --print('dssp: ' .. t['res'] .. t['resn'])
      if t['res3'] then t['res'] = utils.res1[t['res3']] end
      res = self:getResidue(t['res'], tonumber(t['resn']))
   else
      --print(t[1])
      local akl = utils.splitAtomKey(t[1])
      res = self:getResidue( akl[2], akl[1] )
   end

   res:load(t)
   --print(res['res'] .. res['resn'] .. ' ' .. self:countDihedra())

   
end

--- implement assign 'prev' and 'next' for passed Residues[] index and residue
function Chain:setPrevNext(i,r)
   if i > 1 then r['prev'] = self['residues'][i-1] end
   if i < #self['residues'] then r['next'] = self['residues'][i+1] end
end

--- assign 'prev' and 'next' for Residues in self['residues']; trigger Residue to link its dihedra (create map of keys to dihedra, link dihedron atoms into backbone and sidechain tables)
function Chain:linkResidues()
   for i,r in ipairs(self['residues']) do
      self:setPrevNext(i,r)
      r:linkDihedra()
   end
end

--- query each Residue for its count of hedra
-- @return sum of hedra for Residues in chain
function Chain:countHedra()
   local c=0
   for i,v in ipairs(self['residues']) do
      --print(i,v['res'], v['resn'])
      c = c + v:countHedra()
   end
   return c
end

--- query each Residue for its count of dihedra
-- @return sum of dihedra for Residues in chain
function Chain:countDihedra()
   local c=0
   for i,v in ipairs(self['residues']) do
      --print(v['res'] .. v['resn'] .. ' ' .. v:countDihedra())
      c = c + v:countDihedra()
   end
   return c
end

--- query each Resdiue for a DSSP record
-- @return sum of Residues in chain with DSSP records
function Chain:countDSSPs()
   local c=0
   for i,v in ipairs(self['residues']) do
      if v['dssp'] then c = c + 1 end
   end
   return c
end

--- concatenate 1-letter amino acid codes for Residues in this chain
-- @return chain amino acid sequence as string
function Chain:seqStr()
   local s=''
   for i,v in ipairs(self['residues']) do
      s = s .. v['res']
   end
   return s
end

--- trigger generation of atom coordinates for updated Hedra in Residues in this Chain, then do same for updated Dihedra.
-- <br>
-- Dihedra depend on hedra in adjacent Residues, so must complete hedra first
function Chain:renderDihedra()
   for i,r in ipairs(self['residues']) do r:renderHedra() end
   for i,r in ipairs(self['residues']) do r:renderDihedra() end
end   

--- foreach chain, complete residue, dihedra, hedra data structures from protein space atom coordinates
function Chain:dihedraFromAtoms()
   for i,r in ipairs(self['residues']) do
      self:setPrevNext(i,r)
      r:dihedraFromAtoms()
   end
end

--- for each Residue in Chain, set startPos and then trigger Residue to assemble atoms from its Dihedrons starting with the startPos coordinates for N, CA, C
--<br>
-- first Residue startPos is set from DSSP or other source to match PDB file, subsequent startPos's read from previous residue backbone atom coordinates in protein coordinate space
-- @todo incorporate parallel threads here: divide into n segments for n threads, then assemble segments together at end
-- @see setStartCoords
function Chain:assembleResidues()
   local c=1
   local ndx=0
   for i,r in ipairs(self['residues']) do
      local startPos
      if r['prev'] then
         startPos={}
         local rp = r['prev']
         local akl = r:NCaCKeySplit()
         for ai,ak in ipairs(akl) do
            startPos[ak] = rp['atomCoords'][ak]
            --print('start from previous: ' .. ak .. ' ' .. startPos[ak]:transpose():pretty())
         end
      else
         startPos = self['initNCaC']
      end
      r['atomCoords'] = r:assemble(startPos)

      --s,ndx = r:toPDB('A',ndx,r['atomCoords'])
      --io.write(s)
      --print()
      --c = c+1
      --if c>4 then os.exit() end
      
   end
end

--- trigger each Residue in Chain to generate PDB ATOM records, add TER record at end of chain
-- @return string of PDB format records for ATOMS in Chain, plus TER record
function Chain:toPDB(ndx)
   --print('chain topdb '  .. ndx)
   local s = ''
   local ls
   for i,v in ipairs(self['residues']) do
      ls,ndx = v:toPDB(self['id'],ndx)
      s = s .. ls
   end
   local res = self['residues'][#self['residues']]
   ndx = ndx + 1
   s = s .. string.format('TER   %5d      %3s %1s%4d%54s\n',ndx, utils.res3[ res['res'] ], self['id'], res['resn'], '')

   return s,ndx
end

function Chain:toInternalCoords()
   local s = ''

   if self['initNCaC'] then
      local n = 1
      local r = self['residues'][1]
      local akl =  r:NCaCKeySplit()
      for ai,ak in ipairs(akl) do
         local aks = utils.splitAtomKey(ak)
         s = s .. utils.atomString(n,aks[3],r['res'],self['id'],r['resn'],self['initNCaC'][ak][1][1],self['initNCaC'][ak][2][1],self['initNCaC'][ak][3][1],1.0,0.0)
         n = n+1
      end
   end
   for i,r in ipairs(self['residues']) do
      s = s .. r:toInternalCoords(self['pdbid'], self['id'])
   end
   return s
end

--- if first Residue in Chain has 'dssp' field, set self['initNCaC'] to hold those protein space coordinates to build the rest of the Chain from
-- @see assembleResidues
function Chain:setStartCoords()
   local r = self['residues'][1]
   local akl =  r:NCaCKeySplit()
   local initCoords
   if r['dssp'] then
      initCoords = {}
      for ai,ak in ipairs(akl) do
         initCoords[ak] = r:dsspAtom(ak)
      end
   elseif r['atomCoords'][akl[1]] and r['atomCoords'][akl[2]] and r['atomCoords'][akl[3]] then
      initCoords = {}
      for ai,ak in ipairs(akl) do
         initCoords[ak] = r['atomCoords'][ak]
      end
   else
      print('no start coordinates set')
   end
   self['initNCaC'] = initCoords
end

return chain
