--[[
   chain.lua
   
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

--- protein chain class for rFold
--
-- @classmod Chain

local utils = require 'rfold.utils'
local chemdata = require 'rfold.chemdata'
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
   
   if not o['initNCaC'] then
      o['initNCaC'] = {}
   end
   
   return o
end

--- find or initialise Residue object for specified amino acid and sequence position
-- @param ires uppercase 1-letter amino acid code
-- @param iresn sequence position of residue in this chain
-- @return Residue object, minimally initialised and referenced from self['residues'][iresn]
function Chain:getResidue( ires, iresn )
   if 'X' == ires then return nil end   -- dssp missing/disordered residue
   if (not self['residues'][iresn]) or (not self['residues'][iresn]['ordered'])  then   -- no residue or only placeholder at this position
      self['residues'][iresn] = residue.new{ resn = iresn, res = ires, ordered = true }
      -- lua needs arrays starting at 1 with no gaps, so fill in as needed here with ' ' instead of residue opbject
      if not self['firstPos'] then
         self['firstPos'] = iresn
         if 1 < iresn then   
            for i=1,iresn-1 do self['residues'][i] = {} end
         end
      elseif not self['residues'][iresn-1] then -- there is a gap before this position (possibly this is atom record indicating start after chain break)
         for i=self['firstPos'],iresn-1 do
            if not self['residues'][i] then self['residues'][i] = {} end
         end
      end
   elseif ires ~= self['residues'][iresn]['res'] then
      print(iresn,ires,self['residues'][iresn]['resn'],self['residues'][iresn]['res'])
      assert(nil,'stored residue ' .. self['residues'][iresn]['res'] .. iresn .. ' does not match new residue ' .. ires .. iresn)
   end
   return self['residues'][iresn]
end

--- find or create Residue for passed data table, pass table to it for loading
-- @param t table created by file parser, passed here as callback
function Chain:load(t)
   local res

   --for k,v in pairs(t) do print('c: ' .. k,v) end
   --print()
   
   if t['resn'] then -- dssp record or pdb atom record
      --print('dssp: ' .. t['res'] .. t['resn'])
      if t['res3'] then t['res'] = chemdata.res1[t['res3']] end
      --print(t['chn'],t['res'],t['resn'])
      if (t['res'] and not t['icode']) and ((not t['altloc']) or ('A' == t['altloc'])) then     -- standard residue, no insdertion codes, only first altloc if present
         --print(t['res'],t['resn'],t['atom'],t['altloc'])
         res = self:getResidue(t['res'], tonumber(t['resn']))
      end
   elseif t['seqres'] then -- pdb seqres record
      if self['seqres'] then
         self['seqres'] = self['seqres'] .. t['seqres']
      else
         self['seqres'] = t['seqres']
      end
   else  -- hedron or dihedron record
      --print(t[1])
      local akl = utils.splitAtomKey(t[1])
      res = self:getResidue( akl[2], akl[1] )
   end

   if res then   -- might be missing
      res:load(t)
      --print(res['res'] .. res['resn'] .. ' ' .. self:countDihedra())
   end

   
end

--- custom iterator for residue array to skip missing residues and start at sequence positions other than 1
-- @param a array of residues to iterate in numeric order
-- @param i initial value
local function ordResIter(a,i)
   local r={}
   while true do
      i = i+1
      local r = a[i]
      if r and r['ordered'] then return i,r
      elseif not r then return nil end
   end
end

--- initialiser for Chain:oResIter()
-- @param a array of residues to iterate in numeric order
function Chain:orderedResidues()
   return ordResIter, self['residues'], self['firstPos']-1
end
   
--- custom iterator for residue array to process all residue positions (including disordered) and start at sequence positions other than 1
-- @param a array of residues to iterate in numeric order
-- @param i initial value
local function allResIter(a,i)
   local r={}
   while true do
      i = i+1
      local r = a[i]
      if r then return i,r 
      elseif not r then return nil end
   end
end

--- initialiser for Chain:aResIter()
-- @param a array of residues to iterate in numeric order
function Chain:allResidues()
   return allResIter, self['residues'], self['firstPos']-1
end
   

--- implement assign 'prev' and 'next' for passed Residues[] index and residue
function Chain:setPrevNext(i,r)
   if i > self['firstPos'] then r['prev'] = self['residues'][i-1] end
   if i < #self['residues'] then r['next'] = self['residues'][i+1] end
end

--- assign 'prev' and 'next' for Residues in self['residues']; trigger Residue to link its dihedra (create map of keys to dihedra, link dihedron atoms into backbone and sidechain tables)
function Chain:linkResidues()
   for i,r in self:orderedResidues() do
      self:setPrevNext(i,r)
      r:linkDihedra()
   end
end

--- query each Residue for its count of hedra
-- @return sum of hedra for Residues in chain
function Chain:countHedra()
   local c=0
   for i,r in self:orderedResidues() do
      --print(i,v['res'], v['resn'])
      c = c + r:countHedra()
   end
   return c
end

--- query each Residue for its count of dihedra
-- @return sum of dihedra for Residues in chain
function Chain:countDihedra()
   local c=0
   for i,r in self:orderedResidues() do
      --print(r['res'] .. r['resn'] .. ' ' .. r:countDihedra())
      c = c + r:countDihedra()
   end
   return c
end

--- query each Resdiue for a DSSP record
-- @return sum of Residues in chain with DSSP records
function Chain:countDSSPs()
   local c=0
   for i,r in self:orderedResidues() do
      if r['dssp'] then c = c + 1 end
   end
   return c
end

function Chain:countInitNCaCs()
   local c=0
   if {} ~= self['initNCaC'] then
      local n = 1
      for i,r in self:orderedResidues() do
         if self['initNCaC'][i] then c = c + 1 end
      end
   end
   return c
end

            
--- delete protein space atom coordinate data from this chain
function Chain:clearAtomCoords()
   for i,r in self:orderedResidues() do
      r:clearAtomCoords()
   end
end
   
--- delete internal coordinate data from this chain
function Chain:clearInternalCoords()
   for i,r in self:orderedResidues() do
      r:clearInternalCoords()
   end
end

--- concatenate 1-letter amino acid codes for Residues in this chain
-- @return chain amino acid sequence as string
function Chain:seqStr()
   local s=''
   for i,r in self:orderedResidues() do
      s = s .. r['res']
   end
   return s
end

--- trigger generation of atom coordinates for updated Hedra in Residues in this Chain, then do same for updated Dihedra.
-- <br>
-- Dihedra depend on hedra in adjacent Residues, so must complete hedra first
function Chain:renderDihedra()
   for i,r in self:orderedResidues() do r:renderHedra() end
   for i,r in self:orderedResidues() do r:renderDihedra() end
end   

--- foreach chain, complete residue, dihedra, hedra data structures from protein space atom coordinates
function Chain:dihedraFromAtoms()
   for i,r in self:orderedResidues() do
      self:setPrevNext(i,r)  -- set 'prev' and 'next' for this residue
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
   for i,r in self:orderedResidues() do
      local startPos
      if r['prev'] and r['prev']['ordered'] then  -- if sequential, i.e. no chain breaks
         startPos={}
         local rp = r['prev']
         local akl = r:NCaCKeySplit()
         for ai,ak in ipairs(akl) do
            startPos[ak] = rp['atomCoords'][ak]
            --print('start from previous: ' .. ak .. ' ' .. startPos[ak]:transpose():pretty())
         end
      else
         startPos = self['initNCaC'][r['resn']]
      end
      r['atomCoords'] = r:assemble(startPos)

      --s,ndx = r:writePDB('A',ndx,r['atomCoords'])
      --io.write(s)
      --print()
      --c = c+1
      --if c>4 then os.exit() end
      
   end
end

--- trigger each Residue in Chain to generate PDB ATOM records, add TER record at end of chain
-- @return string of PDB format records for ATOMS in Chain, plus TER record
function Chain:writePDB(ndx)
   --print('chain topdb '  .. ndx)
   local s = ''
   local ls
   local lastRes
   for i,r in self:orderedResidues() do
      ls,ndx = r:writePDB(self['id'],ndx)
      s = s .. ls
      lastRes=i
   end
   local res = self['residues'][lastRes]
   ndx = ndx + 1
   s = s .. string.format('TER   %5d      %3s %1s%4d%54s\n',ndx, chemdata.res3[ res['res'] ], self['id'], res['resn'], '')

   return s,ndx
end

function Chain:writeInternalCoords()
   local s = ''

   if {} ~= self['initNCaC'] then
      local n = 1
      for i,r in self:orderedResidues() do
         if self['initNCaC'][i] then
            local initNCaC = self['initNCaC'][i]
            local akl =  r:NCaCKeySplit()
            for ai,ak in ipairs(akl) do
               local aks = utils.splitAtomKey(ak)
               s = s .. utils.atomString(n,aks[3],r['res'],self['id'],r['resn'],initNCaC[ak][1][1],initNCaC[ak][2][1],initNCaC[ak][3][1],1.0,0.0)
               n = n+1
            end
         end
      end
   end
   for i,r in self:orderedResidues() do
      s = s .. r:writeInternalCoords(self['pdbid'], self['id'])
   end
   return s
end

--- write chain data to rfold database
-- @param rfpg open database handle
-- @param pdb_no chain ID in pdb_chain table
-- @param update optional flag, if false silently skip if [pdb_no, res_ndx] entry exists already in residue table 
function Chain:writeDb(rfpg, pdb_no, update)

   for i,r in self:orderedResidues() do

      local res_id = rfpg.Q("select res_id from residue where pdb_no=" .. pdb_no .. " and res_ndx = " .. i )
      if not res_id then
         res_id = rfpg.Q( "insert into residue (res_ndx,pdb_no) values (" .. i .. "," .. pdb_no .. ") returning res_id" )
      elseif not update then
         goto skipResidue  -- can't happen? because existing chain would have been skipped
      end
      res_id = res_id[1]

      local qc = "update residue set (res, res_pos"
      local qv = " = ('" .. r['res'] .. "','" .. r['resn'] .. "'"
      if r['prev'] then
         qc = qc .. ", prev_res"
         qv = qv .. ",'" .. r['prev']['res'] .. "'"
      end
      if r['next'] then
         qc = qc .. ", next_res"
         qv = qv .. ",'" .. r['next']['res'] .. "'"
      end
      
      rfpg.Qcur( qc .. ')' .. qv .. ') where res_id = ' .. res_id )

      if ({} ~= self['initNCaC'] and self['initNCaC'][i]) then  -- have init coords for chain segment, need to store
         local akl =  r:NCaCKeySplit()
         local initCoords = self['initNCaC'][i]
         for ai,ak in ipairs(akl) do
            local q
            local v = initCoords[ak]
            local haveRow = rfpg.Q("select 1 from atom_coordinates where res_id=" .. i .. " and atom = '" .. ak .. "'")
            if (not haveRow) then
               q = 'insert into atom_coordinates (res_id, atom, x, y, z) values (' .. i .. ",'" .. ak .. "',"
            elseif update then
               q = 'update atom_coordinates set (x, y, z) = ('
            end
            if q then -- ((not haveRow) or update)
               for j=1,3 do q = q .. (j>1 and ',' or '')  .. v[j][1] end
               q = q .. ')'
               if (haveRow) then  -- update
                  q = q .. ' where res_id=' .. i .. " and atom= '" .. ak .."'"
               end
               rfpg.Qcur(q)
            end
         end
      end

      r:writeDb(rfpg, res_id, update)
      
      ::skipResidue::
   end
   

   --[[
   below is all to just list available data at this point
   need to add 3x3 chain inital coords to db tables
   
   local s=''
   
   if {} ~= self['initNCaC'] then
      local n = 1
      for i,r in self:orderedResidues() do
         if self['initNCaC'][i] then
            local initNCaC = self['initNCaC'][i]
            local akl =  r:NCaCKeySplit()
            for ai,ak in ipairs(akl) do
               local aks = utils.splitAtomKey(ak)
               s = s .. utils.atomString(n,aks[3],r['res'],self['id'],r['resn'],initNCaC[ak][1][1],initNCaC[ak][2][1],initNCaC[ak][3][1],1.0,0.0)
               n = n+1
            end
         end
      end
   end

   for i,r in self:orderedResidues() do
      s = s .. r:writeInternalCoords(self['pdbid'], self['id'])
   end
   print(s)
   print('----------------------------------------------------')

   return s
   --]]
end


--- if first Residue in Chain has 'dssp' field, set self['initNCaC'] to hold those protein space coordinates to build the rest of the Chain from
-- @see assembleResidues
function Chain:setStartCoords()
   local c=0
   for i,r in self:orderedResidues() do  -- find any residues with no previous or previous is disordered
      if not (self['residues'][i-1] and self['residues'][i-1]['ordered']) then
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
         self['initNCaC'][i] = initCoords
         if initCoords then c=c+1 end
      end
   end
   if utils.warn then io.stderr:write('set ' .. c .. ' start coords\n') end
end

--- report ordered (with coordinates) residues in this chain
-- @return count of ordered residues
function Chain:countResidues()
   local o,a = 0,0
   for i,r in self:allResidues() do
      a=a+1
      if r['ordered'] then o=o+1 end
   end
   return o,a
end

function Chain:printInfo()
   for k,v in self:orderedResidues() do
      v:printInfo()
   end
end

return chain
