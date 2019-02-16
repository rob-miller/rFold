--[[
   chain.lua
   
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

--- protein chain class for rFold
--
-- @classmod Chain

local utils = require 'rfold.utils'
local chemdata = require 'rfold.chemdata'
local residue = require 'rfold.residue'
local geom3d = require 'rfold.geom3d'

local chain = {}  -- module
local Chain = {}  -- class table

------------------------------------------------------------------------------------------------
-- Chain
------------------------------------------------------------------------------------------------
-- class properties

--- array of residue ordered and disordered objects indexed by sequence position
-- @field residues array of residue objects

--- array of atomCoordinates (4x1 float matrices) for initial N,Ca,C residues of all ordered segments of this chain indexed by atomKeys (resPos res atom)
-- @field initNCaC

--- chain ID from PDB file
-- @field id char

--- 4-position PDB identifier
-- @field pdbid string

--- signed integer indicating sequence position of first residue in chain
-- @field firstPos integer

--- PDB file SEQRES data for this chain
-- @field seqres string

--- Chain class object initialiser (not a class method)
-- @param o table with field 'id' = chain ID; ' ' or '' is val
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

--- query database for info about specified chain (not a class method)
-- @param rfpg open database handle
-- @param pdbid 4 character PDB code
-- @param chainid 1 character chain identiifier
-- @return report string for printing
function chain.reportDbChain(rfpg,pdbid,chainid)
   local s=' ' .. chainid
   local pdbno = rfpg.Q("select pdb_no from pdb_chain where pdbid='" .. pdbid .. "' and chain='" .. chainid .. "'")
   if pdbno then
      pdbno = pdbno[1]
      local rTotal = rfpg.Q("select count(*) from residue where residue.pdb_no=" .. pdbno )[1]
      local rOrdered = rfpg.Q("select count(*) from residue where residue.pdb_no=" .. pdbno .. ' and ordered is true' )[1]
      s = s .. ' ordered residues: ' .. rOrdered .. ' total residues: ' .. rTotal 
      local d = rfpg.Q("select count(dssp.*) from dssp,residue where residue.pdb_no=" .. pdbno .. " and dssp.res_id = residue.res_id")[1]
      s = s .. ' DSSPs: ' .. d
      local dh = rfpg.Q("select count(dihedral.*) from dihedral,residue where residue.pdb_no=" .. pdbno .. " and dihedral.res_id = residue.res_id")[1]
      s = s .. ' dihedra: ' .. dh 
      --local h = rfpg.Q("select count(angle.*) from angle,residue where residue.pdb_no=" .. pdbno .. " and angle.res_id=residue.res_id")[1]
      local h = rfpg.Q('select count(*) from (select distinct angle.key from angle, residue where residue.pdb_no=' .. pdbno .. ' and angle.res_id=residue.res_id) as c')[1]
      s = s .. ' hedra: ' .. h

      --[[
         select count(*) from (
   select distinct angle.key 
    from angle, dihedral, residue 
    where residue.pdb_no=6 and dihedral.res_id=residue.res_id and angle.dihedral_id = dihedral.dihedral_id
    ) as c;
      --]]
      --local ac = rfpg.Q("select count(atom_coordinates.*) from atom_coordinates,residue where residue.pdb_no=" .. pdbno .. " and atom_coordinates.res_id = residue.res_id")[1]
      --s = s .. ac .. ' atom coordinates ' 
      local cs = rfpg.Q("select count(atom_coordinates.*) from atom_coordinates,residue where residue.pdb_no=" .. pdbno .. " and atom_coordinates.res_id = residue.res_id and init is true")[1] /3
      s = s .. ' initial coordinate sets: ' .. cs
      local fp = rfpg.Q('select res_pos from residue where pdb_no=' .. pdbno .. ' order by res_ndx limit 1')[1]
      s = s .. ' firstPos: ' .. fp
   end

   return s
end


-- class methods follow


--- find or initialise Residue object (ordered by definition) for specified amino acid and sequence position
-- will set firstPos (first residue in this chain) if appropriate (if not set already)
-- else any missing residues between firstPos and this one will be set to {} (disordered)
-- @param ires uppercase 1-letter amino acid code
-- @param iresn sequence position of residue in this chain
-- @param disordered optional boolean to indicate disordered residue (default is ordered)
-- @return Residue object, minimally initialised and referenced from self['residues'][iresn]
function Chain:getResidue( ires, iresn, disordered )
   if 'X' == ires then return nil end   -- dssp missing/disordered residue
   if disordered then 
      self['residues'][iresn] = {}
      if not self['firstPos'] then self['firstPos'] = iresn end
   else
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
   return ordResIter, self['residues'],self['firstPos']-1
end

--- count ordered residues in chain according to orderedResidues()
-- @return count of ordered residues
function Chain:countOrderedResidues()
   local c=0
   for i,r in self:orderedResidues() do
      c=c+1
   end
   return c
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
   if not self['firstPos'] then return nil end
   return allResIter, self['residues'], self['firstPos']-1
end

--- count all residues in chain according to allResidues()
-- @return count of all residues
function Chain:countAllResidues()
   local c=0
   for i,r in self:allResidues() do
      c=c+1
   end
   return c
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
   if nil ~= next(self['initNCaC']) then
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
-- @see Chain:setStartCoords()
-- @see Chain:writeSCAD()
-- @param start optional int sequence position start of range
-- @param fin optional int sequence position end of range
function Chain:assembleResidues(start,fin)
   local c=1
   local ndx=0
   for i,r in self:orderedResidues() do
      if (not start or start <= i) and (not fin or fin >= i) then      
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
end

--- trigger each Residue in Chain to generate PDB ATOM records, add TER record at end of chain
-- @param ndx counter for lines in pdb output file
-- @param range optional begin:end filter for residues to print
-- @return string of PDB format records for ATOMS in Chain, plus TER record
function Chain:writePDB(ndx,range)
   --print('chain topdb '  .. ndx)
   local s = ''
   local ls
   local lastRes
   local start
   local fin
   if range then
      start,fin = range:match('(%d+):(%d+)')
      start = tonumber(start)
      fin = tonumber(fin)
   end
   for i,r in self:orderedResidues() do
      if (not start or start <= i) and (not fin or fin >= i) then
         ls,ndx = r:writePDB(self['id'],ndx)
         s = s .. ls
         lastRes=i
      end
   end
   local res = self['residues'][lastRes]
   ndx = ndx + 1
   s = s .. string.format('TER   %5d      %3s %1s%4d%54s\n',ndx, chemdata.res3[ res['res'] ], self['id'], res['resn'], '')

   return s,ndx
end

--- generate PIC hedra and dihedra records for each residue in chain
-- @param range optional begin:end filter for residues to print
-- @return string of PIC format records for residues in Chain
function Chain:writeInternalCoords(range)
   local s = ''

   if nil ~= next(self['initNCaC']) then
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

   local start
   local fin
   if range then
      start,fin = range:match('(%d+):(%d+)')
      start = tonumber(start)
      fin = tonumber(fin)
      print('start',start,'fin',fin)
   end

   for i,r in self:orderedResidues() do
      if (not start or start <= i) and (not fin or fin >= i) then
         s = s .. r:writeInternalCoords(self['pdbid'], self['id'])
      end
   end
   return s
end

--- write chain data to rfold database
-- @param rfpg open database handle
-- @param pdb_no chain ID in pdb_chain table
-- @param update optional flag, if false silently skip if [pdb_no, res_ndx] entry exists already in residue table 
function Chain:writeDb(rfpg, pdb_no, update)

   local ndx=0
   for i,r in self:allResidues() do   -- note i = res_pos = varchar (from pdb so could be weird)
      --print(r:tostring())
      ndx = ndx+1
      local res_id = rfpg.Q("select res_id from residue where pdb_no=" .. pdb_no .. " and res_pos = '" .. i .. "'" )
      if not res_id then
         res_id = rfpg.Q( "insert into residue (res_ndx,pdb_no) values (" .. ndx .. ',' .. pdb_no .. ") returning res_id" )
      elseif not update then
         goto skipResidue  -- can't happen? because existing chain would have been skipped
      end
      res_id = res_id[1]

      if r['ordered'] then  -- allResidues() iterator above includes disordered residues
         local qc = "update residue set (res, res_pos, ordered"
         local qv = " = ('" .. r['res'] .. "','" .. r['resn'] .. "'" .. (r['ordered'] and ', true' or ', NULL')
         if r['prev'] and r['prev']['res'] then
            qc = qc .. ", prev_res"
            qv = qv .. ",'" .. r['prev']['res'] .. "'"
         end
         if r['next'] and r['next']['res'] then
            qc = qc .. ", next_res"
            qv = qv .. ",'" .. r['next']['res'] .. "'"
         end
      
         rfpg.Qcur( qc .. ')' .. qv .. ') where res_id = ' .. res_id )
      end
      
      if (nil ~= next(self['initNCaC']) and self['initNCaC'][i]) then  -- have init coords for chain segment, need to store
         local akl =  r:NCaCKeySplit()
         local initCoords = self['initNCaC'][i]
         for ai,ak in ipairs(akl) do
            local q
            local v = initCoords[ak]
            local haveRow = rfpg.Q("select 1 from atom_coordinates where res_id=" .. res_id .. " and atom = '" .. ak .. "'")
            if (not haveRow) then
               q = 'insert into atom_coordinates (res_id, atom, init, x, y, z) values (' .. res_id .. ",'" .. ak .. "', true,"
            elseif update then
               q = 'update atom_coordinates set (init, x, y, z) = ( true,'
            end
            if q then -- ((not haveRow) or update)
               for j=1,3 do q = q .. (j>1 and ',' or '')  .. v[j][1] end
               q = q .. ')'
               if (haveRow) then  -- update
                  q = q .. ' where res_id=' .. res_id .. " and atom= '" .. ak .."'"
               end
               rfpg.Qcur(q)
               
            end
         end
      end

      if r['ordered'] then r:writeDb(rfpg, res_id, update) end
      
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

--- populate chain object from database using supplied pdb and chain id 
-- @param rfpg open database handle
function Chain:dbLoad(rfpg)
   local pdbno = self['pdb_no'] and self['pdb_no'] or rfpg.Q("select pdb_no from pdb_chain where pdbid='" .. self['pdbid'] .. "' and chain= '" .. self['id'] .. "'")[1]
   --print(self['pdbid'],self['id'],pdbno)
   local cur = rfpg.Qcur('select res_id, res, res_pos, ordered from residue where pdb_no = ' .. pdbno .. ' order by res_ndx')
   local resData = cur:fetch({},'a')
   while resData do
      local resPos = tonumber(resData['res_pos'])
      --print(resPos, resData['res'], resData['ordered'], resData['res_id'])
      if resPos then   -- if not disordered
         self['residues'][resPos] = self:getResidue(resData['res'],resPos,(not resData['ordered']))
         self['residues'][resPos]['res_id'] = resData['res_id']
         
         local cur2 =    rfpg.Qcur('select atom, x, y, z from atom_coordinates where res_id = ' .. resData['res_id'] .. ' and init is true')
         local initCoords = {}
         
         local atomData = cur2:fetch({},'a')
         
         while atomData do
            local atom = geom3d.get41mtx()
            atom[1][1] = tonumber(atomData['x'])
            atom[2][1] = tonumber(atomData['y'])
            atom[3][1] = tonumber(atomData['z'])
            initCoords[ atomData['atom'] ] = atom
            --print(atomData['atom'])
            atomData = cur2:fetch({},'a')
         end

         if nil ~= next(initCoords) then
            self['initNCaC'][resPos] = initCoords
         end
      end
      resData = cur:fetch({},'a')
   end

   for i,r in self:orderedResidues() do
      r:dbLoad(rfpg)
   end
   --[[
   for i,r in ipairs(self['residues']) do
      print(i)
      if next(r) then r:dbLoad(rfpg) end
   end
   --]]
   
   self:linkResidues()
end   


--- query database for info about this chain
-- @param rfpg open database handle
-- @return report string for printing
function Chain:reportDb(rfpg)
   return chain.reportDbChain(rfpg,self['pdbid'],self['id'])
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

local function writeSCADdihed(d,transformations,hedraNdx)
   local s = ''
   s = s .. '    //   ' .. d['key'] .. '[ ' .. d['h1key'] .. '  --  ' .. d['h2key'] .. '\n'
   s = s .. '  [ ' .. d['dihedral1'] .. ', ' .. hedraNdx[d['h1key']]-1 .. ', ' .. hedraNdx[d['h2key']]-1 .. ', ' .. (d['reverse'] and '1' or '0') .. ', '   -- openSCAD arrays start from 0
   s = s .. geom3d.writeSCAD(transformations[d['key3']]) .. ' ]'
   return s
   --[[
   print('backbone dk:',dk)
   --for x,y in pairs(d) do print(x,y) end
   print('key3:',d['key3'],'key32:',d['key32'])
   print('key3 transformations:')
   print(transformations[ d['key3'] ]:transpose():pretty())
   print('dihed1:',d['dihedral1'])
   print('hndx:',hedraNdx[ d['key3'] ] .. ', ' .. hedraNdx[ d['key32'] ])
   --]]
end   

--- generate string of OpenSCAD commands to render chain
-- @param range start:finish filter for residues to output
-- @param backboneOnly optional boolean default false, only output for backbone atoms
-- @return string of OpenSCAD data and command lines
function Chain:writeSCAD(range,backboneOnly)
   local start
   local fin
   local s=''
   if range then
      start,fin = range:match('(%d+):(%d+)')
      start = tonumber(start)
      fin = tonumber(fin)
      --print('start',start,'fin',fin)
   end

   -- collect all hedra in hash to eliminate redundant references
   local hedra = {}
   for i,r in self:orderedResidues() do
      if (not start or start <= i) and (not fin or fin >= i) then
         local resHedra = r['hedra']
         for k,h in pairs(resHedra) do
            hedra[k] = h
         end
      end
   end
-- [[
   -- generate openSCAD array of hedra plus index table so can refer to each by array index
   local hedraNdx = {}
   local hedraSCAD = ''
   local ndx=0
   for k,h in utils.pairsByKeys(hedra, residue.keysort) do
      ndx = ndx + 1
      hedraNdx[k] = ndx
      local atoms = {}
      for i,k in pairs(utils.splitKey(k))  do
         local a = utils.splitAtomKey(k)
         atoms[i] = a[3]:sub(1,1)  -- get just the element e.g. N, C, O
      end
      hedraSCAD = hedraSCAD .. '    // ' .. h['key'] .. '\n'
      hedraSCAD = hedraSCAD .. string.format(' [ %9.5f, %9.5f, %9.5f, "%s", "%s", "%s" ],\n', h['len1'], h['angle2'], h['len3'], atoms[1], atoms[2], atoms[3])
   end
--]]

   -- for x,y in pairs(hedraNdx) do print(x,y) end
   
   -- generate openSCAD array of dihedra grouped by residue, each dihedron with residue-space transformation matrix

   local resNdx = {}
   local resSCAD = ''
   ndx = 0
   for i,r in self:orderedResidues() do
      if (not start or start <= i) and (not fin or fin >= i) then
         ndx = ndx + 1
         resNdx[r['resn']] = ndx
         resSCAD = resSCAD .. ' [  // ' .. r['res'] .. r['resn'] .. ' backbone \n'
         local transformations = r:assemble(nil,true)  -- assemble with no start position, return transformation matrices
         --for x,y in pairs(transformations) do print(x, y:transpose():pretty()) end
         --for x,y in pairs(transformations) do print(x) end
         local first=true
         for dk,d in pairs(r['dihedra']) do
            --print('dk:',dk)
            if d:isBackbone() and hedraNdx[d['h2key']] then
               if first then first=false else resSCAD = resSCAD .. ',\n' end
               resSCAD = resSCAD .. writeSCADdihed(d,transformations,hedraNdx)
            end
         end
         if not backboneOnly then
            resSCAD = resSCAD .. ',\n    // ' .. r['res'] .. r['resn'] .. ' sidechain \n'
            first = true
            for dk,d in pairs(r['dihedra']) do
               --print('dk:',dk)
               if (not d:isBackbone()) and hedraNdx[d['h2key']] then
                  if first then first=false else resSCAD = resSCAD .. ',\n' end
                  resSCAD = resSCAD .. writeSCADdihed(d,transformations,hedraNdx)
               end
            end
         end

         resSCAD = resSCAD .. '\n ],\n'

      end
   end

   resSCAD = resSCAD:sub(1,-3) -- lose last ','

   self:clearAtomCoords()   -- all chain fragments start at default dihedron coordinates
   self['initNCaC'] = {}    -- ignore world offset for PDB file
   self:assembleResidues(start,fin)   -- generate atom coordinates starting from origin (default init pos for first dihedron) - any chain breaks handled manually

   local chainSCAD = ''
   local first = true
   for i,r in self:orderedResidues() do
      if (not start or start <= i) and (not fin or fin >= i) then
         local mt, mtr
         local akl = r:NCaCKeySplit()
         if r['prev'] and r['prev']['ordered'] then
            local atomCoords = r['prev']['atomCoords']
            mt, mtr = geom3d.coordSpace( atomCoords[akl[1]], atomCoords[akl[2]], atomCoords[akl[3]], true ) -- get transforms to, from coord space of NCaC for this residue in world coords
         else
            mtr = geom3d.get44mtx()
         end

         if first then first=false else chainSCAD = chainSCAD .. ',\n' end
         chainSCAD = chainSCAD .. ' [ ' .. resNdx[r['resn']]-1 .. ',   // ' .. r['res'] .. r['resn'] .. '\n'
         chainSCAD = chainSCAD .. '   ' .. geom3d.writeSCAD(mtr) .. ' ]'
      end
   end

   chainSCAD = chainSCAD .. '\n'
      
   s = s .. 'hedra = [\n' .. hedraSCAD:sub(1,-3) .. '\n];\n\n'
   s = s .. 'residues = [\n' .. resSCAD .. '\n];\n'
   s = s .. 'chains = [\n' .. chainSCAD .. '\n];\n'
   
   return s
end

return chain
