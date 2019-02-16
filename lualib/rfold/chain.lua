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
      s = s .. ' firstPos: ' .. (fp or 'nil')
   end

   return s
end


-- class methods follow


--- find or initialise Residue object (ordered by definition) for specified amino acid and sequence position
-- will set firstPos (first residue in this chain) if appropriate (if not set already)
-- else any missing residues between firstPos and this one will be set to {} (disordered)
-- @param ires uppercase 1-letter amino acid code
-- @param iresn sequence position of residue in this chain
-- @param isDisordered optional boolean to indicate disordered residue (default is ordered)
-- @param noCheckRes optional boolean to disable checking previous residue for ordered state according to atomCoords
-- @return Residue object, minimally initialised and referenced from self['residues'][iresn]
function Chain:getResidue( ires, iresn, isDisordered, noCheckRes )
   --if 'X' == ires then return nil end   -- dssp missing/disordered residue [ 5.v.17 ??? should be '!' ... some case on db load ? ]
   --[[
   if disordered then 
      self['residues'][iresn] = {}
      if not self['firstPos'] then self['firstPos'] = iresn end
   else
   --]]
   --print(iresn,ires)
   if (ires and (not self['residues'][iresn]) or (not self['residues'][iresn]['ordered']))  then   -- have input residue and currently no residue or only placeholder at this position
         self['residues'][iresn] = residue.new{ resn = iresn, res = ires, ordered = (not isDisordered) }
         -- lua needs arrays starting at 1 with no gaps, so fill in as needed here with ' ' instead of residue opbject
         if not self['firstPos'] then
            self['firstPos'] = iresn
            --if 1 < iresn then   
            --   for i=1,iresn-1 do self['residues'][i] = {} end
            --end
         elseif not self['residues'][iresn-1] then -- there is a gap before this position (possibly this is atom record indicating start after chain break)
            for i=self['firstPos'],iresn-1 do
               if not self['residues'][i] then self['residues'][i] = {} end
            end
         else
            --print(ires, iresn)
            --for x,y in ipairs(self['residues'][iresn-1]) do print(x,y) end
            if (not noCheckRes) and self['residues'][iresn-1]['ordered'] then
               self['residues'][iresn-1]:checkResidue()  -- have finished previous residue, check that it at least has complete backbone to be considered ordered
            end
         end
   elseif ires and ires ~= self['residues'][iresn]['res'] then   -- have input residue but does not match what we have already
         print(iresn,ires,self['residues'][iresn]['resn'],self['residues'][iresn]['res'])
         -- assert(nil,'stored residue ' .. self['residues'][iresn]['res'] .. iresn .. ' does not match new residue ' .. ires .. iresn)
         print('stored residue ' .. self['residues'][iresn]['res'] .. iresn .. ' does not match new residue ' .. ires .. iresn)
         return nil
   end
   --end
   --print('got ' .. self['residues'][iresn]['resn'] .. ' ' .. self['residues'][iresn]['res'], self['residues'][iresn]['ordered'])
   return self['residues'][iresn]
end

--- find or create Residue for passed data table, pass table to it for loading
-- NB: callback for parseProteinData, line-by-line processing
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
         if not ('X' == t['res'] and t['ss']) then        -- 
            res = self:getResidue(t['res'], tonumber(t['resn']))
         end

         if t['TER'] then
            if not res then
               res = self:getResidue(nil, tonumber(t['resn']))  -- handle 4bhr faulty TER line for chain A
            end
            if res then
               res:checkResidue()    -- check for terminal residue is ordered (has N-CA-C)
               res = nil             -- no further processing for this residue below
            end
         end
         
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

   if res then   -- getResidue can return nil so might be missing
      res:load(t)
      --print('Chain:load ', res['res'] .. res['resn'] .. ' ' .. self:countDihedra())
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

--- generate di/hedra keys for all residues
function Chain:initEmpty()
   for i,r in self:allResidues() do
      r:initEmpty()
   end
   self['residues'][#self['residues']]:initEmpty(true)  -- set last one to terminal residue
   --print('term ', #self['residues'])
   --os.exit()
end

--- scale atom coordinates by parameter
-- @param scale float scale multiplier
function Chain:scaleAtomCoords(scale)
   for i,r in self:allResidues() do
      r:scaleAtomCoords(scale)
   end
end

--- clear initial coordinates for this chain
function Chain:clearInitNCaC()
   self['initNCaC'] = {}
end

--- populate Residue hedra and dihedra values with average results using rFold db residue selection query
-- add hedron and dihedron tags sd, min, max
-- @param rfpg open database handle
-- @param resSelector rFold database query returning res_id, e.g. "select res_id from dssp where struc='H' and struc2=' X S+   '" limits to residues inside alpha helices
function Chain:getDbStats(rfpg, resSelector)
   for i,r in self:allResidues() do
      r:getDbStats(rfpg, resSelector)
   end
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

--- query each Residue for its count of hedra - note there are redundant hedra between residues
-- @return sum of hedra for Residues in chain
function Chain:countHedra()
   local c=0
   for i,r in self:orderedResidues() do
      --print(i,r['resn'], r['res'])
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
   --print('lastRes',lastRes)
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

      -- if r['ordered'] then  -- allResidues() iterator above includes disordered residues
      if r['res'] and r['resn'] then
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
         local q = qc .. ')' .. qv .. ') where res_id = ' .. res_id
         --print(q)
         rfpg.Qcur( q )
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
   local lastPos
   while resData do
      local resPos = tonumber(resData['res_pos'])
      --print(resPos, resData['res'], resData['ordered'], resData['res_id'])
      if resPos then   -- if not disordered or (firstPos and possibly disordered)
         if lastPos and (lastPos < resPos-1) then
            for i = lastPos+1, resPos-1 do
               self['residues'][i] = {}
               --print('insert at ' .. i)
            end
         end
         lastPos = resPos
         local res = self:getResidue(resData['res'],resPos,(not resData['ordered']), true)
         assert(res,'getResidue failed for database load of ' .. resData['res'] .. ' ' .. resPos .. (resData['ordered'] and 'ordered' or 'disordered'))  -- shoudl never happen as database only has clean data
         self['residues'][resPos] = res
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

         --print('Chain:dbLoad resPos ' .. resPos .. ' : ' .. self['residues'][resPos]:tostring() .. ( self['initNCaC'][resPos] and ' has init coords' or '' ) )

      end
      resData = cur:fetch({},'a')
   end

   --[[
   print()
   for i,r in self:allResidues() do
      --print('Chain:dbLoad all ' .. i .. ' : ' .. ( self['residues'][i] and (self['residues'][i]:tostring() .. ( self['initNCaC'][i] and ' has init coords' or '' )) or 'nil' ) )
      io.write('Chain:dbLoad all ' .. i .. ' ' .. type(self['residues'][i]) .. ' : ' )
      for x,y in pairs(self['residues'][i]) do io.write(x .. ' ') end
      print()
   end
   
   print()
   --]]
   
   for i,r in self:orderedResidues() do
      --print('dbLoad ', i, r['res'])
      r:dbLoad(rfpg)
      --print('Chain:dbLoad ordered ' .. i .. ' : ' .. r:tostring())
   end
   
   self:linkResidues()
end   


--- query database for info about this chain
-- @param rfpg open database handle
-- @return report string for printing
function Chain:reportDb(rfpg)
   return chain.reportDbChain(rfpg,self['pdbid'],self['id'])
end


--- if first Residue in Chain has 'dssp' field or has NCaC atoms, set self['initNCaC'] to hold those protein space coordinates to build the rest of the Chain from
-- @see assembleResidues
function Chain:setStartCoords()
   local c=0
   for i,r in self:orderedResidues() do  -- find any residues with no previous or previous is disordered
      --print('setStartCoords: ', i,r['resn'],r['res'],self['residues'][i-1])
      if not (self['residues'][i-1] and self['residues'][i-1]['ordered']) then
         -- print('setStartCoords for ndx ' .. i)
         local akl =  r:NCaCKeySplit()
         local initCoords
         if r['dssp'] then
            initCoords = {}
            for ai,ak in ipairs(akl) do
               initCoords[ak] = r:dsspAtom(ak)
               --print('dssp startCoords: ' .. ak)
            end
         elseif r['atomCoords'][akl[1]] and r['atomCoords'][akl[2]] and r['atomCoords'][akl[3]] then
            initCoords = {}
            for ai,ak in ipairs(akl) do
               initCoords[ak] = r['atomCoords'][ak]
               --print('pdb startCoords: ' .. ak)
            end
            --print('setAtomCoords for ' .. r['resn'] .. r['res'])
         else
            if utils.warn then
               print('no start coordinates set for ' .. r['resn'] .. r['res'] .. ' ' .. akl[1] .. ' ' .. akl[2] .. ' ' .. akl[3] )
            end
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

local function indent(spaces)
   s=''
   for i=1,spaces do s = s .. ' ' end
   return s
end

local function writeSCADdihed(d,transformations,hedraNdx)
   local s = ''
   s = s .. '[ ' .. d['dihedral1'] .. ', ' .. hedraNdx[d['h1key']] .. ', ' .. hedraNdx[d['h2key']] .. ', ' .. (d['reverse'] and '1' or '0') .. ', '   
   s = s .. '    //   ' .. d['key'] .. ' [ ' .. d['h1key'] .. '  --  ' .. d['h2key'] .. ' ] ' .. (d['reverse'] and '(reversed) ' or '') .. (d['name'] and d['name'] or '') .. '\n'
   s = s .. indent(8) .. geom3d.writeSCAD( transformations[d['key3']]  ) .. ' ]'

   --[[
   print('backbone dk:',dk)
   --for x,y in pairs(d) do print(x,y) end
   print('key3:',d['key3'],'key32:',d['key32'])
   print('key3 transformations:')
   print(transformations[ d['key3'] ]:transpose():pretty())
   print('dihed1:',d['dihedral1'])
   print('hndx:',( hedraNdx[ d['key3'] ] and hedraNdx[ d['key3'] ] or 'nil' )  .. ', ' .. ( hedraNdx[ d['key32'] ] and hedraNdx[ d['key32'] ] or 'nil' ))
   --]]

   return s
end   

--- generate string of OpenSCAD format arrays 
-- @param scale slicer units per angstrom (usually millimetres)
-- @param range start:finish filter for residues to output
-- @param backboneOnly optional boolean default false, only output for backbone atoms
-- @return string of OpenSCAD data and command lines
function Chain:writeSCAD(scale,range,backboneOnly)
   local start
   local fin

   if range then
      start,fin = range:match('(%d+):(%d+)')
      start = tonumber(start)
      fin = tonumber(fin)
      --print('start',start,'fin',fin)
   end

   -- collect all hedra in hash to eliminate redundant references
   -- generate map from residue, di/hedra_class to di/hedra_key
   local hedra = {}
   local hedraClass = {}
   for i,r in self:orderedResidues() do
      if (not start or start <= i) and (not fin or fin >= i) then
         if not hedraClass[i] then hedraClass[i] = {} end
         local resHedra = r['hedra']
         for k,h in pairs(resHedra) do
            hedra[k] = h
            hedraClass[i][h['class']] = k
            h['len1'] = h['len1'] -- * scale
            h['len3'] = h['len3'] -- * scale
            h:initPos()
         end
         local resDihedra = r['dihedra']
         for k,d in pairs(resDihedra) do
            hedraClass[i][d['class']] = k
            d:initPos()
         end
      end
   end

   -- generate openSCAD array of hedra plus index table so can refer to each by array index
   local hedraNdx = {}
   local hedraSCAD = ''
   local ndx=0
   for k,h in utils.pairsByKeys(hedra, residue.keysort) do
      hedraNdx[k] = ndx
      ndx = ndx + 1
      local atoms = {}
      for i,k in pairs(utils.splitKey(k))  do
         local a = utils.splitAtomKey(k)
         --atoms[i] = a[3]:sub(1,1)  -- get just the element e.g. N, C, O
         atoms[i] = chemdata.residue_atom_bond_state['X'][a[3]] or chemdata.residue_atom_bond_state[a[2]][a[3]]  
      end
      hedraSCAD = hedraSCAD .. indent(5) .. string.format('[ %9.5f, %9.5f, %9.5f, "%s", "%s", "%s" ], ', h['len1'], h['angle2'], h['len3'], atoms[1], atoms[2], atoms[3])
      hedraSCAD = hedraSCAD .. '    // ' .. h['key'] .. '\n'
   end

   hedraSCAD =  hedraSCAD:sub(1,-2);   -- trim last \n but leaving ',' as harder to trim
    
   -- for x,y in pairs(hedraNdx) do print(x,y) end
   
   -- generate openSCAD array of dihedra grouped by residue, each dihedron with residue-space transformation matrix
   -- generate index tables by residue and dihedron key

   local resNdx = {}
   local dihedraNdx = {}
   local resSCAD = ''
   
   ndx = 0
   for i,r in self:orderedResidues() do
      if (not start or start <= i) and (not fin or fin >= i) then
         resNdx[r['resn']] = ndx

         resSCAD = resSCAD .. indent(5) .. '[ // ' .. ndx .. ': ' .. r['resn'] .. r['res'] .. ' backbone \n'

         ndx = ndx + 1       -- openSCAD arrays index from 0 so post-increment

         local transformations = r:assemble(nil,true)  -- assemble with no start position, return transformation matrices

         --for x,y in pairs(transformations) do print(x, y:transpose():pretty()) end
         --for x,y in pairs(transformations) do print(x) end
         --print()
         
         local ndx2 = 0
         local first=true
         for dk,d in pairs(r['dihedra']) do
            --print('dk:',dk)
            if d:isBackbone() and hedraNdx[d['h2key']] then
               if first then first=false else resSCAD = resSCAD .. ',\n' end
               resSCAD = resSCAD .. indent(6) .. writeSCADdihed(d,transformations,hedraNdx)
               dihedraNdx[dk] = ndx2
               ndx2 = ndx2+1
            end
         end
         if not backboneOnly then
            resSCAD = resSCAD .. ',\n' .. indent(7) .. '// ' .. r['resn'] .. r['res'] .. ' sidechain \n'
            first = true
            for dk,d in pairs(r['dihedra']) do
               --print('dk:',dk)
               if (not d:isBackbone()) and hedraNdx[d['h2key']] then
                  if first then first=false else resSCAD = resSCAD .. ',\n' end
                  resSCAD = resSCAD .. indent(6) .. writeSCADdihed(d,transformations,hedraNdx)
                  dihedraNdx[dk] = ndx2
                  ndx2 = ndx2+1
               end
            end
         end

         resSCAD = resSCAD .. '\n' .. indent(5) .. '],\n'

      end
   end

   resSCAD = resSCAD:sub(1,-3) -- lose last ','

   --[[ -- amide set data redundant and not used
      
   -- generate openSCAD array of indexes to
   -- 1N-1Ca-1C hedron            (alpha carbon hedron)
   -- 1Ca-1C-2N-2Ca dihedron      (omega dihedral)
   -- 1N-1Ca-1C-1O, 2C-2Ca-2N-2H  (dihedrals covering rest of amide plane, will not be fully rendered)

   local amideSCAD = ''
   for i,r in self:orderedResidues() do
      if (not start or start <= i) and (not fin or fin >= i) then
         if hedraClass[ i ][ 'CACNCA' ] and hedraClass[i+1] then

            amideSCAD = amideSCAD .. indent(5) .. '[ ' .. resNdx[r['resn'] ]
            amideSCAD = amideSCAD .. ', ' .. hedraNdx[ hedraClass[ i ][ 'NCAC' ] ] 
            amideSCAD = amideSCAD .. ', ' .. dihedraNdx[ hedraClass[ i ][ 'CACNCA' ] ]
            amideSCAD = amideSCAD .. ', ' .. dihedraNdx[ hedraClass[ i ][ 'NCACO' ] ]

            if ( hedraClass[ i+1 ][ 'CCANH' ] ) then  -- XXXX TODO: FIX: ALERT: Need amide protein for H bonds in secondary structure
               amideSCAD = amideSCAD .. ', ' .. dihedraNdx[ hedraClass[ i+1 ][ 'CCANH' ] ] 
            end
            
            amideSCAD = amideSCAD ..  ',   // xxx ' .. r['resn'] .. r['res'] .. ' ' .. hedraClass[ i ][ 'NCAC' ] .. ' ' .. hedraClass[ i ][ 'CACNCA' ] .. ' ' .. hedraClass[ i ][ 'NCACO' ] .. ' '
            amideSCAD = amideSCAD ..  ( hedraClass[ i+1 ][ 'CCANH' ] and hedraClass[ i+1 ][ 'CCANH' ] or 'no CCANH' ) .. '\n'

            if ( hedraClass[ i+1 ][ 'CCANH' ] ) then
               -- need transformation to add next residue to r, with r at origin
               local atomCoords = r:assemble(nil)
               atomCoords = r['next']:assemble(atomCoords)
               local akl = utils.splitKey(hedraClass[ i+1 ][ 'CCANH' ])
               local mt, mtr = geom3d.coordSpace( atomCoords[akl[1] ], atomCoords[akl[2] ], atomCoords[akl[3] ], true )
               --local mt, mtr = geom3d.coordSpace( atomCoords[akl[3] ], atomCoords[akl[2] ], atomCoords[akl[1] ], true )
               amideSCAD = amideSCAD .. indent(6) .. geom3d.writeSCAD(mtr)
            end
            
            amideSCAD = amideSCAD .. indent(5) .. '],\n'

         end
      end
   end
   amideSCAD = amideSCAD:sub(1,-3) -- lose last ','

   --]]
   
   --self:clearAtomCoords()   -- all chain fragments start at default dihedron coordinates
   --self['initNCaC'] = {}    -- ignore world offset for PDB file
   --self:assembleResidues(start,fin)   -- generate atom coordinates starting from origin (default init pos for first dihedron) - any chain breaks handled manually

   local chainSCAD = ''
   local first = true
   for i,r in self:orderedResidues() do
      if (not start or start <= i) and (not fin or fin >= i) then
         local mt, mtr
         local akl = r:NCaCKeySplit()

         if r['prev'] and r['prev']['ordered'] then
            local atomCoords = r['prev']['atomCoords'] -- r['atomCoords'] -- r['prev']['atomCoords'] -- ... no difference ?????
            mt, mtr = geom3d.coordSpace( atomCoords[akl[1]], atomCoords[akl[2]], atomCoords[akl[3]], true ) -- get transforms to, from coord space of NCaC for this residue in world coords
         else
            mtr = geom3d.get44mtx()
         end

         if first then first=false else chainSCAD = chainSCAD .. ',\n' end
         chainSCAD = chainSCAD .. indent(5) .. '[ ' .. resNdx[r['resn']] .. ',   "' .. r['resn'] .. r['res'] .. '",\n'
         chainSCAD = chainSCAD .. indent(6) .. geom3d.writeSCAD(mtr) .. ' ]'
      end
   end

   --chainSCAD = chainSCAD .. '\n'

   local s= indent(3) .. '"' .. self['id'] .. '",  // chain ID\n'   
   s = s .. indent(3) .. '[  // hedra\n' .. hedraSCAD .. '\n' .. indent(3) .. '],\n'
   s = s .. indent(3) .. '[  // residue array of dihedra\n' .. resSCAD .. '\n' .. indent(3) .. '],\n'
   -- s = s .. indent(3) .. '[  // amideSet - backbone only selection\n' .. amideSCAD .. '\n' .. indent(3) .. '],\n'
   s = s .. indent(3) .. '[  // chain - world transform for each residue\n' .. chainSCAD .. '\n' .. indent(3) .. ']\n'
   
   return s
end

return chain
