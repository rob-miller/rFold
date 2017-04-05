--[[
   protein.lua
   
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

--- protein structure class for rFold 
--
-- @classmod Protein

local utils = require 'rfold.utils'
local chain = require 'rfold.chain'

local protein = {}  -- module

local Protein = {}  -- class table

protein.proteins = {} -- list of loaded proteins to find by pdbid

------------------------------------------------------------------------------------------------
-- Protein
------------------------------------------------------------------------------------------------
-- class properties

--- 4 position PDB identifier for this protein
-- @field id string

--- array of chain objects indexed by chain ID
-- @field chains array of chain objects 

--- array of chain IDs indexed by integer order of chains in PDB file
-- @field chainOrder array of char chain IDs

--- pdb HEADER record if supplied
-- @field header string

--- pdb TITLE record if supplied
-- @field title string

--- pdb TITLE record if supplied
-- @field title string

--- pdb COMPND record if supplied
-- @field compnd string 

--- location of file this protein data was read from
-- @field filename string


--- Protein class object initialiser (not a class method)
-- @param o table with field 'id' = Brookhaven / RCSB PDB alphanumeric ID or similar
-- @return minimally initialised Protein object, object added to proteins list
function protein.new (o)
   --for k,v in pairs(o) do print(k,v) end
   assert((o and o['id']),'protein.new() called without id')
   setmetatable(o, { __index = Protein })

   if not o['chains'] then
      o['chains'] = {}
   end
   if not o['chainOrder'] then
      o['chainOrder'] = {}
   end

   if not o['header'] then
      o['header'] = string.format('HEADER    %-40s%9s   %s','DE NOVO PROTEIN',os.date("%d-%b-%Y"),o['id'])
      --print(o['header'])
   end
   protein.proteins[o['id']] = o
   
   return o
end

--- find the Protein with this id or new Protein if it does not exist (not a class method)
-- @param id Brookhaven / RCSB PDB alphanumeric ID or similar
-- @return already existing or minimally initialised protein object
function protein.get(id)
   if protein.proteins[id] then
      return protein.proteins[id]
   else 
      return protein.new{ ['id'] = id }
   end
end

--- generate a new key for the data stored in protein.proteins for the passed id key.
-- <br> new key of form <id>-<n> where n is number uch that <id>-<n> is available
-- <br> enables re-loading file for line-by-line compare to generated daa for test (probably need prot:setStartCoords())
-- @param id Brookhaven / RCSB PDB alphanumeric ID or similar
-- @return new key or nil if no entry exists for passed id
function protein.stashId(id)
   if not protein.proteins[id] then return nil end
   local n=1
   local s1 = id .. '-' .. n
   while protein.proteins[s1] do
      n = n+1
      s1 = id .. '-' .. n
   end
   protein.proteins[s1] = protein.proteins[id]
   protein.proteins[id] = nil
   return s1
end

--- change the key for an entry in the protein.proteins table
-- @param old old key
-- @param new new key
function protein.rename(old,new)
   assert(protein.proteins[old], 'no protein with id ' .. old .. ' loaded')
   protein.proteins[new] = protein.proteins[old]
   protein.drop(old)
end

--- remove an entry from the protein.proteins table so garbage collector can free memory
-- @param id key to drop from table
function protein.drop(id)
   if id then protein.proteins[id] = nil end
end

--- remove all entries from the protein.proteins table so garbage collector can free memory
function protein.clean()
   protein.proteins = {}
end

                      
--- add the passed table data to the protein specified (pdbid) in the table (not a class method - may create new Protein)
-- @param t table with field 'pdbid' and preferably chain identifier field 'chn', passed here as callback from file parser
function protein.load(t)
   if not t['pdbid'] then
      for k,v in pairs(t) do print(k,v) end
      assert(nil,'protein.load - no pdbid for ' .. tostring(t))
   end

   local thisProtein = protein.get(t['pdbid'])

   --for k,v in pairs(t) do print('p: ' .. k,v) end

   if t['header'] then
      thisProtein['header'] = t['header'] --  .. t['pdbid']
   elseif t['title'] then
      if not thisProtein['title'] then thisProtein['title'] = t['title'] end
   elseif t['compnd'] then
      thisProtein['compnd'] = t['compnd']
   elseif t['filename'] then
      thisProtein['filename'] = t['filename']
   else
      local chn = t['chn'] or ''
      if not thisProtein['chains'][chn] then
         thisProtein['chains'][chn] = chain.new({ ['id'] = chn, ['pdbid'] = t['pdbid'] })
         thisProtein['chainOrder'][#thisProtein['chainOrder']+1] = chn
      end
      thisProtein['chains'][chn]:load(t)
   end
end

--- copy DSSP data from src protein object to dst protein object (not a class method)
-- will fall over if protein objects do not otherwise match
-- @param src Protein object with DSSP data
-- @param dst Protein object to get DSSP data
function protein.copyDSSP(src,dst)
   for cid,c in pairs(src['chains']) do
      for resn,res in pairs(c['residues']) do
         if res['dssp'] then
            dst['chains'][cid]['residues'][resn]['dssp'] = res['dssp']
         end
      end
   end
end

---- class methods follow ----

--- compute atomCoords for already loaded internal coordinates; wrapper grouping several Protein: class methods
function Protein:internalToAtomCoords()
   self:linkResidues()              -- set prev, next for each residue in each chain; create tables from residue to dihedron positions according to PDB backbone, sidechain atoms
   self:renderDihedra()             -- generate hedron and dihedron space coordinates for hedron angle, lengths and dihedron angle
   
   self:assembleResidues()          -- assemble overlapping dihedra of each residue to complete atomCoords list, copy coordinates for i+1 N, CA, C to next residue 
   
   --return self:PDB()              -- output PDB format text
end


--- return table of component counts for chains in protein similar to tostring()
-- @return table of [chain id][residue count, dssp count, dihedron count, hedron count]
function Protein:stats()
   local rt = {}
   for k,v in pairs(self['chains']) do
      rt[k] = { #v['residues'], v:countDSSPs(), v:countDihedra(), v:countHedra() }
   end
   return rt
end

--- sum dihedron counts for chains
-- @return total number of dihedra loaded for protein
function Protein:countDihedra()
   local dc = 0
   for k,v in pairs(self['chains']) do
      dc = dc + v:countDihedra()
   end
   return dc
end

--- sum dssp counts for chains
-- @return total number of dssp entries loaded for protein
function Protein:countDSSPs()
   local dc = 0
   for k,v in pairs(self['chains']) do
      dc = dc + v:countDSSPs()
   end
   return dc
end

--- generate sequence strings for chains in this protein
-- @return chain sequence strings, 1n separated
function Protein:tosequence()
   local s=''
   for k,v in ipairs(self['chainOrder']) do
      local chain = self['chains'][v]
      s = s .. self.id .. ' ' .. v .. ' ' .. chain:seqStr() .. '\n'
   end
   return s
end

--- generate descriptive string for Protein : id, sequence, counts of residues, DSSP entries, hedra and dihedra
-- @return descriptive string
function Protein:tostring()
   local s = self['id'] .. ' chains: ' .. table.getn(self['chainOrder']) .. '\n'
   for k,v in ipairs(self['chainOrder']) do
      local chain = self['chains'][v]
      if next(chain) then -- possibly empty if chain specified for load
         s = s .. ' ' .. v  
         local cs = chain:countInitNCaCs()
         local ordRes = chain:countOrderedResidues()
         local allRes = chain:countAllResidues()
         s = s .. '  ordered residues: ' .. ordRes .. ' total residues: ' .. allRes .. ' DSSPs: '  .. chain:countDSSPs() .. ' dihedra: ' .. chain:countDihedra() .. ' hedra: ' .. chain:countHedra() ..  ' initial coordinate sets: ' .. cs .. ' firstPos: ' .. chain['firstPos'] .. '\n'
      end
   end
   return s
end

--- generate PDB format text data, adding HEADER, CAVEAT and END records to Chain:writePDB() result
--  note Residue structures need dihedron data structures to write complete chain
-- @param noRemark boolean do not write REMARK records indicating this is RFOLD output
-- @param range start:fin filter for residues to output 
-- @return string containing PDB format text (complete)
function Protein:writePDB(noRemark,range)
   --print('protein topdb ')
   local s= ''
   if self['header'] then s = s .. self['header'] ..'\n' end
   if self['title'] then s = s .. self['title'] ..'\n' end
   if not noRemark then
      s = s .. string.format('%-80s\n','REMARK   5')
      s = s .. string.format('%-80s\n','REMARK   5 WARNING')
      s = s .. string.format('%-80s\n','REMARK   5 ' .. self['id'] .. ':   ' .. 'RFOLD OUTPUT ' .. os.date("%d-%b-%Y %H:%M:%S") )
   end
   
   local ndx=0
   local ls
   for k,v in ipairs(self['chainOrder']) do
      local chain = self['chains'][v]
      ls, ndx = chain:writePDB(ndx,range)
      s = s .. ls
   end
   s = s .. string.format('%-80s','END')
   return s
end

--- generate text data specifying hedra with length, angle, length ; dihedra with angle ; relevant PDB records for annotation and start coordinates
-- @param noTitle do not write TITLE line if true; DSSP does not output so compare fails
-- @param range start:fin filter for residues to output 
-- @return string containing structure specification in internal coordinates (complete)
function Protein:writeInternalCoords(noTitle,range)
   local s = ''
   if self['header'] then s = s .. self['header'] ..'\n' end
   if (not noTitle) and self['title'] then s = s .. self['title'] ..'\n' end
   for k,v in ipairs(self['chainOrder']) do
      local chain = self['chains'][v]
      s = s .. chain:writeInternalCoords(range)
   end
   return s
end

--- write protein data to rfold database
-- @param rfpg open database handle
-- @param update optional flag, if false (default at protein level only) silently skip if [pdbid,chain,filename] entry exists already in pdb_chain
-- @param src string to put in pdb_chain:src column
function Protein:writeDb(rfpg,update,src)
   update = update or false
   src = src or false

   --[[
      -- takes longer (47s for update code vs 1:15 for cascade delete then load for 1mud)
   if update then
      rfpg.Qcur("delete from pdb_chain where pdbid='" .. self['id'] .. "'")
      update = false
   end
   --]]

   --print('pdbid', self['id'], update)
   local first=true
   local tct = "','"
   local maxChain = rfpg.Q("select max(chain_order) from pdb_chain where pdbid='" .. self['id'] .. "' and filename = '" .. self['filename'] .. "'")
   if maxChain then maxChain = maxChain[1] end
   if not maxChain then maxChain=0 end
   --print('maxChain=', maxChain, self['id'], self['filename'])
   for k,v in ipairs(self['chainOrder']) do
      local pdb_no = rfpg.Q("select pdb_no from pdb_chain where pdbid='" .. self['id'] .. "' and chain = '" .. v .. "' and filename = '" .. self['filename'] .. "'")
      if not pdb_no then 
         pdb_no  = rfpg.Q("insert into pdb_chain (pdbid, chain, filename, chain_order) values ('".. self['id'] .. tct .. v .. tct .. self['filename'] .. "'," .. k+maxChain .. ") returning pdb_no")[1]
         if src then
            rfpg.Qcur("update pdb_chain set (src) = ('" .. src .. "') where pdb_no = " .. pdb_no )
         end
      elseif not update then
         goto skipChain
      else
         pdb_no = pdb_no[1]
         rfpg.Qcur("update pdb_chain set (chain_order) = (" .. k .. ") where pdb_no = " .. pdb_no )
         if src then
            rfpg.Qcur("update pdb_chain set (src) = ('" .. src .. "') where pdb_no = " .. pdb_no )
         end
      end

      --print('pdbno', pdb_no, 'pdbid', self['id'], 'chain', v)
      if first then
         if self['header'] then rfpg.Qcur("update pdb_chain set (header) = ('".. self['header'] .. "') where pdb_no = " .. pdb_no) end
         if self['title'] then rfpg.Qcur("update pdb_chain set (title) = ('".. self['title'] .. "') where pdb_no = " .. pdb_no) end
         first = false
      end

      local chain = self['chains'][v]
      
      local seqStr = chain:seqStr()   -- not needed?
      local seqRes = chain['seqres'] or ''

      rfpg.Qcur("update pdb_chain set (sequence, seqres) = ('" .. seqStr .. tct .. seqRes.. "') where pdb_no = " .. pdb_no)

      self['chains'][v]:writeDb(rfpg, pdb_no, update)

      ::skipChain::
   end
   
end

--- query database for info about this protein
-- @param rfpg open database handle
-- @return report string for printing
function Protein:reportDb(rfpg,s)
   local s = 'rFold DB: ' .. self['id'] .. ' '
   local n = rfpg.Q("select distinct pdb_no from pdb_chain where pdbid='" .. self['id'] .. "' order by pdb_no")
   local c = rfpg.Q("select count(*) from pdb_chain where pdbid='" .. self['id'] .. "'")[1]
   s = s .. 'chains: ' .. c .. ' pdb_no: '
   for i,j in ipairs(n) do
      s = s .. j .. ' '
   end
   s = s .. '\n'

   for k,v in ipairs(self['chainOrder']) do
      local chain = self['chains'][v]
      s = s .. ' ' .. chain:reportDb(rfpg) .. '\n'
   end
   
   return s
end

--- populate protein object from database using supplied pdb id
-- @param rfpg open database handle
function Protein:dbLoad(rfpg)
   local cur
   if self['chain'] then
      cur = rfpg.Qcur("select pdb_no, chain, filename, seqres, 1 as chain_order, header, title from pdb_chain where pdbid='" .. self['id'] .. "' and chain='" .. self['chain'] .. "'")
   else
      cur = rfpg.Qcur("select pdb_no, chain, filename, seqres, chain_order, header, title from pdb_chain where pdbid='" .. self['id'] .. "' order by chain_order")
   end
   local chainData = cur:fetch({},'a')
   local ndx=1
   while chainData do
      --print(chainData['filename'], chainData['chain'])
      self['chainOrder'][tonumber(chainData['chain_order'])] = chainData['chain']
      self['filename'] = chainData['filename']
      if chainData['header'] then self['header'] = chainData['header'] end
      if chainData['title'] then self['title'] = chainData['title'] end
      self['chains'][chainData['chain']] = chain.new({ ['pdbid'] = self['id'],['id'] = chainData['chain'],['seqres'] = chainData['seqres'],['pdb_no'] = chainData['pdb_no'] })
      chainData = cur:fetch({},'a')
   end
   
   for k,v in ipairs(self['chainOrder']) do
      self['chains'][v]:dbLoad(rfpg)
   end
   
end

--- trigger linking of residues, dihedrons and hedrons within chains
function Protein:linkResidues()
   for k,v in pairs(self['chains']) do
      --for x,y in pairs(v) do print(k,x,y) end
      v:linkResidues()
   end
end

--- trigger update of initial atom coordinates (dihedron coordinate space) in each dihedron according to dihedron angle
function Protein:renderDihedra()
   for k,v in pairs(self['chains']) do
      v:renderDihedra()
   end
end

--- delete protein space atom coordinate data throughout protein object
function Protein:clearAtomCoords()
   for k,v in pairs(self['chains']) do
      v:clearAtomCoords()
   end
end
   
--- delete internal coordinate data throughout protein object
function Protein:clearInternalCoords()
   for k,v in pairs(self['chains']) do
      v:clearInternalCoords()
   end
end
   
--- inverse of linkResidues() and renderDihedra() : from input atom coordinates, build complete hedra and dihedra 
function Protein:atomsToInternalCoords()
   for k,v in pairs(self['chains']) do
      v:dihedraFromAtoms()
   end
end

--- trigger calculation/update of residue['backbone'] and residue['sidechain'] atom coordinates (protein coordinate space) by sequentially assembling dihedrons
function Protein:assembleResidues()
   for k,v in pairs(self['chains']) do
      v:assembleResidues()
   end
end

--- set N, CA, C coordinates (protein coordinate space) for first residue in each chain from DSSP data for residue if available
function Protein:setStartCoords()
   for k,v in pairs(self['chains']) do
      v:setStartCoords()
   end
end

--- generate string indicating chain IDs and residues in each
-- @return string
function Protein:report()
   local s=''
   for k,v in ipairs(self['chainOrder']) do
      local chain = self['chains'][v]
      local ordered,all = chain:countResidues()
      s = s .. 'chain ' .. v .. ' ' .. ordered ..'/' .. all .. ' residues  start= ' .. chain['firstPos'] 
   end
   return s
end
      
function Protein:printInfo()
   for k,v in pairs(self['chains']) do
      v:printInfo()
   end
end

   
return protein
