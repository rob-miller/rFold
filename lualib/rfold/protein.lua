--[[
   protein.lua
   
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
   protein.proteins[id] = nil
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
      thisProtein['title'] = t['title']
   elseif t['filename'] then
      thisProtein['filename'] = t['filename']
   else
      local chn = t['chn'] or ''
      if not thisProtein['chains'][chn] then
         thisProtein['chains'][chn] = chain.new{ ['id'] = chn, ['pdbid'] = t['pdbid'] }
         thisProtein['chainOrder'][#thisProtein['chainOrder']+1] = chn
      end
      thisProtein['chains'][chn]:load(t)
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

--- generate descriptive string for Protein : id, sequence, counts of residues, DSSP entries, hedra and dihedra
-- @return descriptive string
function Protein:tostring()
   local s = ''
   for k,v in pairs(self['chains']) do
      s = self.id .. ' ' .. k .. ' ' .. v:seqStr() .. '\n'
      s = s .. '  ' .. #v['residues'] .. ' residues '  .. v:countDSSPs() .. ' dssps ' .. v:countDihedra() .. ' dihedra ' .. v:countHedra() ..  ' hedra '
   end
   return s
end

--- generate PDB format text data, adding HEADER, CAVEAT and END records to Chain:writePDB() result
--  note Residue structures need dihedron data structures to write complete chain
-- @return string containing PDB format text (complete)
function Protein:writePDB(noRemark)
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
      ls, ndx = chain:writePDB(ndx)
      s = s .. ls
   end
   s = s .. string.format('%-80s','END')
   return s
end

--- generate text data specifying hedra with length, angle, length ; dihedra with angle ; relevant PDB records for annotation and start coordinates
-- @return string containing structure specification in internal coordinates (complete)
function Protein:writeInternalCoords()
   local s = ''
   if self['header'] then s = s .. self['header'] ..'\n' end
   if self['title'] then s = s .. self['title'] ..'\n' end
   for k,v in ipairs(self['chainOrder']) do
      local chain = self['chains'][v]
      s = s .. chain:writeInternalCoords()
   end
   return s
end

--- write protein data to rfold database
-- @param rfpg open database handle
-- @param update optional flag, if false (default at protein level only) silently skip if [pdbid,chain,filename] entry exists already in pdb_chain 
function Protein:writeDb(rfpg,update)
   update = update or false
   local first=true
   local tct = "','"
   for k,v in ipairs(self['chainOrder']) do
      local pdb_no = rfpg.Q("select pdb_no from pdb_chain where pdbid='" .. self['id'] .. "' and chain = '" .. v .. "' and filename = '" .. self['filename'] .. "'")

      if not pdb_no then 
         pdb_no  = rfpg.Q("insert into pdb_chain (pdbid, chain, filename, chain_order) values ('".. self['id'] .. tct .. v .. tct .. self['filename'] .. "'," .. k .. ") returning pdb_no")
      elseif not update then
         goto skipChain
      else
         rfpg.Qcur("update pdb_chain set (chain_order) = (" .. k .. ") where pdb_no = " .. pdb_no[1] ) 
      end
      pdb_no = pdb_no[1]

      if first then
         if self['header'] then rfpg.Qcur("update pdb_chain set (header) = ('".. self['header'] .. "') where pdb_no = " .. pdb_no) end
         if self['title'] then rfpg.Qcur("update pdb_chain set (title) = ('".. self['title'] .. "') where pdb_no = " .. pdb_no) end
         first = false
      end

      local chain = self['chains'][v]
      
      local seqStr = chain:seqStr()
      local seqRes = chain['seqres'] or ''

      rfpg.Qcur("update pdb_chain set (sequence, seqres) = ('" .. seqStr .. tct .. seqRes.. "') where pdb_no = " .. pdb_no)

      self['chains'][v]:writeDb(rfpg, pdb_no, update)

      ::skipChain::
   end
   
end


--- trigger linking of residues, dihedrons and hedrons within chains
function Protein:linkResidues()
   for k,v in pairs(self['chains']) do
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
      s = s .. 'chain ' .. v .. ' ' .. ordered ..'/' .. all .. ' residues  '
   end
   return s
end
      
function Protein:printInfo()
   for k,v in pairs(self['chains']) do
      v:printInfo()
   end
end

   
return protein
