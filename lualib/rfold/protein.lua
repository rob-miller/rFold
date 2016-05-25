
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

--- add the passed table data to the protein specified (pdbid) in the table (not a class method - may create new Protein)
-- @param t table with field 'pdbid' and preferably chain identifier field 'chn', passed here as callback from file parser
function protein.load(t)
   if not t['pdbid'] then
      for k,v in pairs(t) do print(k,v) end
      assert(nil,'protein.load - no pdbid for ' .. tostring(t))
   end

   local thisProtein = protein.get(t['pdbid'])

   if t['header'] then
      thisProtein['header'] = t['header'] --  .. t['pdbid']
   elseif t['title'] then
      thisProtein['title'] = t['title']
   else
      local chn = t['chn'] or ''
      if not thisProtein['chains'][chn] then
         thisProtein['chains'][chn] = chain.new{ ['id'] = chn, ['pdbid'] = t['pdbid'] }
         thisProtein['chainOrder'][#thisProtein['chainOrder']+1] = chn
      end
      thisProtein['chains'][chn]:load(t)
   end
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

--- generate descriptive string for Protein : id, sequence, counts of residues, DSSP entries, hedra and dihedra
-- @return descriptive string
function Protein:tostring()
   local s = ''
   for k,v in pairs(self['chains']) do
      s = self.id .. ' ' .. k .. ' ' .. v:seqStr() .. '\n'
      s = s .. '  ' .. #v['residues'] .. ' residues '  .. v:countDSSPs() .. ' dssps ' .. v:countDihedra() .. ' dihedra' .. v:countHedra() ..  ' hedra '
   end
   return s
end

--- generate PDB format text data, adding HEADER, CAVEAT and END records to Chain:toPDB() result
-- @return string containing PDB format text (complete)
function Protein:toPDB()
   --print('protein topdb ')   
   local s= ''
   if self['header'] then s = s .. self['header'] ..'\n' end
   if self['title'] then s = s .. self['title'] ..'\n' end
   s = s .. string.format('%-80s\n','REMARK   5')
   s = s .. string.format('%-80s\n','REMARK   5 WARNING')
   s = s .. string.format('%-80s\n','REMARK   5 ' .. self['id'] .. ':   ' .. 'RFOLD OUTPUT ' .. os.date("%d-%b-%Y %H:%M:%S") )

   local ndx=0
   local ls
   for k,v in ipairs(self['chainOrder']) do
      local chain = self['chains'][v]
      ls, ndx = chain:toPDB(ndx)
      s = s .. ls
   end
   s = s .. string.format('%-80s','END')
   return s
end

--- generate text data specifying hedra with length, angle, length ; dihedra with angle ; relevant PDB records for annotation and start coordinates
-- @return string containing structure specification in internal coordinates (complete)
function Protein:toInternalCoords()
   local s = ''
   if self['header'] then s = s .. self['header'] ..'\n' end
   if self['title'] then s = s .. self['title'] ..'\n' end
   for k,v in ipairs(self['chainOrder']) do
      local chain = self['chains'][v]
      s = s .. chain:toInternalCoords()
   end
   return s
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

--- inverse of linkResidues() and renderDihedra() : from input atom coordinates, build complete hedra and dihedra 
function Protein:dihedraFromAtoms()
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


   
return protein
