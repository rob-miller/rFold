
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
function protein.load (t)
   if not t['pdbid'] then
      for k,v in pairs(t) do print(k,v) end
      assert(nil,'protein.load - no pdbid for ' .. tostring(t))
   end

   local thisProtein = protein.get(t['pdbid'])

   if t['header'] then
      thisProtein['header'] = t['header'] .. t['pdbid']
   else
      local chn = t['chn'] or ''
      if not thisProtein['chains'][chn] then
         thisProtein['chains'][chn] = chain.new{ ['id'] = chn }
         thisProtein['chainOrder'][#thisProtein['chainOrder']+1] = chn
      end
      thisProtein['chains'][chn]:load(t)
   end
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
   local s= self['header'] ..'\n'
   s = s .. string.format('%-80s\n','REMARK   5')
   s = s .. string.format('%-80s\n','REMARK   5 WARNING')
   s = s .. 'REMARK   5 ' .. self['id'] .. ':   ' .. 'RFOLD OUTPUT ' .. os.date("%d-%b-%Y %H:%M:%S") .. '\n'

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

--- trigger linking of residues, dihedrons and hedrons within chains
function Protein:linkResidues()
   for k,v in pairs(self['chains']) do
      v:linkResidues()
   end
end

--- trigger update of initial atom coordinates (dihedron coordinate space) in each dihedron according to dihedron angle
function Protein:renderDihedrons()
   for k,v in pairs(self['chains']) do
      v:renderDihedrons()
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
