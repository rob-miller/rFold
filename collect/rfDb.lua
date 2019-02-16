#!/usr/bin/env luajit

--[[
   rfDb.lua
   
Copyright 2017 Robert T. Miller

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

--- rfDb.lua [ options ] [<PDBid[chainId]>] ...
-- 
-- Query rFold database
--
-- options:
-- 
-- -l                                    : list pdb ids, chain ids loaded in database
-- -w = pic|gpdb|<filename(.pic|.gpdb)>  : instead of printing, write output to filename with '.pic'.
--                                        or '.gpdb' (generated PDB) extension 
-- -f=<list of PDB codes>                : read PDB codes and optional chain IDs from specified file
-- 

--- not really finished!!  Just a test fixture for pulling proteins from database



local parsers = require 'rfold.parsers'
local protein = require "rfold.protein"
local utils = require 'rfold.utils'
local chain = require 'rfold.chain'
local residue = require 'rfold.residue'
local chemdata = require 'rfold.chemdata'

local rfpg = require 'rfold.rfPostgresql'  -- autocommit false by default
--rfpg.dbReconnect(0) -- no autocommit


local args = parsers.parseCmdLine(
   '[ <pdbid[chain id]> ] ...',
   '\n retrieve internal coordinates for specified pdbid from rFold database with optional processing.\n',
   {
      ['l'] = 'list:  list available PDBids in database',
      ['s'] = 'statistics: report position statistics over +/- 1 std dev for rs subset below',
      ['a'] = 'average: generate internal coordinates for average ALA residue over rs subset below'
   },
   {
--      ['w'] = 'pic|gpdb|<filename(.pic|.gpdb)> : write specified format to pdbid.<format> or <filename>',
      ['f'] = '<input file> : process IDs listed in <input file> (1 per line) followed by any on command line',
      ['r'] = '<firstResidue:lastResidue> : only output residues within the range (positions) specified ',
      ['rs'] = '"\\\"<residue selector string>\\\"" : rFold db query returning res_id e.g. [-rs="\\\"select res_id from dssp where struc=' .. "'H'" .. " and struc2=' X S+   '" .. '\\\""] - see rFold schema (and watch quotes!)'
   }
)

if (args['a'] or args['s']) and not args['rs'] then
   print("options 'a' and 's' require option 'rs'")
   os.exit()
end

local resCount=0
local resDb
local resSelectQry = 'select res_id from temp_selected_residues'
local remarkString

--resSelectQry = args['rs'][1]

if (args['rs']) then
   rfpg.Qcur('begin')  -- keep the temp table for this session
   rfpg.Qcur('create temp table temp_selected_residues (res_id bigint) on commit drop')
   rfpg.Qcur('insert into temp_selected_residues ' .. args['rs'][1])
   --rfpg.Qcur('with t as ( ' .. args['rs'][1] .. ' ) select * into temp_selected_residues from t') 
   local resCur = rfpg.Qcur('select distinct res from residue where res_id in ( ' .. resSelectQry .. ' )')
   local resRow = resCur:fetch({})
   while resRow do
      resCount = resCount+1
      resDb = resRow[1]
      resRow = resCur:fetch(resRow)
   end
   --rfpg.Qcur('rollback')
   if 1 ~= resCount then resDb = '_' end
  -- print('resCount = ', resCount, 'res= ', resDb)
   if 0 == resCount then
      print('residue selector "' .. args['rs'][1] .. "' gets no results")
      os.exit()
   end
   resCount = rfpg.Q('select count(*) from temp_selected_residues')[1]
   remarkString = utils.remarkString('RFOLD AVERAGE OVER ' .. resCount .. ' RESIDUES')
end


if args['s'] then
   print(remarkString)
   --print(args['rs'][1])
   --print('resDb',resDb)
   local res = residue.new({['res'] = resDb, ['resn'] = 2})
   res:initEmpty()
   res:getDbStats(rfpg,resSelectQry)
   print(res:writeInternalCoords('AVG','A', true))
   
end

if args['a'] then
   local pdbid='AVG'
   local p = protein.get( pdbid )
   local chn = 'A'

   -- generate a 2 residue empty chain, residue 1 wildcard, residue 2 resDb 
   p['chains'][chn] = chain.new({ ['id'] = chn, ['pdbid'] = pdbid })
   p['chainOrder'][#p['chainOrder']+1] = chn
   p['remarks'] = { remarkString }
   --p['chains'][chn]:getResidue(resDb, 1, false, true) 
   for i=1,8 do p['chains'][chn]:getResidue(( i==1 and '_' or resDb), i, false, true) end

   p:initEmpty()
   p:getDbStats(rfpg,resSelectQry)
   io.write(p:writeInternalCoords())
   
end


if args['l'] then
   local cur = rfpg.Qcur('select pdb_no, pdbid, chain, filename from pdb_chain order by pdbid, chain, filename')
   local pdb = cur:fetch({},'a')
   while pdb  do
      print(pdb['pdb_no'] .. ' ' .. pdb['pdbid'] .. ' ' .. pdb['filename'] .. ' : ' .. chain.reportDbChain(rfpg,pdb['pdbid'],pdb['chain']))
      pdb = cur:fetch({},'a')
   end
end



local toProcess={}
if args['f'] then
   for i,f in ipairs(args['f']) do
      local fh=io.open(f)
      for line in fh:lines() do
         table.insert(toProcess, line)
      end
   end
end

for i,a in ipairs(args) do table.insert(toProcess, a) end 

for i,a in ipairs(toProcess) do
   if a:match('^IDs%s') then toProcess[i]='' end     -- header line for Dunbrack list : 'IDs         length Exptl.  resolution  R-factor FreeRvalue'
   if a:match('^#') then toProcess[i]='' end         -- comment line starts with '#'
   if a:ematch('^(%d%w%w%w)(%w?)%s+') then           -- looks like pdbcode with optional chain, read as compressed file from PDB_repository_base
      local pdbid, chn = _1, _2
      toProcess[i] = pdbid .. (chn and chn or '') 
   end
end


--[[

   DSSP strcuture info:
   
   STRUCTURE	This is a complex column containing multiple sub columns. The first column contains a letter indicating the secondary structure assigned to this residue. Valid values are:
   Code	Description
   H	Alpha helix
   B	Beta bridge
   E	Strand
   G	Helix-3
   I	Helix-5
   T	Turn
   S	Bend
   What follows are three column indicating for each of the three helix types (3, 4 and 5) whether this residue is a candidate in forming this helix. A > character indicates it starts a helix, a number indicates it is inside such a helix and a < character means it ends the helix.

   The next column contains a S character if this residue is a possible bend.

   Then there's a column indicating the chirality and this can either be positive or negative (i.e. the alpha torsion is either positive or negative).

   The last two column contain beta bridge labels. Lower case here means parallel bridge and thus upper case means anti parallel.

   ----------------------
   
   And then there is 1 extra column that I picked up!

   canonical helix interior is ' X S+   '
   
--]]


for i,arg in ipairs(toProcess) do
   if '' ~= arg then

      if 'quit' == arg then os.exit() end    -- so can insert 'quit' in input file of PDB IDs for testing and debugging

      local pdb,chain = arg:match('^(%d%w%w%w)(%w?)%s*$')
      if chain == '' then chain = nil end


      print('pdb= ' .. pdb .. ' chain= ' .. (chain and chain or ''))
      local prot = protein.get(pdb)
      --print(prot:tostring())
      prot:dbLoad(rfpg)
      print('loaded:')
      print(prot:tostring())

      print(prot:writeInternalCoords(nil,args['r'] and args['r'][1] or nil))   -- needs to be the i'th file being processed....
      
      --prot:internalToAtomCoords()
      --print(prot:writePDB(true,args['r'][1]))
      
         --[[
      local pdbid = parsers.parseProteinData(file, function (t) protein.load(t) end, chain)
      local prot = protein.get(pdbid)
      --print(pdbid,prot,prot:countDihedra())
      --print(prot:tostring())


      --- need dihedra and (optionally not) DSSP data
      -- have either PDB coords (pdb file), dihedra data only (pic file), or dihedra and DSSP data (rtm mkdssp output file)
      -- want to test for data consistency each way, essentially repeating buildProtein.lua test mode

      local dsspCount = prot:countDSSPs()  -- useful to know below
      local coordsInternal
      local coords3D
         -- ]]
      --end
   end
end

