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

local rfpg = require 'rfold.rfPostgresql'  -- autocommit false by default
rfpg.dbReconnect(1)


local args = parsers.parseCmdLine(
   '[ <pdbid[chain id]> ] ...',
   '\n retrieve (optionally convert) coordinates for specified pdbid from rFold database.\n',
   {
      ['l'] = 'list:  list available PDBids in database'
   },
   {
      ['w'] = 'pic|gpdb|<filename(.pic|.gpdb)> : write specified format to pdbid.<format> or <filename>',
      ['f'] = '<input file> : process IDs listed in <input file> (1 per line) followed by any on command line',
      ['r'] = '<firstResidue:lastResidue> : only output residues within the range specified '
   }
)


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

      print(prot:writeInternalCoords(nil,args['r'][1]))   -- needs to be the i'th file being processed....
      
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

