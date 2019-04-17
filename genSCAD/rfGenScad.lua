#!/usr/bin/env luajit

--[[
   rfGenSCAD.lua
   
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

--- rfGenSCAD.lua [ options ] [<PDBid[chainId]>] ...
-- 
-- Generate OpenSCAD output from rFold data
--
-- options:
-- 
-- -l                                    : list pdb ids, chain ids loaded in database
-- -f=<list of PDB codes>                : read PDB codes and optional chain IDs from specified file
-- -r=<start>:<finish>                   : only get specified range of residue positions from target protein
-- 

--- work in progress


local parsers = require 'rfold.parsers'
local protein = require "rfold.protein"
local utils = require 'rfold.utils'
local chain = require 'rfold.chain'
local chemdata = require 'rfold.chemdata'

--local rfpg = require 'rfold.rfPostgresql'  -- autocommit false by default
--rfpg.dbReconnect(1)

local scaleDefault=2.0 --10.0
--datafileDefault='rFold.scad'

local args = parsers.parseCmdLine(
   '[ <pdbid[chain id]>  | <pdb file> | <pic file> ] ...',
   '\n retrieve (optionally convert) coordinates for specified pdbid from rFold database, output as scad file for 3D printing\n',
   {
      ['l'] = 'list:  list available PDBids in database',
      ['b'] = 'backbone only: do not generate sidechains',
      ['wpa'] = 'write peptideOffsets.scad file'
   },
   {
      ['s'] = '<scale> : slicer units per angstrom (usually millimetres). default=' .. scaleDefault,
      --['d'] = '<datafile.scad> : coordinate data filename (only written if input structure data).  default=' .. datafileDefault,

      --['w'] = '<filename.scad> : write commands to <filename.scad>  default=stdout.',
      ['r'] = '<firstResidue:lastResidue> : only output residues within the range specified ',
      ['f'] = '<input file> : process IDs listed in <input file> (1 per line) followed by any on command line'
   }
)

local scale = args['s'] and args['s'][1] or scaleDefault
local datafile = args['w'] and args['w'][1] or datafileDefault

if args['wpa'] then
   utils.writeFile('peptideOffsets.scad', protein.scadOffsets())
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


for i,arg in ipairs(toProcess) do
   if '' ~= arg then
      local file, chain, pdbid
      
      if 'quit' == arg then os.exit() end    -- so can insert 'quit' in input file of PDB IDs for testing and debugging

      pdbid,chain = arg:match('^(%d%w%w%w)(%w?)%s*$')   -- try pdb id for db query first

      if not pdbid then
         file,chain = arg:match('^(%S+)%s?(%w?)%s*$')   -- needs space at end if doing on command line, e.g. '7RSAA '
      end
      
      if chain == '' then chain = nil end
      
      if not pdbid then
         -- load current file
         pdbid = parsers.parseProteinData(file, function (t) protein.load(t) end, chain)
      end
      
      if not pdbid then goto continue end   -- give up

      local prot = protein.get(pdbid)

      if not file then 
         prot:dbLoad(rfpg)
      end

      if prot:countDihedra() > 0 then     -- need atom coords to scale
         prot:internalToAtomCoords()
         prot:clearInternalCoords()
      end

      prot:scaleAtomCoords(scale)

      prot:atomsToInternalCoords()
      prot:clearInitNCaC()
      prot:clearAtomCoords()
      prot:internalToAtomCoords()
      
      --print('loaded:')
      --print(prot:tostring())

      local scadData =  prot:writeSCAD(scale,(args['r'] and args['r'][1] or nil), (args['b'] and args['b'] or nil))

      utils.writeFile(pdbid .. 'coords.scad', 'protein = ' .. scadData .. ';\n')

      
      --print('chain(chains,residues,hedra);')

      --print('translate([10,10,' .. (chemdata.covalent_radii['Nres']/2) * scale .. '])')
      --print('  rotate([0,90,0])')
      --print('    rotate([0,0,residues[0][2][0]])')
      --print('      amideChain(amides,chains,residues,hedra);')

      
      -- print(prot:writeInternalCoords(nil,args['r'][1]))   -- 'r' needs to be the i'th file being processed....
      
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

      ::continue::

   end
end



--[[

   OpenSCAD data format:
   
   chain:
   [ res_id = <chain><position>,
     res_code (G,A,V,..),
     render_option,
     [ residue transform matrix in world ],
     res_ndx
   ] ...

   residue:
   [
     [ dangle,
       h1_ndx,
       h2_ndx,
       reverse,
       [ dihedral transform matrix in residue ],
     ] ...
   ] ...
   
   hedra:
   [ len1,
     angle2,
     len3,
     atom1 (N,C,S,...),
     atom2,
     atom3,
     render_option
   ] ...

--]]
