--[[
   parsers.lua
   
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

--- file parsers:
-- command line parser, PDB file / (modified) DSSP output parser
--
-- not a class but works better with ldoc if we say it is
-- @classmod parsers

local utils = require 'rfold.utils'
local chemdata = require 'rfold.chemdata'

local parsers = {}

local function usage(argFormat, info, valid_flags, valid_params)
   print( arg[0] .. ':\n')
   if argFormat then
      if (valid_flags or valid_params) then
         print(' usage: ' .. arg[0] .. '  [options ... ] ' .. argFormat )
      else
         print(' usage: ' .. arg[0] .. argFormat )
      end
   else
      if (valid_flags or valid_params) then
         print(' usage: ' .. arg[0] .. ' [ options ... ] ' )
      else
         print(' usage: ' .. arg[0] )
      end      
   end
   
   if (info) then 
      print(info)
   end
   
   if valid_flags then
      print(' valid flags: ')
      for f,v in pairs(valid_flags) do io.write( '   -' .. f .. ' ' .. v ..'\n') end
      print()
   end
   if valid_params then
      print(' valid parameters: ')
      for p,v in pairs(valid_params) do io.write( '     -' .. p .. '=' .. v .. '\n') end
      print()
   end
end

--- parse specified flags and parameters from command line, pass remaining args through
--
-- ./parseCmdLine.lua -p1=hello -f1 -f1 fee fie fo fum
-- @param argFormat format of command line argments for help e.g. '[ <pdbid[chain id]> | <pdb file> | <pic file> ] ...'
-- @param info descriptive text for what the program does
-- @param valid_flags { ['f1'] = "flag1 effect message", ['f2'] = "flag2 message" }
-- @param valid_params { ['p1'] = "parameter1 explanation ...", ['p2'] = "parameter2 explanation ..." }
-- @return table with 'param' values as specified, 'flag' values counting number of times flag seen, and remaining args in sequential numbered slots 
function parsers.parseCmdLine (argFormat, info, valid_flags, valid_params)
   local args = {}

   if (argFormat and (0 == #arg)) then
      usage(argFormat, info, valid_flags, valid_params)
      os.exit()
   end
   for i=1,#arg do
      if string.match(arg[i],'^-') then
         if string.match(arg[i],'=') then
            local key, val
            key, val = string.match(arg[i],'^-(%S+)=(%S+)$')
            if not val then
               key, val = string.match(arg[i],'^-(%S+)="(.+)"$')
            end
            if not (valid_params and valid_params[ key ]) then
               usage(argFormat, info, valid_flags, valid_params)
               print('parameter -' .. (key and key or 'nil') .. ' value ' .. (val and val or 'nil') .. ' not recognised in fragment: ' .. arg[i])
               os.exit()
            end
            if not args[key] then args[key]={} end
            args[ key ][#args[key]+1] = val
         else
            flag = string.match(arg[i],'^-(%S+)$')
            if not (valid_flags and valid_flags[ flag ]) then
               usage(argFormat, info, valid_flags, valid_params)
               print('flag -' .. flag .. ' not recognised')
               os.exit()
            end
            args[flag] = (args[flag] or 0)+1
         end
      else
         args[#args +1] = arg[i]
      end
   end

   return args
end


--- parse PDB file, internal coordinates (.pic) file, or rfold modified DSSP output file
-- @param infile (file or string) rfold modified DSSP output or PDB file or lines with hedron or dihedron specifications
-- @param callback (optional) call with table containing parsed fields for DSSP data lines, PDB ATOM etc. lines, hedron length-angle-length lines, or dihedron angle lines
-- @param chain only process line if specified chain matches input (where relevant)
-- @param infname (optional) specify source for infile=string data
-- @return pdbid if callback specified, else sequential list of tables as for 'callback' above in order read
function parsers.parseProteinData (infile, callback, chain, infname)
   local rdt = {}
   rdt['triples']={}
   rdt['quads']={}
   rdt['dssp']={}
   local pdbid='xxxxx'
   local fh
   assert(infile,'parseProteinData: no input file')

   local err
   local stringArray
   local fname=nil
   if infile:match('%.gz$') then
      fh,err = io.popen('zcat ' .. infile)
      fname = infile:match('(%S+).gz$')
   elseif infile:match('.[\r\n].') then
      stringArray = utils.split_newlines(infile)
      --print('stringArray count',#stringArray)
      fname = infname or 'text_string'
   else
      fh,err = io.open(infile)
      fname = infile
   end
   if err then
      print(infile,err)
      return nil
   end
   if fh then
      stringArray = {}
      for line in fh:lines() do
         stringArray[#stringArray+1] = line  -- hard on memory usage but need flexibility to call with string
      end
   end
   for key,line in ipairs(stringArray) do
      --print(line)
      if line:match('^HEADER ') then  -- dssp or pdb or pic
         pdbid = line:match('(%S+)%s*%.?$')
         --print(pdbid)
         local header = line:gsub('%.$','')       -- wipe the '.' dssp puts at end
         header = header:gsub('%s*$','')           -- wipe all trailing space
         header = string.format('%-80s\n', header)   -- reformat to space-filled 80 chars long
         
         --local header = line:match('^(.+' .. pdbid .. ').+$') .. '              ' 
         --print(pdbid,header)
         if callback then
            callback({['pdbid'] = pdbid, ['header'] = header })
         else
            rdt['pdbid'] = pdbid
            rdt['header'] = header
         end
         --print(pdbid)
      elseif line:match('^TITLE') then  -- pdb only -- dssp does not pass
         local title = line
         title = title:gsub("'",'-PRIME')
         title = title:gsub('%*$','')           -- wipe all trailing space
         title = string.format('%-80s\n', title)   -- reformat to space-filled 80 chars long
         if callback then
            callback({['pdbid'] = pdbid, ['title'] = title })
         else
            if not rdt['title'] then rdt['title'] = title end
         end
         --print(pdbid)
      elseif line:match('^COMPND') then
         local compnd = line
         if callback then
            callback({['pdbid'] = pdbid, ['compnd'] = compnd })
         else
            rdt['compnd'] = compnd
         end
         --print(pdbid)
      elseif line:ematch('^'..pdbid..'%s+(%a?)%s*(-?[%w_]+):(-?[%w_]+):(-?[%w_]+)%s+(%S+)%s+(%S+)%s+(%S+)%s*$')  then  -- hedron spec   (.pic)
         local trip = {}
         trip['pdbid'],trip['chn'],trip[1],trip[2],trip[3],trip['len1'],trip['angle2'],trip['len3'] = pdbid, _1, _2, _3, _4, tonumber(_5), tonumber(_6), tonumber(_7)
         if not chain or chain == trip['chn'] then
            if callback then
               callback(trip)
            else
               rdt['triples'][#rdt['triples']+1] = trip
            end
         end
         --print(pdbid,chn,a1,a2,a3,v1,v2,v3)
      elseif line:ematch('^'..pdbid..'%s+(%a?)%s*(-?[%w_]+):(-?[%w_]+):(-?[%w_]+):(-?[%w_]+)%s+(%S+)%s*$')  then   -- dihedron spec     (.pic)
         local quad={}
         quad['pdbid'],quad['chn'],quad[1],quad[2],quad[3],quad[4],quad['dihedral1'] = pdbid, _1, _2, _3, _4, _5, tonumber(_6)
         if not chain or chain == quad['chn'] then
            if callback then
               callback(quad)
            else
               rdt['quads'][#rdt['quads']+1] = quad
            end
         end
         --print(pdbid,chn,a1,a2,a3,a4,v1)
      elseif line:ematch('^%s+(%d+)%s+(-?%d+)%s(%a?)%s(%a)%s') then    -- modified dssp output
         --print(_1,_2,_3,_4)
         local dsp = {}
         local ndx
         dsp['ndx'], dsp['resn'], dsp['chn'], dsp['res'] = tonumber(_1),tonumber(_2),_3,_4
         if (dsp['res'] == dsp['res']:lower()) then    -- cys partner encoding by dssp
            dsp['res'] = 'C'
         end
         --print('dssp ',dsp['ndx'], dsp['resn'], dsp['chn'], dsp['res'])
         dsp['ss'] = line:sub(17,17)
         dsp['ss2'] = line:sub(19,26)
         dsp['acc'] = tonumber(line:sub(36,38))
         local nums = line:sub(106)
         --print('nums', nums)
         assert(nums:ematch('^%s*(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s*$'), 'failed to parse dssp: '..line)

         dsp['n'] = {}
         dsp['ca'] = {}
         dsp['c'] = {}
         dsp['cb'] = {}
         dsp['o'] = {}

         dsp['pdbid'] = pdbid
         dsp['phi'], dsp['psi'], dsp['omg'], dsp['n']['x'], dsp['n']['y'], dsp['n']['z'] = tonumber(_1),tonumber(_2),tonumber(_3),tonumber(_4),tonumber(_5),tonumber(_6)
         dsp['ca']['x'], dsp['ca']['y'], dsp['ca']['z'], dsp['c']['x'], dsp['c']['y'], dsp['c']['z'] = tonumber(_7),tonumber(_8),tonumber(_9),tonumber(_10),tonumber(_11),tonumber(_12)
         dsp['cb']['x'], dsp['cb']['y'], dsp['cb']['z'],dsp['o']['x'], dsp['o']['y'], dsp['o']['z'] = tonumber(_13),tonumber(_14),tonumber(_15),tonumber(_16),tonumber(_17),tonumber(_18)
         dsp['bca'] = tonumber(_19)
         if not chain or chain == dsp['chn'] then
            if callback then
               callback(dsp)
            else
               rdt['dssp'][#rdt['dssp']+1] = dsp
            end
         end
         --print('dssp:', dsp['ndx'], dsp['resn'], dsp['chn'], dsp['res'], '.'.. (dsp['ss'] and dsp['ss'] or 'no_ss') ..'.','>'.. (dsp['ss2'] and dsp['ss2'] or 'no_ss2') ..'<',dsp['phi'],dsp['psi'],dsp['omg'],dsp['n']['z'],dsp['c']['z'],dsp['bca'])
--[[
         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N
ATOM   2606  CD  GLN A 324    -147.432-100.309  -1.110  1.00  0.00           C
--]]
         --                ATOM   2606  CD  GLN A 324    -147.432-100.309  -1.110  1.00  0.00           C
      elseif line:ematch('^ATOM%s+(%d+)%s*%w+') then
         --print(line)
         local pdb = {}
         pdb['ndx'] = tonumber(_1)
         pdb['atom'] = line:sub(13,16):gsub('%s','')
         pdb['altloc'] = line:sub(17,17);
         if ' '==pdb['altloc'] then pdb['altloc']=nil end
         
         pdb['res3'] = line:sub(18,20)
         pdb['chn'] = line:sub(22,22)
         pdb['resn'] = tonumber(line:sub(23,26))
         pdb['icode'] = line:sub(27,27)
         if ' '==pdb['icode'] then pdb['icode']=nil end
         
         pdb['x'],pdb['y'],pdb['z'] = tonumber(line:sub(31,38)),tonumber(line:sub(39,46)),tonumber(line:sub(47,54))
         pdb['occ'],pdb['tempfact'] = tonumber(line:sub(55,60)),tonumber(line:sub(61,66))
         pdb['elem'] = line:sub(77,78):gsub('%s','')
         
         pdb['pdbid'] = pdbid

         --for k,v in pairs(pdb) do print(k,v) end
         --print()

         if (not chain) or (chain == pdb['chn']) then
            if callback then
               callback(pdb)
            else
               rdt['pdb'][#rdt['pdb']+1] = pdb
            end
         end
         -- lines to ignore:
         -- dssp / pdb
      elseif line:match('^--$') then
      elseif line:match('^==== Secondary') then
      elseif line:match('^REFERENCE ') then
      elseif line:match('^SOURCE') then
      elseif line:match('^AUTHOR') then
      elseif line:match(' TOTAL NUMBER OF ') then
      elseif line:match('ACCESSIBLE SURFACE') then
      elseif line:match(' HISTOGRAMS ') then
      elseif line:match(' PER ') then
      elseif line:match('^%s+#') then
         -- dssp
      elseif line:match('^%s+%d+%s+%!%s+') then
         -- dssp chain break
      elseif line:match('^%s+%d+%s+!%*%s+') then
         -- dssp chain termination
      elseif line:match('^SEQADV') then
      elseif line:match('^LINK') then
      elseif line:match('^KEYWDS') then
      elseif line:match('^EXPDTA') then
      elseif line:match('^REVDAT') then
      elseif line:match('^JRNL') then
      elseif line:match('^REMARK') then
      elseif line:match('^DBREF') then
      elseif line:match('^SEQRES') then
         local sr = {}
         sr['chn'] = line:sub(12,12)
         sr['pdbid'] = pdbid
         local s = ''
         local psr = line:sub(20)
         for res3 in psr:gmatch("%a+") do
            local r1 = chemdata.res1[res3]
            if not r1 then
               -- io.stderr:write("SEQRES code '" .. res3 .."' not recognised, file: " .. infile .. "\n")
               r1 = 'x'
               -- os.exit()
            end
            s = s .. r1
         end
         sr['seqres'] = s
         if (not chain) or (chain == sr['chn']) then
            if callback then
               callback(sr)
            else
               if not rdt['seqres'] then
                  rdt['seqres'] = { sr['chn'] }
               elseif not rdt['seqres'][ sr['chn'] ] then
                  rdt['seqres'][ sr['chn'] ] = s
               else
                  rdt['seqres'][ sr['chn'] ] = rdt['seqres'][ sr['chn'] ] .. s
               end
            end
         end
      elseif line:match('^HET ') then
      elseif line:match('^HETNAM') then
      elseif line:match('^HETSYN') then
      elseif line:match('^FORMUL') then
      elseif line:match('^HELIX') then
      elseif line:match('^SHEET') then
      elseif line:match('^SSBOND') then
      elseif line:match('^CISPEP') then
      elseif line:match('^SITE') then
      elseif line:match('^CRYST') then
      elseif line:match('^ORIGX') then
      elseif line:match('^MTRIX') then
      elseif line:match('^SCALE') then
      elseif line:match('^ANISOU') then
      elseif line:ematch('^TER%s+(%d+)') then
         --print(line)
         local ter = { ['TER'] = true }
         ter['ndx'] = tonumber(_1)
         ter['res3'] = line:sub(18,20)
         ter['chn'] = line:sub(22,22)
         ter['resn'] = tonumber(line:sub(23,26))
         ter['pdbid'] = pdbid
         
         if (not chain) or (chain == ter['chn']) then
            if callback then
               callback(ter)
            else
               rdt['pdb'][#rdt['pdb']+1] = ter
            end
         end
         
      elseif line:match('^HETATM') then    -- todo: add support
      elseif line:match('^MODRES') then
      elseif line:match('^MODEL') then
      elseif line:match('^NUMMDL') then         
      elseif line:match('^CAVEAT') then
      elseif line:match('^SPRSDE') then
      elseif line:match('^CONECT') then
      elseif line:match('^MASTER') then
      elseif line:match('^END') then
         
      else
         io.stderr:write('failed to parse: >' .. line .. '<\n')
      end
      
   end
   io.close()

   if 'xxxxx'==pdbid then
      print('read problem on ' .. fname)
      return nil
   end
   
   if callback then
      callback({['pdbid'] = pdbid, ['filename'] = fname})
      return pdbid
   end
   rdt['filename'] = fname
   return rdt
end

--  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO   KAPPA  ALPHA   PHI    PSI    OMG    X-N    Y-N    Z-N    X-CA   Y-CA   Z-CA   X-C    Y-C    Z-C    X-CB   Y-CB   Z-CB   B-CA
--   54   54 A V  H >< S+     0   0    0     -4,-1.5     3,-1.2    -3,-0.2     4,-0.3   0.878  111.0   52.5  -74.9  -41.5  179.2   31.0   27.8   17.9   31.5   27.1   16.7   31.8   28.0   15.5   30.5   25.9   16.4    8.9
--   84   84 A a  E     +BC  44  99A   0    -40,-2.9   -40,-2.1    -2,-0.4     2,-0.4  -0.990   10.7  178.2 -115.7  124.1  178.6   29.2    9.4   15.9   28.3    8.5   15.5   28.9    7.6   14.4   27.0    9.1   15.0    9.6

--1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
--         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2      
return parsers
