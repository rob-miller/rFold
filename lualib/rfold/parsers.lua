
--- file parsers:
-- modified DSSP output parser, command line parser, PDB file parser
-- @module parsers

local parsers = {}

local function usage(valid_flags, valid_params)
   print(arg[0] .. ': command line error')
   if valid_flags then
      io.write(' valid flags: ')
      for f,v in pairs(valid_flags) do io.write( '-' .. f .. ' ') end
      print()
   end
   if valid_params then
      io.write(' valid parameters: ')
      for p,v in pairs(valid_params) do io.write( '-' .. p .. '=  ') end
      print()
   end
end

--- parse specified flags and parameters from command line, pass remaining args through
--
-- ./parseCmdLine.lua -p1=hello -f1 -f1 fee fie fo fum
-- @param valid_flags { ['f1'] = true, ['f2'] = true }
-- @param valid_params { ['p1'] = true, ['p2'] = true }
-- @return table with 'param' values as specified, 'flag' values counting number of times flag seen, and remaining args in sequential numbered slots 
function parsers.parseCmdLine (valid_flags, valid_params)
   local args = {}

   for i=1,#arg do
      if string.match(arg[i],'^-') then
         if string.match(arg[i],'=') then
            local key, val
            key, val = string.match(arg[i],'^-(%S+)=(%S+)$')
            if not (valid_params and valid_params[ key ]) then
               usage(valid_flags, valid_params)
               assert(nil, 'parameter -' .. key .. ' value ' .. val .. ' not recognised')
            end
            args[ key ] = val
         else
            flag = string.match(arg[i],'^-(%S+)$')
            if not (valid_flags and valid_flags[ flag ]) then
               usage(valid_flags, valid_params)
               assert(nil, 'flag -' .. flag .. ' not recognised')
            end
            args[flag] = (args[flag] or 0)+1
         end
      else
         args[#args +1] = arg[i]
      end
   end

   return args
end


--- parse rfold modified DSSP output file
-- @param infile DSSP output
-- @param callback (optional) call with table containing parsed fields for DSSP data lines, hedron length-angle-length lines, or dihedron angle lines
-- @return pdbid if callback specified, else sequential list of tables as for 'callback' above in order read
function parsers.prd (infile, callback)
   local rdt = {}
   rdt['triples']={}
   rdt['quads']={}
   rdt['dssp']={}
   local pdbid='xxxxx'
   if infile then io.input(infile) end
   for line in io.lines() do
      --print(line)
      if line:match('^HEADER ') then
         pdbid = line:match('(%S+)%s+%.$')
         local header = line:match('^(.+)'..pdbid..'.+$')
         if callback then
            callback({['pdbid'] = pdbid, ['header'] = header })
         else
            rdt['pdbid'] = pdbid
         end
         --print(pdbid)
      elseif line:ematch('^'..pdbid..'%s+(%a?)%s*(%w+):(%w+):(%w+)%s+(%S+)%s+(%S+)%s+(%S+)%s*$')  then
         local trip = {}
         trip['pdbid'],trip['chn'],trip[1],trip[2],trip[3],trip['len1'],trip['angle2'],trip['len3'] = pdbid, _1, _2, _3, _4, tonumber(_5), tonumber(_6), tonumber(_7)
         if callback then
            callback(trip)
         else
            rdt['triples'][#rdt['triples']+1] = trip
         end
         --print(pdbid,chn,a1,a2,a3,v1,v2,v3)
      elseif line:ematch('^'..pdbid..'%s+(%a?)%s*(%w+):(%w+):(%w+):(%w+)%s+(%S+)%s*$')  then
         local quad={}
         quad['pdbid'],quad['chn'],quad[1],quad[2],quad[3],quad[4],quad['dihedral1'] = pdbid, _1, _2, _3, _4, _5, tonumber(_6)
         if callback then
            callback(quad)
         else
            rdt['quads'][#rdt['quads']+1] = quad
         end
         --print(pdbid,chn,a1,a2,a3,a4,v1)
      elseif line:ematch('^%s+(%d+)%s+(%d+)%s(%a?)%s(%a)%s') then
         local dsp = {}
         local ndx
         dsp['ndx'], dsp['resn'], dsp['chn'], dsp['res'] = tonumber(_1),tonumber(_2),_3,_4
         if (dsp['res'] == dsp['res']:lower()) then    -- cys partner encoding by dssp
            dsp['res'] = 'C'
         end
         local ss = line:sub(17,17)
         local ss2 = line:sub(19,26)
         local nums = line:sub(106)
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
         if callback then
            callback(dsp)
         else
            rdt['dssp'][#rdt['dssp']+1] = dsp
         end
         --print(ndx, resn, chn, res, '.'..ss..'.','>'..ss2..'<',phi,psi,omg,zn,zc,bca)
      elseif line:match('^--$') then
      elseif line:match('^==== Secondary') then
      elseif line:match('^REFERENCE ') then
      elseif line:match('^COMPND') then
      elseif line:match('^SOURCE') then
      elseif line:match('^AUTHOR') then
      elseif line:match(' TOTAL NUMBER OF ') then
      elseif line:match('ACCESSIBLE SURFACE') then
      elseif line:match(' HISTOGRAMS ') then
      elseif line:match(' PER ') then
      elseif line:match('^%s+#') then
      else
         print('reading dssp -- failed to parse: ' .. line)
      end
      
   end
   io.close()
   if callback then
      return pdbid
   end
   return rdt
end

--  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO   KAPPA  ALPHA   PHI    PSI    OMG    X-N    Y-N    Z-N    X-CA   Y-CA   Z-CA   X-C    Y-C    Z-C    X-CB   Y-CB   Z-CB   B-CA
--   54   54 A V  H >< S+     0   0    0     -4,-1.5     3,-1.2    -3,-0.2     4,-0.3   0.878  111.0   52.5  -74.9  -41.5  179.2   31.0   27.8   17.9   31.5   27.1   16.7   31.8   28.0   15.5   30.5   25.9   16.4    8.9
--   84   84 A a  E     +BC  44  99A   0    -40,-2.9   -40,-2.1    -2,-0.4     2,-0.4  -0.990   10.7  178.2 -115.7  124.1  178.6   29.2    9.4   15.9   28.3    8.5   15.5   28.9    7.6   14.4   27.0    9.1   15.0    9.6

--1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
--         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2      
return parsers
