--[[
   utils.lua
   
Copyright 2016, 2017 Robert T. Miller

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

--- utility functions:
-- modified extended string:match(); atom token,hedron and dihedron key manipulators
-- not a class but works better with ldoc if we say it is
-- @classmod utils

local utils = {}
local chemdata = require 'rfold.chemdata'

------- global flags
utils.warn=nil

--- extend string.match to test and capture patterns - goes into string: namespace
--
-- <br>`https://inspired-lua.org/index.php/2013/05/extend-string-match-to-test-and-capture-patterns/`
-- @param pat pattern same as for string.match()
-- @return same as string.match() but sets global variables _1, _2, ... for captures on success
function string:ematch(pat)
   local matches = {string.match(self, pat)}    -- call the original match to do the work
   for i = 1, #matches do                 -- #matches == 0 if no matches
      _G["_" .. i] = matches[i]          -- assign captures to global variables
   end
   return unpack(matches)                 -- return original results
end

--- split_newlines
-- http://lua-users.org/wiki/EnhancedFileLines
--
-- Splits string s into array of lines, returning the result.
-- New-line character sequences ("\n", "\r\n", "\r"),
-- if any, are NOT included at the ends of the lines.
-- @param s string containing newline characters
-- @return array of lines 
function utils.split_newlines(s)
  local ts = {}
  local posa = 1
  while 1 do
    local pos, chars = s:match('()([\r\n].?)', posa)
    if pos then
       if 1 < pos then           -- only if not newline only
          local npos = pos -1    -- drop newline
          local line = s:sub(posa, npos)
          ts[#ts+1] = line
       end
       if chars == '\r\n' then pos = pos + 1 end  -- skip past cmplete newline cahr
       posa = pos + 1
    else
      local line = s:sub(posa)
      if line ~= '' then ts[#ts+1] = line end
      break      
    end
  end
  return ts
end

--- execute line-by-line compare of 2 passed strings
-- @param s1 string containing newline characters
-- @param s2 string containing newline characters
-- @return false if not perfect match, else number of lines compared
function utils.lineByLineCompare(s1,s2)
   local s1t, s2t = {}, {}
   for line in s1:gmatch("[^\r\n]+") do s1t[#s1t+1] = line end
   for line in s2:gmatch("[^\r\n]+") do s2t[#s2t+1] = line end

   local ecnt = 2
   for i,lin in ipairs(s1t) do
      if lin ~= s2t[i] then
         print()
         print(lin)
         print(s2t[i])
         ecnt = ecnt-1
         if 0>ecnt then return false end
      --else
         --print(lin,s2t[i])         
      end
   end
   if #s1t ~= #s2t then
      print(' different linecounts.')
      return false
   end
   --return true
   return #s1t
end


--- generate a string key for hedon or dihedron consisting of atom tokens separated by ':'s
--
-- @param atoms ... tokens of the form (sequence postion)(residue)(atom string) e.g. 224PCB = C-beta of Proline at sequence position 224
-- @return string key e.g. 1MN:1MCA:1MC:2QN
function utils.genKey(...)
   local key = ''
   for i, s in ipairs{...} do
      if '' ~= key then key = key .. ':' end
      key = key .. s
   end
   return key
end

--- split 3- (hedron) or 4- (dihedron) element string key into atom token constituents
--
-- @param k string key to split e.g. 2ECA:2EC:3TN:3TCA
-- @return table of sequential fields e.g. { '2ECA', '2EC', '3TN', '3TCA' }
function utils.splitKey(k)
   if k:ematch('^(-?%d+_?%w+):(-?%d+_?%w+):(-?%d+_?%w+):(-?%d+_?%w+)$') then return { _1, _2, _3, _4 }
   elseif k:ematch('^(-?%d+_?%w+):(-?%d+_?%w+):(-?%d+_?%w+)$') then return { _1, _2, _3 }
   else assert(nil,'utils.splitKey fail on '..k) end
end

--- split atom token of form (sequence postion)(residue)(atom string) into constituents
-- @param k atom token key to split
-- @return table of constituents in order [1] sequence postion [2] residue [3] atom
function utils.splitAtomKey(k)
   if k:ematch('^(-?%d+)(%a)(%w+)$') then return { tonumber(_1), _2, _3 }
   elseif k:ematch('^(-?%d+)(_)(%w+)$') then return { tonumber(_1), _2, _3 }
   else assert(nil,'utils.splitAtomKey fail on '..k) end
end

--- convert table of atom-name-only (from chemdata.sidechains) to atomKey
-- @param base <seqpos><residue>
-- @param t table of atom names in ipairs from chemdata.sidechains e.g. { 'CA', 'CB', 'CG1' }
-- @return new table with atomKeys
function utils.addBase(base,t)
   local nt = {}
   for i,v in ipairs(t) do
      nt[i] = base .. t[i]
   end
   return nt
end

--- round number to value presented for %9.5lf in C
-- @param num the initial value
-- @return the value rounded to 5 decimal places
function utils.setAccuracy95(num)
   local r = tonumber( string.format("%9.5f",num) )
   if -0.00000 == r then r=0.00000 end
   return r
end

--- round number to value presented for %8.3lf in C
-- @param num the initial value
-- @return the value rounded to 3 decimal places
function utils.setAccuracy83(num)
   local r = tonumber( string.format("%8.3f",num) )
   if -0.000 == r then r=0.000 end
   return r 
end

--- sort table by keys, from pil 1st ed 19.3
-- @param t table of key, value pairs
-- @param f optional comparison function
-- @return iterator for use as 'for name, line in pairsByKeys(lines) do'
function utils.pairsByKeys (t, f)
   local a = {}
   for n in pairs(t) do table.insert(a, n) end
   if f then table.sort(a, f) end
   local i = 0      -- iterator variable
   local iter = function ()   -- iterator function
      i = i + 1
      if a[i] == nil then return nil
      else return a[i], t[a[i]]
      end
   end
   return iter
end

--- generate a PDB ATOM record
-- @param ndx index of record
-- @param name atom NAME (CA)
-- @param res 1-letter residue code
-- @param chain chain ID
-- @param resn residue number 
-- @param ax atom X coordinate
-- @param ay atom Y coordinate
-- @param az atom Z coordinate
-- @param occ occupancy
-- @param tempFact temperature factor (B factor)
-- @return string with text and newline for one ATOM record
function utils.atomString(ndx,name,res,chain,resn,ax,ay,az,occ,tempFact)
   return string.format('%6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n',
                        'ATOM  ', ndx, name,' ', chemdata.res3[res], chain, resn,' ', ax, ay, az, occ, tempFact, '    ', string.sub(name,1,1),'  ')
end

--- open a named file, write the passed string, close the file
-- @param fname name of file to open
-- @param s string to write
function utils.writeFile(fname,s)
   print('open ' .. fname)
   local file,err = io.open(fname,'w')
   assert(file,err)
   file:write(s)
   file:close()
end

--- get our hostname
-- @return string with hostname
function utils.getHostname()
    local f = io.popen ("/bin/hostname")
    local hostname = f:read("*a") or ""
    f:close()
    hostname =string.gsub(hostname, "\n$", "")
    return hostname
end

--- sort a pair of variables
-- @param a1 first thing
-- @param a2 second thing
-- @return table with [1] and [2] set to a1/a2 according to lua table.sort()
function utils.orderPair(a1, a2)
   local t = { a1, a2 }
   table.sort(t)
   return t
end
   
return utils

