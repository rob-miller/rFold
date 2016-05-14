
--- utility functions:
-- modified extended string:match(); atom token,hedron and dihedron key manipulators
-- @module utils

local utils = {}

--- residue single letter to 3-letter name conversion
-- @table res3
utils.res3 = { G = 'GLY', A = 'ALA', V = 'VAL', L = 'LEU', I = 'ILE', M = 'MET', F = 'PHE', P = 'PRO', S = 'SER', T = 'THR',
               C = 'CYS', N = 'ASN', Q = 'GLN', Y = 'TYR', W = 'TRP', D = 'ASP', E = 'GLU', H = 'HIS', K = 'LYS', R = 'ARG' }

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
-- @param k string key to split
-- @return table of sequential fields
function utils.splitKey(k)
   if k:ematch('^(%w+):(%w+):(%w+):(%w+)$') then return { _1, _2, _3, _4 }
   elseif k:ematch('^(%w+):(%w+):(%w+)$') then return { _1, _2, _3 }
   else assert(nil,'utils.splitKey fail on '..k) end
end

--- split atom token of form (sequence postion)(residue)(atom string) into constituents
-- @param k atom token key to split
-- @return table of constituents in order [1] sequence postion [2] residue [3] atom
function utils.splitAtomKey(k)
   if k:ematch('^(%d+)(%a)(%w+)$') then return { tonumber(_1), _2, _3 }
   else assert(nil,'utils.splitAtomKey fail on '..k) end
end


return utils

