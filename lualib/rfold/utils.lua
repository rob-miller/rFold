
--- utility functions:
-- modified extended string:match(); atom token,hedron and dihedron key manipulators
-- @module utils

local utils = {}

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

--- convert table of atom-name-only (from utils.sidechains) to atomKey
-- @param base <seqpos><residue>
-- @param t table of atom names in ipairs from utils.sidechains e.g. { 'CA', 'CB', 'CG1' }
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
   return tonumber( string.format("%9.5f",num) )
end

--- round number to value presented for %8.3lf in C
-- @param num the initial value
-- @return the value rounded to 3 decimal places
function utils.setAccuracy83(num)
   return tonumber( string.format("%8.3f",num) )
end

--- sort table by keys, from pil 1st ed 19.3
-- @param t table of key, value pairs
-- @param f optional comparison function
-- @return iterator for use as 'for name, line in pairsByKeys(lines) do'
function utils.pairsByKeys (t, f)
   local a = {}
   for n in pairs(t) do table.insert(a, n) end
   table.sort(a, f)
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
                        'ATOM  ', ndx, name,' ', utils.res3[res], chain, resn,' ', ax, ay, az, occ, tempFact, '    ', string.sub(name,1,1),'  ')
end

--- residue single letter to 3-letter name conversion
-- @table res3
utils.res3 = { G = 'GLY', A = 'ALA', V = 'VAL', L = 'LEU', I = 'ILE', M = 'MET', F = 'PHE', P = 'PRO', S = 'SER', T = 'THR',
               C = 'CYS', N = 'ASN', Q = 'GLN', Y = 'TYR', W = 'TRP', D = 'ASP', E = 'GLU', H = 'HIS', K = 'LYS', R = 'ARG' }

--- residue 3-letter name to single letter conversion
-- @table res1
utils.res1 = { ['GLY'] = 'G', ['ALA'] = 'A', ['VAL'] = 'V', ['LEU'] = 'L', ['ILE'] = 'I', ['MET'] = 'M', ['PHE'] = 'F', ['PRO'] = 'P', ['SER'] = 'S', ['THR'] = 'T',
               ['CYS'] = 'C', ['ASN'] = 'N', ['GLN'] = 'Q', ['TYR'] = 'Y', ['TRP'] = 'W', ['ASP'] = 'D', ['GLU'] = 'E', ['HIS'] = 'H', ['LYS'] = 'K', ['ARG'] = 'R' }


--- per residue sidechain hedra and dihedra definitions, in order of output for internal coordinates spcification file
-- @table sidechains
utils.sidechains = {
   ['V'] = {
      { 'CA', 'CB', 'CG1' },
      { 'N', 'CA', 'CB', 'CG1' }, -- chi1
      { 'CA', 'CB', 'CG2' },
      { 'N', 'CA', 'CB', 'CG2' },
      ['chi1'] = 2
   },
   ['L'] = {
      { 'CA', 'CB', 'CG' },
      { 'N', 'CA', 'CB', 'CG' }, -- chi1
      { 'CB', 'CG', 'CD1' },
      { 'CA', 'CB', 'CG', 'CD1' }, -- chi2
      { 'CB', 'CG', 'CD2' },
      { 'CA', 'CB', 'CG', 'CD2' },
      ['chi1'] = 2, ['chi2'] = 4
   },
   ['I'] = {
      { 'CA', 'CB', 'CG1' },
      { 'N', 'CA', 'CB', 'CG1' },   -- chi1
      { 'CB', 'CG1', 'CD1' },
      { 'CA', 'CB', 'CG1', 'CD1' },   -- chi2
      { 'CA', 'CB', 'CG2' },
      { 'N', 'CA', 'CB', 'CG2' },
      ['chi1'] = 2, ['chi2'] = 4
   },
   ['M'] = {
      { 'CA', 'CB', 'CG' },
      { 'N', 'CA', 'CB', 'CG' },   -- chi1
      { 'CB', 'CG', 'SD' },
      { 'CA', 'CB', 'CG', 'SD' },   -- chi2
      { 'CG', 'SD', 'CE' },
      { 'CB', 'CG', 'SD', 'CE' },   -- chi3
      ['chi1'] = 2, ['chi2'] = 4, ['chi3'] = 6
   },
   ['F'] = {
      { 'CA', 'CB', 'CG' },
      { 'N', 'CA', 'CB', 'CG' },   -- chi1
      { 'CB', 'CG', 'CD1' },
      { 'CA', 'CB', 'CG', 'CD1' },   -- chi2
      { 'CG', 'CD1', 'CE1' },
      { 'CB', 'CG', 'CD1', 'CE1' },
      { 'CD1', 'CE1', 'CZ' },
      { 'CG', 'CD1', 'CE1', 'CZ' },
      { 'CB', 'CG', 'CD2' },
      { 'CA', 'CB', 'CG', 'CD2' },
      { 'CG', 'CD2', 'CE2' },
      { 'CB', 'CG', 'CD2', 'CE2' },
      ['chi1'] = 2, ['chi2'] = 4
   },
   ['P'] = {
      { 'CA', 'CB', 'CG' },
      { 'N', 'CA', 'CB', 'CG' },   -- chi1
      { 'CB', 'CG', 'CD' },
      { 'CA', 'CB', 'CG', 'CD' },   -- chi2
      ['chi1'] = 2, ['chi2'] = 4
   },
   ['S'] = {
      { 'CA', 'CB','OG' },
      { 'N', 'CA', 'CB','OG' },   -- chi1
      ['chi1'] = 2
   },
   ['T'] = {
      { 'CA', 'CB', 'OG1' },
      { 'N', 'CA', 'CB', 'OG1' },   -- chi1
      { 'CB', 'OG1', 'CG2' },
      { 'CA', 'CB', 'OG1', 'CG2' },
      ['chi1'] = 2
   },
   ['C'] = {
      { 'CA', 'CB', 'SG' },
      { 'N', 'CA', 'CB', 'SG' },   -- chi1
      ['chi1'] = 2
   },
   ['N'] = {
      { 'CA', 'CB', 'CG' },
      { 'N', 'CA', 'CB', 'CG' },   -- chi1
      { 'CB', 'CG', 'OD1' },
      { 'CA', 'CB', 'CG', 'OD1' },   -- chi2
      { 'CB', 'CG', 'ND2' },
      { 'CA', 'CB', 'CG', 'ND2' },
      ['chi1'] = 2, ['chi2'] = 4
   },
   ['Q'] = {
      { 'CA', 'CB', 'CG' },
      { 'N', 'CA', 'CB', 'CG' },   -- chi1
      { 'CB', 'CG', 'CD' },
      { 'CA', 'CB', 'CG', 'CD' },   -- chi2
      { 'CG', 'CD', 'OE1' },
      { 'CB', 'CG', 'CD', 'OE1' },   -- chi3
      { 'CG', 'CD', 'NE2' },
      { 'CB', 'CG', 'CD', 'NE2' },
      ['chi1'] = 2, ['chi2'] = 4, ['chi3'] = 6
   },
   ['Y'] = {
      { 'CA', 'CB', 'CG' },
      { 'N', 'CA', 'CB', 'CG' },   -- chi1
      { 'CB', 'CG', 'CD1' },
      { 'CA', 'CB', 'CG', 'CD1' },   -- chi2
      { 'CG', 'CD1', 'CE1' },
      { 'CB', 'CG', 'CD1', 'CE1' },
      { 'CD1', 'CE1', 'CZ' },
      { 'CG', 'CD1', 'CE1', 'CZ' },
      { 'CE1', 'CZ', 'OH' },
      { 'CD1', 'CE1', 'CZ', 'OH' },
      { 'CB', 'CG', 'CD2' },
      { 'CA', 'CB', 'CG', 'CD2' },
      { 'CG', 'CD2', 'CE2' },
      { 'CB', 'CG', 'CD2', 'CE2' },
      ['chi1'] = 2, ['chi2'] = 4
   },
   ['W'] = {
      { 'CA', 'CB', 'CG' },
      { 'N', 'CA', 'CB', 'CG' },   -- chi1
      { 'CB', 'CG', 'CD1' },
      { 'CA', 'CB', 'CG', 'CD1' },   -- chi2
      { 'CG', 'CD1','NE1' },
      { 'CB', 'CG', 'CD1','NE1' },
      { 'CB', 'CG', 'CD2' },
      { 'CA', 'CB', 'CG', 'CD2' },
      { 'CG', 'CD2','CE2' },
      { 'CB', 'CG', 'CD2','CE2' },
      { 'CD2','CE2','CZ2' },
      { 'CG', 'CD2','CE2','CZ2' },
      { 'CE2','CZ2','CH2' },
      { 'CD2','CE2','CZ2','CH2' },
      { 'CG', 'CD2','CE3' },
      { 'CB', 'CG', 'CD2','CE3' },
      { 'CD2','CE3','CZ3' },
      { 'CG','CD2','CE3','CZ3' },
      ['chi1'] = 2, ['chi2'] = 4
   },
   ['D'] = {
      { 'CA', 'CB', 'CG' },
      { 'N', 'CA', 'CB', 'CG' },   -- chi1
      { 'CB', 'CG', 'OD1' },
      { 'CA', 'CB', 'CG', 'OD1' },   -- chi2
      { 'CB', 'CG', 'OD2' },
      { 'CA', 'CB', 'CG', 'OD2' },
      ['chi1'] = 2, ['chi2'] = 4
   },
   ['E'] = {
      { 'CA', 'CB', 'CG' },
      { 'N', 'CA', 'CB', 'CG' },   -- chi1
      { 'CB', 'CG', 'CD' },
      { 'CA', 'CB', 'CG', 'CD' },   -- chi2
      { 'CG', 'CD', 'OE1' },
      { 'CB', 'CG', 'CD', 'OE1' },   -- chi3
      { 'CG', 'CD', 'OE2' },
      { 'CB', 'CG', 'CD', 'OE2' },
      ['chi1'] = 2, ['chi2'] = 4, ['chi3'] = 6
   },
   ['H'] = {
      { 'CA', 'CB', 'CG' },
      { 'N', 'CA', 'CB', 'CG' },   -- chi1
      { 'CB', 'CG', 'ND1' },
      { 'CA', 'CB', 'CG', 'ND1' },   -- chi2
      { 'CG', 'ND1','CE1' },
      { 'CB', 'CG', 'ND1','CE1' },
      { 'CB', 'CG', 'CD2' },
      { 'CA', 'CB', 'CG', 'CD2' },
      { 'CG', 'CD2','NE2' },
      { 'CB', 'CG', 'CD2','NE2' },
      ['chi1'] = 2, ['chi2'] = 4
   },
   ['K'] = {
      { 'CA', 'CB', 'CG' },
      { 'N', 'CA', 'CB', 'CG' },   -- chi1
      { 'CB', 'CG', 'CD' },
      { 'CA', 'CB', 'CG', 'CD' },   -- chi2
      { 'CG', 'CD', 'CE' },
      { 'CB', 'CG', 'CD', 'CE' },   -- chi3
      { 'CD', 'CE', 'NZ' },
      { 'CG', 'CD', 'CE', 'NZ' },   -- chi4
      ['chi1'] = 2, ['chi2'] = 4, ['chi3'] = 6, ['chi4'] = 8
   },
   ['R'] = {
      { 'CA', 'CB', 'CG' },
      { 'N', 'CA', 'CB', 'CG' },   -- chi1
      { 'CB', 'CG', 'CD' },
      { 'CA', 'CB', 'CG', 'CD' },   -- chi2
      { 'CG', 'CD', 'NE' },
      { 'CB', 'CG', 'CD', 'NE' },   -- chi3
      { 'CD', 'NE', 'CZ' },
      { 'CG', 'CD', 'NE', 'CZ' },   -- chi4
      { 'NE', 'CZ', 'NH1' },
      { 'CD', 'NE', 'CZ', 'NH1' },   -- chi5
      { 'NE', 'CZ', 'NH2' },
      { 'CD', 'NE', 'CZ', 'NH2' },   -- chi5
      ['chi1'] = 2, ['chi2'] = 4, ['chi3'] = 6, ['chi4'] = 8, ['chi5'] = 10
   }
}
   
      
   
   
return utils

