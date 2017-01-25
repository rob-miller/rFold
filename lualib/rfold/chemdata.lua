--[[
   chemdata.lua
   
Copyright 2016 Robert T. Miller

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

-- chemistry data including symbols, relationships and constants

local chemdata = {}

function chemdata.covalent_volume(atomClass) 
   return (4/3)*math.pi*chemdata.covalent_radii[atomClass]
end

   



---------------------------------------------- tables -------------------------------------------------------------------------------------------------------------
-- @section tables

--- residue single letter to 3-letter name conversion
-- @table res3
chemdata.res3 = { G = 'GLY', A = 'ALA', V = 'VAL', L = 'LEU', I = 'ILE', M = 'MET', F = 'PHE', P = 'PRO', S = 'SER', T = 'THR',
               C = 'CYS', N = 'ASN', Q = 'GLN', Y = 'TYR', W = 'TRP', D = 'ASP', E = 'GLU', H = 'HIS', K = 'LYS', R = 'ARG',
               X = 'UNK' }

--- residue 3-letter name to single letter conversion
-- @table res1
chemdata.res1 = { ['GLY'] = 'G', ['ALA'] = 'A', ['VAL'] = 'V', ['LEU'] = 'L', ['ILE'] = 'I', ['MET'] = 'M', ['PHE'] = 'F', ['PRO'] = 'P', ['SER'] = 'S', ['THR'] = 'T',
               ['CYS'] = 'C', ['ASN'] = 'N', ['GLN'] = 'Q', ['TYR'] = 'Y', ['TRP'] = 'W', ['ASP'] = 'D', ['GLU'] = 'E', ['HIS'] = 'H', ['LYS'] = 'K', ['ARG'] = 'R',
               ['UNK'] = 'X' }


--- per residue sidechain hedra and dihedra definitions, in order of output for internal coordinates spcification file
-- @table sidechains
chemdata.sidechains = {
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

-- covalent radii from Heyrovska, Raji : 'Atomic Structures of all the Twenty Essential Amino Acids and a Tripeptide, with Bond Lengths as Sums of Atomic Covalent Radii'
-- https://arxiv.org/pdf/0804.2488.pdf

-- adding Ores between Osb and Odb for Asp and Glu, Nres between Nsb and Ndb for Arg, as PDB does not specify

-- http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/formuleAA/ for naming of individual atoms

chemdata.covalent_radii = {
   ['Csb'] = 0.77, ['Cres'] = 0.72, ['Cdb'] = 0.67, ['Osb'] = 0.67, ['Ores'] = 0.635, ['Odb'] = 0.60, ['Nsb'] = 0.70, ['Nres'] = 0.66, ['Ndb'] = 0.62, ['Hsb'] = 0.37, ['Ssb'] = 1.04
}

chemdata.atomic_weight = { ['C'] = 12.0107, ['O'] = 15.9994, ['N'] = 14.0067, ['H'] = 1.0079, ['S'] = 32.065 }
chemdata.electronegativity = { ['C'] = 2.55, ['O'] = 3.44, ['N'] = 3.04, ['H'] = 2.20, ['S'] = 2.58 }

chemdata.residue_atom_class = {
   ['X'] = { ['N'] = 'Nsb', ['CA'] = 'Csb', ['C'] = 'Cdb', ['O'] = 'Odb', ['OXT'] = 'Osb', ['CB'] = 'Csb' },
   ['V'] = { ['CG1'] = 'Csb', ['CG2'] = 'Csb' },
   ['L'] = { ['CG'] = 'Csb', ['CD1'] = 'Csb', ['CD2'] = 'Csb' },
   ['I'] = { ['CG1'] = 'Csb', ['CG2'] = 'Csb', ['CD1'] = 'Csb' },
   ['M'] = { ['CG'] = 'Csb', ['SD'] = 'Ssb', ['CE'] = 'Csb' },
   ['F'] = { ['CG'] = 'Cdb', ['CD1'] = 'Cres', ['CD2'] = 'Cres', ['CE1'] = 'Cdb', ['CE2'] = 'Cdb', ['CZ'] = 'Cres' },
   ['P'] = { ['CG'] = 'Csb' },
   ['S'] = { ['OG'] = 'Osb' },
   ['T'] = { ['OG1'] = 'Osb', ['CG2'] = 'Csb' },
   ['C'] = { ['SG'] = 'Ssb' },
   ['N'] = { ['CG'] = 'Csb', ['OD1'] = 'Odb', ['ND2'] = 'Ndb' },
   ['Q'] = { ['CG'] = 'Csb', ['CD'] = 'Csb', ['OE1'] = 'Odb', ['NE2'] = 'Ndb' },
   ['Y'] = { ['CG'] = 'Cdb', ['CD1'] = 'Cres', ['CD2'] = 'Cres', ['CE1'] = 'Cdb', ['CE2'] = 'Cdb', ['CZ'] = 'Cres', ['OH'] = 'Osb' },
   ['W'] = { ['CG'] = 'Cdb', ['CD1'] = 'Cdb', ['CD2'] = 'Cres', ['NE1'] = 'Nsb', ['CE2'] = 'Cdb', ['CE3'] = 'Cdb', ['CZ2'] = 'Cres', ['CZ3'] = 'Cres', ['CH2'] = 'Cdb' }, 
   ['D'] = { ['CG'] = 'Csb', ['OD1'] = 'Ores', ['OD2'] = 'Ores' },
   ['E'] = { ['CG'] = 'Csb', ['CD'] = 'Csb', ['OE1'] = 'Ores', ['OE2'] = 'Ores' },
   ['H'] = { ['CG'] = 'Cdb', ['CD1'] = 'Cdb', ['ND1'] = 'Nsb', ['CE1'] = 'Cdb', ['NE2'] = 'Ndb' },
   ['K'] = { ['CG'] = 'Csb', ['CD'] = 'Csb', ['CE'] = 'Csb', ['NZ'] = 'Nsb' },
   ['R'] = { ['CG'] = 'Csb', ['CD'] = 'Csb', ['NE'] = 'Nsb', ['CZ'] = 'Cdb', ['NH1'] = 'Nres', ['NH2'] = 'Nres' }
}

return chemdata
