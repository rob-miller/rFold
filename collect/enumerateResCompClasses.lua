#!/usr/bin/env luajit

require 'firePostgresql'
dbReconnect(1)

local update=true

aa = { 'g', 'a', 'v', 'l', 'i', 'm', 'f', 'p', 's', 't', 'c', 'n', 'q', 'y', 'w', 'd', 'e', 'h', 'k', 'r' }

atom_class = { c1 = 0.77, cr = 0.72, c2 = 0.67, o1 = 0.67, o2 = 0.60, n1 = 0.70, n2 = 0.62, h1 = 0.37, s1 = 1.04 }

periodic_table = { h = 1.008, c = 12.011, n = 14.007, o = 15.99, s = 32.06 }

back_atoms = { n = 'n1', ca = 'c1', c = 'c2', o = 'o2' }
back_bonds = { {'_c','n'}, {'n','ca'}, {'ca','c'}, {'c','o'}, {'c','n_'} }

mc_angle = { {'_ca','_c','_n'}, {'_c','n','ca'}, {'n','ca','c'}, {'ca','c','n_'}, {'ca','c','o'} }

sc_atoms = { g = { hb = 'h1' }, a = { cb = 'c1' }, v = { cb = 'c1', cg1 = 'c1', cg2 = 'c1' }, l = { cb = 'c1', cg = 'c1', cd1 = 'c1', cd2 = 'c1' },
             i = { cb = 'c1', cg1 = 'c1', cg2 = 'c1', cd1 = 'c1' }, m = { cb = 'c1', cg = 'c1', sd = 's1', ce = 'c1'  },
             f = { cb = 'c1', cg = 'c2', cd1 = 'cr', cd2 = 'cr', ce1 = 'c2', ce2 = 'c2', cz = 'cr' }, p = { cb = 'c1', cg = 'c1', cd = 'c1' },
             s = { cb = 'c1', og = 'o1' }, t = { cb = 'c1', og1 = 'o1', cg2 = 'c1' }, c = { cb = 'c1', sg = 's1' },
             n = { cb = 'c1', cg = 'c2', od1 = 'o2', nd2 = 'n1' }, q = { cb = 'c1', cg = 'c1', cd = 'c2', oe1 = 'o2', ne2 = 'n1' },
             y = { cb = 'c1', cg = 'c2', cd1 = 'cr', cd2 = 'cr', ce1 = 'c2', ce2 = 'c2', cz ='cr', oh = 'o1' },
             w = { cb = 'c1', cg = 'c2', cd1 = 'c2', cd2 = 'cr', ne1 = 'n1', ce2 = 'c2', ce3 = 'c2', cz2 = 'cr', cz3 = 'cr', ch2 = 'c2' },
             d = { cb = 'c1', cg = 'c2', od1 = 'o2', od2 = 'o1' }, e = { cb = 'c1', cg = 'c1', cd = 'c2', oe1 = 'o2', oe2 = 'o1' },
             h = { cb = 'c1', cg = 'c2', nd1 = 'n1', cd2 = 'c2', ce1 = 'c2', ne2 = 'n2' }, k = { cb = 'c1', cg = 'c1', cd = 'c1', ce = 'c1', nz = 'n1' },
             r = { cb = 'c1', cg = 'c1', cd = 'c1', ne = 'n1', cz = 'c2', nh1 = 'n2', nh2 = 'n1' } }

-- note protonation states not addressed for d,e,h,k,r atom classes!  could determine based on planarity of relevant atoms ?


sc_bonds = { g = {{'ca','hb'}}, a = {{'ca','cb'}}, v = {{'ca','cb'},{'cb','cg1'},{'cb','cg1'}}, l = {{'ca','cb'},{'cb','cg'},{'cg','cd1'},{'cg','cd2'}},
             i = {{'ca','cb'},{'cb','cg1'},{'cb','cg2'},{'cg1','cd1'}}, m = {{'ca','cb'},{'cb','sd'},{'sd','ce'}},
             f = {{'ca','cb'},{'cb','cg'},{'cg','cd1'},{'cd1','ce1'},{'ce1','cz'},{'cz','ce2'},{'ce2','cd2'},{'cd2','cg'}},
             p = {{'ca','cb'},{'cb','cg'},{'cg','cd'},{'cd','n'}}, s = {{'ca','cb'},{'cb','og'}}, t = {{'ca','cb'},{'cb','og1'},{'cb','cg2'}},
             c = {{'ca','cb'},{'cb','sg'}}, n = {{'ca','cb'},{'cb','cg'},{'cg','od1'},{'cg','nd2'}}, q = {{'ca','cb'},{'cb','cg'},{'cg','cd'},{'cd','oe1'},{'cg','ne2'}},
             y = {{'ca','cb'},{'cb','cg'},{'cg','cd1'},{'cd1','ce1'},{'ce1','cz'},{'cz','ce2'},{'ce2','cd2'},{'cd2','cg'},{'cz','oh'}},
             w = {{'ca','cb'},{'cb','cg'},{'cg','cd1'},{'cd1','ne1'},{'ne1','ce2'},{'ce2','cz2'},{'cz2','ch2'},{'ch2','cz3'},{'cz3','ce3'},{'ce3','cd2'},{'cd2','cg'},{'cd2','ce2'}}
             d = {{'ca','cb'},{'cb','cg'},{'cg','od1'},{'cg','od2'}}, e = {{'ca','cb'},{'cb','cg'},{'cg','cd'},{'cd','oe1'},{'cg','oe2'}},
             h = {{'ca','cb'},{'cb','cg'},{'cg','nd1'},{'nd1','ce1'},{'ce','ne2'},{'ne2','cd2'},{'cd2','cg'}}
             k = {{'ca','cb'},{'cb','cg'},{'cg','cd'},{'cd','ce'},{'ce','nz'}}, r = {{'ca','cb'},{'cb','cg'},{'cg','cd'},{'cd','ne'},{'ne','cz'},{'cz','nh1'},{'cz','nh2'}} }



-- mc_angle = { {'_ca','_c','_n'}, {'_c','n','ca'}, {'n','ca','c'}, {'ca','c','n_'}, {'ca','c','o'}, {'cb','ca','c'}, {'n','ca','cb'} }

mc_dihed = { omg = { '_ca', '_c', 'n', 'ca' }, phi = { '_c', 'n', 'ca', 'c' }, psi = { 'n', 'ca', 'c', 'n_' } }
mc_dihed_2 = { psi2 = { 'n', 'ca', 'c', 'o' }, psi3 = { 'cb', 'ca', 'c', 'o' }, psi2 = { 'n', 'ca', 'cb', 'c' } }
-- psi2, psi3 follow psi but needed to build carbonyl o and c-beta

sc_angle = { 
}
sc_dihed = { v = { [1] = { 'n', 'ca', 'cb', 'cg1' } },
             l = { [1] = { 'n', 'ca', 'cb', 'cg' },  [2] = { 'ca', 'cb', 'cg', 'cd1' } },
             i = { [1] = { 'n', 'ca', 'cb', 'cg1' }, [2] = { 'ca', 'cb', 'cg1', 'cd' } },
             m = { [1] = { 'n', 'ca', 'cb', 'cg' },  [2] = { 'ca', 'cb', 'cg', 'sd' } },
             f = { [1] = { 'n', 'ca', 'cb', 'cg' },  [2] = { 'ca', 'cb', 'cg', 'cd1' } },
             p = { [1] = { 'n', 'ca', 'cb', 'cg' },  [2] = { 'ca', 'cb', 'cg', 'cd' } },
             s = { [1] = { 'n', 'ca', 'cb', 'og' } },
             t = { [1] = { 'n', 'ca', 'cb', 'og1' } },
             c = { [1] = { 'n', 'ca', 'cb', 'sg' } },
             n = { [1] = { 'n', 'ca', 'cb', 'cg' },  [2] = { 'ca', 'cb', 'cg', 'od1' } },
             q = { [1] = { 'n', 'ca', 'cb', 'cg' },  [2] = { 'ca', 'cb', 'cg', 'cd' }, ['3'] = { 'cb', 'cg', 'cd', 'oe1' } },
             y = { [1] = { 'n', 'ca', 'cb', 'cg' },  [2] = { 'ca', 'cb', 'cg', 'cd1' } },
             w = { [1] = { 'n', 'ca', 'cb', 'cg' },  [2] = { 'ca', 'cb', 'cg', 'cd1' } },
             d = { [1] = { 'n', 'ca', 'cb', 'cg' },  [2] = { 'ca', 'cb', 'cg', 'od1' } },
             e = { [1] = { 'n', 'ca', 'cb', 'cg' },  [2] = { 'ca', 'cb', 'cg', 'cd' }, ['3'] = { 'cb', 'cg', 'cd', 'oe1' } },
             h = { [1] = { 'n', 'ca', 'cb', 'cg' },  [2] = { 'ca', 'cb', 'cg', 'nd1' } },
             k = { [1] = { 'n', 'ca', 'cb', 'cg' },  [2] = { 'ca', 'cb', 'cg', 'cd' }, ['3'] = { 'cb', 'cg', 'cd', 'ce' }, [4] = { 'cb', 'cg', 'ce', 'nz' } },
             r = { [1] = { 'n', 'ca', 'cb', 'cg' },  [2] = { 'ca', 'cb', 'cg', 'cd' }, ['3'] = { 'cb', 'cg', 'cd', 'ne' }, [4] = { 'cb', 'cg', 'ne', 'cz' }, [5] = {  'cg', 'ne', 'cz', 'nh1' } } }



for k,v in pairs(periodic_table) do
   s = "insert into periodic_table ( atom, weight ) values ( '".. k .."', ".. v .." ) on conflict(atom) " .. (update and "do update set weight = " .. v or "do nothing")
   --print(s)
   --pgQcur(s)
end

for k,v in pairs(atom_class) do
   s = "insert into atom_class ( class, r_covalent ) values ( '".. k .."', ".. v .." ) on conflict(class) " .. (update and "do update set r_covalent = " .. v or "do nothing")
   --print(s)
   --pgQcur(s)
end

for k,v in pairs(back_atoms) do
   local atom = string.sub(k,1,1)
   s = "insert into atoms ( name, class, atom ) values ( '" .. k .. "','" .. v .. "','" .. atom .. "' ) on conflict(name) " .. (update and "do update set class = '" .. v .. "', atom = '" .. atom .. "'" or "do nothing")
   --print( s )
   --pgQcur(s)
end

for i,r in ipairs(aa) do
   for k,v in pairs(back_atoms) do
      local atom = string.sub(k,1,1)
      s = "insert into atoms ( name, class, atom ) values ( '" .. r..k .. "','" .. v .. "','" .. atom .. "' ) on conflict(name) " .. (update and "do update set class = '" .. v .. "', atom = '" .. atom .. "'" or "do nothing")
      --print( s )
      --pgQcur(s)
   end
   local sc = sc_atoms[r]
   for k,v in pairs(sc) do
      local atom = string.sub(k,1,1)
      s = "insert into atoms ( name, class, atom ) values ( '" .. r..k .. "','" .. v .. "','" .. atom .. "' ) on conflict(name) " .. (update and "do update set class = '" .. v .. "', atom = '" .. atom .. "'" or "do nothing")
      --print( s )
      --pgQcur(s)
   end      
end


--[[
for i,res in pairs(sc_atoms) do
   
   print(i,res)
end
]]--

