#!/usr/bin/env luajit

require 'rtmLualib'
require 'parseCmdLine'
require 'parseRtmDssp'
require 'rtmGeom3d'

local triples
local tripMap = {}
local quads
local quadMap = {}
local dsps

function genKey(a1,a2,a3,a4)
   if a4 then 
      return a1 .. ':' .. a2 .. ':' .. a3 .. ':' .. a4
   end
   return a1 .. ':' .. a2 .. ':' .. a3
end

function tripleInitPos(trip)
   -- generate coordinates for 3 atoms with specified bond lengths and angle between on XZ plane (Y=0 for all atoms)
   -- match aacs : a2 = 0,0,0 ; a3 on Z-axis ; a1 on XZ plane
   -- reverse coordinates swap : a2 = 0,0,0 ; a1 on Z-axis ; a3 on XZ plane
   trip['a1initial'] = get41mtx()   
   trip['a2initial'] = get41mtx()    -- a2 to 0,0,0
   trip['a3initial'] = get41mtx()

   trip['a3initial'][3][1] = trip['len3']   -- a3 is len3 up from a2 on Z axis, X=Y=0

   --local ar = trip['angle2'] * d2r
   local sar = (180.0 - trip['angle2']) * d2r    -- angles which add to 180 are supplementary

   trip['a1initial'][1][1] = math.sin(sar) * trip['len1']      -- a1 X is sin( sar ) * len1
   trip['a1initial'][3][1] = - (math.cos(sar) * trip['len1'])  -- a1 Z is -(cos( sar ) * len1) (assume angle2 always obtuse, so a1 is in -Z

   -- same again but a3 at 0,0,0
   trip['a1initialr'] = get41mtx()   
   trip['a2initialr'] = get41mtx()   -- a2 to 0,0,0
   trip['a3initialr'] = get41mtx()    

   trip['a1initialr'][3][1] = trip['len1']   -- a1r is len1 up from a2 on Z axis, X=Y=0

   --local ar = trip['angle2'] * d2r
   local sar = (180.0 - trip['angle2']) * d2r    -- angles which add to 180 are supplementary

   trip['a3initialr'][1][1] = math.sin(sar) * trip['len3']    -- a3r X is sin( sar ) * len3
   trip['a3initialr'][3][1] = - (math.cos(sar) * trip['len3'])  -- a3r Z is -(cos( sar ) * len3)

   return trip
end

function quadInitPos(quad)
   -- generate coordinates for 4 atoms with specified dihedral, first 3 on XZ plane, 4th in +Z
   local t1key = genKey(quad[1], quad[2], quad[3])
   local trip1 = triples[ tripMap[ t1key ] ]
   if not trip1 then return quadInitPosR(quad) end
   --assert(trip1,'failed to locate 1st triple ' .. t1key .. ' a4= ' .. quad['a4'])

   local t2key = genKey(quad['a2'], quad['a3'], quad['a4'])
   local trip2 = triples[ tripMap[ t2key ] ]
   assert(trip2,'failed to locate 2nd triple ' .. t2key .. ' a1= ' .. quad[1])

   quad['a1initial'] = trip1['a1initial']:copy()
   quad['a2initial'] = trip1['a2initial']:copy()
   quad['a3initial'] = trip1['a3initial']:copy()

   quad['a4initial'] = trip2['a3initialr']:copy()
   quad['a4initial'][3][1] = quad['a4initial'][3][1] * -1                   -- a4 to +Z
   quad['a4initial'][3][1] = quad['a4initial'][3][1] + trip2['len1']        -- trip2 shift up so a2 at 0,0,0
   
   local mrz = genMrz( quad['dihedral1'] * d2r )
   quad['a4initial'] = mrz * quad['a4initial']                              -- dihedral set

   -- print('qip : ' .. genKey(quad['a1'], quad['a2'], quad['a3'],quad['a4']) .. ' ' .. t1key .. ' ' .. t2key .. ' [' ..  quad['a1initial']:transpose():pretty() .. '][' ..  quad['a2initial']:transpose():pretty() .. '][' ..  quad['a3initial']:transpose():pretty() .. '][' ..  quad['a4initial']:transpose():pretty() .. ']')


   -- local mt,mr = coordSpace
   
   return quad
end

function quadInitPosR(quad)
   -- generate coordinates for 4 atoms with specified dihedral, first 3 on XZ plane with Z>=0, 4th in +Z -- triples in reverse order from above
   local t1key = genKey(quad[3], quad[2], quad[1])
   local trip1 = triples[ tripMap[ t1key ] ]
   assert(trip1,'failed to locate 1st triple ' .. t1key .. ' a4= ' .. quad[4])

   local t2key = genKey(quad['a4'], quad['a3'], quad['a2'])
   local trip2 = triples[ tripMap[ t2key ] ]
   assert(trip2,'failed to locate 2nd triple ' .. t2key .. ' a1= ' .. quad[1])

   --print('t1k: ' .. t1key .. '  t2k: ' .. t2key)
   
   
   quad['a1initial'] = trip1['a3initialr']:copy()
   quad['a2initial'] = trip1['a2initialr']:copy()
   quad['a3initial'] = trip1['a1initialr']:copy()

   quad['a4initial'] = trip2['a1initial']:copy()
   
   quad['a4initial'][3][1] = quad['a4initial'][3][1] * -1                   -- a4 to +Z
   quad['a4initial'][3][1] = quad['a4initial'][3][1] + trip2['len3']        -- trip2 shift up so a2 at 0,0,0  -- reverse from above so  use len3
   
   local mrz = genMrz( quad['dihedral1'] * d2r )
   quad['a4initial'] = mrz * quad['a4initial']                               -- dihedral set

   --print('qipR: ' .. genKey(quad['a1'], quad['a2'], quad['a3'],quad['a4']) .. ' ' .. t1key .. ' ' .. t2key .. ' [' ..  quad['a1initial']:pretty() .. '][' ..  quad['a2initial']:pretty() .. '][' ..  quad['a3initial']:pretty() .. '][' ..  quad['a4initial']:pretty() .. ']') 
   return quad
end

--[[
   -- moved to rtmGeom3d
   
function coordSpace(a1,a2,a3, r)
   -- return transform matrix placing a1 on XZ plane, a2 at 0,0,0, a3 at 0,0,+Z
   -- a's need to be 4x1 matrices
   -- if r then also generate reverse transformation matrix to bring back 

   local dbg=nil
   if dbg then
      print(a1:transpose():pretty())
      print(a2:transpose():pretty())
      print(a3:transpose():pretty())
   end
   
   local tm = genMt(-a2[1][1], -a2[2][1], -a2[3][1])  -- translation matrix for a2 to origin
      
   -- now get a3 to Z-axis
   local p3 = get41mtx()   -- direct translation of a3 using a2 to origin
   p3[1][1] = a3[1][1] - a2[1][1]
   p3[2][1] = a3[2][1] - a2[2][1]
   p3[3][1] = a3[3][1] - a2[3][1]

   local sc = getSphericalCoordinates(p3)

   if dbg then
      print(sc:pretty())
   end
   
   local mrz = genMrz(-sc[1][1])  -- rotate translated a3 -theta about Z
   local mry = genMry(-sc[2][1])  -- rotate translated a3 -phi about Y
   
   local mt = mrz * tm
   mt = mry * mt                  -- mt completes a2-a3 on Z-axis, still need to align a1 with XZ plane

   if dbg then
      print((mt * a3):pretty())
   end
   
   p3 = mt * a1                   -- p3 not needed, re-use it here

   local sc2 = getSphericalCoordinates(p3)    -- need theta of translated a1
   local mrz2 = genMrz(-sc2[1][1])            -- rotate a1 -theta about Z to align with X

   mt = mrz2 * mt                 -- mt transforms to specified coordinate space

   if not r then return mt end

   -- generate the reverse transformation

   mrz2 = genMrz(sc2[1][1])            -- rotate a1 theta about Z, reversing alignment with X
   mry  = genMry(sc[2][1])             -- rotate a3 phi about Y
   mrz  = genMrz(sc[1][1])             -- rotate a3 theta about Z
   tm   = genMt(a2[1][1], a2[2][1], a2[3][1])  -- translation matrix for origin to a2

   local mr = mry * mrz2
   mr = mrz * mr
   mr = tm * mr

   return mt,mr
end
--]]

function testDihed(dk,p1,p2,p3,p4)
   local m1 = get41mtx()
   local m2 = get41mtx()
   local m3 = get41mtx()
   local m4 = get41mtx()

   --print( p1['x'] .. ' ' .. p1['y'] .. ' ' .. p1['z'] )
   m1[1][1] = p1['x']; m1[2][1] = p1['y']; m1[3][1] = p1['z']; 
   m2[1][1] = p2['x']; m2[2][1] = p2['y']; m2[3][1] = p2['z']; 
   m3[1][1] = p3['x']; m3[2][1] = p3['y']; m3[3][1] = p3['z']; 
   m4[1][1] = p4['x']; m4[2][1] = p4['y']; m4[3][1] = p4['z']; 

   -- get transform to put p1,p2,p3 into standard coordinate space (and reverse transform)
   
   mt,mtr = coordSpace(m1,m2,m3,true)

   -- transform dssp coords to standard coord space 
   
   local m1t = mt * m1
   local m2t = mt * m2
   local m3t = mt * m3
   local m4t = mt * m4

   -- (reverse) transform standard coord space dihedral atoms to match dssp coords
   local quad = quads[ quadMap[dk] ]

   local d1r = mtr * quad['a1initial']
   local d2r = mtr * quad['a2initial']
   local d3r = mtr * quad['a3initial']
   local d4r = mtr * quad['a4initial']

   -- reverse transform the transformed dssp coords back to original positions
   local m1r = mtr * m1t
   local m2r = mtr * m2t
   local m3r = mtr * m3t
   local m4r = mtr * m4t

   -- get distance between transformed dssp atom 4 and generated dihedral atom 4
   
   local dist = getDistance3d(m4t,quad['a4initial']) 
   print(dk .. ' ' .. dist)
  
   if dist > 1.0e-04 then  -- if not close enough
      print(' [' ..  m1:transpose():pretty() .. '][' ..  m2:transpose():pretty() .. '][' ..  m3:transpose():pretty() .. '][' ..  m4:transpose():pretty() .. ']')   -- print dssp coords 
      if not (matrixAreEqual(m1,m1r) and matrixAreEqual(m2,m2r) and matrixAreEqual(m3,m3r) and matrixAreEqual(m4,m4r)) then             -- if dssp coords transformed back do not match originals, print out
         print('x[' ..  m1r:transpose():pretty() .. '][' ..  m2r:transpose():pretty() .. '][' ..  m3r:transpose():pretty() .. '][' ..  m4r:transpose():pretty() .. ']')
      end
      print(' [' ..  m1t:transpose():pretty() .. '][' ..  m2t:transpose():pretty() .. '][' ..  m3t:transpose():pretty() .. '][' ..  m4t:transpose():pretty() .. ']')   -- print dssp coords in standard system

      -- print generated dihedral coords
      print(dk .. ' : ' .. ' [' ..  quad['a1initial']:transpose():pretty() .. '][' ..  quad['a2initial']:transpose():pretty() .. '][' ..  quad['a3initial']:transpose():pretty() .. '][' ..  quad['a4initial']:transpose():pretty() .. ']')
      print('[' ..  d1r:transpose():pretty() .. '][' ..  d2r:transpose():pretty() .. '][' ..  d3r:transpose():pretty() .. '][' ..  d4r:transpose():pretty() .. ']') -- print generated dihedrals transformed to dssp coordinate space

      print('--')
   end
   
end

function testDssp(i)
--   dsps[i-1],j,dsps[i+1]
--[[
   local rn={}
   local r={}
   rn[1] = d1['resn']
   rn[2] = d2['resn']
   rn[3] = d3['resn']
   r[1] = d1['res']
   r[2] = d2['res']
   r[3] = d3['res']

   print(r[1] .. rn[1] .. '-' .. r[2] .. rn[2] .. '-' .. r[3] .. rn[3])
--]]

   -- omega : d1 ca .. d2 ca                           K1CA:K1C:E2N:E2CA
   -- phi   : d1 c  .. d2 c                            K1C:E2N:E2CA:E2C
   -- psi   : d2 n  .. d3 n       K1N:K1CA:K1C:E2N
   -- a1    : d2 n  .. d2 o       K1N:K1CA:K1C:K1O
   -- a2    : d2 o  .. d2 cb      K1O:K1C:K1CA:K1CB

   local d1
   local d2 = dsps[i]
   local d3

   --print('i= ' .. i)
   if dsps[i-1] then      
      d1 = dsps[i-1]
      local omgk = genKey( d1['res']..d1['resn']..'CA', d1['res']..d1['resn']..'C', d2['res']..d2['resn']..'N',  d2['res']..d2['resn']..'CA' )
      print('omega')
      testDihed(omgk,d1['ca'],d1['c'],d2['n'],d2['ca'])
      print('phi')
      local phik = genKey( d1['res']..d1['resn']..'C', d2['res']..d2['resn']..'N', d2['res']..d2['resn']..'CA',  d2['res']..d2['resn']..'C' )
      testDihed(phik,d1['c'],d2['n'],d2['ca'],d2['c'])
   end

   local a1k = genKey( d2['res']..d2['resn']..'N', d2['res']..d2['resn']..'CA', d2['res']..d2['resn']..'C',  d2['res']..d2['resn']..'O' )
   print('a1')
   testDihed(a1k,d2['n'],d2['ca'],d2['c'],d2['o'])
   local a2k = genKey( d2['res']..d2['resn']..'O', d2['res']..d2['resn']..'C', d2['res']..d2['resn']..'CA',  d2['res']..d2['resn']..'CB' )
   print('a2')
   testDihed(a2k,d2['o'],d2['c'],d2['ca'],d2['cb'])

   if dsps[i+1] then
      d3 = dsps[i+1]
      local psik = genKey( d2['res']..d2['resn']..'N', d2['res']..d2['resn']..'CA', d2['res']..d2['resn']..'C',  d3['res']..d3['resn']..'N' )
      print('psi')
      testDihed(psik,d2['n'],d2['ca'],d2['c'],d3['n'])
   end
end
   
local dsspdat

local args = parseCmdLine();
dsspdat = prd(args[1])

triples = dsspdat['triples']

for i,j in ipairs(triples) do
   triples[i] = tripleInitPos(j)
   tripMap[ genKey(j[1],j[2],j[3]) ] = i
end

quads = dsspdat['quads']

for i,j in ipairs(quads) do
   quads[i] = quadInitPos(j)
   quadMap[ genKey(j[1],j[2],j[3],j[4]) ] = i
end

dsps = dsspdat['dssp']

for i,j in ipairs(dsps) do
   --if i>1 and i<#dsps then
      testDssp(i)
   --end
end


--[[
for i,j in ipairs(triples) do
   print('['.. j['a1initial']:pretty() .. '] [' .. j['a2initial']:pretty() .. '] [' .. j['a3initial']:pretty() .. '] ' .. j['len1'] .. ' ' .. j['angle2'] .. ' '  .. j['len3'] )
   --j['a1initial']:print()
   --print( j['a1initial']:pretty())
end
--]]
