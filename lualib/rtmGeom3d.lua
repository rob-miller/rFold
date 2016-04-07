#!/usr/bin/env luajit

require 'numlua'


-- degrees to radians, radians to degrees constants
d2r = (math.pi/180.0)
r2d = (180.0/math.pi)

function get14mtx()
   local m = matrix.zeros(1,4)
   m[1][4]=1.0
   return m
end

function get41mtx()
   local m = matrix.zeros(4,1)
   m[4][1]=1.0
   return m
end

function get44mtx()
   local m = matrix.zeros(4,4)
   m[4][4]=1.0
   return m
end

function genMt(x,y,z)
   local mt = get44mtx()
   mt[1][1] = 1.0
   mt[2][2] = 1.0
   mt[3][3] = 1.0

   mt[1][4] = x
   mt[2][4] = y
   mt[3][4] = z
   return mt
end


function genMrx(angr)
   local mrx = get44mtx()
   local cosang = math.cos(angr)
   local sinang = math.sin(angr)
   mrx[1][1] = 1.0
   mrx[2][2] = cosang
   mrx[3][3] = cosang
   mrx[2][3] = -sinang
   mrx[3][2] = sinang
   return mrx
end

function genMry(angr)
   local mry = get44mtx()
   local cosang = math.cos(angr)
   local sinang = math.sin(angr)
   mry[1][1] = cosang
   mry[2][2] = 1.0
   mry[3][3] = cosang
   mry[3][1] = -sinang
   mry[1][3] = sinang
   return mry
end

function genMrz(angr)
   local mrz = get44mtx()
   local cosang = math.cos(angr)
   local sinang = math.sin(angr)
   mrz[1][1] = cosang
   mrz[2][2] = cosang
   mrz[3][3] = 1.0
   mrz[1][2] = -sinang
   mrz[2][1] = sinang
   return mrz
end

function getSphericalCoordinates(m)
   -- given m[ x y z ], return m[ theta phi distance ] -- angles in radians
   local retm = get41mtx()
   -- distance
   retm[3][1] = math.sqrt( (m[1][1] * m[1][1]) + (m[2][1] * m[2][1]) + (m[3][1] * m[3][1]) )
   -- theta
   local sign = (m[2][1] < 0.0 and -1.0) or 1.0
   retm[1][1] = (m[1][1] == 0.0 and ((math.pi/2.0) * sign)) or (math.atan( m[2][1]/m[1][1] ))   -- if x = 0 then theta = 90 else theta = atan(y/x)
   if m[1][1] < 0.0 then     -- extend sign to 180 degrees
      retm[1][1] = retm[1][1] + (sign * math.pi)
      if retm[1][1] == -math.pi then retm[1][1] = math.pi end
   end
   -- phi
   if 0.0 == retm[3][1] then
      retm[2][1] = 0.0     -- if dist = 0 then point is origin and phi = 0
   else
      retm[2][1] = math.acos(m[3][1] / retm[3][1])     -- else phi = acos(z/distance)
   end
   
   return retm
end

function getDistance3d(m1,m2)
   local d1 = m1[1][1] - m2[1][1]
   local d2 = m1[2][1] - m2[2][1]
   local d3 = m1[3][1] - m2[3][1]

   return math.sqrt( (d1*d1) + (d2*d2) + (d3*d3) )
end

function nzt(val)
   if math.abs(val) < 1.0e9 then
      return false
   end
   return true
end

function matrixAreEqual(m1,m2)
   local ma = matrix.add(m1,m2,-1)
   --print(ma:pretty())
   if matrix.find(ma,nzt) then
      return false
   end
   return true
end

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
