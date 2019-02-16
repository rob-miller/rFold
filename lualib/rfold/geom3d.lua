--[[
   geom3d.lua
   
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

--- 3D geometry module : rotation,transformation matrix generation and manipulation
--
-- **external dependencies** : numlua (https://github.com/carvalho/numlua) -- needed to pull the lua_number2int fix for this repository
--
-- not a class but works better with ldoc if we say it is
-- @classmod geom3d


require 'numlua'


local geom3d = {}

--- initialises a 1x4 matrix
-- @return a 1x4 matrix with element [1][4]=1.0
function geom3d.get14mtx()
   local m = matrix.zeros(1,4)
   m[1][4]=1.0
   return m
end

--- initialises a 4x1 matrix
-- @return a 4x1 matrix with element [4][1]=1.0
function geom3d.get41mtx()
   local m = matrix.zeros(4,1)
   m[4][1]=1.0
   return m
end

--- initialises a 4x4 matrix
-- @return a 4x4 matrix with diagonal elements = 1.0
function geom3d.get44mtx()
   local m = matrix.zeros(4,4)
   m[1][1] = 1.0
   m[2][2] = 1.0
   m[3][3] = 1.0
   m[4][4] = 1.0
   return m
end

--- generate 3D transformation matrix
-- @param x offset in X dimension
-- @param y offset in Y dimension
-- @param z offset in Z dimension
-- @return a 4x4 3D transformation matrix
function geom3d.genMt(x,y,z)
   local mt = geom3d.get44mtx()
   --mt[1][1] = 1.0
   --mt[2][2] = 1.0
   --mt[3][3] = 1.0

   mt[1][4] = x
   mt[2][4] = y
   mt[3][4] = z
   return mt
end

--- generate 3D X axis rotation matrix
-- @param angr rotation angle in radians
-- @return a 4x4 3D rotation matrix
function geom3d.genMrx(angr)
   local mrx = geom3d.get44mtx()
   local cosang = math.cos(angr)
   local sinang = math.sin(angr)
   --mrx[1][1] = 1.0
   mrx[2][2] = cosang
   mrx[3][3] = cosang
   mrx[2][3] = -sinang
   mrx[3][2] = sinang
   return mrx
end

--- generate 3D Y axis rotation matrix
-- @param angr rotation angle in radians
-- @return a 4x4 3D rotation matrix
function geom3d.genMry(angr)
   local mry = geom3d.get44mtx()
   local cosang = math.cos(angr)
   local sinang = math.sin(angr)
   mry[1][1] = cosang
   --mry[2][2] = 1.0
   mry[3][3] = cosang
   mry[3][1] = -sinang
   mry[1][3] = sinang
   return mry
end

--- generate 3D Z axis rotation matrix
-- @param angr rotation angle in radians
-- @return a 4x4 3D rotation matrix
function geom3d.genMrz(angr)
   local mrz = geom3d.get44mtx()
   local cosang = math.cos(angr)
   local sinang = math.sin(angr)
   mrz[1][1] = cosang
   mrz[2][2] = cosang
   --mrz[3][3] = 1.0
   mrz[1][2] = -sinang
   mrz[2][1] = sinang
   return mrz
end

--- given m[ x y z ], return m[ theta phi distance 1 ] 
-- @param m nx1 (n>=3) matrix specifying X, Y, Z coordinates in first dimension
-- @return 4x1 matrix [1][1] = theta radians, [2][1] = phi radians, [3][1] = distance
function geom3d.getSphericalCoordinates(m)
   local retm = geom3d.get41mtx()
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

--- calculate distance between two points (1D vectors)
-- @param m1 nx1 matrix specifying coordinates for first point, e.g. [1][1] = x, [2][1] = y, [3][1] = z
-- @param m2 nx1 matrix specifying coordinates for second point
-- @return euclidean distance between m1 and m2
function geom3d.getDistance3d(m1,m2)
   local s1, s2 = matrix.size(m1,1), matrix.size(m2,1)
   assert(s1 == s2,'matrices not equal in dimension 1')
   local dsum = 0
   for i=1,s1 do
      local v = m1[i][1] - m2[i][1]
      dsum = dsum + (v * v)
   end
   return math.sqrt(dsum)
end

--- generate string with matrix in openSCAD format
-- @param mtx nxm matrix
-- @return string openSCAD representation of m
function geom3d.writeSCAD(mtx)
   local m = mtx:size(1)
   local n = mtx:size(2)

   local s = '[ '
   for i = 1, m do
      s = s .. '[ '
      for j = 1, n do
         s = s .. mtx[i][j] .. (j<n and ', ' or ' ')
      end
      s = s .. ']' .. (i<m and ', ' or ' ')
   end
   s = s .. ']'
   return s
end

local function nzt(val)
   if math.abs(val) > 1.0e-9 then
      return false
   end
   return true
end

--- calculate angle a1a2a3 from lengths of 3 sides of triangle
-- @param d0 length a1a2
-- @param d1 length a2a3
-- @param d2 length a1a3
-- @return angle a1a2a3 in radians
function geom3d.getAngleS3(d0,d1,d2)
   return math.acos( ((d0*d0) + (d1*d1) - (d2*d2)) / (2*d0*d1) )
end

--- test two matrices for equality
-- @param m1 matrix 1 to compare
-- @param m2 matrix 2 to compare
-- @return boolean true if equal
function geom3d.matrixAreEqual(m1,m2)
   local ma = matrix.add(m1,m2,-1)
   --print(ma:pretty())
   if matrix.find(ma,nzt) then
      return false
   end
   return true
end

--- generate transformation matrix to coordinate space defined by 3 points
-- @param a1 nx1 (n>=3) matrix specifying X,Y,Z coordinates of point on XZ plane orienting the new coordinate space
-- @param a2 nx1 (n>=3) matrix specifying X,Y,Z coordinates of point at origin of the new coordinate space
-- @param a3 nx1 (n>=3) matrix specifying X,Y,Z coordinates of point on Z axis orienting the new coordinate space
-- @param r boolean flag, if true return second matrix which reverses the transformation, bringing points in the new coordinate space back to the original coordinate space 
-- @return 4x4 3D transformation matrix placing a1 on XZ plane, a2 at 0,0,0, a3 at 0,0,+Z ; optionally (r flag) also return matrix implementing the reverse transformation
function geom3d.coordSpace(a1,a2,a3, r)
   -- return transform matrix 
   -- a's need to be 4x1 matrices
   -- if r then also generate reverse transformation matrix to bring back 

   local dbg=nil
   if dbg then
      print(a1:transpose():pretty())
      print(a2:transpose():pretty())
      print(a3:transpose():pretty())
   end
   
   local tm = geom3d.genMt(-a2[1][1], -a2[2][1], -a2[3][1])  -- translation matrix for a2 to origin
      
   -- now get a3 to Z-axis
   local p3 = geom3d.get41mtx()   -- direct translation of a3 using a2 to origin
   p3[1][1] = a3[1][1] - a2[1][1]
   p3[2][1] = a3[2][1] - a2[2][1]
   p3[3][1] = a3[3][1] - a2[3][1]

   local sc = geom3d.getSphericalCoordinates(p3)

   if dbg then
      print(sc:pretty())
   end
   
   local mrz = geom3d.genMrz(-sc[1][1])  -- rotate translated a3 -theta about Z
   local mry = geom3d.genMry(-sc[2][1])  -- rotate translated a3 -phi about Y
   
   local mt = mrz * tm
   mt = mry * mt                  -- mt completes a2-a3 on Z-axis, still need to align a1 with XZ plane

   if dbg then
      print((mt * a3):pretty())
   end
   
   p3 = mt * a1                   -- p3 not needed, re-use it here

   local sc2 = geom3d.getSphericalCoordinates(p3)    -- need theta of translated a1
   local mrz2 = geom3d.genMrz(-sc2[1][1])            -- rotate a1 -theta about Z to align with X

   mt = mrz2 * mt                 -- mt transforms to specified coordinate space

   if not r then return mt end

   -- generate the reverse transformation

   mrz2 = geom3d.genMrz(sc2[1][1])            -- rotate a1 theta about Z, reversing alignment with X
   mry  = geom3d.genMry(sc[2][1])             -- rotate a3 phi about Y
   mrz  = geom3d.genMrz(sc[1][1])             -- rotate a3 theta about Z
   tm   = geom3d.genMt(a2[1][1], a2[2][1], a2[3][1])  -- translation matrix for origin to a2

   local mr = mry * mrz2
   mr = mrz * mr
   mr = tm * mr

   return mt,mr
end

return geom3d
