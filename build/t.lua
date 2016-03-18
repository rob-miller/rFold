#!/usr/bin/env luajit

require 'rtmGeom3d'

local m1 = matrix.zeros(4,1)
local m2 = matrix.zeros(4,1)

m2[3][1]=1

if matrixAreEqual(m1,m2) then
   print('equal')
else
   print('not equal')
end
