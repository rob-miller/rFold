#!/usr/bin/env luajit

require 'rtmlualib'
require 'parseCmdLine'
require 'parseRtmDssp'

require 'numlua'

function genMrx(angr)
   local mrx = matrix.zeros(3,3)
   local cosang = math.cos(angr)
   local sinang = math.sin(angr)
   mrx[1][1] = 1.0
   mrx[2][2] = cosang
   mrx[3][3] = cosang
   mrx[2][3] = sinang
   mrx[3][2] = -sinang
   return mrx
end

function genMry(angr)
   local mry = matrix.zeros(3,3)
   local cosang = math.cos(angr)
   local sinang = math.sin(angr)
   mry[1][1] = cosang
   mry[2][2] = 1.0
   mry[3][3] = cosang
   mry[3][1] = sinang
   mry[1][3] = -sinang
   return mry
end

function genMrz(angr)
   local mrz = matrix.zeros(3,3)
   local cosang = math.cos(angr)
   local sinang = math.sin(angr)
   mrz[1][1] = cosang
   mrz[2][2] = cosang
   mrz[3][3] = 1.0
   mrz[1][2] = sinang
   mrz[2][1] = -sinang
   return mrz
end



local m0 = matrix.zeros(1,3)
m0[1][1] = 1.0

print(m0:pretty())

local ang = 90.0 * d2r
--[[
local mry = genMry(ang)
m0 = m0 * mry
print(m0:pretty())

local mrx = genMrx(ang)
m0 = m0 * mrx
print(m0:pretty())
--]]

local mrz = genMrz(ang)
m0 = m0 * mrz
print(m0:pretty())
