#!/usr/bin/env luajit

function parseCmdLine (valid_flags, valid_params)
   local args = {}

   for i=1,#arg do
      if string.match(arg[i],'^-') then
         if string.match(arg[i],'=') then
            local key, val
            key, val = string.match(arg[i],'^-(%S+)=(%S+)$')
            assert(valid_params[ key ], 'parameter -' .. key .. ' value ' .. val .. ' not recognised')
            args[ key ] = val
         else
            flag = string.match(arg[i],'^-(%S+)$')
            assert(valid_flags[ flag ], 'flag -' .. flag .. ' not recognised')
            args[flag] = (args[flag] or 0)+1
         end
      else
         args[#args +1] = arg[i]
      end
   end

   return args
end

-- ./parseCmdLine.lua -p1=hello -f1 -f1 fee fie fo fum

--[[
local vflags = { ['f1'] = true, ['f2'] = true }
local vparams = { ['p1'] = true, ['p2'] = true }

local args = parseCmdLine(vflags, vparams)

for i,j in pairs(args) do
   print('param/flag/arg ', i, j)
end

--]]

   
