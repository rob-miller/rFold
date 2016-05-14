#!/usr/bin/env luajit

function usage(valid_flags, valid_params)
   print(arg[0] .. ' command line error')
   if valid_flags then
      io.write(' valid flags: ')
      for i,f in ipairs(valid_flags) do io.write( '-' .. f .. ' ') end
      print()
   end
   if valid_params then
      io.write(' valid parameters: ')
      for i,p in ipairs(valid_params) do io.write( '-' .. p .. '=  ') end
      print()
   end
end
   
function parseCmdLine (valid_flags, valid_params)
   local args = {}

   for i=1,#arg do
      if string.match(arg[i],'^-') then
         if string.match(arg[i],'=') then
            local key, val
            key, val = string.match(arg[i],'^-(%S+)=(%S+)$')
            if not (valid_params and valid_params[ key ]) then
               usage(valid_flags, valid_params)
               assert(nil, 'parameter -' .. key .. ' value ' .. val .. ' not recognised')
            end
            args[ key ] = val
         else
            flag = string.match(arg[i],'^-(%S+)$')
            if not (valid_flags and valid_flags[ flag ]) then
               usage(valid_flags, valid_params)
               assert(nil, 'flag -' .. flag .. ' not recognised')
            end
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

   
