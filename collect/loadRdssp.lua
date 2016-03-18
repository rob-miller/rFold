#!/usr/bin/env luajit

require 'rtmlualib'
require 'parseCmdLine'
require 'firePostgresql'

require 'parseRtmDssp'

dbReconnect(1)

local update=true


function rdssp2db ( rdssp )
   for i,j in pairs(rdssp) do
      if ('table' == type(j)) then
         print(i,'table: ' .. table.getn(j) .. ' elements')
      else
         print(i, j)
      end
   end
end




local args = parseCmdLine();
if (0 < table.getn(args)) then
   for i,j in ipairs(args) do
      rdssp2db(prd(args[i]))
   end
else
   rdssp2db(prd())
end
