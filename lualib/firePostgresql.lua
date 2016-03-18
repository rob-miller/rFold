#!/usr/bin/env luajit

function pgQ(sql)
   return assert(pgcon:execute(sql),sql):fetch({})
end

function pgQa(sql)
   return assert(pgcon:execute(sql)):fetch({},'a')
end

function pgQcur(sql)
   return assert(pgcon:execute(sql))
end

--[[
local iterate = function(cur)
   return function() return cur:fetch({},'a') end
end

local cur = pgQcur("select target_id, text from pdbenviron.annotation where text like 'c2cStructMap:%'")

for smRow in iterate(cur) do
]]


function pgas2tx (astr)  -- postgresql array str to table
   print('pgas2t enter: ', astr)
   astr = string.match(astr,'{(.+)}')
   --print('pgas2t remove brackets: ', astr)
   local rslt = {}
   for t in string.gmatch(astr,"[^,]+") do
      print ('t= ' .. t )
      rslt[#rslt +1] = t
   end
   --print('pgas2t output: ', (rslt))
   return rslt
end

function pgas2t (astr)  -- postgresql array str to table
   local rslt = {}
   for t in string.gmatch(string.match(astr,'{(.+)}'),"[^,]+") do rslt[#rslt +1] = t end
   return rslt
end

function pgas2tn (astr)  -- postgresql array str to table of numbers
   local rslt = {}
   for t in string.gmatch(string.match(astr,'{(.+)}'),"[^,]+") do rslt[#rslt +1] = tonumber(t) end
   --[[
   io.write('pgas2tn: [' .. #rslt ..'] { ')
   for t=1,#rslt do
      io.write( rslt[t] .. ' ')
   end
   print('}')
   ]]
   return rslt
end


function dbReconnect(autocommit)
   pg = luasql.postgres()
   pgcon = assert(pg:connect("rFold","postgres","postgres","localhost"))
   if autocommit then
      pgcon:setautocommit(true)
   else 
      pgcon:setautocommit(false)
   end
end

function dbClose()
   if pgcon then pgcon:close() end
end

luasql = require "luasql.postgres"
--pg = assert (luasql.postgres())

dbReconnect()


--pgcon = pg:connect("heptares","postgres","postgres","localhost")
