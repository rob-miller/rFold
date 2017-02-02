--[[
   rfPostgresql.lua
   
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

--- Postgresql database access functions:

luasql = require "luasql.postgres"

local rfpg = {}

-- modify this line to change database access parameters:
rfpg.db, rfpg.user, rfpg.pass, rfpg.host = 'rFold','postgres','postgres','localhost'

rfpg.dbg=nil

function rfpg.Q(sql)
   if (rfpg.dbg) then print(sql) end
   return assert(rfpg.con:execute(sql),sql):fetch({})
end

function rfpg.Qa(sql)
   if (rfpg.dbg) then print(sql) end
   return assert(rfpg.con:execute(sql)):fetch({},'a')
end

function rfpg.Qcur(sql)
   if (rfpg.dbg) then print(sql) end
   return assert(rfpg.con:execute(sql))
end

--[[
local iterate = function(cur)
   return function() return cur:fetch({},'a') end
end

local cur = rfpg.Qcur("select target_id, text from pdbenviron.annotation where text like 'c2cStructMap:%'")

for smRow in iterate(cur) do
]]


function rfpg.as2tx (astr)  -- postgresql array str to table
   print('rfpg.as2t enter: ', astr)
   astr = string.match(astr,'{(.+)}')
   --print('rfpg.as2t remove brackets: ', astr)
   local rslt = {}
   for t in string.gmatch(astr,"[^,]+") do
      print ('t= ' .. t )
      rslt[#rslt +1] = t
   end
   --print('rfpg.as2t output: ', (rslt))
   return rslt
end

function rfpg.as2t (astr)  -- postgresql array str to table
   local rslt = {}
   for t in string.gmatch(string.match(astr,'{(.+)}'),"[^,]+") do rslt[#rslt +1] = t end
   return rslt
end

function rfpg.as2tn (astr)  -- postgresql array str to table of numbers
   local rslt = {}
   for t in string.gmatch(string.match(astr,'{(.+)}'),"[^,]+") do rslt[#rslt +1] = tonumber(t) end
   --[[
   io.write('rfpg.as2tn: [' .. #rslt ..'] { ')
   for t=1,#rslt do
      io.write( rslt[t] .. ' ')
   end
   print('}')
   ]]
   return rslt
end

function rfpg.pgArrayOfStrings(...) -- return variable number of strings as postgres array syntax
   local rstr = '{'
   local start=1
   for i,v in ipairs({...}) do
      rstr = rstr .. (start and '"' or ',"') .. v .. '"'
      start = nil
   end
   return rstr ..'}'
end

   





---------------------------------

function rfpg.dbReconnect(autocommit)
   rfpg.pg = luasql.postgres()
   rfpg.con = assert(rfpg.pg:connect(rfpg.db,rfpg.user,rfpg.pass,rfpg.host))
   if autocommit then
      rfpg.con:setautocommit(true)
   else 
      rfpg.con:setautocommit(false)
   end
end

function rfpg.dbClose()
   if rfpg.con then rfpg.con:close() end
end

rfpg.dbReconnect()

return rfpg

