#!/usr/bin/env luajit

function filterPdbLine(line)
      if line:match('^HEADER ') then  
      elseif line:match('^TITLE ') then  -- dssp does not pass but will grab from pdb
      elseif line:match('^TER') then
      elseif line:match('^END') then
      elseif line:match('^ATOM') then
         line = line:sub(1,56) .. '1.00  0.00           ' .. line:sub(-3,-2) .. ' '   -- replace occupancy, temperature factor, segment identifier and charge on atom
      else
         return
      end

      io.write(line .. '\n')
end

if #arg > 0 then
   for i=1,#arg do
      io.input(arg[i])
      for line in io.lines() do filterPdbLine(line) end
   end
else
   for line in io.lines() do filterPdbLine(line) end
end
