#!/usr/bin/env luajit

-- extend string.match to test and capture patterns
-- https://inspired-lua.org/index.php/2013/05/extend-string-match-to-test-and-capture-patterns/

do
    local smatch = string.match     -- keep the original definition of string.match
 
    -- String matching function
    -- Same results as string:match but as a side effect
    -- places the captures in global variables _1, _2, ...
    function string:match(pat)
        local matches = {smatch(self, pat)}    -- call the original match to do the work
        for i = 1, #matches do                 -- #matches == 0 if no matches
            _G["_" .. i] = matches[i]          -- assign captures to global variables
        end
        return unpack(matches)                 -- return original results
    end
end

