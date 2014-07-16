pro oploterr, x, y=y, er=er, low=low, up=up, line=line, col=col

;   print, up
;   print, low
   n = n_elements(x)
;   print, n
   if not defined(line) then line=0
   if not defined(col) then col=0
   if defined(y) and n_elements(x) ne n_elements(y) then stop, "X and Y don't have the same length"

   if (not defined(er)) and (not defined(low)) and (not defined(up)) then stop, 'Input error'

   if (not defined(low)) then low = y-er
   if (not defined(up)) then up = y+er

   for i=0,n-1 do oplot, [x[i],x[i]],[low[i],up[i]], col=col, line=line
end

