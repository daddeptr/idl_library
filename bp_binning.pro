function bp_binning, cls, bin_file, verbose=verbose

   nl = n_elements(cls)
   readcol, bin_file, fl, ll
   nb = n_elements(ll)

   if keyword_set(verbose) then print, fl[nb-1], ll[nb-1]

   lmax = min( [nl-1, ll[nb-1]] )
   if keyword_set(verbose) then print, lmax

   bok = where(ll ge lmax)
   bmax = nb-1
   if ( bok[0] gt 0 ) then bmax = bok[0]

;## bin_cls = fltarr(nb)
   bin_cls = fltarr( bmax+1 )
;## for ib=0,bmax do bin_cls[ib] = total( cls[fl[ib]:ll[ib]] ) / (ll[ib]-fl[ib]+1)
   for ib=0,bmax do begin
       if ll[ib] le lmax then begin
           if keyword_set(verbose) then print, fl[ib], ll[ib]
           ndl = ll[ib]-fl[ib]+1
           bin_cls[ib] = total( cls[fl[ib]:ll[ib]] ) / ndl
       endif else begin
           if keyword_set(verbose) then print, fl[ib], lmax
           ndl = lmax-fl[ib]+1
           bin_cls[ib] = total( cls[fl[ib]:*] ) / ndl
       endelse
       if keyword_set(verbose) then print, bin_cls[ib]
   endfor
   return, bin_cls

end
