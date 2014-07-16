function xf_binning, cls, bin_file, verbose=verbose, lcen=lcen, errors=errors

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
   lcen = fltarr( bmax+1 )
   if keyword_set(errors) then cls = cls[*]^2
   bin_cls = fltarr( bmax+1 )
;## for ib=0,bmax do bin_cls[ib] = total( cls[fl[ib]:ll[ib]] ) / (ll[ib]-fl[ib]+1)
   for ib=0,bmax do begin
       if ll[ib] le lmax then begin
           if keyword_set(verbose) then print, fl[ib], ll[ib]
           ndl = ll[ib]-fl[ib]+1
           wlol = fltarr(ndl) + 1.
           wl = findgen(ndl) + fl[ib]
           wl = (wl+0.5)/(wl+1.)
           norm = total( wl * wlol)
           wlol /= norm
           bin_cls[ib] = total( cls[fl[ib]:ll[ib]] * wlol * (findgen(ndl)+fl[ib]) * (findgen(ndl)+fl[ib]+0.5)) / 2. / !pi
           lcen[ib] = (fl[ib]+ll[ib])/2.
       endif else begin
           if keyword_set(verbose) then print, fl[ib], lmax
           ndl = lmax-fl[ib]+1
           wlol = fltarr(ndl) + 1.
           wl = findgen(ndl) + fl[ib]
           wl = (wl+0.5)/(wl+1.)
           norm = total( wl * wlol)
           wlol /= norm
           bin_cls[ib] = total( cls[fl[ib]:*] * wlol * (findgen(ndl)+fl[ib]) * (findgen(ndl)+fl[ib]+0.5) ) /2. / !pi
           lcen[ib] = (fl[ib]+lmax)/2.
       endelse
       if keyword_set(verbose) then print, bin_cls[ib]
   endfor

   if keyword_set(errors) then bin_cls = sqrt(bin_cls[*])
   return, bin_cls

end
