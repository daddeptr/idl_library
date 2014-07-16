function deconvolve_kernel, pcls, kernel=kernel, fkernel=fkernel, norm=norm, write_kern=write_kern, inv_fkernel=inv_fkernel, idipmon_rem=idipmon_rem, write_cls=write_cls

   True = 1b
   False = 0b

   ncl = n_elements( pcls[0,*] )
   print, ' - ncl = ', ncl

   if (not keyword_set(norm)) and (keyword_set(kernel) or keyword_set(fkernel))then norm = 1./(4.*!pi)

   if keyword_set(kernel) and keyword_set(fkernel) then STOP, ' - insert EITHER kernel OR  kernelfile!'

   if keyword_set(fkernel) then begin
       kernel = mrdfits(fkernel, 0, hdr)
       klmax = n_elements( kernel[*,0] ) - 1

       nfield = n_elements( kernel[0,*] ) / klmax

       dipmon_rem = 0
       kern = kernel[0:klmax,0:klmax]
       pkern = kernel[0:klmax,(klmax+1):(klmax+1)+klmax]
       mkern = kernel[0:klmax,2*(klmax+1):2*(klmax+1)+klmax]
       xkern = kernel[0:klmax,3*(klmax+1):3*(klmax+1)+klmax]

       if total(kernel[*,0]) eq 0. then begin
           print, ' - removing monopole block...'
           kern = kern[1:klmax,1:klmax]
           pkern = pkern[1:klmax,1:klmax]
           mkern = mkern[1:klmax,1:klmax]
           xkern = xkern[1:klmax,1:klmax]
           dipmon_rem = 1
       endif
       if total(kernel[*,0]) eq 0. then begin
           print, ' - removing dipole block...'
           kern = kern[1:klmax-1,1:klmax-1]
           pkern = pkern[1:klmax-1,1:klmax-1]
           mkern = mkern[1:klmax-1,1:klmax-1]
           xkern = xkern[1:klmax-1,1:klmax-1]
           dipmon_rem = 2
       endif

       help, kern
;stop
;##       if not keyword_set(inv_fkernel) then begin

           print, ' - kernel normalization = ', norm
           kern = kern * norm
           if ncl gt 1 then begin
               pkern = pkern * norm * 0.5
               mkern = mkern * norm * 0.5
               xkern = xkern * norm * 0.5
           endif

           print, ' - inverting kernel...'
           ttkern = kern
           if ncl gt 1 then begin
;           eekern = 
           endif
           km1 = invert( transpose( ttkern ), err, /double )
           print, ' - inverting kernel...Done'
           print, ' - inversion status:', err
           
           if keyword_set(write_kern) then begin
               print, ' - saving inverse kernel...'+write_kern
               hdr = [hdr[0:n_elements(hdr)-2],'DIPMON_REM = '+strtrim(string(dipmon_rem),2),'END']
;##           print, hdr
               mwrfits, km1, write_kern, hdr
               print, ' - saving inverse kernel...Done'
           endif
;##       endif
       endif

   if keyword_set(inv_fkernel) then begin
       print, ' - reading in inv-kernel:'+inv_fkernel
       km1 = mrdfits(inv_fkernel, 0)
;##       print, hdr[n_elements(hdr)-2]
;##       tmp = strsplit(hdr[n_elements(hdr)-1],'=',/extract)
       if not keyword_set(idipmon_rem) then dipmon_rem = 2 else dipmon_rem = idipmon_rem
       print, ' - dipmon_rem = ',dipmon_rem
       klmax = n_elements( km1[*,0] ) - 1 + dipmon_rem
   endif

   lmax = min([n_elements(pcls)-1,klmax])
   print, ' - lmax = ', lmax
;##   cls = km1[2:lmax,2:lmax] ## pcls[2:lmax]
   cls = km1[0:lmax-dipmon_rem,0:lmax-dipmon_rem] ## pcls[dipmon_rem:lmax]

   ocls = fltarr(lmax+1) ;pcls * 0.
   ocls[2:lmax] = cls
   if keyword_set(write_cls) then begin
       print, ' - saving cls...'+write_cls
       cl2fits, ocls, write_cls
       print, ' - saving cls...Done'
   endif
   return, ocls

       if False then begin 
           fsky_i = intarr(n_elements(hdr))
           w2_i = intarr(n_elements(hdr))
           totwind_i = intarr(n_elements(hdr))

           for i=0,n_elements(hdr)-1 do fsky_i[i] = strcmp(hdr[i],'FSKY',4)
           fsky_i = where(fsky_i eq True)
           
           for i=0,n_elements(hdr)-1 do w2_i[i] = strcmp(hdr[i],'W2',2)
           w2_i = where(w2_i eq True)

           for i=0,n_elements(hdr)-1 do totwind_i[i] = strcmp(hdr[i],'TOTWIND',7)
           totwind_i = where(totwind_i eq True)
       
           fsky = strsplit( hdr[fsky_i[0]], '=', /extract)
;       print, fsky
           fsky = strsplit( fsky[1], '/', /extract)
;       print, fsky
           fsky = double( fsky[0] )
           print, ' - fsky = ', fsky

           w2 = strsplit( hdr[w2_i[0]], '=', /extract)
;       print, w2
           w2 = strsplit( w2[1], '/', /extract)
;       print, w2
           w2 = double( fsky[0] )
           print, ' - w2 = ', w2
       
           totwind = strsplit( hdr[totwind_i[0]], '=', /extract)
;       print, totwind
           totwind = strsplit( totwind[1], '/', /extract)
;       print, totwind
           totwind = double( totwind[0] )
           print, ' - totwind = ', totwind
       
           norm = fsky * w2 / totwind
       endif
       norm = 1.d0/ (4.d0*!dpi)
;stop

end
