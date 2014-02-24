function deconvolve_kernel, pcls, kernel=kernel, fkernel=fkernel, norm=norm, write_kern=write_kern

   True = 1b
   False = 0b

   if (not keyword_set(norm)) and keyword_set(kernel) then norm = 1./(4.*!pi)
   if keyword_set(kernel) and keyword_set(fkernel) then STOP, ' - insert EITHER kernel OR  kernelfile!'
   if keyword_set(fkernel) then begin
       kernel = mrdfits(fkernel, 0, hdr)
       klmax = n_elements( kernel[*,0] ) - 1
       if total(kernel[*,0]) eq 0. then begin
           print, ' - removing monopole block...'
           kernel = kernel[1:klmax,1:klmax]
       endif
       if total(kernel[*,0]) eq 0. then begin
           print, ' - removing dipole block...'
           kernel = kernel[1:klmax-1,1:klmax-1]
       endif

;##       kernel = kernel[2:klmax,2:klmax]
       help, kernel

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
   endif

;stop

   print, ' - kernel normalization = ', norm
   kernel = kernel * norm
   print, ' - inverting kernel...'
   km1 = invert( transpose( kernel ), err, /double )
   print, ' - inverting kernel...Done'
   print, ' - inversion status:', err

   cls = km1 ## pcls[2:*]

   ocls = pcls * 0.
   ocls[2:*] = cls
   return, ocls

end
