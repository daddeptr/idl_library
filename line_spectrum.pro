function line_spectrum, nu, approx_sed, nu_ref=nu_ref

   if ( not keyword_set(nu_ref) ) then nu_ref = 101.2d9

   if (1b) then begin
       x = [ approx_sed[*,0], nu_ref ]
       y = [ approx_sed[*,1], 1.d0 ]
       ix = sort(x)
       
       x = x[ix]
       y = y[ix]
       sed = interpol(y,x,nu)
   endif else begin
       help, approx_sed

       x = [ approx_sed[*,0], nu_ref ]
       y = [ approx_sed[*,1], 1.d0 ]
       ix = sort(x)
       
       x = x[ix]
       y = y[ix]
       
       iref = where(x eq nu_ref)
       
       y2 = spl_init( x, y )
       
       sed = SPL_INTERP(X, Y, Y2, nu)
       
       sed = sed / sed[iref]
       print, iref
   endelse
;stop
   return, sed

END
