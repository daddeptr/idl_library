; 
; Smoothing to lower resolution to derive weights
;
   True = 1b
   False = 0b

   isynfast, 'data/camb_97283428_scalcls_uK.fits', 'cmb_ns2048.fits', nside=2048, nlmax=4096, fwhm_arcmin=7.03, simul_type=2, iseed=-10
   mapfile = 'cmb_ns2048.fits'

   smthv   = 15.
   Nside  = 2048
   dNside = 1024

   gb = gaussbeam(smthv, 5000l)
   
   wl143i = gaussbeam(7.03, 4000)
   wl143q = gaussbeam(7.03, 4000)
   wl143u = gaussbeam(7.03, 4000)

   hwf = healpixwindow( Nside,3 )
   hwf = double( hwf )

   dhwf = healpixwindow( dNside,3 )
   dhwf = double( dhwf )

   lmax = min( [n_elements( wl143i )-1, n_elements(dhwf[*,0])-1] ) 
   print, ' - lmax = ', lmax

   ewl1 = fltarr( lmax+1, 3 ) + 1.
   ewl1[2:lmax,0:2] = [ [ gb[2:lmax] / wl143i[2:lmax] / hwf[2:lmax,0] * dhwf[2:lmax,0] ], $
                        [ gb[2:lmax] / wl143q[2:lmax] / hwf[2:lmax,1] * dhwf[2:lmax,1] ], $
                        [ gb[2:lmax] / wl143u[2:lmax] / hwf[2:lmax,1] * dhwf[2:lmax,1] ] ]

   loadct, 39
   window, 10 & plot, ewl1[*,0], /ylog, chars=1.5, yr=[1.e-16,1]
   oplot, gb
   oplot, gaussbeam(sqrt(15.^2-7.03^2), 5000), col=70
   oplot, wl143i, line=2, col=245
   oplot, wl143i*ewl1, line=2, col=205

   print, ' - anafasting...'
   ianafast, mapfile, cls1, simul_type=2, alm1_out='tmpalms1.fits', nlmax=lmax, /silent, /double
   fits2alm, index, alms1, 'tmpalms1.fits', 'ALL'

;stop
   print, ' - smoothing alms...'
   talms1 = alms1 * 0.
   tcls1 = cls1 * 0.
   
   for i=0l, lmax do begin
       li=i^2 + i + 1
       lf=i^2 + 2*i + 1
       ili=where(index eq li)
       ilf=where(index eq lf)
       
       tcls1[i,0] = cls1[i,0] * ewl1[i,0]^2
       tcls1[i,1] = cls1[i,1] * ewl1[i,1]^2
       tcls1[i,2] = cls1[i,2] * ewl1[i,2]^2
       tcls1[i,3] = cls1[i,3] * ewl1[i,0] * ewl1[i,1]
       tcls1[i,4] = cls1[i,4] * ewl1[i,0] * ewl1[i,2]
       tcls1[i,5] = cls1[i,5] * ewl1[i,1] * ewl1[i,2]
       
       talms1[ili:ilf,0,0] = alms1[ili:ilf,0,0] * ewl1[i,0]
       talms1[ili:ilf,1,0] = alms1[ili:ilf,1,0] * ewl1[i,0]
       talms1[ili:ilf,0,1] = alms1[ili:ilf,0,1] * ewl1[i,1]
       talms1[ili:ilf,1,1] = alms1[ili:ilf,1,1] * ewl1[i,1]
       talms1[ili:ilf,0,2] = alms1[ili:ilf,0,2] * ewl1[i,2]
       talms1[ili:ilf,1,2] = alms1[ili:ilf,1,2] * ewl1[i,2]
   endfor
   print, ' - smoothing alms. Done'
    
;print, total(alms1)   
;stop
   alm2fits, index, talms1, 'tmpalms1.fits'
   isynfast, tcls1, mapfile+'.b'+string(smth,format='(f4.1)')+'.v2', alm_in='tmpalms1.fits', apply_window=0, simul_type=2, nside=dNside, nlmax=lmax, /double

   mollview, mapfile + '.b'+string(smth,format='(f4.1)')+'.v2', /asinh, win=11, px=650

   stop, ' --- End of Program ---'
end
