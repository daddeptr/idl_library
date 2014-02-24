; 
; Smoothing to lower resolution to derive weights
;
pro change_map_resolution, filename, dNside=dNside, smthv=smthv, map_beamfile=map_beamfile, map_reso=map_reso, tonly=tonly, root=root, silent=silent, no_highwin=no_highwin, bw_lim=bw_lim, $
                           beam_file=beam_file
   True = 1b
   False = 0b

   same_beam = False

   spawn, 'rm tmpalm*.fits tmpbeam*.fits'
   if keyword_set(silent) then silent=1 else silent=False

   if keyword_set(smthv) then begin
       if ( (n_elements(smthv) eq 1) and (not keyword_set(tonly)) ) then begin
           smthv = reform( [smthv, smthv, smthv] )
           same_beam = True
       endif
       gb = fltarr(4*dNside+1,3) + 1.
       for i=0,n_elements(smthv)-1 do begin
           print, ' - beam value = ', smthv[i]
           if (smthv[i] gt 0.) then gb[*,i] = gaussbeam(smthv[i], 4*dNside)
       endfor
   endif

   cosine_apo = fltarr(4*dNside+1)+1.
   cosine_apo[0:2*dNside] = 1.
   indx = lindgen(4*dNside+1)
   cosine_apo[2*dNside+1:3*dNside] = (1.+cos(!pi*(indx[2*dNside+1:3*dNside]-2*dNside)/dNside) )/2.
   cosine_apo[3*dNside+1:*] = 0.
   
;plot, cosine_apo
;stop

   mapfile = filename
   if not keyword_set(root) then maproot = filename else maproot = root
   if not keyword_set(dNside) then stop, ' - dNside not provided'
   if ( (not keyword_set(smthv)) and (not keyword_set(beam_file)) ) then stop, ' - smthv/beam not provided'
   print, ' - mapfile = ', mapfile
   print, ' - dNside  = ', dNside
   if keyword_set(smthv) then print, ' - smthv   = ', smthv
   if keyword_set(beam_file) then begin
       print, ' - beam_file   = ', beam_file
       fits2bl, gb, beam_file
   endif
   print, ' - gb      = '
   help, gb
   read_fits_map, filename, m, nside=Nside

   if (not keyword_set(map_reso)) then begin
       if (not keyword_set(tonly)) then readcol, map_beamfile, l, wl143i, wl143q, wl143u, format='f,f,f,f' else readcol, map_beamfile, l, wl143i, format='f,f'
   endif

   if keyword_set(map_reso) then begin
       wl143i = gaussbeam(map_reso[0],4*dNside)
       if not keyword_set(tonly) then begin
           wl143q = gaussbeam(map_reso[1],4*dNside)
           wl143u = gaussbeam(map_reso[2],4*dNside)
       endif
   endif 

   wl143i = double( wl143i[where(wl143i gt 0.)] )
   if not keyword_set(tonly) then wl143q = double( wl143q[where(wl143q gt 0.)] )
   if not keyword_set(tonly) then wl143u = double( wl143u[where(wl143u gt 0.)] )

   hwf = healpixwindow( Nside, 3 )
   hwf = double( hwf )

   dhwf = healpixwindow( dNside, 3 )
   dhwf = double( dhwf )

   lmax = min( [n_elements( wl143i )-1, n_elements(dhwf[*,0])-1, 4*dNside ] ) 
   print, ' - lmax = ', lmax

   if not keyword_set(tonly) then begin
       ewl1 = fltarr( lmax+1, 3 ) + 1.
       if keyword_set(no_highwin) then begin
           ewl1[2:lmax,0:2] = [ [ gb[2:lmax,0] / wl143i[2:lmax] * dhwf[2:lmax,0] ], $
                                [ gb[2:lmax,1] / wl143q[2:lmax] * dhwf[2:lmax,1] ], $
                                [ gb[2:lmax,2] / wl143u[2:lmax] * dhwf[2:lmax,1] ] ]
       endif else begin
           ewl1[2:lmax,0:2] = [ [ gb[2:lmax,0] / wl143i[2:lmax] / hwf[2:lmax,0] * dhwf[2:lmax,0] ], $
                                [ gb[2:lmax,1] / wl143q[2:lmax] / hwf[2:lmax,1] * dhwf[2:lmax,1] ], $
                                [ gb[2:lmax,2] / wl143u[2:lmax] / hwf[2:lmax,1] * dhwf[2:lmax,1] ] ]
       endelse
       if keyword_set(bw_lim) then ewl1[0:lmax,0:2] = ewl1[0:lmax,0:2] * [ [cosine_apo[0:lmax]], [cosine_apo[0:lmax]], [cosine_apo[0:lmax]] ]
   endif else begin
       ewl1 = fltarr( lmax+1 ) + 1.
       ewl1[2:lmax] = [ gb[2:lmax] / wl143i[2:lmax] / hwf[2:lmax,0] * dhwf[2:lmax,0] ]
       if keyword_set(bw_lim) then ewl1[0:lmax] = ewl1[0:lmax] * cosine_apo[0:lmax]
   endelse       
   ewl2 = ewl1
       
   bl2fits, ewl1, 'tmpbeam.fits'

   print, ' - anafasting...'
   if not keyword_set(tonly) then begin
       ianafast, mapfile, cls1, simul_type=2, alm1_out='tmpalms1.fits', nlmax=lmax, silent=silent, /double
       fits2alm, index, alms1, 'tmpalms1.fits', 'ALL'
   endif else begin
       ianafast, mapfile, cls1, simul_type=1, alm1_out='tmpalms1.fits', nlmax=lmax, silent=silent, /double
       fits2alm, index, alms1, 'tmpalms1.fits'
   endelse

   print, ' - smoothing alms...'
   talms1 = alms1 * 0.
   tcls1 = cls1 * 0.

   for i=0l, lmax do begin
       li=i^2 + i + 1
       lf=i^2 + 2*i + 1
       ili=where(index eq li)
       ilf=where(index eq lf)
       
       tcls1[i,0] = cls1[i,0] * ewl1[i,0]^2
       if not keyword_set(tonly) then begin
           tcls1[i,1] = cls1[i,1] * ewl1[i,1]^2
           tcls1[i,2] = cls1[i,2] * ewl1[i,2]^2
           tcls1[i,3] = cls1[i,3] * ewl1[i,0] * ewl1[i,1]
           tcls1[i,4] = cls1[i,4] * ewl1[i,0] * ewl1[i,2]
           tcls1[i,5] = cls1[i,5] * ewl1[i,1] * ewl1[i,2]
       endif
       
       talms1[ili:ilf,0,0] = alms1[ili:ilf,0,0] * ewl1[i,0]
       talms1[ili:ilf,1,0] = alms1[ili:ilf,1,0] * ewl1[i,0]
       if not keyword_set(tonly) then begin
           talms1[ili:ilf,0,1] = alms1[ili:ilf,0,1] * ewl1[i,1]
           talms1[ili:ilf,1,1] = alms1[ili:ilf,1,1] * ewl1[i,1]
           talms1[ili:ilf,0,2] = alms1[ili:ilf,0,2] * ewl1[i,2]
           talms1[ili:ilf,1,2] = alms1[ili:ilf,1,2] * ewl1[i,2]
       endif
   endfor
   print, ' - smoothing alms. Done'
   
   alm2fits, index, talms1, 'tmpalms1.fits'

   if not keyword_set(tonly) then isynfast, tcls1, maproot+'.ns'+string(dNside,format='(i4.4)')+'.b'+string(reform(smthv[0]),format='(f5.1)'), $
     alm_in='tmpalms1.fits', apply_window=0, simul_type=2, nside=dNside, nlmax=lmax, /double, silent=silent else $
     isynfast, tcls1, maproot+'.ns'+string(dNside,format='(i4.4)') + '.b'+string(reform(smthv[0]),format='(f5.1)'), alm_in='tmpalms1.fits', apply_window=0, simul_type=1, nside=dNside, nlmax=lmax, /double, silent=silent

   mollview, maproot+'.ns'+string(dNside,format='(i4.4)') + '.b'+string(reform(smthv[0]),format='(f5.1)'), /asinh, win=11, px=650

   if same_beam then smthv = reform( smthv[0] )
   print, ' >>> change_map_resolution: --- End of Program ---'

end
