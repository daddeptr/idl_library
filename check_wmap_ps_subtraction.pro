   true = 1b
   false = 0b

   sfreq=['K','Ka','Q','V','W']

   nfreq=n_elements(sfreq)

   fwhm = [0.88, 0.66, 0.51, 0.35, 0.22] * 60.

   beam = 60.

   ns = 512l

   ff = 0
   lf = 4

   do_wmap = true
   do_pws = false
   do_smh = false

   do_smth = false
   do_stamps = false
   do_dg = true

   if (do_stamps) then begin
   for ifreq=ff,lf do begin
       print, sfreq[ifreq]

;;        spawn, 'mv wmap_p__ca__'+sfreq[ifreq]+'.txt wmap_ps_cat_'+sfreq[ifreq]+'.txt'
;;        spawn, 'mv wmap_p__ca__'+sfreq[ifreq]+'_ps_map.fits wmap_ps_cat_'+sfreq[ifreq]+'_ps_map.fits'
       map_infile='wmap7_raw_band_'+sfreq[ifreq]+'_ring_uK.fits'

       if (do_wmap) then begin
          catalog_infile = 'wmap_ps_cat_'+sfreq[ifreq]+'.fits'
          psmap_infile = 'wmap_ps_cat_'+sfreq[ifreq]+'_ps_map.fits'
          rootfile = 'wmap7_ps-removed'
       endif

       if (do_pws) then begin
           catalog_infile='/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/dx7/ns2048/dx7_catalogs/DX7_PwS_'+sfreq[ifreq]+'.fits'
           psmap_infile = 'DX7_PwS_'+sfreq[ifreq]+'_ps_map.fits'
           rootfile = 'dx7_PwS-ps_removed'
       endif
       
       if (do_smh) then begin
           catalog_infile='/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/dx7/ns2048/dx7_catalogs/DX7_MHWv2_'+sfreq[ifreq]+'.fits'
           psmap_infile = 'DX7_MHWv2_'+sfreq[ifreq]+'_ps_map.fits'
           rootfile = 'dx7_MHWv2-ps_removed'
       endif

       read_fits_map, map_infile, map
       read_fits_map, psmap_infile, psmap

       if (false) then begin
       if (false) then begin
           smica_file = 'smica_maps/subtracted_map_'+sfreq[ifreq]+'.fits'
           read_fits_map, smica_file, smica, nside=smicans, order=smicaord
           if ( (smicans ne 2048) and false) then begin
               ud_grade, smica, smicad, order_in=smicaord, order_out='ring', nside_out=2048
               smica = smicad
           endif
           smica = smica*1.e6
           zeros= where(smica eq 0.)
           remove_dipole, smica, nside=2048, ordering='ring', /onlymonopole, gal_cut=60
           smica[zeros] = 0.
           write_fits_map, 'smica_maps/my_subtracted_map_'+sfreq[ifreq]+'.fits', smica, /ring, units='!7l!8K CMB'
       endif else begin
           smica_file = 'smica_maps/my_subtracted_map_'+sfreq[ifreq]+'.fits'
           read_fits_map, smica_file, smica, nside=smicans, order=smicaord
           if (false) then ismoothing, smica_file, 'smica_maps/my_smooth_map_'+sfreq[ifreq]+'_ns2048.fits', fwhm_arcmin=sqrt(beam^2-fwhm[ifreq]^2)
        endelse
        endif

       if (do_smth and false) then ismoothing, map_infile, 'dx7_smooth_map_'+sfreq[ifreq]+'_ns2048.fits', fwhm_arcmin=sqrt(beam^2-fwhm[ifreq]^2)
       if (do_smth and true) then ismoothing, map_infile, 'wmap7_smooth_raw_band_'+sfreq[ifreq]+'.fits', fwhm_arcmin=sqrt(beam^2-fwhm[ifreq]^2)

       clean = map-psmap
       write_fits_map, rootfile+'_map_'+sfreq[ifreq]+'.fits', clean, /ring, units='!7l!8K CMB'

;;        remove_dipole, clean, nside=ns, ordering='ring', /onlymonopole, gal_cut=50

       if (do_smth) then ismoothing, rootfile+'_map_'+sfreq[ifreq]+'.fits', rootfile+'_smooth_map_'+sfreq[ifreq]+'.fits', fwhm_arcmin=sqrt(beam^2-fwhm[ifreq]^2)

       read_fits_map, rootfile+'_smooth_map_'+sfreq[ifreq]+'.fits', smthmap

       read_fits_s, catalog_infile, hdr, table

       theta   = table.GLAT
       phi     = table.GLON

       flux    = table.FLUX
       fluxerr = table.FLUX_ERR

       s2n = flux / fluxerr
    
       gsource = where( s2n ge 2.5)
       isource = reverse(sort(flux[gsource]))
       gsource = gsource[isource]

       theta = theta[gsource]
       phi = phi[gsource]

       nsource = n_elements(gsource)

       frfile = ''
       for isource=0,9 do begin
;           gnomview, 'wmap7_smooth_raw_band_'+sfreq[ifreq]+'.fits', chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Smoothed In Map', win=0, grat=[1,1], min=-400, max=400
;           gnomview, map, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Input Map', win=1, grat=[1,1], min=-400, max=400
;           gnomview, psmap, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='PS Map', win=3, grat=[1,1], min=-400, max=400
;           gnomview, clean, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='PS removed Map', win=2, grat=[1,1], min=-400, max=400
;           gnomview, smthmap, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Smoothed Out Map', win=4, grat=[1,1], min=-400, max=400
; ---
           if ((isource lt 20) and true) then begin
               gnomview, 'wmap7_smooth_raw_band_'+sfreq[ifreq]+'.fits', chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Smoothed In Map', grat=[1,1], win=-1, png='pics/src_'+sfreq[ifreq]+'_insmth_'+string(isource,format='(i2.2)')+'.png', min=-400, max=400
               gnomview, map, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Input Map', grat=[1,1], win=-1, png='pics/src_'+sfreq[ifreq]+'_in_'+string(isource,format='(i2.2)')+'.png', min=-400, max=400
               gnomview, psmap, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='PS Map', grat=[1,1], win=-1, png='pics/src_'+sfreq[ifreq]+'_ps_'+string(isource,format='(i2.2)')+'.png', min=-400, max=400
               gnomview, clean, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='PS removed Map', grat=[1,1], win=-1, png='pics/src_'+sfreq[ifreq]+'_out_'+string(isource,format='(i2.2)')+'.png', min=-400, max=400
               gnomview, rootfile+'_smooth_map_'+sfreq[ifreq]+'.fits', chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Smoothed Map', grat=[1,1], win=-1, png='pics/src_'+sfreq[ifreq]+'_smth_'+string(isource,format='(i2.2)')+'.png', min=-400, max=400
               
               file = 'pics/src_'+sfreq[ifreq]+'_insmth_'+string(isource,format='(i2.2)')+'.png ' + $
                 ' pics/src_'+sfreq[ifreq]+'_in_'+string(isource,format='(i2.2)')+'.png' + $
                 ' pics/src_'+sfreq[ifreq]+'_ps_'+string(isource,format='(i2.2)')+'.png' + $
                 ' pics/src_'+sfreq[ifreq]+'_out_'+string(isource,format='(i2.2)')+'.png' + $
                 ' pics/src_'+sfreq[ifreq]+'_smth_'+string(isource,format='(i2.2)')+'.png'

               spawn, 'montage '+file+ ' -adjoin -tile 5x -geometry 200x200 -background none pics/src_'+sfreq[ifreq]+'_strip_'+string(isource,format='(i2.2)')+'.png'

               frfile = frfile + 'pics/src_'+sfreq[ifreq]+'_strip_'+string(isource,format='(i2.2)')+'.png '
;               spawn, 'display pics/src_'+sfreq[ifreq]+'_strip_'+string(isource,format='(i2.2)')+'.png &'
           endif
;stop
       endfor
       spawn, 'montage ' + frfile + ' -adjoin -tile x10 -geometry +1+1 -background none pics/src_'+sfreq[ifreq]+'_1-10.png'

    endfor
    endif

   if (do_dg and false) then for ifreq=ff,lf do ud_grade,'wmap7_ps-removed_smooth_map_'+sfreq[ifreq]+'.fits', 'wmap7_ps-removed_smooth_map_'+sfreq[ifreq]+'_ns128.fits', nside_out=128

   if (do_dg and true) then for ifreq=ff,lf do ud_grade,'wmap7_smooth_raw_band_'+sfreq[ifreq]+'.fits', 'wmap7_smooth_map_'+sfreq[ifreq]+'_ns128.fits', nside_out=128


end
