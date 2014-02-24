   true = 1b
   false = 0b

   sfreq=['030','044','070','100','143','217','353', '545', '857']

   nfreq=n_elements(sfreq)

   fwhm = [32.65, 27.00, 13.01, 9.94,   7.04,    4.66,    4.41,    4.47,    4.23]

   beam = 40.

   cal0 = [0.023113417d9, $
           0.056961569d9, $
           0.13454397d9, $
           0.241487d9, $
           0.369344d9, $
           0.481127d9, $
           0.288331d9, $
           0.0582363d9, $
           0.00223400d9 ]
 
   ff = 0
   lf = 0

   do_pws = false
   do_smh = false
   do_dx8 = true

   do_smth = false
   do_plt = true

   for ifreq=ff,lf do begin
       print, sfreq[ifreq]

;;        map_infile='/global/homes/d/dpietrob/myscratch/dx7_delivery_analysis/input/maps/mondip_subtracted_map/dx7_Imap_'+sfreq[ifreq]+'_f.fits'
;;        map_infile='/global/homes/d/dpietrob/myscratch/dx7_delivery_analysis/input/maps/dx7_Imap'+sfreq[ifreq]+'GHz_ns2048_uK_f.fits'
       map_infile='/global/scratch/sd/dpietrob/dx8/maps/ns2048/dx8_Imap_'+sfreq[ifreq]+'GHz_ns2048_uK_fs.fits'

       if (do_dx8) then begin
           catalog_infile='/global/scratch/sd/dpietrob/dx8/input/catalogues/hfi/DX8_MHW_'+sfreq[ifreq]+'.fits'
           psmap_infile = 'dx8_smhw_'+sfreq[ifreq]+'_ps_map.fits'
           rootfile = 'dx8_ps_removed'
           fits2cl, chisq_vec, 'dx8_smhw_v3_'+sfreq[ifreq]+'_chisq_vec.fits'
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
           if (false) then begin
               smica_file = 'smica_maps/my_subtracted_map_'+sfreq[ifreq]+'.fits'
               read_fits_map, smica_file, smica, nside=smicans, order=smicaord
               if (false) then ismoothing, smica_file, 'smica_maps/my_smooth_map_'+sfreq[ifreq]+'_ns2048.fits', fwhm_arcmin=sqrt(beam^2-fwhm[ifreq]^2)
           endif
       endelse

       if (do_smth) then ismoothing, map_infile, 'dx7_smooth_map_'+sfreq[ifreq]+'_ns2048.fits', fwhm_arcmin=sqrt(beam^2-fwhm[ifreq]^2)

       clean = map-psmap
       write_fits_map, rootfile+'_map_'+sfreq[ifreq]+'_ns2048.fits', clean, /ring, units='!7l!8K CMB'

       remove_dipole, clean, nside=2048, ordering='ring', /onlymonopole, gal_cut=50

       if (do_smth) then ismoothing, rootfile+'_map_'+sfreq[ifreq]+'_ns2048.fits', rootfile+'_smooth_map_'+sfreq[ifreq]+'_ns2048.fits', fwhm_arcmin=sqrt(beam^2-fwhm[ifreq]^2)

       read_fits_map, rootfile+'_smooth_map_'+sfreq[ifreq]+'_ns2048.fits', smthmap

       if (false) then begin
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
       endif else begin

           all_source = n_elements( where(chisq_vec[*,2] gt 0.) )

           igood = where( (chisq_vec[*,2] gt 0.) and (chisq_vec[*,2] lt 1.1) and (chisq_vec[*,3] eq 0) )
           phi = chisq_vec[igood,0]
           theta = chisq_vec[igood,1] 
           nsource = n_elements(igood)

           print, ' Good fit: ', nsource, '/', all_source

       endelse
       frfile = ''
       if (do_plt) then begin
           fsrc = 0
           lsrc = 20

           for isource=fsrc,lsrc do begin
               print, ' Source characterization:', isource, '/', nsource
               print, chisq_vec[igood[isource], *]

               thetar = (90.-theta[isource])/!radeg
               phir = phi[isource]/!radeg
               ang2pix_ring, 2048, thetar,phir,ipix

;           gnomview, 'dx7_smooth_map_'+sfreq[ifreq]+'_ns2048.fits', chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Smoothed In Map', win=0, grat=[1,1], min=-400, max=400
               cc = map
               cc[ipix] = -1.6375e30
               gnomview, cc, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Input Map '+sfreq[ifreq], win=1, grat=[1,1], /log, reso=0.75, min=10, max=1.e5

               cc = psmap
               cc[ipix] = -1.6375e30
               gnomview, cc, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='PS Map '+sfreq[ifreq], win=3, grat=[1,1], /log, reso=0.75, min=10, max=1.e5

               cc = clean
               cc[ipix] = -1.6375e30
               gnomview, cc, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='PS removed Map '+sfreq[ifreq], win=2, grat=[1,1], /log, reso=0.75, min=10, max=1.e5
;           gnomview, smica, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Smica PS removed Map', win=5, grat=[1,1], min=-400, max=400
;           gnomview, smthmap, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Smoothed Out Map', win=4, grat=[1,1], min=-400, max=400
; ---
               if (isource lt 20) then begin
                   gnomview, 'dx7_smooth_map_'+sfreq[ifreq]+'_ns2048.fits', chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Smoothed In Map', grat=[1,1], win=-1, png='pics/src_'+sfreq[ifreq]+'_insmth_'+string(isource,format='(i2.2)')+'.png', /log, min=10, max=1.e5
                   gnomview, map, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Input Map', grat=[1,1], win=-1, png='pics/src_'+sfreq[ifreq]+'_in_'+string(isource,format='(i2.2)')+'.png', /log, min=10, max=1.e5
                   gnomview, psmap, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='PS Map', grat=[1,1], win=-1, png='pics/src_'+sfreq[ifreq]+'_ps_'+string(isource,format='(i2.2)')+'.png', /log, min=10, max=1.e5
                   gnomview, clean, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='PS removed Map', grat=[1,1], win=-1, png='pics/src_'+sfreq[ifreq]+'_out_'+string(isource,format='(i2.2)')+'.png', /log, min=10, max=1.e5
                   gnomview, smica, chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Smica PS removed Map', grat=[1,1], win=-1, png='pics/src_'+sfreq[ifreq]+'_smica_'+string(isource,format='(i2.2)')+'.png', /log, min=10, max=1.e5
                   gnomview, rootfile+'_smooth_map_'+sfreq[ifreq]+'_ns2048.fits', chars=1.5, rot=[phi[isource], theta[isource], 0], tit='Smoothed Map', grat=[1,1], win=-1, png='pics/src_'+sfreq[ifreq]+'_smth_'+string(isource,format='(i2.2)')+'.png', /log, min=10, max=1.e5
               
                   file = 'pics/src_'+sfreq[ifreq]+'_insmth_'+string(isource,format='(i2.2)')+'.png ' + $
                     ' pics/src_'+sfreq[ifreq]+'_in_'+string(isource,format='(i2.2)')+'.png' + $
                     ' pics/src_'+sfreq[ifreq]+'_ps_'+string(isource,format='(i2.2)')+'.png' + $
                     ' pics/src_'+sfreq[ifreq]+'_out_'+string(isource,format='(i2.2)')+'.png' + $
                     ' pics/src_'+sfreq[ifreq]+'_smth_'+string(isource,format='(i2.2)')+'.png'

;                     ' pics/src_'+sfreq[ifreq]+'_smica_'+string(isource,format='(i2.2)')+'.png' + $
;                     ' pics/src_'+sfreq[ifreq]+'_smth_'+string(isource,format='(i2.2)')+'.png'
;                   spawn, 'montage '+file+ ' -adjoin -tile 6x -geometry 250x250 -background none pics/src_'+sfreq[ifreq]+'_strip_'+string(isource,format='(i2.2)')+'.png'

                   spawn, 'montage '+file+ ' -adjoin -tile 5x -geometry 250x250 -background none pics/src_'+sfreq[ifreq]+'_strip_'+string(isource,format='(i2.2)')+'.png'
                   

                   frfile = frfile + 'pics/src_'+sfreq[ifreq]+'_strip_'+string(isource,format='(i2.2)')+'.png '
;               spawn, 'display pics/src_'+sfreq[ifreq]+'_strip_'+string(isource,format='(i2.2)')+'.png &'
               endif

;stop
           endfor
           spawn, 'montage ' + frfile + ' -adjoin -tile x10 -geometry +1+1 -background none pics/src_'+sfreq[ifreq]+'_1-10.png'
           
       endif

   endfor


end
