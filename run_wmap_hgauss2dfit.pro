pro run_wmap_hgauss2dfit, cat=cat, stp=stp, plt=plt, ff=ff, lf=lf

   true = 1b
   false = 0b

   if (not keyword_set(cat) ) then cat = 'pws'
   if (not keyword_set(stp) ) then stp = false
   if (not keyword_set(plt) ) then plt = 'no'
   if (not keyword_set(ff) )  then ff = 0
   if (not keyword_set(lf) ) then lf = 6

   sfreq=['K','Ka','Q','V','W']

   nfreq=n_elements(sfreq)

   fwhm = [0.88, 0.66, 0.51, 0.35, 0.22 ] * 60. ; in arcmin

   sigma0 = [1.456, 1.490, 2.197, 3.137, 6.549] ; in mK

   cal0 = [0.023113417d9, $
           0.056961569d9, $
           0.13454397d9, $
           0.241487d9, $
           0.369344d9, $
           0.481127d9, $
           0.288331d9, $
           0.0582363d9, $
           0.00223400d9 ]
 
   cal0 = cal0[*] * 0. + 1.

;   ff = 0
;   lf = 6

   do_pws = false
   do_smh = false

   if (cat eq 'pws') then do_pws = true
   if (cat eq 'smh') then do_smh = true
   if (cat eq 'wmap') then do_wmp = true

   do_all = true

   if (do_wmp and false) then begin
       cfile = '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/wmap/wmap_ptsrc_catalog_p4_7yr_v4.txt'
       readcol, cfile, ra, dec, glon, glat, flux1, flux2, flux3, flux4, flux5, fluxerr1, fluxerr2, fluxerr3, fluxerr4, fluxerr5, format='f, f, f, f, f, f, f, f, f, f, f, f, f, f'

       nsrc = n_elements(flux1)

       flux = fltarr(nfreq, nsrc)
       fluxerr = fltarr(nfreq, nsrc)
; ---
       flux[0,*] = flux1
       flux[1,*] = flux2
       flux[2,*] = flux3
       flux[3,*] = flux4
       flux[4,*] = flux5
; ---
       fluxerr[0,*] = fluxerr1
       fluxerr[1,*] = fluxerr2
       fluxerr[2,*] = fluxerr3
       fluxerr[3,*] = fluxerr4
       fluxerr[4,*] = fluxerr5

       for ifreq=0, nfreq-1 do begin
           gsrc = where( flux[ifreq,*] ne -9.9 )

           data = { RA : ra[gsrc], $
                    DEC : dec[gsrc], $
                    GLON : glon[gsrc], $
                    GLAT : glat[gsrc], $
                    FLUX : flux[ifreq, gsrc], $
                    FLUX_ERR : fluxerr[ifreq, gsrc], $
                    EXTENDED : ra[gsrc]*0 }

           mwrfits, data, 'wmap_ps_cat_'+sfreq[ifreq]+'.fits', /create
       endfor

   endif

   if (do_wmp and false) then begin
       for ifreq=0,nfreq-1 do begin
           read_fits_map, '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/wmap/combined_freq/raw_maps/wmap_band_iqumap_r9_7yr_'+sfreq[ifreq]+'_v4.fits', map
           help, map
           ud_grade, map, mapr, order_in='nest', order_out='ring'
           map = mapr
;           mollview, map[*,0]
           write_fits_map, 'wmap7_raw_band_'+sfreq[ifreq]+'_ring_uK.fits', map[*,0]*1.e3, /ring, units='!7l!8K CMB'
           write_fits_map, 'wmap7_raw_RMS_band_'+sfreq[ifreq]+'_ring_uK.fits', sqrt( (sigma0[ifreq]*1.e3)^2 / map[*,3] ), /ring, units='!7l!8K CMB'
           mollview, 'wmap7_raw_band_'+sfreq[ifreq]+'_ring_uK.fits', /hist, win=1
           mollview, 'wmap7_raw_RMS_band_'+sfreq[ifreq]+'_ring_uK.fits', win=2
       endfor
   endif
;stop

   for ifreq=ff,lf do begin
       print, sfreq[ifreq]
       if (do_pws) then cfile = '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/dx7/ns2048/dx7_catalogs/DX7_PwS_'+sfreq[ifreq]+'.fits'
       if (do_smh) then cfile = '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/dx7/ns2048/dx7_catalogs/DX7_MHWv2_'+sfreq[ifreq]+'.fits'
       if (do_wmp) then cfile = '/global/homes/d/dpietrob/myscratch/dx7_delivery_analysis/sources_Sergi/wmap_ps_cat_'+sfreq[ifreq]+'.fits'

       mfile = 'wmap7_raw_band_'+sfreq[ifreq]+'_ring_uK.fits'
       vfile = 'wmap7_raw_RMS_band_'+sfreq[ifreq]+'_ring_uK.fits'
       solution = dp_hgauss2dfit(robust=1b, map_infile=mfile, catalog_infile=cfile, varmap_infile=vfile, save=1b, r_fwhm=1.75, first_source=0, input_fwhm=fwhm[ifreq], map_calibration=cal0[ifreq], plot=plt, do_all_cases=do_all, do_stop=stp)
       
    
   endfor

end
