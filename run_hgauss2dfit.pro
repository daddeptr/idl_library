pro run_hgauss2dfit, cat=cat, stp=stp, plt=plt, ff=ff, lf=lf, tag=tag, fit_area_size=fit_area_size, fsrc=fsrc

   true = 1b
   false = 0b

   if (not keyword_set(cat) ) then cat = 'dx8'
   if (not keyword_set(stp) ) then stp = false
   if (not keyword_set(plt) ) then plt = 'no'
   if (not keyword_set(fsrc) )  then fsrc = 0
   if (not keyword_set(ff) )  then ff = 0
   if (not keyword_set(lf) ) then lf = 6
   if (not keyword_set(tag) ) then tag='test'
   if (not keyword_set(fit_area_size) ) then fit_area_size=2.

   sfreq=['030','044','070','100','143','217','353', '545', '857']

   nfreq=n_elements(sfreq)

   fwhm = [32.65, 27.00, 13.01, 9.94,   7.04,    4.66,    4.41,    4.47,    4.23]

   cal0 = [0.023113417d9, $
           0.056961569d9, $
           0.13454397d9, $
           0.241487d9, $
           0.369344d9, $
           0.481127d9, $
           0.288331d9, $
           0.0582363d9, $
           0.00223400d9 ]
 
   do_pws = false
   do_smh = false
   do_dx8 = false

   if (cat eq 'pws') then do_pws = true
   if (cat eq 'smh') then do_smh = true
   if (cat eq 'dx8') then do_dx8 = true

   do_all = true

   for ifreq=ff,lf do begin
       if (do_pws) then cfile = '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/dx7/ns2048/dx7_catalogs/DX7_PwS_'+sfreq[ifreq]+'.fits'
       if (do_smh) then cfile = '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/dx7/ns2048/dx7_catalogs/DX7_MHWv2_'+sfreq[ifreq]+'.fits'
       if (do_dx8) then cfile = '/global/homes/d/dpietrob/myscratch/dx8/input/catalogues/hfi/DX8_MHW_'+sfreq[ifreq]+'.fits'

;       mfile = '/global/homes/d/dpietrob/myscratch/dx7_delivery_analysis/input/maps/mondip_subtracted_map/dx7_Imap_'+sfreq[ifreq]+'_f.fits'
;       vfile = '/global/homes/d/dpietrob/myscratch/dx7_delivery_analysis/input/dx7_Irms'+sfreq[ifreq]+'GHz_ns2048_uK_f.fits'

;       mfile = '/global/homes/d/dpietrob/myscratch/dx7_delivery_analysis/input/maps/dx7_Imap'+sfreq[ifreq]+'GHz_ns2048_uK_f.fits'
;       vfile = '/global/homes/d/dpietrob/myscratch/dx7_delivery_analysis/input/maps/dx7_Irms'+sfreq[ifreq]+'GHz_ns2048_uK_f.fits'

;       mfile = '/global/homes/d/dpietrob/myscratch/dx7_delivery_analysis/input/dx7_Imap'+sfreq[ifreq]+'_f.fits'
;       vfile = '/global/homes/d/dpietrob/myscratch/dx7_delivery_analysis/input/dx7_Irms'+sfreq[ifreq]+'GHz_ns2048_uK_f.fits'

; --- dx8 
       mfile = '/global/homes/d/dpietrob/myscratch/dx8/maps/ns2048/dx8_Imap_'+sfreq[ifreq]+'GHz_ns2048_uK_fs.fits'
       vfile = '/global/homes/d/dpietrob/myscratch/dx8/maps/ns2048/dx8_Irms_'+sfreq[ifreq]+'GHz_ns2048_uK_fs.fits'

       solution = dp_hgauss2dfit(is_rms=true, robust=true, map_infile=mfile, $ 
                                 catalog_infile=cfile, $ ; coord=[[184.5],[-5.7]], $ ; 
                                 varmap_infile=vfile, save=1b, r_fwhm=fit_area_size, first_source=fsrc, input_fwhm=fwhm[ifreq], map_calibration=cal0[ifreq], plot=plt, do_all_cases=do_all, do_stop=stp, tag=tag+'_'+sfreq[ifreq])
       
    
   endfor

end
