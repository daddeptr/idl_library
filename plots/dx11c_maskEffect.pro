   True = 1b
   False = 0b

;## ------
   A = FIndGen(16) * (!PI*2/16.)
   UserSym, cos(A), sin(A), /fill
;## ---
   freq = ['143', 'recal_217']
   nfreq = n_elements(freq)
   mask = ['30','40','50','60','70','80','90']
   nmask = n_elements(mask)

   !p.multi=[0,3,1]
   window, 0, xsize=1300, ysize=600
   file = xfdir + 'dx11c_SDR_yr_1-2_IQUmap_recal_217_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const_dl50_fisherWins_xfcl_l3000_combined_30_x_combined_30_l4000.newdat'
   s = extract_xfaster_newdat( file, lcen=lc)
   specs = fltarr(nfreq, nmask, n_elements(s))
   for ifreq=0,nfreq-1 do begin
       file = xfdir + 'dx11c_SDR_yr_1-2_IQUmap_'+freq[ifreq]+'_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const_dl50_fisherWins_xfcl_l3000_combined_30_x_combined_30_l4000.newdat'
       s = extract_xfaster_newdat( file, lcen=lc)
       specs[ifreq,0,*] = s
       plot, lc, s/s, chars=3, xtit='!6!8l!6', ytit='!6!8D!dl!u!6fsky!n!8/D!dl!u!630', yr=[0.95,1.05], xr=[0,2500], tit=freq[ifreq]+' spectra ratio'
       for imask=1,nmask-1 do begin
           file = xfdir + 'dx11c_SDR_yr_1-2_IQUmap_'+freq[ifreq]+'_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const_dl50_fisherWins_xfcl_l3000_combined_'+mask[imask]+'_x_combined_'+mask[imask]+'_l4000.newdat'
           sf=extract_xfaster_newdat( file, lcen=lc)
           specs[ifreq,imask,*] = sf
           oplot, lc, sf/s, col=40*imask, psym=-8
       endfor
       legend, mask[1:*]+'/30', psym=lonarr(nmask-1)+8, col=(lindgen(nmask-1)+1)*40, /bottom
   endfor
   
   plot, lc, s/s, chars=3, xtit='!6!8l!6', ytit='!6!8D!dl!u!6fsky,143!n!8/D!dl!u!6fsky,217', yr=[0.975,1.025], xr=[0,2500], tit='143-to-217 spectrum ratio'
   for imask=0,nmask-1 do oplot, lc, specs[0,imask,*]/specs[1,imask,*], col=40*imask, psym=-8
   legend, 'f!dsky!n='+mask, psym=lonarr(nmask)+8, col=lindgen(nmask)*40
   write_png, 'dx11c_maskEffect.png', tvrd(/true)
end
