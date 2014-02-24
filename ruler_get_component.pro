   True = 1b
   False = 0b

   sfreq = ['030','044','070','100','143','217','353']
   Nfreq = n_elements(sfreq)

   comp = ['cmb', 'co', 'thermal_dust', 'ff']
   Ncomp = 1

   root = 'chains/glss/pix_7b/'

   Nside = 2048l
   Npix = nside2npix(Nside)

   for icomp=0,Ncomp-1 do begin
       solution = fltarr(Npix)
       for ifreq=0, Nfreq-1 do begin
           read_fits_map, '../ns2048/dx9_Imap_'+sfreq[ifreq]+'_ns2048_uK.fits', map
           read_fits_map, root+'avrg_weight_'+comp[icomp]+'_'+sfreq[ifreq]+'_ns2048.fits', weight
           solution = solution + map*weight
       endfor
       mollview, solution, min=-300, max=300
   endfor

   stop

end
