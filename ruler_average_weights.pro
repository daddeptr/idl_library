   True = 1b
   False = 0b

   sfreq = ['030','044','070','100','143','217','353']
   Nfreq = n_elements(sfreq)

   comp = ['cmb', 'co', 'thermal_dust', 'ff']
   Ncomp = 1

   root = 'chains/glss/pix_7b/'

   Nside = 128l
   Npix = nside2npix(Nside)

   for icomp=0,Ncomp-1 do begin
       for ifreq=0, Nfreq-1 do begin
           spawn, 'ls '+root+'*weights*'+comp[icomp]+'*_'+sfreq[ifreq]+'.fits', files
;           print, files
           Nfiles = n_elements(files)
           print, Nfiles
;stop
           ave_weight = dblarr(Npix)
           for ifile=0,Nfiles-1 do begin
               read_fits_map, files[ifile], map
               ave_weight = ave_weight + map
           endfor
           ave_weight = ave_weight / Nfiles
           print, 'Saving '+root+'avrg_weight_'+comp[icomp]+'_'+sfreq[ifreq]+'.fits...'
           write_fits_map, root+'avrg_weight_'+comp[icomp]+'_'+sfreq[ifreq]+'.fits', ave_weight, /ring, units='[dimensionless]'
           print, 'Upgrading...'
           hres_weight = upgrade_index_ring(ave_weight, 2048l)
           write_fits_map, root+'avrg_weight_'+comp[icomp]+'_'+sfreq[ifreq]+'_ns2048.fits', hres_weight, /ring, units='[dimensionless]'
           mollview, root+'avrg_weight_'+comp[icomp]+'_'+sfreq[ifreq]+'_ns2048.fits', px=500
;stop
       endfor
   endfor

   stop

end
