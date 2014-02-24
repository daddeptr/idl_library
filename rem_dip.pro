   pro rem_dip, gsig=gsig, ns=ns, smoo=smoo

      freq=['030', '044', '070', '100', '143', '217','353']

      nfreq = n_elements(freq)

      ss = strtrim(string(long(gsig)),2)
      s_ns = strtrim(string(ns),2)
      ssmoo = strtrim(string(long(smoo)),2)

      start_freq = 0 
      for ifreq = start_freq, nfreq-1 do begin

         read_fits_map,'ns'+s_ns+'/madam_map_'+freq[ifreq]+'GHz_TQU_smth'+ssmoo+'_uK'+ss+'_ns'+s_ns+'.fits',map
;         mollview, map
         temp = map[*,0]
         remove_dipole, temp, ordering='ring', nside=ns, gal_cut=30., /onlymonopole
;         mollview, temp
         write_fits_map, 'ns'+s_ns+'/noDip_madam_map_'+freq[ifreq]+'GHz_I_smth'+ssmoo+'_uK'+ss+'_ns'+s_ns+'.fits', temp, /ring, units='!7l!17K Thermodynamic'
         mollview, 'ns'+s_ns+'/noDip_madam_map_'+freq[ifreq]+'GHz_I_smth'+ssmoo+'_uK'+ss+'_ns'+s_ns+'.fits'
;stop
      endfor

   end
