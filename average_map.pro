function average_map, root, std=std

   spawn, 'ls '+root+'*.fits', files
   nfiles = n_elements(files)

   print, 'Reading files: ', nfiles
   read_fits_map, files[0], m, nside=dpns, order=dpord
   avg = dblarr(12l*dpns^2)
   if ( keyword_set(std) ) then avg2 = avg
   for ifile=0, nfiles-1 do begin
       read_fits_map, files[ifile], m, nside=rns, order=rord
       if (rns ne dpns) then stop, ' Different Nside', files[i]
       if (rord ne dpord) then stop, ' Different Ordering', files[i]
       avg = avg + m
       if ( keyword_set(std) ) then avg2 = avg2 + m^2
   endfor

   avg = avg / nfiles
   if ( keyword_set(std) ) then begin
       avg2 = avg2 / nfiles
       std = sqrt(avg2-avg^2)
   endif

return, avg
end
