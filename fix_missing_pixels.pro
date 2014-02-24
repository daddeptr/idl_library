True = 1b
False = 0b

dataset = 'dx8'
tags = ['cmb', 'co', 'dust', 'lowfreq', 'sample']

ntags = n_elements(tags)

Nside = 2048l
Npix = 12l * Nside^2
mask = fltarr(Npix)
mask[*] = 1.

for itag=0,ntags-1 do begin
   print, tags[itag]
   spawn, 'ls '+dataset+'/*'+tags[itag]+'*.fits', files
   nfiles = n_elements(files)
   for ifile=0,nfiles-1 do begin
      print, files[ifile]
      read_fits_map, files[ifile], map, nside=dpns, order=dpord
      if (dpns eq Nside) then begin
         bp = finite( map, /nan)
         ibp = where( bp eq True )
         help, ibp
         if (ibp[0] ne -1) then begin
            map[ibp] = 0.
            mask[ibp] = 0.	
            write_fits_map, files[ifile], map, order=dpord, units='!7l!8K'
            mollview, files[ifile], min=-300, max=300, chars=1.5, px=1000, win=1
         endif
      endif else begin
         print, 'Different Nside', dpns, ': ', files[ifile]
      endelse
   endfor
endfor

write_fits_map, 'dx8_commander-ruler_v4_mask.fits', mask, /ring

STOP

END
