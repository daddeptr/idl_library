function ruler_read_maps_in, channels, $
                     Nside = Nside, $
                     verbosity=verbosity

; Ruler: linear solution in pixel space of the component separation
; problem
;
; D.Pietrobon Jul 2011
;
; Reading maps in
 
   False = 0b
   True  = 1b

   if (verbosity gt 0) then loud = True

   code = ' --> Ruler <--'
   temp_folder = '/global/scratch/sd/dpietrob/rubbish/'

   nfreq = n_elements(channels)
   
   Npix = 12l * Nside^2

; ----------------------------------------------------------------------
   Nfreq = n_elements(channels)
   freq = dblarr(Nfreq)
   a2t = dblarr(Nfreq)
   print, 'Frequency Channels:'
   print, ' Number of frequencies used: ', Nfreq

   tmpmap = { nametag:'', $
              cfreq: 0.d0, $
              map: dblarr(Npix), $
              rms: dblarr(Npix), $
              a2t: 0.d0, $
              monopole: 0.d0, $
              dipole: dblarr(3) }

   channel_maps = replicate(tmpmap, Nfreq)

   for ifreq=0,nfreq-1 do begin
       freq[ifreq] = channels[ifreq].cfreq
       a2t[ifreq] = conversionfactor(freq[ifreq], /antenna2thermo)

       channel_maps[ifreq].nametag = channels[ifreq].nametag
       channel_maps[ifreq].cfreq = channels[ifreq].cfreq
       channel_maps[ifreq].a2t = conversionfactor(freq[ifreq], /antenna2thermo)
       channel_maps[ifreq].monopole = channels[ifreq].monopole
       channel_maps[ifreq].dipole = channels[ifreq].dipole

       if (loud) then print, ifreq, '- ', freq[ifreq]
; ------ reading map
       if (loud) then print, 'reading ', channels[ifreq].map_filename 
       read_fits_map, channels[ifreq].map_filename, map, nside=mapnside, order=maporder
       if (mapnside ne Nside) then STOP, 'Nside mismatch: ', channels[ifreq].nametag, mapnside
       if (maporder eq 'NESTED') then begin
           if (loud) then print, 'reordering ', channels[ifreq].nametag
           map = reorder(map, /n2r)
       endif
       if (channel_maps[ifreq].monopole ne 0.d0) then map = map-channel_maps[ifreq].monopole 
       if ( (abs(channel_maps[ifreq].dipole[0]) gt 1.d-11) and (abs(channel_maps[ifreq].dipole[1]) gt 1.d-11) and (abs(channel_maps[ifreq].dipole[2]) gt 1.d-11) ) then begin
           dipole = make_dipole(Nside, channel_maps[ifreq].dipole)
           if (verbosity gt 2) then mollview, dipole, win=ifreq, tit=channels[ifreq].nametag
           map = map - dipole
           dipole = 0.
       endif
       channel_maps[ifreq].map = map
       map = 0.

; ------ reading rms
       if (loud) then print, 'reading ', channels[ifreq].rms_filename 
       read_fits_map, channels[ifreq].rms_filename, rms, nside=mapnside, order=maporder
       if (mapnside ne Nside) then STOP, 'Nside mismatch: ', channels[ifreq].nametag, mapnside
       if (maporder eq 'NESTED') then begin
           if (loud) then print, 'reordering ', channels[ifreq].nametag
           rms = reorder(rms, /n2r)
       endif
       channel_maps[ifreq].rms = rms
       rms = 0.
       
       if (loud) then print, freq[ifreq], a2t[ifreq]

   endfor ;frequency loop

return, channel_maps

end

