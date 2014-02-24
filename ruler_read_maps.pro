function ruler_read_maps, channels, $
                          Nside = Nside, $
                          verbosity=verbosity

; Ruler: linear solution in pixel space of the component separation
; problem
;
; Reading maps in
 
   False = 0b
   True  = 1b

   if (not keyword_set(verbosity)) then verbosity = 0
   if (verbosity gt 0) then loud = True

   code = ' --> Ruler <--'

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

       if (loud) then print, ifreq, '- ', freq[ifreq]
; ------ reading map
       if (loud) then print, 'reading ', channels[ifreq].map_filename 
       read_fits_map, channels[ifreq].map_filename, map, nside=mapnside, order=maporder
       if (mapnside ne Nside) then STOP, 'Nside mismatch: ', channels[ifreq].nametag, mapnside
       if (maporder eq 'NESTED') then begin
           if (loud) then print, 'reordering ', channels[ifreq].nametag
           map = reorder(map, /n2r)
       endif

       channel_maps[ifreq].map = map

; ------ reading rms
       if (loud) then print, 'reading ', channels[ifreq].rms_filename 
       read_fits_map, channels[ifreq].rms_filename, rms, nside=mapnside, order=maporder
       if (mapnside ne Nside) then STOP, 'Nside mismatch: ', channels[ifreq].nametag, mapnside
       if (maporder eq 'NESTED') then begin
           if (loud) then print, 'reordering ', channels[ifreq].nametag
           rms = reorder(rms, /n2r)
       endif
       channel_maps[ifreq].rms = rms
       
       if (loud) then print, freq[ifreq], a2t[ifreq]

   endfor ;frequency loop

return, channel_maps

end
