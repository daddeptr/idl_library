   function make_dipole, nside, direction=direction, silent=silent, theta=theta, phi=phi, amplitude=amplitude

      if (keyword_set(direction) and (keyword_set(theta) or keyword_set(phi) or (keyword_set(amplitude))) ) then stop, 'Inconsistent inputs' 
      if ( not keyword_set(direction) and (keyword_set(theta) and keyword_set(phi) and (keyword_set(amplitude))) ) then begin
          theta = (90.-theta)/!radeg
          phi = phi / !radeg
          ang2vec, theta, phi, direction
          direction = direction * amplitude
      endif

      if not keyword_set(silent) then print, direction
 

      npix = 12l*nside^2

      dipole = fltarr(npix)
      
;      for ipix=0l,npix-1 do begin
         pix2vec_ring, nside, lindgen(npix), vec
;help, direction ## vec
;stop
         dipole = direction[0] * vec[*,0] + direction[1] * vec[*,1] + direction[2] * vec[*,2]
;      endfor
      
      return, dipole

   end
