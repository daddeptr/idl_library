   function make_disc, radius, nside, direction=direction, radians=radians

      if (not keyword_set(direction)) then direction = [1,0,0]

      nside = long(nside)
      print, radius, nside 
      npix = 12l*nside^2

      mask = fltarr(npix)
;##      ipix = lindgen(npix)

      if (not keyword_set(radians)) then radius = radius / !radeg

      query_disc, nside, direction, radius, listpix

;##      pix2ang_ring, nside, ipix, theta, phi	

;##      theta_deg = 90. - theta * !radeg	

;##      bp = where(abs(theta_deg) lt angle_deg )


      mask[listpix] = 1.
;##      mask[bp] = 0.	

      return ,mask

   end
