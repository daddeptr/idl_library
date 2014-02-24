   function make_sky_cut, angle_deg, nside

      nside = long(nside)
      print, angle_deg, nside 
      npix = 12l*nside^2

      mask = fltarr(npix)
      ipix = lindgen(npix)

      pix2ang_ring, nside, ipix, theta, phi	

      theta_deg = 90. - theta * !radeg	

      bp = where(abs(theta_deg) lt angle_deg )


      mask[*] = 1.
      mask[bp] = 0.	

      return ,mask

   end
