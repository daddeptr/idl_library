function make_ps_map, catalog_file, nside_out

   !path=!path+':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'
   !path=!path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'

   true = 1b
   false = 0b

   npix = 12l*nside_out^2
   psmap = fltarr(npix)

   print, ' - Reading '+catalog_file
   table = mrdfits(catalog_file,1)

   flux = table.FLUX
;   ra = table.RA
;   dec = table.DEC
   glat = table.GLAT
   glon = table.GLON
;   extended = table.EXTENDED

;   s2n = table.FLUX/table.FLUX_ERR

   theta = glat
   phi = glon
;   print, 'theta ',minmax(theta)
;   print, 'phi ',minmax(phi)

   phi = phi/!radeg
   theta = !pi/2 - theta/!radeg
      
;   print, 'theta ',minmax(theta)
;   print, 'phi ',minmax(phi)
;stop
;##   ang2vec, theta, phi, vec
   ang2pix_ring, nside_out, theta, phi, ipix

   psmap[ipix] = psmap[ipix] + flux

   return, psmap

end
