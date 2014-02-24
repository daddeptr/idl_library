function upgrade_map, map_in, ordering_in, nside_up

   True  = 1b
   False = 0b

   sz = size(map_in)
   npix_in = sz[1]
   nside_in = npix2nside(npix_in)

   ud_grade, map_in, map_up, nside_out=nside_up, order_in=ordering_in, order_out='ring'

   write_fits_map, 'tmp_map_up.fits', map_up, /ring
   ianafast, map_up, 'tmp_ind_up_cls.fits', alm1_out='tmp_ind_up_alms.fits', simul_type=1, iter=1, nlmax=4*nside_in, /silent, /ring
   isynfast, 'tmp_ind_up_cls.fits', map, nside=nside_up, alm_in='tmp_ind_up_alms.fits', nlmax=4*nside_in, simul_type=1, /silent
   spawn, 'rm tmp_ind_up_cls.fits tmp_ind_up_alms.fits tmp_map_up.fits'

   if ( (ordering_in eq 'nest') or (ordering_in eq 'NEST') ) then map_out = reorder(map, in='ring', out=ordering_in) else map_out=map

   return, map_out

end

