function upgrade_index_ring, lowres_index, highres_nside, fwhm

   if not keyword_set(fwhm) then fwhm = 0.

   true = 1b
   false = 0b

   trash = '/global/scratch/sd/dpietrob/rubbish/'

   low_Nside = npix2nside(n_elements(lowres_index[*,0]))

   ud_grade, lowres_index, udmap, nside_out=highres_Nside, order_in='ring'
;##   ianafast, udmap, cls, alm1_out=root+'tmp_ind_up_alms.fits', simul_type=1, nlmax=4*low_Nside, /silent, /ring
;##   isynfast, cls, synmap, nside=highres_Nside, alm_in=root+'tmp_ind_up_alms.fits', nlmax=4*low_Nside, simul_type=1, /silent
; Smart Loris'

   if (fwhm gt 0.) then gb = gaussbeam(fwhm, 2*highres_nside) else gb=fltarr(2*highres_nside+1)+1.
   hwl = healpixwindow(low_Nside)
   hwh = healpixwindow(highres_Nside)
   bl = hwh / hwl * gb
   bl2fits, bl, trash + 'tmp_beam.fits'
   ismoothing, udmap, synmap, fwhm_arcmin=0., simul_type=1, nlmax=4*low_Nside, /silent, /ring, beam_file=trash + 'tmp_beam.fits'

return, synmap

end
