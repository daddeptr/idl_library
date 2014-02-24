!path=!path+'/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'
!path=!path+'/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'

true = 1b
false = 0b

freq = [30,44,70,100,143,217,353]
sfreq= string(freq,format='(i3.3)')
nfreq = n_elements(freq)

Nside = [1024l,1024l,1024l,2048l,2048l,2048l,2048l]
Npix = 12l*Nside^2 

beam = [33., 23., 12.7143, 9.64461, 7.11419, 4.72230, 4.5648]
; to avoid ringing
;beam = [33., 23., 12.7143, 10., 10., 10., 10.]

root = '/project/projectdirs/planck/user/graca/CSI/ffp3V4PwS201/'

;mNside = 2048l
;mNpix = 12l*mNside^2

;mask_tot = fltarr(mNpix)
;mask_tot[*] = 1.

ff=0
lf=6

if (false) then begin
for ifreq=ff,lf do begin
   print, freq[ifreq]

   if (false) then begin
      cf = conversionfactor(freq[ifreq],/brightness2muKthermo)
      print, root + 'pws201_'+sfreq[ifreq]+'_2_cat.csv'
      readcol, root + 'pws201_'+sfreq[ifreq]+'_2_cat.csv', theta, phi, flux, flux_err, stn

      print, ' conversion: ', cf

     help, flux
      gsource = where(stn gt 4.)
      help, gsource
   
      ns = Nside[ifreq]

      npx = 12l * ns^2
      ipix = lindgen(npx)
      pix2ang_ring, ns, ipix, ptheta, pphi
      pix2vec_ring, ns, ipix, pvec

;   ang2pix_ring, mNside, theta, phi, ipix
;mask_tot[ipix] = 0.

;mollview, mask_tot
;mask = fltarr(Npix[ifreq])
;mask[*] = 1.
;mask[ipix] = 0. 

      ang2vec, theta, phi, vec

;   ang2pix_ring, ns, theta, phi, ipring
;   np = n_elements(vec[*,0])
      ngsource = n_elements(gsource)

;;    for ip=0,np-1 do begin
;;       query_disc, mnside, vec[ip,0:*], 3.*beam[ifreq]/60./!radeg,listpix
;;       mask_tot[listpix] = 0.
;;    endfor

      mask = fltarr(Npix[ifreq])
      parea = 4.*!pi / Npix[ifreq]

      beam_rad = beam[ifreq]/2. / 60. * !pi/180.
 
      barea = (beam[ifreq]/2./60.*!pi/180.)^2

      radius = 5. * !pi/180

      pix2vec_ring, ns, 0, v0
      cdot = reform( transpose(pvec[*,*]) ## reform(v0) ) 
;      thetas = acos(cdot)
      query_disc, ns, v0, radius, listpix
;      g = where(cdot gt 0.)
      expo = cdot * 0.
      expo[listpix] = exp( -(1.-cdot[listpix]^2)/(2.*beam_rad^2) ) ;/ (2.*!pi*beam_rad^2)

      norm = total(expo)
      print, norm
      print, 1./norm

      expo = expo / norm

;      m = cdot*0.
;      m[0] = 1.
;      ismoothing, m, ms, fwhm_arcmin=beam[ifreq], /ring, /silent

;      gnomview, expo/ms, rot=[0,90,0], reso=.5, min=-50, max=50
      gnomview, expo, rot=[0,90,0], reso=.5                                                                                                                
;      print, total(ms)
;stop
      for ip=0,ngsource-1 do begin
;      help, pvec
;      help, vec[gsource[ip]]
         query_disc, ns, vec[gsource[ip],*], radius, listpix
         cdot = transpose(pvec[listpix,*]) ## reform(transpose(vec[gsource[ip],*]))
;      help, cdot
         cdot = reform(cdot)
;      help, cdot
;      cdot[where(cdot lt 0.)] = 0
;         temp=mask
;         temp[listpix]=cdot
;         mollview, temp
;      gnomview, cdot, rot=vec[gsource[ip],*]
;         thetas = acos(cdot)
;      x = 
;      gpix = where()
;      dtheta = theta[gsource[ip]]-ptheta
;      mollview, dtheta
;      mask = mask + flux[gsource[ip]] * cf * 1.e-9/parea * exp(-thetas^2/(2.*beam_rad^2))/(2.*!pi*beam_rad^2)
         mask[listpix] = mask[listpix] + flux[gsource[ip]] * cf * 1.e-9/parea * exp(-(1.-cdot^2)/(2.*beam_rad^2)) / norm ;total(exp(-(1.-cdot^2)/(2.*beam_rad^2))); norm ;* (4.)
;      mollview, mask
;      gnomview, mask, rot=[phi[gsource[ip]], theta[gsource[ip]], 0]*180./!pi

;      stop
      endfor                    ;
; mollview, mask
;   ismoothing, mask, smask, fwhm_arcmin=beam[ifreq], /ring, /silent, nlmax=6000
;; mollview, smask
   write_fits_map,'ffp3b_ps_map_'+sfreq[ifreq]+'.fits', mask, /ring

endif

   gnomview, 'ffp3b_ps_map_'+sfreq[ifreq]+'.fits', max=500, rot=[70,30,0], reso=.5, min=-200
   gnomview, 'ffp3-b_Imap_'+sfreq[ifreq]+'.fits', max=500, rot=[70,30,0], reso=.5, min=-200
   read_fits_map,'ffp3b_ps_map_'+sfreq[ifreq]+'.fits', p
   read_fits_map, 'ffp3-b_Imap_'+sfreq[ifreq]+'.fits', m
   gnomview, m-p, max=500, min=-200, rot=[70,30,0], reso=.5, tit='Point Sources removed'

 stop
endfor
endif
;; mollview, mask_tot
;stop
if (true) then begin
   mask_tot=fltarr(12l*2048l^2)
   mask_tot[*] = 1.
   for ifreq=ff,lf do begin
      print, sfreq[ifreq]
      read_fits_map,'/global/scratch/sd/dpietrob/ffp3b/ns2048/ffp3b_ps_map_'+sfreq[ifreq]+'.fits', map, nside=lns, ordering=ord
      mapu = map
      print, lns, ' ', ord
      if (lns ne 2048) then ud_grade, map, mapu, nside_out=2048, order_in=ord, order_out='ring'
;      gnomview, mapu, rot=[70,30,0], /asinh, chars=1.5
      help, mapu
      bi = where(mapu[*,0] gt 1.)
      help, bi
      mask_tot[bi] = 0. 
;      gnomview, mask_tot, rot=[70,30,0], chars=1.5
;      mollview, mask_tot, rot=[70,30,0], chars=1.5

   endfor
;;    mollview, mask_tot, /asinh, tit='ps mask'
;   gnomview, mask_tot, rot=[70,30,0], chars=1.5
;stop
;   cc = mask_tot*0.
;   cc[where(abs(mask_tot) gt 0)] = 0.
;   cc[where(mask_tot eq 0)] = 1.
   write_fits_map,'ffp3b_ps_mask.fits',mask_tot,/ring
   mollview, 'ffp3b_ps_mask.fits'
   read_fits_map,'mask_0.80_ns2048.fits', m80
   write_fits_map,'ffp3b_mask_80ps.fits', m80*mask_tot, /ring
   mollview, 'ffp3b_mask_80ps.fits'
endif

end
