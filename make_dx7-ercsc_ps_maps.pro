pro make_dx7_ercsc_ps_maps, ff, lf

!path=!path+':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'
!path=!path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'

true = 1b
false = 0b

freq = [30,44,70,100,143,217,353, 545, 857]
sfreq= string(freq,format='(i3.3)')
nfreq = n_elements(freq)
pfreq=[28.4, 44.1, 70.3, 101.2, 142.6, 221.9, 360.6, 557.6, 857]

ps_nsdie=8192l

Nside = [1024l, 1024l, 1024l, 2048l, 2048l, 2048l, 2048l, 2048l, 2048l]
Nside[*] = 2048l
Npix = 12l*Nside^2 

beam = [33., 23., 12.7143, 9.64461, 7.11419, 4.72230, 4.5648]

;root = '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/dx7/ns2048/ercsc_catalogs/'
;file_name = 'ERCSC_f'
;tag = 'dx7_ercsc_ps_v2'
root = '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/dx7/ns2048/dx7_catalogs/'
file_name = 'DX7_PwS_'
tag = 'dx7_PwS_ps'



;mNside = 2048l
;mNpix = 12l*mNside^2

;mask_tot = fltarr(mNpix)
;mask_tot[*] = 1.

;ff=6
;lf=6

make_mask = false
get_sources = true

radec2gal = 1

if (get_sources) then begin
for ifreq=ff,lf do begin
   print, freq[ifreq]

   if (true) then begin
      cf = conversionfactor(freq[ifreq],/brightness2muKthermo)
      print, ' conversion: ', cf

      infile = root + file_name + sfreq[ifreq]+'.fits'
      print, infile
      read_fits_s, infile, hdr, table

      flux = table.FLUX
      ra = table.RA
      dec = table.DEC
      glat = table.GLAT
      glon = table.GLON
      extended = table.EXTENDED

      s2n = table.FLUX/table.FLUX_ERR

      help, flux

      gsource = where( (extended eq 0) and (s2n ge 5.) )
gsource=gsource[0]
      help, gsource
   
      print, 'ra ',minmax(ra[gsource])
      print, 'dec ',minmax(dec[gsource])

;      euler, ra[gsource], dec[gsource], phi, theta, radec2gal
      theta = glat
      phi = glon
      print, 'theta ',minmax(theta)
      print, 'phi ',minmax(phi)

      phi = phi/!radeg
      theta = !pi/2 - theta/!radeg
      
      print, 'theta ',minmax(theta)
      print, 'phi ',minmax(phi)

      ang2vec, theta, phi, vec

      ns = 4096l ;2*nside[ifreq]
      npx = 12l * ns^2
      ipix = lindgen(npx)
      pix2ang_ring, ns, ipix, ptheta, pphi
      pix2vec_ring, ns, ipix, pvec

      ngsource = n_elements(gsource)

      mask = fltarr( npx )
      parea = 4.*!pi / npx

      beam_rad = beam[ifreq]/2. / 60. * !pi/180.
 
      barea = (beam[ifreq]/2./60.*!pi/180.)^2

      radius = 2.5*beam[ifreq]/2. * !pi/180.

      pix2vec_ring, ns, 0, v0
      cdot = reform( transpose(pvec[*,*]) ## reform(v0) ) 
      thetas = acos(cdot)
      query_disc, ns, v0, radius, listpix
;      g = where(cdot gt 0.)
;;       expo = cdot * 0.
;;       expo[listpix] = exp( -(1.-cdot[listpix]^2)/(2.*beam_rad^2) ) ;/ (2.*!pi*beam_rad^2)
      expo = cdot * 0.
      expo[listpix] = exp( -(thetas[listpix]^2)/(2.*beam_rad^2) ) ;/ (2.*!pi*beam_rad^2)

      norm = total(expo)
      print, norm
      print, 1./norm

      expo = expo / norm
;      gnomview, expo, reso=1.5, chars=1.5, grat=[10,1], rot=[0,90,0]
;stop
      for ip=0,ngsource-1 do begin

         query_disc, ns, vec[ip,*], radius, listpix
         cdot = transpose(pvec[listpix,*]) ## reform(transpose(vec[ip,*]))

         cdot = reform(cdot)
         thetas = acos(cdot)
;      mask = mask + flux[gsource[ip]] * cf * 1.e-9/parea * exp(-thetas^2/(2.*beam_rad^2))/(2.*!pi*beam_rad^2)
;;          mask[listpix] = mask[listpix] + flux[gsource[ip]] * cf * 1.e-9/parea * exp(-(1.-cdot^2)/(2.*beam_rad^2)) / norm ;total(exp(-(1.-cdot^2)/(2.*beam_rad^2))); norm ;* (4.)
         mask[listpix] = mask[listpix] + flux[gsource[ip]] * cf * 1.e-9/parea * exp(-(thetas^2)/(2.*beam_rad^2)) / norm ;total(exp(-(1.-cdot^2)/(2.*beam_rad^2))); norm ;* (4.)

;;          mask[listpix] = mask[listpix] + flux[gsource[ip]] * cf * exp(-(1.-cdot^2)/(2.*beam_rad^2)) / norm ;total(exp(-(1.-cdot^2)/(2.*beam_rad^2))); norm ;* (4.)

gnomview, mask, chars=1.5, min=-550, max=550, rot=[glon[gsource[ip]],glat[gsource[ip]],0], grat=[5,5]
stop
      endfor                    ;

   ud_grade, mask, maskd, nside_out=Nside[ifreq], order_in='ring', order_out='ring'
   mask = 0.
   write_fits_map,'single_src_'+tag+'_map_'+sfreq[ifreq]+'.fits', maskd, /ring
   maskd = 0.

endif

   gnomview, tag+'_map_'+sfreq[ifreq]+'.fits', rot=[70,45,0], reso=1.5, max=500, min=-200, grat=[5,5]
   gnomview, '/global/scratch/sd/dpietrob/dx7_delivery_analysis/input/dx7_Imap'+sfreq[ifreq]+'_f.fits', rot=[70,45,0], reso=1.5, grat=[5,5], min=-400, max=400
   read_fits_map,tag+'_map_'+sfreq[ifreq]+'.fits', p
   read_fits_map, '/global/scratch/sd/dpietrob/dx7_delivery_analysis/input/dx7_Imap'+sfreq[ifreq]+'_f.fits', m
   gnomview, m-p, rot=[70,45,0], reso=1.5, tit='Point Sources removed', grat=[5,5], min=-400, max=400
   write_fits_map, tag+'_ps-removed_map_'+sfreq[ifreq]+'.fits', m-p, /ring, units='!7l!8K CMB'
;   mollview, tag+'_ps-removed_map_'+sfreq[ifreq]+'.fits', chars=1.5, px=1200, /asinh

   p=0.
   m=0.

 stop
endfor
endif
;; mollview, mask_tot


;stop

if (make_mask) then begin
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
