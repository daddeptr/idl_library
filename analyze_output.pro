   true  = 1b
   false = 0b

   dataset = ['f', '1h', '2h']
   nb = ['2b', '3b', '4b', '2b_5arc']

   fb = 3
   lb = 3

   dx = 'dx7'

   mfile = '../dx7_mask_common.fits'

   tag = '-CO_tgh'

   if (true) then begin
      for ib=fb, lb do begin
          if (true) then begin
;;       read_fits_map, nb[ib]+'_f/fg_amp_map_no02_c0001_k00001.fits',m
              read_fits_map, nb[ib]+'_1h'+tag+'/fg_amp_map_no02_c0001_k00001.fits',m1
              read_fits_map, nb[ib]+'_2h'+tag+'/fg_amp_map_no02_c0001_k00001.fits',m2
              s = (m1+m2)/2.
              d = (m1-m2)/2.
              write_fits_map, dx+'_'+nb[ib]+tag+'_s.fits', float(s), /ring, units='!7l!8K CMB'
              write_fits_map, dx+'_'+nb[ib]+tag+'_d.fits', float(d), /ring, units='!7l!8K CMB'
              write_fits_map, dx+'_'+nb[ib]+tag+'_1h.fits', float(m1), /ring, units='!7l!8K CMB'
              write_fits_map, dx+'_'+nb[ib]+tag+'_2h.fits', float(m2), /ring, units='!7l!8K CMB'
          endif

          if (true) then begin
              ianafast, dx+'_'+nb[ib]+tag+'_s.fits', dx+'_'+nb[ib]+tag+'_cls.fits',maskfile=mfile, regression=1, iter=1
              ianafast, dx+'_'+nb[ib]+tag+'_d.fits', dx+'_'+nb[ib]+tag+'_nls.fits',maskfile=mfile, regression=1, iter=1
              ianafast, dx+'_'+nb[ib]+tag+'_1h.fits', dx+'_'+nb[ib]+tag+'_xcls.fits',maskfile=mfile, regression=1, iter=1, map2_in=dx+'_'+nb[ib]+tag+'_2h.fits'
          endif
      endfor
   endif
stop
   if (false) then begin

      fits2cl, cls143, 'dx7_143GHz_cls.fits'
      fits2cl, nls143, 'dx7_143GHz_nls.fits'

      mollview,randomn(-1,12), window=-1
      loadct,39
      !p.color=0
      !p.background=255
      
      cols = [245, 70, 210, 100]

      l = findgen(4097)
      ll = l*(l+1)/2./!pi
      window,1,xsize=1344, ysize=840
      plot_oi, l, smooth(cls143*ll,20), thick=2, chars=1.5, xr=[1,3000], xs=1
      oplot, l, smooth(nls143*ll,20), thick=2

      for ib=0,n_elements(nb)-1 do begin
         fits2cl, cls, dx+'_'+nb[ib]+tag+'_cls.fits'
         fits2cl, nls, dx+'_'+nb[ib]+tag+'_nls.fits'
         oplot, l, smooth(cls*ll,20), col=cols[ib], thick=2
         oplot, l, smooth(nls*ll,20), col=cols[ib], thick=2
      endfor

      fits2cl, smi, 'cls_smica.fits'
      smi = smi * 1.e12
      oplot, l, smooth(smi*ll,20), col=40, thick=2

      legend, ['143GHz', nb,'Smica'], col=[0, cols,40], line=[0,cols*0,0], chars=1.5, thick=2

   endif

   if (false) then begin
       for ib=0,n_elements(nb)-1 do begin
           mollview, 'dx7_'+nb[ib]+tag+'_s.fits', /no_monopole, gal_cut=30, min=-300, max=300, chars=1.5, png='../../ns128/png/cmb_'+nb[ib]+tag+'.png', window=-1
           mollview, 'dx7_'+nb[ib]+tag+'_d.fits', /no_monopole, gal_cut=30, min=-150, max=150, chars=1.5, png='../../ns128/png/cmb_rms_'+nb[ib]+tag+'.png', window=-1
       endfor
   endif
 
   if (false) then begin
      read_fits_map,'../../ns128/chains/cls_recCMB_smth60_ns128_RING.fits', cmb
      read_fits_map,'../../ns128/dx6_mask_0.99_ns128.fits', mask
      bp = where(mask eq 0)
      for ib=0,0 do begin
;         gb=gaussbeam(60,4096)
;         fits2cl,bl,'dx6_'+nb[ib]+'_effBeam_v2.fits'
;         cl2fits,gb/bl,'beam.fits'
;         ismoothing, 'dx7_'+nb[ib]+'-CO_s.fits', 'smth.fits', beam_file='beam.fits'
;         ud_grade, 'smth.fits', '../../ns128/cmb_'+nb[ib]+'-CO_downg.fits', nside_out=128, order_out='ring'
         read_fits_map,'../../ns128/cmb_'+nb[ib]+'-CO_downg.fits',band
         ccc = band
         remove_dipole, ccc, mask, nside=128, ordering='ring', /onlymonopole, gal_cut=30
         ccc[bp] = -1.6375e30
         mollview, ccc, min=-300, max=300, chars=1.5, tit=nb[ib]+' 60 arcmin', png='../../ns128/png/cmb_60arcmin_'+nb[ib]+'.png', window=-1

         ccc = band-cmb
         remove_dipole, ccc, mask, nside=128, ordering='ring', gal_cut=30
         ccc[bp] = -1.6375e30
         mollview, ccc, min=-50, max=50, chars=1.5, tit=nb[ib]+' 60 arcmin - CMB Posterior Average', png='../../ns128/png/h-lres_'+nb[ib]+'.png', window=-1
      endfor


;      read_fits_map,'dx7_2b-CO_s.fits', m2
;      read_fits_map,'dx7_3b-CO_s.fits', m3
;      read_fits_map,'dx7_4b-CO_s.fits', m4
;      read_fits_map,'dx7_2b_5arc-CO_s.fits', m25
      
;      mollview, m3-m2, min=-300, max=300, chars=1.5, units='!7l!8K CMB', tit='3b - 2b'
;      mollview, m3-m4, min=-300, max=300, chars=1.5, units='!7l!8K CMB', tit='3b - 4b'
;      mollview, m3-m25, min=-300, max=300, chars=1.5, units='!7l!8K CMB', tit='3b - 2b_5arc'
;      mollview, m4-m2, min=-300, max=300, chars=1.5, units='!7l!8K CMB', tit='4b - 2b'



   endif
end
