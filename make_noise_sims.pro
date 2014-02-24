;## pro make_noise_sims,
;##
;## Procedure:
;## 1) Making half-difference map or weighted difference map
;## 2) Normalize the hd map by its RMS
;## 3) Compute the power spectrum of that map
;## 3a) Smoothing it: only if no mask is applied, otherwise
;mask-induced coupling is going to show up
;## 4) Generate simulations for that power spectrum
;## 5) Multiply these sims by the RMS or the hd map


; Generalized to generic DX adding lfi/hfi_dx_tag
; 
; Polarization added in maps only (no variance)
; DX10
;
; Only Delta_dx9 is used for the first delivery. 353 is still the old one.
; Most updated and probably correct one

;##   !path = !path + ':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'

   True = 1b
   False = 0b

   do_moll = True
   do_pol  = True

   fstf = 4 
   lstf = 4

   freq=[30, 44, 70, 100, 143, 217, 353, 545, 857]
   sfreq = string(freq, format='(i3.3)')

   cv_MJy2Kcmb = 1.d0 / [1.,1.,1.,1.,1.,1.,1., 57.976641, 2.2690745] ; only for 545 and 857 channels

   nfreq = n_elements(freq)

   dx_version = 'dx11_pre'
   pre_tag = 'DetsetMap';'DetsetMap'; 'SkyMap'
   lfi_dx_tag1 = '_1024_DX10'
   lfi_dx_tag2 = '_1024_DX10'
   hfi_dx_tag1 = '-ds1_2048_v60';'_2048_v60'
   hfi_dx_tag2 = '-ds2_2048_v60';'_2048_v60'
   mission = 'full' ;'full'
   mission_tag = '_ds_1-2'
   dir = '/project/projectdirs/planck/data/mission/DPC_maps/'+dx_version+'/'
   hfi_dir = dir+'hfi/extra/standard_dpc_gains/FREQ/'
   lfi_dir = dir+'../dx10/lfi/'
 ;  tag = '_ds2'

   nsim = 50

   beam_dir = '/global/scratch/sd/paganol/run_FEBeCoP/output_dx9/transFn/Bl/'
   real_beam = [32.65, 27.00, 13.01, 9.94,   7.04,    4.66,    4.41,    4.47,    4.23]

   beams = 40.
   eff_beams = sqrt( float(beams)^2-real_beam^2 )

   spawn, 'mkdir -p -m g+xr /global/scratch2/sd/dpietrob/'+dx_version
   spawn, 'mkdir -p -m g+xr /global/scratch2/sd/dpietrob/'+dx_version+'/noise'

   bad_value = -1.63750e+30

   for ifreq=fstf,lstf do begin
if False then begin
       spawn, 'mkdir -p -m g+xr /global/scratch2/sd/dpietrob/'+dx_version+'/noise/'+sfreq[ifreq]
       out_dir = '/global/scratch2/sd/dpietrob/'+dx_version+'/noise/'+sfreq[ifreq]+'/'
       print, sfreq[ifreq]
       if (ifreq lt 3) then begin
           print, ''+lfi_dir+'LFI_'+pre_tag+'_'+sfreq[ifreq]+lfi_dx_tag1+'_'+mission+'.fits'
           spawn, 'ls '+lfi_dir+'LFI_'+pre_tag+'_'+sfreq[ifreq]+lfi_dx_tag1+'_'+mission+'.fits', files1
           print, ''+lfi_dir+'LFI_'+pre_tag+'_'+sfreq[ifreq]+lfi_dx_tag2+'_'+mission+'.fits'
           spawn, 'ls '+lfi_dir+'LFI_'+pre_tag+'_'+sfreq[ifreq]+lfi_dx_tag2+'_'+mission+'.fits', files2
       endif else begin
           print, ''+hfi_dir+'HFI_'+pre_tag+'_'+sfreq[ifreq]+hfi_dx_tag1+'_'+mission+'.fits'
           spawn, 'ls '+hfi_dir+'HFI_'+pre_tag+'_'+sfreq[ifreq]+hfi_dx_tag1+'_'+mission+'.fits', files1
           print, ''+hfi_dir+'HFI_'+pre_tag+'_'+sfreq[ifreq]+hfi_dx_tag2+'_'+mission+'.fits'
           spawn, 'ls '+hfi_dir+'HFI_'+pre_tag+'_'+sfreq[ifreq]+hfi_dx_tag2+'_'+mission+'.fits', files2
       endelse

       help, files1, files2
       print, files1, files2

       read_fits_map, files1, map1, hdr, xhdr, nside=dpns1, order=dpord1
       read_fits_map, files2, map2, hdr, xhdr, nside=dpns2, order=dpord2
       print, ' reordering...'
       if (strupcase(dpord1) eq 'RING') then map1 = reorder(map1, in=dpord, out='NESTED')
       if (strupcase(dpord2) eq 'RING') then map2 = reorder(map2, in=dpord, out='NESTED')

       tbp1 = where(map1[*,0] eq bad_value)
       pbp1 = where( (map1[*,1] eq bad_value) or (map1[*,2] eq bad_value) )
;       gp1 = where(map1[*,0:2] ne bad_value)
       tbp2 = where(map2[*,0] eq bad_value)
       pbp2 = where( (map2[*,1] eq bad_value) or (map2[*,2] eq bad_value) )
;       gp2 = where(map2[*,0:2] ne bad_value)

;print, "WARNING: to be refind taking into account T-P correlation: variance for now."
;stop
       d = ( map1[*,0:2] - map2[*,0:2] ) / 2. * 1.e6

       var = [ [(map1[*,4]+map2[*,4])], [(map1[*,7]+map2[*,7])], [(map1[*,9]+map2[*,9])] ] / 4.
       rms = sqrt( var ) * 1.e6

       d = d / rms

       if tbp1[0] ne -1  then d[tbp1,0] = bad_value
       if tbp2[0] ne -1  then d[tbp2,0] = bad_value

       if pbp1[0] ne -1  then d[pbp1,1:2] = bad_value
       if pbp2[0] ne -1  then d[pbp2,1:2] = bad_value

       if (tbp1[0] ne -1) then begin
           print, ' :DP - Missing Pixels found in T...1'
           help, tbp1
           n_bp = n_elements(tbp1)
           for ipix=0l,n_bp-1 do begin
               pix2vec_ring, dpns1, tbp1[ipix], vec0
               query_disc, dpns1, vec0, 2.*real_beam[ifreq]/60./!radeg, listp, /nest
               gp = where(d[listp,0] ne bad_value)
               listp = listp[gp]
               d[tbp1[ipix],0] = median(d[listp,0])
               rms[tbp1[ipix],0] = median(rms[listp,0])
           endfor
       endif

       if (tbp2[0] ne -1) then begin
           print, ' :DP - Missing Pixels found in T...2'
           help, tbp2
           n_bp = n_elements(tbp2)
           for ipix=0l,n_bp-1 do begin
               pix2vec_ring, dpns2, tbp2[ipix], vec0
               query_disc, dpns2, vec0, 2.*real_beam[ifreq]/60./!radeg, listp, /nest
               gp = where(d[listp,0] ne bad_value)
               listp = listp[gp]
               d[tbp2[ipix],0] = median(d[listp,0])
               rms[tbp2[ipix],0] = median(rms[listp,0])
           endfor
       endif

       if (pbp1[0] ne -1) then begin
           print, ' :DP - Missing Pixels found in P...1'
           help, tbp1
           n_bp = n_elements(pbp1)
           for ipix=0l,n_bp-1 do begin
               pix2vec_ring, dpns1, pbp1[ipix], vec0
               query_disc, dpns1, vec0, 2.*real_beam[ifreq]/60./!radeg, listp, /nest
               gp = where( (d[listp,1] ne bad_value) and (d[listp,2] ne bad_value))
               listp = listp[gp]
               d[pbp1[ipix],1] = median(d[listp,1])
               d[pbp1[ipix],2] = median(d[listp,2])
               rms[pbp1[ipix],1] = median(rms[listp,1])
               rms[pbp1[ipix],2] = median(rms[listp,2])
           endfor
       endif

       if (pbp2[0] ne -1) then begin
           print, ' :DP - Missing Pixels found in P...2'
           help, pbp2
           n_bp = n_elements(pbp1)
           for ipix=0l,n_bp-1 do begin
               pix2vec_ring, dpns2, pbp2[ipix], vec0
               query_disc, dpns2, vec0, 2.*real_beam[ifreq]/60./!radeg, listp, /nest
               gp = where( (d[listp,1] ne bad_value) and (d[listp,2] ne bad_value))
               listp = listp[gp]
               d[pbp2[ipix],1] = median(d[listp,1])
               d[pbp2[ipix],2] = median(d[listp,2])
               rms[pbp2[ipix],1] = median(rms[listp,1])
               rms[pbp2[ipix],2] = median(rms[listp,2])
           endfor
       endif

       rrms = reorder(rms, in='nest', out='ring')

       clfile = out_dir+dx_version+mission_tag+'_noise_'+sfreq[ifreq]+'_ns'+string(dpns1,format='(i4.4)')+'_uK_cls.fits'
       ianafast, d, cls, /nest, simul_type=2
       cl2fits, cls, clfile
endif
       fits2cl, nls, clfile

;       mollview, d, win=0, px=650
;       loadct, 39
;       !p.color=0
;       !p.background=255

       pwf = healpixwindow(dpns1)
       snls = nls
       l = findgen(4097)
       ll = l*(l+1)/2./!pi
       il = lindgen(4095)+2

;## EXP
;       pars = [ [ 1.e-7, 12., -5, 1.e-7, 100., -1.15, 1.e-7, 3500., 3700., 0.], $
;                [ 1.e-7, 12., -5., 1.e-7, 100., -1.15, 1.07e-7, 3500., 3700., 0.], $
;                [ 0., 12., -5., 1.e-7, 200., -1.15, 1.07e-7, 3500., 3700., 0.], $
;                [ 1.e-7, 12., -4., 0., 100., -1.15, -1.e-7, 3500., 3700., 0.], $
;                [ -5.e-7, 15., -4., 1.e-7, 100., -1.15, -4.5e-7, 25., 25., 0.], $
;                [ -1.e-6, 10., -6., 1.5e-7, 80., -2.15, 0., 3500., 3700., 0.] ]
                
;## COS
       pars = [ [ 1.e-7, 3000., 0., 1.e-7, 12., -5., 1.e-7, 100., -1.15, 0.], $
                [ 1.e-7, 3000., 0., 1.e-7, 12., -5., 1.e-7, 100., -1.15, 0.], $
                [ 1.e-7, 3000., 0., 1.e-7, 12., -5., 1.e-7, 100., -1.15, 0.], $
                [ 1.e-7, 3000., 0., 1.e-7, 12., -5., 1.e-7, 100., -1.15, 0.], $
                [ 1.e-7, 3000., 0., 1.e-7, 12., -5., 1.e-7, 100., -1.15, 0.], $
                [ 1.e-7, 3000., 0., 1.e-7, 12., -5., 1.e-7, 100., -1.15, 0.] ]

                
       sigma_pars = pars*0.
       fitpars = [ [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 0], $
                   [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 0], $
                   [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 0], $
                   [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 0], $
                   [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 0], $
                   [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 0] ]


       tol = 1.e-9

       !p.multi=[0,3,1]
       for i=0,5 do begin
           print, 'icl = ', i
           nls_err = sqrt(2./(2*l+1)) * nls[*,i]
           weights = 1. / ( nls_err )
           bl = bp_binning(l, 'data/bins/ctp/CTP_bin_TT')
           bll = bp_binning(ll, 'data/bins/ctp/CTP_bin_TT')
           y = bp_binning(nls[*,i], 'data/bins/ctp/CTP_bin_TT')
           yll = bp_binning(nls[*,i]*ll, 'data/bins/ctp/CTP_bin_TT')
           er = sqrt( bp_binning( nls_err^2, 'data/bins/ctp/CTP_bin_TT') )
;          erll = sqrt( bp_binning( nls_err^2*ll^2, 'data/bins/ctp/CTP_bin_TT') )
;           y = bp_binning(nls[*,i]*ll*(2*l+1), 'data/bins/ctp/CTP_bin_TT')
;           er = sqrt( bp_binning( nls_err^2*ll^2*(2*l+1)^2, 'data/bins/ctp/CTP_bin_TT') )
           w = 1./er^2 * 8. ;sqrt( bp_binning( weights, 'data/bins/ctp/CTP_bin_TT') )
;           wll = 1./erll^2 * 4. ;sqrt( bp_binning( weights, 'data/bins/ctp/CTP_bin_TT') )
           window, 0, xsize=720*1.8, ysize=450*1.25
           plot_oi, bl, (y), xr=[10,4000], tit=i+1, ytit='!6C!dl!n', chars=3.
           errplot, bl, (y-er), (y+er)
           ipars = reform( pars[*,i] )
           ifitpars = reform( fitpars[*,i] )
;           fit = curvefit( l[il], nls[il,i], weights[il], ipars, sigma_pars, function_name='noise_func', tol=tol, status=err, chisq=chi2, yerror=yr, fita=ifitpars, /double)
;##           fit = curvefit( bl, y, w, ipars, isig, function_name='noise_func', tol=tol, status=err, chisq=chi2, yerror=yr, fita=ifitpars, /double, itmax=100 )
           fit = curvefit( bl, y, w, ipars, isig, function_name='noise_func_cos', tol=tol, status=err, chisq=chi2, yerror=yr, fita=ifitpars, /double, itmax=100 )
;           fit = curvefit( bl, yll, wll, ipars, sigma_pars, function_name='noise_func', tol=tol, status=err, chisq=chi2, yerror=yr, fita=ifitpars, /double, iter=50 )
           pars[*,i] = ipars
           sigma_pars[*,i] = isig
           oplot, bl, fit, col=245
;           oplot, l[il], nls[il,i]*ll[il], col=200, psym=3
;##           noise_func, l[il], ipars, f
           noise_func_cos, l[il], ipars, f
           snls[il,i] = f
           oplot, l[il], f, col=100
           print, pars[*,i]
           print, sigma_pars[*,i]
           print, err
           print, chi2, total( (f-nls[il,i])^2/weights^2 ) / 4095, total( (fit-y)^2/er^2 ) / n_elements(y)
;           print, yr

           plot_oo, bl, yll, xr=[1,4000], tit=i+1, psym=2, ytit='!6D!dl!n', chars=3.
           errplot, bl, (yll+erll), (yll+erll)
           oplot, l[il], f*ll[il], col=100

           plot_oi, bl, (fit*bll-yll)/fit/bll, xr=[1,4200], ytit='!6D!dl!n residuals', chars=3.
;           errplot, bl, y-er, y+er
;           oplot, bl, fit, col=245
;           oplot, bl, fit-y, psym=2
;           errplot, bl, fit-y-er, fit-y+er
           oplot, bl, bl*0, lin=2, col=245
;           oplot, l[il], nls[il,i], col=200, psym=3
;           oplot, l[il], f, col=100
           stop
       endfor

;       isim = 1
stop
       for isim=1l,nsim do begin
           print, isim
;           bl2fits, [[findgen(8193)+1.],[findgen(8193)+1.]], 'tmpwindow.fits'
           isynfast, snls, sim, simul_type=2, nside=dpns1, nlmax=2*dpns1, iseed=-isim
           sim = sim * rrms

;       ianafast, nsim, ols, /ring, simul_type=2, /show_cl
;           mollview, sim, win=1
           mapfile = out_dir+dx_version+mission_tag+'_noise_'+sfreq[ifreq]+'_ns'+string(dpns1,format='(i4.4)')+'_uK_'+string(isim, format='(i5.5)')+'.fits'
           if (not do_pol) then write_fits_map, mapfile, sim, /ring, units='!7l!6K CMB' else write_tqu, mapfile, sim, /ring, units='!7l!6K CMB'
           if (do_moll) then begin
               mymoll, mapfile, min=-50, max=50, px=500, win=ifreq, imap=0
               mymoll, mapfile, min=-50, max=50, px=500, win=ifreq+1, imap=1
           endif
       endfor
stop                   
   endfor   

   print, ' --- End of program --- '
stop

end
