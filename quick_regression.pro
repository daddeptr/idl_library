;; pro make_plots

!path=!path+':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'
!path=!path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'

true  = 1b
false = 0b
mollview, randomn(-1,12), window=-1
loadct, 39
!p.background=255
!p.color=0
!p.charsize=2.5
!p.thick = 2

sbeam = '60'
;##tag = 'paper_cls'
tag = 'paper_pix_wnl'
figs_dir = 'pics/'

;#tag = 'outcmd_w_nless'

do_smth = False
minus_cmb = false
do_png = False

;run_tag = 'cls_wmap_tot_dust'
;run_tag = 'test_nless'
;##run_tag = 'test_w_nless'
;##run_tag = 'pix_wmap_standard_nless'
run_tag = 'cls_wmap_standard' ;'check'

restore,'chains/'+run_tag+'.res.sav'
mondip = total(foreg,1) / nsample

print, ' Monopole and Dipole channel by channel:'
print, mondip

print, 'Rms'
mdrms = mondip*0.
for i=0,nsample-1 do mdrms = mdrms + (foreg[i,*,*]-mondip) ^2
mdrms = mdrms / nsample
mdrms = sqrt(mdrms)
print, mdrms
print, '---------------------'

freq=[22.8, 33., 41., 61., 94.]
sfreq = ['K', 'Ka', 'Q', 'V', 'W']
nfreq = n_elements(freq)
ns = 128l
sns = strtrim(string(ns),2)
npix = 12l*ns^2
ncl = 350l
indx = lindgen(nfreq)
cols = 200 + indx*10

real_beam=[0.88, 0.66,  0.51,   0.35,   0.22] * 60.

; ------ WMAP power laws regression ------------------------------------
if (true) then begin

   restore,'chains/'+run_tag+'.res.sav'
   if (do_smth) then begin
       ismoothing, cmb, scmb, fwhm_arcmin=60, /ring
       cmb = scmb
   endif

   read_fits_map,'wmap7_temperature_analysis_mask_ns128.fits',mask
   bp = where(mask eq 0.)
   gp = where(mask eq 1.)
   ngp = n_elements(gp)

   print, float(n_elements(where(glb_chi2 gt 15)))/n_elements(glb_chi2)

   ; regressing Foreground @ 23 GHz
   label = ['Low-freq. Comp. @ 23 GHz','Dust @ 94 GHz','Synchrotron @ 408 MHz']

   if (true) then begin
       tamp = amp
       tamp2 = amp2
       tamp[*,0] = amp[*,2]
       tamp[*,1] = amp[*,0]
       tamp[*,2] = amp[*,1]
       amp = tamp
       
       tamp2[*,0] = amp2[*,2]
       tamp2[*,1] = amp2[*,0]
       tamp2[*,2] = amp2[*,1]
       amp2 = tamp2
   endif

   test_lab = ['K','W']

   for iamp=1,namp[0]-1 do begin
;##   for iamp=0,namp[0]-1 do begin
       map = amp[*,iamp]
       rms = sqrt(amp2[*,iamp]-amp[*,iamp]^2)

;##   rms[*] = 1.

;   tag = 'incmdK'
;   read_fits_map,'map_'+test_lab[iamp]+'_v2.fits',map
;   read_fits_map,'rms_'+test_lab[iamp]+'_v2.fits',rms
; ---
;   tag = 'inwK'
;   read_fits_map,'nless_map_'+test_lab[iamp]+'.fits',map
;   read_fits_map,'nless_rms_'+test_lab[iamp]+'.fits',rms

       cc = map
       cc[bp] = -1.6375e30
       mollview, cc, chars=1.5, units='!7l!8K CMB', px=1000, min=-200, max=400, tit='!17Map @ 23 GHz', win=0
       if (do_png) then mollview, cc, chars=1.5, units='!7l!8K CMB', min=-200, max=400, ps=figs_dir+'reg_chk_'+tag+'.eps', window=-1, tit='!17Map @ 23 GHz'
       
       cc = rms
       cc[bp] = -1.6375e30
       mollview, cc, chars=1.5, units='!7l!8K CMB', px=1000, win=1
       if (do_png) then mollview, cc, chars=1.5, units='!7l!8K CMB', ps=figs_dir+'reg_chk_'+tag+'_rms.eps', window=-1

       nfg = 3
       A = dblarr(nfg)
       oneOsigma = dblarr(nfg, nfg)
       sigma = oneOsigma
       a_err = dblarr(nfg)

       temp = dblarr(npix,nfg)
       read_fits_map, 'lambda_haslam_60arcmin_ns128_ring.fits', t
       temp[*,0] = t
       read_fits_map, 'lambda_halpha_60arcmin_ns128_ring.fits', t
       temp[*,1] = t
       read_fits_map, 'lambda_fds_94GHz_uKantenna_60arcmin_ns128_ring.fits', t
       temp[*,2] = t

       T = dblarr(nfg, nfg)
       B = dblarr(nfg)
       for ifg=0,nfg-1 do begin
;##           B[ifg] = total( map[gp]*temp[gp,ifg] / rms[gp]^2 ) / ngp
           B[ifg] = total( (map[gp]-mean(map[gp]))*(temp[gp,ifg]-mean(temp[gp,ifg])) / rms[gp]^2 ) / ngp
           for jfg=ifg,nfg-1 do begin
;##               T[ifg,jfg] = total( temp[gp,ifg]*temp[gp,jfg] / rms[gp]^2 ) / ngp
               T[ifg,jfg] = total( (temp[gp,ifg]-mean(temp[gp,ifg]))*(temp[gp,jfg]-mean(temp[gp,jfg])) / rms[gp]^2 ) / ngp

               T[jfg,ifg] = T[ifg,jfg]

               oneOsigma[ifg,jfg] = 2. * T[ifg,jfg]
               oneOsigma[jfg,ifg] = oneOsigma[ifg,jfg]
           endfor
       endfor

       Tm1 = invert(T, /double, status)
       print, 'Inversion status: ', status
       
       A = Tm1 ## B
       
       sigma[*,*] = invert(reform(oneOsigma[*,*]),/double, status)
       print, 'Inversion status sigma: ', status
       for kfg = 0, nfg-1 do begin
           if (sigma[kfg,kfg] gt 0.) then begin
               a_err[kfg] = sqrt(sigma[kfg,kfg])
           endif else begin
               print, 'Sigma Error'
           endelse
       endfor

       print, ' ============================================================'
       print, 'Coefficients: Haslam, f-f, Dust'

       print, 'A:     ', A
       print, 'A_err: ', A_err
       print, ' ============================================================'
   
       fit = map*0.
       for ifg=0,nfg-1 do fit=fit + temp[*,ifg]*A[ifg]
       fit[bp] = -1.6375e30
       mollview, fit, chars=1.5, tit='!17Template Linear Combination: WMAP7', min=-200, max=400, px=1000, win=2
;   mollview, fit, chars=1.5, tit='!17Template Linear Combination: WMAP7', min=-200, max=400, ps=figs_dir+'wmap-pl_fit_fg_'+strtrim(string(iamp+1),2)+'.eps', window=-1

       if not (minus_cmb) then res = map-fit
       if (minus_cmb) then res = map-fit-cmb
       res[bp] = -1.6375e30
       mollview, res, chars=1.5, tit='!17Residual Foregrounds', min=-100, max=100, units='!7l!8K Antenna',px=1000, win=3
       if (do_png) then mollview, res, chars=1.5, tit='!17Residual Foregrounds', min=-200, max=200, units='!7l!8K Antenna', ps=figs_dir+'reg_chk_'+strtrim(string(iamp+1),2)+'_'+tag+'.eps',window=-1

       if (iamp eq 0) then begin
           hazepix=-1
           thres = 3.
           sig = res/rms
           hazepix = where((sig) gt thres)
           if (hazepix[0] lt 0) then begin
               print, 'No haze'
               stop
           endif
           spuriouspix = where(-(sig) gt thres)
           
           disc = sig*0.-1.6375e30
           query_disc, ns, [1.,0.,0.], !pi/5, listpix
           disc[listpix] = res[listpix]
           centerhazepix = where(disc/rms gt thres)
           mollview, disc, chars=1.5, px=1000
       endif
       
       if not (minus_cmb) then hchisq=((map-fit)/rms)^2
       if (minus_cmb) then hchisq=((map-fit-scmb)/rms)^2
       mollview, hchisq, chars=1.5, px=1000, tit='!7v!u2!8', max=16, win=4
       if (do_png) then mollview, hchisq, chars=1.5, tit='!7v!u2!8', max=16, ps=figs_dir+'reg_chk_chisq_'+strtrim(string(iamp+1),2)+'_'+tag+'.eps'
       print, ' Fit chisq: ', mean(hchisq[gp])
       
       loadct, 39
       !p.background=255
       !p.color=0
       !p.charsize=2.5
       !p.thick = 2
       window, 2+iamp, xsize=1344, ysize=840
       plot, map[gp], fit[gp], psym=3, ytit='!17Fit Amplitude', xtit='!17Commander Posterior Mean', tit='!17'+label[iamp]+': WMAP7'
       oplot, fit[gp], fit[gp], psym=3, col=245
       oplot, map[hazepix], fit[hazepix], col=70, psym=1
       oplot, map[spuriouspix], fit[spuriouspix], col=150, psym=1
       oplot, map[centerhazepix], fit[centerhazepix], col=210, psym=2
       legend, ['!17Positive 3 !7r!17 residuals', '!17Negative 3 !7r!17 residuals', '!17Haze'], col=[70,150,210], psym=[1,1,2], chars=2
       
       if not (minus_cmb) then print, ' Regression coefficient: ', mean( (map[gp] * fit[gp])/sqrt(map[gp]^2 * fit[gp]^2) ), correlate(map[gp], fit[gp])
       if (minus_cmb) then print, ' Regression coefficient: ', mean( (map[gp]-scmb[gp]) * fit[gp]/sqrt((map[gp]-scmb[gp])^2 * fit[gp]^2) ), correlate(map[gp]-scmb[gp], fit[gp])
stop
; --- Plotting
      set_plot, 'ps'
      igp = lindgen(ngp/10)
      igp = igp * 10
      device, file=figs_dir+'wmap-pl_scatter_'+strtrim(string(iamp),2)+'.eps', /col, bits=8
      plot, map[gp[igp]], fit[gp[igp]], psym=3, ytit='!17Fit Amplitude', xtit='!17Commander Posterior Mean', tit='!17'+label[iamp]+': WMAP7'
      oplot, fit[gp[igp]], fit[gp[igp]], psym=3, col=245
      oplot, map[hazepix], fit[hazepix], col=70, psym=1
      oplot, map[spuriouspix], fit[spuriouspix], col=150, psym=1
      oplot, map[centerhazepix], fit[centerhazepix], col=210, psym=2
      legend, ['!17Positive 3 !7r!17 residuals', '!17Negative 3 !7r!17 residuals', '!17Haze'], col=[70,150,210], psym=[1,1,2], chars=2
      device, /close
      set_plot, 'x'
;stop
; ---
   endfor
   STOP
endif


print, ' --- End of Procedure ---'
STOP

end
