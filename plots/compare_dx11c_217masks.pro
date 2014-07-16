pro double_power_law, x, A, F, pder

;## Amplitudes
   A[0] = abs( A[0] )
   A[2] = abs( A[2] )
;## INdices
   A[1] = min( [-0.125, A[1]] )
   A[1] = max( [-3.5, A[1]] )
   A[3] = min( [-3.55, A[3]] )
   A[3] = max( [-10, A[3]] )

   A = double( A )

   xterm = exp(double(x))^(A[3]-A[1])

   F = alog(A[0]) + double(x)*A[1] + alog( 1. + A[2]/A[0]*xterm )

   den = ( 1. + A[2]/A[0]*xterm )
   pder = [ [x*0.+1/A[0]-A[2]/A[0]^2*xterm/den], [double(x) - A[2]/A[0]*xterm*x/den], [xterm/den/A[0]], [A[2]/A[0]*xterm*x/den] ]

end

pro double_power_law_lin, x, A, F, pder

; Amplitudes
   A[0] = abs( A[0] )
   A[2] = abs( A[2] )
;## INdices
   A[1] = min( [ -0.125, A[1]] )
   A[1] = max( [-3.0, A[1]] )
   A[3] = min( [ -3.5, A[3]] )
   A[3] = max( [-10, A[3]] )

   A = double( A )
   x = double( x )

   F = A[0]*x^A[1] + A[2]*x^A[3]

   pder = [ $
            [x^A[1]], $
            [A[0] * x^A[1] * alog(x)], $
            [x^A[3]], $
            [A[2] * x^A[3] * alog(x)] $
          ]
end

pro double_power_law_lin_norm, x, A, F, pder

; Amplitudes
   A[0] = abs( A[0] )
   A[2] = abs( A[2] )
;## INdices
   A[1] = min( [ -0.125, A[1]] )
   A[1] = max( [-3.0, A[1]] )
   A[3] = min( [ -3.5, A[3]] )
   A[3] = max( [-10, A[3]] )

   A = double( A )
   x = double( x )

   F = A[0]*(x/3000.)^A[1] + A[2]*(x/100.)^A[3]

   pder = [ $
            [(x/3000.)^A[1]], $
            [A[0] * (x/3000.)^A[1] * alog(x/3000.)], $
            [(x/100.)^A[3]], $
            [A[2] * (x/100.)^A[3] * alog(x/100.)] $
          ]
end

pro triple_power_law, x, A, F, pder

; Amplitudes
   A[0] = abs( A[0] )
   A[2] = abs( A[2] )
   A[4] = abs( A[4] )
; Indices
;   A[1] = max([-5.5, A[1]] )
;   A[1] = min([-4.1, A[1]] )
;##   A[3] = max([-2.55, A[3]] )
;##   A[3] = min([-2.35, A[3]] )
   A[1] = max([-6.5, A[1]] )
   A[1] = min([-3.1, A[1]] )
;   A[3] = max([-2.7, A[3]] )
   A[3] = max([-3.0, A[3]] )
   A[3] = min([-2.0, A[3]] )
; For CMB channels
;##   A[5] = max([0, A[5]] )
;##   A[5] = min([6, A[5]] )
; For 545
;## ;   A[5] = max([-1.05, A[5]] )
;## ;   A[5] = min([-0.55, A[5]] )
   A[5] = max([-1.15, A[5]] )
   A[5] = min([-0.15, A[5]] )

   A = double( A )
;print, 'A=', A
;## print, double(exp(x))
   xterm1 = exp(double(x))^(A[3]-A[1])
   xterm2 = exp(double(x))^(A[5]-A[1])
;print, 'exp=',xterm
   F = alog(A[0]) + double(x)*A[1] + alog( 1. + A[2]/A[0]*xterm1 + A[4]/A[0]*xterm2 )
;print, 'F=',F
   den = ( 1. + A[2]/A[0]*xterm1 + A[4]/A[0]*xterm2 )
   pder = [ $
            [x*0.+1/A[0]-(A[2]/A[0]^2*xterm1 + A[4]/A[0]^2*xterm2) / den], $
            [double(x) * (1. + A[2]/A[0]*xterm1 + A[4]/A[0]*xterm2) / den], $
            [xterm1/den/A[0]], $
            [A[2]/A[0]*xterm1*x/den], $
            [xterm2/den/A[0]], $
            [A[4]/A[0]*xterm2*x/den] $
          ]
;print, 'der=',pder
end

pro triple_power_law_lin, x, A, F, pder

; Amplitudes
   A[0] = abs( A[0] )
   A[2] = abs( A[2] )
   A[4] = abs( A[4] )
; Indices
;   A[1] = max([-5.5, A[1]] )
;   A[1] = min([-4.1, A[1]] )
;##   A[3] = max([-2.55, A[3]] )
;##   A[3] = min([-2.35, A[3]] )
   A[1] = max([-6.5, A[1]] )
   A[1] = min([-3.5, A[1]] )
   A[3] = max([-3.0, A[3]] )
   A[3] = min([-1.5, A[3]] )
; For CMB channels
;##   A[5] = max([0, A[5]] )
;##   A[5] = min([6, A[5]] )
; For 545
;##;   A[5] = max([-1.05, A[5]] )
;##;   A[5] = min([-0.55, A[5]] )
   A[5] = max([-1.15, A[5]] )
   A[5] = min([-0.15, A[5]] )

   A = double( A )
   x = double( x )
;print, 'A=', A
;## print, double(exp(x))
;   xterm1 = exp(double(x))^(A[3]-A[1])
;   xterm2 = exp(double(x))^(A[5]-A[1])
;print, 'exp=',xterm
   F = A[0]*x^A[1] + A[2]*x^A[3] + A[4]*x^A[5]

   pder = [ $
            [x^A[1]], $
            [A[0] * x^A[1] * alog(x)], $
            [x^A[3]], $
            [A[2] * x^A[3] * alog(x)], $
            [x^A[5]], $
            [A[4] * x^A[5] * alog(x)] $
          ]
;print, 'der=',pder
end


; Filter 545 ell< 500 and remake spectra for different masks
; Spectrum with PS mask only

;   !path=!path+':/global/scratch2/sd/dpietrob/Software/src/pro/myastrolib/pro/'
;   !path=!path+':/global/scratch2/sd/dpietrob/Software/src/pro/'

   True = 1b
   False = 0b

   do_png           = True

   do_545spectraHigh= False
   do_545spectra    = False
   do_545spectraFix = False
   do_545diff       = False
   do_545diffFix    = False

;   do_217diff      = False
;   do_143diff      = False
;   do_100diff      = False

   do_217diffLin   = False
   do_143diffLin   = False
   do_100diffLin   = False

;   do_217cleanDiff = False
;   do_143cleanDiff = False
;   do_100cleanDiff = False

   do_217cleanHighEll = False
   do_217highEll      = False

   do_217cleanDiffLin = True
   do_143cleanDiffLin = False
   do_100cleanDiffLin = False

   do_217residuals = False
   do_217maskRes   = False

   scaling = 0.0081
   scaling = scaling / (1.-scaling)

   mollview, findgen(12), win=0, px=450
   loadct, 39
   !p.color=0
   !p.background=255

   spawn, 'ls /global/scratch2/sd/dpietrob/Software/XFaster/data/mask/dx11c/combined_*.fits', maskf
   nmask = n_elements(maskf)
;##    for i=0,nmask-1 do mollview, maskf[i], px=600, win=10+i

   if do_545spectra then begin

;       fwhm = 10.
;       arcmin2rad = !pi / 180. / 60.
;       fwhm2sig = 1./ sqrt(8.*alog(2.)) 
;       sig = fwhm  * arcmin2rad * fwhm2sig 
;       gb = gaussbeam( fwhm,6000 )
       ell = findgen(6001)
       Lc = 500.
       sig = 100.
;## ;       f = exp(-0.5*(ell-Lc)*(ell-Lc+1.) * sig^2)
       f = exp(-0.5*(ell-Lc)^2 / sig^2)
       f[where(ell lt Lc)] = 1.
       f = 1.-f
;       plot, ell, f, chars=1.5
;       oplot, ell, f, line=2, col=245
       mf = 1.-f
;       cl2fits, f, '/global/homes/d/dpietrob/myscratch/Software/XFaster/uscomp/xfaster/high-pass_filter.fits'
;       cl2fits, 1.-f, '/global/homes/d/dpietrob/myscratch/Software/XFaster/uscomp/xfaster/low-pass_filter.fits'
;stop
       if False then begin
           print, 'high-pass year-1...'
           ismoothing, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year1.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year1_high-ell.fits', beam_file='/global/homes/d/dpietrob/myscratch/Software/XFaster/uscomp/xfaster/high-pass_filter.fits', simul_type=1, nlmax=6000
           print, 'low-pass year-1...'
           ismoothing, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year1.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year1_low-ell.fits', beam_file='/global/homes/d/dpietrob/myscratch/Software/XFaster/uscomp/xfaster/low-pass_filter.fits', simul_type=1, nlmax=6000
           mymoll, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year1_low-ell.fits', /ash10
           mymoll, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year1_high-ell.fits', /ash10, win=1
;
           print, 'high-pass year-2...'
           ismoothing, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year2.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year2_high-ell.fits', beam_file='/global/homes/d/dpietrob/myscratch/Software/XFaster/uscomp/xfaster/high-pass_filter.fits', simul_type=1, nlmax=6000
;
           print, 'low-pass year-2...'
           ismoothing, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year2.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year2_low-ell.fits', beam_file='/global/homes/d/dpietrob/myscratch/Software/XFaster/uscomp/xfaster/low-pass_filter.fits', simul_type=1, nlmax=6000
;
           print, 'x-anafast high...'
           ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year1_high-ell.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_high-ell_xyr_cls.fits', map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year2_high-ell.fits'
;
           print, 'x-anafast low...'
           ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year1_low-ell.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_low-ell_xyr_cls.fits', map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year2_low-ell.fits'
;
           print, 'x-anafast high PSmask...'
           ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year1_high-ell.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_high-ell_xyr_cls_PSmask.fits', map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year2_high-ell.fits', maskfile='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/masks/mask_ptsrc_apo15_100_143_217.fits', regression=2
           print, 'x-anafast high COmask...'
           ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year1_high-ell.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_high-ell_xyr_cls_COmask.fits', map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year2_high-ell.fits', maskfile='masks/mask_CO3.fits', regression=2
;
           print, 'full PSmask...'
           ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year1.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_xyr_cls_PSmask.fits', map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year2.fits', maskfile='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/masks/mask_ptsrc_apo15_100_143_217.fits', regression=2
       endif

       fits2cl, cls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_high-ell_xyr_cls.fits'
       readcol,'beams/dx11c_wl_545.dat', il, wl
       hpw = healpixwindow(2048)
       loadct, 39
       !p.color=0
       !p.background=255
       window, 0, xsize=800, ysize=800
       ll = ell*(ell+1.)/2./!pi * 0 +1.

       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT'
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT'
       bell = bp_binning(ell,binsf)

       cls = cls / wl^2 / hpw^2 * scaling^2
       plot_oo, ell, cls*ll*0, chars=2, xr=[10,30000], yr=[1.e-8,10000], ys=1, xs=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6DX11c 545 GHz spectrum (rescaled by 0.0081)'

;##       xyouts, 1500, 2, '!6Full-Sky', chars=1.5

;##       fits2cl, pscls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_xyr_cls_PSmask.fits'
;##       mpscls = deconvolve_kernel(pscls, inv_fkernel='xfaster/mask_ptsrc_apo15_100_143_217_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_xyr_Mlike-cls_PSmask.fits')
       fits2cl, mpscls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_xyr_Mlike-cls_PSmask.fits'
       mpscls = mpscls / wl^2 / hpw^2 * scaling^2
       oplot, bell, bp_binning(mpscls * ll,binsf)
;
;##       fits2cl, pscls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_high-ell_xyr_cls_PSmask.fits'
;##       mpscls = deconvolve_kernel(pscls, fkernel='xfaster/mask_ptsrc_apo15_100_143_217_x_mask_ptsrc_apo15_100_143_217_kernel_l4000_v1.95.fits', write_kern='xfaster/mask_ptsrc_apo15_100_143_217_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_high-ell_xyr_Mlike-cls_PSmask.fits')
;       fits2cl, mpscls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_high-ell_xyr_Mlike-cls_PSmask.fits'
;       mpscls = mpscls / wl^2 / hpw^2 * scaling^2
;       oplot, bell, bp_binning(mpscls * ll,binsf)
;       xyouts, 1500, .1, '!6Ptsrc mask', chars=1.5
;       L1 = 600
;       L2 = 3800
;       x = alog( ell[L1:L2]/ell[3000] )
;       y = alog( mpscls[L1:L2] )
;       pars = linfit(x, y, yfit=fit, sigma=pars_err)
;       oplot, ell[L1:L2], exp(fit)*ll[L1:L2], col=245
;       xyouts, 450, 20, '!8l!6=['+strtrim(string(L1),2)+','+strtrim(string(L2),2)+']', col=245, chars=2
;       xyouts, 450, 8, '!6(A!d3000!n,n)=('+string(exp(pars[0]),format='(f6.3)')+','+string(pars[1],format='(f6.3)')+')', col=245, chars=2

       A = [ 5, -2.53, 1.e-5, 0.65 ]
       A3 = [ 1., -4.2, 5, -2.50, 0.08, -0.9 ]
       n1 = 0.
       n2 = 0.
       n3 = 0.
       n3w = 0.
       Amp3 = 0.
       w = 0.
       cnt = 0
       for i=0,nmask-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, ' - --> '+maskf[i]
           if False then begin
; HIGH-ELL
               ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year1_high-ell.fits', $
                 '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_high-ell_xyr_cls_fsky'+fsky+'.fits', $
                 map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year2_high-ell.fits', maskfile=maskf[i], regression=2
               fits2cl, cls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_high-ell_xyr_cls_fsky'+fsky+'.fits'
               mcls = deconvolve_kernel( cls, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_high-ell_xyr_Mlike-cls_fsky'+fsky+'.fits')
; LOW-ELL
               ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year1_low-ell.fits', $
                 '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_low-ell_xyr_cls_fsky'+fsky+'.fits', $
                 map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_year2_low-ell.fits', maskfile=maskf[i], regression=2
               fits2cl, cls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_low-ell_xyr_cls_fsky'+fsky+'.fits'
;                   oplot, ell, cls / wl^2 / hpw^2 * scaling^2, col=40*(i)
               mcls = deconvolve_kernel( cls, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_low-ell_xyr_Mlike-cls_fsky'+fsky+'.fits')
; FULL
               ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year1.fits', $
                     '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_xyr_cls_fsky'+fsky+'.fits', regression=2, $
                 maskfile=maskf[i], map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year2.fits'
               fits2cl, cls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_xyr_cls_fsky'+fsky+'.fits'
               mcls = deconvolve_kernel( cls, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_xyr_Mlike-cls_fsky'+fsky+'.fits')                   
           endif
;##               fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_high-ell_xyr_Mlike-cls_fsky'+fsky+'.fits'
;##               mcls = mcls / wl^2 / hpw^2 * scaling^2
;##               oplot, bell, bp_binning(mcls * ll,binsf), col=40*(i)
               ;oplot, ell, mcls*ll, col=40*(i), ns=5
;
;##               fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_low-ell_xyr_Mlike-cls_fsky'+fsky+'.fits'
;##               mcls = mcls / wl^2 / hpw^2 * scaling^2
;##               oplot, bell, bp_binning(mcls * ll * mf,binsf), col=40*(i)
;               oplot, ell, mcls*ll, col=40*(i), ns=5
;
           fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_xyr_Mlike-cls_fsky'+fsky+'.fits'
           mcls = mcls / wl^2 / hpw^2 * scaling^2
           bmcls = bp_binning(mcls,binsf)
;##           oplot, bell, bmcls, col=40*(i)

           L1 = 50
           print, min(abs(bell-L1), imn), imn
           L2 = 3000
           print, min(abs(bell-L2), imx), imx
;##               x = alog( ell[L1:L2]/ell[3000] )
           x = bell[imn:imx]/3000.
;##               y = alog( mcls[L1:L2] )
           y = bmcls[imn:imx]
           err = sqrt(2./(2.*x*3000.+1.))
           indx = where( bmcls[imn:imx] gt 0.)
;##               help, indx
;##               print, imx-imn+1
           x = alog( x[indx] )
           y = alog( y[indx] )
           err = err[indx]

           oplot, exp(x)*3000., exp(y), psym=4, col=40*i
           errplot, exp(x)*3000., exp(y)*(1.-err), exp(y)*(1.+err), col=40*i;, psym=4
;stop
;               A = [ 9.e2, -2.25, .1, -0.95 ]
;               double_power_law, x, A, chk, xxx
;               oplot, x, chk, thick=2
;               print,  A
;               err = sqrt(2./(2.*bell+1.))*abs(mcls)
;               range = L1+lindgen(L2-L1+1)
;               print, minmax(range)
;               weights = 1./err[range[indx]]^2 * 0. + 1.
;               errplot, ell, mcls-err, mcls+err
;               errplot, x, y-err, y+err
;##               oplot, x, y-err;, y+err
;##               oplot, x, y+err
           weights = 1./( err*y )^2 ;* 0. + 1.
;           amask = [1,1,1,1]
;           yfit = curvefit( x, y, weights, A, sigma, function_name='double_power_law', chisq=c2, itmax=1000, fita=amask)
;           print, A
;           print, sigma
;           print, c2
;               oplot, exp(x), exp(yfit), col=40*(i), line=2, thick=2

           amask3 = [1,1,1,1,1,1]
           yfit3 = curvefit( x, y, weights, A3, sigma3, function_name='triple_power_law', chisq=c2, itmax=1000000, fita=amask3, iter=it, tol=1.e-6)
           print, A3
           print, sigma3
           print, c2
           print, it
           oplot, exp(x)*3000., exp(yfit3), col=40*(i), line=5, thick=2

           leg_step = abs(alog(1.e-4)-alog(1.e-7))/8
;               xyouts, 1050, exp(alog(1.e-7)+i*leg_step), '!6F!dsky!n='+fsky+':(n1,n2)='+string(A[1],format='(f5.2)')+','+string(A[3],format='(f6.3)'), chars=1.5, col=40*(i)
;##               xyouts, 1050, exp(alog(1.e-7)+i*leg_step), '!6F!dsky!n='+fsky+':(A!d2!n,n!d2!n)='+string(A[2],format='(e10.4)')+','+string(A[3],format='(f6.3)'), chars=1.5, col=40*(i)
           xyouts, 550, exp(alog(5.e-0)+(i+1)*leg_step), '!6F!dsky!n='+fsky+':(n!d1!n,n!d2!n,n!d3!n)=('+string(A3[1],format='(f5.2)')+','+string(A3[3],format='(f5.2)')+','+string(A3[5],format='(f5.2)')+')', chars=1.5, col=40*(i)
           xyouts, 20, exp(alog(1.e-8)+(i+1)*leg_step), '!6(A!d1!n,A!d2!n,A!d3!n)=('+string(A3[0],format='(e9.2)')+','+string(A3[2],format='(e9.2)')+','+string(A3[4],format='(e9.2)')+')', chars=1.5, col=40*(i)

;stop
           n1 += A3[1] ;*(float(fsky))
           n2 += A3[3] ;*(float(fsky))
           n3 += A3[5] ;*(100-float(fsky))
           cnt += 1
           Amp3 += A3[4] * (100-float(fsky))
           n3w += A3[5] *(100-float(fsky))
           w += (100-float(fsky))
;           print, ' - Amp ratio -->', A3[0], A3[2]/A3[0], A3[4]/A3[0]
       endfor
       n1 /= cnt
       n2 /= cnt
       n3 /= cnt
       Amp3 /= w
       n3w /= w
       print, ' - Average spectral indices = ', n1, n2, n3
       xyouts, 800, 3, '!6<n!d1!n,n!d2!n,n!d3!n>=('+string(n1,format='(f5.2)')+','+string(n2,format='(f5.2)')+','+string(n3,format='(f5.2)')+')', chars=1.5
       xyouts, 800, 1, '!6<A!d3!n>!dw!n='+string(Amp3,format='(e9.2)')+';', chars=1.5
       xyouts, 5000, 1, '!6<n!d3!n>!dw!n='+string(n3w,format='(f5.2)'), chars=1.5
       if do_png then write_png,'fg545_3pl_L'+string(L1,format='(i2.2)')+'.png',tvrd(/true)
       stop
   endif   

   if do_545spectraHigh then begin

       ell = findgen(6001)

       fits2cl, cls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_high-ell_xyr_cls.fits'
       readcol,'beams/dx11c_wl_545.dat', il, wl
       hpw = healpixwindow(2048)
       loadct, 39
       !p.color=0
       !p.background=255
       !p.multi = 0
       window, 0, xsize=800, ysize=800
       ll = ell*(ell+1.)/2./!pi * 0 +1.

       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT'
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT'
       bell = bp_binning(ell,binsf)

       cls = cls / wl^2 / hpw^2 * scaling^2
       plot_oo, ell, cls*ll*0, chars=2, xr=[10,30000], yr=[1.e1,1.e5], ys=1, xs=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6DX11c 545 GHz spectrum (rescaled by 0.0081)'

       fits2cl, mpscls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_xyr_Mlike-cls_PSmask.fits'
       mpscls = mpscls / wl^2 / hpw^2 * scaling^2
;##       oplot, bell, bp_binning(mpscls * ll,binsf)*bell*(bell+1)/2./!pi
;
       A3 = [ 1., -4.18, 5, -2.39, 5.6e-5, -0.44 ]
       amask3 = [1,1,1,1,1,1]
       Amp3 = 0.
       Amp3w = 0.
       n1 = 0.
       n2 = 0.
       cnt = 0
       w = 0.
       for i=0,nmask-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, ' - --> '+maskf[i]
;
           fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_xyr_Mlike-cls_fsky'+fsky+'.fits'
           mcls = mcls / wl^2 / hpw^2 * scaling^2
           bmcls = bp_binning(mcls,binsf)
;##           oplot, bell, bmcls, col=40*(i)

           L1 = 30
           print, min(abs(bell-L1), imn), imn
           L2 = 3000
           print, min(abs(bell-L2), imx), imx
;##               x = alog( ell[L1:L2]/ell[3000] )
           x = bell[imn:imx]/3000.
;##               y = alog( mcls[L1:L2] )
           y = bmcls[imn:imx]
           err = sqrt(2./(2.*x*3000.+1.))
           indx = where( bmcls[imn:imx] gt 0.)
;##               help, indx
;##               print, imx-imn+1
           x = alog( x[indx] )
           y = alog( y[indx] )
           err = err[indx]

           oplot, exp(x)*3000., exp(y)*bell*(bell+1)/2./!pi, psym=4, col=40*i
           errplot, exp(x)*3000., exp(y)*(1.-err)*bell*(bell+1)/2./!pi, exp(y)*(1.+err)*bell*(bell+1)/2./!pi, col=40*i;, psym=4
;stop
           weights = 1./(err*abs(y))^2 ;* 0. + 1.

           yfit3 = curvefit( x, y, weights, A3, sigma3, function_name='triple_power_law', chisq=c2, itmax=100000, fita=amask3, tol=1.e-6, iter=it)
           print, A3
           print, sigma3
           print, c2
           print, it
           oplot, exp(x)*3000., exp(yfit3)*bell*(bell+1)/2./!pi, col=40*(i), line=5, thick=2

           leg_step = abs(alog(1.e4)-alog(1.e3))/8

;##           xyouts, 550, exp(alog(1.e3)+(i+1)*leg_step), '!6F!dsky!n='+fsky+':(n!d1!n,n!d2!n)=('+string(A3[1],format='(f5.2)')+','+string(A3[3],format='(f5.2)')+')', chars=1.5, col=40*(i)
;##           xyouts, 20, exp(alog(1.e3)+(i+1)*leg_step), '!6(A!d1!n,A!d2!n)=('+string(A3[0],format='(e9.2)')+','+string(A3[2],format='(e9.2)')+')', chars=1.5, col=40*(i)
           xyouts, 500, exp(alog(2500)+(i+1)*leg_step), '!6(A!d3!n,n!d3!n)=('+string(A3[4],format='(e9.2)')+','+string(A3[5],format='(f7.3)')+')', chars=1.5, col=40*(i)

;stop
           n1 += A3[1]
           n2 += A3[3]
;##           Amp3w += A3[4]*(100.-float(fsky))
           cnt += 1
           print , (100.-float(fsky))
           w += (100.-float(fsky))
;           print, ' - Amp ratio -->', A3[0], A3[2]/A3[0], A3[4]/A3[0]
       endfor
;##       Amp3 /= cnt
;##       Amp3w /= w
       n1 /= cnt
       n2 /= cnt
;##       print, ' - Average point source amplitude (3rd power law) = ', Amp3, Amp3w
;##       xyouts, 2000, exp(alog(1.e-8)+(i)*leg_step), '!6<A!d3!n> = ('+string(Amp3,format='(e9.2)')+')', chars=1.5
;##       xyouts, 2000, exp(alog(1.e-8)+(i-1)*leg_step), '!6<A!d3!n>!dw!n = ('+string(Amp3w,format='(e9.2)')+')', chars=1.5
       xyouts, 1000, 3, '!6<n!d1!n,n!d2!n!n>=('+string(n1,format='(f5.2)')+','+string(n2,format='(f5.2)')+')', chars=1.5
       if do_png then write_png,'fg545_3pl_highEll_L'+string(L1,format='(i2.2)')+'.png',tvrd(/true)
       stop
   endif   

   if do_545spectraFix then begin

       ell = findgen(6001)

       fits2cl, cls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_high-ell_xyr_cls.fits'
       readcol,'beams/dx11c_wl_545.dat', il, wl
       hpw = healpixwindow(2048)
       loadct, 39
       !p.color=0
       !p.background=255
       !p.multi=0
       window, 0, xsize=800, ysize=800
       ll = ell*(ell+1.)/2./!pi * 0 +1.

       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT'
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT'
       bell = bp_binning(ell,binsf)

       cls = cls / wl^2 / hpw^2 * scaling^2
       plot_oo, ell, cls*ll*0, chars=2, xr=[10,30000], yr=[1.e-8,10000], ys=1, xs=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6DX11c 545 GHz spectrum (rescaled by 0.0081)'

       fits2cl, mpscls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_xyr_Mlike-cls_PSmask.fits'
       mpscls = mpscls / wl^2 / hpw^2 * scaling^2
       oplot, bell, bp_binning(mpscls * ll,binsf)
;
       A3 = [ 1., -4.18, 5, -2.39, 5.6e-5, -0.44 ]
       amask3 = [1,1,1,1,0,0]
       Amp3 = 0.
       Amp3w = 0.
       n1 = 0.
       n2 = 0.
       cnt = 0
       w = 0.
       fg_fit = fltarr(nmask,6)
       for i=0,nmask-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, ' - --> '+maskf[i]
;
           fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_xyr_Mlike-cls_fsky'+fsky+'.fits'
           mcls = mcls / wl^2 / hpw^2 * scaling^2
           bmcls = bp_binning(mcls,binsf)
;##           oplot, bell, bmcls, col=40*(i)

           L1 = 50
           print, min(abs(bell-L1), imn), imn
           L2 = 3000
           print, min(abs(bell-L2), imx), imx
;##               x = alog( ell[L1:L2]/ell[3000] )
           x = bell[imn:imx]/3000.
;##               y = alog( mcls[L1:L2] )
           y = bmcls[imn:imx]
           err = sqrt(2./(2.*x*3000.+1.))
           indx = where( bmcls[imn:imx] gt 0.)
;##               help, indx
;##               print, imx-imn+1
           x = alog( x[indx] )
           y = alog( y[indx] )
           err = err[indx]

           oplot, exp(x)*3000., exp(y), psym=4, col=40*i
           errplot, exp(x)*3000., exp(y)*(1.-err), exp(y)*(1.+err), col=40*i;, psym=4
;stop
           weights = 1./(err*abs(y))^2 ;* 0. + 1.

           yfit3 = curvefit( x, y, weights, A3, sigma3, function_name='triple_power_law', chisq=c2, itmax=100000, fita=amask3, tol=1.e-6, iter=it)
           print, A3
           fg_fit[i,*] = A3
           print, sigma3
           print, c2
           print, it
           oplot, exp(x)*3000., exp(yfit3), col=40*(i), line=5, thick=2

           leg_step = abs(alog(1.e-4)-alog(1.e-7))/8

           xyouts, 550, exp(alog(5.e-0)+(i+1)*leg_step), '!6F!dsky!n='+fsky+':(n!d1!n,n!d2!n)=('+string(A3[1],format='(f5.2)')+','+string(A3[3],format='(f5.2)')+')', chars=1.5, col=40*(i)
           xyouts, 20, exp(alog(1.e-8)+(i+1)*leg_step), '!6(A!d1!n,A!d2!n)=('+string(A3[0],format='(e9.2)')+','+string(A3[2],format='(e9.2)')+')', chars=1.5, col=40*(i)

;stop
           n1 += A3[1]
           n2 += A3[3]
;##           Amp3w += A3[4]*(100.-float(fsky))
           cnt += 1
           print , (100.-float(fsky))
           w += (100.-float(fsky))
;           print, ' - Amp ratio -->', A3[0], A3[2]/A3[0], A3[4]/A3[0]
       endfor
;##       Amp3 /= cnt
;##       Amp3w /= w
       n1 /= cnt
       n2 /= cnt
;##       print, ' - Average point source amplitude (3rd power law) = ', Amp3, Amp3w
;##       xyouts, 2000, exp(alog(1.e-8)+(i)*leg_step), '!6<A!d3!n> = ('+string(Amp3,format='(e9.2)')+')', chars=1.5
;##       xyouts, 2000, exp(alog(1.e-8)+(i-1)*leg_step), '!6<A!d3!n>!dw!n = ('+string(Amp3w,format='(e9.2)')+')', chars=1.5
       xyouts, 1000, 3, '!6<n!d1!n,n!d2!n!n>=('+string(n1,format='(f5.2)')+','+string(n2,format='(f5.2)')+')', chars=1.5
       if do_png then write_png,'fg545_3pl_fix_L'+string(L1,format='(i2.2)')+'.png',tvrd(/true)
       stop
   endif   

   if do_545diff then begin

       fits2cl, cls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_xyr_cls_PSmask.fits'
       readcol,'beams/dx11c_wl_545.dat', il, wl
       hpw = healpixwindow(2048)
       loadct, 39
       !p.color=0
       !p.background=255
       window, 0, xsize=800, ysize=800
       ell = findgen( 4001 )
       ll = ell*(ell+1.)/2./!pi * 0 +1.
       cls = cls / wl^2 / hpw^2 * scaling^2
       plot_oo, ell, cls*ll*0, chars=2, xr=[10,30000], yr=[1.e-8,10000], ys=1, xs=1
;##       xyouts, 1500, 2, '!6Full-Sky', chars=1.5

       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT'
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT'
       bell = bp_binning(ell,binsf)
       cls_ar = fltarr(n_elements(cls),nmask)
       bcls_ar = fltarr(n_elements(bell),nmask)

       A = [ 5, -2.55, 1.e-5, -0.65 ]
       for i=0,nmask-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, maskf[i]
;
           fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_xyr_Mlike-cls_fsky'+fsky+'.fits'
           mcls = mcls / wl^2 / hpw^2 * scaling^2
           cls_ar[*,i] = mcls
           bmcls = bp_binning(mcls,binsf)
           bcls_ar[*,i] = bmcls
           oplot, bell, bmcls, col=40*(i)
       endfor

       L1 = 30
       L2 = 3000
       print, min(abs(bell-L1),imn), min(abs(bell-L2), imx)
       plot, ell, cls*0, chars=2, xr=[50,4000], yr=[1.e-8,10000], /xlog, /ylog, tit='545 spectrum differences', ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6'
       leg_step = (alog(1.e3)-alog(1.e0) )/6
       A = [1.e-6, -2.53, 1.e-4, -3]
       n1 = 0.
       n2 = 0.
       n1w = 0.
       n2w = 0.
       w = 0.
       cnt = 0
       SigAvg = A*0.
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (bcls_ar[*,i]-bcls_ar[*,0]), col=40*i
           x = alog(bell[imn:imx]/3000)
           y = (bcls_ar[*,i]-bcls_ar[*,0])
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(bcls_ar[imn:imx,i]^2+bcls_ar[imn:imx,0]^2)/y[imn:imx]
           y = y[imn:imx]
           indx = where(y gt 0.)
           x = x[indx]
           y = alog(y[indx])
           err = err[indx]
           err *= y
;           oplot, exp(x)*3000, exp(y), psym=5
           errplot, exp(x)*3000, exp(y)*(1.-err/y), exp(y)*(1.+err/y), col=40*i
;           pars = linfit( x, y, yfit=fit, /double, covar=pars_cov, sigma=pars_err, measure_errors=err)
;           print, pars
;           print, pars_err
;           oplot, exp(x)*3000, exp(fit), col=40*i, thick=2

           weights = 1./err^2
           amask = [1,1,1,1]
           yfit = curvefit( x, y, weights, A, sigma, function_name='double_power_law', chisq=c2, itmax=1000, fita=amask, iter=it, tol=1.e-6)
           print, A
           print, sigma
           SigAvg += sigma^2
           print, it
           print, ' ------ '
           oplot, exp(x)*3000, exp(yfit), col=40*i, thick=2, line=5
;##           xyouts, 200, exp(1.e0+leg_step*i), '!6F!dsky!n='+fsky+'-30: (A!d3000!n,n)=('+string(exp(pars[0]),format='(f10.7)')+','+string(pars[1],format='(f5.2)')+')', col=40*i, chars=1.5
           xyouts, 20, exp(alog(1.e-8)+i*leg_step), '!6(n!d1!n,n!d2!n)='+string(A[1],format='(f7.3)')+','+string(A[3],format='(f7.3)'), chars=1.5, col=40*(i)
;##               idx += pars[1]
           n1 += max([A[1],A[3]])
           n2 += min([A[1],A[3]])
           n1w += max([A[1],A[3]]) * (float(fsky)) 
           n2w += min([A[1],A[3]]) * (float(fsky)) 
           w += (float(fsky))
           cnt += 1
       endfor
       n1 /= cnt
       n2 /= cnt
       n1w /= w
       n2w /= w
       SigAvg /= cnt
       SigAvg = sqrt(SigAvg)
       print, SigAvg
       print, ' - Average spectral indices = ',n1, n2
       xyouts, 200, 1000, '!6<(n!d1!n,n!d2!n)>=('+string(n1,format='(f5.2)')+','+string(n2,format='(f6.3)')+')', chars=1.5
       xyouts, 200, 250, '!6<(n!d1!n,n!d2!n)>!dw!n=('+string(n1w,format='(f5.2)')+','+string(n2w,format='(f6.3)')+')', chars=1.5
       if do_png then write_png,'fg545diff_v2_'+string(L1,format='(i2.2)')+'.png', tvrd(/true)
stop
   endif   

   if do_545diffFix then begin

       fits2cl, cls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_xyr_cls_PSmask.fits'
       readcol,'beams/dx11c_wl_545.dat', il, wl
       hpw = healpixwindow(2048)
       loadct, 39
       !p.color=0
       !p.background=255
       window, 0, xsize=800, ysize=800
       ell = findgen( 4001 )
       ll = ell*(ell+1.)/2./!pi * 0 +1.
       cls = cls / wl^2 / hpw^2 * scaling^2
       plot_oo, ell, cls*ll*0, chars=2, xr=[10,30000], yr=[1.e-8,10000], ys=1, xs=1
;##       xyouts, 1500, 2, '!6Full-Sky', chars=1.5

       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT'
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT'
       bell = bp_binning(ell,binsf)
       cls_ar = fltarr(n_elements(cls),nmask)
       bcls_ar = fltarr(n_elements(bell),nmask)

       for i=0,nmask-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, maskf[i]
;
           fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_xyr_Mlike-cls_fsky'+fsky+'.fits'
           mcls = mcls / wl^2 / hpw^2 * scaling^2
           cls_ar[*,i] = mcls
           bmcls = bp_binning(mcls,binsf)
           bcls_ar[*,i] = bmcls
           oplot, bell, bmcls, col=40*(i)
       endfor

       L1 = 30
       L2 = 3000
       print, min(abs(bell-L1),imn), min(abs(bell-L2), imx)
       plot, ell, cls*0, chars=2, xr=[50,4000], yr=[1.e-8,10000], /xlog, /ylog, tit='545 spectrum differences: indices fixed to average', ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6'
       leg_step = (alog(1.e3)-alog(1.e0) )/6
       if L1 eq 30 then A = [1.e-6, -2.50, 1.e-10, -5.3] ; L1=30
       if L1 eq 50 then A = [1.e-6, -2.50, 1.e-10, -4.9] ; L1=50
       amask = [1,0,1,0]
       n1 = 0.
       n2 = 0.
       cnt = 0
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (bcls_ar[*,i]-bcls_ar[*,0]), col=40*i
           x = alog(bell[imn:imx]/3000)
           y = (bcls_ar[*,i]-bcls_ar[*,0])
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(bcls_ar[imn:imx,i]^2+bcls_ar[imn:imx,0]^2)/y[imn:imx]
           y = y[imn:imx]
           indx = where(y gt 0.)
           x = x[indx]
           y = alog(y[indx])
           err = err[indx]
           err *= y
;           oplot, exp(x)*3000, exp(y), psym=5
           errplot, exp(x)*3000, exp(y)*(1.-err/y), exp(y)*(1.+err/y), col=40*i

           weights = 1./err^2
           yfit = curvefit( x, y, weights, A, sigma, function_name='double_power_law', chisq=c2, itmax=500000, tol=1.e-6, fita=amask)
           print, A
           print, sigma
           print, c2
           print, ' ------ '
           oplot, exp(x)*3000, exp(yfit), col=40*i, thick=2, line=5
;##           xyouts, 200, exp(1.e0+leg_step*i), '!6F!dsky!n='+fsky+'-30: (A!d3000!n,n)=('+string(exp(pars[0]),format='(f10.7)')+','+string(pars[1],format='(f5.2)')+')', col=40*i, chars=1.5
           xyouts, 20, exp(alog(1.e-8)+i*leg_step), '!6(A!d1!n,A!d2!n)='+string(A[0],format='(e9.2)')+','+string(A[2],format='(e9.2)'), chars=1.5, col=40*(i)
;##               idx += pars[1]
;           n1 += max([A[1],A[3]])
;           n2 += min([A[1],A[3]])
;           cnt += 1
       endfor
;       n1 /= cnt
;       n2 /= cnt
;       print, ' - Average spectral indices = ',n1, n2
       xyouts, 200, 1000, '!6<(n1,n2)>=('+string(A[1],format='(f6.2)')+','+string(A[3],format='(f6.2)')+')', chars=1.5
       if do_png then write_png,'fg545diff_v2_fix.png', tvrd(/true)
stop
   endif   

;## ------ DOT REMEMBER
   if False then begin
;   read_fits_map, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year1.fits', m1, order=ord
;   m1 = m1[*,0]
;   tmp = m1
;   remove_dipole, tmp, fltarr(12l*2048l^2)+1., ordering=ord, nside=2048
;   write_fits_map, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545_templates_year1.fits', m1, /nest
;   read_fits_map, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year1.fits', m2, order=ord
;   m2 = m2[*,0]
;   remove_dipole, m2[*,0], ordering=ord, nside=2048
;   write_fits_map, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545_templates_year2.fits', m2, /nest
;   ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545_templates_year1.fits', $
;     '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_templates_xyr_cls.fits', map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545_templates_year2.fits';, iter=3

;##   ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year1.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_xyr_cls.fits', map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year2.fits';, iter=3

;##   mymoll, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year1.fits', file2='/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year2.fits', win=25
       fits2cl, clfs, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_xyr_cls.fits'
;##   fits2cl, clfs, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_545_templates_xyr_cls.fits'
       clfs = clfs[*,0]/wl^2/hpw^2*scaling^2
       bclfs = xf_binning( clfs,'/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT', lcen=l )
       mollview, findgen(12), win=0, px=800
       loadct, 39
       !p.color=0
       !p.background=255
       window, 0 
       plot, l, bclfs, chars=1.5, psym=-4, /ylog, /xlog, xr=[1,4000], yr=[1,1.e7]
;##   plot, l, bclfs/bclfs, chars=1.5, psym=-4, yr=[0.5, 2], /ylog, /xlog
       readcol,'beams/dx11c_wl_545.dat', il, wl
       hpw = healpixwindow(2048)

       oplot, il, clfs*(il*(il+1.)/2./!pi), line=2
;stop
;##   ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_SDR_545-353_templates_yrs.fits','/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrs_cls.fits';, iter=3
;##   fits2cl, xxx, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrs_cls.fits'
;##   oplot, xxx, col=245
;##   ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_SDR_545-353_templates_yrd.fits','/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrd_nls.fits';, iter=3
;##   fits2cl, yyy, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrd_nls.fits'
;##   oplot, yyy
;##   oplot, xxx-yyy, col=70
;stop
       if False then begin
           for i=0,nmask-1 do begin
               fsky = string(30+i*10,format='(i2.2)')
               print, maskf[i]
               if True then begin
;##               mymoll, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_SDR_545-353_templates_yrs.fits', /ash10, maskfile=maskf[i], win=20
                   if False then begin
                       ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_SDR_545-353_templates_yrs.fits','/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrs_cls_fsky'+fsky+'.fits', maskfile=maskf[i]
                       ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_SDR_545-353_templates_yrd.fits','/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrs_nls_fsky'+fsky+'.fits', maskfile=maskf[i]
                       fits2cl, cls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrs_cls_fsky'+fsky+'.fits'
                       fits2cl, nls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrs_nls_fsky'+fsky+'.fits'
                       oplot, (cls[*,0]-nls[*,0])/wl^2/hpw^2, col=40*(i+1)
                       icls = (cls[*,0]-nls[*,0])/wl^2/hpw^2
                   endif else begin
;##                   ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year1.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_xyr_cls_fsky_'+fsky+'.fits', map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year2.fits', maskfile=maskf[i], regression=2
                       fits2cl, icls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_xyr_cls_fsky_'+fsky+'.fits'
                   endelse
;##               mcls = deconvolve_kernel( (cls[*,0]-nls[*,0]), fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_x_combined_'+fsky+'_kernel_l4000_v1.95.fits', write_kern='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrs_Mlike-cls_fsky'+fsky+'.fits')
                   mcls = deconvolve_kernel( icls, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrs_Mlike-cls_fsky'+fsky+'.fits')
                   fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrs_Mlike-cls_fsky'+fsky+'.fits'
                   mcls = mcls / wl^2 / hpw^2 * scaling^2
                   bmcls = xf_binning(mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT')
                   oplot, l, bmcls, col=40*(i+1), psym=-5
               endif
               fits2cl, cl, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrs_Mlike-cls_fsky'+fsky+'.fits'
               bcl = xf_binning(cl,'/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT')
               print, total(bcl), total(bclfs)
               wset, 0 & oplot, l, bcl/bclfs, col=40*(i+1)
           endfor
           stop
       endif
   endif   

;## ------ 100 RAW LINEAR
   if do_100diffLin then begin
       readcol,'beams/dx11c_wl_100.dat', il, wl100
       hpw = healpixwindow(2048)
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT'
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT'
       ell = findgen(4001)
       bell = bp_binning(ell, binsf)
       cls_ar = fltarr(n_elements(bell),nmask)
       !p.multi=[0,2,1]
       window, 0, xsize=1200, ysize=600
       plot, /nodata, [10,2500], [1.e2,1.e4], /xlog, chars=2, /ylog, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6DX11c 100 GHz yr1 x yr2', xs=1
       leg_step = abs( alog(1000)-alog(200))/6
       for i=0,nmask-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, maskf[i]
           if False then begin
               ianafast, '/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_100_year1_uK.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_100_xyr_cls_fsky_'+fsky+'.fits', map2_in='/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_100_year2_uK.fits', maskfile=maskf[i], regression=2
               fits2cl, icls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_100_xyr_cls_fsky_'+fsky+'.fits'
               mcls = deconvolve_kernel( icls, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_100_yrs_Mlike-cls_fsky'+fsky+'.fits')
           endif
           fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_100_yrs_Mlike-cls_fsky'+fsky+'.fits'
           mcls = mcls / wl100^2 / hpw^2
           bmcls = bp_binning(mcls, binsf)
           cls_ar[*,i] = bmcls
           oplot, bell, bmcls*bell*(bell+1)/2./!pi, col=40*(i)
           xyouts, 20, exp(alog(150)+leg_step*i), '!6F!dsky!n='+fsky, chars=1.5, col=40*(i)
       endfor

       L1 = 30
       L2 = 1500
       print, min(abs(bell-L1),imn), min(abs(bell-L2), imx)
       plot, /nodata, [10,2500], [-750, 4000], /xlog, chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0.
       amask = [1,0,1,1]
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
;           oplot, bell, -(cls_ar[*,i]-cls_ar[*,0]), col=40*i, line=2
           x = ( bell[imn:imx] )
           y = ( cls_ar[*,i]-cls_ar[*,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
           y = y[imn:imx]
           errplot, x, (y-err)*bell*(bell+1)/2./!pi, (y+err)*bell*(bell+1)/2./!pi, col=40*i
           weights = 1./err^2 ;* 0. + 1.
           A = [1.e-3, -2.5, 1., -4.2]

           yfit3 = curvefit( x, y, weights, A, sigma3, function_name='double_power_law_lin_norm', chisq=c2, itmax=1000, fita=amask, iter=it, tol=1.e-6, /double)
           print, A
;##           xyouts, 20, 2000+200*i, '!6F!dsky!n='+fsky+'-30:(n1,n2)=('+string(A[1],format='(f5.2)')+','+string(A[3],format='(f6.3)')+')', chars=1.5, col=40*(i)
           xyouts, 20, 3600, '!6(A!d1!n,n!d1!n,A!d2!n,n!d2!n)', chars=1.5
           xyouts, 70, 800+400*i, '!6+/-('+string(sigma3[0],format='(e9.2)')+','+string(sigma3[1],format='(f7.3)')+','+string(sigma3[2],format='(e9.2)')+','+string(sigma3[3],format='(f7.3)')+')', chars=1.125, col=40*(i)
           xyouts, 20, 1000+400*i, '!6F!dsky!n='+fsky+'-30:('+string(A[0],format='(e9.2)')+','+string(A[1],format='(f7.3)')+','+string(A[2],format='(e9.2)')+','+string(A[3],format='(f7.3)')+')', chars=1.25, col=40*(i)
           print, sigma3
           oplot, x, yfit3*bell*(bell+1)/2./!pi, col=40*i, thick=2, line=5
           print, ' ------ '
       endfor
       if do_png then write_png, 'fg100r_lin.png', tvrd(/true)

       !p.multi=0
       window, 1, xsize=600, ysize=600
       plot, /nodata, [10,2500], [-750, 1000], /xlog, chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0.
       amask = [1,1,1,1]
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
;           oplot, bell, -(cls_ar[*,i]-cls_ar[*,0]), col=40*i, line=2
           x = ( bell[imn:imx] )
           y = ( cls_ar[*,i]-cls_ar[*,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
           y = y[imn:imx]
           errplot, x, (y-err)*bell*(bell+1)/2./!pi, (y+err)*bell*(bell+1)/2./!pi, col=40*i
           weights = 1./err^2 ;* 0. + 1.
           A = [1.e-3, -2.5, 1., -4.2]

           yfit3 = curvefit( x, y, weights, A, sigma3, function_name='double_power_law_lin_norm', chisq=c2, itmax=1000, fita=amask, iter=it, tol=1.e-6, /double)
           oplot, x, yfit3*bell*(bell+1)/2./!pi, col=40*i, thick=2, line=5
       endfor
       if do_png then write_png, 'fg100r_lin_zoom.png', tvrd(/true)
       stop
   endif   

;## ------ 143 RAW LINEAR
   if do_143diffLin then begin
       readcol,'beams/dx11c_wl_143.dat', il, wl143
       hpw = healpixwindow(2048)
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT'
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT'
       ell = findgen(4001)
       bell = bp_binning(ell, binsf)
       cls_ar = fltarr(n_elements(bell),nmask)
       !p.multi=[0,2,1]
       window, 0, xsize=1200, ysize=600
       plot, /nodata, [10,3000], [1.e2,1.e4], /xlog, chars=2, /ylog, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6DX11c 143 GHz yr1 x yr2', xs=1
       leg_step = abs( alog(1000)-alog(200))/6
       for i=0,nmask-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, maskf[i]
           if False then begin
               ianafast, '/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_143_year1_uK.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_143_xyr_cls_fsky_'+fsky+'.fits', map2_in='/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_143_year2_uK.fits', maskfile=maskf[i], regression=2
               fits2cl, icls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_143_xyr_cls_fsky_'+fsky+'.fits'
               mcls = deconvolve_kernel( icls, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_143_yrs_Mlike-cls_fsky'+fsky+'.fits')
           endif
           fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_143_yrs_Mlike-cls_fsky'+fsky+'.fits'
           mcls = mcls / wl143^2 / hpw^2
           bmcls = bp_binning(mcls, binsf)
           cls_ar[*,i] = bmcls
           oplot, bell, bmcls*bell*(bell+1)/2./!pi, col=40*(i)
           xyouts, 20, exp(alog(150)+leg_step*i), '!6F!dsky!n='+fsky, chars=1.5, col=40*(i)
       endfor

       L1 = 30
       L2 = 2500
       print, min(abs(bell-L1),imn), min(abs(bell-L2), imx)
       amask = [1,0,1,1]
       plot, /nodata, [10,3000], [-750, 4000], /xlog, chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0.
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
;           oplot, bell, -(cls_ar[*,i]-cls_ar[*,0]), col=40*i, line=2
           x = ( bell[imn:imx] )
           y = ( cls_ar[*,i]-cls_ar[*,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
           y = y[imn:imx]
           errplot, x, (y-err)*bell*(bell+1)/2./!pi, (y+err)*bell*(bell+1)/2./!pi, col=40*i
           weights = 1./err^2 ;* 0. + 1.
           A = [1.e-5, -2.5, 1.e-5, -4.2]

           yfit3 = curvefit( x, y, weights, A, sigma3, function_name='double_power_law_lin_norm', chisq=c2, itmax=10000l, fita=amask, iter=it, tol=1.e-6, /double)
           print, A
           print, sigma3
           print, c2
           print, it
;##           xyouts, 20, 2000+200*i, '!6F!dsky!n='+fsky+'-30:(n1,n2)=('+string(A[1],format='(f5.2)')+','+string(A[3],format='(f6.3)')+')', chars=1.5, col=40*(i)
           xyouts, 40, 3600, '!6(A!d1!n,n!d1!n,A!d2!n,n!d2!n)', chars=1.5
           xyouts, 35, 800+400.*i, '!6      +/-('+string(sigma3[0],format='(e9.2)')+','+string(sigma3[1],format='(f5.2)')+','+string(sigma3[2],format='(e9.2)')+','+string(sigma3[3],format='(f6.3)')+')', chars=1.25, col=40*(i)
           xyouts, 35, 1000+400*i, '!6'+fsky+'-30:('+string(A[0],format='(e9.2)')+','+string(A[1],format='(f5.2)')+','+string(A[2],format='(e9.2)')+','+string(A[3],format='(f6.3)')+')', chars=1.5, col=40*(i)
           oplot, x, yfit3*bell*(bell+1)/2./!pi, col=40*i, thick=2, line=5
           print, ' ------ '
       endfor
       if do_png then write_png, 'fg143r_lin.png', tvrd(/true)

       !p.multi=0
       window, 1, xsize=600, ysize=600
       plot, /nodata, [10,3000], [-750, 1000], /xlog, chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0.
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
;           oplot, bell, -(cls_ar[*,i]-cls_ar[*,0]), col=40*i, line=2
           x = ( bell[imn:imx] )
           y = ( cls_ar[*,i]-cls_ar[*,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
           y = y[imn:imx]
           errplot, x, (y-err)*bell*(bell+1)/2./!pi, (y+err)*bell*(bell+1)/2./!pi, col=40*i
           weights = 1./err^2 ;* 0. + 1.

           A = [1.e-5, -2.5, 1.e-5, -4.2]
           yfit3 = curvefit( x, y, weights, A, sigma3, function_name='double_power_law_lin_norm', chisq=c2, itmax=1000, fita=amask, iter=it, tol=1.e-6, /double)
           oplot, x, yfit3*bell*(bell+1)/2./!pi, col=40*i, thick=2, line=5
       endfor
       if do_png then write_png, 'fg143r_lin_zoom.png', tvrd(/true)
       stop
   endif   

;## ------ 217 RAW LINEAR
   if do_217diffLin then begin
       readcol,'beams/dx11c_wl_217.dat', il, wl217
       hpw = healpixwindow(2048)
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT'
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT'
       ell = findgen(4001)
       bell = bp_binning(ell, binsf)
       cls_ar = fltarr(n_elements(bell),nmask)
       !p.multi=[0,2,1]
       window, 0, xsize=1200, ysize=600
       plot, /nodata, [10,3500], [1.e2,1.e4], /xlog, chars=2, /ylog, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6DX11c 217 GHz yr1 x yr2', xs=1
       leg_step = abs( alog(1000)-alog(200))/6
       clean_217spe = fltarr(n_elements(bell),nmask)
       for i=0,nmask-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, maskf[i]
           if False then begin
               ianafast, '/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_217_year1_uK.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_xyr_cls_fsky_'+fsky+'.fits', map2_in='/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_217_year2_uK.fits', maskfile=maskf[i], regression=2
               fits2cl, icls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_xyr_cls_fsky_'+fsky+'.fits'
               mcls = deconvolve_kernel( icls, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_yrs_Mlike-cls_fsky'+fsky+'.fits')
           endif
           fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_yrs_Mlike-cls_fsky'+fsky+'.fits'
           mcls = mcls / wl217^2 / hpw^2
           triple_power_law_lin, ell/3000., fg_fit[i,*], F, xxx
           bmcls = bp_binning(mcls, binsf)
           fg = bp_binning(F, binsf)
           cls_ar[*,i] = bmcls
           oplot, bell, bmcls*bell*(bell+1)/2./!pi, col=40*(i)
;##           oplot, bell, fg*bell*(bell+1)/2./!pi, col=40*(i), psym=3
;##           oplot, bell, (bmcls-fg)*bell*(bell+1)/2./!pi, col=40*(i)
;##           clean_217spe[*,i] = (bmcls-fg)*bell*(bell+1)/2./!pi
           xyouts, 20, exp(alog(150)+leg_step*i), '!6F!dsky!n='+fsky, chars=1.5, col=40*(i)
       endfor

;##stop
       L1 = 50
       L2 = 3000
       print, min(abs(bell-L1),imn), min(abs(bell-L2), imx)
;##       A = [1.e-3, -2.5, 1., -4.2]
       amask = [1,0,1,1]
       plot, /nodata, [10,3500], [-750, 40000], /xlog, chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0.
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
;           oplot, bell, -(cls_ar[*,i]-cls_ar[*,0]), col=40*i, line=2
           x = ( bell[imn:imx] )
           y = ( cls_ar[imn:imx,i]-cls_ar[imn:imx,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
;##           y = y[imn:imx]
           errplot, x, (y-err)*x*(x+1)/2./!pi, (y+err)*x*(x+1)/2./!pi, col=40*i
           weights = 1./err^2 ;* 0. + 1.

           A = [1.e-6, -2.5, 1.e-5, -4.2]
           yfit3 = curvefit( x, y, weights, A, sigma3, function_name='double_power_law_lin_norm', chisq=c2, itmax=50000l, fita=amask, iter=it, tol=1.e-7, /double)
           print, A
           print, sigma3
           print, it
           print, c2
;##           print, yfit3
;##           xyouts, 20, 2000+200*i, '!6F!dsky!n='+fsky+'-30:(n1,n2)=('+string(A[1],format='(f5.2)')+','+string(A[3],format='(f6.3)')+')', chars=1.5, col=40*(i)
           xyouts, 40, 37000, '!6(A!d1!n,n!d1!n,A!d2!n,n!d2!n)', chars=1.5
           xyouts, 35, 8000+4000.*i, '!6      +/-('+string(sigma3[0],format='(e9.2)')+','+string(sigma3[1],format='(f5.2)')+','+string(sigma3[2],format='(e9.2)')+','+string(sigma3[3],format='(f6.3)')+')', chars=1.25, col=40*(i)
           xyouts, 35, 10000+4000.*i, '!6'+fsky+'-30:('+string(A[0],format='(e9.2)')+','+string(A[1],format='(f5.2)')+','+string(A[2],format='(e9.2)')+','+string(A[3],format='(f6.3)')+')', chars=1.5, col=40*(i)
           oplot, x, yfit3*x*(x+1)/2./!pi, col=40*i, thick=2, line=5
           print, ' ------ '
       endfor
       if do_png then write_png, 'fg217r_lin.png', tvrd(/true)

       !p.multi=0
       window, 1, xsize=600, ysize=600
       plot, /nodata, [10,3500], [-750, 4000], /xlog, chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0.
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
;           oplot, bell, -(cls_ar[*,i]-cls_ar[*,0]), col=40*i, line=2
           x = ( bell[imn:imx] )
           y = ( cls_ar[imn:imx,i]-cls_ar[imn:imx,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
;##           y = y[imn:imx]
           errplot, x, (y-err)*x*(x+1)/2./!pi, (y+err)*x*(x+1)/2./!pi, col=40*i
           weights = 1./err^2 ;* 0. + 1.

           A = [1.e-10, -2.5, 1.e-8, -4.2]
           yfit3 = curvefit( x, y, weights, A, sigma3, function_name='double_power_law_lin_norm', chisq=c2, itmax=50000l, fita=amask, iter=it, tol=1.e-6, /double)
           oplot, x, yfit3*x*(x+1)/2./!pi, col=40*i, thick=2, line=5
       endfor
       if do_png then write_png, 'fg217r_lin_zoom.png', tvrd(/true)

       stop
   endif   

;## ------ 217 RAW HighEll LINEAR
   if do_217highEll then begin
       readcol,'beams/dx11c_wl_217.dat', il, wl217
       hpw = healpixwindow(2048)
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT'
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT'
       ell = findgen(4001)
       bell = bp_binning(ell, binsf)
       cls_ar = fltarr(n_elements(bell),nmask)
       !p.multi=[0,2,1]
       window, 0, xsize=1200, ysize=600
       plot, /nodata, [10,3500], [1.e2,1.e4], /xlog, chars=2, /ylog, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6DX11c 217 GHz yr1 x yr2', xs=1
       leg_step = abs( alog(1000)-alog(200))/6
       for i=0,nmask-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, maskf[i]
           if False then begin
               ianafast, '/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_217_year1_uK.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_xyr_cls_fsky_'+fsky+'.fits', map2_in='/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_217_year2_uK.fits', maskfile=maskf[i], regression=2
               fits2cl, icls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_xyr_cls_fsky_'+fsky+'.fits'
               mcls = deconvolve_kernel( icls, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_yrs_Mlike-cls_fsky'+fsky+'.fits')
           endif
           fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_yrs_Mlike-cls_fsky'+fsky+'.fits'
           mcls = mcls / wl217^2 / hpw^2
           bmcls = bp_binning(mcls, binsf)
           cls_ar[*,i] = bmcls
           oplot, bell, bmcls*bell*(bell+1)/2./!pi, col=40*(i)
           xyouts, 20, exp(alog(150)+leg_step*i), '!6F!dsky!n='+fsky, chars=1.5, col=40*(i)
       endfor

       L1 = 30
       L2 = 3000
       print, min(abs(bell-L1),imn), min(abs(bell-L2), imx)
;##       A = [1.e-3, -2.5, 1., -4.2]
       amask = [1,1,1,1]
       plot, /nodata, [10,3500], [-50, 350], chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0.
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
;           oplot, bell, -(cls_ar[*,i]-cls_ar[*,0]), col=40*i, line=2
           x = ( bell[imn:imx] )
           y = ( cls_ar[*,i]-cls_ar[*,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
           y = y[imn:imx]
;##           errplot, x, (y-err)*bell*(bell+1)/2./!pi, (y+err)*bell*(bell+1)/2./!pi, col=40*i, psym=3
           oplot, x, (y-err)*bell*(bell+1)/2./!pi, col=40*i, psym=3
           oplot, x, (y+err)*bell*(bell+1)/2./!pi, col=40*i, psym=3
           weights = 1./err^2 ;* 0. + 1.

           A = [1.e-6, -2.5, 1.e-5, -4.2]
           yfit3 = curvefit( x, y, weights, A, sigma3, function_name='double_power_law_lin_norm', chisq=c2, itmax=50000l, fita=amask, iter=it, tol=1.e-7, /double)
           print, A
           print, sigma3
           print, it
           print, c2
;##           print, yfit3
;##           xyouts, 20, 2000+200*i, '!6F!dsky!n='+fsky+'-30:(n1,n2)=('+string(A[1],format='(f5.2)')+','+string(A[3],format='(f6.3)')+')', chars=1.5, col=40*(i)
           xyouts, 40, 37000, '!6(A!d1!n,n!d1!n,A!d2!n,n!d2!n)', chars=1.5
           xyouts, 35, 8000+4000.*i, '!6      +/-('+string(sigma3[0],format='(e9.2)')+','+string(sigma3[1],format='(f5.2)')+','+string(sigma3[2],format='(e9.2)')+','+string(sigma3[3],format='(f6.3)')+')', chars=1.25, col=40*(i)
           xyouts, 35, 10000+4000.*i, '!6'+fsky+'-30:('+string(A[0],format='(e9.2)')+','+string(A[1],format='(f5.2)')+','+string(A[2],format='(e9.2)')+','+string(A[3],format='(f6.3)')+')', chars=1.5, col=40*(i)
;##           oplot, x, yfit3*bell*(bell+1)/2./!pi, col=40*i, thick=2, line=5
           oplot, x, (y-yfit3)*bell*(bell+1)/2./!pi, col=40*i, thick=2, line=5
           print, ' ------ '
       endfor
       if do_png then write_png, 'fg217r_highEll_lin.png', tvrd(/true)
       stop
   endif   

;## ------ 217 CLEAN High-ell
   if do_217cleanHighEll then begin
       readcol,'beams/dx11c_wl_217.dat', il, wl217
       hpw = healpixwindow(2048)
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT'
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT'
       ell = findgen(4001)
       bell = bp_binning(ell, binsf)
       cls_ar = fltarr(n_elements(bell),nmask)
       !p.multi=[0,2,1]
       window, 0, xsize=1200, ysize=600
       plot, /nodata, [10,3500], [100, 10000], /xlog, chars=2, ys=1, /ylog, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6DX11c 217 GHz Undusted yr1 x yr2', xs=1
       leg_step = abs( alog(1000)-alog(200))/6
       for i=0,nmask-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, maskf[i]
           if False then begin
               ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/maps/dx11c_SDR_yr_1-2_IQUmap_217_extMask_545_coTP_undusted_split_ns2048_uK_hr1.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_undusted_xyr_cls_fsky_'+fsky+'.fits', map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/maps/dx11c_SDR_yr_1-2_IQUmap_217_extMask_545_coTP_undusted_split_ns2048_uK_hr2.fits', maskfile=maskf[i], regression=2
               fits2cl, icls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_undusted_xyr_cls_fsky_'+fsky+'.fits'
               mcls = deconvolve_kernel( icls, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_undusted_yrs_Mlike-cls_fsky'+fsky+'.fits')
           endif
           fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_undusted_yrs_Mlike-cls_fsky'+fsky+'.fits'
           mcls = mcls / wl217^2 / hpw^2
           bmcls = bp_binning(mcls, binsf)
           cls_ar[*,i] = bmcls
           oplot, bell, bmcls*bell*(bell+1)/2./!pi, col=40*(i)
           xyouts, 20, exp(alog(200)+i*leg_step), '!6F!dsky!n='+fsky, chars=1.5, col=40*(i)
       endfor

       L1 = 30
       L2 = 3000
       print, min(abs(bell-L1),imn), min(abs(bell-L2), imx)
       plot, /nodata, [10,3500], [-50, 100], chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
           x = ( bell[imn:imx] )
           y = ( cls_ar[*,i]-cls_ar[*,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
           y = y[imn:imx]
;##           errplot, x, (y-err)*bell*(bell+1)/2./!pi, (y+err)*bell*(bell+1)/2./!pi, col=40*i
           errplot, x, (y-err)*bell*(bell+1)/2./!pi, col=40*i, psym=3
           errplot, x, (y+err)*bell*(bell+1)/2./!pi, col=40*i, psym=3
           
           weights = 1./(err)^2
           a = [1.e-10, -2.5, 1.e-8, -6.2]
           amask = [1,0,1,1]
           yfit = curvefit( x, y, weights, A, sigma, function_name='double_power_law_lin_norm', chisq=c2, itmax=1000, fita=amask, iter=it, tol=1.e-6, /double)
           print, A
           print, sigma
           print, it
           print, c2
           oplot, x, yfit*bell*(bell+1)/2./!pi, col=40*i, thick=2, line=2
           print, ' ------ '
           xyouts, 20, 3600, '!6(A!d1!n,n!d1!n,A!d2!n,n!d2!n)', chars=1.5
           xyouts, 35, 800+400.*i, '!6  +/-('+string(sigma[0],format='(e9.2)')+','+string(sigma[1],format='(f7.3)')+','+string(sigma[2],format='(e9.2)')+','+string(sigma[3],format='(f7.3)')+')', chars=1.25, col=40*(i)
           xyouts, 20, 1000+400*i, '!6'+fsky+'-30:('+string(A[0],format='(e9.2)')+','+string(A[1],format='(f7.3)')+','+string(A[2],format='(e9.2)')+','+string(A[3],format='(f7.3)')+')', chars=1.5, col=40*(i)
;           xyouts, 20, 800+150*(nmask), '!6(A!d1!n,n!d1!n,A!d2!n,n!d2!n)', chars=1.5
;           xyouts, 20, 800+150*i, '!6'+fsky+'-30:('+string(A[0],format='(e9.2)')+','+string(A[1],format='(f5.2)')+','+string(A[2],format='(e9.2)')+','+string(A[3],format='(f6.3)')+')', chars=1.5, col=40*(i)
;##               idx += pars[1]
;           idx += A[1]
;           cnt += 1
       endfor
       if do_png then write_png, 'fg217u_highEll_lin.png', tvrd(/true)
       stop
   endif   

;## ------ 217 CLEAN
   if do_217cleanDiffLin then begin
       readcol,'beams/dx11c_wl_217.dat', il, wl217
       hpw = healpixwindow(2048)
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT'
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT'
       ell = findgen(4001)
       bell = bp_binning(ell, binsf)
       cls_ar = fltarr(n_elements(bell),nmask)
       !p.multi=[0,2,1]
       window, 0, xsize=1200, ysize=600
       plot, /nodata, [10,3500], [100, 10000], /xlog, chars=2, ys=1, /ylog, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6DX11c 217 GHz Undusted yr1 x yr2', xs=1
       leg_step = abs( alog(1000)-alog(200))/6
       for i=0,nmask-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, maskf[i]
           if False then begin
               ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/maps/dx11c_SDR_yr_1-2_IQUmap_217_extMask_545_coTP_undusted_split_ns2048_uK_hr1.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_undusted_xyr_cls_fsky_'+fsky+'.fits', map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/maps/dx11c_SDR_yr_1-2_IQUmap_217_extMask_545_coTP_undusted_split_ns2048_uK_hr2.fits', maskfile=maskf[i], regression=2
               fits2cl, icls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_undusted_xyr_cls_fsky_'+fsky+'.fits'
               mcls = deconvolve_kernel( icls, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_undusted_yrs_Mlike-cls_fsky'+fsky+'.fits')
           endif
           fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_undusted_yrs_Mlike-cls_fsky'+fsky+'.fits'
           mcls = mcls / wl217^2 / hpw^2
           bmcls = bp_binning(mcls, binsf)
           cls_ar[*,i] = bmcls
           oplot, bell, bmcls*bell*(bell+1)/2./!pi, col=40*(i)
;##           oplot, bell, clean_217spe[*,i], col=40*(i)
           xyouts, 20, exp(alog(200)+i*leg_step), '!6F!dsky!n='+fsky, chars=1.5, col=40*(i)
       endfor

;##stop
       L1 = 50
       L2 = 3000
       print, min(abs(bell-L1),imn), min(abs(bell-L2), imx)
       plot, /nodata, [10,3500], [-750, 4000], /xlog, chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
           x = ( bell[imn:imx] )
           y = ( cls_ar[imn:imx,i]-cls_ar[imn:imx,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
;##           y = y[imn:imx]
           errplot, x, (y-err)*x*(x+1)/2./!pi, (y+err)*x*(x+1)/2./!pi, col=40*i
           
           weights = 1./(err)^2
;##           a = [1.e-10, -2.5, 1.e-8, -6.2]
;##           amask = [1,0,1,1]
           a = [1.e-10, -2.5, 1.e-20, 0]
           amask = [1,1,0,0]
           yfit = curvefit( x, y, weights, A, sigma, function_name='double_power_law_lin_norm', chisq=c2, itmax=1000, fita=amask, iter=it, tol=1.e-6, /double)
           print, A
           print, sigma
           print, it
           print, c2
           oplot, x, yfit*x*(x+1)/2./!pi, col=40*i, thick=2, line=2
           print, ' ------ '
           xyouts, 20, 3600, '!6(A!d1!n,n!d1!n,A!d2!n,n!d2!n)', chars=1.5
           xyouts, 35, 800+400.*i, '!6  +/-('+string(sigma[0],format='(e9.2)')+','+string(sigma[1],format='(f7.3)')+','+string(sigma[2],format='(e9.2)')+','+string(sigma[3],format='(f7.3)')+')', chars=1.25, col=40*(i)
           xyouts, 20, 1000+400*i, '!6'+fsky+'-30:('+string(A[0],format='(e9.2)')+','+string(A[1],format='(f7.3)')+','+string(A[2],format='(e9.2)')+','+string(A[3],format='(f7.3)')+')', chars=1.5, col=40*(i)
;           xyouts, 20, 800+150*(nmask), '!6(A!d1!n,n!d1!n,A!d2!n,n!d2!n)', chars=1.5
;           xyouts, 20, 800+150*i, '!6'+fsky+'-30:('+string(A[0],format='(e9.2)')+','+string(A[1],format='(f5.2)')+','+string(A[2],format='(e9.2)')+','+string(A[3],format='(f6.3)')+')', chars=1.5, col=40*(i)
;##               idx += pars[1]
;           idx += A[1]
;           cnt += 1
       endfor
       if do_png then write_png, 'fg217u_lin.png', tvrd(/true)

       !p.multi=0
       window, 1, xsize=600, ysize=600
       plot, /nodata, [10,3500], [-750, 1000], /xlog, chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
           x = ( bell[imn:imx] )
           y = ( cls_ar[imn:imx,i]-cls_ar[imn:imx,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
;##           y = y[imn:imx]
           errplot, x, (y-err)*x*(x+1)/2./!pi, (y+err)*x*(x+1)/2./!pi, col=40*i
           
           weights = 1./(err)^2
;##           a = [1.e-10, -2.5, 1., -6.2]
           yfit = curvefit( x, y, weights, A, sigma, function_name='double_power_law_lin_norm', chisq=c2, itmax=1000, fita=amask, iter=it, tol=1.e-6, /double)
           oplot, x, yfit*x*(x+1)/2./!pi, col=40*i, thick=2, line=2
       endfor
       if do_png then write_png, 'fg217u_lin_zoom.png', tvrd(/true)
stop
   endif   


;## ------ 143 CLEAN
   if do_143cleanDiffLin then begin
       readcol,'beams/dx11c_wl_143.dat', il, wl143
       hpw = healpixwindow(2048)
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT'
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT'
       ell = findgen(4001)
       bell = bp_binning(ell, binsf)
       cls_ar = fltarr(n_elements(bell),nmask)
       !p.multi=[0,2,1]
       window, 0, xsize=1200, ysize=600
       plot, /nodata, [10,3000], [100, 10000], /xlog, chars=2, ys=1, /ylog, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6DX11c 143 GHz Undusted yr1 x yr2', xs=1
       leg_step = abs( alog(1000)-alog(200))/6
       for i=0,nmask-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, maskf[i]
           if False then begin
               ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/maps/dx11c_SDR_yr_1-2_IQUmap_143_extMask_545_coTP_undusted_split_ns2048_uK_hr1.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_143_undusted_xyr_cls_fsky_'+fsky+'.fits', map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/maps/dx11c_SDR_yr_1-2_IQUmap_143_extMask_545_coTP_undusted_split_ns2048_uK_hr2.fits', maskfile=maskf[i], regression=2
               fits2cl, icls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_143_undusted_xyr_cls_fsky_'+fsky+'.fits'
               mcls = deconvolve_kernel( icls, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_143_undusted_yrs_Mlike-cls_fsky'+fsky+'.fits')
           endif
           fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_143_undusted_yrs_Mlike-cls_fsky'+fsky+'.fits'
           mcls = mcls / wl143^2 / hpw^2
           bmcls = bp_binning(mcls, binsf)
           cls_ar[*,i] = bmcls
           oplot, bell, bmcls*bell*(bell+1)/2./!pi, col=40*(i)
           xyouts, 20, exp(alog(200)+i*leg_step), '!6F!dsky!n='+fsky, chars=1.5, col=40*(i)
       endfor

       L1 = 30
       L2 = 2500
       print, min(abs(bell-L1),imn), min(abs(bell-L2), imx)
       amask = [1,0,1,1]
       plot, /nodata, [10,3000], [-750, 4000], /xlog, chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
           x = ( bell[imn:imx] )
           y = ( cls_ar[*,i]-cls_ar[*,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
           y = y[imn:imx]
           errplot, x, (y-err)*bell*(bell+1)/2./!pi, (y+err)*bell*(bell+1)/2./!pi, col=40*i
           
           weights = 1./(err)^2
           a = [1., -2.5, 1., -4.2]
           yfit = curvefit( x, y, weights, A, sigma, function_name='double_power_law_lin_norm', chisq=c2, itmax=1000, fita=amask, iter=it, tol=1.e-6)
           print, A
           print, sigma
           print, c2
           print, it
           oplot, x, yfit*bell*(bell+1)/2./!pi, col=40*i, thick=2, line=2
           print, ' ------ '
           xyouts, 20, 3600, '!6(A!d1!n,n!d1!n,A!d2!n,n!d2!n)', chars=1.5
           xyouts, 20, 800+400.*i, '!6      +/-('+string(sigma3[0],format='(e9.2)')+','+string(sigma3[1],format='(f5.2)')+','+string(sigma3[2],format='(e9.2)')+','+string(sigma3[3],format='(f6.3)')+')', chars=1.25, col=40*(i)
           xyouts, 20, 1000+400*i, '!6'+fsky+'-30:('+string(A[0],format='(e9.2)')+','+string(A[1],format='(f5.2)')+','+string(A[2],format='(e9.2)')+','+string(A[3],format='(f6.3)')+')', chars=1.5, col=40*(i)
       endfor
       if do_png then write_png, 'fg143u_lin.png', tvrd(/true)

       !p.multi=0
       window, 1, xsize=600, ysize=600
       plot, /nodata, [10,3000], [-750, 1000], /xlog, chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
           x = ( bell[imn:imx] )
           y = ( cls_ar[*,i]-cls_ar[*,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
           y = y[imn:imx]
           errplot, x, (y-err)*bell*(bell+1)/2./!pi, (y+err)*bell*(bell+1)/2./!pi, col=40*i
           
           weights = 1./(err)^2
           a = [1., -2.5, 1., -4.2]
           yfit = curvefit( x, y, weights, A, sigma, function_name='double_power_law_lin_norm', chisq=c2, itmax=1000, fita=amask, iter=it, tol=1.e-6)
           oplot, x, yfit*bell*(bell+1)/2./!pi, col=40*i, thick=2, line=2
       endfor
       if do_png then write_png, 'fg143u_lin_zoom.png', tvrd(/true)
stop
   endif   


;## ------ 100 CLEAN
   if do_100cleanDiffLin then begin
       readcol,'beams/dx11c_wl_100.dat', il, wl100
       hpw = healpixwindow(2048)
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT'
       binsf = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT'
       ell = findgen(4001)
       bell = bp_binning(ell, binsf)
       cls_ar = fltarr(n_elements(bell),nmask)
       !p.multi=[0,2,1]
       window, 0, xsize=1200, ysize=600
       plot, /nodata, [10,2500], [100, 10000], /xlog, chars=2, ys=1, /ylog, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6DX11c 100 GHz Undusted yr1 x yr2', xs=1
       leg_step = abs( alog(1000)-alog(200))/6
       for i=0,nmask-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, maskf[i]
           if False then begin
               ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/maps/dx11c_SDR_yr_1-2_IQUmap_100_extMask_545_coTP_undusted_split_ns2048_uK_hr1.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_100_undusted_xyr_cls_fsky_'+fsky+'.fits', map2_in='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/maps/dx11c_SDR_yr_1-2_IQUmap_100_extMask_545_coTP_undusted_split_ns2048_uK_hr2.fits', maskfile=maskf[i], regression=2
               fits2cl, icls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_100_undusted_xyr_cls_fsky_'+fsky+'.fits'
               mcls = deconvolve_kernel( icls, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_100_undusted_yrs_Mlike-cls_fsky'+fsky+'.fits')
           endif
           fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_100_undusted_yrs_Mlike-cls_fsky'+fsky+'.fits'
           mcls = mcls / wl100^2 / hpw^2
           bmcls = bp_binning(mcls, binsf)
           cls_ar[*,i] = bmcls
           oplot, bell, bmcls*bell*(bell+1)/2./!pi, col=40*(i)
           xyouts, 20, exp(alog(200)+i*leg_step), '!6F!dsky!n='+fsky, chars=1.5, col=40*(i)
       endfor

       L1 = 30
       L2 = 1500
       print, min(abs(bell-L1),imn), min(abs(bell-L2), imx)

       plot, /nodata, [10,2500], [-750, 4000], /xlog, chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0
       amask = [1,0,1,1]
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
           x = ( bell[imn:imx] )
           y = ( cls_ar[*,i]-cls_ar[*,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
           y = y[imn:imx]
           errplot, x, (y-err)*bell*(bell+1)/2./!pi, (y+err)*bell*(bell+1)/2./!pi, col=40*i
           
           weights = 1./(err)^2
           a = [1., -2.5, 1., -4.2]
           yfit = curvefit( x, y, weights, A, sigma, function_name='double_power_law_lin_norm', chisq=c2, itmax=1000, fita=amask, iter=it, tol=1.e-6, /double)
           print, A
           print, sigma
           print, c2
           print, it
           oplot, x, yfit*bell*(bell+1)/2./!pi, col=40*i, thick=2, line=2
           print, ' ------ '
           xyouts, 20, 3600, '!6(A!d1!n,n!d1!n,A!d2!n,n!d2!n)', chars=1.5
           xyouts, 70, 800+400*i, '!6+/-('+string(sigma3[0],format='(e9.2)')+','+string(sigma3[1],format='(f7.3)')+','+string(sigma3[2],format='(e9.2)')+','+string(sigma3[3],format='(f7.3)')+')', chars=1.125, col=40*(i)
           xyouts, 20, 1000+400*i, '!6F!dsky!n='+fsky+'-30:('+string(A[0],format='(e9.2)')+','+string(A[1],format='(f7.3)')+','+string(A[2],format='(e9.2)')+','+string(A[3],format='(f7.3)')+')', chars=1.25, col=40*(i)
       endfor
       if do_png then write_png, 'fg100u_lin.png', tvrd(/true)

       !p.multi=0
       window, 1, xsize=600, ysize=600
       plot, /nodata, [10,2500], [-750, 1000], /xlog, chars=2, ys=1, ytit='!8C!dl!n!6 [!7l!6K!u2!n]', xtit='!8l!6', tit='!6Mask difference fit', xs=1
       oplot, bell, bell*0
       amask = [1,0,1,1]
       for i=1,nmask-1 do begin
           fsky = string(30+10*i, format='(i2.2)')
           oplot, bell, (cls_ar[*,i]-cls_ar[*,0])*bell*(bell+1)/2./!pi, col=40*i
           x = ( bell[imn:imx] )
           y = ( cls_ar[*,i]-cls_ar[*,0] )
           err = sqrt(2./(2.*bell[imn:imx]+1))*sqrt(cls_ar[imn:imx,i]^2+cls_ar[imn:imx,0]^2)
           y = y[imn:imx]
           errplot, x, (y-err)*bell*(bell+1)/2./!pi, (y+err)*bell*(bell+1)/2./!pi, col=40*i
           
           weights = 1./(err)^2
           a = [1., -2.5, 1., -4.2]
           yfit = curvefit( x, y, weights, A, sigma, function_name='double_power_law_lin_norm', chisq=c2, itmax=1000, fita=amask, iter=it, tol=1.e-6, /double)
           oplot, x, yfit*bell*(bell+1)/2./!pi, col=40*i, thick=2, line=2
       endfor
       if do_png then write_png, 'fg100u_lin_zoom.png', tvrd(/true)
stop
   endif   


;## ------ 217 Raw spectra as function of the mask
   if do_217maskRes then begin
       readcol,'beams/dx11c_wl_217.dat', il, wl217
       hpw = healpixwindow(2048)
       spawn, 'ls /global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/*217*combined_*.newdat', files 
       print, files
       nfiles = n_elements(files)
       cl30 = extract_xfaster_newdat(files[0], lcen=lc, cler=er30, btcl=th)
;##       mymoll, '/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_217_year1_uK.fits', min=-300, max=300, win=3
;##       mymoll, '/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_217_year2_uK.fits', min=-300, max=300, win=4
;##       mymoll, '/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_217_year1_uK.fits', file2='/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_217_year2_uK.fits', win=5, min=-100, max=100
;stop
;##       ianafast, '/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_217_year1_uK.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_xyr_cls.fits', map2_in='/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_217_year2_uK.fits
       fits2cl, clfs, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_xyr_cls.fits'
       !p.multi=0
       if True then begin
           window, 0, xsize=800, ysize=800
           plot, il, clfs/wl217^2/hpw^2 * (il*(il+1)/2./!pi), chars=2, /ylog, yr=[1e1,1.e4], ys=1, /xlog, xr=[100,10000], tit='545 (x 0.0081) and 217 Raw maps'
;       oplot, lc, cl30,psym=-4
;       errplot, lc, cl30-er30, cl30+er30

           print, min(abs(lc-1000), imin)
           print, min(abs(lc-100), imin2)
       
           leg_step = (alog(9000.)-alog(1000))/8

           for i=0,nmask-1 do begin
               fsky = string(30+i*10,format='(i2.2)')
               print, maskf[i]
               if True then begin
;##                   ianafast, '/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_217_year1_uK.fits', '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_xyr_cls_fsky_'+fsky+'.fits', map2_in='/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_solarDipoleRemoved_IQUmap_217_year2_uK.fits', maskfile=maskf[i], regression=2
                   fits2cl, icls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_xyr_cls_fsky_'+fsky+'.fits'
                   mcls = deconvolve_kernel( icls, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_yrs_Mlike-cls_fsky'+fsky+'.fits')
               endif
               fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_yrs_Mlike-cls_fsky'+fsky+'.fits'
               mcls = mcls / wl217^2 / hpw^2
               bmcls = xf_binning(mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT')
               oplot, lc, bmcls, col=40*(i), thick=2
     
               fits2cl, cl, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrs_Mlike-cls_fsky'+fsky+'.fits'
               bcl = xf_binning(cl,'/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT')
               print, total(bcl), total(bclfs)
               oplot, lc, bcl, col=40*(i)
;               oplot, l, bmcls-bcl, col=30*(i+1), psym=-4

               cl = extract_xfaster_newdat(files[i], lcen=lc, cler=er)
;               oplot, lc, cl, col=30*(i+1)
;               errplot, lc, cl-er, cl+er, col=30*(i+1)

               y = alog( bcl[imin:n_elements(lc)-1] )
               x = alog( lc[imin:n_elements(lc)-1]/lc[n_elements(lc)-1] )
               pars = linfit( x, y, yfit=fit, /double, covar=pars_cov, sigma=pars_err, measure_errors=er[imin:n_elements(lc)-1])
               xyouts, 1750, exp(alog(1000)+i*leg_step), '!6f!dsky!n'+fsky+':(A,n)=('+string(exp(pars[0]),format='(f6.2)')+','+string(pars[1],format='(f6.3)')+')', col=40*(i), chars=1.5
;## Low-ell portion of the spectrum 
;               y = alog( bcl[imin2:imin-1] )
;               x = alog( lc[imin2:imin-1] )
;               pars = linfit( x, y, yfit=fit, /double, covar=pars_cov, sigma=pars_err, measure_errors=er[imin2:imin-1])
;               xyouts, 10050, exp(alog(1000)+i*leg_step), '!6(A,n)=('+string(exp(pars[0]),format='(f6.1)')+','+string(pars[1],format='(f8.3)')+')', col=30*(i+1), chars=1.5
           endfor
           xyouts, 500, 7000, '!6217 raw map', chars=2
           xyouts, 500, 20, '!6545 raw map', chars=2
stop
       endif

       if True then begin
           window, 1, xsize=800, ysize=800
           plot, lc, th/th, line=2, chars=2, yr=[0.8, 1.5], xr=[100,3000], xs=1, ys=1;clfs/wl217^2/hpw^2 * (il*(il+1)/2./!pi), chars=2, yr=[1e1,600], ys=1, xr=[1000,4000], xs=1, /ylog
           oplot, lc, cl30/th, psym=-4
           errplot, lc, (cl30-er30)/th, (cl30+er30)/th
       
           print, min(abs(cl-1000), imin)
           print, min(abs(cl-100), imin2)
           
           leg_step = (alog(9000.)-alog(1000))/8

           for i=0,nmask-1 do begin
               fsky = string(30+i*10,format='(i2.2)')
               print, maskf[i]
               fits2cl, mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_yrs_Mlike-cls_fsky'+fsky+'.fits'
               bmcls = xf_binning(mcls, '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT')
               oplot, l, bmcls/th, col=30*(i+1), line=2
;
               fits2cl, cl, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_templates_yrs_Mlike-cls_fsky'+fsky+'.fits'
               bcl = xf_binning(cl,'/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT')
               print, total(bcl), total(bclfs)
;               oplot, l, bcl, col=30*(i+1)
               oplot, l, (bmcls-bcl)/th, col=30*(i+1), psym=5

               cl = extract_xfaster_newdat(files[i], lcen=lc, cler=er)
               oplot, lc, cl/th, col=30*(i+1), psym=4
;               errplot, lc, cl-er, cl+er, col=30*(i+1)

;               y = alog( bcl[imin:n_elements(lc)-1] )
;               x = alog( lc[imin:n_elements(lc)-1] )
;               pars = linfit( x, y, yfit=fit, /double, covar=pars_cov, sigma=pars_err, measure_errors=er[imin:n_elements(lc)-1])
;               xyouts, 2050, exp(alog(1000)+i*leg_step), '!6(A,n)=('+string(exp(pars[0]),format='(f6.1)')+','+string(pars[1],format='(f8.3)')+')', col=30*(i+1), chars=1.5
;## Low-ell portion of the spectrum 
;               y = alog( bcl[imin2:imin-1] )
;               x = alog( lc[imin2:imin-1] )
;               pars = linfit( x, y, yfit=fit, /double, covar=pars_cov, sigma=pars_err, measure_errors=er[imin2:imin-1])
;               xyouts, 10050, exp(alog(1000)+i*leg_step), '!6(A,n)=('+string(exp(pars[0]),format='(f6.1)')+','+string(pars[1],format='(f8.3)')+')', col=30*(i+1), chars=1.5

        endfor
    endif

       stop
   endif

   if do_217residuals then begin
;## ------ 217 GHz residuals
       spawn, 'ls /global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/*217*combined_*.newdat', files 
       print, files
       nfiles = n_elements(files)
       cl30 = extract_xfaster_newdat(files[0], lcen=lc, cler=er30)
       print, min(abs(lc-150),imin)
       print, min(abs(lc-2500),imax)
       onePowerLaw = False
       twoPowerLaws = True
       mollview, findgen(12), win=0, px=800
       loadct, 39
       !p.color=0
       !p.background=255
       !p.multi=[0,2,2]
       window, 0, xsize=1200, ysize=850
       plot, lc, cl30*0, chars=1.5, xtit='!8l', ytit='!7D!8D!dl!n', yr=[-70.0,170.0], ys=1, xr=[1,4500], xs=1
       lowest = 0.
       cls = fltarr(n_elements(lc), nfiles)
       Mcls = fltarr(n_elements(lc), nfiles)
       clser = fltarr(n_elements(lc), nfiles)
       cls[*,0] = cl30
       readcol,'beams/dx11c_wl_217.dat', l, wl
       hpw = healpixwindow(2048)
       for i=1,nfiles-1 do begin
           fsky = string(30+i*10,format='(i2.2)')
           print, fsky
           cl = extract_xfaster_newdat(files[i], cler=er)
           fits2cl, icls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_yr_1-2_IQUmap_217_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_cls_combined_'+fsky+'_x_combined_'+fsky+'.fits'
           fits2cl, inls, '/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_yr_1-2_IQUmap_217_extMask_545_coTP_undusted_split_ns2048_uK_hrhd_1_cls_combined_'+fsky+'_x_combined_'+fsky+'.fits'
;##       mcl = deconvolve_kernel( (icls[*,0]-inls[*,0])/hpw^2/wl^2, inv_fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_invTkernel_l4000.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_yrs_Mlike-cls_fsky'+fsky+'.fits')
;##       mcl = deconvolve_kernel( (icls[*,0]-inls[*,0])/hpw^2/wl^2, fkernel='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/combined_'+fsky+'_x_combined_'+fsky+'_kernel_l4000_v1.95.fits', write_cls='/global/scratch2/sd/dpietrob/Software/XFaster/uscomp/xfaster/dx11c_SDR_217_yrs_Mlike-cls_fsky'+fsky+'.fits')
;##       plot_oo, mcl
;##       bmcl = xf_binning( mcl,'/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT' )
;##       print, total(bmcl-cl)
           cls[*,i] = cl
           clser[*,i] = er
           oplot, lc[imin:imax], cl[imin:imax]-cl30[imin:imax], col=40*i, ns=8;, psym=-4
;##       oplot, lc, bmcl-cl30, col=40*i, psym=-5
;       errplot, lc, (cl-cl30)-er, (cl-cl30)+er,col=40*i
           lowest = min( [lowest,min(cl[imin:imax]-cl30[imin:imax]) ] )
           xyouts, 2850, 5+i*8, '!6'+string(30+i*10,format='(1i4)')+'% - 30% '+string(sqrt(mean((cl[imin:imax]-cl30[imin:imax])^2)),format='(f6.1)'), col=40*i, chars=1.5
;stop
       endfor

       res = fltarr(imax-imin+1, nfiles)

       plot, lc, cl30*0, chars=1.5, xtit='!8l', ytit='!7D!8D!dl!n', yr=[-40.0,40.0], ys=1, xr=[1,5000], xs=1
       for i=1,nfiles-1 do begin
           cl = cls[*,i]
           er = clser[*,i]
;       cl = extract_xfaster_newdat(files[i], cler=er)
; ------ Linear
           y = cl[imin:imax]-cl30[imin:imax]
           x = lc[imin:imax]/lc[imax]
; ------ Log
           y = alog(cl[imin:imax]-cl30[imin:imax]-lowest*1.1)
           x = alog( lc[imin:imax]/lc[imax] )
; ------ Linear
;           oplot, lc[imin:imax], y, col=40*i
; ------ Log
           oplot, lc[imin:imax], exp(y)+lowest*1.1, col=40*i
; ------ Single Power Law
           pars = linfit( x, y, yfit=fit, /double, covar=pars_cov, sigma=pars_err, measure_errors=sqrt(er[imin:imax]^2+er30[imin:imax]^2), chisq=chi2 ) 
           print, chi2
; Double Power Law
           A = [0.5, -0.5, 2, 0.4 ]
           weights = 1. / (er[imin:imax]^2+er30[imin:imax]^2)
           yfit = curvefit( x, y, weights, A, sigA, chisq=chi2, function_name='double_power_law',itmax=1000 ) 
           print, A
           print, sigA
           print, chi2
           print, ' ------ '
;##           est = exp( pars[0]+x*pars[1] )+lowest*1.1
; ------ Linear
;           res[*,i] = fit
; ------ Log
           if onePowerLaw then res[*,i] = exp(fit)+lowest*1.1
; Double power law
           if twoPowerLaws then res[*,i] = exp(yfit)+lowest*1.1
;       print, pars
; ------ Linear
           oplot, lc[imin:imax], res[*,i], col=40*i, thick=2
; ------ Single power Law
           if onePowerLaw then begin
               xyouts, 3050, -30+i*10, '!6(A,n)=('+string( pars[0],format='(f6.2)')+','+string(pars[1],format='(f6.2)')+')', col=40*i, chars=1.5
               xyouts, 3050, -35+i*10, '!6     +/-('+string( pars_err[0],format='(f6.2)')+','+string(pars_err[1],format='(f6.2)')+')', col=40*i, chars=1
           endif
; ------ Double Power Law
           if twoPowerLaws then begin
               xyouts, 3050, -30+i*10, '!6(n!d!n,n!d2!n)=('+string( A[1],format='(f6.2)')+','+string(A[3],format='(f6.2)')+')', col=40*i, chars=1.5
;##               xyouts, 3050, -35+i*10, '!6     +/-('+string( pars_err[0],format='(f6.2)')+','+string(pars_err[1],format='(f6.2)')+')', col=40*i, chars=1
           endif
;
;       xyouts, 3050, -30+i*10, '!6(A,n)=('+string(exp(pars[0]),format='(f6.2)')+','+string(pars[1],format='(f6.3)')+')', col=40*i, chars=1.5
;       xyouts, 3050, -35+i*10, '!6     +/-('+string(exp(pars[0])*pars_err[0]/pars[0],format='(f10.5)')+','+string(pars_err[1],format='(f6.3)')+')', col=40*i, chars=1
;       print, exp(pars[0])*pars_err[0]/pars[0]
;stop
       endfor

       plot, lc, cl30*0, chars=1.5, xtit='!8l', ytit='!7D!8D!dl!n', yr=[-70.0,70.0], ys=1, xr=[1,4500], xs=1
       for i=1,nfiles-1 do begin
           cl = cls[imin:imax,i]-cls[imin:imax,0]
           oplot, lc[imin:imax], cl-res[*,i], col=40*i
;       oplot, lc[imin:imax], res[*,i], col=40*i
           xyouts, 2850, 5+i*8, '!6'+string(30+i*10,format='(1i4)')+'% - 30% '+string(sqrt(mean((cl-res[*,i])^2)),format='(f6.1)'), col=40*i, chars=1.5
       endfor

       plot, lc, cl30*0, chars=1.5, xtit='!8l', ytit='!7D!8D!dl!n', yr=[-70.0,170.0], ys=1, xr=[1,4500], xs=1
       for i=1,nfiles-1 do begin
           cl = cls[imin:imax,i]-cls[imin:imax,0]
           oplot, lc[imin:imax], cl-res[*,i], col=40*i, ns=8
;       oplot, lc[imin:imax], res[*,i], col=40*i
           xyouts, 2850, 5+i*8, '!6'+string(30+i*10,format='(1i4)')+'% - 30% '+string(sqrt(mean((cl-res[*,i])^2)),format='(f6.1)'), col=40*i, chars=1.5
       endfor


stop
       plot, /nodata, [0,1], [0,1], col=255
       xyouts, -0.075,0.9,'Study of the foreground residuals in 217 clean map', chars=2.
       if onePowerLaw then xyouts, -0.05,0.75,'One power law', chars=1.5
       if twoPowerLaws then xyouts, -0.05,0.75,'Two power laws', chars=1.5
;       xyouts, -0.05,0.65,'of the DX11c 217 clean map for different masks,', chars=1.5
;       xyouts, -0.05,0.55,'differentiated against the largest one: 30% of the sky retained.', chars=1.5
;       xyouts, -0.05,0.45,'The top-right panel shows a power law fit to the residuals.', chars=1.5
;       xyouts, -0.05,0.35,'The bottom-left the residuals after fit subtraction.', chars=1.5
;       xyouts, -0.05,0.25,'The fit is performed in l=[300,2500]. The trend is possibly real', chars=1.5
;       xyouts, -0.05,0.15,'but the parameters are well within the error bars. Error bars', chars=1.5
;       xyouts, -0.05,0.05,'were computed assuming no error on the spectra. Accounting for the', chars=1.5
;       xyouts, -0.05,-0.05,'spectrum uncertainty makes the error on (A,n)=(40,5).', chars=1.5
       
       !p.multi = 0

   endif

   
end
