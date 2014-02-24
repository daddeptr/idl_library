 function td_sed, nu, emis, Td, nu0

    k_B     = 1.3806503d-23
    h       = 6.626068d-34
    c       = 2.99792458d8

    x = h / (k_B*Td)
    td=(exp(x*nu0)-1.d0) / (exp(x*nu)-1.d0) * (nu / nu0)^(emis+1.d0)

    return, td
 end
; ---
 function sd_sed, nu, alpha, b, nu0

    sd = exp( -((nu-nu0)/b)^2/2. + alpha )

    return, sd
end
; ---
function check_priors, x, xmn, xmx
   np = n_elements(x)
   passed = intarr(np)
   ok = 0b
   for ix=0,np-1 do if ( (x[ix] ge xmn[ix]) and (x[ix] le xmx[ix]) ) then passed[ix]=1
;print, x, xmn, xmx
;print, passed
   if (total(passed) eq np) then ok = 1b
;stop
   return, ok
end
; ---
function regress_map, map, a_err=a_err, do_res=do_res

   if not keyword_set(do_res) then do_res = 0b

   nfg = 3
   ns = 128l
   npix = 12l*ns^2
   temp = dblarr(npix,nfg)

;;    mollview, map, /hist, chars=1.5

   read_fits_map, 'lambda_fds_94GHz_uKantenna_60arcmin_ns128_ring.fits', fds
   read_fits_map, 'lambda_halpha_60arcmin_ns128_ring.fits', halpha
   read_fits_map, 'lambda_haslam_60arcmin_ns128_ring.fits', haslam

   read_fits_map, 'masks/gibbs_mask_10_ns128.fits', mask, nside=mns
   if (mns ne ns) then stop, 'Nside mismatch'
   gp = where(mask gt 0.)


   temp[*,0] = haslam
   temp[*,1] = halpha
   temp[*,2] = fds[*,0]

   iamp = 3

   print, ' - :DP Warning: rms = 1'
   rms = map*0 + 1.
   A = dblarr(nfg)
   oneOsigma = dblarr(nfg, nfg)
   sigma = oneOsigma
   a_err = dblarr(nfg)
       
   T = dblarr(nfg, nfg)
   B = dblarr(nfg)
   for ifg=0,nfg-1 do begin
       B[ifg] = correlate(map[gp]/ rms[gp],temp[gp,ifg]/ rms[gp], /covariance)
       for jfg=ifg,nfg-1 do begin
           T[ifg,jfg] = correlate( temp[gp,ifg]/ rms[gp],temp[gp,jfg]/ rms[gp], /covariance)
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

   if (do_res) then begin
       fit = map*0
       for ifg=0,nfg-1 do fit=fit+A[ifg]*temp[*,ifg]
       cc = map-fit
       bp = where(mask eq 0.)
       cc[bp] = -1.6375e30
       mollview, cc, chars=1.5, px=1000, /hist, tit='Regression Residuals', colt=3
   endif

   return, A
end
; ---
pro mc_sed, nstep

   !path=!path+':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'
   !path=!path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'     

   true  = 1b
   false = 0b                                                                                                                                                                                      
   mollview, randomn(-1,12), window=-1
   loadct, 39
   !p.color=0
   !p.background=255
   !p.multi=0

   ix=findgen(1001)
   nu = .408+ix*(900.)/1000
   nu = nu*.1e9 
   
   window,1, xsize=720, ysize=450
   readcol, '../../wmap/ns128/spin_dust_spectrum_v2.txt', wx, wy
   plot_oo, wx, wy, psym=1, xr=[10,1000], yr=[.1,1000]
   oplot, wx, wy, psym=1, col=230 
   
   wfreq=[22.8, 33., 41., 61., 94.]
   wsfreq = ['K', 'Ka', 'Q', 'V', 'W']
   wnfreq = n_elements(wfreq)

   pfreq=[30.,44.,70.,100.,143.,217.,353.,545.]
   psfreq = string(pfreq,format='(i3.3)')
   pnfreq = n_elements(pfreq)
   pfreq=[28.4,44.1,70.3,101.2,142.6,221.9,360.6,557.6]

   ns = 128l
   sns = strtrim(string(ns),2)
   npix = 12l*ns^2

   pindx = npix/2

   cut = make_sky_cut(30.,128l)
   cutp = where(cut gt 0.)

   do_compute_sed = false
   if (do_compute_sed) then begin
       dx = fltarr(13)
       dy = fltarr(13,3)
       de = fltarr(13,3)
       de[*,*] = 1.
       for ifreq=0,7 do begin
           dx[ifreq] = pfreq[ifreq]
           read_fits_map, 'nl_dx7_map_'+psfreq[ifreq]+'.fits', map
;; ;    dy[ifreq] = map[pindx]*conversionfactor(pfreq[ifreq], /thermo2antenna)
;;     dy[ifreq] = mean(map[cutp])*conversionfactor(pfreq[ifreq], /thermo2antenna)
        
           A = regress_map(map*conversionfactor(pfreq[ifreq], /thermo2antenna))
           dy[ifreq,*] = A
        
;;     read_fits_map, 'dx7_rms_'+psfreq[ifreq]+'.fits', map
;; ;    de[ifreq] = map[pindx]*conversionfactor(pfreq[ifreq], /thermo2antenna)
;;     de[ifreq] = sqrt(mean(map[cutp]^2))*conversionfactor(pfreq[ifreq], /thermo2antenna)
           de[ifreq,*] = 1.
       endfor

       for ifreq=0,4 do begin
           dx[8+ifreq] = wfreq[ifreq]
           read_fits_map, 'nl_wmap7_map_'+wsfreq[ifreq]+'.fits', map
;; ;    dy[8+ifreq] = map[pindx]*conversionfactor(wfreq[ifreq], /thermo2antenna)
;;     dy[8+ifreq] = mean(map[cutp])*conversionfactor(wfreq[ifreq], /thermo2antenna)

           A = regress_map(map*conversionfactor(pfreq[ifreq], /thermo2antenna))
           dy[8+ifreq,*] = A
;;     read_fits_map, 'wmap7_rms_'+wsfreq[ifreq]+'.fits', map
;; ;    de[8+ifreq] = map[pindx]*conversionfactor(wfreq[ifreq], /thermo2antenna)
;;     de[8+ifreq] = sqrt(mean(map[cutp]^2))*conversionfactor(wfreq[ifreq], /thermo2antenna)
           de[8+ifreq,*] = 1.
       endfor
       
       id = sort(dx)
       dx = dx[id]
       dy = dy[id,*]
       de = de[id,*]
       
       openw,1,'dust_SED.txt'
       for i=0,n_elements(dx)-1 do printf, 1, dx[i], dy[i,0], dy[i,1], dy[i,2]
       close,1
   endif
   spawn, 'head dust_SED.txt'
   
   readcol,'dust_SED.txt', dx, dy0, dy1, dy2
   
   dy = dy2
   de = dy2 * 0.+1.

   plot, dx, dy, psym=4, /xlog, /ylog, xr=[10,1000], yr=[.001, 1000], chars=1.5
   oplot, dx, dy, psym=4, col=245
   errplot, dx, dy-de, dy+de, col=245
; oplot, nu, (sd_sed(nu,-0.6,14, 20.)+td_sed(nu,1.8, 18, 545.) )*wy[5], col=210
;stop                                                                                                                                                                                      
; nstep = 2                                                                                                                                                                         
   p = [[.21, .0, -5, 5], $
        [12, 3., 5, 40], $
        [1.2, 0.3, 0.5, 2.2], $
        [18, 3., 5, 35]]       ;, $
;      [1, 0.1, 0., 50.]]                                                                                                                                                             

   np = n_elements(p[0,*])
          
   nu1 = 20.e9
   nu2 = 353.e9
   loglike = 0.
   oldlike = loglike
   oldp = p[0,*]
   accepted = 0l
   
   chp = fltarr(nstep,np+1)
   num=fltarr(nstep,np+1)
   chlike = fltarr(nstep)

   sy = (sd_sed(nu, p[0,0], p[0,1], nu1) + td_sed(nu, p[0,2], p[0,3], nu2) )
; pv = min(abs(nu-wx[5]),im)
; aa = wy[5]*sy[im]/sy[im]^2
; oplot, nu, sy * aa, line=2
   pv = min(abs(nu-dx[7]*1.e9),im)
   aa = dy[7]*sy[im]/sy[im]^2
   oplot, nu, sy * aa, line=2
 
   xdata = wx[1:*]
   ydata = wy[1:*]

   xdata = dx[*]
   ydata = dy[*]

   ndata = n_elements(xdata)

   stop

   for isamp=0, nstep-1 do begin
       bin = nstep / 10
     if ( (isamp/bin)*bin eq isamp) then print, isamp, ' / ', nstep, oldlike
     dp = randomn(-2*(isamp+1),np)
     num[isamp,0:np-1] = dp
     nwp = oldp + dp*p[1,*]
     if (isamp eq 0) then nwp = oldp
;     print, nwp
;     sed = sd_sed(nu,nwp[0], nwp[1], nu1) + td_sed(nu, nwp[2], nwp[3], nu2)
;     oplot, nu, sed, psym=3
;     print, nwp
;     print, p[2,*]
;     print, p[3,*]
;     print, check_priors(nwp,p[2,*], p[3,*] )
     if ( check_priors(nwp,p[2,*], p[3,*] ) ) then begin 
;         print, 'Priors OK...'
;stop
         sedp = sd_sed(xdata,nwp[0], nwp[1], nu1) + td_sed(xdata, nwp[2], nwp[3], nu2)
         amp = total(ydata[ndata-1]*sedp[ndata-1]) / total(sedp[ndata-1]^2)
;         newchisq = total( ( (ydata-sedp*amp)/0.25)^2 )
         newchisq = total( ( (ydata-sedp*amp)/de)^2 )
;     print, newchisq
         newlike = exp(-newchisq/2.)
;     print, oldlike, newlike
         likeratio = newlike/oldlike
;         print, 'likeratio', likeratio
         if (likeratio gt 1.) then begin
;             print, 'accepted: 1st', nwp
;             print, oldlike, newlike
             chlike[isamp] = newlike
             oldlike = newlike
             oldp = nwp
             oplot, xdata, sedp*amp, psym=3, col=210
             accepted = accepted + 1
             chp[isamp,0:np-1] = nwp
             chp[isamp,np] = amp
         endif

         if (likeratio le 1.) then begin
             t = randomu(-isamp-1)
             num[isamp,np] = t
             if (likeratio gt t) then begin
;                 print, 'accepted: 2nd', nwp
;                 print, oldlike, newlike
                 chlike[isamp] = newlike
                 oldlike = newlike
                 oldp = nwp
                 oplot, xdata, sedp*amp, psym=3, col=210
                 accepted = accepted + 1
                 chp[isamp,0:np-1] = nwp
                 chp[isamp,np] = amp
             endif
         endif
;     print, ' --- '
     endif
 endfor
; oplot, wx, wy, psym=1, col=245
;  oplot, wx, wy, col=245, thick=1.5

print, float(accepted)/nstep   

nbin=18
fitp = fltarr(np+1, 2, nbin)
peakp = fltarr(np+1)
!p.multi=[0,2,3]
for ip=0,np do begin
    gp = where(chp[*,ip] ne 0.d0)
    h=histogram(chp[gp,ip], locations=a, nbins=nbin)
    fitp[ip,0, *] = a
    fitp[ip,1, *] = h
    mm = max(h,im)
    peakp[ip] = a[im]
;    plot, a, h, psym=10
endfor

psed = sd_sed(xdata, peakp[0], peakp[1], nu1) + td_sed(xdata, peakp[2], peakp[3], nu2)
oplot, xdata, psed*peakp[4], col=70, thick=1.5

blike = max(chlike,im)
print, 'Best chisq = ', -2*alog(blike)
print, chp[im,*]
bsed = sd_sed(xdata, chp[im,0], chp[im,1], nu1) + td_sed(xdata, chp[im,2], chp[im,3], nu2)
oplot, xdata, bsed*chp[im,4], col=150, thick=1.5

;oplot, nu, (nu/30.)^(-3)
;oplot, nu, (nu/30.)^(-2.15)
oplot, nu, (sd_sed(nu, chp[im,0], chp[im,1], nu1) + td_sed(nu, chp[im,2], chp[im,3], nu2) )*chp[im,4]

legend, ['!17W7 Empirical SED (Pietrobon 2011)', '!17MCMC Marginal', '!17MCMC Best Fit'], col=[245,70,150], line=[0,0,0], /top, /right, chars=1.2

;!p.multi=0
window, 2, xsize=720*2/3, ysize=450*3*2/3 
for ip=0,np do begin
    plot, fitp[ip,0,*], fitp[ip,1,*], psym=10, chars=1.5
endfor

;window, 3
;for ip=0,np do begin
;    h=histogram(num[*,ip], locations=a, nbins=nbin)
;    plot, a, h, psym=10
;endfor

!p.multi=0
STOP                                                                                                                                                                                  
END                                                                                                                                                                                   
