 !path=!path+':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'
 !path=!path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'

 true  = 1b
 false = 0b

 mollview, randomn(-1,12), window=-1
 loadct, 39
 !p.color=0
 !p.background=255

 ix=findgen(1001)
 nu = .408+ix*(900.)/1000
 b=14
 alpha=-1.
 nu0 = 20.
 nu1 = 545.
 window,1 
 plot_oo, nu, (nu/30.)^(-3), xr=[1,1000], yr=[.001,10]
 readcol, '../../wmap/ns128/spin_dust_spectrum_v2.txt', wx, wy
 oplot, wx, wy, psym=1, col=230
 k_B     = 1.3806503d-23
 h       = 6.626068d-34
 c       = 2.99792458d8
 t=18
 emis=1.8
 x = h / (k_B*Td)

 print, min(abs(nu-wx[5]),im), im

 td=(exp(x*nu1)-1.d0) / (exp(x*nu)-1.d0) * (nu / nu1)^(emis+1.d0)
; oplot, nu, td*wy[4]/td[im]
 sd = exp(-((nu-nu0)/b)^2/2.+alpha)
 oplot, nu, (sd+td)*wy[5]/(sd[im]+td[im]), col=210

 nstep = 1000
 p = [[-0.6,0.1], $
      [14,1], $
      [1.8, 0.3], $
      [18, 3] ]


STOP
END

; =========================================================================
pro check_pix_sed, ipix=ipix,alpha=alpha,beta=beta

!path=!path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'
!path=!path+':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'

freq = [22.8, 33, 41, 61, 94, 30, 44, 70, 100, 143, 217, 353, 545, 0.408]
nfreq = n_elements(freq)
sfreq = ['K', 'Ka', 'Q', 'V', 'W', '030', '044', '070', '100', '143', '217', '353', '545', 'Haslam']

f = fltarr(nfreq)
sed = fltarr(nfreq)
er = fltarr(nfreq)

for ifreq=0,nfreq-1 do begin
   if (ifreq lt 5) then begin
      read_fits_map, 'wmap7_map_'+sfreq[ifreq]+'.fits', map
      read_fits_map, 'wmap7_rms_'+sfreq[ifreq]+'.fits', rms
   endif

   if ( (ifreq ge 5) and (ifreq le 12) ) then begin
      read_fits_map, 'dx7_map_'+sfreq[ifreq]+'.fits', map
      read_fits_map, 'dx7_rms_'+sfreq[ifreq]+'.fits', rms
   endif

   if (ifreq eq 13) then begin
      read_fits_map, 'haslam_ns128.fits', map
      read_fits_map, 'haslam_rms.fits', rms
   endif

   f[ifreq] = freq[ifreq]
   print, freq[ifreq],conversionfactor(freq[ifreq],/thermo2antenna)
   sed[ifreq] = map[ipix,0] * conversionfactor(freq[ifreq],/thermo2antenna)
   er[ifreq] = rms[ipix,0] * conversionfactor(freq[ifreq],/thermo2antenna)
endfor

mollview, randomn(-1,12), window=-1
loadct, 39
!p.color=0
!p.background=255

!p.thick=2
window,1, xsize=720*1.5, ysize=450*1.5
plot_oo, f, sed, psym=4, chars=2., yr=[1,1.e9], xtit='!7m!8 (GHz)', ytit='!17SED'
;oplot, f, sed, col=245, psym=4
;errplot, f, sed-er, sed+er

x=findgen(550)+0.4
x_piv = beta
x_0 = beta
;ps = conversionfactor(x,/antenna2thermo) * ( ( (x/x_piv)^alpha+(x/x_piv)^beta) / ( (x_0/x_piv)^alpha+(x_0/x_piv)^beta) ) * ( ( tanh((x-18.)/2.)+1.) /2.)
;ps = ( ( (x/x_piv)^alpha+(x/x_piv)^beta) / ( (x_0/x_piv)^alpha+(x_0/x_piv)^beta) ) * ( ( 1.*tanh((x-25.)/4.)+1.) /2.)

  X0 = beta;*1.d9
  b = X0 / alpha
  a = X0 * (1.-alog(X0))

  ps = x^alpha * exp(-(x-a)/b) / ( x_piv^alpha * exp(-(x_piv-a)/b) )

mn = min( abs(freq-x_0), im )
norm = sed[im]
oplot, x, ps*norm[0], col=70
oplot, x, (x/22.8)^(-3.05)*sed[0], col=245

oplot, x, (x/30.)^(-2.15)*sed[5], col=100
oplot, x, (x/353.)^(1.5)*sed[11], col=210

legend, ['Synchrotron', 'Free-free', 'Thermal Dust', 'Spinning Dust'], col=[245,100,210,70], line=[0,0,0,0], /top,/right, chars=1.5
stop


end
