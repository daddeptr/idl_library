   scratch_dir = '~dpietrob/Gscratch'

   freq=[30, 44, 70, 100, 143, 217, 353]
   sfreq=['30', '44', '70', '100', '143', '217', '353']
;   sfreq=string(freq, format='(1i3.3)')
   nfreq = n_elements(freq)

;   real_beam = [60.*0.54424700, 60.*0.46527683, 60.*0.21676942, 9.94,   7.04,    4.66,    4.41,    4.47,    4.23]

   real_beam = [32.65, 27.00, 13.01, 9.94,   7.04,    4.66,    4.41,    4.47,    4.23]

true = 1b
false = 0b

dx = 'dx7'
maskfile = 'dx7_mask_common.fits.fits'

if (true) then begin
   print, ' reading files...'

   read_fits_map,dx+'_Imap100GHz_ns2048_uK.fits',m100
   read_fits_map,dx+'_Irms100GHz_ns2048_uK.fits',r100

   read_fits_map,dx+'_Imap143GHz_ns2048_uK.fits',m143
   read_fits_map,dx+'_Irms143GHz_ns2048_uK.fits',r143

   read_fits_map,dx+'_Imap217GHz_ns2048_uK.fits',m217
   read_fits_map,dx+'_Irms217GHz_ns2048_uK.fits',r217

   read_fits_map,dx+'_Imap353GHz_ns2048_uK.fits',m353
   read_fits_map,dx+'_Irms353GHz_ns2048_uK.fits',r353

   readcol,scratch_dir+'/graca/planck/luca/output_nominal/transFn/transFn_GH_100_nominal_symm.txt',bl100, format='x,x,x,x,x,x,x,x,f'
   bl100 = sqrt(bl100)
   bl100[0:1] = 1.

   readcol,scratch_dir+'/graca/planck/luca/output_nominal/transFn/transFn_GH_143_nominal_symm.txt',bl143, format='x,x,x,x,x,x,x,x,f'
   bl143 = sqrt(bl143)
   bl143[0:1] = 1.

   readcol,scratch_dir+'/graca/planck/luca/output_nominal/transFn/transFn_GH_217_nominal_symm.txt',bl217, format='x,x,x,x,x,x,x,x,f'
   bl217 = sqrt(bl217)
   bl217[0:1] = 1.

   readcol,scratch_dir+'/graca/planck/luca/output_nominal/transFn/transFn_GH_353_nominal_symm.txt',bl353, format='x,x,x,x,x,x,x,x,f'
   bl353 = sqrt(bl353)
   bl353[0:1] = 1.

   read_fits_map, maskfile, mask
   gp = where(mask gt 0.)
   bp = where(mask eq 0.)

   bl2 = ( bl217*1./mean(r217^2) + bl353*1./mean(r353^2) ) / ( 1./mean(r217^2) + 1./mean(r353^2) )
   bl3 = ( bl143*1./mean(r143^2) + bl217*1./mean(r217^2) + bl353*1./mean(r353^2) ) / ( 1./mean(r143^2) + 1./mean(r217^2) + 1./mean(r353^2) )
   bl4 = ( bl100*1./mean(r100^2) + bl143*1./mean(r143^2) + bl217*1./mean(r217^2) + bl353*1./mean(r353^2) ) / ( 1./mean(r100^2) + 1./mean(r143^2) + 1./mean(r217^2) + 1./mean(r353^2) )

   cl2fits, bl2, 'chains/'+dx+'_2b_effBeam_v2.fits'
   cl2fits, bl3, 'chains/'+dx+'_3b_effBeam_v2.fits'
   cl2fits, bl4, 'chains/'+dx+'_4b_effBeam_v2.fits'

   stop
endif

if (false) then begin
    m3 = (m143 * 1./r143^2 + m217 * 1./r217^2 + m353 * 1./r353^2 ) / (1./r143^2 + 1./r217^2 + 1./r353^2)

    mollview, m3, chars=1.5, min=-300, max=300
    write_fits_map, dx+'_3b_inv_var_cmb.fits', m3, /ring

    mollview, 'chains/3b/fg_amp_map_no02_c0001_k00001.fits', chars=1.5, min=-300, max=300
    
    ianafast, 'chains/3b/fg_amp_map_no02_c0001_k00001.fits', xxx, maskfile=dx+'_mask_ns2048.fits', regression=2, iter=2

    ianafast, dx+'_3b_inv_var_cmb.fits', yyy, maskfile=dx+'_mask_ns2048.fits', regression=2, iter=2

    l = findgen(4097)
    ll = l*(l+1)/2./!pi

    mollview, randomn(-1,12), window=-1
    loadct,39
    !p.color=0
    !p.background=255

    window, 1 & plot, l, xxx*ll
    oplot, l, yyy*ll, col=245

endif

if (true) then begin
    fits2cl,cl,'chains/'+dx+'_2b_effBeam_v2.fits'
    openw,1,'chains/'+dx+'_2b_effBeam_v2.dat'
    for i=0,4096l do printf,1,i,cl[i]
    close,1

    fits2cl,cl,'chains/'+dx+'_3b_effBeam_v2.fits'
    openw,1,'chains/'+dx+'_3b_effBeam_v2.dat'
    for i=0,4096l do printf,1,i,cl[i]
    close,1

    fits2cl,cl,'chains/'+dx+'_4b_effBeam_v2.fits'
    openw,1,'chains/'+dx+'_4b_effBeam_v2.dat'
    for i=0,4096l do printf,1,i,cl[i]
    close,1
endif

end
