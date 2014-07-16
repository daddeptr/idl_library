

True = 1b
False = 0b

spawn, 'ls ../../mirror/foreground_templates/lambda*.fits', tfiles
print, tfiles
add_smth = sqrt( 60.^2-[6.,6.,60.]^2 )

read_fits_map, '../../mirror/PLA/COM_CompMap_Lfreqfor-commrul_0256_R1.00.fits', lowf, order=dpord, nside=dpns
if strupcase(dpord) ne 'RING' then lowf = reorder(lowf, in=dpord, out='RING')
map = lowf[*,0]
rms = lowf[*,1]
npix = 12l*dpns^2
;write_fits_map, 'low-freq_temp_m.fits', fltarr(npix)+1., /ring
;mollview, map, /hist
;mollview, rms, /hist

;for i=0,2 do begin
;    ismoothing, tfiles[i], 'low-freq_temp_'+strtrim(string(i),2)+'.fits', fwhm_arcmin=add_smth[i], /silent
;    ud_grade, 'low-freq_temp_'+strtrim(string(i),2)+'.fits', 'low-freq_temp_'+strtrim(string(i),2)+'.fits', nside_out=dpns, order_out='ring'
;    mollview, 'low-freq_temp_'+strtrim(string(i),2)+'.fits', /hist
;endfor

;maskf = 'data/mask/dx11/dx9_common_xGal06_QUmask_DX10_2048_T5.0_60_X_ps100-353_mask_ns2048.fits'
;ud_grade, maskf, 'low-freq_mask.fits', nside_out=dpns, order_out='ring'
;read_fits_map, 'low-freq_mask.fits', mask
;mask[ where( mask[*,0] lt 0.75), 0 ] = 0.
;mask[ where( mask[*,0] ge 0.75), 0 ] = 1.
;write_fits_map, 'low-freq_mask.fits', mask[*,0], /ring
;mollview, 'low-freq_mask.fits'

spawn, 'ls low-freq_temp*.fits', tfiles
print, tfiles
;for i=0,n_elements(tfiles)-1 do mollview, tfiles[i], /asinh

res = map_regression( tfiles,rms_map=rms, maskfile='low-freq_mask.fits', inmap=map, fit=fit, a_coeff=a )

l = findgen(513)
ll = l*(l+1)/2./!pi
il = lindgen(511)+2
window, 10 & plot, /nodata, [1,500], [1.e-2,2500], chars=2, /ylog, /xlog

cols=[0,70,100]

for i=0,2 do begin
    ianafast, tfiles[i], cls, maskfile='low-freq_mask.fits', regression=2
    oplot, l[il], cls[il,0]*ll[il]*a[i]^2, col=cols[i]
endfor

;stop
;spawn, 'ls ../../dx11_pre/maps/dx11_pre_Imap_*_full_uK.fits', dxfiles
;tmp = ilc(dxfiles[0:4], maskfile='low-freq_mask.fits', apply_smooth=sqrt(60.^2-[10.,7.03, 5, 5, 5]^2) )
;ud_grade, tmp[*,1], tmpu, nside_out=256, order_in='ring'
;write_fits_map, 'dx11_ilc_ns256_60a.fits', tmpu, /ring
ianafast, tmpu, cls, maskfile='low-freq_mask.fits', regression=2, /ring
oplot, l[il], cls[il,0]*ll[il], col=245
legend, ['FDS','Ha','Haslam', 'dx11-ILC'], col=[[cols],245], lin=0, /bottom, chars=1.5

stop

end
