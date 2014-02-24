true = 1b
false = 0b

do_png = false
do_compute = true

apply_mask = false ; Always, bad performance otherwise
;##mfile = 'dust_ps_mask.fits'
;##mfile = 'ffp4_commomn_mask_ns128.fits'

freq = [30, 44, 70, 100, 143, 217, 353]
sfreq = string(freq, format='(i3.3)')
nfreq = n_elements(freq)

cmb2sz = fltarr(nfreq)
h = 6.626068 * 10.^(-34)
k = 1.3806503 * 10.^(-23)
Tcmb = 2.725
x = (h*freq*10.^9) / (k * Tcmb)
cmb2sz = x * (exp(x)+1.)/(exp(x)-1.)-4.
;print, cmb2sz
;stop

ns = 128l
npix = 12l*ns^2
maps = fltarr(npix, nfreq)

lmax = 3*ns
cls_matrix = fltarr(lmax+1, nfreq, nfreq)
nalms = (lmax+2)*(lmax+1)/2
alms_array = fltarr(nalms, 2, nfreq)

if (do_compute) then begin
    for ifreq = 0, nfreq-1 do begin
        read_fits_map,'ffp4_01a_'+sfreq[ifreq]+'_ns128_uK.fits', temp
        maps[*,ifreq] = temp[*,0]
    endfor

;##for ifreq=0,6 do mollview, maps[*,ifreq], min=-300, max=300, units='!7l!6K CMB', win=ifreq+1 & stop

    window, 20 & plot_io, /nodata, [0,50], [1,10000], chars=1.5
    for ifreq=0, nfreq-1 do begin
        ianafast, maps[*,ifreq], cls, nlmax=3*ns, /ring, alm1_out='tmp_amls1.fits', /silent
        fits2alm, indx, alms, 'tmp_amls1.fits'
        alms_array[*,*,ifreq] = alms
        for jfreq=ifreq,nfreq-1 do begin
;##        print, ifreq, jfreq, icount
            if (not apply_mask) then ianafast, maps[*,ifreq], cls, nlmax=3*ns, map2_in=maps[*,jfreq], /ring, /silent else $
              ianafast, maps[*,ifreq], cls, nlmax=3*ns, map2_in=maps[*,jfreq], /ring, /silent, maskfile=mfile, regression=2
            cls_matrix[*, ifreq, jfreq] = cls
            cls_matrix[*, jfreq, ifreq] = cls
            oplot, cls[0:50]
;##        icount = icount + 1
        endfor
    endfor

    save, filename='hilc.sav', maps, cls_matrix, mask, alms_array, indx
endif else begin
    restore, 'hilc.sav'
endelse

e = fltarr(nfreq)+1.

ilc_alms = fltarr(nalms, 2)

for ell=2l,lmax do begin
    m1 = invert( reform( cls_matrix[ell,*,*] ), info, /double)
    if (info ne 0) then print, ell, info
    norm = reform( e ## ( m1 ## transpose( e ) ) )
    well = reform( m1 ## transpose( e ) ) / norm[0]
    if ( (ell/20)*20 eq ell ) then begin
        print, ell, ' / ', lmax
        print, well
    endif

    low_m = ell^2 + ell + 1
    high_m = ell^2 + 2*ell + 1
    ilow_m = where( indx eq low_m)
    ihigh_m = where( indx eq high_m)
    
;##    print, ell, low_m, high_m, ilow_m, ihigh_m

    for ifreq=0,nfreq-1 do ilc_alms[ilow_m:ihigh_m,*] = ilc_alms[ilow_m:ihigh_m,*] + well[ifreq]*alms_array[ilow_m:ihigh_m,*,ifreq]
endfor

alm2fits, indx, ilc_alms, 'ffp4_01a_ilc_alms.fits'

isynfast, 'camb_ctp3_init_cls.fits', ilc_map, nside=ns, nlmax=lmax, alm_in='ffp4_01a_ilc_alms.fits', simul_type=1
write_fits_map, 'ffp4_01a_hilc_map.fits', ilc_map, /ring, units='!7l!6K CMB'

mollview, 'ffp4_01a_hilc_map.fits', min=-300, max=300, win=10

print, ' --- End of program ---'

stop

end
