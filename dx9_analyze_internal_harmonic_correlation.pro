;## ====================================================================

!path=!path+':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'
!path=!path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'

True = 1b
False = 0b

incomp = ['cmb', 'co','dust', 'ff']
ncomp = n_elements(comp)

sets = ['dx9']

internal_corr = True

dataset = sets[0]

;---------------------------------------------------------------------

if (internal_corr) then begin

    Nside = 2048l
    Npix = 12l*Nside^2
    nComp = 4
    solution = fltarr(Npix, nComp)

    components = ['cmb', 'fg_co_amp', 'fg_dust_amp', 'fg_lowfreq_amp']

    sversion = 'v1.0.128'

    root='/global/scratch/sd/dpietrob/commander-ruler/'+dataset+'/'+sversion+'/dx9_commander-ruler_'+sversion+'_'

    if (False) then begin
    gmask = fltarr(Npix)+1.
    for icomp=0,ncomp-1 do begin
        read_fits_map, root+components[icomp]+'_uK.fits', map
        bp = where( (finite(map, /nan) eq True) or (map eq 0.) or (map eq -1.6375e30) )
        map[bp] = -1.6375e30 ;0.
        gmask[bp] = 0.
        solution[*,icomp] = map
    endfor
    
    lmax = 3000l
    xcls = fltarr(lmax+1, nComp*(Ncomp+1)/2)

    count = 0l
    
;##    window, 1
;##    plot_oi, /nodata, [1,lmax], [-1.e4,1.e4]
;##    l = findgen(lmax+1)
;##    ll = l*(l+1)/2./!pi
;##    il = lindgen(lmax-1)+2
    for icomp=0,nComp-1 do begin
        print, count, icomp
        ianafast, solution[*,icomp], cls, /ring, nlmax=lmax, simul_type=1, regression=2
        xcls[*,count] = cls
;##        oplot, l[il], cls[il]*ll[il], thick=0.5, col=245

;##        ianafast, root+components[icomp]+'_hr1_uK.fits', cls, /ring, nlmax=lmax, simul_type=1, /silent, map2_in=root+components[icomp]+'_hr2_uK.fits', regression=2
;##        oplot, l[il], cls[il]*ll[il], thick=1.5, col=50*(icomp+1)

        count = count+1
        for jcomp=icomp+1,nComp-1 do begin
            print, count, icomp, jcomp
            ianafast, solution[*,icomp], cls, /ring, nlmax=lmax, simul_type=1, map2_in=solution[*,jcomp], regression=2
            xcls[*,count] = cls
            count = count + 1
;##            oplot, l[il], cls[il]*ll[il], col=50*(icomp+1)+5*jcomp
;stop
        endfor
;stop
    endfor

    save, filename='dx9_internal_xcsl.sav', xcls
endif else begin

    restore, 'dx9_internal_xcsl.sav'
    
endelse

    corr = xcls

    corr[*,1] = xcls[*,1] / sqrt( xcls[*,0] * xcls[*,4] )
    corr[*,2] = xcls[*,2] / sqrt( xcls[*,0] * xcls[*,7] )
    corr[*,3] = xcls[*,3] / sqrt( xcls[*,0] * xcls[*,9] ) 

    corr[*,5] = xcls[*,5] / sqrt( xcls[*,4] * xcls[*,7] )
    corr[*,6] = xcls[*,6] / sqrt( xcls[*,4] * xcls[*,9] )

    corr[*,8] = xcls[*,8] / sqrt( xcls[*,7] * xcls[*,9] )

;##    mollview, randomn(-1,12), win=-1
    loadct, 39
    !p.color=0
    !p.background=255
    !p.charsize=1.5
    window, 1
    plot_oi, /nodata, [1,lmax], [-1.2,1.2]
    l = findgen(lmax+1)
    ll = l*(l+1)/2./!pi
    il = lindgen(lmax-1)+2

    oplot, l[il], l[il]*0.

    count = 0
    for icomp=0,ncomp-1 do begin
        count = count + 1
        for jcomp=icomp+1, ncomp-1 do begin
            oplot, l[il], corr[il,count], col=50*(icomp+1)+5*jcomp
            count = count + 1
        endfor
    endfor
stop



    clab = strarr(6)
    icnt = 0
    for i=0,3 do begin
        for j=i+1,3 do begin
            clab[icnt] = components[i]+'-'+components[j]
            out = reorder(correlation_maps[*,icnt], in='nest', out='ring')
            write_fits_map, dataset+'_cross_commander_correlation_'+clab[icnt]+'_'+strtrim(string(low_Nside),2)+'_'+sversion+'.fits', out, /ring
            mollview, out, tit='Correlation '+clab[icnt], min=-1, max=1, chars=1.5, win=icnt+1
            mollview, out, tit='Correlation '+clab[icnt], min=-1, max=1, chars=1.5, png=dataset+'_cross_cr_correlation_'+clab[icnt]+'_'+strtrim(string(low_Nside),2)+'_'+sversion+'.png', win=-1
            icnt = icnt + 1
            print, clab
        endfor
    endfor


endif

;# ---------------------------------------------------------------------

   STOP


END
