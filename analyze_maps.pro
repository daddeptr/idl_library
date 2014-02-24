!path=!path+':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'
!path=!path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'

True = 1b
False = 0b

comp = ['cmb', 'co_amp','dust_amp', 'lowfreq_amp']
ncomp = n_elements(comp)

sets = ['dx8', 'ffp4']

do_anafast    = True
do_30deg      = False
do_noise      = False
do_noise30deg = False
do_857        = False
do_SDF        = True
do_30degSDF   = False
do_30deg857   = False
do_143        = True
do_30deg143   = False

dataset = sets[0]

;## mollview, randomn(-1,12), win=-1
loadct, 39
!p.color=0
!p.background=255

if (False) then begin
    !p.multi=[0,1,2]
    window, 3, xs=720, ys=450*2

    for i=0,1 do begin
        dataset = sets[i]
        
        if (do_anafast) then begin
;##            ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_hr1_uK.fits', dataset+'_xcls.fits', regression=2, maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', map2_in=dataset+'/'+dataset+'_commander-ruler_v3_cmb_hr2_uK.fits'

;##            ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_error_uK.fits', dataset+'_nls.fits', regression=2, maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits'
      
;##            ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cls.fits', regression=2, maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits'
if (dataset eq 'ffp4') then begin
            ianafast, '/global/scratch/sd/dpietrob/magique_outputs/ffp4_7b_mondip_hr1_avrg_cmb.fits', dataset+'_xcls.fits', regression=2, maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', map2_in='/global/scratch/sd/dpietrob/magique_outputs/ffp4_7b_mondip_hr2_avrg_cmb.fits'

read_fits_map, '/global/scratch/sd/dpietrob/magique_outputs/ffp4_7b_mondip_hr1_avrg_cmb.fits', hr1
read_fits_map, '/global/scratch/sd/dpietrob/magique_outputs/ffp4_7b_mondip_hr2_avrg_cmb.fits', hr2
bp=where((hr1 eq -1.6375e30) or (hr2 eq -1.6375e30))
d=(hr1-hr2)2.
d[bp]=0.
write_fits_map,'new_'+dataset+'_7b_mondip_error_avrg_cmb.fits', d, /ring, units='!7l!6K CMB'
            ianafast, 'new_'+dataset+'_7b_mondip_error_avrg_cmb.fits', dataset+'_nls.fits', regression=2, maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits'
      
            ianafast, '/global/scratch/sd/dpietrob/magique_outputs/ffp4_7b_mondip_avrg_cmb.fits', dataset+'_cls.fits', regression=2, maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits'
        endif else begin
            print, 'Not ready'
            stop
        endelse
        endif
   
        if (do_30deg) then begin
            ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_hr1_uK.fits', dataset+'_xcls_30d.fits', regression=2, maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', map2_in=dataset+'/'+dataset+'_commander-ruler_v3_cmb_hr2_uK.fits', theta_cut_deg = 30
            
            ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_error_uK.fits', dataset+'_nls_30d.fits', regression=2, maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', theta_cut_deg = 30
            
            ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cls_30d.fits', regression=2, maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', theta_cut_deg = 30
        endif

        fits2cl, xcls, dataset+'_xcls.fits'
        fits2cl, cls, dataset+'_cls.fits'
        fits2cl, nls, dataset+'_nls.fits'
        
        fits2cl, xcls_30d, dataset+'_xcls_30d.fits'
        fits2cl, cls_30d, dataset+'_cls_30d.fits'
        fits2cl, nls_30d, dataset+'_nls_30d.fits'
   
        l=findgen(4097)
        ll=l*(l+1)/2./!pi
        il=lindgen(4095)+2

        plot_oo, l[il], xcls[il]*ll[il], chars=1.5, xr=[1,4000], xtit='!6l', ytit='!6l(l+1)C!dl!n/2!7p!6', yr=[1,10000], tit='!6pC!dl!n - '+dataset, ns=5
        oplot, l[il], nls[il]*ll[il], col=70, ns=5
        oplot, l[il], cls[il]*ll[il], col=70, ns=5
        oplot, l[il], xcls[il]*ll[il], col=245, ns=5
        oplot, l[il], nls_30d[il]*ll[il], col=100, ns=5
        oplot, l[il], cls_30d[il]*ll[il], col=100, ns=5
        oplot, l[il], xcls_30d[il]*ll[il], col=215, ns=5
  
        legend, ['!6xC!dl!n', '!6C!dl!n', '!6N!dl!n', '!6xC!dl!u30deg!n', '!6C!dl!u30deg!n', '!6N!dl!u30deg!n'], col=[245,70,70,215,100,100], line=lindgen(6)*0., /bottom, /left, chars=1.5
    endfor
endif

if (True) then begin
    cls=fltarr(4097,10)
    cls_30d=fltarr(4097,10)

    nls = fltarr(4097,ncomp)
    nls_30d = fltarr(4097,ncomp)

    count = 0
;    read_fits_map, dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', mask
;    gp = where(mask gt 0.)
;    bp = where(mask eq 0.)

    if (do_857 and (dataset eq 'dx8') ) then begin
;        read_fits_map, '/project/projectdirs/planck/data/mission/DPC_maps/dx8/hfi/FREQ/HFI_857_2048_20120113.fits', m857
;        m857 = m857*1.e6
        ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cmb_X_857_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in=dataset+'_857_map_uK.fits'
        ianafast, dataset+'_857_map_uK.fits', dataset+'_857_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, simul_type=1
        ianafast, dataset+'_857_map_uK.fits', dataset+'_857_X_SDF_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits', simul_type=1
    endif

    if (do_857 and (dataset eq 'ffp4') ) then begin
        ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cmb_X_857_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='/project/projectdirs/planck/data/ffp4/challenge01/map01_857.fits', simul_type=1
        ianafast, '/project/projectdirs/planck/data/ffp4/challenge01/map01_857.fits', dataset+'_857_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, simul_type=1
        ianafast, dataset+'_857_map_uK.fits', dataset+'_857_X_SDF_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits', /ring, simul_type=1
    endif

; ---------
    if (do_30deg857 and (dataset eq 'dx8') ) then begin
;        read_fits_map, '/project/projectdirs/planck/data/mission/DPC_maps/dx8/hfi/FREQ/HFI_857_2048_20120113.fits', m857
;        m857 = m857*1.e6
;        write_fits_map, dataset+'_857_map_uK.fits', m857[*,0], /ring, units='!7l!8K CMB'
        ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cmb_X_857_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in=dataset+'_857_map_uK.fits', theta_cut_deg=30.
        ianafast, dataset+'_857_map_uK.fits', dataset+'_857_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, simul_type=1, theta_cut_deg=30.
        ianafast, dataset+'_857_map_uK.fits', dataset+'_857_X_SDF_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits', simul_type=1, theta_cut_deg=30.
    endif

    if (do_30deg857 and (dataset eq 'ffp4') ) then begin
        ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cmb_X_857_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='/project/projectdirs/planck/data/ffp4/challenge01/map01_857.fits', simul_type=1, theta_cut_deg=30.
        ianafast, '/project/projectdirs/planck/data/ffp4/challenge01/map01_857.fits', dataset+'_857_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, simul_type=1, theta_cut_deg=30.
        ianafast, dataset+'_857_map_uK.fits', dataset+'_857_X_SDF_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits', /ring, simul_type=1, theta_cut_deg=30.
    endif

    fits2cl, xcls857, dataset+'_cmb_X_857_cls.fits'
    fits2cl, cls857, dataset+'_857_cls.fits'
    fits2cl, xcls857_SDF, dataset+'_857_X_SDF_cls.fits'

    fits2cl, xcls857_30d, dataset+'_cmb_X_857_cls_30d.fits'
    fits2cl, cls857_30d, dataset+'_857_cls_30d.fits'
    fits2cl, xcls857_SDF_30d, dataset+'_857_X_SDF_cls_30d.fits'

;# -------------------------------------------
    if (do_143 and (dataset eq 'dx8') ) then begin
        read_fits_map, '/project/projectdirs/planck/data/mission/DPC_maps/dx8/hfi/FREQ/HFI_143_2048_20120113.fits', m143
        m143 = m143*1.e6
        write_fits_map, dataset+'_143_map_uK.fits', m143[*,0], /ring, units='!7l!8K CMB'
        ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cmb_X_143_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in=dataset+'_143_map_uK.fits'
        ianafast, dataset+'_143_map_uK.fits', dataset+'_143_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, simul_type=1
        ianafast, dataset+'_143_map_uK.fits', dataset+'_143_X_SDF_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits', simul_type=1
    endif

    if (do_143 and (dataset eq 'ffp4') ) then begin
        ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cmb_X_143_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='/project/projectdirs/planck/data/ffp4/challenge01/map01_143.fits', simul_type=1
        ianafast, '/project/projectdirs/planck/data/ffp4/challenge01/map01_143.fits', dataset+'_143_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, simul_type=1
        ianafast, dataset+'_143_map_uK.fits', dataset+'_143_X_SDF_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits', /ring, simul_type=1
    endif

; ---------
    if (do_30deg143 and (dataset eq 'dx8') ) then begin
        ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cmb_X_143_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in=dataset+'_143_map_uK.fits', theta_cut_deg=30.
        ianafast, dataset+'_143_map_uK.fits', dataset+'_143_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, simul_type=1, theta_cut_deg=30.
        ianafast, dataset+'_143_map_uK.fits', dataset+'_143_X_SDF_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits', simul_type=1, theta_cut_deg=30.
    endif

    if (do_30deg143 and (dataset eq 'ffp4') ) then begin
        ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cmb_X_143_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='/project/projectdirs/planck/data/ffp4/challenge01/map01_143.fits', simul_type=1, theta_cut_deg=30.
        ianafast, '/project/projectdirs/planck/data/ffp4/challenge01/map01_143.fits', dataset+'_143_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, simul_type=1, theta_cut_deg=30.
        ianafast, dataset+'_143_map_uK.fits', dataset+'_143_X_SDF_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits', /ring, simul_type=1, theta_cut_deg=30.
    endif

    fits2cl, xcls143, dataset+'_cmb_X_143_cls.fits'
    fits2cl, cls143, dataset+'_143_cls.fits'
    fits2cl, xcls143_SDF, dataset+'_143_X_SDF_cls.fits'

    fits2cl, xcls143_30d, dataset+'_cmb_X_143_cls_30d.fits'
    fits2cl, cls143_30d, dataset+'_143_cls_30d.fits'
    fits2cl, xcls143_SDF_30d, dataset+'_143_X_SDF_cls_30d.fits'
;# -------------------------------------------
    if (do_SDF ) then begin
;        ianafast, 'SFD_i100_healpix_1024.fits', 'tempcls.fits', alm1_out='tempalms.fits', nlmax=4096, simul_type=1
;        isynfast, 'tempcls.fits', 'SFD_i100_healpix_2048.fits', alm_in='tempalms.fits', nside=2048, nlmax=4096, simul_type=1
;        mollview, 'SFD_i100_healpix_2048.fits', /hist
        ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cmb_X_SDF_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits'
        ianafast, 'SFD_i100_healpix_2048.fits', dataset+'_SDF_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, simul_type=1
    endif

    if (do_30degSDF ) then begin
        ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cmb_X_SDF_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits', theta_cut_deg=30.
        ianafast, 'SFD_i100_healpix_2048.fits', dataset+'_SDF_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, simul_type=1, theta_cut_deg=30.
    endif

    fits2cl, xclsSDF, dataset+'_cmb_X_SDF_cls.fits'
    fits2cl, clsSDF, dataset+'_SDF_cls.fits'
    fits2cl, xclsSDF_30d, dataset+'_cmb_X_SDF_cls_30d.fits'
    fits2cl, clsSDF_30d, dataset+'_SDF_cls_30d.fits'

;# -------------------------------------------
    for icomp=0,ncomp-1 do begin
;        read_fits_map, dataset+'/'+dataset+'_commander-ruler_v3_'+comp[icomp]+'_uK.fits', map1
        if (do_anafast) then ianafast, dataset+'/'+dataset+'_commander-ruler_v3_'+comp[icomp]+'_uK.fits', dataset+'_'+comp[icomp]+'_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2
        if (do_noise) then ianafast, dataset+'/'+dataset+'_commander-ruler_v3_'+comp[icomp]+'_error_uK.fits', dataset+'_'+comp[icomp]+'_nls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2
; 
        if (do_30deg) then ianafast, dataset+'/'+dataset+'_commander-ruler_v3_'+comp[icomp]+'_uK.fits', dataset+'_'+comp[icomp]+'_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, theta_cut_deg=30.
        if (do_noise30deg) then ianafast, dataset+'/'+dataset+'_commander-ruler_v3_'+comp[icomp]+'_error_uK.fits', dataset+'_'+comp[icomp]+'_nls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, theta_cut_deg=30.

        fits2cl, tcls, dataset+'_'+comp[icomp]+'_cls.fits'
        cls[*,count] = tcls
        fits2cl, tcls, dataset+'_'+comp[icomp]+'_cls_30d.fits'
        cls_30d[*,count] = tcls
        fits2cl, tcls, dataset+'_'+comp[icomp]+'_nls.fits'
        nls[*,icomp] = tcls
        fits2cl, tcls, dataset+'_'+comp[icomp]+'_nls_30d.fits'
        nls_30d[*,icomp] = tcls

        count = count + 1
        for jcomp=icomp+1,ncomp-1 do begin
;            read_fits_map, dataset+'/'+dataset+'_commander-ruler_v3_'+comp[jcomp]+'_uK.fits', map2

;            print, ' - '+comp[icomp]+ ' X '+comp[jcomp] + ' correlation = ', correlate(map1[gp], map2[gp]), correlate(map1[bp], map2[bp]), correlate(map1, map2)

            if (do_anafast) then ianafast, dataset+'/'+dataset+'_commander-ruler_v3_'+comp[icomp]+'_uK.fits', dataset+'_'+comp[icomp]+'_X_'+comp[jcomp]+'_cls.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in=dataset+'/'+dataset+'_commander-ruler_v3_'+comp[jcomp]+'_uK.fits'
            if (do_30deg) then ianafast, dataset+'/'+dataset+'_commander-ruler_v3_'+comp[icomp]+'_uK.fits', dataset+'_'+comp[icomp]+'_X_'+comp[jcomp]+'_cls_30d.fits', maskfile=dataset+'/'+dataset+'_commander-ruler_v3_mask.fits', regression=2, map2_in=dataset+'/'+dataset+'_commander-ruler_v3_'+comp[jcomp]+'_uK.fits', theta_cut_deg=30.
            fits2cl, tcls, dataset+'_'+comp[icomp]+'_X_'+comp[jcomp]+'_cls.fits'
            cls[*,count] = tcls
            fits2cl, tcls, dataset+'_'+comp[icomp]+'_X_'+comp[jcomp]+'_cls_30d.fits'
            cls_30d[*,count] = tcls
            count = count + 1
        endfor
    endfor

    window, 13, xs=720*1.5, ys=450*1.5
    l=findgen(4097)
    il = lindgen(4095)+2
    ll =l*(l+1)/2./!pi
;##     !p.multi=[0,1,2]
    plot_oi, l[il], cls[il,0]*0., xr=[1, 4000], yr=[-1.2,1.2], chars=1.5, xtit='!6l', ytit='!6C!dl!u12!n', tit='!6Correlation Spectra'
    oplot, l[il], cls[il,1]/sqrt(cls[il,0]*cls[il,4]), col=70, ns=5
    oplot, l[il], cls[il,2]/sqrt(cls[il,0]*cls[il,7]), col=245, ns=5
    oplot, l[il], cls[il,3]/sqrt(cls[il,0]*cls[il,9]), col=210, ns=5
;    oplot, l[il], cls[il,1]/sqrt(cls[il,0]*cls[il,4])
    oplot, l[il], cls_30d[il,1]/sqrt(cls_30d[il,0]*cls_30d[il,4]), col=100, ns=5
    oplot, l[il], cls_30d[il,2]/sqrt(cls_30d[il,0]*cls_30d[il,7]), col=225, ns=5
    oplot, l[il], cls_30d[il,3]/sqrt(cls_30d[il,0]*cls_30d[il,9]), col=200, ns=5
;
    oplot, l[il], xcls857/sqrt(cls[il,0]*cls857), col=40, ns=5
    oplot, l[il], xclsSDF/sqrt(cls[il,0]*clsSDF), col=150, ns=5

    oplot, l[il], xcls857_SDF/sqrt(cls857*clsSDF), col=60, ns=5
    oplot, l[il], xcls143_SDF/sqrt(cls143*clsSDF), col=0, ns=5

    oplot, l[il], xcls857_30d/sqrt(cls_30d[il,0]*cls857_30d), col=50, ns=5
    oplot, l[il], xclsSDF_30d/sqrt(cls_30d[il,0]*clsSDF_30d), col=170, ns=5

    oplot, l[il], xcls857_SDF_30d/sqrt(cls857_30d*clsSDF_30d), col=90, ns=5
    oplot, l[il], xcls143_SDF_30d/sqrt(cls143_30d*clsSDF_30d), col=10, ns=5

    legend, ['cmb X co', 'cmb X dust', 'cmb X low-freq', 'cmb X 857GHz', 'cmb X SDFi100', '857 X SDFi100'], col=[70,245,210,40,150,60], line=[0,0,0,0,0,0], chars=1.2, /bottom, /left

    if (False) then begin
        pivot_f = [101.2, 353., 30.]
        cf = conversionfactor(pivot_f,/antenna2thermo)

        plot_oo, l[il], cls[il,0]*ll[il], xr=[1,4000], yr=[.1,1000000], chars=1.5, xtit='!6l!n', ytit='!6l(l+1)C!dl!n/2!7p!6', tit='!6Component Spectra in !7l!6K CMB'
        oplot, l[il], nls[il,0]*ll[il]
        oplot, l[il], cls[il,4]*ll[il]*cf[0]^2, col=245
        oplot, l[il], nls[il,1]*ll[il]*cf[0]^2, col=245
        oplot, l[il], cls[il,7]*ll[il]*cf[1]^2, col=70
        oplot, l[il], nls[il,2]*ll[il]*cf[1]^2, col=70
        oplot, l[il], cls[il,9]*ll[il]*cf[2]^2, col=210
        oplot, l[il], nls[il,3]*ll[il]*cf[2]^2, col=210
;
        oplot, l[il], cls_30d[il,0]*ll[il]
        oplot, l[il], nls_30d[il,0]*ll[il]
        oplot, l[il], cls_30d[il,4]*ll[il]*cf[0]^2, col=225
        oplot, l[il], nls_30d[il,1]*ll[il]*cf[0]^2, col=225
        oplot, l[il], cls_30d[il,7]*ll[il]*cf[1]^2, col=100
        oplot, l[il], nls_30d[il,2]*ll[il]*cf[1]^2, col=100
        oplot, l[il], cls_30d[il,9]*ll[il]*cf[2]^2, col=200
        oplot, l[il], nls_30d[il,3]*ll[il]*cf[2]^2, col=200
        legend, ['!6cmb', '!6co @ 101.2GHz', '!6dust @ 353GHz', '!6low-freq @ 30GHz'], col=[0,245,70,210], line=[0,0,0,0], chars=1.2, /bottom, /right
        !p.multi=0
    endif
stop
endif



 
   STOP

END
