;## ====================================================================

!path=!path+':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'
!path=!path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'

True = 1b
False = 0b

comp = ['cmb', 'co_amp','dust_amp', 'lowfreq_amp']
ncomp = n_elements(comp)

sets = ['dx8', 'ffp4']

mollview, randomn(-1,12), win=-1
loadct, 39
!p.color=0
!p.background=255

do_commander = False
do_gmca    = False
do_smica   = False
do_anafast = False

plot_corr                 = False
plot_pxcl                 = False
corr_maps                 = False
compare_cmd_ruler_maps    = False
smooth_corr               = False
internal_corr             = False
compare_old_new           = False
compute_high_res_chisq_31 = True
compute_high_res_chisq    = False
check_patch               = False
check_spectra             = False
compare_co                = False

avrg_maps                 = False

chname                    = False

dataset = sets[1]

if (chname) then begin
    sversion='v3.1'
    comp = ['cmb','fg_co_amp','fg_dust_amp','fg_lowfreq_amp']
    for ic=1,3 do begin
        read_fits_map, 'ffp4/'+sversion+'/ffp4_commander-ruler_'+sversion+'_'+comp[ic]+'_hr1_uK.fits', map1
        read_fits_map, 'ffp4/'+sversion+'/ffp4_commander-ruler_'+sversion+'_'+comp[ic]+'_hr2_uK.fits', map2

        diff = (map1-map2) /2.

        write_fits_map, 'ffp4/'+sversion+'/ffp4_commander-ruler_'+sversion+'_'+comp[ic]+'_error_uK.fits', diff, /ring, units='!7l!6K'
        mollview, 'ffp4/'+sversion+'/ffp4_commander-ruler_'+sversion+'_'+comp[ic]+'_error_uK.fits'
    endfor
endif

;# --------------------------------------------------------------------

if (compare_cmd_ruler_maps) then begin
    if (False) then begin
        ianafast, 'dx8/dx8_commander-ruler_v3_cmb_uK.fits', 'tempcls.fits', alm1_out='tempalms.fits', nlmax=3000
        fits2alm, index, alms, 'tempalms.fits'
        fits2cl, bl, 'dx8_7b240arcmin_bl.fits'
        for ell=0l,3000 do begin
            fm=ell^2+ell+1
            lm=ell^2+2*ell+1
            fmi=where(index eq fm)
            lmi=where(index eq lm)
;print, fm, lm
;print, fmi, lmi
;stop
            alms[fmi:lmi,0] = alms[fmi:lmi,0] * bl[ell]
            alms[fmi:lmi,1] = alms[fmi:lmi,1] * bl[ell]
        endfor
        alm2fits, index, alms, 'stempalms.fits'
        isynfast, 'tempcls.fits', 'dx8_commander-ruler_v3_cmb_uK_40arcmin.fits', alm_in='stempalms.fits', nlmax=3000, nside=2048
        ud_grade, 'dx8_commander-ruler_v3_cmb_uK_40arcmin.fits', 'dx8_commander-ruler_v3_cmb_uK_40arcmin_ns0256.fits', nside_out=256
        stop
    endif
    read_fits_map, 'dx8_commander-ruler_v3_cmb_uK_40arcmin_ns0256.fits', dgdx8
    restore, '/global/scratch/sd/dpietrob/dx8/output/chains/pix_raw_effreq_scan.res.sav'
    mollview, cmb-dgdx8, min=-15, max=15, tit='!6Commander Mean Posterior - Smoothed Ruler', chars=1.5, units='!7l!8K CMB', px=1000
endif

lat_cuts = [0,10,20,30,40,50,60,70,80]
ncuts = n_elements(lat_cuts)

if (plot_corr) then begin
    if (True) then begin
;##        read_fits_map, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dx8
        read_fits_map, '/global/scratch/sd/loris/CompSep/GLSS/FFP/4/Davide_new/pix_raw_effreq_scan_vtd_5deg_ps/'+dataset+'_7b_mondip_avrg_cmb.fits', dx8
        read_fits_map, '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/gmca/'+dataset+'_gmca_v1_cmb.fits', gmca_dx8
        read_fits_map, '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/smica/'+dataset+'_smica_needlet_cmb.fits', smica_dx8

        gmask = fltarr(12l*2048l^2)
        gmask[*] = 1.
        bp = where(dx8 eq 0.)
        gmca_dx8r = reorder(gmca_dx8, in='nest', out='ring')
        gmcabp = where(gmca_dx8r eq 0.)
        smicabp = where(smica_dx8 eq 0.)
;##    help, where(dx8 eq -1.6375e30)
;##    help, where(gmca_dx8 eq -1.6375e30)
;##    help, where(smica_dx8 eq -1.6375e30)

        gmask[bp] = 0.
        gmask[gmcabp] = 0.
        gmask[smicabp] = 0.

        write_fits_map, dataset+'_global_mask.fits', gmask, /ring
;    mollview, dataset+'_global_mask.fits'
;stop
    endif

    cls = fltarr(4097,3,ncuts)
    cls_gmca = fltarr(4097,3,ncuts)
    cls_smica = fltarr(4097,3,ncuts)
    for icut=0,ncuts-1 do begin
        scut = string(lat_cuts[icut], format='(i2.2)')
        if (do_commander) then begin
            if (do_anafast) then begin
;##                ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cmb_X_SDF_cls_'+scut+'d.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits', theta_cut_deg=lat_cuts[icut], simul_type=1, maskfile=dataset+'_global_mask.fits'

;##                ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_uK.fits', dataset+'_cmb_cls_'+scut+'d.fits', regression=2, theta_cut_deg=lat_cuts[icut], simul_type=1, maskfile=dataset+'_global_mask.fits'

;##                ianafast, 'SFD_i100_healpix_2048.fits', 'SDF_cls_'+scut+'d.fits', regression=2, simul_type=1, theta_cut_deg=lat_cuts[icut], maskfile=dataset+'_global_mask.fits'
                ianafast, '/global/scratch/sd/loris/CompSep/GLSS/FFP/4/Davide_new/pix_raw_effreq_scan_vtd_5deg_ps/'+dataset+'_7b_mondip_avrg_cmb.fits', dataset+'_cmb_X_SDF_cls_'+scut+'d.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits', theta_cut_deg=lat_cuts[icut], simul_type=1, maskfile=dataset+'_global_mask.fits'

                ianafast, '/global/scratch/sd/loris/CompSep/GLSS/FFP/4/Davide_new/pix_raw_effreq_scan_vtd_5deg_ps/'+dataset+'_7b_mondip_avrg_cmb.fits', dataset+'_cmb_cls_'+scut+'d.fits', regression=2, theta_cut_deg=lat_cuts[icut], simul_type=1, maskfile=dataset+'_global_mask.fits'

                ianafast, 'SFD_i100_healpix_2048.fits', 'SDF_cls_'+scut+'d.fits', regression=2, simul_type=1, theta_cut_deg=lat_cuts[icut], maskfile=dataset+'_global_mask.fits'
            endif
        endif
        fits2cl, xcls, dataset+'_cmb_X_SDF_cls_'+scut+'d.fits'
        fits2cl, cmbcls, dataset+'_cmb_cls_'+scut+'d.fits'
        fits2cl, sdfcls, 'SDF_cls_'+scut+'d.fits'

        cls[*,0,icut] = cmbcls
        cls[*,1,icut] = sdfcls
        cls[*,2,icut] = xcls

;# ------ GMCA
        if (do_gmca) then begin
            if (do_anafast) then begin
                ianafast, '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/gmca/'+dataset+'_gmca_v1_cmb.fits', dataset+'_gmca_cmb_X_SDF_cls_'+scut+'d.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits', theta_cut_deg=lat_cuts[icut], simul_type=1, maskfile=dataset+'_global_mask.fits'

                ianafast, '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/gmca/'+dataset+'_gmca_v1_cmb.fits', dataset+'_gmca_cmb_cls_'+scut+'d.fits', regression=2, theta_cut_deg=lat_cuts[icut], simul_type=1, maskfile=dataset+'_global_mask.fits'
            endif
        endif
        fits2cl, xcls, dataset+'_gmca_cmb_X_SDF_cls_'+scut+'d.fits'
        fits2cl, cmbcls, dataset+'_gmca_cmb_cls_'+scut+'d.fits'
        fits2cl, sdfcls, 'SDF_cls_'+scut+'d.fits'

        cls_gmca[*,0,icut] = cmbcls
        cls_gmca[*,1,icut] = sdfcls
        cls_gmca[*,2,icut] = xcls

;# ------ SMICA
        if (do_smica) then begin
            if (do_anafast) then begin
                ianafast, '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/smica/'+dataset+'_smica_needlet_cmb.fits', dataset+'_smica_cmb_X_SDF_cls_'+scut+'d.fits', regression=2, map2_in='SFD_i100_healpix_2048.fits', theta_cut_deg=lat_cuts[icut], simul_type=1, maskfile=dataset+'_global_mask.fits'

                ianafast, '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/smica/'+dataset+'_smica_needlet_cmb.fits', dataset+'_smica_cmb_cls_'+scut+'d.fits', regression=2, theta_cut_deg=lat_cuts[icut], simul_type=1, maskfile=dataset+'_global_mask.fits'
            endif
        endif
        fits2cl, xcls, dataset+'_smica_cmb_X_SDF_cls_'+scut+'d.fits'
        fits2cl, cmbcls, dataset+'_smica_cmb_cls_'+scut+'d.fits'
        fits2cl, sdfcls, 'SDF_cls_'+scut+'d.fits'

        cls_smica[*,0,icut] = cmbcls
        cls_smica[*,1,icut] = sdfcls
        cls_smica[*,2,icut] = xcls

endfor

    window, 13, xs=720*2, ys=450*2
    l=findgen(4097)
    il = lindgen(4095)+2
    ll =l*(l+1)/2./!pi
    !p.multi=[0,3,3]
    for icut=0,ncuts-1 do begin
        plot_oi, l[il], cls[il,0,0]*0., xr=[1, 4000], yr=[-1.,1.], chars=2., xtit='!6l', ytit='!6C!dl!u12!n', ys=1, tit='!6Correlation CMB X SDFi100: '+strtrim(string(lat_cuts[icut]),2)+'deg'
        oplot, l[il], cls[il,2,icut] / sqrt(cls[il,0,icut]*cls[il,1,icut]), col=245, ns=5, thick=1.5
        oplot, l[il], cls_gmca[il,2,icut] / sqrt(cls_gmca[il,0,icut]*cls_gmca[il,1,icut]), col=70, ns=5, thick=1.5
        oplot, l[il], cls_smica[il,2,icut] / sqrt(cls_smica[il,0,icut]*cls_smica[il,1,icut]), col=210, ns=5, thick=1.5
        legend, ['Commander', 'GMCA', 'Need-Smica'], col=[245,70,210], line=[0,0,0], chars=1.
    endfor
    write_png, dataset+'_corr_plot.png', tvrd(/true)
;##        legend, 'cmb X SDFi100 '+strtrim(string(lat_cuts),2)+'deg', col=lindgen(ncuts)*40, line=lindgen(ncuts)*0., chars=1.5, /bottom, /right
    !p.multi=0
endif

;# ---------------------------------------------------------------------
if (plot_pxcl) then begin
    pcls = fltarr(4097,ncuts)
    pcls_gmca = fltarr(4097,ncuts)
    pcls_smica = fltarr(4097,ncuts)
    for icut=0,ncuts-1 do begin
        scut = string(lat_cuts[icut], format='(i2.2)')

        
        if (do_anafast) then ianafast, dataset+'/'+dataset+'_commander-ruler_v3_cmb_hr1_uK.fits', dataset+'_cmb_xcls_'+scut+'d.fits', regression=2, theta_cut_deg=lat_cuts[icut], simul_type=1, map2_in=dataset+'/'+dataset+'_commander-ruler_v3_cmb_hr2_uK.fits', maskfile=dataset+'_global_mask.fits'
        fits2cl, cmbcls, dataset+'_cmb_xcls_'+scut+'d.fits'
        pcls[*,icut] = cmbcls

;# ------ GMCA
        if (do_anafast) then ianafast, '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/gmca/'+dataset+'_gmca_v1_cmb_hr1.fits', dataset+'_gmca_cmb_xcls_'+scut+'d.fits', regression=2, theta_cut_deg=lat_cuts[icut], simul_type=1, map2_in='/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/gmca/'+dataset+'_gmca_v1_cmb_hr2.fits', maskfile=dataset+'_global_mask.fits'
        fits2cl, cmbcls, dataset+'_gmca_cmb_xcls_'+scut+'d.fits'
        pcls_gmca[*,icut] = cmbcls*1.e12

;# ------ SMICA
        if (do_anafast) then ianafast, '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/smica/'+dataset+'_smica_needlet_cmb_hr1.fits', dataset+'_smica_cmb_xcls_'+scut+'d.fits', regression=2, theta_cut_deg=lat_cuts[icut], simul_type=1, map2_in='/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/smica/'+dataset+'_smica_needlet_cmb_hr2.fits', maskfile=dataset+'_global_mask.fits'
        fits2cl, cmbcls, dataset+'_smica_cmb_xcls_'+scut+'d.fits'
        pcls_smica[*,icut] = cmbcls
    endfor

    window, 14, xs=720*2, ys=450*2
    l=findgen(4097)
    il = lindgen(2998)+2
    ll =l*(l+1)/2./!pi
    fits2cl, beam, dataset+'/'+dataset+'_commander-ruler_v3_cmb_beam.fits'
    ;## fits2cl, gmca_beam, '/global/scratch/sd/dpietrob/dx8/maps/gmca/dx8_gmca_v1_cmb_beam.fits'
    gmca_beam=mrdfits('/global/scratch/sd/dpietrob/dx8/maps/gmca/dx8_gmca_v1_cmb_beam.fits',0)
    ;## fits2cl, smica_beam, '/global/scratch/sd/dpietrob/dx8/maps/smica/dx8_smica_needlet_cmb_beam.fits'
    smica_beam = gaussbeam(5., 3000)
    !p.multi=[0,3,3]
    for icut=0,ncuts-1 do begin
        plot_oo, l[il], pcls[il,icut]*ll[il]/beam[il]^2, xr=[1, 4000], yr=[1.e0,1.e4], chars=2., xtit='!6l', ytit='!6l(l+1)C!dl!n/2!7p!6', ys=1, tit='!6CMB Pseudo X-Spectra: '+strtrim(string(lat_cuts[icut]),2)+'deg', ns=5
        oplot, l[il], pcls[il,icut]*ll[il]/beam[il]^2, col=245, ns=5, thick=1.5
        oplot, l[il], pcls_gmca[il,icut]*ll[il]/gmca_beam[il]^2, col=70, ns=5, thick=1.5
        oplot, l[il], pcls_smica[il,icut]*ll[il]/smica_beam[il]^2, col=210, ns=5, thick=1.5
        legend, ['Commander', 'GMCA', 'Need-Smica'], col=[245,70,210], line=[0,0,0], chars=1.
    endfor
;##        legend, 'cmb X SDFi100 '+strtrim(string(lat_cuts),2)+'deg', col=lindgen(ncuts)*40, line=lindgen(ncuts)*0., chars=1.5, /bottom, /right
    !p.multi=0
endif

;# ---------------------------------------------------------------------

if (corr_maps) then begin
    
    sdf_file = 'SFD_i100_healpix_2048.fits'

    iris_file = 'IRIS_nohole_4_2048.fits'

    sversion = 'v3.1'
;##    if (dataset eq 'ffp4') then read_fits_map,'/global/scratch/sd/dpietrob/commander-ruler/ffp4/v4/ffp4_7b_mondip_fs_HK_cmb.fits', dx8, order=dpord 
    if (dataset eq 'ffp4') then read_fits_map,'/global/scratch/sd/dpietrob/commander-ruler/ffp4/'+sversion+'/ffp4_commander-ruler_'+sversion+'_cmb_uK.fits', dx8, order=dpord 

    if (dataset eq 'dx8') then read_fits_map,'/global/scratch/sd/loris/CompSep/GLSS/DX/8/Davide_new/pix_raw_effreq_scan_vtd_5deg_ps/dx8_7b_mondip_fs_avrg_cmb.fits', dx8, order=dpord 

    if ( (dpord eq 'ring') or (dpord eq 'RING') ) then begin
        print, ' - Reordering...'
        dx8r = reorder(dx8, in='ring', out='nest')
        dx8 = dx8r
        dx8r = 0.
    endif

    read_fits_map, '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/gmca/'+dataset+'_gmca_v1_cmb.fits', gmca_dx8, order=dpord
    if ( (dpord eq 'ring') or (dpord eq 'RING') ) then begin
        print, ' - Reordering...'
        gmca_dx8r = reorder(gmca_dx8, in='ring', out='nest')
        gmca_dx8 = gmca_dx8r
        gmca_dx8r = 0.
    endif

    read_fits_map, '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/smica/'+dataset+'_smica_needlet_cmb.fits', smica_dx8, order=dpord
    if ( (dpord eq 'ring') or (dpord eq 'RING') ) then begin
        print, ' - Reordering...'
        smica_dx8r = reorder(smica_dx8, in='ring', out='nest')
        smica_dx8 = smica_dx8r
        smica_dx8r = 0.
    endif
        
    Nside = 2048l
    
    Npix = 12l*Nside^2
    low_Nside = 64l
    low_Npix = 12l*low_Nside^2
    
    print, low_Npix
    
    read_fits_map, sdf_file, sdf, nside=mapns, order=mapord
    if (mapns ne Nside) then stop, 'Nside not matching'
    if ( (mapord ne 'nest') and (mapord ne 'NEST') ) then begin
        print, ' - Reordering: ', mapord, ' --> nest'
        sdfr=reorder(sdf, in=mapord, out='nest')
        sdf=sdfr
        sdfr = 0.
    endif
    
    read_fits_map, iris_file, iris, nside=mapns, order=mapord
    if (mapns ne Nside) then stop, 'Nside not matching'
    if ( (mapord ne 'nest') and (mapord ne 'NEST') ) then begin
        print, ' - Reordering: ', mapord, ' --> nest'
        irisr=reorder(iris, in=mapord, out='nest')
        iris=irisr
        irisr = 0.
    endif
        
    gmask = fltarr(Npix)+1.
    mpix = where( (dx8 eq 0.) or (smica_dx8 eq 0.) or (gmca_dx8 eq 0.) or (dx8 eq -1.6375e30) or (smica_dx8 eq -1.6375e30) or (gmca_dx8 eq -1.6375e30) or (finite(dx8, /nan) eq True) or (finite(sdf, /nan) eq True) or (finite(gmca_dx8, /nan) eq True) or (finite(smica_dx8, /nan) eq True) or (finite(iris, /nan) eq True) or (iris lt 0.) )
    
    gmask[mpix] = 0.
    bp = where(gmask eq 0.)
    gp = where(gmask eq 1.)
    
    correlation_maps = fltarr(low_Npix,3)
    correlation_maps_iris = fltarr(low_Npix,3)
    
    correlation_maps[*,*] = -1.6375e30
    correlation_maps_iris[*,*] = -1.6375e30
    
    skewness_maps = correlation_maps
    kurtosis_maps = correlation_maps
    
    pix2vec_nest, low_Nside, lindgen(low_Npix), vec, vertex

    for i=0, low_Npix-1 do begin
        
        if ( (i/1000)*1000 eq i) then print, i, ' / ', low_Npix, ' --> '+string(float(i)/low_Npix*100,format='(f4.1)')+'%'

        query_polygon, Nside, transpose( reform( vertex[i,*,*]) ), pix_list, inclusive=True, /nested
        
        corr_gp = where(gmask[pix_list] gt 0.)
        
        if (corr_gp[0] ne -1) then begin
            correlation_maps[i,0] = correlate(dx8[pix_list[corr_gp]], sdf[pix_list[corr_gp]])
            correlation_maps[i,1] = correlate(gmca_dx8[pix_list[corr_gp]], sdf[pix_list[corr_gp]])
            correlation_maps[i,2] = correlate(smica_dx8[pix_list[corr_gp]], sdf[pix_list[corr_gp]])
            
            correlation_maps_iris[i,0] = correlate(dx8[pix_list[corr_gp]], iris[pix_list[corr_gp]])
            correlation_maps_iris[i,1] = correlate(gmca_dx8[pix_list[corr_gp]], iris[pix_list[corr_gp]])
            correlation_maps_iris[i,2] = correlate(smica_dx8[pix_list[corr_gp]], iris[pix_list[corr_gp]])
            
            skewness_maps[i,0] = skewness(dx8[pix_list[corr_gp]])
            skewness_maps[i,1] = skewness(gmca_dx8[pix_list[corr_gp]])
            skewness_maps[i,2] = skewness(smica_dx8[pix_list[corr_gp]])
            
            kurtosis_maps[i,0] = kurtosis(dx8[pix_list[corr_gp]])
            kurtosis_maps[i,1] = kurtosis(gmca_dx8[pix_list[corr_gp]])
            kurtosis_maps[i,2] = kurtosis(smica_dx8[pix_list[corr_gp]])
        endif
    endfor

    mollview, correlation_maps[*,0], chars=1.5, min=-1, max=1, tit='!6Commander - SDF Correlation', /nest ;, win=1
    mollview, correlation_maps[*,1], chars=1.5, min=-1, max=1, tit='!6GMCA - SDF Correlation', /nest ;, win=2
    mollview, correlation_maps[*,2], chars=1.5, min=-1, max=1, tit='!6Smica - SDF Correlation', /nest ;, win=3
        
;    mollview, skewness_maps[*,0], chars=1.5, min=-1, max=1, tit='Commander Skewness';, win=4
;    mollview, skewness_maps[*,1], chars=1.5, min=-1, max=1, tit='GMCA Skewness'
;    mollview, skewness_maps[*,2], chars=1.5, min=-1, max=1, tit='Smica Skewness'

;    mollview, kurtosis_maps[*,0], chars=1.5, min=-1, max=1, tit='Commander Kurtosis';, win=5
;    mollview, kurtosis_maps[*,1], chars=1.5, min=-1, max=1, tit='GMCA Kurtosis'
;    mollview, kurtosis_maps[*,2], chars=1.5, min=-1, max=1, tit='Smica Kurtosis'

;stop
    cmapr = reorder(correlation_maps,in='NESTED', out='RING')
    cmapr_iris = reorder(correlation_maps_iris,in='NESTED', out='RING')
    smapr = reorder(skewness_maps,in='NESTED', out='RING')
    kmapr = reorder(kurtosis_maps,in='NESTED', out='RING')

    write_tqu, dataset+'_correlation_maps_sdf_'+strtrim(string(low_Nside),2)+'.fits', cmapr, /ring
    write_tqu, dataset+'_correlation_maps_iris_'+strtrim(string(low_Nside),2)+'.fits', cmapr_iris, /ring
    write_tqu, dataset+'_new_skewness_maps_'+strtrim(string(low_Nside),2)+'.fits', smapr, /ring
    write_tqu, dataset+'_new_kurtosis_maps_'+strtrim(string(low_Nside),2)+'.fits', kmapr, /ring

    lab=['Commander','GMCA','Smica']
    for i=0,2 do begin
        
        mollview, correlation_maps[*,i], min=-1, max=1, tit='!6Correlation '+lab[i]+' - SDF',png='pics/ffp4_sdf_correlation_maps_'+lab[i]+'_'+strtrim(string(low_Nside),2)+'_'+sversion+'.png',win=-1, /nest

        mollview, correlation_maps_iris[*,i], min=-1, max=1, tit='!6Correlation '+lab[i]+' - IRIS',png='pics/ffp4_iris_correlation_maps_'+lab[i]+'_'+strtrim(string(low_Nside),2)+'_'+sversion+'.png',win=-1, /nest

        mollview, kurtosis_maps[*,i], min=-1, max=1, tit='!6Kurtosis '+lab[i],png='pics/ffp4_kurtosis_maps_'+lab[i]+'_'+strtrim(string(low_Nside),2)+'_'+sversion+'.png',win=-1, /nest

        mollview, skewness_maps[*,i], min=-1, max=1, tit='!6Skewness '+lab[i],png='pics/ffp4_skewness_maps_'+lab[i]+'_'+strtrim(string(low_Nside),2)+'_'+sversion+'.png',win=-1, /nest
    endfor
    
endif

;# ---------------------------------------------------------------------

if (smooth_corr) then begin

    low_Nside = 64l

    read_fits_map, dataset+'_correlation_maps_'+strtrim(string(low_Nside),2)+'.fits', cmap
    smth = cmap*0.
    for i=0,2 do begin
        ismoothing, cmap[*,i], smap, fwhm_arcmin=240, /ring
        smth[*,i] = smap
        mollview, smth[*,i], min=-1, max=1
    endfor

    write_tqu, dataset+'_correlation_maps_'+strtrim(string(low_Nside),2)+'_smth60.fits', smth, /ring

endif

;#
;---------------------------------------------------------------------

if (internal_corr) then begin

    Nside = 2048l
    Npix = 12l*Nside^2
    nComp = 4
    solution = fltarr(Npix, nComp)

    components = ['cmb', 'fg_co_amp', 'fg_dust_amp', 'fg_lowfreq_amp']

    sversion = 'v3.1'

;##    if (dataset eq 'ffp4') then root='/global/scratch/sd/dpietrob/commander-ruler/ffp4/v4/ffp4_commander-ruler_v4_cmb_uK.fits'
    if (dataset eq 'ffp4') then root='/global/scratch/sd/dpietrob/commander-ruler/ffp4/'+sversion+'/ffp4_commander-ruler_'+sversion+'_'
    if (dataset eq 'dx8') then root='/global/scratch/sd/loris/CompSep/GLSS/DX/8/Davide_new/pix_raw_effreq_scan_vtd_5deg_ps/dx8_7b_mondip_fs_avrg_'

    gmask = fltarr(Npix)+1.
    for icomp=0,ncomp-1 do begin
        read_fits_map, root+components[icomp]+'_uK.fits', map
        bp = where( finite(map, /nan) eq True )
        map[bp] = 0.
        in = reorder(map, in='ring', out='nest')
        solution[*,icomp] = in
;##        mollview, in, chars=1.5, win=icomp+10, /hist, /nest
        gmask[bp] = 0.
    endfor
    
    low_Nside = 64l
    low_Npix = nside2npix(low_Nside)

    correlation_maps = fltarr(low_Npix,6)
    
    correlation_maps[*,*] = -1.6375e30
    
    pix2vec_nest, low_Nside, lindgen(low_Npix), vec, vertex

    for i=0, low_Npix-1 do begin
        
        if ( (i/1000)*1000 eq i) then print, i, ' / ', low_Npix, ' --> '+string(float(i)/low_Npix*100,format='(f4.1)')+'%'

        query_polygon, Nside, transpose( reform( vertex[i,*,*]) ), pix_list, inclusive=True, /nested
        
        corr_gp = where(gmask[pix_list] gt 0.)
        
        if (corr_gp[0] ne -1) then begin
; cmb-co
            correlation_maps[i,0] = correlate(solution[pix_list[corr_gp],0], solution[pix_list[corr_gp],1])
; cmb-dust
            correlation_maps[i,1] = correlate(solution[pix_list[corr_gp],0], solution[pix_list[corr_gp],2])
; cmb-synch
            correlation_maps[i,2] = correlate(solution[pix_list[corr_gp],0], solution[pix_list[corr_gp],3])
; co-dust
            correlation_maps[i,3] = correlate(solution[pix_list[corr_gp],1], solution[pix_list[corr_gp],2])
; co-synch
            correlation_maps[i,4] = correlate(solution[pix_list[corr_gp],1], solution[pix_list[corr_gp],3])
; dust_synch
            correlation_maps[i,5] = correlate(solution[pix_list[corr_gp],2], solution[pix_list[corr_gp],3])

        endif    

    endfor

    clab = strarr(6)
    icnt = 0
    for i=0,3 do begin
        for j=i+1,3 do begin
            clab[icnt] = components[i]+'-'+components[j]
            out = reorder(correlation_maps[*,icnt], in='nest', out='ring')
            write_fits_map, dataset+'_cross_commander_correlation_'+clab[icnt]+'_'+strtrim(string(low_Nside),2)+'_'+sversion+'.fits', out, /ring
            mollview, out, tit='Correlation '+clab[icnt], min=-1, max=1, chars=1.5, win=icnt+1
            mollview, out, tit='Correlation '+clab[icnt], min=-1, max=1, chars=1.5, png='pics/'+dataset+'_cross_commander_correlation_'+clab[icnt]+'_'+strtrim(string(low_Nside),2)+'_'+sversion+'.png', win=-1
            icnt = icnt + 1
            print, clab
        endfor
    endfor


endif

;# ---------------------------------------------------------------------

if (compare_old_new) then begin
    iris_file = 'IRIS_nohole_4_2048.fits'
    read_fits_map, iris_file, dust
;##    read_fits_map, '/global/scratch/sd/loris/CompSep/GLSS/FFP/4/7b/avrg_maps/ffp4_7b_mondip_avrg_cmb.fits', old
    read_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/ffp4/v3/ffp4_commander-ruler_v3_cmb_uK.fits', old
;##    read_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/ffp4/v4/ffp4_commander-ruler_v4_cmb_uK.fits', new
    read_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/ffp4/v3.1/ffp4_commander-ruler_v3.1_cmb_uK.fits', new

    bp = where( (finite(old, /nan) eq True) or (finite(new, /nan) eq True) or (dust lt 0.) ) 

    gmask = new*0.+1
    gmask[bp] = 0.

    old[bp] = 0.
    new[bp] = 0.
    dust[bp] = 0.

    dustr = reorder(dust, in='ring', out='nest')
    oldr = reorder(old, in='ring', out='nest')
    newr = reorder(new, in='ring', out='nest')

    dust = dustr
    dustr = 0.
    old = oldr
    oldr = 0.
    new = newr
    newr = 0. 

;##    c = correlate_maps_nest(oldr, newr, 64)
    low_Nside = 64l
    low_Npix = nside2npix(low_Nside)

    correlation_maps = fltarr(low_Npix,3)
    correlation_maps[*,*] = -1.6375e30
    
    pix2vec_nest, low_Nside, lindgen(low_Npix), vec, vertex

    for i=0, low_Npix-1 do begin
        if ( (i/1000)*1000 eq i) then print, i, ' / ', low_Npix, ' --> '+string(float(i)/low_Npix*100,format='(f4.1)')+'%'
        query_polygon, Nside, transpose( reform( vertex[i,*,*]) ), pix_list, inclusive=True, /nested
        corr_gp = where(gmask[pix_list] gt 0.)        
        if (corr_gp[0] ne -1) then begin
            correlation_maps[i,0] = correlate(old[pix_list[corr_gp]], new[pix_list[corr_gp]])
            correlation_maps[i,1] = correlate(old[pix_list[corr_gp]], dust[pix_list[corr_gp]])
            correlation_maps[i,2] = correlate(new[pix_list[corr_gp]], dust[pix_list[corr_gp]])
        endif
    endfor
    mollview, correlation_maps[*,0], /nest, chars=1.5, tit='!6v3 - v3.1 Correlation Map', win=21
    mollview, correlation_maps[*,1], min=-1, max=1, /nest, chars=1.5, tit='!6v3 - IRIS Correlation Map', win=22
    mollview, correlation_maps[*,2], min=-1, max=1, /nest, chars=1.5, tit='!6v3.1 - IRIS Correlation Map', win=23
    mollview, correlation_maps[*,1]-correlation_maps[*,2], min=-.25, max=.25, /nest, chars=1.5, tit='!6v3 - v3.1 IRIS Correlation', win=24

    mollview, correlation_maps[*,0], /nest, chars=1.5, tit='!6v3 - v3.1 Correlation Map', win=-1, png='pics/correlation_v3-v3.1.png'
    mollview, correlation_maps[*,1], min=-1, max=1, /nest, chars=1.5, tit='!6v3 - IRIS Correlation Map', win=-1, png='pics/correlation_v3-iris.png'
    mollview, correlation_maps[*,2], min=-1, max=1, /nest, chars=1.5, tit='!6v3.1 - IRIS Correlation Map', win=-1, png='pics/correlation_v3.1-iris.png'
    mollview, correlation_maps[*,1]-correlation_maps[*,2], min=-.25, max=.25, /nest, chars=1.5, tit='!6v3 - v3.1 IRIS Correlation', win=-1, png='pics/correlation_v3_v3.1_diff.png'
endif

;# --------------------------------------------------------------------
if (compute_high_res_chisq_31) then begin
;# Warning: Mean maps only
    freq = [ 28.4, 44.1, 70.3, 142.7, 221.9, 101.2, 360.6 ]
    sfreq = ['030','044','070','143','217','100', '353']
    nfreq = n_elements(freq)

    co_sed = [0, 0, 0, 0, 0.33, 1., 0.109] ;vtd

    restore, '../ffp4/maps/ns0256/chains/pix_raw_effreq_scan_vtd_5deg_ps.res.sav'

    monodip = total(foreg,1)/nsample

    print, monodip

;stop

    Nside = 2048l
    Npix = 12l*Nside^2
    nComp = 4
    solution = fltarr(Npix, nComp)

    components = ['cmb', 'fg_co_amp', 'fg_dust_amp', 'fg_lowfreq_amp']

    if (dataset eq 'ffp4') then root='/global/scratch/sd/dpietrob/commander-ruler/ffp4/v3.1/ffp4_commander-ruler_v3.1_'
    if (dataset eq 'dx8') then root='/global/scratch/sd/loris/CompSep/GLSS/DX/8/Davide_new/pix_raw_effreq_scan_vtd_5deg_ps/dx8_7b_mondip_fs_avrg_'

    gmask = fltarr(Npix)+1.
    print, ' - Reading solutions...'
    for icomp=0,ncomp-1 do begin
        read_fits_map, root+components[icomp]+'_uK.fits', map
        bp = where( (finite(map, /nan) eq True) or (map eq -1.6375e30) )
        map[bp] = 0.
;##        in = reorder(map, in='ring', out='nest')
        solution[*,icomp] = map
;##        mollview, in, chars=1.5, win=icomp+10, /hist, /nest
        gmask[bp] = 0.
    endfor

    print, ' - Reading free-free index...'
    read_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/ffp4/v3.1/ffp4_commander-ruler_v3.1_fg_lowfreq_mean_index.fits', lfc_indx
    lfc_curv = 0.
    read_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/ffp4/v3.1/ffp4_commander-ruler_v3.1_fg_dust_mean_emissivity.fits', dust_emis
    read_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/ffp4/v3.1/ffp4_commander-ruler_v3.1_fg_dust_mean_temperature.fits', dust_temp

    lfc_pivot = 30.
    dust_pivot = 353.

    cf = conversionfactor(freq, /antenna2thermo)

    h = 6.626068 * 10.^(-34)
    k = 1.3806503 * 10.^(-23)
    c = 299792458.

    icmb  = 0
    ico   = 1
    idust = 2
    isync = 3

    chisq = fltarr(Npix)
    
    ff = 0
    lf = 6
    
    for ifreq=ff,lf do begin
        print, sfreq[ifreq]
        read_fits_map, '/global/scratch/sd/loris/CompSep/GLSS/FFP/4/uKInputs/challenge01/map01_'+sfreq[ifreq]+'.fits', inmap
        read_fits_map, '/global/scratch/sd/dpietrob/ffp4/maps/ns2048/ffp4_rms_'+sfreq[ifreq]+'_uK.fits', rms

        dip_vec = reform( monodip[1:3,ifreq] )
        dipole = make_dipole( dip_vec, Nside)
;mollview, dipole
;stop
        x = h*1.d9/k/dust_temp
        bb  = 1. / ( exp(x*freq[ifreq])-1. )
        bb0 = 1. / ( exp(x*dust_pivot)-1. )
        fit = solution[*,icmb] + monodip[0,ifreq] + dipole + cf[ifreq] * (  solution[*,idust] * (freq[ifreq]/dust_pivot)^(dust_emis+1.)*bb/bb0 $
                                              + solution[*,isync]  * (freq[ifreq]/lfc_pivot)^(lfc_indx+lfc_curv*(freq[ifreq]/lfc_pivot)) $
                                              + solution[*,ico] * CO_sed[ifreq] )

;##        mollview, fit, min=-500, max=500, win=3
;        mollview, inmap - fit, min=-300, max=300, win=2, tit='!6Ffp4 '+sfreq[ifreq]+' Residuals', units='!7l!6K CMB'
;        mollview, inmap - fit, min=-300, max=300, win=-1, tit='!6Ffp4 '+sfreq[ifreq]+' Residuals', units='!7l!6K CMB', png='pics/residuals_'+sfreq[ifreq]+'_v3.1.png'

        mollview, (inmap-fit)/rms, min=-3, max=3, win=3, chars=1.5, units='!7r', tit='!6FFP4 '+sfreq[ifreq]+' Residuals'
        mollview, (inmap-fit)/rms, min=-3, max=3, win=-3, chars=1.5, units='!7r', tit='!6FFP4 '+sfreq[ifreq]+' Residuals', png='pics/residuals_'+sfreq[ifreq]+'_v3.1.png'
        
        chisq = chisq + ((inmap-fit)/rms)^2

;stop
    endfor

    mollview, chisq, max=15, tit='!7v!u!62!n High Resolution', chars=1.5
    mollview, chisq, max=15, tit='!7v!u!62!n High Resolution', chars=1.5,win=-1, png='pics/chisq_v3.1.png'

    ud_grade, chisq, chisqdg, nside_out=256, order_in='ring', order_out='ring'

    mollview, chisqdg, max=15, tit='!7v!u!62!n High Resolution (down to 256)', chars=1.5
    mollview, chisqdg, max=15, tit='!7v!u!62!n High Resolution (down to 256)', chars=1.5,win=-1, png='pics/chisq_downg_v3.1.png'

endif

;# --------------------------------------------------------------------
if (compute_high_res_chisq) then begin
    freq = [ 28.44, 44.12, 70.33, 101.28, 142.85, 222.51, 361.46 ]
    sfreq = ['030','044','070','100','143','217','353']
    nfreq = n_elements(freq)

;##    co_sed = [0, 0, 0, 1, 0, 0.720469820083658, 0.116895820776374] ;hke
    co_sed = [0, 0, 0, 1, 0, 0.33, 0.109] ;vtd


    dust_sed = [1.16221264816912, $
                0.306775144929560, $
                0.133584244149504, $
                0.148429290375434, $     
                0.315703299477422, $   
                0.566630763761037, $     
                1.00009001093963 ]

    monodip = [$
                [140.7941,       1.017129,       4.882650,       8.552243], $    
                [94.80920,      8.4179789E-02,   2.850955,       8.657156], $    
                [76.19675,      0.9259039,       3.521682,       8.915731], $    
                [94.19647,       1.939186,       3.858565,       9.081417], $    
                [116.8398,       1.237141,       4.244828,       8.999939], $    
                [120.8427,       1.075501,       4.304761,       8.963346], $    
                [105.0304,       1.515763,       4.079051,       8.923319    ] ]

    monodip = transpose(monodip)

    help, monodip

;stop

    Nside = 2048l
    Npix = 12l*Nside^2
    nComp = 4
    solution = fltarr(Npix, nComp)

    components = ['cmb', 'co', 'dust', 'synch']

    if (dataset eq 'ffp4') then root='/global/scratch/sd/dpietrob/commander-ruler/ffp4/v4/ffp4_7b_mondip_fs_HK_'
    if (dataset eq 'dx8') then root='/global/scratch/sd/loris/CompSep/GLSS/DX/8/Davide_new/pix_raw_effreq_scan_vtd_5deg_ps/dx8_7b_mondip_fs_avrg_'

    gmask = fltarr(Npix)+1.
    print, ' - Reading solutions...'
    for icomp=0,ncomp-1 do begin
        read_fits_map, root+components[icomp]+'.fits', map
        bp = where( finite(map, /nan) eq True )
        map[bp] = 0.
;##        in = reorder(map, in='ring', out='nest')
        solution[*,icomp] = map
;##        mollview, in, chars=1.5, win=icomp+10, /hist, /nest
        gmask[bp] = 0.
    endfor

    print, ' - Reading free-free index...'
    read_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/ffp4/v4/synch_index_ns2048.fits', sync_index
    ff_pivot = 30.

    cf = conversionfactor(freq, /antenna2thermo)

    icmb  = 0
    ico   = 1
    idust = 2
    isync = 3

    chisq = fltarr(Npix)

    for ifreq=0,nfreq-1 do begin
        print, sfreq[ifreq]
        read_fits_map, '/global/scratch/sd/loris/CompSep/GLSS/FFP/4/uKInputs/challenge01/map01_'+sfreq[ifreq]+'.fits', inmap
        read_fits_map, '/global/scratch/sd/dpietrob/ffp4/maps/ns2048/ffp4_rms_'+sfreq[ifreq]+'_uK.fits', rms

        dip_vec = reform( monodip[ifreq, 1:3] )
        dx = dip_vec[0]
        dy = dip_vec[1]
        dz = dip_vec[2]
        dipole = make_dipole( dip_vec, Nside)
;mollview, dipole
;stop
        fit = solution[*,icmb] + monodip[ifreq,0] + dipole + cf[ifreq]*( solution[*,ico]*co_sed[ifreq] + solution[*,idust]*dust_sed[ifreq] + solution[*,isync]*(freq[ifreq]/ff_pivot)^sync_index )
        
;        write_fits_map, ''

;        mollview, fit, min=-300, max=300, win=1

;        mollview, inmap - fit, min=-30, max=30, win=2

;        mollview, (inmap-fit)/rms, min=-3, max=3, win=3
        
        chisq = chisq + ((inmap-fit)/rms)^2

;stop
    endfor

    mollview, chisq, max=15, tit='!7v!u!62!n High Resolution', chars=1.5

    ud_grade, chisq, chisqdg, nside_out=256, order_in='ring', order_out='ring'

    mollview, chisqdg, max=15, tit='!7v!u!62!n High Resolution (down to 256)', chars=1.5

endif

;# --------------------------------------------------------------------

if (check_patch) then begin
    read_fits_map, 'ffp4_correlation_maps_sdf_8.fits', corr_sdf
    read_fits_map, 'ffp4_correlation_maps_iris_8.fits', corr_iris

    read_fits_map, 'ffp4_new_skewness_maps_8.fits', skew_maps

    read_fits_map, 'ffp4_new_kurtosis_maps_8.fits', kurt_maps

    m = fltarr(12l*8l^2)
    hot_pix = [40, 58, 80, 272, 620]
    m[hot_pix] = corr_sdf[hot_pix,0]

;##    mollview, m, min=-1, max=1, win=20, tit='!6Regions'

    pix2ang_ring, 8, hot_pix, theta, phi
    
;[ [60.73, 10.79],[60.73,334.82], [54.63,307.02], [14.19,0.94], [-35.08,316.34] ]

    theta = 90.-theta*!radeg
    phi = phi * !radeg

    Nside = 2048l
    Npix = 12l*Nside^2
    nComp = 4
    solution = fltarr(Npix, nComp)

    components = ['cmb', 'co', 'dust', 'synch']

    if (dataset eq 'ffp4') then root='/global/scratch/sd/dpietrob/commander-ruler/ffp4/v4/ffp4_7b_mondip_fs_HK_'
    if (dataset eq 'dx8') then root='/global/scratch/sd/loris/CompSep/GLSS/DX/8/Davide_new/pix_raw_effreq_scan_vtd_5deg_ps/dx8_7b_mondip_fs_avrg_'

    gmask = fltarr(Npix)+1.
    print, ' - Reading solutions...'
    for icomp=0,ncomp-1 do begin
        read_fits_map, root+components[icomp]+'.fits', map
        bp = where( finite(map, /nan) eq True )
        map[bp] = 0.
        if (icomp eq 0) then remove_dipole, map, nside=2048, ordering='ring', gal_cut=5 
        solution[*,icomp] = map
        gmask[bp] = 0.
    endfor

    read_fits_map, 'SFD_i100_healpix_2048.fits', sdf

    read_fits_map,  'IRIS_nohole_4_2048.fits', iris
    iris[where(iris lt 0)] = -1.6375e30


    icmb  = 0
    ico   = 1
    idust = 2
    isync = 3

    if (True) then begin
        write_fits_map, 'tmp_map_up.fits', solution[*,icmb], /ring
        ianafast, 'tmp_map_up.fits', 'tmp_ind_up_cls.fits', alm1_out='tmp_ind_up_alms.fits', simul_type=1, iter=1, nlmax=2500, /silent
        isynfast, 'tmp_ind_up_cls.fits', map, nside=2048, alm_in='tmp_ind_up_alms.fits', nlmax=2500, simul_type=1, /silent
        spawn, 'rm tmp_ind_up_cls.fits tmp_ind_up_alms.fits tmp_map_up.fits'
        cmb3000 = map
        write_fits_map, dataset+'_hke_cmb_l3000.fits', cmb3000, /ring, units='!7l!6K CMB'
    endif else begin
        read_fits_map, dataset+'_hke_cmb_l3000.fits', cmb3000
    endelse

    read_fits_map, '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/gmca/'+dataset+'_gmca_v1_cmb.fits', gmca_dx8, order=dpord
;print, dpord
;stop
    if ( (dpord eq 'nest') or (dpord eq 'NESTED') ) then begin
        print, ' - Reordering...'
        gmca_dx8r = reorder(gmca_dx8, in='nest', out='ring')
        gmca_dx8 = gmca_dx8r*1.e6
        gmca_dx8r = 0.
        remove_dipole, gmca_dx8, nside=2048, ordering='ring', gal_cut=5        
    endif

    read_fits_map, '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/smica/'+dataset+'_smica_needlet_cmb.fits', smica_dx8, order=dpord
    remove_dipole, smica_dx8, nside=2048, ordering='ring', gal_cut=5        

    for i=0,4 do begin
;        gnomview, corr_sdf[*,0], min=-1, max=1, tit='!6Correlation Commander - SDF', rot=[phi[i], theta[i]], reso=4, win=1, grat=[5,5]
;        gnomview, corr_sdf[*,1], min=-1, max=1, tit='!6Correlation GMCA - SDF', rot=[phi[i], theta[i]], reso=4, win=2, grat=[5,5]
;        gnomview, corr_sdf[*,2], min=-1, max=1, tit='!6Correlation Smica - SDF', rot=[phi[i], theta[i]], reso=4, win=3, grat=[5,5]

;        gnomview, corr_iris[*,0], min=-1, max=1, tit='!6Correlation Commander - IRIS', rot=[phi[i], theta[i]], reso=4, win=4, grat=[5,5]
;        gnomview, corr_iris[*,1], min=-1, max=1, tit='!6Correlation GMCA - IRIS', rot=[phi[i], theta[i]], reso=4, win=5, grat=[5,5]
;        gnomview, corr_iris[*,2], min=-1, max=1, tit='!6Correlation Smica - IRIS', rot=[phi[i], theta[i]], reso=4, win=6, grat=[5,5]

;        gnomview, skew_maps[*,0], min=-1, max=1, tit='!6Commander Skewness', rot=[phi[i], theta[i]], reso=4, win=4, grat=[5,5]
;        gnomview, kurt_maps[*,0], min=-1, max=1, tit='!6Commander Kurtosis', rot=[phi[i], theta[i]], reso=4, win=4, grat=[5,5]

;        gnomview, cmb3000, min=-300, max=300, rot=[phi[i], theta[i]], tit='!6Commander CMB', reso=4, win=7, grat=[5,5]
;        gnomview, gmca_dx8, min=-300, max=300, rot=[phi[i], theta[i]], tit='!6GMCA CMB', reso=4, win=8, grat=[5,5]
;        gnomview, smica_dx8, min=-300, max=300, rot=[phi[i], theta[i]], tit='!6Smica CMB', reso=4, win=9, grat=[5,5]

;        gnomview, cmb3000-gmca_dx8, min=-30, max=30, rot=[phi[i], theta[i]], tit='!6Commander-GMCA', reso=4, win=12, grat=[5,5]
;        gnomview, cmb3000-smica_dx8, min=-30, max=30, rot=[phi[i], theta[i]], tit='!6Commander-Smica', reso=4, win=13, grat=[5,5]
;        gnomview, smica_dx8-gmca_dx8, min=-30, max=30, rot=[phi[i], theta[i]], tit='!6Smica-GMCA', reso=4, win=14, grat=[5,5]

;        gnomview, sdf, rot=[phi[i], theta[i]], tit='!6SDF', reso=4, win=10, grat=[5,5], /hist
;        gnomview, iris, rot=[phi[i], theta[i]], tit='!6IRIS', reso=4, win=11, grat=[5,5], /hist

;# ---
        gnomview, corr_iris[*,0], min=-1, max=1, tit='!6Correlation Commander - IRIS', rot=[phi[i], theta[i]], reso=4, win=-1, grat=[5,5], png='pics/patch_01_corr-iris_'+strtrim(string(i+1),2)+'.png'

        gnomview, cmb3000, min=-300, max=300, rot=[phi[i], theta[i]], tit='!6Commander CMB', reso=4, win=-1, grat=[5,5], png='pics/patch_02_commander_'+strtrim(string(i+1),2)+'.png'
        gnomview, gmca_dx8, min=-300, max=300, rot=[phi[i], theta[i]], tit='!6GMCA CMB', reso=4, win=-1, grat=[5,5], png='pics/patch_03_gmca_'+strtrim(string(i+1),2)+'.png'
        gnomview, smica_dx8, min=-300, max=300, rot=[phi[i], theta[i]], tit='!6Smica CMB', reso=4, win=-1, grat=[5,5], png='pics/patch_04_smica_'+strtrim(string(i+1),2)+'.png'

        gnomview, cmb3000-gmca_dx8, min=-30, max=30, rot=[phi[i], theta[i]], tit='!6Commander-GMCA', reso=4, win=-1, grat=[5,5], png='pics/patch_diff_commander-gmca_'+strtrim(string(i+1),2)+'.png'
        gnomview, cmb3000-smica_dx8, min=-30, max=30, rot=[phi[i], theta[i]], tit='!6Commander-Smica', reso=4, win=-1, grat=[5,5], png='pics/patch_diff_commander-smica_'+strtrim(string(i+1),2)+'.png'
        gnomview, smica_dx8-gmca_dx8, min=-30, max=30, rot=[phi[i], theta[i]], tit='!6Smica-GMCA', reso=4, win=-1, grat=[5,5], png='pics/patch_diff_smica-gmca_'+strtrim(string(i+1),2)+'.png'

        gnomview, sdf, rot=[phi[i], theta[i]], tit='!6SDF', reso=4, win=-1, grat=[5,5], /hist, png='pics/patch_05_sdf_'+strtrim(string(i+1),2)+'.png'
        gnomview, iris, rot=[phi[i], theta[i]], tit='!6IRIS', reso=4, win=-1, grat=[5,5], /hist, png='pics/patch_06_iris_'+strtrim(string(i+1),2)+'.png'

        spawn, 'montage pics/patch_*_'+strtrim(string(i+1),2)+'.png -tile 6x -adjoin -geometry x300\>+1+1 -background none pics/strip_'+strtrim(string(i+1),2)+'.png'

        spawn, 'montage pics/patch_diff_*_'+strtrim(string(i+1),2)+'.png -tile 3x -adjoin -geometry x300\>+1+1 -background none pics/strip_diff_'+strtrim(string(i+1),2)+'.png'

;stop
    endfor

endif

;# --------------------------------------------------------------------

if (check_spectra) then begin
    ianafast, '/global/scratch/sd/dpietrob/commander-ruler/ffp4/v4/ffp4_commander-ruler_v4_cmb_lowres_uK.fits', cls_low, nlmax=750, theta_cut_deg=20, regression=2
    ianafast, '/global/scratch/sd/dpietrob/commander-ruler/ffp4/v4/ffp4_commander-ruler_v4_cmb_uK.fits', cls_high, nlmax=750, theta_cut_deg=20, regression=2

    fits2cl, tf, 'ffp4/v4/ffp4_commander-ruler_v4_cmb_beam.fits'

    window, 3 & plot, cls_low/gaussbeam(40.,750)^2/cls_high*tf^2/healpixwindow(256)^2*healpixwindow(2048)^2, yr=[0.8,1.2]

endif

;# --------------------------------------------------------------------

if (compare_co) then begin

    read_fits_map, 'ffp4/v3.1/ffp4_commander-ruler_v3.1_fg_co_amp_uK.fits', new
    read_fits_map, 'ffp4/v3/ffp4_commander-ruler_v3_co_amp_uK.fits', old
    read_fits_map, 'ffp4/v4/ffp4_commander-ruler_v4_fg_co_amp_uK.fits', hke
    read_fits_map, '../magique_inputs/diffuse/co/co_ampl10.fits', inco

    incf =  0.0041167480
    fact = 2.6e-5
    t2a = conversionfactor(101.28,/thermo2antenna)
    
    inco = inco * incf * fact * t2a

    cut = make_sky_cut(25.,2048)
    gp = where(cut eq 0.)

    window, 22, xsize=720*1.2, ysize=720*1.2
    plot, old[gp[0:*:100]], new[gp[0:*:100]], psym=3, chars=1.5, xtit='!7l!6K Antenna', ytit='!7l!6K Antenna', xr=[0,1000], yr=[0,1000]
    oplot, old[gp[0:*:100]], hke[gp[0:*:100]], psym=3, col=70
    oplot, new[gp[0:*:100]], hke[gp[0:*:100]], psym=3, col=210
    oplot, [0,2.e4], [0,2.e4], col=245
    xyouts, 50,950, '!6Constant T!dd!n CO vs. Varying T!dd!n CO (v3.1)', chars=1.5
    xyouts, 50,900, '!6Constant T!dd!n CO vs. Effective Dust CO (v4)', col=70, chars=1.5
    xyouts, 50,850, '!6Varying T!dd!n CO vs. Effective Dust CO', col=210, chars=1.5
    write_png, 'pics/co_comparison_v3.1.png', tvrd(/true)

    window, 23, xsize=720*1.2, ysize=720*1.2
    plot, inco[gp[0:*:100]], new[gp[0:*:100]], psym=3, chars=1.5, xtit='!6Input CO 1-0 (x 2.3e-5)', ytit='!6Varying T!dd!n CO', xr=[0,1000], yr=[0,1000]
    oplot, [0,2.e4], [0,2.e4], col=245
    write_png, 'pics/co_input-comparison_v3.1.png', tvrd(/true)

endif

;# --------------------------------------------------------------------

if (avrg_maps) then begin
;# From Loris'
;#    fg = 'co'

    afg = ['co','dust', 'ff']
    set = ['hr1', 'hr2']
  
    for ifg=0,2 do begin
        for iset=0,1 do begin

            fg = afg[ifg]

    dir = '/global/scratch/sd/loris/CompSep/GLSS/FFP/4/Davide_new/pix_raw_effreq_scan_vtd_5deg_ps/'

    root = 'ffp4_7b_mondip_'+set[iset]+'_'

    nsample = 34
    nchain  = 1
    ntot    = nchain*nsample
    nfreq = 7
    ns = 2048l
    np = ns*ns*12l

    tot = dblarr(np)
    tot[*] = 0.d0

    for i=1,nchain do begin
        cstring = 'c000'+string(i,'(i1)')  
        cx = dir+'c00'+string(i,'(i1)')  
        for j=1,nsample do begin
            sstring = 'k00'+string(j,'(i3.3)')
            print,cx,' ',sstring
        
            name = cx+'/'+root+cstring+'_'+sstring+'_'+fg+'.fits'
            read_fits_map,name,tmp
            tot = tot+tmp
        end     
    end

    tot = tot/ntot

    good = tot gt -1.d10
    id   = where(good lt 0.5)

    tot[id] = !healpix.bad_value

    name = 'ffp4/v3.1/'+root+'avrg_'+fg+'.fits'
    write_fits_map,name,tot,/ring

endfor
endfor

        
endif


   STOP


END
