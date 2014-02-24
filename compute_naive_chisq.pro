   !path = !path+':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'
   !path = !path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'
   
   False = 0b
   True  = 1b

   dataset = 'dx9'

   if dataset eq 'ffp4' then freq = [28.44, 44.12, 70.33, 101.28, 142.85, 222.51, 361.56]*1.d0 ;46]*1.d0
   if dataset eq 'dx9' then freq = [28.4, 44.1, 70.3, 101.2, 142.6, 221.9, 360.6]*1.d0
   sfreq = ['030', '044', '070', '100', '143', '217', '353']

   wfreq = [23, 33, 41, 61, 94]*1.d0
   wsfreq = ['K','Ka','Q','V','W']
   nwfreq = n_elements(wfreq)

   nfreq = n_elements(freq)
   ns = 128l ;2048l
   npix = 12l*ns^2
   if dataset eq 'ffp4' then  co_spec= [0., 0., 0., 1., 0., 0.65, 0.5]*1.d0
   if dataset eq 'dx9' then co_spec= [0., 0., 0., 1., 0., .674, .259]*1.d0

   do_samples = True

   if (do_samples) then begin
; chi2 for samples in a chain
       if dataset eq 'ffp4' then begin
           root = 'pix_01a_v3'
;##           root = 'pix_01a_v3_ff'
;##           root = 'pix_01a_v3_wmap'
           maskfile = '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/dust_1.5mK_mask_ns0128.fits'
;##           maskfile = '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/dust_0.25mK_mask_ns0128.fits'
           dust_pivot = 361.46
           lowf_pivot = 28.44
           ff_pivot = 70.33
       endif
           
       if dataset eq 'dx9' then begin
;##           root = 'pix_7b'
;##           root = 'pix_7b_2pl_lowf'
;##           root = 'pix_wmap' ; hke-like
           root = 'pix_wmap_ff'
           maskfile = '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/dust_1.5mK_mask_ns0128.fits'
           dust_pivot = 360.6
           lowf_pivot = 28.4
           ff_pivot = 70.3
       endif

       read_fits_map, maskfile, mask
       bp = where(mask eq 0.)

       if dataset eq 'ffp4' then restore, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'.res.sav'
       if dataset eq 'dx9' then restore, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'.res.sav'

       ich = 1
       sch = string(ich, format='(i4.4)')
       isample = 150

       chisq = fltarr(npix, nsample)

    for isamp=isample,isample do begin
        ssamp = string(isamp, format='(i5.5)')
        print, ssamp
    
       foreg = fltarr(npix, nfreq)
       res = fltarr(npix, nfreq)

        read_fits_map, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'/fg_amp_map_no01_c'+sch+'_k'+ssamp+'.fits', dust
        if dataset eq 'dx9' then begin
            read_fits_map, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'/fg_amp_map_no03_c'+sch+'_k'+ssamp+'.fits', co_map
            read_fits_map, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'/fg_amp_map_no02_c'+sch+'_k'+ssamp+'.fits', lowf
        endif
        if dataset eq 'ffp4' then begin
            read_fits_map, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'/fg_amp_map_no02_c'+sch+'_k'+ssamp+'.fits', co_map
            read_fits_map, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'/fg_amp_map_no03_c'+sch+'_k'+ssamp+'.fits', lowf
        endif
        read_fits_map, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'/fg_amp_map_no04_c'+sch+'_k'+ssamp+'.fits', cmb

        if ((root eq 'pix_01a_v3_ff') or (root eq 'pix_7b_2pl_lowf') or (root eq 'pix_wmap_ff')) $
          then read_fits_map, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'/fg_amp_map_no05_c'+sch+'_k'+ssamp+'.fits', ff


        read_fits_map, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'/fg_ind_map_no01_c0001_k'+ssamp+'.fits', emis
        read_fits_map, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'/fg_ind_map_no02_c0001_k'+ssamp+'.fits', temp
        read_fits_map, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'/fg_ind_map_no03_c0001_k'+ssamp+'.fits', beta

        readcol, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'/foreground_c0001_k'+string(isamp-1, format='(i5.5)')+'.dat', m, dx, dy, dz

       cc = cmb
       remove_dipole, cc, mask, nside=ns, ordering='ring'
       mollview, cc, min=-300, max=300, chars=1.5, tit='!6CMB', px=600, png=dataset+'_'+root+'_cmb.png', win=-1
       mollview, cc, min=-300, max=300, chars=1.5, tit='!6CMB', px=600, win=1

       read_fits_map, '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/chains/'+root+'/chisq_c'+sch+'_k'+ssamp+'.fits', glb_chi2
       glb_chi2[where(glb_chi2 eq 0.)] = -1.6375e30
       mollview, glb_chi2, max=21, chars=1.5, tit='!7v!6!u2!n', px=600, png=dataset+'_'+root+'_chisq.png', win=-1
       mollview, glb_chi2, max=21, chars=1.5, tit='!7v!6!u2!n', px=600, win=2

       mollview, dust, min=0, max=400, chars=1.5, units='!7l!6K RJ', tit='!6Thermal Dust Solution', px=600, png=dataset+'_'+root+'_dust.png', win=-1, no_monopole=1, gal_cut=80
       mollview, co_map, min=0, max=400, chars=1.5, units='!7l!6K RJ', tit='!6CO Solution', px=600, png=dataset+'_'+root+'_co.png', win=-1, no_monopole=1, gal_cut=80
       mollview, lowf, min=0, max=400, chars=1.5, units='!7l!6K RJ', tit='!6Low-freq Solution', px=600, png=dataset+'_'+root+'_lowf.png', win=-1, no_monopole=1, gal_cut=80

       mollview, dust, min=0, max=400, chars=1.5, units='!7l!6K RJ', tit='!6Thermal Dust Solution', px=600, win=3, no_monopole=1, gal_cut=80
       mollview, co_map, min=0, max=400, chars=1.5, units='!7l!6K RJ', tit='!6CO Solution', px=600, win=4, no_monopole=1, gal_cut=80
       mollview, lowf, min=0, max=400, chars=1.5, units='!7l!6K RJ', tit='!6Low-freq Solution', px=600, win=5, no_monopole=1, gal_cut=80

       if ((root eq 'pix_01a_v3_ff') or (root eq 'pix_7b_2pl_lowf') or (root eq 'pix_wmap_ff')) then begin
           mollview, ff, min=0, max=400, chars=1.5, units='!7l!6K RJ', tit='!6Ff Solution', px=600, png=dataset+'_'+root+'_f-f.png', win=-1, no_monopole=1, gal_cut=80
           mollview, ff, min=0, max=400, chars=1.5, units='!7l!6K RJ', tit='!6Ff Solution', px=600, win=6, no_monopole=1, gal_cut=80
       endif

       mollview, emis, min=1, max=3, chars=1.5, units='', tit='!6Thermal Dust Emissivity', px=600, png=dataset+'_'+root+'_emis.png', win=-1
       mollview, temp, min=16, max=20, chars=1.5, units='', tit='!6Thermal Dust Temperature', px=600, png=dataset+'_'+root+'_temp.png', win=-1
       mollview, beta, min=-4, max=-2, chars=1.5, units='', tit='!6Low-freq Index', px=600, png=dataset+'_'+root+'_beta.png', win=-1

;stop
        dx = dx * sqrt(3.)
        dy = dy * sqrt(3.)
        dz = dz * sqrt(3.)

        h = 6.626068 * 10.^(-34)
        k = 1.3806503 * 10.^(-23)
        c = 2.99792458d8
        
        for ifreq=0,nfreq-1 do begin

            print, sfreq[ifreq]

            if dataset eq 'dx9' then file = '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/'+dataset+'_Imap_'+sfreq[ifreq]+'_ns128_uK.fits'
            if dataset eq 'ffp4' then file = '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/'+dataset+'_01a_'+sfreq[ifreq]+'_ns128_uK.fits'
            print, file
            read_fits_map, file, map

            if dataset eq 'dx9' then file = '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/'+dataset+'_rms_'+sfreq[ifreq]+'.fits'
            if dataset eq 'ffp4' then file = '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/'+dataset+'_drms_'+sfreq[ifreq]+'.fits'
            print, file
            read_fits_map, file, rms 

            bbb = fltarr(npix)

            x = h*1.d9/k/temp

            bb  = 1. / ( exp(x*freq[ifreq])-1. )

            bb0 = 1. / ( exp(x*dust_pivot)-1. )


;##            print, 'Dipole direction: ',reform([dx[ifreq], dy[ifreq], dz[ifreq]]), total([dx[ifreq], dy[ifreq], dz[ifreq]]*[dx[ifreq], dy[ifreq], dz[ifreq]])

            dipole = make_dipole( ns, [dx[ifreq], dy[ifreq], dz[ifreq]])
mollview, dipole, px=400, tit=sfreq[ifreq], win=15+ifreq
mollview, dipole, px=600, tit='!6Residual Dipole @ '+sfreq[ifreq], win=-1, png=dataset+'_'+root+'_residual_dipole_'+sfreq[ifreq]+'.png'
            foreg[*,ifreq] = m[ifreq] + dipole + conversionfactor(freq[ifreq],/antenna2thermo) $
              * ( dust*(freq[ifreq]/dust_pivot)^(emis+1)*bb/bb0 + $
                  co_spec[ifreq]*co_map + $
                  lowf*(freq[ifreq]/lowf_pivot)^beta )
             
            if ((root eq 'pix_01a_v3_ff') or (root eq 'pix_7b_2pl_lowf') or (root eq 'pix_wmap_ff')) $
              then foreg[*,ifreq] = foreg[*,ifreq] + ff*(freq[ifreq]/ff_pivot)^(-2.15)*conversionfactor(freq[ifreq],/antenna2thermo)

            ccc = map-cmb-foreg[*,ifreq]
            res[*,ifreq] = ccc
;##            remove_dipole, ccc, mask, nside=ns, ordering='ring' ;, /onlymonopole
            mollview, ccc, min=-20, max=20, win=ifreq+1, tit='!6Residuals: '+sfreq[ifreq], chars=1.5, units='!7l!6K CMB', px=600
            mollview, ccc, min=-20, max=20, win=-1, tit='!6Residuals: '+sfreq[ifreq], chars=1.5, units='!7l!6K CMB', px=600, png=dataset+'_'+root+'_residuals_'+sfreq[ifreq]+'.png'
            chisq = chisq + ccc^2/rms^2
             
;##        write_fits_map, 'res_'+sfreq[ifreq]+'_hres_'+root+'.fits', ccc, /ring, units='!7l!8K CMB'
            ccc[bp] = -1.6375e30
            thres = stddev(ccc[where(ccc ne -1.6375e30)])

;##        mollview, ccc, min=-thres, max=thres, chars=1.5, units='!7l!17K CMB', tit='Residuals: '+sfreq[ifreq]+' - '+root, png='../ns128/png/res_'+sfreq[ifreq]+'_hres_'+root+'.png', window=-1
;stop
        endfor

;##    write_fits_map, 'chisq_hres_'+root+'.fits', chisq, /ring
;    chisq[bp] = -1.6375e30
        mollview, chisq, chars=1.5, max=21, tit='!7v!u2!n!8', win=nfreq+2, px=600
;    thres =stddev(chisq[where(chisq ne -1.6375e30)])
;;     mollview, chisq, chars=1.5, max=30, tit='!7v!u2!n!8', png='../ns128/png/chisq_hres_'+root+'.png', window=-1
        print, ' Wmap channels next... Type .c'
        stop

        
        for ifreq=0,nwfreq-1 do begin

            print, sfreq[ifreq]

            if dataset eq 'dx9' then file = '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/wmap7_map_'+wsfreq[ifreq]+'.fits'
            if dataset eq 'ffp4' then file = '/global/scratch/sd/dpietrob/'+dataset+'/maps/ns0128/'+dataset+'_Imap_'+wsfreq[ifreq]+'_ns128_uK.fits'
            print, file
            read_fits_map, file, map

            file = '/global/scratch/sd/dpietrob/dx9/maps/ns0128/wmap7_rms_'+wsfreq[ifreq]+'.fits'
            print, file
            read_fits_map, file, rms 

            bbb = fltarr(npix)

            x = h*1.d9/k/temp

            bb  = 1. / ( exp(x*wfreq[ifreq])-1. )

            bb0 = 1. / ( exp(x*dust_pivot)-1. )

            dipole = make_dipole( ns, [dx[7+ifreq], dy[7+ifreq], dz[7+ifreq]])
mollview, dipole, px=400, tit=wsfreq[ifreq], win=15+ifreq
mollview, dipole, tit='!6Residual Dipole @ '+wsfreq[ifreq], win=-1, png=dataset+'_'+root+'_residual_dipole_'+wsfreq[ifreq]+'.png'
            foreg[*,ifreq] = m[7+ifreq] + dipole + conversionfactor(wfreq[ifreq],/antenna2thermo) $
              * ( dust*(wfreq[ifreq]/dust_pivot)^(emis+1)*bb/bb0 + $
;                  co_spec[ifreq]*co_map + $
                  lowf*(wfreq[ifreq]/lowf_pivot)^beta )
             
            if ((root eq 'pix_01a_v3_ff') or (root eq 'pix_7b_2pl_lowf') or (root eq 'pix_wmap_ff')) $
              then foreg[*,ifreq] = foreg[*,ifreq] + ff*(wfreq[ifreq]/ff_pivot)^(-2.15)*conversionfactor(wfreq[ifreq],/antenna2thermo)

            ccc = map-cmb-foreg[*,ifreq]
            res[*,ifreq] = ccc
;##            remove_dipole, ccc, mask, nside=ns, ordering='ring' ;, /onlymonopole
            mollview, ccc, min=-20, max=20, win=ifreq+1, tit='!6Residuals: '+wsfreq[ifreq], chars=1.5, units='!7l!6K CMB', px=600
            mollview, ccc, min=-20, max=20, win=-1, tit='!6Residuals: '+wsfreq[ifreq], chars=1.5, units='!7l!6K CMB', px=600, png=dataset+'_'+root+'_residuals_'+wsfreq[ifreq]+'.png'
;            chisq = chisq + ccc^2/rms^2
        endfor

        stop
    endfor
endif else begin
; ------ chisq from single maps ----------------------------------------
    dust_pivot = 360.6
    sync_pivot = 28.4
    ff_pivot = 70.3

    read_fits_map, '/global/homes/d/dpietrob/myscratch/'+dataset+'/output/glss/constant_indices/outmap/'+dataset+'_constIndx_fs_dust.fits', dust
    read_fits_map, '/global/homes/d/dpietrob/myscratch/'+dataset+'/output/glss/constant_indices/outmap/'+dataset+'_constIndx_fs_co.fits', co_map
    read_fits_map, '/global/homes/d/dpietrob/myscratch/'+dataset+'/output/glss/constant_indices/outmap/'+dataset+'_constIndx_fs_sync.fits', sync
    read_fits_map, '/global/homes/d/dpietrob/myscratch/'+dataset+'/output/glss/constant_indices/outmap/'+dataset+'_constIndx_fs_ff.fits', ff

    read_fits_map, '/global/homes/d/dpietrob/myscratch/'+dataset+'/output/glss/constant_indices/outmap/'+dataset+'_constIndx_fs_cmb.fits', cmb
    
    emis = 1.6
    temp = 18.
    beta_s = -3
    beta_f = -2.15

    h = 6.626068 * 10.^(-34)
    k = 1.3806503 * 10.^(-23)
    c = 2.99792458d8

    foreg = fltarr(Npix, Nfreq)
    chisq = fltarr(Npix)

    for ifreq=0,nfreq-1 do begin

        print, sfreq[ifreq]

        file = '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/ns2048/'+dataset+'_Imap_'+sfreq[ifreq]+'_ns2048_uK.fits'
        print, file
        read_fits_map, file, map

        file = '/global/homes/d/dpietrob/myscratch/'+dataset+'/maps/ns2048/'+dataset+'_Irms_'+sfreq[ifreq]+'_ns2048_uK.fits'
        print, file
        read_fits_map, file, rms 

        bbb = fltarr(npix)

        x = h*1.d9/k/temp

        bb  = 1. / ( exp(x*freq[ifreq])-1. )

        bb0 = 1. / ( exp(x*dust_pivot)-1. )

        foreg[*,ifreq] = conversionfactor(freq[ifreq],/antenna2thermo) $
          * ( dust*(freq[ifreq]/dust_pivot)^(emis+1)*bb/bb0 + $
              co_spec[ifreq]*co_map + $
              sync*(freq[ifreq]/sync_pivot)^beta_s + $
              ff*(freq[ifreq]/ff_pivot)^beta_f )
        
        ccc = map-cmb-foreg[*,ifreq]
        mollview, ccc, min=-50, max=50, win=ifreq+1, tit='!6Residuals: '+sfreq[ifreq], chars=1.5, units='!7l!6K CMB'

        chisq = chisq + ccc^2/rms^2
    endfor

    mollview, chisq, chars=1.5, max=21, tit='!7v!u2!n!8', win=nfreq+2
    
endelse
STOP

end
