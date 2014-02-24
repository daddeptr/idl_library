pro make_glss_paramfile, runtag, localdir, pbsout

   print, 'restoring '+localdir+runtag+'.par.sav...' 
   restore, localdir+runtag+'.par.sav'

   spawn, 'mkdir -p '+localdir+'glss_'+chains_dir
   spawn, 'mkdir -p '+localdir+'runfiles'

   for i=0,2 do begin
       if (not do_halves) then begin
           lfi_names[i,0] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/dx9_Delta_Imap_'+lfi_tags[i]+'_ns2048_uK.fits'
           lfi_names[i,1] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/dx9_Delta_Irms_'+lfi_tags[i]+'_ns2048_uK.fits'
       endif else begin
           lfi_names[i,0] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/dx9_Delta_Imap_'+lfi_tags[i]+'_ns2048_uK_'+half_tag+'.fits'
           lfi_names[i,1] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/dx9_Delta_Irms_'+lfi_tags[i]+'_ns2048_uK_'+half_tag+'.fits'
       endelse
   endfor

; --- HFI --------------------------------------------------------------
   for i=0,5 do begin
       if (not do_halves) then begin
           hfi_names[i,0] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/dx9_Imap_'+hfi_tags[i]+'_ns2048_uK.fits'
           hfi_names[i,1] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/dx9_Irms_'+hfi_tags[i]+'_ns2048_uK.fits'
       endif else begin
           hfi_names[i,0] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/dx9_Imap_'+hfi_tags[i]+'_ns2048_uK_'+half_tag+'.fits'
           hfi_names[i,1] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/dx9_Irms_'+hfi_tags[i]+'_ns2048_uK_'+half_tag+'.fits'
       endelse
   endfor

; --- 100GHz -----------------------------------------------------------
   nfiles = n_elements(d100_tags)
   for i=0,nfiles-1 do begin
       if (not do_halves) then begin
           d100_names[i,0] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/nominal/bolo/dx9_Imap100_ns2048_uK_'+d100_tags[i]+'.fits'
           d100_names[i,1] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/nominal/bolo/dx9_Irms100_ns2048_uK_'+d100_tags[i]+'.fits'
       endif else begin
           d100_names[i,0] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/nominal/bolo/dx9_Imap100_ns2048_uK_'+d100_tags[i]+'_'+half_tag+'.fits'
           d100_names[i,1] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/nominal/bolo/dx9_Irms100_ns2048_uK_'+d100_tags[i]+'_'+half_tag+'.fits'
       endelse
   endfor

; --- 217GHz -----------------------------------------------------------
   nfiles = n_elements(d217_tags)
   for i=0,nfiles-1 do begin
       if (not do_halves) then begin
           d217_names[i,0] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/nominal/bolo/dx9_Imap217_ns2048_uK_'+d217_tags[i]+'.fits'
           d217_names[i,1] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/nominal/bolo/dx9_Irms217_ns2048_uK_'+d217_tags[i]+'.fits'
       endif else begin
           d217_names[i,0] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/nominal/bolo/dx9_Imap217_ns2048_uK_'+d217_tags[i]+'_'+half_tag+'.fits'
           d217_names[i,1] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/nominal/bolo/dx9_Irms217_ns2048_uK_'+d217_tags[i]+'_'+half_tag+'.fits'
       endelse
   endfor

; --- 353GHz -----------------------------------------------------------
   nfiles = n_elements(d353_tags)
   for i=0,nfiles-1 do begin
       if (not do_halves) then begin
           d353_names[i,0] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/nominal/bolo/dx9_Imap353_ns2048_uK_'+d353_tags[i]+'.fits'
           d353_names[i,1] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/nominal/bolo/dx9_Irms353_ns2048_uK_'+d353_tags[i]+'.fits'
       endif else begin
           d353_names[i,0] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/nominal/bolo/dx9_Imap353_ns2048_uK_'+d353_tags[i]+'_'+half_tag+'.fits'
           d353_names[i,1] = '/global/scratch/sd/dpietrob/dx9/maps/ns2048/nominal/bolo/dx9_Irms353_ns2048_uK_'+d353_tags[i]+'_'+half_tag+'.fits'
       endelse
   endfor

; --- WMAP 7-yr --------------------------------------------------------
   for i=0,4 do begin
       if (not do_halves) then begin
           wmap_names[i,0] = '/global/homes/d/dpietrob/myscratch/dx9/wmap_ns1024/wmap7_'+wmap_tags[i]+'_map_ns2048_uK.fits'
           wmap_names[i,1] = '/global/homes/d/dpietrob/myscratch/dx9/wmap_ns1024/wmap7_'+wmap_tags[i]+'_rms_ns2048_uK.fits'
       endif else begin
           print, '!!! WARNING: WMAP only full dataset'
           wmap_names[i,0] = '/global/homes/d/dpietrob/myscratch/dx9/wmap_ns1024/wmap7_'+wmap_tags[i]+'_map_ns2048_uK.fits'
           wmap_names[i,1] = '/global/homes/d/dpietrob/myscratch/dx9/wmap_ns1024/wmap7_'+wmap_tags[i]+'_rms_ns2048_uK.fits'
       endelse
   endfor

; --- Haslam -----------------------------------------------------------

; ----------------------------------------------------------------------

print, ' --> Components: ', ncomp

print, ' Number of channel selected: ', Nchannels

iname=0
if ilfi[0] ge 0 then for i=0,n_elements(ilfi)-1 do begin
    names[iname,*] = lfi_names[ilfi[i],*]
    print, ' ---> '+names[iname,*]
    iname = iname + 1
endfor
if ihfi[0] ge 0 then for i=0,n_elements(ihfi)-1 do begin
    names[iname,*] = hfi_names[ihfi[i],*]
    print, ' ---> '+names[iname,*]
    iname = iname + 1
endfor
if iwmap[0] ge 0 then for i=0,n_elements(iwmap)-1 do begin
    names[iname,*]=wmap_names[iwmap[i],*]
    print, ' ---> '+names[iname,*]
    iname = iname+1
endfor
if ihaslam[0] ge 0 then begin
    names[iname,*] = haslam_names
    print, ' ---> '+names[iname,*]
    iname = iname+1
endif
if i100[0] ge 0 then for i=0,n_elements(i100)-1 do begin
    names[iname,*]=d100_names[i100[i],*]
    print, ' ---> '+names[iname,*]
    iname = iname+1
endfor
if i217[0] ge 0 then for i=0,n_elements(i217)-1 do begin
    names[iname,*]=d217_names[i217[i],*]
    print, ' ---> '+names[iname,*]
    iname=iname+1
endfor
if i353[0] ge 0 then for i=0,n_elements(i353)-1 do begin
    names[iname,*]=d353_names[i353[i],*]
    print, ' ---> '+names[iname,*]
    iname = iname + 1
endfor

for i=0,nfreq-1 do print, freq[i]

names = names[ifreq,*]

if (not no_sample) then begin
    if (ico1 gt 0) then begin
        print, 'reading, 1-0...'
        spawn, 'wc -l '+localdir+chains_dir+'/fg_tab_c'+string(nchain,format='(i4.4)')+'_no'+string(ico1, format='(i2.2)')+'.dat', lines
        help, lines
;##        print, fix(lines[0])

        nsample = fix(lines[0])-1
        xxx = ' '
        co1_fit = fltarr(nfreq+1, nsample)
        openr,1,localdir+chains_dir+'/fg_tab_c'+string(nchain,format='(i4.4)')+'_no'+string(ico1, format='(i2.2)')+'.dat'
        readf,1,xxx
        readf,1,co1_fit
        close,1
    endif
;print, co1_fit[*,nsample-1]
;stop

    if (ico2 gt 0) then begin
        print, 'reading, 2-1...'
        spawn, 'wc -l '+localdir+chains_dir+'/fg_tab_c'+string(nchain,format='(i4.4)')+'_no'+string(ico2, format='(i2.2)')+'.dat', lines
        help, lines
;##        print, fix(lines[0])

        nsample = fix(lines[0])-1
        xxx = ' '
        co2_fit = fltarr(nfreq+1, nsample)
        openr,1,localdir+chains_dir+'/fg_tab_c'+string(nchain,format='(i4.4)')+'_no'+string(ico2, format='(i2.2)')+'.dat'
        readf,1,xxx
        readf,1,co2_fit
        close,1
    endif

;print, co2_fit[*,nsample-1]
;stop
    if (ico3 gt 0) then begin
        print, 'reading, 3-2...'
        spawn, 'wc -l '+localdir+chains_dir+'/fg_tab_c'+string(nchain,format='(i4.4)')+'_no'+string(ico3, format='(i2.2)')+'.dat', lines
        help, lines
;##        print, fix(lines[0])

        nsample = fix(lines[0])-1
        xxx = ' '
        co3_fit = fltarr(nfreq+1, nsample)
        openr,1,localdir+chains_dir+'/fg_tab_c'+string(nchain,format='(i4.4)')+'_no'+string(ico3, format='(i2.2)')+'.dat'
        readf,1,xxx
        readf,1,co3_fit
        close,1
    endif
;print, co3_fit[*,nsample-1]
;stop
    if ( (idust gt 0) and (do_dust_corr) )then begin
        print, 'reading dust...'
        spawn, 'wc -l '+localdir+chains_dir+'/fg_tab_c'+string(nchain,format='(i4.4)')+'_no'+string(idust, format='(i2.2)')+'.dat', lines
        help, lines
;##        print, fix(lines[0])

        nsample = fix(lines[0])-1
        xxx = ' '
        dust_fit = fltarr(nfreq+1, nsample)
        openr,1,localdir+chains_dir+'/fg_tab_c'+string(nchain,format='(i4.4)')+'_no'+string(idust, format='(i2.2)')+'.dat'
        readf,1,xxx
        readf,1,dust_fit
        close,1
    endif
endif else begin
    fsample = 1
    lsample = 1
    nsample = 1

endelse

   for isample = fsample, min([lsample,nsample]) do begin
;##       print, co_sed
       co_sed_fit = co_sed * 0.
       dust_sed_fit = fltarr( nfreq )
       if (not no_sample) then begin 
           if (ico1 gt 0) then co_sed_fit[*,0] = co1_fit[1:*, isample-1]
           if (ico2 gt 0) then co_sed_fit[*,1] = co2_fit[1:*, isample-1]
           if (ico3 gt 0) then co_sed_fit[*,2] = co3_fit[1:*, isample-1]
       endif else begin
           co_sed_fit = co_sed
       endelse
       if ( (idust gt 0) and (do_dust_corr) )then dust_sed_fit = dust_fit[1:*, isample-1]
       print, ' ----- ', isample
       print, co_sed_fit
       print, ''
       print, dust_sed_fit

       foreg = dblarr(4,nfreq)
       if (not zero_monodipole) then begin
           openr,1,localdir+chains_dir+'/foreground_c'+string(nchain,format='(i4.4)')+'_k'+string(isample, format='(i5.5)')+'.dat'
           readf,1,foreg
           close,1
       endif

       lineco1file = localdir+'glss_'+chains_dir+'/lines_'+runtag+'_CO_100_s'+string(isample,format='(i5.5)')+'.txt'
       lineco2file = localdir+'glss_'+chains_dir+'/lines_'+runtag+'_CO_217_s'+string(isample,format='(i5.5)')+'.txt'
       lineco3file = localdir+'glss_'+chains_dir+'/lines_'+runtag+'_CO_353_s'+string(isample,format='(i5.5)')+'.txt'
       linedustfile = localdir+'glss_'+chains_dir+'/lines_'+runtag+'_dust_s'+string(isample,format='(i5.5)')+'.txt'
       if (not do_halves) then parfile = localdir+'glss_'+chains_dir+'/dp_glss_paramfile_'+runtag+'_s'+string(isample,format='(i5.5)')+'.txt' else parfile = localdir+'glss_'+chains_dir+'/dp_glss_paramfile_'+runtag+'_s'+string(isample,format='(i5.5)')+'_'+half_tag+'.txt'

       openw, 1, parfile

       if (ico1 ne -1) then openw, 2, lineco1file
       if (ico2 ne -1) then openw, 3, lineco2file
       if (ico3 ne -1) then openw, 4, lineco3file
       if ( (idust ne -1) and (do_dust_corr) )then openw,5, linedustfile

       printf,1,"# ------ '+runtag+ '---------------------------------------------"
       if (not do_halves) then printf,1,"output_root = '!"+localdir+"glss_"+chains_dir+"/dx9_CO_ns2048_s"+string(isample,format='(i5.5)')+"_'" else printf,1,"output_root = '!"+localdir+"glss_"+chains_dir+"/dx9_CO_ns2048_s"+string(isample,format='(i5.5)')+"_"+half_tag+"_'"
       printf,1,"nside = 2048"
       printf,1,"t_cmb = 2.725d0" 
       printf,1,"#polarization not supported"
       printf,1,"has_pol        = .false."
       printf,1,"output_weights = .false."
       printf,1,"output_noise   = .false."
       printf,1,"do_sampling    = .false."
       printf,1," "
       printf,1,"seed_01 = 5690" 
       printf,1,"seed_02 = 11195" 
       printf,1,"seed_03 = 6007"
       printf,1,"seed_04 = 5138"
       printf,1," "
       printf,1,"number_of_frequencies = "+string(nfreq,format='(1i2.2)')    
       printf,1,"#must be the same for all the channels. only rms supported for now"            
       printf,1,"noise_format = rms"
       printf,1," "

       if (ico1 ne -1) then printf,2,nfreq
       if (ico2 ne -1) then printf,3,nfreq 
       if (ico3 ne -1) then printf,4,nfreq
       if ( (idust ne -1) and (do_dust_corr) )then printf,5,nfreq

       for i=0,nfreq-1 do begin
           printf,1,"freq_"+string(i+1,format='(1i2.2)')+"  = "+string(freq[i], format='(1f10.3)')
           printf,1,"name_"+string(i+1,format='(1i2.2)')+"  = "+det_tags[i]
           printf,1,"data_"+string(i+1,format='(1i2.2)')+"  = '"+names[i,0]+"'"
           printf,1,"noise_"+string(i+1,format='(1i2.2)')+" = '"+names[i,1]+"'"
           printf,1,"monopole_"+string(i+1,format='(1i2.2)')+" = "+string(foreg[0,i], format='(1f15.5)')
           printf,1,"dipole_x_"+string(i+1,format='(1i2.2)')+" = "+string(foreg[1,i], format='(1f15.5)')
           printf,1,"dipole_y_"+string(i+1,format='(1i2.2)')+" = "+string(foreg[2,i], format='(1f15.5)')
           printf,1,"dipole_z_"+string(i+1,format='(1i2.2)')+" = "+string(foreg[3,i], format='(1f15.5)')
           printf,1," "
; ---
;##       printf,5,string(freq[i],format='(1f10.3)'), ' 1.0 F'
           if (ico1 ne -1) then printf,2,det_tags[i], string(co_sed_fit[i,0],format='(1f15.5)') ;, ' '+sampleCO100[i]+' '+det_tags[i]
           if (ico2 ne -1) then printf,3,det_tags[i], string(co_sed_fit[i,1],format='(1f15.5)') ;, ' '+sampleCO217[i]+' '+det_tags[i]
           if (ico3 ne -1) then printf,4,det_tags[i], string(co_sed_fit[i,2],format='(1f15.5)') ;, ' '+sampleCO353[i]+' '+det_tags[i]
           if ( (idust ne -1) and (do_dust_corr) )then printf,5,det_tags[i], string(dust_sed_fit[i],format='(1f15.5)') ;, ' '+sampleCO353[i]+' '+det_tags[i]
; ---
       endfor
       printf,1,"# ------ Foreground ----------------------------------------------------"
       printf,1,"number_of_components = "+strtrim(string(ncomp),2)
       printf,1," "
       if (icmb gt 0) then begin
           printf,1,"fg_type_"+string(icmb, format='(1i2.2)')+" = cmb"
           printf,1,"fg_name_"+string(icmb, format='(1i2.2)')+" = cmb"
           printf,1," "
       endif
       if ( (idust gt 0) and (not do_dust_corr) )then begin
           printf,1,"fg_type_"+string(idust, format='(1i2.2)')+" = gray_body"
           printf,1,"fg_name_"+string(idust, format='(1i2.2)')+" = dust"
           printf,1,"nu_pivot_"+string(idust, format='(1i2.2)')+" = "+string(dust_ref,format='(f10.3)')
           printf,1,"emissivity_map_"+string(idust, format='(1i2.2)')+"  = '"+localdir+chains_dir+"/fg_ind_map_no01_c"+string(nchain, format='(i4.4)')+"_k"+string(isample, format='(i5.5)')+".fits'"
           printf,1,"temperature_map_"+string(idust, format='(1i2.2)')+" = '"+localdir+chains_dir+"/fg_ind_map_no02_c"+string(nchain, format='(i4.4)')+"_k"+string(isample, format='(i5.5)')+".fits'"
           printf,1," "
       endif
       if ( (idust gt 0) and (do_dust_corr) )then begin
           printf,1,"fg_type_"+string(idust, format='(1i2.2)')+" = modulated_gray_body"
           printf,1,"fg_name_"+string(idust, format='(1i2.2)')+" = dust"
           printf,1,"nu_pivot_"+string(idust, format='(1i2.2)')+" = "+string(dust_ref,format='(f10.3)')
           printf,1,"emissivity_map_"+string(idust, format='(1i2.2)')+"  = '"+localdir+chains_dir+"/fg_ind_map_no01_c"+string(nchain, format='(i4.4)')+"_k"+string(isample, format='(i5.5)')+".fits'"
           printf,1,"temperature_map_"+string(idust, format='(1i2.2)')+" = '"+localdir+chains_dir+"/fg_ind_map_no02_c"+string(nchain, format='(i4.4)')+"_k"+string(isample, format='(i5.5)')+".fits'"
           printf,1,"lines_file_"+string(idust, format='(1i2.2)')+" = '"+linedustfile+"'"
           printf,1," "
       endif
       if (isync gt 0) then begin
           printf,1,"fg_type_"+string(isync, format='(1i2.2)')+" = power_law"
           printf,1,"fg_name_"+string(isync, format='(1i2.2)')+" = sync"
           printf,1,"nu_pivot_"+string(isync, format='(1i2.2)')+" = "+string(sync_ref,format='(f10.3)')
           printf,1,"index_map_"+string(isync, format='(1i2.2)')+" = '"+localdir+chains_dir+"/fg_ind_map_no03_c"+string(nchain, format='(i4.4)')+"_k"+string(isample, format='(i5.5)')+".fits'"
           printf,1," "
       endif
       if (ico1 gt 0) then begin
           printf,1,"fg_type_"+string(ico1, format='(1i2.2)')+"   = line"
           printf,1,"fg_name_"+string(ico1, format='(1i2.2)')+"   = co10"
           printf,1,"lines_file_"+string(ico1, format='(1i2.2)')+" = '"+lineco1file+"'"
;##       printf,1,"lines_file_"+string(ico1, format='(1i2.2)')+" = '"+localdir+"lines_dx9_CO_100_spencer.txt'"
           printf,1," "
       endif
       if (ico2 gt 0) then begin
           printf,1,"fg_type_"+string(ico2, format='(1i2.2)')+"   = line"
           printf,1,"fg_name_"+string(ico2, format='(1i2.2)')+"   = co21"
           printf,1,"lines_file_"+string(ico2, format='(1i2.2)')+" = '"+lineco2file+"'"
;##       printf,1,"lines_file_"+string(ico2, format='(1i2.2)')+" = '"+localdir+"lines_dx9_CO_217_spencer.txt'"
           printf,1," "
       endif
       if (ico3 gt 0) then begin
           printf,1,"fg_type_"+string(ico3, format='(1i2.2)')+"   = line"
           printf,1,"fg_name_"+string(ico3, format='(1i2.2)')+"   = co32"
           printf,1,"lines_file_"+string(ico3, format='(1i2.2)')+" = '"+lineco3file+"'"
;##       printf,1,"lines_file_"+string(ico3, format='(1i2.2)')+" = '"+localdir+"lines_dx9_CO_353_spencer.txt'"
       endif

       close, /all


;   print, ' runtag = ', runtag
;   print, ' localdir = ', localdir

;##       openw, 1, localdir+runtag+'_glss_s'+string(isample, format='(i5.5)')+'_run.pbs'
       if (not do_halves) then pbsout = localdir+'runfiles/'+runtag+'_glss_run_s'+string(isample, format='(i5.5)')+'.pbs' else pbsout = localdir+'runfiles/'+runtag+'_glss_run_s'+string(isample, format='(i5.5)')+'_'+half_tag+'.pbs'
       openw, 1, pbsout
       printf, 1, '#PBS -S /bin/tcsh'
       printf, 1, '#PBS -q serial'
       printf, 1, '### #PBS -q usplanck'
       printf, 1, '#PBS -l walltime=00:20:00'
       printf, 1, '#PBS -l pvmem=16GB'
       printf, 1, '### #PBS -A usplanck'
       printf, 1, '#PBS -j eo'
       printf, 1, '#PBS -e eg.$PBS_JOBID.err'
       printf, 1, '#PBS -N glss.0001'
       printf, 1, '#PBS -m e'
       printf, 1, ''
       printf, 1, 'cd '+localdir
       printf, 1, ''
       printf, 1, 'module swap pgi intel ; module swap openmpi openmpi-intel ; module load mkl ; echo Intel loaded'
       printf, 1, 'setenv HEALPIX /project/projectdirs/cmb/modules/carver/gnu/cmb/2.5.1/healpix_2.15a-2.5.1/'
       printf, 1, 'setenv OMP_NUM_THREADS 1'
       printf, 1, ''
       printf, 1, 'time /global/scratch/sd/dpietrob/Tools/glss_dust_corr.x.carver '+parfile
       close, 1

       if (run_pbs) then begin
           print, ' submitting '+pbsout
           spawn, 'qsub '+pbsout
       endif
   endfor

;##   stop, ' make_glss_paramfile: END'

end
