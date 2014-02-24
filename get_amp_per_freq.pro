   True = 1b
   False = 0b


   if (True) then begin
       print, ' Copying files...'
       dataset = 'ffp4'

       if (dataset eq 'dx8') then begin
           dir = '/global/scratch/sd/loris/CompSep/GLSS/DX/8/7b/avrg_maps/'
       endif
       if (dataset eq 'ffp4') then begin
           dir = '/global/scratch/sd/loris/CompSep/GLSS/FFP/4/7b/avrg_maps/'
       endif           

       comps = ['cmb','co','dust','lowfreq']
       ncomps = n_elements(comps)

;# ---------------------------------------------------------------------
       if (False) then begin
           spawn, 'ls '+dir+'*hr1*', files
           print, files
           for ifile=0,n_elements(files)-1 do begin
               if (ifile eq 0) then spawn, 'cp '+files[ifile]+' /global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_'+comps[ifile]+'_hr1_uK.fits'
               if (ifile ne 0) then spawn, 'cp '+files[ifile]+' /global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_'+comps[ifile]+'_amp_hr1_uK.fits'
           endfor
       endif

;# ---------------------------------------------------------------------
       if (False) then begin
           spawn, 'ls '+dir+'*hr2*', files
           print, files
           for ifile=0,n_elements(files)-1 do begin
               if (ifile eq 0) then spawn, 'cp '+files[ifile]+' /global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_'+comps[ifile]+'_hr2_uK.fits'
               if (ifile ne 0) then spawn, 'cp '+files[ifile]+' /global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_'+comps[ifile]+'_amp_hr2_uK.fits'
           endfor
stop
       endif

;# ---------------------------------------------------------------------
       if (True) then begin
           for ifile=0, ncomps-1 do begin
               if (ifile eq 0) then begin
                   read_fits_map, '/global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_'+comps[ifile]+'_hr1_uK.fits', hr1
                   read_fits_map, '/global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_'+comps[ifile]+'_hr2_uK.fits', hr2
               endif
               if (ifile ne 0) then begin
                   read_fits_map, '/global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_'+comps[ifile]+'_amp_hr1_uK.fits', hr1
                   read_fits_map, '/global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_'+comps[ifile]+'_amp_hr2_uK.fits', hr2
               endif

               error = (hr1-hr2)/2.
               if (ifile eq 0) then write_fits_map, '/global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_'+comps[ifile]+'_error_uK.fits', error, /ring, units='!7l!8K' 
               if (ifile ne 0) then write_fits_map, '/global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_'+comps[ifile]+'_amp_error_uK.fits', error, /ring, units='!7l!8K' 
               mollview, error
           endfor
       endif

;# ---------------------------------------------------------------------
       if (False) then spawn, 'cp /global/scratch/sd/loris/CompSep/GLSS/FFP/4/7b/tf/'+dataset+'_7b_tf.fits /global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_cmb_beam.fits'

;# ---------------------------------------------------------------------
       if (False) then begin
           spawn, 'ls '+dir+dataset+'*.fits', files

           print, files
           for ifile = 0, n_elements(files)-1 do begin
               if (ifile eq 0) then spawn, 'cp '+files[ifile]+ ' /global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_'+comps[ifile]+'_uK.fits' 
               if (ifile ne 0) then spawn, 'cp '+files[ifile]+ ' /global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_'+comps[ifile]+'_amp_uK.fits' 
           endfor
       endif

;# ---------------------------------------------------------------------
       if (False) then begin
           if (dataset eq 'dx8') then begin
               dir = '/global/scratch/sd/loris/CompSep/GLSS/DX/8/7b/'
               root = 'fs_samples'
           endif
           if (dataset eq 'ffp4') then begin
               dir = '/global/scratch/sd/loris/CompSep/GLSS/FFP/4/7b/'
               root = 'samples'
           endif
           spawn, 'cp '+dir+'/'+dataset+'_7b_mondip_'+root+'_cmb_cmb.fits /global/scratch/sd/dpietrob/commander-ruler/' + dataset+'/'+dataset+'_commander-ruler_v3_sample_covariance_cmb_cmb.fits'
           spawn, 'cp '+dir+'/'+dataset+'_7b_mondip_'+root+'_co_cmb.fits /global/scratch/sd/dpietrob/commander-ruler/' + dataset+'/'+dataset+'_commander-ruler_v3_sample_covariance_co_cmb.fits'
           spawn, 'cp '+dir+'/'+dataset+'_7b_mondip_'+root+'_dust_cmb.fits /global/scratch/sd/dpietrob/commander-ruler/' + dataset+'/'+dataset+'_commander-ruler_v3_sample_covariance_dust_cmb.fits'
           spawn, 'cp '+dir+'/'+dataset+'_7b_mondip_'+root+'_ff_cmb.fits /global/scratch/sd/dpietrob/commander-ruler/' + dataset+'/'+dataset+'_commander-ruler_v3_sample_covariance_lowfreq_cmb.fits'
           spawn, 'cp '+dir+'/'+dataset+'_7b_mondip_'+root+'_co_co.fits /global/scratch/sd/dpietrob/commander-ruler/' + dataset+'/'+dataset+'_commander-ruler_v3_sample_covariance_co_co.fits'
           spawn, 'cp '+dir+'/'+dataset+'_7b_mondip_'+root+'_dust_co.fits /global/scratch/sd/dpietrob/commander-ruler/' + dataset+'/'+dataset+'_commander-ruler_v3_sample_covariance_dust_co.fits'
           spawn, 'cp '+dir+'/'+dataset+'_7b_mondip_'+root+'_dust_dust.fits /global/scratch/sd/dpietrob/commander-ruler/' + dataset+'/'+dataset+'_commander-ruler_v3_sample_covariance_dust_dust.fits'
           spawn, 'cp '+dir+'/'+dataset+'_7b_mondip_'+root+'_dust_ff.fits /global/scratch/sd/dpietrob/commander-ruler/' + dataset+'/'+dataset+'_commander-ruler_v3_sample_covariance_dust_lowfreq.fits'
           spawn, 'cp '+dir+'/'+dataset+'_7b_mondip_'+root+'_ff_ff.fits /global/scratch/sd/dpietrob/commander-ruler/' + dataset+'/'+dataset+'_commander-ruler_v3_sample_covariance_lowfreq_lowfreq.fits'
           spawn, 'cp '+dir+'/'+dataset+'_7b_mondip_'+root+'_ff_co.fits /global/scratch/sd/dpietrob/commander-ruler/' + dataset+'/'+dataset+'_commander-ruler_v3_sample_covariance_lowfreq_co.fits'

           for icomp=0,3 do begin
               read_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/' + dataset+'/'+dataset+'_commander-ruler_v3_sample_covariance_'+comps[icomp]+'_'+comps[icomp]+'.fits', map
               if (icomp eq 0) then write_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/' + dataset+'/'+dataset+'_commander-ruler_v3_'+comps[icomp]+'_rms_uK.fits', sqrt(map), /ring, units='!7l!8K'
               if (icomp ne 0) then write_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/' + dataset+'/'+dataset+'_commander-ruler_v3_'+comps[icomp]+'_amp_rms_uK.fits', sqrt(map), /ring, units='!7l!8K'
           endfor
stop
       endif

stop
   endif

;# ---------------------------------------------------------------------
   if (False) then begin
       dataset = 'dx8'
       spawn, 'ls commander-ruler/' + dataset+'/', files
       components = ['dust', 'co', 'lowfreq']
       ncomp = 3
       for icomp=0, ncomp-1 do begin
           print, components[icomp]
           for ifile=0,n_elements(files)-1 do begin
               if ( strmatch(files[ifile], components[icomp]+'_amp') ) then begin
                   cut = strsplit( files[ifile], components[icomp] )
                   print, cut
                   stop
               endif
           endfor
       endfor
       stop
   endif


;# ---------------------------------------------------------------------
   if (False) then begin
       components = ['cmb', 'co_amp', 'dust_amp', 'lowfreq_amp']
       ncomp = 4
       dataset = 'ffp4'
       dir = '/global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'
       for icomp=0,ncomp-1 do begin
           read_fits_map, dir+dataset+'_commander-ruler_v3_'+components[icomp]+'_hr1_uK.fits', h1 
           read_fits_map, dir+dataset+'_commander-ruler_v3_'+components[icomp]+'_hr2_uK.fits', h2
           d = (h1-h2)/2.
           write_fits_map, dir+dataset+'_commander-ruler_v3_'+components[icomp]+'_error_uK.fits', d, /ring, units='!7l!8K'
           mollview, dir+dataset+'_commander-ruler_v3_'+components[icomp]+'_error_uK.fits', /hist
       endfor
   endif

;# ---------------------------------------------------------------------
   if (False) then begin
       print, ' Mean low-frequency component maps'
       ff_pivot = 30.
       nchain = 4
;##       nsample = 45
       dataset = 'ffp4'
       print, dataset
       root = 'pix_raw_effreq_scan_v5m'
       print, root
       if (dataset eq 'dx8') then dir = '/global/scratch/sd/loris/CompSep/GLSS/DX/8/7b/'
       if (dataset eq 'ffp4') then dir = '/global/scratch/sd/loris/CompSep/GLSS/FFP/4/7b/'
       freq = [28.4, 44.1, 70.3, 101.2, 142.6, 221.9, 360.6]
       sfreq = ['030','044','070','100','143','217','353']
       nfreq = n_elements(freq)
       Nside = 2048l
       Npix = 12l*Nside^2
       amp_maps = dblarr(Npix,nfreq)
       amp_maps2 = dblarr(Npix,nfreq)
       tot_samp = 0l
       for ich=1,nchain do begin
           print, ich
           schain = string(ich,format='(i3.3)')
           if (dataset eq 'dx8') then spawn, 'ls '+ dir + 'c'+schain+'/'+dataset+'_7b_mondip_fs_c0'+schain+'*_ff.fits', files
           if (dataset eq 'ffp4') then spawn, 'ls '+ dir + 'c'+schain+'/'+dataset+'_7b_mondip_c0'+schain+'*_ff.fits', files
           nsample = n_elements(files)
           print, ' Nsample = ', nsample
;stop
           for isamp=1,nsample do begin
               tot_samp = tot_samp + 1
               ssamp = string(isamp, format='(i5.5)')
               if (dataset eq 'dx8') then file = dir + 'c'+schain + '/'+dataset+'_7b_mondip_fs_c0'+schain+'_k'+ssamp+'_ff.fits'
               if (dataset eq 'ffp4') then file = dir + 'c'+schain + '/'+dataset+'_7b_mondip_c0'+schain+'_k'+ssamp+'_ff.fits'
;               file = files[isamp]
               read_fits_map, file, amp
               file = '/global/scratch/sd/dpietrob/'+dataset+'/input/high_res_indices/'+root+'/lfc_indx_map_c0'+schain+'_k'+ssamp+'_ns2048.fits'
               read_fits_map, file, indx
               for ifreq=0, nfreq-1 do begin
                   map = amp * (freq[ifreq]/ff_pivot)^indx
                   amp_maps[*,ifreq] = amp_maps[*,ifreq] + map
                   amp_maps2[*,ifreq] = amp_maps2[*,ifreq] + map^2
;##                   mollview, map, /hist, tit=sfreq[ifreq], win=ifreq
               endfor
               if ( (tot_samp/20)*20 eq tot_samp ) then begin
                   print, tot_samp, '/180 saving...'
                   for ifreq=0,nfreq-1 do begin
                       print, sfreq[ifreq]
                       write_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_lowfreq_amp_'+sfreq[ifreq]+'_uK.fits', amp_maps[*,ifreq]/tot_samp, /ring, units='!7l!8K Antenna'
                       write_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_lowfreq_rms_'+sfreq[ifreq]+'_uK.fits', sqrt(amp_maps2[*,ifreq]/tot_samp-(amp_maps[*,ifreq]/tot_samp)^2), /ring, units='!7l!8K Antenna'
                   endfor
;stop
               endif
           endfor
       endfor
       print, tot_samp, ' Final saving...'
       for ifreq=0,nfreq-1 do begin
           print, sfreq[ifreq]
           write_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_lowfreq_amp_'+sfreq[ifreq]+'_uK.fits', amp_maps[*,ifreq]/tot_samp, /ring, units='!7l!8K Antenna'
           write_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_lowfreq_rms_'+sfreq[ifreq]+'_uK.fits', sqrt(amp_maps2[*,ifreq]/tot_samp-(amp_maps[*,ifreq]/tot_samp)^2), /ring, units='!7l!8K Antenna'
       endfor
       
   endif

;# ---------------------------------------------------------------------
   if (False) then begin
       print, ' Mean dust map.'
       dust_pivot = 353.
       nchain = 4
;       nsample = 45
       dataset = 'ffp4'
       print, dataset
       root = 'pix_raw_effreq_scan_v5m'
       print, root
       if (dataset eq 'dx8') then dir = '/global/scratch/sd/loris/CompSep/GLSS/DX/8/7b/'
       if (dataset eq 'ffp4') then dir = '/global/scratch/sd/loris/CompSep/GLSS/FFP/4/7b/'
       freq = [28.4, 44.1, 70.3, 101.2, 142.6, 221.9, 360.6]
       sfreq = ['030','044','070','100','143','217','353']
       nfreq = n_elements(freq)
       Nside = 2048l
       Npix = 12l*Nside^2
       amp_maps = dblarr(Npix,nfreq)
       amp_maps2 = dblarr(Npix,nfreq)
       tot_samp = 0l
       read_fits_map, '/global/scratch/sd/dpietrob/'+dataset+'/input/high_res_indices/dust_temp_map_ns2048.fits', dust_temp

       h = 6.626068 * 10.^(-34)
       k = 1.3806503 * 10.^(-23)
       c = 299792458.

       for ich=1,nchain do begin
           print, ich
           schain = string(ich,format='(i3.3)')
           spawn, 'ls '+ dir + 'c'+schain+'/ffp4_7b_mondip_c0'+schain+'*_dust.fits', files
           nsample = n_elements(files)
           print, ' Nsample = ', nsample
           for isamp=1,nsample do begin
               tot_samp = tot_samp + 1
               schain = string(ich,format='(i3.3)')
               ssamp = string(isamp, format='(i5.5)')
               if (dataset eq 'dx8') then file = dir + 'c'+schain + '/'+dataset+'_7b_mondip_fs_c0'+schain+'_k'+ssamp+'_dust.fits'
               if (dataset eq 'ffp4') then file = dir + 'c'+schain + '/'+dataset+'_7b_mondip_c0'+schain+'_k'+ssamp+'_dust.fits'
               read_fits_map, file, amp
               file = '/global/scratch/sd/dpietrob/'+dataset+'/input/high_res_indices/'+root+'/dust_emis_map_c0'+schain+'_k'+ssamp+'_ns2048.fits'
               read_fits_map, file, indx
               for ifreq=0, nfreq-1 do begin

                   x = h*1.d9/k/dust_temp

                   bb  = 1. / ( exp(x*freq[ifreq])-1. )
                   bb0 = 1. / ( exp(x*dust_pivot)-1. )

                   map = amp * (freq[ifreq]/dust_pivot)^(1+indx)*bb/bb0

                   amp_maps[*,ifreq] = amp_maps[*,ifreq] + map
                   amp_maps2[*,ifreq] = amp_maps2[*,ifreq] + map^2
;##                   mollview, map, /hist, tit=sfreq[ifreq], win=ifreq
               endfor
               if ( (tot_samp/20)*20 eq tot_samp ) then begin
                   print, tot_samp, '/ 180 --> saving...'
                   for ifreq=0,nfreq-1 do begin
                       print, sfreq[ifreq]
                       write_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_dust_amp_'+sfreq[ifreq]+'_uK.fits', amp_maps[*,ifreq]/tot_samp, /ring, units='!7l!8K Antenna'
                       write_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_dust_rms_'+sfreq[ifreq]+'_uK.fits', sqrt(amp_maps2[*,ifreq]/tot_samp-(amp_maps[*,ifreq]/tot_samp)^2), /ring, units='!7l!8K Antenna'
                   endfor
;stop
               endif
;               stop
           endfor
       endfor
       print, tot_samp, ' Final saving...'
       for ifreq=0,nfreq-1 do begin
           print, sfreq[ifreq]
           write_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_dust_amp_'+sfreq[ifreq]+'_uK.fits', amp_maps[*,ifreq]/tot_samp, /ring, units='!7l!8K Antenna'
           write_fits_map, '/global/scratch/sd/dpietrob/commander-ruler/'+dataset+'/'+dataset+'_commander-ruler_v3_dust_rms_'+sfreq[ifreq]+'_uK.fits', sqrt(amp_maps2[*,ifreq]/tot_samp-(amp_maps[*,ifreq]/tot_samp)^2), /ring, units='!7l!8K Antenna'
       endfor
       
   endif

STOP

END
