   True = 1b
   False = 0b

   ;# ---------------------------------------------------------------------
   if (False) then begin
       indir = '/global/scratch/sd/loris/CompSep/GLSS/FFP/4/7b_rerun/avrg_maps/'
       dataset = 'ffp4'
       version = 'v3.2.128'
       outdir = '/global/scratch/sd/dpietrob/commander-ruler/'+dataset+'/'+version+'/'

       components = ['cmb', 'fg_dust_amp', 'fg_co_amp', 'fg_lowfreq_amp']
       incomp = ['cmb', 'dust', 'co', 'ff' ]
       ncomp = 4
       set = ['', '_hr1', '_hr2'] ;['fs','1h','2h']
       outset = ['', '_hr1', '_hr2']
       for icomp=0,0 do begin
           print, components[icomp]
           for iset=0,2 do begin
               infile = indir+'ffp4_7b_mondip'+set[iset]+'_avrg_'+incomp[icomp]+'.fits'
               outfile = outdir+'ffp4_commander-ruler_'+version+'_'+components[icomp]+outset[iset]+'_uK.fits'
               spawn, 'ls '+infile
               spawn, 'ls '+outfile
;stop
               spawn, 'cp '+infile+' '+outfile

           endfor
       endfor
       stop
   endif


;# ---------------------------------------------------------------------
   if (True) then begin
       indir = '/global/scratch/sd/loris/CompSep/GLSS/FFP/4/7b_rerun/components/'
       dataset = 'ffp4'
       version = 'v3.2.128'
       outdir = '/global/scratch/sd/dpietrob/commander-ruler/'+dataset+'/'+version+'/'

       spawn, 'mkdir -p '+outdir+'components/'
       outdir = outdir+'components/'

       components = ['cmb', 'fg_dust_amp', 'fg_co_amp', 'fg_lowfreq_amp']
       comp_maps = ['cmb_scalar', 'compact2', 'diffuse', 'noise']
       incomp = ['cmb', 'dust', 'CO', 'lowfreq' ]
       ncomp = 4
       set = ['', '_hr1', '_hr2'] ;['fs','1h','2h']
       outset = ['', '_hr1', '_hr2']
       for icomp=0,0 do begin
           print, components[icomp]
           for iset=0,2 do begin
               for icompmap=0,3 do begin
                   infile = indir+incomp[icomp]+'/ffp4_7b_'+comp_maps[icompmap]+'_to_'+incomp[icomp]+set[iset]+'_map.fits'
                   outfile = outdir+'ffp4_commander-ruler_'+version+'_'+comp_maps[icompmap]+'_to_'+components[icomp]+outset[iset]+'_uK.fits'
                   print, ' -->'
                   spawn, 'ls '+infile
                   spawn, 'ls '+outfile
;stop
                   spawn, 'cp '+infile+' '+outfile
               endfor
           endfor
       endfor
       stop
   endif


;# ---------------------------------------------------------------------
   if (False) then begin
       components = ['cmb', 'fg_co_amp', 'fg_dust_amp', 'fg_lowfreq_amp']
       ncomp = 1
       dataset = 'ffp4'
       version = 'v3.2.128'
       dir = '/global/homes/d/dpietrob/myscratch/commander-ruler/'+dataset+'/'+version+'/'
       for icomp=0,ncomp-1 do begin
           read_fits_map, dir+dataset+'_commander-ruler_'+version+'_'+components[icomp]+'_hr1_uK.fits', h1 
           read_fits_map, dir+dataset+'_commander-ruler_'+version+'_'+components[icomp]+'_hr2_uK.fits', h2
           d = (h1-h2)/2.
           write_fits_map, dir+dataset+'_commander-ruler_'+version+'_'+components[icomp]+'_error_uK.fits', d, /ring, units='!7l!8K'
           mollview, dir+dataset+'_commander-ruler_'+version+'_'+components[icomp]+'_error_uK.fits', /hist
       endfor
   endif

;# ---------------------------------------------------------------------
   if (False) then begin
       print, ' Mean low-frequency component maps'
       ff_pivot = 30.
       nchain = 4
       nsample = 45
       dataset = 'ffp4'
       print, dataset
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
           for isamp=1,nsample do begin
               tot_samp = tot_samp + 1
               schain = string(ich,format='(i3.3)')
               ssamp = string(isamp, format='(i5.5)')
               if (dataset eq 'dx8') then file = dir + 'c'+schain + '/'+dataset+'_7b_mondip_fs_c0'+schain+'_k'+ssamp+'_ff.fits'
               if (dataset eq 'ffp4') then file = dir + 'c'+schain + '/'+dataset+'_7b_mondip_c0'+schain+'_k'+ssamp+'_ff.fits'
               read_fits_map, file, amp
               file = '/global/scratch/sd/dpietrob/'+dataset+'/input/high_res_indices/pix_raw_effreq_scan_v5/lfc_indx_map_c0'+schain+'_k'+ssamp+'_ns2048.fits'
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
       nsample = 45
       dataset = 'ffp4'
       print, dataset
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
           for isamp=1,nsample do begin
               tot_samp = tot_samp + 1
               schain = string(ich,format='(i3.3)')
               ssamp = string(isamp, format='(i5.5)')
               if (dataset eq 'dx8') then file = dir + 'c'+schain + '/'+dataset+'_7b_mondip_fs_c0'+schain+'_k'+ssamp+'_dust.fits'
               if (dataset eq 'ffp4') then file = dir + 'c'+schain + '/'+dataset+'_7b_mondip_c0'+schain+'_k'+ssamp+'_dust.fits'
               read_fits_map, file, amp
               file = '/global/scratch/sd/dpietrob/'+dataset+'/input/high_res_indices/pix_raw_effreq_scan_v5/dust_emis_map_c0'+schain+'_k'+ssamp+'_ns2048.fits'
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
