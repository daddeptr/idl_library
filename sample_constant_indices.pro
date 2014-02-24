   True = 1b
   False = 0b

   root  = 'pix_01a_v2'
   sfreq = ['030','044','070','100','143','217','353']
   freq  = [28.44, 44.12, 70.33, 101.28, 142.85, 222.51, 361.46]
   nfreq = n_elements(freq)
   mfile = 'ffp4_01a_'+sfreq+'_ns128_uK.fits'
   rfile = 'ffp4_01a_srms_'+sfreq+'_v3.fits'

   Nsamples = 2l

   Nside = 128l
   Npix = 12l*Nside^2
   Ncomp = 4

   icmb = 3

   nofile = ''

; ----------------------------------------------------------------------
   tmpmap = {nametag:nofile, $
             cfreq: 0.d0, $
             map_filename: nofile, $
             rms_filename: nofile, $
             band_pass: nofile }
   
   channels = replicate(tmpmap, nfreq) 

;##   readcol, 'chains/'+root+'/foreground_c0001_k00001.dat', m, dx, dy, dz
;##   for i=0,6 do print,m[i],dx[i],dy[i],dz[i]

   for ifreq=0,nfreq-1 do begin
       channels[ifreq].cfreq = freq[ifreq]
       channels[ifreq].nametag = sfreq[ifreq]
       channels[ifreq].map_filename = mfile[ifreq]
       channels[ifreq].rms_filename = rfile[ifreq]
;##       channels[ifreq].monopole = m[ifreq]
;##       channels[ifreq].dipole = [dx[ifreq], dy[ifreq], dz[ifreq]]
   endfor

; ----------------------------------------------------------------------
   index_pars = { mean:0.d0, sigma:0.d0, lower:0.d0, upper:0.d0}

   tmpcomp = { nametag:nofile, $
               type: nofile, $
               nu_ref:0.d0, $
               indx1_pars: index_pars, $
               indx1_filename: nofile, $
               indx1_current: 0.d0, $
               indx1_constant: True, $
               indx2_pars: index_pars, $
               indx2_constant: True, $
               indx2_filename: nofile, $
               indx2_current: 0.d0, $
               spectrum_filename: nofile }
   
   foregrounds = replicate(tmpcomp,ncomp)

   foregrounds[0].type = 'thermal_dust'
   foregrounds[0].nametag = 'thermal_dust'
   foregrounds[0].nu_ref = 361.46d0
   foregrounds[0].indx1_pars.mean = 1.5
   foregrounds[0].indx1_pars.sigma = 0.3
   foregrounds[0].indx1_pars.lower = 0.5
   foregrounds[0].indx1_pars.upper = 4.0
   foregrounds[0].indx1_current = foregrounds[0].indx1_pars.mean
   foregrounds[0].indx2_pars.mean = 18.
   foregrounds[0].indx2_pars.sigma = 3.
   foregrounds[0].indx2_pars.lower = 5.
   foregrounds[0].indx2_pars.upper = 40.
   foregrounds[0].indx2_current = foregrounds[0].indx2_pars.mean

   
   foregrounds[1].type = 'line_spectrum'
   foregrounds[1].nametag = 'co'
   foregrounds[1].nu_ref = 101.28d0
   foregrounds[1].spectrum_filename = 'CO_spectrum_0.6_0.2_ffp4effreq.dat'

   foregrounds[2].type = 'power_law'
   foregrounds[2].nametag = 'low-freq'
   foregrounds[2].nu_ref = 30.d0
   foregrounds[2].indx1_pars.mean = -3.05
   foregrounds[2].indx1_pars.sigma = 0.3
   foregrounds[2].indx1_pars.lower = -5.
   foregrounds[2].indx1_pars.upper = -1.
   foregrounds[0].indx1_current = foregrounds[0].indx1_pars.mean

   foregrounds[3].type = 'cmb'
   foregrounds[3].nametag = 'cmb'
   foregrounds[3].nu_ref = 100.d0

   maps = ruler_read_maps(channels, Nside=Nside, verbosity=2)
;stop
   foreground_fit = ruler_read_foregrounds(foregrounds, maps.cfreq, Nside=Nside, verbosity=2, component_spectral_response=foreground_spectral_response)

help, maps
help, foregrounds
STOP
   for isample=1, Nsample do begin
       
       ; Sample amplitudes:
       result = func_ruler(channels=channels, foreground_components=foregrounds_fit, Nside=Nside, maskfile=maskfile, verbosity=2, do_check='f')
       ; Sample indices
       for icomp=0, Ncomp-1 do begin
           if (icomp ne icmb) then tmp_indx = ruler_sample_indices(channels, foreground_fit, icomp, verbosity=2)
       endfor
;##       result = func_ruler(channels=channels, foreground_components=foregrounds_fit, Nside=Nside, maskfile=maskfile, verbosity=2, do_check='f')

   endfor

stop

end


