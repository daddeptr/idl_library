;## To be run from /global/scratch/sd/dpietrob/dx9/maps/ns0128/

;##   root  = 'cls' ;## 'pix_01a_v2' ;FFP4
   root = 'pix_7b' ;DX9
   sfreq = ['030','044','070','100','143','217','353', '545']

;##   freq  = [28.44, 44.12, 70.33, 101.28, 142.85, 222.51, 361.46]
;##   mfile = 'ffp4_01a_'+sfreq+'_ns128_uK.fits'
;##   rfile = 'ffp4_01a_srms_'+sfreq+'_v3.fits'

   freq  = [28.4, 44.1, 70.3, 101.2, 142.6, 221.9, 360.6, 557.6]
   mfile = 'dx9_Imap_'+sfreq+'_ns128_uK.fits'
   rfile = 'dx9_rms_'+sfreq+'.fits'
   co_file = 'CO_spectrum_0.6_0.2_effreq.txt'

   nfreq = n_elements(freq)

   Nside = 128l
   Npix  = 12l*Nside^2
   Ncomp = 4 ;6

;##   ssamp = '001'
;##   outtag = 'chains/glss/'+root+'/dx9_'+ssamp

   Nsamp = 1l ;201l
   Fsamp = 1l ;4l

   for isamp=Fsamp, Nsamp do begin
       if ( (Fsamp/25)*25 eq Fsamp ) then print, isamp, '/', Nsamp
       ssamp = string(isamp, format='(i3.3)')
       outtag = 'chains/glss/'+root+'/dx9_test8b'+ssamp
; ----------------------------------------------------------------------
       tmpmap = {nametag:'', $
                 cfreq: 0.d0, $
                 map_filename: '', $
                 rms_filename: '', $
                 monopole: 0.d0, $
                 dipole: dblarr(3) }
   
       channels = replicate(tmpmap, nfreq) 

       print, 'Warning: Mono/dipole set to 0!'
       stop
;##   readcol, 'chains/'+root+'/foreground_c0001_k00'+ssamp+'.dat', m, dx, dy, dz
       m = fltarr(Nfreq)
       dx = fltarr(Nfreq)
       dy = fltarr(Nfreq)
       dz = fltarr(Nfreq)
;##   for i=0,6 do print,m[i],dx[i],dy[i],dz[i]

       for ifreq=0,nfreq-1 do begin
           channels[ifreq].cfreq = freq[ifreq]
           channels[ifreq].nametag = sfreq[ifreq]
           channels[ifreq].map_filename = mfile[ifreq]
           channels[ifreq].rms_filename = rfile[ifreq]
           channels[ifreq].monopole = m[ifreq]
           channels[ifreq].dipole = [dx[ifreq], dy[ifreq], dz[ifreq]]
       endfor

; ----------------------------------------------------------------------
       tmpcomp = { nametag:'', $
                   type: '', $
                   nu_ref:0.d0, $
                   indx1_filename: '', $
                   indx2_filename: '', $
                   spectrum_filename: '' }
   
       foregrounds = replicate(tmpcomp,ncomp)

       foregrounds[0].type = 'thermal_dust'
       foregrounds[0].nametag = 'thermal_dust'
       foregrounds[0].nu_ref = 360.6d0
       foregrounds[0].indx1_filename = 'chains/'+root+'/fg_ind_map_no01_c0001_k00'+ssamp+'.fits'
       foregrounds[0].indx2_filename = 'chains/'+root+'/fg_ind_map_no02_c0001_k00'+ssamp+'.fits'
       
       foregrounds[1].type = 'line_spectrum'
       foregrounds[1].nametag = 'co'
       foregrounds[1].nu_ref = 101.2d0
       foregrounds[1].spectrum_filename = co_file

       foregrounds[2].type = 'power_law'
       foregrounds[2].nametag = 'ff'
       foregrounds[2].nu_ref = 28.4d0
       foregrounds[2].indx1_filename = 'chains/'+root+'/fg_ind_map_no03_c0001_k00'+ssamp+'.fits'

       foregrounds[3].type = 'cmb'
       foregrounds[3].nametag = 'cmb'
       foregrounds[3].nu_ref = 100.d0

;   foregrounds[4].type = 'power_law'
;   foregrounds[4].nametag = 'sync'
;   foregrounds[4].nu_ref = 0.408d0
;   foregrounds[4].indx1_filename = 'chains/'+root+'/fg_ind_map_no05_c0001_k00'+ssamp+'.fits'

;   foregrounds[5].type = 'curved_power_law'
;   foregrounds[5].nametag = 'ame'
;   foregrounds[5].nu_ref = 28.4d0
;   foregrounds[5].indx1_filename = 'chains/'+root+'/fg_ind_map_no03_c0001_k00'+ssamp+'.fits'
;   foregrounds[5].indx2_filename = 'chains/'+root+'/fg_ind_map_no04_c0001_k00'+ssamp+'.fits'

       result = func_ruler(channels=channels, foreground_components=foregrounds, Nside=Nside, maskfile='fsky.fits', verbosity=1, do_check='f', tag=outtag)

   endfor

stop

end


