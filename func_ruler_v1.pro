function func_ruler_v1, channels=channels, $
                     Nside = Nside, $
                     foreground_components=foreground_components, $
                     write_maps=write_maps, $
                     tag=tag, $
                     maskfile=maskfile, $
                     first_pix_frac=first_pix_frac, $
                     last_pix_frac=last_pix_frac, $
                     verbosity=verbosity, $
                     do_check=do_check, $
                     const_rms = const_rms

; Ruler: linear solution in pixel space of the component separation
; problem
;
; D.Pietrobon Jul 2011
;
; To be improved and generalized to be flexible on foreground model
; and translated into f90 for speed
;
; chisq = sum_ch (D_ch-M_ch(A_i))^2/N_ch^2
;
; dchisq/dA_i = 0 = \sum_ch 1/N_ch^2 (D_ch-\sum_j A_j*Ind_ch,j)*I_ch,i
;
; M # A = B
; M = \sum_ch Ind[ch,i] * Ind[ch,j] N[ch]^-1
; M = Ind ## N^-1 ## Ind^t
; --> M^-1
;
; A = (cmb, synch, ff, dust)
;
; B = \sum_ch D_ch/N_ch^2 * Ind_ch ; Ind_ch array of spectral responses
; B = \sum D[ch] * Ind[ch,i] * N[ch]^1 
; B = Ind ## N^-1 ## D
;
; A = B # M^-1
;
; A_i = \sum_j B_j * (M^-1)_ij 
;     = \sum_j (\sum_ch D[ch]*Ind[ch,j] / N[ch]^2) * (\sum_chp 1./N[chp]^2 * Ind[chp] * Ind[chp])^-1_ij =
;     = \sum_ch D[ch] (\sum_j Ind[ch]/N[ch]^2)*(\sum_chp 1./N_chp^2*Ind_chp*Ind_chp)^-1_ij 
;
; A = (Ind # N^-1 # Ind^t)^-1 * (Ind # N^-1 D^t) --> W = (Ind # N^-1 # Ind)^1*(Ind # N^-1 ) 
;
; W = (Ind # N^-1 # Ind)^-1 * (Ind # N^-1 )
;
; Den = \sum_Nfreq spectral_response[Nfreq, Ncomp] ## N[Nfreq]^-1 ## transpose( spectral_response[Nfreq, Ncomp] )
; Num = spectral_response[Nfreq, Ncomp] # N^-1
;
; --- Foreground Model definition
 
   False = 0b
   True  = 1b

   if (not keyword_set(maskfile)) then begin
       full_sky = True
   endif else begin
       read_fits_map, maskfile, mask, nside=masknside, order=maskorder
       if ( maskorder eq 'NESTED' ) then mask = reorder(mask, /n2r)
       gpix = where(mask gt 0.)
       ngpix = n_elements(gpix)
       full_sky = False
       print, 'Fsky = ', float(ngpix)/nside2npix(masknside)
   endelse
   if (not keyword_set(first_pix_frac)) then first_pix_frac = 0. 
   if (not keyword_set(last_pix_frac)) then last_pix_frac = 1. 
   if (not keyword_set(verbosity)) then verbosity = 0
   if (verbosity gt 0) then loud = True

   if ( not keyword_set(tag) ) then tag='test_IDL_ruler'

   if (not keyword_set(do_check)) then begin
       print, 'do_check missing'
       do_check = 'f'
       print, 'Set to "f"'
       stop
   endif

   if (do_check eq 't') then check = true
   if (do_check eq 'f') then check = false

   if (not keyword_set(const_rms)) then const_rms = 0

   code = ' --> Ruler <--'
   temp_folder = '/global/scratch/sd/dpietrob/rubbish/'

   nfreq = n_elements(channels)
   
   Npix = 12l * Nside^2

   do_sed_precompute = True
   if (Nside gt 256l) then do_sed_precompute = False

   if (full_sky) then begin
       print, 'First/Last pixel fraction: ', first_pix_frac, last_pix_frac
       fpix = long( (npix)*first_pix_frac )
       lpix = long( (npix)*last_pix_frac )-1
       print, 'First/Last pixel: ', fpix, lpix
       gpix = lindgen(lpix-fpix+1)+fpix
       ngpix = n_elements(gpix)
   endif else begin
       print, 'Masked Applied'
   endelse

; ----------------------------------------------------------------------
   Ncomp = n_elements(foreground_components)
   component_tags = strarr(Ncomp)
;##   fg_nu_ref = dblarr(Ncomp)

   tmpcomp = { nu_ref:0.d0, $
               nametag:'', $
               type:'', $
               indx1: dblarr(Npix), $
               indx2: dblarr(Npix), $
               approx_sed: dblarr(2*Nfreq,2), $
               fg_amp: dblarr(Npix) }

   components = replicate(tmpcomp, Ncomp)

   for icomp=0, ncomp-1 do begin
;##       fg_nu_ref[icomp] = foreground_components[icomp].nu_ref
       component_tags[icomp] = foreground_components[icomp].type
       components[icomp].nu_ref = foreground_components[icomp].nu_ref ;##fg_nu_ref[icomp]
       components[icomp].type = foreground_components[icomp].type
       components[icomp].nametag = foreground_components[icomp].nametag

       if (component_tags[icomp] eq 'thermal_dust') then begin
           read_fits_map, foreground_components[icomp].indx1_filename, indx, nside=dpns, order=dpord
           if (dpord eq 'NESTED') then indx = reorder(indx, in=dpord, out='RING')
           if (dpns ne Nside) then indx = upgrade_index_ring(indx, Nside)
           if (dpns ne Nside) then print, ' Index upgrade performed...',component_tags[icomp]
           components[icomp].indx1 = indx
           indx = 0.

           read_fits_map, foreground_components[icomp].indx2_filename, indx
           if (dpord eq 'NESTED') then indx = reorder(indx, in=dpord, out='RING')
           if (dpns ne Nside) then indx = upgrade_index_ring(indx, Nside)
           if (dpns ne Nside) then print, ' Index upgrade performed...',component_tags[icomp]
           components[icomp].indx2 = indx
           indx = 0.
       endif 

       if (component_tags[icomp] eq 'power_law') then begin
           read_fits_map, foreground_components[icomp].indx1_filename, indx, nside=dpns, order=dpord
           if (dpord eq 'NESTED') then indx = reorder(indx, in=dpord, out='RING')
           if (dpns ne Nside) then indx = upgrade_index_ring(indx, Nside)
           if (dpns ne Nside) then print, ' Index upgrade performed...',component_tags[icomp]
           components[icomp].indx1 = indx
           indx = 0.
       endif 

       if (component_tags[icomp] eq 'curved_power_law') then begin
           read_fits_map, foreground_components[icomp].indx1_filename, indx, nside=dpns, order=dpord
           if (dpord eq 'NESTED') then indx = reorder(indx, in=dpord, out='RING')
           if (dpns ne Nside) then indx = upgrade_index_ring(indx, Nside)
           if (dpns ne Nside) then print, ' Index upgrade performed...',component_tags[icomp]
           components[icomp].indx1 = indx
           indx = 0.

           read_fits_map, foreground_components[icomp].indx2_filename, indx
           if (dpord eq 'NESTED') then indx = reorder(indx, in=dpord, out='RING')
           if (dpns ne Nside) then indx = upgrade_index_ring(indx, Nside)
           if (dpns ne Nside) then print, ' Index upgrade performed...',component_tags[icomp]
           components[icomp].indx2 = indx
           indx = 0.
       endif 

       if (component_tags[icomp] eq 'line_spectrum') then begin
           readcol, foreground_components[icomp].spectrum_filename, f, s
           for inu=0,n_elements(f)-1 do print, f[inu], s[inu]
           tmp_approx_sed = dblarr(n_elements(f),2)
           tmp_approx_sed[*,0] = f*1.d9
           tmp_approx_sed[*,1] = s
           components[icomp].approx_sed[0:n_elements(f)-1,*] = tmp_approx_sed[0:n_elements(f)-1,*]
       endif

;##       if (loud) then print, fg_nu_ref[icomp], component_tags[icomp]
       if (loud) then print, components[icomp].nu_ref, ' ', component_tags[icomp]
   endfor

; ----------------------------------------------------------------------
   Nfreq = n_elements(channels)
   freq = dblarr(Nfreq)
   a2t = dblarr(Nfreq)
   print, 'Frequency Channels:'
   print, ' Number of frequencies used: ', Nfreq

   tmpmap = { nametag:'', $
              cfreq: 0.d0, $
              map: dblarr(Npix), $
              rms: dblarr(Npix), $
              a2t: 0.d0, $
              monopole: 0.d0, $
              dipole: dblarr(3) }

   channel_maps = replicate(tmpmap, Nfreq)

   for ifreq=0,nfreq-1 do begin
       freq[ifreq] = channels[ifreq].cfreq
       a2t[ifreq] = conversionfactor(freq[ifreq], /antenna2thermo)

       channel_maps[ifreq].nametag = channels[ifreq].nametag
       channel_maps[ifreq].cfreq = channels[ifreq].cfreq
       channel_maps[ifreq].a2t = conversionfactor(freq[ifreq], /antenna2thermo)
       channel_maps[ifreq].monopole = channels[ifreq].monopole
       channel_maps[ifreq].dipole = channels[ifreq].dipole

       if (loud) then print, ifreq, '- ', freq[ifreq]
; ------ reading map
       if (loud) then print, 'reading ', channels[ifreq].map_filename 
       read_fits_map, channels[ifreq].map_filename, map, nside=mapnside, order=maporder
       if (mapnside ne Nside) then STOP, 'Nside mismatch: ', channels[ifreq].nametag, mapnside
       if (maporder eq 'NESTED') then begin
           if (loud) then print, 'reordering ', channels[ifreq].nametag
           map = reorder(map, /n2r)
       endif
       if (channel_maps[ifreq].monopole ne 0.d0) then map = map-channel_maps[ifreq].monopole 
       if ( (abs(channel_maps[ifreq].dipole[0]) gt 1.d-11) and (abs(channel_maps[ifreq].dipole[1]) gt 1.d-11) and (abs(channel_maps[ifreq].dipole[2]) gt 1.d-11) ) then begin
           dipole = make_dipole(Nside, channel_maps[ifreq].dipole)
           if (verbosity gt 2) then mollview, dipole, win=ifreq, tit=channels[ifreq].nametag
           map = map - dipole
           dipole = 0.
       endif
       channel_maps[ifreq].map = map
       map = 0.

; ------ reading rms
       if (loud) then print, 'reading ', channels[ifreq].rms_filename 
       read_fits_map, channels[ifreq].rms_filename, rms, nside=mapnside, order=maporder
       if (mapnside ne Nside) then STOP, 'Nside mismatch: ', channels[ifreq].nametag, mapnside
       if (maporder eq 'NESTED') then begin
           if (loud) then print, 'reordering ', channels[ifreq].nametag
           rms = reorder(rms, /n2r)
       endif
       channel_maps[ifreq].rms = rms
       rms = 0.
       
       if (loud) then print, freq[ifreq], a2t[ifreq]

   endfor ;frequency loop

; ----------------------------------------------------------------------

; ------ computing spectral response by frequency
if (do_sed_precompute) then begin
   component_spectral_response = dblarr(Npix, Nfreq, Ncomp)

   for ifreq=0,nfreq-1 do begin
       for icomp=0, Ncomp-1 do begin
           fg_type = components[icomp].type
           case fg_type of
               'thermal_dust': begin
                   emis = components[icomp].indx1
                   temp = components[icomp].indx2
                   nu_ref = components[icomp].nu_ref * 1.d9
                   component_spectral_response[*,ifreq,icomp] = a2t[ifreq] * gray_body(freq[ifreq]*1.d9, emis, temp, nu_ref=nu_ref)
               end

               'power_law': begin
                   beta = components[icomp].indx1
                   nu_ref = components[icomp].nu_ref * 1.d9
                   component_spectral_response[*,ifreq,icomp] = a2t[ifreq] * power_law(freq[ifreq]*1.d9, beta, nu_ref=nu_ref)
               end
               
               'curved_power_law': begin
                   beta = components[icomp].indx1
                   curv = components[icomp].indx2
                   nu_ref = components[icomp].nu_ref * 1.d9
                   component_spectral_response[*,ifreq,icomp] = a2t[ifreq] * power_law(freq[ifreq]*1.d9, beta, nu_ref=nu_ref, curv=curv)
               end
               
               'line_spectrum': begin
                   sed = dblarr(Npix)+1.d0
                   approx_sed = components[icomp].approx_sed
                   nu_ref = components[icomp].nu_ref * 1.d9
                   component_spectral_response[*,ifreq,icomp] = a2t[ifreq] * line_spectrum(freq[ifreq]*1.d9, approx_sed, nu_ref=nu_ref) * sed
               end

               'cmb': begin
                   sed = dblarr(Npix) + 1.d0
                   component_spectral_response[*,ifreq,icomp] = sed
               end
           endcase
           
       endfor ;component spectral response loop
   endfor; frequency loop
endif
;##   weights = fltarr(npix, ncomp, ncomp) ; TO BE IMPROVED, SINCE DOESN'T GIVE LORIS' WEIGHTS
   if (do_sed_precompute) then weights = fltarr(npix, Nfreq, Ncomp) ; IMPROVED
;##   map = channel_maps.map
;##   rms = channel_maps.rms
;##   reduced_data = channel_maps.map / channel_maps.rms^2

; --- System Solution
   print, ' Starting computation...'
   delta_pix = ngpix / 10
   t1 = systime(/seconds)

   if (const_rms eq 1) then begin
       print, ' - Setting RMS = 1.d0...'
       channel_maps.rms[*,*] = 1.d0
   endif

   for ipix = 0l, ngpix-1 do begin

       if ( (ipix/delta_pix)*delta_pix eq ipix) then print, ipix, ' / ', ngpix

;##       reddata = channel_maps.map[gpix[ipix],*] / channel_maps.rms[gpix[ipix],*]^2 ;##reduced_data[gpix[ipix],*]
       reddata = channel_maps[*].map[gpix[ipix]] / channel_maps[*].rms[gpix[ipix]]^2 ;##reduced_data[gpix[ipix],*]

       if (do_sed_precompute) then spectral_response = transpose( reform( component_spectral_response[gpix[ipix], *, *] ) )
       if (not do_sed_precompute) then begin
           spectral_response = dblarr(Nfreq, Ncomp)

           for ifreq=0l,nfreq-1 do begin
               for icomp=0l, Ncomp-1 do begin
                   fg_type = components[icomp].type
                   case fg_type of
                       'thermal_dust': begin
                           emis = components[icomp].indx1[gpix[ipix]]
                           temp = components[icomp].indx2[gpix[ipix]]
                           nu_ref = components[icomp].nu_ref * 1.d9
                           spectral_response[ifreq,icomp] = a2t[ifreq] * gray_body(freq[ifreq]*1.d9, emis, temp, nu_ref=nu_ref)
                       end

                       'power_law': begin
                           beta = components[icomp].indx1[gpix[ipix]]
                           nu_ref = components[icomp].nu_ref * 1.d9
                           spectral_response[ifreq,icomp] = a2t[ifreq] * power_law(freq[ifreq]*1.d9, beta, nu_ref=nu_ref)
                       end
               
                       'curved_power_law': begin
                           beta = components[icomp].indx1[gpix[ipix]]
                           curv = components[icomp].indx2[gpix[ipix]]
                           nu_ref = components[icomp].nu_ref * 1.d9
                           spectral_response[ifreq,icomp] = a2t[ifreq] * power_law(freq[ifreq]*1.d9, beta, nu_ref=nu_ref, curv=curv)
                       end
               
                       'line_spectrum': begin
                           sed = 1.d0
                           approx_sed = components[icomp].approx_sed
                           nu_ref = components[icomp].nu_ref * 1.d9
                           spectral_response[ifreq,icomp] = a2t[ifreq] * line_spectrum(freq[ifreq]*1.d9, approx_sed, nu_ref=nu_ref) * sed
                       end

                       'cmb': begin
                           sed = 1.d0
                           spectral_response[ifreq,icomp] = sed
                       end
                   endcase
           
               endfor           ;component spectral response loop
           endfor               ; frequency loop
           spectral_response = transpose(spectral_response)
       endif

       B = dblarr(ncomp)
       A = dblarr(ncomp, ncomp)

       B = transpose( spectral_response ) ## reform( reddata )
       B = reform( B )

;help, spectral_response
;help, reddata
;help, B

;stop
       for i=0l,ncomp-1 do begin
           for j=i,ncomp-1 do begin
;##               A[i,j] = total( spectral_response[i,*]*spectral_response[j,*] / channel_maps.rms[gpix[ipix],*]^2 )
               A[i,j] = total( spectral_response[i,*]*spectral_response[j,*] / channel_maps[*].rms[gpix[ipix]]^2 )
;##               for ifreq=0,nfreq-1 do begin
;##                   A[i,j] = A[i,j] + spectral_response[i,ifreq]*spectral_response[j,ifreq] / rms[gpix[ipix],ifreq]^2
;##               endfor
               A[j,i] = A[i,j]
           endfor
       endfor
       
       Am1 = la_invert(A,/double, status=info)

; ------ Computing weights
       if (do_sed_precompute) then begin
;help, weights
;help, spectral_response
;help, channel_maps.rms
;help, Am1
;stop
           for ifreq=0l, Nfreq-1 do for icomp=0l,Ncomp-1 do weights[gpix[ipix], ifreq, icomp] = total( reform(spectral_response[*,ifreq]) / channel_maps[ifreq].rms[gpix[ipix]]^2 * reform(Am1[icomp,*]) )
       endif

       if (check) then begin
           print, 'check:'
           print, 'A ', A
           print, ''
           print, 'Am1 ', Am1
           eval = EIGENQL(A, /double, EIGENVECTORS=evec)
           print, ''
           print, 'EVEC ', evec
           print, ''
           print, 'EVAL', eval

           evalm1 = HQR(ELMHES(Am1), /DOUBLE)
           evecm1 = EIGENVEC(Am1, evalm1, /double)
           print, ''
           print, 'EVECm1 ', evecm1
           print, ''
           print, 'EVALm1', evalm1
           stop
       endif

       X = Am1 ## transpose(B)
       X = reform( X )

;       print, X
; ------ Using weights : Check!   
;##       Xx = reform( weights[gpix[ipix], *, *] ) ## reform(map[gpix[ipix],*])
;##       print, reform(Xx)
;##stop

       if (info /= 0) then begin
           print, 'Pixel badly conditioned:'
           mask[gpix[ipix]] = 0
           print, gpix[ipix], info
       endif

       components[*].fg_amp[gpix[ipix]] = X

   endfor ; pixel loop

;cputime, t2, ts
;if (loud) then print, t2-t1, ' Sec'

;   stop
   t2=systime(/seconds)
   print, ' - End of calculation. Elapsed time:', (t2-t1)/60., ' minutes.'
   print, ' - Saving files...'
   for icomp=0,Ncomp-1 do begin
       write_fits_map, tag+'_ruler_components_'+components[icomp].nametag+'.fits', components[icomp].fg_amp, /ring, units='!7l!6K'
       for ifreq=0,Nfreq-1 do begin
           if (do_sed_precompute) then write_fits_map, tag+'_ruler_components_weights_'+components[icomp].nametag+'_'+channel_maps[ifreq].nametag+'.fits', weights[*,ifreq, icomp], /ring, units='!6[Dimensionless]'
           if (do_sed_precompute and (verbosity gt 2)) then mollview, tag+'_ruler_components_weights_'+components[icomp].nametag+'_'+channel_maps[ifreq].nametag+'.fits', /hist, chars=1.5, px=500 
       endfor
       if (verbosity gt 1) then mollview, tag+'_ruler_components_'+components[icomp].nametag+'.fits', chars=1.5, min=-300, max=300, px=500
   endfor
;return

end
