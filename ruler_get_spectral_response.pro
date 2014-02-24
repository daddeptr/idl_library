function ruler_get_spectral_response, channel_maps=channel_maps, $
;                                      Nside = Nside, $
                                      components=components
;                     write_maps=write_maps, $
;                     tag=tag, $
;                     maskfile=maskfile, $
;                     first_pix_frac=first_pix_frac, $
;                     last_pix_frac=last_pix_frac, $
;                     verbosity=verbosity
;                     do_check=do_check, $
;                     const_rms = const_rms

; Ruler: linear solution in pixel space of the component separation
; problem
;
; D.Pietrobon Jul 2011
;
; Get spectral response of components given frequencies
 
   False = 0b
   True  = 1b

;   if (not keyword_set(verbosity)) then verbosity = 0
;   if (verbosity gt 0) then loud = True

   code = ' --> Ruler <--'
   temp_folder = '/global/scratch/sd/dpietrob/rubbish/'

;   do_sed_precompute = True
;   if (Nside gt 256l) then do_sed_precompute = False

   Ncomp = n_elements(components)

   Nfreq = n_elements(channel_maps)

   Npix = n_elements( channel_maps[0].map)

   a2t = channel_maps.a2t

   freq = channel_maps.cfreq

; ------ computing spectral response by frequency
;   if (do_sed_precompute) then begin
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
               
           endfor               ;component spectral response loop
       endfor                   ; frequency loop
;   endif
;##   weights = fltarr(npix, ncomp, ncomp) ; TO BE IMPROVED, SINCE DOESN'T GIVE LORIS' WEIGHTS
return, component_spectral_response

end
