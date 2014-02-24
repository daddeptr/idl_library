function ruler_sample_indices, channel_maps, foreground_components, icomp, $
                                 verbosity = verbosity;, $
;                               Ngibbs = Ngibbs

; Ruler: linear solution in pixel space of the component separation
; problem

 
   False = 0b
   True  = 1b

   if (not keyword_set(verbosity)) then verbosity = 0
   if (verbosity gt 0) then loud = True
   if (not keyword_set(Ngibbs)) then Ngibbs = 1
   
   code = ' --> Ruler <--'
   
   if ( (foreground_components[icomp].indx1_pars.sigma eq 0.) or $
        (foreground_components[icomp].indx1_pars.lower eq foreground_components[icomp].indx1_pars.upper) ) then begin
       if loud then print, ' Sampling not required'
       return
   endif

   Nstep = 100l
   Ncomp = n_elements(foreground_components)

   freq = channel_maps[*].cfreq
   Nfreq = n_elements(channel_maps)
   Npix = n_elements(channel_maps[0].map)
   If loud then print, ' Nside = ', npix2nside(Npix)

   icmb = where(foreground_components.type eq 'cmb')
   if loud then print, ' iCMB = ', icmb

   indices = dblarr(Npix,2)

   sky_model = dblarr(Npix, Nfreq)
; adding components that are not sampled
   component_spectral_response = dblarr(Nfreq, Ncomp)
   for jcomp=0, Ncomp-1 do begin
; ------ Spectral response
;        ##       component_spectral_response = dblarr(Npix, Nfreq,
;        Ncomp) Not necessary for now, since indices are assumed constant
       if (jcomp ne icomp) then begin
           fg_type = components[jcomp].type
           for ifreq=0, Nfreq-1 do begin
               case fg_type of
                   'thermal_dust': begin
                       emis = foreground_components[jcomp].indx1_current
                       temp = foreground_components[jcomp].indx2_current
                       nu_ref = foreground_components[jcomp].nu_ref * 1.d9
                       component_spectral_response[ifreq,jcomp] = a2t[ifreq] * gray_body(freq[ifreq]*1.d9, emis, temp, nu_ref=nu_ref)
                   end
           
                   'power_law': begin
                       beta = foreground_components[jcomp].indx1_current
                       nu_ref = foreground_components[jcomp].nu_ref * 1.d9
                       component_spectral_response[ifreq,jcomp] = a2t[ifreq] * power_law(freq[ifreq]*1.d9, beta, nu_ref=nu_ref)
                   end
                   
                   'line_spectrum': begin
;##                   sed = dblarr(Npix)+1.d0 Not necessary for now,
;since indices are assumed constant
                       sed = 1.d0
                       approx_sed = foreground_components[jcomp].approx_sed
                       nu_ref = foreground_components[jcomp].nu_ref * 1.d9
                       component_spectral_response[ifreq,jcomp] = a2t[ifreq] * line_spectrum(freq[ifreq]*1.d9, approx_sed, nu_ref=nu_ref) * sed
                   end
           
                   'cmb': begin
;##                   sed = dblarr(Npix) + 1.d0 Not necessary for now,
;since indices are assumed constant
                       sed = 1.d0
                       component_spectral_response[ifreq,jcomp] = sed
                   end
               endcase
               sky_model[*,ifreq] = sky_model[*,ifreq] + foreground_components[jcomp].amp[*] * component_spectral_response[ifreq,jcomp]
           endfor               ;Nfreq
       endif ; Check whether sampling component
   endfor ; Ncomp

   endfor

   chisq = fltarr(Nstep)

; ## Sampling THE component   
   fg_type = foreground_components[icomp].type
; ------ Index 1
   if (foreground_components[icomp].indx1_constant) then begin
       lower = foreground_components[icomp].indx1_pars.lower
       upper = foreground_components[icomp].indx1_pars.upper
       indx1 = dblarr(Nstep-1)/Nstep*(upper-lower)+lower
   endif else begin
       stop, 'Spatially varying indices not implemented yet'
   endif
; ------ Index 2
   if (foreground_components[icomp].indx2_constant) then begin
       lower = foreground_components[icomp].indx2_pars.lower
       upper = foreground_components[icomp].indx2_pars.upper
       indx2 = dblarr(Nstep-1)/Nstep*(upper-lower)+lower
   endif else begin
       stop, 'Spatially varying indices not implemented yet'
   endif


; ------ Sampling Index 1
   for istep=0,Nstep-1 do begin
       if (loud) then print, fg_nu_ref[icomp], component_tags[icomp]
; ------ Spectral response
;##               component_spectral_response = dblarr(Npix, Nfreq, Ncomp)
       fg_type = components[icomp].type
       sampled_fg = dblarr(Npix, Nfreq)
       for ifreq=0, Nfreq-1 do begin
           case fg_type of
               'thermal_dust': begin
                   emis = indx1[istep]
                   temp = foreground_components[icomp].indx2_current
                   nu_ref = foreground_components[icomp].nu_ref * 1.d9
                   component_spectral_response[ifreq,icomp] = a2t[ifreq] * gray_body(freq[ifreq]*1.d9, emis, temp, nu_ref=nu_ref)
               end
           
               'power_law': begin
                   beta = indx1[istep]
                   nu_ref = foreground_components[icomp].nu_ref * 1.d9
                   component_spectral_response[ifreq,icomp] = a2t[ifreq] * power_law(freq[ifreq]*1.d9, beta, nu_ref=nu_ref)
               end
               
               'line_spectrum': begin
                   print, ' Warning: no line sampling implemented'
;                   sed = dblarr(Npix)+1.d0
;                   approx_sed = foreground_components[icomp].approx_sed
;                   nu_ref = foreground_components[icomp].nu_ref * 1.d9
;                   component_spectral_response[ifreq,icomp] = a2t[ifreq] * line_spectrum(freq[ifreq]*1.d9, approx_sed, nu_ref=nu_ref) * sed
               end
           endcase
           sampled_fg[*, ifreq] = foreground_components[icomp].amp * component_spectral_response[ifreq,icomp]
       endfor                   ;Nfreq
; Call get_frequency_maps
       chisq[istep] = chisq[istep] + total( (channel_maps[ifreq].map-(sky_model[*,ifreq]+sample_fg[*,ifreq]) )^2/channel_maps[ifreq].rms^2 )
   endfor ; Nstep

   
   plot, chisq
   mn = min(chisq, imn)
   foreground_components[icomp].indx1_current = indx1[imn]

; ------ Sampling Index 2
   for istep=0,Nstep-1 do begin
       if (loud) then print, fg_nu_ref[icomp], component_tags[icomp]
; ------ Spectral response
;##               component_spectral_response = dblarr(Npix, Nfreq, Ncomp)
       fg_type = components[icomp].type
       sampled_fg = dblarr(Npix, Nfreq)
       for ifreq=0, Nfreq-1 do begin
           case fg_type of
               'thermal_dust': begin
                   emis = foreground_components[icomp].indx1_current
                   temp = indx2[istep]
                   nu_ref = foreground_components[icomp].nu_ref * 1.d9
                   component_spectral_response[ifreq,icomp] = a2t[ifreq] * gray_body(freq[ifreq]*1.d9, emis, temp, nu_ref=nu_ref)
               end
           
               'power_law': begin
                   print, ' Warning: Power law has got one index only'
;##                   beta = foreground_components[icomp].indx1_current
;##                   curv = indx2[istep]
;##                   nu_ref = foreground_components[icomp].nu_ref * 1.d9
;##                   component_spectral_response[ifreq,icomp] = a2t[ifreq] * power_law(freq[ifreq]*1.d9, beta, nu_ref=nu_ref, curv=curv)
               end
               
               'line_spectrum': begin
                   print, ' Warning: no line sampling implemented'
;                   sed = dblarr(Npix)+1.d0
;                   approx_sed = foreground_components[icomp].approx_sed
;                   nu_ref = foreground_components[icomp].nu_ref * 1.d9
;                   component_spectral_response[ifreq,icomp] = a2t[ifreq] * line_spectrum(freq[ifreq]*1.d9, approx_sed, nu_ref=nu_ref) * sed
               end
           endcase
           sampled_fg[*, ifreq] = foreground_components[icomp].amp * component_spectral_response[ifreq,icomp]
       endfor                   ;Nfreq
; Call get_frequency_maps
       chisq[istep] = chisq[istep] + total( (channel_maps[ifreq].map-(sky_model[*,ifreq]+sample_fg[*,ifreq]) )^2/channel_maps[ifreq].rms^2 )
   endfor ; Nstep

   
   plot, chisq
   mn = min(chisq, imn)
   foreground_components[icomp].indx2_current = indx2[imn]

stop ; Done until here

;   Npix = 12l * Nside^2

; ----------------------------------------------------------------------
;   Nfreq = n_elements(Freq)
;   a2t = conversionfactor(Freq, /antenna2thermo)

;   Ncomp = n_elements(foreground_components)
;   component_tags = strarr(Ncomp)
;   fg_nu_ref = dblarr(Ncomp)
;
;   index_pars = { mean: 0.d0, sigma:0.d0, lower:0.d0, upper:0.d0 }
;
;   tmpcomp = { nu_ref:0.d0, $
;               nametag:'', $
;               type:'', $
;               indx1: dblarr(Npix), $
;               indx1_pars: index_pars, $
;               indx1_constant: False, $
;               indx2: dblarr(Npix), $
;               indx2_pars: index_pars, $               
;               indx2_constant: False, $
;               approx_sed: dblarr(2*Nfreq,2), $
;               fg_amp: dblarr(Npix) }
;
;   components = replicate(tmpcomp, Ncomp)
;   component_spectral_response = dblarr(Npix, Nfreq, Ncomp)

   for icomp=0, ncomp-1 do begin
       fg_nu_ref[icomp] = foreground_components[icomp].nu_ref
       component_tags[icomp] = foreground_components[icomp].type
       components[icomp].nu_ref = fg_nu_ref[icomp]
       components[icomp].type = foreground_components[icomp].type
       components[icomp].nametag = foreground_components[icomp].nametag
       
       components[icomp].indx1_pars = foreground_components[icomp].indx1_pars
       components[icomp].indx2_pars = foreground_components[icomp].indx2_pars

       components[icomp].indx1_constant = foreground_components[icomp].indx1_constant
       components[icomp].indx2_constant = foreground_components[icomp].indx1_constant

       if (component_tags[icomp] eq 'thermal_dust') then begin
           if ((strlen(foreground_components[icomp].indx1_filename) gt 0) and (strlen(foreground_components[icomp].indx2_filename) gt 0)) then begin
               read_fits_map, foreground_components[icomp].indx1_filename, indx1
               read_fits_map, foreground_components[icomp].indx2_filename, indx2
               components[icomp].indx1 = indx1
               components[icomp].indx2 = indx2
           endif else begin
               components[icomp].indx1[*] = foreground_components[icomp].indx1_pars.mean
               components[icomp].indx2[*] = foreground_components[icomp].indx2_pars.mean
           endelse
       endif 

       if (component_tags[icomp] eq 'power_law') then begin
           if (strlen(foreground_components[icomp].indx1_filename) gt 0) then begin
               read_fits_map, foreground_components[icomp].indx1_filename, indx1
               components[icomp].indx1 = indx1
           endif else begin
               components[icomp].indx1[*] = foreground_components[icomp].indx1_pars.mean
           endelse
       endif 

       if (component_tags[icomp] eq 'line_spectrum') then begin
           readcol, foreground_components[icomp].spectrum_filename, f, s
           for inu=0,n_elements(f)-1 do print, f[inu], s[inu]
           tmp_approx_sed = dblarr(n_elements(f),2)
           tmp_approx_sed[*,0] = f*1.d9
           tmp_approx_sed[*,1] = s
           components[icomp].approx_sed[0:n_elements(f)-1,*] = tmp_approx_sed[0:n_elements(f)-1,*]
       endif

       if (loud) then print, fg_nu_ref[icomp], component_tags[icomp]
; ------ Spectral response
       fg_type = components[icomp].type
       for ifreq=0, Nfreq-1 do begin
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
       endfor
   endfor

; ----------------------------------------------------------------------

return, indices

end
