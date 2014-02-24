function ruler_read_foregrounds, foreground_components, Freq, $
                                 Nside = Nside, $
                                 verbosity = verbosity, $
                                 component_spectral_response = component_spectral_response

; Ruler: linear solution in pixel space of the component separation
; problem
;
; Reading maps in
 
   False = 0b
   True  = 1b

   if (not keyword_set(verbosity)) then verbosity = 0
   if (verbosity gt 0) then loud = True

   code = ' --> Ruler <--'

   Npix = 12l * Nside^2

; ----------------------------------------------------------------------
   Nfreq = n_elements(Freq)
   a2t = conversionfactor(Freq, /antenna2thermo)

   Ncomp = n_elements(foreground_components)
   component_tags = strarr(Ncomp)
   fg_nu_ref = dblarr(Ncomp)

   index_pars = { mean: 0.d0, sigma:0.d0, lower:0.d0, upper:0.d0 }

   tmpcomp = { nu_ref:0.d0, $
               nametag:'', $
               type:'', $
               indx1: dblarr(Npix), $
               indx1_pars: index_pars, $
               indx1_constant: False, $
               indx2: dblarr(Npix), $
               indx2_pars: index_pars, $               
               indx2_constant: False, $
               approx_sed: dblarr(2*Nfreq,2), $
               fg_amp: dblarr(Npix) }

   components = replicate(tmpcomp, Ncomp)
   component_spectral_response = dblarr(Npix, Nfreq, Ncomp)

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

       if (component_tags[icomp] eq 'curved_power_law') then begin
           if ( (strlen(foreground_components[icomp].indx1_filename) gt 0) and (strlen(foreground_components[icomp].indx2_filename) gt 0) ) then begin
               read_fits_map, foreground_components[icomp].indx1_filename, indx1
               components[icomp].indx1 = indx1
               read_fits_map, foreground_components[icomp].indx2_filename, indx2
               components[icomp].indx2 = indx2
           endif else begin
               components[icomp].indx1[*] = foreground_components[icomp].indx1_pars.mean
               components[icomp].indx2[*] = foreground_components[icomp].indx2_pars.mean
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

return, components

end
