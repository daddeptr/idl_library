function ruler_read_fg_components_in, foreground_components, $
                                      Nside=Nside, $
                                      verbosity=verbosity

; Ruler: linear solution in pixel space of the component separation
; problem
;
; D.Pietrobon Jul 2011
;
; Reading foreground components in
;
; --- Foreground Model definition
 
   False = 0b
   True  = 1b

   if (verbosity gt 0) then loud = True

   if (not keyword_set(const_rms)) then const_rms = 0

   code = ' --> Ruler <--'
   temp_folder = '/global/scratch/sd/dpietrob/rubbish/'

;##   nfreq = n_elements(channels)
   
   Npix = 12l * Nside^2

   Ncomp = n_elements(foreground_components)
   component_tags = strarr(Ncomp)

   tmpcomp = { nu_ref:0.d0, $
               nametag:'', $
               type:'', $
               indx1: dblarr(Npix), $
               indx2: dblarr(Npix), $
               approx_sed: dblarr(20,2), $
               fg_amp: dblarr(Npix) }

   components = replicate(tmpcomp, Ncomp)

   for icomp=0, ncomp-1 do begin
       component_tags[icomp] = foreground_components[icomp].type
       components[icomp].nu_ref = foreground_components[icomp].nu_ref
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

       if (loud) then print, components[icomp].nu_ref, ' ', component_tags[icomp]
   endfor

return, components

end

