function func_ruler_v2, channels=channels, $
                     Nside = Nside, $
                     foreground_components=foreground_components, $
                     write_maps=write_maps, $
                     write_weights=write_weights, $
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

   if (not keyword_set(write_maps)) then write_maps = True
   if (not keyword_set(write_weights)) then write_weights = False

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
       print, 'Masked Applied: '+maskfile
   endelse

; --- Setting components
   Ncomp = n_elements(foreground_components)

   components = ruler_read_fg_components_in(foreground_components, Nside=Nside, verbosity=verbosity)

; --- Setting frequency maps
   Nfreq = n_elements(channels)

   print, 'Frequency Channels:'
   print, ' Number of frequencies used: ', Nfreq

   channel_maps = ruler_read_maps_in(channels, Nside=Nside, verbosity=verbosity) ;replicate(tmpmap, Nfreq)

   a2t  = channel_maps[*].a2t
   freq = channel_maps[*].cfreq

; --- computing spectral response by frequency
   if (do_sed_precompute) then component_spectral_response = ruler_get_spectral_response(channel_maps=channel_maps, components=components)

; --- System Solution
   print, ' Starting computation...'
   delta_pix = ngpix / 10
   t1 = systime(/seconds)

   if (const_rms eq 1) then begin
       print, ' - Setting RMS = 1.d0...'
       channel_maps.rms[*,*] = 1.d0
   endif

   if (do_sed_precompute) then components = ruler_solve_glss(channel_maps, components, gpix, component_spectral_response=component_spectral_response, verbosity=verbosity, do_check=do_check, weights=weights) else $
     components = ruler_solve_glss(channel_maps, components, gpix, verbosity=verbosity, do_check=do_check)

; --- Outputing the solution
   print, ' - Saving files...'
   for icomp=0,Ncomp-1 do begin
       if (write_maps) then begin
           write_fits_map, tag+'_ruler_components_'+components[icomp].nametag+'.fits', components[icomp].fg_amp, /ring, units='!7l!6K'
           if (verbosity gt 1) then mollview, tag+'_ruler_components_'+components[icomp].nametag+'.fits', chars=1.5, min=-300, max=300, px=500
       endif
       for ifreq=0,Nfreq-1 do begin
           if (do_sed_precompute and write_weights) then write_fits_map, tag+'_ruler_components_weights_'+components[icomp].nametag+'_'+channel_maps[ifreq].nametag+'.fits', weights[*,ifreq, icomp], /ring, units='!6[Dimensionless]'
           if (do_sed_precompute and write_weights and (verbosity gt 2)) then mollview, tag+'_ruler_components_weights_'+components[icomp].nametag+'_'+channel_maps[ifreq].nametag+'.fits', /hist, chars=1.5, px=500 
       endfor
   endfor
;return

end
