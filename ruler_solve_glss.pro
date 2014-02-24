function ruler_solve_glss, channel_maps, $
                           components, $
                           gpix, $
                           component_spectral_response=component_spectral_response, $
                           verbosity=verbosity, $
                           do_check=do_check, $
                           weights=weights

; Ruler: linear solution in pixel space of the component separation
; problem
;
; D.Pietrobon Jul 2011
;
; GLSS solution
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

   if (not keyword_set(verbosity)) then verbosity = 0
   if (verbosity gt 0) then loud = True

   if (not keyword_set(do_check)) then begin
       print, 'do_check missing'
       do_check = 'f'
       print, 'Set to "f"'
       stop
   endif

   if (do_check eq 't') then check = true
   if (do_check eq 'f') then check = false

   do_sed_precompute = True
   if (not keyword_set(component_spectral_response) ) then do_sed_precompute = False

;   if (not keyword_set(const_rms)) then const_rms = 0

   code = ' --> Ruler <--'
   temp_folder = '/global/scratch/sd/dpietrob/rubbish/'

   Nfreq = n_elements(channel_maps)
   
;   Npix = 12l * Nside^2
   Npix = n_elements(channel_maps[0].map)
   
   Ngpix = n_elements(gpix)

   Ncomp = n_elements(components)

   a2t  = channel_maps[*].a2t
   freq = channel_maps[*].cfreq

   if (do_sed_precompute) then begin
       weights = fltarr(npix, Nfreq, Ncomp) ; IMPROVED
   endif

; --- System Solution
   print, ' Starting computation...'
   delta_pix = ngpix / 10
   t1 = systime(/seconds)

   for ipix = 0l, ngpix-1 do begin

       if ( (ipix/delta_pix)*delta_pix eq ipix) then print, ipix, ' / ', ngpix

       reddata = channel_maps[*].map[gpix[ipix]] / channel_maps[*].rms[gpix[ipix]]^2

       if (do_sed_precompute) then begin
           spectral_response = transpose( reform( component_spectral_response[gpix[ipix], *, *] ) )
       endif else begin
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
       endelse

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

   t2=systime(/seconds)
   print, ' - End of calculation. Elapsed time:', (t2-t1)/60., ' minutes.'

return, components

end
