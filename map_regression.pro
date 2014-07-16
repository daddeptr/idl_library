function map_regression, template_files, $
                         label=label, $
                         fit=fit, $
                         maskfile=maskfile, $
                         A_coeff=A, $
                         a_err=a_err, $
                         do_plot=do_plot, $
                         md_rem=md_rem, $
                         rms_file = rms_file, $
                         rms_map = rms_map, $
                         inmap=inmap, $
                         infile = infile, $
                         cut = cut, $
                         scale = scale, $
                         debug = debug

;##   !path=!path+':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'
;##   !path=!path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'

   if (not keyword_set(label)) then label = '!6Residuals'
   if (not keyword_set(do_plot)) then do_plot = 'no'
   if (not keyword_set(maskfile)) then begin
       maskfile = '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/wmap/ns128/wmap7_temperature_analysis_mask_ns128.fits'
       print, 'WARNING: maskfile not provided. Set to: '+maskfile
   endif else begin
       print, ' Using mask:', maskfile
   endelse
   
   true  = 1b
   false = 0b

   if (keyword_set(infile) and (not keyword_set(inmap))) then begin
       read_fits_map, infile, inmap, order=ord
       if ord ne 'RING' then inmap = reorder(inmap,in=ord, out='RING')
   endif
   sz = size(inmap[*,0])
   npix = sz[1]
   nside = sqrt( float(npix) / 12. )

   read_fits_map, maskfile, mask, order=dpord, nside=dpns
   if (dpord eq 'NESTED') then rmask=reorder(mask[*,0], in=dpord, out='RING') else rmask=mask[*,0]
   mask = rmask[*,0]
   if (dpns ne nside) then begin
       ud_grade, mask, maskd, nside_out=nside, order_in='ring'
       mask = maskd
   endif
   if (keyword_set(cut)) then begin
       c=make_sky_cut(cut,nside)
       mask = mask * c
   endif
;stop

   bp = where(mask[*,0] eq 0.)
   gp = where(mask[*,0] eq 1.)
   ngp = n_elements(gp)
   print, total(mask,/double)/n_elements(mask[*,0]), double( n_elements(where(mask[*,0] eq 1)) )/n_elements(mask[*,0]), double( n_elements(where(mask[*,0] gt 1)) )/n_elements(mask[*,0])
;stop
   map = reform( inmap[*,0] )                  ;amp[*,iamp]
   if (keyword_set(md_rem)) then remove_dipole, map, mask[*,0], nside=nside, ordering='ring'; else remove_dipole, map, mask[*,0], nside=nside, ordering='ring', /onlymonopole
   if not keyword_set(rms_map) then rms = map * 0.+1 else rms = rms_map
   if (keyword_set(rms_file)) then read_fits_map, rms_file, rms 

   cc = map
   cc[bp] = -1.6375e30
   if (do_plot ne 'no') then mollview, cc, chars=1.5, min=-300, max=300, tit='!6Input map', win=1, px=650

   nfg = n_elements(template_files)
   print, 'nfg = ', nfg
   A = dblarr(nfg)
   oneOsigma = dblarr(nfg, nfg)
   sigma = oneOsigma
   a_err = dblarr(nfg)

   temp = dblarr(npix,nfg)
   for ifg=0, nfg-1 do begin
       read_fits_map, template_files[ifg], t, order=dpord
       if (dpord eq 'NESTED') then tr = reorder(t[*,0],in=dpord, out='ring') else tr = t[*,0]
       if (keyword_set(md_rem)) then begin
           cc = reform( tr )
           remove_dipole, cc, mask[*,0], nside=nside, ordering='ring' ; else remove_dipole, tr[*,0], mask[*,0], nside=nside, ordering='ring', /onlymonopole
           temp[*,ifg] = cc
       endif else begin
           temp[*,ifg] = reform( tr )
       endelse
       if (keyword_set(debug)) then mollview, temp[*,ifg], win=30, px=650, /asinh 
   endfor

   T = dblarr(nfg, nfg)
   B = dblarr(nfg)
   for ifg=0,nfg-1 do begin

; --- Changed
       print, ' - Code changed: correlate is not used anymore!'
       if False then begin
;##       B[ifg] = total( (map[gp]-mean(map[gp]))*(temp[gp,ifg]-mean(temp[gp,ifg])) / rms[gp]^2 ) / ngp
           B[ifg] = correlate( map[gp]/rms[gp], temp[gp,ifg]/rms[gp], /covariance, /double )
           for jfg=ifg,nfg-1 do begin
;##           T[ifg,jfg] = total( (temp[gp,ifg]-mean(temp[gp,ifg]))*(temp[gp,jfg]-mean(temp[gp,jfg])) / rms[gp]^2 ) / ngp
               T[ifg,jfg] = correlate( temp[gp,ifg]/rms[gp], temp[gp,jfg]/rms[gp], /covariance, /double )
               T[jfg,ifg] = T[ifg,jfg]
               oneOsigma[ifg,jfg] = 2. * T[ifg,jfg]
               oneOsigma[jfg,ifg] = oneOsigma[ifg,jfg]
           endfor
       endif

       B[ifg] = total( map[gp]/rms[gp] * temp[gp,ifg]/rms[gp] )
       for jfg=ifg,nfg-1 do begin
           T[ifg,jfg] = total( temp[gp,ifg] * temp[gp,jfg] / rms[gp]^2 )
           T[jfg,ifg] = T[ifg,jfg]

           oneOsigma[ifg,jfg] = 2. * T[ifg,jfg]
           oneOsigma[jfg,ifg] = oneOsigma[ifg,jfg]
       endfor
   endfor

   if (keyword_set(debug)) then print, T, B

   Tm1 = invert(T, /double, status)
   print, 'Inversion status: ', status
       
   A = Tm1 ## B
       
   sigma[*,*] = invert(reform(oneOsigma[*,*]),/double, status)
   print, 'Inversion status sigma: ', status
   for kfg = 0, nfg-1 do begin
       if (sigma[kfg,kfg] gt 0.) then begin
           a_err[kfg] = sqrt(sigma[kfg,kfg])
       endif else begin
           print, 'Sigma Error'
       endelse
   endfor

   print, ' ============================================================'
   print, 'Coefficients:'

   for ifg=0,nfg-1 do print, template_files[ifg], a[ifg],' +/-',a_err[ifg]

;   print, 'A:     ', A
;   print, 'A_err: ', A_err
   print, ' ============================================================'
       
   fit = map*0.
   for ifg=0,nfg-1 do fit=fit + temp[*,ifg]*A[ifg]
   cc = fit
   cc[bp] = -1.6375e30
   if (do_plot eq 'plot') then    mollview, cc, chars=1.5, tit='!17Template Linear Combination', min=-min(abs([min(fit),max(fit)])), max=min(abs([min(fit),max(fit)])), win=2, px=650

   res = map-fit
   cc = res
   cc[bp] = -1.6375e30
   if (do_plot eq 'plot') then    mollview, cc, chars=1.5, tit='!17Residual Foregrounds', min=-min(abs([min(fit),max(fit)])), max=min(abs([min(fit),max(fit)])), win=3, px=650

   print, ' --- End of Procedure ---'

return, res

end
