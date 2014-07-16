
   True = 1b
   False = 0b

   init = True
   do_dx11 = True
   do_pdx11 = False

   if init then begin
       mollview, findgen(12l*8^2), win=3, px=450
       loadct, 39
       !p.color=0
       !p.background=255
   endif

   dops = True
   dopars = True


   restore, '/global/scratch2/sd/dpietrob/Software/XFaster/dx11_pre_likelihoods.sav'
   pdx_likes = likes
   pdx_peaks = peaks
   pdx_error = error
   pdx_drun_tags = run_tags
   pdx_nfold = n_elements(pdx_run_tags)

   restore, '/global/scratch2/sd/dpietrob/Software/XFaster/dx11_likelihoods.sav'
   dx_likes = likes
   dx_peaks = peaks
   dx_error = error
   dx_run_tags = run_tags
   dx_nfold = n_elements(dx_run_tags)
   
   colors = [235, 0, 70, 210, 30, 80, 200, 40, 90, 170 ];, 90, 90];, 70, 100, 100, 100, 210, 210, 210] ;150, 150,180, 180, 0, 0, 120, 120]

   if dopars then begin
       nbins = 26

       npars = 11


       !p.multi = [0,3,4]
       if dops then begin
           set_plot, 'ps'
           if do_pdx11 then device, file='dx11p_pars_comparison.eps', /col, bits=8, /landscape
           if do_dx11 then device, file='dx11_pars_comparison.eps', /col, bits=8, /landscape
           if do_pdx11 and do_dx11 then device, file='dx11p-dx11_pars_comparison.eps', /col, bits=8, /landscape
       endif else begin
           window, 3, xsize=1300, ysize=450*1300./720
       endelse

       for ip=0,npars-1 do begin
       
           nxtic = 4
;##       xtic = findgen(nxtic+1)*(max(likes[*,*,ip,*,ival])-min(likes[*,*,ip,*,ival]))/nxtic + min(likes[*,*,ip,*,ival])
           xtic = findgen(nxtic+1)*( p_range[ip,1]-p_range[ip,1])/nxtic + p_range[ip,0]
           xticnames = [' ', 'CamSpec', '100', '143', '217', ' ']

           nytic = 4
;##       ytic = findgen(nytic+1)*(max(likes[*,*,ip,*,ih])-min(likes[*,*,ip,*,ih])) / nytic + min(likes[*,*,ip,*,ih])
;##           ytic = findgen(nytic+1) / nytic + min(likes[*,ip,*,ih])
;##           yticnames = string(findgen(nytic+1)/nytic,format='(f3.1)')
   
           if dops then plot, /nodata, [0,5], [p_range[ip,0], p_range[ip,1]], xs=1, ys=1, ytit=p_lab[ip], chars=1.5, $
             xtickname=xticnames
           if not dops then plot, /nodata, [0,5], [p_range[ip,0], p_range[ip,1]], xs=1, ys=1, ytit=p_lab[ip], chars=2.5, $
             xtickname=xticnames

           oplot, [1,1], [dx_peaks[0,ip], dx_peaks[0,ip]], psym=5, col=colors[0], thick=1.5
           oplot, [1,1], dx_peaks[0,ip]+[-dx_error[0,ip,0], dx_error[0,ip,1]], col=colors[0], thick=1.5

           if do_dx11 then begin
               for i=1,dx_nfold-1 do begin
                   oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, [dx_peaks[i,ip],dx_peaks[i,ip]], psym=5+1, col=colors[i], thick=1.5
                   oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, dx_peaks[i,ip]+[-dx_error[i,ip,0],dx_error[i,ip,1]], col=colors[i], thick=1.5
               endfor
           endif

           if do_pdx11 then begin
               for i=1,dx_nfold-1-3 do begin
                   oplot, 2+[(i+3-1) mod 3, (i+3-1) mod 3]+((i+3)/3)*0.1, [pdx_peaks[i,1,ip],pdx_peaks[i,1,ip]], psym=3+1, col=colors[i+3], thick=1.5
                   oplot, 2+[(i+3-1) mod 3, (i+3-1) mod 3]+((i+3)/3)*0.1, pdx_peaks[i,1,ip]+[-pdx_error[i,1,ip,0],pdx_error[i,1,ip,1]], col=colors[i+3], thick=1.5
               endfor
           endif

       endfor

   plot, /nodata, [0,1], [0,1], col=255, chars=2.5, position=[0.70,0.05, 0.95,0.25]
   oplot,[0.05,0.05],[0.925,0.925], psym=6, thick=2
   xyouts, 0.1, 0.9, ': DX11 T+P'  ;, chars=1.25

   oplot,[0.6,0.6],[0.925,0.925], psym=4, thick=2
   xyouts, 0.65, 0.9, ': preDX11 T+P';, chars=1.25
   for ixxx=0, 3 do begin
       xyouts, 0., 0.65-0.25*ixxx, dx_run_tags[ixxx], col=colors[ixxx]
   endfor
   for ixxx=4, 6 do begin
       xyouts, 0.4, 0.65-0.25*(ixxx-3), dx_run_tags[ixxx], col=colors[ixxx]
   endfor
   for ixxx=7, n_elements(run_tags)-1 do begin
       xyouts, 0.8, 0.65-0.25*(ixxx-6), dx_run_tags[ixxx], col=colors[ixxx]
   endfor

;   stop

   if dops then begin
       device, /close
       set_plot, 'x'
   endif
stop

   endif

;endfor

   stop, ' --- End of Script ---'

end



