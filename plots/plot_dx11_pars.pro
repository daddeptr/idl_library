
   True = 1b
   False = 0b

   get_lik = True
   init = True

   if init then begin
       mollview, findgen(12l*8^2), win=3, px=450
       loadct, 39
       !p.color=0
       !p.background=255
   endif

   dops = False
   dopars = True

   dir = [ $
           '/global/scratch2/sd/hou/us_planck/for_davide/', $
;
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
;
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
;
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/' $
             ]
   files = [ $
            'dx11_hr_1-2_IQUmap_100_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2000_Tmask_x_Pmask_l4000', $
            'dx11_hr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2500_Tmask_x_Pmask_l4000', $
            'dx11_hr_1-2_IQUmap_217_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l3000_Tmask_x_Pmask_l4000', $
;
            'dx11_yr_1-2_IQUmap_100_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2000_Tmask_x_Pmask_l4000', $
            'dx11_yr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2500_Tmask_x_Pmask_l4000', $
            'dx11_yr_1-2_IQUmap_217_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l3000_Tmask_x_Pmask_l4000', $
;
            'dx11_ds_1-2_IQUmap_100_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2000_Tmask_x_Pmask_l4000', $
            'dx11_ds_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2500_Tmask_x_Pmask_l4000', $
            'dx11_ds_1-2_IQUmap_217_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l3000_Tmask_x_Pmask_l4000' $
           ]

   nfiles = n_elements(files)

   selection = strarr(nfiles+1)
   chains = strarr(nfiles+1)

   selection[0] = 'planck_data/planck_lowl_lowLike/'
   chains[0] = 'base_planck_lowl_lowLike'

;##   selection[1:nfiles1] = 'xfaster_cosmomc/runs/run_plg'+date1+'/run_plg'+date1+'_'+strtrim(string(lindgen(nfiles1)+1),2)+'/chains/'
;##   chains[1:nfiles1] = 'run_plg'+date1+'_'+strtrim(string(lindgen(nfiles1)+1),2)
;##   selection[1+nfiles1:*] = 'xfaster_cosmomc/runs/run_plg'+date2+'/run_plg'+date2+'_'+strtrim(string(lindgen(nfiles2)+1),2)+'/chains/'
;##   chains[1+nfiles1:*] = 'run_plg'+date2+'_'+strtrim(string(lindgen(nfiles2)+1),2)

   selection[1:nfiles] = files[*] + '/chains/'
   chains[1:nfiles] = 'xf_chains_TP'
;##   chains[2] = 'xf_chains_TP'
   chains[nfiles] = 'xf_chains'

;   indx = sort(files)
;
;   files = files[indx]
;   selection[1:*] = selection[indx+1]
;   chains[1:*] = chains[indx+1]

   str_select = '*insideM*'

   grun = strmatch( files, str_select )

   grun = where(grun eq True)

   win = 3

;##   tags = ['a','b']

   grun = [0,[grun]+1]
   print, grun

   selection = selection[grun]
   chains = chains[grun]

   print, chains, selection

   folders = dir + selection + chains
   print, folders

   nfold = n_elements(folders)

   run_tags = folders
   run_tags[0] = 'CamSpec (T)'

;   run_tags[1] = 'hr 143GHz';: clean-IN'
   run_tags[1] = 'hr 100GHz';: clean-IN'
   run_tags[2] = 'hr 143GHz';: clean-IN'
   run_tags[3] = 'hr 217GHz';: clean-IN'

   run_tags[4] = 'yr 100GHz';: clean-IN'
   run_tags[5] = 'yr 143GHz';: clean-IN'
   run_tags[6] = 'yr 217GHz';: clean-IN'

   run_tags[7] = 'ds 100GHz';: clean-IN'
   run_tags[8] = 'ds 143GHz';: clean-IN'
   run_tags[9] = 'ds 217GHz';: clean-IN'

   colors = [235, 0, 70, 210, 30, 80, 200, 40, 90, 170 ];, 90, 90];, 70, 100, 100, 100, 210, 210, 210] ;150, 150,180, 180, 0, 0, 120, 120]

   if dopars then begin
   nbins = 26

   npars = 11

;##   ftag = 0
;##   ltag = 1

   p_lab = strarr(npars)
   p_lab[0] = '!7X!6!db!nh!u2!n'
   p_lab[1] = '!7X!6!dc!nh!u2!n'
   p_lab[2] = '!7H!6!dA!n'
   p_lab[3] = '!6n!ds!n'
   p_lab[4] = '!6A!dns!n'
   p_lab[5] = '!6A!dsz!n'
;##   p_lab[6] = '!6A!dps!u143!n'
   p_lab[6] = '!6A!dps!n'
   p_lab[7] = '!6A!dcib!n'
   p_lab[8] = '!6H!d0!n'
   p_lab[9] = '!6A!dee!n'
   p_lab[10] = '!7s!6'

   p_key = strarr(npars)
   p_key[0] = 'omegabh2'
   p_key[1] = 'omegach2'
   p_key[2] = 'theta'
   p_key[3] = 'ns'
   p_key[4] = 'logA'
;##   p_key[5] = 'asz'
   p_key[5] = 'atsz'
;##   p_key[6] = 'aps'
   p_key[6] = 'apsTT'
   p_key[7] = 'acib'
   p_key[8] = 'H0*'
;##   p_key[9] = 'aee'
   p_key[9] = 'apsEE'
   p_key[10] = 'tau'

   p_range = fltarr(npars,2)
   p_range[0,*] = [0.0205,0.0235]
;##   p_range[0,*] = [0.0195,0.0235]
   p_range[1,*] = [0.110,0.1275]
   p_range[2,*] = [1.035,1.045]
   p_range[3,*] = [0.925,1.015]
;##   p_range[4,*] = [2.95,3.20]
   p_range[4,*] = [2.93,3.20]
   p_range[5,*] = [0,15.]
   p_range[6,*] = [0,500]
   p_range[7,*] = [0,100]
   p_range[8,*] = [62,72]
   p_range[9,*] = [0,50]
   p_range[10,*] = [0.04,0.12]

   if get_lik then begin
       likes = fltarr(nfold,npars,nbins,2)
       peaks = fltarr(nfold,npars)
       error = fltarr(nfold,npars,2)
       ival = 0
       ih = 1
       for i=0,nfold-1 do begin
;##           for itag=ftag,ltag do begin
;##               if (i ne 0) then root = folders[i] + tags[itag] else root = folders[i] + '*'
               root = folders[i]

               print, i, root
               plot_cosmomc_chains, root, par_tags=p_key, values=vals, hist=h, nbins=nbins, /norm_max

               likes[i,*,*,ival] = vals
               likes[i,*,*,ih] = h

               for ip=0,npars-1 do begin
                   print, p_key[ip]
                   mx = max(h[ip,*],im)
                   peaks[i,ip] = vals[ip,im]
                   tlik = h[ip,*]
                   print, total( tlik )
                   if total(tlik gt 0.) then begin
                   tval = vals[ip,*]
                   sur = total( (tlik[0:nbins-2]+tlik[1:nbins-1])*(tval[1:nbins-1]-tval[0:nbins-2])/2 )
                   print, sur
                   tlik = tlik / sur
                   thres = 0.
                   step=0.
                   tmx = max(tlik, tim)
;                   stop
                   while thres lt 0.32 do begin
                       step = step + 0.02
                       gp = where( (tlik lt tmx*step) and (tval lt tval[tim]) )
                       print, tmx*step
                       print, gp
;stop
                       if gp[0] ge 0 then thresl = total( (tlik[gp]+tlik[gp+1])*(tval[gp+1]-tval[gp])/2 ) else thresl = 0.
                       gp = where( (tlik lt tmx*step) and (tval ge tval[tim]) )
                       if gp[0] ge 0 then thresu = total( (tlik[gp]+tlik[gp-1])*(tval[gp]-tval[gp-1])/2 ) else thresu = 0.
                       thres = thresl+thresu
                       print, step, thres, thresl, thresu
;stop
                   endwhile
                   gp = where(tlik gt tmx * step)
                   error[i,ip,0] = tval[tim]-tval[gp[0]]
                   error[i,ip,1] = tval[gp[n_elements(gp)-1]]-tval[tim]
                   print, ' - error: ',error[i,ip,0] , error[i,ip,1] 
;                   stop
                   endif
               endfor

;##           endfor
       endfor
   



       print, ' - got likelihoods.'
; setting ps_143 for camb
       peaks[0,5] = 6.8
       error[0,5,*] = 0.

       peaks[0,6] = 63.3
       error[0,6,*] = 10.

       save, filename='dx11_likelihoods.sav', likes, peaks, error, p_range, p_lab, p_key, run_tags
       
   endif else begin

       restore, 'dx11_likelihoods.sav'

   endelse

   !p.multi = [0,3,4]
   if dops then begin
       set_plot, 'ps'
       device, file='dx11_hr_pars_summary.eps', /col, bits=8, /landscape
   endif else begin
       window, win, xsize=1300, ysize=450*1300./720
   endelse

   for ip=0,npars-1 do begin
       
       nxtic = 4
;##       xtic = findgen(nxtic+1)*(max(likes[*,*,ip,*,ival])-min(likes[*,*,ip,*,ival]))/nxtic + min(likes[*,*,ip,*,ival])
       xtic = findgen(nxtic+1)*( p_range[ip,1]-p_range[ip,1])/nxtic + p_range[ip,0]
       xticnames = [' ', 'CamSpec', '100', '143', '217', ' ']

       nytic = 4
;##       ytic = findgen(nytic+1)*(max(likes[*,*,ip,*,ih])-min(likes[*,*,ip,*,ih])) / nytic + min(likes[*,*,ip,*,ih])
       ytic = findgen(nytic+1) / nytic + min(likes[*,ip,*,ih])
       yticnames = string(findgen(nytic+1)/nytic,format='(f3.1)')
   
       if dops then plot, /nodata, [0,5], [p_range[ip,0], p_range[ip,1]], xs=1, ys=1, ytit=p_lab[ip], chars=1.5, $
         xtickname=xticnames
       if not dops then plot, /nodata, [0,5], [p_range[ip,0], p_range[ip,1]], xs=1, ys=1, ytit=p_lab[ip], chars=2.5, $
         xtickname=xticnames

       if ip eq 0 then begin
           ;legend, ['T-P','T'], line=[0,2], /top, /right
       endif
;       for i=0,nfold-1 do begin
           if ip eq 1 then begin
               ;xyouts, 0.105, max(likes[*,*,ip,*,ih])*(0.9-0.075*i),run_tags[i], col=colors[i], chars=1.25
           endif
;##           for itag=ftag,ltag do begin
;               print, ' - processing "'+folders[i]+tags[itag]+'"'
;##               oplot, likes[i,itag,ip,*,ival], likes[i,itag,ip,*,ih], thick=2, line=2-itag*2, col=colors[i]
;               if i ne 0 then begin
                   oplot, [1,1], [peaks[0,ip], peaks[0,ip]], psym=5, col=colors[0], thick=1.5
                   oplot, [1,1], peaks[0,ip]+[-error[0,ip,0], error[0,ip,1]], col=colors[0], thick=1.5
                   for i=1,nfold-1 do begin
                       oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, [peaks[i,ip],peaks[i,ip]], psym=5+1, col=colors[i], thick=1.5
                       oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, peaks[i,ip]+[-error[i,ip,0],error[i,ip,1]], col=colors[i], thick=1.5
                   endfor
;               endif else begin
;                   oplot, peaks[*,itag,ip], psym=2-itag*2, col=colors
;                   oplot, likes[i,itag,ip,*,ival], smooth( likes[i,itag,ip,*,ih], 3, /edge_mirror ), thick=2, col=colors[i]
;               endelse
;##           endfor
;       endfor
       ;oplot, [bf_pars[0,ip],bf_pars[0,ip]], [0,1]*max(likes[*,*,ip,*,ih]), thick=2;, col=230
   endfor

   plot, /nodata, [0,1], [0,1], col=255, chars=2.5, position=[0.70,0.05, 0.95,0.25]
   oplot,[0.05,0.05],[0.925,0.925], psym=6, thick=2
   xyouts, 0.1, 0.9, ': T+P'  ;, chars=1.25

   oplot,[0.5,0.5],[0.925,0.925], psym=5, thick=2
   xyouts, 0.55, 0.9, ': T';, chars=1.25
   for ixxx=0, 3 do begin
       xyouts, 0., 0.65-0.25*ixxx, run_tags[ixxx], col=colors[ixxx]
   endfor
   for ixxx=4, 6 do begin
       xyouts, 0.4, 0.65-0.25*(ixxx-3), run_tags[ixxx], col=colors[ixxx]
   endfor
   for ixxx=7, n_elements(run_tags)-1 do begin
       xyouts, 0.8, 0.65-0.25*(ixxx-6), run_tags[ixxx], col=colors[ixxx]
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



