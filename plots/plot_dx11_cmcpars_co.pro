
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

   n = 17.0 ; the circle will be "created" with 17 data points (vertices)
   theta = findgen(n)/(n-1.0)*360.0*!DtoR ; 
   x = 1.0*sin(theta)
   y = 1.0*cos(theta)
   usersym, x, y

   dops = True
   dopars = True

   dir = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/'

   files = [ $
            'planck_lowl_lowLike', $
;
            'dx11_yr_1-2_IQUmap_143_full_ns2048_uK_hrhs_xfcl_l2500_Tmask_x_Pmask_l4000', $                                                                  ;TP raw
;
            'dx11_yr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2500_Tmask_x_Pmask_l4000', $          ;                     TP ctp-bin default masks
;
            'dx11_yr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_Tmask_x_Pmask_l4000/' + $
             'dx11_yr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_const2_xfcl_l2500_Tmask_x_Pmask_l4000/', $ ;                     TP const2-bin default masks
;
            'dx11_yr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_Tmask_x_Pmask_l4000/' + $
             'dx11_yr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_const2_xfcl_l2500_Tmask_x_Pmask_l4000/', $ ;                     T  const2-bin default masks
;
            'dx11_yr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/' + $
             'dx11_yr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_const2_xfcl_l2500_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/', $ ; TP const2-bin smaller masks
;
            'dx11_yr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/' + $
             'dx11_yr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_const2_xfcl_l2500_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/', $ ; T  const2-bin smaller masks
;
;## Actually it has a smaller sky fraction            'dx11_yr_1-2_IQUmap_143_full_extMask_545_co_undusted_split_ns2048_uK_hrhs_xfcl_l2500_cons-apo_Tmask_x_Pmask_l4000', $
;            'dx11_yr_1-2_IQUmap_143_full_extMask_545_co_undusted_split_ns2048_uK_hrhs_xfcl_l2500_cons-apo_Tmask3_x_Pmask_l4000', $
;            'dx11_yr_1-2_IQUmap_143_full_extMask_545_co_undusted_split_ns2048_uK_hrhs_const2_xfcl_l3000_cons-apo_Tmask3_x_Pmask_l4000', $
;            'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l3000_cons-apo_Tmask3_x_Pmask_l4000', $
;            'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_xfcl_l3000_Tmask_x_Pmask_l4000', $      ; TP ctp-bin
;            'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l3000_Tmask_x_Pmask_l4000', $ ; TP const2-bin
            'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_Tmask_x_Pmask_l4000/'+ $
              'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l2500_Tmask_x_Pmask_l4000/', $                    ; TP const2-bin default masks
;
            'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_Tmask_x_Pmask_l4000/'+ $
              'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l2500_Tmask_x_Pmask_l4000/', $                    ; T  const2-bin default masks
;
            'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/'+ $
              'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l2500_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/', $ ; TP const2-bin smaller masks
;
            'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/'+ $
              'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l2500_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/', $ ; T  const2-bin smaller masks
;
            'dx11_yr_1-2_IQUmap_217_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/'+ $
              'dx11_yr_1-2_IQUmap_217_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/', $ ; T  const2-bin smaller masks
;
            'dx11_yr_1-2_IQUmap_100_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/'+ $
              'dx11_yr_1-2_IQUmap_100_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l2000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/' $ ; T  const2-bin smaller masks
           ]

   nfiles = n_elements(files)

   selection = strarr(nfiles)
   selection = files + '/results/'

   chains = strarr(nfiles)

   chains[0] = 'base_planck_lowl_lowLike'
   chains[1] = 'xf_chains_TP'
   chains[2] = 'xf_chains_TP'
   chains[3] = 'xf_chains_TP_const2'
   chains[4] = 'xf_chains_T_const2'
   chains[5] = 'xf_chains_TP_const2'
   chains[6] = 'xf_chains_T_const2'
   chains[7] = 'xf_chains_TP_const2'
   chains[8] = 'xf_chains_T_const2'
   chains[9] = 'xf_chains_TP_const2'
   chains[10] = 'xf_chains_T_const2'
   chains[11] = 'xf_chains_TP_const2'
   chains[12] = 'xf_chains_T_const2'

   folders = dir + selection + chains
   print, folders

   nfold = n_elements(folders)

   symbs = [5,6,6,6,5,6,5,6,5,6,5,6,6]
   run_tags = folders
   run_tags[0] = 'CamSpec (T)'
   run_tags[1] = 'ctp raw DefM TP';: clean-IN'
   run_tags[2] = 'ctp IN DefM TP';: clean-IN'
   run_tags[3] = 'const2 IN DefM';: clean-IN'
;   run_tags[4] = 'const2 DefM T';: clean-IN'
   run_tags[4] = 'const2 IN SmallM';: clean-IN'
;   run_tags[6] = 'const2 SmallM T';: clean-IN'
;## Smaller sky fraction   run_tags[3] = 'dx-co apoT: 143GHz';: clean-IN'
   run_tags[5] = 'co const2 DefM';: clean-IN'
;   run_tags[8] = 'co const2 DefM T';: clean-IN'
   run_tags[6] = 'co const2 SmallM';: clean-IN'
;
   run_tags[7] = 'co const2 SmallM';: clean-IN'
;   run_tags[10] = 'co const2 SmallM T';: clean-IN'
;##   run_tags[10] = 'dx-co const2 SmallMask T';: clean-IN'

   colors = [235, 40, 70, 90, 90, 210, 210, 200, 200, 120, 120, 180, 80 ];, 90, 90];, 70, 100, 100, 100, 210, 210, 210] ;150, 150,180, 180, 0, 0, 120, 120]

   if dopars then begin

   npars = 11

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
   p_key[8] = 'H0'
;##   p_key[9] = 'aee'
   p_key[9] = 'apsEE'
   p_key[10] = 'tau'

   p_range = fltarr(npars,2)
   p_range[0,*] = [0.02125,0.02305]
;##   p_range[0,*] = [0.0195,0.0235]
   p_range[1,*] = [0.110,0.1275]
   p_range[2,*] = [1.035,1.045]
   p_range[3,*] = [0.925,1.015]
;##   p_range[4,*] = [2.95,3.20]
   p_range[4,*] = [2.93,3.20]
   p_range[5,*] = [0,15.]
   p_range[6,*] = [0,150]
   p_range[7,*] = [0,100]
   p_range[8,*] = [62,72]
   p_range[9,*] = [0,50]
   p_range[10,*] = [0.04,0.12]

   if get_lik then begin
       peaks = fltarr(nfold,npars)
       error = fltarr(nfold,npars,2)
       for i=0,nfold-1 do begin
               root = folders[i]

               print, i, ': '+root
               readcol, root+'.margestats', keys, p, l, u, format='a,f,x,f,f', skipline=3
               
;               if i gt 7 then begin
;                   p_key[5]='asz'
;                   p_key[6]='aps'
;                   p_key[7]='acl'
;                   p_key[9]='aee'
;               endif
               for ip=0,npars-1 do begin
                   indx = where(keys eq p_key[ip])
                   peaks[i,ip] = p[indx]
                   error[i,ip,0] = p[indx]-l[indx]
                   error[i,ip,1] = u[indx]-p[indx]
                   print, p_key[ip], peaks[i,ip], ' - error: ',error[i,ip,0] , error[i,ip,1] 
               endfor
;               if i gt 7 then begin
;                   peaks[i,9] *= 3.*3001/1001.
;                   error[i,9,*] *= 3.*3001/1001.
;               endif
       endfor
   
       print, ' - got likelihoods.'
; setting ps_143 for camb
       peaks[0,5] = 6.8
       error[0,5,*] = 0.

       peaks[0,6] = 63.3
       error[0,6,*] = 10.

   endif

   !p.multi = [0,3,4]
   if dops then begin
       set_plot, 'ps'
       device, file='dx11_yr-co_pars.eps', /col, bits=8, /landscape
   endif else begin
       window, 3, xsize=1300, ysize=450*1300./720
   endelse

   for ip=0,npars-1 do begin
       
;##       nxtic = 4
;##       xtic = findgen(nxtic+1)*( p_range[ip,1]-p_range[ip,1])/nxtic + p_range[ip,0]
;##       xticnames = [' ', 'CamSpec', '100', '143', '217', ' ']
;##       xticnames = [' ', 'CS', '143', '143', '143', '143', '143', '143', '143', '143', '143', '143', ' ', ' ']
       xticnames = strarr(nfold+2) +' '

       if dops then plot, /nodata, [0,nfold+1], [p_range[ip,0], p_range[ip,1]], xs=1, ys=1, ytit=p_lab[ip], chars=1.5, $
         xtickname=xticnames
       if not dops then plot, /nodata, [0,nfold+1], [p_range[ip,0], p_range[ip,1]], xs=1, ys=1, ytit=p_lab[ip], chars=2.5, $
         xtickname=xticnames

       oplot, [1,1], [peaks[0,ip], peaks[0,ip]], psym=5, col=colors[0], thick=2
       oplot, findgen(nfold+2), peaks[0,ip]+fltarr(nfold+2), col=colors[0], line=2
       oplot, [1,1], peaks[0,ip]+[-error[0,ip,0], error[0,ip,1]], col=colors[0], thick=2

       for i=1,nfold-1 do begin
;##           oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, [peaks[i,ip],peaks[i,ip]], psym=5+1, col=colors[i], thick=2
;##           oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, peaks[i,ip]+[-error[i,ip,0],error[i,ip,1]], col=colors[i], thick=2
           oplot, [i+1,i+1], [peaks[i,ip],peaks[i,ip]], psym=symbs[i], col=colors[i], thick=2
           oplot, [i+1,i+1], peaks[i,ip]+[-error[i,ip,0],error[i,ip,1]], col=colors[i], thick=2
       endfor

   endfor

   plot, /nodata, [0,1], [0,1], col=255, chars=2.5, position=[0.68,0.05, 0.95,0.25]
   oplot,[0.05,0.05],[0.925,0.925], psym=6, thick=2
   xyouts, 0.1, 0.9, ': T+P'  ;, chars=1.25

   oplot,[0.5,0.5],[0.925,0.925], psym=5, thick=2
   xyouts, 0.55, 0.9, ': T';, chars=1.25
   xyouts, 0., 0.65, run_tags[0], col=colors[0]
   xyouts, 0.5, 0.65, 'DX11 143 GHz Year maps'
;##   xyouts, 0.5, 0.65, 'DX11'
   for ixxx=1,3 do begin
       xyouts, 0., 0.65-0.25*ixxx, run_tags[ixxx], col=colors[ixxx]
   endfor
   for ixxx=4,6 do begin
       xyouts, 0.65, 0.65-0.25*(ixxx-3), run_tags[ixxx], col=colors[3+(ixxx-3)*2]
   endfor
;   for ixxx=7, n_elements(run_tags)-1 do begin
;       xyouts, 0.8, 0.65-0.25*(ixxx-6), run_tags[ixxx], col=colors[ixxx]
;   endfor

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



