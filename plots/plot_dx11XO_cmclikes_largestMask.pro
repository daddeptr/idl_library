
   True = 1b
   False = 0b

   get_lik = False
   init = True

   dops = False
   dopars = True

   if init then begin
       mollview, findgen(12l*16^2), win=13, px=450
       loadct, 39
       !p.color=0
       !p.background=255
   endif

   n = 17.0 ; the circle will be "created" with 17 data points (vertices)
   theta = findgen(n)/(n-1.0)*360.0*!DtoR ; 
   x = 1.0*sin(theta)
   y = 1.0*cos(theta)
   usersym, x, y

   dir = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/'

   files = [ $
            'planck_lowl_lowLike', $
;
            'ddx9_check', $
;
            'dx11_XO_co-clean_maps', $
;
            'dx11_XO_co-clean_maps', $  
;
            'dx11_XO_co-clean_maps', $  
;
            'dx11_XO_co-clean_maps', $  
;
            'dx11_XO_co-clean_maps', $
;
            'dx11_XO_co-clean_maps', $
;
            'dx11_XO_co-clean_maps', $
;
            'dx11_XO_co-clean_maps', $  
;
            'dx11_XO_co-clean_maps', $  
;;
            'dx11_XO_co-clean_maps', $  
;;
            'predx11_bpc_smica_cmb_hrhs_05a_2048_IQU', $
;
            'predx11_bpc_smica_cmb_hrhs_05a_2048_IQU' $
              ]

   nfiles = n_elements(files)

   selection = strarr(nfiles)
   selection = files + '/results/plot_data/'

;##   chains = strarr(nfiles)

   chains = [ $
              'base_planck_lowl_lowLike', $
              'base_planck_tauprior', $
              'xf_chains_100_T_const2_largestMask_l1000', $
              'xf_chains_143_T_const2_largestMask_l1000', $
              'xf_chains_217_T_const2_largestMask_l1000', $
              'xf_chains_100_T_const2_largestMask_l1500', $
              'xf_chains_143_T_const2_largestMask_l1500', $
              'xf_chains_217_T_const2_largestMask_l1500', $
              'xf_chains_100_T_const2_largestMask', $
              'xf_chains_143_T_const2_largestMask', $
              'xf_chains_217_T_const2_largestMask', $
              'xf_chains_217_T_const2_largestMask_l3000', $
              'xf_chains_T_const2', $
              'xf_chains_TP_const2' $
               ]

   run_tags = [ $
                'CamSpec', $
                'CamSpec (tau-prior)', $
                'XO-100 (T) l1000', $
                'XO-143 (T) l1000', $
                'XO-217 (T) l1000', $
                'XO-100 (T) l1500', $
                'XO-143 (T) l1500', $
                'XO-217 (T) l1500', $
                'XO-100 (T)', $
                'XO-143 (T)', $
                'XO-217 (T)', $
                'XO-217 (T) l3000', $
                'predx11-Smica (T)', $
                'predx11-Smica (TP)' $
              ]

   folders = dir + selection + chains
   print, folders


   symbs = [5,5,5,5,5,5,5,5,5,5,5,5,6]
   ;##run_tags = [folders,'','']
;##   run_tags[*] = ''

   colors = [245, 235, 70, 110, 205, 70, 110, 205, 70, 110, 205, 205, 5, 5 ]

   indx = lindgen(n_elements(folders)-2)
   folders = folders[indx]
   run_tags = run_tags[indx]
   colors = colors[indx]
   symbs = symbs[indx]

   nfold = n_elements(folders)

   if dopars then begin

       npars = 10

       p_lab = strarr(npars)
       p_lab[0] = '!7X!6!db!nh!u2!n'
       p_lab[1] = '!7X!6!dc!nh!u2!n'
       p_lab[2] = '!7H!6!dA!n'
       p_lab[3] = '!6n!ds!n'
;##       p_lab[4] = '!6A!dns!n'
       p_lab[5] = '!6A!dsz!n'
;##   p_lab[6] = '!6A!dps!u143!n'
       p_lab[6] = '!6A!dps!n'
       p_lab[7] = '!6A!dcib!n'
       p_lab[8] = '!6H!d0!n'
       p_lab[9] = '!6A!dee!n'
;##       p_lab[10] = '!7s!6'
;       p_lab[10] = '!6A!dgal!uTT!n'
       p_lab[4] = '!610!u9!nA!ds!nexp(-2!7s!6)'
       
       p_key = strarr(npars)
       p_key[0] = 'omegabh2'
       p_key[1] = 'omegach2'
       p_key[2] = 'theta'
       p_key[3] = 'ns'
;##       p_key[4] = 'logA'
;##   p_key[5] = 'asz'
       p_key[5] = 'atsz'
;##   p_key[6] = 'aps'
       p_key[6] = 'apsTT'
       p_key[7] = 'acib'
       p_key[8] = 'H0'
;##   p_key[9] = 'aee'
       p_key[9] = 'apsEE'
;##       p_key[10] = 'tau'
;       p_key[10] = 'agalTT'
       p_key[4] = 'clamp'

       p_range = fltarr(npars,2)
       p_range[0,*] = [0.020,0.025]
;##   p_range[0,*] = [0.0195,0.0235]
       p_range[1,*] = [0.09,0.15]
       p_range[2,*] = [1.035,1.05]
       p_range[3,*] = [0.9,1.2]
;##   p_range[4,*] = [2.95,3.20]
;##       p_range[4,*] = [3.00,3.15]
       p_range[5,*] = [0,20.]
       p_range[6,*] = [0,500]
       p_range[7,*] = [0,100]
       p_range[8,*] = [60,82]
       p_range[9,*] = [0,60]
;##       p_range[10,*] = [0.04,0.12]
;       p_range[10,*] = [0, 5]
       p_range[4,*] = [1.75,2.]
       
   endif

   !p.multi = [0,3,4]
   if dops then begin
       set_plot, 'ps'
       device, file='dx11_XO_co_pars.eps', /col, bits=8, /landscape
   endif else begin
       window, 13, xsize=1300, ysize=450*1300./720
   endelse

   for ip=0,npars-1 do begin
       
       xticnames = strarr(nfold+2) +' '

       if dops then plot, /nodata, [p_range[ip,0], p_range[ip,1]], [0.,1.01], xs=1, ys=1, ytit=p_lab[ip], chars=1.5
       if not dops then plot, /nodata, [p_range[ip,0], p_range[ip,1]], [0.,1.01], xs=1, ys=1, ytit=p_lab[ip], chars=2.5

;       oplot, [1,1], [peaks[0,ip], peaks[0,ip]], psym=symbs[0], col=colors[0], thick=2
;       oplot, findgen(nfold+2), peaks[1,ip]+fltarr(nfold+2), line=2
;       oplot, [1,1], peaks[0,ip]+[-error[0,ip,0], error[0,ip,1]], col=colors[0], thick=2

       for i=0,nfold-1 do begin
           root = folders[i]
           spawn, 'ls '+root+'_p_'+p_key[ip]+'.dat', files
           if strlen(files) gt 0 then begin
               readcol, root+'_p_'+p_key[ip]+'.dat', x, y, /silent
;##           oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, [peaks[i,ip],peaks[i,ip]], psym=5+1, col=colors[i], thick=2
;##           oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, peaks[i,ip]+[-error[i,ip,0],error[i,ip,1]], col=colors[i], thick=2
               lin=4
               if i ge 2 then lin = 3 - ( (i-2)/3 +1 )
               oplot, x, y, col=colors[i], thick=2, line= lin
;           oplot, [i+1,i+1], peaks[i,ip]+[-error[i,ip,0],error[i,ip,1]], col=colors[i], thick=2
           endif else begin
               print, ' - skip: '+root+'_p_'+p_key[ip]+'.dat'
           endelse
       endfor
       
   endfor

   plot, /nodata, [0,1], [0,1], col=255, chars=2.5, position=[0.38,0.05, 0.95,0.25]
   oplot,[0.05,0.05],[0.925,0.925], psym=6, thick=2
   xyouts, 0.1, 0.9, ': T+P'  ;, chars=1.25

   oplot,[0.5,0.5],[0.925,0.925], psym=5, thick=2
   xyouts, 0.55, 0.9, ': T';, chars=1.25

   xyouts, 0., 0.65, run_tags[0], col=colors[0], chars=0.65
;##   xyouts, 0.5, 0.65, 'New TP-Masks3'
;##   xyouts, 0.5, 0.65, 'DX11'
   for ixxx=1,3 do begin
       xyouts, 0., 0.65-0.25*(ixxx), run_tags[ixxx], col=colors[ixxx], chars=0.65
   endfor
   for ixxx=4,7 do begin
       xyouts, 0.2, 0.9-0.25*(ixxx-3), run_tags[ixxx], col=colors[ixxx], chars=0.65
   endfor
   for ixxx=8,11 do begin
       xyouts, 0.4, 0.9-0.25*(ixxx-7), run_tags[ixxx], col=colors[ixxx], chars=0.65
   endfor
   for ixxx=12,15 do begin
       xyouts, 0.63, 0.9-0.25*(ixxx-11), run_tags[ixxx], col=colors[ixxx], chars=0.65
   endfor
   for ixxx=16,nfold-1 do begin
       xyouts, 0.85, 0.9-0.25*(ixxx-15), run_tags[ixxx], col=colors[ixxx], chars=0.65
   endfor

   if dops then begin
       device, /close
       set_plot, 'x'
   endif
stop


   stop, ' --- End of Script ---'

end



