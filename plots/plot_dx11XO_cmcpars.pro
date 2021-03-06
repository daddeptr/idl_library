
   True = 1b
   False = 0b

   get_lik = True
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
;
            'dx11_XO_co-clean_maps', $
;
            'dx11_XO_co-clean_maps', $
;
            'dx11_XO_co-clean_maps', $
;
            'predx11_bpc_smica_cmb_hrhs_05a_2048_IQU', $
;
            'predx11_bpc_smica_cmb_hrhs_05a_2048_IQU' $
              ]

   nfiles = n_elements(files)

   selection = strarr(nfiles)
   selection = files + '/results/'

;##   chains = strarr(nfiles)

   chains = [ $
              'base_planck_lowl_lowLike', $
              'base_planck_tauprior', $
              'xf_chains_100_T_const2', $
              'xf_chains_100_TP_const2', $
              'xf_chains_143_T_const2', $
              'xf_chains_143_TP_const2', $
              'xf_chains_217_T_const2', $
              'xf_chains_217_TP_const2', $
              'xf_chains_100_T_const2_oldMask', $
              'xf_chains_100_TP_const2_oldMask', $
              'xf_chains_143_T_const2_oldMask', $
              'xf_chains_143_TP_const2_oldMask', $
              'xf_chains_217_T_const2_oldMask', $$
              'xf_chains_217_TP_const2_oldMask', $
              'xf_chains_100_T_const2_oldMask_raw', $
              'xf_chains_143_T_const2_oldMask_raw', $
              'xf_chains_217_T_const2_oldMask_raw', $
              'xf_chains_143_T_const2_mPix', $
              'xf_chains_100_T_const2_oldMask_raw_lmax1500', $
              'xf_chains_100_TP_const2_oldMask_raw_lmax1500', $
              'xf_chains_100_T_const2_master', $
              'xf_chains_143_T_const2_master', $
              'xf_chains_217_T_const2_master', $
              'xf_chains_100_T_const2_Xmaster', $
              'xf_chains_143_T_const2_Xmaster', $
              'xf_chains_217_T_const2_Xmaster', $
              'xf_chains_100_T_const2_Xmaster_dustmpl', $
              'xf_chains_143_T_const2_Xmaster_dustmpl', $
              'xf_chains_217_T_const2_Xmaster_dustmpl', $
              'xf_chains_100_T_const2_Xmaster_dust', $
              'xf_chains_143_T_const2_Xmaster_dust', $
              'xf_chains_217_T_const2_Xmaster_dust', $
              'xf_chains_T_const2', $
              'xf_chains_TP_const2' $
              ]

   run_tags = [ $
                'CamSpec', $
    'CamSpec (tau-prior)', $
    'XO-100 (T)', $
    'XO-100 (TP)', $
    'XO-143 (T)', $
    'XO-143 (TP)', $
    'XO-217 (T)', $
    'XO-217 (TP)', $
    'XO-100 (T) oM', $
    'XO-100 (TP) oM', $
    'XO-143 (T) oM', $
    'XO-143 (TP) oM', $
    'XO-217 (T) oM', $
    'XO-217 (TP) oM', $
    'XO-100 (T) oM raw', $
    'XO-143 (T) oM raw', $
    'XO-217 (T) oM raw', $
    'XO-143 (T) mPix', $
    'XO-100 (T) oM-1500', $
    'XO-100 (TP) oM-1500', $
    'XO-100 (T) master', $
    'XO-143 (T) master', $
    'XO-217 (T) master', $
    'XO-100 (T) Xmaster', $
    'XO-143 (T) Xmaster', $
    'XO-217 (T) Xmaster', $
    'XO-100 (T) Xmaster dustmpl', $
    'XO-143 (T) Xmaster dustmpl', $
    'XO-217 (T) Xmaster dustmpl', $
    'XO-100 (T) Xmaster dust', $
    'XO-143 (T) Xmaster dust', $
    'XO-217 (T) Xmaster dust', $
    'predx11-Smica (T)', $
    'predx11-Smica (TP)' $
              ]

   folders = dir + selection + chains
   print, folders


   symbs = [5,5,5,6,5,6,5,6,5,6,5,6,5,6,5,5,5,5,5,6,5,5,5,5,5,5,5,5,5,5,5,5,5,6]
   ;##run_tags = [folders,'','']
;##   run_tags[*] = ''

   colors = [245, 235, 70, 70, 110, 110, 205, 205, 70, 70, 110, 110, 205, 205, 80, 140, 200, 110, 70, 70, 70, 110, 205, 70,110,205,70,110,205,70,110,205,5, 5 ]


   indx = [0,1,2,4,6,20,21,22,23,24,25,26,27,28,29,30,31,32,33]
   indx = [0,1,2,4,6,20,21,22,26,27,28]
;##   indx = [0,1,2,3,4,5,6,7]

   folders = folders[indx]
   run_tags = run_tags[indx]
   colors = colors[indx]
   symbs = symbs[indx]

   nfold = n_elements(folders)

   if dopars then begin

       npars = 11

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
       p_lab[10] = '!6A!dgal!uTT!n'
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
       p_key[10] = 'agalTT'
       p_key[4] = 'clamp'

       p_range = fltarr(npars,2)
       p_range[0,*] = [0.02125,0.02305]
;##   p_range[0,*] = [0.0195,0.0235]
       p_range[1,*] = [0.110,0.1275]
       p_range[2,*] = [1.035,1.045]
       p_range[3,*] = [0.925,1.015]
;##   p_range[4,*] = [2.95,3.20]
;##       p_range[4,*] = [3.00,3.15]
       p_range[5,*] = [0,15.]
       p_range[6,*] = [0,150]
       p_range[7,*] = [0,75]
       p_range[8,*] = [62,72]
       p_range[9,*] = [0,60]
;##       p_range[10,*] = [0.04,0.12]
       p_range[10,*] = [0, 5]
       p_range[4,*] = [1.75,1.9]
       
       if get_lik then begin
           peaks = fltarr( nfold,npars )
           error = fltarr( nfold,npars,2 )
           for i=0,nfold-1 do begin
               root = folders[i]

               print, i, ': '+root
;##               readcol, root+'.margestats', keys, p, l, u, format='a,f,x,f,f', skipline=3
               readcol, root+'.margestats', keys, p, s, format='a,f,f', skipline=3
               l = p-s
               u = p+s
               for ip=0,npars-1 do begin
                   indx = where(keys eq p_key[ip])
                   if indx[0] ge 0 then begin
                       peaks[i,ip] = p[indx]
                       error[i,ip,0] = p[indx]-l[indx]
                       error[i,ip,1] = u[indx]-p[indx]
                       print, p_key[ip], peaks[i,ip], ' - error: ',error[i,ip,0] , error[i,ip,1] 
                   endif else begin
                       print, ' - key not found!'
                   endelse
               endfor
       endfor
   
       print, ' - got likelihoods.'
; setting ps_143 for camb
;       peaks[0,5] = 6.8
;       error[0,5,*] = 0.

;       peaks[0,6] = 63.3
;       error[0,6,*] = 10.

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

       if dops then plot, /nodata, [0,nfold+2], [p_range[ip,0], p_range[ip,1]], xs=1, ys=1, ytit=p_lab[ip], chars=1.5, $
         xtickname=xticnames
       if not dops then plot, /nodata, [0,nfold+2], [p_range[ip,0], p_range[ip,1]], xs=1, ys=1, ytit=p_lab[ip], chars=2.5, $
         xtickname=xticnames

       oplot, [1,1], [peaks[0,ip], peaks[0,ip]], psym=symbs[0], col=colors[0], thick=2
       oplot, findgen(nfold+2), peaks[1,ip]+fltarr(nfold+2), line=2
       oplot, [1,1], peaks[0,ip]+[-error[0,ip,0], error[0,ip,1]], col=colors[0], thick=2

       for i=1,nfold-1 do begin
;##           oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, [peaks[i,ip],peaks[i,ip]], psym=5+1, col=colors[i], thick=2
;##           oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, peaks[i,ip]+[-error[i,ip,0],error[i,ip,1]], col=colors[i], thick=2
           oplot, [i+1,i+1], [peaks[i,ip],peaks[i,ip]], psym=symbs[i], col=colors[i], thick=2
           oplot, [i+1,i+1], peaks[i,ip]+[-error[i,ip,0],error[i,ip,1]], col=colors[i], thick=2
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

   endif

;endfor

   stop, ' --- End of Script ---'

end



