
   True = 1b
   False = 0b

   get_lik = True
   init = True

   if init then begin
       mollview, findgen(12l*16^2), win=3, px=450
       loadct, 39
       !p.color=0
       !p.background=255
   endif

   n = 17.0 ; the circle will be "created" with 17 data points (vertices)
   theta = findgen(n)/(n-1.0)*360.0*!DtoR ; 
   x = 1.0*sin(theta)
   y = 1.0*cos(theta)
   usersym, x, y

   dops = False
   dopars = True

   dir = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/'

   files = [ $
            'planck_lowl_lowLike', $
;
            'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/'+ $
              'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l2500_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/', $ ; TP const2-bin smaller masks
;
            'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/'+ $
              'dx11_yr_1-2_IQUmap_143_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l2500_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/', $ ; T  const2-bin smaller masks
;
            'smica_Imap_2013_hrhs_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/'+ $
              'smica_Imap_2013_hrhs_const2_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/' $ ; T  const2-bin smaller masks
;
;##            'dx11_yr_1-2_IQUmap_217_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/'+ $
;##              'dx11_yr_1-2_IQUmap_217_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/', $ ; T  const2-bin smaller masks
;
;##            'dx11_yr_1-2_IQUmap_100_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/'+ $
;##              'dx11_yr_1-2_IQUmap_100_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l2000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/' $ ; T  const2-bin smaller masks
           ]

   nfiles = n_elements(files)

   selection = strarr(nfiles)
   selection = files + '/results/'

   chains = strarr(nfiles)

   chains[0] = 'base_planck_lowl_lowLike'
   chains[1] = 'xf_chains_TP_const2'
   chains[2] = 'xf_chains_T_const2'
   chains[3] = 'xf_chains_T_const2'

   sim_chainsTP = 'xf_chains_TP_const2'
   sim_chainsT = 'xf_chains_T_const2'

   folders = dir + selection + chains
   print, folders

   nfold = n_elements(folders)

   symbs = [5,6,5,5,6,5,6,5,6,6]
   run_tags = [folders,'','']
   run_tags[0] = 'CamSpec (T)'
   run_tags[1] = '143 (TP)'
   run_tags[2] = '143 (T)'
   run_tags[3] = 'Smica2013 (T)'
   run_tags[4] = 'FFP7 SN sims (TP)'
   run_tags[5] = 'FFP7 SN sims (T)'

   colors = [235, 90, 90, 70, 210, 200, 200, 120, 120, 180, 80 ];, 90, 90];, 70, 100, 100, 100, 210, 210, 210] ;150, 150,180, 180, 0, 0, 120, 120]


   if dopars then begin

       npars = 11
       nsims = 99
       sim_bfTP = fltarr(nsims, npars) - 1.6375e30
       sim_bfT = fltarr(nsims, npars) - 1.6375e30

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
       p_range[4,*] = [3.00,3.15]
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
               
               for ip=0,npars-1 do begin
                   indx = where(keys eq p_key[ip])
                   peaks[i,ip] = p[indx]
                   error[i,ip,0] = p[indx]-l[indx]
                   error[i,ip,1] = u[indx]-p[indx]
                   print, p_key[ip], peaks[i,ip], ' - error: ',error[i,ip,0] , error[i,ip,1] 
               endfor
       endfor
   
       for isim=0,nsims-1 do begin
           root = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/ffp7_sn/xf_v1.9/cosmomc/results_'+strtrim(string(isim+1),2)+'/'+sim_chainsTP
;
           print, isim, ': '+root
           spawn, 'ls '+root+'.margestats', files
           nfiles = n_elements(files)
;##               if nfiles[0] ne -1 then begin
           if strlen(files) gt 0 then begin
               readcol, root+'.margestats', keys, p, l, u, format='a,f,x,f,f', skipline=3
;
               for ip=0,npars-1 do begin
                   indx = where(keys eq p_key[ip])
;##                   peaks[i,ip] = p[indx]
                   sim_bfTP[isim,ip] = p[indx]
;##                   error[i,ip,0] = p[indx]-l[indx]
;##                   error[i,ip,1] = u[indx]-p[indx]
                   print, p_key[ip], sim_bfTP[isim,ip]
               endfor
           endif else begin
               print, ' - file not found'
           endelse

               root = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/ffp7_sn/xf_v1.9/cosmomc/results_'+strtrim(string(isim+1),2)+'/'+sim_chainsT
               print, isim, ': '+root
               spawn, 'ls '+root+'.margestats', files
               nfiles = n_elements(files)
               if strlen(files) gt 0 then begin
                   readcol, root+'.margestats', keys, p, l, u, format='a,f,x,f,f', skipline=3
                   for ip=0,npars-1 do begin
                       indx = where(keys eq p_key[ip])
                       sim_bfT[isim,ip] = p[indx]
                       print, p_key[ip], sim_bfT[isim,ip]
                   endfor
               endif else begin
                   print, ' - file not found'
               endelse
           endfor

       print, ' - got likelihoods.'
; setting ps_143 for camb
       peaks[0,5] = 6.8
       error[0,5,*] = 0.

       peaks[0,6] = 63.3
       error[0,6,*] = 10.

       mbfTP = fltarr(npars,2)
       gp=where(sim_bfTP[*,0] ne -1.6375e30) 
       for ip=0,npars-1 do mbfTP[ip,0:1] = [ mean(sim_bfTP[gp,ip]), stddev(sim_bfTP[gp,ip]) ]

       mbfT = fltarr(npars,2)
       gp = where(sim_bfT[*,0] ne -1.6375e30) 
       for ip=0,npars-1 do mbfT[ip,0:1] = [ mean(sim_bfT[gp,ip]), stddev(sim_bfT[gp,ip]) ]

   endif

   !p.multi = [0,3,4]
   if dops then begin
       set_plot, 'ps'
       device, file='dx11_143-yr-co_pars_vs_sims.eps', /col, bits=8, /landscape
   endif else begin
       window, 3, xsize=1300, ysize=450*1300./720
   endelse

   for ip=0,npars-1 do begin
       
       xticnames = strarr(nfold+4) +' '

       if dops then plot, /nodata, [0,nfold+4], [p_range[ip,0], p_range[ip,1]], xs=1, ys=1, ytit=p_lab[ip], chars=1.5, $
         xtickname=xticnames
       if not dops then plot, /nodata, [0,nfold+4], [p_range[ip,0], p_range[ip,1]], xs=1, ys=1, ytit=p_lab[ip], chars=2.5, $
         xtickname=xticnames

       oplot, [1,1], [peaks[0,ip], peaks[0,ip]], psym=5, col=colors[0], thick=2
       oplot, findgen(nfold+4), peaks[0,ip]+fltarr(nfold+4), col=colors[0], line=2
       oplot, [1,1], peaks[0,ip]+[-error[0,ip,0], error[0,ip,1]], col=colors[0], thick=2

       for i=1,nfold-1 do begin
;##           oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, [peaks[i,ip],peaks[i,ip]], psym=5+1, col=colors[i], thick=2
;##           oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, peaks[i,ip]+[-error[i,ip,0],error[i,ip,1]], col=colors[i], thick=2
           oplot, [i+1,i+1], [peaks[i,ip],peaks[i,ip]], psym=symbs[i], col=colors[i], thick=2
           oplot, [i+1,i+1], peaks[i,ip]+[-error[i,ip,0],error[i,ip,1]], col=colors[i], thick=2
       endfor
       
       oplot, (i+1)+fltarr(nsims), sim_bfTP[*,ip], col=colors[i], psym=3
       oplot, [i+1,i+1], mbfTP[ip,0]+[-mbfTP[ip,1], mbfTP[ip,1]], col=colors[i], thick=2

       oplot, (i+2)+fltarr(nsims), sim_bfT[*,ip], col=colors[i+1], psym=3
       oplot, [i+2,i+2], mbfT[ip,0]+[-mbfT[ip,1], mbfT[ip,1]], col=colors[i+1], thick=2

   endfor

   plot, /nodata, [0,1], [0,1], col=255, chars=2.5, position=[0.68,0.05, 0.95,0.25]
   oplot,[0.05,0.05],[0.925,0.925], psym=6, thick=2
   xyouts, 0.1, 0.9, ': T+P'  ;, chars=1.25

   oplot,[0.5,0.5],[0.925,0.925], psym=5, thick=2
   xyouts, 0.55, 0.9, ': T';, chars=1.25
   xyouts, 0., 0.65, run_tags[0], col=colors[0]
   xyouts, 0.5, 0.65, 'New TP-Masks3'
;##   xyouts, 0.5, 0.65, 'DX11'
   for ixxx=1,3 do begin
       xyouts, 0., 0.65-0.25*ixxx, run_tags[ixxx], col=colors[ixxx]
   endfor
   for ixxx=4,5 do begin
       xyouts, 0.65, 0.65-0.25*(ixxx-3), run_tags[ixxx], col=colors[ixxx]
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



