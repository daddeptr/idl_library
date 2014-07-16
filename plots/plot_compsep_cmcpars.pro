
   True = 1b
   False = 0b

   get_lik = True
   init = True

   dops = False
   dopars = True

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

   dir = [ $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/CompSep/outputs/dx11_pre/ns1024', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/CompSep/outputs/dx11_pre/ns1024', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/CompSep/outputs/dx11_pre/ns1024', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/CompSep/outputs/dx11_pre/ns1024', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/CompSep/outputs/dx11_pre/ns1024' $
   ]

   files = [ $
             'planck_lowl_lowLike', $
             '','','','','' $
           ]

   nfiles = n_elements(files)

   selection = strarr(nfiles)
   selection = files + '/results/'

   methods = ['commander','nilc','sevem', 'smica','smicaI_smicaP']

   chains = strarr(nfiles)
   chains[0] = 'base_planck_lowl_lowLike'
   chains[1:4] = 'xf_chains_'+methods[0:3]+'_TP_const2'
   chains[5] = 'xf_chains_'+methods[4]+'_const2'

   folders = dir + selection + chains
   print, folders

   nfold = n_elements(folders)

   symbs = lonarr(nfold) + 5
   symbs[*] = 4
   run_tags = [folders ]
   run_tags = [ $
                'CamSpec', $
                'commander', $
                'nilc', $
                'sevem', $
                'smica', $
                'hybrid smica' $
              ]

   colors = [240, 70,90,105,115,135,155,169,200,210,220, 0,40,60,90,105,115,135,155,169,200,210,220 ]


   if dopars then begin

       npars = 10

       p_lab = strarr(npars)
       p_lab[0] = '!7X!6!db!nh!u2!n'
       p_lab[1] = '!7X!6!dc!nh!u2!n'
       p_lab[2] = '!7H!6!dA!n'
       p_lab[3] = '!6n!ds!n'
;##       p_lab[4] = '!6A!dns!n'
       p_lab[4] = '!610!u10!nA!ds!nexp(-2!7s!6)'
       p_lab[5] = '!6A!dsz!n'
;##   p_lab[6] = '!6A!dps!u143!n'
       p_lab[6] = '!6A!dps!n'
       p_lab[7] = '!6A!dcib!n'
       p_lab[8] = '!6H!d0!n'
       p_lab[9] = '!6A!dee!n'
;##       p_lab[10] = '!7s!6'
       
       p_key = strarr(npars)
       p_key[0] = 'omegabh2'
       p_key[1] = 'omegach2'
       p_key[2] = 'theta'
       p_key[3] = 'ns'
;##       p_key[4] = 'logA'
       p_key[4] = 'clamp'
;##   p_key[5] = 'asz'
       p_key[5] = 'atsz'
;##   p_key[6] = 'aps'
       p_key[6] = 'apsTT'
       p_key[7] = 'acib'
       p_key[8] = 'H0'
;##   p_key[9] = 'aee'
       p_key[9] = 'apsEE'
;##       p_key[10] = 'tau'

       p_mod = strarr(npars)
       p_mod[0] = 0.02205
       p_mod[1] = 0.1199
       p_mod[2] = 1.040325
       p_mod[3] = 0.9619
       p_mod[4] = 3.086
       p_mod[5] = 0.
       p_mod[6] = 0.
       p_mod[7] = 0.
       p_mod[8] = 67.3
       p_mod[9] = 0.
;##       p_mod[10] = 0.092

       p_range = fltarr(npars,2)
       p_range[0,*] = [0.02125,0.02305]
;##   p_range[0,*] = [0.0195,0.0235]
       p_range[1,*] = [0.110,0.1275]
       p_range[2,*] = [1.035,1.045]
       p_range[3,*] = [0.925,1.015]
;##   p_range[4,*] = [2.95,3.20]
;##       p_range[4,*] = [3.00,3.15]
       p_range[4,*] = [1.80,1.90]
       p_range[5,*] = [0,15.]
       p_range[6,*] = [0,150]
       p_range[7,*] = [0,50]
       p_range[8,*] = [62,72]
       p_range[9,*] = [0,75]
;##       p_range[10,*] = [0.04,0.12]

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
   
       print, ' - got likelihoods.'

   endif

   !p.multi = [0,3,4]
   if dops then begin
       set_plot, 'ps'
       device, file='dx11_143-yr-co_pars_vs_sims.eps', /col, bits=8, /landscape
   endif else begin
       window, 3, xsize=1300, ysize=450*1300./720
   endelse

   for ip=0,npars-1 do begin
       
       xticnames = strarr(nfold+2) +' '

       if dops then plot, /nodata, [0,nfold+2], [p_range[ip,0], p_range[ip,1]], xs=1, ys=1, ytit=p_lab[ip], chars=1.5, $
         xtickname=xticnames
       if not dops then plot, /nodata, [0,nfold+2], [p_range[ip,0], p_range[ip,1]], xs=1, ys=1, ytit=p_lab[ip], chars=2.5, $
         xtickname=xticnames

       oplot, [1,1], [peaks[0,ip], peaks[0,ip]], psym=symbs[0], col=colors[0], thick=2

       oplot, findgen(nfold+4), p_mod[ip]+fltarr(nfold+4), col=235, line=2
       if ip eq 4 then begin
           oplot, findgen(nfold+4), 0.994*p_mod[ip]+fltarr(nfold+4), col=25, line=3
           xyouts, 1, 3.02, 'scaling=0.994', col=25
       endif
       if ip eq 3 then begin
           oplot, findgen(nfold+4), 1.007*p_mod[ip]+fltarr(nfold+4), col=25, line=3
           xyouts, 1, 0.935, 'delta=0.007', col=25
       endif

       oplot, [1,1], peaks[0,ip]+[-error[0,ip,0], error[0,ip,1]], col=colors[0], thick=2

       for i=1,nfold-1 do begin
;##           oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, [peaks[i,ip],peaks[i,ip]], psym=5+1, col=colors[i], thick=2
;##           oplot, 2+[(i-1) mod 3, (i-1) mod 3]+(i/3)*0.1, peaks[i,ip]+[-error[i,ip,0],error[i,ip,1]], col=colors[i], thick=2
           oplot, [i+1,i+1], [peaks[i,ip],peaks[i,ip]], psym=symbs[i], col=colors[i], thick=2
           oplot, [i+1,i+1], peaks[i,ip]+[-error[i,ip,0],error[i,ip,1]], col=colors[i], thick=2
       endfor
       
   endfor

   plot, /nodata, [0,1], [0,1], col=255, chars=2.5, position=[0.68,0.05, 0.95,0.25]
   oplot,[0.05,0.05],[0.925,0.925], psym=4, thick=2
   xyouts, 0.1, 0.9, ': T egfg=0'  ;, chars=1.25

   oplot,[0.5,0.5],[0.925,0.925], psym=5, thick=2
   xyouts, 0.55, 0.9, ': T';, chars=1.25
   xyouts, 0., 0.65, run_tags[0], col=colors[0]
;##   xyouts, 0.5, 0.65, 'New TP-Masks3'
;##   xyouts, 0.5, 0.65, 'DX11'
   for ixxx=1,3 do begin
       xyouts, 0., 0.65-0.25*ixxx, run_tags[ixxx], col=colors[ixxx]
   endfor
   for ixxx=4,7 do begin
       xyouts, 0.4, 0.65-0.25*(ixxx-4), run_tags[ixxx], col=colors[ixxx]
   endfor
   for ixxx=8,nfold-1 do begin
       xyouts, 0.75, 0.65-0.25*(ixxx-8), run_tags[ixxx], col=colors[ixxx]
   endfor
;   for ixxx=7, n_elements(run_tags)-1 do begin
;       xyouts, 0.8, 0.65-0.25*(ixxx-8), run_tags[ixxx], col=colors[ixxx]
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



