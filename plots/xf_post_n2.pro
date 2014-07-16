;## /global/scratch2/sd/hou/us_planck/for_davide/xfaster_cosmomc/runs/run_plg131024.readme


   True = 1b
   False = 0b

   get_lik = False
   init = False
   if init then begin
       mollview, findgen(12l*4^2), win=2
       loadct, 39
       !p.color=0
       !p.background=255
   endif

   dospec = False
   dopars = True

   dir = '/global/scratch2/sd/hou/us_planck/for_davide/'
;   dir = dir + 'run_plg131024/'

   tags = ['a','b']
   selection = [ $
                 'planck_data/planck_lowl_lowLike/',$
                 'xfaster_cosmomc/runs/run_plg131024/run_plg131024_1/chains/',$
                 'xfaster_cosmomc/runs/run_plg131024/run_plg131024_2/chains/',$
                 'xfaster_cosmomc/runs/run_plg131024/run_plg131024_3/chains/',$
                 'xfaster_cosmomc/runs/run_plg131024/run_plg131024_4/chains/'$
               ]
   chains = [ $
              'base_planck_lowl_lowLike', $
              'run_plg131024_1', $
              'run_plg131024_2', $
              'run_plg131024_3', $
              'run_plg131024_4' $
            ]

   run_tags = ['CAMSPEC (T)', 'Dset-IN extMask', 'Full-IN', 'Year-IN', 'Full-IN extMask']

   folders = dir + selection + chains
   print, folders
   nfold = n_elements(folders)

;##   colors = findgen( nfold) / (nfold-1) * 245
   colors = [213, 70, 100, 40, 245]

   tvlct, 255, 170,   0, 10  ;; camspec
   tvlct, 220,   0,   0, 20  ;; TT
   tvlct,   0,   0, 220, 30  ;; EE
   tvlct,   0, 220,   0, 40  ;; TE
   tvlct, 200,   0, 255, 50  ;; TT+EE
   tvlct,   0, 255, 255, 60  ;; EE+TE

   if dopars then begin
   nbins = 26

   npars = 11

   bf_pars = [ [0.02205,0.00028], $
               [0.1199,0.0027], $
               [1.04131,0.00063], $
               [0.9609,0.0073], $
               [3.089,0.0255], $
               [0,0], $
               [54.,10.], $
               [0,0], $
               [67.3,1.2], $
               [0,0], $
               [0.089,0.002] ]

   p_lab = strarr(npars)
   p_lab[0] = '!7X!6!db!nh!u2!n'
   p_lab[1] = '!7X!6!dc!nh!u2!n'
   p_lab[2] = '!7H!6!dA!n'
   p_lab[3] = '!6n!ds!n'
   p_lab[4] = '!6A!dns!n'
   p_lab[5] = '!6A!dsz!n'
   p_lab[6] = '!6A!dps!u143!n'
   p_lab[7] = '!6A!dcl!n'
   p_lab[8] = '!6H!d0!n'
   p_lab[9] = '!6A!dee!n'
   p_lab[10] = '!7s!6'

   p_key = strarr(npars)
   p_key[0] = 'omegabh2'
   p_key[1] = 'omegach2'
   p_key[2] = 'theta'
   p_key[3] = 'ns'
   p_key[4] = 'logA'
   p_key[5] = 'asz'
   p_key[6] = 'aps'
   p_key[7] = 'acl'
   p_key[8] = 'H0'
   p_key[9] = 'aee'
   p_key[10] = 'tau'

   if get_lik then begin
   likes = fltarr(nfold,2,npars,nbins,2)
   ival = 0
   ih = 1
   for i=0,nfold-1 do begin
       for itag=0,1 do begin
           if (i ne 0) then root = folders[i] + tags[itag] else root = folders[i] + '*'

;           if False then begin
               readcol, root+'.paramnames', par_labs, format='a'
;           print, par_labs
               ipars = []
               for ix=0,n_elements(p_key)-1 do ipars= [ipars,where(par_labs eq p_key[ix])]
               ipars = ipars[where(ipars ge 0)]
               print, ipars
;stop
;           endif
               plot_cosmomc_chains, root, ipars=3+ipars, values=vals, hist=h, nbins=nbins, /norm_max
;           if itag eq 0 then plot_cosmomc_chains, root, ipars=[3,4,5,7,8,9,10,11,17], values=vals, hist=h, nbins=nbins;, /verbose
;           if itag eq 1 then plot_cosmomc_chains, root, ipars=[3,4,5,7,8,9,10,11,18], values=vals, hist=h, nbins=nbins;, /verbose
               help, vals
               likes[i,itag,0:n_elements(ipars)-1,*,ival] = vals
               likes[i,itag,0:n_elements(ipars)-1,*,ih] = h
;stop
           endfor
       endfor
   endif

;##   !p.multi = [0,3,ceil(npars/3)]
   !p.multi = [0,3,4]
;   window, 1, xsize=1300, ysize=450*1300./720
   set_plot, 'ps'
   device, file='mylike.eps', /col, bits=8, /landscape

   for ip=0,npars-1 do begin
       
       nxtic = 3
       xtic = findgen(nxtic+1)*(max(likes[*,*,ip,*,ival])-min(likes[*,*,ip,*,ival]))/nxtic + min(likes[*,*,ip,*,ival])

       nytic = 4
       ytic = findgen(nytic+1)*(max(likes[*,*,ip,*,ih])-min(likes[*,*,ip,*,ih])) / nytic + min(likes[*,*,ip,*,ih])
       yticnames = string(findgen(nytic+1)/nytic,format='(f3.1)')
   
       plot, /nodata, [min(likes[*,*,ip,*,ival]), max(likes[*,*,ip,*,ival])], [min(likes[*,*,ip,*,ih]), max(likes[*,*,ip,*,ih])], xtit=p_lab[ip], ytit='P', $
         xtickv=xtic, xticks=nxtic, yticks=nytic, ytickv=ytic, ytickname=yticnames,chars=1.5
       if ip eq 0 then begin
           ;legend, ['T-P','T'], line=[0,2], /top, /right
       endif
       for i=0,nfold-1 do begin
           if ip eq 1 then begin
               ;xyouts, 0.105, max(likes[*,*,ip,*,ih])*(0.9-0.075*i),run_tags[i], col=colors[i], chars=1.25
           endif
           for itag=0,1 do begin
               print, ' - processing "'+folders[i]+tags[itag]+'"'
;##               oplot, likes[i,itag,ip,*,ival], likes[i,itag,ip,*,ih], thick=2, line=2-itag*2, col=colors[i]
               oplot, likes[i,itag,ip,*,ival], smooth( likes[i,itag,ip,*,ih], 3, /edge_mirror ), thick=2, line=2-itag*2, col=colors[i]
           endfor
       endfor
       ;oplot, [bf_pars[0,ip],bf_pars[0,ip]], [0,1]*max(likes[*,*,ip,*,ih]), thick=2;, col=230
   endfor

   plot, /nodata, [0,1], [0,1], col=255, chars=2.5
   oplot,[0,0.10],[0.95,0.95], col=0, thick=2
   xyouts, 0.25, 0.95, '(T+P)';, chars=1.25
   oplot,[0.5,0.6],[0.95,0.95], col=0, line=2, thick=2
   xyouts, 0.78, 0.95, '(T)';, chars=1.25
   for ixxx=0, n_elements(run_tags)-1 do begin
       oplot ,[0,0.1],[0.75-0.15*ixxx,0.75-0.15*ixxx], col=colors[ixxx], thick=2
       xyouts, 0.25, 0.75-0.15*ixxx, run_tags[ixxx], col=colors[ixxx];, chars=1.25
   endfor
;   stop

   device, /close
   set_plot, 'x'
stop
   endif

   if dospec then begin
       readcol, '/global/scratch2/sd/hou/us_planck/for_davide/xfaster_cosmomc/runs/run_plg131024.readme', readme, format='a'
       ixf = []
       files = []
       for i=0,n_elements( readme)-1 do begin
           line = readme[i]
           field = strsplit(line,'.', /extract)
           if field[n_elements(field)-1] eq 'newdat' then begin
               ixf = [ixf,i]
               field = strsplit(line,'/',/extract)
               files = [files, field[n_elements(field)-1]]
           endif
       endfor
;       print, readme[ixf]
       print, files
       lmax = 2000
       !p.multi = [0,3,2]
       window, 2, xsize=1300, ysize=450*1300./720

       fst_run=3
; ------ TT
       cl = extract_xfaster_newdata_output('outputs/'+files[0], lcen=l, btcl=tcl)
       plot, l, tcl, chars=3, xr=[1,lmax], xs=1, xtit='!6l', ytit='!6D!dl!uTT!n [!7l!6K!u2!n]'  
       oplot, l, tcl, thick=2
       for i=0,n_elements(ixf)-1 do begin
           cl = extract_xfaster_newdata_output('outputs/'+files[i], lcen=l, cler=cler)
           oplot, l, cl, col=colors[i], psym=6
           errplot, l, cl-cler, cl+cler, col=colors[i]
;           print, ''
       endfor
       legend, run_tags, col=colors, psym=6, /top, /right, chars=1
; ------ EE
       cl = extract_xfaster_newdata_output('outputs/'+files[0], lcen=l, btcl=tcl, ncl='EE')
       plot, l, tcl, chars=3, xr=[1,lmax], xs=1, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]'  
       oplot, l, tcl, thick=2
       for i=0,n_elements(ixf)-1 do begin
           cl = extract_xfaster_newdata_output('outputs/'+files[i], lcen=l, cler=cler, ncl='EE')
           oplot, l, cl, col=colors[i], psym=6
           errplot, l, cl-cler, cl+cler, col=colors[i]
;           print, ''
       endfor
;       legend, run_tags, col=colors, psym=6, /top, /right, chars=1.5

; ------ TE
       cl = extract_xfaster_newdata_output('outputs/'+files[0], lcen=l, btcl=tcl, ncl='TE')
       plot, l, tcl, chars=3, xr=[1,lmax], xs=1, xtit='!6l', ytit='!6D!dl!uTE!n [!7l!6K!u2!n]'  
       oplot, l, tcl, thick=2
       for i=0,n_elements(ixf)-1 do begin
           cl = extract_xfaster_newdata_output('outputs/'+files[i], lcen=l, cler=cler, ncl='TE')
           oplot, l, cl, col=colors[i], psym=6
           errplot, l, cl-cler, cl+cler, col=colors[i]
;           print, ''
       endfor
;       legend, run_tags, col=colors, psym=6, /top, /right, chars=1.5

; ------ Residuas
; ------ TT
       cl = extract_xfaster_newdata_output('outputs/'+files[0], lcen=l, btcl=tcl)
       plot, l, tcl*0., chars=3, yr=[-150,150], xr=[1,lmax], xs=1, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]'  
       oplot, l, tcl*0., thick=2
       for i=0,n_elements(ixf)-1 do begin
           cl = extract_xfaster_newdata_output('outputs/'+files[i], lcen=l, cler=cler)
           oplot, l, cl-tcl, col=colors[i], psym=6
           errplot, l, cl-tcl-cler, cl-tcl+cler, col=colors[i]
;           print, ''
       endfor
; ------ EE
       cl = extract_xfaster_newdata_output('outputs/'+files[0], lcen=l, btcl=tcl, ncl='EE')
       plot, l, tcl*0., chars=3, yr=[-4,4], xr=[1,lmax], xs=1, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]'  
       oplot, l, tcl*0., thick=2
       for i=0,n_elements(ixf)-1 do begin
           cl = extract_xfaster_newdata_output('outputs/'+files[i], lcen=l, cler=cler, ncl='EE')
           oplot, l, cl-tcl, col=colors[i], psym=6
           errplot, l, cl-tcl-cler, cl-tcl+cler, col=colors[i]
;           print, ''
       endfor

; ------ TE
       cl = extract_xfaster_newdata_output('outputs/'+files[0], lcen=l, btcl=tcl, ncl='TE')
       plot, l, tcl*0., chars=3, yr=[-20,20], xr=[1,lmax], xs=1, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]'  
       oplot, l, tcl*0., thick=2
       for i=0,n_elements(ixf)-1 do begin
           cl = extract_xfaster_newdata_output('outputs/'+files[i], lcen=l, cler=cler, ncl='TE')
           oplot, l, cl-tcl, col=colors[i], psym=6
           errplot, l, cl-tcl-cler, cl-tcl+cler, col=colors[i]
;           print, ''
       endfor



;help, readme
;print, readme
!p.multi=0
stop
   endif


   stop, ' --- End of Script ---'

end



