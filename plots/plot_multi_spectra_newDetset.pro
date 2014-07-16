
   True = 1b
   False = 0b

   if False then begin
       mollview, findgen(12l*4^2), win=2
       loadct, 39
       !p.color=0
       !p.background=255
   endif

   dospec = True
   dops = True
   xlog=1

   if dospec then begin
       outdir = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/'
       files = outdir + [ $
                          'dx11_pre_ds_1-2_Imap_217_full_flat_ncl4_v1.8_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_ds_1-2_Imap_217_full_extMask_545_undusted_split_flat_ncl4_v1.8_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_ds_1-2_Imap_217_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_ds_1-2_Imap_143_full_flat_ncl4_v1.8_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_ds_1-2_Imap_143_full_extMask_545_undusted_split_flat_ncl4_v1.8_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_ds_1-2_Imap_100_full_flat_ncl4_v1.8_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_ds_1-2_Imap_100_full_extMask_545_undusted_split_flat_ncl4_v1.8_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_ds_1-2_Imap_100_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat' $
                        ]
                        
       run_tags = '!6ds1-2: '+ [ $
                                 '217 '+['raw','clean-OUT','clean-IN'], $
                                 '143 '+['raw','clean-OUT','clean-IN'], $
                                 '100 '+['raw','clean-OUT','clean-IN'] $
                               ]

       print, run_tags
       print, files
       nfiles = n_elements( files )
       lmax = 3000
       !p.multi = [0,3,2]
;       colors = 245/(nfiles-1)*lindgen(nfiles)
       colors = [0,0,0,40,40,40,95,95,95]
       psyms = [5,2,6,5,2,6,5,2,6]
       if dops then begin
           set_plot, 'ps'
           if not xlog then device, file='preDX11_newDetset_spectra.eps', /col, bits=8, /landscape else device, file='preDX11_newDetset_spectra_xlog.eps', /col, bits=8, /landscape
       endif else begin
           window, 2, xsize=1350, ysize=450*1300./720
       endelse

       readcol,'data/planck_base_lensedCls.dat', fl, ftt, fee, fbb, fte

       fst_run=3
; ------ TT
;       cl = extract_xfaster_newdata_output(files[0], lcen=l, btcl=tcl)
       plot, fl, ftt, chars=1.5, xr=[1,lmax], xs=1, ytit='!8D!dl!uTT!n [!7l!8K!u2!n]', position=[0.06,0.45,0.325,0.975], xtickname=strarr(6)+' ', xlog=xlog
;       oplot, fl, ftt, thick=2
       for i=0,nfiles-1 do begin
           cl = extract_xfaster_newdata_output(files[i], lcen=l, cler=cler, btcl=btcl)
           oplot, l, cl, col=colors[i], psym=psyms[i]
           errplot, l, cl-cler, cl+cler, col=colors[i]
;           print, ''
           if not xlog then begin
               xyouts, 750, 5500-300*i, run_tags[i], col=colors[i]
               oplot, [690,690], [5545-300*i, 5545-300*i], psym=psyms[i], col=colors[i]
           endif else begin
               xyouts, 2, 5500-300*i, run_tags[i], col=colors[i]
               oplot, [1.5,1.5], [5545-300*i, 5545-300*i], psym=psyms[i], col=colors[i]
           endelse
       endfor
;       oplot, fl, ftt, thick=2, col=245
       oplot, l, btcl, psym=6, col=245
;##       legend, run_tags, col=colors, psym=6, /top, /right, chars=1
; ------ EE
;       cl = extract_xfaster_newdata_output(files[0], lcen=l, btcl=tcl, ncl='EE')
       plot, fl, fee, chars=1.5, xr=[1,lmax], xs=1, ytit='!8D!dl!uEE!n [!7l!8K!u2!n]', position=[0.3875,0.45,0.6525,0.975], xtickname=strarr(6)+' ', yr=[0,60], xlog=xlog
;       oplot, fl, fee, thick=2
       for i=0,nfiles-1 do begin
           cl = extract_xfaster_newdata_output(files[i], lcen=l, cler=cler, ncl='EE', btcl=btcl)
           oplot, l, cl, col=colors[i], psym=psyms[i]
           errplot, l, cl-cler, cl+cler, col=colors[i]
;           print, ''
       endfor
;       oplot, fl, fee, thick=2, col=245
       oplot, l, btcl, psym=6, col=245
;       legend, run_tags, col=colors, psym=6, /top, /right, chars=1.5

; ------ TE
;       cl = extract_xfaster_newdata_output(files[0], lcen=l, btcl=tcl, ncl='TE')
       plot, fl, fte, chars=1.5, xr=[1,lmax], xs=1, ytit='!8D!dl!uTE!n [!7l!8K!u2!n]', position=[0.715,0.45,0.98,0.975], xtickname=strarr(6)+' ', xlog=xlog
;       oplot, fl, fte, thick=2
       for i=0,nfiles-1 do begin
           cl = extract_xfaster_newdata_output(files[i], lcen=l, cler=cler, ncl='TE', btcl=btcl)
           oplot, l, cl, col=colors[i], psym=psyms[i]
           errplot, l, cl-cler, cl+cler, col=colors[i]
;           print, ''
       endfor
;       oplot, fl, fte, thick=2, col=245
       oplot, l, btcl, psym=6, col=245
;       legend, run_tags, col=colors, psym=6, /top, /right, chars=1.5

; ------ Residuas
; ------ TT
       ;cl = extract_xfaster_newdata_output(files[0], lcen=l, btcl=tcl)
       plot, fl, ftt*0., chars=1.5, yr=[-750,750], xr=[1,lmax], xs=1, xtit='!8l', ytit='!6Residuals [!7l!8K!u2!n]', position=[0.06,0.075,0.325,0.44], xlog=xlog
       oplot, fl, ftt*0., thick=2
       for i=0,nfiles-1 do begin
           cl = extract_xfaster_newdata_output(files[i], lcen=l, cler=cler, btcl=tcl)
           oplot, l, cl-tcl, col=colors[i], psym=psyms[i]
           errplot, l, cl-tcl-cler, cl-tcl+cler, col=colors[i]
;           print, ''
       endfor
       oplot, fl, ftt*0, thick=2, col=245
; ------ EE
;       cl = extract_xfaster_newdata_output(files[0], lcen=l, btcl=tcl, ncl='EE')
       plot, fl, ftt*0., chars=1.5, yr=[-5,5], xr=[1,lmax], xs=1, xtit='!8l', ytit='!6Residuals [!7l!8K!u2!n]', position=[0.3875,0.075,0.6525,0.44], xlog=xlog  
       oplot, fl, ftt*0., thick=2
       for i=0,nfiles-1 do begin
           cl = extract_xfaster_newdata_output(files[i], lcen=l, cler=cler, ncl='EE', btcl=tcl)
           oplot, l, cl-tcl, col=colors[i], psym=psyms[i]
           errplot, l, cl-tcl-cler, cl-tcl+cler, col=colors[i]
;           print, ''
       endfor
       oplot, fl, ftt*0, thick=2, col=245

; ------ TE
;       cl = extract_xfaster_newdata_output(files[0], lcen=l, btcl=tcl, ncl='TE')
       plot, fl, ftt*0., chars=1.5, yr=[-20,20], xr=[1,lmax], xs=1, xtit='!8l', ytit='!6Residuals [!7l!8K!u2!n]', position=[0.715,0.075,0.98,0.44], xlog=xlog  
       oplot, fl, ftt*0., thick=2
       for i=0,nfiles-1 do begin
           cl = extract_xfaster_newdata_output(files[i], lcen=l, cler=cler, ncl='TE', btcl=tcl)
           oplot, l, cl-tcl, col=colors[i], psym=psyms[i]
           errplot, l, cl-tcl-cler, cl-tcl+cler, col=colors[i]
;           print, ''
       endfor
       oplot, fl, ftt*0, thick=2, col=245

       if dops then begin
           device, /close
           set_plot, 'x'
       endif
;       write_png, 'preDX11_detset_spectra.png', tvrd(/true)

;help, readme
;print, readme
!p.multi=0
;stop
   endif


   stop, ' --- End of Script ---'

end



