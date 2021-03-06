;## /global/scratch2/sd/hou/us_planck/for_davide/xfaster_cosmomc/runs/run_plg131024.readme


   dir = '/global/scratch2/sd/hou/us_planck/for_davide/xfaster_cosmomc/runs/run_plg131024/'

   True = 1b
   False = 0b

   if False then begin
       mollview, findgen(12l*4^2), win=2
       loadct, 39
       !p.color=0
       !p.background=255
;       loadct, 74
   endif

   dospec = True

   if dospec then begin
       outdir = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/'
       files = outdir + [ $
                          'dx11_pre_year_1-2_Imap_143_full_flat_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat', $
                          'dx11_pre_year_1-2_Imap_143_full_flat_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_dx11_extended_QUmask_ps100-353_ns2048_l3000.newdat', $
                          'dx11_pre_year_1-2_Imap_143_full_extMask_undusted_insideM_split_flat_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat', $
                          'dx11_pre_year_1-2_Imap_143_full_extMask_undusted_insideM_split_flat_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_dx11_extended_QUmask_ps100-353_ns2048_l3000.newdat', $
                          'dx11_pre_year_1-2_Imap_143_full_extMask_undusted_split_flat_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_dx11_extended_QUmask_ps100-353_ns2048_l3000.newdat' $
                        ]
                        
       run_tags = '!6yr1-2: '+['raw', 'raw extMask','clean-in','clean-in extMask', 'clean-out extMask']

       print, files
       nfiles = n_elements( files )
       lmax = 2000
       !p.multi = 0
;       colors = 245/(nfiles-1)*lindgen(nfiles)
       colors = [0,70,95,40,245,210]
       set_plot,'ps'
       device, file='preDX11_year_spectra_lowell.eps', /col, bits=8, /landscape
;       window, 2, xsize=1350, ysize=450

       readcol,'data/planck_base_lensedCls.dat', fl, ftt, fee, fbb, fte

; ------ TT
;       cl = extract_xfaster_newdata_output(files[0], lcen=l, btcl=tcl)
;       plot, fl, fee, chars=3, xr=[1,lmax], xs=1, ytit='!8D!dl!uEE!n [!7l!8K!u2!n]', position=[0.06,0.1,0.325,0.975], xtit='!8l'
;       oplot, fl, fee, thick=2
;       for i=0,nfiles-1 do begin
;           cl = extract_xfaster_newdata_output(files[i], lcen=l, cler=cler, ncl='EE')
;           oplot, l, cl, col=colors[i], psym=6
;           errplot, l, cl-cler, cl+cler, col=colors[i]
;           print, ''
;           xyouts, 350, 5500-300*i, run_tags[i], col=colors[i], chars=1.25
;       endfor
;       oplot, fl, fee, thick=2
;##       legend, run_tags, col=colors, psym=6, /top, /right, chars=1
; ------ EE
;       cl = extract_xfaster_newdata_output(files[0], lcen=l, btcl=tcl, ncl='EE')
       plot, fl, fee, chars=1.5, xr=[1,250], xs=1, ytit='!8D!dl!u!6EE!n [!7l!8K!u2!n]', xtit='!8l'
       oplot, fl, fee, thick=2
       for i=0,nfiles-1 do begin
           cl = extract_xfaster_newdata_output(files[i], lcen=l, cler=cler, ncl='EE')
           oplot, l, cl, col=colors[i], psym=6
           errplot, l, cl-cler, cl+cler, col=colors[i]
;           print, ''
       endfor
       oplot, fl, fee, thick=2
;       legend, run_tags, col=colors, psym=6, /top, /right, chars=1.5

; ------ TE
;       cl = extract_xfaster_newdata_output(files[0], lcen=l, btcl=tcl, ncl='TE')
;       plot, fl, fee, chars=3, xr=[1,50], xs=1, ytit='!8D!dl!uEE!n [!7l!8K!u2!n]', position=[0.715,0.1,0.98,0.975], xtit='!8l'
;       oplot, fl, fee, thick=2
;       for i=0,nfiles-1 do begin
;           cl = extract_xfaster_newdata_output(files[i], lcen=l, cler=cler, ncl='EE')
;           oplot, l, cl, col=colors[i], psym=6
;           errplot, l, cl-cler, cl+cler, col=colors[i]
;           print, ''
;       endfor
;       oplot, fl, fee, thick=2
;       legend, run_tags, col=colors, psym=6, /top, /right, chars=1.5
       device, /close
       set_plot, 'x'
;       write_png, 'preDX11_year_spectra_lowell.png', tvrd(/true)

;help, readme
;print, readme
!p.multi=0
stop
   endif


   stop, ' --- End of Script ---'

end



