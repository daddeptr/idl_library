
   True = 1b
   False = 0b

   if True then begin
       mollview, findgen(12l*4^2), win=2, px=550
       loadct, 39
       !p.color=0
       !p.background=255
   endif

   docomp = False
   dospec = True
   dops = False
   xlog=0

   if docomp then begin
       files = outdir + [ $
                          'dx11_pre_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
;                          'dx11_pre_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx9-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx9_ds_1-2_Imap_143_nominal_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx9-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx9_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_ffp7-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat' $
;                          'dx9_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx9-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat' $
                        ]
;       print, files
; ------ dDX9f vs dDX9n
       compare_xfaster_spectra, files[2], files[1], otit='!6143 undested-IN: !7D!6DX9 Full - !7D!6DX9 Nominal', leg_tags=['(!7D!6DX9, FFP7beam)','(!7D!6DX9 Nom, !7D!6DX9beam)'], uleg_pos=[1300,5750], bleg_pos=[1400,-60], win=1
       compare_xfaster_spectra, files[2], files[1], otit='!6143 undested-IN: !7D!6DX9 Full - !7D!6DX9 Nominal', leg_tags=['(!7D!6DX9, FFP7beam, preDX11beam)','(!7D!6DX9 Nom, !7D!6DX9beam)'], uleg_pos=[3000,5750], bleg_pos=[3400,-60], /xlog, win=2

; ------ pDX11f vs dDX9n
       compare_xfaster_spectra, files[0], files[1], otit='!6143 undested-IN: preDX11 - !7D!6DX9 Nominal', leg_tags=['(preDX11, preDX11beam)','(!7D!6DX9 Nom, !7D!6DX9beam)'], uleg_pos=[1300,5750], bleg_pos=[1400,-60], win=1
       compare_xfaster_spectra, files[0], files[1], otit='!6143 undested-IN: preDX11 - !7D!6DX9 Nominal', leg_tags=['(preDX11, preDX11beam)','(!7D!6DX9 Nom, !7D!6DX9beam)'], uleg_pos=[3000,5750], bleg_pos=[3400,-60], /xlog, win=2
      
; ------ pDX11f vs dDX9n
       compare_xfaster_spectra, files[0], files[2], otit='!6143 undested-IN: preDX11 - !7D!6DX9 Full', leg_tags=['(preDX11, preDX11beam)','(!7D!6DX9 Ful, FFP7beam)'], uleg_pos=[1300,5750], bleg_pos=[1400,-60], win=1
       compare_xfaster_spectra, files[0], files[2], otit='!6143 undested-IN: preDX11 - !7D!6DX9 Full', leg_tags=['(preDX11, preDX11beam)','(!7D!6DX9 Ful, FFP7beam)'], uleg_pos=[3000,5750], bleg_pos=[3400,-60], /xlog, win=2
      
stop
endif

   if dospec then begin
       outdir = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/'
       files = outdir + [ $
;                          'dx11_pre_ds_1-2_Imap_143_full_extMask_545_undusted_split_flat_ncl4_v1.8_dx11-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
;                          'dx11_pre_ds_1-2_Imap_143_full_extMask_545_undusted_split_flat_ncl4_v1.8_dx9-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx9-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
;
;                          'dx9_ds_1-2_Imap_143_nominal_extMask_545_undusted_split_flat_ncl4_v1.8_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx9_ds_1-2_Imap_143_nominal_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx9-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
;
;                          'dx9_ds_1-2_Imap_143_full_extMask_545_undusted_split_flat_ncl4_v1.8_dx9-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx9_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_ffp7-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
;                          'dx9_ds_1-2_Imap_143_full_extMask_545_undusted_split_flat_ncl4_v1.8_ffp7-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx9_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx9-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat' $
                        ]
                        
       run_tags =  '!6' + [ $
;                                 '217 '+['raw','clean-OUT','clean-IN'], $
                                 'dx11_pre Full dx11-beam '+['clean-IN'], $
                                 'dx11_pre Full dx9-beam '+['clean-IN'], $
                                 'dx9 Nominal dx9-beam '+['clean-IN'], $
                                 'dx9 Full ffp7-beam '+['clean-IN'], $
                                 'dx9 Full dx9-beam '+['clean-IN'] $
;                                 '100 '+['raw','clean-OUT','clean-IN'] $
                               ]

       print, run_tags
       print, files
       nfiles = n_elements( files )
       lmax = 2500
;       colors = 245/(nfiles-1)*lindgen(nfiles)
       colors = [ $
                  0, $
                  40, $
                  70, $
                  210, $
                  230 $
                ]
       psyms = [ $
                 1, $
                 2, $
                 4, $
                 5, $
                 6 $
                 ]
       if dops then begin
           set_plot, 'ps'
           device, file='preDX11_vs_DX9_detset_ratio.eps', /col, bits=8, /landscape;, xsize=13500, ysize=5000
           !p.multi = [0,3,1]
       endif else begin
           !p.multi = [0,3,1]
           window, 2, xsize=1000, ysize=500
       endelse

       readcol, 'data/beams/dx9/dx9_wl_143_ds_1-2.dat', bl, wld9
       readcol, 'data/beams/dx11_pre/dx11_pre_wl_143_ds_1-2.dat', bl, wld11
       readcol, 'data/beams/ffp7/ffp7_wl_143_ds_1-2.dat', bl, wlf7

       cltag = ['TT', 'EE', 'TE']

       for i=0,0 do begin

           cl1 = extract_xfaster_newdata_output( files[0], lcen=l, cler=cler, btcl=btcl, ncl=cltag[i] )
           cl2 = extract_xfaster_newdata_output( files[1], lcen=l, cler=cler, btcl=btcl, ncl=cltag[i] )

           cl3 = extract_xfaster_newdata_output( files[2], lcen=l, cler=cler, btcl=btcl, ncl=cltag[i] )

           cl4 = extract_xfaster_newdata_output( files[3], lcen=l, cler=cler, btcl=btcl, ncl=cltag[i] )
           cl5 = extract_xfaster_newdata_output( files[4], lcen=l, cler=cler, btcl=btcl, ncl=cltag[i] )

; ------ TT
;       cl = extract_xfaster_newdata_output(files[0], lcen=l, btcl=tcl)
;##           if dops then plot, /nodata, [0,lmax],[0.95,1.05], chars=1.5, xs=1, ys=1, xtit='!8l',ytit='!8D!dl!6!u'+cltag[i]+'!n / D!dl!udx9!n', tit=cltag[i], position=[0.32*i+0.075,0.15,(i+1)+0.32,0.75] else $
;##             plot, /nodata, [0,lmax],[0.95,1.05], chars=3, xs=1, ys=1, xtit='!8l',ytit='!8D!dl!6!u'+cltag[i]+'!n / D!dl!udx9!n',tit=cltag[i], position=[0.32*i+0.075,0.15,(i+1)*0.32,0.75]
;## Dropping EE and TE since the beam was the same
           if dops then plot, /nodata, [0,lmax],[0.95,1.05], chars=1.5, xs=1, ys=1, xtit='!8l',ytit='!8D!dl!6!n Ratio', tit=cltag[i], position=[0.075,0.15,0.64,0.85] else $
             plot, /nodata, [0,lmax],[0.95,1.05], chars=3, xs=1, ys=1, xtit='!8l',ytit='!8D!dl!6!n Ratio',tit=cltag[i], position=[0.075,0.15,0.64,0.85]
           oplot, l, wld11*0.+1., line=1;, col=245
           oplot, l, cl1/cl2
           oplot, l, cl5/cl3, col=150
           oplot, l, cl4/cl5, col=70
           oplot, l, cl1/cl3, col=230
           oplot, l, cl1/cl4, col=210
;##           oplot, l, cl1/cl5, col=0
           oplot, bl, (wld9/wld11)^2, psym=4, ns=25
           oplot, bl, (wld9/wld11)^2, psym=4, ns=20, col=230
           oplot, bl, (wld9/wlf7)^2,  psym=4, ns=20, col=70
           oplot, bl, (wlf7/wld11)^2, psym=4, ns=20, col=210
       endfor

       plot, /nodata, [0,1], [0,1], col=255, position=[0.03,0.9,0.98,0.995]
       xyouts, 0.15, 0.75, '------ 143 undusted-IN: detset1 - detset2 ------', chars=1.5

       oplot, [0.05,0.2], [0.35,0.35], line=0
       xyouts, 0.25, 0.3, '!8D!dl!n!6 ratio;', chars=1.25
;##       oplot, [0.4,0.55], [0.35,0.35], psym=4
       oplot, 0.4+findgen(11)*0.15/10, 0.35+fltarr(11), psym=4
       xyouts, 0.6, 0.3, '!8W!dl!u2!n!6 ratio.', chars=1.25

;##       oplot, [0.015,0.015],[0.4,0.4], psym=psyms[0], col=colors[0]
;##       xyouts, 0.022, 0.35, run_tags[0], col=colors[0], chars=1.

       plot, /nodata, [0,1], [0,1], col=255, position=[0.65,0.15,0.99,0.85]
       xyouts, 0.05, 0.8, '(dx11, dx11-beam)/(dx11, dx9-beam)', chars=1.
       xyouts, 0.05, 0.7, '(dx11, dx11-beam)/(dx9-nom, dx9-beam)', chars=1., col=230
       xyouts, 0.05, 0.6, '(dx9-full, ffp7-beam)/(dx9-full, dx9-beam)', chars=1., col=70
       xyouts, 0.05, 0.5, '(dx11-full, dx11-beam)/(dx9-full, ffp7-beam)', chars=1., col=210
       xyouts, 0.05, 0.4, '(dx9-full, dx9-beam)/(dx9-nom, dx9-beam)', chars=1., col=150

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



