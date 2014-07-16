
   True = 1b
   False = 0b

   if True then begin
       mollview, findgen(12l*4^2), win=2
       loadct, 39
       !p.color=0
       !p.background=255
   endif

   dospec = True
   dops   = False
   xlog   = False

   if dospec then begin
       outdir = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/'
       files = outdir + [ $
                          'dx11_pre_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_cls_Tmask_X_Pmask.fits', $
                          'dx9_ds_1-2_Imap_143_nominal_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_cls_Tmask_X_Pmask.fits', $
                          'dx9_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_cls_Tmask_X_Pmask.fits' $
                        ]
                        
       noisefiles = outdir + [ $
                          'dx11_pre_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhd_cls_Tmask_X_Pmask.fits', $
                          'dx9_ds_1-2_Imap_143_nominal_extMask_545_undusted_insideM_split_ns2048_uK_hrhd_cls_Tmask_X_Pmask.fits', $
                          'dx9_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhd_cls_Tmask_X_Pmask.fits' $
                        ]
                        
       run_tags =  '!6' + [ $
                                 'dx11_pre Full '+['clean-IN'], $
                                 'dx9 Nominal '+['clean-IN'], $
                                 'dx9 Full '+['clean-IN'] $
                               ]

       print, run_tags
       print, files
;stop
       nfiles = n_elements( files )
       lmax = 3000
;       colors = 245/(nfiles-1)*lindgen(nfiles)
       colors = [ $
                  0, $
                  40, $
                  70, $
                  100, $
                  210 $
                ]
       psyms = [ $
                 5, $
                 5, $
                 5, $
                 5, $
                 5 $
                 ]

; ------ First plot
       if dops then begin
           set_plot, 'ps'
           device, file='beamStudy_preDX11_vs_DX9_detset_pCl-ratio.eps', /col, bits=8, /landscape;, xsize=13500, ysize=5000
           !p.multi = [0,4,1]
       endif else begin
;           !p.multi = [0,4,1]
;           window, 2, xsize=1350, ysize=500
           !p.multi = [0,3,1]
           window, 3, xsize=720.*500/420, ysize=500
       endelse

       readcol, 'data/beams/dx9/dx9_wl_143_ds_1-2.dat', bl, wld9
       readcol, 'data/beams/dx11_pre/dx11_pre_wl_143_ds_1-2.dat', bl, wld11
       readcol, 'data/beams/ffp7/ffp7_wl_143_ds_1-2.dat', bl, wlf7

       cltag = ['TT', 'EE', 'BB', 'TE']

       l = findgen(5001)
       ll = l*(l+1)/2./!pi

       ii=[0];,1,3]
       for i=0,n_elements(ii)-1 do begin

           fits2cl, cl1, files[0]
           fits2cl, cl2, files[1]
           fits2cl, cl3, files[2]

           fits2cl, nl1, noisefiles[0]
           fits2cl, nl2, noisefiles[1]
           fits2cl, nl3, noisefiles[2]

; ------ TT
;       cl = extract_xfaster_newdata_output(files[0], lcen=l, btcl=tcl)
           if dops then plot, /nodata, [0,lmax],[0.4,2.1], chars=2, xs=1, ys=1, xtit='!8l',ytit='!8D!dl!6!u'+cltag[i]+'!n ratio', tit=cltag[ii[i]], position=[0.48*i+0.075,0.15,(i+1)*0.48,0.75] else $
;             plot, /nodata, [0,lmax],[0.4,2.1], chars=3, xs=1, ys=1, xtit='!8l',ytit='!8D!dl!6!u'+cltag[i]+'!n / D!dl!udx9!n',tit=cltag[ii[i]], position=[0.32*i+0.075,0.15,(i+1)*0.32,0.75]
             plot, /nodata, [0,lmax],[0.4,2.1], chars=2, xs=1, ys=1, xtit='!8l',ytit='!8pC!dl!6!u'+cltag[i]+'!n ratio',tit=cltag[ii[i]], position=[0.48*i+0.075,0.15,(i+1)*0.48,0.75]
           oplot, l, wld11*0.+1., line=3;, col=245
           oplot, l, cl2[*,ii[i]]/cl1[*,ii[i]], ns=20
           oplot, l, cl3[*,ii[i]]/cl1[*,ii[i]], col=230, ns=20
           oplot, l, cl2[*,ii[i]]/cl3[*,ii[i]], col=210, ns=20

           oplot, l, nl2[*,ii[i]]/nl1[*,ii[i]], ns=20, line=1
           oplot, l, nl3[*,ii[i]]/nl1[*,ii[i]], col=230, ns=20, line=1
           oplot, l, nl2[*,ii[i]]/nl3[*,ii[i]], col=210, ns=20, line=1
       endfor

print, i

       plot, /nodata, [0,lmax],[1,3000], chars=2, xs=1, ys=1, xtit='!8l',ytit='!8pC!dl!6!u'+cltag[i-1]+'!n',tit=cltag[ii[i-1]], position=[0.48*i+0.09,0.15,(i+1)*0.48,0.75], /ylog
;       oplot, l, wld11*0.+1., line=3, col=245
       oplot, l, cl1[*,ii[i-1]]*ll, ns=20
       oplot, l, cl2[*,ii[i-1]]*ll, col=70, ns=20
       oplot, l, cl3[*,ii[i-1]]*ll, col=210, ns=20

       oplot, l, nl1[*,ii[i-1]]*ll, ns=20
       oplot, l, nl2[*,ii[i-1]]*ll, col=70, ns=20
       oplot, l, nl3[*,ii[i-1]]*ll, col=210, ns=20

;       ianafast, '../../dx11_pre/maps/dx11_pre_ds1_Imap_143_full_uK.fits', dx11, maskfile='data/mask/dx11/Tmask.fits', regression=2, map2_in='../../dx11_pre/maps/dx11_pre_ds2_Imap_143_full_uK.fits'
;       ianafast, '../../dx9/maps/dx9_ds1_Imap_143_nominal_uK.fits', dx9n, maskfile='data/mask/dx11/Tmask.fits', regression=2, map2_in='../../dx9/maps/dx9_ds2_Imap_143_nominal_uK.fits'
;       ianafast, '../../dx9/maps/dx9_ds1_Imap_143_full_uK.fits', dx9f, maskfile='data/mask/dx11/Tmask.fits', regression=2, map2_in='../../dx9/maps/dx9_ds2_Imap_143_full_uK.fits'
;       cl2fits, dx11, 'dx11_pre_ds1_x_ds2_Imap_143_full_uK_pCls.fits'
;       cl2fits, dx9n, 'dx9_ds1_x_ds2_Imap_143_nominal_uK_pCls.fits'
;       cl2fits, dx9f, 'dx9_ds1_x_ds2_Imap_143_full_uK_pCls.fits'

       fits2cl, dx11, 'dx11_pre_ds1_x_ds2_Imap_143_full_uK_pCls.fits'
       fits2cl, dx9n, 'dx9_ds1_x_ds2_Imap_143_nominal_uK_pCls.fits'
       fits2cl, dx9f, 'dx9_ds1_x_ds2_Imap_143_full_uK_pCls.fits'

       oplot, l, dx11*ll, ns=20
       oplot, l, dx9n*ll, col=70, ns=20
       oplot, l, dx9f*ll, col=210, ns=20


       plot, /nodata, [0,1], [0,1], col=255, position=[0.05,0.8,0.95,1]
       xyouts, 0., 0.8, 'dx9-Nominal/dx11', chars=1.5
       xyouts, 0.3, 0.8, 'Signal', chars=1.5
       oplot, [0.4,0.45], [0.8,0.8]
       xyouts, 0., 0.5, 'dx9-Full/dx11', chars=1.5, col=230
       xyouts, 0.3, 0.5, 'Noise', chars=1.5
       oplot, [0.4,0.45], [0.5,0.5], line=1
       xyouts, 0., 0.2, 'dx9-Nominal/dx9-Full', chars=1.5, col=210

       xyouts, 0.55, 0.8, 'dx11', chars=1.5
       xyouts, 0.55, 0.5, 'dx9-Nominal', chars=1.5, col=70
       xyouts, 0.55, 0.2, 'dx9-Full', chars=1.5, col=210
;       oplot, [0.65,0.7], [0.85,0.85], line=0
;       oplot, [0.65,0.7], [0.55,0.55], line=2
;       xyouts, 0.725, 0.8, '!8D!dl!n!6 ratio', chars=1.5
;       xyouts, 0.725, 0.5, '!8W!dl!u2!n!6 ratio', chars=1.5

       if dops then begin
           device, /close
           set_plot, 'x'
       endif
;       write_png, 'preDX11_detset_spectra.png', tvrd(/true)

;help, readme
;print, readme
!p.multi=0
;stop




; ------ Second plot
       if dops then begin
           set_plot, 'ps'
           device, file='beamStudy_preDX11_vs_DX9_detset_pCl-ratio.eps', /col, bits=8, /landscape;, xsize=13500, ysize=5000
           !p.multi = [0,4,1]
       endif else begin
;           !p.multi = [0,4,1]
;           window, 2, xsize=1350, ysize=500
           !p.multi = [0,3,1]
           window, 3, xsize=720.*500/420, ysize=500
       endelse

       readcol, 'data/beams/dx9/dx9_wl_143_ds_1-2.dat', bl, wld9
       readcol, 'data/beams/dx11_pre/dx11_pre_wl_143_ds_1-2.dat', bl, wld11
       readcol, 'data/beams/ffp7/ffp7_wl_143_ds_1-2.dat', bl, wlf7

       cltag = ['TT', 'EE', 'BB', 'TE']

       l = findgen(5001)
       ll = l*(l+1)/2./!pi

i=0
       plot, /nodata, [0,1250],[100,3000], chars=2, xs=1, ys=1, xtit='!8l',ytit='!8pC!dl!6!u'+cltag[i-1]+'!n',tit=cltag[ii[i-1]], position=[0.48*i+0.09,0.15,(i+1)*0.48,0.75], /ylog
;       oplot, l, wld11*0.+1., line=3, col=245
;       oplot, l, cl1[*,ii[i-1]]*ll, ns=20
;       oplot, l, cl2[*,ii[i-1]]*ll, col=70, ns=20
;       oplot, l, cl3[*,ii[i-1]]*ll, col=210, ns=20

;       oplot, l, nl1[*,ii[i-1]]*ll, ns=20
;       oplot, l, nl2[*,ii[i-1]]*ll, col=70, ns=20
;       oplot, l, nl3[*,ii[i-1]]*ll, col=210, ns=20

;       ianafast, '../../dx11_pre/maps/dx11_pre_ds1_Imap_143_full_uK.fits', dx11, maskfile='data/mask/dx11/Tmask.fits', regression=2, map2_in='../../dx11_pre/maps/dx11_pre_ds2_Imap_143_full_uK.fits'
;       ianafast, '../../dx9/maps/dx9_ds1_Imap_143_nominal_uK.fits', dx9n, maskfile='data/mask/dx11/Tmask.fits', regression=2, map2_in='../../dx9/maps/dx9_ds2_Imap_143_nominal_uK.fits'
;       ianafast, '../../dx9/maps/dx9_ds1_Imap_143_full_uK.fits', dx9f, maskfile='data/mask/dx11/Tmask.fits', regression=2, map2_in='../../dx9/maps/dx9_ds2_Imap_143_full_uK.fits'
;       cl2fits, dx11, 'dx11_pre_ds1_x_ds2_Imap_143_full_uK_pCls.fits'
;       cl2fits, dx9n, 'dx9_ds1_x_ds2_Imap_143_nominal_uK_pCls.fits'
;       cl2fits, dx9f, 'dx9_ds1_x_ds2_Imap_143_full_uK_pCls.fits'

       fits2cl, dx11, 'dx11_pre_ds1_x_ds2_Imap_143_full_uK_pCls.fits'
       fits2cl, dx9n, 'dx9_ds1_x_ds2_Imap_143_nominal_uK_pCls.fits'
       fits2cl, dx9f, 'dx9_ds1_x_ds2_Imap_143_full_uK_pCls.fits'

       oplot, l, dx11*ll, ns=20
       oplot, l, dx9n*ll, col=70, ns=20
       oplot, l, dx9f*ll, col=210, ns=20

i=1
       plot, /nodata, [350,850],[500,1200], chars=2, xs=1, ys=1, xtit='!8l',ytit='!8pC!dl!6!u'+cltag[i-1]+'!n',tit=cltag[ii[i-1]], position=[0.48*i+0.09,0.15,(i+1)*0.48,0.75]
;       oplot, l, wld11*0.+1., line=3, col=245
;       oplot, l, cl1[*,ii[i-1]]*ll, ns=20
;       oplot, l, cl2[*,ii[i-1]]*ll, col=70, ns=20
;       oplot, l, cl3[*,ii[i-1]]*ll, col=210, ns=20

;       oplot, l, nl1[*,ii[i-1]]*ll, ns=20
;       oplot, l, nl2[*,ii[i-1]]*ll, col=70, ns=20
;       oplot, l, nl3[*,ii[i-1]]*ll, col=210, ns=20

;       ianafast, '../../dx11_pre/maps/dx11_pre_ds1_Imap_143_full_uK.fits', dx11, maskfile='data/mask/dx11/Tmask.fits', regression=2, map2_in='../../dx11_pre/maps/dx11_pre_ds2_Imap_143_full_uK.fits'
;       ianafast, '../../dx9/maps/dx9_ds1_Imap_143_nominal_uK.fits', dx9n, maskfile='data/mask/dx11/Tmask.fits', regression=2, map2_in='../../dx9/maps/dx9_ds2_Imap_143_nominal_uK.fits'
;       ianafast, '../../dx9/maps/dx9_ds1_Imap_143_full_uK.fits', dx9f, maskfile='data/mask/dx11/Tmask.fits', regression=2, map2_in='../../dx9/maps/dx9_ds2_Imap_143_full_uK.fits'
;       cl2fits, dx11, 'dx11_pre_ds1_x_ds2_Imap_143_full_uK_pCls.fits'
;       cl2fits, dx9n, 'dx9_ds1_x_ds2_Imap_143_nominal_uK_pCls.fits'
;       cl2fits, dx9f, 'dx9_ds1_x_ds2_Imap_143_full_uK_pCls.fits'

       fits2cl, dx11, 'dx11_pre_ds1_x_ds2_Imap_143_full_uK_pCls.fits'
       fits2cl, dx9n, 'dx9_ds1_x_ds2_Imap_143_nominal_uK_pCls.fits'
       fits2cl, dx9f, 'dx9_ds1_x_ds2_Imap_143_full_uK_pCls.fits'

       oplot, l, dx11*ll, ns=20
       oplot, l, dx9n*ll, col=70, ns=20
       oplot, l, dx9f*ll, col=210, ns=20




       plot, /nodata, [0,1], [0,1], col=255, position=[0.05,0.8,0.95,1]
       xyouts, 0., 0.8, 'dx9-Nominal/dx11', chars=1.5
       xyouts, 0.3, 0.8, 'Signal', chars=1.5
       oplot, [0.4,0.45], [0.8,0.8]
       xyouts, 0., 0.5, 'dx9-Full/dx11', chars=1.5, col=230
       xyouts, 0.3, 0.5, 'Noise', chars=1.5
       oplot, [0.4,0.45], [0.5,0.5], line=1
       xyouts, 0., 0.2, 'dx9-Nominal/dx9-Full', chars=1.5, col=210

       xyouts, 0.55, 0.8, 'dx11', chars=1.5
       xyouts, 0.55, 0.5, 'dx9-Nominal', chars=1.5, col=70
       xyouts, 0.55, 0.2, 'dx9-Full', chars=1.5, col=210
;       oplot, [0.65,0.7], [0.85,0.85], line=0
;       oplot, [0.65,0.7], [0.55,0.55], line=2
;       xyouts, 0.725, 0.8, '!8D!dl!n!6 ratio', chars=1.5
;       xyouts, 0.725, 0.5, '!8W!dl!u2!n!6 ratio', chars=1.5

       if dops then begin
           device, /close
           set_plot, 'x'
       endif
;       write_png, 'preDX11_detset_spectra.png', tvrd(/true)

;help, readme
;print, readme
       !p.multi=[0,2,1]
;stop

       plot, /nodata, [0,1750],[0.989,1.011], chars=1.5, xs=1, ys=1, xtit='!8l',ytit='!8pC!dl!6!u'+cltag[i-1]+'!n Ratio',tit=cltag[ii[i-1]], position=[0.1,0.1,0.95,0.85]
       fits2cl, dx11, 'dx11_pre_ds1_x_ds2_Imap_143_full_uK_pCls.fits'
       fits2cl, dx9n, 'dx9_ds1_x_ds2_Imap_143_nominal_uK_pCls.fits'
       fits2cl, dx9f, 'dx9_ds1_x_ds2_Imap_143_full_uK_pCls.fits'
       oplot, l, dx11/dx9f, ns=20, col=245
       oplot, l, dx9n/dx9f, col=70, ns=20
       oplot, l, dx9f/dx9f, line=2
       
       plot, /nodata, [0,1],[0,1], chars=1, xs=1, ys=1, position=[0.1,0.9,0.95,1], col=255
       xyouts, 0,0.5, '!6preDX11!uX!n / dx9!dFull!uX!n', col=245, chars=1.5
       xyouts, 0.5,0.5, '!6dx9!dNominal!uX!n / dx9!dFull!uX!n', col=70, chars=1.5

   endif


   stop, ' --- End of Script ---'

end



