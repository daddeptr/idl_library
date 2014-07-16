   True  = 1b
   False = 0b

   do_tt = True
   do_ee = True
   wrt2m = True

   if do_tt then begin
       readcol, 'outputs/dx11c_coTP-clean_maps/dx11c_143_egfg-sub_xfcls_TT.txt',lc,cls,err,egfg,clssub
       readcol, 'outputs/dx11c_coTP-clean_maps/dx11c_217_egfg-sub_xfcls_TT.txt',lc2,cls2,err2,egfg2,clssub2

       ex = extract_xfaster_newdat( 'outputs/dx11c_coTP-clean_maps/dx11c_yr_1-2_IQUmap_143_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_xfcl_l3000_cons-apo_Tmask3_mPix_x_cons-apo_Pmask3_mPix_l4000.newdat', btcl=btcl, binning='data/bins/const/const2' )
       
       plik = mrdfits( 'data/maps/dx11c/plik/plik_of_TT_143x143.fits', 1, hdr )
;##       hilli = mrdfits( 'data/maps/dx11c/hillipop/cross_ForGraca_2_3.fits',1, hdr2 )
       errhilli = mrdfits( 'data/maps/dx11c/hillipop/cross_ForGraca_2_3.fits',1, hdr2 )
       hilli = mrdfits( 'data/maps/dx11c/hillipop/noForegrounds_ForGraca_143x143.fits',1, hdr2 )

       print, n_elements(plik[0,*])
       print, n_elements(plik[1,*])

       fplik = fltarr(50+n_elements(plik[0,*]))
       fplik[50:*] = plik[1,*]
       l = findgen(3001)
       ll = l*(l+1)/2./!pi
       bfplik = xf_binning( fplik/ll, 'data/bins/const/const2_TT')

;
       fhilli = fltarr(2+n_elements(hilli.(1)))
       ferrhilli = fltarr(2+n_elements(errhilli.(2)))
       fhilli[2:*] = hilli.(1)
       ferrhilli[2:*] = errhilli.(2)
       ferrhilli = ( ferrhilli * 1.e12 / ll )^2
;##       bfhilli = xf_binning( fhilli*1.e12/ll, 'data/bins/const/const2_TT')
       bfhilli = xf_binning( fhilli/ll, 'data/bins/const/const2_TT')
       bferrhilli = bp_binning( ferrhilli, 'data/bins/const/const2_TT' )
       bferrhilli = sqrt( bferrhilli[*]/30. ) * lc * (lc+1)/2./!pi 
       
       meancls = (clssub+bfplik+bfhilli)/3

       !p.multi=[0,1,2]
       window,0,xsize=900, ysize=900
       plot, lc[2:*], btcl[2:*], chars=1.5, yr=[-500,6500], ys=1, xr=[1,2050], xs=1, tit='DX11c 143 GHz'
       oplot, lc[2:*], btcl[2:*], col=245
       oplot, lc[2:*], clssub[2:*], psym=-4, col=70, thick=2
;       oplot, lc[2:*], clssub2[2:*], psym=-2, col=210, thick=2
;       oplot, lc[2:*], meancls[2:*]+err[2:*], line=2
;       oplot, lc[2:*], meancls[2:*]-err[2:*], line=2
       oplot, lc[2:*], bfplik[2:*], psym=-6, col=210, thick=2
       oplot, lc[2:*], bfhilli[2:*], psym=-5, col=120, thick=2
;##       legend, ['XFaster','XFaster 217','Plik','Hillipop'],psym=[4,2,6,5],col=[70,210,100,180], /top, /right, chars=2.5, thick=2
       legend, ['XFaster','Plik','Hillipop'],psym=[4,6,5],col=[70,210,120], /top, /right, chars=2.5, thick=2
       
       if wrt2m then btcl = meancls
       if wrt2m then plot, lc[2:*], btcl[2:*]*0., chars=1.5, yr=[-250,250], ys=1, xr=[1,2050], xs=1, tit='Difference wrt mean spectrum' else $
         plot, lc[2:*], btcl[2:*]*0., chars=1.5, yr=[-150,550], ys=1, xr=[1,2050], xs=1, tit='Difference wrt Planck 2013'
       oplot, lc[2:*], btcl[2:*]*0., col=245
;       errplot, lc[2:*], meancls[2:*]-btcl[2:*]-err[2:*], meancls[2:*]-btcl[2:*]+err[2:*], thick=2
;##       oplot, lc[2:*], meancls[2:*]-btcl[2:*]+err[2:*], line=2
;##       oplot, lc[2:*], meancls[2:*]-btcl[2:*]-err[2:*], line=2
       oplot, lc[2:*], clssub[2:*]-btcl[2:*], psym=4, col=70, thick=2
       errplot, lc[2:*], clssub[2:*]-btcl[2:*]-err[2:*], clssub[2:*]-btcl[2:*]+err[2:*], col=70
;       oplot, lc[2:*], clssub[2:*]-btcl[2:*], col=70

       oplot, lc[2:*], bfplik[2:*]-btcl[2:*], psym=6, col=210, thick=2
;       oplot, lc[2:*], bfplik[2:*]-btcl[2:*], col=100

       oplot, lc[2:*]*1.01, bfhilli[2:*]-btcl[2:*], psym=5, col=120, thick=2
       errplot, lc[2:*]*1.01, bfhilli[2:*]-btcl[2:*]-bferrhilli, bfhilli[2:*]-btcl[2:*]+bferrhilli, col=120
;       oplot, lc[2:*], bfhilli[2:*]-btcl[2:*], col=180

;       oplot, lc[2:*], clssub2[2:*]-btcl[2:*], psym=2, col=210, thick=2
;       oplot, lc[2:*], clssub2[2:*]-btcl[2:*], col=210

       if wrt2m then write_png,'dx11c_method_comp_TT_wrtMean.png',tvrd(/true) else write_png,'dx11c_method_comp_TT_wrt2013.png',tvrd(/true)
   endif

stop

   if do_ee then begin
       readcol, 'outputs/dx11c_coTP-clean_maps/dx11c_143_egfg-sub_xfcls_EE.txt',lc,cls,err,egfg,clssub
       readcol, 'outputs/dx11c_coTP-clean_maps/dx11c_217_egfg-sub_xfcls_EE.txt',lc2,cls2,err2,egfg2,clssub2

       ex = extract_xfaster_newdat( 'outputs/dx11c_coTP-clean_maps/dx11c_yr_1-2_IQUmap_143_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_xfcl_l3000_cons-apo_Tmask3_mPix_x_cons-apo_Pmask3_mPix_l4000.newdat', btcl=btcl, binning='data/bins/const/const2', ncl=2 )

       l = findgen(3001)
       ll = l*(l+1)/2./!pi

;   plik = mrdfits( 'data/maps/dx11c/plik/plik_of_TT_143x143.fits', 1, hdr )
;   fplik = fltarr(50+n_elements(plik[0,*]))
;   fplik[50:*] = plik[1,*]
;   bfplik = xf_binning( fplik/ll, 'data/bins/const/const2_EE')

;
;##       hilli = mrdfits( 'data/maps/dx11c/hillipop/noForegrounds_ForGraca_143x143.fits',1, hdr2 )
       errhilli = mrdfits( 'data/maps/dx11c/hillipop/cross_ForGraca_2_3.fits',2, hdr2 )

       fhilli = fltarr(2+n_elements(errhilli.(1)))
       fhilli[2:*] = errhilli.(1)
       bfhilli = xf_binning( fhilli*1.e12/ll, 'data/bins/const/const2_EE')

       ferrhilli = fltarr(2+n_elements(errhilli.(2)))
       ferrhilli[2:*] = errhilli.(2)
       ferrhilli = ( ferrhilli * 1.e12 / ll )^2
;##       bfhilli = xf_binning( fhilli*1.e12/ll, 'data/bins/const/const2_TT')
       bferrhilli = bp_binning( ferrhilli, 'data/bins/const/const2_EE' )
       bferrhilli = sqrt( bferrhilli[*]/60. ) * lc * (lc+1)/2./!pi 

       !p.multi=[0,1,2]
       window,0,xsize=900, ysize=900
       plot, lc[2:*], btcl[2:*], chars=1.5, yr=[-5,50], ys=1, xr=[1,1550], xs=1, tit='DX11c 143 GHz'
       oplot, lc[2:*], btcl[2:*], col=245
       oplot, lc[2:*], clssub[2:*], psym=-4, col=70, thick=2
;##       oplot, lc[2:*], clssub2[2:*], psym=-6, col=210, thick=2
;   oplot, lc[2:*], bfplik[2:*], psym=-6, col=100, thick=2
       oplot, lc[2:*], bfhilli[2:*], psym=-5, col=120, thick=2
;##       legend, ['XFaster','XFaster 217','Hillipop'],psym=[4,6,5],col=[70,210,180], /top, /left, chars=2.5, thick=2
       legend, ['XFaster','Hillipop'],psym=[4,5],col=[70,120], /top, /left, chars=2.5, thick=2
   
       meancls = (clssub+bfhilli)/2
       if wrt2m then btcl = meancls

       if wrt2m then plot, lc[2:*], btcl[2:*]*0., chars=1.5, yr=[-5,5], ys=1, xr=[1,2550], xs=1, tit='Difference wrt mean spectrum' else $
         plot, lc[2:*], btcl[2:*]*0., chars=1.5, yr=[-5,5], ys=1, xr=[1,2550], xs=1, tit='Difference wrt Planck 2013'
       oplot, lc[2:*], btcl[2:*]*0., col=245
;       errplot, lc[2:*], meancls[2:*]-btcl[2:*]-err[2:*], meancls[2:*]-btcl[2:*]+err[2:*], thick=2
;   oplot, lc[2:*], meancls[2:*]-btcl[2:*]+err[2:*], line=2
;   oplot, lc[2:*], meancls[2:*]-btcl[2:*]-err[2:*], line=2
       oplot, lc[2:*], clssub[2:*]-btcl[2:*], psym=-4, col=70, thick=2
       errplot, lc[2:*], clssub[2:*]-btcl[2:*]-err[2:*], clssub[2:*]-btcl[2:*]+err[2:*], col=70
;   oplot, lc[2:*], clssub2[2:*], psym=-6
;       oplot, lc[2:*], bfplik[2:*]-btcl[2:*], psym=-6, col=100, thick=2
       oplot, lc[2:*], bfhilli[2:*]-btcl[2:*], psym=-5, col=120, thick=2
       errplot, lc[2:*], bfhilli[2:*]-btcl[2:*]-bferrhilli[2:*], bfhilli[2:*]-btcl[2:*]+bferrhilli[2:*], col=120
       if wrt2m then write_png,'dx11c_method_comp_EE_wrtMean.png',tvrd(/true) else write_png,'dx11c_method_comp_EE_wrt2013.png',tvrd(/true) 
   endif

stop

end
