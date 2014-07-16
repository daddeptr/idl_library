 f='xfaster/dx11c_SDR_545-353_templates_yrs_ctp_fisherWins_xfcl_l3000_Tmask_badpixmasked_x_Pmask_badpixmasked_l4000.newdat'                                      
 c=extract_xfaster_newdat(f,lcen=lc)                                                                                                                          

 mollview, findgen(12), win=0, px=800
 loadct, 39
 !p.color=0
 !p.background=255

 plot, lc,c, chars=1.5, /xlog, xr=[10,2500], /ylog, xs=1, xtit='!8l', ytit='!6D!d!8l!n', tit='545 GHz XFaster spectrum'
 oplot, lc, c, psym=4
 a = (lc/lc[16])^(-1)*c[16]
 b = (lc/lc[110])^(0.01)*c[110]*0.7
 d = (lc/lc[185])^(4)*c[185]*0.8

 oplot, lc, a, col=245
 oplot, lc, b, col=210
 oplot, lc, d, col=70
;
 oplot, lc, a + b + d, col=40

 xyouts, 30, 7.e7, '!8D!dl!n=c!d1!nl!6!u-1!n', col=245, chars=2
 xyouts, 30, 5.e7, '!8D!dl!n=c!d2!nl!6!u-0.01!n', col=210, chars=2
 xyouts, 30, 3.4e7, '!8D!dl!n=c!d3!nl!6!u4!n', col=70, chars=2

end
