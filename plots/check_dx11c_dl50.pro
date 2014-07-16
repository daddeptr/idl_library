   True = 1b
   False = 0b
   
   xfdir = '/global/scratch2/sd/dpietrob/Software/XFaster/'
   bin = 'dl50'

   freq = [ '100', '143', 'recal_217', 'recal_217' ]
   nf = n_elements(freq)

   root = 'xf_chains_'+freq+'_TP_const_dl50_fisherWins_l2000_Gfg'

   nwdf = xfdir+'uscomp/xfaster/dx11c_SDR_yr_1-2_IQUmap_'+freq+'_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const_dl50_fisherWins_xfcl_l3000_Tmask_badpixmasked_x_Pmask_badpixmasked_l4000.newdat'
   nwdf[3] = xfdir+'uscomp/xfaster/dx11c_SDR_yr_1-2_IQUmap_'+freq[3]+'_full_ns2048_uK_hrhs_const_dl50_fisherWins_xfcl_l3000_Tmask_badpixmasked_x_Pmask_badpixmasked_l4000.newdat'

   bff = xfdir+'uscomp/xfaster/xf_chains_'+freq+'_TP_const_dl50_fisherWins_l2000_Gfg_bestfit.txt'
   fgbff = xfdir+'uscomp/xfaster/xf_chains_'+freq+'_TP_const_dl50_fisherWins_l2000_Gfg_fg_bestfit.txt'

   ll = findgen(5000)+2
   ll = ll*(ll+1.)/2./!pi

   mollview, findgen(12), win=0, px=450
   loadct, 39
   !p.color=0
   !p.background=255
   cols = [0,245,70,210]
   window, 0, xsize=800, ysize=800
   plot, /nodata, [1,2500],[-500,6500], xs=1, ys=1, chars=1.5

   clsar = fltarr(61,nf)
   tar = fltarr(61,nf)
   ear = fltarr(61,nf)
   tear = fltarr(61,nf)


   for i=0,nf-1 do begin
       cls = extract_xfaster_newdat( nwdf[i], lcen=lc, cler=err )
       readcol, bff[i], l, t, te, e
       readcol, fgbff[i], l, tfg, tefg, efg
       t = t/ll
       e = e/ll
       te = te/ll
       
       tfg = tfg/ll
       efg = efg/ll
       tefg = tefg/ll
       
       btfg = xf_binning(reform([0,1,[tfg]]),xfdir+'data/bins/const/const_dl50_TT')
       befg = xf_binning(reform([0,1,efg]),xfdir+'data/bins/const/const_dl50_EE')
       btefg = xf_binning(reform([0,1,tefg]),xfdir+'data/bins/const/const_dl50_TE')

       clsar[*,i] = cls
       tar[*,i] = btfg[0:60]
       ear[*,i] = befg[0:60]
       tear[*,i] = btefg[0:60]

   endfor

   for i=0,nf-1 do begin
       oplot, lc[1:*], clsar[1:*,i], col=cols[i], psym=-4-i
       oplot, lc[1:*], clsar[1:*,i]-tar[1:*,i], col=cols[i], psym=-4-i, thick=2
       oplot, lc[1:*], tar[1:*,i], col=cols[i], line=2 ;, thick=2

   endfor

   ref = ( clsar[*,2] - tar[*,2] )

   window, 1, xsize=800, ysize=800
   plot, /nodata, [1,2500],[-.05,.05], xs=1, ys=1, chars=1.5
   for i=0,nf-1 do begin
       oplot, lc[1:*], (clsar[1:*,i]-tar[1:*,i]-ref[1:*])/ref[1:*], col=cols[i], psym=-4-i, thick=2
       oplot, lc[1:*], (clsar[1:*,i]-tar[1:*,i]-ref[1:*])/ref[1:*], col=cols[i], psym=-4-i, thick=2
   endfor

stop

   if bin eq 'dl50' then binning = 'data/bins/const/const_dl50'
   if bin eq 'ctp' then binning = 'data/bins/ctp/CTP_bin'

;##   ianafast, xfdir + 'uscomp/maps/dx11c_SDR_yr_1-2_IQUmap_143_extMask_545_coTP_undusted_split_ns2048_uK_hr1.fits','dx11c_SDR_yr1_x_yr2_cls_Tmask_badpixmasked.fits', map2_in=xfdir+'uscomp/maps/dx11c_SDR_yr_1-2_IQUmap_143_extMask_545_coTP_undusted_split_ns2048_uK_hr2.fits', simul_type=1, maskfile=xfdir+'uscomp/masks/output/Tmask_badpixmasked.fits', regression=2, /silent

;##   ianafast, xfdir + 'data/maps/dx11c/dx11c_SDR_yr_1-2_IQUmap_recal_217_extMask_545_coTP_undusted_split_ns2048_uK_hr1.fits','dx11c_SDR_yr1_x_yr2_recal_217_cls_Tmask_badpixmasked.fits', map2_in=xfdir+'data/maps/dx11c/dx11c_SDR_yr_1-2_IQUmap_217_extMask_545_coTP_undusted_split_ns2048_uK_hr2.fits', simul_type=1, maskfile=xfdir+'uscomp/masks/output/Tmask_badpixmasked.fits', regression=2

;## ;   spawn, 'mv dx11c_SDR_yr1_x_yr2_cls_Tmask_badpixmasked.fits dx11c_SDR_yr1_x_yr2_143_cls_Tmask_badpixmasked.fits '
;## ;   spawn, 'mv dx11c_SDR_yr1_x_yr2_m-like_cls_Tmask_badpixmasked.fits dx11c_SDR_yr1_x_yr2_143_m-like_cls_Tmask_badpixmasked.fits'
;##   fits2cl, cls, 'dx11c_SDR_yr1_x_yr2_143_cls_Tmask_badpixmasked.fits'
   fits2cl, cls, 'dx11c_SDR_yr1_x_yr2_recal_217_cls_Tmask_badpixmasked.fits'
   
;##   kcls = deconvolve_kernel(cls,inv_fkernel=xfdir+'uscomp/xfaster/Tmask_badpixmasked_x_Pmask_badpixmasked_inv-kernel_l4000_v1.9.fits', write_cls=xfdir+'uscomp/xfaster/dx11c_SDR_yr1_x_yr2_m-like_cls_Tmask_badpixmasked.fits')

;##   kcls = deconvolve_kernel(cls,inv_fkernel=xfdir+'uscomp/xfaster/Tmask_badpixmasked_x_Pmask_badpixmasked_inv-kernel_l4000_v1.9.fits', write_cls=xfdir+'uscomp/xfaster/dx11c_SDR_yr1_x_yr2_217_m-like_cls_Tmask_badpixmasked.fits')

   fits2cl, kcls, xfdir+'uscomp/xfaster/dx11c_SDR_yr1_x_yr2_143_m-like_cls_Tmask_badpixmasked.fits'
   readcol,xfdir+'data/beams/dx11c/dx11c_wl_143.dat', l, bl
   hpw = healpixwindow(2048)
   kcls = kcls / bl^2 / hpw^2
   bkcls143 = xf_binning(kcls,xfdir+binning+'_TT')

   fits2cl, kcls, xfdir+'uscomp/xfaster/dx11c_SDR_yr1_x_yr2_217_m-like_cls_Tmask_badpixmasked.fits'
   readcol,xfdir+'data/beams/dx11c/dx11c_wl_217.dat', l, bl
   hpw = healpixwindow(2048)
   kcls = kcls / bl^2 / hpw^2
   bkcls217 = xf_binning(kcls,xfdir+binning+'_TT')

   if bin eq 'dl50' then begin
       f1 = 'dx11c_SDR_yr_1-2_IQUmap_143_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const_dl50_fisherWins_xfcl_l3000_Tmask_badpixmasked_x_Pmask_badpixmasked_l4000.newdat'
       bcls143 = extract_xfaster_newdat(f1, lcen=lc)

       f2 = 'dx11c_SDR_yr_1-2_IQUmap_recal_217_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const_dl50_fisherWins_xfcl_l3000_Tmask_badpixmasked_x_Pmask_badpixmasked_l4000.newdat'
       bcls217 = extract_xfaster_newdat(f2, lcen=lc, btcl=btcl, binning=xfdir+binning)
   endif

   if bin eq 'ctp' then begin
       f1 = 'dx11c_SDR_yr_1-2_IQUmap_143_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_ctp_fisherWins_xfcl_l3000_Tmask_badpixmasked_x_Pmask_badpixmasked_l4000.newdat'
       bcls143 = extract_xfaster_newdat(f1, lcen=lc)

       f2 = 'dx11c_SDR_yr_1-2_IQUmap_recal_217_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_ctp_fisherWins_xfcl_l3000_Tmask_badpixmasked_x_Pmask_badpixmasked_l4000.newdat'
       bcls217 = extract_xfaster_newdat(f2, lcen=lc, btcl=btcl, binning=xfdir+binning)
   endif

   !p.multi=[0,1,3]
   window,0,ysize=800, xsize=800
   plot, lc, bcls143, chars=2, yr=[-150, 6500], ys=1, xtit='!8l', ytit='!6D!dl!n'
   oplot, lc, bcls143, col=215, thick=2
   oplot, lc, bkcls143, col=245, thick=2
   oplot, lc, bcls217, col=100, thick=2
   oplot, lc, bkcls217, col=70, thick=2
   legend, ['XF143','Ml143','XF217','Ml217'],col=[215,245,100,70],psym=[4,4,4,4], thick=2, chars=1.25, /right

   plot, lc, lc*0, line=2, yr=[-115,115], chars=2, ys=1, tit='!6XFaster vs kernel-inversion', xtit='!8l', ytit='!7D!6D!dl!n'
   oplot, lc, bcls143-bkcls143, psym=-4, thick=2
   oplot, lc, bcls217-bkcls217, psym=-4, col=70, thick=2
   oplot, lc, btcl*0.025-100, col=245 
   legend, ['XF143-Ml143','XF217-Ml217','2013 Planck (x0.025)'], col=[0,70,245], chars=1.25, psym=[4,4,4], thick=2

   plot, lc, lc*0, line=2, yr=[-.05,.05], chars=2, ys=1, tit='!6XFaster vs kernel-inversion', xtit='!8l', ytit='!7D!6D!dl!n/D!dl!n'
   oplot, lc, (bcls143-bkcls143)/btcl, psym=-4, thick=2
   oplot, lc, (bcls217-bkcls217)/btcl, psym=-4, col=70, thick=2
;##   oplot, lc, btcl*0.025-100, col=245 
   legend, ['XF143-Ml143','XF217-Ml217'], col=[0,70], chars=1.25, psym=[4,4], thick=2
   write_png,'xf_dx11c_check_143_'+bin+'.png', tvrd(/true)

; python sequence to get bestfit:
; root='xf_chains_143_TP_const_dl50_fisherWins_l2000_Gfg'
; dir='/global/scratch2/sd/marius/us2/runs/uscomp/xfaster'
; pars=g.read_pars('',filename='results/'+root+'.likestats',dir=dir)
; g.write_pars(pars,root='',dir=dir,filetag='_143_TP_const_dl50_fisherWins_l2000_Gfg')
; g.run_camb('',cambfile='xf_camb_143_TP_const_dl50_fisherWins_l2000_Gfg.ini', dir=dir) 
; g.compute_bestfit('',dir=dir,cambparfile='xf_camb_143_TP_const_dl50_fisherWins_l2000_Gfg.ini',clsfile='xf_chains_143_TP_const_dl50_fisherWins_l2000_Gfg_lensedCls.dat',filetag='xf_chains_143_TP_const_dl50_fisherWins_l2000_Gfg_')
  
   


end
