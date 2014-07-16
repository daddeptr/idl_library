   True = 1b
   False = 0b

   lmax = 3000

   ave_cls = fltarr(4097,6)
   ave_nls = fltarr(4097,6)
   nsims = 100

   if False then begin
       for isim=0,nsims-1 do begin
           fits2cl, cls, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/ffp7/sims/m143/pcls/ffp7_cmb_sim_143_ns2048_uK_hrhs_cls_'+strtrim(string(isim+1),2)+'.fits'
           fits2cl, nls, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/ffp7/sims/m143/pcls/ffp7_cmb_sim_143_ns2048_uK_hrhd_cls_'+strtrim(string(isim+1),2)+'.fits'
           ave_cls += cls
           ave_nls += nls
       endfor

       ave_cls = ave_cls / nsims
       ave_nls = ave_nls / nsims

       for i=0,5 do ave_cls[*,i] = (ave_cls[*,i]-ave_nls[*,i]) / gaussbeam(7.03,4096)^2 / healpixwindow(2048)^2
       ave_file = '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/ffp7/sims/m143/pcls/ffp7_cmb_sim_143_ns2048_uK_hrhs_ave_cls.fits'
       cl2fits, ave_cls, ave_file
       l=findgen(4097)
       ll=l*(l+1)/2./!pi
       ttt = ave_cls[*,0] * ll
       tee = ave_cls[*,1] * ll
       tte = ave_cls[*,3] * ll
   endif else begin

;##       readcol, '/global/scratch2/sd/dpietrob/Software/XFaster/data/camb_97283428_scalcls.dat', l, ttt, tee, tte, format='f,f,f,f'
;##       ave_file = '/global/scratch2/sd/dpietrob/Software/XFaster/data/camb_97283428_scalcls_uK.fits'
;       readcol, '/global/scratch2/sd/dpietrob/Software/XFaster/data/planck_base_lensedCls.dat', l, ttt, tee, tte, format='f,f,f,x,f'
       readcol, '/global/scratch2/sd/dpietrob/Software/XFaster/data/camb_97283428_lensedcls.dat', l, ttt, tee, tte, format='f,f,f,x,f'
       ave_file = '/global/scratch2/sd/dpietrob/Software/XFaster/data/camb_97283428_lensedcls_uK.fits'
   endelse

   dir = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/mock_cmb_white_noise_ns2048_uK_hrhs_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000/'

   root = 'mock_lenscmb_white_noise_b7.03_ns2048_uK_hrhs_1_const2_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'
   froot = 'mock_lenscmb_white_noise_b7.03_ns2048_uK_hrhs_1_const2_aveNls_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'
;##   froot = 'mock_lenscmb_white_noise_b7.03_ns2048_uK_hrhs_1_const2_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'

   file = dir+root+'/'+froot+'.newdat'
   mollview, findgen(12), px=450, win=5
   loadct, 39
   !p.color=0
   !p.background=255
   !p.multi=[0,3,2]
   window, 5, xsize=1200, ysize=650

;### ------- TT
   tt = extract_xfaster_newdat( file, lcen=ltt, binning='/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2', btcl=btttcl, tclfile=ave_file, cler=ertt )

;##   fits2cl, xxx, ave_file
;##   btttcl2 = xf_binning(xxx[*,0],'/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT')
;##   l=findgen(4001)
;##   ll=l*(l+1)/2./!pi
;##   btttcl3 = bp_binning(xxx[*,0]*ll,'/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT')

   cler = ertt^2


   plot, ltt, tt, chars=1.5, ytit='!8D!dl!uTT!n', psym=3, yr=[-550,6550], ys=1, position=[0.06,0.45,0.325,0.975], xtickname=strarr(7)+' ', xr=[1,lmax+50], xs=1
;##   oplot, l, ttt, col=245
   oplot, ltt, btttcl, thick=2, col=245

;stop
   ttar = dblarr(n_elements(tt),nsims)
   ttar[*,0] = tt
   mtt = tt

   for i=2l,nsims do begin
       root = 'mock_lenscmb_white_noise_b7.03_ns2048_uK_hrhs_'+strtrim(string(i),2)+'_const2_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'
       froot = 'mock_lenscmb_white_noise_b7.03_ns2048_uK_hrhs_'+strtrim(string(i),2)+'_const2_aveNls_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'
;##       froot = 'mock_lenscmb_white_noise_b7.03_ns2048_uK_hrhs_'+strtrim(string(i),2)+'_const2_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'
       file = dir+root+'/'+froot+'.newdat'
;##       print, file
       tt = extract_xfaster_newdat( file, cler=ertt )
       cler += ertt^2
       ttar[*,i-1] = tt
       oplot, ltt, tt, psym=3
       mtt += tt
   endfor
   cler /= nsims
   cler = sqrt(cler)
   mtt /= nsims
   mdtt = mtt*0.
   for i=0l,n_elements(tt)-1 do mdtt[i] = median(ttar[i,*])
   oplot, ltt, mtt, col=210, thick=2
   oplot, ltt, mdtt, col=90, thick=2, line=2
   xyouts, 1250, 5500, '!6Planck 2013', col=245, chars=1.5
   xyouts, 1250, 5000, '!6143 mock mean', col=210, chars=1.5
   xyouts, 1250, 4500, '!6143 mock median', col=90, chars=1.5

;33 ------- EE
   ee = extract_xfaster_newdat( file, lcen=lee, binning='/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2', btcl=bteecl, ncl=2, tclfile=ave_file )
   eear = fltarr(n_elements(ee),nsims)
   eear[*,0] = ee
   mee = ee
   plot, lee, ee, chars=1.5, ytit='!8D!dl!uEE!n', psym=3, yr=[-5.5,50.5], ys=1, position=[0.3875,0.45,0.6525,0.975], xtickname=strarr(7)+' ', xr=[1,lmax+50], xs=1
   ;##oplot, l, tee, col=245, thick=2
   oplot, lee, bteecl, col=245, thick=2

   for i=2,nsims do begin
       root = 'mock_lenscmb_white_noise_b7.03_ns2048_uK_hrhs_'+strtrim(string(i),2)+'_const2_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'
       ;## froot = 'mock_lenscmb_white_noise_b7.03_ns2048_uK_hrhs_'+strtrim(string(i),2)+'_const2_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'
       froot = 'mock_lenscmb_white_noise_b7.03_ns2048_uK_hrhs_'+strtrim(string(i),2)+'_const2_aveNls_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'
       file = dir+root+'/'+froot+'.newdat'
       ee = extract_xfaster_newdat( file, ncl=2 )
       eear[*,i-1] = ee
       oplot, lee, ee, psym=3
       mee += ee
   endfor
   mee /= nsims
   mdee = mee*0.
   for i=0,n_elements(ee)-1 do mdee[i] = median(eear[i,*])
   oplot, lee, mee, col=210, thick=2
   oplot, lee, mdee, col=90, thick=2, line=2

;## -------- TE
   te = extract_xfaster_newdat( file, lcen=lte, binning='/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2', btcl=bttecl, ncl=4, tclfile=ave_file )
   mte = te
   tear = fltarr(n_elements(te),nsims)
   tear[*,0] = te
   plot, lte, te, chars=1.5, ytit='!8D!dl!uTE!n', psym=3, yr=[-155,155], ys=1, position=[0.715,0.45,0.98,0.975], xtickname=strarr(7)+' ', xr=[1,lmax+50], xs=1
   ;## oplot, l, tte, col=245, thick=2
   oplot, lte, bttecl, col=245, thick=2

   for i=2,nsims do begin
       root = 'mock_lenscmb_white_noise_b7.03_ns2048_uK_hrhs_'+strtrim(string(i),2)+'_const2_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'
       froot = 'mock_lenscmb_white_noise_b7.03_ns2048_uK_hrhs_'+strtrim(string(i),2)+'_const2_aveNls_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'
       file = dir+root+'/'+froot+'.newdat'
       te = extract_xfaster_newdat( file, ncl=4 )
       oplot, lte, te, psym=3
       tear[*,i-1] = te
       mte += te
   endfor
   mte /= nsims
   mdte = mte*0.
   for i=0,n_elements(te)-1 do mdte[i] = median(tear[i,*])
   oplot, lte, mte, col=210, thick=2
   oplot, lte, mdte, col=90, thick=2, line=2


;## ------ Residuals

;   fits2cl, tnls, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/mock/mock_noise_tnls_143.fits'
;   l = findgen(4096)
;   ll= l*(l+1)/2./!pi
;   bnltt = bp_binning(tnls[*,0]*ll/gaussbeam(7.03,4096)^2,'/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT')
;   bnlee = bp_binning(tnls[*,1]*ll/gaussbeam(7.03,4096)^2,'/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_EE')

;##   plot, ltt, (mtt)/btttcl, chars=1.5, yr=[.945,1.055], position=[0.06,0.075,0.325,0.44], ys=1, psym=-6, ytit='!6TT: <!8D!dl!n!6>-!8D!dl!uTh!n'
   plot, ltt, (mtt)/btttcl, chars=1.5, yr=[.975,1.025], position=[0.06,0.075,0.325,0.44], ys=1, psym=-6, ytit='!6TT: <!8D!dl!n!6>-!8D!dl!uTh!n', xr=[1,lmax+50], xs=1
   oplot, ltt, mtt*0+1, line=2
   oplot, ltt, (mdtt)/btttcl, col=90, psym=-5
;##   oplot, lbtt, (mdtt)/btttcl3, col=90, psym=-5
   xyouts, 550, 0.95, '!6Mean-Bestfit', chars=1.5
   xyouts, 550, 0.925, '!6Median-Bestfit', col=90, chars=1.5
   perr = mtt/btttcl
   oplot, ltt, 1.-(cler)/btttcl, line=3
   oplot, ltt, 1.+(cler)/btttcl, line=3
;##   result = linfit(ltt[2:66], perr[2:66],yfit=fit, chisq=c2)
   result = linfit(ltt[2:66], perr[2:66],yfit=fit, chisq=c2, measure_error=cler[2:66]/btttcl[2:66])
   print, result
   oplot, ltt[2:66], fit, col=240, thick=2
;##   xyouts, 200, 1.15, '<C!dl!n/model>', chars=1.25
   xyouts, 200, 1.11, 'Linear fit: A='+string(result[0],format='(f6.3)')+' b='+string(result[1],format='(e9.2)'), chars=1.25, col=240
   result = linfit(alog(ltt[2:66]), alog(perr[2:66]),yfit=fit, chisq=c2, measure_error=alog(cler[2:66]/btttcl[2:66]) )
   print, result
;stop
   plot, lee, (mee)/bteecl, chars=1.5, yr=[.85,1.15], position=[0.3875,0.075,0.6525,0.44], ys=1, psym=-6, ytit='!6EE: <!8D!dl!n!6>-!8D!dl!uTh!n', xr=[1,lmax+50], xs=1
   oplot, lee, mee*0+1, line=2
   oplot, lee, (mdee)/bteecl, col=90, psym=-5
;##   oplot, lee, mee/bnlee, col=210

   plot, lte, (mte)/bttecl, chars=1.5, yr=[0.85,1.15], position=[0.715,0.075,0.98,0.44], ys=1, psym=-6, ytit='!6TE: <!8D!dl!n!6>-!8D!dl!uTh!n', xr=[1,lmax+50], xs=1
   oplot, lte, mte*0+1, line=2
   oplot, lte, (mdte)/bttecl, col=90, psym=-5


end
