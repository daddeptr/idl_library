   True = 1b
   False = 0b

   ave_file = '/global/scratch2/sd/dpietrob/Software/XFaster/data/planck_lcdm_cl_uK_xf1.e-3.fits'


   if False then begin
       ave_file = '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/ffp7/sims/m143/pcls/ffp7_cmb_sim_143_ns2048_uK_hrhs_ave_cls.fits'
       if False then begin
           ave_cls = fltarr(4097,6)
           ave_nls = fltarr(4097,6)
           nsims = 99
           for isim=0,nsims-1 do begin
               fits2cl, cls, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/ffp7/sims/m143/pcls/ffp7_cmb_sim_143_ns2048_uK_hrhs_cls_'+strtrim(string(isim+1),2)+'.fits'
               fits2cl, nls, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/ffp7/sims/m143/pcls/ffp7_cmb_sim_143_ns2048_uK_hrhd_cls_'+strtrim(string(isim+1),2)+'.fits'
               ave_cls += cls
               ave_nls += nls
           endfor
           ave_cls = ave_cls / nsims
           ave_nls = ave_nls / nsims

           readcol,'/global/scratch2/sd/dpietrob/Software/XFaster/data/beams/ffp7/ffp7_wl_143.dat', l, blt, blq, blu
           bl=[ [blt], [blq], [blu] ]
           hpw = healpixwindow(2048,3)
           lmax = min( [n_elements(l), n_elements(hpw[*,0]),n_elements(ave_cls[*,0]) ] )
           for i=0,5 do begin
               if i lt 3 then ave_cls[0:lmax-1,i] = (ave_cls[0:lmax-1,i]-ave_nls[0:lmax-1,i]) / bl[0:lmax-1,i]^2 / hpw[0:lmax-1,i]^2
               if i eq 3 then ave_cls[0:lmax-1,i] = (ave_cls[0:lmax-1,i]-ave_nls[0:lmax-1,i]) / bl[0:lmax-1,0]/bl[0:lmax-1,1] / hpw[0:lmax-1,0]/hpw[0:lmax-1,1]
               if i eq 4 then ave_cls[0:lmax-1,i] = (ave_cls[0:lmax-1,i]-ave_nls[0:lmax-1,i]) / bl[0:lmax-1,1]/bl[0:lmax-1,2] / hpw[0:lmax-1,0]/hpw[0:lmax-1,2]
               if i eq 5 then ave_cls[0:lmax-1,i] = (ave_cls[0:lmax-1,i]-ave_nls[0:lmax-1,i]) / bl[0:lmax-1,1]/bl[0:lmax-1,2] / hpw[0:lmax-1,2]/hpw[0:lmax-1,2]
           endfor
           cl2fits, ave_cls, ave_file
       endif
       fits2cl, ave_cls, ave_file
       l=findgen(4097)
       ll=l*(l+1)/2./!pi
       att = ave_cls[*,0] * ll
       aee = ave_cls[*,1] * ll 
       ate = ave_cls[*,3] * ll

       readcol, '/global/scratch2/sd/dpietrob/Software/XFaster/data/planck_base_lensedCls.dat', l, ttt, tee, tte, format='f,f,f,x,f'
       ;## readcol, '/global/scratch2/sd/dpietrob/Software/XFaster/data/camb_97283428_scalcls.dat', l, ttt, tee, tte, format='f,f,f,f'
;       readcol,
;       '/global/scratch2/sd/dpietrob/Software/XFaster/data/planck_base_lensedCls.dat', l, ttt, tee, tte, format='f,f,f,x,f'
       ;## readcol,'/global/scratch2/sd/dpietrob/Software/XFaster/data/camb_97283428_lensedcls.dat',lc, tttc, teec, ttec, format='f,f,f,x,f'
       ;## ave_file ='/global/scratch2/sd/dpietrob/Software/XFaster/data/camb_97283428_lensedcls_uK.fits'
 endif

 nsims = 99
 !p.multi=0

 fits2cl, ave_cls, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/ffp7/sims/m143/pcls/ffp7_cmb_sim_143_ns2048_uK_hrhs_ave_cls.fits'

 l=findgen(4097)
 ll=l*(l+1)/2./!pi
 att = ave_cls[*,0] * ll
 aee = ave_cls[*,1] * ll 
 ate = ave_cls[*,3] * ll
 
 readcol, '/global/scratch2/sd/dpietrob/Software/XFaster/data/planck_base_lensedCls.dat', l, ttt, tee, tte, format='f,f,f,x,f'
; readcol, '/global/scratch2/sd/dpietrob/Software/XFaster/data/planck_base_lensedCls.dat', l, ttt, tee, tte, format='f,f,f,x,f'
; window, 0
;   plot, l, ttt/att[2:*], yr=[0.9,1.1], chars=1.5, ns=30
;   oplot, l, l*0.+1., line=2
;   oplot, gaussbeam(0.25,4096)^2, col=245
;stop

 dir = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/ffp7_sn/xf_v1.9/'

   ;## root = 'ffp7_cmb_sim_143_flat_ns2048_uK_hrhs_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'
 root = 'ffp7_cmb_sim_217_flat_ns2048_uK_hrhs_ffp7-beam_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'
   ;## root = 'ffp7_cmb_sim_143_flat_ns2048_uK_hrhs_7.07-beam_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'

 fsim = 1

   file = dir+root+'_'+strtrim(string(fsim),2)+'.newdat'
   mollview, findgen(12), px=450, win=4
   loadct, 39
   !p.color=0
   !p.background=255
   !p.multi=[0,3,2]
   window, 4, xsize=1200, ysize=650

;### ------- TT
   tt = extract_xfaster_newdat( file, lcen=ltt, binning='/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2', btcl=btttcl, tclfile=ave_file )
   plot, ltt, tt, chars=1.5, ytit='!8D!dl!uTT!n', psym=3, yr=[-550,6550], ys=1, position=[0.06,0.45,0.325,0.975], xtickname=strarr(7)+' '
   oplot, l, ttt, col=245, thick=2

   ttar = fltarr(n_elements(tt),nsims)
   ttar[*,0] = tt
   mtt = tt

   for i=fsim+1,nsims do begin
       file = dir+root+'_'+strtrim(string(i),2)+'.newdat'
;##       print, file
       tt = extract_xfaster_newdat( file )
       ttar[*,i-1] = tt
       oplot, ltt, tt, psym=3
       mtt += tt
   endfor
   mtt /= (nsims-fsim+1)
   mdtt = mtt*0.
   for i=0l,n_elements(tt)-1 do mdtt[i] = median(ttar[i,fsim-1:*])
   oplot, ltt, mtt, col=210, thick=2
   oplot, ltt, mdtt, col=90, thick=2, line=2
   xyouts, 1250, 5500, '!6Planck 2013', col=245, chars=1.5
   xyouts, 1250, 5000, '!6143 ffp7 mean', col=210, chars=1.5
   xyouts, 1250, 4500, '!6143 ffp7 median', col=90, chars=1.5

;33 ------- EE
   ee = extract_xfaster_newdat( file, lcen=lee, binning='/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2', btcl=bteecl, ncl=2, tclfile=ave_file )
   eear = fltarr(n_elements(ee),nsims)
   eear[*,0] = ee
   mee = ee
   plot, lee, ee, chars=1.5, ytit='!8D!dl!uEE!n', psym=3, yr=[-5.5,50.5], ys=1, xr=[0,3000], position=[0.3875,0.45,0.6525,0.975], xtickname=strarr(7)+' '
   oplot, l, tee, col=245, thick=2

   for i=fsim+1,nsims do begin
       file = dir+root+'_'+strtrim(string(i),2)+'.newdat'
       ee = extract_xfaster_newdat( file, ncl=2 )
       eear[*,i-1] = ee
       oplot, lee, ee, psym=3
       mee += ee
   endfor
   mee /= (nsims-fsim+1)
   mdee = mee*0.
   for i=0,n_elements(ee)-1 do mdee[i] = median(eear[i,fsim-1:*])
   oplot, lee, mee, col=210, thick=2
   oplot, lee, mdee, col=90, thick=2, line=2

;## -------- TE
   te = extract_xfaster_newdat( file, lcen=lte, binning='/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2', btcl=bttecl, ncl=4, tclfile=ave_file )
   mte = te
   tear = fltarr(n_elements(te),nsims)
   tear[*,0] = te
   plot, lte, te, chars=1.5, ytit='!8D!dl!uTE!n', psym=3, yr=[-155,155], ys=1, position=[0.715,0.45,0.98,0.975], xtickname=strarr(7)+' '
   oplot, l, tte, col=245, thick=2

   for i=fsim+1,nsims do begin
       file = dir+root+'_'+strtrim(string(i),2)+'.newdat'
       te = extract_xfaster_newdat( file, ncl=4 )
       oplot, lte, te, psym=3
       tear[*,i-1] = te
       mte += te
   endfor
   mte /= (nsims-fsim+1)
   mdte = mte*0.
   for i=0,n_elements(te)-1 do mdte[i] = median(tear[i,fsim-1:*])
   oplot, lte, mte, col=210, thick=2
   oplot, lte, mdte, col=90, thick=2, line=2


;## ------ Residuals

;##   plot, ltt, (mtt-btttcl)/btttcl, chars=1.5, yr=[-0.055,0.055], position=[0.06,0.075,0.325,0.44], ys=1, psym=-6, ytit='!6TT: <!8D!dl!n!6>-!8D!dl!uTh!n'
;##   oplot, ltt, mtt*0, line=2
;##   oplot, ltt, (mdtt-btttcl)/btttcl, col=90, psym=-5
;   pwf = healpixwindow(2048)
;   pwf = gaussbeam(0.25,4096)
;   bpwftt = bp_binning(pwf[*,0],'/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TT')
;   bpwfee = bp_binning(pwf[*,0],'/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_EE')
;   bpwfte = bp_binning(pwf[*,0],'/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const2_TE')
;##   bpwf[*] = 1.

   plot, ltt, (mtt)/btttcl, chars=1.5, yr=[0.945,1.055], position=[0.06,0.075,0.325,0.44], ys=1, psym=-6, ytit='!6TT: <!8D!dl!n!6>-!8D!dl!uTh!n'
   oplot, ltt, mtt*0+1, line=2
   oplot, ltt, (mdtt)/btttcl, col=90, psym=-5
   xyouts, 1050, .95, '!6Mean-Bestfit', chars=1.5
   xyouts, 1050, .925, '!6Median-Bestfit', col=90, chars=1.5
   xyouts, 1050, .905, '!6PWF (not squared!!)', col=210, chars=1.5
;   oplot, ltt, 1./bpwftt^2, col=210

   plot, lee, (mee)/bteecl, chars=1.5, yr=[0.945,1.055], position=[0.3875,0.075,0.6525,0.44], xr=[0,3000], ys=1, psym=-6, ytit='!6EE: <!8D!dl!n!6>-!8D!dl!uTh!n'
   oplot, lee, mee*0+1, line=2
   oplot, lee, (mdee)/bteecl, col=90, psym=-5
;   oplot, ltt, 1./bpwfee^2, col=210

   plot, lte, (mte)/bttecl, chars=1.5, yr=[0.85,1.15], position=[0.715,0.075,0.98,0.44], ys=1, psym=-6, ytit='!6TE: <!8D!dl!n!6>-!8D!dl!uTh!n'
   oplot, lte, mte*0+1, line=2
   oplot, lte, (mdte)/bttecl, col=90, psym=-5
;   oplot, ltt, 1./bpwfte^2, col=210


end
