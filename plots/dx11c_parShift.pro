
   True = 1b
   False = 0b

   do_marge = True
   
   mollview, findgen(12), px=450, win=0
   loadct, 39
   !p.color=0
   !p.background=255

;## ------
   A = FIndGen(16) * (!PI*2/16.) 
   UserSym, cos(A), sin(A), /fill 
;## ---

   xfdir = '/global/scratch2/sd/dpietrob/Software/XFaster/xfaster'
   outdir = xfdir+'/../plots/'
   filename = outdir+'dx11c_parShift_v3.png'
   
   roots = [ $
    'xf_chains_217_undust_TP_const_dl50_fisherWins_l1000_Gfg_m30_OK', $
    'xf_chains_143_undust_TP_const_dl50_fisherWins_l1000_Gfg_m30_OK', $
    'xf_chains_217_undust_TP_const_dl50_fisherWins_l2000_Gfg_m30_OK', $
    'xf_chains_143_undust_TP_const_dl50_fisherWins_l2000_Gfg_m30_OK', $
    'xf_chains_217_undust_TP_const_dl50_fisherWins_l1000_Gfg_m70_OK', $
    'xf_chains_143_undust_TP_const_dl50_fisherWins_l1000_Gfg_m70_OK', $
    'xf_chains_217_undust_TP_const_dl50_fisherWins_l2000_Gfg_m70_OK', $
    'xf_chains_143_undust_TP_const_dl50_fisherWins_l2000_Gfg_m70_OK', $
    'xf_chains_217_undust_T_const_dl50_fisherWins_l1000_Gfg_m70_OK', $
    'xf_chains_143_undust_T_const_dl50_fisherWins_l1000_Gfg_m70_OK', $
    'xf_chains_217_undust_T_const_dl50_fisherWins_l2000_Gfg_m70_OK', $
    'xf_chains_143_undust_T_const_dl50_fisherWins_l2000_Gfg_m70_OK', $
    'xf_chains_217_undust_T_const_dl50_fisherWins_l2000_Gfg_m70_OK_notSZ' $
    ]

   tags = ['217', $
           '143', $
           '217', $
           '143' $
   ]

   bin = 'dl50' ;'ctp' ; 'dl50'

   if bin eq 'dl50' then binning = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const_dl50'

   nroots = n_elements(roots)
   pname = ['omegabh2', 'omegach2', 'theta', 'ns', 'H0', 'clamp', 'atSZ','akSZ','apsTT', 'acib', 'apsTE', 'apsEE', 'nCIB', 'agalTT', 'agalTE', 'agalEE', 'ngalTT', 'ngalTE', 'ngalEE' ]
   pindx = [0,1,2,4,23,28,[lindgen(13)+6]]
   tindx = [0,1,2,4,17,22,[lindgen(7)+6]]
   tipos = [0,1,2,3,4,5,6,7,8,9,12,13,16]

; no-tSZ
   stindx = [0,1,2,4,16,21,[lindgen(6)+6]]
   stipos = [0,1,2,3,4,5,7,8,9,12,13,16]

   npars = n_elements(pindx)

;   i_ombh2 = 0
;   i_omch2 = 1
;   i_theta = 2
;   i_ns = 3
;   i_H0 = 4
;   i_amp = 5
;   iplt = [i_ombh2,i_omch2,i_theta,i_ns,i_H0,i_amp]
   iplt = lindgen(npars)

   prange = [ $
              [0.020,0.0245], $
              [0.095,0.14], $
              [1.0345,1.047], $
              [0.9,1.05], $
              [55,80], $
              [1.775,1.92], $
              [0,12], $ ;atsk
              [0,30], $ ;aksz
              [0,100], $ ;apstt
              [0,40], $ ;acib
              [0,35], $ ;apste
              [0,16], $ ;apsee
              [0.45,1.25], $ ;ncib
              [0,30], $ ;agaltt
              [0,20], $ ;agalte
              [0,25], $ ;agalee
              [0,4.5], $ ;ngaltt
              [0,4.5], $ ;ngalte
              [0,4.5] $ ;ngalee
              ]

   pars = fltarr(nroots,npars,4)
   pars[*,*,*] = 0./0.
   lpars = fltarr(nroots,npars,4)
   i_peak = 0
   i_low = 2
   i_up = 3
   i_std = 1

   for i=0,nroots-1 do begin
       if do_marge then begin
       readcol, xfdir+'/results/'+roots[i]+'.margestats',name,peak,std,low1,up1,format='a,f,f,f,f', skipline=3
       ;##readcol, xfdir+'/results/'+roots[i]+'.margestats',name,peak,std,up1,low1,format='a,f,f,x,f,x,f', skipline=3
       if i lt 8 then begin
           pars[i,*,i_peak] = peak[pindx]
           pars[i,*,i_std]  = std[pindx]
           pars[i,*,i_low]  = low1[pindx]
           pars[i,*,i_up]   = up1[pindx]
       endif
       if i ge 8 and i lt 12 then begin
           pars[i,tipos,i_peak] = peak[tindx]
           pars[i,tipos,i_std]  = std[tindx]
           pars[i,tipos,i_low]  = low1[tindx]
           pars[i,tipos,i_up]   = up1[tindx]
       endif
       if i eq 12 then begin
           pars[i,stipos,i_peak] = peak[stindx]
           pars[i,stipos,i_std]  = std[stindx]
           pars[i,stipos,i_low]  = low1[stindx]
           pars[i,stipos,i_up]   = up1[stindx]
       endif
   endif else begin
       readcol, xfdir+'/results/'+roots[i]+'.likestats',name,peak,up1,low1,format='a,f,x,f,f', skipline=3
       if i lt 8 then begin
           pars[i,*,i_peak] = peak[pindx]
;           pars[i,*,i_std]  = std[pindx]
           pars[i,*,i_low]  = low1[pindx]
           pars[i,*,i_up]   = up1[pindx]
       endif
       if i ge 8 and i lt 12 then begin
           pars[i,tipos,i_peak] = peak[tindx]
;           pars[i,tipos,i_std]  = std[tindx]
           pars[i,tipos,i_low]  = low1[tindx]
           pars[i,tipos,i_up]   = up1[tindx]
       endif
       if i eq 12 then begin
           pars[i,stipos,i_peak] = peak[stindx]
;           pars[i,stipos,i_std]  = std[stindx]
           pars[i,stipos,i_low]  = low1[stindx]
           pars[i,stipos,i_up]   = up1[stindx]
       endif
endelse
;
;##       readcol, xfdir+'/results/'+roots[i]+'.likestats',name,peak,low1,up1,format='a,f,f,f', skipline=1
;##       lpars[i,pindx,i_peak] = peak[pindx]
;##       lpars[i,pindx,i_low] = low1[pindx]
;##       lpars[i,pindx,i_up] = up1[pindx]
;
;stop
   endfor

   !p.multi=[0,ceil(sqrt(npars)),floor(sqrt(npars))]
   window, 0, xsize=1200, ysize=800
;##   plot, lindgen(nroots)+1, pars[*,i_ombh2,i_peak], xtickname=[' ',tags],xticks=nroots, psym=8, xticklayout=4, yr=[0.020,0.0235], ys=1, chars=1.5

   xv31 = [1,2]
   xv32 = [3,4]
   xv71 = [1.1,2.1]
   xv72 = [3.1,4.1]
   xvt  = [1.2,2.2]
   xvt2  = [3.2,4.2]

;   ypos = [0.0198,0.103,1.034,0.89,54,1.77]

   for i=0,npars-1 do begin
       print, i, pname[i]
       if total( finite( pars[0:1,iplt[i],i_peak], /nan ) ) eq 0 then $
         plot, xv31, pars[0:1,iplt[i],i_peak], xtickname=[' ',tags,' '], chars=2.5, xr=[0,5], xs=1, yr=prange[*,i], ys=1, ytit=pname[i], psym=8, xticks=5, xmargin=[8,1], ymargin=[2,1]
       oploterr, xv31, low=pars[0:1,iplt[i],i_low], up=pars[0:1,iplt[i],i_up]

       if total( finite( pars[4:5,iplt[i],i_peak], /nan ) ) eq 0 then $
         oplot, xv71, pars[4:5,iplt[i],i_peak], psym=8, col=70
       oploterr, xv71, low=pars[4:5,iplt[i],i_low], up=pars[4:5,iplt[i],i_up], col=70

       if total( finite( pars[2:3,iplt[i],i_peak], /nan ) ) eq 0 then $
         oplot, xv32, pars[2:3,iplt[i],i_peak], psym=6, thick=2
       oploterr, xv32, low=pars[2:3,iplt[i],i_low], up=pars[2:3,iplt[i],i_up]

       if total( finite( pars[6:7,iplt[i],i_peak], /nan ) ) eq 0 then $
         oplot, xv72, pars[6:7,iplt[i],i_peak], psym=6, col=70, thick=2
       oploterr, xv72, low=pars[6:7,iplt[i],i_low], up=pars[6:7,iplt[i],i_up], col=70

       if total( finite( pars[8:9,iplt[i],i_peak], /nan ) ) eq 0 then $
         oplot, xvt, pars[8:9,iplt[i],i_peak], psym=8, col=105, thick=2
       oploterr, xvt, low=pars[8:9,iplt[i],i_low], up=pars[8:9,iplt[i],i_up], col=105

       if total( finite( pars[10:11,iplt[i],i_peak], /nan ) ) eq 0 then $
         oplot, xvt2, pars[10:11,iplt[i],i_peak], psym=6, col=105, thick=2
       oploterr, xvt2, low=pars[10:11,iplt[i],i_low], up=pars[10:11,iplt[i],i_up], col=105

       if total( finite( pars[12,iplt[i],i_peak], /nan ) ) eq 0 then $
         oplot, [3.3,3.3], [pars[12,iplt[i],i_peak],pars[12,iplt[i],i_peak]], psym=6, col=165, thick=2
       oploterr, [3.3], low=pars[12,iplt[i],i_low], up=pars[12,iplt[i],i_up], col=165

   endfor
   plot, /nodata, [0,1], [0,1], col=255
   legend, ['m30 TP','m70 TP','m70 T','m70 T no-tSZ'], psym=[8,8,8,8], col=[0,70,105,165], /top, /right, chars=1.25
   legend, ['l=1000','l=2000'], psym=[8,6], thick=[1,2], chars=1.25
;   oplot, (lindgen(nroots)+1.05), lpars[*,i_ombh2,i_peak], psym=8, col=70
;   oploterr, (lindgen(nroots)+1.05), lpars[*,i_ombh2,i_peak], low=lpars[*,i_ombh2,i_low], up=lpars[*,i_ombh2,i_up], col=70

   write_png, filename, tvrd(/true)

stop
;##   freq = ['100','143','recal_217']
   freq = ['100','143','recal_217']
   nwdfu = xfdir+'dx11c_SDR_yr_1-2_IQUmap_'+freq+'_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const_dl50_fisherWins_xfcl_l3000_combined_70_x_combined_70_l4000.newdat'
;## ------
;##   nwdf = xfdir+'dx11c_SDR_yr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs_const_dl50_fisherWins_xfcl_l3000_Tmask_badpixmasked_x_Pmask_badpixmasked_l4000.newdat'

   bffreq = ['100','143','217']
   bffile = xfdir+'xf_chains_'+bffreq+'_undust_TP_const_dl50_fisherWins_l2000_Gfg_m70_OK_bestfit.txt'
   fgbffile = xfdir+'xf_chains_'+bffreq+'_undust_TP_const_dl50_fisherWins_l2000_Gfg_m70_OK_fg_bestfit.txt'
   cmbbffile = xfdir+'xf_chains_'+bffreq+'_undust_TP_const_dl50_fisherWins_l2000_Gfg_m70_OK_lensedCls.dat'

   ncl = [1,2,4]
   cltag = ['_TT','_EE','_TE']

   readcol, cmbbffile[2], bfl, bftt, bfee, xx, bfte
   bfll = bfl*(bfl+1)/2./!pi
   bftt = bftt/bfll
   bfee = bfee/bfll
   bfte = bfte/bfll

;##   lout = fltarr(100)
   ttout = fltarr(3,100)
   eeout = fltarr(3,100)
   teout = fltarr(3,100)
   tterout = fltarr(3,100)
   eeerout = fltarr(3,100)
   teerout = fltarr(3,100)

   for i=0,2 do begin

       if i eq 0 then btcl = xf_binning(reform([0,1,bftt]),binning+cltag[i])
       if i eq 1 then btcl = xf_binning(reform([0,1,bfee]),binning+cltag[i])
       if i eq 2 then btcl = xf_binning(reform([0,1,bfte]),binning+cltag[i])

       if bin eq 'dl50' then begin
;##           f1 = nwdf[1]
           f1u = nwdfu[1]
;##           bcls143 = extract_xfaster_newdat(f1, lcen=lc, ncl=ncl[i])
           bcls143u = extract_xfaster_newdat(f1u, ncl=ncl[i], cler=er1, lcen=lc)

           lout = lc

           readcol, fgbffile[1], l, tt, ee, te
           tt = reform([0,0,tt])
           l = reform([0,1,l])
           ll = l*(l+1.)/2./!pi
           if i eq 0 then bfgtt1 = xf_binning( tt/ll, binning+cltag[i])
           if i eq 1 then bfgtt1 = xf_binning( ee/ll, binning+cltag[i])
           if i eq 2 then bfgtt1 = xf_binning( te/ll, binning+cltag[i])

;##           f2 = nwdf[2]
           f2u = nwdfu[2]
;##           bcls217 = extract_xfaster_newdat(f2, lcen=lc, binning=binning, ncl=ncl[i] )
           bcls217u = extract_xfaster_newdat(f2u, ncl=ncl[i], cler=er2)
           readcol, fgbffile[2], l, tt, ee, te
           if i eq 0 then bfgtt2 = xf_binning( tt/ll, binning+cltag[i])
           if i eq 1 then bfgtt2 = xf_binning( ee/ll, binning+cltag[i])
           if i eq 2 then bfgtt2 = xf_binning( te/ll, binning+cltag[i])
           
;##           f3 = nwdf[0]
           f3u = nwdfu[0]
;##           bcls100 = extract_xfaster_newdat(f3, lcen=lc, binning=binning, ncl=ncl[i])
           bcls100u = extract_xfaster_newdat(f3u, ncl=ncl[i], cler=er3)
           readcol, fgbffile[0], l, tt, ee, te
           if i eq 0 then bfgtt3 = xf_binning( tt/ll, binning+cltag[i])
           if i eq 1 then bfgtt3 = xf_binning( ee/ll, binning+cltag[i])
           if i eq 2 then bfgtt3 = xf_binning( te/ll, binning+cltag[i])

       endif

       !p.multi=[0,1,2]
       window,i,ysize=900, xsize=900
       if i eq 0 then th = bftt*bfll
       if i eq 1 then th = bfee*bfll
       if i eq 2 then th = bfte*bfll
       if i lt 2 then begin
           if i eq 0 then plot, bfl, th, chars=1.75, ytit='!8D!dl!6!uTT!n [!7l!6K!u2!n]', line=2, position=[0.125,0.35,0.95,0.975], xtickv=findgen(6)*500, xticks=5, xtickname=strarr(6)+' ', xr=[0,2500], xs=1
           if i eq 1 then plot, bfl, th, chars=1.75, ytit='!8D!dl!6!uEE!n [!7l!6K!u2!n]', line=2, position=[0.125,0.35,0.95,0.975], xtickv=findgen(4)*500, xticks=4, xtickname=strarr(5)+' ', xr=[0,2000], xs=1
           oplot, bfl, th, col=210
;   oplot, lc, bcls143u, thick=2, psym=-4, col=245
;   oplot, lc, bcls217u, thick=2, psym=-4
;   oplot, lc, bcls100u, thick=2, psym=-4, col=70
           oplot, lc, (bcls143u-bfgtt1), psym=8, thick=3, col=245
           oplot, lc*1.01, (bcls217u-bfgtt2), psym=8, thick=3
           oplot, lc*0.99, (bcls100u-bfgtt3), psym=8, thick=3, col=70
;
           for x=0,n_elements(lc)-1 do begin
               oplot, [lc[x],lc[x]],[(bcls143u[x]-bfgtt1[x]-er1[x]),(bcls143u[x]-bfgtt1[x]+er1[x])], col=245
               oplot, [lc[x],lc[x]]*1.01,[(bcls217u[x]-bfgtt2[x]-er2[x]),(bcls217u[x]-bfgtt2[x]+er2[x])]
               oplot, [lc[x],lc[x]]*0.99,[(bcls100u[x]-bfgtt3[x]-er3[x]),(bcls100u[x]-bfgtt3[x]+er3[x])], col=70
           endfor
;##           errplot, lc, (bcls143u-bfgtt1-er1), (bcls143u-bfgtt1+er1), psym=0, col=245
;##           errplot, lc, (bcls217u-bfgtt2-er2), (bcls217u-bfgtt2+er2), thick=0
;##           errplot, lc, (bcls100u-bfgtt3-er3), (bcls100u-bfgtt3+er3), thick=0, col=70
       endif else begin
; ------
           plot, bfl, th/bfl, chars=1.75, ytit='!8D!dl!6!uTE!8!n/l!6 [!7l!6K!u2!n]', line=2, position=[0.125,0.35,0.95,0.975], xtickv=findgen(6)*500, xticks=5, xtickname=strarr(6)+' ', xr=[0,2500], xs=1, yr=[-0.5, 0.5], ys=1
           oplot, bfl, th/bfl, col=210
;       oplot, lc, bcls143u/lc, thick=2, psym=-4, col=245
;       oplot, lc, bcls217u/lc, thick=2, psym=-4
;       oplot, lc, bcls100u/lc, thick=2, psym=-4, col=70
           oplot, lc, (bcls143u-bfgtt1)/lc, psym=8, thick=3, col=245
           oplot, lc*1.01, (bcls217u-bfgtt2)/lc, psym=8, thick=3
           oplot, lc*0.99, (bcls100u-bfgtt3)/lc, psym=8, thick=3, col=70
;
           for x=0,n_elements(lc)-1 do begin
               oplot, [lc[x],lc[x]], [(bcls143u[x]-bfgtt1[x]-er1[x]),(bcls143u[x]-bfgtt1[x]+er1[x])]/lc[x], col=245
               oplot, [lc[x],lc[x]]*1.01, [(bcls217u[x]-bfgtt2[x]-er2[x]),(bcls217u[x]-bfgtt2[x]+er2[x])]/lc[x]
               oplot, [lc[x],lc[x]]*0.99, [(bcls100u[x]-bfgtt3[x]-er3[x]),(bcls100u[x]-bfgtt3[x]+er3[x])]/lc[x], col=70
           endfor
;##           errplot, lc, (bcls143u-bfgtt1-er1)/lc, (bcls143u-bfgtt1+er1)/lc, col=245
;##           errplot, lc, (bcls217u-bfgtt2-er2)/lc, (bcls217u-bfgtt2+er2)/lc
;##           errplot, lc, (bcls100u-bfgtt3-er3)/lc, (bcls100u-bfgtt3+er3)/lc, col=70

       endelse

       if i eq 0 then begin
           ttout[1,0:n_elements(lc)-1] = (bcls143u-bfgtt1)
           ttout[2,0:n_elements(lc)-1] = (bcls217u-bfgtt2)
           ttout[0,0:n_elements(lc)-1] = (bcls100u-bfgtt3)
           tterout[1,0:n_elements(lc)-1] = er1
           tterout[2,0:n_elements(lc)-1] = er2
           tterout[0,0:n_elements(lc)-1] = er3
       endif

       if i eq 1 then begin
           eeout[1,0:n_elements(lc)-1] = (bcls143u-bfgtt1)
           eeout[2,0:n_elements(lc)-1] = (bcls217u-bfgtt2)
           eeout[0,0:n_elements(lc)-1] = (bcls100u-bfgtt3)
           eeerout[1,0:n_elements(lc)-1] = er1
           eeerout[2,0:n_elements(lc)-1] = er2
           eeerout[0,0:n_elements(lc)-1] = er3
       endif

       if i eq 2 then begin
           teout[1,0:n_elements(lc)-1] = (bcls143u-bfgtt1)
           teout[2,0:n_elements(lc)-1] = (bcls217u-bfgtt2)
           teout[0,0:n_elements(lc)-1] = (bcls100u-bfgtt3)
           teerout[1,0:n_elements(lc)-1] = er1
           teerout[2,0:n_elements(lc)-1] = er2
           teerout[0,0:n_elements(lc)-1] = er3
       endif

       if i ne 1 then legend, ['XF217','XF143','XF100','217-best fit'],col=[0,245,70,210],psym=[8,8,8,3], thick=3, chars=1.5, /right
       if i eq 1 then legend, ['XF217','XF143','XF100','217-best fit'],col=[0,245,70,210],psym=[8,8,8,3], thick=3, chars=1.5
;##       md = ( (bcls143u-bfgtt1) + (bcls217u-bfgtt2) + (bcls100u-bfgtt3) ) / 3
       md = btcl 
;##       tmp = bcls143u-bfgtt1-md
       if i eq 0 then mx = 60
       if i eq 1 then mx = 10
       if i eq 2 then mx = 0.04
; ------ Second panel
       if i lt 2 then begin
           if i eq 0 then plot, lc, btcl*0., yr=[-mx,mx], line=2, chars=1.75, xtit='!8l!6', ytit='!7D!8D!dl!6!uTT!n [!7l!6K!u2!n]',position=[0.125,0.1,0.95,0.35], xtickv=findgen(6)*500, xticks=5, xr=[0,2500], xs=1
           if i eq 1 then plot, lc, btcl*0., yr=[-mx,mx], line=2, chars=1.75, xtit='!8l!6', ytit='!7D!8D!dl!6!uEE!n [!7l!6K!u2!n]',position=[0.125,0.1,0.95,0.35], xtickv=findgen(5)*500, xticks=4, xr=[0,2000], xs=1
           oplot, lc, btcl*0., line=2, col=210
           oplot, lc, bcls143u-bfgtt1-md, psym=8, thick=3, col=245
           oplot, lc*1.01, bcls217u-bfgtt2-md, psym=8, thick=3
           oplot, lc*.99, bcls100u-bfgtt3-md, psym=8, thick=3, col=70
           for x=0,n_elements(lc)-1 do begin
               oplot, reform([lc[x],lc[x]]), reform([(bcls143u[x]-bfgtt1[x]-md[x]-er1[x]),(bcls143u[x]-bfgtt1[x]-md[x]+er1[x])]), col=245
               oplot, reform([lc[x],lc[x]])*1.01,reform([(bcls217u[x]-bfgtt2[x]-md[x]-er2[x]),(bcls217u[x]-bfgtt2[x]-md[x]+er2[x])])
               oplot, reform([lc[x],lc[x]])*0.99,reform([(bcls100u[x]-bfgtt3[x]-md[x]-er3[x]),(bcls100u[x]-bfgtt3[x]-md[x]+er3[x])]), col=70
           endfor

;##           errplot, lc, bcls143u-bfgtt1-md-er1, bcls143u-bfgtt1-md+er1, col=245
;##           errplot, lc*1.01, bcls217u-bfgtt2-md-er2, bcls217u-bfgtt2-md+er2
;##           errplot, lc*.99, bcls100u-bfgtt3-md-er3, bcls100u-bfgtt3-md+er3, col=70
       endif else begin
           plot, lc, btcl*0., yr=[-mx,mx], line=2, chars=1.75, xtit='!8l!6', ytit='!7D!8D!dl!6!uTE!8!n/l!6 [!7l!6K!u2!n]',position=[0.125,0.1,0.95,0.35], xtickv=findgen(6)*500, xticks=5, xr=[0,2500], xs=1
           oplot, lc, btcl*0., line=2, col=210
           oplot, lc, (bcls143u-bfgtt1-md)/lc, psym=8, thick=3, col=245
           oplot, lc*1.01, (bcls217u-bfgtt2-md)/lc, psym=8, thick=3
           oplot, lc*.99, (bcls100u-bfgtt3-md)/lc, psym=8, thick=3, col=70
           for x=0,n_elements(lc)-1 do begin
               oplot, [lc[x],lc[x]],[(bcls143u[x]-bfgtt1[x]-md[x]-er1[x]),(bcls143u[x]-bfgtt1[x]-md[x]+er1[x])]/lc[x], col=245
               oplot, [lc[x],lc[x]]*1.01,[(bcls217u[x]-bfgtt2[x]-md[x]-er2[x]),(bcls217u[x]-bfgtt2[x]-md[x]+er2[x])]/lc[x]
               oplot, [lc[x],lc[x]]*0.99,[(bcls100u[x]-bfgtt3[x]-md[x]-er3[x]),(bcls100u[x]-bfgtt3[x]-md[x]+er3[x])]/lc[x], col=70
           endfor
;##           errplot, lc, (bcls143u-bfgtt1-md-er1)/lc, (bcls143u-bfgtt1-md+er1)/lc, col=245
;##           errplot, lc*1.01, (bcls217u-bfgtt2-md-er2)/lc, (bcls217u-bfgtt2-md+er2)/lc
;##           errplot, lc*.99, (bcls100u-bfgtt3-md-er3)/lc, (bcls100u-bfgtt3-md+er3)/lc, col=70
       endelse
;##   oplot, lc, btcl*0.025-100, col=245 
;   legend, ['XF143-Fg-P2013','XF217-Fg-P2013'], col=[0,70,245], chars=1.25, psym=[4,4,4], thick=2, /right

       write_png,'xf_dx11c_subEGFG_'+bin+cltag[i]+'.png', tvrd(/true)

;       stop
   endfor

   close,/all
   for i=0,2 do begin
       openw, 1, 'xf_dx11c_subEGFG_'+bin+'_'+bffreq[i]+'.dat'
       printf,1,'# l cl_tt er_tt, cl_ee er_ee cl_te er_te'
       for l=0,n_elements(lc)-1 do printf,1,lc[l],ttout[i,l],tterout[i,l],eeout[i,l],eeerout[i,l],teout[i,l],teerout[i,l], format='(1i12,6f15.8)'
       close,1
;##       spawn, 'head xf_dx11c_subEGFG_'+bin+'_'+bffreq[i]+'.dat'
;##       spawn, 'tail xf_dx11c_subEGFG_'+bin+'_'+bffreq[i]+'.dat'
   endfor

   readcol, '../data/planck_base_lensedCls.dat', a,b,c,x,d
   window, 3, xsize=900, ysize=900
   plot, bfl, bftt*bfll, chars=2, ytit='!8D!dl!u!6TT!n [!7l!6K!u2!n]', position=[0.125,0.35,0.95,0.975], xtickv=findgen(6)*500, xticks=5, xtickname=strarr(6)+' ', xr=[0,2500], xs=1
   oplot, bfl, bftt*bfll, col=70, thick=2
   oplot, a, b, col=245, thick=2
   legend, ['Bestfit 217','Planck 2013'], col=[70,245], line=[0,0], thick=2, /right, chars=1.5

   plot, a, bftt*ll/b, yr=[0.9, 1.1], chars=1.75, xtit='!8l!6', ytit='!8D!dl!6!uTT,217!n / !8D!dl!6!uTT,2013!n',position=[0.125,0.1,0.95,0.35], xtickv=findgen(6)*500, xticks=5, xr=[0,2500], xs=1
   oplot, a, a*0.+1, col=245, line=2
   write_png, 'xf_dx11c_theory_TT.png', tvrd(/true)

   window, 4, xsize=900, ysize=900
   plot, bfl, bfee*bfll, chars=2, ytit='!8D!dl!u!6EE!n [!7l!6K!u2!n]', position=[0.125,0.35,0.95,0.975], xtickv=findgen(6)*500, xticks=5, xtickname=strarr(6)+' ', xr=[0,2500], xs=1
   oplot, bfl, bfee*bfll, col=70, thick=2
   oplot, a, c, col=245, thick=2
   legend, ['Bestfit 217','Planck 2013'], col=[70,245], line=[0,0], thick=2, /right, chars=1.5

   plot, a, bfee*ll/c, yr=[0.9, 1.1], chars=1.75, xtit='!8l!6', ytit='!8D!dl!6!uEE,217!n / !8D!dl!6!uEE,2013!n',position=[0.125,0.1,0.95,0.35], xtickv=findgen(6)*500, xticks=5, xr=[0,2500], xs=1
   oplot, a, a*0.+1, col=245, line=2
   write_png, 'xf_dx11c_theory_EE.png', tvrd(/true)
   
   window, 5, xsize=900, ysize=900
   plot, bfl, bfte*bfll, chars=2, ytit='!8D!dl!u!6TE!n [!7l!6K!u2!n]', position=[0.125,0.35,0.95,0.975], xtickv=findgen(6)*500, xticks=5, xtickname=strarr(6)+' ', xr=[0,2500], xs=1
   oplot, bfl, bfte*bfll, col=70, thick=2
   oplot, a, d, col=245, thick=2
   legend, ['Bestfit 217','Planck 2013'], col=[70,245], line=[0,0], thick=2, /right, chars=1.5

   plot, a, bftt*ll/b, yr=[0.9, 1.1], chars=1.75, xtit='!8l!6', ytit='!8D!dl!6!uTT,217!n / !8D!dl!6!uTT,2013!n',position=[0.125,0.1,0.95,0.35], xtickv=findgen(6)*500, xticks=5, xr=[0,2500], xs=1
   oplot, a, a*0.+1, col=245, line=2
   write_png, 'xf_dx11c_theory_TE.png', tvrd(/true)



end
