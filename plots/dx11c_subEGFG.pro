; python sequence to get bestfit:
; import get_bestfit_pars as g
; ### root='143_raw_T_const_dl50_fisherWins_l2000_Gfg'
; root='143_undust_TP_const_dl50_fisherWins_l2000_Gfg_m70'
; dir='/global/scratch2/sd/dpietrob/Software/XFaster/xfaster'
; pars=g.read_pars('',filename='results/xf_chains_'+root+'.likestats',dir=dir)
; g.write_pars(pars,root='',dir=dir,filetag='_'+root)
; g.run_camb('',cambfile='xf_camb_'+root+'.ini', dir=dir) 
; g.compute_bestfit('',dir=dir,cambparfile='xf_camb_'+root+'.ini',clsfile='xf_chains_'+root+'_lensedCls.dat',filetag='xf_chains_'+root+'_')
  
   True = 1b
   False = 0b
   
   do_100 = False

   mollview, findgen(12), px=450, win=0
   loadct, 39
   !p.color=0
   !p.background=255

;## ------
   A = FIndGen(16) * (!PI*2/16.) 
   UserSym, cos(A), sin(A), /fill 
;## ---

   xfdir = '/global/scratch2/sd/dpietrob/Software/XFaster/xfaster/'
   bin = 'dl50' ;'ctp' ; 'dl50'

   if bin eq 'dl50' then binning = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/const/const_dl50'
   if bin eq 'ctp' then binning = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin'

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
       ; bfm for each channel
       ; 100
       readcol, cmbbffile[0], bfl, bftt, bfee, xx, bfte
       bfll = bfl*(bfl+1)/2./!pi
       bftt = bftt/bfll
       bfee = bfee/bfll
       bfte = bfte/bfll
       if i eq 0 then btcl3 = xf_binning(reform([0,1,bftt]),binning+cltag[i])
       if i eq 1 then btcl3 = xf_binning(reform([0,1,bfee]),binning+cltag[i])
       if i eq 2 then btcl3 = xf_binning(reform([0,1,bfte]),binning+cltag[i])
       ; 143
       readcol, cmbbffile[1], bfl, bftt, bfee, xx, bfte
       bfll = bfl*(bfl+1)/2./!pi
       bftt = bftt/bfll
       bfee = bfee/bfll
       bfte = bfte/bfll
       if i eq 0 then btcl1 = xf_binning(reform([0,1,bftt]),binning+cltag[i])
       if i eq 1 then btcl1 = xf_binning(reform([0,1,bfee]),binning+cltag[i])
       if i eq 2 then btcl1 = xf_binning(reform([0,1,bfte]),binning+cltag[i])
       ; 100
       readcol, cmbbffile[2], bfl, bftt, bfee, xx, bfte
       bfll = bfl*(bfl+1)/2./!pi
       bftt = bftt/bfll
       bfee = bfee/bfll
       bfte = bfte/bfll
       if i eq 0 then btcl2 = xf_binning(reform([0,1,bftt]),binning+cltag[i])
       if i eq 1 then btcl2 = xf_binning(reform([0,1,bfee]),binning+cltag[i])
       if i eq 2 then btcl2 = xf_binning(reform([0,1,bfte]),binning+cltag[i])

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
           if do_100 then oplot, lc*0.99, (bcls100u-bfgtt3), psym=8, thick=3, col=70
;
           for x=0,n_elements(lc)-1 do begin
               oplot, [lc[x],lc[x]],[(bcls143u[x]-bfgtt1[x]-er1[x]),(bcls143u[x]-bfgtt1[x]+er1[x])], col=245
               oplot, [lc[x],lc[x]]*1.01,[(bcls217u[x]-bfgtt2[x]-er2[x]),(bcls217u[x]-bfgtt2[x]+er2[x])]
               if do_100 then oplot, [lc[x],lc[x]]*0.99,[(bcls100u[x]-bfgtt3[x]-er3[x]),(bcls100u[x]-bfgtt3[x]+er3[x])], col=70
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
           if do_100 then oplot, lc*0.99, (bcls100u-bfgtt3)/lc, psym=8, thick=3, col=70
;
           for x=0,n_elements(lc)-1 do begin
               oplot, [lc[x],lc[x]], [(bcls143u[x]-bfgtt1[x]-er1[x]),(bcls143u[x]-bfgtt1[x]+er1[x])]/lc[x], col=245
               oplot, [lc[x],lc[x]]*1.01, [(bcls217u[x]-bfgtt2[x]-er2[x]),(bcls217u[x]-bfgtt2[x]+er2[x])]/lc[x]
               if do_100 then oplot, [lc[x],lc[x]]*0.99, [(bcls100u[x]-bfgtt3[x]-er3[x]),(bcls100u[x]-bfgtt3[x]+er3[x])]/lc[x], col=70
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
       if i eq 0 then mx = 40 ;60
       if i eq 1 then mx = 10
       if i eq 2 then mx = 0.02 ;0.04
; ------ Second panel
       if i lt 2 then begin
           if i eq 0 then plot, lc, btcl*0., yr=[-mx,mx], line=2, chars=1.75, xtit='!8l!6', ytit='!7D!8D!dl!6!uTT!n [!7l!6K!u2!n]',position=[0.125,0.1,0.95,0.35], xtickv=findgen(6)*500, xticks=5, xr=[0,2500], xs=1
           if i eq 1 then plot, lc, btcl*0., yr=[-mx,mx], line=2, chars=1.75, xtit='!8l!6', ytit='!7D!8D!dl!6!uEE!n [!7l!6K!u2!n]',position=[0.125,0.1,0.95,0.35], xtickv=findgen(5)*500, xticks=4, xr=[0,2000], xs=1
           oplot, lc, btcl*0., line=2, col=210
           if False then begin
               ; wrt 217 bfm
               oplot, lc, bcls143u-bfgtt1-md, psym=8, thick=3, col=245
               oplot, lc*1.01, bcls217u-bfgtt2-md, psym=8, thick=3
               if do_100 then oplot, lc*.99, bcls100u-bfgtt3-md, psym=8, thick=3, col=70
               for x=0,n_elements(lc)-1 do begin
                   oplot, reform([lc[x],lc[x]]), reform([(bcls143u[x]-bfgtt1[x]-md[x]-er1[x]),(bcls143u[x]-bfgtt1[x]-md[x]+er1[x])]), col=245
                   oplot, reform([lc[x],lc[x]])*1.01,reform([(bcls217u[x]-bfgtt2[x]-md[x]-er2[x]),(bcls217u[x]-bfgtt2[x]-md[x]+er2[x])])
                   if do_100 then oplot, reform([lc[x],lc[x]])*0.99,reform([(bcls100u[x]-bfgtt3[x]-md[x]-er3[x]),(bcls100u[x]-bfgtt3[x]-md[x]+er3[x])]), col=70
               endfor
           endif else begin
               oplot, lc, bcls143u-bfgtt1-btcl1, psym=8, thick=3, col=245
               oplot, lc*1.01, bcls217u-bfgtt2-btcl2, psym=8, thick=3
               if do_100 then oplot, lc*.99, bcls100u-bfgtt3-btcl3, psym=8, thick=3, col=70
               for x=0,n_elements(lc)-1 do begin
                   oplot, reform([lc[x],lc[x]]), reform([(bcls143u[x]-bfgtt1[x]-btcl1[x]-er1[x]),(bcls143u[x]-bfgtt1[x]-btcl1[x]+er1[x])]), col=245
                   oplot, reform([lc[x],lc[x]])*1.01,reform([(bcls217u[x]-bfgtt2[x]-btcl[x]-er2[x]),(bcls217u[x]-bfgtt2[x]-btcl2[x]+er2[x])])
                   if do_100 then oplot, reform([lc[x],lc[x]])*0.99,reform([(bcls100u[x]-bfgtt3[x]-btcl3[x]-er3[x]),(bcls100u[x]-bfgtt3[x]-btcl3[x]+er3[x])]), col=70
           endfor
           endelse

;##           errplot, lc, bcls143u-bfgtt1-md-er1, bcls143u-bfgtt1-md+er1, col=245
;##           errplot, lc*1.01, bcls217u-bfgtt2-md-er2, bcls217u-bfgtt2-md+er2
;##           errplot, lc*.99, bcls100u-bfgtt3-md-er3, bcls100u-bfgtt3-md+er3, col=70
       endif else begin
           plot, lc, btcl*0., yr=[-mx,mx], line=2, chars=1.75, xtit='!8l!6', ytit='!7D!8D!dl!6!uTE!8!n/l!6 [!7l!6K!u2!n]',position=[0.125,0.1,0.95,0.35], xtickv=findgen(6)*500, xticks=5, xr=[0,2500], xs=1
           oplot, lc, btcl*0., line=2, col=210
           if False then begin
               ; wrt 217 bfm
               oplot, lc, (bcls143u-bfgtt1-md)/lc, psym=8, thick=3, col=245
               oplot, lc*1.01, (bcls217u-bfgtt2-md)/lc, psym=8, thick=3
               if do_100 then oplot, lc*.99, (bcls100u-bfgtt3-md)/lc, psym=8, thick=3, col=70
               for x=0,n_elements(lc)-1 do begin
                   oplot, [lc[x],lc[x]],[(bcls143u[x]-bfgtt1[x]-md[x]-er1[x]),(bcls143u[x]-bfgtt1[x]-md[x]+er1[x])]/lc[x], col=245
                   oplot, [lc[x],lc[x]]*1.01,[(bcls217u[x]-bfgtt2[x]-md[x]-er2[x]),(bcls217u[x]-bfgtt2[x]-md[x]+er2[x])]/lc[x]
                   if do_100 then oplot, [lc[x],lc[x]]*0.99,[(bcls100u[x]-bfgtt3[x]-md[x]-er3[x]),(bcls100u[x]-bfgtt3[x]-md[x]+er3[x])]/lc[x], col=70
               endfor
           endif else begin
               oplot, lc, (bcls143u-bfgtt1-btcl1)/lc, psym=8, thick=3, col=245
               oplot, lc*1.01, (bcls217u-bfgtt2-btcl2)/lc, psym=8, thick=3
               if do_100 then oplot, lc*.99, (bcls100u-bfgtt3-btcl3)/lc, psym=8, thick=3, col=70
               for x=0,n_elements(lc)-1 do begin
                   oplot, [lc[x],lc[x]],[(bcls143u[x]-bfgtt1[x]-btcl1[x]-er1[x]),(bcls143u[x]-bfgtt1[x]-btcl1[x]+er1[x])]/lc[x], col=245
                   oplot, [lc[x],lc[x]]*1.01,[(bcls217u[x]-bfgtt2[x]-btcl2[x]-er2[x]),(bcls217u[x]-bfgtt2[x]-btcl2[x]+er2[x])]/lc[x]
                   if do_100 then oplot, [lc[x],lc[x]]*0.99,[(bcls100u[x]-bfgtt3[x]-btcl3[x]-er3[x]),(bcls100u[x]-bfgtt3[x]-btcl3[x]+er3[x])]/lc[x], col=70
               endfor
           endelse
;##           errplot, lc, (bcls143u-bfgtt1-md-er1)/lc, (bcls143u-bfgtt1-md+er1)/lc, col=245
;##           errplot, lc*1.01, (bcls217u-bfgtt2-md-er2)/lc, (bcls217u-bfgtt2-md+er2)/lc
;##           errplot, lc*.99, (bcls100u-bfgtt3-md-er3)/lc, (bcls100u-bfgtt3-md+er3)/lc, col=70
       endelse
;##   oplot, lc, btcl*0.025-100, col=245 
;   legend, ['XF143-Fg-P2013','XF217-Fg-P2013'], col=[0,70,245], chars=1.25, psym=[4,4,4], thick=2, /right

       write_png,'xf_dx11c_subEGFG_'+bin+cltag[i]+'_wrtSelfBFM.png', tvrd(/true)

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

   if False then begin
       ; BFM differences
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
   endif


end
