pro xf_post, init=init, lmax=lmax, win=win

    if keyword_set(init) then begin
        mollview, findgen(12l*4l^2), px=650, win=win
        loadct, 39
        !p.color=0
        !p.background=255
    endif
    if not keyword_set(lmax) then lmax=2000
    if not keyword_set(win) then win=1

    xfdir = '/global/scratch2/sd/dpietrob/Software/XFaster/'
    fits2cl, tcl, xfdir+'data/planck_lcdm_cl_uK_xf1.e-3.fits'
    tcl[*,4] = 0.
    tcl[*,5] = 0.
    l=findgen(n_elements(tcl[*,0]))
    ll=l*(l+1)/2./!pi
    il = findgen(lmax)+2

    if (0b) then begin
;  -issue EE- 
;    1) overplot 143 for full mission uncleaned and cleaned for
;    intermediate ell and 
    !p.multi = [0,1,2]
    window, win, xsize=720, ysize=450*1.5,tit='EE rise: preDX11 143 GHz'
    plot, l[il], tcl[il,1]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,lmax], yr=[0,60],tit='EE rise: preDX11 143 GHz'
; ---
    file='outputs/dx11_pre_Imap_143_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eeraw=extract_xfaster_newdata_output(file, lcen=lee, cler=eeer, res=eeres, ncl='EE') 
    oplot, lee,eeraw, psym=4
    errplot, lee,eeraw-eeer, eeraw+eeer
; ---
    file='outputs/dx11_pre_Imap_143_undusted_insideM_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eein=extract_xfaster_newdata_output(file, lcen=leein, cler=eeerin, res=eeresin, ncl='EE') 
    oplot, leein,eein, psym=6, col=70
    errplot, leein,eein-eeerin, eein+eeerin, col=70
; ---
    oplot, l[il], tcl[il,1]*ll[il], col=245
    xyouts, 50, 55, 'preDX11 143GHz raw', chars=1.5
    xyouts, 50, 50, 'preDX11 143GHz undusted-IN', chars=1.5, col=70
    plot, l[il], tcl[il,1]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,250]
; ---
    oplot, lee,eeraw, psym=4
    errplot, lee,eeraw-eeer, eeraw+eeer
; ---
    oplot, leein,eein, psym=6, col=70
    errplot, leein,eein-eeerin, eein+eeerin, col=70
; ---
    oplot, l[il], tcl[il,1]*ll[il], col=245

;    2) with hrhd and BB spectra undusted for all range of ell.
    !p.multi = [0,3,1]
    title = 'EE rise: diagnostic using BB spectrum as EE noise bias'
    window, win+2, xsize=720*1.9, ysize=450, tit=title
    plot, l[il], tcl[il,1]*ll[il], chars=2.5, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,lmax], yr=[0,70], ys=1;,tit=title
; ---
    file='outputs/dx11_pre_Imap_143_undusted_insideM_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eein=extract_xfaster_newdata_output(file, lcen=leein, cler=eeerin, res=eeresin, ncl='EE') 
    oplot, leein,eein, psym=6, col=70
    errplot, leein,eein-eeerin, eein+eeerin, col=70
; ---
    file='outputs/dx11_pre_Imap_143_undusted_insideM_ns2048_uK_hrhs_xfcl-fls_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eefls=extract_xfaster_newdata_output(file, lcen=leefls, cler=eeerfls, res=eeresfls, ncl='EE') 
    oplot, leefls,eefls, psym=4
    errplot, leefls,eefls-eeerfls, eefls+eeerfls
    xyouts, 50, 65, 'preDX11 143GHz undusted-IN', chars=1.25, col=70
    xyouts, 50, 61, 'preDX11 143GHz undusted-IN BB-debiased', chars=1.25
    oplot, l[il], tcl[il,1]*ll[il], col=245
 ; ---
    plot, l[il], tcl[il,1]*ll[il], chars=2.5, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,250], yr=[0,3],tit=title
    oplot, leein,eein, psym=6, col=70
    errplot, leein,eein-eeerin, eein+eeerin, col=70
; ---
    oplot, l[il], tcl[il,1]*ll[il], col=245
    oplot, leefls,eefls, psym=4
    errplot, leefls,eefls-eeerfls, eefls+eeerfls
; ---
    plot, l[il], tcl[il,1]*ll[il], chars=2.5, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,50], yr=[0,1];,tit=title
    oplot, leein,eein, psym=6, col=70
    errplot, leein,eein-eeerin, eein+eeerin, col=70
; ---
    oplot, l[il], tcl[il,1]*ll[il], col=245
    oplot, leefls,eefls, psym=4
    errplot, leefls,eefls-eeerfls, eefls+eeerfls
; ---

;   3) full mission vs yearly maps  uncleaned and cleaned for inter ell
    title='preDX11 143 GHz: full mission VS yearly maps'
    !p.multi=[0,1,2]
    window, win+4, xsize=720, ysize=450*1.8,tit=title
    plot, l[il], tcl[il,1]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,lmax], yr=[0,70], ys=1,tit=title
    oplot, l[il], tcl[il,1]*ll[il], col=245
; ---
    file='outputs/dx11_pre_Imap_143_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eeraw=extract_xfaster_newdata_output(file, lcen=lee, cler=eeer, res=eeres, ncl='EE') 
    oplot, lee,eeraw, psym=4
    errplot, lee,eeraw-eeer, eeraw+eeer
; ---
    file='outputs/dx11_pre_Imap_143_undusted_insideM_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eein=extract_xfaster_newdata_output(file, lcen=leein, cler=eeerin, res=eeresin, ncl='EE') 
    oplot, leein,eein, psym=6, col=70
    errplot, leein,eein-eeerin, eein+eeerin, col=70
; ---
    file='outputs/dx11_pre_Imap_143_year_1-2_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eey=extract_xfaster_newdata_output(file, lcen=leey, cler=eeery, res=eeresy, ncl='EE') 
    oplot, leey,eey, psym=6, col=150
    errplot, leey,eey-eeery, eey+eeery, col=150
; ---
    file='outputs/dx11_pre_Imap_143_year_1-2_undusted_insideM_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eeyin=extract_xfaster_newdata_output(file, lcen=leeyin, cler=eeeryin, res=eeresyin, ncl='EE') 
    oplot, leeyin,eeyin, psym=6, col=100
    errplot, leeyin,eeyin-eeeryin, eeyin+eeeryin, col=100

    xyouts, 100, 65, '143GHz raw', chars=1.5
    xyouts, 100, 61, '143GHz undusted-IN', chars=1.5, col=70
    xyouts, 100, 57, '143GHz yearly raw', chars=1.5, col=150
    xyouts, 100, 53, '143GHz yearly undusted-IN', chars=1.5, col=100

    plot, l[il], tcl[il,1]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,250], yr=[0,3], ys=1,tit=title
    oplot, l[il], tcl[il,1]*ll[il], col=245
; ---
    oplot, lee,eeraw, psym=4
    errplot, lee,eeraw-eeer, eeraw+eeer
; ---
    oplot, leein,eein, psym=6, col=70
    errplot, leein,eein-eeerin, eein+eeerin, col=70
; ---
    oplot, leey,eey, psym=6, col=150
    errplot, leey,eey-eeery, eey+eeery, col=150
; ---
    oplot, leeyin,eeyin, psym=6, col=100
    errplot, leeyin,eeyin-eeeryin, eeyin+eeeryin, col=100
stop
;   4) and with hrhd and bb for all ell
    title='preDX11 143 GHz: full missio VS yearly maps -- BB diagnosis'
    !p.multi=[0,3,1]
    window, win+4, xsize=720*1.9, ysize=450,tit=title
    plot, l[il], tcl[il,1]*ll[il], chars=2, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,lmax], yr=[0,70], ys=1
    oplot, l[il], tcl[il,1]*ll[il], col=245
; ---
    file='outputs/dx11_pre_Imap_143_undusted_insideM_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eein=extract_xfaster_newdata_output(file, lcen=leein, cler=eeerin, res=eeresin, ncl='EE') 
    oplot, leein,eein, psym=6
    errplot, leein,eein-eeerin, eein+eeerin
; ---
    file='outputs/dx11_pre_Imap_143_undusted_insideM_ns2048_uK_hrhs_xfcl-fls_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eeinbb=extract_xfaster_newdata_output(file, lcen=leeinbb, cler=eeerinbb, res=eeresinbb, ncl='EE') 
    oplot, leeinbb, eeinbb, psym=6, col=70
    errplot, leeinbb,eeinbb-eeerinbb, eeinbb+eeerinbb, col=70
; ---
    file='outputs/dx11_pre_Imap_143_year_1-2_undusted_insideM_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eeyin=extract_xfaster_newdata_output(file, lcen=leeyin, cler=eeeryin, res=eeresyin, ncl='EE') 
    oplot, leeyin,eeyin, psym=6, col=150
    errplot, leeyin,eeyin-eeeryin, eeyin+eeeryin, col=150
; ---
    file='outputs/dx11_pre_Imap_143_year_1-2_undusted_insideM_ns2048_uK_hrhs_xfcl-fls_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eeyinbb=extract_xfaster_newdata_output(file, lcen=leeyinbb, cler=eeeryinbb, res=eeresyinbb, ncl='EE') 
    oplot, leeyinbb, eeyinbb, psym=6, col=100
    errplot, leeyinbb, eeyinbb-eeeryinbb, eeyinbb+eeeryinbb, col=100

;    xyouts, 100, 65, '143GHz undusted-IN', chars=1.5
;    xyouts, 100, 61, '143GHz undusted-IN BB-debiased', chars=1.5, col=70
;    xyouts, 100, 57, '143GHz yearly undusted-IN', chars=1.5, col=150
;    xyouts, 100, 53, '143GHz yearly undusted-IN BB debiased', chars=1.5, col=100
; ------
    plot, l[il], tcl[il,1]*ll[il], chars=2, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,250], yr=[0,3], ys=1,tit=title
    oplot, l[il], tcl[il,1]*ll[il], col=245
; ---
    oplot, leein,eein, psym=6
    errplot, leein,eein-eeerin, eein+eeerin
; ---
    oplot, leeinbb, eeinbb, psym=6, col=70
    errplot, leeinbb,eeinbb-eeerinbb, eeinbb+eeerinbb, col=70
; ---
    oplot, leeyin,eeyin, psym=6, col=150
    errplot, leeyin,eeyin-eeeryin, eeyin+eeeryin, col=150
; ---
    oplot, leeyinbb, eeyinbb, psym=6, col=100
    errplot, leeyinbb, eeyinbb-eeeryinbb, eeyinbb+eeeryinbb, col=100
    xyouts, 25, 2.8, '143GHz undusted-IN', chars=1.5
    xyouts, 25, 2.6, '143GHz undusted-IN BB-debiased', chars=1.5, col=70
    xyouts, 25, 2.4, '143GHz yearly undusted-IN', chars=1.5, col=150
    xyouts, 25, 2.2, '143GHz yearly undusted-IN BB debiased', chars=1.5, col=100
; ------
    plot, l[il], tcl[il,1]*ll[il], chars=2, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,50], yr=[0,1], ys=1
    oplot, l[il], tcl[il,1]*ll[il], col=245
; ---
    oplot, leein,eein, psym=6
    errplot, leein,eein-eeerin, eein+eeerin
; ---
    oplot, leeinbb, eeinbb, psym=6, col=70
    errplot, leeinbb,eeinbb-eeerinbb, eeinbb+eeerinbb, col=70
; ---
    oplot, leeyin,eeyin, psym=6, col=150
    errplot, leeyin,eeyin-eeeryin, eeyin+eeeryin, col=150
; ---
    oplot, leeyinbb, eeyinbb, psym=6, col=100
    errplot, leeyinbb, eeyinbb-eeeryinbb, eeyinbb+eeeryinbb, col=100

;   5) yearly vs commander cmb map low ell
    !p.multi=0
    window, win+4, xsize=720, ysize=450,tit=title
    plot, l[il], tcl[il,1]*ll[il], chars=1.5, xtit='!8l', ytit='!8D!dl!uEE!n [!7l!8K!u2!n]', xr=[1,250], yr=[0,2], ys=1
    oplot, l[il], tcl[il,1]*ll[il], col=245
; ---
    file='outputs/dx11_pre_Imap_143_year_1-2_undusted_insideM_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eeyin=extract_xfaster_newdata_output(file, lcen=leeyin, cler=eeeryin, res=eeresyin, ncl='EE') 
    oplot, leeyin,eeyin, psym=6
    errplot, leeyin,eeyin-eeeryin, eeyin+eeeryin
; ---
    file='outputs/dx11_pre_Imap_143_year_1-2_undusted_insideM_ns2048_uK_hrhs_xfcl-fls_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eeyinbb=extract_xfaster_newdata_output(file, lcen=leeyinbb, cler=eeeryinbb, res=eeresyinbb, ncl='EE') 
    oplot, leeyinbb, eeyinbb, psym=6, col=70
    errplot, leeyinbb, eeyinbb-eeeryinbb, eeyinbb+eeeryinbb, col=70
; ---
    file='outputs/dx11_pre_commander_year_map_flat_ns256_uK_hrhs_xfcl_l500_dx9_common_xGal06_mask_ns256_X_QUmask_DX10_256_T5.0_60_X_ps100-353_l500.newdat'
    eecmd=extract_xfaster_newdata_output(file, lcen=leecmd, cler=eeercmd, res=eerescmd, ncl='EE') 
    oplot, leecmd, eecmd, psym=6, col=100
    errplot, leecmd, eecmd-eeercmd, eecmd+eeercmd, col=100
    xyouts, 25, 1.8, '!6143GHz yearly undusted-IN', chars=1.5
    xyouts, 25, 1.65, '!6143GHz yearly undusted-IN BB-debiased', chars=1.5, col=70
    xyouts, 25, 1.5, '!6Commander yearly', chars=1.5, col=100
; ---
    file='outputs/dx11_pre_year_1-2_Imap_143_full_extMask_undusted_insideM_split_flat_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eeyinbb=extract_xfaster_newdata_output(file, lcen=leeyinbb, cler=eeeryinbb, res=eeresyinbb, ncl='EE') 
    oplot, leeyinbb, eeyinbb, psym=6, col=70
    errplot, leeyinbb, eeyinbb-eeeryinbb, eeyinbb+eeeryinbb, col=240

stop
endif
;   6)
    !p.multi=0
    title = '100GHz TE anomaly'
    set_plot, 'ps'
    device, file='preDX11_100GHz_TE_anomaly.eps', /col, bits=8
;    window, win+4, xsize=720, ysize=450,tit=title
    plot, l[il], tcl[il,3]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,lmax], yr=[-150,150], ys=1, tit=title
    oplot, l[il], tcl[il,3]*ll[il], col=245
; ---
    file='outputs/dx11_pre_Imap_100_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    te=extract_xfaster_newdata_output(file, lcen=lte, cler=teer, res=teres, ncl='TE') 
    oplot, lte, te, psym=6
    errplot, lte, te-teer, te+teer
; ---
    file='outputs/dx11_pre_Imap_100_undusted_insideM_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    tein=extract_xfaster_newdata_output(file, lcen=ltein, cler=teerin, res=teresin, ncl='TE') 
    oplot, ltein, tein, psym=6, col=70
    errplot, ltein, tein-teerin, te+teerin, col=70
; ---
    file='outputs/dx11_pre_Imap_100_year1-2_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    tey=extract_xfaster_newdata_output(file, lcen=ltey, cler=teery, res=teresy, ncl='TE') 
    oplot, ltey, tey, psym=6, col=100
    errplot, ltey, tey-teery, tey+teery, col=100
; ---
    file='outputs/dx11_pre_Imap_100_year_1-2_undusted_insideM_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    teyin=extract_xfaster_newdata_output(file, lcen=lteyin, cler=teeryin, res=teresyin, ncl='TE') 
    oplot, lteyin, teyin, psym=6, col=210
    errplot, lteyin, teyin-teeryin, teyin+teeryin, col=210

    xyouts, 1050, 130, 'preDX11 100GHz'
    xyouts, 1050, 115, 'preDX11 100GHz undusted', col=70
    xyouts, 1050, 100, 'preDX11 100GHz yearly', col=100
    xyouts, 1050, 85, 'preDX11 100GHz yearly undusted', col=210
;    endif
    device, /close
    set_plot, 'x'
stop
;   7) All results
    !p.multi=0
    window, win+4, xsize=720, ysize=450,tit=title
    plot, l[il], tcl[il,1]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,250], yr=[0,2], ys=1
    oplot, l[il], tcl[il,1]*ll[il], col=245
; ---
    file='outputs/dx11_pre_Imap_143_undusted_insideM_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eein=extract_xfaster_newdata_output(file, lcen=leein, cler=eeerin, res=eeresin, ncl='EE') 
    oplot, leein,eein, psym=6
    errplot, leein,eein-eeerin, eein+eeerin
; ---
    file='outputs/dx11_pre_Imap_143_undusted_insideM_ns2048_uK_hrhs_xfcl-fls_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eeinbb=extract_xfaster_newdata_output(file, lcen=leeinbb, cler=eeerinbb, res=eeresinbb, ncl='EE') 
    oplot, leeinbb,eeinbb, psym=4, col=210
    errplot, leeinbb,eeinbb-eeerinbb, eeinbb+eeerinbb, col=210
; ---
    file='outputs/dx11_pre_Imap_143_year_1-2_undusted_insideM_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eeyin=extract_xfaster_newdata_output(file, lcen=leeyin, cler=eeeryin, res=eeresyin, ncl='EE') 
    oplot, leeyin,eeyin, psym=6, col=70
    errplot, leeyin,eeyin-eeeryin, eeyin+eeeryin, col=70
; ---
    file='outputs/dx11_pre_Imap_143_year_1-2_undusted_insideM_ns2048_uK_hrhs_xfcl-fls_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    eeinybb=extract_xfaster_newdata_output(file, lcen=leeinybb, cler=eeerinybb, res=eeresinybb, ncl='EE') 
    oplot, leeinybb,eeinybb, psym=4, col=150
    errplot, leeinybb,eeinybb-eeerinybb, eeinybb+eeerinybb, col=150
; ---
    file='outputs/dx11_pre_commander_year_map_flat_ns256_uK_hrhs_xfcl_l500_dx9_common_xGal06_mask_ns256_X_QUmask_DX10_256_T5.0_60_X_ps100-353_l500.newdat'
    eecmd=extract_xfaster_newdata_output(file, lcen=leecmd, cler=eeercmd, res=eerescmd, ncl='EE') 
    oplot, leecmd, eecmd, psym=6, col=100
    errplot, leecmd, eecmd-eeercmd, eecmd+eeercmd, col=100

    xyouts, 25, 1.8, '143GHz undusted-IN', chars=1.5
    xyouts, 25, 1.65, '143GHz undusted-IN BB-debiased', chars=1.5, col=210
    xyouts, 25, 1.5, '143GHz yearly undusted-IN', chars=1.5, col=70
    xyouts, 25, 1.35, '143GHz yearly undusted-IN BB-debiased', chars=1.5, col=150
    xyouts, 25, 1.2, 'Commander yearly', chars=1.5, col=100

stop
end
