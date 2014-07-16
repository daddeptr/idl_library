pro plot_hybrid_cmd_yr, init=init, lmax=lmax, win=win

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

    !p.multi=[0,3,2]
;    window, win+4, xsize=1300, ysize=450*1300/720,tit=title
    set_plot, 'ps'
    device, file='preDX11_hybrid_cmd_xf-yr.eps', /landscape, /col, bit=8
    plot, l[il], tcl[il,0]*ll[il], chars=1.5, ytit='!8D!dl!u!6TT!n [!7l!8K!u2!n]', xr=[1,2000], yr=[0,6500], ys=1, position=[0.06,0.45,0.325,0.975], xtickname=strarr(6)+' '
    oplot, l[il], tcl[il,0]*ll[il], col=245
; ---
    file='outputs/dx11_pre_commander_year_map_flat_ns256_uK_hrhs_xfcl_l500_dx9_common_xGal06_mask_ns256_X_QUmask_DX10_256_T5.0_60_X_ps100-353_l500.newdat'
    ttcmd=extract_xfaster_newdata_output(file, lcen=lttcmd, cler=ttercmd, res=ttrescmd, ncl='TT') 
    eecmd=extract_xfaster_newdata_output(file, lcen=leecmd, cler=eeercmd, res=eerescmd, ncl='EE', btcl=beecmd) 
    tecmd=extract_xfaster_newdata_output(file, lcen=ltecmd, cler=teercmd, res=terescmd, ncl='TE') 

    file='outputs/dx11_pre_year_1-2_Imap_143_full_extMask_undusted_insideM_split_flat_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    ttyin = extract_xfaster_newdata_output(file, lcen=lttyin, cler=tteryin, res=ttresyin, ncl='TT') 
    eeyin = extract_xfaster_newdata_output(file, lcen=leeyin, cler=eeeryin, res=eeresyin, ncl='EE') 
    teyin = extract_xfaster_newdata_output(file, lcen=lteyin, cler=teeryin, res=teresyin, ncl='TE') 

    btt = bp_binning(tcl[*,0]*ll, xfdir+'data/bins/ctp/CTP_bin_TT')
    bee = bp_binning(tcl[*,1]*ll, xfdir+'data/bins/ctp/CTP_bin_EE')
    bte = bp_binning(tcl[*,3]*ll, xfdir+'data/bins/ctp/CTP_bin_TE')

    iswitch = where(lttcmd gt 250)
    tt = ttyin
    tt[0:iswitch[0]-1]=ttcmd[0:iswitch[0]-1]
    tter = tteryin
    tter[0:iswitch[0]-1]=ttercmd[0:iswitch[0]-1]

    iswitch = where(leecmd gt 250)
    ee = eeyin
    ee[0:iswitch[0]-1]=eecmd[0:iswitch[0]-1]
    eeer = eeeryin
    eeer[0:iswitch[0]-1]=eeercmd[0:iswitch[0]-1]

    iswitch = where(ltecmd gt 250)
    te = teyin
    te[0:iswitch[0]-1]=tecmd[0:iswitch[0]-1]
    teer = teeryin
    teer[0:iswitch[0]-1]=teercmd[0:iswitch[0]-1]

    oplot, lttyin, tt, psym=6
    errplot, lttyin, tt-tter, tt+tter

; ------
    plot, l[il], tcl[il,1]*ll[il], chars=1.5, ytit='!8D!dl!u!6EE!n [!7l!8K!u2!n]', xr=[1,2000], yr=[0,50], ys=1, position=[0.3875,0.45,0.6525,0.975], xtickname=strarr(6)+' '
    oplot, l[il], tcl[il,1]*ll[il], col=245
    oplot, leeyin, ee, psym=6
    errplot, leeyin, ee-eeer, ee+eeer

; ------
    plot, l[il], tcl[il,3]*ll[il], chars=1.5, ytit='!8D!dl!u!6TE!n [!7l!8K!u2!n]', xr=[1,2000], yr=[-175,175], ys=1, position=[0.715,0.45,0.98,0.975], xtickname=strarr(6)+' '
    oplot, l[il], tcl[il,3]*ll[il], col=245
    oplot, lteyin, te, psym=6
    errplot, lteyin, te-teer, te+teer

;------ Residuals ------
    plot, l[il], tcl[il,0]*0, chars=1.5, ytit='!6Residuals D!dl!u!6TT!n [!7l!8K!u2!n]', xr=[0,2000], yr=[-750,750], xtit='!8l', position=[0.06,0.075,0.325,0.44]
    oplot, l[il], tcl[il,0]*0, col=245
    oplot, lttyin, tt-btt, psym=6
    errplot, lttyin, tt-btt-tter, tt-btt+tter
; ------ EE
    plot, l[il], tcl[il,1]*0, chars=1.5, ytit='!6Residuals D!dl!u!6EE!n [!7l!8K!u2!n]', xr=[0,2000], yr=[-5,5], xtit='!8l', position=[0.3875,0.075,0.6525,0.44]
    oplot, l[il], tcl[il,1]*0, col=245
    oplot, leeyin, ee-bee, psym=6
    errplot, leeyin, ee-bee-eeer, ee-bee+eeer
; ------ TE
    plot, l[il], tcl[il,1]*0, chars=1.5, ytit='!6Residuals D!dl!u!6TE!n [!7l!8K!u2!n]', xr=[0,2000], yr=[-20,20], ys=1, xtit='!8l', position=[0.715,0.075,0.98,0.44]
    oplot, l[il], tcl[il,1]*0, col=245
    oplot, lteyin, te-bte, psym=6
    errplot, lteyin, te-bte-teer, te-bte+teer

    device, /close
    set_plot, 'x'

    stop, ' End End of Program ---'

end
