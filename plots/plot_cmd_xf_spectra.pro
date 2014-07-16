pro plot_cmd_xf_spectra, init=init, lmax=lmax, win=win

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

;   5) yearly vs commander cmb map low ell
    !p.multi=[0,3,2]
;    window, win+4, xsize=1300, ysize=450*1300/720,tit=title
    set_plot, 'ps'
    device, file='preDX11_commander.eps', /col, bits=8, /landscape
    plot, l[il], tcl[il,0]*ll[il], chars=1.5, ytit='!8D!dl!u!6TT!n [!7l!8K!u2!n]', xr=[1,250], yr=[0,6500], ys=1, position=[0.065,0.45,0.33,0.975], xtickname=strarr(6)+' '
    oplot, l[il], tcl[il,0]*ll[il], col=245
; ---
    file='outputs/dx11_pre_commander_year_map_flat_ns256_uK_hrhs_xfcl_l500_dx9_common_xGal06_mask_ns256_X_QUmask_DX10_256_T5.0_60_X_ps100-353_l500.newdat'
    ttcmd=extract_xfaster_newdata_output(file, lcen=lttcmd, cler=ttercmd, res=ttrescmd, ncl='TT') 
    eecmd=extract_xfaster_newdata_output(file, lcen=leecmd, cler=eeercmd, res=eerescmd, ncl='EE', btcl=beecmd) 
    tecmd=extract_xfaster_newdata_output(file, lcen=ltecmd, cler=teercmd, res=terescmd, ncl='TE') 
    oplot, lttcmd, ttcmd, psym=6
    errplot, lttcmd, ttcmd-ttercmd, ttcmd+ttercmd
; ---
    file='outputs/dx11_pre_year_1-2_Imap_143_full_extMask_undusted_insideM_split_flat_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
    ttyin = extract_xfaster_newdata_output(file, lcen=lttyin, cler=tteryin, res=ttresyin, ncl='TT') 
    eeyin = extract_xfaster_newdata_output(file, lcen=leeyin, cler=eeeryin, res=eeresyin, ncl='EE') 
    teyin = extract_xfaster_newdata_output(file, lcen=lteyin, cler=teeryin, res=teresyin, ncl='TE') 
    oplot, lttyin, ttyin, psym=6, col=70
    errplot, lttyin, ttyin-tteryin, ttyin+tteryin, col=70
; ---
; Offset removed from plot
;##    file='/global/scratch2/sd/loris/OD/DX/p11/Commander/spice_cls_dx10pre11_v1_y1_x_y2.fits'
;##    fits2cl, xcls, file
;##    bxtt = bp_binning(xcls[*,0]*ll, xfdir+'data/bins/ctp/CTP_bin_TT')
;##    bxee = bp_binning(xcls[*,1]*ll, xfdir+'data/bins/ctp/CTP_bin_EE')
;##    bxte = bp_binning(xcls[*,3]*ll, xfdir+'data/bins/ctp/CTP_bin_TE')
;##    oplot, lttyin, bxtt, psym=6, col=40
;stop
; ---
    plot, l[il], tcl[il,1]*ll[il], chars=1.5, ytit='!8D!dl!u!6EE!n [!7l!8K!u2!n]', xr=[1,250], yr=[0,2.5], ys=1, position=[0.3925,0.45,0.6575,0.975], xtickname=strarr(6)+' '
    oplot, l[il], tcl[il,1]*ll[il], col=245
; ---
    oplot, leecmd, eecmd, psym=6
    errplot, leecmd, eecmd-eeercmd, eecmd+eeercmd
    oplot, leeyin, eeyin, psym=6, col=70
    errplot, leeyin, eeyin-eeeryin, eeyin+eeeryin, col=70
;##    oplot, leeyin, bxee, psym=6, col=210
    xyouts, 25,2.3, '!6Commander 7-band QU';, chars=1.5
;##    xyouts, 25,2.1, '!6Commander PolSpice X-spectrum', chars=1.5, col=210
    xyouts, 25,2.1, '!6143 yr1-2 IN', col=70
; ------
    plot, l[il], tcl[il,3]*ll[il], chars=1.5, ytit='!8D!dl!u!6TE!n [!7l!8K!u2!n]', xr=[1,250], yr=[-75,75], ys=1, position=[0.72,0.45,0.98,0.975], xtickname=strarr(6)+' '
    oplot, l[il], tcl[il,3]*ll[il], col=245
    oplot, ltecmd, tecmd, psym=6
    errplot, ltecmd, tecmd-teercmd, tecmd+teercmd
    oplot, lteyin, teyin, psym=6, col=70
    errplot, lteyin, teyin-teeryin, teyin+teeryin, col=70
;    oplot, lteyin, bxte, psym=6, col=210
;stop
;------ Residuals ------
    plot, l[il], tcl[il,0]*0, chars=1.5, ytit='!6Residuals !8D!dl!u!6TT!n [!7l!8K!u2!n]', xr=[0,250], yr=[-1000,1000], ys=1, xtit='!8l', position=[0.065,0.075,0.33,0.44]
    oplot, l[il], tcl[il,0]*0, col=245
    oplot, lttyin, ttresyin, psym=6, col=70
    errplot, lttyin, ttresyin-tteryin, ttresyin+tteryin, col=70
    oplot, lttcmd, ttrescmd, psym=6, line=2
;##    errplot, lttcmd, ttrescmd-ttercmd, ttrescmd+ttercmd
; ------ EE
    plot, l[il], tcl[il,1]*0, chars=1.5, ytit='!6Residuals !8D!dl!u!6EE!n [!7l!8K!u2!n]', xr=[0,250], yr=[-1,1], ys=1, xtit='!8l', position=[0.3925,0.075,0.6575,0.44]
    oplot, l[il], tcl[il,1]*0, col=245
    oplot, leeyin, eeresyin, psym=6, col=70
    errplot, leeyin, eeresyin-eeeryin, eeresyin+eeeryin, col=70
    oplot, leecmd, eerescmd, psym=6
    errplot, leecmd, eerescmd-eeercmd, eerescmd+eeercmd
;##    oplot, leeyin, bxee-beecmd, psym=6, col=210

; ------ EE
    plot, l[il], tcl[il,1]*0, chars=1.5, ytit='!6Residuals !8D!dl!u!6TE!n [!7l!8K!u2!n]', xr=[0,250], yr=[-15,15], ys=1, xtit='!8l', position=[0.72,0.075,0.98,0.44]
    oplot, l[il], tcl[il,1]*0, col=245
    oplot, lteyin, teresyin, psym=6, col=70
    errplot, lteyin, teresyin-teeryin, teresyin+teeryin, col=70
    oplot, ltecmd, terescmd, psym=6
    errplot, ltecmd, terescmd-teercmd, terescmd+teercmd

    device, /close
    set_plot, 'x'

    stop, ' End End of Program ---'

end
