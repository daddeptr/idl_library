pro compare_xfaster_spectra, file1, file2, init=init, ncl=ncl, win=win, otit=otit, res=res, lmax=lmax, lmin=lmin, leg_tags=leg_tags, uleg_pos=uleg_pos, bleg_pos=bleg_pos, reso=reso, xlog=xlog, ylog=ylog, tclfile=tclfile, psout=psout, dops=dops, percent=percent, binning=binning, err_bars=err_bars

    print, ' 1- ', file1
    print, ' 2- ', file2

    if keyword_set(xlog) then xlog = 1 else xlog=0
    if keyword_set(ylog) then ylog = 1 else ylog=0

    if not keyword_set(leg_tags) then leg_tags = ['Set-1','Set-2']
    if not keyword_set(otit) then otit='XF spectra comparison: Set-1 - Set-2'
    if not keyword_set(win) then win = 2
    if not keyword_set(reso) then reso = 155
    if keyword_set(init) then begin
        mollview, findgen(12l*4l^2), px=500, win=win
        loadct, 39
        !p.color=0
        !p.background=255
    endif
    if not keyword_set(lmax) then lmax=2000l
    if not keyword_set(lmin) then lmin=1
    if not keyword_set(psout) then psout='xf_spec_dfference.eps'

    xfdir = '/global/scratch2/sd/dpietrob/Software/XFaster/'
;##    fits2cl, tcl, xfdir+'data/lcdm_cl_uK.fits'
    if not keyword_set(tclfile) then tclfile = xfdir + 'data/planck_lcdm_cl_uK_xf1.e-3.fits'
    fits2cl, tcl, tclfile
    tcl[*,4] = 0.
    tcl[*,5] = 0.
    l=findgen(n_elements(tcl[*,0]))
    ll=l*(l+1)/2./!pi

;## NB consistent with my definition of .newdat
;##    cl1 = extract_xfaster_newdata_output(file1,ncl=ncl, lcen=lcen1, cler=cler1)
;##    cl2 = extract_xfaster_newdata_output(file2,ncl=ncl, lcen=lcen2, cler=cler2)
    cl1 = extract_xfaster_newdat(file1,ncl=ncl, lcen=lcen1, cler=cler1)
    cl2 = extract_xfaster_newdat(file2,ncl=ncl, lcen=lcen2, cler=cler2)

    if not keyword_set(binning) then binning = xfdir+'data/bins/ctp/CTP_bin'

    if (ncl eq 1) or (strupcase(ncl) eq 'TT')  then btcl = bp_binning(tcl[*,0]*ll, binning+'_TT')
    if (ncl eq 2) or (strupcase(ncl) eq 'EE')  then btcl = bp_binning(tcl[*,1]*ll, binning+'_EE')
    if (ncl eq 3) or (strupcase(ncl) eq 'BB')  then btcl = bp_binning(tcl[*,2]*ll, binning+'_EE')
    if (ncl eq 4) or (strupcase(ncl) eq 'TE')  then btcl = bp_binning(tcl[*,3]*ll, binning+'_TE')
    if (ncl eq 5) or (strupcase(ncl) eq 'TB')  then btcl = bp_binning(tcl[*,4]*ll, binning+'_TE')
    if (ncl eq 6) or (strupcase(ncl) eq 'EB')  then btcl = bp_binning(tcl[*,5]*ll, binning+'_EE')

    il = lindgen(lmax)+2
    !p.multi = [0,1,2]
    if not keyword_set(dops) then begin
       window, win, xsize=720, ysize=450*1.8
       plot, l[il], tcl[il,ncl-1]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!n [!7l!6K!u2!n]', xr=[lmin,lmax], tit=otit, xs=1, xlog=xlog, ylog=ylog, ys=2;, yr=[-500,6500], ys=1
    endif else begin
       set_plot,'ps'
       device, file=psout+'_ncl'+strtrim(string(ncl),2)+'.eps', /col, bits=8, /landscape
       plot, l[il], tcl[il,ncl-1]*ll[il], xtit='!6l', ytit='!6D!dl!n [!7l!6K!u2!n]', xr=[lmin,lmax], tit=otit, xs=1, xlog=xlog, ylog=ylog;, yr=[-500,6500], ys=1
    endelse
    oplot, l[il], tcl[il,ncl-1]*ll[il], col=245
    oplot, lcen1, btcl, col=245, psym=5
    oplot, lcen1, cl1, psym=4, thick=2
    errplot, lcen1, cl1-cler1, cl1+cler1
    oplot, lcen2, cl2, psym=6, col=70, thick=2
    errplot, lcen2, cl2-cler2, cl2+cler2, col=70, line=2
    if not keyword_set(uleg_pos) then uleg_pos = [1250,5500]
    legend,[leg_tags,'Best fit'], psym=[4,6,5], col=[0,70,245], chars=1.2, pos=uleg_pos

    if not keyword_set(percent) then begin
       if not keyword_set(dops) then $
         plot, l[il], tcl[il,ncl-1]*0., chars=1.5, xtit='!6l', ytit='!7D!6D!dl!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[-reso, reso], ys=1, xs=1, xlog=xlog else $
         plot, l[il], tcl[il,ncl-1]*0., xtit='!6l', ytit='!7D!6D!dl!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[-reso, reso], ys=1, xs=1, xlog=xlog
       oplot, lcen1, cl1-cl2, psym=-4, thick=2
       errplot, lcen1, cl1-cl2-sqrt(cler1^2+cler2^2), cl1-cl2+sqrt(cler1^2+cler2^2), line=1
       if keyword_set(err_bars) then oplot, lcen1, (cler1-cler2), psym=6, col=70
       oplot, l[il], tcl[il,ncl-1]*0., col=245
;       if not keyword_set(bleg_pos) then bleg_pos = [1250,0.75*reso]
;       legend,['Value difference','Error bar difference'], psym=[4,6], col=[0,70], chars=1.2, pos=bleg_pos
    endif else begin
;       if not keyword_set(dops) then $
;         plot, l[il], tcl[il,ncl-1]*0., chars=1.5, xtit='!6l', ytit='!7D!6D!dl!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[-reso, reso], ys=1, xs=1, xlog=xlog, tit='Black: data points difference; blue: data error difference' else $
;         plot, l[il], tcl[il,ncl-1]*0., xtit='!6l', ytit='!7D!6D!dl!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[-reso, reso], ys=1, xs=1, xlog=xlog, tit='Black: data points difference; blue: data error difference'
;       oplot, lcen1, (cl1-cl2)/((cl1+cl2)/2), psym=4
;       errplot, lcen1, (cl1-cl2-sqrt(cler1^2+cler2^2))/((cl1+cl2)/2), (cl1-cl2+sqrt(cler1^2+cler2^2))/((cl1+cl2)/2), line=1
;       oplot, lcen1, (cler1-cler2)/((cl1+cl2)/2), psym=6, col=70
;       oplot, l[il], tcl[il,ncl-1]*0., col=245
       if not keyword_set(dops) then $
         plot, l[il], tcl[il,ncl-1]*0.+1, chars=1.5, xtit='!6l', ytit='!7D!6D!dl!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[1-reso, 1+reso], ys=1, xs=1, xlog=xlog, tit='Black: data points difference; blue: data error difference' else $
         plot, l[il], tcl[il,ncl-1]*0.+1, xtit='!6l', ytit='!7D!6D!dl!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[1-reso, 1.+reso], ys=1, xs=1, xlog=xlog, tit='Black: data points difference; blue: data error difference'
       oplot, lcen1, (cl1/cl2), psym=4
       errplot, lcen1, (cl1-sqrt(cler1^2+cler2^2))/cl2, (cl1+sqrt(cler1^2+cler2^2))/cl2, line=1
       if keyword_set(err_bars) then oplot, lcen1, (cler1/cler2), psym=6, col=70
       oplot, l[il], tcl[il,ncl-1]*0.+1, col=245
;       if not keyword_set(bleg_pos) then bleg_pos = [lmin*1.1,0.75*reso]
;       legend,['Value difference','Error bar difference'], psym=[4,6], col=[0,70], chars=1.2, pos=bleg_pos       
    endelse

    if keyword_set(dops) then begin
       device, /close
       set_plot, 'x'
    endif
!p.multi = 0
;stop
end
