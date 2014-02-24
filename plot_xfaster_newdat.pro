pro plot_xfaster_newdat, file, init=init, pol=pol, win=win, otit=otit, residuals=residuals, xpol=xpol, lmax=lmax, chkbb=chkbb, xlog=xlog, resolution=resolution, tclfile=tclfile, noisefile=noisefile, beams=beams, dops=dops, psout=psout, lmin=lmin, bestfit=bestfit, binning=binning

;##    print, 'Plotting ',n_elements(file), 'files...'

    if not keyword_set(otit) then otit='XF T/P spectra'
    if not keyword_set(resolution) then resolution = 1.
    if not keyword_set(win) then win = 2
    if keyword_set(init) then begin
        mollview, findgen(12l*4l^2), px=650, win=win
        loadct, 39
        !p.color=0
        !p.background=255
    endif
    if not keyword_set(lmin) then lmin=1
    if not keyword_set(lmax) then lmax=2000
    if not keyword_set(beams) then beams=7.03
    if not keyword_set(psout) then psout='xf_spectra.eps'

    xfdir = '/global/scratch2/sd/dpietrob/Software/XFaster/'
    if not keyword_set(binning) then binning=xfdir+'data/bins/ctp/CTP_bin'
;##    fits2cl, tcl, xfdir+'data/lcdm_cl_uK.fits'
    if not keyword_set(tclfile) then begin
;##        fits2cl, tcl, xfdir+'data/camb_97283428_lensedcls_uK.fits'
        fits2cl, tcl, xfdir+'data/planck_lcdm_cl_uK_xf1.e-3.fits'
        if n_elements(tcl[0,*]) lt 6 then begin
            cl=tcl
            tcl = fltarr(n_elements(cl[*,0]),6)
            tcl[*,0:3] = cl
            tcl[*,4] = 0.
            tcl[*,5] = 0.
        endif
;##        tcl[*,4] = 0.
;##        tcl[*,5] = 0.
    endif else begin
        fits2cl, tcl, tclfile
        if n_elements(tcl[0,*]) lt 6 then begin
            cl=tcl
            tcl = fltarr(n_elements(cl[*,0]),6)
            tcl[*,0:3] = cl
            tcl[*,4] = 0.
            tcl[*,5] = 0.
        endif
    endelse
    l=findgen(n_elements(tcl[*,0]))
    ll=l*(l+1)/2./!pi

    if keyword_set(bestfit) then begin
        readcol, bestfit, l, tt, ee, te
    endif
;## NB consistent with my definition of .newdat
;##    if keyword_set(old) then readcol,file,nbtt,nbee,nbbb,nbtb,nbte,nbeb, format='f,f,f,f,f,f',skipline=2,numline=1 else $
;##      readcol,file,nbtt,nbee,nbbb,nbte,nbtb,nbeb, format='f,f,f,f,f,f',skipline=2,numline=1
    head_nlines = 13
    readcol,file,nbtt,nbee,nbbb,nbtb,nbte,nbeb, format='f,f,f,f,f,f',skipline=1,numline=1

    print, ' TT # bin = ', nbtt
    print, ' EE # bin = ', nbee
    print, ' TE # bin = ', nbte
    print, ' BB # bin = ', nbbb
    print, ' TB # bin = ', nbtb
    print, ' EB # bin = ', nbeb

;    print, 6+2*nbtt+2*nbee+2*nbbb

; ------- Binning
;    if not keyword_set(bincl) then begin

    btcltt = bp_binning(tcl[*,0]*ll, binning+'_TT')
    btclee = bp_binning(tcl[*,1]*ll, binning+'_EE')
    btclte = bp_binning(tcl[*,3]*ll, binning+'_TE')
    btclbb = bp_binning(tcl[*,2]*ll, binning+'_EE')
    if keyword_set(xpol) then begin
        btcltb = bp_binning(tcl[*,4]*ll, binning+'_TE')
        btcleb = bp_binning(tcl[*,5]*ll, binning+'_EE')
    endif

    if keyword_set(bestfit) then begin
        x = fltarr(n_elements(tt)+2)
        x[2:*] = tt
        bftt = bp_binning(x, binning+'_TT')
        x[2:*] = ee
        bfee = bp_binning(x, binning+'_EE')
        x[2:*] = te
        bfte = bp_binning(x, binning+'_TE')
    endif

    if keyword_set(noisefile) then begin
        fits2cl, nls, noisefile
    endif else begin
;##        fits2cl, nls, xfdir+'data/mock_noise_tnls.fits'
        nls = tcl * 0.
    endelse
    bnlstt = bp_binning(nls[*,0]*ll/gaussbeam(beams,4096)^2, binning+'_TT')
    bnlsee = bp_binning(nls[*,1]*ll/gaussbeam(beams,4096)^2, binning+'_EE')
    bnlste = bp_binning(nls[*,3]*ll/gaussbeam(beams,4096)^2, binning+'_TE')
    bnlsbb = bp_binning(nls[*,2]*ll/gaussbeam(beams,4096)^2, binning+'_EE')
    bnlstb = bp_binning(nls[*,4]*ll/gaussbeam(beams,4096)^2, binning+'_TE')
    bnlseb = bp_binning(nls[*,5]*ll/gaussbeam(beams,4096)^2, binning+'_EE')

    if (nbtt[0] gt 0) then readcol,file,cltt,cltter,lmntt,lmxtt,format='x,f,f,x,x,f,f', skipline=head_nlines, numline=nbtt[0]
    if (nbee[0] gt 0) then readcol,file,clee,cleeer,lmnee,lmxee,format='x,f,f,x,x,f,f', skipline=head_nlines+2*nbtt[0]+1, numline=nbee[0]
    if (nbbb[0] gt 0) then readcol,file,clbb,clbber,lmnbb,lmxbb,format='x,f,f,x,x,f,f', skipline=head_nlines+2*nbtt[0]+2*nbee[0]+1+1, numline=nbbb[0]
    if (nbte[0] gt 0) then readcol,file,clte,clteer,lmnte,lmxte,format='x,f,f,x,x,f,f', skipline=head_nlines+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+1+1+1, numline=nbte[0]
    if (nbtb[0] gt 0) then readcol,file,cltb,cltber,lmntb,lmxtb,format='x,f,f,x,x,f,f', skipline=head_nlines+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+2*nbte[0]+1+1+1+1, numline=nbtb[0]
    if (nbeb[0] gt 0) then readcol,file,cleb,cleber,lmneb,lmxeb,format='x,f,f,x,x,f,f', skipline=head_nlines+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+2*nbte[0]+2*nbtb[0]+1+1+1+1+1, numline=nbeb[0]


    if (not keyword_set(pol)) then begin
        il = findgen(max(lmxtt))+2
        if not keyword_set(residuals) then begin
            !p.multi = 0
            if not keyword_set(dops) then begin
;                set_plot, 'x'
                window, win, xsize=720, ysize=450,tit=otit 
            endif else begin
                set_plot, 'ps'
                device, file=psout, /col, bits=8, /landscape
            endelse
        endif else begin
            !p.multi=[0,1,2]
;            window, win, xsize=720, ysize=450*1.5,tit=otit
            if not keyword_set(dops) then begin
;                set_plot, 'x'
                window, win, xsize=720, ysize=450*1.5,tit=otit 
            endif else begin
                set_plot, 'ps'
                device, file=psout, /col, bits=8, /landscape
            endelse
        endelse

        if not keyword_set(xlog) then plot, l[il], tcl[il,0]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uTT!n [!7l!6K!u2!n]', xr=[lmin,lmax]
        if keyword_set(xlog) then plot, l[il], tcl[il,0]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uTT!n [!7l!6K!u2!n]', xr=[lmin,lmax], /xlog, xs=1
        oplot, (lmntt+lmxtt)/2,btcltt, psym=4, col=245
;##        oplot, l[2:*], tt, col=210
;##        oplot, (lmntt+lmxtt)/2,bnlstt/sqrt((lmxtt-lmntt+1)), psym=4, col=215
        oplot, (lmntt+lmxtt)/2,cltt, psym=4
        errplot, (lmntt+lmxtt)/2,cltt-cltter, cltt+cltter
        oplot, l[il], tcl[il,0]*ll[il], col=245
        xyouts, max(il), 5500, 'TT Spectrum', chars=1.5

; ------ Residuals
        if (keyword_set(residuals)) then begin
            if keyword_set(bestfit) then btcltt = bftt
            if not keyword_set(xlog) then plot, (lmntt+lmxtt)/2,cltt-btcltt, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=1.5, yr=[-250,250]*resolution, xr=[lmin,lmax]
            if keyword_set(xlog) then plot, (lmntt+lmxtt)/2,cltt-btcltt, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=1.5, yr=[-250,250]*resolution, xr=[lmin,lmax], xs=1, /xlog
            errplot, (lmntt+lmxtt)/2, cltt-btcltt-cltter, cltt-btcltt+cltter
            oplot, (lmntt+lmxtt)/2,cltt*0., thick=0.5, col=245
        endif

    endif else begin
        il = findgen(max(lmxtt))+2
        !p.multi = [0,3,1]
        if (keyword_set(residuals)) then !p.multi=[0,3,2]
;        window, win, xsize=720*1.8, ysize=450*1.25,tit=otit
        if not keyword_set(dops) then begin
;            set_plot, 'x'
            window, win, xsize=720*1.8, ysize=450*1.25,tit=otit
        endif else begin
            set_plot, 'ps'
            device, file=psout, /col, bits=8, /landscape
        endelse
        print, lmntt[0], lmxtt[0]
        if not keyword_set(xlog) then plot, l[il], tcl[il,0]*ll[il], chars=1.5, ytit='!6D!dl!uTT!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[-500,6500], ys=1, position=[0.06,0.45,0.325,0.975], xtickname=strarr(7)+' '
        if keyword_set(xlog) then plot, l[il], tcl[il,0]*ll[il], chars=1.5, ytit='!6D!dl!uTT!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[-500,6500], /xlog, xs=1, ys=1
        oplot, (lmntt+lmxtt)/2,btcltt, psym=4, col=245
        if keyword_set(bestfit) then oplot, l[2:*], tt, col=210
;##        oplot, (lmntt+lmxtt)/2,bnlstt/sqrt((lmxtt-lmntt+1)), psym=4, col=215
        oplot, (lmntt+lmxtt)/2,cltt, psym=4
        errplot, (lmntt+lmxtt)/2,cltt-cltter, cltt+cltter
        oplot, l[il], tcl[il,0]*ll[il], col=245

        print, lmnee[0], lmxee[0]
        il = findgen(max(lmxee))+2
        ll=l*(l+1)/2./!pi
        if not keyword_set(xlog) then plot, l[il], tcl[il,1]*ll[il], chars=1.5, ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[lmin,lmax], position=[0.3875,0.45,0.6525,0.975], xtickname=strarr(7)+' '
        if keyword_set(xlog) then plot, l[il], tcl[il,1]*ll[il], chars=1.5, ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[lmin,lmax], /xlog, xs=1
        oplot, (lmnee+lmxee)/2,btclee, psym=4, col=245
        if keyword_set(bestfit) then oplot, l[2:*], ee, col=210
;##        oplot, (lmnee+lmxee)/2,bnlsee/sqrt((lmxee-lmnee+1)), psym=4, col=215
        oplot, (lmnee+lmxee)/2,clee, psym=4
        errplot, (lmnee+lmxee)/2,clee-cleeer, clee+cleeer
        oplot, l[il], tcl[il,1]*ll[il], col=245

        print, lmnte[0], lmxte[0]
        il = findgen(max(lmxte))+2
        ll=l*(l+1)/2./!pi
        if not keyword_set(xlog) then plot, l[il], tcl[il,3]*ll[il], chars=1.5, ytit='!6D!dl!uTE!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[-200,200], ys=1, position=[0.715,0.45,0.98,0.975], xtickname=strarr(7)+' '
        if keyword_set(xlog) then plot, l[il], tcl[il,3]*ll[il], chars=1.5, ytit='!6D!dl!uTE!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[-200,200], /xlog, xs=1, ys=1
        oplot, (lmnte+lmxte)/2,btclte, psym=4, col=245
        if keyword_set(bestfit) then oplot, l[2:*], te, col=210
        oplot, (lmnte+lmxte)/2,clte, psym=4
        errplot, (lmnte+lmxte)/2,clte-clteer, clte+clteer
        oplot, l[il], tcl[il,3]*ll[il], col=245

; ------ Residuals
        if (keyword_set(residuals)) then begin
            if keyword_set(bestfit) then begin
                btcltt = bftt
                btclee = bfee
                btclte = bfte
            endif
            if not keyword_set(xlog) then plot, (lmntt+lmxtt)/2,cltt-btcltt, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=1.5, yr=[-75,75]*resolution, xr=[lmin,lmax], position=[0.06,0.075,0.325,0.44]
            if keyword_set(xlog) then plot, (lmntt+lmxtt)/2,cltt-btcltt, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=1.5, yr=[-75,75]*resolution, xr=[lmin,lmax], xs=1, /xlog
;##            errplot, (lmntt+lmxtt)/2, -bnlstt*sqrt(2./(2.*((lmntt+lmxtt)/2)+1)/(lmxtt-lmntt+1)), bnlstt*sqrt(2./(2.*((lmntt+lmxtt)/2)+1)/(lmxtt-lmntt+1)), col=215, line=2
            errplot, (lmntt+lmxtt)/2, cltt-btcltt-cltter, cltt-btcltt+cltter
;##            errplot, (lmntt+lmxtt)/2, cltt-btcltt-bnlstt*sqrt(2./(2.*((lmntt+lmxtt)/2)+1)/(lmxtt-lmntt+1)), cltt-btcltt+bnlstt*sqrt(2./(2.*((lmntt+lmxtt)/2)+1)/(lmxtt-lmntt+1)), col=100, line=4
            oplot, (lmntt+lmxtt)/2,cltt*0., thick=0.5, col=245

; --- Res EE
            if not keyword_set(xlog) then plot, (lmnee+lmxee)/2,clee-btclee, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=1.5, xr=[lmin,lmax], yr=[-4,4]*resolution, ys=1, position=[0.3875,0.075,0.6525,0.44] 
            if keyword_set(xlog) then plot, (lmnee+lmxee)/2,clee-btclee, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=1.5, xr=[lmin,lmax], yr=[-4,4]*resolution, ys=1, xs=1, /xlog
;##            errplot, (lmnee+lmxee)/2, -bnlsee*sqrt(2./(2.*((lmnee+lmxee)/2)+1)/(lmxtt-lmntt+1)), bnlsee*sqrt(2./(2.*((lmnee+lmxee)/2)+1)/(lmxtt-lmntt+1)), col=215, line=2
            errplot, (lmnee+lmxee)/2, clee-btclee-cleeer, clee-btclee+cleeer
            oplot, (lmnee+lmxee)/2,clee*0., thick=0.5, col=245

; --- Res TE
;            if not keyword_set(xlog) then plot, (lmnte+lmxte)/2,clte-btclte, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=1.5, xr=[lmin,lmax], yr=[-15,15]*resolution, position=[0.715,0.075,0.98,0.44] 
;            if keyword_set(xlog) then plot, (lmnte+lmxte)/2,clte-btclte, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=1.5, xr=[lmin,lmax], yr=[-15,15]*resolution, xs=1, /xlog
            plot, (lmnte+lmxte)/2, (clte-btclte)/((lmnte+lmxte)/2), psym=4, xtit='!6l', ytit='!6Residuals (!8D!dl!n-Bf)/l!6[!7l!6K!u2!n]', chars=1.5, xr=[lmin,lmax], yr=[-0.05,0.05]*resolution, xs=1, xlog=xlog, position=[0.715,0.075,0.98,0.44]
            errplot, (lmnte+lmxte)/2, (clte-btclte-clteer)/((lmnte+lmxte)/2), (clte-btclte+clteer)/((lmnte+lmxte)/2)
            oplot, (lmnte+lmxte)/2,clte*0., thick=0.5, col=245
        endif

    endelse

    if (keyword_set(xpol)) then begin
        il = findgen(max(lmxtt))+2
        !p.multi = [0,3,2]
;        window, win, xsize=720*1.8, ysize=450*1.25,tit=otit
        if not keyword_set(dops) then begin
;            set_plot, 'x'
            window, win, xsize=720*1.8, ysize=450*1.25,tit=otit
        endif else begin
            set_plot, 'ps'
            device, file=psout, /col, bits=8, /landscape
        endelse
; --- TT
        print, lmntt[0],lmxtt[0]
        if not keyword_set(xlog) then plot, l[il], tcl[il,0]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uTT!n [!7l!6K!u2!n]', xr=[lmin,lmax]
        if keyword_set(xlog) then plot, l[il], tcl[il,0]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uTT!n [!7l!6K!u2!n]', xr=[lmin,lmax], xs=1, /xlog
        oplot, (lmntt+lmxtt)/2,btcltt, psym=4, col=245
        oplot, (lmntt+lmxtt)/2,cltt, psym=4
        errplot, (lmntt+lmxtt)/2,cltt-cltter, cltt+cltter
        oplot, l[il], tcl[il,0]*ll[il], col=245
; --- EE
        print, lmnee[0],lmxee[0]
        il = findgen(max(lmxee))+2
        ll=l*(l+1)/2./!pi
        if not keyword_set(xlog) then plot, l[il], tcl[il,1]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[lmin,lmax]
        if keyword_set(xlog) then plot, l[il], tcl[il,1]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[lmin,lmax], xs=1, /xlog
        oplot, (lmnee+lmxee)/2,btclee, psym=4, col=245
        oplot, (lmnee+lmxee)/2,clee, psym=4
        errplot, (lmnee+lmxee)/2,clee-cleeer, clee+cleeer
        oplot, l[il], tcl[il,1]*ll[il], col=245

; --- BB
        print, lmnbb[0],lmxbb[0]
        if (keyword_set(chkbb)) then begin
            oplot, (lmnee+lmxee)/2,clee-clbb, col=90, psym=6
            legend, ['C!dl!uEE!n-C!dl!uBB!n'], psym=6, col=90, /top, /left, chars=1.8
        endif
        il = findgen(max(lmxbb))+2
        ll=l*(l+1)/2./!pi
        if not keyword_set(xlog) then plot, l[il], tcl[il,2]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uBB!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[0,max(tcl[where(il lt lmax),1]*ll[where(il lt lmax)])]
        if keyword_set(xlog) then plot, l[il], tcl[il,2]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uBB!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[0,max(tcl[where(il lt lmax),1]*ll[where(il lt lmax)])], xs=1, /xlog
        oplot, (lmnbb+lmxbb)/2,btclbb, psym=4, col=245
        oplot, (lmnbb+lmxbb)/2,clbb, psym=4
        errplot, (lmnbb+lmxbb)/2,clbb-clbber, clbb+clbber
        oplot, l[il], tcl[il,2]*ll[il], col=245
        if (keyword_set(chkbb)) then begin
            oplot, (lmnee+lmxee)/2,clee-btclee, col=90, psym=6
            legend, ['C!dl!uEE!n-C!dl!uEE,Th!n'], psym=6, col=90, /top, /right, chars=1.8
        endif
; --- TE
        print, lmnte[0],lmxte[0]
        il = findgen(max(lmxte))+2
        ll=l*(l+1)/2./!pi
        if not keyword_set(xlog) then plot, l[il], tcl[il,3]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uTE!n [!7l!6K!u2!n]', xr=[lmin,lmax]
        if keyword_set(xlog) then plot, l[il], tcl[il,3]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uTE!n [!7l!6K!u2!n]', xr=[lmin,lmax], xs=1, /xlog
        oplot, (lmnte+lmxte)/2,btclte, psym=4, col=245
        oplot, (lmnte+lmxte)/2,clte, psym=4
        errplot, (lmnte+lmxte)/2,clte-clteer, clte+clteer
        oplot, l[il], tcl[il,3]*ll[il], col=245
; --- TB
        print, lmntb[0],lmxtb[0]
        il = findgen(max(lmxtb))+2
        ll=l*(l+1)/2./!pi
        if not keyword_set(xlog) then plot, l[il], tcl[il,4]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uTB!n [!7l!6K!u2!n]', xr=[lmin,lmax],yr=[-50,50]
        if keyword_set(xlog) then plot, l[il], tcl[il,4]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uTB!n [!7l!6K!u2!n]', xr=[lmin,lmax],yr=[-50,50], xs=1, /xlog
        oplot, (lmntb+lmxtb)/2,btcltb, psym=4, col=245
        oplot, (lmntb+lmxtb)/2,cltb, psym=4
        errplot, (lmntb+lmxtb)/2,cltb-cltber, cltb+cltber
        oplot, l[il], tcl[il,4]*ll[il], col=245
; --- EB
        print, lmneb[0],lmxeb[0]
        il = findgen(max(lmxeb))+2
        ll=l*(l+1)/2./!pi
        if not keyword_set(xlog) then plot, l[il], tcl[il,5]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uEB!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[-1,1]
        if keyword_set(xlog) then plot, l[il], tcl[il,5]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uEB!n [!7l!6K!u2!n]', xr=[lmin,lmax], yr=[-1,1], xs=1, /xlog
        oplot, (lmneb+lmxeb)/2,btcleb, psym=4, col=245
        oplot, (lmneb+lmxeb)/2,cleb, psym=4
        errplot, (lmneb+lmxeb)/2,cleb-cleber, cleb+cleber
        oplot, l[il], tcl[il,5]*ll[il], col=245
    endif
    
    if keyword_set(dops) then begin
        device, /close
        set_plot, 'x'
    endif
!p.multi = 0
;stop
end
