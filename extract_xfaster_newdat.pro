function extract_xfaster_newdat, file, ncl=ncl, lcen=lcen, cler=cler, xfdir=xfdir, res=res, btcl=btcl, tclfile=tclfile, slines=slines, binning=binning, covmat=covmat, lmin=lmin, lmax=lmax

    if (not keyword_set(ncl)) then ncl = 1
    if (not keyword_set(xfdir)) then xfdir = '/global/scratch2/sd/dpietrob/Software/XFaster/'
;##    fits2cl, tcl, xfdir+'data/lcdm_cl_uK.fits'
    if (not keyword_set(tclfile)) then tclfile=xfdir+'data/camb_97283428_lensedcls_uK.fits'
    fits2cl, tcl, tclfile
    if n_elements(tcl[0,*]) le 4 then begin
        cls = tcl
        tcl = fltarr(n_elements(cls[*,0]), 6)
        tcl[*,0:3] = cls
        tcl[*,4] = 0.
        tcl[*,5] = 0.
        cls = 0
    endif

    if not keyword_set(slines) then slines = 13
;## NB consistent with my definition of .newdat
;##    if (ncl eq 1) or (strupcase(ncl) eq 'TT') then begin
;##        readcol,file,nbtt, format='f',skipline=1,numline=1
;##        print, ' TT # bin = ', nbtt
;##    endif else begin
;## ;##        if keyword_set(old) then readcol,file,nbtt,nbee,nbbb,nbtb,nbte,nbeb, format='f,f,f,f,f,f',skipline=2,numline=1 else $
;## ;##          readcol,file,nbtt,nbee,nbbb,nbte,nbtb,nbeb, format='f,f,f,f,f,f',skipline=2,numline=1
        readcol,file,nbtt,nbee,nbbb,nbtb,nbte,nbeb, format='f,f,f,f,f,f',skipline=1,numline=1
        
        print, ' TT # bin = ', nbtt
        print, ' EE # bin = ', nbee
        print, ' TE # bin = ', nbte
        print, ' BB # bin = ', nbbb
        print, ' TB # bin = ', nbtb
        print, ' EB # bin = ', nbeb

        nbins = nbtt+nbee+nbbb+nbtb+nbte+nbeb
        print, ' Total # bin = ', nbins

;##    endelse

; ------- Binning
    if not keyword_set(binning) then binning = xfdir+'data/bins/ctp/CTP_bin'
    if (ncl eq 1) or (strupcase(ncl) eq 'TT') then begin
        readcol,file,cl,cler,lmn,lmx,format='x,f,f,x,x,f,f', skipline=slines, numline=nbtt[0]
        btcl = xf_binning(tcl[*,0], binning+'_TT', lcen=lbin)
    endif
    if (ncl eq 2) or (strupcase(ncl) eq 'EE') then begin
        readcol,file,cl,cler,lmn,lmx,format='x,f,f,x,x,f,f', skipline=2*nbtt[0]+slines+1, numline=nbee[0]
        btcl = xf_binning(tcl[*,1], binning+'_EE', lcen=lbin)
    endif
    if (ncl eq 3) or (strupcase(ncl) eq 'BB') then begin
        readcol,file,cl,cler,lmn,lmx,format='x,f,f,x,x,f,f', skipline=2*nbtt[0]+2*nbee[0]+slines+1+1, numline=nbbb[0]
        btcl = xf_binning(tcl[*,2], binning+'_EE', lcen=lbin)
    endif
    if (ncl eq 4) or (strupcase(ncl) eq 'TE') then begin
        readcol,file,cl,cler,lmn,lmx,format='x,f,f,x,x,f,f', skipline=slines+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+1+1+1, numline=nbte[0]
        btcl = xf_binning(tcl[*,3], binning+'_TE', lcen=lbin)
    endif
    if (ncl eq 5) or (strupcase(ncl) eq 'TB') then begin
        readcol,file,cl,cler,lmn,lmx,format='x,f,f,x,x,f,f', skipline=slines+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+2*nbte[0]+1+1+1+1, numline=nbtb[0]
        btcl = xf_binning(tcl[*,4], binning+'_TE', lcen=lbin)
    endif
    if (ncl eq 6) or (strupcase(ncl) eq 'EB') then begin
        readcol,file,cl,cler,lmn,lmx,format='x,f,f,x,x,f,f', skipline=slines+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+2*nbte[0]+2*nbtb[0]+1+1+1+1+1, numline=nbeb[0]
        btcl = xf_binning(tcl[*,5], binning+'_EE', lcen=lbin)
    endif

    lcen = (lmx+lmn)/2.
    res = cl-btcl

    lmin = lmn
    lmax = lmx

    if (nbtb[0] gt 0) or ( nbeb[0] gt 0) then begin
        print, ' covmat meaningful for ncl=4 only. Skipped.'
    endif else begin
        if nbee[0] gt 0 then begin
            skipline = slines+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+2*nbte[0]+1+1+1
            openr,1,file
            tmpstr = ''
            for i=1,skipline do begin
                readf,1,tmpstr
;##                print, tmpstr
            endfor
            covmat = fltarr(nbins,nbins)
;##            help, covmat
            readf,1,covmat
            close,1
;##            help, covmat
;##            print, covmat[0:10,0:10]
        endif else begin
            skipline = slines+2*nbtt[0]
            openr,1,file
            tmpstr = ''
            for i=1,skipline do begin
                readf,1,tmpstr
;##                print, tmpstr
            endfor
            covmat = fltarr(nbins,nbins)
;##            help, covmat
            readf,1,covmat
            close,1
        endelse
    endelse

    return, cl
end
