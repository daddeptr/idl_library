; Simple ILC code
function array_ilc, maps=maps, gw=gw, bw=bw, silent=silent, check=check, double=double, mask=mask

   true = 1b
   false = 0b

   nfreq = n_elements(maps[0,*])
   npix = n_elements(maps[*,0])
;   ns = npix2nside(npix)
;stop
   if not keyword_set(mask) then mask = dblarr(Npix)+1
   for ifreq = 0, nfreq-1 do begin
       temp = maps[*,ifreq]
       bp = where( (temp eq -1.6375e30) or (finite(temp, /nan) eq 1) )
       if (bp[0] ne -1) then temp[bp] = 0.
       if (bp[0] ne -1) then mask[bp] = 0.
       maps[*,ifreq] = temp[*,0]
   endfor
    
;    if (keyword_set(md_rem)) then begin
;        print, 'removing mono/dipole...'
;        for ifreq=0,nfreq-1 do begin
;            temp = maps[*,ifreq]
;            remove_dipole, temp, mask, ordering='ring', nside=ns
;        endfor
;    endif

    gpix = where( mask gt 0. )
    ngpix = n_elements(gpix)

    bpix = where( mask eq 0. )
    nbpix = n_elements(bpix)

    gR = dblarr(nfreq, nfreq)
    bR = dblarr(nfreq, nfreq)
    print, ' - :DP - ILC: computing correlation matrix...'

    for ifreq = 0, nfreq-1 do begin
        if not ( keyword_set(silent) ) then print, ifreq
        for jfreq = ifreq, nfreq-1 do begin
            gR[ifreq, jfreq] = correlate( maps[gpix,ifreq], maps[gpix,jfreq], /covariance )
            gR[jfreq, ifreq] = gR[ifreq, jfreq]
            
            if (bpix[0] ne -1) then begin        
                bR[ifreq, jfreq] = correlate( maps[bpix,ifreq], maps[bpix,jfreq], /covariance )
                bR[jfreq, ifreq] = bR[ifreq, jfreq]
            endif
        endfor
    endfor
; ------


    gRm1 = invert(gR, /double, status)
    print, status
    if (bpix[0] ne -1) then begin
        bRm1 = invert(bR, /double, status)
        print, status
    endif
    a = findgen(nfreq) * 0. + 1.

    gw = dblarr(nfreq)
    gw = total(gRm1,2) / total(gRm1)
    if (bpix[0] ne -1) then bw = total(bRm1,2) / total(bRm1)

    print, ' - gw: ', gw
    if (bpix[0] ne -1) then print, ' - bw: ', bw

    ilc = fltarr(npix,3)
    if (keyword_set(double)) then ilc = dblarr(Npix,3)

    for ifreq = 0, nfreq-1 do begin
        ilc[*,0] = ilc[*,0] + maps[*,ifreq] * gw[ifreq] 
        if (bpix[0] ne -1) then ilc[*,1] = ilc[*,1] + maps[*,ifreq] * bw[ifreq]
    endfor

    ilc[gpix,2] = ilc[gpix,0]
    if (bpix[0] ne -1) then ilc[bpix,2] = ilc[bpix,1]

    print, 'STDDEV ILC gp = ', stddev(ilc[gpix,0])
    if (bpix[0] ne -1) then print, 'STDDEV ILC bp = ', stddev(ilc[bpix,1])

    print, ' >>> array_ilc: --- End of Program ---'
    return, ilc
end
