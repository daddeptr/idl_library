; Simple ILC code
function ilc, mapfiles, maskfile=maskfile, gal_cut=gal_cut, gw=gw, bw=bw, silent=silent, check=check, $
              show_mask=show_mask, double=double, md_rem=md_rem, debug=debug, apply_smooth=apply_smooth, $
              mn=mn, mx=mx, imap=imap

  if ( not (keyword_set(maskfile)) and (not keyword_set(gal_cut)) ) then print, ' --> ILC performed on the whole sky'
  if not keyword_set(mn) then mn=-300
  if not keyword_set(mx) then mx=300
  if not keyword_set(imap) then imap=0

  true = 1b
  false = 0b

  do_png = false

  maptype = typename( mapfiles[0] )
  if not keyword_set( silent ) then print, maptype

;cmb2sz = fltarr(nfreq)
;h = 6.626068 * 10.^(-34)
;k = 1.3806503 * 10.^(-23)
;Tcmb = 2.725
;x = (h*freq*10.^9) / (k * Tcmb)
;cmb2sz = x * (exp(x)+1.)/(exp(x)-1.)-4.
;print, cmb2sz
;stop

    if ( strupcase( maptype ) eq 'STRING' ) then begin
        nfreq = n_elements(mapfiles)
        
        read_fits_map, mapfiles[0], map, nside=ns
        map = 0.
        npix = 12l*ns^2

        maps = fltarr(Npix,nfreq)
        if (keyword_set(double)) then maps = dblarr(Npix, nfreq)
;        mask = dblarr(Npix)+1
        print, 'reading maps...'
;        pmask = mask * 0.+1
        for ifreq = 0, nfreq-1 do begin
            print, 'map ', ifreq+1, ' of', nfreq, ' :'+mapfiles[ifreq]
            read_fits_map, mapfiles[ifreq], temp, order=dpord
;            bp = where( (temp eq -1.6375e30) or (finite(temp, /nan) eq 1) )
;            if (bp[0] ne -1) then temp[bp] = 0.
;            if (bp[0] ne -1) then pmask[bp] = 0.
            if (strupcase(dpord) eq 'NESTED') then begin
                print, ' - reordering '+mapfiles[ifreq], ' ', dpord, ' RING'
                temp = reorder(temp, in=dpord, out='ring') 
            endif
            if (imap lt 3) then begin
                maps[*,ifreq] = temp[*,imap] 
            endif else begin 
                bp = where( (temp[*,1] eq -1.6375e30) or (temp[*,2] eq -1.6375e30) )
                maps[*,ifreq] = sqrt( double(temp[*,1])^2+double(temp[*,2])^2 )
                maps[bp,ifreq] = -1.6375e30
            endelse
            if (keyword_set(check)) then mollview, maps[*,ifreq], /hist, px=600, win=ifreq
        endfor
    endif else begin
        print, ' - input maps are arrays: RING order assumed!'
        maps = mapfiles
        mapfiles = ''
        nfreq = n_elements( maps[0,*] )
        npix = n_elements(maps[*,0])
        ns = npix2nside(npix)
;stop
;        mask = dblarr(Npix)+1
    endelse

    if keyword_set(apply_smooth) then begin
        for ifreq =0, nfreq-1 do begin
            if keyword_set(check) then print, ' - smoothing map ', ifreq+1
            temp = maps[*,ifreq]
            ismoothing, temp, smth, /ring, fwhm_arcmin=apply_smooth[ifreq], simul_type=1
            maps[*,ifreq] = smth
            if keyword_set(check) then mollview, maps[*,ifreq], /hist, win=ifreq
        endfor
    endif
    mask = fltarr(Npix)+1
    if (keyword_set(double)) then mask = dblarr(Npix)+1
    pmask = mask * 0.+1
    for ifreq = 0, nfreq-1 do begin
;        print, 'map ', ifreq+1, ' of', nfreq
;        read_fits_map, mapfile[ifreq], temp, order=dpord
        temp = maps[*,ifreq]
        bp = where( (temp eq -1.6375e30) or (finite(temp, /nan) eq 1) )
        if (bp[0] ne -1) then temp[bp] = 0.
        if (bp[0] ne -1) then pmask[bp] = 0.
        maps[*,ifreq] = temp[*,0]
        if (keyword_set(check)) then mollview, temp, px=600, win=ifreq, min=-300, max=300. 
    endfor
    
    if ( keyword_set(maskfile) ) then begin
        if (strupcase(typename(maskfile)) eq 'STRING') then begin
            print, 'reading '+maskfile+'...'
            read_fits_map, maskfile, gmask, nside=maskns, order=dpord
            if ( strupcase(dpord) eq 'NESTED') then begin
                print, 'Reordering mask...'
                gmask=reorder(gmask, in=dpord, out='ring')
            endif
        endif else begin
            gmask = maskfile
            maskns = long( sqrt(n_elements(gmask)/12) )
        endelse
        mask = mask * gmask
        if (maskns ne ns) then begin
            print, 'Mismatch maps mask Nside', ns, maskns, '. Upgrading...'
            ud_grade, mask, masku, nside_out=ns, order_in='RING'
            mask = masku
            mask[where(mask ge 0.5)] = 1.
            mask[where(mask lt 0.5)] = 0.
        endif
    endif else begin
        mask = fltarr(npix) + 1.
        if (keyword_set(double)) then mask = dblarr(Npix)+1
    endelse
    
    if (keyword_set(gal_cut)) then mask = mask * make_sky_cut(gal_cut,ns) 
    if (keyword_set(show_mask) or keyword_set(check)) then mollview, mask, px=600, win=0, chars=1.5, tit='!6Mask used'
    
    if not keyword_set(silent) then print, ' - fsky = ', total(mask)/n_elements(mask)
;for ifreq = 0, nfreq-1 do maps[*,ifreq] = maps[*,ifreq] * cmb2sz[ifreq]

    if (keyword_set(md_rem)) then begin
        print, 'removing mono/dipole...'
        for ifreq=0,nfreq-1 do begin
            temp = maps[*,ifreq]
            remove_dipole, temp, mask, ordering='ring', nside=ns
            maps[*,ifreq] = temp
        endfor
    endif

    gpix = where( (mask[*,0] gt 0.) and (pmask gt 0.) )
    ngpix = n_elements(gpix)
    bpix = where(mask[*,0] eq 0.)
    nbpix = n_elements(bpix)


    gR = dblarr(nfreq, nfreq)
    bR = dblarr(nfreq, nfreq)
    print, ' - :DP - ILC: computing correlation matrix...'
;------ Pedestrian
if (False) then begin
    ave = dblarr(nfreq,2)

    for ifreq = 0, nfreq-1 do begin
        gavei = mean(maps[gpix, ifreq])
        ave[ifreq,0] = gavei
        if (bpix[0] ne -1) then begin
            bavei = mean(maps[bpix, ifreq])
            ave[ifreq,1] = bavei
        endif
        for jfreq = ifreq, nfreq-1 do begin
            avej = mean(maps[gpix,jfreq])
            gR[ifreq, jfreq] = total( (maps[gpix,ifreq]-gavei) * (maps[gpix,jfreq]-avej) ) / ngpix
;        gR[ifreq, jfreq] = total( (maps[gpix,ifreq]) * (maps[gpix,jfreq]) ) / ngpix
            gR[jfreq, ifreq] = gR[ifreq, jfreq]
            
            if (bpix[0] ne -1) then begin        
                avej = mean(maps[bpix,jfreq])
                bR[ifreq, jfreq] = total( (maps[bpix,ifreq]-bavei) * (maps[bpix,jfreq]-avej) ) / nbpix
;        bR[ifreq, jfreq] = total( (maps[bpix,ifreq]) * (maps[bpix,jfreq]) ) / nbpix
                bR[jfreq, ifreq] = bR[ifreq, jfreq]
            endif
        endfor
    endfor
endif
;------
   for ifreq = 0, nfreq-1 do begin
       if not ( keyword_set(silent) ) then print, ' - freq: ', ifreq+1
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

   if keyword_set( check ) then print, gR;, bR
   if ( keyword_set( check ) and (bpix[0] ne -1)) then print, bR

   gRm1 = invert(gR, /double, status)

   if keyword_set( check ) then print, gRm1;, bR

   print, ' - inversion status (gp): ', status
   if (bpix[0] ne -1) then begin
       bRm1 = invert(bR, /double, status)
       print, ' - inversion status (bp): ', status
   endif

   if ( keyword_set( check ) and (bpix[0] ne -1)) then print, bRm1
   a = findgen(nfreq) * 0. + 1.

   gw = dblarr(nfreq)
   gw = total(gRm1,2) / total(gRm1)
   if (bpix[0] ne -1) then bw = total(bRm1,2) / total(bRm1)

   print, ' - gw: ', gw
   if (bpix[0] ne -1) then print, ' - bw: ', bw

   ilc = fltarr(npix,3)
   if (keyword_set(double)) then ilc = dblarr(Npix, 3)

   for ifreq = 0, nfreq-1 do begin
;##	ilc[gpix] = ilc[gpix] + maps[gpix,ifreq] * gw[ifreq] 
;##        ilc[bpix] = ilc[bpix] + maps[bpix,ifreq] * bw[ifreq]
       ilc[*,0] = ilc[*,0] + maps[*,ifreq] * gw[ifreq] 
       if (bpix[0] ne -1) then ilc[*,1] = ilc[*,1] + maps[*,ifreq] * bw[ifreq]
   endfor

   ilc[gpix,2] = ilc[gpix,0]
   ilc[bpix,2] = ilc[bpix,1]

   print, 'STDDEV ILC gp (fs,out,in)       = ', stddev(ilc[*,0]), stddev(ilc[gpix,0]), stddev(ilc[bpix,0])
   print, 'STDDEV ILC bp (fs,out,in)       = ', stddev(ilc[*,1]), stddev(ilc[gpix,1]), stddev(ilc[bpix,1])
   print, 'STDDEV ILC combined (fs,out,in) = ', stddev(ilc[*,2]), stddev(ilc[gpix,2]), stddev(ilc[bpix,2])

;write_fits_map, 'ffp4_01a_ILCout.fits', ilc[*,0], /ring, units='!7l!6K'
;write_fits_map, 'ffp4_01a_ILCin.fits', ilc[*,1], /ring, units='!7l!6K'


   if (not keyword_set(silent)) then mollview, ilc[*,0], chars=1.5, tit='!6ILC: weights outside mask. N!dside!n='+string(ns,format='(i4.4)'), grat=[10,10], px=650, min=mn, max=mx ;, no_monopole=true, gal_cut=40, min=-300, max=300
   if (not keyword_set(silent) and (bpix[0] ne -1) ) then mollview, ilc[*,1], chars=1.5, tit='!6ILC: weights inside mask. N!dside!n='+string(ns,format='(i4.4)'), grat=[10,10], px=650, min=mn, max=mx ;, no_monopole=true, gal_cut=40, min=-300, max=300
   if (not keyword_set(silent) and (bpix[0] ne -1) ) then mollview, ilc[*,2], chars=1.5, tit='!6ILC: combined', grat=[10,10], px=650, min=mn, max=mx ;, no_monopole=true, gal_cut=40, min=-300, max=300
   if (not keyword_set(silent) and (bpix[0] ne -1) ) then mollview, ilc[*,0]-ilc[*,1], chars=1.5, tit='!6ILC Difference', grat=[10,10], px=650, min=mn/3, max=mx/3 ;, no_monopole=true, gal_cut=40, min=-15, max=15

   if (do_png) then begin
       restore, 'chains/pix_01a_v2.res.sav'
       mollview, cmb, min=-300, max=300, chars=1.5, win=4, tit='!6Commander', no_monopole=true, gal_cut=40
       mollview, cmb-ilc[*,0], min=-30, max=30, chars=1.5, win=5, tit='!6Commander-ILC!dout!n', no_monopole=true, gal_cut=40
       mollview, cmb-ilc[*,1], min=-30, max=30, chars=1.5, win=6, tit='!6Commander-ILC!din!n', no_monopole=true, gal_cut=40
       
       read_fits_map, 'ffp4_scalar_cmb_ns128_60arcmin_uK.fits', inp
       mollview, cmb-inp, min=-30, max=30, chars=1.5, win=7, tit='!6Commander-Input', no_monopole=true, gal_cut=40
       mollview, ilc[*,0]-inp, min=-30, max=30, chars=1.5, win=8, tit='!6ILC!dout!n-Input', no_monopole=true, gal_cut=40
       mollview, ilc[*,1]-inp, min=-30, max=30, chars=1.5, win=9, tit='!6ILC!din!n-Input', no_monopole=true, gal_cut=40
       
       mollview, ilc[*,0], min=-300, max=300, chars=1.5, win=-1, tit='!6ILC: weights outside mask', no_monopole=true, gal_cut=40, png='ffp4_01a_ILCout.png'
       mollview, ilc[*,1], min=-300, max=300, chars=1.5, win=-2, tit='!6ILC: weights inside mask', no_monopole=true, gal_cut=40, png='ffp4_01a_ILCin.png'
       mollview, ilc[*,0]-ilc[*,1], min=-30, max=30, chars=1.5, win=-3, tit='!6ILC Difference', no_monopole=true, gal_cut=40, png='ffp4_01a_ILCout-in.png'
       
       mollview, cmb, min=-300, max=300, chars=1.5, win=-4, tit='!6Commander', no_monopole=true, gal_cut=40, png='ffp4_01a_CMD.png'
       mollview, cmb-ilc[*,0], min=-30, max=30, chars=1.5, win=-5, tit='!6Commander-ILC!dout!n', no_monopole=true, gal_cut=40, png='ffp4_01a_CMD-ILCout.png'
       mollview, cmb-ilc[*,1], min=-30, max=30, chars=1.5, win=-6, tit='!6Commander-ILC!din!n', no_monopole=true, gal_cut=40, png='ffp4_01a_CMD-ILCin.png'

       mollview, cmb-inp, min=-30, max=30, chars=1.5, win=-7, tit='!6Commander-Input', no_monopole=true, gal_cut=40, png='ffp4_01a_CMD-INP.png'
       mollview, ilc[*,0]-inp, min=-30, max=30, chars=1.5, win=-8, tit='!6ILC!dout!n-Input', no_monopole=true, gal_cut=40, png='ffp4_01a_ILCout-INP.png'
       mollview, ilc[*,1]-inp, min=-30, max=30, chars=1.5, win=-9, tit='!6ILC!din!n-Input', no_monopole=true, gal_cut=40, png='ffp4_01a_ILCin-INP.png'
   endif

   if keyword_set(debug) then stop, 'DEBUG mode:'
   print, ' --- End of Program ---'

   return, ilc
;stop

; matrix multiplication failures
gw = grm1##a 
n = reform(a##(grm1##a))
gw[*] = gw[*] / n[0]

print, gw

bw = brm1##a
n = reform(a##(brm1##a))
bw[*] = bw[*] / n[0]


stop


end
