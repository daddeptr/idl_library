
pro find_indices_peak, root=root, $
                       nchain=nchain, $
                       fst_ch=fst_ch, $
                       nforeg=nforeg, $
                       tag=tag, $
                       do_read=do_read, $
                       nband=nband, $
                       cls_samp=cls_samp, $
                       nsd=nsd, $
                       nsample=nsample, $
                       gsig=gsig, $
                       ncl=ncl, $
                       beams=beams, $
                       do_pol=do_pol, $
                       maskfile=maskfile, $
                       fst_samp=fst_samp, $
                       fg_pix=fg_pix, $
                       pl=pl, $
                       icmb=icmb, $
                       noise=noise, $
                       wmap=wmap


   close, /all
   loadct, 39
   set_plot, 'x'
   !p.multi = 0
   !p.background=255
   !p.color=0
   do_residuals = 0b
   if not ( keyword_set(wmap) ) then wmap = 'wmap'
   do_wmap = 0b
   if (wmap eq 'wmap') then do_wmap      = 1b
   print, ' do_wmap = ', do_wmap
   print, " type '.c' to continue"
;   stop

   do_corr      = 0b

   ns = long(nsd)
   nside = ns
   npix = 12l*ns^2
   s_ns = strtrim(string(ns),2)

   if not (keyword_set(beams)) then beams = 60.
   s_beams = strtrim(string(long(beams)),2)

    if not (keyword_set(root)) then begin
        print, 'root not defined'
        stop
    endif
    if not (keyword_set(tag)) then begin
        print, 'tag not defined'
        stop
    endif        

   if not (keyword_set(maskfile)) then maskfile='../sky_cut_30_ns'+s_ns+'_ring.fits'
   if not (keyword_set(nchain)) then nchain = 1 ;10l ;; 4l
   if not (keyword_set(fst_ch)) then fst_ch = 1
   if not (keyword_set(do_read)) then do_read = 'read' ;10l ;; 4l
   if not (keyword_set(nforeg)) then nforeg = 0
   if not (keyword_set(nband))  then nband = 7 
   if not (keyword_set(cls_samp))  then cls_samp = 'on'
   if not (keyword_set(gsig)) then gsig = 3.
;;    if not (keyword_set(nsample))  then nsample = -1l

   if not (keyword_set(ncl)) then ncl = 40l ;180l ;32l
   if not (keyword_set(do_pol)) then do_pol = 'no'
   if not (keyword_set(fst_samp)) then fst_samp = 1l
   if not (keyword_set(fg_pix)) then fg_pix = 't'
   if not (keyword_set(pl)) then pl = 't'
   if not (keyword_set(icmb)) then icmb = '0'
   if not (keyword_set(noise)) then noise = 'f'
   
; ------ Printing parameters                                                                                                                         
   openw, 1, 'find_indices_peak_'+tag+'_runs.txt'
   printf, 1, "cls_stats, maskfile="+ strtrim(maskfile,2)+ $
     ", wmap="    + strtrim(wmap)                        + $
     ", nsd="     + strtrim(string(nsd),2)               + $
     ", tag="     + strtrim(tag)                         + $
     ", root="    + strtrim(root)                        + $
     ", nchain="  + strtrim(string(nchain),2)            + $
     ", fst_ch="  + strtrim(string(fst_ch),2)            + $
     ", do_read=" + strtrim(do_read)                     + $
     ", nforeg="  + strtrim(string(nforeg),2)            + $
     ", nband="   + strtrim(string(nband),2)             + $
     ", cls_samp="+ strtrim(cls_samp)                    + $
     ", gsig="    + strtrim(string(gsig),2)              + $
     ", nsample=" + strtrim(string(nsample),2)           + $
     ", ncl="     + strtrim(string(ncl),2)               + $
     ", do_pol="  + strtrim(do_pol,2)                    + $
     ", fst_samp="+ strtrim(string(fst_samp),2)          + $
     ", fg_pix="  + strtrim(fg_pix,2)                    + $
     ", pl="      + strtrim(pl,2)                        + $
     ", icmb="    + strtrim(string(icmb),2)              + $
     ", noise="   + strtrim(string(noise),2)             + $
     ", beams="   + strtrim(string(beams),2)
   close, 1
   
; ------                                                                                                                                             

   if (noise eq 't') then begin
      do_noise = 1b
   endif else begin
      do_noise = 0b
   endelse

;   if (do_pol eq 'no') then begin
;       do_pol = 0b
;       nspec = 1
;   endif

   print, do_pol

   if (do_pol eq 'pol') then begin
       do_pol = 1b
       nspec = 6
   endif else begin  
       do_pol = 0b
       nspec = 1
   endelse

   if (do_read eq 'read') then begin
   if (fg_pix eq 't') then begin
       fg_pix = 1b
       spawn, 'ls -l '+root+'fg_amp_map_no*_c0001_k00001.fits | wc -l', namp
       namp = reform(long(namp))
       print, ' amplitude maps # = ', namp
       spawn, 'ls -l '+root+'fg_ind_map_no*_c0001_k00001.fits | wc -l', nind
       nind = reform(long(nind))
       print, ' index maps # = ', nind
       if not (do_pol) then amp = fltarr(npix, namp[0])
       if (do_pol) then amp = fltarr(npix, 3, namp[0])
       amp2 = amp
;;        amp1 = fltarr(npix)
;;        amp2 = fltarr(npix)
;;        amp3 = fltarr(npix)
       ind  = 0.
       ind2 = 0.
       if (nind gt 0) then begin
          ind = fltarr(npix, nind[0])
          ind2 = ind
       endif
;;        ind1 = fltarr(npix)
;;        ind2 = fltarr(npix)
;;        ind3 = fltarr(npix)
;;        ind4 = fltarr(npix)
   endif else begin
       fg_pix = 0b
   endelse
   endif

   if (pl eq 't') then begin
       pl = 1b
   endif else begin
       pl = 0b
   endelse

   if (nforeg ge 0) then nforeg = nforeg + 4

   print, ' nforeg : ', nforeg
   print, ' nchain : ', nchain
   print, ' fst_ch : ', fst_ch
   print, ' do_read: ', do_read
   print, ' cls_samp: ', cls_samp
   print, ' ncl:      ', ncl
   print, ' gsig:     ', gsig
   print, ' beams:    ', beams
   print, ' do_residuals', do_residuals
   print, ' do_wmap     ', do_wmap
;;    print, ' do_corr     ', do_corr
   print, ' do_pol      ', do_pol
   print, ' nband:    ', nband, ' <== !!! ==> '
   print, ' icmb:      ', icmb

;;    freq        = ['30', '44', '70', '100', '143', '217', '353']
;;    freq        = ['30', '44', '70', '100', '143', '217', '353', '33', '23', '41', '61', '94']
   freq        = [30, 44, 70, 100, 143, 217, 353, 33, 23, 41, 61, 94]
   freq = string(freq, format='(i3.3)')

   nfreq       = n_elements(freq)
   col_array   = 30 + lindgen(nfreq)*250/nfreq ;; [200, 205, 210, 215, 220, 225, 230]
   line_array  = col_array * 0
   spectra     = ['TT', 'EE', 'BB', 'TE', 'TB', 'EB']

   read_fits_map, maskfile, mask

   gpix = where(mask[*,0] gt 0.)
   ngpix = n_elements(gpix)
   bpix = where(mask[*,0] eq 0.)
   bad_value = -1.6375e30
   fsky = float(ngpix) / n_elements(mask[*,0])
   print, " - fsky = ", fsky

   if (do_read eq 'read') then begin
       nsampch = lonarr(nchain)
       if (nsample lt 0l) then begin
           nsample = 0l
           iich = 0l
           for ich=fst_ch,fst_ch+nchain-1 do begin
               sch = 'c'+string(ich,format='(1i4.4)')
               spawn, 'ls -ltrh '+root+'chisq_'+sch+'_k* | wc -l', nsamp
               nsampch[iich] = long(nsamp)
               nsample = nsample + nsampch[iich]
               iich = iich + 1l
           endfor
       endif else if (nsample gt 0l) then begin 
           nsampch[*] = nsample / nchain
       endif
   
   print, ' nsample:      ', nsample
   print, ' nsampch:      ', nsampch
   print, ' do_residuals: ', do_residuals

   if (min(nsampch) lt fst_samp) then begin
       print, " Not enough samples..."
       stop
   endif

   print, nsample
   nsample = nsample - nchain*(fst_samp-1)
   print, fst_samp, nchain
   print, nsample

   ind_arr = fltarr(npix, reform(nind[0]), nsample)

;;    cl=fltarr(nsample,nspec,2 * (ncl+1))

   l = findgen(ncl+1)

   if (nforeg gt 0) then foreg = fltarr(nsample, nforeg, nband)

   chi2 = dblarr(nsample)

   if (cls_samp eq 'on') then begin
      if not (do_pol) then cmb = dblarr(npix)
      if (do_pol) then cmb = dblarr(npix,3) 
      cmb2 = cmb
;;       glb_chi2 = cmb * 0.
   endif else begin
      if not (do_pol) then cmb = dblarr(npix)
      if (do_pol) then cmb = dblarr(npix,3)
      cmb2 = cmb
;;       glb_chi2 = reform( amp[*,0,*] ) * 0.
   endelse

;;    freq = ['030', '044', '070', '100', '143', '217', '353']
;;    freq        = ['30', '44', '70', '100', '143', '217', '353', '33', '23', '41', '61', '94']

;;    if (do_residuals) then residuals = fltarr(npix,nband,nsample)

   iich = -1l
   tot_isamp = -1l

   warn = 1b

   if (do_noise) then begin
       if (cls_samp eq 'off') then read_fits_map, root+'fg_amp_map_no'+string(icmb+1,format='(i2.2)')+'_c'+string(fst_ch,format='(i4.4)')+'_k'+string(fst_samp,format='(i5.5)')+'.fits', ref_samp
       if (cls_samp eq 'on') then read_fits_map, root+'s_i_'+sch+'_'+ssam+'.fits', ref_samp 
       nl = fltarr(nsample,ncl+1)
   endif

   for ich=fst_ch, fst_ch+nchain-1 do begin
      sch = 'c'+string(ich,format='(1i4.4)')
      iich = iich + 1

      for isam=fst_samp, nsampch[iich] do begin
;         cl_tmp = fltarr(nspec+1, 2*(ncl+1))
;         if (nforeg gt 0) then fg = fltarr(nforeg, nband)
         tot_isamp = tot_isamp+1
         ts = tot_isamp / 100
         if ( (ts*100) eq tot_isamp ) then begin
             print, tot_isamp
             if not (do_residuals) then save, filename='find_indices_peak_'+tag+'.res.sav', foreg, cmb, cmb2, chi2, nsample, nsampch, nside, tot_isamp, amp, amp2, ind, ind2, namp, nind, ncl, fst_samp, nchain, fst_ch, ind_arr
             if (do_residuals) then save, filename='find_indices_peak_'+tag+'.res.sav', foreg, cmb, cmb2, chi2, residuals, nsample, nsampch, nside, tot_isamp, amp, amp2, ind, ind2, namp, nind, ncl, fst_samp,nchain,fst_ch, ind_arr
         endif

         ssam = 'k'+string(isam,format='(1i5.5)')

;         read_fits_map, root+'chisq_'+sch+'_'+ssam+'.fits', chisq
;
;         glb_chi2 = glb_chi2 + chisq

         if (cls_samp eq 'on') then begin
;             openr,1,root+'cl_'+sch+'_'+ssam+'.dat'
;             readf,1,cl_tmp
;             close, 1
             
;             cl[tot_isamp,*,*] = reform(cl_tmp[1:*,*])

             read_fits_map, root+'s_i_'+sch+'_'+ssam+'.fits', cmbmap
             cmb = cmb + cmbmap
             cmb2 = cmb2 + cmbmap^2

         endif

;;          chi2[tot_isamp] = total(chisq) / ngpix / nband

         if (cls_samp eq 'off') then begin
             if (warn) then begin
                 print, " WARNING: cls_off polarization not implemented yet!!!"
                 warn = 0b
             endif
             
             read_fits_map, root+'fg_amp_map_no'+string(icmb+1, format='(i2.2)')+'_'+sch+'_'+ssam+'.fits', cmbmap
             cmb = cmb + cmbmap
             cmb2 = cmb2 + cmbmap^2
         endif

;         if (do_noise)  then begin
;             nmap = (ref_samp - cmbmap) / sqrt(2.)
;             ianafast, nmap, nls, maskfile=maskfile, /ring, nlmax=ncl, /silent
;             nl[tot_isamp,*] = nls
;         endif


;         if (nforeg gt 0) then begin
;             openr, 1, root+'foreground_'+sch+'_'+ssam+'.dat'
;             readf,1,fg
;             close,1
;             foreg[tot_isamp,*,*] = fg
;         endif

         if (fg_pix) then begin
;             for iamp = 1, namp[0] do begin
;                read_fits_map, root+'fg_amp_map_no'+string(iamp,format='(i2.2)')+'_'+sch+'_'+ssam+'.fits', map
;                amp[*,iamp-1] = amp[*,iamp-1] + map
;                amp2[*,iamp-1] = amp2[*,iamp-1] + map^2
;             endfor
             
             if (nind gt 0) then begin
                for iind =1, nind[0] do begin
                   read_fits_map, root+'fg_ind_map_no'+string(iind,format='(i2.2)')+'_'+sch+'_'+ssam+'.fits', map
                   ind[*,iind-1] = ind[*,iind-1] + map
                   ind2[*,iind-1] = ind2[*,iind-1] + map^2
                   ind_arr[*,iind-1,tot_isamp] = map
                endfor
             endif

         endif

;         if (do_residuals) then begin
;             for iband=0,nband-1 do begin
;                 read_fits_map, '../map_freefree_0'+freq[iband]+'.00GHz.smth.ns'+s_ns+'.fits', ff
;                 read_fits_map, '../map_synchrotron_0'+freq[iband]+'.00GHz.smth.ns'+s_ns+'.fits', ss
;                 read_fits_map, '../map_thermal_dust_0'+freq[iband]+'.00GHz.smth.ns'+s_ns+'.fits', dd
;                 read_fits_map, '../noDip_dx4_'+freq[iband]+'GHz_I_smth600_uK9_ns'+s_ns+'.fits', imap
;                 dipole = make_dipole(reform(fg[1:3,iband]),ns)
;                 residuals[*,iband,tot_isamp] = (imap - (map + fg[0,iband] + dipole + fg[4,iband]*ff + fg[5,iband]*ss + fg[6,iband]*dd))              
;             endfor
;         endif


     endfor
 endfor

 print, nsample

 if (cls_samp eq 'on') then begin
    cmb = cmb / (nsample)
    cmb2 = cmb2 / (nsample)
 endif
;;    glb_chi2 = glb_chi2 / nsample

   if (fg_pix) then begin
;       amp = amp / nsample
;       amp2 = amp2 / nsample
;       ampstd=sqrt(amp2-amp^2)

       ind = ind / nsample
       ind2 = ind2 / nsample

       indstd=sqrt(ind2-ind^2)
   endif


   if not (do_residuals) then save, filename='find_indices_peak_'+tag+'.res.sav', foreg, cmb, cl, chi2, cmb2, nsample, nsampch, nside, amp, amp2, ind, ind2, namp, nind, glb_chi2, ncl, mask, nl,fst_samp, nchain,fst_ch, ind_arr
   if (do_residuals) then save, filename='find_indices_peak_'+tag+'.res.sav', foreg, cmb, cl, chi2, cmb2, residuals, nsample, nsampch, nside, amp, amp2, ind, ind2, namp, nind, glb_chi2, ncl, mask, nl,fst_samp, nchain,fst_ch, ind_arr
   endif

   if (do_read ne 'read') then begin
        print, "restoring..."
;;         freq = ['030', '044', '070', '100', '143', '217', '353']
        freq        = ['030', '044', '070', '100', '143', '217', '353', '033', '023', '041', '061', '094']
        nfreq = n_elements(freq)
	restore, filename='find_indices_peak_'+tag+'.res.sav'
        npix = 12l*nside^2
        map = fltarr(npix)
        fg_pix = 0b
        if (fg_pix eq 't') then fg_pix = 1b
   endif

   peak_ind = ind*0.

   print, 'Looping'
   window, 3 & plot, /nodata, [1,3], [0,100], chars=1.5
   for ipix=0l,npix-1 do begin
       if ( (ipix/10000l)*10000l eq ipix) then print, ipix
       for iind=0,nind[0]-1 do begin
           h = histogram(ind_arr[ipix,iind,*], nbins=7, locations=a)
           m = max(h,im)
;print, m, im
           if ( ( (ipix/5000l)*5000l eq ipix) and (iind eq 0) ) then oplot, a, h, psym=10
           peak_ind[ipix,iind] = a[im]
;stop
       endfor
   endfor


   write_fits_map,tag+'_peak_emissivity.fits', peak_ind[*,0],/ring
   write_fits_map,tag+'_peak_t.fits', peak_ind[*,1],/ring
   write_fits_map,tag+'_peak_beta.fits', peak_ind[*,2],/ring

   for iind=0,nind[0]-1 do begin
       mollview, ind[*,iind], tit='Mean'
       mollview, peak_ind[*,iind], tit='Peak'
   endfor
   print,"END OF NEW PIECE"

stop
   if (do_noise) then begin
      avenls = total(nl,1)/nsample
      cl2fits, avenls, tag+'_avenls.fits'
   endif

   mean_md = total(foreg,1) / nsample

   !p.background=255
   !p.color=0
   !p.multi = 0
   !p.charsize=1.5
   window, 11 & plot, chi2, yr=[0.9*min(chi2),1.1*max(chi2)], ys=1, ytit='!7v!u2!n', xtit='!17Sample'
   chihst =  histogram(chi2, nbins=50, locations=achi);, min=mean(chi2)-4.*stddev(chi2), max=mean(chi2)+4.*stddev(chi2) )
   window, 12 & plot, achi, chihst/total(chihst), ys=2, xtit='!7v!u2!n',ytit='!17P', psym=10;, xtit='!17Sample'

   set_plot, 'ps'
   device, file=tag+'_chi2.eps',/col, bits=8
   plot, chi2, yr=[0.9*min(chi2),1.1*max(chi2)], ys=1, ytit='!7v!u2!n', xtit='!17Sample'
   device, /close
; ---
   device, file=tag+'_chi2_hist.eps',/col, bits=8
   plot, achi, chihst/total(chihst), ys=2, xtit='!7v!u2!n',ytit='!17P', psym=10
   device, /close
   set_plot, 'x'

   
   mxchi2 = max(chihst,im)
   mxchi2 = achi[im]
   mnchi2 = achi[0]
   print, mxchi2, mnchi2
   print, n_elements(where(chi2 gt 2*mxchi2-mnchi2)) / nchain
;
;   cmb = cmb / (nsample)
;   cmb2 = cmb2 / (nsample)
;   glb_chi2 = glb_chi2 / nsample
;
;   if (fg_pix) then begin
;       amp1 = amp1 / nsample
;       amp2 = amp2 / nsample
;       ind1 = ind1 / nsample
;       ind2 = ind2 / nsample
;       ind3 = ind3 / nsample
;   endif

   if (do_residuals) then begin 
       res = total(residuals,3) / (nsample-nchain*fst_samp)
       mollview, res[*,3]*mask, min=-100, max=100
       if (nband eq 5) then resw = create_struct('HDR','Residuals','GHz '+freq[0],reform(residuals[*,0]),'GHz '+freq[1],reform(residuals[*,1]),'GHz '+freq[2],reform(residuals[*,2]),'GHz '+freq[3],reform(residuals[*,3]),'GHz '+freq[4],reform(residuals[*,4]) )
       if (nband eq 6) then resw = create_struct('HDR','Residuals','GHz '+freq[0],reform(residuals[*,0]),'GHz '+freq[1],reform(residuals[*,1]),'GHz '+freq[2],reform(residuals[*,2]),'GHz '+freq[3],reform(residuals[*,3]),'GHz '+freq[4],reform(residuals[*,4]),'GHz '+freq[5],reform(residuals[5]) )
;stop
   endif


   if (cls_samp eq 'on') then begin
      cmbv = sqrt(cmb2-cmb^2)

      if not (do_pol) then begin
         write_fits_map, tag+'_recCMB_ns'+s_ns+'_RING.fits',cmb, /ring,units='!7l!17K'
         write_fits_map, tag+'_rmsCMB_ns'+s_ns+'_RING.fits',cmbv, /ring,units='!7l!17K'
      endif

      if (do_pol) then begin
         write_tqu, tag+'_recCMB_ns'+s_ns+'_RING.fits',cmb, /ring,units='!7l!17K'
         write_tqu, tag+'_rmsCMB_ns'+s_ns+'_RING.fits',cmbv, /ring,units='!7l!17K'
      endif
   endif

   if (cls_samp eq 'off') then begin

      if (1b) then begin
         if (do_pol) then begin
            cmb = reform( amp[*,0,icmb] )
            cmb2 = reform( amp2[*,0,icmb] )
         endif
         if not (do_pol) then begin
            cmb = reform( amp[*,icmb] )
            cmb2 = reform( amp2[*,icmb] )
         endif

         cmbv = sqrt(cmb2-cmb^2)

         extr = max([max(cmb*mask),-min(cmb*mask)])
         if (extr gt 200) then extr = 200.
         remove_dipole, cmb, mask[*,0], ordering='ring', nside=ns
         write_fits_map, tag+'_recCMB_ns'+s_ns+'_RING.fits',cmb, /ring,units='!7l!17K'
         write_fits_map, tag+'_rmsCMB_ns'+s_ns+'_RING.fits',cmbv, /ring,units='!7l!17K'
         ccc = cmb
         ccc[bpix] = bad_value
         mollview, ccc, min=-extr, max=extr, chars=1.5, tit='!17Reconstructed CMB: '+tag
         mollview, ccc, min=-250, max=250, chars=1.5, tit='!17dx5 CMB', units='!7l!17K thermo', ps=tag+'_dx5_cmb.eps'

         ccc = cmbv
         ccc[bpix] = bad_value
         mollview, ccc, min=min(cmbv[where(mask gt 0.)]*mask[where(mask gt 0.)]), chars=1.5, tit='!17Reconstructed CMB sigma: '+tag
      endif else begin
         print, ' cls_samp = off: CMB map not saved!'
      endelse
   endif


   if (cls_samp eq 'on') then begin
       print, " smoothing signal to FWHM ", beams
       ismoothing,tag+'_recCMB_ns'+s_ns+'_RING.fits',tag+'_recCMB_smth'+s_beams+'_ns'+s_ns+'_RING.fits', fwhm_arc=beams, /silent
       ismoothing,tag+'_rmsCMB_ns'+s_ns+'_RING.fits',tag+'_rmsCMB_smth'+s_beams+'_ns'+s_ns+'_RING.fits', fwhm_arc=beams, /silent

       read_fits_map, tag+'_recCMB_smth'+s_beams+'_ns'+s_ns+'_RING.fits', smcmb
       read_fits_map, tag+'_rmsCMB_smth'+s_beams+'_ns'+s_ns+'_RING.fits', smcmbv
       extr = max([max(smcmb*mask),-min(smcmb*mask)])
       if (extr gt 300) then extr = 300.

       ccc = smcmb
       ccc[bpix] = bad_value
       mollview, ccc, min=-extr, max=extr, chars=1.5, tit='!17Reconstructed CMB: '+tag

       ccc = smcmbv
       ccc[bpix] = bad_value
       mollview, ccc, min=min(smcmbv[where(mask gt 0.)]*mask[where(mask gt 0.)]), chars=1.5, tit='!17Reconstructed CMB sigma: '+tag

if (do_wmap) then begin
       tempfile = '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/wmap/wmap7_ilc_smth'+s_beams+'_ns'+s_ns+'_ring.fits'
       read_fits_map,tempfile,ilc
       ilc = ilc*1.e3
       ccc = ilc
       ccc[bpix] = bad_value
       mollview, ccc, chars=1.5, tit='!17WMAP7 ILC', min=-200, max=200, units='!7l!17K', ps='wilc.eps'
;;        mollview, ccc, chars=1.5, tit='!17WMAP7 ILC', min=-200, max=200, units='!7l!17K', png='wilc.png', window=-1

       ccc = smcmb
       ccc[bpix] = bad_value
       mollview, ccc, chars=1.5, tit='!17CMD Mean map', min=-200, max=200, units='!7l!17K', ps=tag+'_cmd_mean.eps'
;;        mollview, ccc, chars=1.5, tit='!17CMD Mean map', min=-200, max=200, units='!7l!17K', png=tag+'_cmd_mean.png', window=-1

       ccc = ilc-smcmb
       ccc[bpix] = bad_value
       mollview, ccc, chars=1.5, tit='!17WMAP7 ILC - CMB', min=-40, max=40, units='!7l!17K', ps=tag+'_ilc-cmd.eps'
;;        mollview, ccc, chars=1.5, tit='!17WMAP7 ILC - CMB', min=-40, max=40, units='!7l!17K', png=tag+'_ilc-cmd.png', window=-1
       mollview, ccc, chars=1.5, tit='!17WMAP7 ILC - CMB', min=-40, max=40, units='!7l!17K'
       ianafast, tempfile,ilc_cl, maskfile=maskfile, /silent, regression=2

       ianafast, tag+'_recCMB_smth'+s_beams+'_ns'+s_ns+'_RING.fits', cmb_cl, maskfile=maskfile, /silent, regression=2
print, ' needlet ILC comparison commented'
;;        read_fits_map,'needlet_ILC_dpc_ns128_smth90_ring.fits', needilc
;;        ccc = needilc*1.e6-smcmb
;;        ccc[bpix] = bad_value
;;        mollview, ccc, min=-50, max=50, chars=1.5, units='!7l!17K',tit='!17Needlet ILC - COMMANDER CMB', ps=tag+'_needilc-cmd.eps'
endif else begin 
;;        tempfile = '../map_cmb_smth'+s_beams+'_ns'+s_ns+'_ring.fits'
;;        read_fits_map,tempfile,ilc
;;        ilc = ilc * 1.e6
;;        ccc = ilc
;;        ccc[bpix] = bad_value
;;        mollview, ccc, chars=1.5, tit='!17Input map', min=-200, max=200, units='!7l!17K', ps='inmap.eps'
;; 
;;        ccc = smcmb
;;        ccc[bpix] = bad_value
;;        mollview, ccc, chars=1.5, tit='!17CMD Mean map', min=-200, max=200, units='!7l!17K', ps=tag+'_cmd_mean.eps'
;; 
;;        ccc = ilc-smcmb
;;        ccc[bpix] = bad_value
;;        mollview, ccc, chars=1.5, tit='!17Inpu - CMB', min=-40, max=40, units='!7l!17K', ps=tag+'_inp-cmd.eps'
;;        mollview, ccc, chars=1.5, tit='!17Input - CMB', min=-40, max=40, units='!7l!17K'
;; 
;;        ianafast, tempfile,ilc_cl, maskfile=maskfile, /silent

endelse

   endif

   cc = glb_chi2
   cc[bpix] = bad_value
   mollview, cc, chars=1.5, tit='!7v!u2!17', max=5*nband
   mollview, cc, chars=1.5, tit='!7v!u2!17', max=5*nband, ps=tag+'_glb_chi2.eps'

   if (fg_pix) then begin
       for iamp=1,namp[0] do begin 
           mollview, amp[*,iamp-1], chars=1.5, tit='!17Amp '+string(iamp, format='(i2.2)'), max=300, min=-300
       endfor
       for iind=1,nind[0] do begin
           mollview, ind[*,iind-1], chars=1.5, tit='!17Ind '+string(iind, format='(i2.2)')
       endfor



;;        read_fits_map,'thermal_dust_353GHz_ns128_smth90.fits',dd
;;        hfc = dd * conversionfactor(353.,/thermo2antenna)
;;        mollview, hfc, chars=1.5, /hist, tit='!17FFP3 Thermal Dust Model @ 353GHz', ps='ffp3_thermal_dust_353.eps'

;;        read_fits_map,'freefree_030GHz_ns128_smth90.fits',ff
;;        read_fits_map,'synchrotron_030GHz_ns128_smth90.fits',ss
;;        lfc = (ff+ss)*conversionfactor(30.,/brightness2muKthermo)
;;        mollview, lfc, chars=1.5, /hist, tit='!17FFP3 Synch & f-f Model @ 30GHz', ps='ffp3_sync-ff_030.eps'
   endif
if not (do_wmap) then begin
    if (0b) then begin
       read_fits_map,'../map_cmb_smth'+s_beams+'_ns'+s_ns+'_ring.fits', inp
       inp = inp * 1.e6
       cc = inp
       cc[bpix] = bad_value
       mollview, cc, tit='!17Input CMB', chars=1.5, min=-200, max=200, ps=tag+'_inpmap.eps'


       if (cls_samp eq 'on') then cc = smcmb
       if (cls_samp eq 'off') then cc = cmb
       cc[bpix]= bad_value
       mollview, cc, tit='!17Recovered CMB', chars=1.5, min=-200, max=200, ps=tag+'_cmbmap.eps'

       if (cls_samp eq 'on') then t = inp - smcmb
       if (cls_samp eq 'off') then t = inp - cmb
       t[bpix] = bad_value
;;        remove_dipole, t, nside=ns, ordering='ring'
       mollview, t, chars=1.5, min=-max(inp)/10, max=max(inp)/10, units='!7l!17K thermo', tit='!17Input CMB - Commander CMB', ps=tag+'_CMB_residuals.eps'
       mollview, t, chars=1.5, min=-max(inp)/10, max=max(inp)/10, units='!7l!17K thermo', tit='!17Input CMB - Commander CMB'
print, ' Foreground modeling...'
stop
       fg_model = fltarr(npix, nband)
       inmaps = fltarr(npix, nband)
       convf = conversionfactor(freq, /antenna2thermo)
       print, convf

       h = 6.626068 * 10.^(-34)
       k = 1.3806503 * 10.^(-23)
       c = 3.d8

       

       mon = total(foreg[*,0,*], 1) / nsample
;;      help, mon
       dip = total(foreg[*,1:3,*], 1) / nsample 
;;        help, dip

       print, ' !!! WARINING: template reading removed!!!'
;;        read_fits_map, '../Halpha_90_ns'+s_ns+'_ring.fits', halpha
;;        read_fits_map, '../haslam408_dsds_90_ns'+s_ns+'_ring.fits', haslam
;;        read_fits_map, '../fds_94_90_ns'+s_ns+'_ring.fits', fds

       if not (fg_pix) then begin
           f_amp = total(foreg[*,4,*], 1) / nsample
           s_amp = total(foreg[*,5,*], 1) / nsample
           d_amp = total(foreg[*,6,*], 1) / nsample
       endif

       for ifreq = 0, nband-1 do begin
           print, ' --- '
           print, freq[ifreq]

           if not (fg_pix) then fg_model[0:*, ifreq] = f_amp[0,ifreq] * halpha + s_amp[0,ifreq]*haslam + d_amp[0,ifreq]*fds

           if (fg_pix) then begin
               if (pl) then fg_model[0:*,ifreq] = convf[ifreq] * ( amp1*(float(freq[ifreq])/353.)^ind1 + amp2*(float(freq[ifreq])/30.)^ind2 )
               if not (pl) then begin
                   x = h*freq[ifreq]*1.d9/k/ind3
;;                    mollview, x
                   bb = 1. / (exp(x)-1)
;;                    mollview, bb
                   fg_model[0:*,ifreq] = convf[ifreq] * (amp2*(float(freq[ifreq])/359.919)^(ind2)*bb + amp2*(float(freq[ifreq])/28.45733)^(ind1) + amp3*(float(freq[ifreq])/30.)^(-2.15));;; + mon[0,ifreq] + make_dipole(dip[*,ifreq], ns)
               endif
           endif
           cc = fg_model[*,ifreq]
           cc[bpix] = bad_value

           mollview, cc, chars=1.5, tit='!17Foreg Model '+freq[ifreq], units='!7l!17K thermo', /asinh, min=-100, max=40000, colt=3, ps=tag+'_fg_model_'+freq[ifreq]+'.eps'
;;            mollview, cc, chars=1.5, tit='!17Foreg Model '+freq[ifreq], units='!7l!17K thermo', /asinh, min=-100, max=40000, colt=3

;;            t = fg_model[*,ifreq]
;;            vec0 = reform( mean_md[1:3,ifreq] )
;;            dip = make_dipole(vec0,ns)
;;            mollview, dip+mean_md[0,ifreq], chars=1.5, tit='Dipole: '+freq[ifreq]

;;            remove_dipole, t, mask[*,0], nside=ns, ordering='ring', gal_cut=30
;;            mollview, t, /hist
;stop
;
           read_fits_map, '../ffp3_'+freq[ifreq]+'.fits', map
           inmaps[*,ifreq] = map[*,0] ; --- Temperature olny
           cc = inmaps[*,ifreq] - inp
           cc[bpix] = bad_value
           mollview, cc, chars=1.5, tit='!17Input foreground '+freq[ifreq], units='!7l!17K thermo', /asinh, min=-100, max=40000, colt=3, ps=tag+'_input_fg_'+freq[ifreq]+'.eps'
;;            read_fits_map, '../clean_map_'+freq[ifreq]+'_uK_ns128_90_ring.fits', map

           
           if (cls_samp eq 'on') then cc = inmaps[*,ifreq] - fg_model[*,ifreq] - smcmb - mon[0,ifreq] - make_dipole(reform(dip[*,ifreq]), ns, silent='f')
           if (cls_samp eq 'off') then cc = inmaps[*,ifreq] - fg_model[*,ifreq] - cmb - mon[0,ifreq] - make_dipole(reform(dip[*,ifreq]), ns, silent='f')
           cc[bpix] = bad_value
           print, ' monopole: ', mon[0,ifreq]
           print, ' dipole:   ', dip[*,ifreq]
;;            remove_dipole, cc, nside=ns, ordering='ring'
           mollview, cc, chars=1.5, tit='Input Fg - Cmd model: '+freq[ifreq], units='!7l!17K thermo', min=-40, max=40, ps=tag+'_fg_residuals_'+freq[ifreq]+'.eps'
           mollview, cc, chars=1.5, tit='Input Fg - Cmd model: '+freq[ifreq], units='!7l!17K thermo', /hist;, min=-40, max=40
stop
           cc = inmaps[*,ifreq]-inp[*,0]
           cc[bpix] = bad_value
           mollview, cc, chars=1.5, tit='Input Fg: '+freq[ifreq], units='!7l!17K thermo', /hist, min=0, max=500                      
stop
           cc = fg_model[*,ifreq] - mon[0,ifreq] - make_dipole(reform(dip[*,ifreq]), ns, silent='f')
           cc[bpix] = bad_value
           mollview, cc, chars=1.5, tit='Cmd model: '+freq[ifreq], units='!7l!17K thermo', /hist , min=0, max=500
print, ' Next frequency'
stop
;           c = map - t - cmb - dip - mean_md[0, ifreq]
;           remove_dipole, c, nside=ns, ordering='ring', gal_cut=30
;           mollview, c*mask, chars=1.5, /hist, tit=freq[ifreq]

;;            cc = map - fg_model[*,ifreq]
;;            cc[bpix] = bad_value
;;            mollview, cc, chars=1.5, tit='CMB: '+freq[ifreq], units='!7l!17K thermo'
;
; print, 'Are you sure you want to rewrite clean_map? '
;;            write_fits_map,'../clean_map_'+freq[ifreq]+'_uK_ns128_90_ring.fits', map-t, /ring, units='!7l!17K'

       endfor

       window, 22
       !p.background=255
       !p.color=0
       !p.charsize=1.5
       plot, (inmaps[*,6]-inp)*mask[*,0], fg_model[*,6]*mask[*,0], psym=1, xtit='!17Input pixel value', ytit='!17Commander Fg model', xr=[-1000, 10000], yr=[-1000, 10000]
       oplot, (inmaps[*,0]-inp)*mask[*,0], fg_model[*,0]*mask[*,0], psym=2, col=245
       oplot, fg_model[*,6]*mask[*,0] , fg_model[*,6]*mask[*,0], col=70
       legend,['30 GHz', '353 GHz'], line=[0,0], col=[245, 0]

;;        set_plot, 'ps'
;;        device, file=tag+'_fg_input_VS_output.eps', /col, bits=8
;;        plot, (inmaps[*,6]-inp)*mask[*,0], fg_model[*,6]*mask[*,0], psym=1, xtit='!17Input pixel value', ytit='!17Commander Fg model', xr=[-1000, 10000], yr=[-1000, 10000]
;;        oplot, (inmaps[*,0]-inp)*mask[*,0], fg_model[*,0]*mask[*,0], psym=2, col=245
;;        oplot, fg_model[*,6]*mask[*,0] , fg_model[*,6]*mask[*,0], col=70
;;        legend,['30 GHz', '353 GHz'], line=[0,0], col=[245, 0]
;;        device, /close
;;        set_plot, 'x'
stop
;                                                                                                                                       

endif

endif else begin
    print, ' ==> WARNING: input comparison skipped <== '
endelse


   fits2cl, pwf, '../pixel_window_n'+string(ns,format='(1i4.4)')+'.fits'

;;    ifreq = 2 ;70GHz
;;    gnoise = gsig*randomn(-4+ifreq,npix)
;;    print, " computing example noise cls..."
;;    ianafast, gnoise, ncls, nlmax= ncl, ordering='ring', /silent

   gb = gaussbeam(beams,ncl)
   l = lindgen(ncl+1)
   ll = l*(l+1)/2./!pi	

   bins = max([5, nsample / 50 / nband]) ;75 ;13 ;25 
   hist=fltarr(bins,nband)
   famp = fltarr(bins, nband)
   !p.background=255
   !p.color=0

   if (nforeg gt 0) then begin 

       if (nforeg gt 4) then begin
           for i=0,nband-1 do begin
               hist[*,i]=histogram(foreg[*,4,i],nbins=bins,locations=a)
;;       hist[*,i]=histogram(foreg[0,*,i,4],nbins=bins,locations=a)
               famp[*,i] = a
           endfor
           window, 1 & plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist)/total(hist[*,0])], xtit='!17A!dff!n', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
           for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
;   freq = ['30', '44', '70', '100', '143', '217', '353']
           legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1]
           set_plot, 'ps'
           device, file=tag+'_Famp.eps',/col, bits=8
           plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist)/total(hist[*,0])], xtit='!17A!dff!n', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
           for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
;   freq = ['30', '44', '70', '100', '143', '217', '353']
;   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
           legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1]
           device, /close
           set_plot, 'x'

           !p.background=255
           !p.color=0
           if (nforeg gt 5) then begin
               famp = fltarr(bins, nband)
               for i=0,nband-1 do begin
                   hist[*,i]=histogram(foreg[*,5,i],nbins=bins,locations=a)
;;       hist[*,i]=histogram(foreg[0,*,i,5],nbins=bins,locations=a)
                   famp[*,i] = a
               endfor
               window, 2 & plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist)/total(hist[*,0])], xtit='!17A!dsynch!n', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
               for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
;   freq = ['30', '44', '70', '100', '143', '217', '353']
;   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
               legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1]
               set_plot, 'ps'
               device, file=tag+'_Samp.eps',/col, bits=8
               plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist)/total(hist[*,0])], xtit='!17A!dsynch!n', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
               for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
;   freq = ['30', '44', '70', '100', '143', '217', '353']
;   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
               legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1]
               device, /close
               set_plot, 'x'
           endif

           !p.background=255
           !p.color=0
           if (nforeg gt 6) then begin
               famp = fltarr(bins, nband)
               for i=0,nband-1 do begin
                   hist[*,i]=histogram(foreg[*,6,i],nbins=bins,locations=a)
;;       hist[*,i]=histogram(foreg[0,*,i,6],nbins=bins,locations=a)
                   famp[*,i] = a
               endfor

               window, 3 & plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17A!dtd!n', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
               for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
;   freq = ['30', '44', '70', '100', '143', '217', '353']
;   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
               legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1]
               set_plot, 'ps'
               device, file=tag+'_Damp.eps',/col, bits=8
               plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17A!dtd!n', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
               for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
;   freq = ['30', '44', '70', '100', '143', '217', '353']
;   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
               legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1]
               device, /close
           endif

           set_plot, 'x'
       endif

       !p.background=255
       !p.color=0
       for i=0,nband-1 do begin
           hist[*,i]=histogram(foreg[*,0,i],nbins=bins,locations=a)
;;       hist[*,i]=histogram(foreg[0,*,i,0],nbins=bins,locations=a)
           famp[*,i] = a
       endfor
       window, 4 & plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17Monopole', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
       for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
;   freq = ['30', '44', '70', '100', '143', '217', '353']
;   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
       legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1]

       set_plot, 'ps'
       device, file=tag+'_monopole.eps',/col, bits=8
       plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17Monopole', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
       for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10   freq = ['30', '44', '70', '100', '143', '217', '353']
;   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
       legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1]
       device, /close
       set_plot, 'x'

; --- Dipole: X ---
       !p.background=255
       !p.color=0
       for i=0,nband-1 do begin
           hist[*,i]=histogram(foreg[*,1,i],nbins=bins,locations=a)
;;       hist[*,i]=histogram(foreg[0,*,i,0],nbins=bins,locations=a)
           famp[*,i] = a
       endfor
       window, 5 & plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17Dipole X', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
       for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
       legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1], chars=1.2
       
       set_plot, 'ps'
       device, file=tag+'_dipole_X.eps',/col, bits=8
       plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17Dipole X', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
       for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
       legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1], chars=1.2
       device, /close
       set_plot, 'x'

; --- Dipole: Y ---
       !p.background=255
       !p.color=0
       for i=0,nband-1 do begin
           hist[*,i]=histogram(foreg[*,2,i],nbins=bins,locations=a)
;;       hist[*,i]=histogram(foreg[0,*,i,0],nbins=bins,locations=a)
           famp[*,i] = a
       endfor
       window, 6 & plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17Dipole Y', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
       for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
       legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1], chars=1.2

       set_plot, 'ps'
       device, file=tag+'_dipole_Y.eps',/col, bits=8
       plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17Dipole Y', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
       for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
       legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1], chars=1.2
       device, /close
       set_plot, 'x'

; --- Dipole: Z ---
       !p.background=255
       !p.color=0
       for i=0,nband-1 do begin
           hist[*,i]=histogram(foreg[*,3,i],nbins=bins,locations=a)
;;       hist[*,i]=histogram(foreg[0,*,i,0],nbins=bins,locations=a)
           famp[*,i] = a
       endfor
       window, 7 & plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17Dipole Z', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
       for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
       legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1], chars=1.2
       
       set_plot, 'ps'
       device, file=tag+'_dipole_Z.eps',/col, bits=8
       plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17Dipole Z', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
       for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=col_array[i],psym=10
       legend, [freq[0:nband-1]+' GHz'], col=col_array[0:nband-1], line=line_array[0:nband-1], chars=1.2
       device, /close
       set_plot, 'x'


       !p.background=255
       !p.color=0
;stop
       if (nforeg gt 5) then begin
           window, 23 &  plot, foreg[*,4,0],foreg[*,5,0] ,psym=3, xr=[min(foreg[*,4,0]),max(foreg[*,4,0])], yr=[min(foreg[*,5,0]),max(foreg[*,5,0])],xtit='!17Free-Free', ytit='!17Synch' 
;; for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,*,4,iba],foreg[ich-1,*,5,iba],col=90+10*(ich+iba),psym=3
; oplot, foreg[0,*,4,1],foreg[0,*,6,1],col=100,psym=3
; oplot, foreg[0,*,4,3],foreg[0,*,6,3],col=255*2l^8l,psym=3
; oplot, foreg[0,*,4,2],foreg[0,*,6,2],col=250,psym=3

           if (nforeg gt 6) then begin
               window, 24 & plot, foreg[*,4,0],foreg[*,6,0] ,psym=3, xr=[min(foreg[*,4,0]),max(foreg[*,4,0])], yr=[min(foreg[*,6,0]),max(foreg[*,6,0])],xtit='!17Free-Free', ytit='!17Dust'
;; for iba=0,nband-1 do oplot, foreg[0:nsample/2-1,4,iba],foreg[0:nsample/2-1,6,iba],col=130,psym=3
;; for iba=0,nband-1 do oplot, foreg[nsample/2:*,4,iba],foreg[nsample/2:*,6,iba],col=90,psym=3
;; for ich=fst_ch, fst_ch+nchain-1 do for iba=0, nband-1 do oplot, foreg[ich-1,*,4,iba],foreg[ich-1,*,6,iba],col=90+10*(ich+iba),psym=3
;oplot, foreg[0,*,4,1],foreg[0,*,6,1],col=100,psym=3
;oplot, foreg[0,*,4,3],foreg[0,*,6,3],col=255*2l^8l,psym=3
;oplot, foreg[0,*,4,2],foreg[0,*,6,2],col=250,psym=3


               window, 25 & plot, foreg[*,5,0],foreg[*,6,0] ,psym=3, xr=[min(foreg[*,5,0]),max(foreg[*,5,0])], yr=[min(foreg[*,6,0]),max(foreg[*,6,0])],xtit='!17Synch', ytit='!17Dust' 
;; for iba=0,nband-1 do oplot, foreg[0:nsample/2-1,5,iba],foreg[0:nsample/2-1,6,iba],col=130,psym=3
;; for iba=0,nband-1 do oplot, foreg[nsample/2:*,5,iba],foreg[nsample/2:*,6,iba],col=90,psym=3
;; for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,*,5,iba],foreg[ich-1,*,6,iba],col=90+10*(ich+iba),psym=3
;oplot, foreg[0,*,5,1],foreg[0,*,6,1],col=100,psym=3
;oplot, foreg[0,*,5,3],foreg[0,*,6,3],col=255*2l^8l,psym=3
;oplot, foreg[0,*,5,2],foreg[0,*,6,2],col=250,psym=3

           endif

           set_plot, 'ps'
           device, file=tag+'_scatter_1.eps',/col, bits=8
           plot, foreg[*,4,0],foreg[*,5,0] ,psym=3, xr=[min(foreg[*,4,0]),max(foreg[*,4,0])], yr=[min(foreg[*,5,0]),max(foreg[*,5,0])]
           for iba=0,nband-1 do oplot, foreg[0:nsample/2-1,4,iba],foreg[0:nsample/2-1,5,iba],col=130,psym=3
           for iba=0,nband-1 do oplot, foreg[nsample/2:*,4,iba],foreg[nsample/2:*,5,iba],col=90,psym=3
;plot, foreg[0,*,4,0],foreg[0,*,5,0] ,psym=3, xr=[0.8,2.], yr=[-0.5,1.5]
;for ich=fst_ch, fst_ch+nchain-1 do oplot, foreg[0,*,4,ich],foreg[0,*,5,ich],col=90+10*ich,psym=3
;oplot, foreg[0,*,4,1],foreg[0,*,5,1],col=90,psym=3
;oplot, foreg[0,*,4,3],foreg[0,*,5,3],col=220,psym=3
;oplot, foreg[0,*,4,2],foreg[0,*,5,2],col=250,psym=3
           device, /close


           if (nforeg gt 6) then begin
               device, file=tag+'_scatter_2.eps',/col, bits=8
               plot, foreg[*,4,0],foreg[*,6,0] ,psym=3,xr=[min(foreg[*,4,0]),max(foreg[*,4,0])], yr=[min(foreg[*,6,0]),max(foreg[*,6,0])] 
               for iba=0,nband-1 do oplot, foreg[0:nsample/2-1,4,iba],foreg[0:nsample/2-1,6,iba],col=130,psym=3
               for iba=0,nband-1 do oplot, foreg[nsample/2:*,4,iba],foreg[nsample/2:*,6,iba],col=90,psym=3
;plot, foreg[0,*,4,0],foreg[0,*,6,0] ,psym=3, xr=[0.8,2.], yr=[-0.5,1.5]
;oplot, foreg[0,*,4,1],foreg[0,*,6,1],col=90,psym=3
;oplot, foreg[0,*,4,3],foreg[0,*,6,3],col=220,psym=3
;oplot, foreg[0,*,4,2],foreg[0,*,6,2],col=250,psym=3
               device, /close

               device, file=tag+'_scatter_3.eps',/col, bits=8
               plot, foreg[*,5,0],foreg[*,6,0] ,psym=3,xr=[min(foreg[*,5,0]),max(foreg[*,5,0])], yr=[min(foreg[*,6,0]),max(foreg[*,6,0])]
               for iba=0,nband-1 do oplot, foreg[0:nsample/2-1,5,iba],foreg[0:nsample/2-1,6,iba],col=130,psym=3
               for iba=0,nband-1 do oplot, foreg[nsample/2:*,5,iba],foreg[nsample/2:*,6,iba],col=90,psym=3
;plot, foreg[0,*,5,0],foreg[0,*,6,0] ,psym=3, xr=[-0.8,2.], yr=[-0.5,1.5]
;oplot, foreg[0,*,5,1],foreg[0,*,6,1],col=120,psym=3
;oplot, foreg[0,*,5,3],foreg[0,*,6,3],col=220,psym=3
;oplot, foreg[0,*,5,2],foreg[0,*,6,2],col=250,psym=3
               device, /close
           endif
           
           set_plot, 'x'

       endif

   endif


; ------ Cls Histrograms computation ------

   print, " Computing histograms..."

   !p.background=255
   !p.color=0

   ttcl = fltarr( nsample, 2*(ncl+1) )
   ttcl = reform(cl[*, 0, *])

   maxbins = nsample/20 ;100 ;500
;   htt = fltarr(maxbins,2*(ncl+1))
;   tta = fltarr(maxbins,2*(ncl+1))
   htt = fltarr(maxbins, nspec, 2*(ncl+1))
   tta = fltarr(maxbins, nspec, 2*(ncl+1))

   err = fltarr(nspec, 2*(ncl+1), 2)
   clstats = fltarr(2*(ncl+1),4)
   clmode = fltarr(nspec, 2 * (ncl+1) )

   low = fltarr(2*(ncl+1))
   up  = fltarr(2*(ncl+1))

   ncorr = fix(nsample / 10l) ;1000l
   corr = fltarr(nspec, 2*(ncl+1), ncorr)
   chi2_corr = fltarr(ncorr)

   if (cls_samp eq 'on') then begin

       chatty = 0b

       for ispec=0,nspec-1 do begin

           ttcl = reform(cl[*, ispec, *])
;           print, total(ttcl[*,0:ncl])
;           stop
           if (total(ttcl[*,0:ncl]) ne 0.) then begin
               print, ' '
               print, ' spectrum: ', spectra[ispec]
               
               for icl=0, 2*(ncl+1)-1 do begin

                   is = sort(ttcl[*,icl])

                   ilow = long(.16 * nsample)
                   err[ispec,icl,0] = ttcl[is[ilow],icl]

                   iup = long(.84 * nsample)
                   err[ispec,icl,1] = ttcl[is[iup],icl]

                   bins = (err[ispec,icl,1]-err[ispec,icl,0])/10

                   ilow     = long(.005 * nsample)
                   low[icl] = ttcl[is[ilow],icl]
                   iup      = long(.995 * nsample)
                   up[icl]  = ttcl[is[iup],icl]
        
                   clstats[icl,*] = moment(ttcl[*,icl])

                   temph = histogram(ttcl[*,icl], locations=a, nbins=min([nsample/25+1/nband,maxbins]) );, max = min( [clstats[icl,0]+4.*sqrt(clstats[icl,1]), up[icl]]), min=max( reform([0., clstats[icl,0]-4.*sqrt(clstats[icl,1]), low[icl] ]) ) )

                   if (0b) then begin
                       if (bins eq 0.) then temph = histogram(ttcl[*,icl], locations=a, nbins=min([nsample/25,maxbins]), max = min( [clstats[icl,0]+4.*sqrt(clstats[icl,1]), up[icl]]), min=max( reform([0., clstats[icl,0]-4.*sqrt(clstats[icl,1]), low[icl] ]) ) )

                       if (bins gt 0.) then begin
                       
                           bbb = ( min( [clstats[icl,0]+4.*sqrt(clstats[icl,1]), up[icl]])-max( reform([clstats[icl,0]-4.*sqrt(clstats[icl,1]), low[icl] ]) ) ) / bins
                           temph = histogram(ttcl[*,icl], locations=a, nbins=min([bbb,maxbins]) , max = min( [clstats[icl,0]+4.*sqrt(clstats[icl,1]), up[icl]]), min=max( reform([clstats[icl,0]-4.*sqrt(clstats[icl,1]), low[icl] ]) ) )
                       
                       endif
                   endif

                   temph = smooth(temph,5)

                   nb = n_elements(temph)


;       htt[0:nb-1, icl] = temph ;htt[0:nb[1]-1, icl] = temph
;       tta[0:nb-1, icl] = a ;tta[0:nb[1]-1, icl] = a
           
                   htt[0:nb-1, ispec, icl] = temph 
                   tta[0:nb-1, ispec, icl] = a 

; error bars with cumulants, but if not two-tail distribution doesn't work
;;        frac = htt[0:nb[1]-1,icl] / total(htt[0:nb[1],icl])
;;         cum = frac[*] * 0.
;;         for ia = 0,nb[1]-1 do cum[ia]=total(frac[0:ia])
;;         print, frac
;;         print, cum
;; 
;;         ilow = min(abs(cum-.16),i)
;;         err[icl, 0] = tta[i,icl]
;;         if (i eq 0) then err[icl,0] = 0.
;;         iup  = min(abs(cum-.84),i)
;;         err[icl, 1] = tta[i,icl]

                   frac = htt[0:nb-1,ispec, icl] / max(htt[0:nb-1,ispec,icl]) ;frac = htt[0:nb[1]-1,icl] / max(htt[0:nb[1],icl])
                   
                   imax = max(frac,i)
                   if (chatty) then print, imax, i
                   clmode[ispec,icl] = tta[i,ispec,icl]

                   err[ispec, icl,0] = 0.
                   if (i ne 0) then begin
                       l32 = min( abs(frac[0:i-1l] -.32), i32 )
                       err[ispec, icl, 0] = tta[i32, ispec, icl]
                       if (i32 eq 0) then err[ispec, icl,0] = 0. 
                       if (chatty) then print, err[ispec, icl, 0], i32
                   endif

                   l32 = min( abs(frac[i+1l:*] -.32), i32 )
                   err[ispec, icl, 1] = tta[i32+i, ispec, icl]
                   if (chatty) then print, err[ispec, icl,1], i32+i
;stop
                   if (chatty) then print, icl, min(tta[*,ispec,icl]), err[ispec,icl,0], clmode[ispec,icl], err[ispec,icl,1], max(tta[*,ispec,icl])
       
                   if (chatty) then print, " "

;stop
;       is = sort(tta[*,icl] * htt[*,icl]/nsample)
;       ilow = long(.16 * bins)
;       err[icl,0] = tta[is[ilow],icl]
;       iup = long(.84 * bins)
 ;      err[icl,1] = tta[is[iup],icl]
;
    
                   if (cls_samp eq 'on' and do_corr) then for b=0l,ncorr-1 do corr[ispec,icl,b] = correlate(cl[0:nsample-ncorr,ispec,icl],cl[0+b:nsample-ncorr+b,ispec,icl] )

               endfor ; loop over cls
           endif ; check if cls for that field are vanishing
       endfor ; loop over nspec 
   endif ; cls on

print, 'Correlation not printed...'
;;    for b=0l, ncorr-1 do chi2_corr[b] = correlate(chi2[0:nsample-ncorr],chi2[0+b:nsample-ncorr+b] )
;;    window, 8 & plot, chi2_corr, ytit='!7v!u2 !n!17Corr', xtit='!17N!dsample!n',yr=[-0.2,1.2], xr=[-10,300], xs=1, ys=1
;;    write_png,tag+'_chi2_correlation.png',tvrd(/true)

   gb = gaussbeam(beams,ncl)

   !p.multi = 0
   !p.charsize=1.5

   if (do_wmap) then readcol, '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/wmap/wmap_tt_spectrum_7yr_v4p1.txt', wmapl, wmapcl, wmaperr
   if not (do_wmap) then begin
       print, " anafast for ffp3 input map..."
       if not (do_pol) then ianafast, '../../map_cmb.fits', ctp3cls, nlmax=ncl, /silent, regression=2
       if (do_pol) then ianafast, '../../map_cmb.fits', ctp3cls, nlmax=ncl, simul_type=2, /silent, regression=2
;       large_gb=gaussbeam(600,ncl)
       large_wl=healpixwindow(2048)
       ctp3cls[*,*] = ctp3cls[*,*] * 1.e12 ;/ large_wl^2

; --- Healpix order TT EE BB TE TB EB
; --- Commander order TT TE TB EE EB BB
      if (do_pol) then begin 
       comm_cls = ctp3cls*0.
       comm_cls[*,0] = ctp3cls[*,0]
       comm_cls[*,1] = ctp3cls[*,3]
       comm_cls[*,2] = ctp3cls[*,4]
       comm_cls[*,3] = ctp3cls[*,1]
       comm_cls[*,4] = ctp3cls[*,5]
       comm_cls[*,5] = ctp3cls[*,2]

       ctp3cls = comm_cls
endif
   endif

; ------ plotting histrograms (if computed) ------
   if (cls_samp eq 'on') then begin
       for ispec=0,nspec-1 do begin

           if (total(tta[*,ispec,*]) ne 0.) then begin

               print, spectra[ispec]
               window, 16+2*ispec, tit='Spectrum '+spectra[ispec] & 
;               !p.multi = [0,3,3] 
               !p.multi = [0,4,4]
;               lplot = [2,3,4,6,8,10,12,14,16,18,20,24,28,32,36,40] ;,62]

;               lplot = lonarr(9) 
               lplot = lonarr(16)
;               for il=0,8 do lplot[il] = 4+il*(ncl-2)/9
;               for il=0,15 do lplot[il] = 4+il*(ncl-2)/16l
;;                for il=0,15 do lplot[il] = 4+il*(90-2)/16l
; --- To plot sigma ell of the sky
               max_ncl = 250
               for il=0,15 do lplot[il] = 4+il*(max_ncl-4)/16l
   
               print, lplot

	       slplot = string(lplot, format='(1i3.3)')
               !p.background=255
               !p.color=0

               for il=0,n_elements(lplot)-1 do begin
if (1b) then begin
x = tta[*, ispec, lplot[il]]
y = htt[*, ispec, lplot[il]] / total(htt[*, ispec, lplot[il]])
ym = max(y,im)
thres = 0.01
lx = where(y[0:im] lt ym*thres)
if (lx[0] eq -1) then lx = 0
if (lx[0] ne -1) then lx = lx[n_elements(lx)-1]
ux = where(y[im:*] lt ym*thres)
ux = ux[0]
print, x[lx[0]],x[im+ux[0]]
                   plot, tta[*, ispec, lplot[il]], htt[*, ispec, lplot[il]] / total(htt[*, ispec, lplot[il]]), ytit='!17P', xtit='!17'+strtrim(string(lplot[il]),2), psym=10, xs=1, xr=[x[lx[0]],x[im+ux[0]]]
;;                    oplot, tta[*, ispec, ncl+lplot[il]], htt[*, ispec, ncl+lplot[il]] / total(htt[*, ispec, ncl+lplot[il]])
                   oplot, tta[*, ispec, lplot[il]], smooth(htt[*, ispec, lplot[il]] / total(htt[*, ispec, lplot[il]]),5), col=70, thick=1.5;, ytit='!17P', xtit='!17C!d'+strtrim(string(lplot[il]),2), psym=10;, chars=2.5;, xr=[3*err[ispec,lplot[il],0], 3*err[ispec,lplot[il],1]] ;xr=[low[lplot[il]], up[lplot[il]] ]

endif
if (0b) then begin
   readcol,'cls_2pl_single_TT_l0'+slplot[il]+'_0'+slplot[il]+'.dat',x,y
                   plot, x, y, ytit='!17P', xtit='!17C!d'+slplot[il], psym=10, chars=2.5
endif

                   if (do_wmap) then oplot, reform(lplot[*]*0.+wmapcl[lplot[il]-2]), reform(lplot[*]*0.+findgen(n_elements(lplot))/n_elements(lplot)),col=245, line=2

                   if not (do_wmap) then oplot, reform(lplot[*]*0.+ctp3cls[lplot[il],ispec] *lplot[il]*(lplot[il]+1)/2./!pi), reform(lplot[*]*0.+findgen(n_elements(lplot))/n_elements(lplot)),col=245, line=3

               endfor ; loop over l_plot

;;                write_png,'hke_'+tag+'_ClsHist.png',tvrd(/true)
               write_png,tag+'_ClsHist.png',tvrd(/true)

               if (do_corr) then begin
                   window, 9+2*ispec, tit='Correlation '+spectra[ispec] & 
                   for il=0,n_elements(lplot)-1 do plot, corr[ispec,lplot[il],*], ytit='!17C!dl !nCorr', xtit='!17N!dsample!n', xr=[-10,300],xs=1,yr=[-0.2,1.2],ys=1
                   write_png,tag+'_cls_correlation.png',tvrd(/true)
               end

           endif ; check if histograms are vanishing

       endfor ; loop over spec

   endif

;;    set_plot, 'ps'
;;    device, file=tag+'_cls.eps', /col, bits=8, /encapsulated    
;;    for il=0,n_elements(lplot)-1 do begin
;;        plot, tta[*, lplot[il]],htt[*, lplot[il]]/total(htt[*, lplot[il]]), ytit='!17P', xtit='!17C!d'+strtrim(string(lplot[il]),2), psym=10
;;         oplot, reform(lplot[*]*0.+wmapcl[lplot[il]-2]), reform(lplot[*]*0.+findgen(n_elements(lplot))/n_elements(lplot)),col=245, line=2
;;    endfor
;; device, /close
;; set_plot, 'x'


   save, filename=tag+'_histCls.res.sav', htt, tta, corr, chi2_corr

   if not (do_pol) then !p.multi=0
   if (do_pol) then !p.multi=[0,2,3] 

   if (cls_samp eq 'on') then begin

   gb = gaussbeam(beams,ncl)
   l = lindgen(ncl+1)
   ll = l*(l+1)/2./!pi

   window, 30, tit='CMD Spectra'
   for ispec=0, nspec-1 do begin

       if (do_wmap) then begin
           fits2cl, lcdm_cl, '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/camb_85760952_scalcls.fits'
           plot, wmapl[0:ncl], wmapcl[0:ncl,ispec], xtit='!17l', ytit='!17l(l+1)C!dl!n/2!7p!17 !7l!17K!u2',yr=[-100.,6000], background=255, color=0, ys=1,xr=[2,lplot[n_elements(lplot)-1]],xs=1, line=2
           oplot, wmapl[0:ncl]*(1.03), wmapcl[0:ncl,ispec], col=245, thick=2;, line=2
           errplot, wmapl[0:ncl]*(1.03), wmapcl[0:ncl,ispec]-wmaperr[0:ncl],wmapcl[0:ncl,ispec]+wmaperr[0:ncl], col=245
;;            oplot, l, lcdm_cl[*,0]*ll*1.e12, line=4, thick=2
           legend, ['!17White Noise C!dl!n', '7yr WMAP', '!17C!dl!n mode'], col=[0,245,70], line=[3,2,0], chars=1
       endif
       if not (do_wmap) then begin
           plot, l[2:*], ctp3cls[2:*,ispec]*ll[2:*], xtit='!17l', ytit='!17l(l+1)C!dl!n/2!7p!17 !7l!17K!u2', background=255, color=0, ys=2;,yr=[min(clmode[ispec,2:*]),max(clmode[ispec,2:*])]
           oplot, l[2:*], ctp3cls[2:*,ispec]*ll[2:*], col=245
           legend, ['!17White Noise C!dl!n', '!17FS C!dl', '!17C!dl!n mode'], col=[0,245,70], line=[3,0,4], chars=1
       endif
;;        oplot, l, ncls/gb^2*ll, line=3
       oplot, l, l[*]*0.+gsig^2*4.*!pi/npix/gb^2*ll, line=3

       oplot, l[2:*], clmode[ispec,2:*], col=70, thick=2
       oplot, l[2:*], clmode[ispec,2:*], col=70, psym=1
       errplot, l[2:*], (err[ispec,2:*,0]), (err[ispec,2:*,1]), col=70


       if not (do_pol) then begin
           write_png,tag+'_Cls.png',tvrd(/true)
       if not (do_wmap) then begin
           window, 29
           plot, l[2:*], ctp3cls[2:*,ispec]*ll[2:*]-clmode[ispec,2:*], xtit='!17l', ytit='!17Input - Commander l(l+1)C!dl!n/2!7p!17 !7l!17K!u2', background=255, color=0, ys=2, psym=4, xr=[2,lplot[n_elements(lplot)-1]]
           oplot, l[2:*], l[2:*]*0., line=2
           errplot, l[2:*], ctp3cls[2:*,ispec]*ll[2:*]-clmode[ispec,2:*]-(clmode[ispec,2:*]-err[ispec,2:*,0]), ctp3cls[2:*,ispec]*ll[2:*]-clmode[ispec,2:*]+(err[ispec,2:*,1]-clmode[ispec,2:*])
           write_png,tag+'_Cls_diff.png', tvrd(/true)
       endif
       if (do_wmap) then begin
           if (1b) then begin
               window, 29
               plot, wmapl, wmapcl[0:ncl,ispec] - clmode[ispec,2:ncl], xr=[2,lplot[n_elements(lplot)-1]], xs=1, yr=[-1500,1500], ys=1, xtit='!17l',psym=4
               oplot, wmapl, wmapcl[0:ncl,ispec] - clmode[ispec,2:ncl], col=245, thick=2, psym=4
               oplot, wmapl, wmapl*0.+0.
               errplot, l[2:*], (wmapcl[0:ncl,ispec] - clmode[ispec,2:ncl])-(clmode[ispec,2:ncl] - err[ispec,2:*,0] ), (wmapcl[0:ncl,ispec] - clmode[ispec,2:ncl])+(err[ispec,2:*,1]-clmode[ispec,2:ncl])
               legend,['WMAP7 - CMD'],col=[245],line=0, pos=[10,1400], chars=1.2, psym=4, thick=2
               
               write_png,tag+'_Cls_diff.png',tvrd(/true)
               if (0b) then begin
                   plot, l[2:*], lcdm_cl[2:*]*ll[2:*]*1.e12*gb^2*fsky, ys=2, chars=1.5, xtit='!17l',ytit='l(l+1)C!dl!n/2!7p!17'
                   oplot, l[2:*], ilc_cl[2:*]*ll[2:*], col=70
                   oplot, l[2:*], cmb_cl[2:*]*ll[2:*], col=245
                   legend, ['WMAP7 ILC', 'Mean CMB'], col=[70, 245], line=[0,0], pos=[160,950]
                   write_png,tag+'_ilc-cmb_cls_comp.png',tvrd(/true)
               endif
           endif
       endif

       endif
;       for icl=0,nsample-1 do oplot, cl[icl,ispec,0:ncl],psym=3

   endfor ; loop over nspec

    if (do_pol) then write_png,tag+'_convCls_2.png',tvrd(/true)

    data = fltarr(ncl+1,nspec*3)
    for icl=0,ncl do begin
        for ispec=0,nspec-1 do data[icl,3*ispec:3*(ispec+1)-1] = reform( [clmode[ispec,icl], err[ispec,icl,0], err[ispec,icl,1]] )
    endfor

    openw, 1,tag+'_cls.dat'
    for icl=0,ncl do printf, 1, icl, data[icl,*], format='(1i20,'+string(nspec*3)+'f20.3)'
    close, 1

    endif

   if (cls_samp eq 'off') then begin

   if (do_wmap) then readcol, '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/wmap/wmap_tt_spectrum_7yr_v4p1.txt', wmapl, wmapcl


   gb = gaussbeam(beams,ncl)
   l = lindgen(ncl+1)
   ll = l*(l+1)/2./!pi

   ianafast, tag+'_recCMB_ns'+s_ns+'_RING.fits', cls, nlmax=ncl, maskfile=maskfile, /silent, regression=2

   if (do_wmap) then begin
       window, 30 & plot, wmapl[0:*], wmapcl, xtit='!17l', ytit='!17l(l+1)C!dl!n/2!7p!17 !7l!17K!u2!n',yr=[-500.,7000], background=255, color=0, xr=[2,ncl]
       oplot, wmapl, wmapcl, col=245
       legend, ['!17Reg.W.Noise C!dl!n', '7yr WMAP', '!17Mean map C!dl!n'], col=[0,245,70], line=[3,0,4]
   endif
   if not (do_wmap) then begin
       window, 30 & plot, l[2:*], ctp3cls[2:*]*ll[2:*], xtit='!17l', ytit='!17l(l+1)C!dl!n/2!7p!17 !7l!17K!u2!n',yr=[-500.,6000], background=255, color=0, xr=[2,ncl]
       oplot, l[2:*], ctp3cls[2:*]*ll[2:*], col=245
       legend, ['!17Reg.W.Noise C!dl!n', '!17FS C!dl', '!17Mean map C!dl!n'], col=[0,245,70], line=[3,0,4]
   endif
;;    oplot, l, ncls/gb^2*ll, line=3
   oplot, l, l[*]*0.+gsig^2*4.*!pi/npix/gb^2*ll

   oplot, l[2:*], cls[2:*]*ll[2:*]/gb[2:*]^2/fsky, col=70, line=4
    write_png,tag+'_convCls_2.png',tvrd(/true)
   endif 


   if (do_noise) then begin
       oplot, l, avenls*ll/gb^2/fsky
       if (cls_samp eq 'off') then oplot, (cls-avenls)*ll/fsky/gb^2
   endif

print, " --- End of Code ---"

stop

end
