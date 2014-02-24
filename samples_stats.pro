
pro samples_stats, root=root, fst_ch=fst_ch, nforeg=nforeg, tag=tag, do_read=do_read, nband=nband, cls_samp=cls_samp, nsample=nsample, do_pol=do_pol, maskfile=maskfile, fst_samp=fst_samp, icmb=icmb, do_corr=do_corr, corr_ipix=corr_ipix


   close, /all
   loadct, 39
   set_plot, 'x'
   mollview, randomn(-1,12), window=-1
   !p.multi = 0
   !p.background=255
   !p.color=0

   if not (keyword_set(beams)) then beams = 40.
   s_beams = strtrim(string(long(beams)),2)

    if not (keyword_set(root)) then begin
        print, 'root not defined'
        stop
    endif
    if not (keyword_set(tag)) then begin
        print, 'tag not defined'
        stop
    endif        

;   if not (keyword_set(nchain)) then nchain = 1 ;10l ;; 4l
   if not (keyword_set(fst_ch)) then fst_ch = 1
   if not (keyword_set(do_read)) then do_read = 'read' ;10l ;; 4l
   if not (keyword_set(nforeg)) then nforeg = 0

   if not (keyword_set(cls_samp))  then cls_samp = 'off'

   if not (keyword_set(do_pol)) then do_pol = 'no'
   if not (keyword_set(fst_samp)) then fst_samp = 1l

   if not (keyword_set(icmb)) then icmb = 3

   if not (keyword_set(do_corr)) then do_corr = 'corr'

   print, do_pol

   if (do_pol eq 'pol') then begin
       do_pol = 1b
       nspec = 6
   endif else begin  
       do_pol = 0b
       nspec = 1
   endelse

   print, ' Restoring '+tag+'.res.sav to load mean and rms'
   restore, tag+'.res.sav'
   nchains = n_elements(nsampch)
   if (cls_samp eq 'off') then begin
       cmbvar = amp2[*,icmb]-amp[*,icmb]^2
       avecmb = amp[*,icmb]
   endif
   if (cls_samp eq 'on') then begin
       cmbvar = cmb2-cmb^2
       avecmb = cmb
   endif

   ampvar = amp2-amp^2
   indvar = ind2-ind^2
   
   ns = nside
   npix = 12l*ns^2
   s_ns = strtrim(string(ns),2)

   if (do_corr eq 'corr') then begin
       if not (keyword_set(corr_ipix)) then corr_ipix = [10, npix/2, npix-1]
       do_corr = 1b
       print, corr_ipix
       ncpix = n_elements(corr_ipix)
       corrmap = fltarr(npix, ncpix)

       cmbdist = fltarr(nsample, ncpix)
       ampdist = fltarr(nsample, namp, ncpix)
       inddist = fltarr(nsample, nind, ncpix)
   endif

;   if (do_read eq 'read') then begin
;   if (fg_pix eq 't') then begin
;       fg_pix = 1b
;       spawn, 'ls -l '+root+'fg_amp_map_no*_c0001_k00001.fits | wc -l', namp
;       namp = reform(long(namp))
;       print, ' amplitude maps # = ', namp
;       spawn, 'ls -l '+root+'fg_ind_map_no*_c0001_k00001.fits | wc -l', nind
;       nind = reform(long(nind))
;       print, ' index maps # = ', nind
;       if not (do_pol) then amp = fltarr(npix, namp[0])
;       if (do_pol) then amp = fltarr(npix, 3, namp[0])

       amp3 = amp * 0.   
       amp4 = amp * 0.
       ind3 = ind * 0.
       ind4 = ind * 0.

;       if (nind gt 0) then begin
;          ind = fltarr(npix, nind[0])
;          ind2 = ind
;       endif
;   endif else begin
;       fg_pix = 0b
;   endelse
;   endif

;   if (pl eq 't') then begin
;       pl = 1b
;   endif else begin
;       pl = 0b
;   endelse

;   if (nforeg ge 0) then nforeg = nforeg + 4

;   print, ' nforeg : ', nforeg
   print, ' nchain : ', nchain
   print, ' fst_ch : ', fst_ch
   print, ' do_read: ', do_read
   print, ' cls_samp: ', cls_samp
;   print, ' ncl:      ', ncl
;   print, ' gsig:     ', gsig
;   print, ' beams:    ', beams
;   print, ' do_residuals', do_residuals
;   print, ' do_wmap     ', do_wmap
   print, ' do_corr     ', do_corr
   print, ' do_pol      ', do_pol
;   print, ' nband:    ', nband, ' <== !!! ==> '
   print, ' icmb:      ', icmb

   freq        = [30, 44, 70, 100, 143, 217, 353, 33, 23, 41, 61, 94]
   freq = string(freq, format='(i3.3)')

   nfreq       = n_elements(freq)
   col_array   = 40 * lindgen(nfreq); *250/nfreq ;; [200, 205, 210, 215, 220, 225, 230]
   line_array  = col_array * 0
   spectra     = ['TT', 'EE', 'BB', 'TE', 'TB', 'EB']

;   read_fits_map, maskfile, mask

   gpix = where(mask[*,0] gt 0.)
   ngpix = n_elements(gpix)
   bpix = where(mask[*,0] eq 0.)
   bad_value = -1.6375e30
   fsky = float(ngpix) / n_elements(mask[*,0])
   print, " - fsky = ", fsky

   if (do_read eq 'read') then begin
;       nsampch = lonarr(nchain)
;       if (nsample lt 0l) then begin
;           nsample = 0l
;           iich = 0l
;           for ich=fst_ch,fst_ch+nchain-1 do begin
;               sch = 'c'+string(ich,format='(1i4.4)')
;               spawn, 'ls -ltrh '+root+'chisq_'+sch+'_k* | wc -l', nsamp
;               nsampch[iich] = long(nsamp)
;               nsample = nsample + nsampch[iich]
;               iich = iich + 1l
;           endfor
;       endif else if (nsample gt 0l) then begin 
;           nsampch[*] = nsample / nchain
;       endif
   
   print, ' nsample:      ', nsample
   print, ' nsampch:      ', nsampch

   if (min(nsampch) lt fst_samp) then begin
       print, " Not enough samples..."
       stop
   endif

   print, nsample
   nsample = nsample - nchain*(fst_samp-1)
   print, fst_samp, nchain
   print, nsample

   cl=fltarr(nsample,nspec,2 * (ncl+1))

   l = findgen(ncl+1)

   if (nforeg gt 0) then foreg = fltarr(nsample, nforeg, nband)

   chi2 = dblarr(nsample)

   if (cls_samp eq 'on') then begin
      if not (do_pol) then cmb = dblarr(npix)
      if (do_pol) then cmb = dblarr(npix,3) 
      cmb2 = cmb
      cmb3 = cmb
      cmb4 = cmb
      glb_chi2 = cmb * 0.
   endif else begin
      if not (do_pol) then cmb = dblarr(npix)
      if (do_pol) then cmb = dblarr(npix,3)
      cmb2 = cmb
      cmb3 = cmb
      cmb4 = cmb
      glb_chi2 = reform( amp[*,0,*] ) * 0.
   endelse

;   if (do_residuals) then residuals = fltarr(npix,nband,nsample)

   iich = -1l
   tot_isamp = -1l

   warn = 1b

;   if (do_noise) then begin
;       if (cls_samp eq 'off') then read_fits_map, root+'fg_amp_map_no'+string(icmb+1,format='(i2.2)')+'_c'+string(fst_ch,format='(i4.4)')+'_k'+string(fst_samp,format='(i5.5)')+'.fits', ref_samp
;       if (cls_samp eq 'on') then read_fits_map, root+'s_i_'+sch+'_'+ssam+'.fits', ref_samp 
;       nl = fltarr(nsample,ncl+1)
;   endif

   for ich=fst_ch, fst_ch+nchain-1 do begin
      sch = 'c'+string(ich,format='(1i4.4)')
      iich = iich + 1

      for isam=fst_samp, nsampch[iich] do begin
         cl_tmp = fltarr(nspec+1, 2*(ncl+1))
         if (nforeg gt 0) then fg = fltarr(nforeg, nband)
         tot_isamp = tot_isamp+1
         ts = tot_isamp / 100
         if ( (ts*100) eq tot_isamp ) then begin
             print, tot_isamp
             save, filename=tag+'.samples_stats.sav', cmb, cmb2, cmb3, cmb4, cmbvar, avecmb, corrmap, nsample, nsampch, nside, fsky, tot_isamp, amp, amp2, ind, ind2, namp, nind, mask, amp3, amp4, ind3, ind4, ampdist, inddist, cmbdist
         endif

         ssam = 'k'+string(isam,format='(1i5.5)')

;         read_fits_map, root+'chisq_'+sch+'_'+ssam+'.fits', chisq
;         glb_chi2 = glb_chi2 + chisq

         if (cls_samp eq 'on') then begin
;             openr,1,root+'cl_'+sch+'_'+ssam+'.dat'
;             readf,1,cl_tmp
;             close, 1
;             
;             cl[tot_isamp,*,*] = reform(cl_tmp[1:*,*])
;
             read_fits_map, root+'s_i_'+sch+'_'+ssam+'.fits', cmbmap
             cmb3 = cmb3 + (cmbmap-avecmb)^3
             cmb4 = cmb4 + (cmbmap-avecmb)^4
         endif

;         chi2[tot_isamp] = total(chisq) / ngpix / nband

         if (cls_samp eq 'off') then begin
             if (warn) then begin
                 print, " WARNING: cls_off polarization not implemented yet!!!"
                 warn = 0b
             endif
             
             read_fits_map, root+'fg_amp_map_no'+string(icmb+1, format='(i2.2)')+'_'+sch+'_'+ssam+'.fits', cmbmap
             cmb3 = cmb3 + (cmbmap-avecmb)^3
             cmb4 = cmb4 + (cmbmap-avecmb)^4
         endif

         if (do_corr) then begin
             temp = cmbmap - avecmb
             for icpix=0,ncpix-1 do begin
                 corrmap[*,icpix] = corrmap[*,icpix] + temp*temp[corr_ipix[icpix]]
                 cmbdist[tot_isamp,icpix] = cmbmap[corr_ipix[icpix]]
             endfor
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

;         if (fg_pix) then begin
             for iamp = 1, namp[0] do begin
                read_fits_map, root+'fg_amp_map_no'+string(iamp,format='(i2.2)')+'_'+sch+'_'+ssam+'.fits', map
                amp3[*,iamp-1] = amp3[*,iamp-1] + (map-amp[*,iamp-1])^3
                amp4[*,iamp-1] = amp4[*,iamp-1] + (map-amp[*,iamp-1])^4
                for icpix=0,ncpix-1 do ampdist[tot_isamp,iamp-1,icpix] = map[corr_ipix[icpix]]
             endfor
             
             if (nind gt 0) then begin
                for iind =1, nind[0] do begin
                   read_fits_map, root+'fg_ind_map_no'+string(iind,format='(i2.2)')+'_'+sch+'_'+ssam+'.fits', map
                   ind3[*,iind-1] = ind3[*,iind-1] + (map-ind[*,iind-1])^3
                   ind4[*,iind-1] = ind4[*,iind-1] + (map-ind[*,iind-1])^4
                   for icpix=0,ncpix-1 do inddist[tot_isamp,iind-1,icpix] = map[corr_ipix[icpix]]
                endfor
             endif
;         endif

;         if (do_residuals) then begin
;             for iband=0,nband-1 do begin
;
;                 read_fits_map, '../map_freefree_0'+freq[iband]+'.00GHz.smth.ns'+s_ns+'.fits', ff
;                 read_fits_map, '../map_synchrotron_0'+freq[iband]+'.00GHz.smth.ns'+s_ns+'.fits', ss
;                 read_fits_map, '../map_thermal_dust_0'+freq[iband]+'.00GHz.smth.ns'+s_ns+'.fits', dd
;                 read_fits_map, '../noDip_dx4_'+freq[iband]+'GHz_I_smth600_uK9_ns'+s_ns+'.fits', imap
;                 dipole = make_dipole(reform(fg[1:3,iband]),ns)
;                 residuals[*,iband,tot_isamp] = (imap - (map + fg[0,iband] + dipole + fg[4,iband]*ff + fg[5,iband]*ss + fg[6,iband]*dd))
;              
;             endfor
;         endif


     endfor
 endfor

 print, nsample

 cmb3 = cmb3 / nsample
 cmb4 = cmb4 / nsample

 amp3 = amp3 / nsample
 amp4 = amp4 / nsample

 ind3 = ind3 / nsample
 ind4 = ind4 / nsample

 corrmap = corrmap / nsample

; cmb3n = cmb3 / sqrt(cmbvar)^3
; cmb4n = cmb4 / cmbvar^2 - 3.
;
; corrmapn = corrmap * 0.
; for icpix=0,ncpix-1 do corrmapn[*,icpix] = corrmap[icpix] / (sqrt(cmbvar*cmbvar[corr_ipix[icpix]]))

; if (cls_samp eq 'on') then begin
;    cmb = cmb / (nsample)
;    cmb2 = cmb2 / (nsample)
; endif
;   glb_chi2 = glb_chi2 / nsample
;
;   if (fg_pix) then begin
;       amp = amp / nsample
;       amp2 = amp2 / nsample
;
;       ampstd=sqrt(amp2-amp^2)
;
;       ind = ind / nsample
;       ind2 = ind2 / nsample
;
;       indstd=sqrt(ind2-ind^2)
;   endif


   save, filename=tag+'.samples_stats.sav', cmb, cmb2, cmb3, cmb4, cmbvar, avecmb, corrmap, nsample, nsampch, nside, fsky, amp, amp2, ind, ind2, namp, nind, mask, amp3, amp4, ind3, ind4, ampdist, inddist, cmbdist

   endif

   if (do_read ne 'read') then begin
        print, "restoring..."
;;         freq = ['030', '044', '070', '100', '143', '217', '353']
        freq        = ['030', '044', '070', '100', '143', '217', '353', '033', '023', '041', '061', '094']
        nfreq = n_elements(freq)
	restore, filename=tag+'.samples_stats.sav'
        npix = 12l*nside^2
        map = fltarr(npix)
   endif

 cmb3n = cmb3 / sqrt(cmbvar)^3
 cmb4n = cmb4 / cmbvar^2 - 3.

 amp3n = amp3 / sqrt(ampvar)^3
 amp4n = amp4 / ampvar^2 - 3.

 ind3n = ind3 / sqrt(indvar)^3
 ind4n = ind4 / indvar^2 - 3.

 corrmapn = corrmap * 0.

 for icpix=0,ncpix-1 do corrmapn[*,icpix] = corrmap[*,icpix] / sqrt(cmbvar*cmbvar[corr_ipix[icpix]])

   ccc = cmb3n
   if (bpix[0] gt 0) then ccc[bpix] = bad_value	
   mollview, ccc, chars=1.5, min=-.5, max=.5, tit='CMB map, skewness' 
;   mollview, cmb3n, chars=1.5, min=-1, max=1, tit='CMB map, skewness', png=tag+'_CMBskewness.png', window=-1

   ccc = cmb4n
   if (bpix[0] gt 0) then ccc[bpix] = bad_value
   mollview, ccc, chars=1.5, min=-.5, max=.5, tit='CMB map kurtosis'
;   mollview, cmb4n, chars=1.5, min=-1, max=1, tit='CMB map kurtosis', png=tag+'_CMBkurtosis.png', window=-1

   sAmp = ['Dust','Sync','CO']
   for iamp=1,namp[0] do begin
       ccc = amp3n[*,iamp-1]
       if (bpix[0] gt 0) then ccc[bpix] = bad_value	
       mollview, ccc, chars=1.5, min=-.5, max=.5, tit=sAmp[iamp-1]+' Amp skewness'

       ccc = amp4n[*,iamp-1]
       if (bpix[0] gt 0) then ccc[bpix] = bad_value
       mollview, ccc, chars=1.5, min=-.5, max=.5, tit=sAmp[iamp-1]+' Amp kurtosis'
   endfor
;                        
   sInd = ['Emissivity','Td','Beta']                                                    
   for iind=1,nind[0] do begin
       ccc = ind3n[*,iind-1]
       if (bpix[0] gt 0) then ccc[bpix] = bad_value	
       mollview, ccc, chars=1.5, min=-1, max=1, tit=sInd[iind-1]+' Ind skewness'

       ccc = ind4n[*,iind-1]
       if (bpix[0] gt 0) then ccc[bpix] = bad_value
       mollview, ccc, chars=1.5, min=-1, max=1, tit=sInd[iind-1]+' Ind kurtosis'
   endfor

   if (do_corr) then begin
       for icpix=0,ncpix-1 do mollview, corrmapn[*,icpix], chars=1.5, min=-.5, max=.5, tit='!8Correlation map: No '+string(corr_ipix[icpix],2);, /asinh
;       for icpix=0,ncpix-1 do mollview, corrmapn[*,icpix], chars=1.5, min=-1, max=1, tit='Correlation map: '+corr_ipix[icpix], png=tag+'_corrmap_'+strtrim(string(corr_ipix[icpix]),2)+'.png'
   endif

test = randomn(-1,nsample)
print, moment(test)

nb = 8
cd = fltarr(nb,ncpix)
acd = fltarr(nb, ncpix)
window,1 & plot, /nodata, [-4,4], [0,0.5], chars=2, xtit='!7s!8', ytit='!8P'
for icpix=0,ncpix-1 do begin
    h = (cmbdist[*,icpix]-avecmb[corr_ipix[icpix]]) / sqrt(cmbvar[corr_ipix[icpix]])
    cd[*,icpix] = histogram(h, nbins=nb, locations=a)
    acd[*,icpix] = a
    oplot, acd[*,icpix], cd[*,icpix]/nsample, col=col_array[icpix], thick=2, psym=10
endfor
legend,'Pix '+string(corr_ipix), col=col_array[0:ncpix-1], line=lindgen(ncpix)*0., chars=2

print, " --- End of Code ---"

stop

end
