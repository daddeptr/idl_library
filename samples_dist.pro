
pro samples_dist, root=root, fst_ch=fst_ch, nforeg=nforeg, tag=tag, do_read=do_read, nband=nband, cls_samp=cls_samp, nsample=nsample, do_pol=do_pol, maskfile=maskfile, fst_samp=fst_samp, icmb=icmb, do_corr=do_corr, corr_ipix=corr_ipix, indx=indx, do_amp=do_amp, do_ind=do_ind


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

   print, ' nchain : ', nchain
   print, ' fst_ch : ', fst_ch
   print, ' do_read: ', do_read
   print, ' cls_samp: ', cls_samp
   print, ' do_corr     ', do_corr
   print, ' do_pol      ', do_pol
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

       dist = fltarr(nsample, npix)

       iich = -1l
       tot_isamp = -1l
   
       warn = 1b

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

               if (do_amp) then begin
;             for iamp = 1, namp[0] do begin
                   iamp = indx
                   read_fits_map, root+'fg_amp_map_no'+string(iamp,format='(i2.2)')+'_'+sch+'_'+ssam+'.fits', map
                   dist[tot_isamp, *] = map[*,0]
;             endfor
               endif
             
               if (do_ind) then begin
;                for iind =1, nind[0] do begin
                   iind = indx
                   read_fits_map, root+'fg_ind_map_no'+string(iind,format='(i2.2)')+'_'+sch+'_'+ssam+'.fits', map
                   dist[tot_isamp, *] = map[*,0]
               endif
;         endif

     endfor
 endfor

 print, nsample
   save, filename=tag+'.samples_stats.sav', cmb, cmb2, cmbvar, avecmb, corrmap, nsample, nsampch, nside, fsky, amp, amp2, ind, ind2, namp, nind, mask, ampdist, inddist, cmbdist, dist
stop


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
