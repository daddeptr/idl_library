pro residuals, root=root, nchain=nchain, fst_ch=fst_ch, nforeg=nforeg, tag=tag, do_read=do_read, nband=nband

   close, /all
   loadct, 39
   set_plot, 'x'
   !p.multi = 0

;;    root = '0002/'
;;    root='/project/projectdirs/planck/user/graca/COMMANDER/runs/ctp3/ctp3_commander_inputs/chains/'
   if not (keyword_set(root)) then root = '7b_DipSamp/'
   if not (keyword_set(tag)) then tag = '7b_DipSamp'
;   if not (keyword_set(nsample)) then nsample = 3294l ;1000l ;; 416l
   if not (keyword_set(nchain)) then nchain = 1 ;10l ;; 4l
   if not (keyword_set(fst_ch)) then fst_ch = 4
   if not (keyword_set(do_read)) then do_read = 'read' ;10l ;; 4l
   if not (keyword_set(nforeg)) then nforeg = 3
   if not (keyword_set(nband))  then nband = 7

   ncl = 40l ;32l
   nspec = 1
;   nband = 7
   nforeg = nforeg + 4

   print, ' nforeg : ', nforeg
   print, ' nchain : ', nchain
   print, ' fst_ch : ', fst_ch
   print, ' do_read: ', do_read

;   read_fits_map, '../ext_ctp3_mask_ns16_RING.fits', mask
   read_fits_map, '../sky_cut_20_ns16_ring.fits', mask
   fsky = float(n_elements(where(mask[*,0] ne 0.))) / n_elements(mask[*,0])
   gpix = where(mask[*,0] ne 0.)
   ngpix = n_elements(gpix)

   if (do_read eq 'read') then begin
   nsample = 0l
   iich = 0l
   nsampch = lonarr(nchain)
   for ich=fst_ch,fst_ch+nchain-1 do begin
      if (ich lt 10) then sch = 'c000'+strtrim(string(ich),2)
      if (ich lt 100 and ich ge 10) then sch = 'c00'+strtrim(string(ich),2)
      spawn, 'ls -ltrh '+root+'s_i_'+sch+'_k* | wc -l', nsamp
      nsampch[iich] = long(nsamp)
      nsample = nsample + nsampch[iich]
      iich = iich + 1l
   endfor

   print, ' nsample: ', nsample
   print, nsampch
;;    cl=fltarr(nchain,nsample,6,2 * (ncl+1))
   cl=fltarr(nchain,nsample,nspec,2 * (ncl+1))

   l = findgen(ncl+1)

;;   if (do_read eq 'read') then begin

   foreg = fltarr(nchain, nsample, nforeg, nband)

; --- cls
;   for ich=fst_ch,fst_ch+nchain-1 do for isam=1l,nsample do begin
;       cl_tmp = fltarr(nspec+1, 2*(ncl+1))
;      if (isam lt 10) then ssam = 'k0000'+strtrim(string(isam),2)
;      if (isam lt 100 and isam ge 10) then ssam = 'k000'+strtrim(string(isam),2)
;      if (isam lt 1000 and isam ge 100) then ssam = 'k00'+strtrim(string(isam),2)
;
;      if (ich lt 10) then sch = 'c000'+strtrim(string(ich),2)
;      if (ich lt 100 and ich ge 10) then sch = 'c00'+strtrim(string(ich),2)
;      openr,1,root+'cl_'+sch+'_'+ssam+'.dat'
;      readf,1,cl_tmp
;;      print, cl_tmp[0,0:89]
;      close, 1
;      cl[ich-fst_ch,isam-1,*,*] = reform(cl_tmp[1:*,*])
;
;   endfor

; --- chi2 from the map
;   read_fits_map, '../ext_ctp3_mask_ns16_RING.fits', mask
   read_fits_map, '../sky_cut_20_ns16_ring.fits', mask
   fsky = float(n_elements(where(mask[*,0] ne 0.))) / n_elements(mask[*,0])
   gpix = where(mask[*,0] ne 0.)
   ngpix = n_elements(gpix)
   chi2 = dblarr(nchain,nsample)
   cmb = mask[*,0] * 0.
   cmb2 = cmb

   freq = ['030', '044', '070', '100', '143', '217']
   nside = 16l
   npix = 12l*nside^2
   residuals = fltarr(npix,nband)

;print, do_read
;stop
;   if (do_read eq 'read') then begin

;   for ich=fst_ch,fst_ch+nchain-1 do begin
;      if (ich lt 10) then sch = 'c000'+strtrim(string(ich),2)
;      if (ich lt 100 and ich ge 10) then sch = 'c00'+strtrim(string(ich),2)
;	spawn, 'ls -ltrh '+root+'s_i_'+sch+'_k* | wc -l', nsampch

   iich = -1l
   for ich=fst_ch,fst_ch+nchain-1 do begin
      if (ich lt 10) then sch = 'c000'+strtrim(string(ich),2)
      if (ich lt 100 and ich ge 10) then sch = 'c00'+strtrim(string(ich),2)
;      spawn, 'ls -ltrh '+root+'s_i_'+sch+'_k* | wc -l', nsampch
;      nsampch = long(nsampch)
      iich = iich + 1
      for isam=1l,nsampch[iich] do begin
       cl_tmp = fltarr(nspec+1, 2*(ncl+1))
       fg = fltarr(nforeg, nband)
      if (isam lt 10) then ssam = 'k0000'+strtrim(string(isam),2)
      if (isam lt 100 and isam ge 10) then ssam = 'k000'+strtrim(string(isam),2)
      if (isam lt 1000 and isam ge 100) then ssam = 'k00'+strtrim(string(isam),2)
      if (isam lt 10000 and isam ge 1000) then ssam = 'k0'+strtrim(string(isam),2)
      if (isam lt 100000 and isam ge 10000) then ssam = 'k'+strtrim(string(isam),2)

      if (ich lt 10) then sch = 'c000'+strtrim(string(ich),2)
      if (ich lt 100 and ich ge 10) then sch = 'c00'+strtrim(string(ich),2)
      read_fits_map, root+'chisq_'+sch+'_'+ssam+'.fits', chisq
;      readf,1,cl_tmp
;      print, cl_tmp[0,0:89]
;      close, 1
      openr,1,root+'cl_'+sch+'_'+ssam+'.dat'
      readf,1,cl_tmp
;      print, cl_tmp[0,0:89]
      close, 1
      cl[ich-fst_ch,isam-1,*,*] = reform(cl_tmp[1:*,*])
      chi2[ich-fst_ch,isam-1] = reform(total(chisq,/double)/double(ngpix)/nband)
;mollview, (chisq[*]*mask[*,0])/ngpix
;print, ich, isam, chi2[ich-fst_ch,isam-1]
;stop
      read_fits_map, root+'s_i_'+sch+'_'+ssam+'.fits', map
      cmb = cmb + map
      cmb2 = cmb2 + map^2
;mollview, map
;mollview, cmb
;mollview, cmb2
;stop
      openr, 1, root+'foreground_'+sch+'_'+ssam+'.dat'
      readf,1,fg
      close,1
      foreg[ich-fst_ch,isam-1,*,*] = fg

      for iband=0,nband-1 do begin

         read_fits_map, '../map_freefree_0'+freq[iband]+'.00GHz.smth.ns16.fits', ff
         read_fits_map, '../map_synchrotron_0'+freq[iband]+'.00GHz.smth.ns16.fits', ss
         read_fits_map, '../map_thermal_dust_0'+freq[iband]+'.00GHz.smth.ns16.fits', dd
         read_fits_map, '../noDip_madam_map_'+freq[iband]+'GHz_I_smth600_uK10_ns16.fits', imap
         dipole = make_dipole(reform(fg[1:3,iband]),nside)
         residuals[*,iband] = residuals[*,iband] + (imap - (map + fg[0,iband] + dipole + fg[4,iband]*ff + fg[5,iband]*ss + fg[6,iband]*dd))
        
      endfor
;mollview, residuals
;stop
   endfor
   endfor

   print, nsample

   save, filename=tag+'.sav', foreg, cmb, cl, chi2, cmb2, residuals, nsample, nsampch
   endif

   if (do_read ne 'read') then begin
        print, "restoring..."
        freq = ['030', '044', '070', '100', '143', '217']
        nfreq = n_elements(freq)
	restore, filename=tag+'.sav'
   endif

   !p.multi = 0
   !p.charsize=1.5
   window, 11 & plot, chi2, yr=[0.9*min(chi2),1.1*max(chi2)], ys=1, ytit='!7v!u2!n', xtit='!17Sample'
   chihst =  histogram(chi2, nbins=50, locations=achi);, min=mean(chi2)-4.*stddev(chi2), max=mean(chi2)+4.*stddev(chi2) )
   window, 12 & plot, achi, chihst/total(chihst), ys=2, xtit='!7v!u2!n',ytit='!17P', psym=10;, xtit='!17Sample'
   set_plot, 'ps'
   device, file=tag+'_chi2.eps',/col, bits=8
   plot, chi2, yr=[0.9*min(chi2),1.1*max(chi2)], ys=1, ytit='!7v!u2!n', xtit='!17Sample'
   device, /close
   set_plot, 'x'

   cmb = cmb / nsample
   cmb2 = cmb2 / nsample

   residuals = residuals / nsample
mollview, residuals*mask
;stop
   cmbv = sqrt(cmb2-cmb^2)

   write_fits_map, tag+'_recCMB_ns16_RING.fits',cmb, /ring,units='!7l!17K'
   write_fits_map, tag+'_varCMB_ns16_RING.fits',cmbv, /ring,units='!7l!17K'

;   write_fits_map, tag+'_residuals_ns16_RING.fits',residuals, /ring,units='!7l!17K'
   if (nband eq 5) then resw = create_struct('HDR','Residuals','GHz '+freq[0],reform(residuals[*,0]),'GHz '+freq[1],reform(residuals[*,1]),'GHz '+freq[2],reform(residuals[*,2]),'GHz '+freq[3],reform(residuals[*,3]),'GHz '+freq[4],reform(residuals[*,4]) )
   if (nband eq 6) then resw = create_struct('HDR','Residuals','GHz '+freq[0],reform(residuals[*,0]),'GHz '+freq[1],reform(residuals[*,1]),'GHz '+freq[2],reform(residuals[*,2]),'GHz '+freq[3],reform(residuals[*,3]),'GHz '+freq[4],reform(residuals[*,4]),'GHz '+freq[5],reform(residuals[5]) )

;   write_fits_sb, tag+'_residuals_ns16_RING.fits', 0, resw, ordering='ring', nside=16l;, coord='G', /ring, /nothealpix

   ismoothing,tag+'_recCMB_ns16_RING.fits',tag+'_recCMB_smth600_ns16_RING.fits', fwhm_arc=600
   ismoothing,tag+'_varCMB_ns16_RING.fits',tag+'_varCMB_smth600_ns16_RING.fits', fwhm_arc=600

;   mollview, tag+'recCMB_smth600_ns16_RING.fits', min=-100, max=100, chars=1.5, tit='!17Reconstructed CMB: '+tag
;   mollview, tag+'varCMB_smth600_ns16_RING.fits', chars=1.5, tit='!17STDDEV CMB: '+tag

   read_fits_map, tag+'_recCMB_smth600_ns16_RING.fits', smcmb
   read_fits_map, tag+'_varCMB_smth600_ns16_RING.fits', smcmbv
   mollview, smcmb*mask, min=-100, max=100, chars=1.5, tit='!17Reconstructed CMB: '+tag
   mollview, smcmbv*mask, min=min(smcmbv[where(mask gt 0.)]*mask[where(mask gt 0.)]), chars=1.5, tit='!17Reconstructed CMB: '+tag

;   mollview, tag+'recCMB_smth600_ns16_RING.fits', min=-100, max=100, chars=1.5, tit='!17Reconstructed CMB: '+tag,ps=tag+'_recCMB.eps'
;   mollview, cmb, min=-100, max=100, tit='!17'+tag+' Reconstructed CMB', units='!7l!17K', chars=1.5
;   mollview, cmb*mask[*,0], min=-100, max=100, tit='!17'+tag+' Reconstructed CMB', units='!7l!17K', chars=1.5, ps=tag+'_recCMB.eps',pxs=1600
;   mollview, '../input/camb_cmb_smth_uK_ns16.fits', min=-100, max=100, tit='!17'+tag+' Input CMB', units='!7l!17K', chars=1.5
;   mollview, '../input/camb_cmb_smth_uK_ns16.fits', min=-100, max=100, tit='!17Original CMB', units='!7l!17K', chars=1.5, ps='camb_origCMB_ns16.eps',pxs=1600
   mollview, 'wmap_ilc_7ys_smth600_ns16.fits', min=-100, max=100, factor=1.e3, tit='!17WMAP 7yr', units='!7l!17K', chars=1.5;,pxs=1600
   mollview, 'wmap_ilc_7ys_smth600_ns16.fits', min=-100, max=100, factor=1.e3, tit='!17WMAP 7yr', units='!7l!17K', chars=1.5, ps='wmap7.eps',pxs=1600

   fits2cl, tcl, '../../../input/camb_sim/camb_85760952_scalcls.fits'
   tcl = tcl * 1.e12
;;    ianafast, '../input/camb_85760952_map_ns16.fits', scl, nlmax = 64
;;    ianafast, '../input/camb_cmb_smth_uK_ns16.fits', bscl, nlmax = 64, alm1_out='alms.fits'
   ianafast, '/project/projectdirs/planck/user/iodwyer/CompSep/maps/simmaps/nside16/map_cmb.smth.ns16.fits', bscl, nlmax = 64, alm1_out='alms.fits'
;   bscl = bscl * 1.e12

;fits2alm,index,alms,'alms.fits'
;help, alms
;gb = gaussbeam(600,ncl)
;for il=1,ncl do begin
;l1=il^2+il+1
;l2=il^2+2*il+1
;il1=where(index eq l1)
;il2=where(index eq l2)
;alms[il1:il2,*]=alms[il1:il2,*]/gb[il]
;endfor
;isynfast, bscl, ocmb, alm_in=alms, nside=nside, nlmax=64
;stop
   fits2cl, pwf, '/project/projectdirs/cmb/modules/pdsf/software_gnu/cmb/2.4.0/healpix_2.11c_20090219-2.4.0/data/pixel_window_n0016.fits'
;;    ianafast, '../input/camb_ctp3_CompSep_070GHz_noise_ns16_RING.fits', cl70n, nlmax = 64
   ianafast, '../../../new_ctp3_CompSep_070GHz_noise_ns16_RING.fits', cl70n, nlmax = 64
   ifreq = 2 ;70GHz
   gsig = 2.
   npix = 12l*16l^2
   gnoise = gsig*randomn(-4+ifreq,npix)
   ianafast, gnoise, ncls, nlmax=64, ordering='ring'
;   ianafast, cmb, rcl, /ring, maskfile='../ctp3_mask_ns16_RING.fits'
;   ianafast, cmb, fsrcl, /ring
;   ianafast, '/project/projectdirs/planck/user/iodwyer/CompSep/maps/simmaps/nside16/map_cmb.smth.ns16.fits', tcl,maskfile='../ctp3_mask_ns16_RING.fits'
   gb = gaussbeam(600,ncl)
   l = lindgen(ncl+1)
   ll = l*(l+1)/2./!pi	
;   window, 3 & plot_io, l, tcl*gb^2*pwf^2, xtit='!17l', ytit='!17C!dl!n'
;   oplot, l, bscl*1.e12/gb^2/pwf^2, col=245, line=4
;   oplot, l, bscl*1.e12, col=245, line=4
;   oplot, l, gb^2
;   oplot, l, pwf^2
;   oplot, l, cl70n, col=90, line=2
;   oplot, l, ncls, line=3
;   legend, ['!17Input C!dl!n', '!17FS C!dl'], col=[0,245], line=[0,3], pos=[3,1.e-6]
;   oplot, total(cl[0,*,0,*],2)/nsample/ll
;
;   set_plot, 'ps'
;   device, file=tag+'_spectra.eps',/col,bits=8
;   plot_io, l, tcl*gb^2*pwf^2, xtit='!17l', ytit='!17C!dl!n'
;   oplot, l, bscl*1.e12/gb^2/pwf^2, col=245, line=4
;   oplot, l, bscl, col=245, line=4
;   oplot, l, gb^2
;   oplot, l, pwf^2
;   oplot, l, cl70n, col=90, line=2
;   oplot, l, ncls, line=3
;   legend, ['!17Input C!dl!n', '!17FS C!dl'], col=[0,245], line=[0,3], pos=[3,1.e-6]
;   oplot, total(cl[0,*,0,*],2)/nsample/ll

;plot_io, l, tcl*gb^2*pwf^2, xtit='!17l', ytit='!17C!dl!n'
;   oplot, l, scl*1.e12, col=245, line=3
;   oplot, l, gb^2
;   oplot, l, pwf^2
;   oplot, l, cl70n, col=90, line=2
;   oplot, l, ncls, line=3
;   legend, ['!17Input C!dl!n', '!17FS C!dl'], col=[0,245], line=[0,3]
;   oplot, total(cl[0,*,*,0],2)/nsample
;;    plot, l, tcl*ll/gb^2*1.e12, xtit='!17l', ytit='!17l(l+1)C!dl!n/2!7p!17', tit='!17'+tag+' spectra'
;;    oplot, l, rcl*ll, col=245, line=3
;;    legend, ['!17Input C!dl!n', '!17Rec C!dl'], col=[0,245], line=[0,3]
;   device, /close
;   set_plot, 'x'

   bins=13 ;25 
   hist=fltarr(bins,nband)
   famp = fltarr(bins, nband)

   if (nforeg gt 4) then begin
   for i=0,nband-1 do begin
      hist[*,i]=histogram(foreg[0,*,4,i],nbins=bins,locations=a)
;;       hist[*,i]=histogram(foreg[0,*,i,4],nbins=bins,locations=a)
      famp[*,i] = a
   endfor
   window, 2 & plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist)/total(hist[*,0])], xtit='!17A!dff!n', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=200+i*5,psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
   freq = ['30', '44', '70', '100', '143', '217', '353']
   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
   set_plot, 'ps'
   device, file=tag+'_Famp.eps',/col, bits=8
   plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist)/total(hist[*,0])], xtit='!17A!dff!n', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=200+i*5,psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
   freq = ['30', '44', '70', '100', '143', '217', '353']
   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
   device, /close
   set_plot, 'x'

if (nforeg gt 5) then begin
   famp = fltarr(bins, nband)
   for i=0,nband-1 do begin
      hist[*,i]=histogram(foreg[0,*,5,i],nbins=bins,locations=a)
;;       hist[*,i]=histogram(foreg[0,*,i,5],nbins=bins,locations=a)
      famp[*,i] = a
   endfor
   window, 7 & plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist)/total(hist[*,0])], xtit='!17A!dsynch!n', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=200+i*5,psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
   freq = ['30', '44', '70', '100', '143', '217', '353']
   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
   set_plot, 'ps'
   device, file=tag+'_Samp.eps',/col, bits=8
   plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist)/total(hist[*,0])], xtit='!17A!dsynch!n', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=200+i*5,psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
   freq = ['30', '44', '70', '100', '143', '217', '353']
   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
   device, /close
   set_plot, 'x'
endif

if (nforeg gt 6) then begin
   famp = fltarr(bins, nband)
   for i=0,nband-1 do begin
      hist[*,i]=histogram(foreg[0,*,6,i],nbins=bins,locations=a)
;;       hist[*,i]=histogram(foreg[0,*,i,6],nbins=bins,locations=a)
      famp[*,i] = a
   endfor

   window, 5 & plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17A!dtd!n', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=200+i*5,psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
   freq = ['30', '44', '70', '100', '143', '217', '353']
   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
   set_plot, 'ps'
   device, file=tag+'_Damp.eps',/col, bits=8
   plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17A!dtd!n', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=200+i*5,psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
   freq = ['30', '44', '70', '100', '143', '217', '353']
   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
   device, /close
endif

   set_plot, 'x'
   endif

   for i=0,nband-1 do begin
      hist[*,i]=histogram(foreg[0,*,0,i],nbins=bins,locations=a)
;;       hist[*,i]=histogram(foreg[0,*,i,0],nbins=bins,locations=a)
      famp[*,i] = a
   endfor
   window, 4 & plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17Monopole', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=200+i*5,psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10
   freq = ['30', '44', '70', '100', '143', '217', '353']
   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]

   set_plot, 'ps'
   device, file=tag+'_monopole.eps',/col, bits=8
   plot, famp[*,0],hist[*,0]/total(hist[*,0]), psym=10, yr=[0,max(hist[*,0])/total(hist[*,0])], xtit='!17Monopole', ytit='!17P',xr=[min(famp),max(famp)], chars=1.5
   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/total(hist[*,i]),col=200+i*5,psym=10
;   for i=0,nband-1 do oplot, famp[*,i],hist[*,i]/max(hist[*,i]),col=200+i*5,psym=10   freq = ['30', '44', '70', '100', '143', '217', '353']
   legend, [freq+' GHz'], col=[200, 205, 210, 215, 220, 225, 230], line=[0,0,0,0,0,0,0]
   device, /close
   set_plot, 'x'


;stop
   if (nforeg gt 5) then begin
window, 23 &  plot, foreg[0,*,4,0],foreg[0,*,5,0] ,psym=3, xr=[min(foreg[*,*,4,0]),max(foreg[*,*,4,0])], yr=[min(foreg[*,*,5,0]),max(foreg[*,*,5,0])],xtit='!17Free-Free', ytit='!17Synch' 
for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,0:nsampch[ich-1]/2-1,4,iba],foreg[ich-1,0:nsampch[ich-1]/2-1,5,iba],col=130,psym=3
for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,nsampch[ich-1]/2:*,4,iba],foreg[ich-1,nsampch[ich-1]/2:*,5,iba],col=90,psym=3
;; for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,*,4,iba],foreg[ich-1,*,5,iba],col=90+10*(ich+iba),psym=3
; oplot, foreg[0,*,4,1],foreg[0,*,6,1],col=100,psym=3
; oplot, foreg[0,*,4,3],foreg[0,*,6,3],col=255*2l^8l,psym=3
; oplot, foreg[0,*,4,2],foreg[0,*,6,2],col=250,psym=3

if (nforeg gt 6) then begin
window, 24 & plot, foreg[0,*,4,0],foreg[0,*,6,0] ,psym=3, xr=[min(foreg[*,*,4,0]),max(foreg[*,*,4,0])], yr=[min(foreg[*,*,6,0]),max(foreg[*,*,6,0])],xtit='!17Free-Free', ytit='!17Dust'
for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,0:nsampch[ich-1]/2-1,4,iba],foreg[ich-1,0:nsampch[ich-1]/2-1,6,iba],col=130,psym=3
for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,nsampch[ich-1]/2:*,4,iba],foreg[ich-1,nsampch[ich-1]/2:*,6,iba],col=90,psym=3
;; for ich=fst_ch, fst_ch+nchain-1 do for iba=0, nband-1 do oplot, foreg[ich-1,*,4,iba],foreg[ich-1,*,6,iba],col=90+10*(ich+iba),psym=3
;oplot, foreg[0,*,4,1],foreg[0,*,6,1],col=100,psym=3
;oplot, foreg[0,*,4,3],foreg[0,*,6,3],col=255*2l^8l,psym=3
;oplot, foreg[0,*,4,2],foreg[0,*,6,2],col=250,psym=3


window, 25 & plot, foreg[0,*,5,0],foreg[0,*,6,0] ,psym=3, xr=[min(foreg[*,*,5,0]),max(foreg[*,*,5,0])], yr=[min(foreg[*,*,6,0]),max(foreg[*,*,6,0])],xtit='!17Synch', ytit='!17Dust' 
for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,0:nsampch[ich-1]/2-1,5,iba],foreg[ich-1,0:nsampch[ich-1]/2-1,6,iba],col=130,psym=3
for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,nsampch[ich-1]/2:*,5,iba],foreg[ich-1,nsampch[ich-1]/2:*,6,iba],col=90,psym=3
;; for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,*,5,iba],foreg[ich-1,*,6,iba],col=90+10*(ich+iba),psym=3
;oplot, foreg[0,*,5,1],foreg[0,*,6,1],col=100,psym=3
;oplot, foreg[0,*,5,3],foreg[0,*,6,3],col=255*2l^8l,psym=3
;oplot, foreg[0,*,5,2],foreg[0,*,6,2],col=250,psym=3

endif

set_plot, 'ps'
device, file=tag+'_scatter_1.eps',/col, bits=8
plot, foreg[0,*,4,0],foreg[0,*,5,0] ,psym=3, xr=[min([foreg[0,*,4,0],foreg[0,*,5,0]]),max([foreg[0,*,4,0],foreg[0,*,5,0]])], yr=[min([foreg[0,*,4,0],foreg[0,*,5,0]]),max([foreg[0,*,4,0],foreg[0,*,5,0]])]
for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,0:nsampch[ich-1]/2-1,4,iba],foreg[ich-1,0:nsampch[ich-1]/2-1,5,iba],col=130,psym=3
for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,nsampch[ich-1]/2:*,4,iba],foreg[ich-1,nsampch[ich-1]/2:*,5,iba],col=90,psym=3
;plot, foreg[0,*,4,0],foreg[0,*,5,0] ,psym=3, xr=[0.8,2.], yr=[-0.5,1.5]
;for ich=fst_ch, fst_ch+nchain-1 do oplot, foreg[0,*,4,ich],foreg[0,*,5,ich],col=90+10*ich,psym=3
;oplot, foreg[0,*,4,1],foreg[0,*,5,1],col=90,psym=3
;oplot, foreg[0,*,4,3],foreg[0,*,5,3],col=220,psym=3
;oplot, foreg[0,*,4,2],foreg[0,*,5,2],col=250,psym=3
device, /close


if (nforeg gt 6) then begin
device, file=tag+'_scatter_2.eps',/col, bits=8
plot, foreg[0,*,4,0],foreg[0,*,6,0] ,psym=3, xr=[-1.,3.], yr=[-1.,3.]
for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,0:nsampch[ich-1]/2-1,4,iba],foreg[ich-1,0:nsampch[ich-1]/2-1,6,iba],col=130,psym=3
for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,nsampch[ich-1]/2:*,4,iba],foreg[ich-1,nsampch[ich-1]/2:*,6,iba],col=90,psym=3
;plot, foreg[0,*,4,0],foreg[0,*,6,0] ,psym=3, xr=[0.8,2.], yr=[-0.5,1.5]
;oplot, foreg[0,*,4,1],foreg[0,*,6,1],col=90,psym=3
;oplot, foreg[0,*,4,3],foreg[0,*,6,3],col=220,psym=3
;oplot, foreg[0,*,4,2],foreg[0,*,6,2],col=250,psym=3
device, /close

device, file=tag+'_scatter_3.eps',/col, bits=8
plot, foreg[0,*,5,0],foreg[0,*,6,0] ,psym=3, xr=[-1.,3.], yr=[-1.,3]
for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,0:nsampch[ich-1]/2-1,5,iba],foreg[ich-1,0:nsampch[ich-1]/2-1,6,iba],col=130,psym=3
for ich=fst_ch, fst_ch+nchain-1 do for iba=0,nband-1 do oplot, foreg[ich-1,nsampch[ich-1]/2:*,5,iba],foreg[ich-1,nsampch[ich-1]/2:*,6,iba],col=90,psym=3
;plot, foreg[0,*,5,0],foreg[0,*,6,0] ,psym=3, xr=[-0.8,2.], yr=[-0.5,1.5]
;oplot, foreg[0,*,5,1],foreg[0,*,6,1],col=120,psym=3
;oplot, foreg[0,*,5,3],foreg[0,*,6,3],col=220,psym=3
;oplot, foreg[0,*,5,2],foreg[0,*,6,2],col=250,psym=3
device, /close
endif

 set_plot, 'x'

   endif


;;   clen = min([200,nsample])
;;   corr=dblarr(clen,nchain,nspec,ncl)
;;   ave = total(cl,2) / nsample
;;	
;;   for n=0,clen-2 do begin
;;      indx = lindgen(nsample-n)
;;      for ich=fst_ch,fst_ch+nchain-1 do for kcl=0,nspec-1 do for icl=2,ncl-1 do corr[n,ich-fst_ch,kcl,icl] = mean( (cl[ich-fst_ch,indx,kcl,icl]-ave[ich-fst_ch,kcl,icl])*(cl[ich-fst_ch,indx+n,kcl,icl]-ave[ich-fst_ch,kcl,icl]) ) / variance(cl[ich-fst_ch,indx,kcl,icl])
;;   endfor
;;
;;  window, 27& plot, corr[*,0,0,4]
;;
   ttcl = fltarr(nsample, 2*(ncl+1))
   i = 0l
   for ichain=0,nchain-1 do for isample=0l,nsampch[ichain]-1 do begin
       ttcl[i,*] = cl[ichain, isample, 0, *]
       i = i + 1
   endfor

   bins = 50
   htt = fltarr(bins,2*(ncl+1))
   tta = fltarr(bins,2*(ncl+1))
   err = fltarr(2*(ncl+1))
   for icl=0, 2*(ncl+1)-1 do begin
       err[icl] = stddev(ttcl[*,icl])
       htt[*,icl] = histogram(ttcl[*,icl], locations=a, nbins=bins, max=mean(ttcl[*,icl])+4.*stddev(ttcl[*,icl]),min=mean(ttcl[*,icl])-4.*stddev(ttcl[*,icl]))
;;        htt[*,icl] = histogram(ttcl[*,icl], locations=a, nbins=bins)
       tta[*, icl] = a
   endfor

   gb = gaussbeam(600,ncl)

   !p.multi = 0
   !p.charsize=1.5
   window, 6 & 
   !p.multi = [0,4,4]
   lplot = [2,3,4,6,8,10,12,14,16,18,20,24,28,32];,42,52,62]
   fscl = bscl*1.e12/gb^2/pwf^2

   readcol, '../../wmap_tt_spectrum_7yr_v4p1.txt', wmapl, wmapcl 

   for il=0,n_elements(lplot)-1 do begin
       plot, tta[*, lplot[il]],htt[*, lplot[il]]/total(htt[*, lplot[il]]), ytit='!17P', xtit='!17C!d'+strtrim(string(lplot[il]),2), psym=10, chars=1.5
;;        oplot, reform(lplot[*]*0.+fsrcl[il]*lplot[il]*(lplot[il]+1)/2./!pi), reform(lplot[*]*0.+findgen(n_elements(lplot))/n_elements(lplot) ), col=245
;        oplot, reform(lplot[*]*0.+rcl[il]*lplot[il]*(lplot[il]+1)/2./!pi), reform(lplot[*]*0.+findgen(n_elements(lplot))/n_elements(lplot) ), col=90
;        oplot, reform(lplot[*]*0.+fscl[lplot[il]]*lplot[il]*(lplot[il]+1)/2./!pi), reform(lplot[*]*0.+findgen(n_elements(lplot))/n_elements(lplot)),col=245, line=2
        oplot, reform(lplot[*]*0.+wmapcl[lplot[il]-2]), reform(lplot[*]*0.+findgen(n_elements(lplot))/n_elements(lplot)),col=245, line=2
;       oplot, tta[*, ncl+1+lplot[il]],htt[*, ncl+1+lplot[il]]/total(htt[*, ncl+1+lplot[il]]), col=240, psym=10
   endfor
   write_png,tag+'_ctp3_cls.png',tvrd(/true)

   set_plot, 'ps'
   device, file=tag+'_ctp3_cls.eps', /col, bits=8, /encapsulated    
   for il=0,n_elements(lplot)-1 do begin
       plot, tta[*, lplot[il]],htt[*, lplot[il]]/total(htt[*, lplot[il]]), ytit='!17P', xtit='!17C!d'+strtrim(string(lplot[il]),2), psym=10
;;        oplot, reform(lplot[*]*0.+fsrcl[il]*lplot[il]*(lplot[il]+1)/2./!pi), reform(lplot[*]*0.+findgen(n_elements(lplot))/n_elements(lplot)),col=245
;        oplot, reform(lplot[*]*0.+rcl[il]*lplot[il]*(lplot[il]+1)/2./!pi), reform(lplot[*]*0.+findgen(n_elements(lplot))/n_elements(lplot) ), col=90
        oplot, reform(lplot[*]*0.+fscl[lplot[il]]*lplot[il]*(lplot[il]+1)/2./!pi), reform(lplot[*]*0.+findgen(n_elements(lplot))/n_elements(lplot)),col=245, line=2
;       oplot, tta[*, ncl+1+lplot[il]],htt[*, ncl+1+lplot[il]]/total(htt[*, ncl+1+lplot[il]]), col=240, psym=10
   endfor
device, /close
set_plot, 'x'

save, filename=tag+'_histCls.sav', htt, tta

!p.multi=0


   gb = gaussbeam(600,ncl)
   l = lindgen(ncl+1)
   ll = l*(l+1)/2./!pi
;;    window, 31 & plot_io, l, tcl*gb^2*pwf^2, xtit='!17l', ytit='!17C!dl!n'
   window, 31 & plot_io, wmapl[0:ncl-2], wmapcl/wmapl/(wmapl+1)*2.*!pi, xtit='!17l', ytit='!17C!dl!n'
   oplot, wmapl, wmapcl/wmapl/(wmapl+1)*2.*!pi, thick=2, col=210
;   oplot, l, tcl*gb^2*pwf^2, col=0
;   oplot, l, bscl*1.e12/gb^2/pwf^2, col=245, line=4
;   oplot, l, bscl*1.e12, col=245, line=4
;   oplot, l, gb^2
;   oplot, l, pwf^2
;   oplot, l, cl70n, col=90, line=2
   oplot, l, ncls, line=3
   legend, ['!17Mean C!dl!n', '!17FS C!dl', '!17C!dl!n mode', '7yr WMAP'], col=[0,245,70,210], line=[0,3,2,0], pos=[3,1.e-6]
   oplot, total(total(cl[*,*,0,*],2),1)/nsample/ll, line=1


ianafast, '../../wmap_ilc_7ys_smth600_ns16.fits', wmapclsmoo, nlmax=ncl
oplot, l[2:*], wmapclsmoo[2:*]*1.e6 / gb[2:*]^2 / pwf[2:*]^2, thick=2, col=245
oplot, l, l*0.+10.^2*4.*!pi/(12l*16l^2)/gb^2

   clmode = fltarr(ncl+1)
   for i=0l,ncl do clmode[i] = tta[where(htt[*,i] eq max(htt[*,i])), i]
   oplot, l, clmode/ll, col=70, line=2
   errplot, (clmode-err)/ll, (clmode+err)/ll

oplot, l, l*0.+10.^2*4.*!pi/(12l*16l^2)

;;    oplot, l, (clmode-mean(clmode[25:*]) )/ll, line=2

   write_png,tag+'_convCls.png',tvrd(/true)


   window, 30 & plot_io, wmapl[0:ncl-2], wmapcl[0:ncl-2], xtit='!17l', ytit='!17l(l+1)C!dl!n/2!7p!17'
   oplot, wmapl[0:ncl-2], wmapcl[0:ncl-2], thick=2, col=210
;   oplot, l, tcl*gb^2*pwf^2, col=0
;   oplot, l, bscl*1.e12/gb^2/pwf^2, col=245, line=4
;   oplot, l, bscl*1.e12, col=245, line=4
;   oplot, l, gb^2
;   oplot, l, pwf^2
;   oplot, l, cl70n, col=90, line=2
   oplot, l, ncls, line=3
   legend, ['!17Mean C!dl!n', '!17C!dl!n mode', '7yr WMAP'], col=[0,70,210], line=[0,2,0];, pos=[3,1.e-6]
   oplot, l, total(total(cl[*,*,0,*],2),1)/nsample, line=1
   errplot, clmode-err, clmode+err
;clmode = fltarr(ncl+1)
;  for i=0,ncl do clmode[i] = tta[where(htt[*,i] eq max(htt[*,i])), i]
   oplot, l, clmode, col=70, line=2
   write_png,tag+'_convCls_2.png',tvrd(/true)


print, " Old piece of code from hereafter"
print, " --- End of Code ---"

stop
   ch_ave = total(cl, 2) / nsample

   window, 0 & plot, l[0:89], ch_ave[0,0,0:89]
   oplot, l[0:89], ch_ave[1,0,0:89], col=200
   oplot, l, ch_ave[2,0,0:89], col=210
   oplot, l, ch_ave[3,0,0:89], col=220

   oplot, l, ch_ave[0,0,90:*], line=2
   oplot, l, ch_ave[1,0,90:*], col=200, line=2
   oplot, l, ch_ave[2,0,90:*], col=210, line=2
   oplot, l, ch_ave[3,0,90:*], col=220, line=2

   cl_ave = total(ch_ave,1) / nchain

   oplot, l, cl_ave[0,0:89], thick=2
   oplot, l, cl_ave[0,90:*], thick=2, col=200

   
     nside = 64l
     npix = 12l*nside^2
     tot_chi2 = fltarr(nchain,nsample)

     for ichain=1, nchain do begin
        if (ichain lt 10) then sch='000'+strtrim(string(ichain),2)
        if (ichain ge 10 and ichain lt 100) then sch='00'+strtrim(string(ichain),2)
        for isample=1, nsample do begin
           if (isample lt 10) then ssam = '0000'+strtrim(string(isample),2)
           if (isample ge 10 and isample lt 100) then ssam = '000'+strtrim(string(isample),2)
           if (isample ge 100 and isample lt 1000) then ssam = '00'+strtrim(string(isample),2)
           if (isample ge 1000 and isample lt 10000) then ssam = '0'+strtrim(string(isample),2)
           if (isample ge 10000 and isample lt 100000) then ssam = strtrim(string(isample),2)
           file = root+'chisq_c'+sch+'_k'+ssam+'.fits'
           read_fits_map, file, chi2
           tot_chi2[ichain-1,isample-1] = total(chi2) / npix
        endfor
     endfor

     window, 1 & plot, /nodata, [0,nsamples], [0, 2]
     for ichain=0, nchain-1 do oplot, tot_chi2[ichain,*], col=200+ichain

 




stop

end
