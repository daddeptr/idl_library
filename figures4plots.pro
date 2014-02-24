   mollview, randomn(-1,12), window=-1
   loadct, 39
   !p.color=0
   !p.background=255

   freq=[30,44,70,100,143,217,353]

   sfreq = string(freq, format='(i3.3)')

   nfreq = n_elements(freq)


; ------
if (1b) then begin
restore,'chains/cls_hp.res.sav'
readcol,'../ns2048/chains/xfaster/output/xfaster_dx6_5arcmin_TT',x,y,er
bin_str='/project/projectdirs/planck/user/paganol/xfaster/inputs/binning/CTP_bin_TT'
readcol,'../../wmap/wmap_binned_tt_spectrum_7yr_v4p1.txt',wl,wcl,wer,format='f,x,x,f,f'
st=fltarr(4,502)
for it=0l,501 do st[*,it] = moment(cl[*,0,it])
cmd_bin = bp_binning( st[0,251:*], bin_str ) 
cmd_er = bp_binning( sqrt( st[1,251:*] ), bin_str )
window, 1, xsize=1344, ysize=840  & plot_oo, x,y,chars=1.5, xtit='!8l',ytit='!8l(l+1)C!dl!n/2!7p!8', xr=[1,250], yr=[10,6500], thick=2, ys=1, xs=1
oplot, wl, wcl,col=160,thick=4
errplot, x*0.99,y-er, y+er, thick=2
oplot, x[0:50], cmd_bin[0:50], col=245, thick=2
errplot, x[0:50]*1.01, cmd_bin[0:50]-cmd_er[0:50],cmd_bin[0:50]+cmd_er[0:50], col=245, thick=2
;; oplot, wl, wcl, col=160, thick=3
;errplot,wl*1.02, wcl-wer, wcl+wer, col=160, thick=2
legend,['Xfaster PS', 'Commander PS','WMAP PS'], col=[0,245,160], line=[0,0,0], pos=[2, 5000], chars=1.5
write_png, 'png/ps_comp.png', tvrd(/true)
stop
endif

if (0b) then begin
;fits2cl,nls, '../ns2048/chains/dcls_5arcmin_apod.fits'
readcol,'../ns2048/chains/xfaster/output/xfaster_dx6_5arcmin_w7bin_TT',xw,yw,erw
readcol,'../../wmap/wmap_binned_tt_spectrum_7yr_v4p1.txt',wl,wcl,wer,format='f,x,x,f,f'
window, 3, xsize=1344,ysize=840 & plot, xw,yw, chars=2,thick=2, xtit='!8l',ytit='!8l(l+1)C!dl!n/2!7p!8',xr=[1,1200],xs=1
errplot,xw*0.995,yw-erw,yw+erw, thick=2
oplot, wl, wcl, col=245, thick=2
errplot,wl*1.005, wcl-wer, wcl+wer, col=245, thick=2
legend, ['DX6 Xfaster PS', 'WMAP 7yr PS'], col=[0,245], line=[0,0], pos=[800, 5800], chars=1.5
;write_png,'png/dx6_vsw7.png',tvrd(/true)
stop

window, 4, xsize=1344,ysize=840 & plot_oi, xw,yw, chars=2,thick=2, xtit='!8l',ytit='!8l(l+1)C!dl!n/2!7p!8',xr=[1,1200],xs=1
errplot,xw*0.99,yw-erw,yw+erw, thick=2
oplot, wl, wcl, col=245, thick=2
errplot,wl*1.01, wcl-wer, wcl+wer, col=245, thick=2
legend, ['DX6 Xfaster PS', 'WMAP 7yr PS'], col=[0,245], line=[0,0], pos=[2, 5800], chars=1.5
write_png,'png/dx6_vsw7_log.png',tvrd(/true)
stop
endif

; ------

; ------
   if (0b) then begin
   restore, 'chains/cls_hp_99.res.sav'
   bad_val = -1.6375e30
   bp = where(mask eq 0.)
   gp = where(mask eq 1.)
   ccc = cmb
   ccc[bp] = bad_val
   mollview, ccc, chars=1.5, min=-300, max=300, tit='Posterior Average of CMB', png=figs_dir+'cmb_0.99.png'
endif
; ------


   root = 'cls_hp'
   cls_samp = 'on'
   figs_dir = './png/'

   restore, 'chains/'+root+'.res.sav'

   bad_val = -1.6375e30
   bp = where(mask eq 0.)
   gp = where(mask eq 1.)

    sind = ['!17Dust Emissivity', '!17Dust Temperature', '!17Low Frequency Spectral Index']
    samp = ['!17Foreground Amplitude @ 353 GHz', '!17Foreground Amplitude @ 30GHz', 'CO Amplitude @ 100GHz']

    fgind = ['dust_emis', 'dust_t', 'beta']
    fgamp = ['dust_amp', 'sync_amp', 'co_amp']
    
    if (0b) then begin
        ccc= cmb
        ccc[bp] = bad_val
        mollview, ccc, chars=1.5, tit='!17CMB Posterior Average', units='!7l!8K CMB',png=figs_dir+'cmb.png'

        ccc= sqrt(cmb2-cmb^2)
        ccc[bp] = bad_val

        mollview, ccc, chars=1.5, tit='!17CMB RMS', units='!7l!8K CMB',png=figs_dir+'cmb_rms.png'
    endif

    if (0b) then begin
        read_fits_map, 'chains/cls_hp_recCMB_smth90_ns128_RING.fits', smcmb
        for ifreq=0,nfreq-1 do begin
            read_fits_map, 'dx6_map_'+sfreq[ifreq]+'.fits', map
            mollview, map-smcmb, chars=1.5,units='!7l!8K CMB', tit='!17Foregrounds @ '+sfreq[ifreq]+' GHz', /hist, png=figs_dir+'fg_'+sfreq[ifreq]+'.png',window=-1
        endfor
    endif


;    ccc= glb_chi2
;    ccc[bp] = bad_val
;    mollview, ccc, chars=1.5, tit='!8Mean !7v!8!u2', png=figs_dir+'chisq.png'

if (0b) then begin
    for ifreq=0,nfreq-1 do begin
        mollview, 'dx6_map_'+sfreq[ifreq]+'.fits', chars=1.5, units='!7l!8K CMB', tit='!17Sky @ '+sfreq[ifreq]+' GHz - 90 arcmin',/asinh, png=figs_dir+'map_'+sfreq[ifreq]+'.png', window=-1
        mollview, 'dx6_rms_'+sfreq[ifreq]+'.fits', /asinh, chars=1.5,units='!7l!8K CMB', tit='!17RMS @ '+sfreq[ifreq]+' GHz - 90 arcmin', png=figs_dir+'rms_'+sfreq[ifreq]+'.png', window=-1
    endfor
endif

if (0b) then begin
fits2cl,nls, '../ns2048/chains/dcls_5arcmin_apod.fits'
l=findgen(2048)
ll=l*(l+1)/2./!pi
fsky = 0.75
readcol,'../ns2048/chains/xfaster/output/xfaster_dx6_5arcmin_TT',x,y,er
readcol,'../../wmap/wmap_binned_tt_spectrum_7yr_v4p1.txt', wl, wcl, wer, format='f,x,x,f,f'

window, 1, xsize=1344, ysize=840  & plot, x,y,chars=1.5, xtit='!8l',ytit='!8l(l+1)C!dl!n/2!7p!8', xr=[1,1600], yr=[0,6500], thick=2, ys=1, xs=1
errplot, x,y-er, y+er
oplot, wl, wcl, col=245, thick=2
oplot, nls*ll/fsky, col=70, thick=2
legend, ['DX6 Xfaster PS', 'WMAP 7yr PS','CMB Map Noise Power'], col=[0,245,70], line=[0,0,0], pos=[1100, 6200], chars=1.5
write_png,figs_dir+'dx6_ps.png', tvrd(/true)

window, 2, xsize=1344, ysize=840 & plot_oi, x,y,chars=1.5, xtit='!8l',ytit='!8l(l+1)C!dl!n/2!7p!8', xr=[1,1600], yr=[0,6500], thick=2, ys=1, xs=1
errplot, x,y-er, y+er
oplot, wl, wcl, col=245, thick=2
oplot, nls*ll/fsky, col=70, thick=2
legend, ['DX6 Xfaster PS', 'WMAP 7yr PS','CMB Map Noise Power'], col=[0,245,70], line=[0,0,0], pos=[2, 6200], chars=1.5
write_png,figs_dir+'dx6_ps_log.png', tvrd(/true)

;window, 3 & plot, wl, wcl,chars=1.5, xtit='!8l',ytit='!8l(l+1)C!dl!n/2!7p!8', xr=[1,1600], yr=[0,6500], thick=2, ys=1, xs=1
;errplot, wl, wcl-wer, wcl+wer, thick=2
;oplot, x, y, col=245, thick=2
;legend, ['DX6 Xfaster PS', 'WMAP 7yr PS'], col=[245,0], line=[0,0], pos=[1100, 6200], chars=1.5
;write_png,figs_dir+'w7_ps.png', tvrd(/true)
;
;window, 4 & plot_oi, wl, wcl,chars=1.5, xtit='!8l',ytit='!8l(l+1)C!dl!n/2!7p!8', xr=[1,1600], yr=[0,6500], thick=2, ys=1, xs=1
;errplot, wl, wcl-wer, wcl+wer, thick=2
;oplot, x, y, col=245, thick=2
;legend, ['DX6 Xfaster PS', 'WMAP 7yr PS'], col=[245,0], line=[0,0], pos=[2, 6200], chars=1.5
;write_png,figs_dir+'w7_ps_log.png', tvrd(/true)

endif
stop


   if (0b) then begin
      ilabs = [' ','K',' ']
    for iind=0,nind[0]-1 do begin
        ccc = ind[*,iind]
        mn = moment(ccc[gp])
        ccc[bp] = bad_val
        mollview, ccc, chars=1.5, min=mn[0]-3.*sqrt(mn[1]), max=mn[0]+3.*sqrt(mn[1]), tit=sind[iind], units=ilabs[iind], png=figs_dir+fgind[iind]+'.png', window=-1
       shist = histogram(ccc[gp], nbins=50, locations=a)
       window,0 & plot, a,shist, chars=1.5, ytit='!17P', psym=10, thick=2, tit=sind[iind]
       write_png,figs_dir+fgind[iind]+'_hist.png', tvrd(/true)
;stop
    endfor

endif

   if (0b) then begin

    for iamp=0,namp[0]-2 do begin
        ccc=amp[*,iamp]
        ccc[bp] = bad_val
        mollview, ccc, chars=1.5, min=0, max=500,units='!7l!8K !17Antenna', tit=samp[iamp], png=figs_dir+fgamp[iamp]+'.png', window=-1
       shist = histogram(ccc[gp], nbins=50, locations=a)
       window,1 & plot, a,shist, chars=1.5, ytit='!17P', psym=10, thick=2, tit=samp[iamp]
       write_png,figs_dir+fgamp[iamp]+'_hist.png', tvrd(/true)
;stop
    endfor

endif

; ------ Skewness and Kurtosis
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

; ------

        if (0b) then begin
           ccc = sqrt(cmbvar)
           ccc[bp] = bad_val
           mollview, ccc, chars=1.5, tit='CMB RMS', png=figs_dir+'cmb_rms.png', window=-1
        endif

   if (0b) then begin
   for iamp=1,namp[0]-1 do begin
       ccc = sqrt(ampvar[*,iamp-1])
       ccc[bp] = bad_val
       mollview, ccc, chars=1.5, tit=sAmp[iamp-1]+' RMS', png=figs_dir+fgAmp[iamp]+'_rms.png', window=-1
    endfor
    endif

   if (0b) then begin
      ilabs = [' ','K',' ']
   for iind=1,nind[0] do begin
       ccc = sqrt(indvar[*,iind-1])
       ccc[bp] = bad_val
       mollview, ccc, chars=1.5, tit=sInd[iind-1]+' RMS', units=ilabs[iind-1], png=figs_dir+fgInd[iind-1]+'_rms.png', window=-1
       if (iind eq nind[0]) then mollview, ccc, chars=1.5, min=.25, tit=sInd[iind-1]+' RMS', png=figs_dir+fgInd[iind-1]+'_rms.png', window=-1
    endfor
    endif

stop
; ------
stop

if (0b) then begin
        s_ns = strtrim(string(nside),2)

        npix = 12l*nside^2

        print, "restoring..."
        restore, filename='chains/'+root+'.samples_stats.sav'
;   nchains = n_elements(nsampch)
        cmb3n = cmb3 / sqrt(cmbvar)^3
        cmb4n = cmb4 / cmbvar^2 - 3.

        amp3n = amp3 / sqrt(ampvar)^3
        amp4n = amp4 / ampvar^2 - 3.

        ind3n = ind3 / sqrt(indvar)^3
        ind4n = ind4 / indvar^2 - 3.

        corrmapn = corrmap * 0.

        corr_ipix = [10,npix/2,npix-1]
        ncpix = n_elements(corr_ipix)
        for icpix=0,ncpix-1 do corrmapn[*,icpix] = corrmap[*,icpix] / sqrt(cmbvar*cmbvar[corr_ipix[icpix]])

        ccc = cmb3n
        ccc[bp] = bad_val
;        mollview, ccc, chars=1.5, min=-.5, max=.5, tit='!17CMB map skewness'
   mollview, cmb3n, chars=1.5, min=-.5, max=.5, tit='!17CMB Map Skewness', png=figs_dir+'CMBskewness.png', window=-1

        ccc = cmb4n
        ccc[bp] = bad_val
;        mollview, ccc, chars=1.5, min=-.5, max=.5, tit='CMB map kurtosis'
   mollview, ccc, chars=1.5, min=-.5, max=.5, tit='CMB Map Kurtosis', png=figs_dir+'CMBkurtosis.png', window=-1


;;     sAmp = ['Dust','Sync','CO']
   for iamp=1,namp[0]-1 do begin
       ccc = amp3n[*,iamp-1]
       ccc[bp] = bad_val
;       mollview, ccc, chars=1.5, min=-.5, max=.5, tit=sAmp[iamp-1]+' Skewness'
       mollview, ccc, chars=1.5, min=-.5, max=.5, tit=sAmp[iamp-1]+' Skewness', png=figs_dir+fgAmp[iamp]+'_skew.png', window=-1
       shist = histogram(ccc[gp], nbins=50, locations=a)
       window,0 & plot, a,shist, chars=1.5, ytit='!17P', psym=10, thick=2, tit='!17Skewess Histogram: '+sAmp[iamp-1]
       write_png,figs_dir+fgamp[iamp-1]+'_skew_hist.png', tvrd(/true)

       ccc = amp4n[*,iamp-1]
       ccc[bp] = bad_val
;       mollview, ccc, chars=1.5, min=-.5, max=.5, tit=sAmp[iamp-1]+' Kurtosis'
       mollview, ccc, chars=1.5, min=-.5, max=.5, tit=sAmp[iamp-1]+' Kurtosis', png=figs_dir+fgAmp[iamp]+'_kurt.png', window=-1
       shist = histogram(ccc[gp], nbins=50, locations=a)
       window,1 & plot, a,shist, chars=1.5, ytit='!17P', psym=10, thick=2, tit='!17Kurtosis Histogram: '+sAmp[iamp-1]
       write_png,figs_dir+fgamp[iamp-1]+'_kurt_hist.png', tvrd(/true)
;stop
   endfor
                        


;;    sInd = ['Emissivity','Td','Beta']
   for iind=1,nind[0] do begin
       ccc = ind3n[*,iind-1]
       ccc[bp] = bad_val
;       mollview, ccc, chars=1.5, min=-1, max=1, tit=sInd[iind-1]+' Skewness'
       mollview, ccc, chars=1.5, min=-1, max=1, tit=sInd[iind-1]+' Skewness', png=figs_dir+fgInd[iind-1]+'_skew.png', window=-1

       shist = histogram(ccc[gp], nbins=50, locations=a)
       window,0 & plot, a, shist, chars=1.5, ytit='!17P', psym=10, thick=2, tit='!17Skewess Histogram: '+sind[iind-1]
       write_png,figs_dir+fgind[iind-1]+'_skew_hist.png', tvrd(/true)

       ccc = ind4n[*,iind-1]
       ccc[bp] = bad_val
;       mollview, ccc, chars=1.5, min=-1, max=1, tit=sInd[iind-1]+' Kurtosis'
       mollview, ccc, chars=1.5, min=-1, max=1, tit=sInd[iind-1]+' Kurtosis', png=figs_dir+fgInd[iind-1]+'_kurt.png', window=-1

       shist = histogram(ccc[gp], nbins=50, locations=a)
       window,1 & plot, a,shist, chars=1.5, ytit='!17P', psym=10, thick=2, tit='!17Kurtosis Histogram: '+sind[iind-1]
       write_png,figs_dir+fgind[iind-1]+'_kurt_hist.png', tvrd(/true)
;stop
   endfor

   rot_pos = fltarr(ncpix,3)
   rot_pos[0,*] = [0,90,0] 
   rot_pos[1,*]= [0,-180,0]
   rot_pos[2,*]= [0,-90,0]

   grat_pos = fltarr(ncpix,2)
   grat_pos[0,*] = [15,1]
   grat_pos[1,*]= [1,1]
   grat_pos[2,*]= [15,1]


       for icpix=0,ncpix-1 do begin
           mollview, corrmapn[*,icpix], chars=1.5, min=-.5, max=.5, tit='!17Correlation map: pix No '+strtrim(string(corr_ipix[icpix]),2), rot=reform(rot_pos[icpix,*]), png=figs_dir+'corr_mv_'+strtrim(string(corr_ipix[icpix]),2)+'.png', window=-1
;          mollview, corrmapn[*,icpix], chars=1.5, min=-1, max=1, tit='Correlation map: '+corr_ipix[icpix], png=fig_dir+'corrmap_'+strtrim(string(corr_ipix[icpix]),2)+'.png'

           gnomview, corrmapn[*,icpix], chars=1.5, min=-.5, max=.5, rot=reform(rot_pos[icpix,*]), reso=2.5, grat=reform(grat_pos[icpix,*]), /asinh, tit='!17Zoom', png=figs_dir+'corr_gv_'+strtrim(string(corr_ipix[icpix]),2)+'.png', window=-1
       endfor

col_array = [0,70,245]

nb = 25
cd = fltarr(nb,ncpix)
acd = fltarr(nb, ncpix)
window,1 & plot, /nodata, [-4,4], [0,0.2], chars=2, xtit='!7r!8', ytit='!8P'
for icpix=0,ncpix-1 do begin
    h = (cmbdist[*,icpix]-avecmb[corr_ipix[icpix]]) / sqrt(cmbvar[corr_ipix[icpix]])
    cd[*,icpix] = histogram(h, nbins=nb, locations=a)
    acd[*,icpix] = a
    oplot, acd[*,icpix], cd[*,icpix]/nsample, col=col_array[icpix], thick=2, psym=10
endfor
 legend,'Pix '+string(corr_ipix), col=col_array[0:ncpix-1], line=lindgen(ncpix)*0., chars=1.2
 write_png,figs_dir+'pix_hist.png', tvrd(/true)



endif

 print, " ---End of Code ---"



stop
end

