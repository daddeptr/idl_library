; WMAP 7yr

do_sim = 1b

sfreq = ['Ka','K','Q','V','W']
nfreq = n_elements(sfreq)

sig_mK = [1.437, 1.470, 2.197, 3.137, 6.549]
sig_uK = sig_mK * 1.e3

down_beam = 60.
beam_arcmin = 60. * [0.88, 0.66, 0.51, 0.35, 0.22]

eff_beam = sqrt(down_beam^2-beam_arcmin^2)

print, eff_beam

root = '../combined_freq/wmap_band_iqumap_r9_7yr_'

Nsim = 5000

hNside = 512l
hNpix = 12l*hNside^2

Nside = 128l
Npix = 12l*Nside^2

Nobs = fltarr(hnpix,nfreq)
for ifreq=0,nfreq-1 do begin
   read_fits_map,root+sfreq[ifreq]+'_v4.fits',map, hdr, xhdr
   Nobs[*,ifreq] = map[*,3]
endfor

iseed = -1

ave = fltarr(Npix, nfreq)
ave2 = fltarr(Npix, nfreq)


if (do_sim) then begin

for isim=0l,Nsim-1 do begin
;   print, isim
   for ifreq=0,nfreq-1 do begin
      noise = randomn(iseed,hNpix) * sig_uK[ifreq] / sqrt(Nobs[*,ifreq])
      ismoothing, noise, smoo_noise, /nest, fwhm_arcmin=down_beam, /silent
      ud_grade, smoo_noise, low_noise, order_in='nest', order_out='ring', nside_out=Nside
      ave[*,ifreq] = ave[*,ifreq] + low_noise
      ave2[*,ifreq] = ave2[*,ifreq] + low_noise^2
   endfor
   if ( (isim/50)*50 eq isim) then begin
print, 'OK ', isim
	tave = ave / (isim+1)
	tave2 = ave2 / (isim+1)

rms = sqrt(tave2 - tave^2)
for ifreq=0,nfreq-1 do write_fits_map,'noise_sims/w7_rms_smoo'+strtrim(string(down_beam),2)+'_Nsim'+strtrim(string(Nsim),2)+'_ns'+strtrim(string(Nside),2)+'_'+sfreq[ifreq]+'_uK_ring.fits', rms[*,ifreq], /ring, units='!7l!17K CMB'

   endif
endfor

ave = ave / Nsim
ave2 = ave2 / Nsim

rms = sqrt(ave2 - ave^2)

for ifreq=0,nfreq-1 do write_fits_map,'noise_sims/w7_rms_smoo'+strtrim(string(down_beam),2)+'_Nsim'+strtrim(string(Nsim),2)+'_ns'+strtrim(string(Nside),2)+'_'+sfreq[ifreq]+'_uK_ring.fits', rms[*,ifreq], /ring, units='!7l!17K CMB'

endif

for ifreq=0,nfreq-1 do mollview, 'noise_sims/w7_rms_smoo'+strtrim(string(down_beam),2)+'_Nsim'+strtrim(string(Nsim),2)+'_ns'+strtrim(string(Nside),2)+'_'+sfreq[ifreq]+'_uK_ring.fits', units='!7l!17K CMB', tit='RMS '+sfreq[ifreq]+' - '+strtrim(string(long(down_beam)),2)+' arcmin', chars=1.5, ps='rms_'+sfreq[ifreq]+'_'+strtrim(string(long(down_beam)),2)+'.eps', window=-1

end
