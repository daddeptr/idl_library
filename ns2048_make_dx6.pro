; High resolution version of ns128 script:
; Only two HFI frequencies are retrieved according to the setup agreed
; (mixing matrix derived at low resolution only)
;
; !!! - WARNING - !!!
; Missing pixel filled in with median filter value
;
; Modified to handle dx6 dataset
; HFI beams as dx4 ?!?!
;
; Specific for ns=128 only; 
; Uses median filter to handle missing pixels.
; Simple version: up/down grading maps to nside=1024 to run commander.
; LFI: ready
; HFI:  one down
;
;
;; 
;pro make_dx4, ns=ns, beams=beams;, gsig=gsig

;   if not (keyword_set(ns)) or not (keyword_set(beams)) then begin
;      print, 'usage: make_rdata, ns=, beams='
;      stop
;   endif
;

   do_lfi = 0b
   do_hfi = 1b

   do_fs = 1b
   do_1h = 1b
   do_2h = 1b

   do_5arcmin = 0b

   ff = 3 
   lf = 3

   ns = 2048l
;;    beams = 90.
;;    gsig = 3.

   freq=[30, 44, 70, 100, 143, 217, 353]
   sfreq=['30', '44', '70', '100', '143', '217', '353']
;   sfreq=string(freq, format='(1i3.3)')
   nfreq = n_elements(freq)

;   real_beam = [60.*0.54424700, 60.*0.46527683, 60.*0.21676942, 9.94,   7.04,    4.66,    4.41,    4.47,    4.23]

   real_beam = [32.65, 27.00, 13.01, 9.94,   7.04,    4.66,    4.41,    4.47,    4.23]

;;    s_beams = strtrim(string(long(beams)),2)
;;    act_beam = sqrt(float(beams)^2-real_beam^2)

;;    eff_beams = sqrt(float(beams)^2-real_beam^2)

;;    print, eff_beams

;;    dir = '/project/projectdirs/planck/user/planck/dx4/'

   dir = '/project/projectdirs/planck/data/mission/DPC_maps/dx6/'

   out_dir = './'

   npix = 12l*ns^2
   
   bad_value = -1.63750e+30

   if (do_lfi) then begin

   for ifreq=0, 2 do begin
       print, sfreq[ifreq]

       infile = dir+'lfi/LFI_'+sfreq[ifreq]+'_1024_20110304.fits'
       print, ' reading: ', infile
       read_fits_map, infile, map, hdr, xhdr

       bp = where(map[*,0] eq bad_value)
       gp = where(map[*,0] ne bad_value)
       map[gp,0] = map[gp,0] * 1.e6
       map[gp,4] = map[gp,4] * 1.e12
       sz = size(map)
       print, sz

       n_bp = n_elements(bp)
       if (bp[0] ne -1) then begin
           for ipix=0l,n_bp-1 do begin
               pix2vec_nest, npix2nside(sz[1]), bp[ipix], vec0
               query_disc, npix2nside(sz[1]), vec0, 10.*sqrt(4.*!pi / sz[1]), listp, /nest
;                 help, listp
               gp = where(map[listp,0] ne bad_value)
;                 help, gp
;stop
               listp = listp[gp]

;;                print, ' '
;;                print, median(map[listp,0])
               map[bp[ipix],0] = median(map[listp,0])
               map[bp[ipix],4] = median(map[listp,4])
           endfor
       endif

       ud_grade, map, map1024,order_in='nest', order_out='ring'

       write_fits_map,'dx6_tempmap.fits',map1024[*,0],/ring
       
       ismoothing, 'dx6_tempmap.fits', 'dx6_smthmap.fits', fwhm_arcmin=eff_beams[ifreq], nlmax=min([4*ns,4000l])

       ud_grade,'dx6_smthmap.fits',out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK.fits', nside_out=ns

       mollview, out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK.fits', /hist

   endfor
   endif


   if (do_hfi) then begin

print, "-- HFI ---"
;stop
   if (do_fs) then begin
   for ifreq=ff,lf do begin
       print, sfreq[ifreq]

       infile = dir+'/hfi/HFI_'+sfreq[ifreq]+'_2048_20110411.fits'
       print, ' reading: ', infile
       read_fits_map, infile, map
       bp = where(map[*,0] eq bad_value)
       gp = where(map[*,0] ne bad_value)
       map[gp,0] = map[gp,0] * 1.e6
       map[gp,2] = map[gp,2] * 1.e12

       sz = size(map)

       n_bp = n_elements(bp)
       if (bp[0] ne -1) then begin
           for ipix=0l,n_bp-1 do begin
               pix2vec_ring, npix2nside(sz[1]), bp[ipix], vec0
               query_disc, npix2nside(sz[1]), vec0, 10.*sqrt(4.*!pi / sz[1]), listp

               gp = where(map[listp,0] ne bad_value)

               listp = listp[gp]

               map[bp[ipix],0] = median(map[listp,0])
               map[bp[ipix],2] = median(map[listp,2])
           endfor
       endif

       write_fits_map, out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK.fits', map[*,0], /ring, units='!7l!17K CMB'

      write_fits_map, out_dir+'dx6_Irms'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK.fits', sqrt(map[*,2]), /ring, units='!7l!17K CMB'

       mollview, out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK.fits',/hist
       mollview, out_dir+'dx6_Irms'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK.fits',/hist

    endfor
   endif
 
   if (do_1h) then begin
   for ifreq=ff,lf do begin
       print, sfreq[ifreq]

       infile = dir+'/hfi/HFI_'+sfreq[ifreq]+'_2048_20110412_half_1.fits'
       print, ' reading: ', infile
       read_fits_map, infile, map
       bp = where(map[*,0] eq bad_value)
       gp = where(map[*,0] ne bad_value)
       map[gp,0] = map[gp,0] * 1.e6
       map[gp,2] = map[gp,2] * 1.e12

       sz = size(map)

       n_bp = n_elements(bp)
       if (bp[0] ne -1) then begin
           for ipix=0l,n_bp-1 do begin
               pix2vec_ring, npix2nside(sz[1]), bp[ipix], vec0
               query_disc, npix2nside(sz[1]), vec0, 10.*sqrt(4.*!pi / sz[1]), listp

               gp = where(map[listp,0] ne bad_value)

               listp = listp[gp]

               map[bp[ipix],0] = median(map[listp,0])
               map[bp[ipix],2] = median(map[listp,2])
            endfor
        endif

       write_fits_map, out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_1h.fits', map[*,0], /ring, units='!7l!17K CMB'

      write_fits_map, out_dir+'dx6_Irms'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_1h.fits', sqrt(map[*,2]), /ring, units='!7l!17K CMB'

       mollview, out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_1h.fits',/hist
       mollview, out_dir+'dx6_Irms'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_1h.fits',/hist

    endfor
endif

   if (do_2h) then begin
   for ifreq=ff,lf do begin
       print, sfreq[ifreq]

       infile = dir+'/hfi/HFI_'+sfreq[ifreq]+'_2048_20110412_half_2.fits'
       print, ' reading: ', infile
       read_fits_map, infile, map
       bp = where(map[*,0] eq bad_value)
       gp = where(map[*,0] ne bad_value)
       map[gp,0] = map[gp,0] * 1.e6
       map[gp,2] = map[gp,2] * 1.e12

       sz = size(map)

       n_bp = n_elements(bp)
       if (bp[0] ne -1) then begin
           for ipix=0l,n_bp-1 do begin
               pix2vec_ring, npix2nside(sz[1]), bp[ipix], vec0
               query_disc, npix2nside(sz[1]), vec0, 10.*sqrt(4.*!pi / sz[1]), listp

               gp = where(map[listp,0] ne bad_value)

               listp = listp[gp]

               map[bp[ipix],0] = median(map[listp,0])
               map[bp[ipix],2] = median(map[listp,2])
            endfor
        endif

       write_fits_map, out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_2h.fits', map[*,0], /ring, units='!7l!17K CMB'

      write_fits_map, out_dir+'dx6_Irms'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_2h.fits', sqrt(map[*,2]), /ring, units='!7l!17K CMB'

       mollview, out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_2h.fits',/hist
       mollview, out_dir+'dx6_Irms'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_2h.fits',/hist

    endfor




   endif



   endif

; ---------------------------------------------------------------------

   if (do_5arcmin) then begin

      for ifreq=ff, lf do begin

         ismoothing, out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_2h.fits',out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_5arcmin.fits',fwhm_arcmin=sqrt(5.^2-real_beam[ifreq]^2)

         ismoothing, out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_1h.fits',out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_1h_5arcmin.fits',fwhm_arcmin=sqrt(5.^2-real_beam[ifreq]^2)

         ismoothing, out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_2h.fits',out_dir+'dx6_Imap'+string(freq[ifreq],format='(1i3.3)')+'GHz_ns'+strtrim(string(ns),2)+'_uK_2h_5arcmin.fits',fwhm_arcmin=sqrt(5.^2-real_beam[ifreq]^2)


      endfor

   endif	


   print, ' --- End of program --- '
stop

end
