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
   do_hfi = 0b
   do_smoothing = 0b
   do_plotting = 1b

   ns = 2048l
;;    beams = 90.
;;    gsig = 3.

   freq=[30, 44, 70, 100, 143, 217, 353]
   sfreq=['30', '44', '70', '100', '143', '217', '353']
;   sfreq=string(freq, format='(1i3.3)')
   nfreq = n_elements(freq)

   real_beam = [32.65, 27.00, 13.01, 9.94,   7.04,    4.66,    4.41,    4.47,    4.23]


   dir = '/project/projectdirs/planck/data/mission/DPC_maps/dx6/'

   out_dir = './'

   npix = 12l*ns^2
   
   bad_value = -1.63750e+30

   read_fits_map, 'chains/dx6_cmb_5arcmin.fits', cmb

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

;;        write_fits_map,'dx6_tempmap.fits',map1024[*,0],/ring
       
       if (real_beam[ifreq] gt 5.) then begin
           ismoothing, cmb, cmbsm, /ring, fwhm_arcmin=sqrt(real_beam[ifreq]^2-5.^2), nlmax=min([4*ns,4000l])
       endif else begin
           cmbudg = cmb
       endelse

       ud_grade, cmbsm, cmbudg, nside_out=npix2nside(sz[1]), order_in='ring', order_out='ring'	
       fg = map1024-cmbudg
;       mollview, fg, /hist, chars=1.5

       write_fits_map, out_dir+'fg_'+sfreq[ifreq]+'.fits',fg, /ring

       mollview, out_dir+'fg_'+sfreq[ifreq]+'.fits', /hist

   endfor
   endif


   if (do_hfi) then begin

print, "-- HFI ---"
;stop
ff = 4
   for ifreq=ff,nfreq-1 do begin
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

       map1024 = map[*,0]

       if (real_beam[ifreq] gt 5.) then begin
           ismoothing, cmb, cmbsm, /ring, fwhm_arcmin=sqrt(real_beam[ifreq]^2-5.^2), nlmax=min([4*ns,4000l])
       endif else begin
           cmbudg = cmb
       endelse

       ud_grade, cmbsm, cmbudg, nside_out=npix2nside(sz[1]), order_in='ring', order_out='ring'
       fg = map1024-cmbudg
;       mollview, fg, /hist, chars=1.5

       write_fits_map, out_dir+'fg_'+sfreq[ifreq]+'.fits',fg, /ring

       mollview, out_dir+'fg_'+sfreq[ifreq]+'.fits', /hist


    endfor

 
   endif	

   if (do_smoothing) then  begin
       for ifreq=0,nfreq-1 do begin
           ismoothing, 'fg_'+sfreq[ifreq]+'.fits','smth.fits', fwhm_arcmin=sqrt(60.^2-real_beam[ifreq]^2)
           ud_grade, 'smth.fits','fg_'+sfreq[ifreq]+'_60arcmin_ns512.fits',nside_out=512, order_out='ring'
       endfor
   endif

   if (do_plotting) then begin
       fig_dir = '../ns128/png/'
       for ifreq=0,nfreq-1 do begin
           mollview, 'fg_'+sfreq[ifreq]+'.fits', min=0, chars=1.5, /hist, tit='!17Foregrounds @ '+sfreq[ifreq]+' GHz', units='!7l!8K CMB', png=fig_dir+'fg_'+sfreq[ifreq]+'.png',window=-1
           mollview, 'fg_'+sfreq[ifreq]+'_60arcmin_ns512.fits', min=0, chars=1.5, /hist, tit='!17Foregrounds @ '+sfreq[ifreq]+' GHz - 1 deg', units='!7l!8K CMB', png=fig_dir+'fg_'+sfreq[ifreq]+'_1deg.png', window=-1
       endfor
   endif

   print, ' --- End of program --- '
stop

end
