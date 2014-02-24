; Only Delta_dx9 is used for the first delivery. 353 is still the old one.
; Most updated and probably correct one

   !path = !path + ':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'

   True = 1b
   False = 0b

   do_moll = False

   do_lfi   = False
   do_hfi   = True

   do_berkeley      = False
   do_Delta         = True;False;True
   do_MJyResca_only = True;False;True
   do_new_zodi      = True
   do_mirror        = False

   do_var   = False

   if (do_berkeley or do_new_zodi) then do_var = False

   tag = 'Delta_noZodi' ;'Delta_new_zodi' ;'berkeley';'Delta_MJyResca_only';'' ;'Delta'

   fstf = 7 
   lstf = 8

;##   lmax  = 3500l
   ns    = 128l
   beams = 60.

   freq=[30, 44, 70, 100, 143, 217, 353, 545, 857]
   sfreq = string(freq, format='(i3.3)')

   cv_MJy2Kcmb = 1.d0 / [1.,1.,1.,1.,1.,1.,1., 57.976641, 2.2690745] ; only for 545 and 857 channels

   nfreq = n_elements(freq)

   dx_version = 'dx9'

   beam_dir = '/global/scratch/sd/paganol/run_FEBeCoP/output_dx9/transFn/Bl/'

   real_beam = [32.65, 27.00, 13.01, 9.94,   7.04,    4.66,    4.41,    4.47,    4.23]

   gbl = gaussbeam(beams, 4l*2048l)

   eff_beams = sqrt( float(beams)^2-real_beam^2 )

   dir = '/project/projectdirs/planck/data/mission/DPC_maps/dx9/'
   if (do_mirror) then dir = '/global/homes/d/dpietrob/myscratch/mirror/'

   out_dir = '/global/scratch/sd/dpietrob/'+dx_version+'/delta_maps/'+dx_version+'/'
   spawn, 'mkdir -p -m g+xr /global/scratch/sd/dpietrob/'+dx_version+'/delta_maps'
   spawn, 'mkdir -p -m g+xr /global/scratch/sd/dpietrob/'+dx_version+'/delta_maps/'+dx_version

   npix = 12l*ns^2
   
   bad_value = -1.63750e+30

   if (do_lfi) then begin

       if (not do_Delta) then spawn, 'ls '+dir+'lfi/LFI*9_nominal.fits ', files else $
         files = '/project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/lfi/'+ $
         ['LFI_30_1024_20120914_nominal_1s.fits', $
          'LFI_44_1024_20120914_nominal_1s.fits', $
          'LFI_70_1024_20120912_nominal_1s.fits']

       if (do_berkeley) then files = '/project/projectdirs/planck/data/mission/dx9/us_maps/map_dx9rc1_'+sfreq[0:2]+'_nom_01s.fits'
       help, files
       print, files

       if (do_mirror) then spawn, 'ls '+ dir + 'lfi_delta/*.fits', files
       help, files
       print, files
;stop

       for ifreq=fstf,min([lstf,2]) do begin
           print, sfreq[ifreq]

           infile = files[ifreq]
           print, ' reading: ', infile
           read_fits_map, infile, map, hdr, xhdr, nside=dpns, order=dpord

           print, ' reordering...'
           map = reorder(map, in=dpord, out='RING')

           bp = where(map[*,0] eq bad_value)
           gp = where(map[*,0] ne bad_value)
           map[gp,0] = map[gp,0] * 1.e6
           if (do_var) then map[gp,4] = map[gp,4] * 1.e12
           sz = size(map)
           print, sz

           n_bp = n_elements(bp)
           if (bp[0] ne -1) then begin
               print, ' :DP - Missing Pixels found...'
               help, bp
               for ipix=0l,n_bp-1 do begin
                   pix2vec_ring, npix2nside(sz[1]), bp[ipix], vec0
                   query_disc, npix2nside(sz[1]), vec0, 2.*real_beam[ifreq]/60./!radeg, listp, /nest

                   gp = where(map[listp,0] ne bad_value)

                   listp = listp[gp]

                   map[bp[ipix],0] = median(map[listp,0])
                   if (do_var) then map[bp[ipix],4] = median(map[listp,4])
               endfor
           endif

; ------ Not stored because we are going to use Belen's ones

           map1024 = map[*,0]
           if (do_moll) then mollview, map1024, win=1, /hist, tit=files[ifreq], px=500
           if (do_var and do_moll) then mollview, map[*,4], win=2, tit=files[ifreq]+' RMS', px=500
           
           print, ' harmonic space downgrading:'
           print, ' anafasting...'

           fits2cl, egbl, beam_dir+'Bl_GB_cent_'+sfreq[ifreq]+'_dx9_CMB.fits'
;print, egbl[0:10]
           wl = healpixwindow(dpns)
           hwl = healpixwindow(ns)
;print, wl[0:10]
           tf = hwl[*,0]*0.
           tf[0:1] = 1.d0
           tf[2:n_elements(tf)-1] = gbl[2:n_elements(tf)-1]/(egbl[2:n_elements(tf)-1]*wl[2:n_elements(tf)-1]) ;* hwl[2:*]
;print, tf[0:10]
           cl2fits, tf, 'tmpbeam.fits'
           help, egbl, gbl, wl, tf
           if (do_moll) then begin
               window, 0 & plot_io, tf, chars=2
               oplot, egbl
               oplot, gbl, col=245
           endif
;stop
           ismoothing, map1024, smth, /ring, fwhm_arcmin=0.d0, beam_file='tmpbeam.fits', /silent
           ud_grade, smth, mapud, nside_out=ns, order_in='ring'
;##           mapfile = out_dir+dx_version+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
;##           if (do_Delta) then mapfile = out_dir+dx_version+'_'+tag+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
           mapfile = out_dir+dx_version+'_'+tag+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
           write_fits_map, mapfile, mapud, /ring, units='!7l!6K CMB'
                  
           if (do_moll) then mollview, mapfile, min=-300, max=300, px=500, win=ifreq
;##           ianafast, out_dir+dx_version+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits', tt, /show_cl, iter=2, /silent
;stop
           if (do_var) then begin
               var = map[*,4]
               ismoothing, sqrt(var), smth, /ring, fwhm_arcmin=0.d0, beam_file='tmpbeam.fits', /silent
               ud_grade, smth, varud, nside_out=ns, order_in='ring'
;##               varfile = out_dir+dx_version+'_Irms_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
;##               if (do_Delta) then varfile = out_dir+dx_version+'_'+tag+'_Irms_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
               varfile = out_dir+dx_version+'_'+tag+'_Irms_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
               write_fits_map, varfile, varud, /ring, units='!7l!6K CMB'
               if (do_moll) then mollview, varfile, px=500, win=ifreq
           endif
       endfor
   endif

;stop
   if (do_hfi) then begin

       print, "-- HFI ---"

       for ifreq = max([fstf,3]),lstf do begin
           print, sfreq[ifreq]
;##       infile = files[ifreq]
           if (not do_delta) then begin
               if (ifreq lt 8) then spawn, 'ls '+dir+'hfi/official/HFI_'+sfreq[ifreq]+'_2048_20120611_nominal.fits', infile else $
                 spawn, 'ls '+dir+'hfi/official/HFI_'+sfreq[ifreq]+'_2048_20120612_nominal.fits', infile
           endif else begin
               print, 'HFI delta'
               if (do_new_zodi) then begin
                   print, 'new_zodi'
                   spawn, 'ls /global/scratch/sd/loris/OD/DX/FromMagique/dx9/nominal/hfi_delta_nozodi/hfi_'+sfreq[ifreq]+'_nominal_nozodi.fits', infile
               endif else begin
                   if (ifreq lt 6) then spawn, 'ls '+dir+'hfi/official/HFI_'+sfreq[ifreq]+'_2048_20120611_nominal.fits', infile
                   if (ifreq eq 6) then spawn, 'ls /project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/hfi/noZodi/FREQ/HFI_'+sfreq[ifreq]+'_2048_20121128_nominal_noZodi.fits', infile
                   if ( (ifreq gt 6) and (not do_MJyResca_only) ) then spawn, 'ls /project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/hfi/MJyResca_noZodi/FREQ/HFI_'+sfreq[ifreq]+'_2048_20121129_nominal_MJyResca_noZodi.fits', infile
                   if ( (ifreq gt 6) and (do_MJyResca_only) ) then spawn, 'ls /project/projectdirs/planck/data/mission/DPC_maps/dx9_delta/hfi/MJyResca/FREQ/HFI_'+sfreq[ifreq]+'_2048_20121128_nominal_MJyResca.fits', infile
               endelse
           endelse
               
           if (do_berkeley) then infile = '/project/projectdirs/planck/data/mission/dx9/us_maps/map_dx9rc1_'+sfreq[ifreq]+'_nom_60s.fits'

;##           if (do_mirror) then infile =  dir + 'hfi_delta/HFI_'+sfreq[ifreq]+'_2048_20121208_nominal.fits'
           if (do_mirror) then infile =  dir + 'hfi_delta_withzodi/hfi_'+sfreq[ifreq]+'_nominal_withzodi.fits'

           print, ' reading: ', infile
           read_fits_map, infile, map, nside=dpns, order=dpord
           print, dpord
           map = reorder(map, in=dpord, out='RING')

           bp = where(map[*,0] eq bad_value)
           gp = where(map[*,0] ne bad_value)
           map[gp,0] = map[gp,0] * 1.e6
           if (ifreq gt 6) then begin
               if (do_delta or do_mirror) then map[gp,0] = map[gp,0] * cv_MJy2Kcmb[ifreq]
           endif
;##           map[gp,2] = map[gp,2] * 1.e12
           if (do_var) then begin
               if (ifreq lt 7) then  map[gp,4] = map[gp,4] * 1.e12 else  map[gp,2] = map[gp,2] * 1.e12

               if (ifreq gt 6) then begin
                   if (do_delta or do_mirror) then map[gp,2] = map[gp,2] * cv_MJy2Kcmb[ifreq]^2
               endif
;       map[bp,*] = 0.
           endif
           sz = size(map)

           n_bp = n_elements(bp)
           if (bp[0] ne -1) then begin
               print, ' :DP - Missing Pixels Found:'
               help, bp
               for ipix=0l,n_bp-1 do begin
                   pix2vec_ring, npix2nside(sz[1]), bp[ipix], vec0
                   query_disc, npix2nside(sz[1]), vec0, 2.*real_beam[ifreq]/60./!radeg, listp
;               help, listp
                   gp = where(map[listp,0] ne bad_value)
;               help, gp
;stop
                   listp = listp[gp]
;               help, listp
;;                print,' '
;;                print, median(map[listp,0])
                   map[bp[ipix],0] = median(map[listp,0])
                   if (do_var) then begin
                       if (ifreq lt 7) then map[bp[ipix],4] = median(map[listp,4]) else  map[bp[ipix],2] = median(map[listp,2])
                   endif
               endfor
           endif

           map1024 = map[*,0]
           if (do_moll) then mollview, map1024, /hist, tit=infile, win=1, px=500
           if (do_var) then begin
               if (ifreq lt 7) then var = map[*,4] else var = map[*,2]
               if (do_moll) then mollview, var, tit=infile+' RMS', win=2, px=500
           endif

           fits2cl, egbl, beam_dir+'Bl_BS_Mars12_cent_'+sfreq[ifreq]+'_dx9_CMB.fits'
;print, egbl[0:10]
;                                                                                 
           wl = healpixwindow(dpns)
           hwl = healpixwindow(ns)
;print, wl[0:10]
;                                                                                 
           tf = hwl[*,0] * 0.
           tf[0:1] = 1.d0
           tf[2:n_elements(tf)-1] = gbl[2:n_elements(tf)-1]/(egbl[2:n_elements(tf)-1]*wl[2:n_elements(tf)-1]) ;* hwl[2:*]
;print, tf[0:10]
;                                                                                 
           cl2fits, tf, 'tmpbeam.fits'
           help, egbl, gbl, wl, tf
           if (do_moll) then window, 0 & plot_io, tf
;stop
;                                                                                 
           ismoothing, map1024, mapsmth, /ring, fwhm_arcmin=0.d0, beam_file='tmpbeam.fits', /silent
           ud_grade, mapsmth, mapud, nside_out=ns,order_in='ring'

;##           mapfile = out_dir+dx_version+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
;##           if (do_delta) then mapfile = out_dir+dx_version+'_'+tag+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
           mapfile = out_dir+dx_version+'_'+tag+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
           write_fits_map, mapfile, mapud, /ring
           if (do_moll) then mollview, mapfile, min=-300, max=300, px=500, win=ifreq
;           ianafast, out_dir+dx_version+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits', tt, /show_cl, /silent
;stop
           if (do_var) then begin
               ismoothing, sqrt(var), smth, /ring, fwhm_arcmin=0.d0, beam_file='tmpbeam.fits', /silent
               ud_grade, smth, varud, nside_out=ns, order_in='ring'
; ---
;##               mollview, sqrt(var), win=10, px=500, tit='sqrt(var)'
;##               mollview, smth, win=11, px=500, tit='smth'
;##               mollview, varud, win=12, px=500, tit='ud_grade smth'
; ---
;##               varfile = out_dir+dx_version+'_Irms_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
;##               if (do_delta) then varfile = out_dir+dx_version+'_'+tag+'_Irms_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
               varfile = out_dir+dx_version+'_'+tag+'_Irms_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
               write_fits_map, varfile, varud, /ring, units='!7l!6K CMB'
               if (do_moll) then mollview, varfile, px=500, win=ifreq
           endif
       endfor
       
   endif

   print, ' removing *tmp*'
   spawn, 'rm *tmp*'
   spawn, 'chgrp planck '+out_dir+dx_version+'*_Imap_*_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
   spawn, 'chmod g+r '+out_dir+dx_version+'*_Imap_*_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
   

   print, ' --- End of program --- '
stop

end
