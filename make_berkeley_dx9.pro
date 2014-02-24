; Most updated and probably correct one

   !path = !path + ':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'

   true = 1b
   false = 0b

   do_lfi = True
   do_Delta = True
   do_hfi = True

   fstf = 2
   lstf = 8

   lmax  = 3750l
   ns    = 1024l
   beams = 15.

   freq=[30, 44, 70, 100, 143, 217, 353, 545, 857]
   sfreq = string(freq, format='(i3.3)')

   nfreq = n_elements(freq)

   dx_version = 'dx9'

;##   hfi_TrFn_dir = '/project/projectdirs/planck/data/febecop/nominal/hfi/bls/fits/'
;##   lfi_TrFn_dir = '/project/projectdirs/planck/data/febecop/nominal/lfi/bls/fits/'

   beam_dir = '/global/scratch/sd/paganol/run_FEBeCoP/output_dx9/transFn/Bl/'

   real_beam = [32.65, 27.00, 13.01, 9.94,   7.04,    4.66,    4.41,    4.47,    4.23]

   gbl = gaussbeam(beams, 4l*2048l)

   eff_beams = sqrt( float(beams)^2-real_beam^2 )

;;    dir = '/project/projectdirs/planck/user/planck/dx4/'
;;   dir = '/project/projectdirs/planck/data/mission/DPC_maps/dx7/'
;##   dir = '/project/projectdirs/planck/data/mission/DPC_maps/dx8/'
   dir = '/project/projectdirs/planck/data/mission/DPC_maps/dx9/'

   out_dir = '/global/scratch/sd/dpietrob/'+dx_version+'/maps/ns'+string(ns, format='(i4.4)')+'/'
   spawn, 'mkdir -p /global/scratch/sd/dpietrob/'+dx_version+'/maps'
   spawn, 'mkdir -p /global/scratch/sd/dpietrob/'+dx_version+'/maps/ns'+string(ns, format='(i4.4)')+'/'

   npix = 12l*ns^2
   
   bad_value = -1.63750e+30

   if (do_lfi) then begin

       spawn, 'ls '+dir+'lfi/LFI*9_nominal.fits ', files
       if (do_Delta) then files = ['/global/scratch/sd/planck/user/zonca/DX9_Delta/LFI_30_1024_20120914_nominal_1s.fits', $
                                  '/global/scratch/sd/planck/user/zonca/DX9_Delta/LFI_44_1024_20120914_nominal_1s.fits', $
                                  '/global/scratch/sd/planck/user/zonca/DX9_Delta/LFI_70_1024_20120912_nominal_1s.fits']
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
           map[gp,4] = map[gp,4] * 1.e12
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
                   map[bp[ipix],4] = median(map[listp,4])
               endfor
           endif

; ------ Not stored because we are going to use Belen's ones
;##           write_fits_map, out_dir+dx_version+'_variance_'+sfreq[ifreq]+'_nominal_ns2048_uK.fits', map[*,4], /ring

           map1024 = map[*,0]
           mollview, map1024, chars=1.5, win=1, /hist, tit=files[ifreq], px=500
           
           var = map[*,4]
           ud_grade, var, varud, nside_out=ns, order_in='ring'
           varfile = out_dir+dx_version+'_Irms_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits'
           if (do_Delta) then varfile = out_dir+dx_version+'_Delta_Irms_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits'
           write_fits_map, varfile, sqrt(varud), /ring, units='!7l!6K CMB'
           mollview, varfile, px=500, win=ifreq

           print, ' harmonic space downgrading:'
           print, ' anafasting...'

           fits2cl, egbl, beam_dir+'Bl_GB_cent_'+sfreq[ifreq]+'_dx9_CMB.fits'
;print, egbl[0:10]
           wl = healpixwindow(dpns)
           hwl = healpixwindow(ns)
;print, wl[0:10]
           tf = hwl[*,0]*0.
           tf[0:1] = 1.d0
           tf[2:n_elements(tf)-1] = gbl[2:n_elements(tf)-1]/(egbl[2:n_elements(tf)-1]*wl[2:n_elements(tf)-1]) * hwl[2:*]
;print, tf[0:10]
           cl2fits, tf, 'tmpbeam.fits'
           help, egbl, gbl, wl, tf
           window, 0 & plot_io, tf
;stop
           ismoothing, map1024, smth, /ring, fwhm_arcmin=0.d0, beam_file='tmpbeam.fits', /silent
           ud_grade, smth, mapud, nside_out=ns, order_in='ring'
           mapfile = out_dir+dx_version+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_h.fits'
           if (do_Delta) then mapfile = out_dir+dx_version+'_Delta_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_h.fits'
           write_fits_map, mapfile, mapud, /ring, units='!7l!6K CMB'
                  
           mollview, mapfile, chars=1.5, min=-500, max=500, px=500, win=ifreq
;##           ianafast, out_dir+dx_version+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits', tt, /show_cl, iter=2, /silent
;stop
       endfor
   endif

;stop
   if (do_hfi) then begin

       print, "-- HFI ---"

       for ifreq = max([fstf,3]),lstf do begin
           print, sfreq[ifreq]
;##       infile = files[ifreq]
           if (ifreq lt 8) then spawn, 'ls '+dir+'hfi/official/HFI_'+sfreq[ifreq]+'_2048_20120611_nominal.fits', infile else $
             spawn, 'ls '+dir+'hfi/official/HFI_'+sfreq[ifreq]+'_2048_20120612_nominal.fits', infile

           print, ' reading: ', infile
           read_fits_map, infile, map, nside=dpns, order=dpord
           print, dpord
           map = reorder(map, in=dpord, out='RING')

           bp = where(map[*,0] eq bad_value)
           gp = where(map[*,0] ne bad_value)
           map[gp,0] = map[gp,0] * 1.e6
;##           map[gp,2] = map[gp,2] * 1.e12
           if (ifreq lt 7) then  map[gp,4] = map[gp,4] * 1.e12 else  map[gp,2] = map[gp,2] * 1.e12
;       map[bp,*] = 0.

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
                   if (ifreq lt 7) then map[bp[ipix],4] = median(map[listp,4]) else  map[bp[ipix],2] = median(map[listp,2])
               endfor
           endif

           if (ifreq lt 7) then var = map[*,4] else var = map[*,2]
           ud_grade, var, varud, nside_out=ns, order_in='ring'
           varfile = out_dir+dx_version+'_Irms_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits'
;           if (do_Delta) then mapfile = out_dir+dx_version+'_Delta_Irms_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits'
           write_fits_map, varfile, sqrt(varud), /ring, units='!7l!6K CMB'
           mollview, varfile, px=500, win=ifreq

           map1024 = map[*,0]
           mollview, map1024, /hist, tit=infile, win=2, px=500

           fits2cl, egbl, beam_dir+'Bl_BS_Mars12_cent_'+sfreq[ifreq]+'_dx9_CMB.fits'
;print, egbl[0:10]
;                                                                                 
           wl = healpixwindow(dpns)
           hwl = healpixwindow(ns)
;print, wl[0:10]
;                                                                                 
           tf = hwl[*,0] * 0.
           tf[0:1] = 1.d0
           tf[2:n_elements(tf)-1] = gbl[2:n_elements(tf)-1]/(egbl[2:n_elements(tf)-1]*wl[2:n_elements(tf)-1]) * hwl[2:*]
;print, tf[0:10]
;                                                                                 
           cl2fits, tf, 'tmpbeam.fits'
           help, egbl, gbl, wl, tf
           window, 0 & plot_io, tf
;stop
;                                                                                 
           ismoothing, map1024, mapsmth, /ring, fwhm_arcmin=0.d0, beam_file='tmpbeam.fits', /silent

           ud_grade, mapsmth,mapud, nside_out=ns,order_in='ring'

           write_fits_map, out_dir+dx_version+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_h.fits', mapud, /ring
           mollview, out_dir+dx_version+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_h.fits', min=-500, max=500, px=500, win=ifreq
;           ianafast, out_dir+dx_version+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits', tt, /show_cl, /silent
;stop
       endfor
       
   endif

   print, ' removing *tmp*'
   spawn, 'remove *tmp*'

   print, ' --- End of program --- '
stop

end
