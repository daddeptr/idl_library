; Generalized to FFP
; Generalized to generic DX adding lfi/hfi_dx_tag
; 
; Polarization added in maps only (no variance)
; DX10
;
; Only Delta_dx9 is used for the first delivery. 353 is still the old one.
; Most updated and probably correct one

;##   !path = !path + ':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'

   True = 1b
   False = 0b

   do_moll = False

;   do_lfi  = False
;   do_hfi  = True

   do_var  = False
   do_pol  = True

   fstf = 3 
   lstf = 8

;##   lmax  = 3500l
   ns    = -1l
   if (ns lt 0) then do_smooth = False
   print, 'do_smooth = ', do_smooth

   beams = 60.

   freq=[30, 44, 70, 100, 143, 217, 353, 545, 857]
   sfreq = string(freq, format='(i3.3)')

   cv_MJy2Kcmb = 1.d0 / [1.,1.,1.,1.,1.,1.,1., 57.976641, 2.2690745] ; only for 545 and 857 channels

   nfreq = n_elements(freq)

   dx_version = 'ffp7'
   tag = '';'nobpm_';'cmb_'
   lfi_dx_tag = 'ffp7_total_bpm' ;'ffp7_total_nobpm' ;'ffp7_cmb';
;##   hfi_dx_tag = 'ffp7_cmb';'ffp7_total_bpm'
   mission = 'full_map'
   dir = '/project/projectdirs/planck/data/'+dx_version+'/'
   hfi_dir = dir+'maps/total/';cmb/';
   lfi_dir = dir+'maps/total/';cmb/';

   beam_dir = '/global/scratch/sd/paganol/run_FEBeCoP/output_dx9/transFn/Bl/'
   real_beam = [32.65, 27.00, 13.01, 9.94,   7.04,    4.66,    4.41,    4.47,    4.23]
   gbl = gaussbeam(beams, 4l*2048l)
   eff_beams = sqrt( float(beams)^2-real_beam^2 )

   spawn, 'mkdir -p -m g+xr /global/scratch2/sd/dpietrob/'+dx_version
   spawn, 'mkdir -p -m g+xr /global/scratch2/sd/dpietrob/'+dx_version+'/maps'
;##   spawn, 'mkdir -p -m g+xr /global/scratch2/sd/dpietrob/'+dx_version+'/maps/ns2048'
;##   out_dir = '/global/scratch2/sd/dpietrob/'+dx_version+'/maps/ns2048/'
   out_dir = '/global/scratch2/sd/dpietrob/'+dx_version+'/maps/'

   npix = 12l*ns^2
   
   bad_value = -1.63750e+30

   for ifreq=fstf,lstf do begin
       print, sfreq[ifreq]
       if ifreq gt 6 then do_pol = False

       fstr = lfi_dir+lfi_dx_tag+'_'+sfreq[ifreq]+'_'+mission+'*.fits'
       print, fstr
       spawn, 'ls '+fstr, files
       help, files
       print, files

       print, n_elements(files)
       for ifile=0,n_elements(files)-1 do begin
           infile = files[ifile]
           print, ' reading: ', infile

           tags = strsplit(infile,'/',/extract)
           tags = tags[n_elements(tags)-1]
           tags = strsplit(tags,'.',/extract)
;           print, tags
           tags = tags[0]
           imission = strpos(tags,mission)
           mission_tag = strmid(tags,imission)
           print, mission_tag
;stop
           read_fits_map, infile, map, hdr, xhdr, nside=dpns, order=dpord

           print, ' reordering...'
           if ( (dpord eq 'RING') or (dpord eq 'ring') ) then map = reorder(map, in=dpord, out='NESTED')

           bp = where(map[*,0] eq bad_value)
           gp = where(map[*,0] ne bad_value)
           map[gp,0] = map[gp,0] * 1.e6
           if (do_var) then map[gp,4] = map[gp,4] * 1.e12
           sz = size(map)
           print, sz

           n_bp = n_elements(bp)
           if (bp[0] ne -1) then begin
               print, ' :DP - Missing Pixels found in T...'
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

           if (do_pol) then begin
               for ip=1,2 do begin
                   bp = where(map[*,ip] eq bad_value)
                   gp = where(map[*,ip] ne bad_value)
                   map[gp,ip] = map[gp,ip] * 1.e6
                   
                   n_bp = n_elements(bp)
                   if (bp[0] ne -1) then begin
                       print, ' :DP - Missing Pixels found in P...'
                       help, bp
                       for ipix=0l,n_bp-1 do begin
                           pix2vec_ring, npix2nside(sz[1]), bp[ipix], vec0
                           query_disc, npix2nside(sz[1]), vec0, 2.*real_beam[ifreq]/60./!radeg, listp, /nest
                           gp = where(map[listp,ip] ne bad_value)
                           listp = listp[gp]
                           map[bp[ipix],ip] = median(map[listp,ip])
                       endfor
                   endif
               endfor
           endif
; ------ Not stored because we are going to use Belen's ones

           if (not do_pol) then map1024 = map[*,0] else map1024 = map[*,0:2]
           
           if (do_moll) then mollview, map1024, win=1, /asinh, tit=files[ifile], px=500, /nest
           if (do_var and do_moll) then mollview, map[*,4], win=2, tit=files[ifile]+' RMS', px=500, /nest
           
           if do_smooth then begin
               if do_pol then stop, 'Polarization not implemented'
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
               ismoothing, map1024, smth, /nest, fwhm_arcmin=0.d0, beam_file='tmpbeam.fits', /silent
               ud_grade, smth, mapud, nside_out=ns, order_in='nest', order_out='ring'
;##           mapfile = out_dir+dx_version+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
;##           if (do_Delta) then mapfile = out_dir+dx_version+'_'+tag+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
               mapfile = out_dir+dx_version+'_Imap_'+tag+'map_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
               write_fits_map, mapfile, mapud, /ring, units='!7l!6K CMB'
                   
               if (do_moll) then mollview, mapfile, min=-300, max=300, px=500, win=ifreq
;##           ianafast, out_dir+dx_version+'_Imap_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits', tt, /show_cl, iter=2, /silent
;stop
               if (do_var) then begin
                   var = map[*,4]
                   ismoothing, sqrt(var), smth, /nest, fwhm_arcmin=0.d0, beam_file='tmpbeam.fits', /silent
                   ud_grade, smth, varud, nside_out=ns, order_in='nest', order_out='ring'
;##               varfile = out_dir+dx_version+'_Irms_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
;##               if (do_Delta) then varfile = out_dir+dx_version+'_'+tag+'_Irms_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
                   varfile = out_dir+dx_version+tag+'_Irms_'+sfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK_'+string(fix(beams),format='(i2.2)')+'a.fits'
                   write_fits_map, varfile, varud, /ring, units='!7l!6K CMB'
                   if (do_moll) then mollview, varfile, px=500, win=ifreq
               endif
               
           endif else begin
               print, 'No processing...'
               mapfile = out_dir+dx_version+'_IQU'+tag+'map_'+sfreq[ifreq]+'_'+mission_tag+'_uK.fits'
               if (not do_pol) then write_fits_map, mapfile, map1024, /nest, units='!7l!6K CMB' else write_tqu, mapfile, map1024, /nest, units='!7l!6K CMB'
;##               if (do_moll) then mollview, mapfile, min=-300, max=300, px=500, win=ifreq
               if (do_moll) then mollview, mapfile, /asinh, px=500, win=ifreq
               
               if (do_var) then begin
                   varfile = out_dir+dx_version+tag+'_Irms_'+sfreq[ifreq]+'_'+mission_tag+'_uK.fits'
                   write_fits_map, varfile, map[*,4], /nest, units='!7l!6K CMB'
               endif
           endelse
;               stop
       endfor
   endfor

   print, ' removing *tmp*'
   spawn, 'rm *tmp*'
   spawn, 'chgrp planck '+ out_dir+dx_version+'_IQU'+tag+'map_*_'+mission_tag+'_uK.fits'
   spawn, 'chmod g+r '+ out_dir+dx_version+'_IQU'+tag+'map_*_'+mission_tag+'_uK.fits'
   

   print, ' --- End of program --- '
stop

end
