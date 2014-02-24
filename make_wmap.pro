; on Magique3
; /space/ashdown/dr3/outputs/commander
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
   true  = 1b
   false = 0b

   do_map = true
   do_adjrms = false
   do_raw = true
   do_clean = false

   ns = 128l
   npix = 12l * ns^2
   beams = 60.

   freq=[30, 44, 70, 100, 143, 217, 353, 545, 857]
   sfreq = string(freq, format='(i3.3)')

   nfreq = n_elements(freq)

   ff = 1
   lf = 4
   df = 1

   hfi_TrFn_dir = '/project/projectdirs/planck/data/febecop/nominal/hfi/bls/fits/'
   lfi_TrFn_dir = '/project/projectdirs/planck/data/febecop/nominal/lfi/bls/fits/'

;   TrFn_dir = '/global/scratch/sd/paganol/run_FEBeCoP/output_ffp4/transFn/bls'
   TrFn_dir = '/global/homes/d/dpietrob/myscratch/ffp4/input/quickbeam/freqs/'

;   real_beam = [32.65, 27.00, 13.01, 9.94,   7.04,    4.66,    4.41,    4.47,    4.23]
; LFI --> http://hfilfi.planck.fr/index.php/WG/Febecop10 (left panel)
; HFI --> http://hfilfi.planck.fr/index.php/WG/Febecop8 (right panel)
   real_beam = [32.908, 27.019, 13.268, 9.544,   7.141,    4.877,    4.681,    4.653,    4.348]

   gbl = gaussbeam(beams, 4l*2048l)

   eff_beams = sqrt( float(beams)^2-real_beam^2 )

   dx = 'dx9'

   dir = '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/wmap/'

   out_dir = '/global/homes/d/dpietrob/myscratch/'+dx+'/maps/ns'+string(ns,format='(i4.4)')+'/'

   npix = 12l*ns^2
   
   bad_value = -1.63750e+30

   wsfreq = ['K','Ka','Q','V','W']
   real_beam = 60.*[0.88, 0.66, 0.51, 0.35, 0.22]
   sig0 = [1.437, 1.470, 2.197, 3.137, 6.549]
   sig0 = sig0 * 1.e3

   if (do_map) then begin
       for ifreq=ff,lf,df do begin

           if (do_raw) then begin
               infile = dir+'combined_freq/raw_maps/wmap_band_iqumap_r9_7yr_'+wsfreq[ifreq]+'_v4.fits'
               tag = 'raw'
               print, tag+' maps...'
           endif
           if (do_clean) then begin
               if ff eq 0 then return
               infile = dir+'combined_freq/clean_maps/wmap_band_forered_iqumap_r9_7yr_'+wsfreq[ifreq]+'_v4.fits'
               tag = 'clean'
               print, tag+' maps...'
           endif
           print, ' reading: ', infile
           read_fits_map, infile, map, hdr, xhdr, nside=dpns, order=dpord
;    help, map       
;stop
           bp = where(map[*,0] eq bad_value)
           gp = where(map[*,0] ne bad_value)
           map[gp,0] = map[gp,0] * 1.e3
           map[gp,3] = map[gp,3] ;* 1.e3
           
           sz = size(map)
           print, sz
               
           n_bp = n_elements(bp)
           if (bp[0] ne -1) then begin
               print, ' :DP - Missing Pixels found...'
               help, bp
               print, real_beam[ifreq], ' --> ', real_beam[ifreq]/60./!radeg
               for ipix=0l,n_bp-1 do begin
                   pix2vec_nest, npix2nside(sz[1]), bp[ipix], vec0
                   query_disc, npix2nside(sz[1]), vec0, real_beam[ifreq]/60./!radeg, listp, /nest
                   gp = where(map[listp,0] ne bad_value)
                   
                   listp = listp[gp]

                   map[bp[ipix],0] = median(map[listp,0])
                   map[bp[ipix],3] = median(map[listp,3])
               endfor
           endif
           
           if (dpord eq 'NESTED') then begin
               print, ' reordering...'
               ud_grade, map, map1024,order_in=dpord, order_out='ring'
           endif
;               mollview, map1024, chars=1.5, win=1, min=-500, max=500, tit=wsfreq[ifreq], px=500
               
           fits2cl, egbl, dir+'beams/wmap_'+wsfreq[ifreq]+'_beam_7yr.fits'

           wl = healpixwindow(dpns)

           tf = gbl/egbl        ;/wl

           cl2fits, tf, 'gbeam.fits'
               
           ismoothing, map1024[*,0], smth, /ring, fwhm_arcmin=0., beam_file='gbeam.fits', nlmax=3750, simul_type=1, /silent
           ud_grade, smth, mapud, nside_out=ns, order_in='ring'
           write_fits_map, out_dir+'wmap7_Imap_'+tag+'_'+wsfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits', mapud, /ring, units='!7l!6K'
           mollview, out_dir+'wmap7_Imap_'+tag+'_'+wsfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits', win=ifreq, px=500, min=-300, max=300
           
           ud_grade, sig0[ifreq]^2/map1024[*,3], varud, nside_out=ns, order_in='ring'
           write_fits_map, out_dir+'wmap7_Irms_'+tag+'_'+wsfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits', sqrt(varud), /ring, units='!7l!6K'
           mollview, out_dir+'wmap7_Irms_'+tag+'_'+wsfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits', px=500
       endfor
   endif

   if (do_adjrms) then begin
       for ifreq=ff,lf,df do begin

           if (do_raw) then begin
               infile = dir+'combined_freq/raw_maps/wmap_band_iqumap_r9_7yr_'+wsfreq[ifreq]+'_v4.fits'
               tag = 'raw'
               print, tag+' maps...'
           endif
           if (do_clean) then begin
               if ff eq 0 then return
               infile = dir+'combined_freq/clean_maps/wmap_band_forered_iqumap_r9_7yr_'+wsfreq[ifreq]+'_v4.fits'
               tag = 'clean'
               print, tag+' maps...'
           endif
           print, ' reading: ', infile
           read_fits_map, infile, map, hdr, xhdr, nside=dpns, order=dpord
           
           bp = where(map[*,0] eq bad_value)
           gp = where(map[*,0] ne bad_value)
           map[gp,0] = map[gp,0] * 1.e3
           map[gp,3] = map[gp,3] ;* 1.e3
           
           sz = size(map)
           print, sz
               
           n_bp = n_elements(bp)
           if (bp[0] ne -1) then begin
               print, ' :DP - Missing Pixels found...'
               help, bp
               print, real_beam[ifreq], ' --> ', real_beam[ifreq]/60./!radeg
               for ipix=0l,n_bp-1 do begin
                   pix2vec_nest, npix2nside(sz[1]), bp[ipix], vec0
                   query_disc, npix2nside(sz[1]), vec0, real_beam[ifreq]/60./!radeg, listp, /nest
                   gp = where(map[listp,0] ne bad_value)
                   
                   listp = listp[gp]

                   map[bp[ipix],0] = median(map[listp,0])
                   map[bp[ipix],3] = median(map[listp,3])
               endfor
           endif
           
           if (dpord eq 'NESTED') then begin
               print, ' reordering...'
               ud_grade, map, map1024,order_in=dpord, order_out='ring'
           endif
;               mollview, map1024, chars=1.5, win=1, min=-500, max=500, tit=wsfreq[ifreq], px=500
               
           fits2cl, egbl, dir+'beams/wmap_'+wsfreq[ifreq]+'_beam_7yr.fits'

           wl = healpixwindow(dpns)

           tf = gbl/egbl        ;/wl                                                                                                                   

           cl2fits, tf, 'gbeam.fits'

           ave = dblarr(npix)
           ave2 = dblarr(npix)

           nsims = 10

           for isim=0,nsims-1 do begin
               print, isim
               noise = randomn(-(1+isim),12l*dpns^2) * sig0[ifreq] / sqrt(map1024[*,3])
               ismoothing, noise, smth, /ring, fwhm_arcmin=0., beam_file='gbeam.fits', simul_type=1, /silent

               ud_grade, smth, mapud, nside_out=ns, order_in='ring'
;               mollview, mapud, px=500
               ave = ave + mapud
               ave2 = ave2 + mapud^2
           endfor

           ave = ave / nsims
           ave2 = ave2 / nsims

           std = sqrt(ave2 - ave^2)
           gp = where(finite(std,/nan) eq 0.)
           help, gp
           read_fits_map, out_dir+'wmap7_Irms_'+tag+'_'+wsfreq[ifreq]+'_ns'+strtrim(string(ns),2)+'_uK.fits', rms
           scaling = mean( std[gp]/rms[gp] )
           print, wsfreq[ifreq], scaling

           write_fits_map, out_dir+'wmap7_adjrms_'+tag+'_'+wsfreq[ifreq]+'.fits', rms*scaling, /ring, units='!7l!6K'
           mollview, out_dir+'wmap7_adjrms_'+tag+'_'+wsfreq[ifreq]+'.fits', px=500
;stop
       endfor
   endif
   
   print, ' removing *temp*'
   spawn, 'remove *temp* gbeam.fits'
   
   print, ' --- End of program --- '
   stop
   
end
