; changed 'all' --> 'e2e2'

; readcol,'/project/projectdirs/planck/data/ffp3/mission/fpdb_ffp3.txt',ch,nu,fwhm,format='a,x,x,x,x,x,f,x,x,x,x,x,x,x,x,x,f'
; l70=where(fwhm gt 12 and fwhm lt 14)
; print, mean(fwhm[l70])

     true = 1b
     false = 0b

     freq = [30, 44, 70, 100, 143, 217, 353]
     convf = conversionfactor(freq,/antenna2thermo)
     print, freq
     print, convf

     sfreq = string(freq, format='(i3.3)')

     nfreq = n_elements(freq)

     smoo = 15.
     strsmoo = strtrim( string(long(smoo)),2 )

     beam = [33., 23., 12.7143, 9.64461, 7.11419, 4.72230, 4.5648]

     bad_value = -1.63750e+30

     ns = [512l, 1024l, 1024l, 2048l, 2048l, 2048l, 2048l ]
     npix = 12 * ns^2
     sns = strtrim(string(ns),2)

     print,  'Down smoothing: ', smoo

     if (0b) then begin
        plot_io, /nodata, [1,4000], [1.e-8,1.e4], chars=1.5
;;      for ifreq=2, nfreq-1 do begin
        for ifreq=4, nfreq-1 do begin
           print, sfreq[ifreq]

           infile = '/project/projectdirs/planck/data/ffp3/outputs/maps/ffp3.springtide'+sns[ifreq]+'.'+sfreq[ifreq]+'.e2e2.fits'
           print, 'Reading '+infile
           print, freq[ifreq]
           print, ns[ifreq]
           print, npix[ifreq]
           print, beam[ifreq]

           read_fits_map,infile,map, hdr, xhdr

           sz = size(map)
           print, sz

           bp = where(map[*,0] eq bad_value)
           gp = where(map[*,0] ne bad_value)
           help, bp
;           if (total(bp) gt 0) then stop

           temp = fltarr(12l*ns[ifreq]^2)
           var = temp

           temp[gp] = map[gp,0] * 1.e6  * convf[ifreq]
           var[gp]  = map[gp,3] * 1.e12 * convf[ifreq]^2

           if (bp[0] ne -1) then begin
              temp[bp] = bad_value
              var[bp] = bad_value
           endif
        
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
;                 help, listp
;                 print, ' '
;                 print, median(map[listp,0])
                 temp[bp[ipix]] = median(map[listp,0]) * 1.e6
                 var[bp[ipix]] = median(map[listp,3]) * 1.e12
              endfor
          endif
;; mollview, temp, /nest, /hist
;; stop         
           ianafast, temp, cls, /nest
           oplot, cls
          ismoothing, temp, tempsm, fwhm_arcmin=sqrt(smoo^2-beam[ifreq]^2), /nest, iter=1
;mollview, tempsm,/nest
          ianafast, tempsm, cls, /nest
          oplot, cls

         ud_grade, tempsm, map, order_out='ring', order_in='nest'
;mollview, map
;;            ud_grade, temp, map, order_out='ring', order_in='nest'
;                            
;stop
;;           gsig = sqrt(smoo^2-beam[ifreq]^2)
;;           bl = gaussbeam(gsig, 2*nside)^2

; --- RMS computed from simulations --> ffp3forcmd_noise.pro
           ud_grade, sqrt(var), rms, order_out='ring', order_in='nest'

;;           ismoothing, map, mapsm, fwhm_arcmin=gsig, /ring

          write_fits_map, 'ffp3_Imap_'+sfreq[ifreq]+'_smth'+strsmoo+'.fits', map, /ring, units='!7l!17K'
          mollview, 'ffp3_Imap_'+sfreq[ifreq]+'_smth'+strsmoo+'.fits', chars=1.5, /hist
;;           write_fits_map,'ffp3_Irms_'+sfreq[ifreq]+'_smoo.fits', rms, /ring, units='!7l!17K'
; ---
;;            write_fits_map, 'ffp3_Imap_'+sfreq[ifreq]+'.fits', map, /ring, units='!7l!17K'
;           mollview, 'ffp3_Imap_'+sfreq[ifreq]+'.fits', chars=1.5, /hist
           write_fits_map,'ffp3_Irms_'+sfreq[ifreq]+'.fits', rms, /ring, units='!7l!17K'                                                                        
;           mollview, 'ffp3_Irms_'+sfreq[ifreq]+'.fits'

;stop
      endfor

  endif

; ----------------------------------------------------------------------

  if (0b) then begin
         nsim = 100

;;          out_dir = '/global/scratch/sd/dpietrob/ffp3/ns2048/mc_noise/'
         out_dir = './'

     for isim=nsim, nsim do begin

        ssim=string(isim,format='(i5.5)')
                                                                                                                               
        for ifreq=4, nfreq-1 do begin
           print, sfreq[ifreq]

           infile = '/project/projectdirs/planck/data/ffp3/outputs/mc_noise/ffp3.nmc.'+ssim+'.'+sfreq[ifreq]+'.e2e2.2048.fits'
           print, 'Reading '+infile
           print, freq[ifreq]
           print, ns[ifreq]
           print, npix[ifreq]
           print, beam[ifreq]

           read_fits_map,infile,map, hdr, xhdr

           sz = size(map)
           print, sz

           bp = where(map[*,0] eq bad_value)
           gp = where(map[*,0] ne bad_value)
           help, bp
;           if (total(bp) gt 0) then stop

           temp = fltarr(12l*ns[ifreq]^2)
;;          var = temp

           temp[gp] = map[gp,0] * 1.e6  * convf[ifreq]
;;          var[gp]  = map[gp,3] * 1.e12 * convf[ifreq]^2

           if (bp[0] ne -1) then begin
              temp[bp] = bad_value
;;              var[bp] = bad_value
           endif

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
;                 help, listp
;                 print, ' '
;                 print, median(map[listp,0])
                 temp[bp[ipix]] = median(temp[listp])
;;                  var[bp[ipix]] = median(map[listp,3]) * 1.e12
              endfor
           endif

;;            ianafast, temp, cls, /nest
;;            oplot, cls
;;            ismoothing, temp, tempsm, fwhm_arcmin=sqrt(smoo^2-beam[ifreq]^2), /nest, iter=1
;mollview, tempsm,/nest
;                                                                                            
;;            ianafast, tempsm, cls, /nest
;;            oplot, cls
;;            ud_grade, tempsm, map, order_out='ring', order_in='nest'
           ud_grade, temp, map, order_out='ring', order_in='nest'
         
; --- RMS computed from simulations --> ffp3forcmd_noise.pro
;;          ud_grade, sqrt(var), rms, order_out='ring', order_in='nest'
           write_fits_map, out_dir+'ffp3_noise_'+ssim+'_'+sfreq[ifreq]+'.fits', map, /ring, units='!7l!17K'
           mollview, out_dir+'ffp3_noise_'+ssim+'_'+sfreq[ifreq]+'.fits', chars=1.5, /hist
;;            write_fits_map,'ffp3_Irms_'+sfreq[ifreq]+'.fits', rms, /ring, units='!7l!17K'
;;            write_fits_map, out_dir+'ffp3_noise_'+ssim+'_'+sfreq[ifreq]+'_smth'+strsmoo+'.fits', map, /ring, units='!7l!17K'
;;            mollview, out_dir+'ffp3_noise_'+ssim+'_'+sfreq[ifreq]+'_smth'+strsmoo+'.fits', chars=1.5, /hist



;stop

        endfor

     endfor

  end
  
; ----------------------------------------------------------------------

  if (0b) then begin
      init_healpix
      mollview, randomn(-1,12), window=-1
      loadct, 39
      !p.color=0
      !p.background=255

      cls = fltarr(1501,nfreq)

      l = findgen(1501)
      ll = l*(l+1)/2./!pi
      window, 1
      plot_oo, /nodata, [1,1500], [1,100000], xs=1, ys=1, chars=1.5

      for ifreq=2,nfreq-1 do begin
          ianafast, 'ffp3_Imap_'+sfreq[ifreq]+'_smth'+strsmoo+'.fits', tmpcl ;, nlmax=1500, theta_cut_deg=30, regression=2, iter=2, /double
          cls[*,ifreq] = tmpcl[0:1500]
          oplot, l, tmpcl*ll*2., col=200+5*(ifreq-2)
      endfor

  endif



  if (0b) then begin

      for ifreq=2,nfreq-1 do begin
          read_fits_map,'4ffp3_Irms_'+sfreq[ifreq]+'_1000sims_.fits', rms
          rms = rms * convf[ifreq]
          write_fits_map, 'ffp3_Irms_'+sfreq[ifreq]+'_1000sims.fits', rms, /ring, units='!7l!17K cmb'
          mollview,  'ffp3_Irms_'+sfreq[ifreq]+'_1000sims.fits', chars=1.5

      endfor

  endif


 if (0b) then begin

      for ifreq=2,nfreq-1 do begin
          ud_grade,'ffp3_Imap_'+sfreq[ifreq]+'_smth30.fits','ffp3_Imap_'+sfreq[ifreq]+'_smth30_ns512.fits', nside_out=512, order_out='ring'
          read_fits_map,'ffp3_Irms_'+sfreq[ifreq]+'_1000sims.fits', rms
          ud_grade, rms^2, rms512, nside_out=512, order_in='ring', order_out='ring'
          write_fits_map, 'ffp3_Irms_'+sfreq[ifreq]+'_ns512.fits', sqrt(rms512), /ring, units='!7l!17K cmb'
          mollview,  'ffp3_Irms_'+sfreq[ifreq]+'_ns512.fits', chars=1.5

      endfor

  endif

  if (0b) then begin

      ns = 512l
      npix = 12l*ns^2
      scaling = [1,1,1,1,1,1,5]

      mollview, randomn(-1,12), window=-1
      loadct, 39
      !p.color=0
      !p.background=255
      window,1 & plot_oo, /nodata, [1,1024], [1.e-6,1.e4], chars=1.5, xs=1, ys=1

      for ifreq=2,nfreq-1 do begin
          read_fits_map, 'ffp3_Imap_'+sfreq[ifreq]+'_smth30_ns512.fits', map
          read_fits_map, 'noise_sims/ffp3_Irms_'+sfreq[ifreq]+'_ns512_1000sims.fits', rms

          noise = randomn(-(ifreq+1),npix) * rms * scaling[ifreq]
          map = map + noise

          write_fits_map, 'ffp3_map_'+sfreq[ifreq]+'_smth30_ns512_Nadd.fits', map, /ring
          write_fits_map, 'ffp3_rms_'+sfreq[ifreq]+'.fits', rms*scaling[ifreq], /ring
          ianafast, 'ffp3_map_'+sfreq[ifreq]+'_smth30_ns512_Nadd.fits', cls
          oplot, cls, col=40*ifreq
          oplot, cls*0+total(rms^2)/npix*4.*!pi/npix*scaling[ifreq]^2, col=40*ifreq

      endfor

  endif

  if (true) then begin
     print, " Computation of the effective beam based on combination of Gaussians"
     
     read_fits_map, 'mask_80_ps.fits', mask
     gp = where(mask gt 0.)
     bp = where(mask eq 0.)

;;      mollview, mask
     read_fits_map, 'ffp3_Irms_143.fits', r143
     read_fits_map, 'ffp3_Irms_217.fits', r217
     read_fits_map, 'ffp3_Irms_353.fits', r353

     bl143 = gaussbeam(beam[4],4096)
     bl217 = gaussbeam(beam[5],4096)
     bl353 = gaussbeam(beam[6],4096)

     b3 = ( bl143 * 1./mean( r143[gp] ) + bl217 * 1./mean( r217[gp] ) + bl353 * 1./mean( r353[gp] ) ) / ( 1./mean( r143[gp] ) + 1./mean( r217[gp] ) + 1./mean( r353[gp] ) )
     b2 = ( bl217 * 1./mean( r217[gp] ) + bl353 * 1./mean( r353[gp] ) ) / ( 1./mean( r217[gp] ) + 1./mean( r353[gp] ) )

;;      mollview, randomn(-1,12), window=-1
     loadct, 39
     !p.color=0
     !p.background=255
;;      window, 1 & plot, bl143
;;      oplot, bl143, col=200
;;      oplot, bl217, col=245
;;      oplot, bl353, col=70

;;      oplot, b3, thick=2
;;      oplot, b2, line=2, thick=2

     cl2fits, b3, 'ffp3_3b_effBeam.fits'
     openw,1,'ffp3_3b_effBeam.dat'
     for i=0l,4096 do printf,1,i,b3[i]
     close,1

     cl2fits, b2, 'ffp3_2b_effBeam.fits'
     openw,1,'ffp3_2b_effBeam.dat'
     for i=0l,4096 do printf,1,i,b2[i]
     close,1


  endif


end
