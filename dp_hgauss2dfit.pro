FUNCTION DP_HGAUSS2DFIT, ROBUST = ROBUST, DETECT = DETECT, BAND = BAND, $
                      R_FWHM = R_FWHM, R_MAX = R_MAX, OBJECT = OBJECT, $
                      COORD = COORD, ECLIPTIC = ECLIPTIC, $
                      CATALOG = CATALOG, L_CUT = L_CUT, B_CUT = B_CUT, $
                      MAP_FILE = MAP_FILE, VERSION = VERSION, END_RING = END_RING, LOCAL = LOCAL, $
                      MADAM = MADAM, HALF = HALF, SURVEY = SURVEY, $
                      PLOT = PLOT, NESTED = NESTED, QUIET=QUIET, $
                      ONLY_T = ONLY_T, $
                      SCALAR = SCALAR, LINEAR = LINEAR, PARABOLIC = PARABOLIC, $
                      SAVE = SAVE, OUTFILE = OUTFILE, DO_STOP = DO_STOP, JPEG = JPEG, $
                      catalog_infile=catalog_infile, map_infile=map_infile, varmap_infile=varmap_infile, first_source=first_source, $
                      input_fwhm=input_fwhm, map_calibration=map_calibration, last_source=last_source, do_all_cases=do_all_cases, $
                         is_rms=is_rms, tag=tag ;, $
; , DIR = DIR, ERROR = ERROR, DATA = DATA, 

; --- :DP
  true = 1b
  false = 0b
  undef = -1b

  !path=!path+':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'
  !path=!path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'
; ---

; 2 dimensional Gaussian fit in a genuine Healpix format
; For typical sources, the data are accessed via the cropped maps. However, to save space in the 
; disk, IF THE KEYWORD COORD IS USED, THE MAPS WILL BE OPENNED AND the data for the fit will not be saved.
; One may give by hand a DIR and an OBJNAME ( array if POL is considered )
; If not, read_info_map will select and OPTION:
; 1.- Official maps through VERSION and END_RING
; 2.- CROP: set DIR = 'crop'
; 3.- Madam: MADAM = [ str1, str2, ... ]
; 4.- Own computer
; 5.- MAP_INFILE: allows for a local FITS file. Include all the path.
; Initiating
; By default background is a scalar

  If ( (NOT Keyword_set( LINEAR )) AND (NOT Keyword_set( PARABOLIC )) ) then SCALAR = 1
  If ( NOT Keyword_set( R_FWHM ) ) then R_FWHM = 2.5d0 ; related with R_MAX (arcmin)
  If ( NOT Keyword_set( ROBUST ) ) then ROBUST = true
  If ( NOT Keyword_set( VERSION ) ) then VERSION = 'D.Pietrobon - Modified from S.Hildebrandt'
  If ( (NOT Keyword_set( OBJECT )) AND (NOT Keyword_set( COORD )) AND (NOT keyword_set(catalog_infile)) ) then OBJECT = 'Crab'

  If ( NOT Keyword_set( BAND ) and (not keyword_set(map_infile)) ) then BAND = [ 30, 44, 70, 100, 143, 217, 353, 545, 857 ]

; --- :DP
  if ( not keyword_set(do_all_cases) ) then do_all_cases = false

  if (do_all_cases) then begin
      PARABOLIC = true
      N_PAR = 12
  endif else begin
      if ( keyword_set(SCALAR) ) then print, SCALAR
      if ( keyword_set(LINEAR) ) then print, LINEAR
      if ( keyword_set(PARABOLIC) ) then print, PARABOLIC
      If ( Keyword_set( SCALAR ) ) then N_PAR = 7 else If ( Keyword_set( LINEAR ) ) then N_PAR = 9 else if ( Keyword_set( PARABOLIC ) ) then N_PAR = 12 else begin
          print, ' - :DP No procedure selected.'
          print, '       SCALAR set default'
          PARABOLIC = false
          LINEAR    = true
          SCALAR    = false
      endelse
      print, 'N_PAR = ', N_PAR
  endelse

  if ( not keyword_set(is_rms) ) then is_rms = false

  if ( not keyword_set(only_t) ) then only_t = true

  if ( keyword_set(only_t) and (only_t eq false) ) then stop, ' - :DP WARNING: Polarization not checked!!!'

  if ( not keyword_set(first_source)) then first_source = 0l

  if ( not keyword_set(map_calibration) ) then map_calibration = 1.

  if ( keyword_set(map_infile) ) then begin
      print, ' - :DP map_infile = ', strtrim(map_infile) 
      if ( not keyword_set(input_fwhm) ) then stop, ' - :DP FWHM of the channel not set'
      if ( keyword_set(varmap_infile) ) then print, ' - :DP varmap_infile = ', strtrim(varmap_infile) else if ( keyword_set(map_infile) ) then print, ' - :DP map_infile provided, variance not!'
      if ( not keyword_set(varmap_infile) ) then varmap_infile = undef
;      band = ['From File']
  endif else begin
      map_infile = undef
  endelse

  if ( not keyword_set(save) ) then save=false

  if ( not keyword_set(outfile) and save) then outfile = 'test.sav' else outfile=''

  If ( Keyword_set(COORD) and (not keyword_set(catalog_infile)) ) then begin
      tag = tag + '_coord'
      object = 'Coordinate:'
      N_SOURCES   = n_elements( reform( coord[ *, 0 ] ) )
      N_SOURCES_2 = n_elements( reform( coord[ *, 1 ] ) )
      help, n_sources

      If (n_sources ne n_sources_2) then stop,'Number of sources and coordinates not compliant'
      ALL_L_OBJECT = reform( coord[ *, 0 ] ) 
      ALL_B_OBJECT = reform( coord[ *, 1 ] )
      If (Keyword_set( ECLIPTIC ) ) then begin
          euler, all_l_object, all_b_object, all_l_object_ecl, all_b_object_ecl, 6 ; Galactic to Ecliptic
          ALL_L_OBJECT = all_l_object_ecl
          ALL_B_OBJECT = all_b_object_ecl
      Endif
  Endif else begin                         ; COORD
      COORD = false
      print, 'COORD = ', COORD
  endelse

  if ( keyword_set(catalog_infile) and (not keyword_set(coord)) ) then begin
      if (keyword_set(coord)) then stop, ' - :DP Inconsitency: both COORD and CATALOG_INFILE.'
      radec2gal = 1
      print, 'CATALOG_INFILE keyword set.'
      print, 'reading source catalog from: ', strtrim(catalog_infile,2)

      read_fits_s, catalog_infile, hdr, table

      flux = table.FLUX
      fluxerr = table.FLUX_ERR
      ra = table.RA
      dec =table.DEC
      extended = table.EXTENDED

      read_fits_map, map_infile, xxx, nside=mapns
      
      if (n_elements(BAND) gt 1) then stop, ' -:DP Inconsistency BAND and INFILE'

; catalog driven
;      cf = conversionfactor(float(band),/brightness2muKthermo)
;      print, ' conversion: ', cf
;      mapns=1024l
;      npix = 12l*mapns^2
;      parea = 4.*!pi / npix
;      cflux = flux[gsource] * cf * 1.e-9/parea

      s2n = flux / fluxerr

;      print, minmax(ra)
;      print, minmax(dec)

;      help, flux

      gsource = where( (extended eq 0) and (s2n ge 5.) )
      print, ' - :DP Thresholding done on EXTENDED and S2N...'
      help, gsource
;gsource = gsource[0:0]

      isource = reverse(sort(flux[gsource]))
      
      ngsource = n_elements(gsource)

      gsource = gsource[isource]

;      print, 'ra ',  minmax( ra[gsource] )
;      print, 'dec ', minmax( dec[gsource] )

; DX7 alteady in Galactic
;      euler, ra[gsource], dec[gsource], phi, theta, radec2gal
      theta = table.GLAT[gsource] ;ra[gsource]
      phi = table.GLON[gsource] ;dec[gsource]

      print, 'theta ',minmax(theta)
      print, 'phi ',minmax(phi)

      coord = [[phi], [theta]]
      help, coord
      ALL_L_OBJECT = reform( coord[ *, 0 ] )
      ALL_B_OBJECT = reform( coord[ *, 1 ] )
      N_SOURCES = ngsource

      file_arr = strsplit(catalog_infile,'/', /extract)
      cfile = strtrim( file_arr[n_elements(file_arr)-1], 2 )
      print, cfile
;stop
      cat_outfile = strjoin(strsplit( cfile, '.fits', /extract ),'_')
      openw, 1, cat_outfile + '.txt'
;      openw,1,'dx7_catalog_5stnr_'+band[0]+'.txt'
      printf, 1, n_elements(gsource) 
      for i=0,n_elements(gsource)-1 do printf,1,strtrim(string(phi[i]),2)+':'+strtrim(string(theta[i]),2),', ', phi[i],', ', theta[i],', ', 15.,', ', flux[gsource[i]], $
        format='(a30,a2,f10.4,a2,f10.4,a2,f10.4,a2,f12.4)'
      close,1
      print, ' - :DP Catalog written in format...'
      object = string(lindgen(ngsource))
;help, object
;      stop

  endif
; ---

  on_error, 1

  init_healpix
  mollview, randomn(-1,12), window=-1
  loadct, 39
  !p.color=0
  !p.background=255

  if (strlen(strtrim(map_infile)) le 1) then begin
      print, ' :DP - No map_infile supplied!'
      FWHM = [ 36d0, 23d0, 14d0, 10d0, 7d0, 5d0, 5d0, 5d0, 5d0 ] ; arcmin
      BAND_TMP = [ '30','44', '70', '100', '143', '217', '353', '545', '857' ]
      N_BAND = n_elements( band )
      print, ' - :DP n_band = ', n_band
  endif else begin
      file_arr = strsplit(map_infile,'/',/extract)
      BAND = file_arr[n_elements(file_arr)-1]
      N_BAND = 1
      FWHM = input_fwhm
      BAND_TMP = undef
      CAL_0 = map_calibration
  endelse

  If ( NOT Keyword_set( R_MAX ) ) then R_MAX = r_fwhm * fwhm  ; arcmin ( HFI enough for point sources) 

; --- :DP ----------------------------------------------------------------------

  if ( not keyword_set(DETECT) )          then DETECT = undef
  if ( not keyword_set(ECLIPTIC) )      then ECLIPTIC = undef
  if ( not keyword_set(CATALOG) )       then CATALOG = undef
  if ( not keyword_set(L_CUT) )         then L_CUT = undef
  if ( not keyword_set(B_CUT) )         then B_CUT = undef
  if ( not keyword_set(MAP_FILE) )      then MAP_FILE = undef
  if ( not keyword_set(END_RING) )      then END_RING = undef
  if ( not keyword_set(LOCAL) )         then LOCAL = undef
  if ( not keyword_set(VERSION) )       then VERSION = undef
  if ( not keyword_set(MADAM) )         then MADAM = undef
  if ( not keyword_set(HALF) )          then HALF = undef
  if ( not keyword_set(SURVEY) )        then SURVEY = undef
  if ( not keyword_set(NESTED) )        then NESTED = undef
  if ( not keyword_set(SCALAR) )        then SCALAR = undef
  if ( not keyword_set(LINEAR) )        then LINEAR = undef
  if ( not keyword_set(PARABOLIC) )     then PARABOLIC = undef
  if ( not keyword_set(SAVE) )          then SAVE = undef
  if ( not keyword_set(JPEG) )          then JPEG = undef
  if ( not keyword_set(varmap_infile) ) then varmap_infile = undef
  if ( not keyword_set(catalog_infile) ) then catalog_infile = undef
  if ( not keyword_set(QUIET) )          then QUIET = true
  if ( not keyword_set(PLOT) )          then PLOT = 'test'
  if ( not keyword_set(DO_STOP) )          then DO_STOP = false
  
;  print, 'HGAUSS2DFIT Variable:'
;  print, 'ROBUST = ', ROBUST
;      print, ' DATA = ', DATA 
      print, 'DETECT = ', DETECT 
      print, 'BAND = ', BAND 
      print, 'R_FWHM = ', R_FWHM 
      print, 'R_MAX = ', R_MAX 
      print, 'OBJECT = ', size( OBJECT ) 
      print, 'COORD = ', size( COORD ) 
      print, 'ECLIPTIC = ', ECLIPTIC 
      print, 'CATALOG = ', CATALOG  
      print, 'L_CUT = ', L_CUT 
      print, 'B_CUT = ', B_CUT 
      print, 'MAP_FILE = ', MAP_FILE 
      print, 'VERSION = ', VERSION 
      print, 'END_RING = ', END_RING 
      print, 'LOCAL = ', LOCAL 
      print, 'MADAM = ', MADAM 
      print, 'HALF = ', HALF 
;      print, 'DIR = ', DIR 
      print, 'SURVEY = ', SURVEY 
      print, 'PLOT = ', PLOT 
      print, 'NESTED = ', NESTED 
      print, 'QUIET = ', QUIET 
;    'ERROR = ', ERROR 
      print, 'ONLY_T = ', ONLY_T 
      print, 'SCALAR = ', SCALAR 
      print, 'LINEAR = ', LINEAR 
      print, 'PARABOLIC = ', PARABOLIC 
      print, 'SAVE = ', SAVE 
      print, 'OUTFILE = ', OUTFILE 
      print, 'DO_STOP = ', DO_STOP 
      print, 'JPEG = ', JPEG 
      print, 'catalog_infile = ', strtrim(catalog_infile,2) 
      print, 'map_infile = ', strtrim(map_infile,2) 
      print, 'varmap_infile = ', strtrim(varmap_infile,2) 
      print, 'do_all_cases = ', DO_ALL_CASES
;stop
  print, ' --- :DP ---'
  print, ''
; ----------------------------------------------------------------------

; Opening the maps, and then lopping over objects
; Defining the set of maps to be considered at each band
  N_MAPS = 4
; --- :DP
;  If Keyword_set( MAP_FILE ) then N_MAPS = 1
;  If ( Keyword_set( MAP_FILE ) or keyword_set(map_infile) ) then N_MAPS = 1
  If ( ONLY_T ) then N_MAPS = 1
  print, ' - :DP N_MAPS = ', N_MAPS
; ---
; To be re-defined when Calatog is used ( see below )
  PAR_FIT = dblarr( n_sources, n_band, n_maps, n_par + 2 ) ; simply because IDL does not concatenate different dimensions autmotically
                                                        ; the last two dimensions do not play any role but for returning the result
                                                        ; in a matrix with par_fit_err together i
                                                        ; (if changed, change PAR_FIT after RES)
  PAR_FIT_ERR = dblarr( n_sources, n_band, n_maps, n_par + 2 ) ; Chi^2 and DoF added to the error estimates of each parameter
; Some generic quantities. Re-defined when usinga catalog (diff # of sources per band)
  L_MAX_VEC = dblarr( n_sources, n_band )
  COB_MAX_VEC = l_max_vec
  FWHM_AVG_ARCMIN_VEC = dblarr( n_sources, n_band, n_maps )
  FWHM_AVG_ARCMIN_ERR_VEC = fwhm_avg_arcmin_vec
  ELIP_VEC = fwhm_avg_arcmin_vec
  ELIP_ERR_VEC = fwhm_avg_arcmin_vec
  DELTA_L_VEC = fwhm_avg_arcmin_vec
  DELTA_B_VEC = fwhm_avg_arcmin_vec
  FLUX_JY_VEC = fwhm_avg_arcmin_vec
  FLUX_JY_ERR_VEC = fwhm_avg_arcmin_vec

; LOOP OVER BANDS/DETECT
  I_1 = 0
  I_2 = n_band - 1
  For i_band = i_1, i_2 do begin 
      if ( strlen(strtrim(map_infile)) le 1 ) then begin
          print, ' map_infile undefined', map_infile

      endif else begin
; --- :DP --------------------------------------------------------------
  ; MAP_INFILE reading file. MAP_INFILE is a string containing path and fits_name
          print, ' - :DP Reading map from: ', strtrim(map_infile,2)
          read_fits_map, map_infile, map_fits, nside=dpns, ordering=dporder
          n_side=dpns
          npix = 12l*n_side^2
          detect_tmp = band[i_band]
;      mollview, map_infile, chars=1.5, min=-400, max=400
          print, ' - :DP RING order assumed...'
          print, dporder
          if (dporder eq 'nested') then begin
              print, ' - :DP Reordering NEST --> RING'
              ud_grade, map_fits, map_ring, order_in=dporder, order_out='ring'
              map_fits = map_ring
          endif
          bpix = where(map_fits eq -1.6375e30)
          if (bpix[0] ne -1) then map_fits[bpix] = 0.
          bpix = finite( map_fits, /nan )
          if (bpix[0] ne -1)then map_fits[bpix] = 0.
;          help, map_fits
;          ismoothing, map_fits, smth_map, fwhm_arcmin=sqrt(40.^2-fwhm^2), /ring, /silent
;          mollview, smth_map, chars=1.5, min=-500, max=500, win=3
;          mollview, map_fits, chars=1.5, min=-500, max=500, win=4
          map_stokes = dblarr(n_maps, 12l*dpns^2)
          help, map_stokes
          if ( strlen(strtrim(varmap_infile)) ge 2 ) then read_fits_map, varmap_infile, var_fits, nside=dpns, ordering=dporder else var_fits = map_fits*0.+1.
          bpix = where(var_fits eq -1.6375e30)
          if (bpix[0] ne -1)then begin
              print, ' WARNING bad pixels found'
              var_fits[bpix] = 1.6375e30
          endif
          bpix = finite( var_fits, /nan )
          if (bpix[0] ne -1)then begin
              print, ' WARNING NAN pixels found'
              var_fits[bpix] = 1.6375e30
          endif

          if (is_rms) then var_fits = var_fits^2
          if (dporder eq 'nested') then begin
              print, ' - :DP Reordering NEST --> RING'
              ud_grade, var_fits, var_ring, order_in=dporder, order_out='ring'
              var_fits = var_ring
          endif

          map_stokes_var = map_stokes
;stop
; ---
          map_stokes[ 0, * ] = reform( map_fits[ *, 0 ] )
          map_stokes_var[ 0, * ] = reform( var_fits[ *, 0 ] )
          if (not only_t) then begin
              map_stokes[ 1, * ] = reform( map_fits[ *, 1 ] )
              map_stokes[ 2, * ] = reform( map_fits[ *, 2 ] )
              map_stokes[ 3, * ] = reform( sqrt( ( map_fits[ *, 1 ] )^2 + ( map_fits[ *, 2 ] )^2 ) )
              undefine, map_fits
          endif

          q_band = 0
          
      Endelse
; ----------------------------------------------------------------------
; Eliminating Healpix bad values from the 1 (4) maps
      I_FINAL = n_maps - 1
      For i_tmp = 0, i_final do begin 
          Q_TMP = where( abs( reform( map_stokes[ i_tmp, * ] ) ) gt 1d30 )
          If q_tmp[ 0 ] ne -1 then map_stokes[ i_tmp, q_tmp ] = 0d0
      Endfor

      ps_map = fltarr(12l*dpns^2)

      if ( (last_source eq -1) ) then last_source = n_sources-1

  ;  Source by source (faster for several sources, to read the map only once for each band and then iterate over sources)
      print, ' :DP - Starting Loop over Sources...'
      
      chisq_vec = fltarr(n_sources, 6)

      For source = first_source, last_source do begin
          if ( (source/50)*50 eq source ) then begin
              print, 'Saving map, Source ',source,'/',n_sources-1
              write_fits_map, tag+'_ps_map.fits', ps_map, /ring, units='!7l!8K CMB'
              ;mollview, tag+'_ps_map.fits', win=25, max=750, chars=1.5
              cl2fits, chisq_vec, tag+'_chisq_vec.fits'
;stop
          endif
          L_OBJECT = all_l_object[ source ]
          B_OBJECT = all_b_object[ source ]

          chisq_vec[source,0] = l_object
          chisq_vec[source,1] = b_object

          print,''
          print,'SOURCE # ', str( 1 + source ) +'/' + str( n_sources ), l_object, b_object
          print,''
          
;;          If ( ( ECLIPTIC eq undef) ) then print, 'INPUT COORDINATES FOR THE OBJECT (GAL): ', str( [l_object, b_object] ) else $
;;            print, 'INPUT COORDINATES FOR THE OBJECT (ECL): ', str( [l_object, b_object] )

  ; Redefining the region

          ang2pix_ring, n_side, !dpi/2d0 - b_object /180d0 *!dpi, l_object * !dpi/180d0, I_PIX_SRC

          TAIL_PLOT = [ ' I ', ' U ', ' Q ', ' P ' ]
          For s = 0, n_maps - 1 do begin
              MAP_TMP = reform( map_stokes[ s, * ] )
              RMS_TMP = reform( sqrt(map_stokes_var[ s, * ]) )

  ; Finding the maximum in Intensity (Pol is noisier)
              If (s eq 0) then begin 
                  SIGN_TMP = sign( map_tmp[ i_pix_src[0] ] )
                  MAP_TMP = sign_tmp * map_tmp
                  pix2vec_ring, n_side, i_pix_src, VEC_SRC

                  query_disc, n_side, vec_src,  ( 0.5 * fwhm[ q_band ] ) / 60d0 / 180d0 * !dpi, LIST_PIX_MAX 
;;                   query_disc, n_side, vec_src,  ( 1. * fwhm[ q_band ] ) / 60d0 / 180d0 * !dpi, LIST_PIX_MAX 

                  print, ' MAX VALUE nearby maximum: ', max( map_tmp[ list_pix_max ], Q_MAX), Q_MAX
                  print, ' MIN VALUE nearby maximum: ', min( map_tmp[ list_pix_max ], Q_MIN), Q_MIN
                  
;stop
                  If ( (q_max[ 0 ] eq -1) OR (n_elements( q_max ) gt 2)) then GOTO, NO_OBSERVED
                  I_PIX_MAX = list_pix_max[ q_max[ 0 ] ]
                  I_PIX_MIN = list_pix_max[ q_min[ 0 ] ]
                  pix2vec_ring, n_side, i_pix_max, VEC_MAX
              Endif             ; s eq 0
; Creating the array of data
              TAIL_NESTED='RING'

              query_disc, n_side, reform( vec_max ), ( r_max[ q_band ] )[ 0 ] / 60d0 / 180d0 * !dpi, LIST_L_MAX 

              SOURCE_DATA = reform( map_stokes[ s, list_l_max ] )
              SOURCE_DATA_ERR = sqrt( reform( map_stokes_var[ s, list_l_max ] ) ) 
              pix2ang_ring, n_side, i_pix_max, COB_MAX, L_MAX
              l_max_vec[ source, i_band ] = l_max[ 0 ]
              cob_max_vec[ source, i_band ] = cob_max[ 0 ]

  ; Special case of polarization intensity (I need to incorporate COV(Q,U)
              If (s eq 3) then begin
                  SOURCE_DATA  = sqrt( reform( map_stokes[ 1, list_l_max ] )^2d0 + reform( map_stokes[ 2, list_l_max ] )^2d0 )
                  SOURCE_DATA_ERR = abs( reform( map_stokes[ 1, list_l_max ] ) ) / source_data * sqrt( reform( map_stokes_var[ 1, list_l_max ] ) )+ $
                    abs( reform( map_stokes[ 1, list_l_max ] ) ) / source_data * sqrt( reform( map_stokes_var[ 2, list_l_max ] ) )
              Endif             ; Err in s=3 (pol int)
; 2D Gaussian Healpix fit
              PARINFO = replicate( { LIMITED: [ 0, 0 ], LIMITS: [0d0, 0d0], VALUE: 0d, FIXED: 0 }, n_par )
; Better leave room for a background of different sign that the source ( Q or U may be < 0, while >0 in I, and necessarily in P )
; PARAMETERS
; 0: Background, constant
; 1: Amplitude of the 2D Gaussian
; 2: FWHM_LONGITUDE
; 3: FWHM_LATITUDE
; 4: CENTER_LONGITUDE
; 5: CENTER_LATITUDE
; 6: TILT ELLIPSE 

; BACKGROUND (un-comment to make it positive -not in general!/ possibility of imposing it only for P)
; If s eq 3 then begin
;  parinfo[ 0 ].limited[ 0 ] = 1
;  parinfo[ 0 ].limits[ 0 ]  = 0d0 
; Endif

; AMPLITUDE 2D GAUSSIAN
; If s eq 3 then begin
; --- :DP
              parinfo[ 1 ].limited[ 1 ] = 1b
;              parinfo[ 1 ].limits[ 1 ]  = max( [map_tmp[ i_pix_max ], map_tmp[i_pix_max]-min([map_tmp[i_pix_min],-500]) ] ) ; max value is the maximum of the region (not excellent but to avoid troubles with extended sources)
              parinfo[ 1 ].limits[ 1 ]  = map_tmp[ i_pix_max ] - min( [map_tmp[i_pix_min], -400] )
              parinfo[ 1 ].limited[ 0 ] = 1b
              parinfo[ 1 ].limits[ 0 ]  = 0.d0

;              parinfo[1].fixed = true
 
              print, ' Amplitude range = ---> ', [parinfo[ 1 ].limits[ 0 ], parinfo[ 1 ].limits[ 1 ]]

; Endif

; FWHM_LONGITUDE
              parinfo[ 2 ].limited[ 0 ] = 1
              parinfo[ 2 ].limits[ 0 ] = fwhm[ q_band ]/ 60d0 / 180d0 * !dpi*0.5d0 ;0d0
              parinfo[ 2 ].limited[ 1 ] = 1
              parinfo[ 2 ].limits[ 1 ] = fwhm[ q_band ]/ 60d0 / 180d0 * !dpi*1.5d0
; FWHM_LATITUDE
              parinfo[ 3 ].limited[ 0 ] = 1
              parinfo[ 3 ].limits[ 0 ] = fwhm[ q_band ]/ 60d0 / 180d0 * !dpi*0.5d0 ;0d0
              parinfo[ 3 ].limited[ 1 ] = 1
              parinfo[ 3 ].limits[ 1 ] = fwhm[ q_band ]/ 60d0 / 180d0 * !dpi*1.5d0
; TILT (angle wrt X-Y axes)
              parinfo[ 6 ].limited[ 0 ] = 1
              parinfo[ 6 ].limited[ 1 ] = 1
              parinfo[ 6 ].limits[ 0 ] = 0.d0
              parinfo[ 6 ].limits[ 1 ] = !dpi/2.d0 ;6.28318d0 ; 2d0 * !dpi
;              parinfo[ 6 ].fixed = true
; First values
              A_0 = dblarr( n_par )
;              A_0[ 0 : 1 ] = map_tmp[ i_pix_max ] # [1d0/100d0, 1d0 ] ; It can be improved 
              A_0[ 0 : 1 ] = [ mean(map_tmp[list_l_max]), map_tmp[i_pix_max]*0.9 ]
; rad.
              A_0[ 2 : 6 ] = [ fwhm[ q_band ] / 60d0 / 180d0 * !dpi, fwhm[ q_band ] / 60d0 / 180d0 * !dpi, $
                               0.d0, 0.d0,  !dpi/8.d0 ]

;print, 'do_all_cases', do_all_cases
;stop
if (do_all_cases) then begin
; --- scalar
;    print, 'scalar...'
    npar = 7
    
              DATA_FIT_scal = { ROBUST : robust, $
                           DATA : source_data, $
                           ERR_DATA : source_data_err, $
                           LIST_PIX : list_l_max, $
                           L_MAX : l_max, $
                           COB_MAX : cob_max, $
                           N_SIDE : n_side, $
                           NESTED : false, $
                           N_PAR : npar }
             
              RES_scal= mpfit( 'dp_funct_fit_hgauss2d', A_0, FUNCTARGS = data_fit_scal, PARINFO  = parinfo, $
                          PERROR = perror, FASTNORM = FASTNORM,  AUTODERIVATIVE = AUTODERIVATIVE, $
                          MAXITER = MAXITER, BESTNORM = CHI2_scal, $
                          COVAR = COVAR, QUIET = QUIET, /NOCATCH, DOF= DOF, $
                          ERRMSG= ERRMSG, STATUS = STATUS, NITER = 200 )

    tres = res_scal
    tres[2] = tres[2] * 60. * 180. / !pi
    tres[3] = tres[3] * 60. * 180. / !pi
    tres[4] = tres[4] * 60. * 180. / !pi
    tres[5] = tres[5] * 60. * 180. / !pi
    tres[6] = tres[6] * 180. / !pi
    
    amp_scal = tRES[1]
    angle_scal = tRES[6]
    if (not quiet) then begin
        print, ' --- --- ---'
        print, ' Resulting parameters: Scalar'
        print, ' background offset: ', tRES[0]
        print, ' Gauss Amplitude: ', tRES[1]
        print, ' y-fwhm arcmin: ', tRES[2]
        print, ' x-fwhm arcmin: ', tRES[3]
        print, ' peak center x: ', tRES[4]
        print, ' peak center y: ', tRES[5]
        print, ' Rotation angle: ', tRES[6]
    endif
; --- linear
;    print, 'linear...'
    npar = 9          
              DATA_FIT_lin = { ROBUST : robust, $
                           DATA : source_data, $
                           ERR_DATA : source_data_err, $
                           LIST_PIX : list_l_max, $
                           L_MAX : l_max, $
                           COB_MAX : cob_max, $
                           N_SIDE : n_side, $
                           NESTED : false, $
                           N_PAR : npar }
             
              RES_lin= mpfit( 'dp_funct_fit_hgauss2d', A_0, FUNCTARGS = data_fit_lin, PARINFO  = parinfo, $
                          PERROR = perror, FASTNORM = FASTNORM,  AUTODERIVATIVE = AUTODERIVATIVE, $
                          MAXITER = MAXITER, BESTNORM = CHI2_lin, $
                          COVAR = COVAR, QUIET = QUIET, /NOCATCH, DOF= DOF, $
                          ERRMSG= ERRMSG, STATUS = STATUS, NITER = 200 )

    tres = res_lin
    tres[2] = tres[2] * 60. * 180. / !pi
    tres[3] = tres[3] * 60. * 180. / !pi
    tres[4] = tres[4] * 60. * 180. / !pi
    tres[5] = tres[5] * 60. * 180. / !pi
    tres[6] = tres[6] * 180. / !pi
    amp_lin = tRES[1]
    angle_lin = tRES[6]
    if (not quiet) then begin
        print, ' --- --- ---'
        print, ' Resulting parameters: Linear'
        print, ' background offset: ', tRES[0]
        print, ' Gauss Amplitude: ', tRES[1]
        print, ' y-fwhm arcmin: ', tRES[2]
        print, ' x-fwhm arcmin: ', tRES[3]
        print, ' peak center x: ', tRES[4]
        print, ' peak center y: ', tRES[5]
        print, ' Rotation angle: ', tRES[6]
    endif
; --- parabolic
;    print, 'parabolic...'
    npar = 12
              DATA_FIT_par = { ROBUST : robust, $
                           DATA : source_data, $
                           ERR_DATA : source_data_err, $
                           LIST_PIX : list_l_max, $
                           L_MAX : l_max, $
                           COB_MAX : cob_max, $
                           N_SIDE : n_side, $
                           NESTED : false, $
                           N_PAR : npar }
             
              RES_par= mpfit( 'dp_funct_fit_hgauss2d', A_0, FUNCTARGS = data_fit_par, PARINFO  = parinfo, $
                          PERROR = perror, FASTNORM = FASTNORM,  AUTODERIVATIVE = AUTODERIVATIVE, $
                          MAXITER = MAXITER, BESTNORM = CHI2_par, $
                          COVAR = COVAR, QUIET = QUIET, /NOCATCH, DOF= DOF, $
                          ERRMSG= ERRMSG, STATUS = STATUS, NITER = 200 )
    
    tres = res_par
    tres[2] = tres[2] * 60. * 180. / !pi
    tres[3] = tres[3] * 60. * 180. / !pi
    tres[4] = tres[4] * 60. * 180. / !pi
    tres[5] = tres[5] * 60. * 180. / !pi
    tres[6] = tres[6] * 180. / !pi
    amp_par = tRES[1]
    angle_par = tRES[6]
    if (not quiet) then begin
        print, ' --- --- ---'
        print, ' Resulting parameters: Parabolic'
        print, ' background offset: ', tRES[0]
        print, ' Gauss Amplitude: ', tRES[1]
        print, ' y-fwhm arcmin: ', tRES[2]
        print, ' x-fwhm arcmin: ', tRES[3]
        print, ' peak center x: ', tRES[4]
        print, ' peak center y: ', tRES[5]
        print, ' Rotation angle: ', tRES[6]
    endif

    chi2_solution = [chi2_scal, chi2_lin, chi2_par]
    amp_solution = [amp_scal, amp_lin, amp_par]

    if (not quiet) then begin
        print, ' chisq/ampl for different cases:'
        print, 'Scalar:    ', chi2_solution[0], ' --> ', amp_solution[0]
        print, 'Linear:    ', chi2_solution[1], ' --> ', amp_solution[1]
        print, 'Parabolic: ', chi2_solution[2], ' --> ', amp_solution[2]
    endif
 
    print, 'Minimum: ', min(chi2_solution, im)
     
    chi2 = chi2_solution[im]

    chisq_vec[source,2] = chi2 / dof

    if ( (amp_scal eq amp_lin) and (amp_lin eq amp_par) ) then chisq_vec[source,3] = 1b
    if ( (angle_scal eq angle_lin) and (angle_lin eq angle_par) ) then chisq_vec[source,3] = 2b

;    print, chisq_vec[source, *]

    if (im eq 0) then begin
        res = res_scal
        data_fit = data_fit_scal 
        best_fit_tag = 'Scalar'
    endif
    if (im eq 1) then begin
        res = res_lin
        data_fit = data_fit_lin 
        best_fit_tag = 'Linear'
    endif
    if (im eq 2) then begin
        res = res_par
        data_fit = data_fit_par 
        best_fit_tag = 'Parabolic'
    endif

    print, ' - :DP Best Fit: ', best_fit_tag

endif else begin
    print, 'Single fit performed...'
                           DATA_FIT = { ROBUST : robust, $
                           DATA : source_data, $
                           ERR_DATA : source_data_err, $
                           LIST_PIX : list_l_max, $
                           L_MAX : l_max, $
                           COB_MAX : cob_max, $
                           N_SIDE : n_side, $
                           NESTED : false, $
                           N_PAR : n_par }
             
              RES= mpfit( 'dp_funct_fit_hgauss2d', A_0, FUNCTARGS = data_fit, PARINFO  = parinfo, $
                          PERROR = perror, FASTNORM = FASTNORM,  AUTODERIVATIVE = AUTODERIVATIVE, $
                          MAXITER = MAXITER, BESTNORM = CHI2, $
                          COVAR = COVAR, QUIET = QUIET, /NOCATCH, DOF= DOF, $
                          ERRMSG= ERRMSG, STATUS = STATUS, NITER = 200 )

              if (scalar) then best_fit_tag = 'Scalar'
              if (linear) then best_fit_tag = 'Linear'
              if (parabolic) then best_fit_tag = 'Parabolic'

endelse
              PAR_FIT[ source, i_band, s, 0 : n_par - 1 ] = res ; last 2 are fiducial
              PAR_FIT_ERR[ source, i_band, s, * ] = [ perror, chi2, dof ]
; Analyzing the results

              if (map_infile eq undef) then begin
                  Case detect_tmp of
                      '30'  : CAL_0 = 0.023113417d9
                      '44'  : CAL_0 = 0.056961569d9
                      '70'  : CAL_0 = 0.13454397d9
                      '100' : CAL_0 = 0.241487d9
                      '143' : CAL_0 = 0.369344d9
                      '217' : CAL_0 = 0.481127d9
                      '353' : CAL_0 = 0.288331d9
                      '545' : CAL_0 = 0.0582363d9
                      '857' : CAL_0 = 0.00223400d9
                  Endcase
              endif
              FWHM_AVG = sqrt( res[ 2 ] * res[ 3 ] ) 
              FWHM_AVG_ARCMIN = fwhm_avg * 10800d0 / !dpi ; arcmin
              FWHM_AVG_ARCMIN_VEC[ source, i_band, s ] = fwhm_avg_arcmin
              ELIP = max( [ res[ 2 ], res[ 3 ] ] ) / min( [ res[ 2 ], res[ 3 ] ] )
              ELIP_VEC[ source, i_band, s ] = elip
              L_SOURCE_FIT = l_max + res[ 4 ]
              B_SOURCE_FIT = !dpi / 2d0 - cob_max - res[ 5 ]
              DELTA_L = 60d0 * 60d0 * ( l_object - 180d0 / !dpi * l_source_fit ) ; arcsec
              DELTA_L_VEC[ source, i_band, s ] = delta_l
              DELTA_B = 60d0 * 60d0 * ( b_object - 180d0 / !dpi * b_source_fit ) ; arcsec 
              DELTA_B_VEC[ source, i_band, s ] = delta_b
;print, 'before multiplying'
              FLUX_JY = !dpi / 4d0 / alog( 2d0 ) * res[ 1 ] * fwhm_avg^2d0 * cal_0 ; it is integrated, so no 1/sr
;stop
              FLUX_JY_VEC[ source, i_band, s ] = flux_jy 
; Computing relative errors
              ERR_REL_FLUX = perror[ 2 ] / res[ 2 ] + perror[ 3 ] / res[ 3 ] + perror[ 1 ] / res[ 1 ]
              ERR_REL_ELIP = perror[ 2 ] / res[ 2 ] + perror[ 3 ] / res[ 3 ]
              ERR_REL_FWHM = 0.5d0 * err_rel_elip
; Computing absolute values
              FLUX_JY_ERR_VEC[ source, i_band, s ] = err_rel_flux * flux_jy ; Jy
              FWHM_AVG_ARCMIN_ERR_VEC[ source, i_band, s ] = err_rel_fwhm * fwhm_avg_arcmin ; arcmin
              ELIP_ERR_VEC[ source, i_band, s ] = err_rel_elip * elip ; adimensional

;              If Keyword_set( PLOT ) then if ( PLOT eq 'test' ) then begin
              F_FIT = DP_FUNCT_FIT_HGAUSS2D( res, N_PAR = n_par, ROBUST =ROBUST, $
                                             DATA = DATA_FIT.data, ERR_DATA = data_fit.ERR_DATA, $
                                             LIST_PIX = data_fit.LIST_PIX, HGAUSS2D = HGAUSS2D, $
                                             L_MAX = data_fit.L_MAX, COB_MAX = data_fit.COB_MAX, $
                                             L_CENTR = L_CENTR, COB_CENTR = COB_CENTR, $
                                             N_SIDE = data_fit.N_SIDE, $
                                             src_model=src_map )

              TITLE_WINDOW = best_fit_tag

;              If (n_par eq 7)  then TITLE_WINDOW = 'SCALAR'
;              If (n_par eq 9)  then TITLE_WINDOW = 'LINEAR'
;              If (n_par eq 12) then TITLE_WINDOW = 'PARABOLIC'
;;               If NOT Keyword_set( JPEG ) then  window, n_sources * source + s, xsize = 800, ysize = 800, title = title_window, retain = 2 else  window,  source mod 31, xsize = 800, ysize = 800, title = title_window, retain = 2
; --- :DP
              help, data_fit.data, hgauss2d
;                  stop
;                  ps_map[list_l_max] = hgauss2d
              if ( (chisq_vec[source,2] lt 1.1) and (chisq_vec[source,3] eq 0) ) then begin
                  ps_map[list_l_max] = ps_map[list_l_max] + src_map*res[1] ; res[1] amplitude of the source

; --- correcting input map to remove better other sources
                  map_stokes[s, list_l_max] = map_stokes[s, list_l_max] - ps_map[list_l_max]
              endif else begin
                  ps_map[list_l_max] = -1.6375e30
              endelse
              

;                  gnomview, reform(map_stokes[s,*]), rot=[l_object, b_object, 0], grat=[1,1], chars=1.5, /hist, reso=0.6, tit='Clean Input Map', win=31
; ---
              if (plot eq 'test') then begin
;                  cc = map_tmp
;                  cc[i_pix_src] = -1.6375e30
;                  cc[i_pix_max] = -1.6375e30
;                  gnomview, cc, rot=[l_object, b_object, 0], grat=[1,1], chars=1.5, min=-500, max=500, reso=0.4, tit='Input Map', win=30

;                  cc = rms_tmp
;                  cc[i_pix_src] = -1.6375e30
;                  cc[i_pix_max] = -1.6375e30
;                  gnomview, cc, rot=[l_object, b_object, 0], grat=[1,1], chars=1.5, min=-50, max=50, reso=0.4, tit='Input RMS Map', win=27

                  cc = ps_map*0-1.6375e30
                  cc[list_l_max] = map_tmp[list_l_max]
                  cc[i_pix_src] = -1.6375e30
                  cc[i_pix_max] = -1.6375e30
                  gnomview, cc, rot=[l_object, b_object, 0], chars=1.5, grat=[1,1], min=-500, max=500, reso=0.4, tit='Input Map: Fit Region', win=31

                  cc = ps_map*0-1.6375e30
                  cc[list_l_max] = ps_map[list_l_max]
                  cc[i_pix_src] = -1.6375e30
                  cc[i_pix_max] = -1.6375e30
                  gnomview, cc, rot=[l_object, b_object, 0], chars=1.5, grat=[1,1], min=-500, max=500, reso=0.4, tit='PS Solution', win=29

;                  cc = map_tmp-ps_map
;                  cc[i_pix_src] = -1.6375e30
;                  cc[i_pix_max] = -1.6375e30
;                  gnomview, cc, rot=[l_object, b_object, 0], chars=1.5, grat=[1,1], min=-500, max=500, reso=0.4, tit='Difference', win=28

;                  cc = ps_map*0-1.6375e30
;                  cc[list_l_max] = 1.-ps_map[list_l_max]/map_tmp[list_l_max]
;                  gnomview, cc, rot=[l_object, b_object, 0], chars=1.5, grat=[1,1], min=-1.2, max=1.2, reso=0.4, tit='Relative Difference', win=27

                  cc = ps_map*0-1.6375e30
                  cc[list_l_max]=map_tmp[list_l_max]-ps_map[list_l_max]
                  cc[i_pix_max] = -1.6375e30
                  gnomview, cc, rot=[l_object, b_object, 0], chars=1.5, grat=[1,1], min=-500, max=500, reso=0.4, tit='Difference: Fit Region', win=26
; ---
                  If (JPEG eq undef) then  window, 0, xsize = 720*1.5, ysize = 450*1.5, title = title_window, retain = 2 else  window,  source mod 31, xsize = 800, ysize = 800, title = title_window, retain = 2
              
; --- :DP
                  print, ' - :DP create_ctb commented???'
;              create_ctb
                  !p.multi = [ 0, 2, 2 ]
                  END_PLOT = [ ' I ', ' Q ', ' U ', ' P ' ]
  ; 1st plot
                  plot, hgauss2d-src_map*res[1], /nodata, xtitle = ' SAMPLE ', ytitle = 'Fit-Source', $
                    charsize = 1.5, /xs, ys=2, title = object[ source ] + ' @ ' + detect_tmp + ': ' + $
                    end_plot[ s ]
                  oplot, (hgauss2d-src_map*res[1]), col = 30, psym = 3

                  legend, [ 'DATA', 'FIT' ], psym = [ 6, 6 ], col = [ 30, 104 ]

; --- :DP
                  plot, (data_fit.data-src_map*res[1]), /nodata, xtitle = 'SAMPLE', ytitle = 'Data-Fit', $
                    charsize = 1.5, /xs, ys=2, title = object[ source ] + ' @ ' + detect_tmp + ': ' + $
                    end_plot[ s ]
                  oplot, (data_fit.data-src_map*res[1]), col = 70, psym = 3
                  oplot, (hgauss2d-src_map*res[1]), col = 230, psym = 3
                  legend, [ 'DATA-PS', 'Background' ], psym = [ 7, 7 ], col = [ 70, 230 ]
; ---
  ; 2nd plot
                  L_CENTR_PLOT = ( l_centr - res[ 4 ] ) * 180d0 / !dpi
                  plot, l_centr_plot * 60d0, data_fit.data, /nodata, xtitle='ANGULAR DISTANCE (arcmin )', $
                    ytitle = ' T_THERMO ( !7l!8K CMB)', /xs, ys=2, charsize = 1.5;, yr=[-1000, 15000]
                  oplot, l_centr_plot * 60d0, data_fit.data, col = 0, psym = 3
                  oplot, l_centr_plot * 60d0, hgauss2d, psym = 3, col = 104 
                  oplot, l_centr_plot * 60d0, src_map*res[1], psym = 3, col = 75 
                  oplot, l_centr_plot * 60d0, src_map*0.+res[0], psym = 3, col = 245 
                  legend, [ 'DATA', 'FIT', 'PS map', 'Offset' ], psym = [ 7, 7, 7, 7 ], col = [ 0, 104, 75, 245 ]
                  errplot, l_centr_plot * 60d0, data_fit.data-data_fit.err_data, data_fit.data+data_fit.err_data
  ; 3rd plot
                  B_CENTR_PLOT = ( cob_centr - res[ 5 ] ) * 180d0 / !dpi
                  plot, b_centr_plot * 60d0, data_fit.data, /nodata, xtitle='ANGULAR DISTANCE (arcmin )', $
                    ytitle = ' T_THERMO ( !7l!8K CMB)', /xs, ys=2, charsize = 1.5;, yr=[-1000, 15000]
                  oplot, b_centr_plot * 60d0, data_fit.data, col = 0, psym = 7
                  oplot, b_centr_plot * 60d0, hgauss2d, psym = 7, col = 104
                  oplot, b_centr_plot * 60d0, src_map*res[1], psym = 7, col = 75 
                  oplot, b_centr_plot * 60d0, src_map*0.+res[0], psym = 7, col = 245 
                  legend, [ 'DATA', 'FIT', 'PS map', 'Offset' ], psym = [ 7, 7, 7, 7 ], col = [ 0, 104, 75, 245 ]

;                  IF (JPEG ne undef) then begin
;                      PATH_JPEG = read_path() +'polar/sources/hg2d/jpeg/'
;                      TAIL_JPEG = [ '_I', '_Q', '_U','_P' ]
;                      FILENAME_JPEG = 'hg2d_fit_' + object[ source ] + '_' + detect_tmp + tail_jpeg[ s ]  + '.jpeg'
;                      write_jpeg, path_jpeg + filename_jpeg, TVRD(/TRUE), /TRUE
;                  Endif
              Endif

;endif

; Calibrating in Jy
; Official values from http://wiki.planck.fr/index.php/Proc/DataStatus
; Writing them in Jy/sr/K_thermo
              L_PIX = sqrt( 4d0 * !dpi / 12d0 / n_side^2d0 )
; --- :DP
;(MJy/sr IRAS)/KCMB
; GHz              Uc             UcErr
; 100        243.86045      0.19353327
; 143        371.60990      0.054413923
; 217        483.09579      0.0070682290
; 353        287.41872      0.0032819714
; 545        58.011943      0.0020943308
; 857        2.2383275      0.00057251453

              if (map_infile eq undef) then begin
                  Case detect_tmp of
;                  '30'  : CAL_0 = 1d0 ; to be updated
;                  '44'  : CAL_0 = 1d0 ; to be updated
;                  '70'  : CAL_0 = 1d0 ; to be updated
                      '30'  : CAL_0 = 0.023113417d9
                      '44'  : CAL_0 = 0.056961569d9
                      '70'  : CAL_0 = 0.13454397d9
                      '100' : CAL_0 = 0.241487d9
                      '143' : CAL_0 = 0.369344d9
                      '217' : CAL_0 = 0.481127d9
                      '353' : CAL_0 = 0.288331d9
                      '545' : CAL_0 = 0.0582363d9
                      '857' : CAL_0 = 0.00223400d9
                  Endcase
              endif
              SR_PIX = l_pix^2d0
              K_2_JY = cal_0 * sr_pix
;              If s eq 0 then print,'UNIT CONVERSION FROM K_THERMO TO Jy: ', str( k_2_jy)
; Integrated flux
;;INT_FLUX[source, i_band, s ] = res_fit_min[ 0 ] * k_2_jy * sign_tmp
;;INT_FLUX_ERR[ source, i_band, s ] = pcerror[ 0 ] * k_2_jy * sign_tmp ; res_fit_err[ 0 ] * k_2_jy * sign_tmp
; Expected value
;;INT_FLUX_TOT[ source, i_band, s ] = p_accum_cnt * k_2_jy * sign_tmp 
              print,'BAND: ', band[ q_band ]
          Endfor                ; Stokes (I, Q, U )
;          print,'Returning the integrated flux (T, Q, U & P) in Jy ...'
          if (strlen(strtrim(catalog_infile)) ge 1) then begin
;              print, ' - Source Flux = ', strtrim(string(flux_jy_vec[ source, i_band, * ]*1.e-3),2) + ' mJy vs PwS Flux = ' + strtrim(string(flux[source]),2)+' +/- '+strtrim(string(fluxerr[source]),2) + ' --> ' + strtrim(string(flux_jy_vec[ source, i_band, * ]*1.e-3/flux[source]),2)
              chisq_vec[source,4] = flux_jy_vec[source,i_band,0]*1.e-3
              chisq_vec[source,5] = flux[source]
              print, reform( chisq_vec[source,*] )
          endif else begin
              print, ' - Source Flux = ', strtrim(string(flux_jy_vec[ source, i_band, * ]*1.e-3),2) + ' mJy'
          endelse
NO_OBSERVED:
          if q_max[ 0 ] eq -1 then print,'Source not observed'
; --- :DP
          print, ' --- 0 --- 0 --- 0 --- 0 --- 0 --- 0 --- '
if (do_stop) then begin
    tres = res
    tres[2] = tres[2] * 60. * 180. / !pi
    tres[3] = tres[3] * 60. * 180. / !pi
    tres[4] = tres[4] * 60. * 180. / !pi
    tres[5] = tres[5] * 60. * 180. / !pi
    tres[6] = tres[6] * 180. / !pi
    print, ' Resulting parameters:'
    print, ' background offset: ', tRES[0]
    print, ' Gauss Amplitude: ', tRES[1]
    print, ' y-fwhm arcmin: ', tRES[2]
    print, ' x-fwhm arcmin: ', tRES[3]
    print, ' peak center x: ', tRES[4]
    print, ' peak center y: ', tRES[5]
    print, ' Rotation angle: ', tRES[6]
    if (n_elements(tRES) eq 7) then print, 'Scalar'
    if (n_elements(tRES) eq 9) then print, 'Linear'
    if (n_elements(tRES) eq 12) then print, 'Parabolic'
    print, ' Chisq = ', chi2, ' for ', dof, ' dof --> red_chi2 = ', chi2 / (dof-n_par)
    print, ' --- Source Loop ---'
    stop
endif

      Endfor                    ;  sources
      print, 'Saving map, Source ',source,'/',n_sources-1
      write_fits_map, tag+'_ps_map.fits', ps_map, /ring, units='!7l!8K CMB'
;      mollview, tag+'_ps_map.fits', win=25, max=750, chars=1.5
      cl2fits, chisq_vec, tag+'_chisq_vec.fits'

      bpix = where(ps_map eq -1.6375e30)
      clean = map_fits-ps_map
      clean[bpix] = -1.6375e30
      write_fits_map, tag+'_ps-removed_map.fits', clean, /ring, units='!7l!8K CMB'

      mollview, tag+'_ps-removed_map.fits', chars=1.5, win=1

  Endfor                        ; bands

;RESULT = [ [int_flux], [int_flux_err] ]
  If Keyword_set( R_FIX ) then RESULT = [ [int_flux], [int_flux_tot] , [int_flux_err] ] 


;; Writing the results
;HEADER = [ 'Coordinates: Galatic', ' Flux units: Jy' ] 
;DATA = {BANDS: band, $
;        L_COORD: all_l_object, $
;        B_COORD: all_b_object, $
;        FLUX_I: result[ *, *, 0 ], $
;        FLUX_Q: result[ *, *, 1 ], $
;        FLUX_U: result[ *, *, 2 ], $
;        FLUX_P: result[ *, *, 3 ]}

;  If Keyword_set( SAVE ) then begin
  If ( SAVE ) then begin
      if (false) then begin
          PATH = read_path()
          TAIL_BG = ''
          If ( ROBUST ) then TAIL_RB = '_robust' else TAIL_RB = ''
          If Keyword_set( LINEAR ) then TAIL_BG = '_linear'
          If Keyword_set( PARABOLIC ) then TAIL_BG = '_parabolic'
          If Keyword_set( SURVEY ) then TAIL_SURVEY = '_survey' + str( survey ) else TAIL_SURVEY = ''
          END_TAIL_HALF = [ 'first', 'second' ]
          If Keyword_set( HALF ) then TAIL_HALF = '_' + end_tail_half[ half -1 ] + '_half' else TAIL_HALF = ''
          If NOT Keyword_set( FILE ) then FILE_HG2D = path + 'polar/sources/hg2d/h2gd_test' +  tail_survey + $
            tail_half + tail_bg + tail_rb + '.sav'
;          FILE_HG2D = path + file + tail_survey + $
;            tail_half + tail_bg + tail_rb + '.sav' else $
          If Keyword_set( DETECT ) AND Keyword_set( CATALOG )  then FILE_HG2D = path + 'polar/sources/hg2d/h2gd_' + catalog + '_' + detect + tail_survey + tail_bg + tail_rb + '.sav'
          OBJECT_INFO = ''
          If Keyword_set( COORD ) then OBJECT_INFO = coord
          If Keyword_set( OBJECT ) then OBJECT_INFO = object
          If Keyword_set( CATALOG ) then OBJECT_INFO = name_source
      endif
      FILE_HG2D = outfile
      save, cob_max_vec, delta_b_vec, delta_l_vec, elip_vec, flux_jy_vec, fwhm_avg_arcmin_vec, l_max_vec, $
        object_info, par_fit, par_fit_err, flux_jy_err_vec, fwhm_avg_arcmin_err_vec, elip_err_vec, $
        FILENAME = file_hg2d
      print,'SAVED FILE FROM HG2FLUX FITIN: ', file_hg2d
  Endif

  If n_elements( res ) gt 0 then begin
      print,'INITIAL PARAMS: ', str( a_0 )
      print,'FINAL PARAMS: ', str( res )
      print,'FWHM (arcmin): ', str( [ res[ 2 ], res[ 3 ] ] * 10800d0 / !dpi )
      L_SOURCE_FIT = l_max + res[ 4 ]
      B_SOURCE_FIT = !dpi / 2d0 - cob_max - res[ 5 ]
      print,'POSITION OF THE SOURCE. (L, B): ' + str( l_source_fit * 180d0 / !dpi ) + ', ' + $
        str( b_source_fit * 180d0 / !dpi )
      print,'POSITION IN J2000-NED. (L, B): ' + str( l_object ) + ', ' + str( b_object )
      print,'DIFFERENCE IN ARCMIN. (L, B): ' + str( ( l_source_fit * 180d0 / !dpi - l_object ) * 60d0 ) + ', ' + $
        str( ( b_source_fit * 180d0 / !dpi - b_object ) * 60d0 )
  ; Writing them in Jy/sr/K_thermo
      SIZE_PIX = sqrt( 4d0 * !dpi / 12d0 / n_side^2d0 ) ; rad
      if (map_infile eq undef) then begin
          Case detect_tmp of
              '30'  : CAL_0 = 0.023113417d9
              '44'  : CAL_0 = 0.056961569d9
              '70'  : CAL_0 = 0.13454397d9
              '100' : CAL_0 = 0.241487d9
              '143' : CAL_0 = 0.369344d9
              '217' : CAL_0 = 0.481127d9  
              '353' : CAL_0 = 0.288331d9
              '545' : CAL_0 = 0.0582363d9
              '857' : CAL_0 = 0.00223400d9
          Endcase
      endif
      SR_PIX = size_pix^2d0
      K_2_JY = cal_0 * sr_pix
  ;print,' FLUX (HEALPix GAUSSIAN2D): ', str( !dpi * res[ 1 ] * res[ 2 ] * res[ 3 ] / (8d0 * alog( 2d0 ) ) * k_2_jy )
  ;PCERROR = perror * sqrt( chi2 / dof ) 
  ;print,'ERR FINAL PARAMS: ', str( pcerror )
      F_INI = DP_FUNCT_FIT_HGAUSS2D( A_0, ROBUST = ROBUST, N_PAR = n_par, $
                                     DATA = DATA_FIT.data, ERR_DATA = data_fit.ERR_DATA, $
                                     LIST_PIX = data_fit.LIST_PIX, HGAUSS2D = HGAUSS2D, $
                                     L_MAX = data_fit.L_MAX, COB_MAX = data_fit.COB_MAX, $
                                     N_SIDE = data_fit.N_SIDE)
      
      F_FIT = DP_FUNCT_FIT_HGAUSS2D( res, ROBUST =ROBUST, N_PAR = n_par, $
                                     DATA = DATA_FIT.data, ERR_DATA = data_fit.ERR_DATA, $
                                     LIST_PIX = data_fit.LIST_PIX, HGAUSS2D = HGAUSS2D, $
                                     L_MAX = data_fit.L_MAX, COB_MAX = data_fit.COB_MAX, $
                                     N_SIDE = data_fit.N_SIDE )
      
  Endif
  print,'END'
;If Keyword_set( DO_STOP ) then stop
  If n_elements( res ) gt 0 then return, [ [par_fit ], [par_fit_err ] ] else return,'NOTHING TO FIT'

END
