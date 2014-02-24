; Updated removing features: only RIMO supported
; Adapted from hfi one
;+
; NAME:
;
;   lfi_read_bandpass
;
; PURPOSE:
;
;   Read LFI bandpass info from official .save files
;
; CATEGORY:
;

FUNCTION lfi_read_bandpass, rimofile, bp_info=bp_info

; -- Definitions
   CONST = { C: 2.9979246d+08, H: 6.6260755d-34, HBAR: 1.0545727d-34, K: 1.3806580d-23, TCMB: 2.725d0 }

   N_DETECT = 22                ; Planck LFI

   N_MEASUREMENTS = 500l

   N_FREQ = 3

   BP_INFO = replicate( { NAME: '', FREQ: dblarr( n_measurements), TRANS: dblarr( n_measurements) }, n_detect ) 
   BP_FREQ = replicate( { NAME: '', FREQ: dblarr( n_measurements), TRANS: dblarr( n_measurements), CFREQ: 0.d0 }, n_freq ) 

;   print, ' - :DP FFP6'
   RIMO_FITS = rimofile ;'/project/projectdirs/planck/data/ffp6/mission/HFI_RIMO_FFP6-v1.fits'
   fitsshift = 3

   FITS_RIMO = mrdfits( rimo_fits, 1, HDR_RIMO , /silent)

   If strpos( hdr_rimo[ 4 ], '22' ) eq -1 then stop,'Is it an LFI RIMO file?'

   BOLO_ID_TMP = []
   for i=0,n_detect-1 do bolo_id_tmp = [BOLO_ID_TMP, fits_rimo[i].(0)]

   For i = 0, n_detect-1 do begin
       NAME_RIMO_TMP = strmid( fits_rimo[ i ].detector, 0, 3 ) + '_' + strmid( fits_rimo[ i ].detector, 4  )
       Q_NAME = where( strmid( bolo_id_tmp, 3) eq name_rimo_tmp )
       BP_INFO[ i ].NAME = bolo_id_tmp[ q_name[ 0 ] ] ; completing the name of the detector (e.g. 100-1a in the RIMO)
   Endfor
  ; Joining the information of the bandpasses in the bp_info structure
   For i = 0, n_detect - 1 do begin
      
       FITS_RIMO = mrdfits( rimo_fits, fitsshift + i ,/silent)

       f = fits_rimo.(0) * 1.d9
       b = fits_rimo.(1)
       cnt = 0

       ff = []
       bb = []
       for j=0,n_elements(f)-2 do begin
           if (f[j] ne f[j+1]) then begin
               ff = [ ff,f[j+1] ]
               bb = [ bb,b[j+1] ]
               cnt = cnt + 1
           endif
       endfor

       nf = n_elements(ff)

       BP_INFO[ i ].FREQ  = ff ;fits_rimo.wavenumber * 1d2 * const.c ; wave number given in cm^-1
       BP_INFO[ i ].TRANS = bb ;fits_rimo.transmission ; values of the transmission (normalized later)
           
   Endfor
  
   For i = 0, n_freq - 1 do begin
      
       FITS_RIMO = mrdfits( rimo_fits, fitsshift + n_detect + i ,/silent)
           

       f = fits_rimo.(0) * 1.d9
       b = fits_rimo.(1)
       cnt = 0

       ff = []
       bb = []
       for j=0,n_elements(f)-2 do begin
           if (f[j] ne f[j+1]) then begin
               ff = [ ff,f[j+1] ]
               bb = [ bb,b[j+1] ]
               cnt = cnt + 1
           endif
       endfor

       nf = n_elements(ff)

       BP_FREQ[ i ].FREQ  = ff ;fits_rimo.wavenumber * 1d2 * const.c ; wave number given in cm^-1
       BP_FREQ[ i ].TRANS = bb ;fits_rimo.transmission ; values of the transmission (normalized later)
       gf = where(BP_FREQ[ i ].FREQ gt 0.)
       BP_FREQ[ i ].CFREQ = int_tabulated(BP_FREQ[ i ].FREQ[gf], BP_FREQ[ i ].FREQ[gf]*BP_FREQ[ i ].TRANS[gf], /double) / int_tabulated(BP_FREQ[ i ].FREQ[gf], BP_FREQ[ i ].TRANS[gf], /double)
       
   Endfor

   print, ' HFI central frequencies: ', BP_FREQ.CFREQ * 1.d-9

  return, bp_freq
END
