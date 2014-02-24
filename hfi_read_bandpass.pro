; Updated removing features: only RIMO supported
;+
; NAME:
;
;   hfi_read_bandpass
;
; PURPOSE:
;
;   Read HFI bandpass info from official .save files
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;   bp = hfi_read_bandpass(version) + Keywords 
;
; INPUTS:
;
;   version :  HFI bandpass version 
;              v101
;              v201
;
; OPTIONAL INPUTS:
;
;   None
;
; KEYWORD PARAMETERS:
;  
;   IN_PATH:  input path for save datasets. Default same as scripts (change it for CVS).
;   IMO:  If an IMO label is to be read. By default, it is 2_55 for v101 and 2_57 for v201
;   RIMO: If a RIMO fits file is to read
; 
; OUTPUTS:
;
;   
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;    
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;    Planck collaboration only.
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;   bp = hfi_read_bandpass('v101', /imo)
;   bp = hfi_read_bandpass('v101', /rimo)
;
; MODIFICATION HISTORY:
;
;   Created by S. Hildebrandt and J.F. Macias-Perez, Oct, 2010
;-

FUNCTION hfi_read_bandpass, rimofile, bp_info=bp_info

; -- Definitions
   CONST = { C: 2.9979246d+08, H: 6.6260755d-34, HBAR: 1.0545727d-34, K: 1.3806580d-23, TCMB: 2.725d0 }

   N_DETECT = 52                ; Planck HFI

   N_MEASUREMENTS = 17000l

   N_FREQ = 6

   BP_INFO = replicate( { NAME: '', FREQ: dblarr( n_measurements), TRANS: dblarr( n_measurements), CFREQ: 0.d0 }, n_detect ) 
   BP_FREQ = replicate( { NAME: '', FREQ: dblarr( n_measurements), TRANS: dblarr( n_measurements), CFREQ: 0.d0 }, n_freq ) 

;   print, ' - :DP FFP6'
   RIMO_FITS = rimofile ;'/project/projectdirs/planck/data/ffp6/mission/HFI_RIMO_FFP6-v1.fits'
   fitsshift = 3

   FITS_RIMO = mrdfits( rimo_fits, 1, HDR_RIMO , /silent)
   If strpos( hdr_rimo[ 4 ], '52' ) eq -1 then stop,'Is it an HFI RIMO file?'
   BOLO_ID_TMP = [ '00_100_1a', '01_100_1b', '20_100_2a', '21_100_2b', '40_100_3a', '41_100_3b', $            
                       '80_100_4a', '81_100_4b', '02_143_1a', '03_143_1b', '30_143_2a', '31_143_2b', $
                       '50_143_3a', '51_143_3b', '82_143_4a', '83_143_4b', '10_143_5',  '42_143_6',  $            
                       '60_143_7',  '70_143_8',  '11_217_5a', '12_217_5b', '43_217_6a', '44_217_6b', $
                       '61_217_7a', '62_217_7b', '71_217_8a', '72_217_8b', '04_217_1',  '22_217_2',  $
                       '52_217_3',  '84_217_4',  '23_353_3a', '24_353_3b', '32_353_4a', '33_353_4b', $
                       '53_353_5a', '54_353_5b', '63_353_6a', '64_353_6b', '05_353_1',  '13_353_2',  $
                       '45_353_7',  '85_353_8',  '14_545_1',  '34_545_2',  '55_545_3',  '73_545_4',  $
                       '25_857_1',  '35_857_2',  '65_857_3',  '74_857_4' ]
   For i = 0, n_detect-1 do begin
       NAME_RIMO_TMP = strmid( fits_rimo[ i ].detector, 0, 3 ) + '_' + strmid( fits_rimo[ i ].detector, 4  )
       Q_NAME = where( strmid( bolo_id_tmp, 3) eq name_rimo_tmp )
       BP_INFO[ i ].NAME = bolo_id_tmp[ q_name[ 0 ] ] ; completing the name of the detector (e.g. 100-1a in the RIMO)
   Endfor
  ; Joining the information of the bandpasses in the bp_info structure
   For i = 0, n_detect - 1 do begin
      
       FITS_RIMO = mrdfits( rimo_fits, fitsshift + i ,/silent)

       f = fits_rimo.wavenumber * 1d2 * const.c
       b = fits_rimo.transmission
       cnt = 0

;help, where(f gt 0.)

;       for j=0,n_elements(f)-2 do begin
;           if (f[j] eq f[j+1]) then print, j 
;       endfor
;stop
       ff = f[where(f gt 0.)];[]
       bb = b[where(f gt 0.)];[]
;       for j=0,n_elements(f)-2 do begin
;           if (f[j] ne f[j+1]) then begin
;               ff = [ ff,f[j+1] ]
;               bb = [ bb,b[j+1] ]
;               cnt = cnt + 1
;           endif
;       endfor

;       nf = n_elements(ff)

       BP_INFO[ i ].FREQ  = ff ;fits_rimo.wavenumber * 1d2 * const.c ; wave number given in cm^-1
       BP_INFO[ i ].TRANS = bb ;fits_rimo.transmission ; values of the transmission (normalized later)
       gf = where(BP_INFO[ i ].FREQ gt 0.)
;help, gf
;stop
;##       BP_INFO[ i ].CFREQ = int_tabulated(BP_INFO[ i ].FREQ[gf], BP_INFO[ i ].FREQ[gf]*BP_INFO[ i ].TRANS[gf], /double) / int_tabulated(BP_INFO[ i ].FREQ[gf], BP_INFO[ i ].TRANS[gf], /double)
       BP_INFO[ i ].CFREQ = int_tabulated(ff, ff*bb, /double) / int_tabulated(ff, bb, /double)
           
   Endfor
  
   For i = 0, n_freq - 1 do begin
      
       FITS_RIMO = mrdfits( rimo_fits, fitsshift + n_detect + i ,/silent)
           

       f = fits_rimo.wavenumber * 1d2 * const.c
       b = fits_rimo.transmission
       cnt = 0

       ff = f[where(f gt 0.)];[]
       bb = b[where(f gt 0.)];[]
;       for j=0,n_elements(f)-2 do begin
;           if (f[j] ne f[j+1]) then begin
;               ff = [ ff,f[j+1] ]
;               bb = [ bb,b[j+1] ]
;               cnt = cnt + 1
;           endif
;       endfor

;       nf = n_elements(ff)

       BP_FREQ[ i ].FREQ  = ff ;fits_rimo.wavenumber * 1d2 * const.c ; wave number given in cm^-1
       BP_FREQ[ i ].TRANS = bb ;fits_rimo.transmission ; values of the transmission (normalized later)
       gf = where(BP_FREQ[ i ].FREQ gt 0.)
;##       BP_FREQ[ i ].CFREQ = int_tabulated(BP_FREQ[ i ].FREQ[gf], BP_FREQ[ i ].FREQ[gf]*BP_FREQ[ i ].TRANS[gf], /double) / int_tabulated(BP_FREQ[ i ].FREQ[gf], BP_FREQ[ i ].TRANS[gf], /double)
       BP_FREQ[ i ].CFREQ = int_tabulated(ff, ff*bb, /double) / int_tabulated(ff, bb, /double)
       
   Endfor

   print, ' HFI central frequencies: ', BP_FREQ.CFREQ * 1.d-9

   if (keyword_set(BP_INFO)) then print, ' HFI detector central frequencies', BP_INFO.CFREQ * 1.d-9

  return, bp_freq
END
