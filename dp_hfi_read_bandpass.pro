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

FUNCTION dp_hfi_read_bandpass, VERSION, IN_PATH =  IN_PATH, IMO = IMO, RIMO = RIMO, PATH_RIMO = PATH_RIMO

; -- Definitions
   CONST = { C: 2.9979246d+08, H: 6.6260755d-34, HBAR: 1.0545727d-34, K: 1.3806580d-23, TCMB: 2.725d0 }

   N_DETECT = 52                ; Planck HFI
   Case VERSION of
       'v101': N_MEASUREMENTS =  17000l
       'v201': N_MEASUREMENTS = 17000l
       'dx8': N_MEASUREMENTS = 17000l
       'ffp6': N_MEASUREMENTS = 17000l
   Endcase

   BP_INFO = replicate( { NAME: '', FREQ: dblarr( n_measurements), TRANS: dblarr( n_measurements) }, n_detect ) 
   If Keyword_set( RIMO ) then begin
       If NOT Keyword_set( PATH_RIMO ) then PATH_RIMO =''
       Case VERSION of
           'v101': BEGIN
               RIMO_FITS = 'HFI-RIMO-20101025_2_55.fits'
               fitsshift = 2
           END
  
           'v201': BEGIN
               RIMO_FITS = 'HFI-RIMO-20101025_2_57.fits'
               fitsshift = 3
           END 
           'dx8' : BEGIN
               print, ' - :DP DX8'
               RIMO_FITS = '/project/projectdirs/planck/data/mission/DPC_maps/dx8/hfi/HFI-RIMO-20120124.fits'
               fitsshift = 3
           END
           'ffp6' : BEGIN
               print, ' - :DP FFP6'
               RIMO_FITS = '/project/projectdirs/planck/data/ffp6/mission/HFI_RIMO_FFP6-v1.fits'
               fitsshift = 3
           END
       Endcase
       FITS_RIMO = mrdfits( path_rimo + rimo_fits, 1, HDR_RIMO , /silent)
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
       For i = 0, 51 do begin
           NAME_RIMO_TMP = strmid( fits_rimo[ i ].detector, 0, 3 ) + '_' + strmid( fits_rimo[ i ].detector, 4  )
           Q_NAME = where( strmid( bolo_id_tmp, 3) eq name_rimo_tmp )
           BP_INFO[ i ].NAME = bolo_id_tmp[ q_name[ 0 ] ] ; completing the name of the detector (e.g. 100-1a in the RIMO)
       Endfor
  ; Joining the information of the bandpasses in the bp_info structure
       For i = 0, n_detect - 1 do begin
      
           FITS_RIMO = mrdfits( path_rimo + rimo_fits, fitsshift + i ,/silent)
           
           BP_INFO[ i ].FREQ  = fits_rimo.wavenumber * 1d2 * const.c ; wave number given in cm^-1
           BP_INFO[ i ].TRANS = fits_rimo.transmission ; values of the transmission (normalized later)
           
       Endfor
   ENDIF 
  
; IMO 

  If Keyword_set( IMO ) then begin
    Case VERSION of
    'v101': IMO_LBL = '2_55'
    'v201': IMO_LBL = '2_57'
    Endcase
  ImoFileMD = '/data/dmc/MISS01/METADATA'
  ImoFile = ImoFileMD + '%lbl:IMO_' + imo_lbl
  ; Open Imo
  IMO_GROUP = PIOOpenIMOFile(ImoFile, 'r')
  print,'IMO_GROUP: ',imo_group
    IF imo_group LE 0 THEN pioerrmess, imo_group
  ; Joining the information of the bandpasses in the bp_info structure
  G_TYPE = 'PIOSTRING'
  TAIL_X = 'SpectralResp:SpecTransmissions:VectX'
  CHANNEL_TMP = [ '100', '143', '217', '353', '545', '857' ]
  N_DETECT = 0
    For i = 0, 5 do begin
    BOLO_LIST = get_hfibolo_list( CHANNEL = channel_tmp[ i ] )
    N_LIST = n_elements( bolo_list )
      For j = 0, n_list - 1 do begin
      BP_INFO[ n_detect + j ].NAME = bolo_list[ j ]
      G_NAME_X = 'IMO:HFI:DET:Phot_Pixel Name="' + bolo_list[ j ]  + '":' + tail_x
      DUMMY  = PIOGetValue( G_VALUE, gerror, g_type, gunit, gcomment, $
                            g_name_x, imo_group)
      Q = strpos( g_value, ':' , /reverse_s)
      G_VALUE_READ = strmid( g_value, 0, q )
      WN = pioread( g_value_read ) ; Wave number in cm^-1
      BP_INFO[ n_detect + j ].FREQ = 2.9979246e+10 * wn ; freq in Hz
      TAIL_Y = 'SpectralResp:SpecTransmissions:VectY'
      G_NAME_Y = 'IMO:HFI:DET:Phot_Pixel Name="' + bolo_list[ j ]  + '":' + tail_y
      DUMMY  = PIOGetValue( G_VALUE, gerror, g_type, gunit, gcomment, $
                            g_name_y, imo_group)
      Q = strpos( g_value, ':' , /reverse_s)
      G_VALUE_READ = strmid( g_value, 0, q )
      BP_INFO[ n_detect + j ].TRANS = pioread( g_value_read )
      Endfor
    N_DETECT = n_detect + n_list
    Endfor
  Endif 


  return, bp_info
END
