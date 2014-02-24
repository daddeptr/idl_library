PRO OBJECT_LIST, OBJECT = OBJECT, L_OBJECT, B_OBJECT, LIST = LIST, ECLIPTIC = ECLIPTIC, M3 = M3

; A way to have main regions together
; List will give the list of objects. ASCII file!
FILE_LIST_OBJECT = read_path( /bp, M3 = M3 ) + 'sources/list/list_sources_galactic.txt'
readcol, file_list_object, NAME_SOURCE, A, B, FORMAT = 'A,D,D', /silent
  If Keyword_set( LIST ) then begin
  L_OBJECT = NAME_SOURCE ; Simply call object_list,/list,LIST_OBJECT
  GOTO, END_Q2
  Endif
; Loop over different objects
N_SOURCES = n_elements( object )
L_OBJECT = strarr( n_sources )
B_OBJECT = l_object
  For s = 0, n_sources - 1 do begin
  OBJECT_TMP = object[ s ]
; For planets rings are read directly with objet2ring
NAME_PLANET = [ 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'NO_OBJECT' ]
Q_PLANET = where( object_tmp eq name_planet )
  If q_planet[ 0 ] ne -1 then begin
  L_OBJECT[ s ] = 'planet' ; label for read_data to select the BigPlanets group
  GOTO, END_Q
  Endif
Q_SOURCE_3 = where( object_tmp eq name_source )
  If q_source_3[ 0 ] eq -1 then begin
  print,'OBJECT_LIST. Source: ' + object_tmp + ' not recognized.' 
  L_OBJECT[ s ] = 'NO_OBJECT'
  GOTO, END_Q
  Endif else  print,'OBJECT_LIST: Considering ', object_tmp
; Coordinates from SIMBAD J2000 (LMC: 280.4652, -32.8884)
; Perseus: Feature D in Tibbs et al. 2010: 03:42:48, +32:39:50 (RA, DEC J2000)
; RA: 55.7000 deg, DEC = 32.663889 deg
; LON = 159.85324 deg, LAT = -17.654644 deg
L_OBJECT[ s ] = ( a[ q_source_3 ] )[ 0 ]
B_OBJECT[ s ] = ( b[ q_source_3 ] )[ 0 ]
TAIL_COORD = 'GALACTIC'
  If Keyword_set( ECLIPTIC ) then begin
  TAIL_COORD = 'ECLIPTIC'
  euler, l_object[ s ], b_object[ s ], l_object_ecliptic, b_object_ecliptic, 6 
  L_OBJECT[ s ] = l_object_ecliptic
  B_OBJECT[ s ] = b_object_ecliptic
  Endif
print, tail_coord + ' COORDINATES: '
print, 'L: ', str( l_object[ s ] )
print, 'B: ', str( b_object[ s ] )
END_Q:
  Endfor ; s (different objects)
; If it is a scalar, force it
  If n_elements( l_object ) eq 1 then begin
  L_OBJECT = l_object[ 0 ]
  B_OBJECT = b_object[ 0 ]
  Endif
END_Q2:
END
