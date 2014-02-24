function rotate_map, Nside, map, angle, nested=nested, radians=radians

   False = 0b
   True  = 1b
   if (not keyword_set(nested)) then nested = False
   if (not keyword_set(radians)) then radians = False

   if (nested) then mapr = reorder(map, in='NESTED', out='RING') else mapr = map
   if (radians) then angle = angle/!radeg

   sz = size( map )
   Npix = sz[1]
   Ns = npix2nside(sz[1])

   if (Ns ne Nside) then stop,'Nside /= Ns'

   ipix = lindgen(npix)
   pix2vec_ring, Nside, ipix, pixvec
;       angle = [1.5,1.5,1.5]
   matrix = euler_matrix_new(angle[0],angle[1],angle[2], deg=True)
   rotvec = rotate_coord(pixvec, euler_matrix=matrix)

   vec2pix_ring, Nside, rotvec, rotpix
   
   rot_map = fltarr( Npix )
   rot_map = mapr[rotpix]

   if (nested) then begin
       rot_mapn = reorder(rot_map, in='ring', out='nest')
       rot_map = rot_mapn
   endif
return, rot_map

END
