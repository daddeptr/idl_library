function fill_map_gaps, map, bad_val=bad_val, radius_arcmin=radius_arcmin, nest=nest, silent=silent
   if not keyword_set(bad_val) then bad_val = -1.6375e30
   if not keyword_set(radius_arcmin) then radius_arcmin = 60.

   if not keyword_set(nest) then mapn=reorder(map, in='ring', out='nest') else mapn = map

   if not keyword_set(silent) then print, ' - fill_map_gaps: Filling using ', string(radius_arcmin, format='(f5.2)')+' arcmin radius.'

   map = mapn
   mapn = 0.

   bp = where(map[*,0] eq bad_val)
   gp = where(map[*,0] ne bad_val)

   sz = size(map)
;   print, sz

   n_bp = n_elements(bp)

   if (bp[0] ne -1) then begin
       if not keyword_set(silent) then print, ' - fill_map_gaps: Bad pixels found...', bad_val
       if not keyword_set(silent) then help, bp
       delta = n_bp / 10
       for ipix=0l,n_bp-1 do begin

           if ( (not keyword_set(silent)) and ( (ipix/delta)*delta eq ipix) ) then print, ' - fill_map_gaps: ......'+string( float(ipix)/n_bp*100., format='(f5.1)')+'%'
           pix2vec_ring, npix2nside(sz[1]), bp[ipix], vec0
           query_disc, npix2nside(sz[1]), vec0, radius_arcmin/60./!radeg, listp, /nest
           
           gp = where(map[listp,0] ne bad_val)
           
           listp = listp[gp]
           
           map[bp[ipix],0] = median(map[listp,0])
       endfor
   endif

   mapr = reorder(map, in='nest', out='ring')
   return, mapr

end
