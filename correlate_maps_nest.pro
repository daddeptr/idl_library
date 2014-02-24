function correlate_maps_nest, map1_in, map2_in, correlation_nside, mask=mask

   True = 1b
   False = 0b

   help, map1_in, map2_in

   sz1 = size(map1_in)
   ns1 = npix2nside(sz1[1])

   sz2 = size(map2_in)
   ns2 = npix2nside(sz2[1])

;print, ns1, ns2
;stop

    if (ns1 ne ns2) then stop, 'Map size inconsistent: ', ns1, ' vs. ', ns2 else Nside=ns1

    Npix = 12l*Nside^2
    if not keyword_set( mask ) then gmask = fltarr(Npix)+1.
    bp = where( (gmask eq 0.) or (map1_in eq -1.6375e30) or (map2_in eq -1.6375e30) or (finite(map1_in, /nan) eq True) or (finite(map2_in, /nan) eq True) )

    low_Nside = correlation_nside
    low_Npix = 12l*low_Nside^2

    correlation_maps = fltarr(low_Npix)
    correlation_maps[*] = -1.6375e30

    pix2vec_nest, low_Nside, lindgen(low_Npix), vec, vertex

    for i=0, low_Npix-1 do begin
        if ( (i/500)*500 eq i) then print, i, ' / ', low_Npix, ' --> '+string(float(i)/low_Npix*100,format='(f4.1)')+'%'

        query_polygon, Nside, transpose( reform( vertex[i,*,*]) ), pix_list, inclusive=True, /nested

        corr_gp = where(gmask[pix_list] gt 0.)

        if (corr_gp[0] ne -1) then correlation_maps[i] = correlate(map1_in[pix_list[corr_gp]], map2_in[pix_list[corr_gp]])
    endfor

    mollview, correlation_maps, chars=1.5, min=-1, max=1, tit='!6Correlation Map', /nest

    return, correlation_maps
end
