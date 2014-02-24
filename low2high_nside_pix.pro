function low2high_nside_pix, high_nside, low_nside, order
   print, ' - Mapping '+strtrim(string(high_nside),2)+' to '+strtrim(string(low_nside),2)
   print, ' - Ordering: ', order

    Npix     = 12l * high_Nside^2
    low_Npix = 12l * low_Nside^2

    ipix = lindgen(Npix)
    if (order eq 'ring') then begin
        pix2vec_ring, high_Nside, ipix, vec
        vec2pix_ring, low_Nside, vec, low_ipix
    endif

    if (order eq 'nest') then begin
        pix2vec_nest, high_Nside, ipix, vec
        vec2pix_nest, low_Nside, vec, low_ipix
    endif

    print, 'Done...'
   return, low_ipix
end
