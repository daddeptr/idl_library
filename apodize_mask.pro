function apodize_mask, mask, nest=nest, radius_arcmin=radius_arcmin, cos=cos, gauss=gauss, verbose=verbose, write=write, root=root

   if not keyword_set(root) then root = 'apodized_mask'

   if not keyword_set(radius_arcmin) then begin
       radius_arcmin = 90.
       print, ' - default radius [arcmin] = ', radius_arcmin
   endif else begin
       print, ' - radius [arcmin] = ', radius_arcmin
   endelse
   if keyword_set(gauss) then stop, ' Gaussian filter not implemented yet'

   mask = float( mask )
   npix = n_elements( mask )
   nside = long( sqrt( long(npix) / 12l ) )
   print, ' - nside = ', nside

   if not keyword_set(nest) then begin
       if keyword_set(verbose) then print, ' - reordering...'
       mask = reorder(mask, in='ring', out='nest')
       print, ' - reordered'
       if keyword_set(verbose) then mollview, mask, /nest, win=19, px=650
   endif

   mpix = where( mask eq 0. )
   nmpix = n_elements( mpix )
   if keyword_set(verbose) then print, ' - vanishing pixels: ', nmpix

   pix2ang_nest, nside, mpix, theta, phi
   pix2vec_nest, nside, mpix, vec0
   pix2vec_nest, nside, lindgen( npix ), vec
   
   if nmpix[0] gt 0 then begin
       for ipix=0l, nmpix-1 do begin
           neighbours_nest, nside, mpix[ipix], listp
;           help, listp
;stop
           if total( mask[listp] gt 0.) then begin
;               print, ipix, nmpix
               if ( (ipix/fix(nmpix/10))*fix(nmpix/10) eq ipix) then print, float(ipix)/nmpix, '% completed'
               query_disc, nside, vec0[ipix,*], radius_arcmin/60./!radeg, listp, nlist, /nest
               dot_prod = transpose( vec0[ipix,*]) ## vec[listp,*]
               alpha = acos( dot_prod )
;##               apodization = sin( alpha/max(alpha) * !dpi/2. )
               apodization = ( 1.+cos( !pi*(1.-alpha/max(alpha)) ) ) /2. ;sin( alpha/max(alpha) * !dpi/2. )
;           help, dot_prod
;           print, dot_prod[0]
;mx = mask*0.
;mx[listp]=1.
;mx[mpix[ipix]]=-1.6375e30
;gnomview, mx, /nest, rot=[45,1], reso=0.5, win=20
;ma = mask*0.
;##               mask[listp] = mask[listp] * apodization
               mask[listp] = mask[listp] * apodization
               mask[mpix[ipix]] = 0.
;           ma = mask
;           ma[mpix[ipix]]=-1.6375e30
               if ( ( (ipix/fix(nmpix/10))*fix(nmpix/10) eq ipix) and keyword_set(verbose) )then gnomview, mask, /nest, rot=[phi[ipix],!pi/2-theta[ipix]]*!radeg, reso=0.5, win=21
;##stop
           endif else begin
;##               print, ' - skipping...???', total( mask[listp] )
           endelse
       endfor
   endif else begin
       print, 'zero-value pixels not found'
   endelse
   
   if keyword_set(verbose) then mollview, mask, /nest, win=21  

   if keyword_set(write) then write_fits_map, root+'_apoFullCos'+strtrim(string(radius_arcmin),2)+'_nest.fits', mask, /nest
   if not keyword_set(nest) then begin
       print, ' - returning ring order...'
       mask=reorder(mask, in='nest', out='ring')
   endif

   return, mask

end
