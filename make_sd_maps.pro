;## CHECK ORDERING
pro make_sd_maps, file1,file2,root

   print, 'Read map 1: '+file1
   read_fits_map, file1, m1, order=dpord1, nside=dpns1
   print, dpord1

   print, 'Read map 2: '+file2
   read_fits_map, file2, m2, order=dpord2, nside=dpns2
   print, dpord1

   if (dpns1 ne dpns2) then stop, 'Different Nside'
   ns = string(dpns1, format='(i4)')
   if (dpord1 ne dpord2) then begin
       print, 'reordering...'
       ud_grade, m2,m2o, order_in=dpord2, order_out=dpord1       
       m2o = m2
       mollview, m1, ordering=dpord1, px=500
       mollview, m2, ordering=dpord1, px=500
   endif

   sz1 = size(m1)
   sz2 = size(m2)
;##    print, sz1, sz2
   if (sz1[2] ne sz2[2]) then stop, 'Check map dimensions'
   if (sz1[2] eq 1) then begin
       print, sz1[2]
       write_fits_map, root+'_ns'+ns+'_uK_hrhs.fits',(m1+m2)/2., order=dpord1, units='!7l!6K'
       write_fits_map, root+'_ns'+ns+'_uK_hrhd.fits',(m1-m2)/2., order=dpord1, units='!7l!6K'
   endif else begin
       print, sz1[2]
       print, 'Writing maps...'
       print, root+'_ns'+ns+'_uK_hrhs.fits'
       write_tqu, root+'_ns'+ns+'_uK_hrhs.fits', float( (m1+m2)/2. ), order=dpord1, units='!7l!6K'

       print, root+'_ns'+ns+'_uK_hrhd.fits'
       write_tqu, root+'_ns'+ns+'_uK_hrhd.fits', float( (m1-m2)/2. ), order=dpord1, units='!7l!6K'
   endelse
   mollview, root+'_ns'+ns+'_uK_hrhs.fits', win=10, px=600, min=-300, max=300
   mollview, root+'_ns'+ns+'_uK_hrhd.fits', win=11, px=600, min=-100, max=100

end

