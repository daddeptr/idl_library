
   True = 1b
   False = 0b

   win = 10
   freq = '217'
   xflmax = '3000'

   mollview, findgen(12*64l^2), win=win, px=450
   loadct, 39
   !p.color=0
   !p.background=255

   dir = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/'

;   mymoll, 'data/mask/dx11/Tmask.fits', win=1
;   mymoll, 'data/mask/dx11/cons-apo_Tmask.fits', win=2
;   mymoll, 'data/mask/dx11/cons-apo_Tmask2.fits', win=2

   read_fits_map, 'data/mask/dx11/Tmask.fits', Tmask
   read_fits_map, 'data/mask/dx11/cons-apo_Tmask.fits', coTmask
   read_fits_map, 'data/mask/dx11/cons-apo_Tmask2.fits', coTmask2
   read_fits_map, 'data/mask/dx11/cons-apo_Tmask3.fits', coTmask3

   print, total(Tmask)/n_elements(Tmask)
   print, total(coTmask)/n_elements(coTmask)
   print, total(coTmask2)/n_elements(coTmask2)
   print, total(coTmask3)/n_elements(coTmask3)

stop

   splits = [ $
          'dx11_hr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs_xfcl_l'+xflmax+'_Tmask_x_Pmask_l4000', $
          'dx11_hr_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l'+xflmax+'_Tmask_x_Pmask_l4000', $
          'dx11_ds_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs_xfcl_l'+xflmax+'_Tmask_x_Pmask_l4000', $
          'dx11_ds_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l'+xflmax+'_Tmask_x_Pmask_l4000', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs_xfcl_l'+xflmax+'_Tmask_x_Pmask_l4000', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l'+xflmax+'_Tmask_x_Pmask_l4000', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_co_undusted_split_ns2048_uK_hrhs_xfcl_l'+xflmax+'_cons-apo_Tmask_x_Pmask_l4000' $
            ]
         
   maps = [ $
          'dx11_hr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs', $
          'dx11_hr_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs', $
          'dx11_ds_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs', $
          'dx11_ds_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_co_undusted_split_ns2048_uK_hrhs' $
            ]
         
   noises = [ $
          'dx11_hr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhd', $
          'dx11_hr_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhd', $
          'dx11_ds_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhd', $
          'dx11_ds_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhd', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhd', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhd', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_co_undusted_split_ns2048_uK_hrhd' $
            ]
         
   masks = [ $
          '_xfcl_l'+xflmax+'_Tmask_x_Pmask_l4000', $
          '_xfcl_l'+xflmax+'_Tmask_x_Pmask_l4000', $
          '_xfcl_l'+xflmax+'_Tmask_x_Pmask_l4000', $
          '_xfcl_l'+xflmax+'_Tmask_x_Pmask_l4000', $
          '_xfcl_l'+xflmax+'_Tmask_x_Pmask_l4000', $
          '_xfcl_l'+xflmax+'_Tmask_x_Pmask_l4000', $
          '_xfcl_l'+xflmax+'_cons-apo_Tmask_x_Pmask_l4000' $
            ]
         
   nfiles = [ $
          'dx11_hr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhd_1_cls_Tmask_x_Pmask.fits', $
          'dx11_hr_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhd_1_cls_Tmask_x_Pmask.fits', $
          'dx11_ds_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhd_1_cls_Tmask_x_Pmask.fits', $
          'dx11_ds_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhd_1_cls_Tmask_x_Pmask.fits', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhd_1_cls_Tmask_x_Pmask.fits', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhd_1_cls_Tmask_x_Pmask.fits', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_co_undusted_split_ns2048_uK_hrhd_1_cls_cons-apo_Tmask_x_Pmask.fits' $
            ]

   sfiles = [ $
          'dx11_hr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs_cls_Tmask_x_Pmask.fits', $
          'dx11_hr_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_cls_Tmask_x_Pmask.fits', $
          'dx11_ds_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs_cls_Tmask_x_Pmask.fits', $
          'dx11_ds_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_cls_Tmask_x_Pmask.fits', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs_cls_Tmask_x_Pmask.fits', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_cls_Tmask_x_Pmask.fits', $
          'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_co_undusted_split_ns2048_uK_hrhs_cls_cons-apo_Tmask_x_Pmask.fits' $
            ]


   tags = 'DX11 '+freq+': '+['hr-raw', 'hr-clean','ds-raw', 'ds-clean', 'yr-raw', 'yr-clean', 'yr-clean-co']

   newdat = dir+splits+'/'+splits+'.newdat'

   nsets = n_elements(splits)

   weights1 = [ $
               [1,0], $
               [1.00199304, -0.00199304], $
               [1,0], $
               [1.00197605, -0.00197605], $
               [1,0], $
               [1.00199282, -0.00199282], $
               [1.00257562, -0.00257562] $
             ]

   weights2 = [ $
               [1,0], $
               [1.00199301, -0.00199301], $
               [1,0], $
               [1.00197096, -0.00197096], $
               [1,0], $
               [1.00199329, -0.00199329], $
               [1.00257574, -0.00257574] $
             ]

   window, win, xsize=720*1.25, ysize=450*1.25
   plot, /nodata, [1,5000], [10, 6500], chars=1.5, ys=1, /ylog, /xlog
   l=findgen(5001)
   ll=l*(l+1)/2./!pi
   lns = [0,2,0,2,0,2,2]
   for iset=1,nsets-1,2 do begin
;##       read_fits_map, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11/'+maps[iset]+'.fits', map
;##       print, ' - minmax = ', minmax(map)
       cls = 0.
       nls = 0.
;       ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11/'+maps[iset]+'.fits', dir+splits[iset]+'/'+maps[iset]+'_cls.fits', simul_type=2, nlmax=5000, /silent
;       ianafast, '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11/'+noises[iset]+'.fits', dir+splits[iset]+'/'+noises[iset]+'_cls.fits', simul_type=2, nlmax=5000, /silent
;##       fits2cl, cls, dir+splits[iset]+'/'+maps[iset]+'_cls.fits'
;##       fits2cl, nls, dir+splits[iset]+'/'+noises[iset]+'_cls.fits'
       fits2cl, cls, dir+splits[iset]+'/'+sfiles[iset]
       fits2cl, nls, dir+splits[iset]+'/'+nfiles[iset]
;       oplot, l, cls*ll/gaussbeam(7.03,5000)^2, col=40*(iset+1)
;       oplot, l, nls*ll/gaussbeam(7.03,5000)^2, col=40*(iset+1)
       oplot, l, cls[*,0]*ll, col=60*((iset)/2), line=lns[iset], ns=50, thick=2
       oplot, l, nls[*,0]*ll, col=60*((iset)/2), line=lns[iset], ns=50, thick=2
;##       oplot, l, cls[*,0]/nls[*,0], col=60*((iset)/2), line=lns[iset], ns=50, thick=2
       if total(cls) eq 0. then stop, 'ERROR somewhere: '+maps[iset]
   endfor
   iset -= 1
       fits2cl, cls, dir+splits[iset]+'/'+sfiles[iset]
       fits2cl, nls, dir+splits[iset]+'/'+nfiles[iset]
       oplot, l, cls[*,0]*ll, col=60*((iset)/2), line=lns[iset], ns=50, thick=2
       oplot, l, nls[*,0]*ll, col=60*((iset)/2), line=lns[iset], ns=50, thick=2
   legend, tags+', ILC w='+string((weights1[0,*]))+string((weights2[0,*])), col=((lindgen(nsets))/2)*60, psym=4, /bottom, /right, chars=1.

   stop

   for iset=0,nsets-1 do begin
       for jset=iset+1, nsets-1 do begin
           compare_xfaster_spectra, newdat[iset], newdat[jset], leg_tags=[tags[iset], tags[jset]], win=jset^2+iset+1, ylog=1, otit='Difference: '+tags[iset]+' - '+tags[jset], lmax=long(xflmax)
stop
       endfor
   endfor
   

end
