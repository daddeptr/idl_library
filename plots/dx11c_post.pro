   True = 1b
   False = 0b

   do_maps = False
   do_spec = True

   freq = ['100','143','217']

   if do_maps then begin
       fmaps_raw = '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_yr_1-2_IQUmap_'+freq+'_ns2048_uK_hrhs.fits'
       fmaps_clean = '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_yr_1-2_IQUmap_'+freq+'_extMask_545_coTP_undusted_split_ns2048_uK_hrhs.fits'

;##   fmaps_rawn = '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_yr_1-2_IQUmap_'+freq+'_ns2048_uK_hrhd.fits'
;##   for i=0,2 do spawn, 'mv /global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11c/dx11c_yr_1-2_IQUmap_'+freq[i]+'_full_ns2048_uK_hrhd.fits '+fmaps_rawn[i]

       beams = [9.64, 7.03, 4.83]
       smth = 180
       maskf = [ $
                 '/global/scratch2/sd/dpietrob/Software/XFaster/data/mask/dx11c/cons-apo_Tmask3_mPix.fits', $
                 '/global/scratch2/sd/dpietrob/Software/XFaster/data/mask/dx11c/cons-apo_Pmask3_mPix.fits', $
                 '/global/scratch2/sd/dpietrob/Software/XFaster/data/mask/dx11c/cons-apo_Pmask3_mPix.fits' $
               ]

;read_fits_map, fmaps_clean[2], m, order=ord
;if ord ne 'RING' then m=reorder(m,in=ord,out='ring')
;write_fits_map, 'monopole.fits', fltarr(12l*2048l^2)+1., /ring

if False then begin
    tmp01 = map_regression(['monopole.fits','data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year1.fits.b100'],maskfile='uscomp/masks/mask_CO3.fits',infile='/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_IQUmap_100_year1_uK.fits', A_coeff=A01)
    tmp02 = map_regression(['monopole.fits','data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year2.fits.b100'],maskfile='uscomp/masks/mask_CO3.fits',infile='/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_IQUmap_100_year2_uK.fits', A_coeff=A02)
    write_fits_map, 'data/maps/dx11c/dx11c_SDR_100_yrs_Treg.fits', (tmp01/(1.-A01[1])+tmp02/(1.-A02[1]))/2., /ring
    write_fits_map, 'data/maps/dx11c/dx11c_SDR_100_yrd_Treg.fits', (tmp01/(1.-A01[1])-tmp02/(1.-A02[1]))/2., /ring

    print, (A01/(1.-A01[1]) + A02/(1.-A02[1]))/2.
    
    tmp11 = map_regression(['monopole.fits','data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year1.fits.b143'],maskfile='uscomp/masks/mask_CO3.fits',infile='/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_IQUmap_143_year1_uK.fits', A_coeff=A11)
    tmp12 = map_regression(['monopole.fits','data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year2.fits.b143'],maskfile='uscomp/masks/mask_CO3.fits',infile='/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_IQUmap_143_year2_uK.fits', A_coeff=A12)
    write_fits_map, 'data/maps/dx11c/dx11c_SDR_143_yrs_Treg.fits', (tmp11/(1.-A11[1])+tmp12/(1.-A12[1]))/2., /ring
    write_fits_map, 'data/maps/dx11c/dx11c_SDR_143_yrd_Treg.fits', (tmp11/(1.-A11[1])-tmp12/(1.-A12[1]))/2., /ring

    print, (A11/(1.-A11[1]) + A12/(1.-A12[1]))/2.
    
    tmp21 = map_regression(['monopole.fits','data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year1.fits.b217'],maskfile='uscomp/masks/mask_CO3.fits',infile='/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_IQUmap_217_year1_uK.fits', A_coeff=A21)
    tmp22 = map_regression(['monopole.fits','data/maps/dx11c/dx11c_solarDipoleRemoved_545-353_templates_year2.fits.b217'],maskfile='uscomp/masks/mask_CO3.fits',infile='/global/scratch2/sd/dpietrob/dx11c/maps/dx11c_IQUmap_217_year2_uK.fits', A_coeff=A22)
    write_fits_map, 'data/maps/dx11c/dx11c_SDR_217_yrs_Treg.fits', (tmp21/(1.-A21[1])+tmp22/(1.-A22[1]))/2., /ring
    write_fits_map, 'data/maps/dx11c/dx11c_SDR_217_yrd_Treg.fits', (tmp21/(1.-A21[1])-tmp22/(1.-A22[1]))/2., /ring

    print, (A21/(1.-A21[1]) + A22/(1.-A22[1]))/2.
endif
   fmaps_Treg = 'data/maps/dx11c/dx11c_SDR_'+freq+'_yrs_Treg.fits'
   newdatf = '/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_yr_1-2_IQUmap_'+freq+'_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_xfcl_l3000_cons-apo_Tmask3_mPix_x_cons-apo_Pmask3_mPix_l4000.newdat'
   newdatf_largeMask = '/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_yr_1-2_IQUmap_'+freq+'_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_xfcl_l3000_pla_2013_largestMask_x_cons-apo_Pmask3_mPix_l4000.newdat'
;   newdatf_raw = '/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_yr_1-2_IQUmap_'+freq+'_ns2048_uK_hrhs_xfcl_l3000_cons-apo_Tmask3_mPix_x_cons-apo_Pmask3_mPix_l4000.newdat'
   newdatf_raw_largeMask = '/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_yr_1-2_IQUmap_'+freq+'_ns2048_uK_hrhs_xfcl_l3000_pla_2013_largestMask_x_cons-apo_Pmask3_mPix_l4000.newdat'
;   mymoll, fmaps_raw[2], maskfile=maskf[0], /ash10, px=1000, units='!6asinh10(!7l!6K)', tit=' ', /no_dipole, png='dx11c_217_rawmap.png', min=-300, max=300
;   mymoll, fmaps_clean[2], maskfile=maskf[0], /ash10, px=1000, units='!6asinh10(!7l!6K)', tit=' ', /no_dipole, png='dx11c_217_cleanmap.png', min=-300, max=300
;   mymoll, fmaps_clean[2], maskfile=maskf[0], px=1000, units='!6asinh10(!7l!6K)', tit='217-100', /no_dipole, file2=fmaps_clean[0], win=1, /ash10, png='dx11c_217-100_cleanmap.png'
;   mymoll, fmaps_clean[2], maskfile=maskf[0], px=1000, units='!6asinh10(!7l!6K)', tit='217-143', /no_dipole, file2=fmaps_clean[1], win=1, /ash10, png='dx11c_217-143_cleanmap.png'
;   mymoll, fmaps_clean[1], maskfile=maskf[0], px=1000, units='!6asinh10(!7l!6K)', tit='143-100', /no_dipole, file2=fmaps_clean[0], win=1, /ash10, png='dx11c_143-100_cleanmap.png'
   mymoll, fmaps_Treg[2], maskfile=maskf[0], px=1000, units='!6asinh10(!7l!6K)', tit='217-100', /no_dipole, file2=fmaps_Treg[0], win=1, /ash10, png='dx11c_217-100_Tregmap.png'
   mymoll, fmaps_Treg[2], maskfile=maskf[0], px=1000, units='!6asinh10(!7l!6K)', tit='217-143', /no_dipole, file2=fmaps_Treg[1], win=1, /ash10, png='dx11c_217-143_Tregmap.png'
   mymoll, fmaps_Treg[1], maskfile=maskf[0], px=1000, units='!6asinh10(!7l!6K)', tit='143-100', /no_dipole, file2=fmaps_Treg[0], win=1, /ash10, png='dx11c_143-100_Tregmap.png'
;   mymoll, fmaps_Treg[2], maskfile=maskf[0], /ash10, px=1000, units='!6asinh10(!7l!6K)', tit='217-143', /no_monopole, file2=fmaps_Treg[1], min=-30, max=30, dosmooth=60., win=2, png='dx11c_217-143_Tregmap.png'
;   mymoll, fmaps_Treg[1], maskfile=maskf[0], /ash10, px=1000, units='!6asinh10(!7l!6K)', tit='143-100', /no_monopole, file2=fmaps_Treg[0], min=-30, max=30, dosmooth=60., win=3, png='dx11c_143-100_Tregmap.png'

;   compare_xfaster_spectra, newdatf_raw_largeMask[2], newdatf[2], /init, reso=150, binning='data/bins/const/const2',leg_tags='dx11c 217 GHz:'+[' Raw on largest 2013 Mask',' Clean on small mask'], uleg_pos=[1300,5800], lmax=3000, win=3;, /dops, psout='dx11c_217_maskcomp'

;   compare_xfaster_spectra, newdatf_raw_largeMask[2], newdatf_largeMask[2], /init, reso=150, binning='data/bins/const/const2',leg_tags='dx11c 217 GHz:'+[' Raw',' Clean']+' on 2013 largest mask', uleg_pos=[1200,5800], lmax=3000, /dops, psout='dx11c_217_maskcomp', otit=' '

stop
   field = ['I','Q','U']
;##   if do_maps then begin
       if False then begin
           for i=0,2 do begin
               for j=0,2 do mymoll, fmaps_raw[i], /ash10,tit='DX11c: '+freq[i]+' Raw map '+field[j], units='!6Asinh10(!7l!6K)', win=-(i+1), px=800, imap=j, png='dx11c_'+freq[i]+'_raw_'+field[j]+'map.png'
;
               for j=0,2 do mymoll, fmaps_clean[i], /ash10,tit='DX11c: '+freq[i]+' Clean map '+field[j], units='!6Asinh10(!7l!6K)', win=-(i+1), px=800, imap=j, png='dx11c_'+freq[i]+'_clean_'+field[j]+'map.png'
           endfor
           for j=0,2 do spawn,'convert -delay 100 dx11c_*_raw_'+field[j]+'map.png -loop 0 dx11c_raw_'+field[j]+'maps.gif'
           for j=0,2 do spawn,'convert -delay 100 dx11c_*_clean_'+field[j]+'map.png -loop 0 dx11c_clean_'+field[j]+'maps.gif'

           mollview, '/global/scratch2/sd/dpietrob/Software/XFaster/data/mask/dx11_pre/hfi_co3_mask_x_P_ns2048.fits', chars=1.5, px=800,tit='Cleaning mask: T and P', grat=[20,20], png='dx11c_cleaning_mask.png'
           mollview, '/global/scratch2/sd/dpietrob/Software/XFaster/data/mask/dx11c/cons-apo_Tmask3_mPix.fits', chars=1.5, px=800,tit='Power spectrum mask: T', grat=[20,20], png='dx11c_T_mask.png'
           mollview, '/global/scratch2/sd/dpietrob/Software/XFaster/data/mask/dx11c/cons-apo_Pmask3_mPix.fits', chars=1.5, px=800,tit='Power spectrum mask: P', grat=[20,20], png='dx11c_P_mask.png'

           j=0
           for i=0,2 do mymoll, fmaps_clean[i], /ash10,tit='DX11c: '+freq[i]+' Clean map '+field[j], units='!6Asinh10(!7l!6K)', win=-(i+1), px=800, imap=j, maskfile=maskf[j], /no_dipole, png='dx11c_'+freq[i]+'_clean_'+field[j]+'map_masked.png'
           spawn,'convert -delay 100 dx11c_*_clean_'+field[j]+'map_masked.png -loop 0 dx11c_clean_'+field[j]+'maps_masked.gif'
           
           for i=0,2 do mymoll, fmaps_clean[i], /ash10,tit='DX11c: '+freq[i]+' Clean map '+field[j]+' '+strtrim(string(smth),2)+'arcmin', units='!6Asinh10(!7l!6K)', win=(i+1), px=800, imap=j, maskfile=maskf[j], /no_dipole, dosmooth=sqrt(smth^2-beams[i]^2), png='dx11c_'+freq[i]+'_clean_'+field[j]+'map_masked_smooth.png'
           spawn,'convert -delay 100 dx11c_*_clean_'+field[j]+'map_masked_smooth.png -loop 0 dx11c_clean_'+field[j]+'maps_masked_smooth.gif'

       j=1
       for i=0,2 do mymoll, fmaps_clean[i], /ash10,tit='DX11c: '+freq[i]+' Clean map '+field[j]+' '+strtrim(string(smth),2)+'arcmin', units='!6Asinh10(!7l!6K)', win=(i+1), px=800, imap=j, maskfile=maskf[j], /no_dipole, dosmooth=sqrt(smth^2-beams[i]^2), png='dx11c_'+freq[i]+'_clean_'+field[j]+'map_masked_smooth.png'
       spawn,'convert -delay 100 dx11c_*_clean_'+field[j]+'map_masked_smooth.png -loop 0 dx11c_clean_'+field[j]+'maps_masked_smooth.gif'

       j=2
       for i=0,2 do mymoll, fmaps_clean[i], /ash10,tit='DX11c: '+freq[i]+' Clean map '+field[j]+' '+strtrim(string(smth),2)+'arcmin', units='!6Asinh10(!7l!6K)', win=(i+1), px=800, imap=j, maskfile=maskf[j], /no_dipole, dosmooth=sqrt(smth^2-beams[i]^2), png='dx11c_'+freq[i]+'_clean_'+field[j]+'map_masked_smooth.png'
       spawn,'convert -delay 100 dx11c_*_clean_'+field[j]+'map_masked_smooth.png -loop 0 dx11c_clean_'+field[j]+'maps_masked_smooth.gif'

       mymoll, 'data/maps/dx11_pre/dx11_pre_hr1_545-353_templates.fits', tit='Dust I template: 545', /ash10, units='!6Asinh10(!7l!6K)', win=-3, px=800, imap=0, png='predx11_545_Idust_template.png'
       mymoll, 'data/maps/dx11_pre/dx11_pre_hr1_545-353_templates.fits', tit='Dust Q template: 353', /ash10, units='!6Asinh10(!7l!6K)', win=-4, px=800, imap=1, png='predx11_353_Qdust_template.png'
       mymoll, 'data/maps/dx11_pre/dx11_pre_hr1_545-353_templates.fits', tit='Dust U template: 353', /ash10, units='!6Asinh10(!7l!6K)', win=-5, px=800, imap=2, png='predx11_353_Udust_template.png'
       endif
   endif

   if do_spec then begin
       newdatf = '/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_yr_1-2_IQUmap_'+freq+'_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_xfcl_l3000_cons-apo_Tmask3_mPix_x_cons-apo_Pmask3_mPix_l4000.newdat'
       newdatf_largeMask = '/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_yr_1-2_IQUmap_'+freq+'_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_xfcl_l3000_pla_2013_largestMask_x_cons-apo_Pmask3_mPix_l4000.newdat'
       newdatf_raw = '/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_yr_1-2_IQUmap_'+freq+'_ns2048_uK_hrhs_xfcl_l3000_cons-apo_Tmask3_mPix_x_cons-apo_Pmask3_mPix_l4000.newdat'
       newdatf_raw_largeMask = '/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_yr_1-2_IQUmap_'+freq+'_ns2048_uK_hrhs_xfcl_l3000_pla_2013_largestMask_x_cons-apo_Pmask3_mPix_l4000.newdat'
       lr = [450, 2250, 2300]
;##       for i=0,2 do plot_xfaster_newdat, newdatf[i], lmax=lr[i], win=i+1, /pol, /resi, /init, binning='data/bins/const/const2', /dops, psout='dx11c_xf-spec_'+freq[i]+'.eps'
;##       for i=2,0,-1 do for j=i-1,0,-1 do compare_xfaster_spectra, newdatf[i], newdatf[j], lmax=lr[i], win=i^2+j+1, /init, binning='data/bins/const/const2', leg_tags='dx11c '+[freq[i],freq[j]]+' clean map', otit='XFaster spectra comparison: TT';, /dops, psout='dx11c_xf-spec_comparison_'+freq[i]+'-'+freq[j]
       if False then begin
           bint = 'const2'
           binp = '/global/homes/d/dpietrob/myscratch/Software/XFaster/data/bins/const/const2'
           newdatf = '/global/homes/d/dpietrob/myscratch/Software/XFaster/uscomp/xfaster/dx11c_SDR_yr_1-2_IQUmap_'+freq+'_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_'+bint+'_fisherWins_xfcl_l3000_Tmask_badpixmasked_x_Pmask_badpixmasked_l4000.newdat'
           newdatf[2] = '/global/homes/d/dpietrob/myscratch/Software/XFaster/uscomp/xfaster/dx11c_SDR_yr_1-2_IQUmap_recal_'+freq[2]+'_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_'+bint+'_fisherWins_xfcl_l3000_Tmask_badpixmasked_x_Pmask_badpixmasked_l4000.newdat'
           for i=2,1,-1 do for j=i-1,1,-1 do compare_xfaster_spectra, newdatf[i], newdatf[j], lmax=2300, win=i^2+j+1, /init, binning=binp, leg_tags='dx11c '+[freq[i],freq[j]]+' clean map', otit='XFaster spectra comparison: TT', /dops, psout='dx11c_xf-spec_comparison_'+freq[i]+'-'+freq[j]+'_v2',uleg_pos=[500,5500]
;stop
           for i=2,1,-1 do for j=i-1,1,-1 do compare_xfaster_spectra, newdatf[i], newdatf[j], lmax=1000, win=i^2+j+1, /init, binning=binp, leg_tags='dx11c '+[freq[i],freq[j]]+' clean map', otit='XFaster spectra comparison: EE', ncl=2, reso=2.05, /dops, psout='dx11c_xf-spec_comparison_'+freq[i]+'-'+freq[j]+'_v2'
;stop
           for i=2,1,-1 do for j=i-1,1,-1 do compare_xfaster_spectra, newdatf[i], newdatf[j], lmax=1600, win=i^2+j+1, /init, binning=binp, leg_tags='dx11c '+[freq[i],freq[j]]+' clean map', otit='XFaster spectra comparison: TE', ncl=4, reso=10.15, /dops, psout='dx11c_xf-spec_comparison_'+freq[i]+'-'+freq[j]+'_v2'
           stop
       endif

       if False then begin
           r = [2.5,2.5,2.5]
           for i=2,0,-1 do for j=i-1,0,-1 do compare_xfaster_spectra, newdatf[i], newdatf[j], lmax=1000, win=i^2+j+1, /init, binning='data/bins/const/const2', leg_tags='dx11c '+[freq[i],freq[j]]+' clean map', otit='XFaster spectra comparison: EE' ;, reso=r[i], ncl=2, uleg_pos=[50,45], /dops, psout='dx11c_xf-spec_comparison_'+freq[i]+'-'+freq[j]
           r = [10.5,10.5,10.5]
           for i=2,0,-1 do for j=i-1,0,-1 do compare_xfaster_spectra, newdatf[i], newdatf[j], lmax=1600, win=i^2+j+1, /init, binning='data/bins/const/const2', leg_tags='dx11c '+[freq[i],freq[j]]+' clean map', otit='XFaster spectra comparison: TE' ;, reso=r[i], ncl=4, uleg_pos=[950,145], /dops, psout='dx11c_xf-spec_comparison_'+freq[i]+'-'+freq[j]
           stop
       endif
;       r = [155,155,350]
;##       for i=2,0,-1 do for j=i-1,0,-1 do compare_xfaster_spectra, newdatf_raw[i], newdatf_raw[j], lmax=lr[i], win=i^2+j+1, /init, binning='data/bins/const/const2', leg_tags='dx11c '+[freq[i],freq[j]]+' raw map', otit='XFaster spectra comparison: raw maps', reso=r[i], /dops, psout='dx11c_xf-spec_raw-comparison_'+freq[i]+'-'+freq[j]
;## --- EE
;       r = [2.5,2.5,2.5]
;       for i=2,0,-1 do for j=i-1,0,-1 do compare_xfaster_spectra, newdatf_raw[i], newdatf_raw[j], lmax=1000, win=i^2+j+1, /init, binning='data/bins/const/const2', leg_tags='dx11c '+[freq[i],freq[j]]+' raw map', otit='XFaster spectra comparison: raw maps', reso=r[i], ncl=2;, /dops, psout='dx11c_xf-spec_raw-comparison_'+freq[i]+'-'+freq[j]
;stop
;##       for i=0,2 do compare_xfaster_spectra, newdatf_raw[i], newdatf[i], lmax=lr[i], win=i+1, /init, binning='data/bins/const/const2', leg_tags='dx11c '+freq[i]+[' raw map',' clean map'], otit='XFaster spectra comparison: raw vs clean', reso=r[i], /dops, psout='dx11c_xf-spec_raw-clean_comparison_'+freq[i]
;stop       
;##       for i=0,2 do compare_xfaster_spectra, newdatf[i], newdatf_largeMask[i], lmax=lr[i], win=i+1, /init, binning='data/bins/const/const2', leg_tags='dx11c '+freq[i]+['',' large Mask'], otit='XFaster spectra comparison: mask dependence', reso=100, /dops, psout='dx11c_xf-spec_mask_comparison_'+freq[i]
;##       r = [100,100,350]
;##       for i=0,2 do compare_xfaster_spectra, newdatf_raw[i], newdatf_raw_largeMask[i], lmax=lr[i], win=i+1, /init, binning='data/bins/const/const2', leg_tags='dx11c raw map '+freq[i]+['',' large Mask'], otit='XFaster spectra comparison: mask dependence', reso=r[i], /dops, psout='dx11c_xf-spec_mask_comparison_raw_'+freq[i]
       r = [100,100,350]
       bestf = '/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/' + [ $
                  ['xf_chains_100_TP_const2_l1750_lensedCls.dat','100_TP_const2_l1750_fg_bestfit.txt'], $
                  ['xf_chains_143_TP_const2_l2000_lensedCls.dat','143_TP_const2_l2000_fg_bestfit.txt'], $
                  ['xf_chains_217_TP_const2_l2500_lensedCls.dat','217_TP_const2_l2500_fg_bestfit.txt'] $
                                                                                                      ]
; EE ------------------------------------------------------------------
;143
       xfcls = extract_xfaster_newdat( newdatf[1], lcen=lc, ncl=2, cler=er, covmat=cvmt )
       readcol, bestf[1,1], x, x, y
       l = findgen(5001)
       ll = l*(l+1)/2./!pi
       bfg = xf_binning(y/ll,'/global/homes/d/dpietrob/myscratch/Software/XFaster/data/bins/const/const2_EE')
       openw,1,'/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_143_egfg-sub_xfcls_EE.txt'
       printf,1,'# lcen cls cls_err egfg cls-egfg'
       for i=0,n_elements(lc)-1 do printf,1,lc[i],xfcls[i],er[i],bfg[i],xfcls[i]-bfg[i],format='(i12,4f15.8)'
       close,1
       cvmtl=n_elements(cvmt[0,*])
;stop
       openw,1,'/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_143_xf_covmat_'+strtrim(string(cvmtl),2)+'_x_'+strtrim(string(cvmtl),2)+'.txt'
       for i=0,cvmtl-1 do printf,1,cvmt[i,*],format='('+strtrim(string(cvmtl),2)+'f15.8)'
       close,1

;217
       xfcls2 = extract_xfaster_newdat( newdatf[2], lcen=lc, ncl=2, cler=er2, covmat=cvmt2 )
       readcol, bestf[1,2], x, x, y
       ll = l*(l+1)/2./!pi
       bfg2 = xf_binning(y/ll,'/global/homes/d/dpietrob/myscratch/Software/XFaster/data/bins/const/const2_EE')
       openw,1,'/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_217_egfg-sub_xfcls_EE.txt'
       printf,1,'# lcen cls cls_err egfg cls-egfg'
       for i=0,n_elements(lc)-1 do printf,1,lc[i],xfcls2[i],er2[i],bfg2[i],xfcls2[i]-bfg2[i],format='(i12,4f15.8)'
       close,1
       cvmtl=n_elements(cvmt2[0,*])
       openw,1,'/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_217_xf_covmat_'+strtrim(string(cvmtl),2)+'_x_'+strtrim(string(cvmtl),2)+'.txt'
       for i=0,cvmtl-1 do printf,1,cvmt2[i,*],format='('+strtrim(string(cvmtl),2)+'f15.8)'
       close,1
;stop
       if False then for i=0,2 do plot_xfaster_newdat, newdatf[i], lmax=lr[i], win=i+1, /init, /pol, /resi, binning='data/bins/const/const2', bestfit=bestf[*,i], /dops, psout='dx11c_xf-spec_residuals_'+freq[i]+'.eps'

       !p.multi=[0,1,2]
       mollview, findgen(12), win=0
       loadct, 39
       !p.color=0
       !p.background=255
       window, 0, xsize=800, ysize=800
       plot, lc, (xfcls-bfg), thick=2, psym=-4, chars=1.5, xr=[75,1150], yr=[-5,50], xs=1, tit='!6EE', xtit='!8l', ytit='!6D!dl!n'
       oplot, lc, (xfcls2-bfg2), thick=2, psym=-4, col=70
       legend, 'DX11c '+['143 yrs','217 yrs'], col=[0,70], psym=[4,4], chars=1.5
       plot, lc, (xfcls-bfg)-(xfcls2-bfg2), thick=2, psym=-4, chars=1.5, xr=[75,1150], yr=[-5,5], xs=1, ys=1, xtit='!8l', ytit='!7D!6!dl!n'
       oplot, lc, lc*0., col=245
       oplot, lc, er, line=2
       oplot, lc, -er, line=2
       xyouts, 100,4,'!6143 GHz - 217 GHz', chars=1.5

; TE
;143
       xfcls = extract_xfaster_newdat( newdatf[1], lcen=lc, ncl=4, cler=er )
       readcol, bestf[1,1], x, x, x, y
       l = findgen(5001)
       ll = l*(l+1)/2./!pi
       bfg = xf_binning(y/ll,'/global/homes/d/dpietrob/myscratch/Software/XFaster/data/bins/const/const2_TE')
       openw,1,'/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_143_egfg-sub_xfcls_TE.txt'
       printf,1,'# lcen cls cls_err egfg cls-egfg'
       for i=0,n_elements(lc)-1 do printf,1,lc[i],xfcls[i],er[i],bfg[i],xfcls[i]-bfg[i],format='(i12,4f15.8)'
       close,1

;217
       xfcls2 = extract_xfaster_newdat( newdatf[2], lcen=lc, ncl=4, cler=er2 )
       readcol, bestf[1,2], x, x, x, y
       ll = l*(l+1)/2./!pi
       bfg2 = xf_binning(y/ll,'/global/homes/d/dpietrob/myscratch/Software/XFaster/data/bins/const/const2_TE')
       openw,1,'/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_217_egfg-sub_xfcls_TE.txt'
       printf,1,'# lcen cls cls_err egfg cls-egfg'
       for i=0,n_elements(lc)-1 do printf,1,lc[i],xfcls2[i],er2[i],bfg2[i],xfcls2[i]-bfg2[i],format='(i12,4f15.8)'
       close,1

       window, 1, xsize=800, ysize=800
       plot, lc, (xfcls-bfg), thick=2, psym=-4, chars=1.5, xr=[75,1750], yr=[-150,150], xs=1, tit='!6TE', xtit='!8l', ytit='!6D!dl!n'
       oplot, lc, (xfcls2-bfg2), thick=2, psym=-4, col=70
       legend, 'DX11c '+['143 yrs','217 yrs'], col=[0,70], psym=[4,4], chars=1.5, /right
       plot, lc, (xfcls-bfg)-(xfcls2-bfg2), thick=2, psym=-4, chars=1.5, xr=[75,1750], yr=[-15,15], xs=1, ys=1, xtit='!8l', ytit='!7D!6!dl!n'
       oplot, lc, lc*0., col=245
       oplot, lc, er, line=2
       oplot, lc, -er, line=2
       xyouts, 100,13,'!6143 GHz - 217 GHz', chars=1.5

; TT
;143
       xfcls = extract_xfaster_newdat( newdatf[1], lcen=lc, ncl=1, cler=er )
       readcol, bestf[1,1], x, y
       l = findgen(5001)
       ll = l*(l+1)/2./!pi
       bfg = xf_binning(y/ll,'/global/homes/d/dpietrob/myscratch/Software/XFaster/data/bins/const/const2_TT')
       openw,1,'/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_143_egfg-sub_xfcls_TT.txt'
       printf,1,'# lcen cls cls_err egfg cls-egfg'
       for i=0,n_elements(lc)-1 do printf,1,lc[i],xfcls[i],er[i],bfg[i],xfcls[i]-bfg[i],format='(i12,4f15.8)'
       close,1

;217
       xfcls2 = extract_xfaster_newdat( newdatf[2], lcen=lc, cler=er2 )
       readcol, bestf[1,2], x, y
       ll = l*(l+1)/2./!pi
       bfg2 = xf_binning(y/ll,'/global/homes/d/dpietrob/myscratch/Software/XFaster/data/bins/const/const2_TT')
       openw,1,'/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_coTP-clean_maps/dx11c_217_egfg-sub_xfcls_TT.txt'
       printf,1,'# lcen cls cls_err egfg cls-egfg'
       for i=0,n_elements(lc)-1 do printf,1,lc[i],xfcls2[i],er2[i],bfg2[i],xfcls2[i]-bfg2[i],format='(i12,4f15.8)'
       close,1

       window, 2, xsize=800, ysize=800
       plot, lc, (xfcls-bfg), thick=2, psym=-4, chars=1.5, xr=[75,2250], yr=[-150,6500], xs=1, tit='!6TT', xtit='!8l', ytit='!6D!dl!n', ys=1
       oplot, lc, (xfcls2-bfg2), thick=2, psym=-4, col=70
       legend, 'DX11c '+['143 yrs','217 yrs'], col=[0,70], psym=[4,4], chars=1.5, /right
       plot, lc, (xfcls-bfg)-(xfcls2-bfg2), thick=2, psym=-4, chars=1.5, xr=[75,2250], yr=[-200,200], xs=1, ys=1, xtit='!8l', ytit='!7D!6!dl!n'
       oplot, lc, lc*0., col=245
       oplot, lc, er, line=2
       oplot, lc, -er, line=2
       xyouts, 100,170,'!6143 GHz - 217 GHz', chars=1.5

       !p.multi=0

   endif


   stop, 'EndEnd of script -'
   do_newdat = False
   do_xnewdat = False
   do_anafast = False

   dir = '/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11c_co-clean_maps/'
   cleanfiles = 'dx11_XO_yr_1-2_IQUmap_'+freq+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000.newdat'
   files = 'dx11_XO_yr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000.newdat'

   sfiles = 'dx11_XO_yr_1-2_IQUmap_'+freq+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_cls_cons-apo_Tmask3_x_cons-apo_Pmask3.fits'
   nfiles = 'dx11_XO_yr_1-2_IQUmap_'+freq+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhd_1_cls_cons-apo_Tmask3_x_cons-apo_Pmask3.fits'

   tmaster = fltarr(4001,3)
   xtmaster = fltarr(4001,3)
   btmaster = fltarr(100,3)
   bxtmaster = fltarr(100,3)
   btt = fltarr(100,3)

   scaling = 0.8

   mp1 = '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11_XO/dx11_XO_v2_yr_1-2_IQUmap_'+freq+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hr1.fits'
   mp2 = '/global/scratch2/sd/dpietrob/Software/XFaster/data/maps/dx11_XO/dx11_XO_v2_yr_1-2_IQUmap_'+freq+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hr2.fits'
   tmask = 'data/mask/dx11/cons-apo_Tmask3.fits'

   pwf = healpixwindow(2048)
   for i=0,2 do begin
       if do_anafast then ianafast, mp1[i], dir+'dx11_XO_Icls_'+freq[i]+'_year1_x_year2_extMask_545_coTP_undusted_split_ns2048_uK_cons-apo_Tmask3.fits', map2_in=mp2[i], maskfile=tmask, regression=2, simul_type=1
       fits2cl, xcls, dir+'dx11_XO_Icls_'+freq[i]+'_year1_x_year2_extMask_545_coTP_undusted_split_ns2048_uK_cons-apo_Tmask3.fits'

       tt = extract_xfaster_newdat(dir+cleanfiles[i], lcen=lc, cler=er, btcl=tcl, binning='data/bins/const/const2')
       btt[*,i] = tt

       readcol, '/global/homes/d/dpietrob/myscratch/Software/XFaster/data/beams/dx11/dx11_wl_'+freq[i]+'.dat', l, bl
       fits2cl, cls, dir+sfiles[i]
       fits2cl, nls, dir+nfiles[i]
       kern = mrdfits( '/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/cons-apo_Tmask3_x_cons-apo_Pmask3_kernel_l4000_v1.9_inverse.fits', 0 )
       help, kern
;       pcls = deconvolve_kernel(cls[*,0]-nls[*,0], fkernel='/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/')
       lmax = n_elements(kern[*,0])+1
       print, lmax
       pcls = kern[0:lmax-2, 0:lmax-2] ## reform( cls[2:lmax,0]-nls[2:lmax,0] )
       tmaster[2:lmax,i] = pcls
       bpcls = xf_binning(tmaster[*,i]/bl^2/pwf^2,'/global/homes/d/dpietrob/myscratch/Software/XFaster/data/bins/const/const2_TT')
       btmaster[0:99,i] = bpcls[0:99]

       pcls = kern[0:lmax-2, 0:lmax-2] ## reform( xcls[2:lmax] )
       xtmaster[2:lmax,i] = pcls
       bpcls = xf_binning(xtmaster[*,i]/bl^2/pwf^2,'/global/homes/d/dpietrob/myscratch/Software/XFaster/data/bins/const/const2_TT')
       bxtmaster[0:99,i] = bpcls[0:99]

       !p.multi=[0,1,2]
       window, i+1, xsize=720*scaling, ysize=450*2*scaling
       plot, lc, tt, chars=1.5, psym=4, yr=[-500, 6500], ys=1, tit='DX11_XO: '+freq[i],xtit='!8l', ytit='!8D!dl!n [!7l!8K!u2!n]'
       oplot, lc, btmaster[*,i], col=245, psym=-4
       oplot, lc, bxtmaster[*,i], col=245, psym=-5
       plot, lc, tt-btmaster[*,i], yr=[-150,150], ys=2, chars=1.5,xtit='!8l',ytit='!6Difference: XFaster-Master [!7l!8K!u2!n]'
       oplot, lc, tt-bxtmaster[*,i]
       oplot, lc, btmaster[*,i]-bxtmaster[*,i], col=160
       oplot, lc, lc*0., line=2
       oplot, lc, er, line=2
       oplot, lc, -er, line=2
       res = linfit(lc[2:66],tt[2:66]-btmaster[2:66,i],yfit=fit)
       xyouts, 400, 170, string(res[0],format='(f8.4)'), col=70
       xyouts, 400, 150, string(res[1],format='(f8.4)'), col=70
       oplot, lc[2:66], fit, col=70
;stop

       if do_newdat then begin
           ofile = dir+'dx11_XO_yr_1-2_IQUmap_'+freq[i]+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_mastercl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000.newdat'
           spawn, 'cp '+dir+cleanfiles[i]+' '+ofile
           for ibin=1,100 do begin
               readcol, ofile, ib, cl, er, er, off, lmn, lmx, tag, format='i,f,f,f,f,i,i,a', skipline=12+ibin, numline=1, /silent
;##           print, ib, cl, er, er, off, lmn, lmx, tag
               oline = '' + $
                 string(ib,format='(i12)') + $
                 string(btmaster[ibin-1,i],format='(e12.4)') + $
                 string(er,format='(e12.4)') + $
                 string(er,format='(e12.4)') + $
                 string(off,format='(e12.4)') + $
                 string(lmn,format='(i12)') + $
                 string(lmx,format='(i12)') + $
                 ' '+tag
;##           print, oline
;##           print, 'sed "'+strtrim(string(ibin+13),2)+'c\ '+oline+'" -i '+ofile
;##           stop
               spawn, 'sed "'+strtrim(string(ibin+13),2)+'c\ '+oline+'" -i '+ofile
           endfor
       endif

       if do_xnewdat then begin
           ofile = dir+'dx11_XO_yr_1-2_IQUmap_'+freq[i]+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_masterXcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000.newdat'
           spawn, 'cp '+dir+cleanfiles[i]+' '+ofile
           for ibin=1,100 do begin
               readcol, ofile, ib, cl, er, er, off, lmn, lmx, tag, format='i,f,f,f,f,i,i,a', skipline=12+ibin, numline=1, /silent
;##           print, ib, cl, er, er, off, lmn, lmx, tag
               oline = '' + $
                 string(ib,format='(i12)') + $
                 string(bxtmaster[ibin-1,i],format='(e12.4)') + $
                 string(er,format='(e12.4)') + $
                 string(er,format='(e12.4)') + $
                 string(off,format='(e12.4)') + $
                 string(lmn,format='(i12)') + $
                 string(lmx,format='(i12)') + $
                 ' '+tag
;##           print, oline
;##           print, 'sed "'+strtrim(string(ibin+13),2)+'c\ '+oline+'" -i '+ofile
;##           stop
               spawn, 'sed "'+strtrim(string(ibin+13),2)+'c\ '+oline+'" -i '+ofile
           endfor
       endif
   endfor


   window, 4, xsize=720*scaling, ysize=450*2*scaling
   plot, lc, tcl, psym=4, xr=[0,2000]
   oplot, lc, tcl, psym=-4, col=245
   for i=0,2 do oplot, lc, btmaster[*,i], col=50*i, psym=-4
   for i=0,2 do oplot, lc, btt[*,i], col=50*i, psym=-5

   plot, lc, lc*0., line=2, yr=[-160,210], xr=[0,2000]
   for i=0,2 do oplot, lc, btmaster[*,i]-tcl, col=50*i, psym=4
   for i=0,2 do oplot, lc*1.01, btt[*,i]-tcl, col=50*i, line=4
   for i=2,0,-1 do for j=i-1,0,-1 do oplot, lc, btmaster[*,i]-btmaster[*,j], psym=-4, col=240-50*(i+j)
   for i=2,0,-1 do for j=i-1,0,-1 do xyouts, 1500-400*(i+j), -150, freq[i]+'-'+freq[j],col=240-50*(i+j)
   
   !p.multi=0
   window, 5
   indx = 42
   plot, lc, bxtmaster[*,0]/btmaster[*,0], psym=-4, yr=[0.99,1.01], chars=1.5
   oplot, lc, lc*0+1, line=2
   for i=0,2 do  oplot, lc, bxtmaster[*,i]/btmaster[*,i], psym=-4, col=50*i
   for i=0,2 do  oplot, lc, lc*0+mean(bxtmaster[0:indx,i]/btmaster[0:indx,i]), line=4, col=50*i
   for i=0,2 do  xyouts, 200+i*500, 1.0085, freq[i], col=50*i, chars=1.5
   for i=0,2 do  xyouts, 200+i*500, 1.007, string(mean(bxtmaster[0:indx,i]/btmaster[0:indx,i]),format='(f9.6)'), col=50*i, chars=1.5
   xyouts, 200, 1.005, '!6Recalibration up to l='+string(lc[indx]), chars=1.25

stop
   r = [60,80,375]



;   for i=0,2 do compare_xfaster_spectra, dir+files[i], dir+cleanfiles[i], /init, binning='data/bins/const/const2', lmax=3000, leg_tags=['raw','clean']+'-'+freq[i], reso=r[i], otit='DX11_XO: raw maps vs clean maps', win=i+10, /dops, psout='dx11xo_cleaning_'+freq[i]

;   for i=2,0,-1 do for j=i-1,0,-1 do compare_xfaster_spectra, dir+files[i], dir+files[j], /init, binning='data/bins/const/const2', lmax=3000, leg_tags=['raw','raw']+'-'+[freq[i], freq[j]], reso=r[i], otit='DX11_XO: raw maps comparison', win=i^2+j+1, /dops, psout='dx11xo_raw_'+freq[i]+'_'+freq[j]


;   for i=2,0,-1 do for j=i-1,0,-1 do compare_xfaster_spectra, dir+cleanfiles[i], dir+cleanfiles[j], /init, binning='data/bins/const/const2', lmax=3000, leg_tags=['clean','clean']+'-'+[freq[i], freq[j]], reso=180, otit='DX11_XO: clean maps comparison', win=i^2+j+1, /dops, psout='dx11xo_clean_'+freq[i]+'_'+freq[j]


end
