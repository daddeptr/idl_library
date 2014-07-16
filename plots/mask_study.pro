freq = ['100','143','217']

spawn, 'ls outputs/dx11_XO_co-clean_maps/dx11_XO_yr_1-2_IQUmap_'+freq+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_xfcl_l3000_Tmask_x_Pmask_l4000.newdat', cof_oldM
spawn, 'ls outputs/dx11_XO_co-clean_maps/dx11_XO_yr_1-2_IQUmap_'+freq+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000.newdat', cof_newM
spawn, 'ls outputs/dx11_XO_co-clean_maps/dx11_XO_yr_1-2_IQUmap_'+freq+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_xfcl_l3000_cons-apo_Tmask3_mPix_x_cons-apo_Pmask3_mPix_l4000.newdat', cof_pixM
;##
spawn, 'ls outputs/dx11_XO_co-clean_maps/dx11_XO_yr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs_xfcl_l3000_Tmask_x_Pmask_l4000.newdat', raw_oldM
spawn, 'ls outputs/dx11_XO_co-clean_maps/dx11_XO_yr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000.newdat', raw_newM
spawn, 'ls outputs/dx11_XO_co-clean_maps/dx11_XO_yr_1-2_IQUmap_'+freq+'_full_ns2048_uK_hrhs_xfcl_l3000_cons-apo_Tmask3_mPix_x_cons-apo_Pmask3_mPix_l4000.newdat', raw_pixM

compare_xfaster_spectra, raw_pixM[0], raw_newM[0], binning='data/bins/const/const2', reso=350, lmax=3000, leg_tags=['raw pix mask','raw new mask'], win=2
compare_xfaster_spectra, raw_oldM[0], raw_newM[0], binning='data/bins/const/const2', reso=350, lmax=3000, leg_tags=['raw old mask','raw new mask'], win=1


end
