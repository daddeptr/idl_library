True = 1b
False = 0b

tmaskf = 'data/mask/dx11/cons-apo_Tmask3.fits'
pmaskf = 'data/mask/dx11/cons-apo_Pmask3.fits'

freq = ['100','143','217']
lmx = [2000,2500,3000]
nmrt = 'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'
fnm= 'outputs/'+nmrt+'/'+nmrt+'.newdat'

for ifreq = 0, 2 do begin
;##ort = 'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_co_undusted_split_ns2048_uK_hrhs_const2_xfcl_l3000_cons-apo_Tmask3_x_Pmask_l4000'
;##nrt = 'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l3000_cons-apo_Tmask3_x_Pmask_l4000'
;##vnrt = 'dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs_const2_xfcl_l3000_cons-apo_Tmask3_x_cons-apo_Pmask3_l4000'

    plot_xfaster_newdat, fnm[ifreq], /pol, /init, /resi, win=ifreq+1, binning='data/bins/const/const2', lmax=lmx[ifreq]
    for jfreq=ifreq+1,2 do begin
        compare_xfaster_spectra, fnm[ifreq], fnm[jfreq], win=4, reso=75, binning='data/bins/const/const2'
        compare_xfaster_spectra, fnm[ifreq], fnm[jfreq], win=5, ncl=2, reso=10, binning='data/bins/const/const2'
stop
    endfor
stop
endfor
stop
;mymoll, 'data/maps/dx11/dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_co_undusted_split_ns2048_uK_hrhs.fits', $
;  file2='data/maps/dx11/dx11_yr_1-2_IQUmap_'+freq+'_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs.fits', min=-2, max=2, maskfile=pmaskf, imap=3, win=1, tit='P difference', dosmooth=sqrt(300.^2-7.^2)

;mymoll, 'data/maps/dx11/dx11_yr_1-2_IQUmap_217_full_extMask_545_co_undusted_split_ns2048_uK_hrhs.fits', $
;  file2='data/maps/dx11/dx11_yr_1-2_IQUmap_217_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs.fits', min=-15, max=15, maskfile=pmaskf, imap=3, win=1, tit='P difference', dosmooth=sqrt(300.^2-5.^2)

;mymoll, 'data/maps/dx11/dx11_yr_1-2_IQUmap_217_full_extMask_545_co_undusted_split_ns2048_uK_hrhs.fits', min=0, max=150, maskfile=pmaskf, imap=3, win=2, tit='T difference', dosmooth=sqrt(300.^2-5.^2)
;mymoll, 'data/maps/dx11/dx11_yr_1-2_IQUmap_217_full_extMask_545_coTP_undusted_split_ns2048_uK_hrhs.fits', min=0, max=150, maskfile=pmaskf, imap=3, win=3, tit='TP difference', dosmooth=sqrt(300.^2-5.^2)

;fnm='outputs/'+nmrt+'/'+nmrt+'.newdat'
;fn='outputs/'+nrt+'/'+nrt+'.newdat'
;fvn='outputs/'+vnrt+'/'+vnrt+'.newdat'
;fo='outputs/'+ort+'/'+ort+'.newdat'

;##compare_xfaster_spectra, fnm, fo, /init, win=1, reso=75, lmax=3000
;compare_xfaster_spectra, fnm, fo, /init, win=2, ncl=2, reso=15, lmax=3000, leg_tags=['co-apoPmask','co-Pmask']
;compare_xfaster_spectra, fnm, fo, /init, win=3, ncl=4, reso=20, lmax=3000, leg_tags=['co-apoPmask','co-Pmask']

;compare_xfaster_spectra, fn, fvn, /init, win=3, ncl=1, reso=50, lmax=3000, leg_tags=['coTP','coTP Pmask']
;compare_xfaster_spectra, fn, fvn, /init, win=4, ncl=2, reso=10, lmax=3000, leg_tags=['coTP','coTP Pmask']
;compare_xfaster_spectra, fn, fvn, /init, win=5, ncl=4, reso=20, lmax=3000, leg_tags=['coTP','coTP Pmask']


end
