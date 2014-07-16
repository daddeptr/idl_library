
   True = 1b
   False = 0b

   if True then begin
       mollview, findgen(12l*4^2), win=2, px=550
       loadct, 39
       !p.color=0
       !p.background=255
   endif

   dospec = True
   dops   = True
   xlog   = False

   if dospec then begin
       outdir = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/'
       dfiles = outdir + [ $
                          'dx11_pre_ds_1-2_Imap_100_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_ds_1-2_Imap_217_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat' $
                        ]
                        
       yfiles = outdir + [ $
                          'dx11_pre_year_1-2_Imap_100_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_year_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', $
                          'dx11_pre_year_1-2_Imap_217_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat' $
                        ]
                        
       run_tags = [ $
                                 '100 dx11_pre Full dx11-beam ', $
                                 '143 dx11_pre Full dx11-beam ', $
                                 '217 dx11_pre Full dx11-beam ' $
                               ];+['clean-IN']

       print, run_tags
       print, yfiles, dfiles

       lmax=[1550,2050,2550]

       fileroot='detset_year_pse_diff_'+['100','143','217']

       for i=0,2 do begin
           compare_xfaster_spectra, dfiles[i], yfiles[i], leg_tags=['Detset 1-2', 'Year 1-2'], otit=run_tags[i], win=1, ncl=1, lmax=lmax[i], lmin=50, xlog=xlog, dops=dops, fileout=fileroot[i]
           compare_xfaster_spectra, dfiles[i], yfiles[i], leg_tags=['Detset 1-2', 'Year 1-2'], otit=run_tags[i], win=2, ncl=2, reso=2, lmax=lmax[i], lmin=50, xlog=xlog, dops=dops, fileout=fileroot[i]
           compare_xfaster_spectra, dfiles[i], yfiles[i], leg_tags=['Detset 1-2', 'Year 1-2'], otit=run_tags[i], win=4, ncl=4, reso=10, lmax=lmax[i], lmin=50, xlog=xlog, dops=dops, fileout=fileroot[i]
;           stop
       endfor


 ;      nfiles = n_elements( files )
       lmax = 3000
;       !p.multi = [0,3,3]
;       colors = 245/(nfiles-1)*lindgen(nfiles)
       colors = [ $
                  0, $
                  70, $
                  210, $
                  230 $
                ]
       psyms = [ $
                 1, $
                 2, $
                 4, $
                 5, $
                 6 $
                 ]

   endif
   stop, ' --- End of Script ---'

end



