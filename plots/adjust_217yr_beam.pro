   d1 = extract_xfaster_newdata_output('../outputs/dx11_pre_ds_1-2_Imap_217_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', ncl=1)
   d2 = extract_xfaster_newdata_output('../outputs/dx11_pre_year_1-2_Imap_217_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000.newdat', ncl=1, lcen=lc)

   mollview, findgen(12), win=0, px=500
   loadct, 39
   !p.color=0
   !p.background=255

   set_plot, 'ps'
   device, file='yr217beam_correction.eps', /col, bits=8
   plot, lc, d2/d1, chars=1.5, yr=[0.95,1.05], /xlog, xtit='!8l', ytit='!6'
   oplot, gaussbeam(0.5,3000)^2, col=245
   xyouts, 15, 1.05, '!6C!dl!u217!dyr!n /C!dl!u217!dds!n', chars=1.5
   xyouts, 15, 1.04, '!6Gaussian beam 0.5 arcmin', chars=1.5, col=245
   device, /close
   set_plot, 'x'

end
