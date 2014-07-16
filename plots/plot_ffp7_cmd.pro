   True = 1b
   False = 0b

   mollview, findgen(12), win=0, px=450
   loadct, 39
   !p.color=0
   !p.background=255

   dir = '/global/scratch2/sd/dpietrob/Software/XFaster/'

   fits2cl, cls, dir+'ffp7_cmb_000_cls.fits'

   cls = cls * 1.e12
   bcls = xf_binning(cls[0:3500], dir+'data/bins/ctp/CTP_bin_TT', lcen = lc)

   fits2cl, spicls, '/global/scratch2/sd/loris/FFP/7/Commander/temp/spice_cls_commander_ffp7_temp_n2048_05arc_v004_hr1_x_hr2_cmb.fits'
   bspicls = xf_binning(spicls[0:3500], dir+'data/bins/ctp/CTP_bin_TT' )

   cmdfile = dir+'outputs/commander_ffp7_temp_n2048_05arc_v004_hrs_cmb_ctp_xfcl_l3000_combined_60_x_combined_60_l4000.newdat'
;##   cmdfile = dir+'outputs/commander_ffp7_temp_n2048_05arc_v004_hrs_cmb_ctp_xfcl_l3000_temp_mask_fsky0.8_src70-353_2048_x_temp_mask_fsky0.8_src70-353_2048_l4000.newdat'
;##   cmdfile = dir+'outputs/commander_ffp7_temp_n2048_05arc_v004_hrs_cmb_ctp_xfcl_l2500_temp_mask_fsky0.8_src70-353_2048_x_temp_mask_fsky0.8_src70-353_2048_l4000.newdat'
   
   cmdcl = extract_xfaster_newdat( cmdfile, cler=er ) 

   window, 0, xsize=1200, ysize=650
   !p.multi = [0,2,1]
   plot, lc, bcls, chars=2, xtit='!8l', ytit='!8D!dl!n [!7l!^K!u2!n]', xr=[0,3000]
   oplot, lc, bspicls, psym=-6, col=75
   oplot, lc, cmdcl, psym=-4
   errplot, lc, cmdcl-er, cmdcl+er
   oplot, lc, bcls, thick=2, col=245

   plot, lc, cmdcl-bcls, chars=2, xtit='!8l', ytit='!7D!8!dl!n [!7l!^K!u2!n]', xr=[1,3000], psym=4, yr=[-1.e3, 1.e3], ys=1
   oplot, lc, bspicls-bcls, col=75, psym=6
   oplot, lc, cmdcl*0, line=2
;   errplot, lc, cmdcl-er, cmdcl+er
;   oplot, lc, bcls, thick=2, col=245

end
