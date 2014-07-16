
   True = 1b
   False = 0b

   init = True

   if init then begin
       mollview, findgen(12l*8^2), win=1, px=450
       loadct, 39
       !p.color=0
       !p.background=255
   endif

   dops = True

   dir = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/'

   dfiles = [ $
            'dx11_hr_1-2_IQUmap_100_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2000_Tmask_x_Pmask_l4000', $
            'dx11_hr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2500_Tmask_x_Pmask_l4000', $
            'dx11_hr_1-2_IQUmap_217_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l3000_Tmask_x_Pmask_l4000', $
;#
            'dx11_yr_1-2_IQUmap_100_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2000_Tmask_x_Pmask_l4000', $
            'dx11_yr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2500_Tmask_x_Pmask_l4000', $
            'dx11_yr_1-2_IQUmap_217_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l3000_Tmask_x_Pmask_l4000', $
;#
            'dx11_ds_1-2_IQUmap_100_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2000_Tmask_x_Pmask_l4000', $
            'dx11_ds_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2500_Tmask_x_Pmask_l4000', $
            'dx11_ds_1-2_IQUmap_217_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l3000_Tmask_x_Pmask_l4000' $
           ]

   dx_files = dfiles + '/' + dfiles + '.newdat'

   pdfiles = [ $
            'dx11_pre_year_1-2_Imap_100_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000', $
            'dx11_pre_year_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000', $
            'dx11_pre_year_1-2_Imap_217_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000', $
;#
            'dx11_pre_ds_1-2_Imap_100_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000', $
            'dx11_pre_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000', $
            'dx11_pre_ds_1-2_Imap_217_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000' $
           ] ;+ '.newdat'
   
   pdx_files = pdfiles + '/' + pdfiles + '.newdat'

   dx_beams = '/global/scratch2/sd/dpietrob/Software/XFaster/data/beams/dx11/dx11_wl_'+ $
     ['100','143','217','100','143','217','100_ds_1-2','143_ds_1-2', '217_ds_1-2'] + '.dat'

   pdx_beams = '/global/scratch2/sd/dpietrob/Software/XFaster/data/beams/dx11_pre/dx11_pre_wl_'+ $
     ['100','143','217','100_ds_1-2','143_ds_1-2', '217_ds_1-2'] + '.dat'

   nfiles = n_elements(dfiles)

   run_tags = strarr(nfiles)
   run_tags[0] = 'hr-100GHz';: clean-IN'
   run_tags[1] = 'hr-143GHz';: clean-IN'
   run_tags[2] = 'hr-217GHz';: clean-IN'

   run_tags[3] = 'yr-100GHz';: clean-IN'
   run_tags[4] = 'yr-143GHz';: clean-IN'
   run_tags[5] = 'yr-217GHz';: clean-IN'

   run_tags[6] = 'ds-100GHz';: clean-IN'
   run_tags[7] = 'ds-143GHz';: clean-IN'
   run_tags[8] = 'ds-217GHz';: clean-IN'

   ftags = ['100','143','217']
   lmx = [2000,2500,3000]
   lmn = 1
   resolution = 1 

   xfdir = '/global/scratch2/sd/dpietrob/Software/XFaster/'
   fits2cl, tcl, xfdir+'data/planck_lcdm_cl_uK_xf1.e-3.fits'
   tcl[*,4] = 0.
   tcl[*,5] = 0.

   l = findgen(n_elements(tcl[*,0]))
   ll = l*(l+1)/2./!pi
   for i=0,5 do tcl[*,i] = tcl[*,i] * ll

   lstr = [800, 1000, 1250]*0+250

   !p.multi=[0,3,1]
;   !p.multi = 0

   if not dops then window, 1, xsize=720*1.8, ysize=450*1.2
   if dops then begin
       set_plot, 'ps'
       device, file='dx11p-dx11_xf-spectra_ratio.eps', /col, /landscape, bits=8
   endif

   calf = [1.013, 1.008, 1.015]
   for i=0,2 do begin

       indx = i

;## --- TT
       dfile = dir+dx_files[indx+3]
       print, dfile
       dbclyrtt = extract_xfaster_newdat( dfile, lcen=bltt, btcl=tbcltt )
;       dbclyree = extract_xfaster_newdat( dfile, lcen=blee, btcl=tbclee, ncl=2 )
;       dbclyrte = extract_xfaster_newdat( dfile, lcen=blte, btcl=tbclte, ncl=4 )

       pdfile = dir+pdx_files[indx]
       print, pdfile
       pdbclyrtt = extract_xfaster_newdat( pdfile )
;       pdbclyree = extract_xfaster_newdat( pdfile, ncl=2 )
;       pdbclyrte = extract_xfaster_newdat( pdfile, ncl=4 )

       dfile = dir+dx_files[indx+6]
       print, dfile
       dbcldstt = extract_xfaster_newdat( dfile )
;       dbcldsee = extract_xfaster_newdat( dfile, ncl=2 )
;       dbcldste = extract_xfaster_newdat( dfile, ncl=4 )

       pdfile = dir+pdx_files[indx+3]
       print, pdfile
       pdbcldstt = extract_xfaster_newdat( pdfile )
;       pdbcldsee = extract_xfaster_newdat( pdfile, ncl=2 )
;       pdbcldste = extract_xfaster_newdat( pdfile, ncl=4 )

;## --- Residuals

       plot, bltt, dbclyrtt/pdbclyrtt, psym=4, xtit='!6l', ytit='!6Ratio DX11-preDX11 (TT)', chars=1.5, xr=[lmn,lmx[i]], ys=1, yr=[0.95,1.05];, position=[0.06,0.075,0.325,0.94]
       oplot, bltt, dbclyrtt/pdbclyrtt, psym=4, col=70
       oplot, bltt, dbclyrtt/pdbclyrtt, col=70
       oplot, bltt, dbclyrtt*0+1, col=245, line=2
       oplot, bltt, dbclyrtt*0+calf[i], line=3

       oplot, bltt, dbcldstt/pdbcldstt, psym=5, col=210
       oplot, bltt, dbcldstt/pdbcldstt, col=210

       xyouts, lstr[indx], 1.04, run_tags[indx+3], col=70
       xyouts, lstr[indx], 1.035, run_tags[indx+6], col=210
       xyouts, lstr[indx], 0.98, 'Recalibration: '+string(calf[i],format='(f5.3)')

       readcol, dx_beams[indx+3], l, dwlyr
       readcol, pdx_beams[indx], l, pdwlyr
       
       byr = bp_binning( (pdwlyr/dwlyr)^2, '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT')

       oplot, bltt, byr, col=70

       readcol, dx_beams[indx+6], l, dwlds
       readcol, pdx_beams[indx+3], l, pdwlds
       bds = bp_binning( (pdwlds/dwlds)^2, '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/CTP_bin_TT' )

       oplot, bltt, bds, col=210

;## --- EE
       if False then begin

       plot, blee, dbclyree-pdbclyree, psym=4, xtit='!6l', ytit='!6Difference DX11-preDX11 [!7l!6K!u2!n]', chars=1.5, xr=[lmn,lmx[i]], position=[0.3875,0.075,0.6525,0.94], ys=1, yr=[-1.25,3]
       oplot, blee, dbclyree-pdbclyree, psym=4, col=70
       oplot, blee, dbclyree-pdbclyree, col=70
       oplot, blee, dbclyree*0, col=245, thick=2

       oplot, blee, dbcldsee-pdbcldsee, psym=5, col=210
       oplot, blee, dbcldsee-pdbcldsee, col=210

       xyouts, lstr[indx], 400, run_tags[indx+3]+': DX11-preDX11', col=70
       xyouts, lstr[indx], 350, run_tags[indx+6]+': DX11-preDX11', col=210

;## --- TE

       plot, blte, dbclyrte-pdbclyrte, psym=4, xtit='!6l', ytit='!6Difference DX11-preDX11 [!7l!6K!u2!n]', chars=1.5, xr=[lmn,lmx[i]], position=[0.715,0.075,0.98,0.94], ys=1, yr=[-10,10]
       oplot, blte, dbclyrte-pdbclyrte, psym=4, col=70
       oplot, blte, dbclyrte-pdbclyrte, col=70
       oplot, blte, dbclyrte*0, col=245, thick=2

       oplot, blte, dbcldste-pdbcldste, psym=5, col=210
       oplot, blte, dbcldste-pdbcldste, col=210

       xyouts, lstr[indx], 400, run_tags[indx+3]+': DX11-preDX11', col=70
       xyouts, lstr[indx], 350, run_tags[indx+6]+': DX11-preDX11', col=210

endif
;stop

   endfor

   if dops then begin
       device, /close
       set_plot, 'x'
   endif


   stop, ' --- End of Script ---'

end



