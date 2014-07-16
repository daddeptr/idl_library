
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

   dx_files = dfiles + '/' + dfiles +'.newdat'

   pdfiles = [ $
            'dx11_pre_year_1-2_Imap_100_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000', $
            'dx11_pre_year_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000', $
            'dx11_pre_year_1-2_Imap_217_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000', $
;#
            'dx11_pre_ds_1-2_Imap_100_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000', $
            'dx11_pre_ds_1-2_Imap_143_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000', $
            'dx11_pre_ds_1-2_Imap_217_full_extMask_545_undusted_insideM_split_flat_ncl4_v1.8_dx11_pre-beam_ns2048_uK_hrhs_xfcl_l3000_Tmask_X_Pmask_l4000' $
           ]
   
   pdx_files = pdfiles + '/' + pdfiles +'.newdat'

;print, dx_files
;print, ''
;print, pdx_files
;stop

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

   lstr = [500, 1000, 1250]

   !p.multi=[0,4,1]

   for i=0,2 do begin
       if not dops then window, 1, xsize=720*1.8, ysize=450*1.2
       if dops then begin
           set_plot, 'ps'
           device, file='dx11p-dx11_xf-spectra_comparison2theory-'+ftags[i]+'.eps', /col, /landscape, bits=8
       endif

       indx = i

;## --- TT
       dfile = dir+dx_files[indx+3]
       print, dfile
       dbclyrtt = extract_xfaster_newdat( dfile, lcen=bltt, btcl=tbcltt )
       dbclyree = extract_xfaster_newdat( dfile, lcen=blee, btcl=tbclee, ncl=2 )
       dbclyrte = extract_xfaster_newdat( dfile, lcen=blte, btcl=tbclte, ncl=4 )

       pdfile = dir+pdx_files[indx]
       print, pdfile
       pdbclyrtt = extract_xfaster_newdat( pdfile )
       pdbclyree = extract_xfaster_newdat( pdfile, ncl=2 )
       pdbclyrte = extract_xfaster_newdat( pdfile, ncl=4 )

       dfile = dir+dx_files[indx+6]
       print, dfile
       dbcldstt = extract_xfaster_newdat( dfile )
       dbcldsee = extract_xfaster_newdat( dfile, ncl=2 )
       dbcldste = extract_xfaster_newdat( dfile, ncl=4 )

       pdfile = dir+pdx_files[indx+3]
       print, pdfile
       pdbcldstt = extract_xfaster_newdat( pdfile )
       pdbcldsee = extract_xfaster_newdat( pdfile, ncl=2 )
       pdbcldste = extract_xfaster_newdat( pdfile, ncl=4 )

;## --- Residuals

       plot, bltt, dbclyrtt-tbcltt, psym=4, xtit='!6l', ytit='!6TT: Data - Model!d2013!n [!7l!6K!u2!n]', chars=1.5, xr=[lmn,lmx[i]], position=[0.06,0.075,0.325,0.90], ys=1, yr=[-265,265]
       oplot, bltt, dbclyrtt-tbcltt, psym=-4, col=70
;       oplot, bltt, dbclyrtt-tbcltt, col=70
       oplot, bltt, pdbclyrtt-tbcltt, psym=-5, col=90
;       oplot, bltt, pdbclyrtt-tbcltt, col=90

       oplot, bltt, dbclyrtt*0, col=245, thick=2

       oplot, bltt, dbcldstt-tbcltt, psym=-4, col=210
;       oplot, bltt, dbcldstt-tbcltt, col=210
       oplot, bltt, pdbcldstt-tbcltt, psym=-5, col=200
;       oplot, bltt, pdbcldstt-tbcltt, col=200

       plot, /nodata, [0,1], [0,1], position=[0.06,0.92,0.94,0.99], ys=1, yr=[-265,265], col=255
       xyouts, 0., 0.7, 'DX11: '+run_tags[indx+3], col=70, chars=1.25
       xyouts, 0.25, 0.7, 'preDX11: '+run_tags[indx+3], col=90, chars=1.25
       xyouts, 0.5, 0.7, 'DX11: '+run_tags[indx+6], col=210, chars=1.25
       xyouts, 0.75, 0.7, 'preDX11: '+run_tags[indx+6], col=200, chars=1.25

;stop
;       readcol,'/global/scratch2/sd/dpietrob/Software/XFaster/data/beams/dx11/dx11_wl_100.dat', l, dwl
;       readcol,'/global/scratch2/sd/dpietrob/Software/XFaster/data/beams/dx11_pre/dx11_pre_wl_100.dat', l, pdwl
;       oplot, l, (pdwl/dwl)^2, col=70

;## --- EE

       plot, blee, dbclyree-tbclee, psym=4, xtit='!6l', ytit='!6EE: Data - Model!d2013!n [!7l!6K!u2!n]', chars=1.5, xr=[lmn,lmx[i]], position=[0.3875,0.075,0.6525,0.90], ys=1, yr=[-5,5], col=0
       oplot, blee, dbclyree-tbclee, psym=-4, col=70
;       oplot, blee, dbclyree-tbclee, col=70

       oplot, blee, pdbclyree-tbclee, psym=-5, col=90
;       oplot, blee, pdbclyree-tbclee, col=90

       oplot, blee, dbclyree*0, col=245, thick=2

       oplot, blee, dbcldsee-tbclee, psym=-4, col=210
;       oplot, blee, dbcldsee-tbclee, col=210

       oplot, blee, pdbcldsee-tbclee, psym=-5, col=200
;       oplot, blee, pdbcldsee-tbclee, col=200

;       xyouts, lstr[indx], 400, run_tags[indx+3]+': DX11-preDX11', col=70
;       xyouts, lstr[indx], 350, run_tags[indx+6]+': DX11-preDX11', col=210


;stop


;## --- TE

       plot, blte, dbclyrte-tbclte, psym=4, xtit='!6l', ytit='!6TE: Data - Model!d2013!n [!7l!6K!u2!n]', chars=1.5, xr=[lmn,lmx[i]], position=[0.715,0.075,0.98,0.90], ys=1, yr=[-15,15]
       oplot, blte, dbclyrte-tbclte, psym=-4, col=70
;       oplot, blte, dbclyrte-tbclte, col=70

       oplot, blte, pdbclyrte-tbclte, psym=-5, col=90
;       oplot, blte, pdbclyrte-tbclte, col=90
       oplot, blte, dbclyrte*0, col=245, thick=2

       oplot, blte, dbcldste-tbclte, psym=-4, col=210
;       oplot, blte, dbcldste-tbclte, col=210

       oplot, blte, pdbcldste-tbclte, psym=-5, col=200
;       oplot, blte, pdbcldste-tbclte, col=200

;       xyouts, lstr[indx], 400, run_tags[indx+3]+': DX11-preDX11', col=70
;       xyouts, lstr[indx], 350, run_tags[indx+6]+': DX11-preDX11', col=210


;stop

       if dops then begin
           device, /close
           set_plot, 'x'
       endif

   endfor

   stop, ' --- End of Script ---'

end



