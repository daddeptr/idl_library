
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

   dir = [ $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
;#
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
;#
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', $
           '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/' $
             ]
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

   run_tags = 'DX11 '+run_tags

   ftags = ['100','143','217']
   lmx = [2000,2500,3000]
   lmn = 1
   resolution = 1 
   !p.multi=[0,3,2]

   xfdir = '/global/scratch2/sd/dpietrob/Software/XFaster/'
   fits2cl, tcl, xfdir+'data/planck_lcdm_cl_uK_xf1.e-3.fits'
   tcl[*,4] = 0.
   tcl[*,5] = 0.

   l = findgen(n_elements(tcl[*,0]))
   ll = l*(l+1)/2./!pi
   for i=0,5 do tcl[*,i] = tcl[*,i] * ll

   lstr = [800, 1000, 1250]

   for i=0,2 do begin
       if not dops then window, 1, xsize=720*1.8, ysize=450*1.8
       if dops then begin
           set_plot, 'ps'
           device, file='dx11_xf-spectra-'+ftags[i]+'.eps', /col, /landscape, bits=8
       endif

       indx = i
;## --- TT
       file = dir[indx]+dfiles[indx]+'/'+dfiles[indx]+'.newdat'
       bclhrtt = extract_xfaster_newdat( file, lcen=bltt, btcl=tbcltt )
       plot, l, tcl[*,0], xr=[lmn,lmx[i]], yr=[-500,6500], ys=1, position=[0.06,0.45,0.325,0.975], xtickname=strarr(6)+' ', ytit='!8D!dl!uTT!n !6[!7l!6K!u2!n]', chars=1.5
       oplot, l, tcl[*,0], col=245, thick=2
       oplot, bltt, bclhrtt, psym=5
       xyouts, lstr[i], 6000, run_tags[indx]
       xyouts, lstr[i], 5500, run_tags[indx+3], col=70
       xyouts, lstr[i], 5000, run_tags[indx+6], col=210

       file = dir[indx+3]+dfiles[indx+3]+'/'+dfiles[indx+3]+'.newdat'
       bclyrtt = extract_xfaster_newdat( file )
       oplot, bltt, bclyrtt, psym=5, col=70

       file = dir[indx+6]+dfiles[indx+6]+'/'+dfiles[indx+6]+'.newdat'
       bcldstt = extract_xfaster_newdat( file )
       oplot, bltt, bcldstt, psym=5, col=210

;## --- EE
       file = dir[indx]+dfiles[indx]+'/'+dfiles[indx]+'.newdat'
       bclhree = extract_xfaster_newdat( file, lcen=blee, btcl=tbclee, ncl=2 )
       plot, l, tcl[*,1], xr=[lmn,lmx[i]], yr=[-5,50], ys=1, position=[0.3875,0.45,0.6525,0.975], xtickname=strarr(6)+' ', ytit='!8D!dl!uEE!n !6[!7l!6K!u2!n]', chars=1.5
       oplot, l, tcl[*,1], col=245, thick=2
       oplot, blee, bclhree, psym=5
       oplot, blee, bclhree

       file = dir[indx+3]+dfiles[indx+3]+'/'+dfiles[indx+3]+'.newdat'
       bclyree = extract_xfaster_newdat( file, ncl=2 )
       oplot, blee, bclyree, psym=5, col=70
       oplot, blee, bclyree, col=70

       file = dir[indx+6]+dfiles[indx+6]+'/'+dfiles[indx+6]+'.newdat'
       bcldsee = extract_xfaster_newdat( file, ncl=2 )
       oplot, blee, bcldsee, psym=5, col=210
       oplot, blee, bcldsee, col=210

;## --- TE
       file = dir[indx]+dfiles[indx]+'/'+dfiles[indx]+'.newdat'
       bclhrte = extract_xfaster_newdat( file, lcen=blte, btcl=tbclte, ncl=4 )
       plot, l, tcl[*,3], xr=[lmn,lmx[i]], yr=[-150,150], ys=1, position=[0.715,0.45,0.98,0.975], xtickname=strarr(6)+' ', ytit='!8D!dl!uTE!n !6[!7l!6K!u2!n]', chars=1.5
       oplot, l, tcl[*,3], col=245, thick=2
       oplot, blte, bclhrte, psym=5
       oplot, blte, bclhrte

       file = dir[indx+3]+dfiles[indx+3]+'/'+dfiles[indx+3]+'.newdat'
       bclyrte = extract_xfaster_newdat( file, ncl=4 )
       oplot, blte, bclyrte, psym=5, col=70
       oplot, blte, bclyrte, col=70

       file = dir[indx+6]+dfiles[indx+6]+'/'+dfiles[indx+6]+'.newdat'
       bcldste = extract_xfaster_newdat( file, ncl=4 )
       oplot, blte, bcldste, psym=5, col=210
       oplot, blte, bcldste, col=210

;## --- Residuals

       plot, bltt, bclhrtt-tbcltt, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=1.5, yr=[-250,250]*resolution, xr=[lmn,lmx[i]], position=[0.06,0.075,0.325,0.44]
       oplot, bltt, bclhrtt*0, col=245, thick=2
       oplot, bltt, bclhrtt-tbcltt
       oplot, bltt, bclyrtt-tbcltt, col=70, psym=4
       oplot, bltt, bclyrtt-tbcltt, col=70
       oplot, bltt, bcldstt-tbcltt, col=210, psym=4
       oplot, bltt, bcldstt-tbcltt, col=210

       plot, blee, bclhree-tbclee, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=1.5, xr=[lmn,lmx[i]], yr=[-5,5]*resolution, ys=1, position=[0.3875,0.075,0.6525,0.44]
       oplot, blee, bclhree*0, col=245, thick=2
       oplot, blee, bclhree-tbclee
       oplot, blee, bclyree-tbclee, col=70, psym=4
       oplot, blee, bclyree-tbclee, col=70
       oplot, blee, bcldsee-tbclee, col=210, psym=4
       oplot, blee, bcldsee-tbclee, col=210

       plot, blte, bclhrte-tbclte, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=1.5, xr=[lmn,lmx[i]], yr=[-20,20]*resolution, position=[0.715,0.075,0.98,0.44]
       oplot, blte, bclhrte*0, col=245, thick=2
       oplot, blte, bclhrte-tbclte
       oplot, blte, bclyrte-tbclte, col=70, psym=4
       oplot, blte, bclyrte-tbclte, col=70
       oplot, blte, bcldste-tbclte, col=210, psym=4
       oplot, blte, bcldste-tbclte, col=210
;stop

       if dops then begin
           device, /close
           set_plot, 'x'
       endif

   endfor

   stop, ' --- End of Script ---'

end



