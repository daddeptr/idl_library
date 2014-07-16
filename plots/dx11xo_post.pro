   True = 1b
   False = 0b

   do_newdat = False
   do_xnewdat = False
   do_anafast = False

   dir = '/global/homes/d/dpietrob/myscratch/Software/XFaster/outputs/dx11_XO_co-clean_maps/'
   freq = ['100','143','217']
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
