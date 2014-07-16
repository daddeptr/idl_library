True = 1b
False = 0b

rt = [ $
      'dx11_yr_1-2_IQUmap_100_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2000_Tmask_x_Pmask_l4000', $
      'dx11_yr_1-2_IQUmap_143_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l2500_Tmask_x_Pmask_l4000', $
      'dx11_yr_1-2_IQUmap_217_full_extMask_545_undusted_insideM_split_ns2048_uK_hrhs_xfcl_l3000_Tmask_x_Pmask_l4000' $
]

dir = '/global/scratch2/sd/dpietrob/Software/XFaster/outputs/'

cmb = fltarr(167,3)
bcl = fltarr(167,3)
fgbcl = fltarr(167,3)

mollview, findgen(12), px=400, win=10
!p.multi= 0
loadct, 39
!p.color=0
!p.background=255

!p.multi=[0,2,1]
window, 10
cols=[0,70,210]
plot_io, /nodata, [1,3250], [10,10000], chars=1.5, thick=2, xs=1, ys=1, xtit='!8l', ytit='!8D!dl!n [!7l!8K!u2!n]'
xyouts, 2000,7000, '!6Planck 2013', col=245
xyouts, 2000,4500, '!6DX11 100', col=cols[0]
xyouts, 2000,3000, '!6DX11 143', col=cols[1]
xyouts, 2000,2000, '!6DX11 217', col=cols[2]

for i=0,2 do begin
    file = dir+rt[i]+'/'+rt[i]+'.newdat'
    tt = extract_xfaster_newdat( file, lcen=cl, btcl=t13 ) 
    oplot, cl, tt, col=cols[i], line=3
    oplot, cl, t13, col=245, thick=2

    readcol, dir+rt[i]+'/'+'fg_bestfit.txt', l, fgtt
    xfgtt = reform( [0,0,[fgtt]] )
    xl = reform( [0,0,[l]] )
    bfgtt = bp_binning(xfgtt, 'data/bins/ctp/CTP_bin_TT')
    bl = bp_binning(xl, 'data/bins/ctp/CTP_bin_TT')

    readcol, dir+rt[i]+'/'+'bestfit.txt', l, bftt
    xtt = reform( [0,0,[bftt]] )
;    xl = reform( [0,0,[l]] )
    bbftt = bp_binning(xtt, 'data/bins/ctp/CTP_bin_TT')
;    bl = bp_binning(xl, 'data/bins/ctp/CTP_bin')

    ctt = bbftt - bfgtt
    oplot, cl, ctt, thick=2, col=cols[i]
    oplot, cl, bfgtt, line=2, col=cols[i]

    cmb[*,i] = ctt[0:166]
    fgbcl[*,i] = bfgtt[0:166]
    bcl[*,i] = bbftt[0:166]
endfor

plot_oo, /nodata, [1,3250], [10,10000], chars=1.5, thick=2, xs=1, ys=1, xtit='!8l', ytit='!8D!dl!n [!7l!8K!u2!n]'
for i=0,2 do begin
    file = dir+rt[i]+'/'+rt[i]+'.newdat'
    tt = extract_xfaster_newdat( file, lcen=cl, btcl=t13 ) 
    oplot, cl, tt, col=cols[i], line=3
    oplot, cl, t13, col=245, thick=2

    readcol, dir+rt[i]+'/'+'fg_bestfit.txt', l, fgtt
    xfgtt = reform( [0,0,[fgtt]] )
    xl = reform( [0,0,[l]] )
    bfgtt = bp_binning(xfgtt, 'data/bins/ctp/CTP_bin_TT')
    bl = bp_binning(xl, 'data/bins/ctp/CTP_bin_TT')

    readcol, dir+rt[i]+'/'+'bestfit.txt', l, bftt
    xtt = reform( [0,0,[bftt]] )
;    xl = reform( [0,0,[l]] )
    bbftt = bp_binning(xtt, 'data/bins/ctp/CTP_bin_TT')
;    bl = bp_binning(xl, 'data/bins/ctp/CTP_bin')

    ctt = bbftt - bfgtt
    oplot, cl, ctt, thick=2, col=cols[i]
    oplot, cl, bfgtt, line=2, col=cols[i]

    cmb[*,i] = ctt[0:166]
    fgbcl[*,i] = bfgtt[0:166]
    bcl[*,i] = bbftt[0:166]
endfor

!p.multi=0
window, 12
plot, cl, cmb[*,0]/t13, chars=1.5, thick=2, xs=1, ys=1, xtit='!8l', ytit='!8D!dl!n / Bf!u2013!n', xr=[1,1550], yr=[0.85, 1.15]
oplot, cl, cmb[*,1]/t13, thick=2, col=cols[1]
oplot, cl, cmb[*,2]/t13, thick=2, col=cols[2]
oplot, cl, cmb[*,2]*0.+1., thick=2, col=245, line=2
xyouts, 800, 0.91, '!6DX11 100/Bfit!u2013!n', chars=1.5
xyouts, 800, 0.885, '!6DX11 143/Bfit!u2013!n', col=70, chars=1.5
xyouts, 800, 0.86, '!6DX11 217/Bfit!u2013!n', col=210, chars=1.5

!p.multi=0
window, 11
plot, cl, cmb[*,0]/cmb[*,2], chars=1.5, thick=2, xs=1, ys=1, xtit='!8l', ytit='!8D!dl!n / D!dl!u217!n', xr=[1,1550], yr=[0.85, 1.15]
oplot, cl, cmb[*,1]/cmb[*,2], thick=2, col=90
oplot, cl, cmb[*,2]*0.+1., thick=2, col=245, line=2
xyouts, 800, 0.925, '!6DX11 100/217', chars=1.5
xyouts, 800, 0.875, '!6DX11 143/217', col=90, chars=1.5
end
