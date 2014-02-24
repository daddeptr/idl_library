; Simple ILC code
function hfi_ilc, tag
;tag = 'delta_dx9'
if (keyword_set(tag)) then tag = 'ffp6'

true = 1b
false = 0b

do_png = false

freq = [100, 143, 217, 353]
;freq = [143, 217, 353]
sfreq = string(freq, format='(i3.3)')
nfreq = n_elements(freq)

cmb2sz = fltarr(nfreq)
h = 6.626068 * 10.^(-34)
k = 1.3806503 * 10.^(-23)
Tcmb = 2.725
x = (h*freq*10.^9) / (k * Tcmb)
cmb2sz = x * (exp(x)+1.)/(exp(x)-1.)-4.
;print, cmb2sz
;stop

ns = 128l ;1024l
npix = 12l*ns^2
maps = fltarr(npix, nfreq)

;read_fits_map, 'chi2_ps_mask_ns128_ring.fits', mask
;read_fits_map, 'ps_mask_ns128_ring.fits', mask
;read_fits_map,'../ns1024/sky_cut_30_ns1024_ring.fits', mask
;mask[*] = 1.
;## mfile = 'dust_ps_mask.fits'
;## mfile = 'ffp4_commomn_mask_ns128.fits'

mfile = '/global/homes/d/dpietrob/myscratch/dx9/maps/ns0128/dust_0.25mK_mask_ns0128.fits'
;mfile = '/global/homes/d/dpietrob/myscratch/dx9/maps/ns0128/dust_0.06mK_mask_ns0128.fits'
;mfile = '/global/homes/d/dpietrob/myscratch/dx9/maps/ns0128/dust_1.5mK_mask_ns0128.fits'

read_fits_map, mfile, mask
mollview, mask, px=600, win=0, chars=1.5, tit='!6Mask used'

;dir = '/global/homes/d/dpietrob/myscratch/dx9/maps/ns0128/'
;
;mapfile = dir+['dx9_Delta_Imap_030_ns128_uK.fits', $
;               'dx9_Delta_Imap_044_ns128_uK.fits', $
;               'dx9_Delta_Imap_070_ns128_uK.fits', $
;               'dx9_Imap_100_ns128_uK.fits', $
;              'dx9_Imap_143_ns128_uK.fits', $
;              'dx9_Imap_217_ns128_uK.fits', $
;              'dx9_Imap_353_ns128_uK.fits']

dir = '/global/homes/d/dpietrob/myscratch/'+tag+'/maps/ns0128/'

mapfile = dir+tag+['_Imap_100_ns128_uK.fits', $
              '_Imap_143_ns128_uK.fits', $
              '_Imap_217_ns128_uK.fits', $
              '_Imap_353_ns128_uK.fits']

for ifreq = 0, nfreq-1 do begin
;    read_fits_map,'dx4_Imap'+sfreq[ifreq]+'GHz_ns128_uK.fits', temp
;    read_fits_map,'noMono_dx4_Imap'+sfreq[ifreq]+'GHz_ns128_uK.fits', temp
;    read_fits_map,'../ns1024/dx4_ugd_Imap_'+sfreq[ifreq]+'GHz_ns1024_uK.fits', temp
;##    if (do_ffp4) then read_fits_map,'ffp4_01a_'+sfreq[ifreq]+'_ns128_uK.fits', temp
    read_fits_map, mapfile[ifreq], temp
    maps[*,ifreq] = temp[*,0]
;    mollview, temp, min=-250, max=250.
    bpix = where(abs(temp) gt 200. )
;    mask[bpix,0] = 0.
endfor

;for ifreq = 0, nfreq-1 do maps[*,ifreq] = maps[*,ifreq] * cmb2sz[ifreq]

gpix = where(mask[*,0] gt 0.)
ngpix = n_elements(gpix)
bpix = where(mask[*,0] eq 0.)
nbpix = n_elements(bpix)


gR = fltarr(nfreq, nfreq)
bR = fltarr(nfreq, nfreq)

ave = fltarr(nfreq,2)

for ifreq = 0, nfreq-1 do begin
    gavei = mean(maps[gpix, ifreq])
    ave[ifreq,0] = gavei
    if (bpix[0] ne -1) then begin
        bavei = mean(maps[bpix, ifreq])
        ave[ifreq,1] = bavei
    endif
    for jfreq = ifreq, nfreq-1 do begin
        avej = mean(maps[gpix,jfreq])
        gR[ifreq, jfreq] = total( (maps[gpix,ifreq]-gavei) * (maps[gpix,jfreq]-avej) ) / ngpix
;        gR[ifreq, jfreq] = total( (maps[gpix,ifreq]) * (maps[gpix,jfreq]) ) / ngpix
        gR[jfreq, ifreq] = gR[ifreq, jfreq]

        if (bpix[0] ne -1) then begin        
            avej = mean(maps[bpix,jfreq])
            bR[ifreq, jfreq] = total( (maps[bpix,ifreq]-bavei) * (maps[bpix,jfreq]-avej) ) / nbpix
;        bR[ifreq, jfreq] = total( (maps[bpix,ifreq]) * (maps[bpix,jfreq]) ) / nbpix
            bR[jfreq, ifreq] = bR[ifreq, jfreq]
        endif
    endfor
endfor

gRm1 = invert(gR, /double, status)
print, status

bRm1 = invert(bR, /double, status)
print, status

a = findgen(nfreq) * 0. + 1.

gw = fltarr(nfreq)
gw = total(gRm1,2) / total(gRm1)
bw = total(bRm1,2) / total(bRm1)

print, ' - gw: ', gw
print, ' - bw: ', bw

hfi_ilc = fltarr(npix,3)

for ifreq = 0, nfreq-1 do begin
;##	ilc[gpix] = ilc[gpix] + maps[gpix,ifreq] * gw[ifreq] 
;##        ilc[bpix] = ilc[bpix] + maps[bpix,ifreq] * bw[ifreq]
	hfi_ilc[*,0] = hfi_ilc[*,0] + maps[*,ifreq] * gw[ifreq] 
        hfi_ilc[*,1] = hfi_ilc[*,1] + maps[*,ifreq] * bw[ifreq]
endfor

hfi_ilc[gpix,2] = hfi_ilc[gpix,0]
hfi_ilc[bpix,2] = hfi_ilc[bpix,1]

;write_fits_map, 'ffp4_01a_ILCout.fits', ilc[*,0], /ring, units='!7l!6K'
;write_fits_map, 'ffp4_01a_ILCin.fits', ilc[*,1], /ring, units='!7l!6K'


mollview, hfi_ilc[*,0], min=-300, max=300, chars=1.5, win=1, tit='!6HFI_ILC: wieghts outside mask';, no_monopole=true, gal_cut=40
mollview, hfi_ilc[*,1], min=-300, max=300, chars=1.5, win=2, tit='!6ILC: wieghts inside mask';, no_monopole=true, gal_cut=40
mollview, hfi_ilc[*,2], min=-300, max=300, chars=1.5, win=3, tit='!6ILC: combined';, no_monopole=true, gal_cut=40
mollview, hfi_ilc[*,0]-hfi_ilc[*,1], min=-30, max=30, chars=1.5, win=4, tit='!6HFI_ILC Difference';, no_monopole=true, gal_cut=40

if (do_png) then begin
    restore, 'chains/pix_01a_v2.res.sav'
    mollview, cmb, min=-300, max=300, chars=1.5, win=4, tit='!6Commander', no_monopole=true, gal_cut=40
    mollview, cmb-hfi_ilc[*,0], min=-30, max=30, chars=1.5, win=5, tit='!6Commander-HFI_ILC!dout!n', no_monopole=true, gal_cut=40
    mollview, cmb-hfi_ilc[*,1], min=-30, max=30, chars=1.5, win=6, tit='!6Commander-HFI_ILC!din!n', no_monopole=true, gal_cut=40
    
    read_fits_map, 'ffp4_scalar_cmb_ns128_60arcmin_uK.fits', inp
    mollview, cmb-inp, min=-30, max=30, chars=1.5, win=7, tit='!6Commander-Input', no_monopole=true, gal_cut=40
    mollview, hfi_ilc[*,0]-inp, min=-30, max=30, chars=1.5, win=8, tit='!6HFI_ILC!dout!n-Input', no_monopole=true, gal_cut=40
    mollview, hfi_ilc[*,1]-inp, min=-30, max=30, chars=1.5, win=9, tit='!6HFI_ILC!din!n-Input', no_monopole=true, gal_cut=40

    mollview, hfi_ilc[*,0], min=-300, max=300, chars=1.5, win=-1, tit='!6HFI_ILC: wieghts outside mask', no_monopole=true, gal_cut=40, png='ffp4_01a_HFI_ILCout.png'
    mollview, hfi_ilc[*,1], min=-300, max=300, chars=1.5, win=-2, tit='!6HFI_ILC: wieghts inside mask', no_monopole=true, gal_cut=40, png='ffp4_01a_HFI_ILCin.png'
    mollview, hfi_ilc[*,0]-hfi_ilc[*,1], min=-30, max=30, chars=1.5, win=-3, tit='!6HFI_ILC Difference', no_monopole=true, gal_cut=40, png='ffp4_01a_HFI_ILCout-in.png'

    mollview, cmb, min=-300, max=300, chars=1.5, win=-4, tit='!6Commander', no_monopole=true, gal_cut=40, png='ffp4_01a_CMD.png'
    mollview, cmb-hfi_ilc[*,0], min=-30, max=30, chars=1.5, win=-5, tit='!6Commander-HFI_ILC!dout!n', no_monopole=true, gal_cut=40, png='ffp4_01a_CMD-HFI_ILCout.png'
    mollview, cmb-hfi_ilc[*,1], min=-30, max=30, chars=1.5, win=-6, tit='!6Commander-HFI_ILC!din!n', no_monopole=true, gal_cut=40, png='ffp4_01a_CMD-HFI_ILCin.png'

    mollview, cmb-inp, min=-30, max=30, chars=1.5, win=-7, tit='!6Commander-Input', no_monopole=true, gal_cut=40, png='ffp4_01a_CMD-INP.png'
    mollview, hfi_ilc[*,0]-inp, min=-30, max=30, chars=1.5, win=-8, tit='!6HFI_ILC!dout!n-Input', no_monopole=true, gal_cut=40, png='ffp4_01a_HFI_ILCout-INP.png'
    mollview, hfi_ilc[*,1]-inp, min=-30, max=30, chars=1.5, win=-9, tit='!6HFI_ILC!din!n-Input', no_monopole=true, gal_cut=40, png='ffp4_01a_HFI_ILCin-INP.png'
endif

print, ' --- End of Program ---'

return, hfi_ilc

;stop

; matrix multiplication failures
gw = grm1##a 
n = reform(a##(grm1##a))
gw[*] = gw[*] / n[0]

print, gw

bw = brm1##a
n = reform(a##(brm1##a))
bw[*] = bw[*] / n[0]


stop


end
