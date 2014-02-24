mollview, randomn(-1,12), win=-1
loadct, 39
!p.color=0
!p.background=255

sfreq = ['030','044','070','100','143','217','353','545','857']
freq = [28.44, 44.12, 70.33, 101.28, 143.85, 222.51, 361.46, 557.51, 864.28]
nfreq = n_elements(freq)-2

cf = conversionfactor(freq, /thermo2antenna)

ns = 128l
npix=12l*ns^2

cut=make_sky_cut(5., ns)
read_fits_map, 'fin_region_disc_mask.fits', cut

gp = where(cut eq 1)
ngp = n_elements(gp)

maps = fltarr(npix, nfreq)

iref = 6

lowf_pivot = freq[0]
dust_pivot = freq[iref]
dust2_pivot = freq[4]

beta = -2.8
emis = 2.6
temp = 16.

emis2 = 1.5
temp2 = 10.

co_spec= [0., 0., 0., 1., 0., .39, .124, 0., 0.]

h = 6.626068 * 10.^(-34)
k = 1.3806503 * 10.^(-23)
c = 2.99792458d8

x = h*1.d9/k/temp
bb  = 1. / ( exp(x*freq)-1. )
bb0 = 1. / ( exp(x*dust_pivot)-1. )

x2 = h*1.d9/k/temp2
bb2  = 1. / ( exp(x2*freq)-1. )
bb02 = 1. / ( exp(x2*dust2_pivot)-1. )

lowf = (freq/lowf_pivot)^beta 
dust = (freq/dust_pivot)^(emis+1)*bb/bb0
dust2 = (freq/dust2_pivot)^(emis2+1)*bb2/bb02

lowf_amp = 1000.
dust_amp = 1000.
dust2_amp = 0.;0.5;1000.
co_amp = 50. 

for ifreq=0, nfreq-1 do begin
    read_fits_map, 'ffp4_01a_'+sfreq[ifreq]+'_ns128_uK.fits', map
    maps[*,ifreq] = map
endfor

window, 10
plot_oo, /nodata, [10,1000], [1.e-2,1.e5], chars=1.5, xtit='!7m', ytit='!6S!d!7m!n!6', ys=1.2, xs=1

sed = fltarr(ngp,nfreq)

model = lowf*lowf_amp + dust*dust_amp + co_spec*co_amp + dust2*dust2_amp

for ipix=0,ngp-1 do begin
    oplot, freq, maps[gp[ipix],*]*cf ;/ model
    sed[ipix,*] = maps[gp[ipix],*]*cf ;/ model
endfor


oplot, freq, model, col=245, thick=1.5
oplot, freq, lowf*lowf_amp, col=210
oplot, freq, dust*dust_amp, col=75
oplot, freq, dust2*dust2_amp, col=105
oplot, freq, co_spec*co_amp, col=45, psym=1

print, total(sed,1)/ngp

STOP

END
