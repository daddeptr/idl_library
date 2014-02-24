 pro plot_cls_input, gsig=gsig, ns=ns, smoo=smoo

    loadct, 39
    freq = ['030', '044', '070', '100', '143', '217', '353']
    nfreq = n_elements(freq)

    ss = strtrim(string(long(gsig)),2)
    s_ns = strtrim(string(ns),2)
    ssmoo = strtrim(string(smoo),2)

;    set_plot, 'ps'
;    device, file='real_data_cls.eps', /col, bits=8
    plot_io, /nodata, [-20, 3*ns], [.01, 100000], xs=2,ys=2, chars=1.5, xtit='!17l', ytit='!17C!dl!n'
    l = findgen(3*ns+1)
    ll = l*(l+1)/2./!pi

;    maskf = 'ctp3_mask_ns16_RING.fits'
;    maskf = 'sky_cut_20_ns16_ring.fits'

    for ifreq=0, nfreq-1 do begin

       ianafast, 'ns'+s_ns+'/noDip_madam_map_'+freq[ifreq]+'GHz_I_smth'+ssmoo+'_uK'+ss+'_ns'+s_ns+'.fits', cls, theta_cut_deg=20., nlmax=3*ns
       oplot, l, cls, col=200+5*ifreq, line=ifreq
    endfor

    npix = 12l*ns^2
    noise = gsig*randomn(-4,npix)
    ianafast, noise, ncls, theta_cut_deg=20., nlmax=3*ns, /ring
    oplot, l, ncls, line=2
    oplot, l, ncls*0.+gsig^2*4.*!pi/npix, line=1
;    oplot, l, ncls*0.+3.^2*4.*!pi/3072l, line=1
    beam = gaussbeam(smoo,3*ns)
    oplot, l, beam^2
    oplot, l, ncls/beam^2, line=3
;    ianafast, noise, ncls, nlmax=40, /ring
;    oplot, l, ncls, line=2

        ianafast, 'internal_sync_templ_030GHz_smth600_ns'+s_ns+'.fits', cls, theta_cut_deg=20., nlmax=3*ns, /ring
    oplot, l, cls, col=70

;    ianafast, '../new_ctp3_CompSep_143GHz_all_regN3_ns16_RING.fits', cls, maskfile=maskf, nlmax=40
;    oplot, l, cls, col=150
;    ianafast, '../new_ctp3_CompSep_143GHz_noise_regN3_ns16_RING.fits', cls, maskfile=maskf, nlmax=40
;    oplot, l, cls, col=150, line=2

        ianafast, 'internal_dust_templ_217GHz_smth600_ns'+s_ns+'.fits', cls, theta_cut_deg=20., nlmax=3*ns, /ring
    oplot, l, cls, col=90
    legend, [freq,'3!7l!17K Noise','Sync 30','Dust 217','ctp3 143'], col=[200+5*lindgen(nfreq), 0, 70,90,150],line=[lindgen(nfreq),2,0,0,0], pos=[-20,160], chars=1. 
;    device, /close
;    set_plot, 'x'

stop
 end
