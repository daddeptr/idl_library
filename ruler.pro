pro ruler, solve_for=solve_for, tag=tag, ichannels=ichannels, dataset=dataset, mapdir=mapdir, maproot=maproot, fpix_frac=fpix_frac, lpix_frac=lpix_frac, mapns=mapns, no_noise=no_noise, rmsdir=rmsdir, filedir=filedir, do_check=do_check, sync_pivot=sync_pivot, dust_pivot=dust_pivot, split_co=split_co, correct=correct

; Ruler: linear solution in pixel space of the component separation
; problem
;
; D.Pietrobon Jul 2011
;
; To be improved and generalized to be flexible on foreground model
; and translated into f90 for speed
;
; chisq = sum_ch (D_ch-M_ch(A_i))^2/N_ch^2
;
; dchisq/dA_i = 0 = \sum_ch 1/N_ch^2 (D_ch-\sum_j A_j*Ind_ch,j)*I_ch,i
;
; M # A = B
;
; A = (cmb, synch, ff, dust)
;
; B = \sum_ch D_ch/N_ch^2 * Ind_ch ; Ind_ch array of spectral
; responses

; --- Foreground Model definition
 
false = 0b
true  = 1b

ruler_mode = 'ruler'
wiener_mode = 'no_wiener'

if not keyword_set(do_check) then begin
   print, 'do_check missing'
   stop
endif

if (do_check eq 't') then check=true
if (do_check eq 'f') then check=false

if (not keyword_set(solve_for)) and (not keyword_set(tag)) and (not keyword_set(ichannels)) and (not keyword_set(dataset)) and (not keyword_set(mapdir)) and (not keyword_set(maproot)) then begin
    print,'Syntax - ruler, solve_for=["cmb","dust","co"], ichannel=[1,2,3], dataset="f/1h/2h", mapdir=, tag=, maproot=, fpix_frac=, lpix_frac, mapns=, split_co='
    return
 endif

;; if ( not keyword_set(ruler_mode)) then ruler_mode = 'ruler'
;; if ( not keyword_set(wiener_mode)) then wiener_mode = 'no_wiener'
if ( not keyword_set(solve_for) ) then solve_for = ['cmb', 'dust', 'co']
if ( not keyword_set(tag)) then tag = '3b'
if ( not keyword_set(ichannels) ) then ichannels = [2, 3, 4]
if ( not keyword_set(dataset)) then dataset = 'f'
if ( not keyword_set(mapdir)) then mapdir = '/global/scratch/sd/dpietrob/dx7_delivery_analysis/input/mondip_subtracted_map/'
if ( not keyword_set(rmsdir)) then rmsdir = '/global/scratch/sd/dpietrob/dx7_delivery_analysis/input/'
if ( not keyword_set(filedir)) then filedir = './'
if ( not keyword_set(maproot)) then maproot = ''
if ( not keyword_set(rmsroot)) then rmsroot = ''
if ( not keyword_set(fpix_frac)) then fpix_frac = 0.
if ( not keyword_set(lpix_frac)) then lpix_frac = 1.
if ( not keyword_set(mapns)) then mapns = 2048l
if ( not keyword_set(no_noise)) then no_noise='f'
if ( not keyword_set(sync_pivot)) then sync_pivot=30.
if ( not keyword_set(dust_pivot)) then dust_pivot=353.
if ( not keyword_set(split_co)) then split_co='dont_split'
if ( not keyword_set(correct)) then correct='f'


do_ruler = 0b
do_wiener = 0b
if (ruler_mode eq 'ruler') then do_ruler = 1b
if (wiener_mode eq 'wiener') then do_wiener = 1b

;; dataset = 'f'
;; filedir = './' 
;; mapdir = 'mondip_subtracted_map/'
;; mapdir = 'co_removed_map/'
;; mapdir = '5arcmin_convolved_map/'
;; rmsdir = './'

if (do_ruler) then begin

code = 'ruler'
temp_folder = '/global/scratch/sd/dpietrob/rubbish/'

cfreq = [30., 44.,70, 100, 143, 217, 353, 545, 857, 23, 33, 41, 61, 94]
; cCO_spec= [1, 0, .6, .2, 0., 0.]
cCO_spec= [0., 0., 0., 1, 0, 0.35, .12, 0., 0., 0, 0, 0, 0, 0] ; new values
; To solve for line ratio
cCO_spec_1_0= [0., 0., 0., 1., 0, 0., .12, 0., 0., 0, 0, 0, 0, 0] ; new values
cCO_spec_2_1= [0., 0., 0., 0., 0, 1., .35, 0., 0., 0, 0, 0, 0, 0] ; new values

real_beam = [32.65, 27.00, 13.01, 9.94,   7.04,    4.66,    4.41,    4.47,    4.23]
pf=(cfreq/30.)^(-3)
;plot_oo, cfreq, pf
;oplot, cfreq, pf*cor, col=245
;oplot, cfreq, (cfreq/30.)^(-2.15)*cor, col=70
;oplot, cfreq, (cfreq/353.)^(1.5)*cor, col=210

;stop

freq = cfreq[ichannels]
sfreq = string(freq,format='(i3.3)')
nfreq = n_elements(freq)
;; ns = 2048l
npix = 12l*mapns^2
h = 6.626068 * 10.^(-34)
k = 1.3806503 * 10.^(-23)
c = 3.d8

print, fpix_frac, lpix_frac

fpix = long( (npix)*fpix_frac )
lpix = long( (npix)*lpix_frac )-1

print, 'f/l pix: ', fpix, lpix

convfact = conversionfactor(freq, /antenna2thermo)
print, convfact

cor = real_beam*0.
l=findgen(4097)
for ifreq=0,6 do begin
    bl=gaussbeam(real_beam[ifreq],4096)
    cor[ifreq] = total( (2*l+1)/2./!pi*bl )
endfor

cor = cor/cor[0]
print, cor

if (correct eq 't') then convfact = convfact / cor

; Moved to parameters
;; sync_pivot = 70.
;; dust_pivot = 353.

if (split_co eq 'dont_split') then co_spec = cCO_spec[ichannels]
if (split_co eq 'split') then begin
    co_spec1 = cCO_spec_1_0[ichannels]
    co_spec2 = cCO_spec_2_1[ichannels]
endif

; --- from pix_2.res.sav              
if (false) then begin
    read_fits_map,'beta_init_tgh.fits', beta
    read_fits_map,'emissivity_init_tgh.fits', emis
    read_fits_map,'t_init_tgh.fits', td
endif
; --- from cls_hke.res.sav
if (true) then begin
   if (mapns eq 2048) then begin
;      read_fits_map,rmsdir+'cls_hke_peak_beta_init_tgh.fits', beta
;      read_fits_map,rmsdir+'cls_hke_peak_emisivity_init_tgh.fits', emis
;      read_fits_map,rmsdir+'cls_hke_peak_t_init_tgh.fits', td
       read_fits_map,rmsdir+'pix_ns256_mean_emissivity_init_harmonics_filter.fits', emis
       read_fits_map,rmsdir+'pix_ns256_mean_t_init_harmonics_filter.fits', td
       read_fits_map,rmsdir+'pix_ns256_mean_beta_init_harmonics_filter.fits', beta
       beta[*] = -3.
   endif

   if (mapns eq 1024) then begin
      read_fits_map,rmsdir+'cls_hke_peak_beta_init_tgh_ns1024.fits', beta
      print, ' !!! - WARNING - !!!'
      print, ' Setting all indices to mean value of the map'
      beta[*] = mean( beta )
      read_fits_map,rmsdir+'cls_hke_peak_emisivity_init_tgh_ns1024.fits', emis
      emis[*] = mean( emis )
      read_fits_map,rmsdir+'cls_hke_peak_t_init_tgh_ns1024.fits', td
      td[*] = mean( td )
   endif

   if (mapns eq 1024) then begin
      read_fits_map,rmsdir+'cls_hke_peak_beta_init_tgh_ns1024.fits', beta
      print, ' !!! - WARNING - !!!'
      print, ' Setting all indices to mean value of the map'
      beta[*] = mean( beta )
      read_fits_map,rmsdir+'cls_hke_peak_emisivity_init_tgh_ns1024.fits', emis
      emis[*] = mean( emis )
      read_fits_map,rmsdir+'cls_hke_peak_t_init_tgh_ns1024.fits', td
      td[*] = mean( td )
      endif
endif
;---

ncomp = n_elements(solve_for) 
map = fltarr(npix,nfreq)
rms = fltarr(npix,nfreq)
components = fltarr(npix, ncomp)

print, ' *** ***'
print, ' - ruler_mode: ', ruler_mode
print, ' - wiener_mode: ', wiener_mode
print, ' - freq: ', freq
print, ' - ncomp: ', ncomp
print, ' - solve_for: ', solve_for
print, ' - tag: ', tag
print, ' - filedir: ', filedir
print, ' - mapdir: ', mapdir
print, ' - maproot: ', maproot
print, ' - rmsdir: ', rmsdir
print, ' - ichannels: --> ', freq
print, ' - dataset: ', dataset
print, ' - do_check: ', do_check, ' --> ', check
print, ' - split_co: ', split_co
print, ' - correct: ', correct

; --- Noise weighted Data
for ifreq=0, nfreq-1 do begin
   print, ' Loading freq '+sfreq[ifreq]
   read_fits_map, mapdir + maproot + '_'+ sfreq[ifreq] + '_' + dataset + '.fits', m, nside=map_ns, order=map_order
;;    read_fits_map,mapdir+'dx7_Imap'+sfreq[ifreq]+'_'+dataset+'.fits', m, nside=map_ns, order=map_order
;;    read_fits_map, mapdir+'sim_map_'+sfreq[ifreq]+'.fits', m, nside=map_ns, order=map_order
;;    read_fits_map,mapdir+'dx7_'+sfreq[ifreq]+'_d.fits', m, nside=map_ns, order=map_order
   if ( (strtrim(map_order,2) eq 'nested') or (strtrim(map_order,2) eq 'NESTED') ) then begin
       print, ' - reordering map...'
       m = reorder(m, /n2r)
   endif
   map[*,ifreq] = m ;; *1.e6
;   mollview, map[*,ifreq], /hist, tit='Map '+sfreq[ifreq]

   if (map_ns ne mapns) then begin
      print, 'map_ns /= ns'
      print, ' --- STOP ---'
      stop
   endif

   if (no_noise eq 'f') then begin
      read_fits_map,rmsdir+'dx7_Irms'+sfreq[ifreq]+'GHz_ns2048_uK_'+dataset+'.fits', m, nside=map_ns, order=map_order

      if ( (strtrim(map_order,2) eq 'nested') or (strtrim(map_order,2) eq 'NESTED') ) then begin
         print, ' - reordering map...'
         m = reorder(m, /n2r)
      endif
      rms[*,ifreq] = m
;      mollview, rms[*,ifreq], /hist, tit='RMS '+sfreq[ifreq]

      if (map_ns ne mapns) then begin
         print, 'map_ns /= ns'
         print, ' --- STOP ---'
         stop
      endif
   endif else begin
      print, ' Noise NOT taken into account'
      rms[*,ifreq] = 1.
   endelse

endfor

mask = intarr(npix)
mask[*] = 1

;; weights = fltarr(npix, ncomp, ncomp)

; --- System Solution
for ipix=fpix, lpix do begin

   if ( (ipix/1000000l)*1000000 eq ipix) then print, ipix, ' / ', npix

   reduced_data = map[ipix,*] / rms[ipix,*]^2
   spectral_response = dblarr(ncomp,nfreq)

   B = dblarr(ncomp)
   A = dblarr(ncomp, ncomp)

   ; to be improved with routines
   for icomp=0,ncomp-1 do begin
; cmb
       if ( solve_for[icomp] eq 'cmb' ) then   spectral_response[icomp,*] = 1. 
;        print, reduced_data
; thermal dust
       if ( solve_for[icomp] eq 'dust' ) then begin
           x = h*1.d9/k/td[ipix]
           bb  = 1. / ( exp(x*freq[*])-1. )
           bb0 = 1. / ( exp(x*dust_pivot)-1. )
           spectral_response[icomp,*] = convfact[*] * (freq[*]/dust_pivot)^(emis[ipix]+1)*bb/bb0 
       endif
; CO
       if ( solve_for[icomp] eq 'co' ) then spectral_response[icomp,*] = convfact[*] * co_spec[*] 
       if ( (split_co eq 'split') and ( solve_for[icomp] eq 'co1' ) ) then spectral_response[icomp,*] = convfact[*] * co_spec1[*] 
       if ( (split_co eq 'split') and ( solve_for[icomp] eq 'co2' ) ) then spectral_response[icomp,*] = convfact[*] * co_spec2[*] 
; sync
       if ( solve_for[icomp] eq 'sync' ) then spectral_response[icomp,*] = convfact[*] * (freq[*]/sync_pivot)^beta[ipix]
   endfor

   B = transpose(spectral_response) ## reduced_data

   for i=0,ncomp-1 do begin
      for j=i,ncomp-1 do begin
	for ifreq=0,nfreq-1 do begin
            A[i,j] = A[i,j] + spectral_response[i,ifreq]*spectral_response[j,ifreq] / rms[ipix,ifreq]^2
        endfor
            A[j,i] = A[i,j]
      endfor
   endfor
   
   Am1 = la_invert(A,/double, status=info)

if (check) then begin
   print, 'check:'
   print, 'A ', A
   print, ''
   print, 'Am1 ', Am1
   eval = EIGENQL(A, /double, EIGENVECTORS=evec)
   print, ''
   print, 'EVEC ', evec
   print, ''
   print, 'EVAL', eval


   evalm1 = HQR(ELMHES(Am1), /DOUBLE)
   evecm1 = EIGENVEC(Am1, evalm1, /double)
   print, ''
   print, 'EVECm1 ', evecm1
   print, ''
   print, 'EVALm1', evalm1
   stop
endif
;;    weights[ipix,*,*] = float( Am1 )

   X = Am1 ## transpose(B)
   
   if (info /= 0) then mask[ipix] = 0

   if ( ((ipix/1000000)*1000000 eq ipix) or (info /= 0) ) then begin
      if (info /= 0) then print, info
;;      help, X
;;       print, X
   endif

   components[ipix,*] = X
;stop
endfor

;cputim, t2, ts

;print, t2-t1, ' Sec'

for icomp=0,ncomp-1 do write_fits_map,filedir+tag+'_ruler_components_'+solve_for[icomp]+'_'+dataset+'.fits', components[*,icomp], /ring, units='!7l!8K'

;; w_struct = create_struct(HDR,['Weights for each component derived by Ruler'],'CMB_')

;; stop

endif ; --- Ruler

if (do_wiener) then begin


    ianafast,'/global/scratch/sd/dpietrob/dx7_delivery_analysis/output/new/dx7_cmb_3b_s.fits','cls.fits', alm1_out='alms.fits'
    fits2cl, cls, 'cls.fits'
    ianafast,'/global/scratch/sd/dpietrob/dx7_delivery_analysis/output/new/dx7_cmb_3b_1h.fits','xls.fits', map2_in='/global/scratch/sd/dpietrob/dx7_delivery_analysis/output/new/dx7_cmb_3b_2h.fits'

fits2cl, xls, 'xls.fits'

fits2alm, index, alm, 'alms.fits'

wl = smooth(xls/cls,25)

alm_wl = alm*0.
for l=0l, n_elements(cls)-1 do begin
    lm2index, l, lindgen(l+1), lind
    alm_wl[lind,0] = alm[lind,0] * wl[l] 
    alm_wl[lind,1] = alm[lind,1] * wl[l]
endfor

alm2fits, index, alm_wl, 'alm_wl.fits'

isynfast, cls, 'wiener_map_3b.fits', alm_in='alm_wl.fits', nside=2048

mollview, 'wiener_map_3b.fits', min=-300, max=300, chars=1.5

end




end
