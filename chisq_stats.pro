function chisq_pdf, dof, x
   k = dof
   pdf = 1./( 2.^(k/2)*gamma(k/2) ) * x^(k/2-1) * exp(-x/2)
return, pdf
end

function cumulative_pdf, dataset, _EXTRA=extra
   ndata = n_elements(dataset)
   dataset = reform(dataset)
   is = sort( dataset )
   sortdata = dataset[ is ]
   pdf = chisq_pdf(7,sortdata)
   cumulative_pdf = fltarr(ndata)
   intg = cumulative_pdf * 0.
   for i=1l,ndata-1 do intg[i] = (sortdata[i]-sortdata[i-1]) * (pdf[i]+pdf[i-1]) / 2.
   norm= total( intg )
   for i=0l,ndata-1 do cumulative_pdf[i] = total( intg[0:i] ) / norm
return, cumulative_pdf
end


pro chisq_stats, root=root, fst_ch=fst_ch, tag=tag, do_read=do_read, nband=nband, nsample=nsample, fst_samp=fst_samp, nchain=nchain, nside=nside

   !path=!path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'
;   print, !path

   false = 0b
   true = 1b

   close, /all
   loadct, 39
   set_plot, 'x'
   mollview, randomn(-1,12), window=-1
   !p.multi = 0
   !p.background=255
   !p.color=0


    if not (keyword_set(root)) then begin
        print, 'root not defined'
        stop
    endif
    if not (keyword_set(tag)) then begin
        print, 'tag not defined'
        stop
    endif        

   if not (keyword_set(nside)) then nside = 128l
   if not (keyword_set(nchain)) then nchain = 1
   if not (keyword_set(fst_ch)) then fst_ch = 1
   if not (keyword_set(do_read)) then do_read = 'read'
   if not (keyword_set(nband)) then nband = 7
   if not (keyword_set(fst_samp)) then fst_samp = 1l

   print, ' ncs    : ', nside
   print, ' nchain : ', nchain
   print, ' fst_ch : ', fst_ch
   print, ' do_read: ', do_read

   npix = 12l*nside^2
   s_nside = strtrim(string(nside),2)
   freq        = [30, 44, 70, 100, 143, 217, 353, 33, 23, 41, 61, 94]
   freq = string(freq, format='(i3.3)')

   nfreq       = n_elements(freq)

   if (do_read eq 'read') then begin
       nsampch = lonarr(nchain)
       if (nsample lt 0l) then begin
           nsample = 0l
           iich = 0l
           for ich=fst_ch,fst_ch+nchain-1 do begin
               sch = 'c'+string(ich,format='(1i4.4)')
               spawn, 'ls -ltrh '+root+'chisq_'+sch+'_k* | wc -l', nsamp
               nsampch[iich] = long(nsamp)
               nsample = nsample + nsampch[iich]
               iich = iich + 1l
           endfor
       endif else if (nsample gt 0l) then begin 
           nsampch[*] = nsample / nchain
       endif
   
       print, ' nsample:      ', nsample
       print, ' nsampch:      ', nsampch

       if (min(nsampch) lt fst_samp) then begin
          print, " Not enough samples..."
          stop
       endif

       print, nsample
       nsample = nsample - nchain*(fst_samp-1)
       print, fst_samp, nchain
       print, nsample

       chi2 = dblarr(nsample)

       glb_chi2 = fltarr(npix)
       chisq_dist = fltarr(npix,nsample)

       iich = -1l
       tot_isamp = -1l

       warn = 1b

       for ich=fst_ch, fst_ch+nchain-1 do begin
          sch = 'c'+string(ich,format='(1i4.4)')
          iich = iich + 1

          for isam=fst_samp, nsampch[iich] do begin
             tot_isamp = tot_isamp+1
             ts = tot_isamp / 100
             if ( (ts*100) eq tot_isamp ) then begin
                print, tot_isamp
                save, filename=tag+'.chisq_stats.sav', nsample, nsampch, nside, tot_isamp, chisq_dist, glb_chi2, chi2
             endif

         ssam = 'k'+string(isam,format='(1i5.5)')

         read_fits_map, root+'chisq_'+sch+'_'+ssam+'.fits', chisq

         glb_chi2 = glb_chi2 + chisq

         chisq_dist[*,tot_isamp] = chisq
         ngpix = n_elements( where(chisq gt 0.) )
         chi2[tot_isamp] = total(chisq) / ngpix / nband

     endfor
 endfor

 print, nsample

 glb_chi2 = glb_chi2 / nsample

 save, filename=tag+'.chisq_stats.sav', nsample, nsampch, nside, chisq_dist, glb_chi2, chi2

   endif

   if (do_read ne 'read') then begin
        print, "restoring..."
;;         freq = ['030', '044', '070', '100', '143', '217', '353']
        freq        = ['030', '044', '070', '100', '143', '217', '353', '033', '023', '041', '061', '094']
        nfreq = n_elements(freq)
	restore, filename=tag+'.chisq_stats.sav'
        npix = 12l*nside^2
   endif

   gp = where(glb_chi2 gt 0.)
   ngp = n_elements(gp)

   bp = where(glb_chi2 eq 0.)

print, moment(glb_chi2)
;mollview, glb_chi2, chars=1.5, max=4*nband, tit='!8Mean !7v!u2!n'
;window,1 & plot, chi2, chars=1.5, xtit='!8Sample', ytit='!7v!u2!n'

peak_chi2 = glb_chi2 * 0.

window, 2 & plot, /nodata, [0,50] ,[0,0.35], chars=1.5, xtit='!7v!u2!n', ytit='!8P', ys=1

chimin = min( [min(chisq_dist), 0] )
chimax = min( [max(chisq_dist),60] )
k = 7
nx = 150
x = findgen(nx)/nx * chimax

pdf = chisq_pdf( k, x ) ;; 1./( 2.^(k/2)*gamma(k/2) ) * x^(k/2-1) * exp(-x/2)

cpdf = cumulative_pdf( x )

plot, x, pdf, yr=[0,1.1], ys=1, xr=[0,60], xs=1
oplot, x, pdf, col=245, thick=2
oplot, x, cpdf, col=70, thick=2
;stop

; thres = where( cpdf gt .99 )
thres = where( cpdf gt .95 )

print, ' - 0.99 chisq', x[thres[0]]

ks_dist = fltarr(npix)

s2n = fltarr(npix)
ks_dist_routine = fltarr(npix)

pix_cpdf = fltarr( npix, nx )

;; indx = lindgen(1000)+1
;; ks_ref = fltarr(nx)
;; for ix=0,nx-1 do ks_ref[ix] = sqrt(2.*!pi)/x[ix] * total(exp(-(2.*indx-1)^2*!pi^2)/(8.*x[ix]^2))
;; for ix=0,nx-1 do ks_ref[ix] = 1. - 2.* total( (-1.)^(indx-1)*exp(-(2.*indx*x)^2))
;; 
;; oplot, x, ks_ref, col=100, thick=2
; stop

;; read_fits_map,'ks_distance_4graca.fits', ks_dist

for ipix=0l,ngp-1 do begin
    if ( (ipix/1000)*1000 eq ipix) then print, ipix
   h = histogram(chisq_dist[gp[ipix],*], locations=c2, binsize=bins, min=0, max=60)
   ; for ibin=0,nbin-1 do pix_cumulative[gp[ipix],ibin] = total( h[0:ibin] ) / nsample

   fnx = x * 0.
   for ix=0l, nx-1 do fnx[ix] = float( n_elements(where(chisq_dist[gp[ipix],*] le x[ix]) ) ) / nsample
   pix_cpdf[gp[ipix]] = fnx
   pix_prob = fnx[ thres[0] ]
   s2n[gp[ipix]] = pix_prob
   

;   oplot, x, fnx

;   d = max( abs(fnx - cpdf) )
;   prob_ks, d, nsample, p
;
;   ks_dist[gp[ipix]] = d
;   s2n[gp[ipix]] = p


;   tot = total(h) * bins
if (false) then begin
   if ( (ipix/1000)*1000 eq ipix) then begin
       oplot, c2, pix_cumulative[ipix,*], col=ipix
       tot = total(h) * bins
       tsamp = tsamp + 1
       oplot, c2, h/tot, psym=10, col=ipix ;, thick=2
;       y2 = spl_init(c2,h)
;       interpolf = spl_interp(c2, h, y2, x)
;       oplot, x, pdf, thick=2, col=245
;print, c2
;print, h
;help, h
;print, minmax(c2)
       
       mean_pdf = mean_pdf + h
       mean_pdf2 = mean_pdf2 + h^2
       ab=h
       ab[where(h gt 0)] = 1
       nsamp = nsamp + ab
   endif
endif

   m=max(h,im)
   peak_chi2[ipix] = c2[im]

;stop
endfor

save, filename=tag+'_pix_cpdf.chisq_stats.sav', pix_cpdf


;ks_dist[bp] = -1.6375e30
s2n[bp] = -1.6375e30

;mollview, ks_dist_routine, chars=1.5
mollview, s2n, chars=1.5
stop

mean_pdf = mean_pdf / nsamp
mean_pdf2 = mean_pdf2 / nsamp
var = mean_pdf2-mean_pdf^2
rms = sqrt(mean_pdf2-mean_pdf^2)

; mean_pdf = mean_pdf / tsamp

norm = (total(mean_pdf)*bins)

oplot, x, pdf, thick=2, col=245,line=2                                                                                                        
oplot, c2, mean_pdf / norm, thick=2, col=70
errplot, c2, (mean_pdf-rms)/norm, (mean_pdf+rms)/norm, col=70, thick=2
errplot, c2, (mean_pdf-1./sqrt(tsamp))/norm, (mean_pdf+1./sqrt(tsamp))/norm, col=40, thick=2
;legend,['!8Peak !7v!u2!n', '!8PDF'], col=[0,245], line=[0,0], /top,/right

mollview, peak_chi2, chars=1.5, max=2*nband, tit='!8Peak !7v!u2!n'
print, moment(peak_chi2)

ks_dist[bp] = -1.6375e30
mollview, ks_dist, chars=1.5, tit='!8KS - Distance'

ks_significance = ks_dist*0.
for ipix=0l,ngp-1 do begin
; print, ks_dist[gp[ipix]], ks_significance[gp[ipix]]
dist = ks_dist[gp[ipix]]
prob_ks, dist, 486, prob
ks_significance[gp[ipix]] = prob
;print, ks_dist[gp[ipix]], ks_significance[gp[ipix]]
;print, '  '
;stop
endfor

ks_significance[bp] = -1.6375e30
mollview, ks_significance, tit='!8KS - Significance'

;ksone, ks_dist, chisq_pdf, D

print, " --- End of Code ---"

stop

end
