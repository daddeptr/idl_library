pro plot_cosmomc_chains, root, par_tags=par_tags, nbins=nbins, init=init, win=win, values=values, hist=hist, doplot=doplot, verbose=verbose, norm_max=norm_max

   if not keyword_set(win) then win=0
   if keyword_set(init) then begin
       mollview, findgen(12l*4^2), win=win
       loadct, 39
       !p.color=0
       !p.background=255
   endif
   if not keyword_set(ipars) then ipars = [3,4,5,7,8,25]
   if not keyword_set(nbins) then nbins = 30

;   ipars = ipars - 1 
;   print, ipars
   npars = n_elements(par_tags)

   print, ' - reading ', root
   print, ' - reading ', par_tags
   readcol, root+'.paramnames', tpars, format='a'
   ipars = intarr( npars )-1
   for i=0,n_elements(par_tags)-1 do ipars[i] = where(tpars eq par_tags[i])
   
   print, tpars[ipars]
   ipars[where(ipars ge 0)] = ipars[where(ipars ge 0)] + 2
;stop
   print, ' - reading parameters ', ipars+1
;   readcol_fmt = strarr( max(ipars) )
;   plen = n_elements(readcol_fmt)
;   readcol_fmt[*] = 'x'
;   readcol_fmt[ ipars-1 ] = 'f'
;   readcol_fmt[0:plen-2] = readcol_fmt[0:plen-2] + ','
;   readcol_fmt = string(readcol_fmt,format='('+strtrim(string(npars),2)+'a)')

   clen = 0l
   ntpars = 0l
   print, root
   spawn, 'ls '+root+'*.txt', cfiles
;print, cfiles
;stop
   nc = n_elements(cfiles)
   samples = 0l
   for ic=0,nc-1 do begin
       print, cfiles[ic]
       spawn, 'wc '+cfiles[ic], finfo
;   spawn, 'wc -l '+file, clen
       finfo = strsplit(finfo, ' ', /extract)
       if keyword_set(verbose) then print, finfo

       rows = long( finfo[0] )
       if ic eq 0 then begin
           cols = long( finfo[1])/long(finfo[0] )
           if cols gt ntpars then ntpars = cols ;  first call
       endif
       if keyword_set(verbose) then print, rows, cols
       samples = samples + rows
;stop
   endfor
   if samples gt clen then clen = samples
   
   if keyword_set(verbose) then print, clen, ntpars
   if keyword_set(verbose) then stop, root

   p_lab = par_tags

   pars = fltarr(ntpars, clen)
   
   samples = 0
   for ifile=0,nc-1 do begin
       spawn, 'wc '+cfiles[ifile], finfo
;   spawn, 'wc -l '+file, clen
       finfo = strsplit(finfo, ' ', /extract)
;   print, finfo
       rows = long( finfo[0] )
       cols = long( finfo[1])/long(finfo[0] )
       tmp = fltarr( cols, rows )
       openr, 1, cfiles[ifile]
       readf, 1, tmp
       close,1
       pars[0:cols-1,samples:samples+rows-1] = tmp
       samples = samples + rows
   endfor

   tmp = 0.

   if keyword_set(doplot) then begin
       !p.multi = [0,3,ceil(npars/3)]
       if not keyword_set(otit) then begin
           otit = strsplit(root,'/',/extract)
           otit = otit[n_elements(otit)-1]
       endif
       window, win, xsize=720*1.25, ysize=450*1.5, tit=otit
   endif

   hist = fltarr(npars, nbins)
   values = fltarr(npars, nbins)
   for i=0,npars-1 do begin
       if (ipars[i] ge 0) then begin
           print, mean( pars[ ipars[i], * ] )
           print, median( pars[ ipars[i], * ] )
           print, stddev( pars[ ipars[i], * ] )
; ## To account for the molteplicity
           h = fltarr(nbins)
           dx = ( max(pars[ ipars[i], * ])-min(pars[ ipars[i], * ]) ) / (nbins+1)
           xint = findgen(nbins+1)*dx + min(pars[ ipars[i], * ])
           for ibin=0,nbins-1 do begin
               ih = where( (pars[ ipars[i], * ] ge xint[ibin]) and (pars[ ipars[i], * ] lt xint[ibin+1]) )
               h[ibin] = total( pars[0,ih] )
           endfor
;       h = histogram( pars[ ipars[i], * ], locations=vals, nbins=nbins )
           vals = fltarr(nbins)
           vals = ( xint[0:nbins-1] + xint[1:nbins] ) / 2.
           h = h / int_tabulated(vals, h )
           if keyword_set(norm_max) then h = h / max(h)
           hist[i,*] = h
           values[i,*] = vals
           nxtic = 3
           xtic = findgen(nxtic+1)*(max(values[i,*])-min(values[i,*]))/nxtic + min(values[i,*])
           nytic = 4
           ytic = findgen(nytic+1)*(max(hist[i,*])-min(hist[i,*])) / nytic + min(hist[i,*])
           yticnames = string(findgen(nytic+1)/nytic,format='(f3.1)')
           if keyword_set(doplot) then plot, vals, h, chars=3, xtit=p_lab[ipars[i]], ytit='P', xr=[min(values[i,*]), max(values[i,*])], yr=[min(hist[i,*]),max(hist[i,*])], xtickv=xtic, xticks=nxtic, yticks=nytic, ytickv=ytic, ytickname=yticnames
           if keyword_set(doplot) then oplot, values[i,*], hist[i,*], thick=2, col=245
       endif
;       stop
   endfor

   if keyword_set(doplot) then !p.multi=0

;stop

   if keyword_set(verbose) then print, values, hist

end
