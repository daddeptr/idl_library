pro mymoll, file1, file2=file2, minv=minv, maxv=maxv, no_monopole=no_monopole, no_dipole=no_dipole, $
            win=win, hist=hist, maskfile=maskfile, grat=grat, px=px, imap=imap, tit=tit, gal_cut=gal_cut, $
            dosmooth=dosmooth, chs=chs, outmap=outmap, units=units, gnom=gnom, pos=pos, nested=nested, ash10=ash10

   ctn = 13
   if keyword_set(ash10) then begin
       FDIR = ''
       CTDIR = '/global/scratch2/sd/dpietrob/Software/src/pro/paperplots-master/IDL/IDL_scripts/'
       CTFILE = 'Planck_CT.tbl' ; The HFI_CT script will create this file for you as needed in the specified CTDIR directory
       HDRDIR = CTDIR
       HDRFILE = 'RGB_Planck_hdr.idl'
       
       planck_ct, CTDIR=CTDIR, CTFILE=CTFILE, /LOAD, /HIGHDR, HDRFILE=HDRFILE
       
       MINV = -3
       MAXV =  7
       
       ctn = 75
   endif

   False = 0b
   True = 1b

   filetype = typename(file1)
;   filetype2 = 'FLOAT'
   if keyword_set(file2) then filetype2 = typename(file2)
   print, filetype

;##   if not keyword_set(filetype) then filetype='FITS'
   if not keyword_set(units) then units='!7l!6K'
   if not keyword_set(chs) then chs=1.
   if keyword_set(minv) then minv=float(minv)
   if keyword_set(maxv) then maxv=float(maxv)
;   if not keyword_set(file2) then file2=''
   if not keyword_set(no_monopole) then no_monopole=False
   if not keyword_set(no_dipole) then no_dipole=False
   if not keyword_set(win) then win=0
   if not keyword_set(hist) then hist=0
   if not keyword_set(grat) then grat=[20,20]
   if not keyword_set(px) then px=650
   if not keyword_set(imap) then imap=0
   if not keyword_set(tit) then tit='Map difference'
   if not keyword_set(gal_cut) then gal_cut=0.

   if strupcase(filetype) eq 'STRING' then begin

       print, file1
       read_fits_map, file1, map1, order=dpord
       if ( (dpord eq 'NESTED') or (dpord eq 'nest') ) then map1=reorder(map1, in='NESTED', out='RING')

       if keyword_set(file2) then begin
           if (strupcase(filetype2) eq "STRING") then begin
               print, file2
               read_fits_map, file2, map2, order=dpord
               if ( (dpord eq 'NESTED') or (dpord eq 'nest') ) then map2=reorder(map2, in='NESTED', out='RING')
           endif
       endif else begin
           map2 = map1*0.
       endelse
       
;##       if imap eq 3 then tmp = (sqrt( (map1[*,1]-map2[*,1])^2+(map1[*,2]-map2[*,2])^2)) else tmp = map1[*,imap]-map2[*,imap]
       if imap eq 3 then tmp = sqrt(map1[*,1]^2+map1[*,2]^2) - sqrt(map2[*,1]^2+map2[*,2]^2) else tmp = map1[*,imap]-map2[*,imap]
       if keyword_set(dosmooth) then begin
           print, ' - smoothing...'
           ismoothing, tmp, tmp, fwhm_arcmin=dosmooth, /ring, /silent
       endif
;       if not keyword_set(minv) then minv=min(tmp)
;       if not keyword_set(maxv) then maxv=max(tmp)
;       print, ' - min/max-v = ', minv, maxv
       if ( keyword_set(maskfile) ) then begin
           read_fits_map, maskfile, mask, nside=dpns, order=dpord
           if ( (dpord eq 'NESTED') or (dpord eq 'nest') ) then mask=reorder(mask, in='NESTED', out='RING')           
;           if n_elements(mask[0,*] eq 1) then mollview, tmp*mask, chars=chs, min=min, max=max, tit=tit, no_monopole=no_monopole, no_dipole=no_dipole, win=win, hist=hist, grat=grat, px=px, gal_cut=gal_cut else $
           if (no_monopole or no_dipole) then begin
              if no_dipole then remove_dipole, tmp, mask, nside=dpns, ordering='RING' else $
                 remove_dipole, tmp, mask, nside=dpns, ordering='RING', /onlymonopole
           endif
           mp = where( mask lt 0.5 )
           if keyword_set(ash10) then begin
               print, ' - scaling using asinh10...'
               tmp = asinh10( tmp[*,0] )
           endif
           tmp[mp] = -1.6375e30
           if (not keyword_set(gnom)) then mollview, tmp, chars=chs, min=minv, max=maxv, tit=tit, no_monopole=no_monopole, no_dipole=no_dipole, win=win, hist=hist, grat=grat, px=px, gal_cut=gal_cut, colt=ctn, units=units else gnomview, tmp, chars=chs, min=minv, max=maxv, tit=tit, win=win, hist=hist, grat=grat, px=px*2/3, rot=pos
       endif else begin
           if keyword_set(ash10) then begin
               print, ' - scaling using asinh10...'
               tmp = asinh10( tmp[*,0] )
           endif
           if (not keyword_set(gnom)) then mollview, tmp, chars=chs, min=minv, max=maxv, tit=tit, no_monopole=no_monopole, no_dipole=no_dipole, win=win, hist=hist, grat=grat, px=px, gal_cut=gal_cut, colt=ctn, units=units else gnomview, tmp, chars=chs, min=minv, max=maxv, tit=tit, win=win, hist=hist, grat=grat, px=px*2/3, rot=pos
       endelse
       if keyword_set(outmap) then outmap = tmp
           
   endif else begin
;## ------ array case
       if keyword_set(nested) then file1 = reorder(file1, in='nested', out='ring')
       if keyword_set(file2) then begin
           if keyword_set(nested) then file2 = reorder(file2, in='nested', out='ring')
;           if imap eq 3 then tmp = (sqrt( (map1[*,1]-map2[*,1])^2+(map1[*,2]-map2[*,2])^2)) else tmp = (file1[*,imap]-file2[*,imap])
           if imap eq 3 then tmp = sqrt(file1[*,1]^2+file1[*,2]^2) - sqrt(file2[*,1]^2+file2[*,2]^2) else tmp = (file1[*,imap]-file2[*,imap])
       endif else begin
           if imap eq 3 then tmp = sqrt(file1[*,1]^2+file1[*,2]^2) else tmp = file1[*,imap]
       endelse

       if keyword_set(dosmooth) then begin
           print, ' - smoothing...'
           ismoothing, tmp, tmp, fwhm_arcmin=dosmooth, /ring, /silent
       endif
;       if not keyword_set(minv) then minv=min(tmp)
;       if not keyword_set(maxv) then maxv=max(tmp)
       if ( keyword_set(maskfile) ) then begin
           read_fits_map, maskfile, mask, nside=dpns, order=dpord
           if ( (dpord eq 'NESTED') or (dpord eq 'nest') ) then mask=reorder(mask, in='NESTED', out='RING')           
;           if n_elements(mask[0,*] eq 1) then mollview, tmp*mask, chars=chs, min=min, max=max, tit=tit, no_monopole=no_monopole, no_dipole=no_dipole, win=win, hist=hist, grat=grat, px=px, gal_cut=gal_cut else $
           mp = where( mask lt 0.5 )
           if keyword_set(ash10) then tmp = asinh10(tmp)
           tmp[mp] = -1.6375e30
       endif
       if keyword_set(ash10) then tmp = asinh10(tmp)
       if not keyword_set(gnom) then mollview, tmp, chars=chs, min=minv, max=maxv, tit=tit, no_monopole=no_monopole, no_dipole=no_dipole, win=win, hist=hist, grat=grat, px=px, gal_cut=gal_cut else $
          gnomview, tmp, chars=chs, min=minv, max=maxv, tit=tit, win=win, hist=hist, grat=grat, px=px*2/3, rot=pos
;         mollview, tmp, chars=chs, min=min, max=max, tit=tit, no_monopole=no_monopole, no_dipole=no_dipole, win=win, hist=hist, grat=grat, px=px, gal_cut=gal_cut
       if keyword_set(outmap) then outmap = tmp
   endelse

return

end
