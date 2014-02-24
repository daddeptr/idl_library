pro my_ls_mollview, map, init=init, window=window, title=title, px=px, no_monopole=no_monopole, no_dipole=no_dipole, gal_cut=gal_cut, png=png, ps=ps, grat=grat, maskfile=maskfile, rot=rot, basic=basic, minv=minv, maxv=maxv

   if (not keyword_set(window)) then window=1
   if (not keyword_set(title))  then title=' '
   if (not keyword_set(px))     then px=650
   if (not keyword_set(minv))     then minv = -300
   if (not keyword_set(maxv))     then maxv = 300
;##   if (not keyword_set(outtag)) then outtag='map'

; --- Location settings
   FDIR = ''
   CTDIR = '/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/paperplots-master/IDL/IDL_scripts/'
   CTFILE = 'Planck_CT.tbl' ; The HFI_CT script will create this file for you as needed in the specified CTDIR directory
   HDRDIR = CTDIR
   HDRFILE = 'RGB_Planck_hdr.idl'

   if (keyword_set(init)) then HFI_CT, CTDIR=CTDIR, CTFILE=CTFILE, /LOAD, /HIGHDR, HDRFILE=HDRFILE

   MINVAL_ = -1.d3
   MAXVAL_ =  1.d7
; ---
   Ngrat = 181d
   out_ = {COORD:'G',RA:DBLARR(Ngrat),DEC:DBLARR(Ngrat), LINESTYLE:0, PSYM:0, SYMSIZE:0} ; the outline structure accepted by mollview.
   Nout = 2 ; the outline is done as two half-curves.
   out = REPLICATE(out_,Nout)
   RA = DBLARR(Ngrat) ; held constant while DEC changes  ;  The RA and DEC are in healpix/mollview notation.
   DEC = DINDGEN(Ngrat) - 90d ; bottom to top, -90 to +90 deg.
                                                               
   out[0].RA  = RA - 180d ;  The half curve at -180 deg.
   out[0].DEC = DEC
   out[0].LINESTYLE=0

   out[1].RA  = RA + 180d ;  The half curve at +180 deg.
   out[1].DEC = DEC
   out[1].LINESTYLE=0

   GR_88  = [60,45] ;  The graticule spacing for the 88 mm figure.                                                                        
   GR_120 = [60,15] ;  The graticule spacing for the 120 mm figure.                                                                       
   GR_180 = [60,15] ;  The graticule spacing for the 180 mm figure.                                                                       

   W_88 = 8.8d                  ; cm 
   W_120 = 12d                  ; cm 
   W_180 = 18d                  ; cm 

;   NmPref = 'PlanckFig_map_columbi1_IDL_HighDR_'
;   NmPrefHFI = 'PlanckFig_map_HighDR_IDL_'
   NmSuf = 'mm'
;   Nm_88  = NmPref + '88' + NmSuf
;   Nm_120 = NmPref + '120'+ NmSuf
;   Nm_180 = NmPref + '180'+ NmSuf
   
   FigRes = 300d ; in dpi.  The paper figures should be at least 300 dpi, according to the A&A author guide.  For Maps, 600 dpi would be better.

   PX_88  =  8.8d/2.54d*FigRes
   PX_120 = 12.0d/2.54d*FigRes
   PX_180 = 18.0d/2.54d*FigRes
   
   !P.FONT = 0
   SZ    = 88d
   SZstr = '88'
;##   PX  = PX_88
   if (not keyword_set(grat)) then GR  = GR_88
   W   =  W_88
;##   Nm  = Nm_88

   if (keyword_set(maskfile)) then begin
       read_fits_map, maskfile, mask, nside=dpns, order=dpord
       if (dpord eq 'NESTED') then mask = reorder(mask, in=dpord, out='ring')
       ns = npix2nside(n_elements(map[*,0]))
       if (ns ne dpns) then begin
           print, 'Downgrading mask...'
           ud_grade, mask, maskd, nside_out=ns, order_in='ring'
           maskd[where(maskd lt 0.5)] = 0.
           maskd[where(maskd ge 0.5)] = 1.
           mask = maskd
       endif
       if (keyword_set(no_monopole)) then remove_dipole, map, mask, nside=ns, ordering='ring', /onlymonopole
       if (keyword_set(no_dipole)) then remove_dipole, map, mask, nside=ns, ordering='ring'
       no_monopole = 0
       no_dipole= 0
       bp = where(mask eq 0.)
       map[bp] = !healpix.bad_value
   endif

   if (not keyword_set(basic)) then LS_mollview, map, COLT=41, CTDIR=CTDIR, CTFILE=CTFILE, MIN=MINVAL_, MAX=MAXVAL_, CHARSIZE=8d/11d, GRATICULE=GR, $
     GLSIZE=1, HXSIZE=W, FLIP=0, GRMIN=[-179, -89], GRMAX=[179,89], GRLS =1, OUTLINE=out, $
     TITLE=title, PXSIZE=PX, CBLBL='u', /CBLIN, /MODASINH, /CBTICKS, /CBTICKLAB, /CBOUT, UNITS='!7l!6K', $
     window=window, no_monopole=no_monopole, no_dipole=no_dipole, gal_cut=gal_cut, rot=rot

   if (keyword_set(basic)) then LS_mollview, map, COLT=41, CTDIR=CTDIR, CTFILE=CTFILE, MIN=MINV, MAX=MAXV, window=window, no_monopole=no_monopole, no_dipole=no_dipole, gal_cut=gal_cut, rot=rot, UNITS='!7l!6K', grat=grat, /modasinh, title=title
 
;CHARSIZE=8d/11d, GRATICULE=GR, $
;     GLSIZE=1, HXSIZE=W, FLIP=0, GRMIN=[-179, -89], GRMAX=[179,89], GRLS =1, OUTLINE=out, $
;     TITLE=title, PXSIZE=PX, CBLBL='u', /CBLIN, /CBTICKS, /CBTICKLAB, /CBOUT, UNITS='!7l!6K', $

   if (keyword_set(png)) then LS_mollview, map, COLT=41, CTDIR=CTDIR, CTFILE=CTFILE, MIN=MINVAL_, MAX=MAXVAL_, CHARSIZE=8d/11d, GRATICULE=GR, $
     GLSIZE=1, HXSIZE=W, FLIP=0, GRMIN=[-179, -89], GRMAX=[179,89], GRLS =1, OUTLINE=out, $
     TITLE=title, PXSIZE=PX, CBLBL='u', /CBLIN, /MODASINH, /CBTICKS, /CBTICKLAB, /CBOUT, UNITS='!7l!6K', $
     window=window, png=png, no_monopole=no_monopole, no_dipole=no_dipole, gal_cut=gal_cut
 
   if (keyword_set(ps)) then LS_mollview, map, COLT=41, CTDIR=CTDIR, CTFILE=CTFILE, MIN=MINVAL_, MAX=MAXVAL_, CHARSIZE=8d/11d, GRATICULE=GR, $
     GLSIZE=1, HXSIZE=W, FLIP=0, GRMIN=[-179, -89], GRMAX=[179,89], GRLS =1, OUTLINE=out, $
     TITLE=title, PXSIZE=PX, CBLBL='u', /CBLIN, /MODASINH, /CBTICKS, /CBTICKLAB, /CBOUT, UNITS='!7l!6K', $
     window=window, ps=ps, no_monopole=no_monopole, no_dipole=no_dipole, gal_cut=gal_cut
   

end
