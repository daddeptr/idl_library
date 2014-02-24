;## ====================================================================

!path=!path+':/project/projectdirs/planck/user/dpietrob/ctp3/CompSep/real_data/pro/'
!path=!path+':/project/projectdirs/planck/user/dpietrob/software/myastrolib/pro/'

True = 1b
False = 0b

incomp = ['cmb', 'co','dust', 'ff']
ncomp = n_elements(comp)

sets = ['dx9']

internal_corr = True

dataset = sets[0]

;---------------------------------------------------------------------

if (internal_corr) then begin

    Nside = 2048l
    Npix = 12l*Nside^2
    nComp = 4
    solution = fltarr(Npix, nComp)

    components = ['cmb', 'fg_co_amp', 'fg_dust_amp', 'fg_lowfreq_amp']

    sversion = 'v1.0.128'

    root='/global/scratch/sd/dpietrob/commander-ruler/'+dataset+'/'+sversion+'/dx9_commander-ruler_'+sversion+'_'

    gmask = fltarr(Npix)+1.
    for icomp=0,ncomp-1 do begin
        read_fits_map, root+components[icomp]+'_uK.fits', map
        bp = where( (finite(map, /nan) eq True) or (map eq 0.) or (map eq -1.6375e30) )
        map[bp] = 0.
        in = reorder(map, in='ring', out='nest')
        solution[*,icomp] = in
;##        mollview, in, chars=1.5, win=icomp+10, /hist, /nest
        gmask[bp] = 0.
    endfor
    
    low_Nside = 64l
    low_Npix = nside2npix(low_Nside)

    correlation_maps = fltarr(low_Npix,6)
    
    correlation_maps[*,*] = -1.6375e30
    
    pix2vec_nest, low_Nside, lindgen(low_Npix), vec, vertex

    for i=0, low_Npix-1 do begin
        
        if ( (i/1000)*1000 eq i) then print, i, ' / ', low_Npix, ' --> '+string(float(i)/low_Npix*100,format='(f4.1)')+'%'

        query_polygon, Nside, transpose( reform( vertex[i,*,*]) ), pix_list, inclusive=True, /nested
        
        corr_gp = where(gmask[pix_list] gt 0.)
        
        if (corr_gp[0] ne -1) then begin
; cmb-co
            correlation_maps[i,0] = correlate(solution[pix_list[corr_gp],0], solution[pix_list[corr_gp],1])
; cmb-dust
            correlation_maps[i,1] = correlate(solution[pix_list[corr_gp],0], solution[pix_list[corr_gp],2])
; cmb-synch
            correlation_maps[i,2] = correlate(solution[pix_list[corr_gp],0], solution[pix_list[corr_gp],3])
; co-dust
            correlation_maps[i,3] = correlate(solution[pix_list[corr_gp],1], solution[pix_list[corr_gp],2])
; co-synch
            correlation_maps[i,4] = correlate(solution[pix_list[corr_gp],1], solution[pix_list[corr_gp],3])
; dust_synch
            correlation_maps[i,5] = correlate(solution[pix_list[corr_gp],2], solution[pix_list[corr_gp],3])

        endif    

    endfor

    clab = strarr(6)
    icnt = 0
    for i=0,3 do begin
        for j=i+1,3 do begin
            clab[icnt] = components[i]+'-'+components[j]
            out = reorder(correlation_maps[*,icnt], in='nest', out='ring')
            write_fits_map, dataset+'_cross_commander_correlation_'+clab[icnt]+'_'+strtrim(string(low_Nside),2)+'_'+sversion+'.fits', out, /ring
            mollview, out, tit='Correlation '+clab[icnt], min=-1, max=1, chars=1.5, win=icnt+1
            mollview, out, tit='Correlation '+clab[icnt], min=-1, max=1, chars=1.5, png=dataset+'_cross_cr_correlation_'+clab[icnt]+'_'+strtrim(string(low_Nside),2)+'_'+sversion+'.png', win=-1
            icnt = icnt + 1
            print, clab
        endfor
    endfor


endif

;# ---------------------------------------------------------------------

   STOP


END
