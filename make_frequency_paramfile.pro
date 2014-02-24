
True = 1b
False = 0b

;## @init_frequencies_CO100_fit_reg
;## @init_frequencies_CO217_fit_reg
@init_frequencies_CO353_fit_reg

;## @init_frequencies_CO217_milca
;## @init_frequencies_CO217_milca_dust

;## @init_frequencies_CO353_milca
;## @init_frequencies_CO353_milca_dust
;## @init_frequencies_CO353_unpol.pro

;## @init_frequencies_haze

;## @init_frequencies_1co
;## @init_frequencies_3co
;## @init_frequencies_check

;## @init_frequencies_pre-run
;## @init_frequencies_run
;## @init_frequencies_test
;## @init_frequencies_CO217_fit
;## @init_frequencies_CO217_fitCO
;## @init_frequencies_CO217_fitDUST

;## @init_frequencies_CO100_fitBOTH
;## @init_frequencies_CO100_fitBOTH_070
;## @init_frequencies_CO100_fitBOTH_070_ff

;## @init_frequencies_CO217_fitBOTH

;## @init_frequencies_CO353
;## @init_frequencies_CO353_fitBOTH
;## @init_frequencies_CO353_fitBOTH_217

use_cmb_cfreq = true ; <--- !!!
                                                                                                        
if (not only_glss) then begin
spawn, 'mkdir '+chains_dir
nprocperband = 4

;## use_bandpasses = '.false.'

; --- LFI
lfi_freq = [ 28.4d0, 44.1d0, 70.3d0 ]
lfi_co =  lfi_freq*0.d0
lfi_tags = ['030','044','070']

lfi_names = strarr(3,2)
lfi_bp = strarr(3)
for i=0,2 do begin
    lfi_names[i,0] = 'dx9_Delta_Imap_'+lfi_tags[i]+'_ns128_uK.fits'
    lfi_names[i,1] = 'dx9_adjrms_'+lfi_tags[i]+'.fits'
    lfi_bp[i] = 'band_pass_'+lfi_tags[i]+'.dat'
endfor

; --- HFI --------------------------------------------------------------
;# hfi_freq = [101.2d0, 143.08d0, 221.9d0, 360.6d0, 559.6d0, 857.d0 ]
hfi_freq = [101.237d0, 142.913d0, 222.141d0, 361.116d0, 557.230d0, 864.356d0 ]
;## hfi_co = dblarr(6)
hfi_co = [1.d0, 0.d0, 0.6d0, 0.08d0, 0.d0, 0.d0]
hfi_tags = ['100','143','217','353','545','857']

hfi_names = strarr(6,2)
hfi_bp = strarr(6)
for i=0,5 do begin
    hfi_names[i,0] = 'dx9_Imap_'+hfi_tags[i]+'_ns128_uK.fits'
    hfi_names[i,1] = 'dx9_adjrms_'+hfi_tags[i]+'.fits'
    hfi_bp[i] = 'band_pass_'+hfi_tags[i]+'.dat'
endfor

; --- 100GHz -----------------------------------------------------------
;## spawn, 'ls /global/scratch/sd/dpietrob/dx9/maps/ns0128/dx9_CO_bp_100*.txt', files
;## spawn, 'ls /global/scratch/sd/dpietrob/dx9/maps/ns0128/dx9_newCO_bp_100_*.txt', files
readcol, '/global/scratch/sd/dpietrob/dx9/maps/ns0128/chains_milca_CO_100.txt.dpc', ff100, d100_co, d100_tags, format='f, f, x, a'
print, ff100, d100_co, d100_tags
ndet100 = n_elements(d100_co)

read_tag = d100_tags
nfiles = n_elements(read_tag)
d100_freq = dblarr(nfiles)

d100_names = strarr(ndet100,2)
d100_bp = strarr(ndet100)
for i=0,nfiles-1 do begin
    d100_freq[i] = ff100[i]
    d100_names[i,0] = 'dx9_Imap_100_ns128_uK_'+d100_tags[i]+'.fits'
    d100_names[i,1] = 'dx9_adjrms_100_ns128_uK_'+d100_tags[i]+'.fits'
    d100_bp[i] = 'dx9_newCO_bp_100_'+read_tag[i]+'.txt'
endfor

; --- 217GHz -----------------------------------------------------------
readcol, '/global/scratch/sd/dpietrob/dx9/maps/ns0128/chains_milca_dust_CO_217.txt.dpc', ff217, d217_co, d217_tags, format='f,f,x,a'
print, ff217, d217_co, d217_tags
ndet217 = n_elements(d217_co)

read_tag = d217_tags
nfiles = n_elements(read_tag)
d217_freq = dblarr(nfiles)

d217_names = strarr(ndet217,2)
d217_bp = strarr(ndet217)
for i=0,nfiles-1 do begin
    d217_freq[i] = ff217[i]
    d217_names[i,0] = 'dx9_Imap_217_ns128_uK_'+d217_tags[i]+'.fits'
    d217_names[i,1] = 'dx9_adjrms_217_ns128_uK_'+d217_tags[i]+'.fits'
    d217_bp[i] = 'dx9_newCO_bp_217_'+read_tag[i]+'.txt'
endfor

; --- 353GHz -----------------------------------------------------------
readcol, '/global/scratch/sd/dpietrob/dx9/maps/ns0128/chains_milca_dust_CO_353.txt.dpc', ff353, d353_co, d353_tags, format='f,f,x,a'
print, ff353, d353_co, d353_tags
ndet353 = n_elements(d353_co)

read_tag = d353_tags
nfiles = n_elements(read_tag)
d353_freq = dblarr(nfiles)

d353_names = strarr(ndet353,2)
d353_bp = strarr(ndet353)
for i=0,nfiles-1 do begin
    d353_freq[i] = ff353[i]
    d353_names[i,0] = 'dx9_Imap_353_ns128_uK_'+d353_tags[i]+'.fits'
    d353_names[i,1] = 'dx9_adjrms_353_ns128_uK_'+d353_tags[i]+'.fits'
    d353_bp[i] = 'dx9_newCO_bp_353_'+read_tag[i]+'.txt'
endfor

; --- WMAP 7-yr --------------------------------------------------------
wmap_freq = [22.570d0, 32.950d0, 40.693d0, 60.942d0, 93.468d0 ]
wmap_co = dblarr(5)
wmap_tags = ['K','Ka','Q','V','W']
wmap_names = strarr(5,2)
wmap_bp = strarr(5)
for i=0,4 do begin
    wmap_names[i,0] = 'wmap7_Imap_raw_'+wmap_tags[i]+'_ns128_uK.fits'
    wmap_names[i,1] = 'wmap7_adjrms_raw_'+wmap_tags[i]+'.fits'
    wmap_bp[i] = 'wmap7_band_pass_'+wmap_tags[i]+'_v2.txt'
endfor

; --- Haslam -----------------------------------------------------------
haslam_freq = [0.408d0]
haslam_co = 0.d0
haslam_tags = ['has']
haslam_names = ['haslam_raw_ns128_uK.fits', 'haslam_raw_rms.fits']
haslam_bp = 'band_pass_has.dat'
;haslam_fixmonopole = 1b
;haslam_fixdipole = 1b
;haslam_sampleCO100 = [0]
;haslam_sampleCO217 = [0]
;haslam_sampleCO353 = [0]

; ----------------------------------------------------------------------
                                                                                                        
co100f = d100_freq[ico1_ref-1]
co217f = d217_freq[ico2_ref-1]
co353f = d353_freq[ico3_ref-1]
print, co100f, co217f, co353f

comp = [icmb, ico1, ico2, ico3, idust, ispin, isync, iff, isz, iame]
ncomp =  n_elements(where(comp ge 0))
print, ' --> Components: ', ncomp

channels = [ilfi, ihfi, i100, i217, i353, ihaslam, iwmap]
nchannels = n_elements(where( channels ge 0))

print, ' Number of channel selected: ', Nchannels

freq = [0.d0]
if ilfi[0] ge 0 then freq = [freq,lfi_freq[ilfi]]
if ihfi[0] ge 0 then freq = [freq,hfi_freq[ihfi]]
if iwmap[0] ge 0 then freq=[freq,wmap_freq[iwmap]]
if ihaslam[0] ge 0 then freq=[freq,haslam_freq[ihaslam]]
if i100[0] ge 0 then freq=[freq,d100_freq[i100]]
if i217[0] ge 0 then freq=[freq,d217_freq[i217]]
if i353[0] ge 0 then freq=[freq,d353_freq[i353]]

freq = freq[1:*]

nfreq = n_elements(freq)
names = strarr(nfreq,2)
bandpasses = strarr(nfreq)

fixmonopole = strarr(nfreq) + '.false.'
fixdipole  = strarr(nfreq) + '.false.'
sampleCO100 = strarr(nfreq) + 'F'
sampleCO217 = strarr(nfreq) + 'F'
sampleCO353 = strarr(nfreq) + 'F'

det_tags = strarr(nfreq)

nco = 3

co_sed = dblarr(nfreq,nco)

iname=0
if ilfi[0] ge 0 then for i=0,n_elements(ilfi)-1 do begin
    names[iname,*] = lfi_names[ilfi[i],*]
    print, ' ---> '+names[iname,*]
    bandpasses[iname] = lfi_bp[ilfi[i]]
    if (lfi_fixmonopole[ilfi[i]]) then fixmonopole[iname] = '.true.'
    if (lfi_fixdipole[ilfi[i]]) then fixdipole[iname] = '.true.'
    if (lfi_sampleCO100[ilfi[i]]) then sampleCO100[iname] = 'T'
    if (lfi_sampleCO217[ilfi[i]]) then sampleCO217[iname] = 'T'
    if (lfi_sampleCO353[ilfi[i]]) then sampleCO353[iname] = 'T'
    co_sed[iname,*] = lfi_co[ilfi[i]]
    det_tags[iname] = lfi_tags[ilfi[i]]
    print, ' ---> '+det_tags[iname]
    iname = iname + 1
endfor
if ihfi[0] ge 0 then for i=0,n_elements(ihfi)-1 do begin
    names[iname,*] = hfi_names[ihfi[i],*]
    print, ' ---> '+names[iname,*]
    bandpasses[iname] = hfi_bp[ihfi[i]]
    if (hfi_fixmonopole[ihfi[i]]) then fixmonopole[iname] = '.true.'
    if (hfi_fixdipole[ihfi[i]]) then fixdipole[iname] = '.true.'
    if (hfi_sampleCO100[ihfi[i]]) then sampleCO100[iname] = 'T'
    if (hfi_sampleCO217[ihfi[i]]) then sampleCO217[iname] = 'T'
    if (hfi_sampleCO353[ihfi[i]]) then sampleCO353[iname] = 'T'
    co_sed[iname,*] = hfi_co[ihfi[i]]
    det_tags[iname] = hfi_tags[ihfi[i]]
    iname = iname + 1
endfor
if iwmap[0] ge 0 then for i=0,n_elements(iwmap)-1 do begin
    names[iname,*]=wmap_names[iwmap[i],*]
    print, ' ---> '+names[iname,*]
    bandpasses[iname] = wmap_bp[iwmap[i]]
    if (wmap_fixmonopole[iwmap[i]]) then fixmonopole[iname] = '.true.'
    if (wmap_fixdipole[iwmap[i]]) then fixdipole[iname] = '.true.'
    if (wmap_sampleCO100[iwmap[i]]) then sampleCO100[iname] = 'T'
    if (wmap_sampleCO217[iwmap[i]]) then sampleCO217[iname] = 'T'
    if (wmap_sampleCO353[iwmap[i]]) then sampleCO353[iname] = 'T'
    co_sed[iname,*] = wmap_co[iwmap[i]]
    det_tags[iname] = wmap_tags[iwmap[i]]
    iname = iname+1
endfor
if ihaslam[0] ge 0 then begin
    names[iname,*] = haslam_names
    print, ' ---> '+names[iname,*]
    bandpasses[iname] = haslam_bp
    if (haslam_fixmonopole) then fixmonopole[iname] = '.true.'
    if (haslam_fixdipole) then fixdipole[iname] = '.true.'
    if (haslam_sampleCO100[ihaslam[0]]) then sampleCO100[iname] = 'T'
    if (haslam_sampleCO217[ihaslam[0]]) then sampleCO217[iname] = 'T'
    if (haslam_sampleCO353[ihaslam[0]]) then sampleCO353[iname] = 'T'
    co_sed[iname,*] = haslam_co
    det_tags[iname] = haslam_tags[ihaslam[0]]
    iname = iname+1
endif
if i100[0] ge 0 then for i=0,n_elements(i100)-1 do begin
    names[iname,*]=d100_names[i100[i],*]
    print, ' ---> '+names[iname,*]
    bandpasses[iname] = d100_bp[i100[i]]
    if (d100_fixmonopole[i100[i]]) then fixmonopole[iname] = '.true.'
    if (d100_fixdipole[i100[i]]) then fixdipole[iname] = '.true.'
    if (d100_sampleCO100[i100[i]]) then sampleCO100[iname] = 'T'
    if (d100_sampleCO217[i100[i]]) then sampleCO217[iname] = 'T'
    if (d100_sampleCO353[i100[i]]) then sampleCO353[iname] = 'T'
    co_sed[iname,0] = d100_co[i100[i]]
    det_tags[iname] = d100_tags[i100[i]]
    iname = iname+1
endfor
if i217[0] ge 0 then for i=0,n_elements(i217)-1 do begin
    names[iname,*]=d217_names[i217[i],*]
    print, ' ---> '+names[iname,*]
    bandpasses[iname] = d217_bp[i217[i]]
    if (d217_fixmonopole[i217[i]]) then fixmonopole[iname] = '.true.'
    if (d217_fixdipole[i217[i]]) then fixdipole[iname] = '.true.'
    if (d217_sampleCO100[i217[i]]) then sampleCO100[iname] = 'T'
    if (d217_sampleCO217[i217[i]]) then sampleCO217[iname] = 'T'
    if (d217_sampleCO353[i217[i]]) then sampleCO353[iname] = 'T'
    co_sed[iname,1] = d217_co[i217[i]]
    det_tags[iname] = d217_tags[i217[i]]
    iname=iname+1
endfor
if i353[0] ge 0 then for i=0,n_elements(i353)-1 do begin
    names[iname,*]=d353_names[i353[i],*]
    print, ' ---> '+names[iname,*]
    bandpasses[iname] = d353_bp[i353[i]]
    if (d353_fixmonopole[i353[i]]) then fixmonopole[iname] = '.true.'
    if (d353_fixdipole[i353[i]]) then fixdipole[iname] = '.true.'
    if (d353_sampleCO100[i353[i]]) then sampleCO100[iname] = 'T'
    if (d353_sampleCO217[i353[i]]) then sampleCO217[iname] = 'T'
    if (d353_sampleCO353[i353[i]]) then sampleCO353[iname] = 'T'
    co_sed[iname,2] = d353_co[i353[i]]
    det_tags[iname] = d353_tags[i353[i]]
    iname = iname + 1
endfor

for i=0,nfreq-1 do print, freq[i];, names[i], bandpasses[i], co_sed[i,*]

ifreq = sort(freq)
freq = freq[ifreq]
names = names[ifreq,*]
bandpasses = bandpasses[ifreq]
fixmonopole = fixmonopole[ifreq]
fixdipole = fixdipole[ifreq]

co_sed = co_sed[ifreq,*]
sampleCO100 = sampleCO100[ifreq]
sampleCO217 = sampleCO217[ifreq]
sampleCO353 = sampleCO353[ifreq]
det_tags = det_tags[ifreq]

save, filename=localdir + runtag+'.par.sav'

if (do_dust_corr) then begin
    dust_corr_file = localdir+runtag+'/fg_tab_c'+string(nchain,format='(i4.4)')+'_no'+string(idust,format='(i2.2)')+'.dat'
    spawn, 'wc -l '+dust_corr_file, lines
    nsample = long(lines[0])-1
    print, nsample
    dust_corr = fltarr(nfreq+1, nsample)
    xxx = ' '
    openr,1,dust_corr_file
    readf,1,xxx
    readf,1,dust_corr
    close,1

    dust_corr = dust_corr[1:*,*]

;##    names[*,0] = names[*,0]+'.mdr'
endif

openw, 1, localdir+'channels_'+runtag+'.tmp'
openw, 5, localdir+runtag+'_totDUST_spectrum.txt'
if (ico1 ne -1) then openw,2, localdir+runtag+'_CO_100.txt'
if (ico2 ne -1) then openw,3, localdir+runtag+'_CO_217.txt'
if (ico3 ne -1) then openw,4, localdir+runtag+'_CO_353.txt'
    printf,1,'# ------ '+runtag+ '---------------------------------------------'
    printf,1,'# Map sets; all data must be given in THERMODYNAMIC temperatures!'
    printf,1,'NUM_REALIZATIONS          = 1'
    printf,1,'NUM_CHAIN_PER_REALIZATION = 1'
    printf,1,'NUM_GROUPS                = 1'
    printf,1,'#NUM_CHAIN                 = 1 #4         # Number of parallel Markov chains to run'
    printf,1,'NUMBAND                   = '+string(nfreq,format='(1i2.2)')
    printf,1,'NUM_PROC_PER_BAND         = '+string(nprocperband,format='(1i2.2)')
    printf,1,'DP_USE_BAND_PASS          = '+use_bandpasses
    printf,1,"CHAIN_DIRECTORY    = '"+chains_dir+"'" ;/";+runtag+"'"
    printf,1,"MASKFILE           = "+maskfile
    printf,1,"CORR_CHISQ_THRESHOLD = 0.1"
    printf,1,"NSIDE_CORR           = 128"
    printf,1,"LMAX_CORR            = 200"
    printf,1,"MASKFILE_CORR        = '"+maskfile+"'"
    printf,1,"OUTPUT_CROSS_CORRELATION_STATS = .false."
    printf,1,"TEMPLATE_AMP_INPUT   = 'none'"
    printf,1,"OUTPUT_MIXING_MATRIX = .false."
    printf,1,"OUTPUT_BAND_CHISQ    = .true."
    printf,1,"NUM_MONO_DIPOLE_BURNIN_STEP = 200"

    for i=0,nfreq-1 do begin
        sfreq = string(i+1,format='(i2.2)')

        print, sfreq+', '+det_tags[i], freq[i], co_sed[i, *]

        printf,1,'# ---------------------------------------------------'
        printf,1,'MAP'+sfreq+'_0001         = '+names[i,0]
        printf,1,'MAP_LOWRES'+sfreq+'       = '+names[i,0]
        printf,1,'NOISE'+sfreq+'            = '+names[i,1]
        printf,1,'NOISE_LOWRES'+sfreq+'     = '+names[i,1]
        printf,1,'BEAM'+sfreq+'             = '+'beam_60arcmin.fits'
        printf,1,'# --- :DP'
        printf,1,'USE_BANDPASS'+sfreq+'     = '+use_bandpasses
        printf,1,'BANDPASS'+sfreq+'         = '+bandpasses[i]
        printf,1,'BAND_PASS_FILE_'+sfreq+'  = '+bandpasses[i]
        printf,1,'A2T'+sfreq+'              = 0.d0'
        printf,1,'FIT_A2T'+sfreq+'              = .false.'
        printf,1,'# ---'
        printf,1,'FREQ_L'+sfreq+'           = '+string(freq[i],format='(1f10.3)')
        printf,1,'FREQ_C'+sfreq+'           = '+string(freq[i],format='(1f10.3)')
        printf,1,'FREQ_U'+sfreq+'           = '+string(freq[i],format='(1f10.3)')
; ---
        if (not do_dust_corr) then printf,5,string(freq[i],format='(1f10.3)'), ' 1.0 F' else printf,5,string(freq[i],format='(1f10.3)'), string(dust_corr[i,nsample-1], format='(1f10.3)'),'  F'
        if (ico1 ne -1) then printf,2,string(freq[i],format='(1f10.3)'), string(co_sed[i,0],format='(1f10.3)'), ' '+sampleCO100[i]+' '+det_tags[i]
        if (ico2 ne -1) then printf,3,string(freq[i],format='(1f10.3)'), string(co_sed[i,1],format='(1f10.3)'), ' '+sampleCO217[i]+' '+det_tags[i]
        if (ico3 ne -1) then printf,4,string(freq[i],format='(1f10.3)'), string(co_sed[i,2],format='(1f10.3)'), ' '+sampleCO353[i]+' '+det_tags[i]
; ---
        printf,1,'FIX_MONOPOLE'+sfreq+'     = '+fixmonopole[i]
        printf,1,'FIX_DIPOLE_X'+sfreq+'     = '+fixdipole[i]
        printf,1,'FIX_DIPOLE_Y'+sfreq+'     = '+fixdipole[i]
        printf,1,'FIX_DIPOLE_Z'+sfreq+'     = '+fixdipole[i]
        printf,1,'FIT_MONO_DIPOLE_BAND'+sfreq+'     = '+sfreq
        if (fixdipole[i] eq '.true.') then printf,1,'FIT_MONO_DIPOLE_IND'+sfreq+'      = 1' else printf,1,'FIT_MONO_DIPOLE_IND'+sfreq+'     = -1'
        if (fixdipole[i] eq '.true.') then printf,1,'SUBTRACT_MONOPOLE_DIPOLE'+sfreq+'      = '+fixdipole[i] else printf,1,'SUBTRACT_MONOPOLE_DIPOLE'+sfreq+'     = '+fixdipole[i]
        if (fixdipole[i] eq '.true.') then printf,1,"MASK_MONOPOLE_DIPOLE"+sfreq+"      = '"+monodipole_maskfile+"'"
        printf,1,''
    endfor
if (false) then begin
    printf,1,'# ----------------------------------------------------------------------'
    printf,1,'# CO template'
    printf,1,"COMP_TYPE'+string(ico1,format='(1i2.2)')+'               = 'tabulated'"
    printf,1,'REFERENCE_FREQUENCY'+string(ico1,format='(1i2.2)')+'     =  '+string(co100f,format='(f10.3)')
    printf,1,"REFERENCE_TEMPLATE'+string(ico1,format='(1i2.2)')+'      =  'none'"
    printf,1,'IMPOSE_FG_ORTHOGONALITY'+string(ico1,format='(1i2.2)')+' = .false.'
    printf,1,'SPECTRUM_FILENAME'+string(ico1,format='(1i2.2)')+'       = '+runtag+'_CO_100.txt'
    printf,1,'# ----------------------------------------------------------------------'
    printf,1,'# CO template'
    printf,1,"COMP_TYPE'+string(ico2,format='(1i2.2)')+'               = 'tabulated'"
    printf,1,'REFERENCE_FREQUENCY'+string(ico2,format='(1i2.2)')+'     =  '+string(co217f,format='(f10.3)')
    printf,1,"REFERENCE_TEMPLATE'+string(ico2,format='(1i2.2)')+'      =  'none'"
    printf,1,'IMPOSE_FG_ORTHOGONALITY'+string(ico2,format='(1i2.2)')+' = .false.'
    printf,1,'SPECTRUM_FILENAME'+string(ico2,format='(1i2.2)')+'       = '+runtag+'_CO_100.txt'
    printf,1,'# ----------------------------------------------------------------------'
    printf,1,'# CO template'
    printf,1,"COMP_TYPE'+string(ico3,format='(1i2.2)')+'               = 'tabulated'"
    printf,1,'REFERENCE_FREQUENCY'+string(ico3,format='(1i2.2)')+'     =  '+string(co353f,format='(f10.3)')
    printf,1,"REFERENCE_TEMPLATE'+string(ico3,format='(1i2.2)')+'      =  'none'"
    printf,1,'IMPOSE_FG_ORTHOGONALITY'+string(ico3,format='(1i2.2)')+' = .false.'
    printf,1,'SPECTRUM_FILENAME'+string(ico3,format='(1i2.2)')+'       = '+runtag+'_CO_100.txt'
    printf,1,'# ----------------------------------------------------------------------'
endif

close, /all

print, ' runtag = ', runtag
print, ' localdir = ', localdir

stop, ' type .c to run make_component_paramfile and make_pbs'

make_component_paramfile, runtag, localdir
make_pbs, runtag, localdir

endif
print, ' GLSS parameter file...'
stop, ' type .c '
;## @make_glss_paramfile
make_glss_paramfile, runtag, localdir, pbsout

;
;FIX_FREE_TEMP03_01 = .true.

end
