pro make_component_paramfile, runtag, localdir

True = 1b
False = 0b

close, /all
print, 'restorig '+localdir+runtag+'.par.sav'
restore, localdir+runtag+'.par.sav'

openw, 1, localdir+"components_"+runtag+".tmp"
printf,1,"SAMPLE_INSIDE_MASK                     = .true."
printf,1,"NUM_FG_SIGNAL_COMPONENTS               = "+strtrim(string(ncomp),2)
printf,1,"NUM_FG_GIBBS_CYCLE                     = 1"
printf,1,"OUTPUT_CMB_FREQUENCY_MAPS         = .false."
printf,1,"# ----------------------------------------------------------------------"
if (icmb ne -1) then begin
    printf,1,"# Necessary to solve for CMB"
    printf,1,"# CMB"
    printf,1,"COMP_TYPE"+string(icmb,format='(1i2.2)')+"                    =  'cmb'"
    printf,1,"REFERENCE_FREQUENCY"+string(icmb,format='(1i2.2)')+"          =  100.d0"
    printf,1,"REFERENCE_TEMPLATE"+string(icmb,format="(1i2.2)")+"           =  'none'"
    printf,1,"IMPOSE_FG_ORTHOGONALITY"+string(icmb,format="(1i2.2)")+"      = .false."
    printf,1,"OUTPUT_FREQUENCY_COMPONENT_MAPS"+string(icmb,format="(1i2.2)")+"         = .false."
    printf,1,"# ----------------------------------------------------------------------"
endif
if (idust ne -1) then begin
    printf,1,"## Thermal dust component"
    printf,1,"COMP_TYPE"+string(idust,format="(1i2.2)")+"                              = 'one-component_dust'"
    printf,1,"REFERENCE_FREQUENCY"+string(idust,format="(1i2.2)")+"                    =  "+string(dust_ref,format='(f10.3)')
    printf,1,"REFERENCE_TEMPLATE"+string(idust,format="(1i2.2)")+"                     =  'none'"
    printf,1,"INITIALIZATION_MODE"+string(idust,format="(1i2.2)")+"                    =  'default'";'input_map'"
    printf,1,"INIT_INDEX_MAP"+string(idust,format="(1i2.2)")+"_01                      =  'e_1.8_init.fits'"
    printf,1,"INIT_INDEX_MAP"+string(idust,format="(1i2.2)")+"_02                      =  '/global/homes/d/dpietrob/myscratch/dx9/maps/ns0128/DR2_dust_temp/tdr2_filled_init.fits'"
    printf,1,"DEFAULT_EMISSIVITY"+string(idust,format="(1i2.2)")+"                     =   "+string(beta,format='(f4.2)');+" #1.6"
    printf,1,"DEFAULT_DUST_TEMP"+string(idust,format="(1i2.2)")+"                      =   "+string(T,format='(f5.2)');+" #17.5"
    printf,1,""
    printf,1,"EMISSIVITY_PRIOR_UNIFORM_LOW_T"+string(idust,format="(1i2.2)")+"         =   1.0"
    printf,1,"EMISSIVITY_PRIOR_UNIFORM_HIGH_T"+string(idust,format="(1i2.2)")+"        =   2.5"
    printf,1,"EMISSIVITY_PRIOR_GAUSSIAN_MEAN_T"+string(idust,format="(1i2.2)")+"       =   "+string(beta,format='(f4.2)');+" #1.6"
    printf,1,"EMISSIVITY_PRIOR_GAUSSIAN_STDDEV_T"+string(idust,format="(1i2.2)")+"     =   0.01"
    printf,1,""
    printf,1,"DUST_TEMP_PRIOR_UNIFORM_LOW_T"+string(idust,format="(1i2.2)")+"          =   6."
    printf,1,"DUST_TEMP_PRIOR_UNIFORM_HIGH_T"+string(idust,format="(1i2.2)")+"         =   35."
    printf,1,"DUST_TEMP_PRIOR_GAUSSIAN_MEAN_T"+string(idust,format="(1i2.2)")+"        =   "+string(T,format='(f5.2)');+" #17.5"
    printf,1,"DUST_TEMP_PRIOR_GAUSSIAN_STDDEV_T"+string(idust,format="(1i2.2)")+"      =   0.01"
    printf,1,""
    printf,1,"CONSTANT_INDEX"+string(idust,format="(1i2.2)")+"_01                      = .false."
    printf,1,"CONSTANT_INDEX"+string(idust,format="(1i2.2)")+"_02                      = .false."
    printf,1,"OUTPUT_FREQUENCY_COMPONENT_MAPS"+string(idust,format="(1i2.2)")+"        = .false."
    printf,1,"" 
    printf,1,"SPECTRUM_FILENAME"+string(idust,format="(1i2.2)")+"                      = '"+runtag+"_totDUST_spectrum.txt'"
    printf,1,""
    printf,1,"IMPOSE_FG_ORTHOGONALITY"+string(idust,format="(1i2.2)")+"                =      .false."
    printf,1,"PRECOMPUTE_GRID"+string(idust,format="(1i2.2)")+"                        =      .true."
    printf,1,"APPLY_JEFFREYS_PRIOR"+string(idust,format="(1i2.2)")+"                   =      .true."
    printf,1,""
    printf,1,"# ----------------------------------------------------------------------"
endif
if (ispin ne -1) then begin
    printf,1,"# Spinning dust component"
    printf,1,"COMP_TYPE"+string(ispin,format="(1i2.2)")+"                          =      'curved_power_law'"
    printf,1,"REFERENCE_FREQUENCY"+string(ispin,format="(1i2.2)")+"                =       22.57"
    printf,1,"REFERENCE_TEMPLATE"+string(ispin,format="(1i2.2)")+"                 =      'none'"
    printf,1,"INITIALIZATION_MODE"+string(ispin,format="(1i2.2)")+"                =      'default'"
    printf,1,"INIT_INDEX_MAP"+string(ispin,format="(1i2.2)")+"_01                  =      'b_init.fits'"
    printf,1,"INIT_INDEX_MAP"+string(ispin,format="(1i2.2)")+"_02                  =      'c_init.fits'"
    printf,1,"DEFAULT_SPECTRAL_INDEX"+string(ispin,format="(1i2.2)")+"             =       -3.5"
    printf,1,"DEFAULT_CURVATURE"+string(ispin,format="(1i2.2)")+"                  =       -1.8"
    printf,1,""
    printf,1,"INDEX_PRIOR_UNIFORM_LOW_T"+string(ispin,format="(1i2.2)")+"          =       -5.15"
    printf,1,"INDEX_PRIOR_UNIFORM_HIGH_T"+string(ispin,format="(1i2.2)")+"         =       -3.3"
    printf,1,"INDEX_PRIOR_GAUSSIAN_MEAN_T"+string(ispin,format="(1i2.2)")+"        =       -3.5"
    printf,1,"INDEX_PRIOR_GAUSSIAN_STDDEV_T"+string(ispin,format="(1i2.2)")+"      =        0.1"
    printf,1,""
    printf,1,"CURVATURE_PRIOR_UNIFORM_LOW_T"+string(ispin,format="(1i2.2)")+"      =       -9.0"
    printf,1,"CURVATURE_PRIOR_UNIFORM_HIGH_T"+string(ispin,format="(1i2.2)")+"     =       -0.5"
    printf,1,"CURVATURE_PRIOR_GAUSSIAN_MEAN_T"+string(ispin,format="(1i2.2)")+"    =       -1.8"
    printf,1,"CURVATURE_PRIOR_GAUSSIAN_STDDEV_T"+string(ispin,format="(1i2.2)")+"  =        0."
    printf,1,""
    printf,1,"CONSTANT_INDEX"+string(ispin,format="(1i2.2)")+"_01 = .false."
    printf,1,"CONSTANT_INDEX"+string(ispin,format="(1i2.2)")+"_02 = .true."
    printf,1,"OUTPUT_FREQUENCY_COMPONENT_MAPS"+string(ispin,format="(1i2.2)")+"         = .false."
    printf,1,""
    printf,1,"IMPOSE_FG_ORTHOGONALITY"+string(ispin,format="(1i2.2)")+"            =      .false."
    printf,1,"PRECOMPUTE_GRID"+string(ispin,format="(1i2.2)")+"                    =      .true."
    printf,1,"APPLY_JEFFREYS_PRIOR"+string(ispin,format="(1i2.2)")+"               =      .true."
    printf,1,""
    printf,1,"# ----------------------------------------------------------------------"
endif
if (iff ne -1) then begin
    printf,1,"# Free-free component"
    printf,1,"COMP_TYPE"+string(iff,format="(1i2.2)")+"                        = 'power_law'"
    printf,1,"REFERENCE_FREQUENCY"+string(iff,format="(1i2.2)")+"              = "+string(ff_ref,format='(f5.2)')
    printf,1,"REFERENCE_TEMPLATE"+string(iff,format="(1i2.2)")+"               =  'none'"
    printf,1,"INITIALIZATION_MODE"+string(iff,format="(1i2.2)")+"              =      'default'"
    printf,1,"INIT_INDEX_MAP"+string(iff,format="(1i2.2)")+"_01                =      'beta_init.fits'"
    printf,1,"DEFAULT_SPECTRAL_INDEX"+string(iff,format="(1i2.2)")+"           =       -2.15"
    printf,1,"DEFAULT_CURVATURE"+string(iff,format="(1i2.2)")+"                =        0."
    printf,1,"INDEX_PRIOR_UNIFORM_LOW_T"+string(iff,format="(1i2.2)")+"        =       -2.35"
    printf,1,"INDEX_PRIOR_UNIFORM_HIGH_T"+string(iff,format="(1i2.2)")+"       =       -1.50"
    printf,1,"INDEX_PRIOR_GAUSSIAN_MEAN_T"+string(iff,format="(1i2.2)")+"      =       -2.15"
    printf,1,"INDEX_PRIOR_GAUSSIAN_STDDEV_T"+string(iff,format="(1i2.2)")+"    =        0.01"
    printf,1,""
    printf,1,"CONSTANT_INDEX"+string(iff,format="(1i2.2)")+"_01                = .false."
    printf,1,"OUTPUT_FREQUENCY_COMPONENT_MAPS"+string(iff,format="(1i2.2)")+"  = .false."
    printf,1,""
    printf,1,"IMPOSE_FG_ORTHOGONALITY"+string(iff,format="(1i2.2)")+"          =      .false."
    printf,1,"PRECOMPUTE_GRID"+string(iff,format="(1i2.2)")+"                  =      .true."
    printf,1,"APPLY_JEFFREYS_PRIOR"+string(iff,format="(1i2.2)")+"             =      .true."
    printf,1,""
    printf,1,"# ----------------------------------------------------------------------"
endif
if (isync ne -1) then begin
    printf,1,"# Synchrotron component"
    printf,1,"COMP_TYPE"+string(isync,format="(1i2.2)")+"                              = 'power_law'"
    printf,1,"REFERENCE_FREQUENCY"+string(isync,format="(1i2.2)")+"                    = "+string(sync_ref,format='(f5.2)')
    printf,1,"REFERENCE_TEMPLATE"+string(isync,format="(1i2.2)")+"                     =  'none'"
    printf,1,"INITIALIZATION_MODE"+string(isync,format="(1i2.2)")+"                    =      'default'"
    printf,1,"INIT_INDEX_MAP"+string(isync,format="(1i2.2)")+"_01              =      'beta_init.fits'"
    printf,1,"DEFAULT_SPECTRAL_INDEX"+string(isync,format="(1i2.2)")+"         =       -3.05"
    printf,1,"DEFAULT_CURVATURE"+string(isync,format="(1i2.2)")+"              =        0."
    printf,1,"INDEX_PRIOR_UNIFORM_LOW_T"+string(isync,format="(1i2.2)")+"      =       -5.25"
    printf,1,"INDEX_PRIOR_UNIFORM_HIGH_T"+string(isync,format="(1i2.2)")+"     =       -2.35"
    printf,1,"INDEX_PRIOR_GAUSSIAN_MEAN_T"+string(isync,format="(1i2.2)")+"    =       -3.05"
    printf,1,"INDEX_PRIOR_GAUSSIAN_STDDEV_T"+string(isync,format="(1i2.2)")+"  =        0.3"
    printf,1,""
    printf,1,"CONSTANT_INDEX"+string(isync,format="(1i2.2)")+"_01              =       .false."
    printf,1,"OUTPUT_FREQUENCY_COMPONENT_MAPS"+string(isync,format="(1i2.2)")+" =      .false."
    printf,1,""
    printf,1,"IMPOSE_FG_ORTHOGONALITY"+string(isync,format="(1i2.2)")+"        =      .false."
    printf,1,"PRECOMPUTE_GRID"+string(isync,format="(1i2.2)")+"                =      .true."
    printf,1,"APPLY_JEFFREYS_PRIOR"+string(isync,format="(1i2.2)")+"           =      .true."
    printf,1,""
    printf,1,"# ----------------------------------------------------------------------"
endif
if (isz ne -1) then begin
    printf,1,"# SZ component"
    printf,1,"COMP_TYPE"+string(isz,format="(1i2.2)")+"                      =       'sz'"
    printf,1,"REFERENCE_FREQUENCY"+string(isz,format="(1i2.2)")+"            =       353."
    printf,1,"IMPOSE_FG_ORTHOGONALITY"+string(isz,format="(1i2.2)")+"        =      .true."
    printf,1,"PRECOMPUTE_GRID"+string(isz,format="(1i2.2)")+"                =      .false."
    printf,1,"APPLY_JEFFREYS_PRIOR"+string(isz,format="(1i2.2)")+"           =      .false."
    printf,1,"# ----------------------------------------------------------------------"
endif
if (ico1 ne -1) then begin
    printf,1,"# CO template"
    printf,1,"COMP_TYPE"+string(ico1,format="(1i2.2)")+"               = 'tabulated'"
    printf,1,"REFERENCE_FREQUENCY"+string(ico1,format="(1i2.2)")+"     =  "+string(co100f,format="(f10.3)")
    printf,1,"REFERENCE_TEMPLATE"+string(ico1,format="(1i2.2)")+"      =  'none'"
    printf,1,"IMPOSE_FG_ORTHOGONALITY"+string(ico1,format="(1i2.2)")+" = .false."
    printf,1,"SPECTRUM_FILENAME"+string(ico1,format="(1i2.2)")+"       = "+runtag+"_CO_100.txt"
    printf,1,"OUTPUT_FREQUENCY_COMPONENT_MAPS"+string(ico1,format="(1i2.2)")+"         = .false."
    printf,1,"# ----------------------------------------------------------------------"
endif
if (ico2 ne -1) then begin
    printf,1,"# CO template"
    printf,1,"COMP_TYPE"+string(ico2,format="(1i2.2)")+"               = 'tabulated'"
    printf,1,"REFERENCE_FREQUENCY"+string(ico2,format="(1i2.2)")+"     =  "+string(co217f,format="(f10.3)")
    printf,1,"REFERENCE_TEMPLATE"+string(ico2,format="(1i2.2)")+"      =  'none'"
    printf,1,"IMPOSE_FG_ORTHOGONALITY"+string(ico2,format="(1i2.2)")+" = .false."
    printf,1,"SPECTRUM_FILENAME"+string(ico2,format="(1i2.2)")+"       = "+runtag+"_CO_217.txt"
    printf,1,"OUTPUT_FREQUENCY_COMPONENT_MAPS"+string(ico2,format="(1i2.2)")+"         = .false."
    printf,1,"# ----------------------------------------------------------------------"
endif
if (ico3 ne -1) then begin
    printf,1,"# CO template"
    printf,1,"COMP_TYPE"+string(ico3,format="(1i2.2)")+"               = 'tabulated'"
    printf,1,"REFERENCE_FREQUENCY"+string(ico3,format="(1i2.2)")+"     =  "+string(co353f,format="(f10.3)")
    printf,1,"REFERENCE_TEMPLATE"+string(ico3,format="(1i2.2)")+"      =  'none'"
    printf,1,"IMPOSE_FG_ORTHOGONALITY"+string(ico3,format="(1i2.2)")+" = .false."
    printf,1,"SPECTRUM_FILENAME"+string(ico3,format="(1i2.2)")+"       = "+runtag+"_CO_353.txt"
    printf,1,"OUTPUT_FREQUENCY_COMPONENT_MAPS"+string(ico3,format="(1i2.2)")+"         = .false."
    printf,1,"# ----------------------------------------------------------------------"
endif
if (iame ne -1) then begin
    printf,1,"# Spinning dust component"
    printf,1,"COMP_TYPE"+string(iame,format="(1i2.2)")+"                          =      'AME'"
    printf,1,"REFERENCE_FREQUENCY"+string(iame,format="(1i2.2)")+"                =       20.00"
    printf,1,"REFERENCE_TEMPLATE"+string(iame,format="(1i2.2)")+"                 =      'none'"
    printf,1,"INITIALIZATION_MODE"+string(iame,format="(1i2.2)")+"                =      'default'"
    printf,1,"INIT_INDEX_MAP"+string(iame,format="(1i2.2)")+"_01                  =      'b_init.fits'"
    printf,1,"INIT_INDEX_MAP"+string(iame,format="(1i2.2)")+"_02                  =      'c_init.fits'"
    printf,1,"DEFAULT_NU_P"+string(iame,format="(1i2.2)")+"             =       20."
    printf,1,"DEFAULT_M_60_"+string(iame,format="(1i2.2)")+"                  =       -1.8"
    printf,1,""
    printf,1,"NU_P_PRIOR_UNIFORM_LOW_T"+string(iame,format="(1i2.2)")+"          =        5"
    printf,1,"NU_P_PRIOR_UNIFORM_HIGH_T"+string(iame,format="(1i2.2)")+"         =        60"
    printf,1,"NU_P_PRIOR_GAUSSIAN_MEAN_T"+string(iame,format="(1i2.2)")+"        =        20"
    printf,1,"NU_P_PRIOR_GAUSSIAN_STDDEV_T"+string(iame,format="(1i2.2)")+"      =        1."
    printf,1,""
    printf,1,"M_60_PRIOR_UNIFORM_LOW_T"+string(iame,format="(1i2.2)")+"      =       -9.0"
    printf,1,"M_60_PRIOR_UNIFORM_HIGH_T"+string(iame,format="(1i2.2)")+"     =       -0.5"
    printf,1,"M_60_PRIOR_GAUSSIAN_MEAN_T"+string(iame,format="(1i2.2)")+"    =       -1.8"
    printf,1,"M_60_PRIOR_GAUSSIAN_STDDEV_T"+string(iame,format="(1i2.2)")+"  =        0.1"
    printf,1,""
    printf,1,"CONSTANT_INDEX"+string(iame,format="(1i2.2)")+"_01 = .false."
    printf,1,"CONSTANT_INDEX"+string(iame,format="(1i2.2)")+"_02 = .false."
    printf,1,"OUTPUT_FREQUENCY_COMPONENT_MAPS"+string(iame,format="(1i2.2)")+"         = .false."
    printf,1,""
    printf,1,"SPECTRUM_FILENAME"+string(iame,format="(1i2.2)")+"                      = 'totDUST_spectrum.txt'"
    printf,1,""
    printf,1,"IMPOSE_FG_ORTHOGONALITY"+string(iame,format="(1i2.2)")+"            =      .false."
    printf,1,"PRECOMPUTE_GRID"+string(iame,format="(1i2.2)")+"                    =      .true."
    printf,1,"APPLY_JEFFREYS_PRIOR"+string(iame,format="(1i2.2)")+"               =      .true."
    printf,1,""
    printf,1,"# ----------------------------------------------------------------------"
endif

close, /all

spawn, "rm param_"+runtag+".txt"
spawn, "cat /global/scratch/sd/dpietrob/dx9/maps/ns0128/param_template.txt "+localdir+"channels_"+runtag+".tmp "+localdir+"components_"+runtag+".tmp > "+localdir+"param_"+runtag+".txt"

stop
;
;FIX_FREE_TEMP03_01 = .true.

end


