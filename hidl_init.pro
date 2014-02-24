   spawn, 'head ~dpietrob/myscratch/Tools/src/pro/hidl_init.pro'

   !path=!path+':/global/scratch/sd/dpietrob/Tools/src/pro/'
   !path=!path+':/global/scratch/sd/dpietrob/Tools/src/myastrolib/pro/'
   !path=!path+':/global/scratch/sd/dpietrob/Tools/src/pro/paperplots-master/IDL/IDL_scripts/'
;.r /global/scratch/sd/dpietrob/Tools/src/myastrolib/pro/legend.pro
;.r HFI_plot.pro

   True = 1b
   False = 0b

   set_init = True

   if set_init then begin
       mollview, randomn(-1,12), px=500, win=0
       loadct, 39
       !p.color=0
       !p.background=255
       
;## --- run HFI_plot first                                                                                                                        
       FDIR = ''
       CTDIR = '/global/scratch/sd/dpietrob/Tools/src/pro/paperplots-master/IDL/IDL_scripts/'
       CTFILE = 'Planck_CT.tbl' ; The HFI_CT script will create this file for you as needed in the specified CTDIR directory                      
       HDRDIR = CTDIR
       HDRFILE = 'RGB_Planck_hdr.idl'

       if set_init then HFI_CT, CTDIR=CTDIR, CTFILE=CTFILE, /LOAD, /HIGHDR, HDRFILE=HDRFILE
;# ---                                                                                                                                            
   endif

