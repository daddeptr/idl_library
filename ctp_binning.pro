function ctp_binning, cls, bin_location=bin_location, bin_root=bin_root, silent=silent

   if not keyword_set( bin_location ) then bin_location = '/global/scratch2/sd/dpietrob/Software/XFaster/data/bins/ctp/'
   if not keyword_set( bin_root ) then bin_root = 'CTP_bin_'

   intype = typename( cls )
   if not keyword_set(silent) then print, intype
   if strupcase(intype) eq "STRING" then fits2cl, incls, cls else incls = cls
   ncl = n_elements(incls[0,*])
   nl = n_elements(incls[*,0])

   

return, bin_cls

end
