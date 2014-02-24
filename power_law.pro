function power_law, nu, beta, nu_ref=nu_ref, curv=curv

   if ( not keyword_set(nu_ref) ) then nu_ref = 30.d0

   true  = 1b
   false = 0b

   if ( not keyword_set(curv) ) then pl_sed = (nu / nu_ref)^beta
   if ( keyword_set(curv) ) then pl_sed = (nu / nu_ref)^(beta+curv*alog(nu/nu_ref))

   return, pl_sed

END
