function gray_body, nu, emis, temp, nu_ref=nu_ref

   if ( not keyword_set(nu_ref) ) then nu_ref = 353.d9

   true  = 1b
   false = 0b

   k_B     = 1.3806503d-23
   h       = 6.626068d-34
   c       = 2.99792458d8

   x = h / (k_B*temp)

   dust_sed = ( exp(x*nu_ref)-1.d0 ) / ( exp(x*nu)-1.d0 ) * (nu / nu_ref)^(emis+1.d0)

;## (exp(x*nu_ref)-1.d0) / (exp(x*nu)-1.d0) * (nu / nu_ref)**(alpha+1.d0)

   return, dust_sed

END
