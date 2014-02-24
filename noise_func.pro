pro noise_func, l, pars, f, pder
; pars : l_piv, indx, const, ampl, l0, dl

   ampl1 = 0
   l_piv = 1
   indx = 2

   ampl2 = 3
   l_piv2 = 4
   indx2 = 5

   ampl3 = 6
   l0 = 7
   dl = 8

   const = 9

   t1 = pars[ampl1]*( l/pars[l_piv] )^pars[indx]
   pder11 = ( l / pars[l_piv] ) ^ pars[indx]
   pder12 = pars[ampl1] * pars[indx] * ( l/pars[l_piv] )^(pars[indx]-1) * l * (-1)/pars[l_piv]^2
   pder13 = pars[ampl1] * ( l/pars[l_piv] )^pars[indx] * alog(l/pars[l_piv])
   pder1 = [[pder11], [pder12], [pder13]]


   t2 = pars[ampl2]*( l/pars[l_piv2] )^pars[indx2]
   pder2 = [ [( l/pars[l_piv2] )^pars[indx2]], [pars[ampl2] * pars[indx2] * (l/pars[l_piv2])^(pars[indx2]-1) * l * (-1)/pars[l_piv2]^2], [pars[ampl2]*( l/pars[l_piv2] )^pars[indx2]*alog(l/pars[l_piv2])] ]
;help, pder1

   t3 = pars[ampl3] * exp( -(l-pars[l0])^2/pars[dl]^2 )
   pder3 = [ [exp(-(l-pars[l0])^2/pars[dl]^2)], [pars[ampl3]*exp(-(l-pars[l0])^2/pars[dl]^2)*2*(l-pars[l0])/pars[dl]^2], [pars[ampl3]*exp(-(l-pars[l0])^2/pars[dl]^2)*(l-pars[l0])^2*2/pars[dl]^3] ]
;help, pder2

   t4 = pars[const]
   pder4 = l*0.+1.

   f =  t1 + t2 + t3 + t4

   if (n_params() ge 4) then $
       pder = [ [pder1],[pder2], [pder3], [pder4] ]

;help, pder
end
