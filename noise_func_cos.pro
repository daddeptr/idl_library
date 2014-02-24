pro noise_func_cos, l, pars, f, pder
; pars : l_piv, indx, const, ampl, l0, dl

   ampl1 = 0
   l_piv = 1
   phi = 2

   ampl2 = 3
   l_piv2 = 4
   indx2 = 5

   ampl3 = 6
   l_piv3 = 7
   indx3 = 8

   const = 9

   f = pars[ampl1] * cos( l/pars[l_piv] + pars[phi] ) + pars[ampl2]*( l/pars[l_piv2] )^pars[indx2] + pars[ampl3]*( l/pars[l_piv3] )^pars[indx3] + pars[const]

   pder1 = cos( l/pars[l_piv] + pars[phi] )
   
   pder2 = pars[ampl1] * sin( l/pars[l_piv] + pars[phi] ) * (l/pars[l_piv]^2)

   pder3 = -pars[ampl1] * sin( l/pars[l_piv] + pars[phi] )

   pder4 = ( l/pars[l_piv2] )^pars[indx2]

   pder5 = pars[ampl2] * pars[indx2] * (l/pars[l_piv2])^(pars[indx2]-1) * l * (-1)/pars[l_piv2]^2

   pder6 = pars[ampl2]*( l/pars[l_piv2] )^pars[indx2]*alog(l/pars[l_piv2])

   pder7 = ( l/pars[l_piv3] )^pars[indx3]

   pder8 = pars[ampl3] * pars[indx3] * (l/pars[l_piv3])^(pars[indx3]-1) * l * (-1)/pars[l_piv3]^2

   pder9 = pars[ampl3]*( l/pars[l_piv3] )^pars[indx3]*alog(l/pars[l_piv3])
   
   pder10 = l*0.+1

; ----------
if (0b) then begin
   f = pars[ampl1] * cos( l/pars[l_piv] + pars[phi] ) * (1. + pars[ampl2]*( l/pars[l_piv2] )^pars[indx2] + pars[ampl3]*( l/pars[l_piv3] )^pars[indx3] ) + pars[const]

   pder1 = cos( l/pars[l_piv] + pars[phi] ) * (1. + pars[ampl2]*( l/pars[l_piv2] )^pars[indx2] + pars[ampl3]*( l/pars[l_piv3] )^pars[indx3] )
   
   pder2 = pars[ampl1] * sin( l/pars[l_piv] + pars[phi] ) * (l/pars[l_piv]^2) * (1. + pars[ampl2]*( l/pars[l_piv2] )^pars[indx2] + pars[ampl3]*( l/pars[l_piv3] )^pars[indx3] )

   pder3 = -pars[ampl1] * sin( l/pars[l_piv] + pars[phi] ) * (1. + pars[ampl2]*( l/pars[l_piv2] )^pars[indx2] + pars[ampl3]*( l/pars[l_piv3] )^pars[indx3] )

   pder4 = pars[ampl1] * cos( l/pars[l_piv] + pars[phi] ) * ( l/pars[l_piv2] )^pars[indx2]

   pder5 = pars[ampl1] * cos( l/pars[l_piv] + pars[phi] ) * pars[ampl2] * pars[indx2] * (l/pars[l_piv2])^(pars[indx2]-1) * l * (-1)/pars[l_piv2]^2

   pder6 = pars[ampl1] * cos( l/pars[l_piv] + pars[phi] ) * pars[ampl2]*( l/pars[l_piv2] )^pars[indx2]*alog(l/pars[l_piv2])

   pder7 = pars[ampl1] * cos( l/pars[l_piv] + pars[phi] ) * ( l/pars[l_piv3] )^pars[indx3]

   pder8 = pars[ampl1] * cos( l/pars[l_piv] + pars[phi] ) * pars[ampl3] * pars[indx3] * (l/pars[l_piv3])^(pars[indx3]-1) * l * (-1)/pars[l_piv3]^2

   pder9 = pars[ampl1] * cos( l/pars[l_piv] + pars[phi] ) * pars[ampl3]*( l/pars[l_piv3] )^pars[indx3]*alog(l/pars[l_piv3])
   
   pder10 = l*0.+1
endif
; ----------------
   if (n_params() ge 4) then $
       pder = [ [pder1],[pder2], [pder3], [pder4], [pder5],[pder6], [pder7], [pder8],[pder9],[pder10] ]

;help, pder
end
