function get_effective_frequency, f, bp, beta

   if not keyword_set(beta) then beta = 1.d0

   f = double(f)
   bp = double(bp)
   beta = double(beta)

;   indx = lindgen(n_elements(f)-1)
;   cfreq = total( (f[indx+1]-f[indx])*(f[indx+1]^beta*bp[indx+1]+f[indx]^beta*bp[indx])/2 ) / $
;      total( (f[indx+1]-f[indx])*(bp[indx+1]+bp[indx])/2 )
   cfreq = int_tabulated(f, bp*f^beta, /double)

   cfreq = cfreq^(1./beta)
    print, '- Effective frequency: ', cfreq

return, cfreq

end


