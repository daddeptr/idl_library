sfreq = ['030','044','070','100','143','217','353','545','857']
nfreq = n_elements(sfreq)

cf = dblarr(nfreq)
cfi = dblarr(nfreq)

;##beta = [1., 1., 1., 1., 1., 1., 1., 1., 1.]
beta = [-3., -2.5, -2.5, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2]

for ifreq=0,nfreq-1 do begin
    readcol, 'band_pass_'+sfreq[ifreq]+'.dat', f, bp, /silent
    indx=lindgen(n_elements(f)-1)
    cf[ifreq] = total( (f[indx+1]-f[indx])*(f[indx+1]*bp[indx+1]+f[indx]*bp[indx])/2 ) / $
      total( (f[indx+1]-f[indx])*(bp[indx+1]+bp[indx])/2 ) 
;# --- Considering beta /= 1.
    cfi[ifreq] = total( (f[indx+1]-f[indx])*(f[indx+1]^beta[ifreq]*bp[indx+1]+f[indx]^beta[ifreq]*bp[indx])/2 ) / $
      total( (f[indx+1]-f[indx])*(bp[indx+1]+bp[indx])/2 ) 
    print, '- Effective frequency: '+sfreq[ifreq], cf[ifreq], cfi[ifreq]^(1./beta[ifreq])
endfor

STOP
end
