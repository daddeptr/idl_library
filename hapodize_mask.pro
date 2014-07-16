pro hapodize_mask, in_maskfile, feedback=feedback, out_maskfile=out_maskfile, sigma=sigma

   openw,1,'apodize_mask.par.tmp'
   printf,1,'feedback='+strtrim(string(feedback),2)
   printf,1,'in_maskfile='+in_maskfile
   if keyword_set(out_maskfile) then printf,1,'out_maskfile='+out_maskfile
   if keyword_set(sigma) then printf,1,'sigma='+sigma
   close,1
   spawn, 'head apodize_mask.par.tmp'
   spawn, '/global/scratch2/sd/dpietrob/Software/src/f90/apodize_mask/apodize_mask.x apodize_mask.par.tmp'

end
