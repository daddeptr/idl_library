pro write_newdat, out_file, band=band,spec_file=spec_file, cov_file=cov_file

   if not defined(spec_file) then spec_file='/global/scratch2/sd/marius/us2/runs/dx11c_allmasks/for_davide/spectrum.txt'
   if not defined(cov_file) then cov_file='/global/scratch2/sd/marius/us2/runs/dx11c_allmasks/for_davide/cov.txt'
   if not defined(band) then band=[12,61]

   readcol,spec_file, ibin, lmn, lmx, dl
   nc = n_elements(ibin)
   cov = fltarr(nc,nc)
   openr,1,cov_file
   readf,1,cov
   close,/all
   
   openw,1,out_file
   tag = strsplit(out_file,'.',/extract)
   printf,1,tag[0]+'_'
   printf,1,nc,0,0,0,0,0,format='(6i10)'
   printf,1,'BAND_SELECTION'
   printf,1,band[0],band[1], format='(2i4)'
   for i=0,4 do printf,1,0,0, format='(2i4)'
   printf,1,'1   1   0'
   printf,1,'0   0   0'
   printf,1,'0   #iliketype'
   printf,1,'TT'
   for i=0,nc-1 do printf,1,i+1,dl[i],sqrt(cov[i,i]),sqrt(cov[i,i]),0.,lmn[i], lmx[i],1,format='(i10,3f15.6,1f5.1,3i10)'
   for i=0,nc-1 do printf,1,fltarr(nc),format='('+strtrim(string(nc))+'f7.3)'
   for i=0,nc-1 do printf,1,cov[i,*],format='('+strtrim(string(nc))+'e13.3)'
   close,1

   

end
