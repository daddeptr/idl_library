pro swap_eenoise_with_bbsignal, sfile, nfile, ofile, smoothing=smoothing

   fits2cl, cls, sfile
   fits2cl, nls, nfile

   fls = nls
   bb = cls[*,2]
   bbs = bb
   if keyword_set(smoothing) then bbs = smooth(bb,smoothing)
   fls[*,1] = bbs

   cl2fits, fls, ofile

end
