function asinh10, x

   x = double(x)
   y = alog10( 0.5d0*(x+sqrt(x^2+4.d0)) )

   return, y
end
