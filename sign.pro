function sign, x
    temx = x*0.
    gx = where(x ne 0.)
    temx[gx] = x[gx] / abs(x[gx])
    return, long(temx)
end

