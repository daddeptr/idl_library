  function ame_spectrum, nu, nu_p, m_60, nu_ref=nu_ref

     if not keyword_set(nu_ref) then nu_ref = 20.d9

     log_S     = -(m_60*alog(nu_p*1.d9)/alog(nu_p/60.d0)+2.d0)*alog(nu)     + m_60*alog(nu)^2     / (2.d0*alog(nu_p/60.d0))
     log_S_ref = -(m_60*alog(nu_p*1.d9)/alog(nu_p/60.d0)+2.d0)*alog(nu_ref) + m_60*alog(nu_ref)^2 / (2.d0*alog(nu_p/60.d0))

     ame_spectrum = exp(log_S-log_S_ref)

     return, ame_spectrum

 end
