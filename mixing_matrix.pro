pro mixing_matrix, solve_for=solve_for, ichannels=ichannels, indices=indices

; Ruler: linear solution in pixel space of the component separation
; problem
;
; D.Pietrobon Jul 2011
;
; To be improved and generalized to be flexible on foreground model
; and translated into f90 for speed
;
; chisq = sum_ch (D_ch-M_ch(A_i))^2/N_ch^2
;
; dchisq/dA_i = 0 = \sum_ch 1/N_ch^2 (D_ch-\sum_j A_j*Ind_ch,j)*I_ch,i
;
; M # A = B
;
; A = (cmb, synch, ff, dust)
;
; B = \sum_ch D_ch/N_ch^2 * Ind_ch ; Ind_ch array of spectral
; responses

; --- Foreground Model definition
 
false = 0b
true  = 1b

if (not keyword_set(solve_for)) and (not keyword_set(ichannels)) then begin
    print,'Syntax - ruler, solve_for=["cmb","dust","co"], ichannel=[1,2,3]'
    return
 endif

if ( not keyword_set(solve_for) ) then solve_for = ['cmb', 'dust', 'co']
if ( not keyword_set(ichannels) ) then ichannels = [2, 3, 4]
if ( not keyword_set(indices) ) then indices = [-3, 1.8, 18]

code = 'ruler'
temp_folder = '/global/scratch/sd/dpietrob/rubbish/'

cfreq = [30, 44, 70, 100, 143, 217, 353, 545, 857, 23, 33, 41, 61, 94]
; cCO_spec= [1, 0, .6, .2, 0., 0.]
cCO_spec= [0., 0., 0., 1, 0, 0.35, .12, 0., 0., 0., 0., 0., 0., 0.] ; new values

freq = cfreq[ichannels]
sfreq = string(freq,format='(i3.3)')
nfreq = n_elements(freq)
h = 6.626068 * 10.^(-34)
k = 1.3806503 * 10.^(-23)
c = 3.d8

convfact = conversionfactor(freq, /antenna2thermo)
print, convfact

sync_pivot = 30.
dust_pivot = 353.

co_spec = cCO_spec[ichannels]

beta = indices[0]
emis = indices[1]
td = indices[2]

ncomp = n_elements(solve_for) 

components = fltarr(ncomp)

rms = fltarr(nfreq)
rms[*] = 1.

; --- System Solution

   spectral_response = dblarr(ncomp,nfreq)

   B = dblarr(ncomp)
   A = dblarr(ncomp, ncomp)

   ; to be improved with routines
   for icomp=0,ncomp-1 do begin
; cmb
       if ( solve_for[icomp] eq 'cmb' ) then   spectral_response[icomp,*] = 1. 
;        print, reduced_data
; thermal dust
       if ( solve_for[icomp] eq 'dust' ) then begin
           x = h*1.d9/k/td
           bb  = 1. / ( exp(x*freq[*])-1. )
           bb0 = 1. / ( exp(x*dust_pivot)-1. )
           spectral_response[icomp,*] = convfact[*] * (freq[*]/dust_pivot)^(emis+1)*bb/bb0 
       endif
; CO
       if ( solve_for[icomp] eq 'co' ) then spectral_response[icomp,*] = convfact[*] * co_spec[*] 
; sync
       if ( solve_for[icomp] eq 'sync' ) then spectral_response[icomp,*] = convfact[*] * (freq[*]/sync_pivot)^beta
   endfor

   for i=0,ncomp-1 do begin
      for j=i,ncomp-1 do begin
	for ifreq=0,nfreq-1 do begin
            A[i,j] = A[i,j] + spectral_response[i,ifreq]*spectral_response[j,ifreq] / rms[ifreq]^2
        endfor
            A[j,i] = A[i,j]
      endfor
   endfor
   
   Am1 = la_invert(A,/double, status=info)

   print, 'check:'
   print, 'A ', A
   print, ''
   print, 'Am1 ', Am1
   eval = EIGENQL(A, /double, EIGENVECTORS=evec)
   print, ''
   print, 'EVEC ', evec
   print, ''
   print, 'EVAL', eval


   evalm1 = HQR(ELMHES(Am1), /DOUBLE)
   evecm1 = EIGENVEC(Am1, evalm1, /double)
   print, ''
   print, 'EVECm1 ', evecm1
   print, ''
   print, 'EVALm1', evalm1
   stop
;;    weights[ipix,*,*] = float( Am1 )

   X = Am1 ## transpose(B)
   
   components[*] = X
;stop


end
