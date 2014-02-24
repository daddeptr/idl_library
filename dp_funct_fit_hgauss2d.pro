FUNCTION DP_FUNCT_FIT_HGAUSS2D, A, ROBUST = ROBUST, DATA = DATA, ERR_DATA = ERR_DATA, $
                            N_PAR = N_PAR, LIST_PIX = LIST_PIX, HGAUSS2D = HGAUSS2D, $
                            L_MAX = L_MAX, COB_MAX = COB_MAX, N_SIDE = N_SIDE, $
                            L_CENTR = L_CENTR, COB_CENTR = COB_CENTR, src_model=src_model

; Need to up-grade to higher n_side in NESTED.
; DATA has: DATA, ERR_DATA, LIST_PIX, L_MAX, COB_MAX, N_SIDE, NESTED
; Better leave room for a background of different sign that the source ( Q or U may be < 0, while >0 in I, and necessarily in P )
; PARAMETERS
; 0: Background, constant
; 1: Amplitude of the 2D Gaussian
; 2: FWHM_LONGITUDE
; 3: FWHM_LATITUDE
; 4: CENTER_LONGITUDE
; 5: CENTER_LATITUDE
; 6: TILT ELLIPSE 

; print, A

   pix2ang_ring, n_side, list_pix, COB_PIX, L_PIX
   N_PIX_FIT = n_elements( l_pix )
   L_CENTR = l_pix - replicate( l_max[ 0 ], n_pix_fit )
   COB_CENTR = cob_pix - replicate( cob_max[ 0 ], n_pix_fit ) 
   
   X_ROT = ( l_centr - A[ 4 ] ) * cos( A[ 6 ] ) - ( cob_centr - A[ 5 ] ) * sin( A[ 6 ] ) ; arcmin
   Y_ROT = ( cob_centr - A[ 5 ] ) * cos( A[ 6 ] ) + ( l_centr - A[ 4 ] ) * sin( A[ 6 ] ) ; arcmin

   FWHM_X_ROT = A[2]            ; arcmin
   FWHM_Y_ROT = A[3]            ; arcmin  

   X_ROT_FWHM = 2d0 * x_rot / fwhm_x_rot
   Y_ROT_FWHM = 2d0 * y_rot / fwhm_y_rot

 ;; Maybe useful for the future.
 ;;Q_SAMPLES = where(abs(x_rot_fwhm) lt 10)
 ;;  If q_samples(0) eq -1 then GOTO,END_SRC_MODEL
 ;;Q_SAMPLES_2 = where(abs(y_rot_fwhm(q_samples)) lt 30) ; should be dependent on the S/N 
 ;;                   ;of the source (For Jupiter is S/N 300 for each circle.
 ;;  If q_samples_2(0) eq -1 then GOTO,END_SRC_MODEL

;;END_SAMPLES = q_samples(q_samples_2)

;;  EXP_X_ROT = - x_rot_fwhm(end_samples) * x_rot_fwhm(end_samples)
;;  EXP_Y_ROT = - y_rot_fwhm(end_samples) * y_rot_fwhm(end_samples)

;; SRC_MODEL(end_samples) = 2d0^(exp_x_rot) * 2d0^(exp_y_rot)

   EXP_X_ROT = - x_rot_fwhm * x_rot_fwhm
   EXP_Y_ROT = - y_rot_fwhm * y_rot_fwhm
   SRC_MODEL = 2d0^(exp_x_rot) * 2d0^(exp_y_rot)
   HGAUSS2D = a[ 0 ] +  a[ 1 ] * src_model
; Including the linear background option
   If (N_PAR eq 9) then HGAUSS2D = hgauss2d + a[ 7 ] * x_rot_fwhm / max( x_rot_fwhm ) + a[ 8 ] * y_rot_fwhm / max( y_rot_fwhm )
;  If N_PAR eq 9 then HGAUSS2D = hgauss2d + a[ 7 ] * x_rot_fwhm + a[ 8 ] * y_rot_fwhm 
; Including the parabolic background
   If (N_PAR eq 12) then HGAUSS2D = hgauss2d + a[ 7 ] * x_rot_fwhm / max( x_rot_fwhm ) + a[ 8 ] * y_rot_fwhm / max( y_rot_fwhm ) - $
     a[ 9 ] * exp_x_rot / max( x_rot_fwhm )^2d0 - $
     a[ 11 ] * exp_y_rot / max( y_rot_fwhm )^2d0 + $
     a[ 10 ] * x_rot_fwhm * y_rot_fwhm / max( x_rot_fwhm ) / max( y_rot_fwhm )
;;    If (N_PAR eq 12) then HGAUSS2D = hgauss2d - a[ 9 ] * exp_x_rot / max( x_rot_fwhm )^2d0 - $
;;      a[ 11 ] * exp_y_rot / max( y_rot_fwhm )^2d0 + $
;;      a[ 10 ] * x_rot_fwhm * y_rot_fwhm / max( x_rot_fwhm ) / max( y_rot_fwhm )
; computing the statistics
   RESULT = ( hgauss2d - data ) / ( err_data )
;   If ( Keyword_set( ROBUST ) ) then RESULT = sqrt( abs( result ) )
;stop, robust
   If ( ROBUST ) then begin
       RESULT = sqrt( abs( result ) )
;       print, 'robust: ', robust
   endif

   return, result

END
