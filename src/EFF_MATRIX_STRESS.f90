SUBROUTINE EFF_MATRIX_STRESS(B_t,B_s1,B_s2,beta,tau_0,sigma_m,sigma_m_eff)
    
	!----------------------------------------------------------------------------------------------------------------------------
	! '''
	! EFF_MATRIX_STRESS: This function calculates the effective matrix state given the decomposed stress vector to the matrix.
	! Fertig failure theory is used for this purpose. 
	! '''
	
    IMPLICIT NONE
	
    ! Define function parameters
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: B_t, B_s1, B_s2, beta, tau_0, sigma_m(6)
	REAL(KIND=KIND(0.D0)), INTENT(OUT) :: sigma_m_eff
	
	! Define internal variables
	REAL(KIND=KIND(0.D0)), PARAMETER :: ONE = 1.D0, FOUR = 4.D0
	REAL(KIND=KIND(0.D0)) :: I_t, I_s1, I_s2, I_h, MACAULAY_I_t, MACAULAY_I_h
    
	!----------------------------------------------------------------------------------------------------------------------------
	! Transversely isotropic stress invariants 
	I_t = ( sigma_m(2) + sigma_m(3) + sqrt( (sigma_m(2)+sigma_m(3))**2-4*(sigma_m(2)*sigma_m(3)-sigma_m(6)) ) )/2
	I_s1 = sigma_m(4)**2 + sigma_m(5)**2
	I_s2 = (ONE/FOUR) * (sigma_m(2)-sigma_m(3))**2 + sigma_m(6)**2
	I_h = sigma_m(2) + sigma_m(3)
	
	! Effective matrix stress
	CALL MACAULAY(I_t, MACAULAY_I_t)
	CALL MACAULAY(-I_h, MACAULAY_I_h)
	
	! ** Debug print **
	!WRITE(*,*) 'MACAULAY_I_t: ', MACAULAY_I_t
	!WRITE(*,*) 'MACAULAY_I_h: ', MACAULAY_I_h
	
	sigma_m_eff = sqrt( (B_t/B_s1) * (MACAULAY_I_t)**2 + 1/(1+((beta/tau_0)* MACAULAY_I_h)) * (I_s1 + ((B_s2/B_s1)*I_s2)) )
	!----------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE EFF_MATRIX_STRESS 

INCLUDE 'MACAULAY.f90'

