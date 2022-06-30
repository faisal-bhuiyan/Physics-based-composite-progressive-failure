SUBROUTINE NCYCLES_TO_FAILURE_CALC(STRAN, DSTRAN, TEMP, NDAMAGE_SUM, NCYCLES_FAIL)
	!----------------------------------------------------------------------------------------------------------------------------
	! '''
	! NCYCLES_TO_FAILURE_CALC: This function calculates the remaining number of cycles to failure from an arbitrary damage state, n.
	! It is used to calculate a dynamic value for cycle jump at each increment of the fatigue step.
	! '''

	! Import MODULEs
	USE KTF_CONSTANTS
	USE MATERIAL_PROPS

	! Get rid of inherent variable definitions
	IMPLICIT NONE

	! Define function parameters
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: STRAN(6),DSTRAN(6),TEMP,NDAMAGE_SUM
	REAL(KIND=KIND(0.D0)), INTENT(OUT) :: NCYCLES_FAIL

	! Define internal variables
	REAL(KIND=KIND(0.D0)) :: STRESS(6)
	REAL(KIND=KIND(0.D0)) :: EFF_STRESS, EFF_STRESS_MIN, sigma_m(6), sigma_m_eff, T
	REAL(KIND=KIND(0.D0)), PARAMETER :: MAX_NCYCLES_FAIL = 1.0e+8, CLOSE_TO_ZERO = 1.0e-4
	!----------------------------------------------------------------------------------------------------------------------------

	! Convert the Temp to Kelvin. This assumes that the zero value is room temp and in Celsius
	T = TEMP + (273.D0 + 22.D0)

    !******************************************************************************************
	! Calculate current effective matrix stress, i.e. EFF_MATRIX_STRESS

	IF (NDAMAGE_SUM .LT. ONE) THEN
		STRESS = MATMUL(C, STRAN + DSTRAN)

		! Calculate effective matrix stress
		sigma_m = MATMUL(Q_m, STRESS) - (TEMP * psi_m)
		! Effective Matrix Stress
		CALL EFF_MATRIX_STRESS(B_t, B_s1, B_s2, beta, tau_0, sigma_m, sigma_m_eff)
		EFF_STRESS = sigma_m_eff

	ELSE
		STRESS = MATMUL(C_mdeg, STRAN + DSTRAN)

		! Calculate effective matrix stress
		sigma_m = MATMUL(Q_m_mdeg, STRESS) - (TEMP * psi_mdeg)
		CALL EFF_MATRIX_STRESS(B_t_mdeg, B_s1_mdeg, B_s2_mdeg, beta, tau_0, sigma_m, sigma_m_eff)
		! Define Effective Stress
		EFF_STRESS = sigma_m_eff
	ENDIF

	! Calculate EFF_STRESS_MIN
	EFF_STRESS_MIN = R * EFF_STRESS

	!******************************************************************************************
	! Calculate number of cycles to failure based on current damage state

	IF (EFF_STRESS_MIN .GT. CLOSE_TO_ZERO) THEN
		IF (lambda .EQ. ONE) THEN
			! Calculate # of cycles to failure for Lambda = 1
			NCYCLES_FAIL = (H/(K*T)**2) * (GAM*F*(EFF_STRESS-EFF_STRESS_MIN)) &
			* exp(U/(K*T)) * (exp((GAM*EFF_STRESS)/(K*T)) - &
			exp((GAM*EFF_STRESS_MIN)/(K*T)))**(-ONE) * log(exp(ONE)- NDAMAGE_SUM * (exp(ONE)-ONE))

		ELSEIF (lambda .EQ. NINE) THEN
			! Calculate # of cycles to failure for Lambda = 9
			NCYCLES_FAIL = (ONE/(ONE-lambda)) * (H/(K*T)**2) * (GAM*F*(EFF_STRESS-EFF_STRESS_MIN)) &
			* exp(U/(K*T)) * (exp((GAM*EFF_STRESS)/(K*T)) - exp((GAM*EFF_STRESS_MIN)/(K*T)))**(-ONE) &
			* ((n_0-NDAMAGE_SUM)**(ONE-lambda) - (ONE-n_0)**(ONE-lambda))
		ENDIF

	ELSEIF (EFF_STRESS_MIN .LT. CLOSE_TO_ZERO .AND. EFF_STRESS_MIN .GT. ZERO) THEN
		NCYCLES_FAIL = MAX_NCYCLES_FAIL

	ELSE
		NCYCLES_FAIL = -ONE
    ENDIF

	! Set max value of NCYCLES_FAIL as MAX_NCYCLES_FAIL
	IF (NCYCLES_FAIL .GT. MAX_NCYCLES_FAIL) THEN
		NCYCLES_FAIL = MAX_NCYCLES_FAIL
	ENDIF

	! Set NCYCLES_FAIL as -1.D0 (required for INDEX_FINDER logic) when it's < 1.0 and NaN
	IF (NCYCLES_FAIL .LT. ONE .OR. NCYCLES_FAIL .NE. NCYCLES_FAIL ) THEN
		NCYCLES_FAIL = -ONE
	ENDIF
	!----------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE NCYCLES_TO_FAILURE_CALC

