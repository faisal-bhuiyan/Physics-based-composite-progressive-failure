SUBROUTINE DAMAGE_INC_FATIGUE(EFF_STRESS, TEMP, CMNAME, KINC, NCYCLES, NDAMAGE_SUM, NFAIL_SUM)
    !----------------------------------------------------------------------------------------------------------------------------
	! '''
	! DAMAGE_INC_FATIGUE: This function calculates the damage increment for a given number of load cycles.
	! '''

	! Import MODULEs
	USE KTF_CONSTANTS
	USE MATERIAL_PROPS

	! Get rid of inherent variable definitions
	IMPLICIT NONE

	! Define function parameters
	INTEGER, INTENT(IN) :: KINC
	CHARACTER(LEN=80), INTENT(IN) :: CMNAME
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: EFF_STRESS, TEMP, NCYCLES
	REAL(KIND=KIND(0.D0)), INTENT(INOUT) :: NDAMAGE_SUM, NFAIL_SUM

	! Define internal variables
	REAL(KIND=KIND(0.D0)) :: NDAMAGE_KTF, T, EFF_STRESS_MIN

	!----------------------------------------------------------------------------------------------------------------------------
	! Convert the Temp to Kelvin. This assumes that the zero value is room temp and in Celsius
	T = TEMP + (273.D0 + 22.D0)

	! Calculate Min Stress
	EFF_STRESS_MIN = R*EFF_STRESS
	!WRITE(*,*) 'EFF_STRESS = ', EFF_STRESS

	!******************************************************************************************
	! Damage state for delamination layer matrix
	IF (CMNAME .EQ. 'MATRIX') THEN
		IF (lambda .EQ. ONE) THEN
			! Calculate damage state for Lambda = 1
			NDAMAGE_KTF = n_0 - (n_0-NDAMAGE_SUM) * (exp(-(((K*T)**2)/H) &
			* exp(-U_MAT/(K*T)) * (NCYCLES/(GAM_MAT*F*(EFF_STRESS-EFF_STRESS_MIN))) &
			* (exp((GAM_MAT*EFF_STRESS)/(K*T)) - exp((GAM_MAT*EFF_STRESS_MIN)/(K*T)))))

		ELSEIF (lambda .EQ. NINE) THEN
			! Calculate damage state for Lambda = 9
			NDAMAGE_KTF = n_0 - ((n_0-NDAMAGE_SUM)**(ONE-lambda) - (ONE-lambda) * &
			 ((K*T)**2/H)*(NCYCLES/(GAM_MAT*(EFF_STRESS-EFF_STRESS_MIN)))* &
			 exp(-U_MAT/(K*T)) * (exp((GAM_MAT*EFF_STRESS)/(K*T)) - &
			 exp((GAM_MAT*EFF_STRESS_MIN)/(K*T))))**(ONE/(ONE-lambda))
		ENDIF

	!******************************************************************************************
	! Damage state for in-situ matrix
	ELSE
		IF (lambda .EQ. ONE) THEN
			! SDVONE: Calculate damage state for Lambda = 1
			NDAMAGE_KTF = n_0 - (n_0-NDAMAGE_SUM) * (exp(-(((K*T)**2)/H) &
			* exp(-U/(K*T)) * (NCYCLES/(GAM*F*(EFF_STRESS-EFF_STRESS_MIN))) &
			* (exp((GAM*EFF_STRESS)/(K*T)) - exp((GAM*EFF_STRESS_MIN)/(K*T)))))


			ELSEIF (lambda .EQ. NINE) THEN
				! SDVONE: Calculate damage state for Lambda = 9
				NDAMAGE_KTF = n_0 - ((n_0-NDAMAGE_SUM)**(ONE-lambda) - (ONE-lambda) * &
				 ((K*T)**2/H)*(NCYCLES/(GAM*(EFF_STRESS-EFF_STRESS_MIN)))* &
				 exp(-U/(K*T)) * (exp((GAM*EFF_STRESS)/(K*T)) - &
				 exp((GAM*EFF_STRESS_MIN)/(K*T))))**(ONE/(ONE-lambda))
		ENDIF
	ENDIF

	! Update total damage/NDAMAGE_SUM
	NDAMAGE_SUM = NDAMAGE_KTF

	! SDV2: Calculate Number of Load Cycles
	NFAIL_SUM = NFAIL_SUM + NCYCLES
	!----------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE DAMAGE_INC_FATIGUE

