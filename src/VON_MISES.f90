SUBROUTINE VON_MISES(STRESS, SMISES)

	!----------------------------------------------------------------------------------------------------------------------------
	! '''
	! VON_MISES: Calculates von Mises stress for a provided stress tensor. Used as effective stress for an
	! isotropic material.
	! '''

	! Get rid of inherent variable definitions
	IMPLICIT NONE

	! Define function parameters
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: STRESS(6)
	REAL(KIND=KIND(0.D0)), INTENT(OUT) :: SMISES

	! Define internal variables
	REAL(KIND=KIND(0.D0)), PARAMETER :: ZERO=0.D0, ONE=1.D0
	REAL(KIND=KIND(0.D0)), PARAMETER :: TWO=2.D0, THREE=3.D0
	REAL(KIND=KIND(0.D0)), PARAMETER :: SIX=6.D0

	!----------------------------------------------------------------------------------------------------------------------------
	! Calculate Von Mises Stress, load scenario: general
	SMISES = (STRESS(1) - STRESS(2))**2 + (STRESS(2)- STRESS(3))**2 &
			+ (STRESS(3) - STRESS(1))**2 + SIX*(STRESS(4)**2 + STRESS(5)**2 &
			+ STRESS(6)**2)

	SMISES = SQRT(SMISES/TWO)
	!----------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE VON_MISES

