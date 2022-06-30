SUBROUTINE RANDOM_WEIBULL(a_shape, b_scale, random_number_weibull)

	!----------------------------------------------------------------------------------------------------------------------------
	! '''
	! RANDOM_WEIBULL: Generates a random variate from the Weibull distribution with
	! probability density:
	!                  a
	!           a-1  -x
	! f(x) = a.x    e
	! '''

	IMPLICIT NONE

	! Function parameters
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: a_shape, b_scale
	REAL(KIND=KIND(0.D0)), INTENT(OUT) :: random_number_weibull

	! Internal variables
	REAL(KIND=KIND(0.D0)) :: x = 1.0

	!----------------------------------------------------------------------------------------------------------------------------
	CALL RANDOM_NUMBER(x)
	random_number_weibull = a_shape*(-LOG(1-(x)))**(1/b_scale)
	!----------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE RANDOM_WEIBULL




















