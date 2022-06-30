SUBROUTINE MACAULAY(input_value, output_value)
	!----------------------------------------------------------------------------------------------------------------------------
	! '''
	! MACAULAY: Performs the MacAulay bracket operation on a term.
	! '''

    IMPLICIT NONE

    !Variable definitions
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: input_value
	REAL(KIND=KIND(0.D0)), INTENT(OUT) :: output_value

	!----------------------------------------------------------------------------------------------------------------------------
	! Definition of MACAULAY bracket
	IF (input_value>0) THEN
		output_value = input_value
	ELSE
		output_value = 0.D0
	ENDIF
	!----------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE

