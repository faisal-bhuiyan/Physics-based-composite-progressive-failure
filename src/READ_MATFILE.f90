1SUBROUTINE  READ_MATERIAL(Material_file, DEG_FACTOR_MATRIX, DEG_FACTOR_DELAM, DEG_FACTOR_FIBER, C, S, C_f, S_f, C_m, S_m, &
	C_mdeg, S_mdeg, C_f_mdeg, S_f_mdeg, C_m_mdeg_ply, S_m_mdeg_ply, C_m_mdeg_delam, S_m_mdeg_delam, &
	C_fdeg, S_fdeg, C_f_fdeg, S_f_fdeg, C_m_fdeg, S_m_fdeg, &
	Q_f, Q_m, Q_f_mdeg, Q_m_mdeg, Q_f_fdeg, Q_m_fdeg, psi_m, psi_mdeg, I)

	!----------------------------------------------------------------------------------------------------------------------------
	! '''
	!  READ_MATERIAL: Reads material data from a text file.
	! Properties read:
	! C, S, C_f, S_f, C_m, S_m
	! C_mdeg, S_mdeg, C_f_mdeg, S_f_mdeg
	! C_m_mdeg_ply, S_m_mdeg_ply, C_m_mdeg_delam, S_m_mdeg_delam
	! C_fdeg, S_fdeg, C_f_fdeg, S_f_fdeg, C_m_fdeg, S_m_fdeg
	! Q_f, Q_m, Q_f_mdeg, Q_m_mdeg, Q_f_fdeg, Q_m_fdeg
	! psi_m, psi_mdeg, I
	! '''

    ! Get rid of inherent variable definitions
	IMPLICIT NONE

	! Define function parameters
    CHARACTER(*), INTENT(IN) :: Material_file
    REAL(KIND=KIND(0.D0)), INTENT(IN) :: DEG_FACTOR_MATRIX, DEG_FACTOR_DELAM, DEG_FACTOR_FIBER
	REAL(KIND=KIND(0.D0)), INTENT(OUT), DIMENSION(6,6) :: C, S, C_f, S_f, C_m, S_m
    REAL(KIND=KIND(0.D0)), INTENT(OUT), DIMENSION(6,6) :: C_mdeg, S_mdeg, C_f_mdeg, S_f_mdeg, &
                                                          C_m_mdeg_ply, S_m_mdeg_ply, C_m_mdeg_delam, S_m_mdeg_delam
	REAL(KIND=KIND(0.D0)), INTENT(OUT), DIMENSION(6,6) :: C_fdeg, S_fdeg, C_f_fdeg, S_f_fdeg, C_m_fdeg, S_m_fdeg
	REAL(KIND=KIND(0.D0)), INTENT(OUT), DIMENSION(6,6) :: Q_f, Q_m, Q_f_mdeg, Q_m_mdeg, Q_f_fdeg, Q_m_fdeg, I
	REAL(KIND=KIND(0.D0)), DIMENSION(6) :: psi_m, psi_mdeg

	!----------------------------------------------------------------------------------------------------------------------------
	OPEN(unit=15, file=Material_file, status='OLD', action='READ')

	READ(15,*)
	READ(15,*) C
	READ(15,*)
	READ(15,*)
	READ(15,*) S
	READ(15,*)
	READ(15,*)
	READ(15,*) C_f
	READ(15,*)
	READ(15,*)
	READ(15,*) S_f
	READ(15,*)
	READ(15,*)
	READ(15,*) C_m
	READ(15,*)
	READ(15,*)
	READ(15,*) S_m
	READ(15,*)
	READ(15,*)
	READ(15,*) C_mdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) S_mdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) C_f_mdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) S_f_mdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) C_m_mdeg_ply
	READ(15,*)
	READ(15,*)
	READ(15,*) S_m_mdeg_ply
	READ(15,*)
	READ(15,*)
	READ(15,*) C_m_mdeg_delam
	READ(15,*)
	READ(15,*)
	READ(15,*) S_m_mdeg_delam
	READ(15,*)
	READ(15,*)
	READ(15,*) C_fdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) S_fdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) C_m_fdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) S_m_fdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) C_f_fdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) S_f_fdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) Q_f
	READ(15,*)
	READ(15,*)
	READ(15,*) Q_m
	READ(15,*)
	READ(15,*)
	READ(15,*) Q_f_mdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) Q_m_mdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) Q_f_fdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) Q_m_fdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) psi_m
	READ(15,*)
	READ(15,*)
	READ(15,*) psi_mdeg
	READ(15,*)
	READ(15,*)
	READ(15,*) I
	CLOSE(15)

	C = 1.0e+5 * C
	S = 1.0e-3 * S
	C_f = 1.0e+5 * C_f
	S_f = 1.0e-3 * S_f
	C_m = 1.0e+3 * C_m
    S_m = 1.0e-3 * S_m

    ! Define matrix-degraded composite layer properties based on DEG_FACTOR_MATRIX
    IF (DEG_FACTOR_MATRIX .EQ. 1.0e-1) THEN
        C_mdeg = 1.0e+5 * C_mdeg
		C_f_mdeg = 1.0e+5 * C_f_mdeg
		C_m_mdeg_ply = 1.0e+2 * C_m_mdeg_ply
		! Multiplying factors might not be correct for compliance
		S_mdeg = 1.0 * S_mdeg
        S_f_mdeg = 1.0e-3 * S_f_mdeg
		S_m_mdeg_ply = 1.0 * S_m_mdeg_ply
    ENDIF

    ! Define degraded delam layer properties based on DEG_FACTOR_DELAM
    IF (DEG_FACTOR_DELAM .EQ. 1.0e-4) THEN
		C_m_mdeg_delam = 1.0e-2 * C_m_mdeg_delam
		! Multiplying factors might not be correct for compliance
        S_m_mdeg_delam = 1.0e+2 * S_m_mdeg_delam

    ELSE IF (DEG_FACTOR_DELAM .EQ. 1.0e-2) THEN
		C_m_mdeg_delam = 1.0 * C_m_mdeg_delam
		! Multiplying factors might not be correct for compliance
        S_m_mdeg_delam = 1.0e+4 * S_m_mdeg_delam
    ENDIF

    ! Define fiber-degraded composite layer properties based on DEG_FACTOR_FIBER
    IF (DEG_FACTOR_FIBER .EQ. 1.0e-4) THEN
        C_fdeg = 1.0e+1 * C_fdeg
		C_m_fdeg = 1.0e-2 * C_m_fdeg
		C_f_fdeg = 1.0e+1 * C_f_fdeg
		! Multiplying factors might not be correct for compliance
		S_fdeg = 1.0 * S_fdeg
        S_m_fdeg = 1.0e+1 * S_m_fdeg
		S_f_fdeg = 1.0e+1 * S_f_fdeg

	ELSE IF (DEG_FACTOR_FIBER .EQ. 1.0e-2) THEN
		C_fdeg = 1.0e+3 * C_fdeg
		C_m_fdeg = 1.0 * C_m_fdeg
		C_f_fdeg = 1.0e+3 * C_f_fdeg
		! Multiplying factors might not be correct for compliance
		S_fdeg = 1.0 * S_fdeg
		S_m_fdeg = 1.0 * S_m_fdeg
		S_f_fdeg = 1.0 * S_f_fdeg
    ENDIF

	Q_f = TRANSPOSE(Q_f)
	Q_m = TRANSPOSE(Q_m)
	Q_f_mdeg = TRANSPOSE(Q_f_mdeg)
	Q_m_mdeg = TRANSPOSE(Q_m_mdeg)
	Q_f_fdeg = TRANSPOSE(Q_f_fdeg)
	Q_m_fdeg = TRANSPOSE(Q_m_fdeg)
	!----------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE

