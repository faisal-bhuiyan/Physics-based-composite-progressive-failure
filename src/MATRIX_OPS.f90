SUBROUTINE GET_STIFFNESS(E, G, NU, C, S)

    !----------------------------------------------------------------------------------------------------------------------------
	! '''
	! GET_STIFFNESS: Defines stiffness matrix C and compliance matrix S based on provided elastic constants.
	! '''

	IMPLICIT NONE
    !Variable definitions
    REAL(KIND=KIND(0.D0)), INTENT(IN) :: E(3), G(3), NU(3)

    REAL(KIND=KIND(0.D0)), INTENT(OUT) :: C(6,6), S(6,6)

    INTEGER :: ErrFlg

    !----------------------------------------------------------------------------------------------------------------------------
    S = 0.D0
    S(1,1) = 1/E(1)
    S(2,2) = 1/E(2)
    S(3,3) = 1/E(3)
    S(4,4) = 1/G(1)
    S(5,5) = 1/G(2)
    S(6,6) = 1/G(3)
    S(1,2) = -NU(1)/E(1)
    S(1,3) = -NU(2)/E(1)
    S(2,3) = -NU(3)/E(2)
    S(2,1) = S(1,2)
    S(3,1) = S(1,3)
    S(3,2) = S(2,3)

    C = S

    CALL FINDInv(S, C, 6, ErrFlg)
	!----------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE

!----------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
!----------------------------------------------------------------------------------------------------------------------------
	! Subroutine to find the inverse of a square matrix
	! Author : Louisda16th a.k.a Ashwith J. Rego
	! Reference : Algorithm has been well explained in:
	! http://math.uww.edu/~mcfarlat/inverse.htm
	! http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html

	IMPLICIT NONE
	!Declarations
	INTEGER, INTENT(IN) :: n
	INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
	REAL(KIND=KIND(0.D0)), INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
	REAL(KIND=KIND(0.D0)), INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix

	LOGICAL :: FLAG = .TRUE.
	INTEGER :: i, j, k, l
	REAL(KIND=KIND(0.D0)) :: m
	REAL(KIND=KIND(0.D0)), DIMENSION(n,2*n) :: augmatrix !augmented matrix
	!----------------------------------------------------------------------------------------------------------------------------

	! Augment input matrix with an identity matrix
	DO i = 1, n
		DO j = 1, 2*n
			IF (j <= n ) THEN
				augmatrix(i,j) = matrix(i,j)
			ELSE IF ((i+n) == j) THEN
				augmatrix(i,j) = 1
			Else
				augmatrix(i,j) = 0
			ENDIF
		END DO
	END DO

	!Reduce augmented matrix to upper traingular form
	DO k =1, n-1
		IF (augmatrix(k,k) == 0) THEN
			FLAG = .FALSE.
			DO i = k+1, n
				IF (augmatrix(i,k) /= 0) THEN
					DO j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					END DO
					FLAG = .TRUE.
					EXIT
				ENDIF
				IF (FLAG .EQV. .FALSE.) THEN
					PRINT*, "Matrix is non - invertible"
					inverse = 0
					errorflag = -1
					return
				ENDIF
			END DO
		ENDIF
		DO j = k+1, n
			m = augmatrix(j,k)/augmatrix(k,k)
			DO i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			END DO
		END DO
	END DO

	!Test for invertibility
	DO i = 1, n
		IF (augmatrix(i,i) == 0) THEN
			PRINT*, "Matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		ENDIF
	END DO

	!Make diagonal elements as 1
	DO i = 1 , n
		m = augmatrix(i,i)
		DO j = i , (2 * n)
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		END DO
	END DO

	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
		DO i =1, k
		m = augmatrix(i,k+1)
			DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			END DO
		END DO
	END DO

	!store answer
	DO i =1, n
		DO j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		END DO
	END DO
	errorflag = 0
	!----------------------------------------------------------------------------------------------------------------------------

END SUBROUTINE FINDinv

