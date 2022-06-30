!DEC$ FREEFORM
!----------------------------------------------------------------------------------------------------------------------------
	! SUBROUTINE UMAT: User defined material routine for progressive composite failure under fatigue loading.

    ! Release notes/ Usage instructions:

	!** Be familiar with the physics-based, multiscale law based on KTF from the paper: Bhuiyan 18. It provides
	!   a comprehensive description of the material properties required by this routine.

	!** Define materials properties in "MODULE KTF_CONSTANTS" and "MODULE MATERIAL_PROPS", in addition to the
	!   constituent and homogenized properties supplied via text file "MAT_FILE". These properties can be
	!   computed using the provided Matlab routines, which require lamina level elastic properties and fatigue
	!   test data.

	!** Read through the usage instruction on individual SUBROUTINES in this file as well as the included functions.

	!** Run subroutine on the "unit tests" provided - especially the one element and the 140 element models to ensure things
	!   are working as intended.

	!** Uncomment the existing WRITE statements as required to print debug statements.

!----------------------------------------------------------------------------------------------------------------------------
! MODULE definitions
!----------------------------------------------------------------------------------------------------------------------------
MODULE KTF_CONSTANTS

	! ** Instructions: Special attention should be paid to ensure consistent units **
	! ** Here, the N-mm-s system is being used **

	! ** Tech debt **: These props should be ideally read in from the materials file along with the other props.

	! Planck's Constant [kJ-s/mol] & Boltzmann Constant [kJ/mol-K]
	REAL(KIND=KIND(0.D0)), PARAMETER :: H = 3.9903e-13, K = 0.008314

	! R (Stress ratio, sigma_max/sigma_min), F (Frequency) and lambda (shape parameter of damage evolution curve)
	REAL(KIND=KIND(0.D0)), PARAMETER :: R = 0.1, F = 10
	REAL(KIND=KIND(0.D0)), PARAMETER :: lambda = 1.D0

	! U (Activation Energy [kJ/mol] and GAM (Activation Volume [kJ/MPa-mol]): intralaminar/in-situ matrix
    REAL(KIND=KIND(0.D0)), PARAMETER :: U = 137.6304,	GAM = 1.1775

	! U (Activation Energy [kJ/mol] and GAM (Activation Volume [kJ/MPa-mol]): interlaminar/delamination layer matrix
    REAL(KIND=KIND(0.D0)), PARAMETER :: U_MAT = 137.6304,	GAM_MAT = 0.7372

	! Constant for KTF based damage calculation: n_0
	REAL(KIND=KIND(0.D0)), PARAMETER :: n_0 = 1.7710

	!******************************************************************************************
	! Variables used for Adaptive Cycle Jumping

	! ** Instructions: These variables should be calculated with an initial run of the UMAT, without the adaptive cycle jumping.
	! ** For example, dimension of UNSORTED_NCYCLES_FAIL_MATRIX and SORTED_NCYCLES_FAIL_MATRIX array is the number of elements in problem.
	REAL(KIND=KIND(0.D0)) :: NCYCLES_INC = 1.0e+0
    INTEGER*8, DIMENSION(1) :: INDEX_FINDER_MATRIX = -1, FAIL_PCT_MATRIX = 1, INDEX_FINDER_FIBER = -1, FAIL_PCT_FIBER = 4
    INTEGER*8, PARAMETER :: BASE_FAILURE_PCT = 1000
    REAL(KIND=KIND(0.D0)), PARAMETER :: MIN_CYCLE_JUMP = 4.0e4
    REAL(KIND=KIND(0.D0)) :: NELEM = 608056, NELEM_LE = 336800, NELEM_UMAT = 268776, NELEM_UMAT_DELAM = 117240

	! Initiate UNSORTED_NCYCLES_FAIL_MATRIX to -1.0, DIMENSION is the number of elements in the model
    REAL(KIND=KIND(0.D0)), DIMENSION(608056) :: UNSORTED_NCYCLES_FAIL_MATRIX = -1.D0, SORTED_NCYCLES_FAIL_MATRIX
	REAL(KIND=KIND(0.D0)), DIMENSION(608056) :: UNSORTED_NCYCLES_FAIL_FIBER = -1.D0, SORTED_NCYCLES_FAIL_FIBER
	REAL(KIND=KIND(0.D0)), PARAMETER :: MAX_NCYCLES_FAIL_FIBER = 1.0e+9

    ! QSORT parameters
	INTEGER*8 :: LEN = 608056, ISIZE = 8
    external COMPAR
    REAL(KIND=KIND(0.D0)) :: MIN_NCYCLES_FAIL_MATRIX, ADAPTIVE_NCYCLES_MATRIX
    REAL(KIND=KIND(0.D0)) :: MIN_NCYCLES_FAIL_FIBER, ADAPTIVE_NCYCLES_FIBER, ADAPTIVE_NCYCLES

    !******************************************************************************************
    ! Variables used for Material Props input and debugging files
    CHARACTER(*), PARAMETER :: MAT_FILE='/scratch/fbhuiyan/Fiber Failure/Properties_Exp1.txt'
    CHARACTER(*), PARAMETER :: LOG_FILE='/scratch/fbhuiyan/Fiber Failure/CycleJump_Mat1000Elem_FF4Elem_DelSD0pt0001_DelGam0pt74.txt'

END MODULE

!----------------------------------------------------------------------------------------------------------------------------
MODULE MATERIAL_PROPS

	! ** Tech debt **: These props should be ideally read in from the materials file along with the other props.

	! Strength properties
    REAL(KIND=KIND(0.D0)), PARAMETER :: B_t = 1.802093586182320e-04, B_s1 = 2.582965325467884e-04, &
        B_s2 = 2.865521091325947e-04, beta = 0.35, tau_0 = 62.221533715312184
	REAL(KIND=KIND(0.D0)), PARAMETER :: B_t_mdeg = 1.8808e-04, B_s1_mdeg = 2.5624e-04, B_s2_mdeg = 4.7118e-04
	REAL(KIND=KIND(0.D0)), PARAMETER :: B_t_fdeg = 1.5457e-04, B_s1_fdeg = 2.2586e-04, B_s2_fdeg = 2.2591e-04

	REAL(KIND=KIND(0.D0)), PARAMETER :: phi_f = 0.685, phi_m = 1-phi_f
	REAL(KIND=KIND(0.D0)), PARAMETER, DIMENSION(6) :: eta_f=0.D0, eta_m=0.D0, eta_c=0.D0
	REAL(KIND=KIND(0.D0)), PARAMETER :: ZERO = 0.D0, ONE = 1.D0, NINE = 9.D0, ONE_K = 1.0e3

	! Composite elastic properties
	REAL(KIND=KIND(0.D0)), PARAMETER, DIMENSION(3) :: E_c = (/1.63586E+05, 9.00065E+03, 9.00065E+03/)
	REAL(KIND=KIND(0.D0)), PARAMETER, DIMENSION(3) :: NU_c = (/3.15770E-01, 3.15770E-01, 4.02876E-01/)
	REAL(KIND=KIND(0.D0)), PARAMETER, DIMENSION(3) :: G_c = (/4.97247E+03, 4.97247E+03, 3.20892E+03/) 	! G_c(3) = E_c(2)/(2*(1+NU_c(3)))

	! Matrix elastic properties
	REAL(KIND=KIND(0.D0)), PARAMETER, DIMENSION(3) :: E_m = (/3.47600E+03, 3.47600E+03, 3.47600E+03/)
	REAL(KIND=KIND(0.D0)), PARAMETER, DIMENSION(3) :: NU_m = (/3.550E-01, 3.550E-01, 3.550E-01/)
    REAL(KIND=KIND(0.D0)), PARAMETER, DIMENSION(3) :: G_m = (/1.28300E+03, 1.28300E+03, 1.28300E+03/)

    ! Material Property Degradation Factors
    REAL(KIND=KIND(0.D0)), PARAMETER :: DEG_FACTOR_MATRIX = 1.0e-1, DEG_FACTOR_DELAM = 1.0e-4, DEG_FACTOR_FIBER = 1.0e-4

	! Based on the elastic properties, stiffness and several other matrices are required to
	! drive the micromechanics code. These matrices were pre-calculated in MATLAB and written
	! in a text file to reduce the computational burden of the UMAT
	!
	! Read in properties from the material file
	REAL(KIND=KIND(0.D0)), DIMENSION(6,6) :: C, S, C_f, S_f, C_m, S_m
	REAL(KIND=KIND(0.D0)), DIMENSION(6,6) :: C_mdeg, S_mdeg, C_f_mdeg, S_f_mdeg, C_m_mdeg_ply, S_m_mdeg_ply, C_m_mdeg_delam, S_m_mdeg_delam
	REAL(KIND=KIND(0.D0)), DIMENSION(6,6) :: C_fdeg, S_fdeg, C_f_fdeg, S_f_fdeg, C_m_fdeg, S_m_fdeg
	REAL(KIND=KIND(0.D0)), DIMENSION(6,6) :: Q_f, Q_m, Q_f_mdeg, Q_m_mdeg, Q_f_fdeg, Q_m_fdeg, I
	REAL(KIND=KIND(0.D0)), DIMENSION(6) :: psi_m, psi_mdeg

	! Parameters A, B required for Weibull distribution of fiber strength
	REAL(KIND=KIND(0.D0)), PARAMETER ::A = 3.9630e+03, B = 28.0751  ! scaled values
	! REAL(KIND=KIND(0.D0)), PARAMETER ::A = 2.973570156341983e+03, B = 0.021065871287020e+03  ! unscaled values
	REAL(KIND=KIND(0.D0)) :: Fiber_Strength

END MODULE

!----------------------------------------------------------------------------------------------------------------------------
! Core of the UMAT
!----------------------------------------------------------------------------------------------------------------------------

SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
	RPL, DDSDDT, DRPLDE, DRPLDT, &
	STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
	NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
	CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, JSTEP, KINC)

	! Import the MODULEs
	USE KTF_CONSTANTS
    USE MATERIAL_PROPS

	! Variable definitions as required by Abaqus
	IMPLICIT NONE			! Each variable must be declared explicitly

	! Variables passed in for information
	INTEGER, INTENT(IN) :: NDI 				! Number of direct stress components at this point
	INTEGER, INTENT(IN) :: NSHR 			! Number of engineering shear stress components at this point
	INTEGER, INTENT(IN) :: NTENS 			! Size of the stress or strain component array (NDI + NSHR)
	INTEGER, INTENT(IN) :: NSTATV 			! Number of solution-dependent state variables that are associated with this material type
	INTEGER, INTENT(IN) :: NPROPS 			! User-defined number of material constants associated with this user material
	INTEGER, INTENT(IN) :: NOEL 			! Element number
	INTEGER, INTENT(IN) :: NPT 				! Integration point number
	INTEGER, INTENT(IN) :: LAYER 			! Layer number (for composite shells and layered solids)
	INTEGER, INTENT(IN) :: KSPT 			! Section point number within the current layer
	INTEGER, INTENT(IN) :: JSTEP(1) 		! Step number
	INTEGER, INTENT(IN) :: KINC 			! Increment number
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: STRAN(NTENS)			! Array containing the total strains at the beginning of the increment
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: DSTRAN(NTENS)			! Array of strain increments
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: TIME(2)				! Value of total time at the beginning of the current increment
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: DTIME					! Time increment
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: TEMP					! Temperature at the start of the increment
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: DTEMP					! Increment of temperature
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: PREDEF(1)				! Array of interpolated values of predefined field variables at this point
																! at the start of the increment, based on the values read in at the nodes
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: DPRED(1)				! Array of increments of predefined field variables
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: PROPS(NPROPS) 			! User-specified array of material constants associated with this user material
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: COORDS(3) 				! An array containing the coordinates of this point
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: DROT(3,3) 				! Rotation increment matrix
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: CELENT 				! Characteristic element length, which is a typical length of a line across an element
																! for a first-order element; it is half of the same typical length for a second-order element
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: DFGRD0(3,3) 			! Array containing the deformation gradient at the beginning of the increment
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: DFGRD1(3,3) 			! Array containing the deformation gradient at the end of the increment
	CHARACTER(LEN=80), INTENT(IN) :: CMNAME						! User-defined material name, left justified

	! Variables to be defined in all situations
	REAL(KIND=KIND(0.D0)), INTENT(OUT) :: STRESS(NTENS) 		! This array is passed in as the stress tensor at the beginning of the increment
																! and must be updated in this routine to be the stress tensor at the end of the increment
	REAL(KIND=KIND(0.D0)), INTENT(OUT) :: STATEV(NSTATV)      	! An array containing the solution-dependent state variables. These are passed in as the values
																! at the beginning of the increment, and must be returned as the values at the end of the increment
	REAL(KIND=KIND(0.D0)), INTENT(OUT) :: DDSDDE(NTENS,NTENS) 	! Jacobian matrix of the constitutive model
	REAL(KIND=KIND(0.D0)), INTENT(OUT) :: SSE 					! Specific elastic strain energy, no effect on the solution, except that they are used for energy output
	REAL(KIND=KIND(0.D0)), INTENT(OUT) :: SPD 					! Specific plastic dissipation, no effect on the solution, except that they are used for energy output
	REAL(KIND=KIND(0.D0)), INTENT(OUT) :: SCD 					! Specific  “creep” dissipation, no effect on the solution, except that they are used for energy output

	! Variables to be defined only in a fully-coupled thermal-stress analysis
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: RPL  					! Volumetric heat generation per unit time at the end of the increment caused by mechanical working of the material
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: DDSDDT(NTENS) 			! Variation of the stress increments with respect to the temperature
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: DRPLDE(NTENS) 			! Variation of RPL with respect to the strain increments
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: DRPLDT 				! Variation of RPL with respect to the temperature

	! Variables that can be updated
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: PNEWDT 				! Ratio of suggested new time increment to the time increment being used (DTIME)

	!******************************************************************************************
	! Internal variables updated as necessary with added features in the UMAT
	REAL(KIND=KIND(0.D0)) :: SMISES
	REAL(KIND=KIND(0.D0)) :: EFF_STRESS, EFF_STRENGTH
	REAL(KIND=KIND(0.D0)) :: sigma_m(6), sigma_f(6), sigma_m_eff
	REAL(KIND=KIND(0.D0)) :: dummy(6,6)

	! Variables used for identifying starting/ending Element number, Integration Pt. number
	INTEGER :: COUNT = 0
	INTEGER :: START_NOEL, START_NPT, MAX_EL, MAX_PT
	LOGICAL :: FIRST_CALL = .TRUE., FIRST_RUN = .TRUE.

	! Variables used for calc. of Miner's rule damage accumulation
	REAL(KIND=KIND(0.D0)) :: S_11_Tension, Ncycles_fiber

	! Text file containing propeties for a particular material
	INTEGER :: FID = 90, ioerr

	!******************************************************************************************
	! Body of the SUBROUTINE UMAT
	!******************************************************************************************

	IF ((NOT(FIRST_CALL)) .AND. (FIRST_RUN)) THEN
		IF ((NOEL .EQ. START_NOEL) .AND. (NPT .EQ. START_NPT)) THEN
			! WRITE(FID,*) '---------FIRST RUN HAS ENDED----------'
			FIRST_RUN = .FALSE.
		ENDIF
	ENDIF

	IF (FIRST_RUN) THEN
		! WRITE(FID,*) 'FIRST_RUN', FIRST_RUN
		! WRITE(FID,*) 'FIRST_CALL', FIRST_CALL
		FIRST_RUN = .FALSE.
		COUNT = COUNT + 1

		IF (FIRST_CALL) THEN
			FIRST_CALL = .FALSE.
			START_NOEL = NOEL
			START_NPT = NPT

			! Determine max element and integration point number
			MAX_EL = NOEL
            MAX_PT = NPT

			! Read off material props from MAT_FILE
			CALL  READ_MATERIAL(MAT_FILE, DEG_FACTOR_MATRIX, DEG_FACTOR_DELAM, DEG_FACTOR_FIBER, C, S, C_f, S_f, C_m, S_m, &
			    C_mdeg, S_mdeg, C_f_mdeg, S_f_mdeg, C_m_mdeg_ply, S_m_mdeg_ply, C_m_mdeg_delam, S_m_mdeg_delam, &
				C_fdeg, S_fdeg, C_f_fdeg, S_f_fdeg, C_m_fdeg, S_m_fdeg, &
				Q_f, Q_m, Q_f_mdeg, Q_m_mdeg, Q_f_fdeg, Q_m_fdeg, psi_m, psi_mdeg, I)

			! Initialize STATEV(1) : damage variable & STATEV(4) : modulus degradation percentage due to damage
			STATEV(1) = 0.0
			STATEV(4) = 1.0

			! Effective Strength
			EFF_STRENGTH = SQRT(ONE / B_s1)

		ELSE
			IF (NOEL .GT. MAX_EL) MAX_EL = NOEL
			IF (NPT .GT. MAX_PT) MAX_PT = NPT
		ENDIF

		! Calculate Stress (undamaged state)
		STRESS = MATMUL(C, STRAN + DSTRAN)
		DDSDDE = C

	ELSE
	! Activate the UMAT at the start of second, i.e. the Fatigue Step
		IF (JSTEP(1) .GT. 1) THEN
			!******************************************************************************************
			! Assign a strength value to each element using Weibull distribution
			IF (KINC .EQ. 1) THEN
				CALL RANDOM_WEIBULL(A, B, Fiber_Strength)
				Fiber_Strength = Fiber_Strength / phi_f
                STATEV(5) = Fiber_Strength

				! ** Debug print **
				! WRITE(FID,*) 'Fiber_Strength: ',Fiber_Strength
			ENDIF
			!******************************************************************************************

			! Matrix damage is < 1.0 -> use unfailed/pristine properties for all constituents
			IF (STATEV(1) .LT. ONE) THEN
				IF (CMNAME .EQ. 'MATRIX') THEN
				! Caculate damage increment for the delamination layer matrix
					! Calculate stress
					STRESS = MATMUL(C_m, STRAN + DSTRAN)
					DDSDDE = C_m

					! Calculate Von Mises stress (i.e. Effective Stress for interlaminar matrix)
					CALL VON_MISES(STRESS, SMISES)
					EFF_STRESS = SMISES
					STATEV(3) = EFF_STRESS

				ELSE
				! Assuming only 2 materials, CMNAME is now the name of the composite material
				! Caculate damage increment for the ply level matrix
					STRESS = MATMUL(C, STRAN + DSTRAN)
					DDSDDE = C

					! Calculate effective matrix stress (i.e. Effective Stress for intralaminar/in-situ matrix)
					sigma_m = MATMUL(Q_m, STRESS) - (TEMP * psi_m)
					CALL EFF_MATRIX_STRESS(B_t, B_s1, B_s2, beta, tau_0, sigma_m, sigma_m_eff)
                    EFF_STRESS = sigma_m_eff

					! Update State Variables
					STATEV(3) = EFF_STRESS
					STATEV(4) = ( EFF_STRESS / ((EFF_STRENGTH - EFF_STRESS) * STATEV(1) + EFF_STRESS) )

					! Safeguard against STATEV(4) becoming so low that a numerical error occurs while calculating it
					IF (STATEV(4) .LT. 1.0E-2) THEN
						STATEV(4) = 1.0E-2
					ENDIF
				ENDIF

				! Calculate Cycles to Failure as well as other STATEVs
				CALL DAMAGE_INC_FATIGUE(EFF_STRESS, TEMP, CMNAME, KINC, NCYCLES_INC, STATEV(1), STATEV(2))

				! Calculate # of cycles to failure from the current state of damage, i.e. STATEV(1)
				CALL NCYCLES_TO_FAILURE_CALC(STRAN, DSTRAN, TEMP, STATEV(1), STATEV(7))
				UNSORTED_NCYCLES_FAIL_MATRIX(NOEL) = STATEV(7)

				! Making sure STATEV(1) does not have a higher value than 1.0
				IF (STATEV(1) .GT. ONE) THEN
					STATEV(1) = ONE
				ENDIF

				!******************************************************************************************
				! Implement fiber failure using Miner's rule for damage accumulation

				! Fiber is not failed yet, calculate damage increment
				IF (STATEV(6) .LT. 1.0) THEN
					sigma_f = MATMUL(Q_f, STRESS) - (TEMP * psi_m)
					Ncycles_fiber = exp((STATEV(5) - sigma_f(1)) / (64.23))
					STATEV(6) = STATEV(6) + (NCYCLES_INC / Ncycles_fiber)

					UNSORTED_NCYCLES_FAIL_FIBER(NOEL) = Ncycles_fiber

					! Set max value of NCYCLES_FAIL as MAX_NCYCLES_FAIL
					IF (UNSORTED_NCYCLES_FAIL_FIBER(NOEL) .GT. MAX_NCYCLES_FAIL_FIBER) THEN
						UNSORTED_NCYCLES_FAIL_FIBER(NOEL) = MAX_NCYCLES_FAIL_FIBER
					ENDIF

					! Set UNSORTED_NCYCLES_FAIL_FIBER(NOEL) as -1.D0 (required for INDEX_FINDER logic)
					! when it's < 1000.0 and NaN
					IF (UNSORTED_NCYCLES_FAIL_FIBER(NOEL) .LT. ONE_K &
					    .OR. UNSORTED_NCYCLES_FAIL_FIBER(NOEL) .NE. UNSORTED_NCYCLES_FAIL_FIBER(NOEL)) THEN
						UNSORTED_NCYCLES_FAIL_FIBER(NOEL) = -ONE
					ENDIF

					IF (STATEV(6) .GT. 0.5) THEN
						STATEV(6) = 1.0
					ENDIF

				! Use fiber degraded props
				ELSE
					DDSDDE = C_fdeg
					STRESS = MATMUL(C_fdeg, STRAN + DSTRAN)

					! Effective Matrix Stress
					sigma_m = MATMUL(Q_m_fdeg, STRESS) - (TEMP * psi_mdeg)
					CALL EFF_MATRIX_STRESS(B_t_fdeg, B_s1_fdeg, B_s2_fdeg, beta, tau_0, sigma_m, sigma_m_eff)
                    EFF_STRESS = sigma_m_eff

                    UNSORTED_NCYCLES_FAIL_FIBER(NOEL) = -ONE

					! Update State Variables
					STATEV(1) = 1.0
					STATEV(2) = STATEV(2)
					STATEV(3) = EFF_STRESS
					STATEV(4) = STATEV(4)
					STATEV(5) = STATEV(5)
					STATEV(6) = 1.0
				ENDIF
				!******************************************************************************************

			! STATEV(1) >= 1.0 -> Degrade matrix properties
			ELSE
				IF (CMNAME .EQ. 'MATRIX') THEN
					! Use Matrix Degraded Props (Delamination layer)
					DDSDDE = C_m_mdeg_delam
                    STRESS = MATMUL(C_m_mdeg_delam, STRAN + DSTRAN)

					! Calculate the Von Mises Stress (i.e. Effective Stress for interlaminar matrix)
					CALL VON_MISES(STRESS, SMISES)
                    EFF_STRESS = SMISES

					! Update State Variables
					STATEV(1) = 1.0
					STATEV(2) = STATEV(2)
					STATEV(3) = EFF_STRESS

				! Assuming only 2 materials, CMNAME is now the composite material
				ELSE
					! Use Matrix Degraded Props (Ply/in-situ matrix)
					DDSDDE = C_mdeg
					STRESS = MATMUL(C_mdeg, STRAN + DSTRAN)

					! Effective Matrix Stress (i.e. Effective Stress for intralaminar/in-situ matrix)
					sigma_m = MATMUL(Q_m_mdeg, STRESS) - (TEMP * psi_mdeg)
					CALL EFF_MATRIX_STRESS(B_t_mdeg, B_s1_mdeg, B_s2_mdeg, beta, tau_0, sigma_m, sigma_m_eff)
					EFF_STRESS = sigma_m_eff

					! Update State Variables
					STATEV(1) = 1.0
					STATEV(2) = STATEV(2) + NCYCLES_INC
					STATEV(3) = EFF_STRESS
					STATEV(4) = STATEV(4)
					STATEV(5) = STATEV(5)
					STATEV(6) = STATEV(6)
					STATEV(7) = 0.0

					! Set "number of cycles to failure" as -1.0 (since matrix already failed) to remove this
					! integration pt/element from INDEX_FINDER_MATRIX calculations
					UNSORTED_NCYCLES_FAIL_MATRIX(NOEL) = - ONE

					!******************************************************************************************
					! Fiber stress from composite stress

					! STATEV(6) < 1.0 -> Use undegraded fiber properties
					IF (STATEV(6) .LT. 1.0) THEN
					! Fiber is not failed yet, calculate damage increment, i.e. STATEV(6)
						sigma_f = MATMUL(Q_f_mdeg, STRESS) - (TEMP * psi_mdeg)
						Ncycles_fiber = exp((STATEV(5) - sigma_f(1)) / (64.23))
						STATEV(5) = STATEV(5)
                        STATEV(6) = STATEV(6) + (NCYCLES_INC / Ncycles_fiber)

						UNSORTED_NCYCLES_FAIL_FIBER(NOEL) = Ncycles_fiber

						! Set max value of NCYCLES_FAIL as MAX_NCYCLES_FAIL
						IF (UNSORTED_NCYCLES_FAIL_FIBER(NOEL) .GT. MAX_NCYCLES_FAIL_FIBER) THEN
							UNSORTED_NCYCLES_FAIL_FIBER(NOEL) = MAX_NCYCLES_FAIL_FIBER
						ENDIF

						! Set UNSORTED_NCYCLES_FAIL_FIBER(NOEL) as -1.D0 (required for INDEX_FINDER logic)
						! when it's < 1000.0 and NaN
						IF (UNSORTED_NCYCLES_FAIL_FIBER(NOEL) .LT. ONE_K &
							.OR. UNSORTED_NCYCLES_FAIL_FIBER(NOEL) .NE. UNSORTED_NCYCLES_FAIL_FIBER(NOEL)) THEN
							UNSORTED_NCYCLES_FAIL_FIBER(NOEL) = -ONE
						ENDIF

						IF (STATEV(6) .GE. 0.5) THEN
							STATEV(6) = 1.0
						ENDIF

					! STATEV(6) >= 1.0 -> Use degraded fiber properties
					ELSE
						DDSDDE = C_fdeg
						STRESS = MATMUL(C_fdeg, STRAN + DSTRAN)

						sigma_m = MATMUL(Q_m_fdeg, STRESS) - (TEMP * psi_mdeg)  ! Effective Matrix Stress
						CALL EFF_MATRIX_STRESS(B_t_fdeg, B_s1_fdeg, B_s2_fdeg, beta, tau_0, sigma_m, sigma_m_eff)
                        EFF_STRESS = sigma_m_eff  ! Define Effective Stress

                        UNSORTED_NCYCLES_FAIL_FIBER(NOEL) = -ONE

						! Update State Variables
						STATEV(1) = 1.0
						STATEV(2) = STATEV(2)
						STATEV(3) = EFF_STRESS
						STATEV(4) = STATEV(4)
						STATEV(5) = STATEV(5)
						STATEV(6) = 1.0
						STATEV(7) = 0.0
					ENDIF
					!******************************************************************************************
				ENDIF
			ENDIF

		ELSE
		! Calculations for the first analysis step (Static Load)
			IF (CMNAME .EQ. 'MATRIX') THEN
			! Caculate stress for the delamination layer (isotropic matrix)
				STRESS = MATMUL(C_m, STRAN + DSTRAN)
				DDSDDE = C_m

			ELSE
			! Caculate stress for the composite ply
				STRESS = MATMUL(C, STRAN + DSTRAN)
				DDSDDE = C
			ENDIF
			STATEV = 0.D0	! State Variables are set to zero for the static load step
		ENDIF
	ENDIF

	! Values that are always defined as zero for this UMAT
	SSE = 0.D0
	SPD = 0.D0
	SCD = 0.D0

    ! Stopping criteria I: End analysis if 2 million load cycles have already been reached
	IF (STATEV(2) .GT. 2.0e+6) THEN
		CALL XIT
	ENDIF

	RETURN

END SUBROUTINE UMAT

!-----------------------------------------------------------------------------------------------------------------------------
! User subroutine UEXTERNALDB
!-----------------------------------------------------------------------------------------------------------------------------
! UEXTERNALDB is called once each at the beginning of the analysis, at the beginning of each increment,
! at the end of each increment, and at the end of the analysis.
! More at: https://help.3ds.com/2018/english/dssimulia_established/SIMACAESUBRefMap/simasub-c-uexternaldb.htm?ContextScope=all#simasub-c-uexternaldb

SUBROUTINE UEXTERNALDB(LOP, LRESTART, TIME, DTIME, KSTEP, KINC, NSTATV, STATEV)

	! Import the MODULEs
	USE KTF_CONSTANTS
    USE MATERIAL_PROPS

	IMPLICIT NONE			! Each variable must be declared explicitly

	! Variables passed in for information
	INTEGER, INTENT(IN) :: NSTATV
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: STATEV(NSTATV)
	! LOP=0 indicates that the user subroutine is being called at the start of the analysis
	! LOP=1 indicates that the user subroutine is being called at the start of the current analysis increment
	! LOP=2 indicates that the user subroutine is being called at the end of the current analysis increment
	! LOP=3 indicates that the user subroutine is being called at the end of the analysis
	INTEGER, INTENT(IN) :: LOP
	INTEGER, INTENT(IN) :: LRESTART
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: TIME(2)				! TIME(1): Value of current step time
																! TIME(2): Value of total time at the beginning of the current increment
	REAL(KIND=KIND(0.D0)), INTENT(IN) :: DTIME					! Time increment
	INTEGER, INTENT(IN) :: KSTEP 								! Step number
	INTEGER, INTENT(IN) :: KINC							        ! Increment number

	! Define variables for internal use
	INTEGER :: FID = 101, ioerr

	!******************************************************************************************
	! Body of the SUBROUTINE UEXTERNALDB
	!******************************************************************************************

	IF (LOP .EQ. 0) THEN
		! At the beginning of analysis, create a log file to write debug data
		OPEN(UNIT = FID, FILE = LOG_FILE, STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ioerr)

	! At the start of each increment, calculate NCYCLES_INC
    ELSEIF (LOP .EQ. 1 .AND. KSTEP .GT. 1 .AND. KINC .GT. 1) THEN
		!******************************************************************************************
        ! Implement adaptive cycle selection for cycle jump based on lifetime of lowest x% elements

        ! Stopping criteria II: End analysis if fiber failure has occured in a pre-defined max percentage of elements
        IF (INDEX_FINDER_FIBER(1) .GT. (NELEM_LE + NELEM_UMAT_DELAM + 0.50 * NELEM_UMAT)) THEN
            CALL XIT
        ENDIF

		! Use a sorting algorithm to find the required minimum values of cycles for matrix & fiber failure
        SORTED_NCYCLES_FAIL_MATRIX = UNSORTED_NCYCLES_FAIL_MATRIX
        SORTED_NCYCLES_FAIL_FIBER = UNSORTED_NCYCLES_FAIL_FIBER

        CALL QSORT(SORTED_NCYCLES_FAIL_MATRIX, LEN, ISIZE, COMPAR)
        CALL QSORT(SORTED_NCYCLES_FAIL_FIBER, LEN, ISIZE, COMPAR)

		! Select the value for cycle jump based on lifetime of lowest x% elements
        MIN_NCYCLES_FAIL_MATRIX = MINVAL(SORTED_NCYCLES_FAIL_MATRIX)
        MIN_NCYCLES_FAIL_FIBER = MINVAL(SORTED_NCYCLES_FAIL_FIBER)

        ! Find the index of the (first positive value)/(#cycles of first unfailed element) in the array
		! "SORTED_NCYCLES_FAIL_MATRIX"
		! FINDLOC: Finds the location of a specified value in an array
		! FLOOR(X): Returns the greatest integer less than or equal to X
		! SIGN: Returns the absolute value of the first argument times the sign of the second argument
		INDEX_FINDER_MATRIX = FINDLOC(FLOOR(SIGN(ONE, SORTED_NCYCLES_FAIL_MATRIX)), VALUE=1)
        ADAPTIVE_NCYCLES_MATRIX = SORTED_NCYCLES_FAIL_MATRIX(INDEX_FINDER_MATRIX(1) + FAIL_PCT_MATRIX(1))

        INDEX_FINDER_FIBER = FINDLOC(FLOOR(SIGN(ONE, SORTED_NCYCLES_FAIL_FIBER)), VALUE=1)
        ADAPTIVE_NCYCLES_FIBER = SORTED_NCYCLES_FAIL_FIBER(INDEX_FINDER_FIBER(1) + FAIL_PCT_FIBER(1))

        ! Implement situation-based choice for FAIL_PCT_MATRIX/cycle jumping, anecdotal information is required for this
        ! to speed up the simulations
        IF (INDEX_FINDER_MATRIX(1) .LT. (NELEM_LE + 0.10 * NELEM_UMAT)) THEN
            FAIL_PCT_MATRIX = 3 * BASE_FAILURE_PCT

            IF (ADAPTIVE_NCYCLES_MATRIX .LT. MIN_CYCLE_JUMP) THEN
                ADAPTIVE_NCYCLES_MATRIX = MIN_CYCLE_JUMP
            ENDIF

        ELSEIF (INDEX_FINDER_MATRIX(1) .GT. (NELEM_LE + 0.10 * NELEM_UMAT) .AND. &
            INDEX_FINDER_MATRIX(1) .LT. (NELEM_LE + 0.30 * NELEM_UMAT)) THEN
            FAIL_PCT_MATRIX = 2 * BASE_FAILURE_PCT

            IF (ADAPTIVE_NCYCLES_MATRIX .LT. MIN_CYCLE_JUMP / 2.0) THEN
                ADAPTIVE_NCYCLES_MATRIX = MIN_CYCLE_JUMP / 1.0   ! 2.0
            ENDIF

        ELSEIF (INDEX_FINDER_MATRIX(1) .GT. (NELEM_LE + 0.30 * NELEM_UMAT) .AND. &
            INDEX_FINDER_MATRIX(1) .LT. (NELEM_LE + 0.50 * NELEM_UMAT)) THEN
            FAIL_PCT_MATRIX = 1 * BASE_FAILURE_PCT

            IF (ADAPTIVE_NCYCLES_MATRIX .LT. MIN_CYCLE_JUMP / 4.0) THEN
                ADAPTIVE_NCYCLES_MATRIX = MIN_CYCLE_JUMP / 1.0    ! 4.0
            ENDIF

        ELSE
            FAIL_PCT_MATRIX = 0.5 * BASE_FAILURE_PCT

            IF (ADAPTIVE_NCYCLES_MATRIX .LT. MIN_CYCLE_JUMP / 10.0) THEN
                ADAPTIVE_NCYCLES_MATRIX = MIN_CYCLE_JUMP / 1.0   ! 10.0
            ENDIF

        ENDIF

        ! Compare the adaptive cycle jumping between matrix and fiber constituents and choose the smaller value
        IF (ADAPTIVE_NCYCLES_MATRIX .LT. ADAPTIVE_NCYCLES_FIBER) THEN
            ADAPTIVE_NCYCLES = ADAPTIVE_NCYCLES_MATRIX
        ELSE
            ADAPTIVE_NCYCLES = ADAPTIVE_NCYCLES_FIBER
        ENDIF

		! Set lowest value of ADAPTIVE_NCYCLES_MATRIX/NCYCLES_INC as MIN_CYCLE_JUMP / 10.0
		IF (ADAPTIVE_NCYCLES .LT. MIN_CYCLE_JUMP / 10.0) THEN
			ADAPTIVE_NCYCLES = MIN_CYCLE_JUMP / 10.0
        ENDIF

        ! Update number of cycles to jump
		NCYCLES_INC = ADAPTIVE_NCYCLES

		! ** Debug print: Write the sorted cycles to failure data into a file for the requested increments **
		! IF (KINC .EQ. 2 .OR. KINC .EQ. 3) THEN
        !     WRITE(101,*) 'UNSORTED_NCYCLES_FAIL_FIBER: ', UNSORTED_NCYCLES_FAIL_FIBER
        !     WRITE(101,*) 'SORTED_NCYCLES_FAIL_FIBER: ', SORTED_NCYCLES_FAIL_FIBER
		! ENDIF

        ! ** Debug print: At all increments
        WRITE(101,*)'KINC: ', KINC
        WRITE(101,*) 'MIN_NCYCLES_FAIL_MATRIX: ', MIN_NCYCLES_FAIL_MATRIX
        WRITE(101,*) 'FAIL_PCT_MATRIX: ', FAIL_PCT_MATRIX
        WRITE(101,*) 'INDEX_FINDER_MATRIX: ', INDEX_FINDER_MATRIX
        WRITE(101,*) 'INDEX_FINDER_FIBER: ', INDEX_FINDER_FIBER
        WRITE(101,*) 'ADAPTIVE_NCYCLES_MATRIX: ', ADAPTIVE_NCYCLES_MATRIX
        WRITE(101,*) 'ADAPTIVE_NCYCLES_FIBER: ', ADAPTIVE_NCYCLES_FIBER
		WRITE(101,*) 'NCYCLES_INC: ', NCYCLES_INC
		WRITE(101,*) '*---------------------------------------------'
		!******************************************************************************************

	! At the end of analysis, close the log file
	ELSEIF (LOP .EQ. 3) THEN
	    CLOSE(101)
	ENDIF

	RETURN

END SUBROUTINE UEXTERNALDB

!-----------------------------------------------------------------------------------------------------------------------------
! External routine COMPAR for the QSORT algorithm to define ordering of elements
INTEGER*2 FUNCTION COMPAR( A, B )
	REAL(KIND=KIND(0.D0)) :: A, B

	IF ( A .LT. B ) COMPAR = -1
	IF ( A .EQ. B ) COMPAR = 0
	IF ( A .GT. B ) COMPAR = 1

	RETURN
END
!-----------------------------------------------------------------------------------------------------------------------------

! Include called functions
INCLUDE 'DAMAGE_INC_FATIGUE.f90'
INCLUDE 'EFF_MATRIX_STRESS.f90'
INCLUDE 'MATRIX_OPS.f90'
INCLUDE 'NCYCLES_TO_FAILURE_CALC.f90'
INCLUDE 'RANDOM_WEIBULL_DST.f90'
INCLUDE 'READ_MATFILE.f90'
INCLUDE 'VON_MISES.f90'

