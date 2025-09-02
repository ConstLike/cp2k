!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Davidson iterative solver for TDDFT excited states.
!>
!> \details This module implements a theory-agnostic Davidson iterative diagonalization
!>          algorithm for computing electronic excited states within Time-Dependent Density
!>          Functional Theory (TDDFT). The implementation supports multiple excitation
!>          theories through a polymorphic interface:
!>          - Spin-flip (SF): Excited states with spin contaminations
!>          - Mixed Reference Spin-Flip (MRSF): Excited states with pure separation
!>            to singlets or triplets
!>
!> \par Algorithm Overview
!>      The Davidson method iteratively builds a Krylov subspace by applying the TDDFT
!>      Hamiltonian (orbital Hessian) to trial vectors and diagonalizing the resulting
!>      subspace matrix. The algorithm follows these key steps:
!>      1. Initialize trial vectors using orbital energy differences
!>      2. Apply TDDFT operators (Hartree-Fock exchange, Coulomb, XC kernel)
!>      3. Build and diagonalize subspace Hamiltonian
!>      4. Form Ritz vectors and compute residuals
!>      5. Check convergence or expand subspace
!>      6. Restart if necessary to maintain efficiency
!>
!> \par Theory Abstraction
!>      The implementation uses object-oriented programming to separate theory-specific
!>      logic from the Davidson algorithm. Each theory implements the tddft_theory_type
!>      interface, providing methods for:
!>      - Initial guess generation
!>      - MO-to-AO transformations
!>      - Operator applications
!>      - AO-to-MO back transformations
!>      - Results analysis and printing
!>
!> \author  Konstantin Komarov (constlike@gmail.com)
!> \date    Aug 2025
! **************************************************************************************************
MODULE qs_tddft_methods
   USE bibliography,                    ONLY: Iannuzzi2005, &
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type, &
                                              tddfpt2_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add, &
                                              cp_fm_trace, &
                                              cp_fm_column_scale
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_types,                     ONLY: cp_fm_set_all, &
                                              cp_fm_to_fm, &
                                              cp_fm_set_element, &
                                              cp_fm_get_element, &
                                              cp_fm_vectorsnorm
   USE cp_log_handling,                 ONLY: cp_get_default_logger, &
                                              cp_logger_get_default_io_unit, &
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_add_iter_level, &
                                              cp_iterate, &
                                              cp_print_key_finished_output, &
                                              cp_print_key_unit_nr, &
                                              cp_rm_iter_level
   USE header,                          ONLY: tddft_header
   USE input_constants,                 ONLY: tddfpt_kernel_full
   USE input_section_types,             ONLY: section_vals_get_subs_vals, &
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush, &
                                              m_walltime
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env, &
                                              qs_environment_type
   USE qs_tddfpt2_methods,              ONLY: kernel_info, &
                                              tddfpt_input
   USE qs_tddfpt2_properties,           ONLY: tddfpt_dipole_operator
   USE qs_tddfpt2_smearing_methods,     ONLY: tddfpt_smeared_occupation
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE qs_tddfpt2_utils,                ONLY: tddfpt_init_mos, &
                                              tddfpt_oecorr, &
                                              tddfpt_release_ground_state_mos
   USE qs_tddft_theory,                 ONLY: tddft_theory_type
   USE qs_tddft_theory_mrsf,            ONLY: tddft_mrsf_type
   USE qs_tddft_theory_sf,              ONLY: tddft_sf_type
   USE qs_tddft_types,                  ONLY: tddft_info_type, &
                                              tddft_initialization, &
                                              tddft_finalization
   USE qs_tddft_utils,                  ONLY: allocate_fm_matrix, &
                                              deallocate_fm_matrix
   USE string_utilities,                ONLY: integer_to_string
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_methods'

   PUBLIC :: tddft

CONTAINS

! **************************************************************************************************
!> \brief Main entry point for TDDFT excited state calculations
!>
!> \details This subroutine orchestrates the complete TDDFT calculation workflow:
!>          - Initializes ground-state molecular orbitals
!>          - Applies orbital energy corrections if needed
!>          - Executes Davidson iterations for multiple theories
!>          - Manages theory selection and iteration over SF, MRSF-S, and MRSF-T
!>
!> \param[inout] qs_env       QuickStep environment containing system information
!> \param[in]    calc_forces  Logical flag for force calculation (currently unused)
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE tddft(qs_env, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft'

      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      TYPE(section_vals_type), POINTER                   :: tddfpt_section
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(tddft_info_type)                              :: tddft_env
      CLASS(tddft_theory_type), ALLOCATABLE              :: theory
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_oep
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      CHARACTER(LEN=6), DIMENSION(3)                     :: theory_to_run
      CHARACTER(LEN=6)                                   :: logo
      INTEGER                                            :: handle, i, log_unit
      LOGICAL                                            :: do_admm, do_exck, do_hfx, do_hfxlr, &
                                                            do_hfxsr
      TYPE(section_vals_type), POINTER                   :: hfxsr_section, lri_section

      CALL timeset(routineN, handle)

      ! Forces are not yet implemented for TDDFT
      IF (calc_forces) THEN
         CPWARN("Force calculation not yet implemented for TDDFT")
      END IF

      ! Initialize logging environment
      NULLIFY (logger)
      logger => cp_get_default_logger()

      ! Get TDDFPT input section
      NULLIFY (tddfpt_section)
      tddfpt_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT")

      ! Parse input and determine calculation settings
      CALL tddfpt_input(qs_env, do_hfx, do_admm, do_exck, &
                        do_hfxsr, do_hfxlr, xc_section, tddfpt_print_section, &
                        lri_section, hfxsr_section)

      ! Open output unit for program banner
      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "PROGRAM_BANNER", &
                                      extension=".tddfptLog")

      ! Extract essential environment components
      CALL get_qs_env(qs_env, &
                      blacs_env=blacs_env, &
                      dft_control=dft_control, &
                      matrix_ks=matrix_ks)

      ! Configure TDDFPT control parameters
      NULLIFY(tddfpt_control)
      tddfpt_control => dft_control%tddfpt2_control
      tddfpt_control%do_hfx = do_hfx
      tddfpt_control%do_admm = do_admm
      tddfpt_control%do_hfxsr = do_hfxsr
      tddfpt_control%hfxsr_primbas = 0
      tddfpt_control%hfxsr_re_int = .TRUE.
      tddfpt_control%do_hfxlr = do_hfxlr
      tddfpt_control%do_exck = do_exck

      ! Cite relevant reference
      ! TODO: Add MRSF-related citations
      CALL cite_reference(Iannuzzi2005)

      ! Print calculation header and kernel information
      CALL tddft_header(log_unit)
      CALL kernel_info(log_unit, dft_control, tddfpt_control, xc_section)

      ! Initialize ground-state molecular orbitals
      NULLIFY (gs_mos)
      CALL tddfpt_init_mos(qs_env, gs_mos, log_unit)

      ! Apply smeared occupation numbers if requested
      IF (tddfpt_control%do_smearing) THEN
         CALL tddfpt_smeared_occupation(qs_env, gs_mos, log_unit)
      END IF

      ! Apply orbital energy corrections to obtain corrected KS-matrix
      CALL tddfpt_oecorr(qs_env, gs_mos, matrix_ks_oep)
      IF (ASSOCIATED(matrix_ks_oep)) matrix_ks => matrix_ks_oep

      ! Check compatibility of LRI approximation with kernel choice
      IF ((tddfpt_control%do_lrigpw) .AND. &
          (tddfpt_control%kernel /= tddfpt_kernel_full)) THEN
         CALL cp_abort(__LOCATION__, "LRI only implemented for full kernel")
      END IF

      ! Define theories to execute
      ! TODO: Make this configurable via input
      theory_to_run = ["    SF", "MRSF-S", "MRSF-T"]

      ! Main theory loop - iterate over all requested theories
      DO i = 1, SIZE(theory_to_run)

         ! Initialize TDDFT environment for current theory
         CALL tddft_initialization(qs_env, tddft_env, gs_mos)

         logo = theory_to_run(i)

         ! Allocate theory-specific handler based on theory type
         SELECT CASE (logo)
            CASE('    SF')
               tddft_env%response_type = 1 ! Spin-flip TDDFT
               ALLOCATE(tddft_sf_type :: theory)
            CASE('MRSF-S')
               tddft_env%response_type = 2 ! Mixed reference spin-flip singlets
               ALLOCATE(tddft_mrsf_type :: theory)
            CASE('MRSF-T')
               tddft_env%response_type = 3 ! Mixed reference spin-flip triplets
               ALLOCATE(tddft_mrsf_type :: theory)
          END SELECT

         ! Execute Davidson solver for current theory
         CALL davidson_solver(tddft_env, theory, qs_env, tddfpt_control, logger, &
                              tddfpt_print_section, gs_mos, &
                              blacs_env, logo)

         ! Clean up theory-specific resources
         DEALLOCATE(theory)

         ! Finalize TDDFT environment
         CALL tddft_finalization(tddft_env)

      END DO

      ! Clean up ground-state molecular orbitals
      DO i = tddft_env%nspins, 1, -1
         CALL tddfpt_release_ground_state_mos(gs_mos(i))
      END DO
      DEALLOCATE (gs_mos)

      CALL timestop(handle)

   END SUBROUTINE tddft

! **************************************************************************************************
!> \brief Davidson iterative solver for TDDFT eigenvalue problems
!>
!> \details Implements the Davidson diagonalization algorithm with the following features:
!>          - Iterative subspace expansion using Krylov vectors
!>          - Automatic restart when subspace becomes too large or stagnates
!>          - Convergence monitoring based on residual norms and eigenvalue changes
!>          - Preconditioned residual vectors for improved convergence
!>
!> \param[inout] tddft_env         TDDFT environment containing Davidson data structures
!> \param[inout] theory_env        Theory-specific handler for operator applications
!> \param        qs_env            QuickStep environment
!> \param[in]    tddfpt_control    TDDFPT control parameters
!> \param[in]    logger            Output logger
!> \param[in]    tddfpt_print_section  Print control section
!> \param[in]    gs_mos            Ground-state molecular orbitals
!> \param[in]    blacs_env         BLACS parallel environment
!> \param[in]    logo              Theory identifier for output
!>
!> \par Algorithm Details
!>      The Davidson method builds a subspace of dimension nVec by applying the TDDFT
!>      Hamiltonian to trial vectors. The subspace Hamiltonian H_ij = <X_i|Sigma_j>
!>      is diagonalized to obtain Ritz values (approximate eigenvalues) and Ritz vectors.
!>      Convergence is assessed by monitoring residual norms ||Sigma_k - omega_k X_k||.
!>
!> \par Restart Strategy
!>      Restarts are triggered when:
!>      - Subspace dimension exceeds maximum (mxvec)
!>      - Convergence stagnates (stall_count > 3)
!>      The restart procedure retains the best Ritz vectors as new trial vectors.
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE davidson_solver(tddft_env, theory_env, qs_env, tddfpt_control, logger, &
                              tddfpt_print_section, gs_mos, blacs_env, logo)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      CLASS(tddft_theory_type), INTENT(INOUT)            :: theory_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt2_control_type), POINTER, INTENT(IN)    :: tddfpt_control
      TYPE(cp_logger_type), POINTER, INTENT(IN)          :: logger
      TYPE(section_vals_type), POINTER, INTENT(IN)       :: tddfpt_print_section
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER, INTENT(IN)                             :: gs_mos
      TYPE(cp_blacs_env_type), POINTER, INTENT(IN)       :: blacs_env
      CHARACTER(LEN=6), INTENT(IN)                       :: logo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'davidson_solver'

      CHARACTER(len=20)                                  :: nStates_str
      INTEGER                                            :: handle, handle2, log_unit, &
                                                            iter, iend, ist, jvec
      REAL(kind=dp)                                      :: t1, t2
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: ostrength

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      ! Initialize theory-specific work arrays and parameters
      CALL theory_env%allocate(qs_env, tddft_env)

      ! Only proceed if iterations are requested
      IF (tddfpt_control%niters > 0) THEN

         ! Open logger
         log_unit = cp_print_key_unit_nr(logger, &
            tddfpt_print_section, "GUESS_VECTORS", extension=".tddfptLog")

         ! Initialize trial subspace with orthonormalized guess vectors
         CALL theory_env%init_guess(qs_env, tddft_env, gs_mos, log_unit)

         ! Close logger
         CALL cp_print_key_finished_output(log_unit, logger, tddfpt_print_section, "GUESS_VECTORS")

         ! Open iteration info logger
         log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "ITERATION_INFO", extension=".tddfptLog")
         IF (log_unit > 0) THEN
            WRITE (log_unit, '(/,1x,79("-"),/," -",17(" "),A,22(" "),"-",/,1x,79("-"))') &
               TRIM(logo)//"-TDDFT WAVEFUNCTION OPTIMIZATION"
            CALL m_flush(log_unit)
         END IF

         CALL cp_add_iter_level(logger%iter_info, "TDDFT_SCF")

         ! Initialize convergence monitoring variables
         tddft_env%max_error_prev = 10.0_dp
         tddft_env%max_delta_eval_prev = 10.0_dp

         ! Outer loop for managing restarts
         davidson_outer_loop: DO
            ist = 1
            iend = tddft_env%nVec
            tddft_env%status(1:2) = .FALSE.
            tddft_env%stall_count = 0

            ! Inner Davidson iterations
            davidson_loop: DO

               CALL timeset(routineN//"_inner_loop", handle2)

               CALL cp_iterate(logger%iter_info, iter_nr_out=iter)

               ! Apply orbital Hessian to new trial vectors
               ! This is the core TDDFT operation: Sigma_j = H * X_j
               DO jvec = ist, iend
                  ! Transform trial vector from MO to AO representation
                  CALL theory_env%mo_to_ao(tddft_env, gs_mos, jvec)

                  ! Apply TDDFT operators (HFX, Coulomb, XC kernel) in AO basis
                  CALL theory_env%apply(qs_env, tddft_env, tddfpt_control)

                  ! Transform resulting sigma vector back to MO representation
                  CALL theory_env%ao_to_mo(tddft_env, gs_mos, jvec)
               END DO

               ! Allocate subspace matrices for current iteration
               CALL allocate_fm_matrix(tddft_env%subspace, tddft_env%nVec, tddft_env%nVec, blacs_env)
               CALL allocate_fm_matrix(tddft_env%subspace_evecs, tddft_env%nVec, tddft_env%nVec, blacs_env)

               ! Build subspace Hamiltonian: H_ij = <X_i|Sigma_j>
               CALL build_subspace(tddft_env)

               ! Diagonalize subspace: H*C = Omega*C
               CALL diagonalize_subspace(tddft_env)

               ! Form Ritz vectors: Z_k = sum_j c_jk X_j
               CALL form_ritz_vectors(tddft_env)

               ! Compute residuals: r_k = Sigma_k - Omega_k Z_k
               CALL compute_residuals(tddft_env)

               ! Clean up subspace matrices
               CALL deallocate_fm_matrix(tddft_env%subspace)
               CALL deallocate_fm_matrix(tddft_env%subspace_evecs)

               ! Log current status of Ritz vectors
               CALL log_ritz_status(tddft_env, log_unit)

               ! Check convergence criteria
               tddft_env%status(1) = is_converged(tddft_env, tddfpt_control)

               IF (tddft_env%status(1)) THEN
                  ! Converged: commit Ritz solution
                  CALL commit_ritz_solution(tddft_env)
               ELSE
                  ! Not converged: continue iterations or restart

                  ! Detect convergence stagnation
                  CALL detect_stall(tddft_env, iter)

                  ! Check if restart is needed
                  tddft_env%status(2) = needs_restart(tddft_env)

                  IF (.NOT. tddft_env%status(2)) THEN
                     ! Continue: precondition and expand subspace

                     ! Apply diagonal preconditioning to residuals
                     CALL apply_preconditioner(tddft_env)

                     ! Orthogonalize residuals and expand Krylov space
                     CALL expand_krylov_space(tddft_env)
                  END IF

                  IF (tddft_env%status(2)) THEN
                     ! Restart: use Ritz vectors as new trial vectors
                     CALL restart_davidson(tddft_env)
                  END IF
               END IF

               ! Print iteration timing and convergence status
               t2 = m_walltime()
               IF (log_unit > 0) THEN
                  WRITE (log_unit, '(T7,I8,T24,F7.1,T40,ES11.4,T66)') iter, t2 - t1, tddft_env%max_error
                  CALL m_flush(log_unit)
               END IF
               t1 = t2

               CALL timestop(handle2)

               ! Exit conditions: converged, max iterations, or restart needed
               IF (tddft_env%status(1) .OR. tddft_env%status(2)) THEN
                  EXIT davidson_loop
               END IF

               ! Update vector indices for next iteration
               ist = tddft_env%novec+1
               iend = tddft_env%nVec

            END DO davidson_loop

            ! Check outer loop exit conditions
            CALL cp_iterate(logger%iter_info, increment=0, iter_nr_out=iter)

            ! Exit if converged or maximum iterations reached
            IF (tddft_env%status(1) .OR. iter >= tddfpt_control%niters) THEN
               EXIT davidson_outer_loop
            END IF

            ! Print restart message
            IF (log_unit > 0) THEN
               WRITE (log_unit, '(1X,25("-"),1X,A,1X,25("-"))') "Restart Davidson iterations"
               CALL m_flush(log_unit)
            END IF

         END DO davidson_outer_loop

         CALL cp_rm_iter_level(logger%iter_info, "TDDFT_SCF")

         ! Print final convergence status
         IF (log_unit > 0) THEN
            CALL integer_to_string(iter, nStates_str)
            IF (tddft_env%status(1)) THEN
               WRITE (log_unit, '(/,1x,79("-"),/," - ",A,/,1x,79("-"))') &
                  TRIM(logo)//"-TDDFT run converged in "//TRIM(nStates_str)//" iteration(s)"
            ELSE
               WRITE (log_unit, '(/,1x,79("-"),/," - ",A,/,1x,79("-"))') &
                  TRIM(logo)//"-TDDFT run did NOT converge after "//TRIM(nStates_str)//" iteration(s)"
            END IF
         END IF

         ! Close logger
         CALL cp_print_key_finished_output(log_unit, logger, tddfpt_print_section, "ITERATION_INFO")
      ELSE
         CALL cp_warn(__LOCATION__, "Skipping TDDFPT wavefunction optimization")
      END IF

      ! Print summary information and analyze results
      log_unit = cp_logger_get_default_io_unit(logger)
      ALLOCATE (ostrength(tddft_env%nstates), source=0.0_dp)

      ! Compute transition dipole operators
      ! TODO: Need to validate it
      CALL tddfpt_dipole_operator(tddft_env%dipole_op_mos_occ, &
                                  tddfpt_control, &
                                  gs_mos, &
                                  qs_env)

      ! Theory-specific analysis and printing
      CALL theory_env%print(qs_env, tddft_env, ostrength, gs_mos, tddfpt_control, log_unit)

      ! Clean up theory-specific arrays
      CALL theory_env%deallocate(qs_env)

      DEALLOCATE (ostrength)

      CALL timestop(handle)

   END SUBROUTINE davidson_solver

! **************************************************************************************************
!> \brief Constructs the subspace Hamiltonian matrix for Davidson diagonalization
!>
!> \details Builds the projected Hamiltonian matrix H_ij = <X_i|Sigma_j> where X_i are
!>          trial vectors and Sigma_j = H*X_j are their images under the TDDFT Hamiltonian.
!>          The matrix is Hermitian due to the self-adjoint nature of the TDDFT operator,
!>          so only the lower triangle is computed and mirrored to the upper triangle.
!>
!> \param[inout] tddft_env  TDDFT environment containing trial and sigma vectors
!>
!> \par Implementation Details
!>      Matrix elements are computed as inner products H_ij = Tr(X_i^T * Sigma_j).
!>      The lower triangular part is calculated explicitly and copied to the upper
!>      triangle by symmetry to reduce computational cost by factor of two.
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE build_subspace(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace'

      REAL(KIND=dp)                                      :: val
      INTEGER                                            :: handle, i, j

      CALL timeset(routineN, handle)

      ! Construct subspace Hamiltonian: H_ij = <X_i | Sigma_j> = Tr(X_i^T * Sigma_j)
      ! Could be DGEMM, but multiplication is between slices of matrix and array of matrices
      ! H[nVec,nVec] = 1.0*(Sigma(nVec)[dim])^T*X(nVec)[dim] + 0.0*H[nVec,nVec]
      DO i = 1, tddft_env%nVec
         DO j = 1, i
            ! Compute matrix element as trace of X_i^T * Sigma_j
            CALL cp_fm_trace(matrix_a=tddft_env%evects_mo(i), &
                             matrix_b=tddft_env%sigma_vecs(j), &
                             trace=val)

            ! Store lower triangle element
            CALL cp_fm_set_element(matrix=tddft_env%subspace, &
                                   irow_global=i, &
                                   icol_global=j, &
                                   alpha=val)

            ! Mirror to upper triangle for i != j
            IF (i/=j) THEN
               CALL cp_fm_set_element(matrix=tddft_env%subspace, &
                                      irow_global=j, &
                                      icol_global=i, &
                                      alpha=val)
            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE build_subspace

! **************************************************************************************************
!> \brief Diagonalizes the subspace Hamiltonian to obtain Ritz values and vectors
!>
!> \details Solves the generalized eigenvalue problem H*C = Omega*C where H is the
!>          subspace Hamiltonian, Omega contains Ritz values (approximate excitation energies),
!>          and C contains eigenVector coefficients for constructing Ritz vectors.
!>
!> \param[inout] tddft_env  TDDFT environment containing subspace matrix and storage
!>                          for eigenvalues and eigenVectors
!>
!> \par Algorithm Details
!>      Uses CP2K's optimal eigensolver selection based on matrix size and parallel
!>      distribution. The eigenVectors are normalized to ensure proper construction
!>      of orthonormal Ritz vectors.
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE diagonalize_subspace(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'diagonalize_subspace'

      INTEGER                                            :: handle, ierr, ivec

      CALL timeset(routineN, handle)

      ! Diagonalize subspace Hamiltonian: H*C = Omega*C
      ierr = 0
      CALL choose_eigv_solver(matrix=tddft_env%subspace, &
                              eigenvalues=tddft_env%ritz_evals(1:tddft_env%nVec), &
                              eigenVectors=tddft_env%subspace_evecs, &
                              info=ierr)

      IF (ierr /= 0) THEN
         CPABORT("Subspace diagonalization failed")
      END IF

      ! Normalize eigenVectors to ensure orthonormal Ritz vectors
      CALL cp_fm_vectorsnorm(tddft_env%subspace_evecs, tddft_env%vector_norm(1:tddft_env%nVec))

      ! Compute normalization factors
      DO ivec = 1, tddft_env%nVec
         IF (tddft_env%vector_norm(ivec) > 1.0E-10_dp) THEN
            tddft_env%vector_scal(ivec) = 1.0_dp / tddft_env%vector_norm(ivec)
         ELSE
            tddft_env%vector_scal(ivec) = 0.0_dp
         END IF
      END DO

      ! Apply normalization to eigenVectors
      CALL cp_fm_column_scale(tddft_env%subspace_evecs, tddft_env%vector_scal(1:tddft_env%nVec))

      CALL timestop(handle)

   END SUBROUTINE diagonalize_subspace

! **************************************************************************************************
!> \brief Forms Ritz vectors from subspace eigenVectors and trial vectors
!>
!> \details Constructs Ritz vectors as linear combinations of trial vectors weighted by
!>          eigenVector coefficients: Z_k = sum_j c_jk X_j. The Ritz vectors represent
!>          the best current approximations to true eigenVectors within the subspace.
!>
!> \param[inout] tddft_env  TDDFT environment containing trial vectors and eigenVector
!>                          coefficients from subspace diagonalization
!>
!> \par Mathematical Background
!>      The Ritz vectors are formed by projecting the exact eigenVectors onto the current
!>      Krylov subspace. They provide improved approximations compared to individual
!>      trial vectors and are used for convergence assessment and restart procedures.
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE form_ritz_vectors(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'form_ritz_vectors'

      REAL(KIND=dp)                                      :: val
      INTEGER                                            :: handle, istate, jvec

      CALL timeset(routineN, handle)

      ! Construct Ritz vectors: Z_k = sum_j c_jk X_j
      DO istate = 1, tddft_env%nstates
         ! Initialize Ritz vector to zero
         CALL cp_fm_set_all(matrix=tddft_env%ritz_vects(istate), alpha=0.0_dp)

         ! Form linear combination of trial vectors
         DO jvec = 1, tddft_env%nVec
             ! Extract coefficient c_jk from eigenVector
             CALL cp_fm_get_element(matrix=tddft_env%subspace_evecs, &
                                    irow_global=jvec, &
                                    icol_global=istate, &
                                    alpha=val)

            ! Accumulate: Z_k += c_jk * X_j
             CALL cp_fm_scale_and_add(&
                     alpha=1.0_dp, matrix_a=tddft_env%ritz_vects(istate), &
                     beta=val, matrix_b=tddft_env%evects_mo(jvec))
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE form_ritz_vectors

! **************************************************************************************************
!> \brief Computes residual vectors for Davidson convergence assessment
!>
!> \details Calculates residuals r_k = Sigma_k - Omega_k Z_k where Sigma_k are the images
!>          of Ritz vectors under the TDDFT Hamiltonian, Omega_k are Ritz values, and Z_k
!>          are Ritz vectors. Small residual norms indicate convergence to true eigenpairs.
!>
!> \param[inout] tddft_env  TDDFT environment containing Ritz vectors, eigenvalues,
!>                          and storage for residual vectors and norms
!>
!> \par Convergence Criterion
!>      The Davidson method converges when ||r_k|| < threshold for all requested states.
!>      The residual represents the error in the eigenvalue equation AZ_k = Omega_k Z_k.
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE compute_residuals(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_residuals'

      REAL(KIND=dp)                                      :: val
      INTEGER                                            :: handle, istate, jvec

      CALL timeset(routineN, handle)

      ! Compute residuals: r_k = Sigma_k - Omega_k Z_k
      DO istate = 1, tddft_env%nstates
         ! Initialize residual vector
         CALL cp_fm_set_all(matrix=tddft_env%residual_vecs(istate), alpha=0.0_dp)

         ! Form Sigma_k = sum_j c_jk Sigma_j (image of Ritz vector)
         DO jvec = 1, tddft_env%nVec
            ! Get eigenVector coefficient
            CALL cp_fm_get_element(matrix=tddft_env%subspace_evecs, &
                                   irow_global=jvec, &
                                   icol_global=istate, &
                                   alpha=val)

            ! Accumulate Sigma_k += c_jk * Sigma_j
            CALL cp_fm_scale_and_add( &
                     alpha=1.0_dp, matrix_a=tddft_env%residual_vecs(istate), &
                     beta=val,   matrix_b=tddft_env%sigma_vecs(jvec))

            ! Subtract Omega_k * c_jk * X_j to form residual
            CALL cp_fm_scale_and_add( &
                     alpha=1.0_dp, matrix_a=tddft_env%residual_vecs(istate), &
                     beta=-val*tddft_env%ritz_evals(istate), &
                     matrix_b=tddft_env%evects_mo(jvec))
         END DO

         ! Compute residual norm ||r_k||
         CALL cp_fm_trace(tddft_env%residual_vecs(istate), &
                          tddft_env%residual_vecs(istate), val)
         tddft_env%residual_norms(istate) = val
      END DO

      ! Track maximum residual for global convergence monitoring
      tddft_env%max_error = MAXVAL(tddft_env%residual_norms(1:tddft_env%nStates))

      CALL timestop(handle)

   END SUBROUTINE compute_residuals

! **************************************************************************************************
!> \brief Logs current status of Ritz values and convergence progress
!>
!> \details Prints a formatted table showing excitation energies and residual norms
!>          for all requested states. This provides real-time monitoring of the
!>          Davidson convergence process.
!>
!> \param[in] tddft_env  TDDFT environment containing current Ritz values and residuals
!> \param[in] log_unit   Output unit for logging (negative values suppress output)
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE log_ritz_status(tddft_env, log_unit)
      TYPE(tddft_info_type), INTENT(IN)                  :: tddft_env
      INTEGER,               INTENT(IN)                  :: log_unit
      INTEGER                                            :: istate

      IF (log_unit <= 0) RETURN

      ! Print table header
      WRITE (log_unit, '(/,4X,A,T14,A,T36,A)') "State", "Exc. energy (eV)", "Convergence (a.u.)"

      ! Print data for each state
      DO istate = 1, tddft_env%nstates
         WRITE (log_unit, '(1X,I8,T12,F14.7,T38,ES11.4)') istate, &
            tddft_env%ritz_evals(istate)*evolt, tddft_env%residual_norms(istate)
      END DO

      WRITE (log_unit, *)
      CALL m_flush(log_unit)

   END SUBROUTINE log_ritz_status

! **************************************************************************************************
!> \brief Assesses convergence based on residual norms and eigenvalue stability
!>
!> \details Implements dual convergence criteria:
!>          1. Residual criterion: max_k ||r_k|| < threshold
!>          2. Eigenvalue criterion: max_k |Omega_k - Omega_k^prev| < threshold/10
!>
!>          The dual approach ensures robust convergence by catching cases where
!>          residual norms may be unstable but eigenvalues have converged.
!>
!> \param[inout] tddft_env      TDDFT environment containing convergence data
!> \param[in]    tddfpt_control Control parameters including convergence threshold
!> \return       converged      .TRUE. if convergence criteria are satisfied
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   FUNCTION is_converged(tddft_env, tddfpt_control) RESULT(converged)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(tddfpt2_control_type), POINTER, INTENT(IN)    :: tddfpt_control
      LOGICAL                                            :: converged

      ! Compute maximum change in eigenvalues from previous iteration
      tddft_env%max_delta_eval = &
         maxval(abs(tddft_env%ritz_evals(1:tddft_env%nstates) &
                  - tddft_env%ritz_evals_prev(1:tddft_env%nstates)))

      ! Apply dual convergence criteria
      ! Primary: residual norm convergence
      ! Secondary: eigenvalue stability (10x tighter threshold)
      converged = (tddft_env%max_error < tddfpt_control%conv) .OR. &
                  (tddft_env%max_delta_eval < tddfpt_control%conv*0.1_dp)

   END FUNCTION is_converged

! **************************************************************************************************
!> \brief Updates working arrays with converged Ritz solution
!>
!> \details Copies converged Ritz vectors to trial vector storage and updates eigenvalues.
!>          This routine is called both upon convergence and before restarts to preserve
!>          the best current approximations.
!>
!> \param[inout] tddft_env  TDDFT environment containing Ritz solution and working arrays
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE commit_ritz_solution(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'commit_ritz_solution'

      INTEGER                                            :: handle, istate

      CALL timeset(routineN, handle)

      ! Transfer Ritz vectors to trial vector storage
      DO istate = 1, tddft_env%nstates
         CALL cp_fm_to_fm(source=tddft_env%ritz_vects(istate), &
                          destination=tddft_env%evects_mo(istate))

         ! Update eigenvalue estimates
         tddft_env%evals(istate) = tddft_env%ritz_evals(istate)
      END DO

      CALL timestop(handle)

   END SUBROUTINE commit_ritz_solution

! **************************************************************************************************
!> \brief Detects convergence stagnation using relative change monitoring
!>
!> \details Monitors convergence progress by tracking relative changes in both residual
!>          norms and eigenvalues between consecutive iterations. Persistent small changes
!>          indicate potential stagnation requiring restart intervention.
!>
!> \param[inout] tddft_env  TDDFT environment containing convergence history
!> \param[in]    iter       Current iteration number
!>
!> \par Stall Detection Logic
!>      Stagnation is detected when both residual and eigenvalue relative changes
!>      fall below 5% for multiple consecutive iterations. After 4 consecutive
!>      stalled iterations, a restart is recommended.
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE detect_stall(tddft_env, iter)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      INTEGER, INTENT(IN)                                :: iter

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'detect_stall'

      REAL(KIND=dp)                                      :: rel_error_change, rel_delta_eval_change
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Calculate relative changes in convergence metrics
      rel_error_change = ABS(tddft_env%max_error_prev &
                           - tddft_env%max_error) / MAX(ABS(tddft_env%max_error_prev), 1.0E-12_dp)
      rel_delta_eval_change = ABS(tddft_env%max_delta_eval_prev &
                           - tddft_env%max_delta_eval) / MAX(ABS(tddft_env%max_delta_eval_prev), 1.0E-12_dp)

      ! Detect stagnation when both metrics show minimal progress
      IF (iter > 1) THEN
         IF (rel_error_change < 5.0E-2 .OR. rel_delta_eval_change < 5.0E-2) THEN
            tddft_env%stall_count = tddft_env%stall_count + 1
         ELSE
            tddft_env%stall_count = 0  ! Reset counter on progress
         END IF
      END IF

      ! Store current values for next iteration comparison
      tddft_env%max_error_prev = tddft_env%max_error
      tddft_env%max_delta_eval_prev = tddft_env%max_delta_eval

      CALL timestop(handle)

   END SUBROUTINE detect_stall

! **************************************************************************************************
!> \brief Determines if Davidson iteration requires restart
!>
!> \details Evaluates restart conditions based on stagnation detection and subspace
!>          size constraints. Restarts help escape local convergence difficulties
!>          and manage memory usage.
!>
!> \param[in] tddft_env  TDDFT environment containing stall and dimension information
!> \return    res        .TRUE. if restart is recommended
!>
!> \par Restart Triggers
!>      - Persistent stagnation (stall_count > 3)
!>      - Maximum subspace dimension reached (nVec >= mxvec)
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   FUNCTION needs_restart(tddft_env) RESULT(res)
      TYPE(tddft_info_type), INTENT(IN)                  :: tddft_env
      LOGICAL                                            :: res

      res = (tddft_env%stall_count > 3) .OR. (tddft_env%nVec >= tddft_env%mxvec)

   END FUNCTION needs_restart

! **************************************************************************************************
!> \brief Applies diagonal preconditioning to residual vectors for enhanced convergence
!>
!> \details Implements preconditioning using the formula r_k^new = r_k / (Omega_k - (e_a - e_i))
!>          where e_a and e_i are virtual and occupied orbital energies. This approximates
!>          the inverse of the orbital Hessian diagonal and accelerates convergence.
!>
!> \param[inout] tddft_env  TDDFT environment containing residuals and orbital energies
!>
!> \par Numerical Stability
!>      Small denominators are regularized to prevent numerical instabilities.
!>      Denominators smaller than 1e-12 are set to 1e-5 while preserving sign.
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE apply_preconditioner(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_preconditioner'

      REAL(KIND=DP)                                      :: denom, scal
      INTEGER                                            :: handle, i, a, ij, istate

      CALL timeset(routineN, handle)

      ! Store eigenvalues for next iteration
      tddft_env%ritz_evals_prev(1:tddft_env%nstates) = tddft_env%ritz_evals(1:tddft_env%nstates)

      ! Apply diagonal preconditioning to all states
      DO istate = 1, tddft_env%nstates
         ij = 0  ! Linear index for (i,a) pairs

         ! Loop over all (occupied, virtual) pairs
         DO a = 1, tddft_env%nvirb
            DO i = 1, tddft_env%nocca
               ij = ij + 1

               ! Compute preconditioner denominator: Omega_k - (mo_e_a - mo_e_i)
               denom = tddft_env%ritz_evals(istate) - tddft_env%e_virt_minus_occ(ij)

               ! Regularize small denominators for numerical stability
               IF (ABS(denom) < 1.e-12_dp) THEN
                  denom = 1.e-5_dp
               ELSE IF (ABS(denom) < 1.e-5_dp) THEN
                  denom = SIGN(1.e-5_dp, denom)  ! Preserve sign
               END IF

               ! Apply preconditioning: r_{ia} = r_{ia} / denom
               CALL cp_fm_get_element(tddft_env%residual_vecs(istate), i, a, scal)
               scal = scal / denom
               CALL cp_fm_set_element(tddft_env%residual_vecs(istate), i, a, scal)
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE apply_preconditioner

! **************************************************************************************************
!> \brief Expands Krylov subspace with orthogonalized, preconditioned residuals
!>
!> \details Performs the following operations to grow the Davidson subspace:
!>          1. Orthogonalize residuals against all existing trial vectors
!>          2. Normalize the orthogonalized residuals
!>          3. Add sufficiently large residuals as new trial vectors
!>          4. Update subspace dimension counters
!>
!> \param[inout] tddft_env  TDDFT environment containing trial vectors and residuals
!>
!> \par Quality Control
!>      Residuals with norm < 1e-9 are rejected to maintain numerical stability.
!>      If no vectors can be added or maximum dimension is reached, a restart is triggered.
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE expand_krylov_space(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'expand_krylov_space'

      REAL(KIND=dp)                                      :: val, scal
      INTEGER                                            :: handle, istate, ivec

      CALL timeset(routineN, handle)

      ! Store current subspace size for operator application tracking
      tddft_env%novec = tddft_env%nVec

      ! Process each residual vector for potential inclusion
      DO istate = 1, tddft_env%nstates
         ! Orthogonalize residual against all existing trial vectors (Gram-Schmidt)
         DO ivec = 1, tddft_env%nVec
            ! Compute projection: <X_j | r_k>
            CALL cp_fm_trace(tddft_env%evects_mo(ivec), &
                             tddft_env%residual_vecs(istate), val)

            ! Remove projection: r_k = r_k - <X_j | r_k> X_j
            CALL cp_fm_scale_and_add( &
                     alpha=1.0_dp, matrix_a=tddft_env%residual_vecs(istate), &
                     beta=-val,     matrix_b=tddft_env%evects_mo(ivec))
         END DO

         ! Compute norm of orthogonalized residual
         CALL cp_fm_trace(tddft_env%residual_vecs(istate), &
                          tddft_env%residual_vecs(istate), val)
         scal = sqrt(val)

         ! Skip vectors with insufficient norm
         ! TODO: maybe replace with some tddfpt thresholds
         IF (scal < 1.e-9_dp) CYCLE

         ! Add normalized residual as new trial vector
         tddft_env%nVec = tddft_env%nVec + 1

         ! Check subspace dimension limits
         IF (tddft_env%nVec == tddft_env%mxvec) THEN
            tddft_env%status(2) = .TRUE.  ! Trigger restart
            CALL timestop(handle)
            RETURN
         END IF

         ! Store normalized residual as new trial vector
         CALL cp_fm_scale_and_add( &
                  alpha=0._dp,     matrix_a=tddft_env%evects_mo(tddft_env%nVec), &
                  beta=1._dp/scal, matrix_b=tddft_env%residual_vecs(istate))
      END DO

      ! Trigger restart if no new vectors were added
      IF ((tddft_env%nVec - tddft_env%novec) == 0) THEN
         tddft_env%status(2) = .TRUE.
      END IF

      CALL timestop(handle)

   END SUBROUTINE expand_krylov_space

! **************************************************************************************************
!> \brief Restarts Davidson iteration by reducing subspace to Ritz vectors
!>
!> \details Implements the Davidson restart procedure by replacing the current trial
!>          vectors with the best Ritz vector approximations. This maintains convergence
!>          progress while reducing memory usage and escaping potential stagnation.
!>
!> \param[inout] tddft_env  TDDFT environment containing current iteration state
!>
!> \par Restart Strategy
!>      The restart preserves the nStates lowest-energy Ritz vectors as the new
!>      trial vector set. All convergence monitoring counters are reset to enable
!>      fresh progress assessment.
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE restart_davidson(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'restart_davidson'

      INTEGER                                            :: istate, handle

      CALL timeset(routineN, handle)

      ! Commit current Ritz solution as new starting point
      CALL commit_ritz_solution(tddft_env)

      ! Reset subspace dimensions to minimal set
      tddft_env%nVec = tddft_env%nstates
      tddft_env%novec = tddft_env%nVec
      tddft_env%stall_count = 0  ! Reset stagnation counter

      ! Clear residual vectors for clean restart
      DO istate = 1, tddft_env%nstates
         CALL cp_fm_set_all(matrix=tddft_env%residual_vecs(istate), alpha=0.0_dp)
      END DO

      CALL timestop(handle)

   END SUBROUTINE restart_davidson

END MODULE qs_tddft_methods
