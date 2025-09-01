!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Davidson iterative solver for TDDFT excited states.
!>          Theory-agnostic implementation. Available theories:
!>          - Spin-flip (SF): Excited states with spin contaminations
!>          - Mixed Reference SF (MRSF): Excited states with pure separation to Singlets or Triplets
!> \author  Konstantin Komarov (constlike@gmail.com) (Luber group)
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
!> \brief ...
!> \param[inout] qs_env ...
!> \param[in] calc_forces ...
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

      IF(calc_forces) THEN
      END IF

      NULLIFY (logger)
      logger => cp_get_default_logger()

      NULLIFY (tddfpt_section)
      tddfpt_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT")

      CALL tddfpt_input(qs_env, do_hfx, do_admm, do_exck, &
                        do_hfxsr, do_hfxlr, xc_section, tddfpt_print_section, &
                        lri_section, hfxsr_section)

      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "PROGRAM_BANNER", &
                                      extension=".tddfptLog")

      CALL get_qs_env(qs_env, &
                      blacs_env=blacs_env, &
                      dft_control=dft_control, &
                      matrix_ks=matrix_ks)
      NULLIFY(tddfpt_control)
      tddfpt_control => dft_control%tddfpt2_control
      tddfpt_control%do_hfx = do_hfx
      tddfpt_control%do_admm = do_admm
      tddfpt_control%do_hfxsr = do_hfxsr
      tddfpt_control%hfxsr_primbas = 0
      tddfpt_control%hfxsr_re_int = .TRUE.
      tddfpt_control%do_hfxlr = do_hfxlr
      tddfpt_control%do_exck = do_exck

      CALL cite_reference(Iannuzzi2005)

      CALL tddft_header(log_unit)
      CALL kernel_info(log_unit, dft_control, tddfpt_control, xc_section)

      NULLIFY (gs_mos)
      CALL tddfpt_init_mos(qs_env, gs_mos, log_unit)

      ! Obtain smeared occupation numbers
      IF (tddfpt_control%do_smearing) THEN
         CALL tddfpt_smeared_occupation(qs_env, gs_mos, log_unit)
      END IF

      ! Obtain corrected KS-matrix
      CALL tddfpt_oecorr(qs_env, gs_mos, matrix_ks_oep)

      IF (ASSOCIATED(matrix_ks_oep)) matrix_ks => matrix_ks_oep

      IF ((tddfpt_control%do_lrigpw) .AND. &
          (tddfpt_control%kernel /= tddfpt_kernel_full)) THEN
         CALL cp_abort(__LOCATION__, "LRI only implemented for full kernel")
      END IF

      theory_to_run = ["    SF", "MRSF-S", "MRSF-T"]
!     IF(tddfpt_control%sf_tddft==1) THEN
!        theory_to_run = ["    SF"]
!     ELSE IF(tddfpt_control%mrsf_type==1) THEN
!        IF (tddfpt_control%mrsf_mult==1) theory_to_run = ["MRSF-S"]
!        IF (tddfpt_control%mrsf_mult==3) theory_to_run = ["MRSF-T"]
!     END IF

      DO i = 1, SIZE(theory_to_run)

         CALL tddft_initialization(qs_env, tddft_env, gs_mos)

         logo = theory_to_run(i)
         SELECT CASE (logo)
            CASE('    SF')
               tddft_env%response_type = 1 ! SF
               ALLOCATE(tddft_sf_type :: theory)
            CASE('MRSF-S')
               tddft_env%response_type = 2 ! MRSF-S
               ALLOCATE(tddft_mrsf_type :: theory)
            CASE('MRSF-T')
               tddft_env%response_type = 3 ! MRSF-T
               ALLOCATE(tddft_mrsf_type :: theory)
          END SELECT

         CALL davidson_solver(tddft_env, theory, qs_env, tddfpt_control, logger, &
                              tddfpt_print_section, gs_mos, &
                              blacs_env, logo)

         DEALLOCATE(theory)

         CALL tddft_finalization(tddft_env)

      END DO

      DO i = tddft_env%nspins, 1, -1
         CALL tddfpt_release_ground_state_mos(gs_mos(i))
      END DO
      DEALLOCATE (gs_mos)

      CALL timestop(handle)

   END SUBROUTINE tddft

! **************************************************************************************************
!> \brief ...
!> \param tddft_env ...
!> \param theory_env ...
!> \param qs_env ...
!> \param tddfpt_control ...
!> \param logger ...
!> \param tddfpt_print_section ...
!> \param gs_mos ...
!> \param blacs_env ...
!> \param logo ...
! **************************************************************************************************
   SUBROUTINE davidson_solver(tddft_env, theory_env, qs_env, tddfpt_control, logger, &
                              tddfpt_print_section, gs_mos, blacs_env, logo)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      CLASS(tddft_theory_type)                           :: theory_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      CHARACTER(LEN=6)                                   :: logo

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'davidson_solver'

      CHARACTER(len=20)                                  :: nStates_str
      INTEGER                                            :: handle, handle2, log_unit, &
                                                            iter, iend, ist, ivec
      REAL(kind=dp)                                      :: t1, t2
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: ostrength

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      ! Initialize theory-specific work arrays and parameters
      CALL theory_env%allocate(qs_env, tddft_env)

      IF (tddfpt_control%niters > 0) THEN

         ! Open logger
         log_unit = cp_print_key_unit_nr(logger, &
            tddfpt_print_section, "GUESS_VECTORS", extension=".tddfptLog")

         ! Initialize trial subspace with orthonormalized guess vectors
         CALL theory_env%init_guess(qs_env, tddft_env, gs_mos, log_unit)

         ! Close logger
         CALL cp_print_key_finished_output(log_unit, logger, tddfpt_print_section, "GUESS_VECTORS")

         ! Open logger
         log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, "ITERATION_INFO", extension=".tddfptLog")
         IF (log_unit > 0) THEN
            WRITE (log_unit, '(/,1x,79("-"),/," -",17(" "),A,22(" "),"-",/,1x,79("-"))') &
               TRIM(logo)//"-TDDFT WAVEFUNCTION OPTIMIZATION"
            CALL m_flush(log_unit)
         END IF

         CALL cp_add_iter_level(logger%iter_info, "TDDFT_SCF")

         tddft_env%max_error_prev = 10.0_dp
         tddft_env%max_delta_eval_prev = 10.0_dp

         ! Outer loop for managing restarts
         davidson_outer_loop: DO
            ist = 1
            iend = tddft_env%nvec
            tddft_env%status(1:2) = .FALSE.
            tddft_env%stall_count = 0

            ! Davidson iterations
            davidson_loop: DO

               CALL timeset(routineN//"_inner_loop", handle2)

               CALL cp_iterate(logger%iter_info, &
                               iter_nr_out=iter)

               ! Apply orbital Hessian to new trial vectors
               DO ivec = ist, iend
                  !> MO->AO transformation for trial vectors
                  CALL theory_env%mo_to_ao(tddft_env, gs_mos, ivec)

                  !> Apply all operators (HFX, Coulomb, XC, ...), build Sigma vectors (=A*X) in AO space
                  CALL theory_env%apply(qs_env, tddft_env, tddfpt_control)

                  !> AO->MO transformation for resulting vectors
                  CALL theory_env%ao_to_mo(tddft_env, gs_mos, ivec)
               END DO

               ! Dynamic allocation of subspace matrix for current size
               CALL allocate_fm_matrix(tddft_env%subspace, tddft_env%nvec, tddft_env%nvec, blacs_env)
               CALL allocate_fm_matrix(tddft_env%subspace_evecs, tddft_env%nvec, tddft_env%nvec, blacs_env)

               ! Construct subspace matrix H_ij = <X_i|Sigma_j>
               CALL build_subspace(tddft_env)

               ! Diagonalize H*C = Omega*C
               CALL diagonalize_subspace(tddft_env)

               ! Form Ritz vectors Z_k = Σ_j c_{jk} X_j
               CALL form_ritz_vectors(tddft_env)

               ! Compute residuals r_k = Sigma_k - Omega_k Z_k
               CALL compute_residuals(tddft_env)

               CALL deallocate_fm_matrix(tddft_env%subspace)
               CALL deallocate_fm_matrix(tddft_env%subspace_evecs)

               ! Print status of trial vectors at current iteration
               CALL log_ritz_status(tddft_env, log_unit)

               ! Check convergence: max||r_k|| < threshold .OR. max|Energy_diff_k| < threshold/10
               tddft_env%status(1) = is_converged(tddft_env, tddfpt_control)

               IF (tddft_env%status(1)) THEN

                  CALL commit_ritz_solution(tddft_env)

               ELSE

                  ! Detect stall based on relative changes in residuals and eigenvalues
                  CALL detect_stall(tddft_env, iter)

                  ! Check if restart needed (stalled or nvec >= max)
                  tddft_env%status(2) = needs_restart(tddft_env)

                  IF (.NOT. tddft_env%status(2)) THEN

                     ! Precondition residuals: r_k = P^{-1} r_k
                     CALL apply_preconditioner(tddft_env)

                     ! Orthogonalize and expand Krylov space
                     CALL expand_krylov_space(tddft_env)

                  END IF

                  IF (tddft_env%status(2)) THEN

                     ! Restart with Ritz vectors as new trial guess vectors
                     CALL restart_davidson(tddft_env)

                  END IF

               END IF

               !> Print iteration time and max_error
               t2 = m_walltime()
               IF (log_unit > 0) THEN
                  WRITE (log_unit, '(T7,I8,T24,F7.1,T40,ES11.4,T66)') iter, t2 - t1, tddft_env%max_error
                  CALL m_flush(log_unit)
               END IF
               t1 = t2

               CALL timestop(handle2)

               ! Check convergence status
               IF (tddft_env%status(1) .OR. tddft_env%status(2)) THEN
                  EXIT    ! Exit davidson_loop
               END IF

               ! Update indices for next iteration
               ist = tddft_env%novec+1
               iend = tddft_env%nvec

            END DO davidson_loop ! Davidson iterations

            ! Check conditions for exiting or restarting
            ! at this point at least one of the following conditions are met:
            ! a) convergence criteria has been achieved;
            ! b) maximum number of iterations has been reached;
            ! c) Davidson iterations must be restarted due to lack of Krylov vectors

            CALL cp_iterate(logger%iter_info, increment=0, iter_nr_out=iter)

            ! Terminate the loop if either (a) or (b) is True
            IF (tddft_env%status(1) .OR. iter >= tddfpt_control%niters) EXIT

            IF (log_unit > 0) THEN
               WRITE (log_unit, '(1X,25("-"),1X,A,1X,25("-"))') "Restart Davidson iterations"
               CALL m_flush(log_unit)
            END IF

         END DO davidson_outer_loop ! Outer loop for managing restarts

         CALL cp_rm_iter_level(logger%iter_info, "TDDFT_SCF")

         ! print convergence summary
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
         CALL cp_print_key_finished_output(log_unit, logger, &
                                           tddfpt_print_section, "ITERATION_INFO")
      ELSE
         CALL cp_warn(__LOCATION__, &
                      "Skipping TDDFPT wavefunction optimization")
      END IF
      ! *** print summary information ***
      log_unit = cp_logger_get_default_io_unit(logger)
      ALLOCATE (ostrength(tddft_env%nstates), source=0.0_dp)

      CALL tddfpt_dipole_operator(tddft_env%dipole_op_mos_occ, &
                                  tddfpt_control, &
                                  gs_mos, &
                                  qs_env)

      CALL theory_env%print(qs_env, tddft_env, ostrength, gs_mos, tddfpt_control, log_unit)

      CALL theory_env%deallocate(qs_env)

      DEALLOCATE (ostrength)

      CALL timestop(handle)

   END SUBROUTINE davidson_solver

! **************************************************************************************************
!> \brief Builds the subspace Hamiltonian matrix for Davidson diagonalization.
!>        Implements H_ij = <X_i|Sigma_j>
!>        where X_i are trial vectors and Sigma_j = A*X_j are their images under
!>        the orbital Hessian operator.
!> \param[inout] tddft_env Davidson solver environment containing trial and sigma vectors
!> \note The subspace matrix is symmetric, so only the lower triangle is computed
!>       and then mirrored to the upper triangle.
! **************************************************************************************************
   SUBROUTINE build_subspace(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_subspace'

      REAL(KIND=dp)                                      :: val
      INTEGER                                            :: handle, i, j

      CALL timeset(routineN, handle)

      ! Build subspace matrix: H_ij = <X_i | Sigma_j>
      ! Could be DGEMM, but multiplication is between slices of matrix and array of matrices
      ! H[nvec,nvec] = 1.0*(Sigma(nvec)[dim])^T*X(nvec)[dim] + 0.0*H[nvec,nvec]
      DO i = 1, tddft_env%nvec
         DO j = 1, i
            CALL cp_fm_trace(matrix_a=tddft_env%evects_mo(i), &
                             matrix_b=tddft_env%sigma_vecs(j), &
                             trace=val)
            ! Store in subspace matrix (symmetric)
            CALL cp_fm_set_element(matrix=tddft_env%subspace, &
                                   irow_global=i, &
                                   icol_global=j, &
                                   alpha=val)
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
!> \brief Diagonalizes the subspace Hamiltonian to obtain Ritz values and vectors.
!>        Implements: Diagonalize H*C = Omega*C
!>        The eigenvalues (Omega) are the Ritz values approximating excitation energies,
!>        and eigenvectors (C) are coefficients for constructing Ritz vectors.
!> \param[inout] tddft_env Davidson solver environment containing subspace matrix
!> \note Uses the optimal eigensolver based on matrix size and parallel distribution.
! **************************************************************************************************
   SUBROUTINE diagonalize_subspace(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'diagonalize_subspace'

      INTEGER                                            :: handle, ierr, ivec

      CALL timeset(routineN, handle)

      ! Diagonalize H*C = Omega*C
      ierr = 0
      CALL choose_eigv_solver(matrix=tddft_env%subspace, &
                              eigenvalues=tddft_env%ritz_evals(1:tddft_env%nvec), &
                              eigenvectors=tddft_env%subspace_evecs, &
                              info=ierr)

      IF (ierr /= 0) THEN
         CPABORT("Subspace diagonalization failed")
      END IF

      ! Normalize eigenvectors
      CALL cp_fm_vectorsnorm(tddft_env%subspace_evecs, tddft_env%vector_norm(1:tddft_env%nvec))
      DO ivec = 1, tddft_env%nvec
         IF (tddft_env%vector_norm(ivec) > 1.0E-10_dp) THEN
            tddft_env%vector_scal(ivec) = 1.0_dp / tddft_env%vector_norm(ivec)
         ELSE
            tddft_env%vector_scal(ivec) = 0.0_dp
         END IF
      END DO

      ! Apply normalization
      CALL cp_fm_column_scale(tddft_env%subspace_evecs, tddft_env%vector_scal(1:tddft_env%nvec))

      CALL timestop(handle)

   END SUBROUTINE diagonalize_subspace

! **************************************************************************************************
!> \brief Forms Ritz vectors from subspace eigenvectors.
!>        Implements Z_k = Σ_j c_{jk} X_j
!>        The Ritz vectors are linear combinations of trial vectors weighted by
!>        the eigenvector coefficients from subspace diagonalization.
!> \param[inout] tddft_env Davidson solver environment containing trial vectors and eigenvectors
!> \note Only computes the nstates lowest Ritz vectors needed for convergence.
! **************************************************************************************************
   SUBROUTINE form_ritz_vectors(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'form_ritz_vectors'

      REAL(KIND=dp)                                      :: val
      INTEGER                                            :: handle, istate, jvec

      CALL timeset(routineN, handle)

      ! Form Ritz vectors: Z_k = Sigma_j c_{jk} X_j
      DO istate = 1, tddft_env%nstates
         CALL cp_fm_set_all(matrix=tddft_env%ritz_vects(istate), alpha=0.0_dp)
         DO jvec = 1, tddft_env%nvec
             CALL cp_fm_get_element(matrix=tddft_env%subspace_evecs, &
                                    irow_global=jvec, &
                                    icol_global=istate, &
                                    alpha=val)
            ! Z_k += c_jk * X_j
             CALL cp_fm_scale_and_add(&
                     alpha=1.0_dp, matrix_a=tddft_env%ritz_vects(istate), &
                     beta=val, matrix_b=tddft_env%evects_mo(jvec))
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE form_ritz_vectors

! **************************************************************************************************
!> \brief Computes residual vectors for convergence assessment.
!>        Implements r_k = Sigma_k - Omega_k Z_k
!>        The residuals measure how well the Ritz pairs approximate true eigenpairs.
!>        Small residual norms indicate convergence.
!> \param[inout] tddft_env Davidson solver environment containing Ritz vectors and values
!> \note Also computes Sigma_k = Sigma_j c_{jk} Sigma_j for efficiency.
! **************************************************************************************************
   SUBROUTINE compute_residuals(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_residuals'

      REAL(KIND=dp)                                      :: val
      INTEGER                                            :: handle, istate, jvec

      CALL timeset(routineN, handle)

      ! Compute residuals: r_k = Sigma_k - Omega_k Z_k
      DO istate = 1, tddft_env%nstates
         ! Initialize residual to zero
         CALL cp_fm_set_all(matrix=tddft_env%residual_vecs(istate), alpha=0.0_dp)

         ! Form Sigma_k = Sigma_j c_{jk} Sigma_j
         DO jvec = 1, tddft_env%nvec
            ! Get coefficient c_{jk}
            CALL cp_fm_get_element(matrix=tddft_env%subspace_evecs, &
                                   irow_global=jvec, &
                                   icol_global=istate, &
                                   alpha=val)

            ! Add c_{jk} * Sigma_j
            CALL cp_fm_scale_and_add( &
                     alpha=1.0_dp, matrix_a=tddft_env%residual_vecs(istate), &
                     beta=val,   matrix_b=tddft_env%sigma_vecs(jvec))

            ! Subtract Omega_k*Z_k
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

      ! Update maximum residual for convergence check
      tddft_env%max_error = MAXVAL(tddft_env%residual_norms(1:tddft_env%nStates))

      CALL timestop(handle)

   END SUBROUTINE compute_residuals

! **************************************************************************************************
!> \brief Logs current Ritz values and residual norms to output.
!> \param[in] tddft_env Davidson solver environment
!> \param[in] log_unit Output unit number
! **************************************************************************************************
   SUBROUTINE log_ritz_status(tddft_env, log_unit)
      TYPE(tddft_info_type), INTENT(IN)                  :: tddft_env
      INTEGER,               INTENT(IN)                  :: log_unit
      INTEGER                                            :: istate

      IF (log_unit <= 0) RETURN

      WRITE (log_unit, '(/,4X,A,T14,A,T36,A)') "State", "Exc. energy (eV)", "Convergence (a.u.)"
      DO istate = 1, tddft_env%nstates
         WRITE (log_unit, '(1X,I8,T12,F14.7,T38,ES11.4)') istate, &
            tddft_env%ritz_evals(istate)*evolt, tddft_env%residual_norms(istate)
      END DO
      WRITE (log_unit, *)
      CALL m_flush(log_unit)

   END SUBROUTINE log_ritz_status

! **************************************************************************************************
!> \brief Checks convergence criteria for Davidson iteration.
!>        Converged if max_k ||r_k|| < threshold .OR. max_k |Omega_k-Omega_old_k| < threshold/10
!> \param[inout] tddft_env Davidson solver environment containing residual norms and eigenvalues
!> \param[in] tddfpt_control Control parameters including convergence threshold
!> \retval[out] converged .TRUE. if converged based on residual OR eigenvalue criteria
!> \note The dual criteria ensure robust convergence: residual norm for accuracy,
!>       eigenvalue stability for unstable residual norm of states.
! **************************************************************************************************
   FUNCTION is_converged(tddft_env, tddfpt_control) RESULT(converged)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(tddfpt2_control_type), POINTER, INTENT(IN)    :: tddfpt_control
      LOGICAL                                            :: converged

      ! Compute max energy difference between this and previous iteration
      tddft_env%max_delta_eval = &
         maxval(abs(tddft_env%ritz_evals(1:tddft_env%nstates) &
                  - tddft_env%ritz_evals_prev(1:tddft_env%nstates)))

      ! Check convergence criteria
      ! max_k ||r_k|| < threshold OR max_k |Omega_k-Omega_old_k| < threshold/10.0
      converged = (tddft_env%max_error < tddfpt_control%conv) .OR. &
                  (tddft_env%max_delta_eval < tddfpt_control%conv*0.1_dp)

   END FUNCTION is_converged

! **************************************************************************************************
!> \brief Commits converged Ritz solution to working arrays.
!>        Updates the trial vectors and eigenvalues with the converged Ritz pairs.
!>        Called both upon convergence and before restarts to preserve best approximations.
!> \param[inout] tddft_env Davidson solver environment
! **************************************************************************************************
   SUBROUTINE commit_ritz_solution(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'commit_ritz_solution'

      INTEGER                                            :: handle, istate

      CALL timeset(routineN, handle)

      ! Update guess vectors to latest Ritz vectors and eigenvalues
      DO istate = 1, tddft_env%nstates
         CALL cp_fm_to_fm(source=tddft_env%ritz_vects(istate), &
                          destination=tddft_env%evects_mo(istate))

         tddft_env%evals(istate) = tddft_env%ritz_evals(istate)
      END DO

      CALL timestop(handle)

   END SUBROUTINE commit_ritz_solution

! **************************************************************************************************
!> \brief Detects iteration stall based on relative change criteria.
!>        Monitors convergence progress by tracking relative changes in residual norms
!>        and eigenvalues between iterations. Small changes indicate potential stagnation.
!> \param[inout] tddft_env Davidson solver environment
!> \param[in] iter Current iteration number
!> \note Stall detection criteria:
!>       - Relative residual change < 5%
!>       - Relative eigenvalue change < 5%
!>       After 4 consecutive stalls, restart is triggered to escape local stagnation.
! **************************************************************************************************
   SUBROUTINE detect_stall(tddft_env, iter)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      INTEGER, INTENT(IN)                                :: iter

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'detect_stall'

      REAL(KIND=dp)                                      :: rel_error_change, rel_delta_eval_change
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Compute relative changes for stall detection
      rel_error_change = ABS(tddft_env%max_error_prev &
                           - tddft_env%max_error) / MAX(ABS(tddft_env%max_error_prev), 1.0E-12_dp)
      rel_delta_eval_change = ABS(tddft_env%max_delta_eval_prev &
                           - tddft_env%max_delta_eval) / MAX(ABS(tddft_env%max_delta_eval_prev), 1.0E-12_dp)

      ! Stall detection: both error and eigenvalue change must be small
      IF (iter > 1) THEN
         IF (rel_error_change < 5.0E-2 .OR. rel_delta_eval_change < 5.0E-2) THEN
            tddft_env%stall_count = tddft_env%stall_count + 1
         ELSE
            tddft_env%stall_count = 0
         END IF
      END IF

      ! Save for next iteration
      tddft_env%max_error_prev = tddft_env%max_error
      tddft_env%max_delta_eval_prev = tddft_env%max_delta_eval

      CALL timestop(handle)

   END SUBROUTINE detect_stall

! **************************************************************************************************
!> \brief Checks if Davidson iteration needs restart.
!>        Checks restart conditions
!> \param[in] tddft_env Davidson solver environment
!> \retval[out] res .TRUE. if restart needed due to stall or reaching max vectors
!> \note Restart is triggered when:
!>       - Stall counter > 3 (persistent stagnation)
!>       - Number of vectors >= maximum allowed
! **************************************************************************************************
   FUNCTION needs_restart(tddft_env) RESULT(res)
      TYPE(tddft_info_type), INTENT(IN)                  :: tddft_env
      LOGICAL                                            :: res

      res = (tddft_env%stall_count > 3) .OR. (tddft_env%nvec >= tddft_env%mxvec)

   END FUNCTION needs_restart

! **************************************************************************************************
!> \brief Applies diagonal preconditioning to Davidson residuals to accelerate convergence.
!>        Implements r_k = P^{-1} r_k
!>        The preconditioner is P^{-1} where P_{ia} = Omega_k - (mo_e_a - mo_e_i)
!> \param[inout] tddft_env Davidson solver environment containing residual vectors
!> \note This preconditioning step improves convergence by scaling residuals based on
!>       the difference between current Ritz values and orbital energy differences.
! **************************************************************************************************
   SUBROUTINE apply_preconditioner(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_preconditioner'

      REAL(KIND=DP)                                      :: denom, scal
      INTEGER                                            :: handle, i, a, ij, istate

      CALL timeset(routineN, handle)

      ! Save for next iteration round:
      tddft_env%ritz_evals_prev(1:tddft_env%nstates) = tddft_env%ritz_evals(1:tddft_env%nstates)

      ! Apply preconditioner P^{-1} to residuals
      DO istate = 1, tddft_env%nstates
         ij = 0
         DO a = 1, tddft_env%nvirb
            DO i = 1, tddft_env%nocca
               ij = ij + 1

               ! Calculate the denominator of the preconditioner
               denom = tddft_env%ritz_evals(istate) - tddft_env%e_virt_minus_occ(ij)

               ! Regularization for numerical stability when denominator approaches zero
               IF (ABS(denom) < 1.e-12_dp) THEN
                  denom = 1.e-5_dp
               ELSE IF (ABS(denom) < 1.e-5_dp) THEN
                  denom = SIGN(1.e-5_dp, denom)  ! Keep sign
               END IF

               ! Apply the preconditioner: r_{ia} = r_{ia} / (Omega_k - (e_a - e_i))
               CALL cp_fm_get_element(tddft_env%residual_vecs(istate), i, a, scal)
               scal = scal / denom
               CALL cp_fm_set_element(tddft_env%residual_vecs(istate), i, a, scal)
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE apply_preconditioner

! **************************************************************************************************
!> \brief Expands Krylov space with preconditioned, orthogonalized residuals.
!>        Implements:
!>        - Orthogonalize residuals against existing trial vectors
!>        - Normalize residuals
!>        - Add to trial vector set
!> \param[inout] tddft_env Davidson solver environment
!> \note Residuals with norm < 1e-9 are skipped to maintain numerical stability.
!>       If no new vectors can be added, triggers restart.
! **************************************************************************************************
   SUBROUTINE expand_krylov_space(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'expand_krylov_space'

      REAL(KIND=dp)                                      :: val, scal
      INTEGER                                            :: handle, istate, ivec

      CALL timeset(routineN, handle)

      ! Save current vector count for operator application range
      tddft_env%novec = tddft_env%nvec

      ! Krylov expansion with orthogonalized residuals
      DO istate = 1, tddft_env%nstates
         ! Orthogonalize residual against all existing trial vectors
         DO ivec = 1, tddft_env%nvec
            CALL cp_fm_trace(tddft_env%evects_mo(ivec), &
                             tddft_env%residual_vecs(istate), val)
            CALL cp_fm_scale_and_add( &
                     alpha=1.0_dp, matrix_a=tddft_env%residual_vecs(istate), &
                     beta=-val,     matrix_b=tddft_env%evects_mo(ivec))
         END DO

         ! Normalize the orthogonalized residual
         CALL cp_fm_trace(tddft_env%residual_vecs(istate), &
                          tddft_env%residual_vecs(istate), val)
         scal = sqrt(val)

         ! Skip if residual norm too small
         IF (scal<1.e-9_dp) CYCLE

         tddft_env%nvec = tddft_env%nvec + 1

         ! Check if reached maximum vectors
         IF (tddft_env%nvec == tddft_env%mxvec) THEN

            tddft_env%status(2) = .TRUE.
            CALL timestop(handle)
            RETURN

         END IF

         ! Add normalized residual as new trial vector
         CALL cp_fm_scale_and_add( &
                  alpha=0._dp,     matrix_a=tddft_env%evects_mo(tddft_env%nvec), &
                  beta=1._dp/scal, matrix_b=tddft_env%residual_vecs(istate))
      END DO

      ! If no new vectors added, trigger restart
      IF ((tddft_env%nvec - tddft_env%novec) == 0) THEN
         tddft_env%status(2) = .TRUE.
      END IF

      CALL timestop(handle)

   END SUBROUTINE expand_krylov_space

! **************************************************************************************************
!> \brief Restarts Davidson iteration with Ritz vectors.
!>        Implements Algorithm 1, Line 20: X_j = Z_j for j ≤ nstates
!>        Reduces the Krylov space to the best current approximations (Ritz vectors)
!>        to continue iteration from a compact, high-quality subspace.
!> \param[inout] tddft_env Davidson solver environment
!> \note Preserves convergence progress while escaping stagnation or memory limits.
! **************************************************************************************************
   SUBROUTINE restart_davidson(tddft_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'restart_davidson'

      INTEGER                                            :: istate, handle

      CALL timeset(routineN, handle)

      ! Prepare to restart Davidson

      CALL commit_ritz_solution(tddft_env)

      tddft_env%nvec = tddft_env%nstates  ! Start with nstates vectors
      tddft_env%novec = tddft_env%nvec  ! Reset to 0 to force operator re-application
      tddft_env%stall_count = 0

      ! Reset residuals
      DO istate = 1, tddft_env%nstates
         CALL cp_fm_set_all(matrix=tddft_env%residual_vecs(istate), alpha=0.0_dp)
      END DO

      CALL timestop(handle)

   END SUBROUTINE restart_davidson

END MODULE qs_tddft_methods
