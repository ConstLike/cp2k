!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddft_types
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create, &
                                              cp_fm_struct_release, &
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create, &
                                              cp_fm_get_info, &
                                              cp_fm_release, &
                                              cp_fm_set_all, &
                                              cp_fm_type
   USE input_constants,                 ONLY: do_potential_coulomb, &
                                              do_potential_long, &
                                              do_potential_mix_cl, &
                                              do_potential_mix_cl_trunc, &
                                              do_potential_short
   USE input_section_types,             ONLY: section_vals_get, &
                                              section_vals_get_subs_vals, &
                                              section_vals_type, &
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env, &
                                              qs_environment_type
   USE qs_tddft_utils,                  ONLY: allocate_fm_matrix, &
                                              deallocate_fm_matrix
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_types'

   PUBLIC :: tddft_info_type
   PUBLIC :: tddft_initialization
   PUBLIC :: tddft_finalization

! **************************************************************************************************
!> \brief General TDDFT buffers and configuration
!>        Holds only what is needed for a Davdiosn iterations,
!>        All data here is private for theory%fields.
!> \author  Konstantin Komarov (constlike@gmail.com)
!> \date    Aug 2025
! **************************************************************************************************
   TYPE tddft_info_type

      !> Type of SCF reference: 1=ROKS/ROHF, 2=UKS/UHF
      INTEGER                                            :: reference_type

      !> Type of SCF reference: 1=SF, 2=MRSF-S, 3=MRSF-T, 4=MRSF-Q
      INTEGER                                            :: response_type

      !> Number of spin channels (1 = restricted, 2 = unrestricted)
      INTEGER                                            :: nSpins

      !> Work fot properties: Dipoles
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: dipole_op_mos_occ

      !--------------------------------------------
      ! DAVIDSON NEEDS
      !--------------------------------------------
      !> Number of states to solve for
      INTEGER                                            :: nstates
      !> Maximum number of states
      INTEGER                                            :: nstates_max !> nocca * nvirb (spin-flip)
      !> Current number of vectors in subspace (Krylov dimension)
      INTEGER                                            :: nvec
      !> The value of nvec at the beginning of the iteration (before adding new vectors),
      !> used to determine how many new vectors are added (nvec - novec).
      INTEGER                                            :: novec
      !> Max number of Krylov vectors per state (from input)
      INTEGER                                            :: maxvec
      !> Maximum number of trial vectors allowed in the Davidson subspace
      INTEGER                                            :: mxvec  !> maxvec * nStates
      !> Number of basis functions (AO)
      INTEGER                                            :: nbf
      !> Number of alpha/beta occupied, virtual orbitals
      INTEGER                                            :: nocca, noccb, nvirb, nvira
      !> Consecutive-iteration stall counter
      INTEGER                                            :: stall_count
      !> Maximum residual norm in current iteration
      REAL(KIND=dp)                                      :: max_error
      !> Maximum change of Ritz eigenvalues in current iteration
      REAL(KIND=dp)                                      :: max_delta_eval
      !> Stored maximum residual norm from previous iteration
      REAL(KIND=dp)                                      :: max_error_prev
      !> Stored maximum eigenvalue change from previous iteration
      REAL(KIND=dp)                                      :: max_delta_eval_prev
      !> Logical flags indicating convergence status:
      !> status(1) = .TRUE. if converged
      !> status(2) = .TRUE. if restart is needed
      LOGICAL, DIMENSION(2)                              :: status

      !> TDDFT excitation eigenvalues (e.g., omega_i)
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: evects
      !> Trial vectors in MO basis: X_k
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: evects_mo
      !> Sigma vectors in MO basis: sigma_k = H X_k
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: sigma_vecs
      !> Ritz vectors in MO basis
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: ritz_vects
      !> Residual vectors for convergence check
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: residual_vecs

      !> Eigenvalues of current subspace (Ritz values)
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: ritz_evals
      !> Eigenvalues from previous iteration (for convergence)
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: ritz_evals_prev
      !> Norms of residual vectors
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: residual_norms

      !> Subspace matrix (Atilde, for diagonalization)
      TYPE(cp_fm_type), POINTER                          :: subspace
      !> Subspace eigenvectors (coefficients for Ritz vectors)
      TYPE(cp_fm_type), POINTER                          :: subspace_evecs

      !> Other working buffers for orthogonalization, normalization, etc.
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: vector_norm
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: vector_scal

      !> Energy difference MO_E_virt - MO_E_occ
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: e_virt_minus_occ

      LOGICAL                                            :: do_hfx
      !> Flag indicating if the functional is range-separated
      LOGICAL                                            :: do_hfx_cam
      REAL(KIND=dp)                                      :: C_hf_alpha
      REAL(KIND=dp)                                      :: C_hf_beta
      REAL(KIND=dp)                                      :: C_hf_omega

   END TYPE tddft_info_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param tddft_env ...
!> \param qs_env ...
!> \param gs_mos ...
! **************************************************************************************************
   SUBROUTINE tddft_initialization(qs_env, tddft_env, gs_mos)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CLASS(tddft_info_type), INTENT(INOUT)              :: tddft_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddft_allocate_data'

      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(section_vals_type), POINTER                   :: hf_section
      INTEGER                                            :: ivec, handle, n_rep_hf, potential_type
      LOGICAL                                            :: explicit

      CALL timeset(routineN, handle)

      ! Get number of HF sections
      hf_section => section_vals_get_subs_vals(section_vals=qs_env%input, &
                                               subsection_name="DFT%XC%HF")
      CALL section_vals_val_get(section_vals=hf_section, &
                                keyword_name="FRACTION", &
                                r_val=tddft_env%C_hf_alpha, &
                                i_rep_section=1)
      CALL section_vals_get(hf_section, n_repetition=n_rep_hf)
      IF (n_rep_hf > 0) THEN
         ! Get FRACTION from the first HF section
         CALL section_vals_val_get(hf_section, "FRACTION", r_val=tddft_env%C_hf_alpha, i_rep_section=1)
         tddft_env%do_hfx = ABS(tddft_env%C_hf_alpha) > 1.0e-10_dp
         IF (tddft_env%do_hfx) THEN
            ! Get POTENTIAL_TYPE, defaults to do_potential_coulomb if not specified
            CALL section_vals_val_get(hf_section, "INTERACTION_POTENTIAL%POTENTIAL_TYPE", &
                                      i_val=potential_type, i_rep_section=1, explicit=explicit)
            IF (.NOT. explicit) potential_type = do_potential_coulomb
            tddft_env%do_hfx_cam = (potential_type == do_potential_mix_cl .OR. &
                                    potential_type == do_potential_mix_cl_trunc .OR. &
                                    potential_type == do_potential_short .OR. &
                                    potential_type == do_potential_long)
         ELSE
            tddft_env%do_hfx_cam = .FALSE.
         END IF
      ELSE
         tddft_env%do_hfx = .FALSE.
         tddft_env%do_hfx_cam = .FALSE.
      END IF

      IF (tddft_env%do_hfx_cam) THEN
         CALL section_vals_val_get(hf_section, "INTERACTION_POTENTIAL%SCALE_COULOMB", &
                                   r_val=tddft_env%C_hf_alpha)
         CALL section_vals_val_get(hf_section, "INTERACTION_POTENTIAL%SCALE_LONGRANGE", &
                                   r_val=tddft_env%C_hf_beta)
         CALL section_vals_val_get(hf_section, "INTERACTION_POTENTIAL%OMEGA", &
                                   r_val=tddft_env%C_hf_omega)
      END IF

      CALL cp_fm_get_info(gs_mos(1)%mos_occ, ncol_global=tddft_env%nocca)
      CALL cp_fm_get_info(gs_mos(2)%mos_occ, ncol_global=tddft_env%noccb)
      CALL cp_fm_get_info(gs_mos(1)%mos_virt, ncol_global=tddft_env%nvira)
      CALL cp_fm_get_info(gs_mos(2)%mos_virt, ncol_global=tddft_env%nvirb)
      CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=tddft_env%nbf)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control)

      ! Set the reference type ROHF/UHF
      tddft_env%reference_type = 0
      IF (dft_control%roks) tddft_env%reference_type = 1
      IF (dft_control%uks) tddft_env%reference_type = 2

      tddft_env%nStates = dft_control%tddfpt2_control%nStates
      tddft_env%nSpins = SIZE(gs_mos, 1)
      tddft_env%nStates_max = tddft_env%nocca*tddft_env%nvirb
      tddft_env%maxvec = dft_control%tddfpt2_control%nkvs
      tddft_env%mxvec = MIN(tddft_env%maxvec*tddft_env%nstates, tddft_env%nStates_max)
      tddft_env%nStates = MIN(tddft_env%nstates, tddft_env%mxvec)
      tddft_env%nvec = MIN(MAX(tddft_env%nStates, 4), tddft_env%mxvec)
      tddft_env%novec = tddft_env%nvec

      IF (tddft_env%nStates > tddft_env%nStates_max) THEN
         CPWARN("NUMBER OF EXCITED STATES COULD LEAD TO PROBLEMS!")
         CPWARN("Experimental: CHANGED nStates TO ITS MAXIMUM VALUE!")
         tddft_env%nStates = tddft_env%nStates_max
      END IF

      CALL get_qs_env(qs_env=qs_env, &
                      blacs_env=blacs_env)

      !--------------------
      ! ALLOCATE  DP
      !--------------------

      ALLOCATE (tddft_env%e_virt_minus_occ(tddft_env%nStates_max), source=0.0_dp)
      ALLOCATE (tddft_env%evals(tddft_env%mxvec), source=0.0_dp)
      ALLOCATE (tddft_env%ritz_evals(tddft_env%mxvec), source=0.0_dp)
      ALLOCATE (tddft_env%ritz_evals_prev(tddft_env%nstates), source=1.0E10_dp)
      ALLOCATE (tddft_env%vector_scal(tddft_env%mxvec), source=0.0_dp)
      ALLOCATE (tddft_env%vector_norm(tddft_env%mxvec), source=0.0_dp)
      ALLOCATE (tddft_env%residual_norms(tddft_env%nstates), source=0.0_dp)

      !--------------------
      ! ALLOCATE  FM
      !--------------------

      !=== evects_mo (mxvec) [nocca, nvirb] ============
      CALL allocate_fm_matrix(matrix=tddft_env%evects_mo, &
                              nmat=tddft_env%mxvec, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      !=== ritz_vects (nStates) [nocca, nvirb] =========
      CALL allocate_fm_matrix(matrix=tddft_env%ritz_vects, &
                              nmat=tddft_env%nStates, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      !=== sigma_vecs (mxvec) [nocca, nvirb] ===========
      CALL allocate_fm_matrix(matrix=tddft_env%sigma_vecs, &
                              nmat=tddft_env%mxvec, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      !=== residual_vecs (nStates) [nocca, nvirb] ======
      CALL allocate_fm_matrix(matrix=tddft_env%residual_vecs, &
                              nmat=tddft_env%nStates, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      !=== evects (1, nVec) [nbf, nocca] ===============
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nocca, &
                               context=blacs_env)
      ALLOCATE (tddft_env%evects(1, tddft_env%nvec))
      DO ivec = 1, tddft_env%nVec
         CALL cp_fm_create(tddft_env%evects(1, ivec), fm_struct)
         CALL cp_fm_set_all(tddft_env%evects(1, ivec), 0.0_dp)
      END DO
      CALL cp_fm_struct_release(fm_struct)
      !=================================================

      CALL timestop(handle)

   END SUBROUTINE tddft_initialization

! **************************************************************************************************
!> \brief ...
!> \param tddft_env ...
! **************************************************************************************************
   SUBROUTINE tddft_finalization(tddft_env)
      CLASS(tddft_info_type), INTENT(INOUT)              :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddft_deallocate_data'

      INTEGER                                            :: ivec, ispin, handle

      CALL timeset(routineN, handle)

      !--------------------
      ! DEALLOCATE  DP
      !--------------------
      IF (ALLOCATED(tddft_env%e_virt_minus_occ)) THEN
         DEALLOCATE (tddft_env%e_virt_minus_occ)
      END IF
      IF (ALLOCATED(tddft_env%evals)) THEN
         DEALLOCATE (tddft_env%evals)
      END IF
      IF (ALLOCATED(tddft_env%ritz_evals)) THEN
         DEALLOCATE (tddft_env%ritz_evals)
      END IF
      IF (ALLOCATED(tddft_env%ritz_evals_prev)) THEN
         DEALLOCATE (tddft_env%ritz_evals_prev)
      END IF
      IF (ALLOCATED(tddft_env%vector_scal)) THEN
         DEALLOCATE (tddft_env%vector_scal)
      END IF
      IF (ALLOCATED(tddft_env%vector_norm)) THEN
         DEALLOCATE (tddft_env%vector_norm)
      END IF
      IF (ALLOCATED(tddft_env%residual_norms)) THEN
         DEALLOCATE (tddft_env%residual_norms)
      END IF

      !--------------------
      ! DEALLOCATE  FM
      !--------------------
      CALL deallocate_fm_matrix(tddft_env%evects_mo)
      CALL deallocate_fm_matrix(tddft_env%ritz_vects)
      CALL deallocate_fm_matrix(tddft_env%sigma_vecs)
      CALL deallocate_fm_matrix(tddft_env%residual_vecs)
      IF (ALLOCATED(tddft_env%evects)) THEN
         DO ivec = 1, SIZE(tddft_env%evects, 2)
            CALL cp_fm_release(tddft_env%evects(1, ivec))
         END DO
         DEALLOCATE (tddft_env%evects)
      END IF
      IF (ALLOCATED(tddft_env%dipole_op_mos_occ)) THEN
         DO ivec = 1, SIZE(tddft_env%dipole_op_mos_occ, 1)
            DO ispin = 1, SIZE(tddft_env%dipole_op_mos_occ, 2)
               CALL cp_fm_release(tddft_env%dipole_op_mos_occ(ivec, ispin))
            END DO
         END DO
         DEALLOCATE (tddft_env%dipole_op_mos_occ)
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddft_finalization

END MODULE qs_tddft_types
