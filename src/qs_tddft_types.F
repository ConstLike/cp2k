!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddft_types
   USE input_constants,                 ONLY: &
        do_hfx_auto_shells, do_potential_coulomb, do_potential_gaussian, do_potential_id, &
        do_potential_long, do_potential_mix_cl, do_potential_mix_cl_trunc, do_potential_mix_lg, &
        do_potential_short, do_potential_truncated, hfx_ri_do_2c_diag, hfx_ri_do_2c_iter
   USE hfx_types,                       ONLY: hfx_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace,&
                                              integrate_v_rspace_one_center
   USE qs_tddfpt2_operators,            ONLY: tddfpt_apply_coulomb,&
                                              tddfpt_apply_xc,&
                                              tddfpt_apply_xc_potential
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho,&
                                              qs_rho_update_tddfpt
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_zero
   USE pw_env_types,                    ONLY: pw_env_get
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_trace,&
                                              cp_fm_transpose,&
                                              cp_fm_column_scale,&
                                              cp_fm_row_scale
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              copy_dbcsr_to_fm,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE cp_dbcsr_api,                    ONLY: dbcsr_type,&
                                              dbcsr_create,&
                                              dbcsr_get_info,&
                                              dbcsr_set,&
                                              dbcsr_release,&
                                              dbcsr_p_type,&
                                              dbcsr_init_p,&
                                              dbcsr_release_p,&
                                              dbcsr_complete_redistribute,&
                                              dbcsr_distribution_release,&
                                              dbcsr_distribution_type,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_deallocate_matrix
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_dbcsr_create_by_dist,&
                                              tddfpt_subgroup_env_type
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_get,&
                                              qs_rho_type
   USE util,                            ONLY: sort
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_p_type,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE admm_methods,                    ONLY: admm_fit_mo_coeffs
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Grimme2013,&
                                              Grimme2016,&
                                              Iannuzzi2005,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_create,&
                                              fm_pool_create_fm,&
                                              fm_pool_release, &
                                              fm_pool_give_back_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type,&
                                              cp_fm_set_element, &
                                              cp_fm_get_element, &
                                              cp_fm_get_submatrix,&
                                              cp_fm_write_formatted,&
                                              cp_fm_write_unformatted,&
                                              cp_fm_vectorsnorm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE exstates_types,                  ONLY: excited_energy_type
   USE header,                          ONLY: tddfpt_header,&
                                              tddfpt_soc_header
   USE hfx_admm_utils,                  ONLY: aux_admm_init
   USE hfx_types,                       ONLY: compare_hfx_sections,&
                                              hfx_create
   USE input_constants,                 ONLY: &
        do_admm_aux_exch_func_none, do_admm_basis_projection, do_admm_exch_scaling_none, &
        do_admm_purify_none, do_potential_truncated, tddfpt_dipole_velocity, tddfpt_kernel_full, &
        tddfpt_kernel_none, tddfpt_kernel_stda,&
                                              no_sf_tddfpt,&
                                              tddfpt_sf_col,&
                                              tddfpt_sf_noncol
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE lri_environment_methods,         ONLY: lri_print_stat
   USE lri_environment_types,           ONLY: lri_density_release,&
                                              lri_env_release
   USE message_passing,                 ONLY: mp_para_env_type
   USE min_basis_set,                   ONLY: create_minbas_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kernel_methods,               ONLY: create_fxc_kernel,&
                                              create_kernel_env
   USE qs_kernel_types,                 ONLY: full_kernel_env_type,&
                                              kernel_env_type,&
                                              release_kernel_env
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: mo_set_type
   USE qs_scf_methods,                  ONLY: eigensolver
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_tddfpt2_assign,               ONLY: assign_state
   USE qs_tddfpt2_densities,            ONLY: tddfpt_construct_aux_fit_density,&
                                              tddfpt_construct_ground_state_orb_density
   USE qs_tddfpt2_eigensolver,          ONLY: tddfpt_davidson_solver,&
                                              tddfpt_orthogonalize_psi1_psi0,&
                                              tddfpt_orthonormalize_psi1_psi1
   USE qs_tddfpt2_forces,               ONLY: tddfpt_forces_main
   USE qs_tddfpt2_fprint,               ONLY: tddfpt_print_forces
   USE qs_tddfpt2_lri_utils,            ONLY: tddfpt2_lri_init
   USE qs_tddfpt2_properties,           ONLY: tddfpt_dipole_operator,&
                                              tddfpt_print_excitation_analysis,&
                                              tddfpt_print_exciton_descriptors,&
                                              tddfpt_print_nto_analysis,&
                                              tddfpt_print_summary
   USE qs_tddfpt2_restart,              ONLY: tddfpt_read_restart,&
                                              tddfpt_write_newtonx_output,&
                                              tddfpt_write_restart
   USE qs_tddfpt2_smearing_methods,     ONLY: tddfpt_smeared_occupation
   USE qs_tddfpt2_soc,                  ONLY: tddfpt_soc
   USE qs_tddfpt2_stda_types,           ONLY: allocate_stda_env,&
                                              deallocate_stda_env,&
                                              stda_env_type,&
                                              stda_init_param
   USE qs_tddfpt2_stda_utils,           ONLY: get_lowdin_mo_coefficients,&
                                              stda_init_matrices
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_sub_env_init,&
                                              tddfpt_sub_env_release,&
                                              tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: hfxsr_create_work_matrices,&
                                              stda_create_work_matrices,&
                                              tddfpt_create_work_matrices,&
                                              tddfpt_ground_state_mos,&
                                              tddfpt_release_work_matrices,&
                                              tddfpt_work_matrices
   USE qs_tddfpt2_utils,                ONLY: &
                                              tddfpt_init_mos,&
                                              tddfpt_oecorr,&
                                              tddfpt_release_ground_state_mos
   USE string_utilities,                ONLY: integer_to_string
   USE xc_write_output,                 ONLY: xc_write

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_types'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

   PUBLIC :: tddft_info_type
   PUBLIC :: tddft_initialization
   PUBLIC :: tddft_finalization
   PUBLIC :: allocate_fm_matrix
   PUBLIC :: deallocate_fm_matrix

! **************************************************************************************************
!> \brief General TDDFT buffers and configuration
!>        Holds only what is needed for a Davdiosn iterations,
!>        All data here is private for theory_handler%fields.
!> \author  Konstantin Komarov (constlike@gmail.com)
!> \date    2025-05-18
! **************************************************************************************************
   TYPE tddft_info_type

      !> Type of SCF reference: 1=ROKS/ROHF, 2=UKS/UHF
      INTEGER :: reference_type

      !> Type of SCF reference: 1=SF, 2=MRSF-S, 3=MRSF-T, 4=MRSF-Q
      INTEGER :: response_type

      !> Number of spin channels (1 = restricted, 2 = unrestricted)
      INTEGER :: nSpins

      !> Fock matrix for Alpha occupied orbitals block [nocca, nocca]
      TYPE(cp_fm_type), POINTER :: Focka_noc_noc => NULL()

      !> Fock matrix for Beta virtual orbitals block [nvirb, nvirb]
      TYPE(cp_fm_type), POINTER :: Fockb_vir_vir => NULL()

      !> Matrix pool for MO-space matrices [nocca, nvirb]
      TYPE(cp_fm_pool_p_type), POINTER :: fm_pool_mo => NULL()

      !> Matrix pools for AO-space matrices [nbf, nbf] per spin channel
      TYPE(cp_fm_pool_p_type), POINTER, DIMENSION(:) :: fm_pool_ao => NULL()

      !> Per-spin AO scratch matrices for Fock transformation [nbf, nbf]
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: work_ao_ao

      !> Per-spin AO Fock matrices (used in transformation to MO)
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: fock

      !> Work fot properties: Dipoles
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :) :: dipole_op_mos_occ

      !> Per-spin AO scratch matrices for Fock transformation [nbf, nbf]
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: dens_diff_ij
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: dens_diff_ab

      !--------------------------------------------
      ! DAVIDSON NEEDS
      !--------------------------------------------
      !> Number of states to solve for
      INTEGER :: nstates
      !> Maximum number of states
      INTEGER :: nstates_max != nocca * nvirb (spin-flip)
      !> Current number of vectors in subspace (Krylov dimension)
      INTEGER :: nvec
      !> The value of nvec at the beginning of the iteration (before adding new vectors),
      !> used to determine how many new vectors are added (nvec - novec).
      INTEGER :: novec
      !> Max number of Krylov vectors per state (from input)
      INTEGER :: maxvec
      !> Maximum number of trial vectors allowed in the Davidson subspace
      INTEGER :: mxvec  != maxvec * nStates
      !> Number of basis functions (AO)
      INTEGER :: nbf
      !> Number of alpha/beta occupied, virtual orbitals
      INTEGER :: nocca, noccb, nvirb, nvira
      !> Consecutive-iteration stall counter
      INTEGER :: stall_count
      !> Maximum residual norm in current iteration
      REAL(KIND=dp) :: max_error
      !> Maximum change of Ritz eigenvalues in current iteration
      REAL(KIND=dp) :: max_delta_eval
      !> Stored maximum residual norm from previous iteration
      REAL(KIND=dp) :: max_error_prev
      !> Stored maximum eigenvalue change from previous iteration
      REAL(KIND=dp) :: max_delta_eval_prev
      !> Logical flags indicating convergence status:
      !> status(1) = .TRUE. if converged
      !> status(2) = .TRUE. if restart is needed
      LOGICAL, DIMENSION(2) :: status

      !> TDDFT excitation eigenvalues (e.g., omega_i)
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: evals
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:,:) :: evects
      !> Trial vectors in MO basis: X_k
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: evects_mo
      !> Sigma vectors in MO basis: sigma_k = H X_k
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: sigma_vecs
      !> Ritz vectors in MO basis
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: ritz_vects
      !> Residual vectors for convergence check
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: residual_vecs

      !> Eigenvalues of current subspace (Ritz values)
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: ritz_evals
      !> Eigenvalues from previous iteration (for convergence)
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: ritz_evals_prev
      !> Norms of residual vectors
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: residual_norms

      !> Subspace matrix (Atilde, for diagonalization)
      TYPE(cp_fm_type), POINTER :: subspace
      !> Subspace eigenvectors (coefficients for Ritz vectors)
      TYPE(cp_fm_type), POINTER :: subspace_evecs

      !> Other working buffers for orthogonalization, normalization, etc.
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: vector_norm
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: vector_scal

      !> Energy difference E_virt - E_occ
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: e_virt_minus_occ

      LOGICAL :: do_hfx
      !> Flag indicating if the functional is range-separated
      LOGICAL :: do_hfx_cam
      REAL(KIND=dp) :: C_hf_alpha
      REAL(KIND=dp) :: C_hf_beta
      REAL(KIND=dp) :: C_hf_omega

   END TYPE tddft_info_type

CONTAINS

   SUBROUTINE allocate_fm_matrix(matrix, nrow, ncol, blacs_env)
      TYPE(cp_fm_type), POINTER, INTENT(OUT)             :: matrix
      INTEGER, INTENT(IN)                                :: nrow, ncol
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'allocate_fm_matrix'

      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=nrow, &
                               ncol_global=ncol, &
                               context=blacs_env)
      ALLOCATE(matrix)
      CALL cp_fm_create(matrix, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(matrix, 0.0_dp)

      CALL timestop(handle)

   END SUBROUTINE allocate_fm_matrix

   SUBROUTINE deallocate_fm_matrix(matrix)
      TYPE(cp_fm_type), POINTER, INTENT(INOUT)           :: matrix

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'deallocate_fm_matrix'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(matrix)) THEN
         CALL cp_fm_release(matrix)
         DEALLOCATE(matrix)
      END IF

      CALL timestop(handle)

   END SUBROUTINE deallocate_fm_matrix

   SUBROUTINE tddft_initialization(tddft_env, qs_env, gs_mos, blacs_env)
      CLASS(tddft_info_type), INTENT(INOUT)              :: tddft_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_allocate_data'

      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(section_vals_type), POINTER                   :: hf_section
      INTEGER                                            :: ivec, ispin, handle
      INTEGER :: n_rep_hf, potential_type
      LOGICAL :: explicit

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s)

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

      If (tddft_env%do_hfx_cam) THEN
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

      ! Set the reference type ROHF/UHF
      tddft_env%reference_type = 0
      IF (dft_control%roks) tddft_env%reference_type = 1
      IF (dft_control%uks) tddft_env%reference_type = 2

      tddft_env%nStates = dft_control%tddfpt2_control%nStates
      tddft_env%nSpins = SIZE(gs_mos, 1)
      tddft_env%nStates_max = tddft_env%nocca * tddft_env%nvirb
      tddft_env%maxvec = dft_control%tddfpt2_control%nkvs
      tddft_env%mxvec = min(tddft_env%maxvec*tddft_env%nstates, tddft_env%nStates_max)
      tddft_env%nStates = min(tddft_env%nstates, tddft_env%mxvec)
      tddft_env%nvec = min(max(tddft_env%nStates, 4), tddft_env%mxvec)
      tddft_env%novec = tddft_env%nvec

      IF (tddft_env%nStates > tddft_env%nStates_max) THEN
         CPWARN("NUMBER OF EXCITED STATES COULD LEAD TO PROBLEMS!")
         CPWARN("Experimental: CHANGED nStates TO ITS MAXIMUM VALUE!")
         tddft_env%nStates = tddft_env%nStates_max
      END IF

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
      !=== focka_noc_noc [nocca, nocca]
      CALL allocate_fm_matrix(matrix=tddft_env%focka_noc_noc, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nocca, &
                              blacs_env=blacs_env)
      !=== fockb_vir_vir [nvirb, nvirb]
      CALL allocate_fm_matrix(matrix=tddft_env%fockb_vir_vir, &
                              nrow=tddft_env%nvirb, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      !--------------------
      ! ALLOCATE  FM POOL
      !--------------------
      !---------- MO POOL [nocca, nvirb] : fm_pool_mo
      ALLOCATE (tddft_env%fm_pool_mo)
      NULLIFY (tddft_env%fm_pool_mo%pool)
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nocca, &
                               ncol_global=tddft_env%nvirb, &
                               context=blacs_env)
      CALL fm_pool_create(tddft_env%fm_pool_mo%pool, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      !---------- END MO POOL [nocca, nvirb]

      ! MO POOL: evects_mo (mxvec)[nocca, nvirb]
      ALLOCATE (tddft_env%evects_mo(tddft_env%mxvec))
      DO ivec = 1, tddft_env%mxvec
         CALL fm_pool_create_fm(tddft_env%fm_pool_mo%pool, &
                                tddft_env%evects_mo(ivec))
      END DO
      ! MO POOL: ritz_vects (nstates)[nocca, nvirb]
      ALLOCATE (tddft_env%ritz_vects(tddft_env%nstates))
      DO ivec = 1, tddft_env%nstates
         CALL fm_pool_create_fm(tddft_env%fm_pool_mo%pool, &
                                tddft_env%ritz_vects(ivec))
      END DO
      ! MO POOL: sigma_vecs (mxvec)[nocca, nvirb]
      ALLOCATE (tddft_env%sigma_vecs(tddft_env%mxvec))
      DO ivec = 1, tddft_env%mxvec
         CALL fm_pool_create_fm(tddft_env%fm_pool_mo%pool, &
                                tddft_env%sigma_vecs(ivec))
      END DO

      ! MO POOL: ritz [nocca, nvirb]
      ALLOCATE(tddft_env%residual_vecs(tddft_env%nstates))
      DO ivec = 1, tddft_env%nstates
         CALL fm_pool_create_fm(tddft_env%fm_pool_mo%pool, tddft_env%residual_vecs(ivec))
      END DO

      !---------- AO POOL [nbf, nbf] : fm_pool_ao
      ALLOCATE (tddft_env%fm_pool_ao(tddft_env%nSpins))
      DO ispin = 1, tddft_env%nSpins
         NULLIFY (tddft_env%fm_pool_ao(ispin)%pool)

         CALL cp_fm_struct_create(fm_struct, &
                                  nrow_global=tddft_env%nbf, &
                                  ncol_global=tddft_env%nbf, &
                                  context=blacs_env)
         CALL fm_pool_create(tddft_env%fm_pool_ao(ispin)%pool, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END DO
      !---------- END AO POOL [nbf, nbf]

      ! AO POOL: Fock (nSpins)[nbf, nbf]
      ALLOCATE (tddft_env%fock(tddft_env%nSpins))
      DO ispin = 1, tddft_env%nSpins
         CALL fm_pool_create_fm(tddft_env%fm_pool_ao(ispin)%pool, &
                                tddft_env%fock(ispin))
      END DO
      ! AO POOL: Fock (nSpins)[nbf, nbf]
      ALLOCATE (tddft_env%work_ao_ao(tddft_env%nSpins))
      DO ispin = 1, tddft_env%nSpins
         CALL fm_pool_create_fm(tddft_env%fm_pool_ao(ispin)%pool, &
                                tddft_env%work_ao_ao(ispin))
      END DO

      ALLOCATE (tddft_env%dens_diff_ij(tddft_env%nStates))
      DO ivec = 1, tddft_env%nStates
         CALL fm_pool_create_fm(tddft_env%fm_pool_ao(1)%pool, &
                                tddft_env%dens_diff_ij(ivec))
      END DO
      ALLOCATE (tddft_env%dens_diff_ab(tddft_env%nStates))
      DO ivec = 1, tddft_env%nStates
         CALL fm_pool_create_fm(tddft_env%fm_pool_ao(1)%pool, &
                                tddft_env%dens_diff_ab(ivec))
      END DO
      CALL cp_fm_struct_create(fm_struct, &
                            nrow_global=tddft_env%nbf, &
                            ncol_global=tddft_env%nocca, &
                            context=blacs_env)
      ALLOCATE(tddft_env%evects(1, tddft_env%nvec))
      DO ivec = 1, tddft_env%nvec
         CALL cp_fm_create(matrix=tddft_env%evects(1, ivec), &
                           matrix_struct=fm_struct)
      END DO
      CALL cp_fm_struct_release(fm_struct)

      !--------------------
      ! ALLOCATE  RHO
      !--------------------

      ! Allocate and set up qs_rho_type
      !
      CALL timestop(handle)

   END SUBROUTINE tddft_initialization

   SUBROUTINE tddft_finalization(this, blacs_env)
      CLASS(tddft_info_type), INTENT(INOUT)              :: this
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_deallocate_data'

      INTEGER                                            :: ivec, ispin, handle

      CALL timeset(routineN, handle)

      blacs_env = blacs_env
      !--------------------
      ! DEALLOCATE  DP
      !--------------------
      IF (ALLOCATED(this%e_virt_minus_occ)) THEN
         DEALLOCATE(this%e_virt_minus_occ)
      END IF
      IF (ALLOCATED(this%evals)) THEN
         DEALLOCATE(this%evals)
      END IF
      IF (ALLOCATED(this%ritz_evals)) THEN
         DEALLOCATE(this%ritz_evals)
      END IF
      IF (ALLOCATED(this%ritz_evals_prev)) THEN
         DEALLOCATE(this%ritz_evals_prev)
      END IF
      IF (ALLOCATED(this%vector_scal)) THEN
         DEALLOCATE (this%vector_scal)
      END IF
      IF (ALLOCATED(this%vector_norm)) THEN
         DEALLOCATE (this%vector_norm)
      END IF
      IF (ALLOCATED(this%residual_norms)) THEN
         DEALLOCATE (this%residual_norms)
      END IF
      !--------------------
      ! DEALLOCATE  FM
      !--------------------
      IF (ASSOCIATED(this%focka_noc_noc)) THEN
         CALL cp_fm_release(this%focka_noc_noc)
         DEALLOCATE(this%focka_noc_noc)
      END IF
      IF (ASSOCIATED(this%fockb_vir_vir)) THEN
         CALL cp_fm_release(this%fockb_vir_vir)
         DEALLOCATE(this%fockb_vir_vir)
      END IF
      !--------------------
      ! DEALLOCATE  FM POOL
      !--------------------
      IF (ALLOCATED(this%fock)) THEN
         DO ispin = 1, SIZE(this%fock, 1)
            CALL fm_pool_give_back_fm(this%fm_pool_ao(ispin)%pool, &
                                      this%fock(ispin))
         END DO
         DEALLOCATE(this%fock)
      END IF
      IF (ALLOCATED(this%work_ao_ao)) THEN
         DO ispin = 1, SIZE(this%work_ao_ao, 1)
            CALL fm_pool_give_back_fm(this%fm_pool_ao(ispin)%pool, &
                                      this%work_ao_ao(ispin))
         END DO
         DEALLOCATE(this%work_ao_ao)
      END IF
      IF (ALLOCATED(this%dens_diff_ab)) THEN
         DO ivec = 1, SIZE(this%dens_diff_ab, 1)
            CALL fm_pool_give_back_fm(this%fm_pool_ao(1)%pool, &
                                      this%dens_diff_ab(ivec))
         END DO
         DEALLOCATE(this%dens_diff_ab)
      END IF
      IF (ALLOCATED(this%dens_diff_ij)) THEN
         DO ivec = 1, SIZE(this%dens_diff_ij, 1)
            CALL fm_pool_give_back_fm(this%fm_pool_ao(1)%pool, &
                                      this%dens_diff_ij(ivec))
         END DO
         DEALLOCATE(this%dens_diff_ij)
      END IF
      IF (ALLOCATED(this%evects_mo)) THEN
         DO ivec = 1, SIZE(this%evects_mo, 1)
            CALL fm_pool_give_back_fm(this%fm_pool_mo%pool, &
                                      this%evects_mo(ivec))
         END DO
         DEALLOCATE(this%evects_mo)
      END IF
      IF (ALLOCATED(this%ritz_vects)) THEN
         DO ivec = 1, SIZE(this%ritz_vects, 1)
            CALL fm_pool_give_back_fm(this%fm_pool_mo%pool, &
                                      this%ritz_vects(ivec))
         END DO
         DEALLOCATE(this%ritz_vects)
      END IF
      IF (ALLOCATED(this%sigma_vecs)) THEN
         DO ivec = 1, SIZE(this%sigma_vecs, 1)
            CALL fm_pool_give_back_fm(this%fm_pool_mo%pool, &
                                      this%sigma_vecs(ivec))
         END DO
         DEALLOCATE(this%sigma_vecs)
      END IF

      IF (ALLOCATED(this%residual_vecs))  THEN
         DO ivec = SIZE(this%residual_vecs),1,-1
            CALL fm_pool_give_back_fm(this%fm_pool_mo%pool, this%residual_vecs(ivec))
         END DO
         DEALLOCATE(this%residual_vecs)
      END IF

      IF (ASSOCIATED(this%fm_pool_mo)) THEN
         CALL fm_pool_release(this%fm_pool_mo%pool)
         DEALLOCATE(this%fm_pool_mo)
      END IF

      IF (ASSOCIATED(this%fm_pool_ao)) THEN
         DO ispin = SIZE(this%fm_pool_ao),1,-1
            CALL fm_pool_release(this%fm_pool_ao(ispin)%pool)
         END DO
         DEALLOCATE(this%fm_pool_ao)
      END IF

      IF (ALLOCATED(this%dipole_op_mos_occ)) THEN
         DO ivec = 1, SIZE(this%dipole_op_mos_occ, 1)
            DO ispin = 1, SIZE(this%dipole_op_mos_occ, 2)
               CALL cp_fm_release(this%dipole_op_mos_occ(ivec, ispin))
            END DO
         END DO
         DEALLOCATE(this%dipole_op_mos_occ)
      END IF
      IF (ALLOCATED(this%evects)) THEN
         DO ivec = 1, SIZE(this%evects, 2)
            CALL cp_fm_release(this%evects(1, ivec))
         END DO
         DEALLOCATE(this%evects)
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddft_finalization

END MODULE qs_tddft_types
