!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddft_types
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
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
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
   PUBLIC :: tddft_type
   PUBLIC :: tddft_davidson_type
   PUBLIC :: tddft_allocate_data
   PUBLIC :: tddft_deallocate_data
   PUBLIC :: allocate_fm_matrix
   PUBLIC :: deallocate_fm_matrix

! **************************************************************************************************
!> \brief Universal TDDFT/Davidson infrastructure for excitation solvers.
!>        Contains all state vectors, residuals, subspace matrices, eigenvalues, etc.
!> \author  Konstantin Komarov (constlike@gmail.com)
!> \date    2025-05-18
! **************************************************************************************************
   TYPE tddft_davidson_type

      !> Number of states to solve for
      INTEGER :: nstates
      !> Current number of vectors in subspace (Krylov dimension)
      INTEGER :: nvec
      !> Max allowed number of vectors in subspace
      INTEGER :: mxvec

      !> Number of basis functions (AO)
      INTEGER :: nbf
      !> Number of alpha/beta occupied, virtual orbitals
      INTEGER :: nocca, noccb, nvirb, nvira
      !> Number of spin channels
      INTEGER :: nspins

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
      TYPE(cp_fm_type) :: subspace
      !> Subspace eigenvectors (coefficients for Ritz vectors)
      TYPE(cp_fm_type) :: subspace_evecs

      !> Other working buffers for orthogonalization, normalization, etc.
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: vector_norm
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: vector_scal

      !> Diagonal MO energy differences (universal for guess, etc)
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: e_virt_minus_occ

   END TYPE tddft_davidson_type

! **************************************************************************************************
!> \brief General TDDFT buffers and configuration
!>        Holds only what is needed for a specific TDDFT theory (SF, MRSF, ...)
!>        All data here is private for theory_handler%fields.
!> \author  Konstantin Komarov (constlike@gmail.com)
!> \date    2025-05-18
! **************************************************************************************************
   TYPE tddft_type

      !> Reference type: 1=ROKS/ROHF, 2=UKS/UHF
      INTEGER :: reference_type

      !> [nbf, nbf] intermediate work matrices (for sym/antisym density etc)
      TYPE(cp_fm_type) :: work1_ao_ao
      TYPE(cp_fm_type) :: work2_ao_ao

      !> [nbf, nbf] AO trial vector in AO form (symmetric, antisymmetric components)
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: evects_ao_ao

      !> [nbf, nbf] DBCSR work matrices for density-like, symmetric/antisymmetric matrices
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:) :: work_rho_ao_ao_symm
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:) :: work_rho_ao_ao_asymm

      !> DBCSR/CP_FM working matrices for HFX, AO sigma
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:) :: work_hmat_symm, work_hmat_asymm

      !> Buffers for AO sigma, AO-to-MO transformation
      TYPE(cp_fm_type) :: hmat_ao_ao
      TYPE(cp_fm_type) :: work_ao_vir
      TYPE(cp_fm_type) :: work_noc_ao
      TYPE(cp_fm_type) :: work_noc_vir

      ! Buffers for alpha and beta Fock(1:2) and Fock submatrices in MO basis
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fock
      TYPE(cp_fm_type), POINTER                          :: Focka_noc_noc => NULL()
      TYPE(cp_fm_type), POINTER                          :: Fockb_vir_vir => NULL()

      !> AO density and grid representations
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:) :: rho_ia_ao
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:) :: A_ia_munu_sub
      TYPE(pw_c1d_gs_type), ALLOCATABLE, DIMENSION(:) :: A_ia_rspace_sub

   END TYPE tddft_type

   TYPE tddft_info_type
      INTEGER                                            :: reference_type = 0 ! 1 = ROKS/ROHF, 2 = UKS/UHF
      INTEGER                                            :: nSpins
      INTEGER                                            :: nbf
      INTEGER                                            :: nocca
      INTEGER                                            :: noccb
      INTEGER                                            :: nvira
      INTEGER                                            :: nvirb
      INTEGER                                            :: nStates
      INTEGER                                            :: nStates_max
      INTEGER                                            :: nvec
      INTEGER                                            :: maxvec
      INTEGER                                            :: mxvec
      INTEGER                                            :: xvec_dim

      !> Energy difference E_virt - E_occ, used in initial guess and post-processing
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_minus_occ
      !> TDDFT excitation eigenvalues (e.g., omega_i)
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals

      ! Fock
      TYPE(cp_fm_type), POINTER                          :: Focka_noc_noc => NULL()
      TYPE(cp_fm_type), POINTER                          :: Fockb_vir_vir => NULL()
      TYPE(cp_fm_type), POINTER                          :: work_ao_vir => NULL()
      TYPE(cp_fm_type), POINTER                          :: work_noc_ao => NULL()

      ! POOL MO [nocca, nvirb]
      TYPE(cp_fm_pool_p_type), POINTER                   :: fm_pool_mo => NULL()
      ! Work
      TYPE(cp_fm_type), POINTER                          :: work_noc_vir => NULL()

      ! evects_mo (mxvec)
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: evects_mo
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: ritz_vects
      ! sigma_vecs (mxvec)
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: sigma_vecs

      ! POOL AO [nbf, nbf]
      TYPE(cp_fm_pool_p_type), POINTER, DIMENSION(:)     :: fm_pool_ao => NULL()
      ! Work
      TYPE(cp_fm_type), POINTER, DIMENSION(:)            :: work_ao_ao => NULL()
      TYPE(cp_fm_type), POINTER                          :: hmat_ao_ao => NULL()
      TYPE(cp_fm_type), POINTER                          :: work1_ao_ao => NULL()
      TYPE(cp_fm_type), POINTER                          :: work2_ao_ao => NULL()
      ! Fock: AO to MO (nSpins)
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fock
      ! evects_ao_ao (nvec)
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: evects_ao_ao
      ! Dypoles dummy
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:,:)      :: dipole_op_mos_occ

      ! DBCSR (nvec)
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: work_rho_ao_ao_symm
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: work_rho_ao_ao_asymm
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: work_hmat_symm
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: work_hmat_asymm
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: A_ia_munu_sub
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: rho_ia_ao

      ! PW
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: A_ia_rspace_sub => NULL()
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_r => NULL()
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_g => NULL()
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: tau_r => NULL()
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: work_v_xc => NULL()
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: work_v_xc_tau => NULL()

      ! RHO
      TYPE(qs_rho_type), POINTER                         :: rho_orb_struct_sub => NULL()

      !> Ground-state density matrix in AO basis (shared across subgroups)
      TYPE(qs_rho_type), POINTER                         :: gs_rho_xc_struct => NULL()
      TYPE(cp_fm_type), POINTER                          :: rho_ia_ao_fm => NULL()

      ! ---------- FOR DAVIDSON ----------
      TYPE(cp_fm_type), POINTER                          :: subspace => NULL()
      TYPE(cp_fm_type), POINTER                          :: subspace_evecs => NULL()
      TYPE(cp_fm_type), POINTER                          :: transformation_matrix => NULL()
      REAL(dp)        , ALLOCATABLE , DIMENSION(:)  :: vector_norm
      REAL(dp)        , ALLOCATABLE , DIMENSION(:)  :: vector_scal
      TYPE(cp_fm_type) , ALLOCATABLE , DIMENSION(:) :: residual_vecs
      REAL(kind=dp)    , ALLOCATABLE , DIMENSION(:) :: residual_norms
      REAL(dp)        , ALLOCATABLE , DIMENSION(:)  :: ritz_evals
      REAL(dp)        , ALLOCATABLE , DIMENSION(:)  :: ritz_evals_prev
      INTEGER                                       :: novec = 0          ! sub-space size
   END TYPE tddft_info_type

CONTAINS

   SUBROUTINE allocate_fm_matrix(matrix, nrow, ncol, sub_env)
      TYPE(cp_fm_type), POINTER, INTENT(OUT)             :: matrix
      INTEGER, INTENT(IN)                                :: nrow, ncol
      TYPE(tddfpt_subgroup_env_type), INTENT(IN)         :: sub_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'allocate_fm_matrix'

      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=nrow, &
                               ncol_global=ncol, &
                               context=sub_env%blacs_env)
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

   SUBROUTINE tddft_allocate_data(tddft_env, qs_env, gs_mos, sub_env)
      CLASS(tddft_info_type), INTENT(INOUT)              :: tddft_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(tddfpt_subgroup_env_type), INTENT(IN)         :: sub_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_allocate_data'

      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: dbcsr_template_hfx
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_hfx
      INTEGER                                            :: ivec, ispin, handle
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s)

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
      tddft_env%xvec_dim = tddft_env%nStates_max
      tddft_env%mxvec = min(tddft_env%maxvec*tddft_env%nstates, tddft_env%xvec_dim)
      tddft_env%nStates = min(tddft_env%nstates, tddft_env%mxvec)
      tddft_env%nvec = min(max(2*tddft_env%nStates, 5), tddft_env%mxvec)
      tddft_env%novec = tddft_env%nvec

      IF (tddft_env%nStates > tddft_env%nStates_max) THEN
         CPWARN("NUMBER OF EXCITED STATES COULD LEAD TO PROBLEMS!")
         CPWARN("Experimental: CHANGED nStates TO ITS MAXIMUM VALUE!")
         tddft_env%nStates = tddft_env%nStates_max
      END IF

      !--------------------
      ! ALLOCATE  DP
      !--------------------
      ALLOCATE (tddft_env%e_virt_minus_occ(tddft_env%xvec_dim), source=0.0_dp)
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
                              sub_env=sub_env)
      !=== fockb_vir_vir [nvirb, nvirb]
      CALL allocate_fm_matrix(matrix=tddft_env%fockb_vir_vir, &
                              nrow=tddft_env%nvirb, &
                              ncol=tddft_env%nvirb, &
                              sub_env=sub_env)
      !=== work_noc_ao [nocca, nbf]
      CALL allocate_fm_matrix(matrix=tddft_env%work_noc_ao, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nbf, &
                              sub_env=sub_env)
      !=== work_ao_vir [nbf, nvirb]
      CALL allocate_fm_matrix(matrix=tddft_env%work_ao_vir, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nvirb, &
                              sub_env=sub_env)
      !=== rho_ia_ao_fm [nbf, nbf]
      CALL allocate_fm_matrix(matrix=tddft_env%rho_ia_ao_fm, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nbf, &
                              sub_env=sub_env)

      ALLOCATE(tddft_env%dipole_op_mos_occ(tddft_env%nstates, tddft_env%nSpins))
      DO ivec = 1, tddft_env%nstates
         DO ispin = 1, tddft_env%nSpins
            CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, matrix_struct=fm_struct)
            CALL cp_fm_create(tddft_env%dipole_op_mos_occ(ivec, ispin), fm_struct)
            CALL cp_fm_set_all(tddft_env%dipole_op_mos_occ(ivec, ispin), 0.0_dp)
            NULLIFY (fm_struct)
         END DO
      END DO
      !--------------------
      ! ALLOCATE  FM POOL
      !--------------------
      !---------- MO POOL [nocca, nvirb] : fm_pool_mo
      ALLOCATE (tddft_env%fm_pool_mo)
      NULLIFY (tddft_env%fm_pool_mo%pool)
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nocca, &
                               ncol_global=tddft_env%nvirb, &
                               context=sub_env%blacs_env)
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
      ! MO POOL: work_noc_vir [nocca, nvirb]
      ALLOCATE(tddft_env%work_noc_vir)
      CALL fm_pool_create_fm(tddft_env%fm_pool_mo%pool, &
                             tddft_env%work_noc_vir)

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
                                  context=sub_env%blacs_env)
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
      ! AO POOL: evects_ao_ao (nspin)
      ALLOCATE (tddft_env%evects_ao_ao(1))
      DO ispin = 1, 1
         CALL fm_pool_create_fm(tddft_env%fm_pool_ao(ispin)%pool, &
                                tddft_env%evects_ao_ao(ispin))
      END DO
      ! AO POOL: Fock (nSpins)[nbf, nbf]
      ALLOCATE (tddft_env%work_ao_ao(tddft_env%nSpins))
      DO ispin = 1, tddft_env%nSpins
         CALL fm_pool_create_fm(tddft_env%fm_pool_ao(ispin)%pool, &
                                tddft_env%work_ao_ao(ispin))
      END DO

      ! AO POOL: Work
      ALLOCATE(tddft_env%work1_ao_ao)
      CALL fm_pool_create_fm(tddft_env%fm_pool_ao(1)%pool, tddft_env%work1_ao_ao)
      ALLOCATE(tddft_env%work2_ao_ao)
      CALL fm_pool_create_fm(tddft_env%fm_pool_ao(1)%pool, tddft_env%work2_ao_ao)
      ALLOCATE (tddft_env%hmat_ao_ao)
      CALL fm_pool_create_fm(tddft_env%fm_pool_ao(1)%pool, tddft_env%hmat_ao_ao)

      !--------------------
      ! ALLOCATE  DBCSR
      !--------------------
      CALL get_qs_env(qs_env, dbcsr_dist=sub_env%dbcsr_dist, sab_orb=sab_hfx)
      dbcsr_template_hfx => matrix_s(1)%matrix

      NULLIFY (tddft_env%work_rho_ao_ao_symm)
      CALL dbcsr_allocate_matrix_set(tddft_env%work_rho_ao_ao_symm, tddft_env%nSpins)
      NULLIFY (tddft_env%work_rho_ao_ao_asymm)
      CALL dbcsr_allocate_matrix_set(tddft_env%work_rho_ao_ao_asymm, tddft_env%nSpins)
      DO ispin = 1, tddft_env%nSpins
         CALL dbcsr_init_p(tddft_env%work_rho_ao_ao_symm(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist( &
                  tddft_env%work_rho_ao_ao_symm(ispin)%matrix, &
                  template=dbcsr_template_hfx, &
                  dbcsr_dist=sub_env%dbcsr_dist, &
                  sab=sab_hfx)

         CALL dbcsr_init_p(tddft_env%work_rho_ao_ao_asymm(ispin)%matrix)
         CALL dbcsr_create(tddft_env%work_rho_ao_ao_asymm(ispin)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric, &
                           template=tddft_env%work_rho_ao_ao_symm(ispin)%matrix)
         CALL dbcsr_complete_redistribute(tddft_env%work_rho_ao_ao_symm(ispin)%matrix, &
                                          tddft_env%work_rho_ao_ao_asymm(ispin)%matrix)
      END DO

      CALL dbcsr_allocate_matrix_set(tddft_env%work_hmat_symm, tddft_env%nSpins)
      CALL dbcsr_allocate_matrix_set(tddft_env%work_hmat_asymm, tddft_env%nSpins)
      DO ispin = 1, tddft_env%nSpins
         CALL dbcsr_init_p(tddft_env%work_hmat_symm(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(tddft_env%work_hmat_symm(ispin)%matrix, &
                                          template=dbcsr_template_hfx, &
                                          dbcsr_dist=sub_env%dbcsr_dist, &
                                          sab=sab_hfx)

         CALL dbcsr_init_p(tddft_env%work_hmat_asymm(ispin)%matrix)
         CALL dbcsr_create(tddft_env%work_hmat_asymm(ispin)%matrix, matrix_type=dbcsr_type_antisymmetric, &
                           template=tddft_env%work_hmat_symm(ispin)%matrix)
         CALL dbcsr_complete_redistribute(tddft_env%work_hmat_symm(ispin)%matrix, &
                                          tddft_env%work_hmat_asymm(ispin)%matrix)
      END DO

      NULLIFY (tddft_env%A_ia_munu_sub)
      CALL dbcsr_allocate_matrix_set(tddft_env%A_ia_munu_sub, 1)
      DO ispin = 1, 1
         CALL dbcsr_init_p(tddft_env%A_ia_munu_sub(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(tddft_env%A_ia_munu_sub(ispin)%matrix, &
                                          template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, &
                                          sab=sub_env%sab_orb)
      END DO

      !--------------------
      ! ALLOCATE  RHO
      !--------------------
      NULLIFY (tddft_env%rho_ia_ao)
      CALL dbcsr_allocate_matrix_set(tddft_env%rho_ia_ao, tddft_env%nSpins)
      DO ispin = 1, tddft_env%nSpins
         CALL dbcsr_init_p(tddft_env%rho_ia_ao(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(matrix=tddft_env%rho_ia_ao(ispin)%matrix, &
                                          template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, &
                                          sab=sub_env%sab_orb)
      END DO

      ! Allocate and set up qs_rho_type
      NULLIFY (tddft_env%rho_orb_struct_sub)
      ALLOCATE (tddft_env%rho_orb_struct_sub)
      CALL qs_rho_create(tddft_env%rho_orb_struct_sub)
      CALL qs_rho_set(tddft_env%rho_orb_struct_sub, rho_ao=tddft_env%rho_ia_ao)
      CALL qs_rho_rebuild(tddft_env%rho_orb_struct_sub, qs_env, &
                          rebuild_ao=.FALSE., &
                          rebuild_grids=.TRUE., &
                          pw_env_external=sub_env%pw_env)

      !--------------------
      ! ALLOCATE  PW
      !--------------------
      CALL pw_env_get(sub_env%pw_env, auxbas_pw_pool=auxbas_pw_pool)
      ALLOCATE (tddft_env%A_ia_rspace_sub(tddft_env%nSpins))
      DO ispin = 1, tddft_env%nspins
         CALL auxbas_pw_pool%create_pw(tddft_env%A_ia_rspace_sub(ispin))
      END DO
      ALLOCATE (tddft_env%rho_r(tddft_env%nspins))
      DO ispin = 1, tddft_env%nspins
         CALL auxbas_pw_pool%create_pw(tddft_env%rho_r(ispin))
      END DO
      ALLOCATE (tddft_env%rho_g(tddft_env%nspins))
      DO ispin = 1, tddft_env%nspins
         CALL auxbas_pw_pool%create_pw(tddft_env%rho_g(ispin))
      END DO
      ALLOCATE (tddft_env%tau_r(tddft_env%nspins))
      DO ispin = 1, tddft_env%nspins
         CALL auxbas_pw_pool%create_pw(tddft_env%tau_r(ispin))
      END DO
      ALLOCATE (tddft_env%work_v_xc(tddft_env%nspins))
      DO ispin = 1, tddft_env%nspins
         CALL auxbas_pw_pool%create_pw(tddft_env%work_v_xc(ispin))
      END DO
      ALLOCATE (tddft_env%work_v_xc_tau(tddft_env%nspins))
      DO ispin = 1, tddft_env%nspins
         CALL auxbas_pw_pool%create_pw(tddft_env%work_v_xc_tau(ispin))
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddft_allocate_data

   SUBROUTINE tddft_deallocate_data(this, sub_env, kernel_env)
      CLASS(tddft_info_type), INTENT(INOUT)              :: this
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(kernel_env_type)                              :: kernel_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_deallocate_data'

      INTEGER                                            :: ivec, ispin, handle
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

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
      IF (ASSOCIATED(this%work_noc_ao)) THEN
         CALL cp_fm_release(this%work_noc_ao)
         DEALLOCATE(this%work_noc_ao)
      END IF
      IF (ASSOCIATED(this%work_ao_vir)) THEN
         CALL cp_fm_release(this%work_ao_vir)
         DEALLOCATE(this%work_ao_vir)
      END IF
      IF (ALLOCATED(this%dipole_op_mos_occ)) THEN
         DO ivec = 1, this%nstates
            DO ispin = 1, this%nSpins
               CALL cp_fm_release(this%dipole_op_mos_occ(ivec, ispin))
            END DO
         END DO
         DEALLOCATE(this%dipole_op_mos_occ)
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
      IF (ASSOCIATED(this%work_noc_vir%matrix_struct)) THEN
         CALL fm_pool_give_back_fm(this%fm_pool_mo%pool, &
                                   this%work_noc_vir)
         DEALLOCATE(this%work_noc_vir)
      END IF

      IF (ALLOCATED(this%evects_ao_ao)) THEN
         DO ivec = 1, SIZE(this%evects_ao_ao, 1)
            CALL fm_pool_give_back_fm(this%fm_pool_mo%pool, &
                                      this%evects_ao_ao(ivec))
         END DO
         DEALLOCATE(this%evects_ao_ao)
      END IF
      IF (ASSOCIATED(this%work1_ao_ao%matrix_struct)) THEN
         CALL fm_pool_give_back_fm(this%fm_pool_mo%pool, &
                                   this%work1_ao_ao)
         DEALLOCATE(this%work1_ao_ao)
      END IF
      IF (ASSOCIATED(this%work2_ao_ao%matrix_struct)) THEN
         CALL fm_pool_give_back_fm(this%fm_pool_mo%pool, &
                                   this%work2_ao_ao)
         DEALLOCATE(this%work2_ao_ao)
      END IF
      IF (ASSOCIATED(this%hmat_ao_ao%matrix_struct)) THEN
         CALL fm_pool_give_back_fm(this%fm_pool_mo%pool, &
                                   this%hmat_ao_ao)
         DEALLOCATE(this%hmat_ao_ao)
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

      IF (ASSOCIATED(this%gs_rho_xc_struct)) THEN
         CALL qs_rho_release(this%gs_rho_xc_struct)
         DEALLOCATE(this%gs_rho_xc_struct)
      END IF

      IF (ASSOCIATED(this%work_rho_ao_ao_symm)) THEN
         DO ivec = 1, SIZE(this%work_rho_ao_ao_symm)
            CALL dbcsr_deallocate_matrix(this%work_rho_ao_ao_symm(ivec)%matrix)
         END DO
         DEALLOCATE(this%work_rho_ao_ao_symm)
      END IF
      IF (ASSOCIATED(this%work_rho_ao_ao_asymm)) THEN
         DO ivec = 1, SIZE(this%work_rho_ao_ao_asymm)
            CALL dbcsr_deallocate_matrix(this%work_rho_ao_ao_asymm(ivec)%matrix)
         END DO
         DEALLOCATE(this%work_rho_ao_ao_asymm)
      END IF
      IF (ASSOCIATED(this%work_hmat_symm)) THEN
         DO ivec = 1, SIZE(this%work_hmat_symm)
            CALL dbcsr_deallocate_matrix(this%work_hmat_symm(ivec)%matrix)
         END DO
         DEALLOCATE(this%work_hmat_symm)
      END IF
      IF (ASSOCIATED(this%work_hmat_asymm)) THEN
         DO ivec = 1, SIZE(this%work_hmat_asymm)
            CALL dbcsr_deallocate_matrix(this%work_hmat_asymm(ivec)%matrix)
         END DO
         DEALLOCATE(this%work_hmat_asymm)
      END IF
      IF (ASSOCIATED(this%rho_ia_ao)) THEN
         DO ivec = 1, SIZE(this%rho_ia_ao)
            CALL dbcsr_deallocate_matrix(this%rho_ia_ao(ivec)%matrix)
         END DO
         DEALLOCATE(this%rho_ia_ao)
      END IF
      IF (ASSOCIATED(this%A_ia_munu_sub)) THEN
         DO ivec = 1, SIZE(this%A_ia_munu_sub)
            CALL dbcsr_deallocate_matrix(this%A_ia_munu_sub(ivec)%matrix)
         END DO
         DEALLOCATE(this%A_ia_munu_sub)
      END IF

      !--------------------
      ! DEALLOCATE  PW
      !--------------------
      CALL pw_env_get(sub_env%pw_env, auxbas_pw_pool=auxbas_pw_pool)
      IF (ASSOCIATED(this%A_ia_rspace_sub)) THEN
         DO ispin = 1, this%nSpins
            CALL auxbas_pw_pool%give_back_pw(this%A_ia_rspace_sub(ispin))
         END DO
         DEALLOCATE(this%A_ia_rspace_sub)
      END IF
      IF (ASSOCIATED(this%rho_r)) THEN
         DO ispin = 1, this%nSpins
            CALL auxbas_pw_pool%give_back_pw(this%rho_r(ispin))
         END DO
         DEALLOCATE(this%rho_r)
      END IF
      IF (ASSOCIATED(this%rho_g)) THEN
         DO ispin = 1, this%nSpins
            CALL auxbas_pw_pool%give_back_pw(this%rho_g(ispin))
         END DO
         DEALLOCATE(this%rho_g)
      END IF
      IF (ASSOCIATED(this%tau_r)) THEN
         DO ispin = 1, this%nSpins
            CALL auxbas_pw_pool%give_back_pw(this%tau_r(ispin))
         END DO
         DEALLOCATE(this%tau_r)
      END IF
      IF (ASSOCIATED(this%work_v_xc)) THEN
         DO ispin = 1, this%nSpins
            CALL auxbas_pw_pool%give_back_pw(this%work_v_xc(ispin))
         END DO
         DEALLOCATE(this%work_v_xc)
      END IF
      IF (ASSOCIATED(this%work_v_xc_tau)) THEN
         DO ispin = 1, this%nSpins
            CALL auxbas_pw_pool%give_back_pw(this%work_v_xc_tau(ispin))
         END DO
         DEALLOCATE(this%work_v_xc_tau)
      END IF

!     IF (ASSOCIATED(this%rho_orb_struct_sub)) THEN
!        CALL qs_rho_release(this%rho_orb_struct_sub)
!        DEALLOCATE(this%rho_orb_struct_sub)
!     END IF

      CALL release_kernel_env(kernel_env%full_kernel)
      NULLIFY (kernel_env%full_kernel)

      CALL timestop(handle)

   END SUBROUTINE tddft_deallocate_data

END MODULE qs_tddft_types
