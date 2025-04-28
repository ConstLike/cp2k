!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddft_methods
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
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
                                              cp_fm_set_element,&
                                              cp_fm_write_formatted
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
   USE machine,                         ONLY: m_flush
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_methods'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: tddft_info_type, davidson_info_type

! **************************************************************************************************
   TYPE tddft_info_type
      INTEGER                                            :: nspins
      INTEGER                                            :: nocca
      INTEGER                                            :: noccb
      INTEGER                                            :: nbf
      INTEGER                                            :: nvecs
      INTEGER                                            :: nstates
      INTEGER                                            :: nstates_max

      !> Energy difference E_virt - E_occ, used in initial guess and post-processing
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: e_virt_minus_occ

      !> TDDFT excitation eigenvalues (e.g., omega_i)
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: evals

      !> MO space trial vectors: shape [nstates], each cp_fm_type of size [nocca x nvirtb]
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:,:)      :: evects_mo

      !> AO space trial vectors: shape [nstates], each cp_fm_type of size [nbf x nbf]
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: evects_ao_ao
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: work_ao_ao_sym
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: work_ao_ao_asym
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fock_mo
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:)      :: work_hmat_symm
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:)      :: work_hmat_asymm
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:)      :: work_rho_ao_ao_symm
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:)      :: work_rho_ao_ao_asymm

      TYPE(cp_fm_type), POINTER                          :: Fock_noc_noc => NULL()
      TYPE(cp_fm_type), POINTER                          :: Fock_vir_vir => NULL()
      TYPE(cp_fm_type), POINTER                          :: work_noc_vir => NULL()

      TYPE(cp_fm_type), POINTER                          :: hmat_ao_ao => NULL()
      TYPE(cp_fm_type), POINTER                          :: fock_ao_ao_alpha => NULL()
      TYPE(cp_fm_type), POINTER                          :: fock_ao_ao_beta => NULL()
      TYPE(cp_fm_type), POINTER                          :: F_alpha_Z_ao_ao => NULL()
      TYPE(cp_fm_type), POINTER                          :: Z_ext_noc_ao => NULL()
      TYPE(cp_fm_type), POINTER                          :: Z_F_beta_noc_ao => NULL()
      TYPE(cp_fm_type), POINTER                          :: work_noc_vir2 => NULL()
      TYPE(cp_fm_type), POINTER                          :: work_noc_ao => NULL()
      TYPE(cp_fm_type), POINTER                          :: work_ao_vir => NULL()

      !> Matrix pools for efficient memory management of MO and AO trial vectors
      TYPE(cp_fm_pool_p_type), POINTER, DIMENSION(:)     :: fm_pool_mo => NULL()
      TYPE(cp_fm_pool_p_type), POINTER, DIMENSION(:)     :: fm_pool_ao => NULL()

      !> Ground-state density matrix in AO basis (shared across subgroups)
      TYPE(qs_rho_type), POINTER                         :: gs_rho_orb_struct => NULL()
      TYPE(qs_rho_type), POINTER                         :: gs_rho_xc_struct => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gs_rho_orb_ao
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gs_rho_xc_ao
      TYPE(cp_fm_type), POINTER                          :: gs_rho_orb_ao_fm => NULL()
      CONTAINS
         PROCEDURE :: init => tddft_initialize_data
         PROCEDURE :: setup => tddft_setup_data
         PROCEDURE :: run => tddft_run
         PROCEDURE :: clean => tddft_finalize_data
   END TYPE tddft_info_type

   TYPE davidson_info_type
      ! Krylov subspace vectors
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :) :: krylov_vects  ! Krylov vectors
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :) :: Aop_krylov   ! A operator applied to Krylov vectors
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :) :: S_krylov     ! S operator applied to Krylov vectors

      ! Ritz vectors and related data
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :) :: ritz_vects      ! Ritz vectors (approximate eigenvectors)
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :) :: Aop_ritz        ! A operator applied to Ritz vectors
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)       :: ritz_evals      ! Ritz values (approximate eigenvalues)
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)       :: ritz_evals_last ! Eigenvalues from the previous iteration

      ! Auxiliary data
      REAL(kind=dp), POINTER, DIMENSION(:, :)        :: Atilde       ! Projected A matrix in Krylov subspace

      ! Algorithm parameters
      INTEGER                                        :: nvects_exist  ! Current number of vectors in subspace
      INTEGER                                        :: nvects_new    ! Number of new vectors added this iteration
      INTEGER                                        :: max_krylov_vects  ! Maximum number of Krylov vectors
      INTEGER                                        :: nstates       ! Number of target states
      INTEGER                                        :: nspins        ! Number of spin channels
      CONTAINS
         PROCEDURE :: init => init_davidson_env
         PROCEDURE :: clean => release_davidson_env
   END TYPE davidson_info_type

! **************************************************************************************************

CONTAINS

   SUBROUTINE init_davidson_env(this, evects, nspins, nstates, &
                                max_krylov_vects, fm_pool)
      CLASS(davidson_info_type)         :: this
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(in) :: evects
      INTEGER, INTENT(in)                           :: nspins, nstates, max_krylov_vects
      TYPE(cp_fm_pool_p_type), DIMENSION(:), INTENT(in) :: fm_pool

      INTEGER :: ispin, istate

      this%nspins = nspins
      ! Set parameters
      this%nstates = nstates
      this%nvects_exist = 0
      this%max_krylov_vects = max_krylov_vects
      this%nvects_new = nstates

      ! Allocate memory for arrays
      ALLOCATE(this%krylov_vects(nspins, max_krylov_vects))
      ALLOCATE(this%Aop_krylov(nspins, max_krylov_vects))
      ALLOCATE(this%S_krylov(nspins, max_krylov_vects))
      ALLOCATE(this%ritz_vects(nspins, nstates))
      ALLOCATE(this%Aop_ritz(nspins, nstates))
      ALLOCATE(this%ritz_evals(nstates))
      ALLOCATE(this%ritz_evals_last(max_krylov_vects))
      ALLOCATE(this%Atilde(1,1))  ! Placeholder for Atilde pointer

      ! Initialize matrices using the matrix pool
      DO istate = 1, nstates
         DO ispin = 1, nspins
            CALL fm_pool_create_fm(fm_pool(ispin)%pool, this%krylov_vects(ispin, istate))
            CALL cp_fm_to_fm(evects(ispin, istate), this%krylov_vects(ispin, istate))
            CALL fm_pool_create_fm(fm_pool(ispin)%pool, this%S_krylov(ispin, istate))
            CALL fm_pool_create_fm(fm_pool(ispin)%pool, this%Aop_krylov(ispin, istate))
            CALL fm_pool_create_fm(fm_pool(ispin)%pool, this%ritz_vects(ispin, istate))
            CALL fm_pool_create_fm(fm_pool(ispin)%pool, this%Aop_ritz(ispin, istate))
         END DO
      END DO
      ! Remaining krylov_vects will be initialized as needed during the algorithm
   END SUBROUTINE init_davidson_env

   SUBROUTINE release_davidson_env(this, fm_pool)
      CLASS(davidson_info_type)         :: this
      TYPE(cp_fm_pool_p_type), DIMENSION(:), INTENT(in) :: fm_pool

      INTEGER :: ispin, istate

      ! Release Krylov subspace vectors
      DO istate = 1, this%max_krylov_vects
         DO ispin = 1, this%nspins
            IF (ASSOCIATED(this%krylov_vects(ispin, istate)%matrix_struct)) THEN
               CALL fm_pool_give_back_fm(fm_pool(ispin)%pool, this%krylov_vects(ispin, istate))
            END IF
            IF (ASSOCIATED(this%Aop_krylov(ispin, istate)%matrix_struct)) THEN
               CALL fm_pool_give_back_fm(fm_pool(ispin)%pool, this%Aop_krylov(ispin, istate))
            END IF
            IF (ASSOCIATED(this%S_krylov(ispin, istate)%matrix_struct)) THEN
               CALL fm_pool_give_back_fm(fm_pool(ispin)%pool, this%S_krylov(ispin, istate))
            END IF
         END DO
      END DO
      DEALLOCATE(this%krylov_vects, this%Aop_krylov, this%S_krylov)

      ! Release Ritz vectors
      DO istate = 1, this%nstates
         DO ispin = 1, this%nspins
            CALL fm_pool_give_back_fm(fm_pool(ispin)%pool, this%ritz_vects(ispin, istate))
            CALL fm_pool_give_back_fm(fm_pool(ispin)%pool, this%Aop_ritz(ispin, istate))
         END DO
      END DO
      DEALLOCATE(this%ritz_vects, this%Aop_ritz)

      ! Release remaining data
      DEALLOCATE(this%ritz_evals, this%ritz_evals_last, this%Atilde)
   END SUBROUTINE release_davidson_env


   ! **************************************************************************************************
   !> \brief Initializes TDDFT matrix containers and internal pools for trial vector storage.
   !> \param this Object of type tddft_info_type
   !> \param qs_env Quickstep environment
   !> \param gs_mos Ground-state molecular orbitals (spin-resolved)
   !> \param nvecs Number of trial vectors to allocate
   ! **************************************************************************************************
   SUBROUTINE tddft_initialize_data(this, qs_env, gs_mos, sub_env)
      CLASS(tddft_info_type), INTENT(INOUT)              :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(tddfpt_subgroup_env_type), INTENT(IN)         :: sub_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_initialize_data'

      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      INTEGER                                            :: nvirtb, ivec, ispin, handle

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s)

      CALL cp_fm_get_info(gs_mos(1)%mos_occ, ncol_global=this%nocca)
      CALL cp_fm_get_info(gs_mos(2)%mos_occ, ncol_global=this%noccb)
      CALL cp_fm_get_info(gs_mos(2)%mos_virt, ncol_global=nvirtb)
      CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=this%nbf)
      this%nvecs = dft_control%tddfpt2_control%nstates
      this%nstates = dft_control%tddfpt2_control%nstates
      this%nspins = SIZE(gs_mos, 1)
      this%nstates_max = this%nocca * nvirtb

      IF (this%nstates > this%nstates_max) THEN
         CPWARN("NUMBER OF EXCITED STATES COULD LEAD TO PROBLEMS!")
         CPWARN("Experimental: CHANGED NSTATES TO ITS MAXIMUM VALUE!")
         this%nstates = this%nstates_max
      END IF

      ! Allocate eigenvalue arrays
      ALLOCATE (this%e_virt_minus_occ(this%nocca * nvirtb))
      ALLOCATE (this%evals(this%nvecs))

      ! Initialize MO matrix pool and allocate MO trial vectors
      ALLOCATE (this%fm_pool_mo(1))
      NULLIFY (this%fm_pool_mo(1)%pool)
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=this%nocca, &
                               ncol_global=nvirtb, &
                               context=sub_env%blacs_env)
      CALL fm_pool_create(this%fm_pool_mo(1)%pool, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ALLOCATE (this%evects_mo(1, this%nvecs))
      DO ivec = 1, this%nvecs
         CALL fm_pool_create_fm(this%fm_pool_mo(1)%pool, &
                                this%evects_mo(1, ivec))
      END DO

      ! Initialize AO matrix pool and allocate AO trial vectors

      ALLOCATE (this%fm_pool_ao(this%nspins))
      DO ispin = 1, this%nspins
         NULLIFY (this%fm_pool_ao(ispin)%pool)
         CALL cp_fm_struct_create(fm_struct, &
                                  nrow_global=this%nbf, &
                                  ncol_global=this%nbf, &
                                  context=sub_env%blacs_env)
         CALL fm_pool_create(this%fm_pool_ao(ispin)%pool, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END DO

      ALLOCATE(this%evects_ao_ao(this%nvecs))
      ALLOCATE(this%work_ao_ao_sym(1))
      ALLOCATE(this%work_ao_ao_asym(1))
      DO ivec = 1, this%nvecs
         CALL fm_pool_create_fm(this%fm_pool_ao(1)%pool, this%evects_ao_ao(ivec))
         CALL fm_pool_create_fm(this%fm_pool_ao(1)%pool, this%work_ao_ao_sym(1))
         CALL fm_pool_create_fm(this%fm_pool_ao(1)%pool, this%work_ao_ao_asym(1))
      END DO

      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=this%nbf, &
                               ncol_global=this%nbf, &
                               context=sub_env%blacs_env)
      ALLOCATE (this%gs_rho_orb_ao_fm)
      CALL cp_fm_create(this%gs_rho_orb_ao_fm, fm_struct)

      ALLOCATE (this%hmat_ao_ao)
      CALL cp_fm_create(this%hmat_ao_ao, fm_struct)

      ALLOCATE (this%fock_ao_ao_alpha)
      CALL cp_fm_create(this%fock_ao_ao_alpha, fm_struct)

      ALLOCATE (this%fock_ao_ao_beta)
      CALL cp_fm_create(this%fock_ao_ao_beta, fm_struct)

      ALLOCATE (this%F_alpha_Z_ao_ao)
      CALL cp_fm_create(this%F_alpha_Z_ao_ao, fm_struct)

      ALLOCATE (this%fock_mo(this%nspins))
      DO ispin = 1, this%nspins
         CALL cp_fm_create(this%fock_mo(ispin), fm_struct)
      END DO

      CALL cp_fm_struct_release(fm_struct)

      ! --- allocate [nocca×nocca] and [nvirtb×nvirtb] ---------------------
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=this%nocca, &
                               ncol_global=this%nocca, &
                               context=sub_env%blacs_env)
      ALLOCATE (this%fock_noc_noc)
      CALL cp_fm_create(this%fock_noc_noc, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=nvirtb, &
                               ncol_global=nvirtb, &
                               context=sub_env%blacs_env)
      ALLOCATE (this%fock_vir_vir)
      CALL cp_fm_create(this%fock_vir_vir, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=this%nbf, &
                               ncol_global=nvirtb, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_ao_vir)
      CALL cp_fm_create(this%work_ao_vir, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=this%nocca, &
                               ncol_global=nvirtb, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_noc_vir)
      CALL cp_fm_create(this%work_noc_vir, fm_struct)
      ALLOCATE(this%work_noc_vir2)
      CALL cp_fm_create(this%work_noc_vir2, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=this%nocca, &
                               ncol_global=this%nbf, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_noc_ao)
      CALL cp_fm_create(this%work_noc_ao, fm_struct)
      ALLOCATE(this%Z_ext_noc_ao)
      CALL cp_fm_create(this%Z_ext_noc_ao, fm_struct)
      ALLOCATE(this%Z_F_beta_noc_ao)
      CALL cp_fm_create(this%Z_F_beta_noc_ao, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ! Allocate DBCSR matrix array
      ALLOCATE(this%gs_rho_orb_ao(this%nspins))
      DO ispin = 1, this%nspins
         CALL dbcsr_init_p(this%gs_rho_orb_ao(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist( &
                  this%gs_rho_orb_ao(ispin)%matrix, &
                  template=matrix_s(1)%matrix, &
                  dbcsr_dist=sub_env%dbcsr_dist, &
                  sab=sub_env%sab_orb)
      END DO

      ! Allocate and set up qs_rho_type
      ALLOCATE(this%gs_rho_orb_struct)
      CALL qs_rho_create(this%gs_rho_orb_struct)
      CALL qs_rho_set(this%gs_rho_orb_struct, rho_ao=this%gs_rho_orb_ao)
      CALL qs_rho_rebuild(this%gs_rho_orb_struct, qs_env, &
                          rebuild_ao=.FALSE., &
                          rebuild_grids=.TRUE., &
                          pw_env_external=sub_env%pw_env)

      IF (dft_control%qs_control%gapw_xc) THEN
         NULLIFY (this%gs_rho_xc_ao)
         CALL dbcsr_allocate_matrix_set(this%gs_rho_xc_ao, this%nspins)
         DO ispin = 1, this%nspins
            CALL dbcsr_init_p(this%gs_rho_xc_ao(ispin)%matrix)
            CALL tddfpt_dbcsr_create_by_dist( &
                     this%gs_rho_xc_ao(ispin)%matrix, &
                     template=matrix_s(1)%matrix, &
                     dbcsr_dist=sub_env%dbcsr_dist, &
                     sab=sub_env%sab_orb)
         END DO
         NULLIFY (this%gs_rho_xc_struct)
         ALLOCATE (this%gs_rho_xc_struct)
         CALL qs_rho_create(this%gs_rho_xc_struct)
         CALL qs_rho_set(this%gs_rho_xc_struct, &
                         rho_ao=this%gs_rho_xc_ao)
         CALL qs_rho_rebuild(this%gs_rho_xc_struct, qs_env, &
                             rebuild_ao=.FALSE., &
                             rebuild_grids=.TRUE., &
                             pw_env_external=sub_env%pw_env)
      END IF

      ALLOCATE(this%work_rho_ao_ao_symm(this%nstates))
      ALLOCATE(this%work_rho_ao_ao_asymm(this%nstates))
      ALLOCATE(this%work_hmat_symm(this%nstates))
      ALLOCATE(this%work_hmat_asymm(this%nstates))

      DO ivec = 1, this%nstates
         CALL dbcsr_init_p(this%work_rho_ao_ao_symm(ivec)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(this%work_rho_ao_ao_symm(ivec)%matrix, &
                                          template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, &
                                          sab=sub_env%sab_orb)
         CALL dbcsr_init_p(this%work_rho_ao_ao_asymm(ivec)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(this%work_rho_ao_ao_asymm(ivec)%matrix, &
                                          template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, &
                                          sab=sub_env%sab_orb)
         CALL dbcsr_init_p(this%work_hmat_symm(ivec)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(this%work_hmat_symm(ivec)%matrix, &
                                          template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, &
                                          sab=sub_env%sab_orb)
         CALL dbcsr_init_p(this%work_hmat_asymm(ivec)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(this%work_hmat_asymm(ivec)%matrix, &
                                          template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, &
                                          sab=sub_env%sab_orb)
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddft_initialize_data

   SUBROUTINE tddft_setup_data(this, qs_env, gs_mos, matrix_ks, logger, xc_section, &
         tddfpt_print_section, tddfpt_control, kernel_env, full_kernel_env, sub_env)
      CLASS(tddft_info_type), INTENT(INOUT)              :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(kernel_env_type)                              :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(full_kernel_env_type), TARGET                 :: full_kernel_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_setup_data'

      INTEGER                                            :: handle, log_unit, ispin
      LOGICAL                                            :: do_exck

      CALL timeset(routineN, handle)

      this%nstates = this%nstates
      CALL tddfpt_construct_ground_state_orb_density( &
               rho_orb_struct=this%gs_rho_orb_struct, &
               rho_xc_struct=this%gs_rho_xc_struct, &
               is_rks_triplets=tddfpt_control%rks_triplets, &
               qs_env=qs_env, &
               sub_env=sub_env, &
               wfm_rho_orb=this%gs_rho_orb_ao_fm)

      ! Full kernel
      CALL create_kernel_env(kernel_env=full_kernel_env, &
                             rho_struct_sub=this%gs_rho_orb_struct, &
                             xc_section=xc_section, &
                             is_rks_triplets=tddfpt_control%rks_triplets, &
                             sub_env=sub_env)
      kernel_env%full_kernel => full_kernel_env
      NULLIFY (kernel_env%admm_kernel)

      ! Fxc from kernel definition
      do_exck = tddfpt_control%do_exck
      kernel_env%full_kernel%do_exck = do_exck
      ! initilize xc kernel
!     IF (do_exck) THEN
!        CALL create_fxc_kernel(this%gs_rho_orb_struct, &
!                               work_matrices%fxc_rspace_sub, &
!                               xc_section, tddfpt_control%rks_triplets, sub_env, qs_env)
!     END IF
!
      ! Transfom matrix_ks from AO to MO basis
      DO ispin = 1, this%nspins
         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, this%fock_mo(ispin))
         ! work_ao_ao_sym = gs_mos(1)%mos_full^T * this%fock_mo
         CALL parallel_gemm('T', 'N', this%nbf, this%nbf, this%nbf, &
                            1.0_dp, gs_mos(ispin)%mos_full, &
                                    this%fock_mo(ispin), &
                            0.0_dp, this%work_ao_ao_sym(1))

         ! this%fock_mo = work_ao_ao_sym * gs_mos(2)%mos_full
         CALL parallel_gemm('N', 'N', this%nbf, this%nbf, this%nbf, &
                            1.0_dp, this%work_ao_ao_sym(1), &
                                    gs_mos(ispin)%mos_full, &
                            0.0_dp, this%fock_mo(ispin))
      END DO

      ! --- fock_occ = F_alpha( occ , occ )  ------------------------------------
      CALL cp_fm_to_fm_submat(this%fock_mo(1), this%fock_noc_noc, &
                              nrow=this%nocca, ncol=this%nocca, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=1, t_firstcol=1)

      ! --- fock_vir = F_beta( vir , vir )  ------------------------------------
      CALL cp_fm_to_fm_submat(this%fock_mo(2), this%fock_vir_vir, &
                              nrow=(this%nbf - this%noccb), &
                              ncol=(this%nbf - this%noccb), &
                              s_firstrow=this%noccb+1, &
                              s_firstcol=this%noccb+1, &
                              t_firstrow=1, t_firstcol=1)

      ! Open logger
      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, &
                                      "GUESS_VECTORS", extension=".tddfptLog")

      CALL tddft_guess_vectors(this%evects_mo, this%evals, this%e_virt_minus_occ, &
                               gs_mos, log_unit, tddfpt_control)
      CALL cp_print_key_finished_output(log_unit, logger, &
                                        tddfpt_print_section, "GUESS_VECTORS")
      ! Close logger

      CALL timestop(handle)

   END SUBROUTINE tddft_setup_data

   ! **************************************************************************************************
   !> \brief Releases all matrices and deallocates memory pools associated with TDDFT trial data.
   !> \param this Object of type tddft_info_type
   ! **************************************************************************************************
   SUBROUTINE tddft_finalize_data(this)
      CLASS(tddft_info_type), INTENT(INOUT)              :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_finalize_data'

      INTEGER                                            :: ivec, ispin, handle

      CALL timeset(routineN, handle)

      IF (ALLOCATED(this%evects_mo)) THEN
         DO ivec = 1, SIZE(this%evects_mo(1,:))
            CALL cp_fm_release(this%evects_mo(1, ivec))
         END DO
         DEALLOCATE(this%evects_mo)
      END IF

      IF (ALLOCATED(this%evects_ao_ao)) THEN
         DO ivec = 1, SIZE(this%evects_ao_ao)
            CALL cp_fm_release(this%evects_ao_ao(ivec))
         END DO
         DEALLOCATE(this%evects_ao_ao)
      END IF
      IF (ALLOCATED(this%work_ao_ao_sym)) THEN
         CALL cp_fm_release(this%work_ao_ao_sym(1))
         DEALLOCATE(this%work_ao_ao_sym)
      END IF
      IF (ALLOCATED(this%work_ao_ao_asym)) THEN
         CALL cp_fm_release(this%work_ao_ao_asym(1))
         DEALLOCATE(this%work_ao_ao_asym)
      END IF
      IF (ASSOCIATED(this%work_ao_vir)) THEN
         CALL cp_fm_release(this%work_ao_vir)
         DEALLOCATE(this%work_ao_vir)
      END IF
      IF (ASSOCIATED(this%work_noc_ao)) THEN
         CALL cp_fm_release(this%work_noc_ao)
         DEALLOCATE(this%work_noc_ao)
      END IF

      IF (ALLOCATED(this%fock_mo)) THEN
         DO ispin = 1, this%nspins
            CALL cp_fm_release(this%fock_mo(ispin))
         END DO
         DEALLOCATE(this%fock_mo)
      END IF

      IF (ASSOCIATED(this%fock_ao_ao_alpha)) THEN
         CALL cp_fm_release(this%fock_ao_ao_alpha)
         DEALLOCATE(this%fock_ao_ao_alpha)
      END IF
      IF (ASSOCIATED(this%fock_ao_ao_beta)) THEN
         CALL cp_fm_release(this%fock_ao_ao_beta)
         DEALLOCATE(this%fock_ao_ao_beta)
      END IF
      IF (ASSOCIATED(this%Z_ext_noc_ao)) THEN
         CALL cp_fm_release(this%Z_ext_noc_ao)
         DEALLOCATE(this%Z_ext_noc_ao)
      END IF
      IF (ASSOCIATED(this%F_alpha_Z_ao_ao)) THEN
         CALL cp_fm_release(this%F_alpha_Z_ao_ao)
         DEALLOCATE(this%F_alpha_Z_ao_ao)
      END IF
      IF (ASSOCIATED(this%Z_F_beta_noc_ao)) THEN
         CALL cp_fm_release(this%Z_F_beta_noc_ao)
         DEALLOCATE(this%Z_F_beta_noc_ao)
      END IF
      IF (ASSOCIATED(this%work_noc_vir)) THEN
         CALL cp_fm_release(this%work_noc_vir)
         DEALLOCATE(this%work_noc_vir)
      END IF
      IF (ASSOCIATED(this%work_noc_vir2)) THEN
         CALL cp_fm_release(this%work_noc_vir2)
         DEALLOCATE(this%work_noc_vir2)
      END IF
      IF (ASSOCIATED(this%fock_noc_noc)) THEN
         CALL cp_fm_release(this%fock_noc_noc)
      END IF
      IF (ASSOCIATED(this%fock_vir_vir)) THEN
         CALL cp_fm_release(this%fock_vir_vir)
      END IF

      IF (ASSOCIATED(this%fm_pool_mo)) THEN
         CALL fm_pool_release(this%fm_pool_mo(1)%pool)
         DEALLOCATE(this%fm_pool_mo)
      END IF

      IF (ASSOCIATED(this%fm_pool_ao)) THEN
         DO ispin = 1, SIZE(this%fm_pool_ao)
            CALL fm_pool_release(this%fm_pool_ao(ispin)%pool)
         END DO
         DEALLOCATE(this%fm_pool_ao)
      END IF

      IF (ALLOCATED(this%e_virt_minus_occ)) &
         DEALLOCATE(this%e_virt_minus_occ)
      IF (ALLOCATED(this%evals))  DEALLOCATE(this%evals)

      IF (ASSOCIATED(this%gs_rho_orb_struct)) THEN
         CALL qs_rho_release(this%gs_rho_orb_struct)
         DEALLOCATE(this%gs_rho_orb_struct)
      END IF

      IF (ASSOCIATED(this%gs_rho_orb_ao)) THEN
         DO ispin = 1, SIZE(this%gs_rho_orb_ao)
            IF (ASSOCIATED(this%gs_rho_orb_ao(ispin)%matrix)) THEN
               CALL dbcsr_release(this%gs_rho_orb_ao(ispin)%matrix)
            END IF
         END DO
         DEALLOCATE(this%gs_rho_orb_ao)
      END IF

      IF (ASSOCIATED(this%gs_rho_orb_ao_fm)) THEN
         CALL cp_fm_release(this%gs_rho_orb_ao_fm)
         DEALLOCATE(this%gs_rho_orb_ao_fm)
      END IF
      IF (ASSOCIATED(this%hmat_ao_ao)) THEN
         CALL cp_fm_release(this%hmat_ao_ao)
         DEALLOCATE(this%hmat_ao_ao)
      END IF

      IF (ASSOCIATED(this%gs_rho_xc_struct)) THEN
         CALL qs_rho_release(this%gs_rho_xc_struct)
         DEALLOCATE(this%gs_rho_xc_struct)
      END IF

      IF (ASSOCIATED(this%gs_rho_xc_ao)) THEN
         DO ispin = 1, SIZE(this%gs_rho_xc_ao)
            IF (ASSOCIATED(this%gs_rho_xc_ao(ispin)%matrix)) THEN
               CALL dbcsr_release(this%gs_rho_xc_ao(ispin)%matrix)
            END IF
         END DO
         DEALLOCATE(this%gs_rho_xc_ao)
      END IF

      IF (ALLOCATED(this%work_rho_ao_ao_symm)) THEN
         DO ivec = 1, SIZE(this%work_rho_ao_ao_symm)
            CALL dbcsr_deallocate_matrix(this%work_rho_ao_ao_symm(ivec)%matrix)
         END DO
         DEALLOCATE(this%work_rho_ao_ao_symm)
      END IF
      IF (ALLOCATED(this%work_rho_ao_ao_asymm)) THEN
         DO ivec = 1, SIZE(this%work_rho_ao_ao_asymm)
            CALL dbcsr_deallocate_matrix(this%work_rho_ao_ao_asymm(ivec)%matrix)
         END DO
         DEALLOCATE(this%work_rho_ao_ao_asymm)
      END IF
      IF (ALLOCATED(this%work_hmat_symm)) THEN
         DO ivec = 1, SIZE(this%work_hmat_symm)
            CALL dbcsr_deallocate_matrix(this%work_hmat_symm(ivec)%matrix)
         END DO
         DEALLOCATE(this%work_hmat_symm)
      END IF
      IF (ALLOCATED(this%work_hmat_asymm)) THEN
         DO ivec = 1, SIZE(this%work_hmat_asymm)
            CALL dbcsr_deallocate_matrix(this%work_hmat_asymm(ivec)%matrix)
         END DO
         DEALLOCATE(this%work_hmat_asymm)
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddft_finalize_data


! **************************************************************************************************
!> \brief The energy calculation has been moved to its own subroutine
!> \param qs_env ...
!> \param nstates ...
!> \param work_matrices ...
!> \param tddfpt_control ...
!> \param logger ...
!> \param tddfpt_print_section ...
!> \param evects ...
!> \param evals ...
!> \param gs_mos ...
!> \param S_evects ...
!> \param matrix_s ...
!> \param kernel_env ...
!> \param matrix_ks ...
!> \param sub_env ...
!> \param xc_section ...
!> \param full_kernel_env ...
!> \param kernel_env_admm_aux ...
! **************************************************************************************************
   SUBROUTINE tddft_run(this, qs_env, tddfpt_control, logger, tddfpt_print_section, &
                        gs_mos, kernel_env, matrix_ks, sub_env)

      CLASS(tddft_info_type), INTENT(INOUT)              :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(kernel_env_type)                              :: kernel_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_run'

      CHARACTER(len=20)                                  :: nstates_str
      INTEGER                                            :: handle, log_unit, &
         energy_unit, iter, niters, nv, iend, ist, iv, ivec, max_krylov_vects, &
         mxvec, nvirtb, i, a, ij
      REAL(KIND=DP)                                      :: conv
      TYPE(davidson_info_type)                           :: davidson
      REAL(dp), POINTER, CONTIGUOUS                      :: evects_data(:,:), aop_data(:,:)

      CALL timeset(routineN, handle)

      gs_mos = gs_mos
      matrix_ks = matrix_ks
      nstates_str = nstates_str
      kernel_env = kernel_env
      logger = logger
      tddfpt_print_section = tddfpt_print_section
      tddfpt_control = tddfpt_control
      qs_env = qs_env
      log_unit = log_unit
      sub_env = sub_env

      niters = tddfpt_control%niters
      IF (niters > 0) THEN
         log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, &
                                         "ITERATION_INFO", extension=".tddfptLog")
         energy_unit = cp_print_key_unit_nr(logger, &
                                            tddfpt_print_section, &
                                            "DETAILED_ENERGY", &
                                            extension=".tddfptLog")

         IF (log_unit > 0) THEN
            WRITE (log_unit, "(1X,A)") "", &
               "-------------------------------------------------------------------------------", &
               "-                     SF-TDDFT WAVEFUNCTION OPTIMIZATION                       -", &
               "-------------------------------------------------------------------------------"

            WRITE (log_unit, '(/,T11,A,T27,A,T40,A,T62,A)') "Step", "Time", "Convergence", "Conv. states"
            WRITE (log_unit, '(1X,79("-"))')
         END IF

         CALL cp_add_iter_level(logger%iter_info, "TDDFT_SCF")

         nvirtb = SIZE(gs_mos(2)%evals_full) - gs_mos(2)%homo
         mxvec = gs_mos(1)%homo * nvirtb
         max_krylov_vects = MIN(MAX(tddfpt_control%nkvs, this%nstates), INT(mxvec))

         CALL davidson%init(evects=this%evects_mo, &
                            nspins=1, &
                            nstates=this%nstates, &
                            max_krylov_vects=max_krylov_vects, &
                            fm_pool=this%fm_pool_mo)

         ist = 1
         iend = this%nstates
         DO iter = 1, niters
            nv = iend-ist+1

            ! MO to AO
            DO ivec = ist, iend
               iv = ivec-ist+1
               CALL cp_fm_set_all(this%evects_ao_ao(iv), 0.0_dp)
               ! Orthogonal_transform
               ! X_ao = C_alpha * X_mo * C_beta^T
               CALL parallel_gemm('N', 'N', this%nbf, this%nbf, gs_mos(1)%homo, &
                                  1.0_dp, gs_mos(1)%mos_occ, &
                                          this%evects_mo(1, ivec), &
                                  0.0_dp, this%work_ao_vir)
               CALL parallel_gemm('N', 'T', this%nbf, this%nbf, gs_mos(2)%homo, &
                                  1.0_dp, this%work_ao_vir, &
                                          gs_mos(2)%mos_virt, &
                                  0.0_dp, this%evects_ao_ao(iv))
               ! Transform to sym/asym
               ! symm
               CALL cp_fm_set_all(this%work_ao_ao_sym(1), 0.0_dp)
               CALL cp_fm_to_fm(this%evects_ao_ao(iv), &
                                this%work_ao_ao_sym(1)) ! X_ao
               CALL parallel_gemm('T', 'N', this%nbf, this%nbf, this%nbf, &
                                  0.5_dp, this%evects_ao_ao(iv), &
                                          this%work_ao_ao_sym(1), &
                                  0.5_dp, this%work_ao_ao_sym(1))  ! (X + X^T)/2

               CALL dbcsr_set(this%work_rho_ao_ao_symm(iv)%matrix, 0.0_dp)
               CALL copy_fm_to_dbcsr(this%work_ao_ao_sym(1), &
                                     this%work_rho_ao_ao_symm(iv)%matrix, &
                                     keep_sparsity=.TRUE.)
               ! asymm
               CALL cp_fm_set_all(this%work_ao_ao_asym(1), 0.0_dp)
               CALL cp_fm_to_fm(this%evects_ao_ao(iv), &
                                this%work_ao_ao_asym(1)) ! X_ao
               CALL parallel_gemm('T', 'N', this%nbf, this%nbf, this%nbf, &
                                 -0.5_dp, this%evects_ao_ao(iv), &
                                          this%work_ao_ao_asym(1), &
                                  0.5_dp, this%work_ao_ao_asym(1))  ! (X - X^T)/2

               CALL dbcsr_set(this%work_rho_ao_ao_asymm(iv)%matrix, 0.0_dp)
               CALL copy_fm_to_dbcsr(this%work_ao_ao_asym(1), &
                                     this%work_rho_ao_ao_asymm(iv)%matrix, &
                                     keep_sparsity=.TRUE.)

               CALL dbcsr_set(this%work_hmat_symm(iv)%matrix, 0.0_dp)
               CALL dbcsr_set(this%work_hmat_asymm(iv)%matrix, 0.0_dp)
            END DO

            CALL tddft_hfx_matrix(this%work_hmat_symm, this%work_rho_ao_ao_symm, qs_env)
            CALL tddft_hfx_matrix(this%work_hmat_asymm, this%work_rho_ao_ao_asymm, qs_env)

            DO ivec = ist, iend
               iv = ivec-ist+1

               CALL copy_dbcsr_to_fm(this%work_hmat_symm(iv)%matrix, this%work_ao_ao_sym(1))
               CALL copy_dbcsr_to_fm(this%work_hmat_asymm(iv)%matrix, this%work_ao_ao_asym(1))

               CALL cp_fm_set_all(this%hmat_ao_ao, 0.0_dp)
               CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=this%work_ao_ao_sym(1), &
                                        beta=0.0_dp, matrix_b=this%hmat_ao_ao)
               CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=this%work_ao_ao_asym(1), &
                                        beta=1.0_dp, matrix_b=this%hmat_ao_ao)

               CALL parallel_gemm('T', 'N', this%nocca, this%nbf, this%nbf, &
                                  1.0_dp, gs_mos(1)%mos_occ, &
                                          this%hmat_ao_ao, &
                                  0.0_dp, this%work_noc_ao)

               CALL parallel_gemm('N', 'N', this%nocca, nvirtb, this%nbf, &
                                  1.0_dp, this%work_noc_ao, &
                                          gs_mos(2)%mos_virt, &
                                  0.0_dp, davidson%Aop_krylov(1, ivec))

               IF (.TRUE.) THEN !ROHF or ROKS
                  CALL copy_dbcsr_to_fm(matrix_ks(1)%matrix, this%fock_ao_ao_alpha)
                  CALL copy_dbcsr_to_fm(matrix_ks(2)%matrix, this%fock_ao_ao_beta)

                  CALL cp_fm_set_all(this%Z_ext_noc_ao, 0.0_dp)
                  CALL cp_fm_to_fm_submat(this%evects_mo(1, ivec), &
                                          this%Z_ext_noc_ao, &
                                          nrow=this%nocca, &
                                          ncol=nvirtb, &
                                          s_firstrow=1, s_firstcol=1, &
                                          t_firstrow=1, t_firstcol=this%noccb + 1)

                  ! --- wrk1 = -fa_occ * abxc2  --------------------------------------
                  CALL parallel_gemm('N','N', this%nocca, nvirtb, this%nocca, &
                                     -1.0_dp, this%fock_noc_noc, &
                                              this%evects_mo(1, ivec), &
                                     0.0_dp,  this%work_noc_vir)

                  ! --- wrk1 +=  abxc2 * fb_vir  -------------------------------------
                  CALL parallel_gemm('N','N', this%nocca, nvirtb, nvirtb, &
                                      1.0_dp, this%evects_mo(1, ivec), &
                                              this%fock_vir_vir, &
                                      1.0_dp, this%work_noc_vir)

                  CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=this%work_noc_vir, &
                                           beta =1.0_dp, matrix_b=davidson%Aop_krylov(1, ivec))

               ELSE ! UHF or UKS
                  CALL cp_fm_get_info(this%evects_mo(1, ivec), local_data=evects_data)
                  CALL cp_fm_get_info(davidson%Aop_krylov(1, ivec), local_data=aop_data)
                  ij = 0
                  DO a = 1, nvirtb
                     DO i = 1, this%nocca
                        ij = ij + 1
                        aop_data(i, a) = aop_data(i, a) + this%e_virt_minus_occ(ij) * evects_data(i, a)
                     END DO
                  END DO
               END IF
            END DO

!           DO ivec = ist, iend
!             iv = ivec-ist+1
!
!             vbp => vb(:,noccb+1:)
!             vap => va(:,1:nocca)
!             call dgemm('t','n',nocca,nbf,nbf, &
!                        1.0_dp,vap,nbf,pao,nbf, &
!                        0.0_dp,scr,nocca)
!              │
!              call dgemm('n','n',nocca,nbf-noccb,nbf, &
!                         1.0_dp,scr,nocca,vbp,nbf, &
!                         0.0_dp,pmo,nocca)
!             IF (roref) THEN
!               call iatogen(evects_mo(:,ivec),abxc,nocca,noccb)
!
!               ! faz
!               call dgemm('n','n',nocca,nbf,nocca, &
!                          1.0_dp,fa,nbf,abxc,nbf, &
!                          0.0_dp,scr1,nbf)
!               ! zfb - faz
!               call dgemm('n','n',nocca,nbf,nbf, &
!                          1.0_dp,abxc,nbf,fb,nbf, &
!                         -1.0_dp,scr1,nbf)
!
!               call sfroesum(scr1,ab2_mo,nocca,noccb,ivec)
!             ELSE
!               call sfesum(mo_energy_a,mo_energy_b,ab2_mo,evects_mo,nocca,noccb,ivec)
!             END IF
!           END DO

            ! ERI contraction is here
            IF (tddfpt_control%do_hfx) THEN
               IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
                  ! full TDDFPT kernel
               ELSE
                  CPABORT("Kernel type not implemented")
               END IF
            END IF


!
            ! Solve lin eqs

            ! Print

!           mxerr = maxval(rnorm)

!           Check convergence
!           converged = mxerr<=cnvtol
!           if (converged) exit

!           No space left for new vectors, exit
!           if (nvec==mxvec) ierr = 1
!           if (ierr/=0) exit

!           call rpanewb(nstates,evects_mo,scr3,novec,nvec,ierr,tamm_dancoff=.true.)

  !         ierr=1 nvec over mxvec: not converged case
!           if (ierr/=0) exit

!           ist = novec+1
!           iend = nvec
!        end do
!        DO
!           ! *** perform Davidson iterations ***
!           conv = tddfpt_davidson_solver( &
!                  evects=evects, &
!                  evals=evals, &
!                  S_evects=S_evects, &
!                  gs_mos=gs_mos, &
!                  tddfpt_control=tddfpt_control, &
!                  matrix_ks=matrix_ks, &
!                  qs_env=qs_env, &
!                  kernel_env=kernel_env, &
!                  sub_env=sub_env, &
!                  logger=logger, &
!                  iter_unit=log_unit, &
!                  energy_unit=energy_unit, &
!                  tddfpt_print_section=tddfpt_print_section, &
!                  work_matrices=work_matrices)
!
!           ! at this point at least one of the following conditions are met:
!           ! a) convergence criteria has been achieved;
!           ! b) maximum number of iterations has been reached;
!           ! c) Davidson iterations must be restarted due to lack of Krylov vectors
!
!           CALL cp_iterate(logger%iter_info, increment=0, iter_nr_out=iter)
!           ! terminate the loop if either (a) or (b) is true ...
!           IF ((conv <= tddfpt_control%conv) .OR. iter >= niters) EXIT
!
!           ! ... otherwise restart Davidson iterations
!           evals = 0.0_dp
!           IF (log_unit > 0) THEN
!              WRITE (log_unit, '(1X,25("-"),1X,A,1X,25("-"))') "Restart Davidson iterations"
!              CALL m_flush(log_unit)
!           END IF
         END DO

         CALL cp_rm_iter_level(logger%iter_info, "TDDFT_SCF")

         ! print convergence summary
         IF (log_unit > 0) THEN
            CALL integer_to_string(iter, nstates_str)
            IF (conv <= tddfpt_control%conv) THEN
               WRITE (log_unit, "(1X,A)") "", &
                  "-------------------------------------------------------------------------------", &
                  "-   TDDFT run converged in "//TRIM(nstates_str)//" iteration(s) ", &
                  "-------------------------------------------------------------------------------"
            ELSE
               WRITE (log_unit, "(1X,A)") "", &
                  "-------------------------------------------------------------------------------", &
                  "-   TDDFT run did NOT converge after "//TRIM(nstates_str)//" iteration(s) ", &
                  "-------------------------------------------------------------------------------"
            END IF
         END IF

         CALL cp_print_key_finished_output(energy_unit, logger, &
                                           tddfpt_print_section, "DETAILED_ENERGY")
         CALL cp_print_key_finished_output(log_unit, logger, &
                                           tddfpt_print_section, "ITERATION_INFO")

         CALL davidson%clean(this%fm_pool_mo)
      ELSE
         CALL cp_warn(__LOCATION__, &
                      "Skipping TDDFPT wavefunction optimization")
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddft_run

   SUBROUTINE tddft_guess_vectors(evects_mo, evals, e_virt_minus_occ, &
                                  gs_mos, log_unit, tddfpt_control)
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(inout)   :: evects_mo
      REAL(kind=dp), DIMENSION(:), INTENT(inout)         :: evals
      REAL(kind=dp), DIMENSION(:)                        :: e_virt_minus_occ
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), INTENT(in)                        :: gs_mos
      INTEGER, INTENT(in)                                :: log_unit
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddft_get_guess'

      INTEGER                                            :: handle, imo_occ, imo_virt, &
                                                            istate, nstates, &
                                                            nstates_selected, istate_orig
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_m_occ

      CALL timeset(routineN, handle)


      nstates = SIZE(evects_mo, 2)

      tddfpt_control = tddfpt_control

      IF (debug_this_module) THEN
         CPASSERT(nstates > 0)
      END IF

      nstates_selected = gs_mos(1)%homo*(SIZE(gs_mos(2)%evals_full)-gs_mos(2)%homo)

      ALLOCATE (inds(nstates_selected))
      ALLOCATE (e_virt_m_occ(nstates_selected))

      istate = 0
      DO imo_virt = gs_mos(2)%homo+1, SIZE(gs_mos(2)%evals_full)
         DO imo_occ = 1, gs_mos(1)%homo
            istate = (imo_virt-gs_mos(2)%homo-1)*gs_mos(1)%homo+imo_occ
            e_virt_minus_occ(istate) = gs_mos(2)%evals_full(imo_virt)-gs_mos(1)%evals_full(imo_occ)
         END DO
      END DO

      e_virt_m_occ(1:nstates_selected) = e_virt_minus_occ(1:nstates_selected)

      CALL sort(e_virt_m_occ, nstates_selected, inds)

      IF (log_unit > 0) THEN
         WRITE (log_unit, "(1X,A)") "", &
            "-------------------------------------------------------------------------------", &
            "-                           SF-TDDFT Initial Guess                             -", &
            "-------------------------------------------------------------------------------"
         WRITE (log_unit, '(T11,A)') "State         Occupied      ->      Virtual          Excitation"
         WRITE (log_unit, '(T11,A)') "number         orbital              orbital          energy (eV)"
         WRITE (log_unit, '(1X,79("-"))')
      END IF

      DO istate = 1, nstates
         evals(istate) = e_virt_m_occ(istate)
         istate_orig = inds(istate)
         imo_occ = MOD(istate_orig - 1, gs_mos(1)%homo) + 1
         imo_virt = (istate_orig - 1) / gs_mos(1)%homo + 1
         CALL cp_fm_set_element(evects_mo(1, istate), imo_occ, imo_virt, 1.0_dp)
         IF (log_unit > 0) &
            WRITE (log_unit, '(T7,I8,T24,I8,T37,A5,T45,I8,T54,A5,T60,2F14.5)') &
            istate, imo_occ, '(alp)', gs_mos(2)%homo + imo_virt, '(bet)', e_virt_m_occ(istate)*evolt, &
            evects_mo(1, istate)%local_data(imo_occ, imo_virt)
      end do

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T7,A,T50,I24)') 'Number of active states:', &
            SIZE(gs_mos(1)%evals_occ)*SIZE(gs_mos(2)%evals_virt)
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      DEALLOCATE (inds, e_virt_m_occ)

      CALL timestop(handle)

   END SUBROUTINE tddft_guess_vectors

END MODULE qs_tddft_methods
