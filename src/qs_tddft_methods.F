!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddft_methods
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_methods'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: tddft_info_type

! **************************************************************************************************
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
      TYPE(cp_fm_type) , ALLOCATABLE , DIMENSION(:) :: ritz_vecs      ! Y_k (MO basis)
      TYPE(cp_fm_type) , ALLOCATABLE , DIMENSION(:) :: ritz_sigma     ! Σ Y_k
      TYPE(cp_fm_type) , ALLOCATABLE , DIMENSION(:) :: residual_vecs  ! pre-conditioned R_k
      REAL(kind=dp)    , ALLOCATABLE , DIMENSION(:) :: residual_norms  ! pre-conditioned R_k
      REAL(dp)        , ALLOCATABLE , DIMENSION(:)  :: ritz_evals      ! ω_k
      INTEGER                                   :: novec = 0          ! sub-space size


      CONTAINS
         PROCEDURE :: init => tddft_initialize_data
         PROCEDURE :: setup => tddft_setup_data
         PROCEDURE :: run => tddft_run
         PROCEDURE :: clean => tddft_finalize_data
   END TYPE tddft_info_type

! **************************************************************************************************

CONTAINS
   SUBROUTINE alloc_fm_mat(mat, nrow, ncol, sub_env)
      TYPE(cp_fm_type), POINTER, INTENT(OUT)             :: mat
      INTEGER, INTENT(IN)                                :: nrow, ncol
      TYPE(tddfpt_subgroup_env_type), INTENT(IN)         :: sub_env

      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=nrow, &
                               ncol_global=ncol, &
                               context=sub_env%blacs_env)
      ALLOCATE(mat)
      CALL cp_fm_create(mat, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(mat, 0.0_dp)

   END SUBROUTINE alloc_fm_mat

   SUBROUTINE dealloc_fm_mat(mat)
      TYPE(cp_fm_type), POINTER, INTENT(INOUT)             :: mat

      IF (ASSOCIATED(mat)) THEN
         CALL cp_fm_release(mat)
         DEALLOCATE(mat)
      END IF
   END SUBROUTINE dealloc_fm_mat

   ! **************************************************************************************************
   !> \brief Initializes TDDFT matrix containers and internal pools for trial vector storage.
   !> \param this Object of type tddft_info_type
   !> \param qs_env Quickstep environment
   !> \param gs_mos Ground-state molecular orbitals (spin-resolved)
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
      TYPE(dbcsr_type), POINTER                          :: dbcsr_template_hfx
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_hfx
      INTEGER                                            :: ivec, ispin, handle
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s)

      CALL cp_fm_get_info(gs_mos(1)%mos_occ, ncol_global=this%nocca)
      CALL cp_fm_get_info(gs_mos(2)%mos_occ, ncol_global=this%noccb)
      CALL cp_fm_get_info(gs_mos(1)%mos_virt, ncol_global=this%nvira)
      CALL cp_fm_get_info(gs_mos(2)%mos_virt, ncol_global=this%nvirb)
      CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=this%nbf)

      ! Set the reference type ROHF/UHF
      this%reference_type = 0
      IF (dft_control%roks) this%reference_type = 1
      IF (dft_control%uks) this%reference_type = 2

      this%nStates = dft_control%tddfpt2_control%nStates
      this%nSpins = SIZE(gs_mos, 1)
      this%nStates_max = this%nocca * this%nvirb
      this%maxvec = dft_control%tddfpt2_control%nkvs
      this%xvec_dim = this%nStates_max
      this%mxvec = min(this%maxvec*this%nstates, this%xvec_dim)
      this%nStates = min(this%nstates, this%mxvec)
      this%nvec = min(max(2*this%nStates, 5), this%mxvec)
      this%novec = this%nvec

      IF (this%nStates > this%nStates_max) THEN
         CPWARN("NUMBER OF EXCITED STATES COULD LEAD TO PROBLEMS!")
         CPWARN("Experimental: CHANGED nStates TO ITS MAXIMUM VALUE!")
         this%nStates = this%nStates_max
      END IF

      !--------------------
      ! ALLOCATE  DP
      !--------------------
      ALLOCATE (this%e_virt_minus_occ(this%xvec_dim))
      ALLOCATE (this%evals(this%mxvec))
      ALLOCATE (this%ritz_evals(this%mxvec), source=0.0_dp)
      ALLOCATE (this%vector_scal(this%mxvec))
      ALLOCATE (this%vector_norm(this%mxvec))
      ALLOCATE (this%residual_norms(this%nstates), source=0.0_dp)

      !--------------------
      ! ALLOCATE  FM
      !--------------------
      !=== focka_noc_noc [nocca, nocca]
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=this%nocca, &
                               ncol_global=this%nocca, &
                               context=sub_env%blacs_env)
      ALLOCATE (this%focka_noc_noc)
      CALL cp_fm_create(this%focka_noc_noc, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      !=== fockb_vir_vir [nvirb, nvirb]
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=this%nvirb, &
                               ncol_global=this%nvirb, &
                               context=sub_env%blacs_env)
      ALLOCATE (this%fockb_vir_vir)
      CALL cp_fm_create(this%fockb_vir_vir, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      !=== work_noc_ao [nocca, nbf]
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=this%nocca, &
                               ncol_global=this%nbf, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_noc_ao)
      CALL cp_fm_create(this%work_noc_ao, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      !=== work_ao_vir [nbf, nvirb]
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=this%nbf, &
                               ncol_global=this%nvirb, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_ao_vir)
      CALL cp_fm_create(this%work_ao_vir, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      !=== rho_ia_ao_fm [nbf, nbf]
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=this%nbf, &
                               ncol_global=this%nbf, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%rho_ia_ao_fm)
      CALL cp_fm_create(this%rho_ia_ao_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      !--------------------
      ! ALLOCATE  FM POOL
      !--------------------
      !---------- MO POOL [nocca, nvirb] : fm_pool_mo
      ALLOCATE (this%fm_pool_mo)
      NULLIFY (this%fm_pool_mo%pool)
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=this%nocca, &
                               ncol_global=this%nvirb, &
                               context=sub_env%blacs_env)
      CALL fm_pool_create(this%fm_pool_mo%pool, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      !---------- END MO POOL [nocca, nvirb]

      ! MO POOL: evects_mo (mxvec)[nocca, nvirb]
      ALLOCATE (this%evects_mo(this%mxvec))
      DO ivec = 1, this%mxvec
         CALL fm_pool_create_fm(this%fm_pool_mo%pool, &
                                this%evects_mo(ivec))
      END DO
      ! MO POOL: sigma_vecs (mxvec)[nocca, nvirb]
      ALLOCATE (this%sigma_vecs(this%mxvec))
      DO ivec = 1, this%mxvec
         CALL fm_pool_create_fm(this%fm_pool_mo%pool, &
                                this%sigma_vecs(ivec))
      END DO
      ! MO POOL: work_noc_vir [nocca, nvirb]
      ALLOCATE(this%work_noc_vir)
      CALL fm_pool_create_fm(this%fm_pool_mo%pool, &
                             this%work_noc_vir)

      ! MO POOL: ritz [nocca, nvirb]
      ALLOCATE(this%ritz_vecs(this%nstates),   &
               this%ritz_sigma(this%nstates),  &
               this%residual_vecs(this%nstates))
      DO ivec = 1, this%nstates
         CALL fm_pool_create_fm(this%fm_pool_mo%pool, this%ritz_vecs(ivec))
         CALL fm_pool_create_fm(this%fm_pool_mo%pool, this%ritz_sigma(ivec))
         CALL fm_pool_create_fm(this%fm_pool_mo%pool, this%residual_vecs(ivec))
      END DO

      !---------- AO POOL [nbf, nbf] : fm_pool_ao
      ALLOCATE (this%fm_pool_ao(this%nSpins))
      DO ispin = 1, this%nSpins
         NULLIFY (this%fm_pool_ao(ispin)%pool)

         CALL cp_fm_struct_create(fm_struct, &
                                  nrow_global=this%nbf, &
                                  ncol_global=this%nbf, &
                                  context=sub_env%blacs_env)
         CALL fm_pool_create(this%fm_pool_ao(ispin)%pool, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END DO
      !---------- END AO POOL [nbf, nbf]

      ! AO POOL: Fock (nSpins)[nbf, nbf]
      ALLOCATE (this%fock(this%nSpins))
      DO ispin = 1, this%nSpins
         CALL fm_pool_create_fm(this%fm_pool_ao(ispin)%pool, &
                                this%fock(ispin))
      END DO
      ! AO POOL: evects_ao_ao (nspin)
      ALLOCATE (this%evects_ao_ao(1))
      DO ispin = 1, 1
         CALL fm_pool_create_fm(this%fm_pool_ao(ispin)%pool, &
                                this%evects_ao_ao(ispin))
      END DO
      ! AO POOL: Fock (nSpins)[nbf, nbf]
      ALLOCATE (this%work_ao_ao(this%nSpins))
      DO ispin = 1, this%nSpins
         CALL fm_pool_create_fm(this%fm_pool_ao(ispin)%pool, &
                                this%work_ao_ao(ispin))
      END DO

      ! AO POOL: Work
      ALLOCATE(this%work1_ao_ao)
      CALL fm_pool_create_fm(this%fm_pool_ao(1)%pool, this%work1_ao_ao)
      ALLOCATE(this%work2_ao_ao)
      CALL fm_pool_create_fm(this%fm_pool_ao(1)%pool, this%work2_ao_ao)
      ALLOCATE (this%hmat_ao_ao)
      CALL fm_pool_create_fm(this%fm_pool_ao(1)%pool, this%hmat_ao_ao)

      !--------------------
      ! ALLOCATE  DBCSR
      !--------------------
      CALL get_qs_env(qs_env, dbcsr_dist=sub_env%dbcsr_dist, sab_orb=sab_hfx)
      dbcsr_template_hfx => matrix_s(1)%matrix

      NULLIFY (this%work_rho_ao_ao_symm)
      CALL dbcsr_allocate_matrix_set(this%work_rho_ao_ao_symm, this%nSpins)
      NULLIFY (this%work_rho_ao_ao_asymm)
      CALL dbcsr_allocate_matrix_set(this%work_rho_ao_ao_asymm, this%nSpins)
      DO ispin = 1, this%nSpins
         CALL dbcsr_init_p(this%work_rho_ao_ao_symm(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist( &
                  this%work_rho_ao_ao_symm(ispin)%matrix, &
                  template=dbcsr_template_hfx, &
                  dbcsr_dist=sub_env%dbcsr_dist, &
                  sab=sab_hfx)

         CALL dbcsr_init_p(this%work_rho_ao_ao_asymm(ispin)%matrix)
         CALL dbcsr_create(this%work_rho_ao_ao_asymm(ispin)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric, &
                           template=this%work_rho_ao_ao_symm(ispin)%matrix)
         CALL dbcsr_complete_redistribute(this%work_rho_ao_ao_symm(ispin)%matrix, &
                                          this%work_rho_ao_ao_asymm(ispin)%matrix)
      END DO

      CALL dbcsr_allocate_matrix_set(this%work_hmat_symm, this%nSpins)
      CALL dbcsr_allocate_matrix_set(this%work_hmat_asymm, this%nSpins)
      DO ispin = 1, this%nSpins
         CALL dbcsr_init_p(this%work_hmat_symm(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(this%work_hmat_symm(ispin)%matrix, &
                                          template=dbcsr_template_hfx, &
                                          dbcsr_dist=sub_env%dbcsr_dist, &
                                          sab=sab_hfx)

         CALL dbcsr_init_p(this%work_hmat_asymm(ispin)%matrix)
         CALL dbcsr_create(this%work_hmat_asymm(ispin)%matrix, matrix_type=dbcsr_type_antisymmetric, &
                           template=this%work_hmat_symm(ispin)%matrix)
         CALL dbcsr_complete_redistribute(this%work_hmat_symm(ispin)%matrix, &
                                          this%work_hmat_asymm(ispin)%matrix)
      END DO

      NULLIFY (this%A_ia_munu_sub)
      CALL dbcsr_allocate_matrix_set(this%A_ia_munu_sub, 1)
      DO ispin = 1, 1
         CALL dbcsr_init_p(this%A_ia_munu_sub(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(this%A_ia_munu_sub(ispin)%matrix, &
                                          template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, &
                                          sab=sub_env%sab_orb)
      END DO

      !--------------------
      ! ALLOCATE  RHO
      !--------------------
      NULLIFY (this%rho_ia_ao)
      CALL dbcsr_allocate_matrix_set(this%rho_ia_ao, this%nSpins)
      DO ispin = 1, this%nSpins
         CALL dbcsr_init_p(this%rho_ia_ao(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(matrix=this%rho_ia_ao(ispin)%matrix, &
                                          template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, &
                                          sab=sub_env%sab_orb)
      END DO

      ! Allocate and set up qs_rho_type
      NULLIFY (this%rho_orb_struct_sub)
      ALLOCATE (this%rho_orb_struct_sub)
      CALL qs_rho_create(this%rho_orb_struct_sub)
      CALL qs_rho_set(this%rho_orb_struct_sub, rho_ao=this%rho_ia_ao)
      CALL qs_rho_rebuild(this%rho_orb_struct_sub, qs_env, &
                          rebuild_ao=.FALSE., &
                          rebuild_grids=.TRUE., &
                          pw_env_external=sub_env%pw_env)

      !--------------------
      ! ALLOCATE  PW
      !--------------------
      CALL pw_env_get(sub_env%pw_env, auxbas_pw_pool=auxbas_pw_pool)
      ALLOCATE (this%A_ia_rspace_sub(this%nSpins))
      DO ispin = 1, this%nspins
         CALL auxbas_pw_pool%create_pw(this%A_ia_rspace_sub(ispin))
      END DO
      ALLOCATE (this%rho_r(this%nspins))
      DO ispin = 1, this%nspins
         CALL auxbas_pw_pool%create_pw(this%rho_r(ispin))
      END DO
      ALLOCATE (this%rho_g(this%nspins))
      DO ispin = 1, this%nspins
         CALL auxbas_pw_pool%create_pw(this%rho_g(ispin))
      END DO
      ALLOCATE (this%tau_r(this%nspins))
      DO ispin = 1, this%nspins
         CALL auxbas_pw_pool%create_pw(this%tau_r(ispin))
      END DO
      ALLOCATE (this%work_v_xc(this%nspins))
      DO ispin = 1, this%nspins
         CALL auxbas_pw_pool%create_pw(this%work_v_xc(ispin))
      END DO
      ALLOCATE (this%work_v_xc_tau(this%nspins))
      DO ispin = 1, this%nspins
         CALL auxbas_pw_pool%create_pw(this%work_v_xc_tau(ispin))
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
      INTEGER                                            :: nrow, ncol

      CALL timeset(routineN, handle)

      this%nStates = this%nStates
      CALL tddfpt_construct_ground_state_orb_density( &
               rho_orb_struct=this%rho_orb_struct_sub, &
               rho_xc_struct=this%gs_rho_xc_struct, &
               is_rks_triplets=tddfpt_control%rks_triplets, &
               qs_env=qs_env, &
               sub_env=sub_env, &
               wfm_rho_orb=this%rho_ia_ao_fm)

      ! Full kernel
      CALL create_kernel_env(kernel_env=full_kernel_env, &
                             rho_struct_sub=this%rho_orb_struct_sub, &
                             xc_section=xc_section, &
                             is_rks_triplets=tddfpt_control%rks_triplets, &
                             sub_env=sub_env)
      kernel_env%full_kernel => full_kernel_env
      NULLIFY (kernel_env%admm_kernel)

      ! Transfom matrix_ks from AO to MO basis
      print *,  "nSpins = ", this%nSpins
      print *,  "nbf = ", this%nbf
      print *,  "nmo (alpha) = ", this%nocca + this%nvira
      print *,  "nmo (beta) = ", this%noccb + this%nvirb
      DO ispin = 1, this%nSpins
         CALL cp_fm_get_info(gs_mos(ispin)%mos_full, nrow_global=nrow, ncol_global=ncol)
         print *,  "gs_mos(ispin)%mos_full: ", nrow, ncol
         CALL cp_fm_get_info(this%fock(ispin), nrow_global=nrow, ncol_global=ncol)
         print *,  "this%fock(ispin): ", nrow, ncol
         CALL cp_fm_get_info(this%work_ao_ao(ispin), nrow_global=nrow, ncol_global=ncol)
         print *,  "this%work_ao_ao: ", nrow, ncol
      END DO

      DO ispin = 1, this%nSpins
         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, this%fock(ispin))
         ! work_ao_ao = gs_mos(1)%mos_full^T * this%fock
         CALL parallel_gemm('T', 'N', this%nbf, this%nbf, this%nbf, &
                            1.0_dp, gs_mos(ispin)%mos_full, &
                                    this%fock(ispin), &
                            0.0_dp, this%work_ao_ao(ispin))

         ! this%fock = work_ao_ao * gs_mos(2)%mos_full
         CALL parallel_gemm('N', 'N', this%nbf, this%nbf, this%nbf, &
                            1.0_dp, this%work_ao_ao(ispin), &
                                    gs_mos(ispin)%mos_full, &
                            0.0_dp, this%fock(ispin))
      END DO

      ! --- focka_occ = F_alpha( occ , occ )  ------------------------------------
      CALL cp_fm_to_fm_submat(this%fock(1), this%focka_noc_noc, &
                              nrow=this%nocca, ncol=this%nocca, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=1, t_firstcol=1)

      ! --- fockb_vir = F_beta( vir , vir )  ------------------------------------
      CALL cp_fm_to_fm_submat(this%fock(2), this%fockb_vir_vir, &
                              nrow=this%nvirb, &
                              ncol=this%nvirb, &
                              s_firstrow=this%noccb+1, &
                              s_firstcol=this%noccb+1, &
                              t_firstrow=1, t_firstcol=1)

      ! Open logger
      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, &
                                      "GUESS_VECTORS", extension=".tddfptLog")

      CALL tddft_guess_vectors(this%evects_mo, this%evals, this%e_virt_minus_occ, &
                               gs_mos, this%nvec, this%nstates, log_unit, tddfpt_control)
      CALL cp_print_key_finished_output(log_unit, logger, &
                                        tddfpt_print_section, "GUESS_VECTORS")
      ! Close logger

      CALL timestop(handle)

   END SUBROUTINE tddft_setup_data

   ! **************************************************************************************************
   !> \brief Releases all matrices and deallocates memory pools associated with TDDFT trial data.
   !> \param this Object of type tddft_info_type
   ! **************************************************************************************************
   SUBROUTINE tddft_finalize_data(this, sub_env, kernel_env)
      CLASS(tddft_info_type), INTENT(INOUT)              :: this
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(kernel_env_type)                              :: kernel_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_finalize_data'

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
      IF (ALLOCATED(this%ritz_vecs))  THEN
         DO ivec = SIZE(this%ritz_vecs),1,-1
            CALL fm_pool_give_back_fm(this%fm_pool_mo%pool,this%ritz_vecs(ivec))
            CALL fm_pool_give_back_fm(this%fm_pool_mo%pool,this%ritz_sigma(ivec))
            CALL fm_pool_give_back_fm(this%fm_pool_mo%pool,this%residual_vecs(ivec))
         END DO
         DEALLOCATE(this%ritz_vecs, this%ritz_sigma, this%residual_vecs)
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

      kernel_env = kernel_env
!     CALL release_kernel_env(kernel_env%full_kernel)

      CALL timestop(handle)

   END SUBROUTINE tddft_finalize_data


! **************************************************************************************************
!> \brief The energy calculation has been moved to its own subroutine
!> \param qs_env ...
!> \param nStates ...
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

      CHARACTER(len=20)                                  :: nStates_str
      INTEGER                                            :: handle, log_unit, &
         energy_unit, iter, niters, iend, ist, ivec, max_krylov_vects, &
         mxvec, i, a, ij, ierr, j, k, iv, ms
      REAL(kind=dp)                                      :: t1, t2
      REAL(KIND=DP)                                      :: conv, val, sig, scal, val2, max_error
      LOGICAL                                            :: is_converged

      CALL timeset(routineN, handle)

      matrix_ks = matrix_ks
      kernel_env = kernel_env

      t1 = m_walltime()

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

            WRITE (log_unit, '(1X,79("-"))')
         END IF

         CALL cp_add_iter_level(logger%iter_info, "TDDFT_SCF")

         mxvec = this%nocca * this%nvirb
         max_krylov_vects = MIN(MAX(tddfpt_control%nkvs, this%nStates), INT(mxvec))

         ist = 1
         iend = this%nvec
         ierr = 0

         DO iter = 1, niters

            ! Transform trial vectors from MO to AO basis
            DO ivec = ist, iend
               ! Orthogonal_transform
               ! X_ao = C_alpha * X_mo * C_beta^T
               ! C[nbf,nvirb] = 1.0*A[nbf,nocca]*B[nocca,nvirb] + 0.0*C[nbf,nvirb]
               CALL parallel_gemm('N', 'N', this%nbf, this%nvirb, this%nocca, &
                                  1.0_dp, gs_mos(1)%mos_occ, &
                                          this%evects_mo(ivec), &
                                  0.0_dp, this%work_ao_vir)
               ! C[nbf,nbf] = 1.0*A[nbf,nvirb]*(B[nbf,nvirb])^T + 0.0*C[nbf,nbf]
               CALL parallel_gemm('N', 'T', this%nbf, this%nbf, this%nvirb, &
                                  1.0_dp, this%work_ao_vir, &
                                          gs_mos(2)%mos_virt, &
                                  0.0_dp, this%evects_ao_ao(1))

               ! Print
!              IF (log_unit > 0) THEN
!                 CALL cp_fm_write_formatted(this%evects_ao_ao(1), log_unit)
!              END IF
               CALL print_mat_ao_ao(this%evects_ao_ao(1), 'evects', logger, ivec)

               ! Transform to sym/asym
               ! symm
               CALL cp_fm_set_all(matrix=this%work1_ao_ao, alpha=0.0_dp)
               ! work1 = 0.5*X
               CALL cp_fm_scale_and_add( &
                         alpha=0.0_dp, matrix_a=this%work1_ao_ao, &
                         beta=0.5_dp,  matrix_b=this%evects_ao_ao(1))
               ! work2 = X^T
               CALL cp_fm_transpose(this%evects_ao_ao(1), this%work2_ao_ao)

               ! work1 += 0.5*work2
               CALL cp_fm_scale_and_add( &
                         alpha=1.0_dp, matrix_a=this%work1_ao_ao, &
                         beta=0.5_dp,  matrix_b=this%work2_ao_ao)

               ! Print
               CALL print_mat_ao_ao(this%work1_ao_ao, 'work_sym', logger, ivec)

               ! Set density-like to 0.0
               CALL dbcsr_set(matrix=this%work_rho_ao_ao_symm(1)%matrix, &
                              alpha=0.0_dp)

               ! FM to DBCSR
               CALL copy_fm_to_dbcsr(fm=this%work1_ao_ao, &
                                     matrix=this%work_rho_ao_ao_symm(1)%matrix, &
                                     keep_sparsity=.TRUE.)
               ! asymm
               CALL cp_fm_set_all(matrix=this%work1_ao_ao, alpha=0.0_dp)

               ! work1 = 0.5*X
               CALL cp_fm_scale_and_add( &
                         alpha=0.0_dp, matrix_a=this%work1_ao_ao, &
                         beta=0.5_dp,  matrix_b=this%evects_ao_ao(1))
               ! work1 -= 0.5*work2
               CALL cp_fm_scale_and_add( &
                         alpha=1.0_dp, matrix_a=this%work1_ao_ao, &
                         beta=-0.5_dp, matrix_b=this%work2_ao_ao)

               ! Print
               CALL print_mat_ao_ao(this%work1_ao_ao, 'work_asym', logger, ivec)

               ! Set density-like to 0.0
               CALL dbcsr_set(matrix=this%work_rho_ao_ao_asymm(1)%matrix, &
                              alpha=0.0_dp)

               ! FM to DBCSR
               CALL copy_fm_to_dbcsr(fm=this%work1_ao_ao, &
                                     matrix=this%work_rho_ao_ao_asymm(1)%matrix, &
                                     keep_sparsity=.TRUE.)

               ! Set sigma to 0.0
               CALL dbcsr_set(matrix=this%work_hmat_symm(1)%matrix, alpha=0.0_dp)
               CALL dbcsr_set(matrix=this%work_hmat_asymm(1)%matrix, alpha=0.0_dp)

               !-----------------------------------------------
               ! ERI contraction: Sigma = ERI * X [DBCSR type]
               !-----------------------------------------------
               qs_env = qs_env
               CALL tddft_hfx_matrix(matrix_ks=this%work_hmat_symm, &
                                     rho_ao=this%work_rho_ao_ao_symm, &
                                     qs_env=qs_env)
               CALL tddft_hfx_matrix(matrix_ks=this%work_hmat_asymm, &
                                     rho_ao=this%work_rho_ao_ao_asymm, &
                                     qs_env=qs_env)
 
               !-----------------------------------------------
               ! Transform back to MO basis
               !-----------------------------------------------
               ! BDCSR to FM
               CALL cp_fm_set_all(matrix=this%work1_ao_ao, alpha=0.0_dp)
               CALL cp_fm_set_all(matrix=this%work2_ao_ao, alpha=0.0_dp)

               CALL copy_dbcsr_to_fm(matrix=this%work_hmat_symm(1)%matrix, &
                                     fm=this%work1_ao_ao)
               CALL copy_dbcsr_to_fm(matrix=this%work_hmat_asymm(1)%matrix, &
                                     fm=this%work2_ao_ao)
               ! Print
               CALL print_mat_ao_ao(this%work1_ao_ao, 'Awork_sym', logger, ivec)
               CALL print_mat_ao_ao(this%work2_ao_ao, 'Awork_asym', logger, ivec)

               CALL cp_fm_set_all(matrix=this%hmat_ao_ao, alpha=0.0_dp)
               ! A = 0.0*A + 1.0*B
               CALL cp_fm_scale_and_add( &
                        alpha=0.0_dp, matrix_a=this%hmat_ao_ao, &
                        beta=1.0_dp,  matrix_b=this%work1_ao_ao)
               ! A = 1.0*A + 1.0*B
               CALL cp_fm_scale_and_add( &
                        alpha=1.0_dp, matrix_a=this%hmat_ao_ao, &
                        beta=1.0_dp,  matrix_b=this%work2_ao_ao)
               ! Print
               CALL print_mat_ao_ao(this%hmat_ao_ao, "hmat_ao_tot", logger, ivec)

               ! C[nocca,nbf] = 1.0*(A[nbf,nocca])^T * B[nbf,nbf] + 0.0*C[nocca,nbf]
               CALL parallel_gemm('T', 'N', this%nocca, this%nbf, this%nbf, &
                                  1.0_dp, gs_mos(1)%mos_occ, &
                                          this%hmat_ao_ao, &
                                  0.0_dp, this%work_noc_ao)

               ! C[nocca,nvirb] = 1.0*A[noc,nbf] * B[nbf,nvirb] + 0.0*C[nocca,nvirb]
               CALL parallel_gemm('N', 'N', this%nocca, this%nvirb, this%nbf, &
                                  1.0_dp, this%work_noc_ao, &
                                          gs_mos(2)%mos_virt, &
                                  0.0_dp, this%sigma_vecs(ivec))
               ! Print
               CALL print_mat_ao_ao(this%sigma_vecs(ivec), "sigma_vecs", logger, ivec)

               SELECT CASE (this%reference_type)
               CASE (1) ! ROKS/ROHF
                  ! Add ROHF corrections using Fock matrix contributions
                  ! Sigma = X*F_b - F_a*X
                  ! C[nocca,nvirb] = -1.0*A[nocca,nocca]*B[nocca,nvirb]+0.0*C[nocca,nvirb]
                  CALL parallel_gemm('N','N', this%nocca, this%nvirb, this%nocca, &
                                     -1.0_dp, this%focka_noc_noc, &
                                              this%evects_mo(ivec), &
                                      0.0_dp, this%work_noc_vir)

                  ! C[nocca,nvirb] = 1.0*A[nocca,nvirb]*B[nvirb,nvirb] + 1.0*C[nocca,nvirb]
                  CALL parallel_gemm('N','N', this%nocca, this%nvirb, this%nvirb, &
                                      1.0_dp, this%evects_mo(ivec), &
                                              this%fockb_vir_vir, &
                                      1.0_dp, this%work_noc_vir)

                  ! Add corrections to sigma_mo
                  ! A = 1.0 * A + 1.0 * B
                  CALL cp_fm_scale_and_add( &
                           alpha=1.0_dp, matrix_a=this%sigma_vecs(ivec), &
                           beta=1.0_dp,  matrix_b=this%work_noc_vir)

               CASE (2) ! UKS/UHF
                  ! Add orbital energy difference contributions for UHF
                  ij = 0
                  DO a = 1, this%nvirb
                     DO i = 1, this%nocca
                        ij = ij + 1

                        CALL cp_fm_get_element(matrix=this%sigma_vecs(ivec), &
                                               irow_global=i, &
                                               icol_global=a, &
                                               alpha=sig)
                        CALL cp_fm_get_element(matrix=this%evects_mo(ivec), &
                                               irow_global=i, &
                                               icol_global=a, &
                                               alpha=val)

                        sig = sig + this%e_virt_minus_occ(ij) * val

                        CALL cp_fm_set_element(matrix=this%sigma_vecs(ivec), &
                                               irow_global=i, &
                                               icol_global=a, &
                                               alpha=sig)
                     END DO
                  END DO
               END SELECT
               ! Print
               CALL print_mat_ao_ao(this%sigma_vecs(ivec), 'sigma_mat', logger, ivec)
            END DO

            ! alloc A[nvec,nvec]
            CALL alloc_fm_mat(this%subspace, this%nvec, this%nvec, sub_env)
            CALL alloc_fm_mat(this%subspace_evecs, this%nvec, this%nvec, sub_env)

            ! (1) build sub-space matrix Atilde = X^T * sigma * X
            ! Could be DGEMM
            ! C[nvec,nvec] = 1.0*(A(nvec)[dim])^T*B(nvec)[dim] + 0.0*C[nvec,nvec]
            DO i = 1, this%nvec
               DO j = 1, i
                  CALL cp_fm_trace(matrix_a=this%evects_mo(i), &
                                   matrix_b=this%sigma_vecs(j), &
                                   trace=val)
                  CALL cp_fm_set_element(matrix=this%subspace, &
                                         irow_global=i, &
                                         icol_global=j, &
                                         alpha=val)
                  IF (i/=j) CALL cp_fm_set_element(matrix=this%subspace, &
                                                   irow_global=j, &
                                                   icol_global=i, &
                                                   alpha=val)
               END DO
            END DO
            ! Print
            CALL print_mat_ao_ao(this%subspace, "subspace", logger, iter)

            ! (2) diagonalise Atilde  → omega_i & coeffs
            ! Diagonalize subspace matrix
            CALL choose_eigv_solver(matrix=this%subspace, &
                                    eigenvalues=this%ritz_evals(1:this%nvec), &
                                    eigenvectors=this%subspace_evecs, &
                                    info=ierr)
            CALL dealloc_fm_mat(this%subspace)
            ! Print
            print *, "KKK| ================ Energy ============== at iter=", iter
            DO ivec = 1, this%nstates
               print *, "KKK| Energy (eV):", ivec, this%ritz_evals(ivec)*evolt
            END DO
            CALL print_mat_ao_ao(this%subspace_evecs, "subspace_evecs", logger, iter)

            ! Normalization
            CALL cp_fm_vectorsnorm(this%subspace_evecs, this%vector_norm(1:this%nvec))
            this%vector_scal(1:this%nvec) = 1.0_dp/this%vector_norm(1:this%nvec)
            DO i = 1, this%nvec
               PRINT *, "Norm of column ", i, ": ", this%vector_norm(i)
            END DO
            CALL cp_fm_column_scale(this%subspace_evecs, this%vector_scal(1:this%nvec))

            ! Print
            CALL print_mat_ao_ao(this%subspace_evecs, "normed subspace_evecs", logger, iter)

!           transformation_matrix = subspace_evecs
            ! Flag(a)
            DO iv = 1, this%nstates
               CALL cp_fm_set_all(matrix=this%residual_vecs(iv), &
                                  alpha=0.0_dp)

               DO j = 1, this%nvec
                  CALL cp_fm_get_element(matrix=this%subspace_evecs, &
                                         irow_global=j, &
                                         icol_global=iv, &
                                         alpha=val)
                  CALL cp_fm_scale_and_add( &
                           alpha=1.0_dp, matrix_a=this%residual_vecs(iv), &
                           beta=val,     matrix_b=this%sigma_vecs(j))
               END DO
            END DO

            ! Print
            DO iv = 1, this%nstates
               CALL print_mat_ao_ao(this%residual_vecs(iv), "residual_vecs", logger, iter)
            END DO

            ! Flag(b)
            DO i = 1, this%nstates
               DO j = 1, this%nvec
                  CALL cp_fm_get_element(matrix=this%subspace_evecs, &
                                         irow_global=j, &
                                         icol_global=i, &
                                         alpha=val)
                   val = -val*this%ritz_evals(i)
                  CALL cp_fm_set_element(matrix=this%subspace_evecs, &
                                         irow_global=j, &
                                         icol_global=i, &
                                         alpha=val)

               END DO
            END DO

            ! Flag(c)
            DO iv = 1, this%nstates
               DO j = 1, this%nvec
                  CALL cp_fm_get_element(matrix=this%subspace_evecs, &
                                         irow_global=j, &
                                         icol_global=iv, &
                                         alpha=val)
                  CALL cp_fm_scale_and_add( &
                           alpha=1.0_dp, matrix_a=this%residual_vecs(iv), &
                           beta=val,     matrix_b=this%evects_mo(j))
               END DO
            END DO
            CALL dealloc_fm_mat(this%subspace_evecs)

            ! Print
            DO iv = 1, this%nstates
               CALL print_mat_ao_ao(this%residual_vecs(iv), "2 residual_vecs", logger, iv)
            END DO

            ! Flag(d)
            DO iv = 1, this%nstates
               CALL cp_fm_trace(this%residual_vecs(iv), &
                                this%residual_vecs(iv), val)
               this%residual_norms(iv) = val
               print *, "dot_prod residual_vecs=", val
            END DO

            DO iv = 1, this%nstates
               ij = 0
               DO a = 1, this%nvirb
                  DO i = 1, this%nocca
                     ij = ij + 1
                     val = this%ritz_evals(iv)-this%e_virt_minus_occ(ij)
                     val2 = abs(val)
                     IF (val2 < 1.e-12_dp) THEN
                        val = 1.e-5_dp
                     ELSE IF (val2< 1.e-5_dp) THEN
                        sig = val2/val
                        val = sig*1.e-5_dp
                     END IF
                     CALL cp_fm_get_element(this%residual_vecs(iv), i, a, scal)
                     scal = scal/val
                     CALL cp_fm_set_element(this%residual_vecs(iv), i, a, scal)
                  END DO
               END DO
            END DO
            IF (log_unit > 0) THEN
               WRITE (log_unit, '(/,4X,A,T14,A,T36,A)') "State", "Exc. energy (eV)", "Convergence (eV)"
               DO iv = 1, this%nstates
                  WRITE (log_unit, '(1X,I8,T12,F14.7,T38,ES11.4)') iv, &
                     this%ritz_evals(iv)*evolt, this%residual_norms(iv)
               END DO
               WRITE (log_unit, *)
               CALL m_flush(log_unit)
            END IF

            max_error = maxval(this%residual_norms)
            t2 = m_walltime()
            IF (log_unit > 0) THEN
               WRITE (log_unit, '(T7,I8,T24,F7.1,T40,ES11.4,T66)') iter, t2 - t1, max_error
               CALL m_flush(log_unit)
            END IF
            t1 = t2

            is_converged = .FALSE.
            is_converged = max_error<=1.e-6_dp
            print *, "is_converged=",is_converged, max_error
            IF (is_converged) EXIT

            if (this%nvec==this%mxvec) ierr = 1
            if (ierr /= 0) EXIT

            this%novec = this%nvec
            ierr = 0

            DO k = 1, this%nstates
               DO iv = 1, this%nvec
                  CALL cp_fm_trace(this%evects_mo(iv), &
                                   this%residual_vecs(k), val)
                  CALL cp_fm_scale_and_add( &
                           alpha=1.0_dp, matrix_a=this%residual_vecs(k), &
                           beta=-val,     matrix_b=this%evects_mo(iv))
               END DO

               CALL cp_fm_trace(this%residual_vecs(k), &
                                this%residual_vecs(k), val)
               scal = sqrt(val)

               IF (scal<1.e-9_dp) CYCLE

               IF (this%nvec == this%mxvec) THEN
                  ierr = 2
                  RETURN
               END IF

               this%nvec = this%nvec + 1
               CALL cp_fm_scale_and_add( &
                        alpha=0._dp,     matrix_a=this%evects_mo(this%nvec), &
                        beta=1._dp/scal, matrix_b=this%residual_vecs(k))
            END DO

            ms = this%nvec - this%novec
            IF (ms == 0) ierr = 3

            IF (ierr /= 0) EXIT

            ist = this%novec+1
            iend = this%nvec

            ! Print
            print *, "ENERGY:",this%ritz_evals(1:this%nstates)

!
            ! ERI contraction is here
            IF (tddfpt_control%do_hfx) THEN
               IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
                  ! full TDDFPT kernel
               ELSE
                  CPABORT("Kernel type not implemented")
               END IF
            END IF


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
            CALL integer_to_string(iter, nStates_str)
            IF (conv <= tddfpt_control%conv) THEN
               WRITE (log_unit, "(1X,A)") "", &
                  "-------------------------------------------------------------------------------", &
                  "-   TDDFT run converged in "//TRIM(nStates_str)//" iteration(s) ", &
                  "-------------------------------------------------------------------------------"
            ELSE
               WRITE (log_unit, "(1X,A)") "", &
                  "-------------------------------------------------------------------------------", &
                  "-   TDDFT run did NOT converge after "//TRIM(nStates_str)//" iteration(s) ", &
                  "-------------------------------------------------------------------------------"
            END IF
         END IF

         CALL cp_print_key_finished_output(energy_unit, logger, &
                                           tddfpt_print_section, "DETAILED_ENERGY")
         CALL cp_print_key_finished_output(log_unit, logger, &
                                           tddfpt_print_section, "ITERATION_INFO")
      ELSE
         CALL cp_warn(__LOCATION__, &
                      "Skipping TDDFPT wavefunction optimization")
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddft_run

   SUBROUTINE print_mat_ao_ao(mat, mat_name, logger, iter)
      TYPE(cp_fm_type), INTENT(IN) :: mat
      TYPE(cp_logger_type), POINTER      :: logger
      INTEGER                             :: iter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_hmat_ao_ao'
      CHARACTER(LEN=*) :: mat_name
      INTEGER :: handle, i, j, unit_nr
      REAL(KIND=dp), DIMENSION(5,5) :: submatrix

      CALL timeset(routineN, handle)

      unit_nr = cp_logger_get_default_io_unit(logger)

      CALL cp_fm_get_submatrix(mat, submatrix, &
                               start_row=1, start_col=1, &
                               n_rows=5, n_cols=5)

      IF (unit_nr > 0) THEN
         WRITE(unit_nr, '(A,A,A,I0,A)') "First 5x5 elements of ",TRIM(mat_name)," at iter ",iter,":"
         DO i = 1, 5
            WRITE(unit_nr, '(5F12.6)') (submatrix(i,j), j=1,5)
         END DO
         WRITE(unit_nr, '(A)') ""
         CALL m_flush(unit_nr)
      END IF

      CALL timestop(handle)
   END SUBROUTINE print_mat_ao_ao

   SUBROUTINE tddft_guess_vectors(evects_mo, evals, e_virt_minus_occ, &
                                  gs_mos, nvec, nStates, log_unit, tddfpt_control)
      TYPE(cp_fm_type), DIMENSION(:), INTENT(inout)   :: evects_mo
      REAL(kind=dp), DIMENSION(:), INTENT(inout)         :: evals
      REAL(kind=dp), DIMENSION(:)                        :: e_virt_minus_occ
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), INTENT(in)                        :: gs_mos
      INTEGER                                            :: nStates, nvec
      INTEGER, INTENT(in)                                :: log_unit
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddft_get_guess'

      INTEGER                                            :: handle, imo_occ, imo_virt, &
                                                            istate, &
                                                            nStates_selected, istate_orig
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_m_occ

      CALL timeset(routineN, handle)


      tddfpt_control = tddfpt_control

      IF (debug_this_module) THEN
         CPASSERT(nStates > 0)
      END IF

      nStates_selected = gs_mos(1)%homo*(SIZE(gs_mos(2)%evals_full)-gs_mos(2)%homo)

      ALLOCATE (inds(nStates_selected))
      ALLOCATE (e_virt_m_occ(nStates_selected))

      istate = 0
      DO imo_virt = gs_mos(2)%homo+1, SIZE(gs_mos(2)%evals_full)
         DO imo_occ = 1, gs_mos(1)%homo
            istate = (imo_virt-gs_mos(2)%homo-1)*gs_mos(1)%homo+imo_occ
            e_virt_minus_occ(istate) = gs_mos(2)%evals_full(imo_virt)-gs_mos(1)%evals_full(imo_occ)
         END DO
      END DO

      e_virt_m_occ(1:nStates_selected) = e_virt_minus_occ(1:nStates_selected)

      CALL sort(e_virt_m_occ, nStates_selected, inds)

      IF (log_unit > 0) THEN
         WRITE (log_unit, "(1X,A)") "", &
            "-------------------------------------------------------------------------------", &
            "-                           SF-TDDFT Initial Guess                             -", &
            "-------------------------------------------------------------------------------"
         WRITE (log_unit, '(T11,A)') "State         Occupied      ->      Virtual          Excitation"
         WRITE (log_unit, '(T11,A)') "number         orbital              orbital          energy (eV)"
         WRITE (log_unit, '(1X,79("-"))')
      END IF

      DO istate = 1, nVec
         evals(istate) = e_virt_m_occ(istate)
         istate_orig = inds(istate)
         imo_occ = MOD(istate_orig - 1, gs_mos(1)%homo) + 1
         imo_virt = (istate_orig - 1) / gs_mos(1)%homo + 1
         CALL cp_fm_set_element(evects_mo(istate), imo_occ, imo_virt, 1.0_dp)
         IF (istate > nStates) CYCLE
         IF (log_unit > 0) &
            WRITE (log_unit, '(T7,I8,T24,I8,T37,A5,T45,I8,T54,A5,T60,2F14.5)') &
            istate, imo_occ, '(alp)', gs_mos(2)%homo + imo_virt, '(bet)', e_virt_m_occ(istate)*evolt, &
            evects_mo(istate)%local_data(imo_occ, imo_virt)
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
