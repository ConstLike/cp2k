!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddft_methods
   USE cp_dbcsr_api,                    ONLY: dbcsr_type,&
                                              dbcsr_create,&
                                              dbcsr_get_info,&
                                              dbcsr_release,&
                                              dbcsr_p_type,&
                                              dbcsr_init_p,&
                                              dbcsr_release_p,&
                                              dbcsr_complete_redistribute,&
                                              dbcsr_distribution_release,&
                                              dbcsr_distribution_type,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_deallocate_matrix
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
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
                                              fm_pool_release
   USE cp_fm_pool_types,                ONLY: fm_pool_create_fm
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

   PUBLIC :: tddft_info_type

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

      !> MO space trial vectors: shape [nstates], each cp_fm_type of size [nocca x nvirb]
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: evects_mo

      !> AO space trial vectors: shape [nstates], each cp_fm_type of size [nbf x nbf]
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: evects_ao

      !> Matrix pools for efficient memory management of MO and AO trial vectors
      TYPE(cp_fm_pool_p_type), POINTER                   :: fm_pool_mo => NULL()
      TYPE(cp_fm_pool_p_type), POINTER                   :: fm_pool_ao => NULL()

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

! **************************************************************************************************

CONTAINS
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
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_mo
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_ao
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
      ALLOCATE (this%fm_pool_mo)
      NULLIFY (this%fm_pool_mo%pool)
      CALL cp_fm_struct_create(fm_struct_mo, &
                               nrow_global=this%nocca, &
                               ncol_global=nvirtb, &
                               context=sub_env%blacs_env)
      CALL fm_pool_create(this%fm_pool_mo%pool, fm_struct_mo)
      CALL cp_fm_struct_release(fm_struct_mo)

      ALLOCATE (this%evects_mo(this%nvecs))
      DO ivec = 1, this%nvecs
         CALL fm_pool_create_fm(this%fm_pool_mo%pool, &
                                this%evects_mo(ivec))
      END DO

      ! Initialize AO matrix pool and allocate AO trial vectors

      ALLOCATE (this%fm_pool_ao)
      NULLIFY (this%fm_pool_ao%pool)
      CALL cp_fm_struct_create(fm_struct_ao, &
                               nrow_global=this%nbf, &
                               ncol_global=this%nbf, &
                               context=sub_env%blacs_env)
      CALL fm_pool_create(this%fm_pool_ao%pool, fm_struct_ao)
      CALL cp_fm_struct_release(fm_struct_ao)

      ALLOCATE(this%evects_ao(this%nvecs))
      DO ivec = 1, this%nvecs
         CALL fm_pool_create_fm(this%fm_pool_ao%pool, this%evects_ao(ivec))
      END DO

      CALL cp_fm_struct_create(fm_struct_ao, &
                               nrow_global=this%nbf, &
                               ncol_global=this%nbf, &
                               context=sub_env%blacs_env)
      ALLOCATE (this%gs_rho_orb_ao_fm)
      CALL cp_fm_create(this%gs_rho_orb_ao_fm, fm_struct_ao)
      CALL cp_fm_struct_release(fm_struct_ao)

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

      CALL timestop(handle)

   END SUBROUTINE tddft_initialize_data

   SUBROUTINE tddft_setup_data(this, qs_env, gs_mos, logger, xc_section, &
         tddfpt_print_section, tddfpt_control, kernel_env, full_kernel_env, sub_env)
      CLASS(tddft_info_type), INTENT(INOUT)              :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(kernel_env_type)                              :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(full_kernel_env_type), TARGET                 :: full_kernel_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_setup_data'

      INTEGER                                            :: handle, log_unit
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

      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, &
                                      "GUESS_VECTORS", extension=".tddfptLog")
      gs_mos = gs_mos

!     CALL tddft_guess_vectors(evects, evals, tddft%evects_mo, gs_mos, log_unit, matrix_s, tddfpt_control)
      CALL cp_print_key_finished_output(log_unit, logger, &
                                        tddfpt_print_section, "GUESS_VECTORS")

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
         DO ivec = 1, SIZE(this%evects_mo)
            CALL cp_fm_release(this%evects_mo(ivec))
         END DO
         DEALLOCATE(this%evects_mo)
      END IF

      IF (ALLOCATED(this%evects_ao)) THEN
         DO ivec = 1, SIZE(this%evects_ao)
            CALL cp_fm_release(this%evects_ao(ivec))
         END DO
         DEALLOCATE(this%evects_ao)
      END IF

      IF (ASSOCIATED(this%fm_pool_mo)) THEN
         CALL fm_pool_release(this%fm_pool_mo%pool)
         DEALLOCATE(this%fm_pool_mo)
      END IF

      IF (ASSOCIATED(this%fm_pool_ao)) THEN
         CALL fm_pool_release(this%fm_pool_ao%pool)
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
      INTEGER                                            :: handle, log_unit

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

      this%nstates = this%nstates


!     niters = tddfpt_control%niters
!     IF (niters > 0) THEN
!        log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, &
!                                        "ITERATION_INFO", extension=".tddfptLog")
!        energy_unit = cp_print_key_unit_nr(logger, &
!                                           tddfpt_print_section, &
!                                           "DETAILED_ENERGY", &
!                                           extension=".tddfptLog")
!
!        IF (log_unit > 0) THEN
!           WRITE (log_unit, "(1X,A)") "", &
!              "-------------------------------------------------------------------------------", &
!              "-                     SF-TDDFT WAVEFUNCTION OPTIMIZATION                       -", &
!              "-------------------------------------------------------------------------------"
!
!           WRITE (log_unit, '(/,T11,A,T27,A,T40,A,T62,A)') "Step", "Time", "Convergence", "Conv. states"
!           WRITE (log_unit, '(1X,79("-"))')
!        END IF
!
!        CALL cp_add_iter_level(logger%iter_info, "TDDFT_SCF")
!
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
!        END DO
!
!        CALL cp_rm_iter_level(logger%iter_info, "TDDFT_SCF")
!
!        ! print convergence summary
!        IF (log_unit > 0) THEN
!           CALL integer_to_string(iter, nstates_str)
!           IF (conv <= tddfpt_control%conv) THEN
!              WRITE (log_unit, "(1X,A)") "", &
!                 "-------------------------------------------------------------------------------", &
!                 "-  TDDFPT run converged in "//TRIM(nstates_str)//" iteration(s) ", &
!                 "-------------------------------------------------------------------------------"
!           ELSE
!              WRITE (log_unit, "(1X,A)") "", &
!                 "-------------------------------------------------------------------------------", &
!                 "-  TDDFPT run did NOT converge after "//TRIM(nstates_str)//" iteration(s) ", &
!                 "-------------------------------------------------------------------------------"
!           END IF
!        END IF
!
!        CALL cp_print_key_finished_output(energy_unit, logger, &
!                                          tddfpt_print_section, "DETAILED_ENERGY")
!        CALL cp_print_key_finished_output(log_unit, logger, &
!                                          tddfpt_print_section, "ITERATION_INFO")
!     ELSE
!        CALL cp_warn(__LOCATION__, &
!                     "Skipping TDDFPT wavefunction optimization")
!     END IF
!
!     DEALLOCATE (e_virt_minus_occ)
!
!     DO istate = 1, nstates
!        CALL cp_fm_release(evects_mo(1, istate))
!     END DO
!     DEALLOCATE (evects_mo)

      CALL timestop(handle)
   END SUBROUTINE tddft_run

   SUBROUTINE tddft_get_guess(evects, evals, evects_mo, gs_mos, log_unit, matrix_s, tddfpt_control)
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(inout)   :: evects
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(inout)   :: evects_mo
      REAL(kind=dp), DIMENSION(:), INTENT(inout)         :: evals
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), INTENT(in)          :: gs_mos
      INTEGER, INTENT(in)                                :: log_unit
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddft_get_guess'

      CHARACTER(len=5)                                   :: spin_label1, spin_label2
      INTEGER                                            :: handle, imo_occ, imo_virt, ind, ispin, &
                                                            istate, jspin, spin1, spin2, nstates, &
                                                            nstates_occ_virt_alpha, &
                                                            nstates_selected
      INTEGER                                            :: nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ_avail, nmo_occ_selected, &
                                                            nmo_virt_selected
      REAL(kind=dp)                                      :: e_occ
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_minus_occ
      TYPE(cp_fm_struct_p_type), DIMENSION(maxspins)     :: fm_struct_evects

      CALL timeset(routineN, handle)

      nstates = SIZE(evects, 2)
      nspins = SIZE(gs_mos, 1)

      IF (debug_this_module) THEN
         CPASSERT(nstates > 0)
         CPASSERT(nspins == 1 .OR. nspins == 2)
      END IF

      DO ispin = 1, nspins
         ! number of occupied orbitals for each spin component
         nmo_occ_avail(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         ! number of occupied and virtual orbitals which can potentially
         ! contribute to the excited states in question.
         nmo_occ_selected(ispin) = MIN(nmo_occ_avail(ispin), nstates)
         nmo_virt_selected(ispin) = MIN(SIZE(gs_mos(ispin)%evals_virt), nstates)

         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, matrix_struct=fm_struct_evects(ispin)%struct)
      END DO

      ! TO DO: the variable 'nstates_selected' should probably be declared as INTEGER(kind=int_8),
      !        however we need a special version of the subroutine sort() in order to do so
      nstates_selected = nmo_occ_selected(1)*nmo_virt_selected(2)

      ALLOCATE (inds(nstates_selected))
      ALLOCATE (e_virt_minus_occ(nstates_selected))

      istate = 0
      ! Guess for spin-flip TDDFT
      DO imo_occ = 1, nmo_occ_selected(1)
         ! Here imo_occ enumerate alpha Occupied orbitals in inverse order (from the last to the first element)
         e_occ = gs_mos(1)%evals_occ(nmo_occ_avail(1) - imo_occ + 1)

         DO imo_virt = 1, nmo_virt_selected(2)
            istate = istate + 1
            e_virt_minus_occ(istate) = gs_mos(2)%evals_virt(imo_virt) - e_occ
         END DO
      END DO

      IF (debug_this_module) THEN
         CPASSERT(istate == nstates_selected)
      END IF

      CALL sort(e_virt_minus_occ, nstates_selected, inds)

      spin1 = 1
      spin2 = 2
      spin_label1 = '(alp)'
      spin_label2 = '(bet)'

      ! Calculate maximum number of spin-flip excitations
      nstates_occ_virt_alpha = nmo_occ_selected(1)*nmo_virt_selected(2)
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
         ! Index of excited state - 1
         ind = inds(istate) - 1

         ! Recover index of occupied MO (imo_occ) and unoccupied MO (imo_virt)
         ! associated to the excited state index (ind+1)
         imo_occ = nmo_occ_avail(spin1) - ind/nmo_virt_selected(spin2)
         imo_virt = MOD(ind, nmo_virt_selected(spin2)) + 1
         ! Assign initial guess for excitation energy
         evals(istate) = e_virt_minus_occ(istate)
         CALL cp_fm_set_element(evects_mo(1, istate), imo_occ, imo_virt, 1.0_dp)

         IF (log_unit > 0) &
            WRITE (log_unit, '(T7,I8,T24,I8,T37,A5,T45,I8,T54,A5,T60,F14.5)') &
            istate, imo_occ, spin_label1, nmo_occ_avail(spin2) + imo_virt, spin_label2, e_virt_minus_occ(istate)*evolt

         DO jspin = 1, SIZE(evects,1)
            ! .NOT. ASSOCIATED(evects(jspin, istate)%matrix_struct))
            CALL cp_fm_create(evects(jspin, istate), fm_struct_evects(jspin)%struct)
            CALL cp_fm_set_all(evects(jspin, istate), 0.0_dp)

            IF (jspin == spin1) &
            ! Half transform excitation vector to ao space:
            ! evects_mi = c_ma*X_ai
               CALL cp_fm_to_fm(gs_mos(spin2)%mos_virt, evects(spin1, istate), &
                                ncol=1, source_start=imo_virt, target_start=imo_occ)
         END DO
      END DO
            CALL tddfpt_print_excitation_analysis( &
               log_unit, &
               evects_mo, &
               evals, &
               gs_mos, &
               matrix_s(1)%matrix, &
               tddfpt_control%spinflip, &
               min_amplitude=tddfpt_control%min_excitation_amplitude)

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T7,A,T50,I24)') 'Number of active states:', &
            SIZE(gs_mos(1)%evals_occ)*SIZE(gs_mos(2)%evals_virt)
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      DEALLOCATE (e_virt_minus_occ)
      DEALLOCATE (inds)

      CALL timestop(handle)

   END SUBROUTINE tddft_guess_vectors

END MODULE qs_tddft_methods
