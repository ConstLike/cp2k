!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddft_methods
   USE qs_tddft_types,                  ONLY: tddft_info_type,&
                                              tddft_allocate_data,&
                                              tddft_deallocate_data,&
                                              allocate_fm_matrix,&
                                              deallocate_fm_matrix
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

   PUBLIC :: tddft_energies

CONTAINS

   SUBROUTINE tddft_energies(qs_env, sub_env, kernel_env, full_kernel_env, &
                             gs_mos, matrix_ks, logger, tddfpt_control, &
                             xc_section, tddfpt_print_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(kernel_env_type)                              :: kernel_env
      TYPE(full_kernel_env_type), TARGET                 :: full_kernel_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_energies'

      TYPE(tddft_info_type)                              :: tddft_env

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL tddft_allocate_data(tddft_env, qs_env, gs_mos, sub_env)

      CALL tddft_setup_data(tddft_env, qs_env, gs_mos, &
                            matrix_ks, logger, xc_section, &
                            tddfpt_print_section, tddfpt_control, &
                            kernel_env, full_kernel_env, sub_env)

      CALL tddft_run(tddft_env, qs_env, tddfpt_control, logger, &
                     tddfpt_print_section, gs_mos, kernel_env, &
                     matrix_ks, sub_env)

      CALL tddft_deallocate_data(tddft_env, sub_env, kernel_env)

      CALL timestop(handle)

   END SUBROUTINE tddft_energies

   SUBROUTINE tddft_setup_data(tddft_env, qs_env, gs_mos, matrix_ks, logger, xc_section, &
         tddfpt_print_section, tddfpt_control, kernel_env, full_kernel_env, sub_env)
      CLASS(tddft_info_type), INTENT(INOUT)              :: tddft_env
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

      tddft_env%nStates = tddft_env%nStates
      CALL tddfpt_construct_ground_state_orb_density( &
               rho_orb_struct=tddft_env%rho_orb_struct_sub, &
               rho_xc_struct=tddft_env%gs_rho_xc_struct, &
               is_rks_triplets=tddfpt_control%rks_triplets, &
               qs_env=qs_env, &
               sub_env=sub_env, &
               wfm_rho_orb=tddft_env%rho_ia_ao_fm)

      ! Full kernel
      CALL create_kernel_env(kernel_env=full_kernel_env, &
                             rho_struct_sub=tddft_env%rho_orb_struct_sub, &
                             xc_section=xc_section, &
                             is_rks_triplets=tddfpt_control%rks_triplets, &
                             sub_env=sub_env)
      kernel_env%full_kernel => full_kernel_env
      NULLIFY (kernel_env%admm_kernel)

      ! Transfom matrix_ks from AO to MO basis
      print *,  "nSpins = ", tddft_env%nSpins
      print *,  "nbf = ", tddft_env%nbf
      print *,  "nmo (alpha) = ", tddft_env%nocca + tddft_env%nvira
      print *,  "nmo (beta) = ", tddft_env%noccb + tddft_env%nvirb
      DO ispin = 1, tddft_env%nSpins
         CALL cp_fm_get_info(gs_mos(ispin)%mos_full, nrow_global=nrow, ncol_global=ncol)
         print *,  "gs_mos(ispin)%mos_full: ", nrow, ncol
         CALL cp_fm_get_info(tddft_env%fock(ispin), nrow_global=nrow, ncol_global=ncol)
         print *,  "tddft_env%fock(ispin): ", nrow, ncol
         CALL cp_fm_get_info(tddft_env%work_ao_ao(ispin), nrow_global=nrow, ncol_global=ncol)
         print *,  "tddft_env%work_ao_ao: ", nrow, ncol
      END DO

      DO ispin = 1, tddft_env%nSpins
         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, tddft_env%fock(ispin))
         ! work_ao_ao = gs_mos(1)%mos_full^T * tddft_env%fock
         CALL parallel_gemm('T', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, gs_mos(ispin)%mos_full, &
                                    tddft_env%fock(ispin), &
                            0.0_dp, tddft_env%work_ao_ao(ispin))

         ! tddft_env%fock = work_ao_ao * gs_mos(2)%mos_full
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, tddft_env%work_ao_ao(ispin), &
                                    gs_mos(ispin)%mos_full, &
                            0.0_dp, tddft_env%fock(ispin))
      END DO

      ! --- focka_occ = F_alpha( occ , occ )  ------------------------------------
      CALL cp_fm_to_fm_submat(tddft_env%fock(1), tddft_env%focka_noc_noc, &
                              nrow=tddft_env%nocca, ncol=tddft_env%nocca, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=1, t_firstcol=1)

      ! --- fockb_vir = F_beta( vir , vir )  ------------------------------------
      CALL cp_fm_to_fm_submat(tddft_env%fock(2), tddft_env%fockb_vir_vir, &
                              nrow=tddft_env%nvirb, &
                              ncol=tddft_env%nvirb, &
                              s_firstrow=tddft_env%noccb+1, &
                              s_firstcol=tddft_env%noccb+1, &
                              t_firstrow=1, t_firstcol=1)

      ! Open logger
      log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, &
                                      "GUESS_VECTORS", extension=".tddfptLog")

      CALL tddft_guess_vectors(tddft_env%evects_mo, tddft_env%evals, tddft_env%e_virt_minus_occ, &
                               gs_mos, tddft_env%nvec, tddft_env%nstates, log_unit, tddfpt_control)
      CALL cp_print_key_finished_output(log_unit, logger, &
                                        tddfpt_print_section, "GUESS_VECTORS")
      ! Close logger

      CALL timestop(handle)

   END SUBROUTINE tddft_setup_data

   SUBROUTINE tddft_run(tddft_env, qs_env, tddfpt_control, logger, tddfpt_print_section, &
                        gs_mos, kernel_env, matrix_ks, sub_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
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

         mxvec = tddft_env%nocca * tddft_env%nvirb
         max_krylov_vects = MIN(MAX(tddfpt_control%nkvs, tddft_env%nStates), INT(mxvec))

         ist = 1
         iend = tddft_env%nvec
         ierr = 0

         DO iter = 1, niters

            ! Transform trial vectors from MO to AO basis
            DO ivec = ist, iend
               ! Orthogonal_transform
               ! X_ao = C_alpha * X_mo * C_beta^T
               ! C[nbf,nvirb] = 1.0*A[nbf,nocca]*B[nocca,nvirb] + 0.0*C[nbf,nvirb]
               CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nvirb, tddft_env%nocca, &
                                  1.0_dp, gs_mos(1)%mos_occ, &
                                          tddft_env%evects_mo(ivec), &
                                  0.0_dp, tddft_env%work_ao_vir)
               ! C[nbf,nbf] = 1.0*A[nbf,nvirb]*(B[nbf,nvirb])^T + 0.0*C[nbf,nbf]
               CALL parallel_gemm('N', 'T', tddft_env%nbf, tddft_env%nbf, tddft_env%nvirb, &
                                  1.0_dp, tddft_env%work_ao_vir, &
                                          gs_mos(2)%mos_virt, &
                                  0.0_dp, tddft_env%evects_ao_ao(1))

               ! Print
               CALL print_mat_ao_ao(tddft_env%evects_ao_ao(1), 'evects', logger, ivec)

               ! Transform to sym/asym
               ! symm
               CALL cp_fm_set_all(matrix=tddft_env%work1_ao_ao, alpha=0.0_dp)
               ! work1 = 0.5*X
               CALL cp_fm_scale_and_add( &
                         alpha=0.0_dp, matrix_a=tddft_env%work1_ao_ao, &
                         beta=0.5_dp,  matrix_b=tddft_env%evects_ao_ao(1))
               ! work2 = X^T
               CALL cp_fm_transpose(tddft_env%evects_ao_ao(1), tddft_env%work2_ao_ao)

               ! work1 += 0.5*work2
               CALL cp_fm_scale_and_add( &
                         alpha=1.0_dp, matrix_a=tddft_env%work1_ao_ao, &
                         beta=0.5_dp,  matrix_b=tddft_env%work2_ao_ao)

               ! Print
               CALL print_mat_ao_ao(tddft_env%work1_ao_ao, 'work_sym', logger, ivec)

               ! Set density-like to 0.0
               CALL dbcsr_set(matrix=tddft_env%work_rho_ao_ao_symm(1)%matrix, &
                              alpha=0.0_dp)

               ! FM to DBCSR
               CALL copy_fm_to_dbcsr(fm=tddft_env%work1_ao_ao, &
                                     matrix=tddft_env%work_rho_ao_ao_symm(1)%matrix, &
                                     keep_sparsity=.TRUE.)
               ! asymm
               CALL cp_fm_set_all(matrix=tddft_env%work1_ao_ao, alpha=0.0_dp)

               ! work1 = 0.5*X
               CALL cp_fm_scale_and_add( &
                         alpha=0.0_dp, matrix_a=tddft_env%work1_ao_ao, &
                         beta=0.5_dp,  matrix_b=tddft_env%evects_ao_ao(1))
               ! work1 -= 0.5*work2
               CALL cp_fm_scale_and_add( &
                         alpha=1.0_dp, matrix_a=tddft_env%work1_ao_ao, &
                         beta=-0.5_dp, matrix_b=tddft_env%work2_ao_ao)

               ! Print
               CALL print_mat_ao_ao(tddft_env%work1_ao_ao, 'work_asym', logger, ivec)

               ! Set density-like to 0.0
               CALL dbcsr_set(matrix=tddft_env%work_rho_ao_ao_asymm(1)%matrix, &
                              alpha=0.0_dp)

               ! FM to DBCSR
               CALL copy_fm_to_dbcsr(fm=tddft_env%work1_ao_ao, &
                                     matrix=tddft_env%work_rho_ao_ao_asymm(1)%matrix, &
                                     keep_sparsity=.TRUE.)

               ! Set sigma to 0.0
               CALL dbcsr_set(matrix=tddft_env%work_hmat_symm(1)%matrix, alpha=0.0_dp)
               CALL dbcsr_set(matrix=tddft_env%work_hmat_asymm(1)%matrix, alpha=0.0_dp)

               !-----------------------------------------------
               ! ERI contraction: Sigma = ERI * X [DBCSR type]
               !-----------------------------------------------
               CALL tddft_hfx_matrix(matrix_ks=tddft_env%work_hmat_symm, &
                                     rho_ao=tddft_env%work_rho_ao_ao_symm, &
                                     qs_env=qs_env)
               CALL tddft_hfx_matrix(matrix_ks=tddft_env%work_hmat_asymm, &
                                     rho_ao=tddft_env%work_rho_ao_ao_asymm, &
                                     qs_env=qs_env)
 
               !-----------------------------------------------
               ! Transform back to MO basis
               !-----------------------------------------------
               ! BDCSR to FM
               CALL cp_fm_set_all(matrix=tddft_env%work1_ao_ao, alpha=0.0_dp)
               CALL cp_fm_set_all(matrix=tddft_env%work2_ao_ao, alpha=0.0_dp)

               CALL copy_dbcsr_to_fm(matrix=tddft_env%work_hmat_symm(1)%matrix, &
                                     fm=tddft_env%work1_ao_ao)
               CALL copy_dbcsr_to_fm(matrix=tddft_env%work_hmat_asymm(1)%matrix, &
                                     fm=tddft_env%work2_ao_ao)
               ! Print
               CALL print_mat_ao_ao(tddft_env%work1_ao_ao, 'Awork_sym', logger, ivec)
               CALL print_mat_ao_ao(tddft_env%work2_ao_ao, 'Awork_asym', logger, ivec)

               CALL cp_fm_set_all(matrix=tddft_env%hmat_ao_ao, alpha=0.0_dp)
               ! A = 0.0*A + 1.0*B
               CALL cp_fm_scale_and_add( &
                        alpha=0.0_dp, matrix_a=tddft_env%hmat_ao_ao, &
                        beta=1.0_dp,  matrix_b=tddft_env%work1_ao_ao)
               ! A = 1.0*A + 1.0*B
               CALL cp_fm_scale_and_add( &
                        alpha=1.0_dp, matrix_a=tddft_env%hmat_ao_ao, &
                        beta=1.0_dp,  matrix_b=tddft_env%work2_ao_ao)
               ! Print
               CALL print_mat_ao_ao(tddft_env%hmat_ao_ao, "hmat_ao_tot", logger, ivec)

               ! C[nocca,nbf] = 1.0*(A[nbf,nocca])^T * B[nbf,nbf] + 0.0*C[nocca,nbf]
               CALL parallel_gemm('T', 'N', tddft_env%nocca, tddft_env%nbf, tddft_env%nbf, &
                                  1.0_dp, gs_mos(1)%mos_occ, &
                                          tddft_env%hmat_ao_ao, &
                                  0.0_dp, tddft_env%work_noc_ao)

               ! C[nocca,nvirb] = 1.0*A[noc,nbf] * B[nbf,nvirb] + 0.0*C[nocca,nvirb]
               CALL parallel_gemm('N', 'N', tddft_env%nocca, tddft_env%nvirb, tddft_env%nbf, &
                                  1.0_dp, tddft_env%work_noc_ao, &
                                          gs_mos(2)%mos_virt, &
                                  0.0_dp, tddft_env%sigma_vecs(ivec))
               ! Print
               CALL print_mat_ao_ao(tddft_env%sigma_vecs(ivec), "sigma_vecs", logger, ivec)

               SELECT CASE (tddft_env%reference_type)
               CASE (1) ! ROKS/ROHF
                  ! Add ROHF corrections using Fock matrix contributions
                  ! Sigma = X*F_b - F_a*X
                  ! C[nocca,nvirb] = -1.0*A[nocca,nocca]*B[nocca,nvirb]+0.0*C[nocca,nvirb]
                  CALL parallel_gemm('N','N', tddft_env%nocca, tddft_env%nvirb, tddft_env%nocca, &
                                     -1.0_dp, tddft_env%focka_noc_noc, &
                                              tddft_env%evects_mo(ivec), &
                                      0.0_dp, tddft_env%work_noc_vir)

                  ! C[nocca,nvirb] = 1.0*A[nocca,nvirb]*B[nvirb,nvirb] + 1.0*C[nocca,nvirb]
                  CALL parallel_gemm('N','N', tddft_env%nocca, tddft_env%nvirb, tddft_env%nvirb, &
                                      1.0_dp, tddft_env%evects_mo(ivec), &
                                              tddft_env%fockb_vir_vir, &
                                      1.0_dp, tddft_env%work_noc_vir)

                  ! Add corrections to sigma_mo
                  ! A = 1.0 * A + 1.0 * B
                  CALL cp_fm_scale_and_add( &
                           alpha=1.0_dp, matrix_a=tddft_env%sigma_vecs(ivec), &
                           beta=1.0_dp,  matrix_b=tddft_env%work_noc_vir)

               CASE (2) ! UKS/UHF
                  ! Add orbital energy difference contributions for UHF
                  ij = 0
                  DO a = 1, tddft_env%nvirb
                     DO i = 1, tddft_env%nocca
                        ij = ij + 1

                        CALL cp_fm_get_element(matrix=tddft_env%sigma_vecs(ivec), &
                                               irow_global=i, &
                                               icol_global=a, &
                                               alpha=sig)
                        CALL cp_fm_get_element(matrix=tddft_env%evects_mo(ivec), &
                                               irow_global=i, &
                                               icol_global=a, &
                                               alpha=val)

                        sig = sig + tddft_env%e_virt_minus_occ(ij) * val

                        CALL cp_fm_set_element(matrix=tddft_env%sigma_vecs(ivec), &
                                               irow_global=i, &
                                               icol_global=a, &
                                               alpha=sig)
                     END DO
                  END DO
               END SELECT
               ! Print
               CALL print_mat_ao_ao(tddft_env%sigma_vecs(ivec), 'sigma_mat', logger, ivec)
            END DO

            ! alloc A[nvec,nvec]
            CALL allocate_fm_matrix(tddft_env%subspace, tddft_env%nvec, tddft_env%nvec, sub_env)
            CALL allocate_fm_matrix(tddft_env%subspace_evecs, tddft_env%nvec, tddft_env%nvec, sub_env)

            ! (1) build sub-space matrix Atilde = X^T * sigma * X
            ! Could be DGEMM
            ! C[nvec,nvec] = 1.0*(A(nvec)[dim])^T*B(nvec)[dim] + 0.0*C[nvec,nvec]
            DO i = 1, tddft_env%nvec
               DO j = 1, i
                  CALL cp_fm_trace(matrix_a=tddft_env%evects_mo(i), &
                                   matrix_b=tddft_env%sigma_vecs(j), &
                                   trace=val)
                  CALL cp_fm_set_element(matrix=tddft_env%subspace, &
                                         irow_global=i, &
                                         icol_global=j, &
                                         alpha=val)
                  IF (i/=j) CALL cp_fm_set_element(matrix=tddft_env%subspace, &
                                                   irow_global=j, &
                                                   icol_global=i, &
                                                   alpha=val)
               END DO
            END DO
            ! Print
            CALL print_mat_ao_ao(tddft_env%subspace, "subspace", logger, iter)

            ! (2) diagonalise Atilde  → omega_i & coeffs
            ! Diagonalize subspace matrix
            CALL choose_eigv_solver(matrix=tddft_env%subspace, &
                                    eigenvalues=tddft_env%ritz_evals(1:tddft_env%nvec), &
                                    eigenvectors=tddft_env%subspace_evecs, &
                                    info=ierr)
            CALL deallocate_fm_matrix(tddft_env%subspace)
            ! Print
            print *, "KKK| ================ Energy ============== at iter=", iter
            DO ivec = 1, tddft_env%nstates
               print *, "KKK| Energy (eV):", ivec, tddft_env%ritz_evals(ivec)*evolt
            END DO
            CALL print_mat_ao_ao(tddft_env%subspace_evecs, "subspace_evecs", logger, iter)

            ! Normalization
            CALL cp_fm_vectorsnorm(tddft_env%subspace_evecs, tddft_env%vector_norm(1:tddft_env%nvec))
            tddft_env%vector_scal(1:tddft_env%nvec) = 1.0_dp/tddft_env%vector_norm(1:tddft_env%nvec)
            DO i = 1, tddft_env%nvec
               PRINT *, "Norm of column ", i, ": ", tddft_env%vector_norm(i)
            END DO
            CALL cp_fm_column_scale(tddft_env%subspace_evecs, tddft_env%vector_scal(1:tddft_env%nvec))

            ! Print
            CALL print_mat_ao_ao(tddft_env%subspace_evecs, "normed subspace_evecs", logger, iter)

!           transformation_matrix = subspace_evecs
            ! Flag(a)
            DO iv = 1, tddft_env%nstates
               CALL cp_fm_set_all(matrix=tddft_env%residual_vecs(iv), &
                                  alpha=0.0_dp)

               DO j = 1, tddft_env%nvec
                  CALL cp_fm_get_element(matrix=tddft_env%subspace_evecs, &
                                         irow_global=j, &
                                         icol_global=iv, &
                                         alpha=val)
                  CALL cp_fm_scale_and_add( &
                           alpha=1.0_dp, matrix_a=tddft_env%residual_vecs(iv), &
                           beta=val,     matrix_b=tddft_env%sigma_vecs(j))
               END DO
            END DO

            ! Print
            DO iv = 1, tddft_env%nstates
               CALL print_mat_ao_ao(tddft_env%residual_vecs(iv), "residual_vecs", logger, iter)
            END DO

            ! Flag(b)
            DO i = 1, tddft_env%nstates
               DO j = 1, tddft_env%nvec
                  CALL cp_fm_get_element(matrix=tddft_env%subspace_evecs, &
                                         irow_global=j, &
                                         icol_global=i, &
                                         alpha=val)
                   val = -val*tddft_env%ritz_evals(i)
                  CALL cp_fm_set_element(matrix=tddft_env%subspace_evecs, &
                                         irow_global=j, &
                                         icol_global=i, &
                                         alpha=val)

               END DO
            END DO

            ! Flag(c)
            DO iv = 1, tddft_env%nstates
               DO j = 1, tddft_env%nvec
                  CALL cp_fm_get_element(matrix=tddft_env%subspace_evecs, &
                                         irow_global=j, &
                                         icol_global=iv, &
                                         alpha=val)
                  CALL cp_fm_scale_and_add( &
                           alpha=1.0_dp, matrix_a=tddft_env%residual_vecs(iv), &
                           beta=val,     matrix_b=tddft_env%evects_mo(j))
               END DO
            END DO
            CALL deallocate_fm_matrix(tddft_env%subspace_evecs)

            ! Print
            DO iv = 1, tddft_env%nstates
               CALL print_mat_ao_ao(tddft_env%residual_vecs(iv), "2 residual_vecs", logger, iv)
            END DO

            ! Flag(d)
            DO iv = 1, tddft_env%nstates
               CALL cp_fm_trace(tddft_env%residual_vecs(iv), &
                                tddft_env%residual_vecs(iv), val)
               tddft_env%residual_norms(iv) = val
            END DO

            DO iv = 1, tddft_env%nstates
               ij = 0
               DO a = 1, tddft_env%nvirb
                  DO i = 1, tddft_env%nocca
                     ij = ij + 1
                     val = tddft_env%ritz_evals(iv)-tddft_env%e_virt_minus_occ(ij)
                     val2 = abs(val)
                     IF (val2 < 1.e-12_dp) THEN
                        val = 1.e-5_dp
                     ELSE IF (val2< 1.e-5_dp) THEN
                        sig = val2/val
                        val = sig*1.e-5_dp
                     END IF
                     CALL cp_fm_get_element(tddft_env%residual_vecs(iv), i, a, scal)
                     scal = scal/val
                     CALL cp_fm_set_element(tddft_env%residual_vecs(iv), i, a, scal)
                  END DO
               END DO
            END DO
            IF (log_unit > 0) THEN
               WRITE (log_unit, '(/,4X,A,T14,A,T36,A)') "State", "Exc. energy (eV)", "Convergence (eV)"
               DO iv = 1, tddft_env%nstates
                  WRITE (log_unit, '(1X,I8,T12,F14.7,T38,ES11.4)') iv, &
                     tddft_env%ritz_evals(iv)*evolt, tddft_env%residual_norms(iv)
               END DO
               WRITE (log_unit, *)
               CALL m_flush(log_unit)
            END IF

            max_error = maxval(tddft_env%residual_norms)
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

            if (tddft_env%nvec==tddft_env%mxvec) ierr = 1
            if (ierr /= 0) EXIT

            tddft_env%novec = tddft_env%nvec
            ierr = 0

            DO k = 1, tddft_env%nstates
               DO iv = 1, tddft_env%nvec
                  CALL cp_fm_trace(tddft_env%evects_mo(iv), &
                                   tddft_env%residual_vecs(k), val)
                  CALL cp_fm_scale_and_add( &
                           alpha=1.0_dp, matrix_a=tddft_env%residual_vecs(k), &
                           beta=-val,     matrix_b=tddft_env%evects_mo(iv))
               END DO

               CALL cp_fm_trace(tddft_env%residual_vecs(k), &
                                tddft_env%residual_vecs(k), val)
               scal = sqrt(val)

               IF (scal<1.e-9_dp) CYCLE

               IF (tddft_env%nvec == tddft_env%mxvec) THEN
                  ierr = 2
                  RETURN
               END IF

               tddft_env%nvec = tddft_env%nvec + 1
               CALL cp_fm_scale_and_add( &
                        alpha=0._dp,     matrix_a=tddft_env%evects_mo(tddft_env%nvec), &
                        beta=1._dp/scal, matrix_b=tddft_env%residual_vecs(k))
            END DO

            ms = tddft_env%nvec - tddft_env%novec
            IF (ms == 0) ierr = 3

            IF (ierr /= 0) EXIT

            ist = tddft_env%novec+1
            iend = tddft_env%nvec

            ! Print
            print *, "ENERGY:",tddft_env%ritz_evals(1:tddft_env%nstates)

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
