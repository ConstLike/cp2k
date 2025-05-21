!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddft_methods
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type,&
                                              mp_request_type
   USE kahan_sum,                       ONLY: accurate_dot_product
   USE input_constants,                 ONLY: tddfpt_dipole_berry,&
                                              tddfpt_dipole_length
   USE qs_tddft_theory_sf,              ONLY: tddft_sf_theory_type
   USE qs_tddft_theory,                 ONLY: tddft_theory_type
   USE qs_tddft_types,                  ONLY: tddft_info_type,&
                                              tddft_initialization,&
                                              tddft_finalization,&
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
                                              cp_fm_frobenius_norm,&
                                              cp_fm_transpose,&
                                              cp_fm_column_scale,&
                                              cp_fm_row_scale
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
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              int_8
   USE lri_environment_methods,         ONLY: lri_print_stat
   USE lri_environment_types,           ONLY: lri_density_release,&
                                              lri_env_release
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
                                              tddfpt_print_exciton_descriptors,&
                                              tddfpt_print_nto_analysis
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

      full_kernel_env = full_kernel_env
      xc_section = xc_section
      kernel_env = kernel_env

      CALL tddft_initialization(tddft_env, qs_env, gs_mos, sub_env)

      CALL tddft_setup_data(tddft_env, qs_env, gs_mos, &
                            matrix_ks, tddfpt_control, sub_env)

      CALL tddft_run(tddft_env, qs_env, tddfpt_control, logger, &
                     tddfpt_print_section, gs_mos, kernel_env, &
                     sub_env)

      CALL tddft_finalization(tddft_env, sub_env, kernel_env)

      CALL timestop(handle)

   END SUBROUTINE tddft_energies

   SUBROUTINE tddft_setup_data(tddft_env, qs_env, gs_mos, matrix_ks, &
                               tddfpt_control, sub_env)
      CLASS(tddft_info_type), INTENT(INOUT)              :: tddft_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_setup_data'

      INTEGER                                            :: handle, ispin

      CALL timeset(routineN, handle)

      CALL tddfpt_construct_ground_state_orb_density( &
               rho_orb_struct=tddft_env%rho_orb_struct_sub, &
               rho_xc_struct=tddft_env%gs_rho_xc_struct, &
               is_rks_triplets=tddfpt_control%rks_triplets, &
               qs_env=qs_env, &
               sub_env=sub_env, &
               wfm_rho_orb=tddft_env%rho_ia_ao_fm)

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


      CALL timestop(handle)

   END SUBROUTINE tddft_setup_data

   SUBROUTINE tddft_run(tddft_env, qs_env, tddfpt_control, logger, tddfpt_print_section, &
                        gs_mos, kernel_env, sub_env)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(kernel_env_type)                              :: kernel_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'tddft_run'

      CHARACTER(len=20)                                  :: nStates_str
      INTEGER                                            :: handle, log_unit, &
         energy_unit, iter, niters, iend, ist, ivec, max_krylov_vects, &
         mxvec, iv
      REAL(kind=dp)                                      :: t1, t2
      CLASS(tddft_theory_type), ALLOCATABLE              :: theory_handler
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: ostrength
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s

      CALL timeset(routineN, handle)

      ALLOCATE(tddft_sf_theory_type :: theory_handler)

      t1 = m_walltime()

      niters = tddfpt_control%niters
      IF (niters > 0) THEN
         ! Open logger
         log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, &
                                         "GUESS_VECTORS", extension=".tddfptLog")

         CALL theory_handler%initialize_guess_vectors(tddft_env, gs_mos, tddfpt_control, log_unit)

         ! Close logger
         CALL cp_print_key_finished_output(log_unit, logger, &
                                           tddfpt_print_section, "GUESS_VECTORS")

         log_unit = cp_print_key_unit_nr(logger, tddfpt_print_section, &
                                         "ITERATION_INFO", extension=".tddfptLog")
         energy_unit = cp_print_key_unit_nr(logger, &
                                            tddfpt_print_section, &
                                            "DETAILED_ENERGY", &
                                            extension=".tddfptLog")

         IF (log_unit > 0) THEN
            WRITE (log_unit, "(1X,A)") "", &
               "-------------------------------------------------------------------------------", &
               "-                     SF-TDDFT WAVEFUNCTION OPTIMIZATION                      -", &
               "-------------------------------------------------------------------------------"

            WRITE (log_unit, '(1X,79("-"))')
         END IF

         CALL cp_add_iter_level(logger%iter_info, "TDDFT_SCF")

         mxvec = tddft_env%nocca * tddft_env%nvirb
         max_krylov_vects = MIN(MAX(tddfpt_control%nkvs, tddft_env%nstates), INT(mxvec))

         tddft_env%max_error_prev = 10.0_dp
         tddft_env%max_delta_eval_prev = 10.0_dp

         ! Outer loop for managing restarts
         DO
            ist = 1
            iend = tddft_env%nvec
            tddft_env%status(1:2) = .FALSE.
            tddft_env%stall_count = 0

            ! Davidson iterations
            DO
               CALL cp_iterate(logger%iter_info, iter_nr_out=iter)

               DO ivec = ist, iend
                  !> MO->AO transformation for trial vectors
                  CALL theory_handler%mo_to_ao_transform(tddft_env, gs_mos, ivec, logger)

                  !> Apply all operators (HFX, XC, ...), build Sigma vectors (=A*X) in AO space
                  CALL theory_handler%apply_operator(qs_env, tddft_env, kernel_env, &
                                                     sub_env, tddfpt_control, ivec, logger)

                  !> AO->MO transformation for resulting vectors
                  CALL theory_handler%ao_to_mo_transform(tddft_env,gs_mos,ivec, logger)
               END DO

               !> Run converger and expand the vectors space if needed
               tddft_env%status = run_converger(tddft_env, sub_env, tddfpt_control, &
                                                log_unit, logger, iter)

               !> Print iteration time and max_error
               t2 = m_walltime()
               IF (log_unit > 0) THEN
                  WRITE (log_unit, '(T7,I8,T24,F7.1,T40,ES11.4,T66)') iter, t2 - t1, tddft_env%max_error
                  CALL m_flush(log_unit)
               END IF
               t1 = t2

               ! Check convergence status
               IF (tddft_env%status(1)) THEN
                  ! Update guess vectors to latest Ritz vectors and eigenvalues
                  DO iv = 1, tddft_env%nstates
                     CALL cp_fm_to_fm(source=tddft_env%ritz_vects(iv), &
                                      destination=tddft_env%evects_mo(iv))

                     tddft_env%evals(iv) = tddft_env%ritz_evals(iv)
                  END DO

                  EXIT    ! Exit Davidson loop
               END IF

               ! Exit to restart Davidson iterations
               IF (tddft_env%status(2)) EXIT

               ist = tddft_env%novec+1
               iend = tddft_env%nvec
            END DO ! Davidson iterations

            ! Check conditions for exiting or restarting
            ! at this point at least one of the following conditions are met:
            ! a) convergence criteria has been achieved;
            ! b) maximum number of iterations has been reached;
            ! c) Davidson iterations must be restarted due to lack of Krylov vectors

            CALL cp_iterate(logger%iter_info, increment=0, iter_nr_out=iter)

            ! Terminate the loop if either (a) or (b) is True
            IF (tddft_env%status(1) .OR. iter >= niters) EXIT

            IF (log_unit > 0) THEN
               WRITE (log_unit, '(1X,25("-"),1X,A,1X,25("-"))') "Restart Davidson iterations"
               CALL m_flush(log_unit)
            END IF
         END DO ! Outer loop for managing restarts

         CALL cp_rm_iter_level(logger%iter_info, "TDDFT_SCF")

         ! print convergence summary
         IF (log_unit > 0) THEN
            CALL integer_to_string(iter, nStates_str)
            IF (tddft_env%status(1)) THEN
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
      ! *** print summary information ***
      log_unit = cp_logger_get_default_io_unit(logger)
      ALLOCATE (ostrength(tddft_env%nstates), source=0.0_dp)
      CALL tddfpt_print_summary(tddft_env= tddft_env, &
                                log_unit=log_unit, &
                                evects=tddft_env%evects_mo, &
                                evals=tddft_env%evals, &
                                ostrength=ostrength, &
                                mult=3, &
                                dipole_form=tddfpt_control%dipole_form)
      matrix_s = matrix_s

      CALL timestop(handle)
   END SUBROUTINE tddft_run
   SUBROUTINE tddfpt_print_excitation_analysis(tddft_env, &
         log_unit, evects, evals, gs_mos, matrix_s, spinflip, min_amplitude)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      INTEGER, INTENT(in)                                :: log_unit
      TYPE(cp_fm_type), DIMENSION(:), INTENT(in)      :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(dbcsr_type), POINTER                          :: matrix_s
      INTEGER                                            :: spinflip
      REAL(kind=dp), INTENT(in)                          :: min_amplitude

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_print_excitation_analysis'

      CHARACTER(len=5)                                   :: spin_label, spin_label2
      INTEGER                                            :: handle, icol, iproc, irow, ispin, spin2, &
                                                            istate, nao, ncols_local, nrows_local, &
                                                            nspins, nstates, state_spin, state_spin2
      INTEGER(kind=int_8)                                :: iexc, imo_occ, imo_virt, ind, nexcs, &
                                                            nexcs_local, nexcs_max_local, &
                                                            nmo_virt_occ, nmo_virt_occ_alpha
      INTEGER(kind=int_8), ALLOCATABLE, DIMENSION(:)     :: inds_local, inds_recv, nexcs_recv
      INTEGER(kind=int_8), DIMENSION(1)                  :: nexcs_send
      INTEGER(kind=int_8), DIMENSION(maxspins)           :: nmo_occ8, nmo_virt8
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ, nmo_virt
      LOGICAL                                            :: do_exc_analysis
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: weights_local, weights_neg_abs_recv, &
                                                            weights_recv
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: S_mos_virt, weights_fm
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(mp_request_type)                              :: send_handler, send_handler2
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: recv_handlers, recv_handlers2

      CALL timeset(routineN, handle)

      nspins = tddft_env%nspins
      nstates = tddft_env%nstates
      do_exc_analysis = min_amplitude < 1.0_dp

      CALL cp_fm_get_info(gs_mos(1)%mos_occ, context=blacs_env, para_env=para_env)
      CALL dbcsr_get_info(matrix_s, nfullrows_total=nao)

      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
         nmo_occ8(ispin) = SIZE(gs_mos(ispin)%evals_occ, kind=int_8)
         nmo_virt8(ispin) = SIZE(gs_mos(ispin)%evals_virt, kind=int_8)
      END DO

      ! *** excitation analysis ***
      IF (do_exc_analysis) THEN
         CPASSERT(log_unit <= 0 .OR. para_env%is_source())
         nmo_virt_occ_alpha = INT(nmo_virt(1), int_8)*INT(nmo_occ(1), int_8)

         IF (log_unit > 0) THEN
            WRITE (log_unit, "(1X,A)") "", &
               "-------------------------------------------------------------------------------", &
               "-                            Excitation analysis                              -", &
               "-------------------------------------------------------------------------------"
            WRITE (log_unit, '(8X,A,T27,A,T49,A,T69,A)') "State", "Occupied", "Virtual", "Excitation"
            WRITE (log_unit, '(8X,A,T28,A,T49,A,T69,A)') "number", "orbital", "orbital", "amplitude"
            WRITE (log_unit, '(1X,79("-"))')

            IF (nspins == 1) THEN
               state_spin = 1
               state_spin2 = 2
               spin_label = '     '
               spin_label2 = '     '
            ELSE IF (spinflip .NE. no_sf_tddfpt) THEN
               state_spin = 1
               state_spin2 = 2
               spin_label = '(alp)'
               spin_label2 = '(bet)'
            END IF
         END IF

         ALLOCATE (S_mos_virt(1), weights_fm(1))
         DO ispin = 1, 1
            IF (spinflip == no_sf_tddfpt) THEN
               spin2 = ispin
            ELSE
               spin2 = 2
            END IF
            CALL cp_fm_get_info(gs_mos(spin2)%mos_virt, matrix_struct=fm_struct)
            CALL cp_fm_create(S_mos_virt(ispin), fm_struct)
            CALL cp_dbcsr_sm_fm_multiply(matrix_s, &
                                         gs_mos(spin2)%mos_virt, &
                                         S_mos_virt(ispin), &
                                         ncol=nmo_virt(spin2), alpha=1.0_dp, beta=0.0_dp)

            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_virt(spin2), ncol_global=nmo_occ(ispin), &
                                     context=blacs_env)
            CALL cp_fm_create(weights_fm(ispin), fm_struct)
            CALL cp_fm_struct_release(fm_struct)
         END DO

         nexcs_max_local = 0
         DO ispin = 1, 1
            CALL cp_fm_get_info(weights_fm(ispin), nrow_local=nrows_local, ncol_local=ncols_local)
            nexcs_max_local = nexcs_max_local + INT(nrows_local, int_8)*INT(ncols_local, int_8)
         END DO

         ALLOCATE (weights_local(nexcs_max_local), inds_local(nexcs_max_local))

         DO istate = 1, nstates
            nexcs_local = 0
            nmo_virt_occ = 0

            ! analyse matrix elements locally and transfer only significant
            ! excitations to the master node for subsequent ordering
            DO ispin = 1, 1
               IF (spinflip == no_sf_tddfpt) THEN
                  spin2 = ispin
               ELSE
                  spin2 = 2
               END IF
               ! compute excitation amplitudes
               CALL parallel_gemm('T', 'N', nmo_virt(spin2), nmo_occ(ispin), nao, 1.0_dp, S_mos_virt(ispin), &
                                  evects(istate), 0.0_dp, weights_fm(ispin))

               CALL cp_fm_get_info(weights_fm(ispin), nrow_local=nrows_local, ncol_local=ncols_local, &
                                   row_indices=row_indices, col_indices=col_indices, local_data=local_data)

               ! locate single excitations with significant amplitudes (>= min_amplitude)
               DO icol = 1, ncols_local
                  DO irow = 1, nrows_local
                     IF (ABS(local_data(irow, icol)) >= min_amplitude) THEN
                        ! number of non-negligible excitations
                        nexcs_local = nexcs_local + 1
                        ! excitation amplitude
                        weights_local(nexcs_local) = local_data(irow, icol)
                        ! index of single excitation (ivirt, iocc, ispin) in compressed form
                        inds_local(nexcs_local) = nmo_virt_occ + INT(row_indices(irow), int_8) + &
                                                  INT(col_indices(icol) - 1, int_8)*nmo_virt8(spin2)
                     END IF
                  END DO
               END DO

               nmo_virt_occ = nmo_virt_occ + nmo_virt8(spin2)*nmo_occ8(ispin)
            END DO

            IF (para_env%is_source()) THEN
               ! master node
               ALLOCATE (nexcs_recv(para_env%num_pe), recv_handlers(para_env%num_pe), recv_handlers2(para_env%num_pe))

               ! collect number of non-negligible excitations from other nodes
               DO iproc = 1, para_env%num_pe
                  IF (iproc - 1 /= para_env%mepos) THEN
                     CALL para_env%irecv(nexcs_recv(iproc:iproc), iproc - 1, recv_handlers(iproc), 0)
                  ELSE
                     nexcs_recv(iproc) = nexcs_local
                  END IF
               END DO

               DO iproc = 1, para_env%num_pe
                  IF (iproc - 1 /= para_env%mepos) &
                     CALL recv_handlers(iproc)%wait()
               END DO

               ! compute total number of non-negligible excitations
               nexcs = 0
               DO iproc = 1, para_env%num_pe
                  nexcs = nexcs + nexcs_recv(iproc)
               END DO

               ! receive indices and amplitudes of selected excitations
               ALLOCATE (weights_recv(nexcs), weights_neg_abs_recv(nexcs))
               ALLOCATE (inds_recv(nexcs), inds(nexcs))

               nmo_virt_occ = 0
               DO iproc = 1, para_env%num_pe
                  IF (nexcs_recv(iproc) > 0) THEN
                     IF (iproc - 1 /= para_env%mepos) THEN
                        ! excitation amplitudes
                        CALL para_env%irecv(weights_recv(nmo_virt_occ + 1:nmo_virt_occ + nexcs_recv(iproc)), &
                                            iproc - 1, recv_handlers(iproc), 1)
                        ! compressed indices
                        CALL para_env%irecv(inds_recv(nmo_virt_occ + 1:nmo_virt_occ + nexcs_recv(iproc)), &
                                            iproc - 1, recv_handlers2(iproc), 2)
                     ELSE
                        ! data on master node
                        weights_recv(nmo_virt_occ + 1:nmo_virt_occ + nexcs_recv(iproc)) = weights_local(1:nexcs_recv(iproc))
                        inds_recv(nmo_virt_occ + 1:nmo_virt_occ + nexcs_recv(iproc)) = inds_local(1:nexcs_recv(iproc))
                     END IF

                     nmo_virt_occ = nmo_virt_occ + nexcs_recv(iproc)
                  END IF
               END DO

               DO iproc = 1, para_env%num_pe
                  IF (iproc - 1 /= para_env%mepos .AND. nexcs_recv(iproc) > 0) THEN
                     CALL recv_handlers(iproc)%wait()
                     CALL recv_handlers2(iproc)%wait()
                  END IF
               END DO

               DEALLOCATE (nexcs_recv, recv_handlers, recv_handlers2)
            ELSE
               ! working node: send the number of selected excited states to the master node
               nexcs_send(1) = nexcs_local
               CALL para_env%isend(nexcs_send, para_env%source, send_handler, 0)
               CALL send_handler%wait()

               IF (nexcs_local > 0) THEN
                  ! send excitation amplitudes
                  CALL para_env%isend(weights_local(1:nexcs_local), para_env%source, send_handler, 1)
                  ! send compressed indices
                  CALL para_env%isend(inds_local(1:nexcs_local), para_env%source, send_handler2, 2)

                  CALL send_handler%wait()
                  CALL send_handler2%wait()
               END IF
            END IF

            ! sort non-negligible excitations on the master node according to their amplitudes,
            ! uncompress indices and print summary information
            IF (para_env%is_source() .AND. log_unit > 0) THEN
               weights_neg_abs_recv(:) = -ABS(weights_recv)
               CALL sort(weights_neg_abs_recv, INT(nexcs), inds)

               WRITE (log_unit, '(T7,I8,F10.5,A)') istate, evals(istate)*evolt, " eV"

               ! This reinitialization is needed to prevent the intel fortran compiler from introduce
               ! a bug when using optimization level 3 flag
               state_spin = 1
               state_spin2 = 1
               IF (spinflip .NE. no_sf_tddfpt) THEN
                  state_spin = 1
                  state_spin2 = 2
               END IF
               DO iexc = 1, nexcs
                  ind = inds_recv(inds(iexc)) - 1
                  IF ((nspins > 1) .AND. (spinflip == no_sf_tddfpt)) THEN
                     IF (ind < nmo_virt_occ_alpha) THEN
                        spin_label = '(alp)'
                        spin_label2 = '(alp)'
                     ELSE
                        state_spin = 2
                        state_spin2 = 2
                        ind = ind - nmo_virt_occ_alpha
                        spin_label = '(bet)'
                        spin_label2 = '(bet)'
                     END IF
                  END IF
                  imo_occ = ind/nmo_virt8(state_spin2) + 1
                  imo_virt = MOD(ind, nmo_virt8(state_spin2)) + 1

                  WRITE (log_unit, '(T27,I8,1X,A5,T48,I8,1X,A5,T70,F9.6)') imo_occ, spin_label, &
                     nmo_occ8(state_spin2) + imo_virt, spin_label2, weights_recv(inds(iexc))
               END DO
            END IF

            ! deallocate temporary arrays
            IF (para_env%is_source()) &
               DEALLOCATE (weights_recv, weights_neg_abs_recv, inds_recv, inds)
         END DO

         DEALLOCATE (weights_local, inds_local)
         IF (log_unit > 0) THEN
            WRITE (log_unit, "(1X,A)") &
               "-------------------------------------------------------------------------------"
         END IF
      END IF

      CALL cp_fm_release(weights_fm)
      CALL cp_fm_release(S_mos_virt)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_print_excitation_analysis

   FUNCTION run_converger(tddft_env, sub_env, tddfpt_control, log_unit, logger, iter) RESULT(status)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      INTEGER                                            :: log_unit
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: iter

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'run_converger'

      REAL(kind=dp)                                      :: val, rel_error_change, rel_delta_eval_change
      REAL(KIND=DP)                                      :: scal
      LOGICAL                                            :: is_converged, need_restart
      INTEGER                                            :: handle, ms, k , iv

      ! Output
      LOGICAL, DIMENSION(2) :: status

      CALL timeset(routineN, handle)

      CALL converger_run(tddft_env, log_unit, sub_env, logger, iter)

      tddft_env%max_error = maxval(tddft_env%residual_norms)
      ! Compute max energy difference between this and previous iteration
      tddft_env%max_delta_eval = maxval(abs(tddft_env%ritz_evals(1:tddft_env%nstates) &
                                - tddft_env%ritz_evals_prev(1:tddft_env%nstates)))
      tddft_env%ritz_evals_prev(1:tddft_env%nstates) = tddft_env%ritz_evals(1:tddft_env%nstates)

      ! ------------------------------------------------------------------------------
      !   Convergence, Stall, and Restart Logic for Davidson Iterations
      ! ------------------------------------------------------------------------------

      ! Check the main convergence criteria
      is_converged = ((tddft_env%max_error < tddfpt_control%conv) .OR. &
                      (tddft_env%max_delta_eval < tddfpt_control%conv*1.e-1_dp))
!     print *, "is_converged=", is_converged, "max_error=", tddft_env%max_error
!     print *, "is_converged=", is_converged, "max_delta_eval=", tddft_env%max_delta_eval
      status(1) = is_converged

      IF (is_converged) THEN
         status(2) = .FALSE.

         CALL timestop(handle)

         RETURN
      END IF

      ! Prepare for next iteration: update previous error
      tddft_env%novec = tddft_env%nvec

      ! --- Stall detection: check if residual error stagnates ---
      ! Compute relative changes for stall detection and logging
      rel_error_change      = &
         ABS(tddft_env%max_error_prev - tddft_env%max_error) / MAX(ABS(tddft_env%max_error_prev), &
         1.0E-12_dp)
      rel_delta_eval_change = &
         ABS(tddft_env%max_delta_eval_prev - tddft_env%max_delta_eval) / MAX(ABS(tddft_env%max_delta_eval_prev), &
         1.0E-12_dp)
      ! Stall detection: both error and eigenvalue change must be small
      IF (iter > 1) THEN
         IF (rel_error_change < 5.0E-2 .OR. rel_delta_eval_change < 5.0E-2) THEN
            tddft_env%stall_count = tddft_env%stall_count + 1
         ELSE
            tddft_env%stall_count = 0
         END IF
         ! Print for monitoring
!        IF (log_unit > 0) THEN
!           WRITE (log_unit, '(A,I5,A,I3,2(A,ES12.5))') ' Iter:', iter, 'stall_count:',tddft_env%stall_count, &
!                                                  '  rel_error_change=', rel_error_change, &
!                                                  '  rel_delta_eval_change=', rel_delta_eval_change
!        END IF
      END IF
      tddft_env%max_error_prev = tddft_env%max_error
      tddft_env%max_delta_eval_prev = tddft_env%max_delta_eval

      need_restart = tddft_env%stall_count >= 4 .OR. tddft_env%nvec >= tddft_env%mxvec
      status(2) = need_restart

      ! --- Handle stall, exhausted Krylov space, or too many iterations ---
      IF (need_restart) THEN
         ! Prepare to restart Davidson
         ! Update guess vectors to latest Ritz vectors and eigenvalues
         DO iv = 1, tddft_env%nstates
            CALL cp_fm_to_fm(source=tddft_env%ritz_vects(iv), &
                             destination=tddft_env%evects_mo(iv))
            tddft_env%evals(iv) = tddft_env%ritz_evals(iv)
         END DO

         tddft_env%nvec = tddft_env%nstates  ! Start with nstates vectors
         tddft_env%novec = tddft_env%nvec
         tddft_env%stall_count = 0

         ! Reset residuals
         DO iv = 1, tddft_env%nstates
            CALL cp_fm_set_all(matrix=tddft_env%residual_vecs(iv), alpha=0.0_dp)
         END DO

         CALL timestop(handle)

         RETURN
      END IF

      ! --- Krylov Expansion (only if not exiting) ---
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
            status(2) = .TRUE.

            CALL timestop(handle)

            RETURN
         END IF

         tddft_env%nvec = tddft_env%nvec + 1
         CALL cp_fm_scale_and_add( &
                  alpha=0._dp,     matrix_a=tddft_env%evects_mo(tddft_env%nvec), &
                  beta=1._dp/scal, matrix_b=tddft_env%residual_vecs(k))
      END DO

      ms = tddft_env%nvec - tddft_env%novec
      IF (ms == 0) status(2) = .TRUE.

      CALL timestop(handle)

   END FUNCTION run_converger



   SUBROUTINE tddfpt_print_summary(tddft_env, log_unit, evects, evals, ostrength, mult, &
                                   dipole_form)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      INTEGER, INTENT(in)                                :: log_unit
      TYPE(cp_fm_type), DIMENSION(:), INTENT(inout)      :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      REAL(kind=dp), DIMENSION(:), INTENT(inout)         :: ostrength
      INTEGER, INTENT(in)                                :: mult
      INTEGER, INTENT(in)                                :: dipole_form

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_print_summary'

      CHARACTER(len=1)                                   :: lsd_str
      CHARACTER(len=20)                                  :: mult_str
      INTEGER                                            :: handle, ideriv, ispin, istate, nspins, &
                                                            nstates
      REAL(kind=dp)                                      :: osc_strength
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: trans_dipoles

      CALL timeset(routineN, handle)

      ideriv = ideriv
      ispin = ispin
      nspins = tddft_env%nSpins
      nstates =tddft_env%nstates
      evects = evects

      IF (nspins > 1) THEN
         lsd_str = 'U'
      ELSE
         lsd_str = 'RO'
      END IF

      ! *** summary header ***
      IF (log_unit > 0) THEN
         CALL integer_to_string(mult, mult_str)
         WRITE (log_unit, '(/,1X,A1,A,1X,A)') lsd_str, "-TDDFPT states of multiplicity", TRIM(mult_str)
         SELECT CASE (dipole_form)
         CASE (tddfpt_dipole_berry)
            WRITE (log_unit, '(1X,A,/)') "Transition dipoles calculated using Berry operator formulation"
         CASE (tddfpt_dipole_length)
            WRITE (log_unit, '(1X,A,/)') "Transition dipoles calculated using length formulation"
         CASE (tddfpt_dipole_velocity)
            WRITE (log_unit, '(1X,A,/)') "Transition dipoles calculated using velocity formulation"
         CASE DEFAULT
            CPABORT("Unimplemented form of the dipole operator")
         END SELECT

         WRITE (log_unit, '(T10,A,T19,A,T37,A,T69,A)') "State", "Excitation", &
            "Transition dipole (a.u.)", "Oscillator"
         WRITE (log_unit, '(T10,A,T19,A,T37,A,T49,A,T61,A,T67,A)') "number", "energy (eV)", &
            "x", "y", "z", "strength (a.u.)"
         WRITE (log_unit, '(T10,72("-"))')
      END IF

      ! transition dipole moment
      ALLOCATE (trans_dipoles(nstates, nderivs, nspins))
      trans_dipoles(:, :, :) = 0.0_dp

      ! nspins == 1 .AND. mult == 3 : spin-flip transitions are forbidden due to symmetry reasons
!     IF (nspins > 1 .OR. mult == 1) THEN
!        DO ispin = 1, nspins
!           DO ideriv = 1, nderivs
!              CALL cp_fm_trace(evects(:), dipole_op_mos_occ(ideriv, ispin), &
!                               trans_dipoles(:, ideriv, ispin))
!           END DO
!        END DO
!
!        IF (nspins == 1) THEN
!           trans_dipoles(:, :, 1) = SQRT(2.0_dp)*trans_dipoles(:, :, 1)
!        ELSE
!           trans_dipoles(:, :, 1) = SQRT(trans_dipoles(:, :, 1)**2 + trans_dipoles(:, :, 2)**2)
!        END IF
!     END IF

      ! *** summary information ***
      DO istate = 1, nstates
         osc_strength = 2.0_dp/3.0_dp*evals(istate)* &
                        accurate_dot_product(trans_dipoles(istate, :, 1), trans_dipoles(istate, :, 1))
         ostrength(istate) = osc_strength
         IF (log_unit > 0) THEN
            WRITE (log_unit, '(1X,A,T9,I7,T19,F11.5,T31,3(1X,ES11.4E2),T69,ES12.5E2)') &
               "TDDFT| ", istate, evals(istate)*evolt, trans_dipoles(istate, 1:nderivs, 1), osc_strength
         END IF
      END DO

      ! punch a checksum for the regs
      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T2,A,E14.6)') 'TDDFT : CheckSum  =', SQRT(SUM(evals**2))
      END IF

      DEALLOCATE (trans_dipoles)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_print_summary

   SUBROUTINE converger_run(tddft_env, log_unit, sub_env, logger, iter)
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      INTEGER                                            :: log_unit
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: iter

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'converger_run'

      INTEGER                                            :: handle, i, j, a, ij, ierr, iv
      REAL(KIND=DP)                                      :: val, sig, scal, val2

      CALL timeset(routineN, handle)

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
      ierr = 0
      CALL choose_eigv_solver(matrix=tddft_env%subspace, &
                              eigenvalues=tddft_env%ritz_evals(1:tddft_env%nvec), &
                              eigenvectors=tddft_env%subspace_evecs, &
                              info=ierr)
      ! Check ierr where?

      CALL deallocate_fm_matrix(tddft_env%subspace)
      ! Print
!     print *, "KKK| ================ Energy ============== at iter=", iter
!     DO iv = 1, tddft_env%nstates
!        print *, "KKK| Energy (eV):", iv, tddft_env%ritz_evals(iv)*evolt
!     END DO
!     CALL print_mat_ao_ao(tddft_env%subspace_evecs, "subspace_evecs", logger, iter)

      ! Normalization
      CALL cp_fm_vectorsnorm(tddft_env%subspace_evecs, tddft_env%vector_norm(1:tddft_env%nvec))
      tddft_env%vector_scal(1:tddft_env%nvec) = 1.0_dp/tddft_env%vector_norm(1:tddft_env%nvec)
      CALL cp_fm_column_scale(tddft_env%subspace_evecs, tddft_env%vector_scal(1:tddft_env%nvec))

      ! Print
      CALL print_mat_ao_ao(tddft_env%subspace_evecs, "normed subspace_evecs", logger, iter)

!     transformation_matrix = subspace_evecs
      DO iv = 1, tddft_env%nstates
         CALL cp_fm_set_all(matrix=tddft_env%ritz_vects(iv), alpha=0.0_dp)
         DO j = 1, tddft_env%nvec
             CALL cp_fm_get_element(matrix=tddft_env%subspace_evecs, &
                                    irow_global=j, &
                                    icol_global=iv, &
                                    alpha=val)
             CALL cp_fm_scale_and_add(&
                     alpha=1.0_dp, matrix_a=tddft_env%ritz_vects(iv), &
                     beta=val, matrix_b=tddft_env%evects_mo(j))
         END DO
      END DO

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

      CALL timestop(handle)
   END SUBROUTINE converger_run

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
      IF (.FALSE.) THEN
         IF (unit_nr > 0) THEN
            WRITE(unit_nr, '(A,A,A,I0,A)') "First 5x5 elements of ",TRIM(mat_name)," at iter ",iter,":"
            DO i = 1, 5
               WRITE(unit_nr, '(5F12.6)') (submatrix(i,j), j=1,5)
            END DO
            WRITE(unit_nr, '(A)') ""
            CALL m_flush(unit_nr)
         END IF
      END IF

      CALL timestop(handle)
   END SUBROUTINE print_mat_ao_ao

END MODULE qs_tddft_methods
