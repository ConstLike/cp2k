!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!
MODULE qs_tddft_theory_mrsf

   USE util,                            ONLY: sort
   USE string_utilities,                ONLY: integer_to_string
   USE input_constants,                 ONLY: &
        do_admm_aux_exch_func_none, do_admm_basis_projection, do_admm_exch_scaling_none, &
        do_admm_purify_none, do_potential_truncated, tddfpt_dipole_velocity, tddfpt_kernel_full, &
        tddfpt_kernel_none, tddfpt_kernel_stda,&
                                              no_sf_tddfpt,&
                                              tddfpt_sf_col,&
                                              tddfpt_sf_noncol
   USE input_constants,                 ONLY: tddfpt_dipole_berry,&
                                              tddfpt_dipole_length
   USE kahan_sum,                       ONLY: accurate_dot_product
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_get,&
                                              qs_rho_type
   USE qs_tddfpt2_densities,            ONLY: tddfpt_construct_aux_fit_density,&
                                              tddfpt_construct_ground_state_orb_density
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_dbcsr_create_by_dist,&
                                              tddfpt_subgroup_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_zero
   USE pw_env_types,                    ONLY: pw_env_get
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_p_type,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE qs_tddft_theory,                 ONLY: tddft_theory_type
   USE qs_tddft_types,                  ONLY: tddft_info_type
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE cp_control_types,                ONLY: tddfpt2_control_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE qs_kernel_types,                 ONLY: kernel_env_type
   USE kinds,                           ONLY: dp, sp
   USE physcon,                         ONLY: evolt
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_gemm,&
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
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE qs_kernel_types,                 ONLY: full_kernel_env_type,&
                                              kernel_env_type,&
                                              release_kernel_env
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_theory_mrsf'

   PUBLIC :: tddft_mrsf_theory_type

   TYPE, EXTENDS(tddft_theory_type) :: tddft_mrsf_theory_type

      !> Hamiltonian matrix in AO basis
      TYPE(cp_fm_type), POINTER, DIMENSION(:) :: hmat_ao_ao => NULL()

      !> AO trial vectors (size is 1)
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: evects_ao_ao

      !> AO–MO transformation buffer [nbf, nvirb]
      TYPE(cp_fm_type), POINTER :: work_ao_vir => NULL()
      !> MO–AO transformation buffer [nocca, nbf]
      TYPE(cp_fm_type), POINTER :: work_noc_ao => NULL()
      !> Work buffer [nocca, nocca]
      TYPE(cp_fm_type), POINTER :: work_noc_noc => NULL()
      !> Work buffer [nbf, nocca]
      TYPE(cp_fm_type), POINTER :: work_ao_noc => NULL()
      !> Work buffer [nbf, nocca]
      TYPE(cp_fm_type), POINTER :: work_vir_vir => NULL()
      !> Scratch matrix [nocca, nvirb], used for Fock response, dipoles, etc.
      TYPE(cp_fm_type), POINTER :: work_noc_vir => NULL()
      !> Scratch matrices used in AO contraction [nbf, nbf] (sym, asym)
      TYPE(cp_fm_type), POINTER :: work1_ao_ao => NULL()
      TYPE(cp_fm_type), POINTER :: work2_ao_ao => NULL()

      TYPE(cp_fm_type), POINTER :: work1_ao_1 => NULL()
      TYPE(cp_fm_type), POINTER :: work2_ao_1 => NULL()
      TYPE(cp_fm_type), POINTER :: work_vir_1 => NULL()
      TYPE(cp_fm_type), POINTER :: work_noc_1 => NULL()
      TYPE(cp_fm_type), POINTER :: work_vir_noc => NULL()
      TYPE(cp_fm_type), POINTER :: work_ao_nocb => NULL()

      !-------------------------------
      ! KERNELs related
      !-------------------------------
      !> [XC kernel] NO NAME
      TYPE(cp_fm_type), POINTER :: rho_ia_ao_fm => NULL()
      TYPE(qs_rho_type), POINTER :: gs_rho_xc_struct => NULL()
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:) :: A_ia_munu_sub
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER :: A_ia_rspace_sub => NULL()
      TYPE(qs_rho_type), POINTER :: rho_orb_struct_sub => NULL()
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:) :: rho_ia_ao
      TYPE(cp_fm_type), POINTER, DIMENSION(:) :: rho_fm_ao => NULL()

      !> [HFX kernel] Symmetric AO density matrices [nSpins]
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:) :: work_rho_ao_ao_symm

      !> [HFX kernel] Antisymmetric AO density matrices [nSpins]
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:) :: work_rho_ao_ao_asymm

      !> [HFX kernel] Symmetric Hamiltonian-like matrices in AO space [nSpins]
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:) :: work_hmat_symm

      !> [HFX kernel] Antisymmetric Hamiltonian-like matrices in AO space [nSpins]
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:) :: work_hmat_asymm

      !> [Coulomb kernel] Density on grid in real-space [nSpins]
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER :: rho_r => NULL()

      !> [Coulomb kernel] Density in G-space [nSpins]
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER :: rho_g => NULL()

      !> [.. kernel] Kinetic energy density (tau) on grid [nSpins]
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER :: tau_r => NULL()

      !> Exchange-correlation potential v_XC(r) [nSpins]
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER :: work_v_xc => NULL()
      !> tau-dependent v_XC components (meta-GGA) [nSpins]
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER :: work_v_xc_tau => NULL()

      ! -- theory-specific fields if needed --
   CONTAINS
      PROCEDURE, PASS(this) :: allocate_theory_data     => mrsf_allocate_theory_data
      PROCEDURE, PASS(this) :: deallocate_theory_data   => mrsf_deallocate_theory_data
      PROCEDURE, PASS(this) :: initialize_guess_vectors => mrsf_initialize_guess_vectors
      PROCEDURE, PASS(this) :: mo_to_ao_transform       => mrsf_mo_to_ao_transform
      PROCEDURE, PASS(this) :: apply_operator           => mrsf_apply_operator
      PROCEDURE, PASS(this) :: ao_to_mo_transform       => mrsf_ao_to_mo_transform
      PROCEDURE, PASS(this) :: print_theory_data        => mrsf_print_theory_data
   END TYPE

CONTAINS

   SUBROUTINE mrsf_allocate_theory_data(this, qs_env, tddft_env, sub_env)
      CLASS(tddft_mrsf_theory_type), INTENT(INOUT) :: this
      TYPE(qs_environment_type), POINTER :: qs_env
      TYPE(tddft_info_type), INTENT(INOUT) :: tddft_env
      TYPE(tddfpt_subgroup_env_type), INTENT(IN) :: sub_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mrsf_allocate_theory_data'

      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_hfx
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      INTEGER                                            :: handle, ivec, ispin, imat

      CALL timeset(routineN, handle)

      SELECT CASE (tddft_env%response_type)
      CASE(2) ! MRSF-S
         tddft_env%nStates = min(tddft_env%nstates, tddft_env%nstates_max-1)
         tddft_env%mxvec = min(tddft_env%maxvec*tddft_env%nstates, &
                               tddft_env%nstates_max-1)
      CASE(3) ! MRSF-T
         tddft_env%nStates = min(tddft_env%nstates, tddft_env%nstates_max-3)
         tddft_env%mxvec = min(tddft_env%maxvec*tddft_env%nstates, &
                               tddft_env%nstates_max-3)
      END SELECT

      !=== hmat_ao_ao [nbf, nbf] =======================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%hmat_ao_ao(7))
      DO imat = 1, 7
         CALL cp_fm_create(this%hmat_ao_ao(imat), fm_struct)
         CALL cp_fm_set_all(this%hmat_ao_ao(imat), 0.0_dp)
      END DO
      CALL cp_fm_struct_release(fm_struct)
      !=================================================
      !=== evects_ao_ao (7) [nbf, nbf] =============
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%evects_ao_ao(7))
      DO imat = 1, 7
         CALL cp_fm_create(this%evects_ao_ao(imat), fm_struct)
         CALL cp_fm_set_all(this%evects_ao_ao(imat), 0.0_dp)
      END DO
      CALL cp_fm_struct_release(fm_struct)
      !=================================================
      !=== work_ao_vir [nbf, nvirb] ====================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nvirb, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_ao_vir)
      CALL cp_fm_create(this%work_ao_vir, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_ao_vir, 0.0_dp)
      !=================================================
      !=== work_noc_ao [nocca, nbf] ====================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nocca, &
                               ncol_global=tddft_env%nbf, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_noc_ao)
      CALL cp_fm_create(this%work_noc_ao, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_noc_ao, 0.0_dp)
      !=================================================
      !=== work_noc_noc [nocca, nocca] =================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nocca, &
                               ncol_global=tddft_env%nocca, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_noc_noc)
      CALL cp_fm_create(this%work_noc_noc, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_noc_noc, 0.0_dp)
      !=================================================
      !=== work_ao_noc [nbf, nocca] ====================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nocca, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_ao_noc)
      CALL cp_fm_create(this%work_ao_noc, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_ao_noc, 0.0_dp)
      !=================================================
      !=== work_vir_vir [nvirb, nvirb] =================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nvirb, &
                               ncol_global=tddft_env%nvirb, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_vir_vir)
      CALL cp_fm_create(this%work_vir_vir, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_vir_vir, 0.0_dp)
      !=================================================
      !=== work_noc_vir [nocca, nvirb] =================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nocca, &
                               ncol_global=tddft_env%nvirb, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_noc_vir)
      CALL cp_fm_create(this%work_noc_vir, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_noc_vir, 0.0_dp)
      !=================================================
      !=== work1_ao_ao [nbf, nbf] ======================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work1_ao_ao)
      CALL cp_fm_create(this%work1_ao_ao, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work1_ao_ao, 0.0_dp)
      !=================================================
      !=== work2_ao_ao [nbf, nbf] ======================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work2_ao_ao)
      CALL cp_fm_create(this%work2_ao_ao, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work2_ao_ao, 0.0_dp)
      !=================================================
      !=== work1_ao_1 [nbf, 1] =========================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=1, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work1_ao_1)
      CALL cp_fm_create(this%work1_ao_1, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work1_ao_1, 0.0_dp)
      !=================================================
      !=== work2_ao_1 [nbf, 1] =========================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=1, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work2_ao_1)
      CALL cp_fm_create(this%work2_ao_1, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work2_ao_1, 0.0_dp)
      !=================================================
      !=== work_vir_1 [nvira, 1] =========================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nvira, &
                               ncol_global=1, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_vir_1)
      CALL cp_fm_create(this%work_vir_1, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_vir_1, 0.0_dp)
      !=================================================
      !=== work_noc_1 [noccb, 1] =========================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%noccb, &
                               ncol_global=1, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_noc_1)
      CALL cp_fm_create(this%work_noc_1, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_noc_1, 0.0_dp)
      !=================================================
      !=== work_vir_noc [nvira, noccb] =========================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nvira, &
                               ncol_global=tddft_env%noccb, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_vir_noc)
      CALL cp_fm_create(this%work_vir_noc, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_vir_noc, 0.0_dp)
      !=================================================
      !=== work_ao_nocb [nbf, noccb] ====================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%noccb, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%work_ao_nocb)
      CALL cp_fm_create(this%work_ao_nocb, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_ao_nocb, 0.0_dp)
      !=================================================
      !=== rho_ia_ao_fm [nbf, nbf] =====================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%rho_ia_ao_fm)
      CALL cp_fm_create(this%rho_ia_ao_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%rho_ia_ao_fm, 0.0_dp)
      !=== rho_fm_ao [nbf, nbf] ========================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%rho_fm_ao(tddft_env%nSpins))
      DO ispin = 1, tddft_env%nSpins
         CALL cp_fm_create(this%rho_fm_ao(ispin), fm_struct)
         CALL cp_fm_set_all(this%rho_fm_ao(ispin), 0.0_dp)
      END DO
      CALL cp_fm_struct_release(fm_struct)
      !=================================================
      !=== dipoles_mo(3) [nocca, nvirb] ================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nocca, &
                               ncol_global=tddft_env%nvirb, &
                               context=sub_env%blacs_env)
      ALLOCATE(this%dipoles_mo(3))
      DO ivec = 1, 3
         CALL cp_fm_create(this%dipoles_mo(ivec), fm_struct)
         CALL cp_fm_set_all(this%dipoles_mo(ivec), 0.0_dp)
      END DO
      CALL cp_fm_struct_release(fm_struct)
      !=================================================

      !--------------------
      ! ALLOCATE  DBCSR
      !--------------------
      CALL get_qs_env(qs_env=qs_env, &
                      matrix_s=matrix_s, &
                      dbcsr_dist=sub_env%dbcsr_dist, &
                      sab_orb=sab_hfx)

      NULLIFY (this%A_ia_munu_sub)
      CALL dbcsr_allocate_matrix_set(this%A_ia_munu_sub, 1)
      DO ispin = 1, 1
         CALL dbcsr_init_p(this%A_ia_munu_sub(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(this%A_ia_munu_sub(ispin)%matrix, &
                                          template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, &
                                          sab=sub_env%sab_orb)
      END DO

      NULLIFY (this%work_rho_ao_ao_symm)
      CALL dbcsr_allocate_matrix_set(this%work_rho_ao_ao_symm, tddft_env%nSpins)
      NULLIFY (this%work_rho_ao_ao_asymm)
      CALL dbcsr_allocate_matrix_set(this%work_rho_ao_ao_asymm, tddft_env%nSpins)
      DO ispin = 1, tddft_env%nSpins
         CALL dbcsr_init_p(this%work_rho_ao_ao_symm(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist( &
                  this%work_rho_ao_ao_symm(ispin)%matrix, &
                  template=matrix_s(1)%matrix, &
                  dbcsr_dist=sub_env%dbcsr_dist, &
                  sab=sab_hfx)

         CALL dbcsr_init_p(this%work_rho_ao_ao_asymm(ispin)%matrix)
         CALL dbcsr_create(this%work_rho_ao_ao_asymm(ispin)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric, &
                           template=this%work_rho_ao_ao_symm(ispin)%matrix)
         CALL dbcsr_complete_redistribute(this%work_rho_ao_ao_symm(ispin)%matrix, &
                                          this%work_rho_ao_ao_asymm(ispin)%matrix)
      END DO

      NULLIFY (this%work_hmat_symm)
      CALL dbcsr_allocate_matrix_set(this%work_hmat_symm, tddft_env%nSpins)
      NULLIFY (this%work_hmat_asymm)
      CALL dbcsr_allocate_matrix_set(this%work_hmat_asymm, tddft_env%nSpins)
      DO ispin = 1, tddft_env%nSpins
         CALL dbcsr_init_p(this%work_hmat_symm(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(this%work_hmat_symm(ispin)%matrix, &
                                          template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, &
                                          sab=sab_hfx)

         CALL dbcsr_init_p(this%work_hmat_asymm(ispin)%matrix)
         CALL dbcsr_create(this%work_hmat_asymm(ispin)%matrix, matrix_type=dbcsr_type_antisymmetric, &
                           template=this%work_hmat_symm(ispin)%matrix)
         CALL dbcsr_complete_redistribute(this%work_hmat_symm(ispin)%matrix, &
                                          this%work_hmat_asymm(ispin)%matrix)
      END DO

      NULLIFY (this%rho_ia_ao)
      CALL dbcsr_allocate_matrix_set(this%rho_ia_ao, tddft_env%nSpins)
      DO ispin = 1, tddft_env%nSpins
         CALL dbcsr_init_p(this%rho_ia_ao(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(matrix=this%rho_ia_ao(ispin)%matrix, &
                                          template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, &
                                          sab=sub_env%sab_orb)
      END DO
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

      ALLOCATE (this%A_ia_rspace_sub(tddft_env%nSpins))
      DO ispin = 1, tddft_env%nSpins
         CALL auxbas_pw_pool%create_pw(this%A_ia_rspace_sub(ispin))
      END DO
      ALLOCATE (this%rho_r(tddft_env%nspins))
      DO ispin = 1, tddft_env%nspins
         CALL auxbas_pw_pool%create_pw(this%rho_r(ispin))
      END DO
      ALLOCATE (this%rho_g(tddft_env%nspins))
      DO ispin = 1, tddft_env%nspins
         CALL auxbas_pw_pool%create_pw(this%rho_g(ispin))
      END DO
      ALLOCATE (this%tau_r(tddft_env%nspins))
      DO ispin = 1, tddft_env%nspins
         CALL auxbas_pw_pool%create_pw(this%tau_r(ispin))
      END DO
      ALLOCATE (this%work_v_xc(tddft_env%nspins))
      DO ispin = 1, tddft_env%nspins
         CALL auxbas_pw_pool%create_pw(this%work_v_xc(ispin))
      END DO
      ALLOCATE (this%work_v_xc_tau(tddft_env%nspins))
      DO ispin = 1, tddft_env%nspins
         CALL auxbas_pw_pool%create_pw(this%work_v_xc_tau(ispin))
      END DO

!     ! Full kernel
!     CALL create_kernel_env(kernel_env=full_kernel_env, &
!                            rho_struct_sub=this%rho_orb_struct_sub, &
!                            xc_section=xc_section, &
!                            is_rks_triplets=tddfpt_control%rks_triplets, &
!                            sub_env=sub_env)
!     kernel_env%full_kernel => full_kernel_env
!     NULLIFY (kernel_env%admm_kernel)

      CALL timestop(handle)

   END SUBROUTINE mrsf_allocate_theory_data

   SUBROUTINE mrsf_deallocate_theory_data(this, sub_env)
      CLASS(tddft_mrsf_theory_type), INTENT(INOUT) :: this
      TYPE(tddfpt_subgroup_env_type), INTENT(IN) :: sub_env

      TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
      INTEGER :: ivec, ispin, imat

      IF (ASSOCIATED(this%hmat_ao_ao)) THEN
         DO imat = 1, SIZE(this%hmat_ao_ao)
            CALL cp_fm_release(this%hmat_ao_ao(imat))
         END DO
         DEALLOCATE(this%hmat_ao_ao)
      END IF
      IF (ALLOCATED(this%evects_ao_ao)) THEN
         DO imat = 1, SIZE(this%evects_ao_ao)
            CALL cp_fm_release(this%evects_ao_ao(imat))
         END DO
         DEALLOCATE(this%evects_ao_ao)
      END IF
      IF (ASSOCIATED(this%work_ao_vir)) THEN
         CALL cp_fm_release(this%work_ao_vir)
         DEALLOCATE(this%work_ao_vir)
      END IF
      IF (ASSOCIATED(this%work_noc_ao)) THEN
         CALL cp_fm_release(this%work_noc_ao)
         DEALLOCATE(this%work_noc_ao)
      END IF
      IF (ASSOCIATED(this%work_noc_noc)) THEN
         CALL cp_fm_release(this%work_noc_noc)
         DEALLOCATE(this%work_noc_noc)
      END IF
      IF (ASSOCIATED(this%work_ao_noc)) THEN
         CALL cp_fm_release(this%work_ao_noc)
         DEALLOCATE(this%work_ao_noc)
      END IF
      IF (ASSOCIATED(this%work_vir_vir)) THEN
         CALL cp_fm_release(this%work_vir_vir)
         DEALLOCATE(this%work_vir_vir)
      END IF
      IF (ASSOCIATED(this%work_noc_vir)) THEN
         CALL cp_fm_release(this%work_noc_vir)
         DEALLOCATE(this%work_noc_vir)
      END IF
      IF (ASSOCIATED(this%work1_ao_ao)) THEN
         CALL cp_fm_release(this%work1_ao_ao)
         DEALLOCATE(this%work1_ao_ao)
      END IF
      IF (ASSOCIATED(this%work2_ao_ao)) THEN
         CALL cp_fm_release(this%work2_ao_ao)
         DEALLOCATE(this%work2_ao_ao)
      END IF
      IF (ASSOCIATED(this%work1_ao_1)) THEN
         CALL cp_fm_release(this%work1_ao_1)
         DEALLOCATE(this%work1_ao_1)
      END IF
      IF (ASSOCIATED(this%work2_ao_1)) THEN
         CALL cp_fm_release(this%work2_ao_1)
         DEALLOCATE(this%work2_ao_1)
      END IF
      IF (ASSOCIATED(this%work_vir_1)) THEN
         CALL cp_fm_release(this%work_vir_1)
         DEALLOCATE(this%work_vir_1)
      END IF
      IF (ASSOCIATED(this%work_noc_1)) THEN
         CALL cp_fm_release(this%work_noc_1)
         DEALLOCATE(this%work_noc_1)
      END IF
      IF (ASSOCIATED(this%work_vir_noc)) THEN
         CALL cp_fm_release(this%work_vir_noc)
         DEALLOCATE(this%work_vir_noc)
      END IF
      IF (ASSOCIATED(this%work_ao_nocb)) THEN
         CALL cp_fm_release(this%work_ao_nocb)
         DEALLOCATE(this%work_ao_nocb)
      END IF
      IF (ASSOCIATED(this%rho_ia_ao_fm)) THEN
         CALL cp_fm_release(this%rho_ia_ao_fm)
         DEALLOCATE(this%rho_ia_ao_fm)
      END IF
      IF (ASSOCIATED(this%rho_fm_ao)) THEN
         DO ispin = 1, SIZE(this%rho_fm_ao)
            CALL cp_fm_release(this%rho_fm_ao(ispin))
         END DO
         DEALLOCATE(this%rho_fm_ao)
      END IF
      IF (ALLOCATED(this%dipoles_mo)) THEN
         DO ivec = 1, SIZE(this%dipoles_mo)
            CALL cp_fm_release(this%dipoles_mo)
         END DO
!        DEALLOCATE(this%dipoles_mo)
      END IF

      IF (ASSOCIATED(this%A_ia_munu_sub)) THEN
         DO ivec = 1, SIZE(this%A_ia_munu_sub)
            CALL dbcsr_deallocate_matrix(this%A_ia_munu_sub(ivec)%matrix)
         END DO
         DEALLOCATE(this%A_ia_munu_sub)
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

      CALL pw_env_get(sub_env%pw_env, auxbas_pw_pool=auxbas_pw_pool)

      IF (ASSOCIATED(this%A_ia_rspace_sub)) THEN
         DO ispin = 1, SIZE(this%A_ia_rspace_sub)
            CALL auxbas_pw_pool%give_back_pw(this%A_ia_rspace_sub(ispin))
         END DO
         DEALLOCATE(this%A_ia_rspace_sub)
      END IF
      IF (ASSOCIATED(this%rho_r)) THEN
         DO ispin = 1, SIZE(this%rho_r)
            CALL auxbas_pw_pool%give_back_pw(this%rho_r(ispin))
         END DO
         DEALLOCATE(this%rho_r)
      END IF
      IF (ASSOCIATED(this%rho_g)) THEN
         DO ispin = 1, SIZE(this%rho_g)
            CALL auxbas_pw_pool%give_back_pw(this%rho_g(ispin))
         END DO
         DEALLOCATE(this%rho_g)
      END IF
      IF (ASSOCIATED(this%tau_r)) THEN
         DO ispin = 1, SIZE(this%tau_r)
            CALL auxbas_pw_pool%give_back_pw(this%tau_r(ispin))
         END DO
         DEALLOCATE(this%tau_r)
      END IF
      IF (ASSOCIATED(this%work_v_xc)) THEN
         DO ispin = 1, SIZE(this%work_v_xc)
            CALL auxbas_pw_pool%give_back_pw(this%work_v_xc(ispin))
         END DO
         DEALLOCATE(this%work_v_xc)
      END IF
      IF (ASSOCIATED(this%work_v_xc_tau)) THEN
         DO ispin = 1, SIZE(this%work_v_xc_tau)
            CALL auxbas_pw_pool%give_back_pw(this%work_v_xc_tau(ispin))
         END DO
         DEALLOCATE(this%work_v_xc_tau)
      END IF
      IF (ASSOCIATED(this%rho_ia_ao)) THEN
         DO ivec = 1, SIZE(this%rho_ia_ao)
            CALL dbcsr_deallocate_matrix(this%rho_ia_ao(ivec)%matrix)
         END DO
         DEALLOCATE(this%rho_ia_ao)
      END IF

   END SUBROUTINE mrsf_deallocate_theory_data

   SUBROUTINE mrsf_initialize_guess_vectors(this, qs_env, sub_env, tddft_env, gs_mos, &
                                          tddfpt_control, log_unit)
      CLASS(tddft_mrsf_theory_type), INTENT(INOUT)         :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(IN)         :: sub_env
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), POINTER                           :: gs_mos
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      INTEGER                                            :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mrsf_initialize_guess_vectors'

      INTEGER                                            :: handle, imo_occ, imo_virt, &
                                                            istate, lr1, lr2, &
                                                            nStates_selected, istate_orig
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      LOGICAL                                             :: do_triplets
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_m_occ

      CALL timeset(routineN, handle)

      CALL tddfpt_construct_ground_state_orb_density( &
               rho_orb_struct=this%rho_orb_struct_sub, &
               rho_xc_struct=this%gs_rho_xc_struct, &
               is_rks_triplets=tddfpt_control%rks_triplets, &
               qs_env=qs_env, &
               sub_env=sub_env, &
               wfm_rho_orb=this%rho_ia_ao_fm)


      IF (.FALSE.) THEN
         print *, this%unused
         NULLIFY(tddfpt_control)
      END IF

      nStates_selected = gs_mos(1)%homo*(SIZE(gs_mos(2)%evals_full)-gs_mos(2)%homo)

      ALLOCATE (inds(nStates_selected))
      ALLOCATE (e_virt_m_occ(nStates_selected))

      lr1 = tddft_env%noccb+1
      lr2 = tddft_env%nocca
      do_triplets = tddft_env%response_type == 3
      istate = 0
      DO imo_virt = gs_mos(2)%homo+1, SIZE(gs_mos(2)%evals_full)
         DO imo_occ = 1, gs_mos(1)%homo
            istate = (imo_virt-gs_mos(2)%homo-1)*gs_mos(1)%homo+imo_occ
            IF (imo_occ==lr1 .AND. imo_virt==lr1) THEN
               tddft_env%e_virt_minus_occ(istate) = (gs_mos(2)%evals_full(lr1)-gs_mos(1)%evals_full(lr1) &
                                                   + gs_mos(2)%evals_full(lr2)-gs_mos(1)%evals_full(lr2))*0.5_dp
               CYCLE
             END IF

             tddft_env%e_virt_minus_occ(istate) = gs_mos(2)%evals_full(imo_virt)-gs_mos(1)%evals_full(imo_occ)

             if (imo_occ==lr2 .and. imo_virt==lr2) THEN
               tddft_env%e_virt_minus_occ(istate) = huge(1.0_dp)
             ELSE IF (do_triplets .and. imo_occ==lr2 .and. imo_virt==lr1) THEN
               tddft_env%e_virt_minus_occ(istate) = huge(1.0_dp)
             ELSE IF (do_triplets .and. imo_occ==lr1 .and. imo_virt==lr2) THEN
               tddft_env%e_virt_minus_occ(istate) = huge(1.0_dp)
             END IF
         END DO
      END DO

      e_virt_m_occ(1:nStates_selected) = tddft_env%e_virt_minus_occ(1:nStates_selected)

      CALL sort(e_virt_m_occ, nStates_selected, inds)

      IF (log_unit > 0) THEN
         WRITE (log_unit, "(1X,A)") "", &
            "-------------------------------------------------------------------------------", &
            "-                         MRSF-TDDFT Initial Guess                            -", &
            "-------------------------------------------------------------------------------"
         WRITE (log_unit, '(T11,A)') "State         Occupied      ->      Virtual          Excitation"
         WRITE (log_unit, '(T11,A)') "number         orbital              orbital          energy (eV)"
         WRITE (log_unit, '(1X,79("-"))')
      END IF

      DO istate = 1, tddft_env%nVec
         tddft_env%evals(istate) = e_virt_m_occ(istate)
         istate_orig = inds(istate)
         imo_occ = MOD(istate_orig - 1, gs_mos(1)%homo) + 1
         imo_virt = (istate_orig - 1) / gs_mos(1)%homo + 1
         CALL cp_fm_set_element(tddft_env%evects_mo(istate), imo_occ, imo_virt, 1.0_dp)
         IF (istate > tddft_env%nStates) CYCLE
         IF (log_unit > 0) &
            WRITE (log_unit, '(T7,I8,T24,I8,T37,A5,T45,I8,T54,A5,T60,2F14.5)') &
            istate, imo_occ, '(adp)', gs_mos(2)%homo + imo_virt, '(adp)', e_virt_m_occ(istate)*evolt, &
            tddft_env%evects_mo(istate)%local_data(imo_occ, imo_virt)
      end do

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T7,A,T50,I24)') 'Number of active states:', &
            SIZE(gs_mos(1)%evals_occ)*SIZE(gs_mos(2)%evals_virt)
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      DEALLOCATE (inds, e_virt_m_occ)

      CALL timestop(handle)

   END SUBROUTINE mrsf_initialize_guess_vectors

   SUBROUTINE mrsf_mo_to_ao_transform(this, tddft_env, gs_mos, ivec, logger)
      CLASS(tddft_mrsf_theory_type), TARGET, INTENT(INOUT) :: this
      TYPE(tddft_info_type), INTENT(INOUT) :: tddft_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), POINTER :: gs_mos
      INTEGER, INTENT(IN) :: ivec
      TYPE(cp_logger_type), POINTER :: logger

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mrsf_mo_to_ao_transform'

      INTEGER :: handle, i, j, m, lr1, lr2, nbf, nocca, noccb, nvira, nvirb
      TYPE(cp_fm_type), POINTER :: ball, bo2v, bo1v, bco1, bco2, o21v, co12
      TYPE(cp_fm_type), POINTER :: va, vb
      REAL(KIND=dp) :: scalar
      INTEGER :: shift

      CALL timeset(routineN, handle)

      logger = logger
      m = m
      j = j
      i = i

      ! Extract dimensions from tddft_env
      nbf = tddft_env%nbf
      nocca = tddft_env%nocca
      noccb = tddft_env%noccb
      nvira = tddft_env%nvira  ! nvira = nbf - nocca
      nvirb = tddft_env%nvirb  ! nvirb = nbf - noccb

      ! Define indices
      lr1 = nocca - 1
      lr2 = nocca
      shift = nocca - noccb

      ! Assign pointers from tddft_env%evects_ao_ao [nbf, nbf]
      bo2v => this%evects_ao_ao(1)
      bo1v => this%evects_ao_ao(2)
      bco1 => this%evects_ao_ao(3)
      bco2 => this%evects_ao_ao(4)
      o21v => this%evects_ao_ao(5)
      co12 => this%evects_ao_ao(6)
      ball => this%evects_ao_ao(7)

      va => gs_mos(1)%mos_full           ! [nbf, nbf]
      vb => gs_mos(2)%mos_full           ! [nbf, nbf]

      ! transform needed row to coulumn, not a transpose.
      ! nocca = noccb + 2
      ! nvira = nvirb - 2
      ! work[1:nvira, 1] => X(1:nocca, 1:nvirb) [nocca, 1+shift:nvirb]
      CALL cp_fm_set_all(this%work_vir_1, 0.0_dp)
      DO i = 1, nvira
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), lr2, i+shift, scalar)
         CALL cp_fm_set_element(this%work_vir_1, i, 1, scalar)
      END DO

      ! tmp (1:nvira, 1) = vb(:,nocca+1:nbf)*bvec(lr2, nocca+1:nbf)
      CALL parallel_gemm('N', 'N', nbf, 1, nvira, &
                         1.0_dp, gs_mos(2)%mos_full, &
                                 this%work_vir_1, &
                         0.0_dp, this%work2_ao_1, &
                         a_first_col=nocca+1)
      ! 1: bo2v += va(:,lr2) * tmp^T
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                         1.0_dp, gs_mos(1)%mos_full, &
                                 this%work2_ao_1, &
                         1.0_dp, this%evects_ao_ao(1), & !
                         a_first_col=lr2)
      ! 5: o21v += tmp * va(:,lr1)^T
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                         1.0_dp, this%work2_ao_1, &
                                 gs_mos(1)%mos_full, &
                         1.0_dp, this%evects_ao_ao(5), &
                         b_first_col=lr1)

      ! work[1:nvira, 1] => X(1:nocca, 1:nvirb) [nocca-1, 1+shift:nvirb]
      CALL cp_fm_set_all(this%work_vir_1, 0.0_dp)
      DO i = 1, nvira
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), lr1, i+shift, scalar)
         CALL cp_fm_set_element(this%work_vir_1, i, 1, scalar)
      END DO

      ! tmp (1:nvira, 1) = vb(:,nocca+1:nbf)*bvec(lr1, nocca+1:nbf)
      CALL parallel_gemm('N', 'N', nbf, 1, nvira, &
                         1.0_dp, gs_mos(2)%mos_full, &
                                 this%work_vir_1, &
                         0.0_dp, this%work2_ao_1, &
                         a_first_col=nocca+1)
      ! 2: bo1v += va(:,lr1) * tmp^T
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                         1.0_dp, gs_mos(1)%mos_full, &
                                 this%work2_ao_1, &
                         1.0_dp, this%evects_ao_ao(2), &
                         a_first_col=lr1)
      ! 5: o21v -= tmp * va(:,lr2)^T
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                        -1.0_dp, this%work2_ao_1, &
                                 gs_mos(1)%mos_full, &
                         1.0_dp, this%evects_ao_ao(5), &
                         b_first_col=lr2)

      ! work[1:noccb, 1] => X(1:nocca, 1:nvirb) [1:noccb, lr1]
      CALL cp_fm_set_all(this%work_noc_1, 0.0_dp)
      DO i = 1, noccb
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), i, lr1, scalar)
         CALL cp_fm_set_element(this%work_noc_1, i, 1, scalar)
      END DO

      ! tmp (1:nvira, 1) = va(:,1:noccb)*bvec(1:noccb, lr1)
      CALL parallel_gemm('N', 'N', nbf, 1, noccb, &
                         1.0_dp, gs_mos(1)%mos_full, &
                                 this%work_noc_1, &
                         0.0_dp, this%work2_ao_1)
      ! 3: bco1 += tmp * vb(:,lr1)
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                         1.0_dp, this%work2_ao_1, &
                                 gs_mos(2)%mos_full, &
                         1.0_dp, this%evects_ao_ao(3), &
                         b_first_col=lr1)
      ! 6: co12 += vb(:,lr2) * tmp1^T
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                         1.0_dp, gs_mos(2)%mos_full, &
                                 this%work2_ao_1, &
                         1.0_dp, this%evects_ao_ao(6), &
                         a_first_col=lr2)

      ! work[1:noccb, 1] => X(1:nocca, 1:nvirb) [1:noccb, lr1]
      CALL cp_fm_set_all(this%work_noc_1, 0.0_dp)
      DO i = 1, noccb
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), i, lr2, scalar)
         CALL cp_fm_set_element(this%work_noc_1, i, 1, scalar)
      END DO

      ! tmp (1:nvira, 1) = va(:,1:noccb)*bvec(1:noccb, lr2)
      CALL parallel_gemm('N', 'N', nbf, 1, noccb, &
                         1.0_dp, gs_mos(1)%mos_full, &
                                 this%work_noc_1, &
                         0.0_dp, this%work2_ao_1)
      ! 4: bco2 += tmp * vb(:,lr2)
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                         1.0_dp, this%work2_ao_1, &
                                 gs_mos(2)%mos_full, &
                         1.0_dp, this%evects_ao_ao(4), &
                         b_first_col=lr2)
      ! 6: co12 -= vb(:,lr1) * tmp1^T
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                        -1.0_dp, gs_mos(2)%mos_full, &
                                 this%work2_ao_1, &
                         1.0_dp, this%evects_ao_ao(6), &
                         a_first_col=lr1)

      ! Sum
      DO i = 1, 4
         CALL cp_fm_scale_and_add( &
               alpha=1.0_dp, matrix_a=this%evects_ao_ao(7), &
               beta=1.0_dp,  matrix_b=this%evects_ao_ao(i))
      END DO

      DO i = 1, nvira
         DO j = 1, noccb
            CALL cp_fm_get_element(tddft_env%evects_mo(ivec), j, i+shift, scalar)
            CALL cp_fm_set_element(this%work_vir_noc, i, j, scalar)
         END DO
      END DO

      ! tmp[nbf,noccb] = A[1:nbf,nocca+1:nbf]*B[1:nvira,noccb]
      CALL parallel_gemm('N', 'N', nbf, noccb, nvira, &
                         1.0_dp, gs_mos(2)%mos_full, &
                                 this%work_vir_noc, &
                         0.0_dp, this%work_ao_nocb, &
                         a_first_col=nocca+1)
      ! ball[nbf, nbf] = A[1:nbf,1:noccb]*B[1:nbf,noccb]
      CALL parallel_gemm('N', 'T', nbf, nbf, noccb, &
                         1.0_dp, gs_mos(1)%mos_full, &
                                 this%work_ao_nocb, &
                         1.0_dp, this%evects_ao_ao(7))


      SELECT CASE (tddft_env%response_type)
      CASE(2) ! MRSF-S
         ! ball += va(:,lr2) * [vb(:,lr1) * bvec(lr2,lr1)]^T
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), lr2, lr1, scalar)
         CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                            scalar, gs_mos(1)%mos_full, &
                                    gs_mos(2)%mos_full, &
                            1.0_dp, this%evects_ao_ao(7), &
                            a_first_col=lr2, &
                            b_first_col=lr1)
         ! ball += va(:,lr1) * [vb(:,lr2) * bvec(lr1,lr2)]^T
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), lr1, lr2, scalar)
         CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                            scalar, gs_mos(1)%mos_full, &
                                    gs_mos(2)%mos_full, &
                            1.0_dp, this%evects_ao_ao(7), &
                            a_first_col=lr1, &
                            b_first_col=lr2)
         ! ball += va(:,lr1) * [vb(:,lr1) * bvec(lr1,lr1) * isqrt2]^T
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), lr1, lr1, scalar)
         scalar = scalar/sqrt(2.0_dp)
         CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                            scalar, gs_mos(1)%mos_full, &
                                    gs_mos(2)%mos_full, &
                            1.0_dp, this%evects_ao_ao(7), &
                            a_first_col=lr1, &
                            b_first_col=lr1)
         ! ball += va(:,lr2) * [-vb(:,lr2) * bvec(lr1,lr1) * isqrt2]^T
         scalar = -scalar
         CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                            scalar, gs_mos(1)%mos_full, &
                                    gs_mos(2)%mos_full, &
                            1.0_dp, this%evects_ao_ao(7), &
                            a_first_col=lr2, &
                            b_first_col=lr2)
      CASE(3) ! MRSF-T
         ! ball += va(:,lr1) * [vb(:,lr1) * bvec(lr1,lr1) * isqrt2]^T
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), lr1, lr1, scalar)
         scalar = scalar/sqrt(2.0_dp)
         CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                            scalar, gs_mos(1)%mos_full, &
                                    gs_mos(2)%mos_full, &
                            1.0_dp, this%evects_ao_ao(7), &
                            a_first_col=lr1, &
                            b_first_col=lr1)
         ! ball += va(:,lr2) * [-vb(:,lr2) * bvec(lr1,lr1) * isqrt2]^T
         CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                            scalar, gs_mos(1)%mos_full, &
                                    gs_mos(2)%mos_full, &
                            1.0_dp, this%evects_ao_ao(7), &
                            a_first_col=lr2, &
                            b_first_col=lr2)
      END SELECT

      CALL timestop(handle)
   END SUBROUTINE mrsf_mo_to_ao_transform

   SUBROUTINE mrsf_apply_operator(this, qs_env, tddft_env, kernel_env, sub_env, &
                                tddfpt_control, ivec, logger)
      CLASS(tddft_mrsf_theory_type), INTENT(INOUT)      :: this
      TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(kernel_env_type)                           :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                  :: sub_env
      TYPE(tddfpt2_control_type), POINTER             :: tddfpt_control
      INTEGER, INTENT(IN)                             :: ivec
      TYPE(cp_logger_type), POINTER                   :: logger
      INTEGER :: imat
      IF (.FALSE. .AND. ivec == 1) THEN
         print *, this%unused
         tddft_env = tddft_env
         NULLIFY(qs_env)
         kernel_env = kernel_env
         sub_env = sub_env
         NULLIFY(tddfpt_control)
      END IF
      !-----------------------------------------------
      ! ERI contraction: Sigma = ERI * X [DBCSR type]
      !-----------------------------------------------
      IF (tddfpt_control%do_hfx) THEN
         DO imat = 7, 7
            !--------------------------------------
            ! Get Symmetric part of X vector
            !--------------------------------------
            CALL cp_fm_set_all(matrix=this%work1_ao_ao, alpha=0.0_dp)
            ! work1 = 0.5*X
            CALL cp_fm_scale_and_add( &
                      alpha=0.0_dp, matrix_a=this%work1_ao_ao, &
                      beta=0.5_dp,  matrix_b=this%evects_ao_ao(imat))
            ! work2 = X^T
            CALL cp_fm_transpose(this%evects_ao_ao(imat), this%work2_ao_ao)

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
            !--------------------------------------
            ! Get Antisymmetric part of X vector
            !--------------------------------------
            CALL cp_fm_set_all(matrix=this%work1_ao_ao, alpha=0.0_dp)

            ! work1 = 0.5*X
            CALL cp_fm_scale_and_add( &
                      alpha=0.0_dp, matrix_a=this%work1_ao_ao, &
                      beta=0.5_dp,  matrix_b=this%evects_ao_ao(imat))
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

            ! Symmetric
            CALL dbcsr_set(matrix=this%work_hmat_symm(1)%matrix, alpha=0.0_dp)
            CALL tddft_hfx_matrix(matrix_ks=this%work_hmat_symm, &
                                  rho_ao=this%work_rho_ao_ao_symm, &
                                  qs_env=qs_env)

            ! Antisymmetric
            CALL dbcsr_set(matrix=this%work_hmat_asymm(1)%matrix, alpha=0.0_dp)
            CALL tddft_hfx_matrix(matrix_ks=this%work_hmat_asymm, &
                                  rho_ao=this%work_rho_ao_ao_asymm, &
                                  qs_env=qs_env)
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

            CALL cp_fm_set_all(matrix=this%hmat_ao_ao(imat), alpha=0.0_dp)
            ! A = 0.0*A + 1.0*B
            CALL cp_fm_scale_and_add( &
                     alpha=0.0_dp, matrix_a=this%hmat_ao_ao(imat), &
                     beta=1.0_dp,  matrix_b=this%work1_ao_ao)
            ! A = 1.0*A + 1.0*B
            CALL cp_fm_scale_and_add( &
                     alpha=1.0_dp, matrix_a=this%hmat_ao_ao(imat), &
                     beta=1.0_dp,  matrix_b=this%work2_ao_ao)
            CALL cp_fm_set_all(matrix=this%work2_ao_ao, alpha=0.0_dp)
            ! Print
            CALL print_mat_ao_ao(this%hmat_ao_ao(imat), "hmat_ao_HF_tot", logger, ivec)
         END DO
      END IF

   END SUBROUTINE mrsf_apply_operator

   SUBROUTINE mrsf_ao_to_mo_transform(this, tddft_env, gs_mos, ivec, logger)
      CLASS(tddft_mrsf_theory_type), TARGET, INTENT(INOUT)    :: this
      TYPE(tddft_info_type),    INTENT(INOUT)         :: tddft_env
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), POINTER                        :: gs_mos
      INTEGER, INTENT(IN)                             :: ivec
      TYPE(cp_logger_type), POINTER                   :: logger
      TYPE(cp_fm_type), POINTER :: agdlr, ado2v, ado1v, adco1, adco2, ao21v, aco12
      INTEGER :: i, j, lr1, lr2, nbf, nocca, noccb, nvira, nvirb, shift
      REAL(KIND=dp) :: scalar, val

      IF (.FALSE. .AND. ivec == 1) THEN
         print *, this%unused
         tddft_env = tddft_env
         gs_mos = gs_mos
         logger = logger
      END IF

      i = i
      j = j

      ! Extract dimensions from tddft_env
      nbf = tddft_env%nbf
      nocca = tddft_env%nocca
      noccb = tddft_env%noccb
      nvira = tddft_env%nvira  ! nvira = nbf - nocca
      nvirb = tddft_env%nvirb  ! nvirb = nbf - noccb

      ! Define indices
      lr1 = nocca - 1
      lr2 = nocca
      shift = nocca - noccb

      ! Assign pointers from tddft_env%evects_ao_ao [nbf, nbf]
      agdlr => this%hmat_ao_ao(7)
      ado2v => this%hmat_ao_ao(1)
      ado1v => this%hmat_ao_ao(2)
      adco1 => this%hmat_ao_ao(3)
      adco2 => this%hmat_ao_ao(4)
      ao21v => this%hmat_ao_ao(5)
      aco12 => this%hmat_ao_ao(6)

      !-----------------------------------------------
      ! Transform back to MO basis
      !-----------------------------------------------
      ! C[nocca,nbf] = 1.0*(A[nbf,nocca])^T * B[nbf,nbf] + 0.0*C[nocca,nbf]
      CALL parallel_gemm('T', 'N', nocca, nbf, nbf, &
                         1.0_dp, gs_mos(1)%mos_occ, &
                                 this%hmat_ao_ao(7), &
                         0.0_dp, this%work_noc_ao)

      ! C[nocca,nvirb] = 1.0*A[noc,nbf] * B[nbf,nvirb] + 0.0*C[nocca,nvirb]
      CALL parallel_gemm('N', 'N', nocca, nvirb, nbf, &
                         1.0_dp, this%work_noc_ao, &
                                 gs_mos(2)%mos_virt, &
                         0.0_dp, tddft_env%sigma_vecs(ivec))
      ! Print
      CALL print_mat_ao_ao(tddft_env%sigma_vecs(ivec), "sigma_vecs", logger, ivec)

      CALL cp_fm_get_element(tddft_env%sigma_vecs(ivec), lr1, shift-1, scalar)
      CALL cp_fm_get_element(tddft_env%sigma_vecs(ivec), lr2, shift, val)
      SELECT CASE (tddft_env%response_type)
      CASE(2) ! MRSF-S
         scalar = (scalar - val)/sqrt(2.0_dp)
      CASE(3) ! MRSF-T
         scalar = (scalar + val)/sqrt(2.0_dp)
      END SELECT
      CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr1, shift-1, scalar)

      CALL parallel_gemm('N', 'N', nbf, 1, nbf, &
                         1.0_dp, this%hmat_ao_ao(2), &  ! ado1v
                                 gs_mos(2)%mos_full, &
                         0.0_dp, this%work1_ao_1, &     ! tmp
                         b_first_col=lr2)
      CALL parallel_gemm('N', 'N', nbf, 1, nbf, &
                         1.0_dp, this%hmat_ao_ao(6), &  ! aco12
                                 gs_mos(2)%mos_full, &
                         1.0_dp, this%work1_ao_1, &     ! tmp
                         b_first_col=lr1)
      CALL parallel_gemm('T', 'N', noccb, 1, nbf, &
                         1.0_dp, gs_mos(2)%mos_occ, &
                                 this%work1_ao_1, &     ! tmp
                         0.0_dp, this%work_noc_1)      ! tmp2[1:noccb, 1]
      DO i = 1, noccb
         CALL cp_fm_get_element(this%work_noc_1, i, 1, scalar)
         CALL cp_fm_get_element(tddft_env%sigma_vecs(ivec), i, lr2-noccb, val)
         val = val + scalar
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), i, lr2-noccb, val)
      END DO

      CALL parallel_gemm('N', 'N', nbf, 1, nbf, &
                         1.0_dp, this%hmat_ao_ao(1), &  ! ado2v
                                 gs_mos(2)%mos_full, &
                         0.0_dp, this%work1_ao_1, &     ! tmp
                         b_first_col=lr1)
      CALL parallel_gemm('N', 'N', nbf, 1, nbf, &
                        -1.0_dp, this%hmat_ao_ao(6), &  ! aco12
                                 gs_mos(2)%mos_full, &
                         1.0_dp, this%work1_ao_1, &     ! tmp
                         b_first_col=lr2)
      CALL parallel_gemm('T', 'N', noccb, 1, nbf, &
                         1.0_dp, gs_mos(2)%mos_occ, &
                                 this%work1_ao_1, &     ! tmp
                         0.0_dp, this%work_noc_1)      ! tmp2[1:noccb, 1]
      DO i = 1, noccb
         CALL cp_fm_get_element(this%work_noc_1, i, 1, scalar)
         CALL cp_fm_get_element(tddft_env%sigma_vecs(ivec), i, lr1-noccb, val)
         val = val + scalar
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), i, lr1-noccb, val)
      END DO

      CALL parallel_gemm('T', 'N', nbf, 1, nbf, &
                         1.0_dp, this%hmat_ao_ao(4), &  ! adco2^T
                                 gs_mos(1)%mos_full, &  ! va
                         0.0_dp, this%work1_ao_1, &     ! tmp
                         b_first_col=lr1)

      CALL parallel_gemm('T', 'N', nbf, 1, nbf, &
                         1.0_dp, this%hmat_ao_ao(5), &  ! ao21v^T
                                 gs_mos(1)%mos_full, &  ! va
                         1.0_dp, this%work1_ao_1, &     ! tmp
                         b_first_col=lr2)

      CALL parallel_gemm('T', 'N', nvira, 1, nbf, &
                         1.0_dp, gs_mos(1)%mos_virt, &  ! va_virt [nbf, nvira]
                                 this%work1_ao_1, &     ! tmp [nbf, 1]
                         0.0_dp, this%work_vir_1)       ! tmp2 [nvira, 1]

      DO i = 1, nvira
         CALL cp_fm_get_element(this%work_vir_1, i, 1, scalar)
         CALL cp_fm_get_element(tddft_env%sigma_vecs(ivec), lr1, i, val)
         val = val + scalar
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr1, i, val)
      END DO

      CALL parallel_gemm('T', 'N', nbf, 1, nbf, &
                         1.0_dp, this%hmat_ao_ao(3), &  ! adco2^T
                                 gs_mos(1)%mos_full, &  ! va
                         0.0_dp, this%work1_ao_1, &     ! tmp
                         b_first_col=lr2)

      CALL parallel_gemm('T', 'N', nbf, 1, nbf, &
                        -1.0_dp, this%hmat_ao_ao(5), &  ! ao21v^T
                                 gs_mos(1)%mos_full, &  ! va
                         1.0_dp, this%work1_ao_1, &     ! tmp
                         b_first_col=lr1)

      CALL parallel_gemm('T', 'N', nvira, 1, nbf, &
                         1.0_dp, gs_mos(1)%mos_virt, &  ! va_virt [nbf, nvira]
                                 this%work1_ao_1, &     ! tmp [nbf, 1]
                         0.0_dp, this%work_vir_1)       ! tmp2 [nvira, 1]

      DO i = 1, nvira
         CALL cp_fm_get_element(this%work_vir_1, i, 1, scalar)
         CALL cp_fm_get_element(tddft_env%sigma_vecs(ivec), lr2, i, val)
         val = val + scalar
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr2, i, val)
      END DO

      SELECT CASE (tddft_env%response_type)
      CASE(2) ! MRSF-S
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr2, shift, 0.0_dp)
      CASE(3) ! MRSF-T
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr1, shift, 0.0_dp)
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr2, shift-1, 0.0_dp)
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr2, shift, 0.0_dp)
      END SELECT

   END SUBROUTINE mrsf_ao_to_mo_transform

   SUBROUTINE mrsf_print_theory_data(this, qs_env, tddft_env, gs_mos, tddfpt_control, log_unit)
      CLASS(tddft_mrsf_theory_type), INTENT(INOUT) :: this
      TYPE(qs_environment_type), POINTER :: qs_env
      TYPE(tddft_info_type), INTENT(INOUT) :: tddft_env
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), POINTER :: gs_mos
      TYPE(tddfpt2_control_type), POINTER  :: tddfpt_control
      INTEGER :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mrsf_print_theory_data'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      this%unused = this%unused
      tddft_env = tddft_env
      qs_env = qs_env
      gs_mos = gs_mos
      tddfpt_control = tddfpt_control
      log_unit = log_unit

      CALL timestop(handle)

   END SUBROUTINE mrsf_print_theory_data

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

END MODULE qs_tddft_theory_mrsf
