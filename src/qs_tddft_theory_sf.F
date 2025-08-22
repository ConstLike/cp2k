!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!
MODULE qs_tddft_theory_sf

   USE qs_energy_types,                 ONLY: qs_energy_type
   USE kinds,                           ONLY: default_path_length
   USE qs_tddfpt2_properties,           ONLY: tddfpt_dipole_operator,&
                                              tddfpt_print_exciton_descriptors,&
                                              tddfpt_print_nto_analysis
   USE pw_env_types,                    ONLY: pw_env_release,&
                                              pw_env_retain,&
                                              pw_env_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
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
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: evolt
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_theory_sf'

   PUBLIC :: tddft_sf_theory_type

   TYPE, EXTENDS(tddft_theory_type) :: tddft_sf_theory_type

      !> Hamiltonian matrix in AO basis
      TYPE(cp_fm_type), POINTER :: hmat_ao_ao => NULL()

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

      !> Per-spin AO Fock matrices (used in transformation to MO)
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: fock

      !> Fock matrix for Alpha occupied orbitals block [nocca, nocca]
      TYPE(cp_fm_type), POINTER :: Focka_noc_noc => NULL()

      !> Fock matrix for Beta virtual orbitals block [nvirb, nvirb]
      TYPE(cp_fm_type), POINTER :: Fockb_vir_vir => NULL()

      !-------------------------------
      ! KERNELs related
      !-------------------------------
      !> [XC kernel] NO NAME
      TYPE(cp_fm_type), POINTER :: rho_ia_ao_fm => NULL()
      TYPE(qs_rho_type), POINTER :: gs_rho_xc_struct => NULL()
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:) :: A_ia_munu_sub
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER :: A_ia_rspace_sub => NULL()
      TYPE(qs_rho_type), POINTER :: rho_orb_struct => NULL()
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

      !> Per-spin AO scratch matrices for Fock transformation [nbf, nbf]
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: dens_diff_ij
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: dens_diff_ab

      ! -- theory-specific fields if needed --
   CONTAINS
      PROCEDURE, PASS(this) :: allocate_theory_data     => sf_allocate_theory_data
      PROCEDURE, PASS(this) :: deallocate_theory_data   => sf_deallocate_theory_data
      PROCEDURE, PASS(this) :: initialize_guess_vectors => sf_initialize_guess_vectors
      PROCEDURE, PASS(this) :: mo_to_ao_transform       => sf_mo_to_ao_transform
      PROCEDURE, PASS(this) :: apply_operator           => sf_apply_operator
      PROCEDURE, PASS(this) :: ao_to_mo_transform       => sf_ao_to_mo_transform
      PROCEDURE, PASS(this) :: print_theory_data        => sf_print_theory_data
   END TYPE

CONTAINS

   SUBROUTINE sf_allocate_theory_data(this, qs_env, tddft_env, blacs_env)
      CLASS(tddft_sf_theory_type), INTENT(INOUT) :: this
      TYPE(qs_environment_type), POINTER :: qs_env
      TYPE(tddft_info_type), INTENT(INOUT) :: tddft_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'sf_allocate_theory_data'

      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(dbcsr_distribution_type), POINTER :: dbcsr_dist
      TYPE(pw_env_type), POINTER :: pw_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
      INTEGER                                            :: handle, ivec, ispin

      CALL timeset(routineN, handle)

      !=== hmat_ao_ao [nbf, nbf] =======================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=blacs_env)
      ALLOCATE(this%hmat_ao_ao)
      CALL cp_fm_create(this%hmat_ao_ao, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%hmat_ao_ao, 0.0_dp)
      !=================================================
      !=== evects_ao_ao (1) [nbf, nbf] =============
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=blacs_env)
      ALLOCATE(this%evects_ao_ao(1))
      DO ivec = 1, 1
         CALL cp_fm_create(this%evects_ao_ao(1), fm_struct)
         CALL cp_fm_set_all(this%evects_ao_ao(1), 0.0_dp)
      END DO
      CALL cp_fm_struct_release(fm_struct)
      !=================================================
      !=== focka_noc_noc [nocca, nocca] ================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nocca, &
                               ncol_global=tddft_env%nocca, &
                               context=blacs_env)
      ALLOCATE(this%focka_noc_noc)
      CALL cp_fm_create(this%focka_noc_noc, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%focka_noc_noc, 0.0_dp)
      !=================================================
      !=== fockb_vir_vir [nvirb, nvirb] ================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nvirb, &
                               ncol_global=tddft_env%nvirb, &
                               context=blacs_env)
      ALLOCATE(this%fockb_vir_vir)
      CALL cp_fm_create(this%fockb_vir_vir, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%fockb_vir_vir, 0.0_dp)
      !=================================================
      !=== fock (nSpins) [nbf, nvirb] ==================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=blacs_env)
      ALLOCATE(this%fock(tddft_env%nSpins))
      DO ispin = 1, tddft_env%nSpins
         CALL cp_fm_create(this%fock(ispin), fm_struct)
         CALL cp_fm_set_all(this%fock(ispin), 0.0_dp)
      END DO
      CALL cp_fm_struct_release(fm_struct)
      !=================================================
      !=== work_ao_vir [nbf, nvirb] ====================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nvirb, &
                               context=blacs_env)
      ALLOCATE(this%work_ao_vir)
      CALL cp_fm_create(this%work_ao_vir, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_ao_vir, 0.0_dp)
      !=================================================
      !=== work_noc_ao [nocca, nbf] ====================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nocca, &
                               ncol_global=tddft_env%nbf, &
                               context=blacs_env)
      ALLOCATE(this%work_noc_ao)
      CALL cp_fm_create(this%work_noc_ao, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_noc_ao, 0.0_dp)
      !=================================================
      !=== work_noc_noc [nocca, nocca] =================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nocca, &
                               ncol_global=tddft_env%nocca, &
                               context=blacs_env)
      ALLOCATE(this%work_noc_noc)
      CALL cp_fm_create(this%work_noc_noc, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_noc_noc, 0.0_dp)
      !=================================================
      !=== work_ao_noc [nbf, nocca] ====================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nocca, &
                               context=blacs_env)
      ALLOCATE(this%work_ao_noc)
      CALL cp_fm_create(this%work_ao_noc, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_ao_noc, 0.0_dp)
      !=================================================
      !=== work_vir_vir [nvirb, nvirb] =================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nvirb, &
                               ncol_global=tddft_env%nvirb, &
                               context=blacs_env)
      ALLOCATE(this%work_vir_vir)
      CALL cp_fm_create(this%work_vir_vir, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_vir_vir, 0.0_dp)
      !=================================================
      !=== work_noc_vir [nocca, nvirb] =================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nocca, &
                               ncol_global=tddft_env%nvirb, &
                               context=blacs_env)
      ALLOCATE(this%work_noc_vir)
      CALL cp_fm_create(this%work_noc_vir, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work_noc_vir, 0.0_dp)
      !=================================================
      !=== work1_ao_ao [nbf, nbf] ======================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=blacs_env)
      ALLOCATE(this%work1_ao_ao)
      CALL cp_fm_create(this%work1_ao_ao, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work1_ao_ao, 0.0_dp)
      !=================================================
      !=== work2_ao_ao [nbf, nbf] ======================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=blacs_env)
      ALLOCATE(this%work2_ao_ao)
      CALL cp_fm_create(this%work2_ao_ao, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%work2_ao_ao, 0.0_dp)
      !=================================================
      !=== rho_ia_ao_fm [nbf, nbf] =====================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=blacs_env)
      ALLOCATE(this%rho_ia_ao_fm)
      CALL cp_fm_create(this%rho_ia_ao_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(this%rho_ia_ao_fm, 0.0_dp)
      !=== rho_fm_ao [nbf, nbf] ========================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=blacs_env)
      ALLOCATE(this%rho_fm_ao(tddft_env%nSpins))
      DO ispin = 1, tddft_env%nSpins
         CALL cp_fm_create(this%rho_fm_ao(ispin), fm_struct)
         CALL cp_fm_set_all(this%rho_fm_ao(ispin), 0.0_dp)
      END DO
      CALL cp_fm_struct_release(fm_struct)
      !=================================================
      !=== dens_diff_ij [nbf, nbf] =====================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=blacs_env)
      ALLOCATE(this%dens_diff_ij(tddft_env%nStates))
      DO ivec = 1, tddft_env%nStates
         CALL cp_fm_create(this%dens_diff_ij(ivec), fm_struct)
         CALL cp_fm_set_all(this%dens_diff_ij(ivec), 0.0_dp)
      END DO
      CALL cp_fm_struct_release(fm_struct)
      !=================================================
      !=== dens_diff_ab [nbf, nbf] =====================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nbf, &
                               context=blacs_env)
      ALLOCATE(this%dens_diff_ab(tddft_env%nStates))
      DO ivec = 1, tddft_env%nStates
         CALL cp_fm_create(this%dens_diff_ab(ivec), fm_struct)
         CALL cp_fm_set_all(this%dens_diff_ab(ivec), 0.0_dp)
      END DO
      CALL cp_fm_struct_release(fm_struct)
      !=================================================
      !=== dipoles_mo(3) [nocca, nvirb] ================
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nocca, &
                               ncol_global=tddft_env%nvirb, &
                               context=blacs_env)
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
                      dbcsr_dist=dbcsr_dist, &
                      pw_env=pw_env, &
                      sab_orb=sab_orb)

      NULLIFY (this%A_ia_munu_sub)
      CALL dbcsr_allocate_matrix_set(this%A_ia_munu_sub, 1)
      DO ispin = 1, 1
         CALL dbcsr_init_p(this%A_ia_munu_sub(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(this%A_ia_munu_sub(ispin)%matrix, &
                                          template=matrix_s(1)%matrix, &
                                          dbcsr_dist=dbcsr_dist, &
                                          sab=sab_orb)
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
                  dbcsr_dist=dbcsr_dist, &
                  sab=sab_orb)

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
                                          dbcsr_dist=dbcsr_dist, &
                                          sab=sab_orb)

         CALL dbcsr_init_p(this%work_hmat_asymm(ispin)%matrix)
         CALL dbcsr_create(this%work_hmat_asymm(ispin)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric, &
                           template=this%work_hmat_symm(ispin)%matrix)
         CALL dbcsr_complete_redistribute(this%work_hmat_symm(ispin)%matrix, &
                                          this%work_hmat_asymm(ispin)%matrix)
      END DO

      NULLIFY (this%rho_orb_struct)
      ALLOCATE (this%rho_orb_struct)
      CALL qs_rho_create(this%rho_orb_struct)
      CALL qs_rho_set(this%rho_orb_struct, rho_ao=this%work_rho_ao_ao_symm)
      CALL qs_rho_rebuild(this%rho_orb_struct, qs_env, &
                          rebuild_ao=.FALSE., &
                          rebuild_grids=.TRUE.)

      !--------------------
      ! ALLOCATE  PW
      !--------------------
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

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
!                            rho_struct_sub=this%rho_orb_struct, &
!                            xc_section=xc_section, &
!                            is_rks_triplets=tddfpt_control%rks_triplets, &
!                            sub_env=sub_env)
!     kernel_env%full_kernel => full_kernel_env
!     NULLIFY (kernel_env%admm_kernel)

      CALL timestop(handle)

   END SUBROUTINE sf_allocate_theory_data

   SUBROUTINE sf_deallocate_theory_data(this, qs_env)
      CLASS(tddft_sf_theory_type), INTENT(INOUT) :: this
      TYPE(qs_environment_type), POINTER :: qs_env

      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
      TYPE(pw_env_type), POINTER :: pw_env
      INTEGER :: ivec, ispin

      CALL get_qs_env(qs_env=qs_env, &
                      pw_env=pw_env, &
                      blacs_env=blacs_env)

      IF (ASSOCIATED(this%hmat_ao_ao)) THEN
         CALL cp_fm_release(this%hmat_ao_ao)
         DEALLOCATE(this%hmat_ao_ao)
      END IF
      IF (ALLOCATED(this%evects_ao_ao)) THEN
         DO ivec = 1, SIZE(this%evects_ao_ao)
            CALL cp_fm_release(this%evects_ao_ao(ivec))
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
      IF (ASSOCIATED(this%rho_ia_ao_fm)) THEN
         CALL cp_fm_release(this%rho_ia_ao_fm)
         DEALLOCATE(this%rho_ia_ao_fm)
      END IF
      IF (ALLOCATED(this%dens_diff_ij)) THEN
         DO ivec = 1, SIZE(this%dens_diff_ij)
            CALL cp_fm_release(this%dens_diff_ij(ivec))
         END DO
         DEALLOCATE(this%dens_diff_ij)
      END IF
      IF (ALLOCATED(this%dens_diff_ab)) THEN
         DO ivec = 1, SIZE(this%dens_diff_ab)
            CALL cp_fm_release(this%dens_diff_ab(ivec))
         END DO
         DEALLOCATE(this%dens_diff_ab)
      END IF
      IF (ASSOCIATED(this%rho_fm_ao)) THEN
         DO ispin = 1, SIZE(this%rho_fm_ao)
            CALL cp_fm_release(this%rho_fm_ao(ispin))
         END DO
         DEALLOCATE(this%rho_fm_ao)
      END IF
      IF (ALLOCATED(this%fock)) THEN
         DO ispin = 1, SIZE(this%fock)
            CALL cp_fm_release(this%fock(ispin))
         END DO
         DEALLOCATE(this%fock)
      END IF
      IF (ASSOCIATED(this%focka_noc_noc)) THEN
         CALL cp_fm_release(this%focka_noc_noc)
         DEALLOCATE(this%Focka_noc_noc)
      END IF
      IF (ASSOCIATED(this%Fockb_vir_vir)) THEN
         CALL cp_fm_release(this%Fockb_vir_vir)
         DEALLOCATE(this%Fockb_vir_vir)
      END IF
      IF (ALLOCATED(this%dipoles_mo)) THEN
         DO ivec = 1, SIZE(this%dipoles_mo)
            CALL cp_fm_release(this%dipoles_mo(ivec))
         END DO
         DEALLOCATE(this%dipoles_mo)
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

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

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

      IF (ASSOCIATED(this%rho_orb_struct)) THEN
!        CALL qs_rho_release(this%rho_orb_struct)
!        DEALLOCATE(this%rho_orb_struct)
         NULLIFY(this%rho_orb_struct)
      END IF

   END SUBROUTINE sf_deallocate_theory_data

   SUBROUTINE sf_initialize_guess_vectors(this, qs_env, tddft_env, gs_mos, &
                                          tddfpt_control, log_unit)
      CLASS(tddft_sf_theory_type), INTENT(INOUT)         :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), POINTER                           :: gs_mos
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      INTEGER                                            :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sf_initialize_guess_vectors'

      INTEGER                                            :: handle, imo_occ, imo_virt, &
                                                            istate, ispin, &
                                                            nStates_selected, istate_orig
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_m_occ
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks

      CALL timeset(routineN, handle)

      IF (.FALSE.) THEN
         print *, this%unused
         qs_env = qs_env
         NULLIFY(tddfpt_control)
      END IF
      CALL get_qs_env(qs_env, &
                      matrix_ks=matrix_ks)

      ! FOCK
      DO ispin = 1, tddft_env%nSpins
         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, this%fock(ispin))
         ! work1_ao_ao = gs_mos(1)%mos_full^T * this%fock
         CALL parallel_gemm('T', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, gs_mos(ispin)%mos_full, &
                                    this%fock(ispin), &
                            0.0_dp, this%work1_ao_ao)

         ! this%fock = work_ao_ao * gs_mos(2)%mos_full
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, this%work1_ao_ao, &
                                    gs_mos(ispin)%mos_full, &
                            0.0_dp, this%fock(ispin))
      END DO

      ! --- focka_occ = F_alpha( occ , occ )  ------------------------------------
      CALL cp_fm_to_fm_submat(this%fock(1), this%focka_noc_noc, &
                              nrow=tddft_env%nocca, ncol=tddft_env%nocca, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=1, t_firstcol=1)

      ! --- fockb_vir = F_beta( vir , vir )  ------------------------------------
      CALL cp_fm_to_fm_submat(this%fock(2), this%fockb_vir_vir, &
                              nrow=tddft_env%nvirb, &
                              ncol=tddft_env%nvirb, &
                              s_firstrow=tddft_env%noccb+1, &
                              s_firstcol=tddft_env%noccb+1, &
                              t_firstrow=1, t_firstcol=1)

      nStates_selected = gs_mos(1)%homo*(SIZE(gs_mos(2)%evals_full)-gs_mos(2)%homo)

      ALLOCATE (inds(nStates_selected))
      ALLOCATE (e_virt_m_occ(nStates_selected))

      istate = 0
      DO imo_virt = gs_mos(2)%homo+1, SIZE(gs_mos(2)%evals_full)
         DO imo_occ = 1, gs_mos(1)%homo
            istate = (imo_virt-gs_mos(2)%homo-1)*gs_mos(1)%homo+imo_occ
            tddft_env%e_virt_minus_occ(istate) = gs_mos(2)%evals_full(imo_virt)-gs_mos(1)%evals_full(imo_occ)
         END DO
      END DO

      e_virt_m_occ(1:nStates_selected) = tddft_env%e_virt_minus_occ(1:nStates_selected)

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

      DO istate = 1, tddft_env%nVec
         tddft_env%evals(istate) = e_virt_m_occ(istate)
         istate_orig = inds(istate)
         imo_occ = MOD(istate_orig - 1, gs_mos(1)%homo) + 1
         imo_virt = (istate_orig - 1) / gs_mos(1)%homo + 1
         CALL cp_fm_set_element(tddft_env%evects_mo(istate), imo_occ, imo_virt, 1.0_dp)
         IF (istate > tddft_env%nStates) CYCLE
         IF (log_unit > 0) &
            WRITE (log_unit, '(T7,I8,T24,I8,T37,A5,T45,I8,T54,A5,T60,F14.5)') &
            istate, imo_occ, '(alp)', gs_mos(2)%homo + imo_virt, '(bet)', e_virt_m_occ(istate)*evolt
      end do

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T7,A,T50,I24)') 'Number of active states:', &
            SIZE(gs_mos(1)%evals_occ)*SIZE(gs_mos(2)%evals_virt)
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      DEALLOCATE (inds, e_virt_m_occ)

      CALL timestop(handle)

   END SUBROUTINE sf_initialize_guess_vectors

   SUBROUTINE sf_mo_to_ao_transform(this, tddft_env, gs_mos, ivec, logger)
      CLASS(tddft_sf_theory_type), TARGET, INTENT(INOUT)      :: this
      TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), POINTER                        :: gs_mos
      INTEGER, INTENT(IN)                             :: ivec
      TYPE(cp_logger_type), POINTER                   :: logger
      IF (.FALSE. .and. ivec == 1) THEN
         print *, this%unused
         tddft_env = tddft_env
         NULLIFY(gs_mos)
      END IF

      !--------------------------------------
      ! Orthogonal transform: X vector to AO
      !--------------------------------------
      ! X_ao = C_alpha * X_mo * C_beta^T
      ! C[nbf,nvirb] = 1.0*A[nbf,nocca]*B[nocca,nvirb] + 0.0*C[nbf,nvirb]
      CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nvirb, tddft_env%nocca, &
                         1.0_dp, gs_mos(1)%mos_occ, &
                                 tddft_env%evects_mo(ivec), &
                         0.0_dp, this%work_ao_vir)
      ! C[nbf,nbf] = 1.0*A[nbf,nvirb]*(B[nbf,nvirb])^T + 0.0*C[nbf,nbf]
      CALL parallel_gemm('N', 'T', tddft_env%nbf, tddft_env%nbf, tddft_env%nvirb, &
                         1.0_dp, this%work_ao_vir, &
                                 gs_mos(2)%mos_virt, &
                         0.0_dp, this%evects_ao_ao(1))

      ! Print
      CALL print_mat_ao_ao(this%evects_ao_ao(1), 'evects', logger, ivec)

   END SUBROUTINE sf_mo_to_ao_transform

   SUBROUTINE sf_apply_operator(this, qs_env, tddft_env, blacs_env, &
                                tddfpt_control, ivec, logger)
      CLASS(tddft_sf_theory_type), INTENT(INOUT)      :: this
      TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(tddfpt2_control_type), POINTER             :: tddfpt_control
      INTEGER, INTENT(IN)                             :: ivec
      TYPE(cp_logger_type), POINTER                   :: logger
      IF (.FALSE. .AND. ivec == 1) THEN
         print *, this%unused
         tddft_env = tddft_env
         NULLIFY(qs_env)
         blacs_env = blacs_env
         NULLIFY(tddfpt_control)
      END IF
      !-----------------------------------------------
      ! ERI contraction: Sigma = ERI * X [DBCSR type]
      !-----------------------------------------------
      IF (tddfpt_control%do_hfx) THEN
         !--------------------------------------
         ! Get Symmetric part of X vector
         !--------------------------------------
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
         !--------------------------------------
         ! Get Antisymmetric part of X vector
         !--------------------------------------
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

         CALL cp_fm_set_all(matrix=this%hmat_ao_ao, alpha=0.0_dp)
         ! A = 0.0*A + 1.0*B
         CALL cp_fm_scale_and_add( &
                  alpha=0.0_dp, matrix_a=this%hmat_ao_ao, &
                  beta=1.0_dp,  matrix_b=this%work1_ao_ao)
         ! A = 1.0*A + 1.0*B
         CALL cp_fm_scale_and_add( &
                  alpha=1.0_dp, matrix_a=this%hmat_ao_ao, &
                  beta=1.0_dp,  matrix_b=this%work2_ao_ao)
         CALL cp_fm_set_all(matrix=this%work2_ao_ao, alpha=0.0_dp)
         ! Print
         CALL print_mat_ao_ao(this%hmat_ao_ao, "hmat_ao_HF_tot", logger, ivec)
      END IF

!     !-----------------------------------------------
!     ! XC Contribution: .. = ..
!     !-----------------------------------------------
!     ! Density set at initial step to work_rho_ao_ao_symm
!     ! Rebuild response density grids with current trial vector
!     CALL qs_rho_update_rho(rho_struct=this%rho_orb_struct, &
!                            qs_env=qs_env)
!     DO i = 1, tddft_env%nSpins
!        CALL pw_zero(this%A_ia_rspace_sub(i))
!     END DO
!     CALL tddfpt_apply_xc(A_ia_rspace=this%A_ia_rspace_sub, &
!                          kernel_env=kernel_env%full_kernel, &
!                          rho_ia_struct=this%rho_orb_struct, &
!                          is_rks_triplets=tddfpt_control%rks_triplets, &
!                          pw_env=pw_env, &
!                          work_v_xc=this%work_v_xc, &
!                          work_v_xc_tau=this%work_v_xc_tau, &
!                          spinflip=tddfpt_control%spinflip) ! 1 collinear, 2=noncollinear
!     CALL pw_scale(this%A_ia_rspace_sub(1), &
!                   this%A_ia_rspace_sub(1)%pw_grid%dvol)
!     ! Integrate XC contribution into the TDDFT operator (add to hmat_ao_ao)
!     CALL dbcsr_set(matrix=this%A_ia_munu_sub(1)%matrix, alpha=0.0_dp)
!     CALL integrate_v_rspace(v_rspace=this%A_ia_rspace_sub(1), &
!                             hmat=this%A_ia_munu_sub(1), &
!                             qs_env=qs_env, &
!                             calculate_forces=.FALSE., &
!                             gapw=.FALSE., &
!                             pw_env_external=pw_env, &
!                             task_list_external=sub_env%task_list_orb)
!     CALL copy_dbcsr_to_fm(matrix=this%A_ia_munu_sub(1)%matrix, &
!                           fm=this%work2_ao_ao)
!     ! Print
!     CALL print_mat_ao_ao(this%work2_ao_ao, 'Awork_xc', logger, ivec)
!     CALL cp_fm_scale_and_add( &
!              alpha=1.0_dp, matrix_a=this%hmat_ao_ao, &
!              beta=1.0_dp,  matrix_b=this%work2_ao_ao)
!     CALL cp_fm_set_all(matrix=this%work2_ao_ao, alpha=0.0_dp)
!     ! Print
!     CALL print_mat_ao_ao(this%hmat_ao_ao, "hmat_ao_HF_XC_tot", logger, ivec)
   END SUBROUTINE sf_apply_operator

   SUBROUTINE sf_ao_to_mo_transform(this, tddft_env, gs_mos, ivec, logger)
      CLASS(tddft_sf_theory_type), TARGET, INTENT(INOUT)      :: this
      TYPE(tddft_info_type),    INTENT(INOUT)         :: tddft_env
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), POINTER                        :: gs_mos
      INTEGER, INTENT(IN)                             :: ivec
      TYPE(cp_logger_type), POINTER                   :: logger

      REAL(KIND=DP)                                      :: val, sig
      INTEGER                                            :: a, i, ij

      !-----------------------------------------------
      ! Transform back to MO basis
      !-----------------------------------------------
      ! C[nocca,nbf] = 1.0*(A[nbf,nocca])^T * B[nbf,nbf] + 0.0*C[nocca,nbf]
      CALL parallel_gemm('T', 'N', tddft_env%nocca, tddft_env%nbf, tddft_env%nbf, &
                         1.0_dp, gs_mos(1)%mos_occ, &
                                 this%hmat_ao_ao, &
                         0.0_dp, this%work_noc_ao)

      ! C[nocca,nvirb] = 1.0*A[noc,nbf] * B[nbf,nvirb] + 0.0*C[nocca,nvirb]
      CALL parallel_gemm('N', 'N', tddft_env%nocca, tddft_env%nvirb, tddft_env%nbf, &
                         1.0_dp, this%work_noc_ao, &
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
                            -1.0_dp, this%focka_noc_noc, &
                                     tddft_env%evects_mo(ivec), &
                             0.0_dp, this%work_noc_vir)

         ! C[nocca,nvirb] = 1.0*A[nocca,nvirb]*B[nvirb,nvirb] + 1.0*C[nocca,nvirb]
         CALL parallel_gemm('N','N', tddft_env%nocca, tddft_env%nvirb, tddft_env%nvirb, &
                             1.0_dp, tddft_env%evects_mo(ivec), &
                                     this%fockb_vir_vir, &
                             1.0_dp, this%work_noc_vir)

         ! Add corrections to sigma_mo
         ! A = 1.0 * A + 1.0 * B
         CALL cp_fm_scale_and_add( &
                  alpha=1.0_dp, matrix_a=tddft_env%sigma_vecs(ivec), &
                  beta=1.0_dp,  matrix_b=this%work_noc_vir)

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
   END SUBROUTINE sf_ao_to_mo_transform

   SUBROUTINE sf_print_theory_data(this, qs_env, tddft_env, gs_mos, ostrength, &
                                   tddfpt_control, log_unit)
      CLASS(tddft_sf_theory_type), INTENT(INOUT) :: this
      TYPE(qs_environment_type), POINTER :: qs_env
      TYPE(tddft_info_type), INTENT(INOUT) :: tddft_env
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), POINTER :: gs_mos
      REAL(kind=dp), DIMENSION(:)                     :: ostrength
      TYPE(tddfpt2_control_type), POINTER  :: tddfpt_control
      TYPE(section_vals_type), POINTER                   :: tddfpt_section
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      INTEGER :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sf_print_theory_data'

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      CHARACTER(len=2)                                   :: lsd_str
      CHARACTER(len=1)                                   :: mult_str
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: trans_dipoles
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: trans_dipoles_exc_exc
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: s2_mat
      REAL(kind=dp)                                      :: osc_strength, val, f, energ, gs_energy
      REAL(KIND=dp) :: n_spin_imbalance, s2_spin, s2_beta_occ, s2_cross, s2_anti, s2_offdiag
      INTEGER                                            :: handle, ideriv, istate, ispin
      INTEGER, ALLOCATABLE :: indices(:), i_indices(:), a_indices(:)
      REAL(KIND=dp), ALLOCATABLE :: amplitudes(:)
      INTEGER, DIMENSION(2)                       :: nmo_occ
      INTEGER :: n_amplitudes, max_amplitudes, i, a, idx, jstate, nao, nspins
      TYPE(qs_energy_type), POINTER                      :: energy
      CHARACTER(LEN=default_path_length)                 :: suffix

      CALL timeset(routineN, handle)

      SELECT CASE (tddft_env%reference_type)
      CASE (1)
         lsd_str = 'RO'
      CASE (2)
         lsd_str = ' U'
      END SELECT

      ! *** summary header ***
      CALL integer_to_string(3, mult_str)
      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,1X,A2,A,1X,A)') lsd_str, "-TDDFT states of multiplicity", TRIM(mult_str)
         SELECT CASE (tddfpt_control%dipole_form)
         CASE (tddfpt_dipole_berry)
            WRITE (log_unit, '(1X,A,/)') "Transition dipoles calculated using Berry operator formulation"
         CASE (tddfpt_dipole_length)
            WRITE (log_unit, '(1X,A,/)') "Transition dipoles calculated using length formulation"
         CASE (tddfpt_dipole_velocity)
            WRITE (log_unit, '(1X,A,/)') "Transition dipoles calculated using velocity formulation"
         CASE DEFAULT
            CPABORT("Unimplemented form of the dipole operator")
         END SELECT

         WRITE (log_unit, '(T3,A,T15,A,T37,A,T69,A)') "State", "Excitation     S^2", &
            "Transition dipole (a.u.)", "Oscillator"
         WRITE (log_unit, '(T3,A,T14,A,T37,A,T49,A,T61,A,T67,A)') "number", "energy (a.u.)", &
            "x", "y", "z", "strength (a.u.)"
         WRITE (log_unit, '(T3,79("-"))')
      END IF

      CALL get_qs_env(qs_env=qs_env, &
                      matrix_s=matrix_s, &
                      energy=energy)
      gs_energy = energy%total

      nspins = SIZE(gs_mos)
      DO ispin = 1, nspins
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, nrow_global=nao, ncol_global=nmo_occ(ispin))
      END DO
      DO ispin = 1, nspins
         CALL parallel_gemm('N', 'T', nao, nao, nmo_occ(ispin), 1.0_dp, &
                            gs_mos(ispin)%mos_occ, gs_mos(ispin)%mos_occ, &
                            0.0_dp, this%rho_fm_ao(ispin))
      END DO

      ! matrix_s to cp_fm uging tmp matrix hmat_ao_ao
      CALL cp_fm_set_all(matrix=this%hmat_ao_ao, alpha=0.0_dp)
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, this%hmat_ao_ao)

      DO ideriv = 1, 3
         CALL parallel_gemm('T', 'N', tddft_env%nocca, tddft_env%nvirb, tddft_env%nbf, &
                            1.0_dp,  tddft_env%dipole_op_mos_occ(ideriv, 1), &
                                     gs_mos(2)%mos_virt, &
                            0.0_dp, this%dipoles_mo(ideriv))
      END DO

      ! transition dipole moment
      ALLOCATE (trans_dipoles(tddft_env%nstates, 3), source=0.0_dp)
      ! s2_mat
      ALLOCATE (s2_mat(tddft_env%nstates), source=0.0_dp)

      DO istate = 1, tddft_env%nStates
         !--------------------------------------
         ! Orthogonal transform: X vector to AO
         !--------------------------------------
         ! C[nbf,nvirb] = 1.0*A[nbf,nocca]*B[nocca,nvirb] + 0.0*C[nbf,nvirb]
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nvirb, tddft_env%nocca, &
                            1.0_dp, gs_mos(1)%mos_occ, &
                                    tddft_env%evects_mo(istate), &
                            0.0_dp, this%work_ao_vir)
         ! C[nbf,nbf] = 1.0*A[nbf,nvirb]*(B[nbf,nvirb])^T + 0.0*C[nbf,nbf]
         CALL parallel_gemm('N', 'T', tddft_env%nbf, tddft_env%nbf, tddft_env%nvirb, &
                            1.0_dp, this%work_ao_vir, &
                                    gs_mos(2)%mos_virt, &
                            0.0_dp, this%evects_ao_ao(1))

         !--------------------------------------------------------------
         ! Compute AO-basis difference density matrices:
         !   dens_diff_ij = AO transform of -X X^T (occupied block)
         !   dens_diff_ab = AO transform of  X^T X (virtual block)
         !   Delta_rho = dens_diff_ab + dens_diff_ij
         !   where X = tddft_env%evects_mo(istate)
         !--------------------------------------------------------------
         ! Compute dens_diff_ij: (hole contribution, occupied block)
         CALL parallel_gemm('N', 'T', tddft_env%nocca, tddft_env%nocca, tddft_env%nvirb, &
                           -1.0_dp, tddft_env%evects_mo(istate), &
                                    tddft_env%evects_mo(istate), &
                            0.0_dp, this%work_noc_noc)
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nocca, tddft_env%nocca, &
                            1.0_dp, gs_mos(1)%mos_occ, &
                                    this%work_noc_noc, &
                            0.0_dp, this%work_ao_noc)
         CALL parallel_gemm('N', 'T', tddft_env%nbf, tddft_env%nbf, tddft_env%nocca, &
                            1.0_dp, this%work_ao_noc, &
                                    gs_mos(1)%mos_occ, &
                            0.0_dp, this%dens_diff_ij(istate))

         ! Compute dens_diff_ab: (particle contribution, virtual block)
         CALL parallel_gemm('T', 'N', tddft_env%nvirb, tddft_env%nvirb, tddft_env%nocca, &
                            1.0_dp, tddft_env%evects_mo(istate), &
                                    tddft_env%evects_mo(istate), &
                            0.0_dp, this%work_vir_vir)
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nvirb, tddft_env%nvirb, &
                            1.0_dp, gs_mos(2)%mos_virt, &
                                    this%work_vir_vir, &
                            0.0_dp, this%work_ao_vir)
         CALL parallel_gemm('N', 'T', tddft_env%nbf, tddft_env%nbf, tddft_env%nvirb, &
                            1.0_dp, this%work_ao_vir, &
                                    gs_mos(2)%mos_virt, &
                            0.0_dp, this%dens_diff_ab(istate))

         !--------------------------------------------------------------
         ! Compute transition dipoles for state istate
         !--------------------------------------------------------------
         DO ideriv = 1, 3
            CALL cp_fm_trace(tddft_env%evects_mo(istate), &
                             this%dipoles_mo(ideriv), &
                             trans_dipoles(istate, ideriv))
         END DO

         n_spin_imbalance = REAL(tddft_env%nocca - tddft_env%noccb, dp)
         s2_spin = 0.25_dp * n_spin_imbalance * (n_spin_imbalance - 2.0_dp)
         s2_beta_occ = REAL(tddft_env%noccb, dp) + 1.0_dp

         !--------------------------------------------------------------
         ! Compute S^2 expectation value for excited state:
         ! Components:
         ! - s2_cross, s2_anti, s2_offdiag (see formulas below)
         !   (using the relaxed density for the alpha/beta blocks)
         !--------------------------------------------------------------
         ! Step 1
!        !--------------------------------------------------------------
!        ! Compute fully relaxed AO-basis density matrix for the excited state:
!        !
!        !   rho_excited = rho_ground + delta_rho
!        !               = [rho_fm_ao(1) + dens_diff_ij] (alpha)
!        !               = [rho_fm_ao(2) + dens_diff_ab] (beta)
!        !
!        ! If orbital relaxation is included, delta_rho should incorporate those corrections.
!        ! In this context, we compute the "unrelaxed" density as ground + linear response.
!        !--------------------------------------------------------------
!        ! Construct alpha excited state density: work1_ao_ao = rho_fm_ao(1) + dens_diff_ij
         CALL cp_fm_scale_and_add(0.0_dp, this%work1_ao_ao, 1.0_dp, this%rho_fm_ao(1))
         CALL cp_fm_scale_and_add(1.0_dp, this%work1_ao_ao, 1.0_dp, this%dens_diff_ij(istate))
!        ! <--- this%work1_ao_ao now contains AO density for alpha (excited state candidate)

         ! scr1 = S * (rho_alpha_excited) * S
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, this%hmat_ao_ao, &
                                    this%work1_ao_ao, &
                            0.0_dp, this%work2_ao_ao)
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, this%work2_ao_ao, &
                                    this%hmat_ao_ao, &
                            0.0_dp, this%work1_ao_ao)

!        ! Construct beta excited state density: work2_ao_ao = rho_fm_ao(2) + dens_diff_ab
         CALL cp_fm_scale_and_add(0.0_dp, this%work2_ao_ao, 1.0_dp, this%rho_fm_ao(2))
         CALL cp_fm_scale_and_add(1.0_dp, this%work2_ao_ao, 1.0_dp, this%dens_diff_ab(istate))
!        ! <--- this%work2_ao_ao now contains AO density for beta (excited state candidate)

         ! s2_cross = -tr[(rho_beta_excited) * scr1]
         CALL cp_fm_trace(this%work2_ao_ao, this%work1_ao_ao, s2_cross)
         s2_cross = -s2_cross

         ! scr1 = S * dens_diff_ij * S [ta stands for dens_diff_ij, that is, the excited-state difference (alpha block)]
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, this%hmat_ao_ao, &
                                    this%dens_diff_ij(istate), &
                            0.0_dp, this%work1_ao_ao)
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, this%work1_ao_ao, &
                                    this%hmat_ao_ao, &
                            0.0_dp, this%work2_ao_ao)

         ! s2_anti = -tr[dens_diff_ab * scr1]
         CALL cp_fm_trace(this%dens_diff_ab(istate), this%work2_ao_ao, s2_anti)
         s2_anti = -s2_anti

         ! s2_offdiag = (1/2) * tr[symmetrized(evects_ao_ao) * Smat]
         CALL cp_fm_set_all(matrix=this%work1_ao_ao, alpha=0.0_dp)
         CALL cp_fm_transpose(matrix=this%evects_ao_ao(1), matrixt=this%work1_ao_ao)
         CALL cp_fm_scale_and_add(0.0_dp, this%work2_ao_ao, 0.5_dp, this%evects_ao_ao(1))
         CALL cp_fm_scale_and_add(1.0_dp, this%work2_ao_ao, 0.5_dp, this%work1_ao_ao)
         CALL cp_fm_trace(this%work2_ao_ao, this%hmat_ao_ao, s2_offdiag)
         s2_offdiag = s2_offdiag / 2.0_dp

         ! s^2 = s2_spin + s2_beta_occ + s2_cross - s2_anti + s2_offdiag^2
         s2_mat(istate) = s2_spin + s2_beta_occ + s2_cross - s2_anti + s2_offdiag**2
      END DO

      trans_dipoles(:, :) = SQRT(2.0_dp)*trans_dipoles(:, :)

      ! *** summary information ***
      DO istate = 1, tddft_env%nstates
         osc_strength = 2.0_dp/3.0_dp*tddft_env%evals(istate)* &
                        accurate_dot_product(trans_dipoles(istate, 1:3), &
                                             trans_dipoles(istate, 1:3))
         ostrength(istate) = osc_strength
         IF (log_unit > 0) THEN
            WRITE (log_unit, '(1X,I4,1X,F22.14,T30,F4.2,T34,3(1X,ES11.4E2),T70,ES12.4E2)') &
              istate, gs_energy+tddft_env%evals(istate), s2_mat(istate), &
              trans_dipoles(istate, 1:3), osc_strength
         END IF
      END DO
      DEALLOCATE (trans_dipoles, s2_mat)

      ! punch a checksum for the regs
      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T2,A,E14.6)') 'TDDFT : CheckSum  =', &
           SQRT(SUM(tddft_env%evals**2))
      END IF

      max_amplitudes = tddft_env%nocca * tddft_env%nvirb
      ALLOCATE(amplitudes(max_amplitudes))
      ALLOCATE(i_indices(max_amplitudes))
      ALLOCATE(a_indices(max_amplitudes))
      ALLOCATE(indices(max_amplitudes))

      IF (log_unit > 0) THEN
         WRITE (log_unit, "(1X,A)") "", &
            "-------------------------------------------------------------------------------", &
            "-                            Excitation analysis                              -", &
            "-------------------------------------------------------------------------------"
         WRITE (log_unit, '(8X,A,T27,A,T49,A,T69,A)') "State", "Occupied", "Virtual", "Excitation"
         WRITE (log_unit, '(8X,A,T28,A,T49,A,T69,A)') "number", "orbital", "orbital", "amplitude"
         WRITE (log_unit, '(1X,79("-"))')
      END IF

      DO istate = 1, tddft_env%nstates
         IF (log_unit > 0) &
            WRITE (log_unit, '(T7,I8,F10.5,A)') istate, tddft_env%evals(istate)*evolt, " eV"
         n_amplitudes = 0
         DO a = 1, tddft_env%nvirb
            DO i = 1, tddft_env%nocca
               CALL cp_fm_get_element(matrix=tddft_env%evects_mo(istate), &
                                      irow_global=i, &
                                      icol_global=a, &
                                      alpha=val)
               IF (ABS(val) > tddfpt_control%min_excitation_amplitude) THEN
                  n_amplitudes = n_amplitudes + 1
                  amplitudes(n_amplitudes) = -ABS(val)
                  i_indices(n_amplitudes) = i
                  a_indices(n_amplitudes) = a
               END IF
            END DO
         END DO
         IF (n_amplitudes > 0) THEN
            CALL sort(amplitudes(1:n_amplitudes), n_amplitudes, indices)

            DO idx = 1, n_amplitudes
               i = i_indices(indices(idx))
               a = a_indices(indices(idx))
               CALL cp_fm_get_element(matrix=tddft_env%evects_mo(istate), &
                                      irow_global=i, &
                                      icol_global=a, &
                                      alpha=val)
               IF (log_unit > 0) &
                  WRITE (log_unit, '(T27,I8,1X,A5,T48,I8,1X,A5,T70,F9.6)') &
                  i, '(alp)', tddft_env%noccb + a, '(bet)', val
            END DO
         END IF
      END DO
      IF (log_unit > 0) THEN
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      DEALLOCATE(amplitudes, i_indices, a_indices, indices)

      ALLOCATE(trans_dipoles_exc_exc(tddft_env%nstates, tddft_env%nstates, 3), source=0.0_dp)

!     DO istate = 1, tddft_env%nstates
!        DO jstate = istate + 1, tddft_env%nstates
!           ! Calculate difference of excitation vectors: X_I - X_J
!           CALL cp_fm_scale_and_add(0.0_dp, this%work_noc_vir, &
!                                    1.0_dp, tddft_env%evects_mo(istate))
!           CALL cp_fm_scale_and_add(1.0_dp, this%work_noc_vir, &
!                                    -1.0_dp, tddft_env%evects_mo(jstate))
!           ! Calculating transition dipole moments
!           DO ideriv = 1, 3
!              CALL cp_fm_trace(this%work_noc_vir, &
!                               this%dipoles_mo(ideriv), &
!                               trans_dipoles_exc_exc(istate, jstate, ideriv))
!           END DO
!        END DO
!     END DO

      IF (log_unit > 0 .and. .FALSE.) THEN
         WRITE (log_unit, "(1X,A)") "", &
            "-------------------------------------------------------------------------------", &
            "-                  Transition Dipoles Between Excited States                  -", &
            "-------------------------------------------------------------------------------"
         WRITE (log_unit, '(12X,A,6X,A,13X,A)') "Excitation", "Transition dipole, a.u.", "Oscillator"
         WRITE (log_unit, '(1X,A,5X,A,12X,A,11X,A,11X,A,7X,A)') "Transition", "eV", "x", "y", "z", "strength (a.u.)"
         WRITE (log_unit, '(1X,79("-"))')

         DO istate = 1, tddft_env%nstates
            DO jstate = istate + 1, tddft_env%nstates
               energ = tddft_env%evals(jstate) - tddft_env%evals(istate)
               f = 2.0_dp / 3.0_dp * energ * SUM(trans_dipoles_exc_exc(istate, jstate, :)**2)
               WRITE (log_unit, '(3X,I0,1X,"->",1X,I0,T11,F11.5,3X,3(1X,ES11.4E2),2X,ES12.4E2)') &
                  istate, jstate, energ * evolt, trans_dipoles_exc_exc(istate, jstate, 1:3), f
            END DO
         END DO
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      DEALLOCATE(trans_dipoles_exc_exc)
      DO istate = 1, tddft_env%nStates
         CALL cp_fm_set_all(matrix=this%work_noc_vir, alpha=0.0_dp)
         CALL cp_fm_to_fm(source=tddft_env%evects_mo(istate), &
                          destination=this%work_noc_vir)
         CALL parallel_gemm('N', 'T', tddft_env%nbf, tddft_env%nocca, tddft_env%nvirb, &
                         1.0_dp, gs_mos(2)%mos_virt, &
                                 this%work_noc_vir, &
                         0.0_dp, tddft_env%evects(1, istate))

      END DO

      suffix = TRIM("sfl-")
      tddfpt_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT")
      tddfpt_print_section => section_vals_get_subs_vals(tddfpt_section, "PRINT")
      CALL tddfpt_print_nto_analysis(qs_env, &
                                     tddft_env%evects, tddft_env%evals, &
                                     ostrength, &
                                     gs_mos, &
                                     matrix_s(1)%matrix, &
                                     tddfpt_print_section, suffix)

      CALL timestop(handle)

   END SUBROUTINE sf_print_theory_data

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

END MODULE qs_tddft_theory_sf
