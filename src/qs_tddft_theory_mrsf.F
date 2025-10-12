!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Mixed Reference Spin-Flip TDDFT (MRSF-TDDFT) theory implementation
!> \details This module implements the MRSF-TDDFT method for computing electronic excited states
!>          with pure spin separation. MRSF-TDDFT addresses spin contamination issues in
!>          conventional spin-flip TDDFT by using a multi-reference approach that yields pure
!>          spin states (singlets and triplets) from a mixed reference configuration.
!>
!> \par Key features
!>      - Pure spin states without contamination
!>      - Balanced treatment of static and dynamic correlation
!>      - Configurable spin-pair coupling parameters
!>      - Support for range-separated hybrid functionals
!>      - ADMM acceleration for Hartree-Fock exchange kernel
!>      - Comprehensive natural transition orbital analysis
!>
!> \author Konstantin Komarov (constlike@gmail.com) (Luber group)
!> \date Aug 2025
! **************************************************************************************************
MODULE qs_tddft_theory_mrsf

   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_complete_redistribute,&
                                              dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_distribution_type,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type_antisymmetric
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_frobenius_norm,&
                                              cp_fm_norm,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace,&
                                              cp_fm_transpose
   USE cp_fm_types,                     ONLY: cp_fm_create, &
                                              cp_fm_get_element,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_maxabsval,&
                                              cp_fm_release, &
                                              cp_fm_set_all,&
                                              cp_fm_set_element,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE hfx_types,                       ONLY: hfx_type
   USE input_constants,                 ONLY: tddfpt_dipole_berry,&
                                              tddfpt_dipole_length,&
                                              tddfpt_dipole_velocity
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kahan_sum,                       ONLY: accurate_dot_product
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE machine,                         ONLY: m_flush
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_scale,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild,&
                                              qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_tddfpt2_properties,           ONLY: tddfpt_print_nto_analysis
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_dbcsr_create_by_dist
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE qs_tddft_theory,                 ONLY: tddft_theory_type
   USE qs_tddft_types,                  ONLY: tddft_info_type
   USE qs_tddft_utils,                  ONLY: allocate_fm_matrix,&
                                              deallocate_fm_matrix
   USE string_utilities,                ONLY: integer_to_string
   USE task_list_types,                 ONLY: task_list_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_theory_mrsf'

   PUBLIC :: tddft_mrsf_type

   ! Local constants for MRSF calculations
   REAL(Kind=dp), PARAMETER, PRIVATE :: sqrt12 = 1.0_dp/SQRT(2.0_dp)
   REAL(KIND=dp), PARAMETER, PRIVATE :: threshold = 1.0E-10_dp

! **************************************************************************************************
!> \brief MRSF-TDDFT derived type implementing the abstract theory interface
!> \details Contains all working arrays and procedures specific to mixed reference spin-flip
!>          TDDFT calculations. The implementation handles complex multi-reference transformations
!>          with spin-pair coupling parameters.
! **************************************************************************************************
   TYPE, EXTENDS(tddft_theory_type) :: tddft_mrsf_type

      ! Theory-specific parameters
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: spc   ! [7, 2] spin-pair coupling

      ! AO-MO transformation working matrices
      TYPE(cp_fm_type), POINTER                          :: work_nao_vir => NULL()   ! [nbf, nvirb]
      TYPE(cp_fm_type), POINTER                          :: work_noc_nao => NULL()   ! [nocca, nbf]
      TYPE(cp_fm_type), POINTER                          :: work_noc_noc => NULL()   ! [nocca, nocca]
      TYPE(cp_fm_type), POINTER                          :: work_nao_noc => NULL()   ! [nbf, nocca]
      TYPE(cp_fm_type), POINTER                          :: work_vir_vir => NULL()   ! [nvirb, nvirb]
      TYPE(cp_fm_type), POINTER                          :: work_noc_vir => NULL()   ! [nocca, nvirb]
      TYPE(cp_fm_type), POINTER                          :: work2_noc_vir => NULL()  ! [nocca, nvirb]

      ! General AO working matrices
      TYPE(cp_fm_type), POINTER                          :: work1_nao_nao => NULL()  ! [nbf, nbf]
      TYPE(cp_fm_type), POINTER                          :: work2_nao_nao => NULL()  ! [nbf, nbf]
      TYPE(cp_fm_type), POINTER                          :: work1_nao_1  => NULL()   ! [nbf, 1]
      TYPE(cp_fm_type), POINTER                          :: work2_nao_1  => NULL()   ! [nbf, 1]

      ! MRSF-specific working arrays for complex transformations
      TYPE(cp_fm_type), POINTER                          :: work_vira_1 => NULL()    ! [nvira, 1]
      TYPE(cp_fm_type), POINTER                          :: work_nocb_1 => NULL()    ! [noccb, 1]
      TYPE(cp_fm_type), POINTER                          :: work_nocb_vira => NULL() ! [noccb, nvira]
      TYPE(cp_fm_type), POINTER                          :: work_nao_nocb => NULL()  ! [nbf, noccb]

      ! ADMM working matrices
      TYPE(cp_fm_type), POINTER                          :: work_aux_aux => NULL()   ! [nao_aux, nao_aux]
      TYPE(cp_fm_type), POINTER                          :: work_aux_nao => NULL()   ! [nao_aux, nao]
      TYPE(cp_fm_type), POINTER                          :: work_nao_aux => NULL()   ! [nao, nao_aux]

      ! Coulomb kernel working arrays
      TYPE(pw_r3d_rs_type), POINTER                      :: work_v_rspace => NULL()  ! real-space potential
      TYPE(pw_c1d_gs_type), POINTER                      :: work_v_gspace => NULL()  ! reciprocal-space potential
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: hmat_coulomb => NULL()   ! Coulomb intermediate
      TYPE(qs_rho_type), POINTER                         :: rho_nao_struct => NULL() ! Density structure

   CONTAINS
      PROCEDURE, PASS(this) :: allocate   => mrsf_allocate
      PROCEDURE, PASS(this) :: init_guess => mrsf_init_guess
      PROCEDURE, PASS(this) :: mo_to_ao   => mrsf_mo_to_ao
      PROCEDURE, PASS(this) :: apply      => mrsf_apply
      PROCEDURE, PASS(this) :: ao_to_mo   => mrsf_ao_to_mo
      PROCEDURE, PASS(this) :: print      => mrsf_print
      PROCEDURE, PASS(this) :: deallocate => mrsf_deallocate
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief Allocate all working arrays for MRSF-TDDFT calculations
!> \param[inout] this      MRSF theory object
!> \param[in]    qs_env    QuickStep environment
!> \param[inout] tddft_env TDDFT environment containing system dimensions
!> \details Sets up spin-pair coupling parameters and allocates complex transformation arrays
! **************************************************************************************************
   SUBROUTINE mrsf_allocate(this, qs_env, tddft_env)
      CLASS(tddft_mrsf_type), INTENT(INOUT)              :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mrsf_allocate'

      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sab_aux_fit
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_s_aux_fit
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(section_vals_type), POINTER                   :: tddfpt_section
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      LOGICAL                                            :: explicit
      INTEGER                                            :: handle, ispin

      CALL timeset(routineN, handle)

      ! Get input sections and control parameters
      NULLIFY (tddfpt_section)
      tddfpt_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT")

      CALL get_qs_env(qs_env, dft_control=dft_control)
      NULLIFY (tddfpt_control)
      tddfpt_control => dft_control%tddfpt2_control

      ! Initialize spin-pair coupling (SPC) parameter matrix
      ALLOCATE (this%spc(7, 2), source=0.0_dp)

      ! Extract SPC parameters from input or set defaults
      CALL section_vals_val_get(tddfpt_section, "SPIN_PAIR_COUPLING", explicit=explicit)

      IF (explicit) THEN
         IF (tddft_env%do_hfx_cam) THEN
            ! Range-separated functionals: separate short-range and long-range parameters
            this%spc(1:4, 1:2) = tddfpt_control%mrsf_SPC(3)  ! CO-OV coupling
            this%spc(5, 1:2) = tddfpt_control%mrsf_SPC(2)    ! OV-OV coupling
            this%spc(6, 1:2) = tddfpt_control%mrsf_SPC(1)    ! CO-CO coupling
         ELSE
            ! Standard functionals: single parameter set
            this%spc(1:4, 1) = tddfpt_control%mrsf_SPC(3)
            this%spc(5, 1) = tddfpt_control%mrsf_SPC(2)
            this%spc(6, 1) = tddfpt_control%mrsf_SPC(1)
         END IF
      ELSE
         ! Default coupling parameters
         IF (tddft_env%do_hfx_cam) THEN
            this%spc(1:6, 1:2) = 0.5_dp
         ELSE
            this%spc(1:6, 1) = 0.5_dp
         END IF
      END IF

      ! Extract HF exchange scaling parameters
      CALL section_vals_val_get(tddfpt_section, "MRSF_HF", explicit=explicit)
      IF (explicit) THEN
         IF (tddft_env%do_hfx_cam) THEN
            this%spc(7, 1) = tddfpt_control%mrsf_HF(1)  ! Alpha parameter
            this%spc(7, 2) = tddfpt_control%mrsf_HF(2)  ! Beta parameter
         ELSE
            this%spc(7, 1) = tddfpt_control%mrsf_HF(1)
         END IF
      ELSE
         ! Use functional HF fractions
         IF (tddft_env%do_hfx_cam) THEN
            this%spc(7, 1) = tddft_env%C_hf_alpha
            this%spc(7, 2) = tddft_env%C_hf_beta
         ELSE
            this%spc(7, 1) = tddft_env%C_hf_alpha
         END IF
      END IF

      ! Apply theory-specific modifications to subspace dimensions
      SELECT CASE (tddft_env%response_type)
      CASE (2) ! MRSF-S: Singlet states
         tddft_env%nStates = MIN(tddft_env%nstates, tddft_env%nstates_max - 1)
         tddft_env%mxvec = MIN(tddft_env%maxvec*tddft_env%nstates, &
                               tddft_env%nstates_max - 1)

      CASE (3) ! MRSF-T: Triplet states
         tddft_env%nStates = MIN(tddft_env%nstates, tddft_env%nstates_max - 3)
         tddft_env%mxvec = MIN(tddft_env%maxvec*tddft_env%nstates, &
                               tddft_env%nstates_max - 3)
         ! Sign flip for triplet coupling
         this%spc(1:6, 1:2) = -this%spc(1:6, 1:2)
      END SELECT

      !--------------------
      ! ALLOCATE  FM
      !--------------------
      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env)

      !=== hmat_nao_nao (7) [nbf, nbf] ===================
      CALL allocate_fm_matrix(matrix=this%hmat_nao_nao, &
                              nmat=7, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nbf, &
                              blacs_env=blacs_env)

      !=== evects_nao_nao (7) [nbf, nbf] =================
      CALL allocate_fm_matrix(matrix=this%evects_nao_nao, &
                              nmat=7, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nbf, &
                              blacs_env=blacs_env)

      !=== dipoles_mo (3) [nocca, nvirb] ===============
      CALL allocate_fm_matrix(matrix=this%dipoles_mo, &
                              nmat=3, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      !=== fock (nSpins) [nbf, nvirb] ==================
      CALL allocate_fm_matrix(matrix=this%fock, &
                              nmat=tddft_env%nSpins, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nbf, &
                              blacs_env=blacs_env)

      !=== focka_noc_noc [nocca, nocca] ================
      CALL allocate_fm_matrix(matrix=this%focka_noc_noc, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nocca, &
                              blacs_env=blacs_env)

      !=== fockb_vir_vir [nvirb, nvirb] ================
      CALL allocate_fm_matrix(matrix=this%fockb_vir_vir, &
                              nrow=tddft_env%nvirb, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      !=== work_nao_vir [nbf, nvirb] ====================
      CALL allocate_fm_matrix(matrix=this%work_nao_vir, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      !=== work_noc_nao [nocca, nbf] ====================
      CALL allocate_fm_matrix(matrix=this%work_noc_nao, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nbf, &
                              blacs_env=blacs_env)

      !=== work_noc_noc [nocca, nocca] =================
      CALL allocate_fm_matrix(matrix=this%work_noc_noc, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nocca, &
                              blacs_env=blacs_env)

      !=== work_nao_noc [nbf, nocca] ====================
      CALL allocate_fm_matrix(matrix=this%work_nao_noc, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nocca, &
                              blacs_env=blacs_env)

      !=== work_vir_vir [nvirb, nvirb] =================
      CALL allocate_fm_matrix(matrix=this%work_vir_vir, &
                              nrow=tddft_env%nvirb, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      !=== work_noc_vir [nocca, nvirb] =================
      CALL allocate_fm_matrix(matrix=this%work_noc_vir, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      !=== work2_noc_vir [nocca, nvirb] =================
      CALL allocate_fm_matrix(matrix=this%work2_noc_vir, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      !=== work1_nao_nao [nbf, nbf] ======================
      CALL allocate_fm_matrix(matrix=this%work1_nao_nao, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nbf, &
                              blacs_env=blacs_env)

      !=== work2_nao_nao [nbf, nbf] ======================
      CALL allocate_fm_matrix(matrix=this%work2_nao_nao, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nbf, &
                              blacs_env=blacs_env)

      !=== work1_nao_1 [nbf, 1] =========================
      CALL allocate_fm_matrix(matrix=this%work1_nao_1, &
                              nrow=tddft_env%nbf, &
                              ncol=1, &
                              blacs_env=blacs_env)

      !=== work2_nao_1 [nbf, 1] =========================
      CALL allocate_fm_matrix(matrix=this%work2_nao_1, &
                              nrow=tddft_env%nbf, &
                              ncol=1, &
                              blacs_env=blacs_env)

      !=== work_vira_1 [nvira, 1] ======================
      CALL allocate_fm_matrix(matrix=this%work_vira_1, &
                              nrow=tddft_env%nvira, &
                              ncol=1, &
                              blacs_env=blacs_env)

      !=== work_nocb_1 [noccb, 1] ======================
      CALL allocate_fm_matrix(matrix=this%work_nocb_1, &
                              nrow=tddft_env%noccb, &
                              ncol=1, &
                              blacs_env=blacs_env)

      !=== work_nocb_vira [noccb, nvira] ===============
      CALL allocate_fm_matrix(matrix=this%work_nocb_vira, &
                              nrow=tddft_env%noccb, &
                              ncol=tddft_env%nvira, &
                              blacs_env=blacs_env)

      !=== work_nao_nocb [nbf, noccb] ===================
      CALL allocate_fm_matrix(matrix=this%work_nao_nocb, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%noccb, &
                              blacs_env=blacs_env)

      IF (tddft_env%do_admm) THEN
         !=== work_aux_aux [nbf_aux, nbf_aux] ==========
         CALL allocate_fm_matrix(matrix=this%work_aux_aux, &
                                nrow=tddft_env%nbf_aux, &
                                ncol=tddft_env%nbf_aux, &
                                blacs_env=blacs_env)

         !=== work_aux_nao [nbf_aux, nbf] ==============
         CALL allocate_fm_matrix(matrix=this%work_aux_nao, &
                                nrow=tddft_env%nbf_aux, &
                                ncol=tddft_env%nbf, &
                                blacs_env=blacs_env)

         !=== work_nao_aux [nbf, nbf_aux] ==============
         CALL allocate_fm_matrix(matrix=this%work_nao_aux, &
                                nrow=tddft_env%nbf, &
                                ncol=tddft_env%nbf_aux, &
                                blacs_env=blacs_env)
      END IF

      !--------------------
      ! ALLOCATE  DBCSR
      !--------------------

      CALL get_qs_env(qs_env=qs_env, &
                      matrix_s=matrix_s, &
                      dbcsr_dist=dbcsr_dist, &
                      sab_orb=sab_orb)

      NULLIFY (this%hmat_coulomb)
      CALL dbcsr_allocate_matrix_set(this%hmat_coulomb, 1)
      DO ispin = 1, 1
         CALL dbcsr_init_p(this%hmat_coulomb(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(this%hmat_coulomb(ispin)%matrix, &
                                          template=matrix_s(1)%matrix, &
                                          dbcsr_dist=dbcsr_dist, &
                                          sab=sab_orb)
      END DO

      NULLIFY (this%work_rho_nao_nao_symm)
      CALL dbcsr_allocate_matrix_set(this%work_rho_nao_nao_symm, tddft_env%nSpins)
      NULLIFY (this%work_rho_nao_nao_asymm)
      CALL dbcsr_allocate_matrix_set(this%work_rho_nao_nao_asymm, tddft_env%nSpins)
      DO ispin = 1, tddft_env%nSpins
         CALL dbcsr_init_p(this%work_rho_nao_nao_symm(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist( &
            this%work_rho_nao_nao_symm(ispin)%matrix, &
            template=matrix_s(1)%matrix, &
            dbcsr_dist=dbcsr_dist, &
            sab=sab_orb)

         CALL dbcsr_init_p(this%work_rho_nao_nao_asymm(ispin)%matrix)
         CALL dbcsr_create(this%work_rho_nao_nao_asymm(ispin)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric, &
                           template=this%work_rho_nao_nao_symm(ispin)%matrix)
         CALL dbcsr_complete_redistribute(this%work_rho_nao_nao_symm(ispin)%matrix, &
                                          this%work_rho_nao_nao_asymm(ispin)%matrix)
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
         CALL dbcsr_create(this%work_hmat_asymm(ispin)%matrix, matrix_type=dbcsr_type_antisymmetric, &
                           template=this%work_hmat_symm(ispin)%matrix)
         CALL dbcsr_complete_redistribute(this%work_hmat_symm(ispin)%matrix, &
                                          this%work_hmat_asymm(ispin)%matrix)
      END DO

      IF (tddft_env%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         CALL get_admm_env(admm_env, &
                           matrix_s_aux_fit=matrix_s_aux_fit, &
                           sab_aux_fit=sab_aux_fit)

         NULLIFY (this%work_rho_aux_aux_symm)
         CALL dbcsr_allocate_matrix_set(this%work_rho_aux_aux_symm, tddft_env%nSpins)
         NULLIFY (this%work_rho_aux_aux_asymm)
         CALL dbcsr_allocate_matrix_set(this%work_rho_aux_aux_asymm, tddft_env%nSpins)
         DO ispin = 1, tddft_env%nSpins
            CALL dbcsr_init_p(this%work_rho_aux_aux_symm(ispin)%matrix)
            CALL tddfpt_dbcsr_create_by_dist( &
               this%work_rho_aux_aux_symm(ispin)%matrix, &
               template=matrix_s_aux_fit(1)%matrix, &
               dbcsr_dist=dbcsr_dist, &
               sab=sab_aux_fit)

            CALL dbcsr_init_p(this%work_rho_aux_aux_asymm(ispin)%matrix)
            CALL dbcsr_create(this%work_rho_aux_aux_asymm(ispin)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric, &
                              template=this%work_rho_aux_aux_symm(ispin)%matrix)
            CALL dbcsr_complete_redistribute(this%work_rho_aux_aux_symm(ispin)%matrix, &
                                             this%work_rho_aux_aux_asymm(ispin)%matrix)
         END DO

         NULLIFY (this%work_hmat_aux_symm)
         CALL dbcsr_allocate_matrix_set(this%work_hmat_aux_symm, tddft_env%nSpins)
         NULLIFY (this%work_hmat_aux_asymm)
         CALL dbcsr_allocate_matrix_set(this%work_hmat_aux_asymm, tddft_env%nSpins)
         DO ispin = 1, tddft_env%nSpins
            CALL dbcsr_init_p(this%work_hmat_aux_symm(ispin)%matrix)
            CALL tddfpt_dbcsr_create_by_dist(this%work_hmat_aux_symm(ispin)%matrix, &
                                             template=matrix_s_aux_fit(1)%matrix, &
                                             dbcsr_dist=dbcsr_dist, &
                                             sab=sab_aux_fit)

            CALL dbcsr_init_p(this%work_hmat_aux_asymm(ispin)%matrix)
            CALL dbcsr_create(this%work_hmat_aux_asymm(ispin)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric, &
                              template=this%work_hmat_aux_symm(ispin)%matrix)
            CALL dbcsr_complete_redistribute(this%work_hmat_aux_symm(ispin)%matrix, &
                                             this%work_hmat_aux_asymm(ispin)%matrix)
         END DO
      END IF

      !--------------------
      ! ALLOCATE  Rho (Coulomb part)
      !--------------------

      NULLIFY (this%rho_nao_struct)
      ALLOCATE (this%rho_nao_struct)
      CALL qs_rho_create(this%rho_nao_struct)
      CALL qs_rho_set(this%rho_nao_struct, rho_ao=this%work_rho_nao_nao_symm)
      CALL qs_rho_rebuild(this%rho_nao_struct, qs_env, &
                          rebuild_ao=.FALSE., &
                          rebuild_grids=.TRUE.)

      !--------------------
      ! ALLOCATE  PW
      !--------------------

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool)

      ALLOCATE (this%work_v_gspace)
      CALL auxbas_pw_pool%create_pw(this%work_v_gspace)
      ALLOCATE (this%work_v_rspace)
      CALL auxbas_pw_pool%create_pw(this%work_v_rspace)

      CALL timestop(handle)

   END SUBROUTINE mrsf_allocate

! **************************************************************************************************
!> \brief Generate initial guess vectors for MRSF-TDDFT with mixed reference configurations
!> \param[inout] this      MRSF theory object
!> \param[in]    qs_env    QuickStep environment
!> \param[inout] tddft_env TDDFT environment
!> \param[in]    gs_mos    ground state molecular orbitals
!> \param[in]    log_unit  output unit for logging
! **************************************************************************************************
   SUBROUTINE mrsf_init_guess(this, qs_env, tddft_env, gs_mos, log_unit)
      CLASS(tddft_mrsf_type), INTENT(INOUT)              :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), POINTER, INTENT(IN)               :: gs_mos
      INTEGER, INTENT(IN)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mrsf_init_guess'

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_m_occ
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      LOGICAL                                            :: do_triplets
      INTEGER                                            :: handle, imo_occ, imo_virt, &
                                                            istate, lr1, lr2, ispin, &
                                                            nStates_selected, istate_orig

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks)

      ! Transform Fock matrices to MO basis
      DO ispin = 1, tddft_env%nSpins
         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, this%fock(ispin))
         ! work1_nao_nao = gs_mos(1)%mos_full^T * this%fock
         CALL parallel_gemm('T', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, gs_mos(ispin)%mos_full, &
                                    this%fock(ispin), &
                            0.0_dp, this%work1_nao_nao)

         ! this%fock = work_nao_nao * gs_mos(2)%mos_full
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, this%work1_nao_nao, &
                                    gs_mos(ispin)%mos_full, &
                            0.0_dp, this%fock(ispin))
      END DO

      ! Extract Fock matrix blocks
      CALL cp_fm_to_fm_submat(this%fock(1), this%focka_noc_noc, &
                              nrow=tddft_env%nocca, ncol=tddft_env%nocca, &
                              s_firstrow=1, &
                              s_firstcol=1, &
                              t_firstrow=1, t_firstcol=1)
      CALL cp_fm_to_fm_submat(this%fock(2), this%fockb_vir_vir, &
                              nrow=tddft_env%nvirb, &
                              ncol=tddft_env%nvirb, &
                              s_firstrow=tddft_env%noccb + 1, &
                              s_firstcol=tddft_env%noccb + 1, &
                              t_firstrow=1, t_firstcol=1)

      ! Print MRSF parameter information
      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T2,A)') REPEAT("-", 79)
         WRITE (log_unit, '(T2,A)') "-                   MRSF-TDDFT Fitting Parameters                             -"
         WRITE (log_unit, '(T2,A)') REPEAT("-", 79)
         IF (tddft_env%do_hfx_cam) THEN
            WRITE (log_unit, '(T14,A,F6.3)') "Fraction of short-range HF exchange (alpha) ........", &
               this%spc(7, 1)
            WRITE (log_unit, '(T14,A,F6.3)') "Fraction of long-range HF exchange (beta) ..........", &
               this%spc(7, 2)
            WRITE (log_unit, '(T14,A,F6.3)') "Total long-range HF exchange (alpha + beta) ........", &
               this%spc(7, 1) + this%spc(7, 2)
            WRITE (log_unit, '(T14,A,F6.3)') "Range-separation parameter (omega) .................", &
               tddft_env%C_hf_omega
         ELSE
            WRITE (log_unit, '(T14,A,F6.3)') "Fraction of Hartree-Fock exchange ........", &
               this%spc(7, 1)
         END IF
         WRITE (log_unit, '(/,T14,A)') "Spin-Pair Coupling parameters:"
         WRITE (log_unit, '(T31,A,T41,A,T51,A)') "CO-CO  |", "OV-OV  |", "CO-OV"
         WRITE (log_unit, '(T30,F6.3,A,T40,F6.3,A,T50,F6.3)') &
            ABS(this%spc(6, 1)), "  |", ABS(this%spc(5, 1)), "  |", ABS(this%spc(1, 1))
         WRITE (log_unit, '(T2,A)') REPEAT("-", 79)
         WRITE (log_unit, '(T2,A,/,T2,A)') &
            "- Legend: CO = Closed-to-Open,  OV = Open-to-Virtual configurations           -", &
            "-         CO-OV = Coupling between CO and OV configurations                   -"
         WRITE (log_unit, '(T2,A)') REPEAT("-", 79)
      END IF

      ! Create energy-ordered excitation list with MRSF modifications
      nStates_selected = tddft_env%nocca*(tddft_env%nbf - tddft_env%noccb)

      ALLOCATE (inds(nStates_selected))
      ALLOCATE (e_virt_m_occ(nStates_selected))

      lr1 = tddft_env%noccb + 1   ! lowest singly occupied
      lr2 = tddft_env%nocca       ! highest singly occupied
      do_triplets = tddft_env%response_type == 3

      ! Compute excitation energies with MRSF-specific modifications
      istate = 0
      DO imo_virt = tddft_env%noccb + 1, tddft_env%nbf
         DO imo_occ = 1, tddft_env%nocca
            istate = (imo_virt - tddft_env%noccb - 1)*tddft_env%nocca + imo_occ

            ! Special handling for critical lr1-lr1 transition
            IF (imo_occ == lr1 .AND. imo_virt == lr1) THEN
               tddft_env%e_virt_minus_occ(istate) = (gs_mos(2)%evals_full(lr1) - gs_mos(1)%evals_full(lr1) &
                                                     + gs_mos(2)%evals_full(lr2) - gs_mos(1)%evals_full(lr2))*0.5_dp
               CYCLE
            END IF

            ! Standard orbital energy difference
            tddft_env%e_virt_minus_occ(istate) = gs_mos(2)%evals_full(imo_virt) - gs_mos(1)%evals_full(imo_occ)

            ! Exclude problematic transitions to maintain MRSF purity
            IF (imo_occ == lr2 .AND. imo_virt == lr2) THEN
               tddft_env%e_virt_minus_occ(istate) = HUGE(1.0_dp)
            ELSE IF (do_triplets .AND. imo_occ == lr2 .AND. imo_virt == lr1) THEN
               tddft_env%e_virt_minus_occ(istate) = HUGE(1.0_dp)
            ELSE IF (do_triplets .AND. imo_occ == lr1 .AND. imo_virt == lr2) THEN
               tddft_env%e_virt_minus_occ(istate) = HUGE(1.0_dp)
            END IF
         END DO
      END DO

      e_virt_m_occ(1:nStates_selected) = tddft_env%e_virt_minus_occ(1:nStates_selected)

      ! Sort excitations by energy
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

      ! Set up initial guess vectors
      DO istate = 1, tddft_env%nVec
         tddft_env%evals(istate) = e_virt_m_occ(istate)
         istate_orig = inds(istate)
         imo_occ = MOD(istate_orig - 1, tddft_env%nocca) + 1
         imo_virt = (istate_orig - 1)/tddft_env%nocca + 1
         CALL cp_fm_set_element(tddft_env%evects_mo(istate), imo_occ, imo_virt, 1.0_dp)
         IF (istate > tddft_env%nStates) CYCLE
         IF (log_unit > 0) &
            WRITE (log_unit, '(T7,I8,T24,I8,T37,A5,T45,I8,T54,A5,T60,F14.5)') &
            istate, imo_occ, '(adp)', tddft_env%noccb + imo_virt, '(adp)', e_virt_m_occ(istate)*evolt
      END DO

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T7,A,T50,I24)') 'Number of active states:', &
            SIZE(gs_mos(1)%evals_occ)*SIZE(gs_mos(2)%evals_virt)
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      DEALLOCATE (inds, e_virt_m_occ)

      CALL timestop(handle)

   END SUBROUTINE mrsf_init_guess

! **************************************************************************************************
!> \brief Transform MRSF trial vectors from MO to AO basis with multi-reference expansion
!> \param[inout] this      MRSF theory object
!> \param[inout] tddft_env TDDFT environment
!> \param[in]    gs_mos    ground state molecular orbitals
!> \param[in]    ivec      index of vector to transform
!> \details Complex multi-reference transformation handling different coupling types and
!>          reference configurations. Creates multiple AO matrices for different excitation types.
! **************************************************************************************************
   SUBROUTINE mrsf_mo_to_ao(this, tddft_env, gs_mos, ivec)
      CLASS(tddft_mrsf_type), INTENT(INOUT)              :: this
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER, INTENT(IN)                             :: gs_mos
      INTEGER, INTENT(IN)                                :: ivec

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mrsf_mo_to_ao'

      REAL(KIND=dp)                                      :: val
      INTEGER                                            :: handle, i, j, &
         lr1, lr2, nbf, nocca, noccb, nvira, nvirb, shift

      CALL timeset(routineN, handle)

      ! Extract system dimensions
      nbf = tddft_env%nbf
      nocca = tddft_env%nocca
      noccb = tddft_env%noccb
      nvira = tddft_env%nvira  ! nvira = nbf - nocca
      nvirb = tddft_env%nvirb  ! nvirb = nbf - noccb

      ! Define critical orbital indices and shift parameter
      lr1 = nocca - 1  ! lowest singly occupied
      lr2 = nocca      ! highest singly occupied
      shift = nocca - noccb

      ! Initialize all AO transformation matrices
      DO i = 1, 7
         CALL cp_fm_set_all(this%evects_nao_nao(i), 0.0_dp)
      END DO

      ! Transform various excitation types to AO basis
      ! The implementation creates different AO matrices for:
      ! evects_nao_nao(1) - bo2v: closed-shell occupied to open-shell virtual (type 1)
      ! evects_nao_nao(2) - bo1v: closed-shell occupied to open-shell virtual (type 2)
      ! evects_nao_nao(3) - bco1: open-shell configurations (type 1)
      ! evects_nao_nao(4) - bco2: open-shell configurations (type 2)
      ! evects_nao_nao(5) - o21v: cross-coupling terms
      ! evects_nao_nao(6) - co12: cross-coupling terms
      ! evects_nao_nao(7) - ball: combined transformation

      ! Dev note, transform needed row to coulumn, not a transpose.
      ! nocca = noccb + 2
      ! noccb = nocca - 2
      ! nvira = nvirb - 2
      ! nvirb = nvira + 2

      ! work[1:nvira, 1] => X(1:nocca, 1:nvirb) [lr2, 1+shift:nvirb]
      CALL cp_fm_set_all(this%work_vira_1, 0.0_dp)
      DO i = 1, nvira
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), lr2, shift + i, val)
         CALL cp_fm_set_element(this%work_vira_1, i, 1, val)
      END DO
      ! C[1:nbf,1] = 1.0*A[1:nbf,nocca+1:nbf]*B[1:nvira,1] + 0.0*C[1:nbf,1] or
      ! C[1:nbf,1] = 1.0*A[1:nbf,1:nvira]*B[1:nvira,1] + 0.0*C[1:nbf,1]
      CALL parallel_gemm('N', 'N', nbf, 1, nvira, &
                         1.0_dp, gs_mos(1)%mos_virt, &   ! MOs beta (=alpha)
                                 this%work_vira_1, &     ! X[shift+1:nvirb, lr2]
                         0.0_dp, this%work2_nao_1)       ! tmp
      ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,lr2]*(B[1:nbf,1])^T + 1.0*C[1:nbf,1:nbf]
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                         1.0_dp, gs_mos(1)%mos_full, &   ! MOs alpha
                                 this%work2_nao_1, &     ! tmp
                         1.0_dp, this%evects_nao_nao(1), & ! bo2v
                         a_first_col=lr2)
      ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,1]*(B[1:nbf,lr1])^T + 1.0*C[1:nbf,1:nbf]
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                         1.0_dp, this%work2_nao_1, &     ! tmp
                                 gs_mos(1)%mos_full, &   ! MOs alpha
                         1.0_dp, this%evects_nao_nao(5), & ! o21v
                         b_first_col=lr1)

      ! work[1:nvira, 1] => X(1:nocca, 1:nvirb) [nocca-1, 1+shift:nvirb]
      CALL cp_fm_set_all(this%work_vira_1, 0.0_dp)
      DO i = 1, nvira
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), lr1, shift + i, val)
         CALL cp_fm_set_element(this%work_vira_1, i, 1, val)
      END DO
      ! C[1:nbf,1] = 1.0*A[1:nbf,nocca+1:nbf]*B[1:nvira,1] + 0.0*C[1:nbf,1] or
      ! C[1:nbf,1] = 1.0*A[1:nbf,1:nvira]*B[1:nvira,1] + 0.0*C[1:nbf,1]
      CALL parallel_gemm('N', 'N', nbf, 1, nvira, &
                         1.0_dp, gs_mos(1)%mos_virt, &   ! MOs beta (=alpha)
                                 this%work_vira_1, &     ! X[shift+1:nvirb, lr1]
                         0.0_dp, this%work2_nao_1)       ! tmp
      ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,lr1]*(B[1:nbf,1])^T + 1.0*C[1:nbf,1:nbf]
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                         1.0_dp, gs_mos(1)%mos_full, &   ! MOs alpha
                                 this%work2_nao_1, &     ! tmp
                         1.0_dp, this%evects_nao_nao(2), & ! bo1v
                         a_first_col=lr1)
      ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,1]*(B[1:nbf,lr2])^T + 1.0*C[1:nbf,1:nbf]
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                        -1.0_dp, this%work2_nao_1, &     ! tmp
                                 gs_mos(1)%mos_full, &   ! MOs alpha
                         1.0_dp, this%evects_nao_nao(5), & ! o21v
                         b_first_col=lr2)

      ! work[1:noccb, 1] => X(1:nbf, 1:nbf) [1:noccb, lr1] (lr1=noccb+1)
      ! or
      ! work[1:noccb, 1] => X(1:nocca, 1:nvirb) [1:noccb, 1(lr1)]
      CALL cp_fm_set_all(this%work_nocb_1, 0.0_dp)
      DO i = 1, noccb
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), i, 1, val)
         CALL cp_fm_set_element(this%work_nocb_1, i, 1, val)
      END DO
      ! C[1:nbf,1] = 1.0*A[1:nbf,1:noccb]*B[1:noccb,1] + 0.0*C[1:nbf,1]
      CALL parallel_gemm('N', 'N', nbf, 1, noccb, &
                         1.0_dp, gs_mos(2)%mos_occ, &    ! MOs alpha (=beta)
                                 this%work_nocb_1, &     ! X[1:noccb, 1(lr1)]
                         0.0_dp, this%work2_nao_1)
      ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,1]*(B[1:nbf,lr1])^T + 1.0*C[1:nbf,1:nbf]
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                         1.0_dp, this%work2_nao_1, &     ! tmp
                                 gs_mos(2)%mos_full, &   ! MOs beta
                         1.0_dp, this%evects_nao_nao(3), & ! bco1
                         b_first_col=lr1)
      ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,lr2]*(B[1:nbf,1])^T + 1.0*C[1:nbf,1:nbf]
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                         1.0_dp, gs_mos(2)%mos_full, &   ! MOs beta
                                 this%work2_nao_1, &     ! tmp
                         1.0_dp, this%evects_nao_nao(6), & ! co12
                         a_first_col=lr2)

      ! work[1:noccb, 1] => X(1:nbf, 1:nbf) [1:noccb, lr2] (lr2=noccb+2)
      ! or
      ! work[1:noccb, 1] => X(1:nocca, 1:nvirb) [1:noccb, 2(lr2)]
      CALL cp_fm_set_all(this%work_nocb_1, 0.0_dp)
      DO i = 1, noccb
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), i, 2, val)
         CALL cp_fm_set_element(this%work_nocb_1, i, 1, val)
      END DO
      ! C[1:nbf,1] = 1.0*A[1:nbf,1:noccb]*B[1:noccb,1] + 0.0*C[1:nbf,1]
      CALL parallel_gemm('N', 'N', nbf, 1, noccb, &
                         1.0_dp, gs_mos(2)%mos_occ, &    ! MOs alpha (=beta)
                                 this%work_nocb_1, &     ! X[1:noccb, 2(lr2)]
                         0.0_dp, this%work2_nao_1)
      ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,1]*(B[1:nbf,lr2])^T + 1.0*C[1:nbf,1:nbf]
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                         1.0_dp, this%work2_nao_1, &     ! tmp
                                 gs_mos(2)%mos_full, &   ! MOs beta
                         1.0_dp, this%evects_nao_nao(4), & ! bco2
                         b_first_col=lr2)
      ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,lr1]*(B[1:nbf,1])^T + 1.0*C[1:nbf,1:nbf]
      CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                        -1.0_dp, gs_mos(2)%mos_full, &   ! MOs beta
                                 this%work2_nao_1, &     ! tmp
                         1.0_dp, this%evects_nao_nao(6), & ! co12
                         a_first_col=lr1)

      ! Sum first four components into combined matrix
      DO i = 1, 4 ! ball += bo2v + bo1v + bco1 + bco2
         CALL cp_fm_scale_and_add( &
            alpha=1.0_dp, matrix_a=this%evects_nao_nao(7), &
            beta=1.0_dp, matrix_b=this%evects_nao_nao(i))
      END DO

      ! work[1:noccb, 1:nvira] => X(1:nbf, 1:nbf) [1:noccb, nocca+1:nbf]
      ! or
      ! work[1:noccb, 1:nvira] => X(1:nocca, 1:nvirb) [1:noccb, shift+1:nvirb]
      DO i = 1, nvira
         DO j = 1, noccb
            CALL cp_fm_get_element(tddft_env%evects_mo(ivec), j, shift + i, val)
            CALL cp_fm_set_element(this%work_nocb_vira, j, i, val)
         END DO
      END DO
      ! C[1:nbf,1:noccb] = 1.0*A[1:nbf,1:nvira]*(B[1:noccb,1:nvira])^T + 0.0*C[1:nbf,1:noccb]
      CALL parallel_gemm('N', 'T', nbf, noccb, nvira, &
                         1.0_dp, gs_mos(1)%mos_virt, &   ! MOs beta (=alpha)
                                 this%work_nocb_vira, &  ! X[1:noccb,shift+1:nvirb]
                         0.0_dp, this%work_nao_nocb)     ! tmp
      ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,1:noccb]*(B[1:nbf,1:noccb])^T + 1.0*C[1:nbf,1:nbf]
      CALL parallel_gemm('N', 'T', nbf, nbf, noccb, &
                         1.0_dp, gs_mos(2)%mos_occ, &    ! MOs alpha (=beta)
                                 this%work_nao_nocb, &   ! tmp
                         1.0_dp, this%evects_nao_nao(7)) ! ball

      ! Apply theory-specific modifications for singlet/triplet states
      SELECT CASE (tddft_env%response_type)
      CASE (2) ! MRSF-S: Singlet states
         ! val = X(1:nbf,1:nbf)[lr2, lr1] or
         ! val = X(1:nocca,1:nvirb)[lr2, 1(lr1)]
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), lr2, 1, val)
         ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,lr2]*(B[1:nbf,lr1])^T + 1.0*C[1:nbf,1:nbf]
         CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                               val, gs_mos(1)%mos_full, &   ! MOs alpha
                                    gs_mos(2)%mos_full, &   ! MOs beta
                            1.0_dp, this%evects_nao_nao(7), & ! ball
                            a_first_col=lr2, &
                            b_first_col=lr1)
         ! val = X(1:nbf,1:nbf)[lr1, lr2] or
         ! val = X(1:nocca,1:nvirb)[lr1, 2(lr2)]
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), lr1, 2, val)
         ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,lr2]*(B[1:nbf,lr1])^T + 1.0*C[1:nbf,1:nbf]
         CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                               val, gs_mos(1)%mos_full, &   ! MOs alpha
                                    gs_mos(2)%mos_full, &   ! MOs beta
                            1.0_dp, this%evects_nao_nao(7), & ! ball
                            a_first_col=lr1, &
                            b_first_col=lr2)
         ! val = X(1:nbf,1:nbf)[lr1, lr1]/sqrt(2.0) or
         ! val = X(1:nocca,1:nvirb)[lr1, 1(lr1)]/sqrt(2.0)
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), lr1, 1, val)
         val = val/SQRT(2.0_dp)
         ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,lr1]*(B[1:nbf,lr1])^T + 1.0*C[1:nbf,1:nbf]
         CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                               val, gs_mos(1)%mos_full, &   ! MOs alpha
                                    gs_mos(2)%mos_full, &   ! MOs beta
                            1.0_dp, this%evects_nao_nao(7), & ! ball
                            a_first_col=lr1, &
                            b_first_col=lr1)
         ! val = -X(1:nbf,1:nbf)[lr1, lr1]/sqrt(2.0) or
         ! val = -X(1:nocca,1:nvirb)[lr1, 1(lr1)]/sqrt(2.0)
         val = -val
         ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,lr2]*(B[1:nbf,lr2])^T + 1.0*C[1:nbf,1:nbf]
         CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                               val, gs_mos(1)%mos_full, &
                                    gs_mos(2)%mos_full, &
                            1.0_dp, this%evects_nao_nao(7), &
                            a_first_col=lr2, &
                            b_first_col=lr2)

      CASE (3) ! MRSF-T: Triplet states
         ! val = X(1:nbf,1:nbf)[lr1, lr1]/sqrt(2.0) or
         ! val = X(1:nocca,1:nvirb)[lr1, 1(lr1)]/sqrt(2.0)
         CALL cp_fm_get_element(tddft_env%evects_mo(ivec), lr1, 1, val)
         val = val/SQRT(2.0_dp)
         ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,lr1]*(B[1:nbf,lr1])^T + 1.0*C[1:nbf,1:nbf]
         CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                               val, gs_mos(1)%mos_full, &
                                    gs_mos(2)%mos_full, &
                            1.0_dp, this%evects_nao_nao(7), &
                            a_first_col=lr1, &
                            b_first_col=lr1)
         ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,lr2]*(B[1:nbf,lr2])^T + 1.0*C[1:nbf,1:nbf]
         CALL parallel_gemm('N', 'T', nbf, nbf, 1, &
                               val, gs_mos(1)%mos_full, &
                                    gs_mos(2)%mos_full, &
                            1.0_dp, this%evects_nao_nao(7), &
                            a_first_col=lr2, &
                            b_first_col=lr2)
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE mrsf_mo_to_ao

! **************************************************************************************************
!> \brief Apply MRSF-TDDFT operators
!> \param[inout] this            MRSF theory object
!> \param[in]    qs_env          QuickStep environment
!> \param[inout] tddft_env       TDDFT environment
!> \details Applies HFX and Coulomb operators with theory-specific coupling parameters for
!>          different excitation types. Each coupling type uses different scaling factors.
!> \note  The HFX kernel is applied to both symmetric and antisymmetric
!>        components of the density matrix. When ADMM is enabled, the
!>        expensive four-center integrals are evaluated in a smaller
!>        auxiliary basis and then transformed back to the primary basis.
! **************************************************************************************************
   SUBROUTINE mrsf_apply(this, qs_env, tddft_env)
      CLASS(tddft_mrsf_type), INTENT(INOUT)              :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mrsf_apply'

      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_ia_g
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: hfx_data
      TYPE(section_vals_type), POINTER                   :: hf_section
      TYPE(task_list_type), POINTER                      :: task_list_orb
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(admm_type), POINTER                           :: admm_env
      REAL(KIND=dp)                                      :: pair_energy
      LOGICAL                                            :: skip_sym, skip_asym
      INTEGER                                            :: handle, imat, n_rep_hf, i_rep, i_thread

      CALL timeset(routineN, handle)

      !-----------------------------------------------
      ! ERI contraction: Sigma = ERI * X [DBCSR type]
      !-----------------------------------------------
      IF(tddft_env%do_hfx)THEN

         CALL get_qs_env(qs_env=qs_env, &
                         x_data=hfx_data, &
                         pw_env=pw_env, &
                         task_list=task_list_orb)

         IF (tddft_env%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
         END IF

         hf_section => section_vals_get_subs_vals(section_vals=qs_env%input, &
                                                  subsection_name="DFT%XC%HF")
         CALL section_vals_get(hf_section, n_repetition=n_rep_hf)

         ! Loop over 6 different coupling types and the total one
         DO imat = 1, 7
            ! Check if coupling parameter is non-zero
            IF (tddft_env%do_hfx_cam) THEN
               IF (ABS(this%spc(imat, 1)) < threshold .AND. &
                   ABS(this%spc(imat, 2)) < threshold) CYCLE
            ELSE
               IF (ABS(this%spc(imat, 1)) < threshold) CYCLE
            END IF

            ! Set HFX scaling parameters for this coupling type
            IF (tddft_env%do_hfx_cam) THEN
               DO i_rep = 1, n_rep_hf
                  DO i_thread = 1, SIZE(hfx_data, 2)
                     hfx_data(i_rep, i_thread)%potential_parameter%scale_coulomb = this%spc(imat, 1)
                     hfx_data(i_rep, i_thread)%potential_parameter%scale_longrange = this%spc(imat, 2)
                  END DO
               END DO
            ELSE
               DO i_rep = 1, n_rep_hf
                  DO i_thread = 1, SIZE(hfx_data, 2)
                     hfx_data(i_rep, i_thread)%general_parameter%fraction = this%spc(imat, 1)
                  END DO
               END DO
            END IF

            !--------------------------------------
            ! Get Symmetric part of X vector
            !--------------------------------------
            CALL cp_fm_set_all(matrix=this%work1_nao_nao, alpha=0.0_dp)
            ! work1 = 0.5*X
            CALL cp_fm_scale_and_add( &
               alpha=0.0_dp, matrix_a=this%work1_nao_nao, &
               beta=0.5_dp, matrix_b=this%evects_nao_nao(imat))
            ! work2 = X^T
            CALL cp_fm_transpose(this%evects_nao_nao(imat), this%work2_nao_nao)

            ! work1 += 0.5*work2
            CALL cp_fm_scale_and_add( &
               alpha=1.0_dp, matrix_a=this%work1_nao_nao, &
               beta=0.5_dp, matrix_b=this%work2_nao_nao)

            skip_sym = is_zero_matrix_fm(matrix=this%work1_nao_nao)

            IF (tddft_env%do_admm) THEN
               ! Project onto the auxiliary basis
               ! C[nbf_aux,nbf] = 1.0*A[nbf_aux,nbf]*B[nbf,nbf] + 0.0*C[nbf_aux,nbf]
               CALL parallel_gemm('N', 'N', tddft_env%nbf_aux, tddft_env%nbf, tddft_env%nbf, &
                                  1.0_dp, admm_env%A, &
                                          this%work1_nao_nao, &
                                  0.0_dp, this%work_aux_nao)
               ! C[nbf_aux,nbf_aux] = 1.0*A[nbf_aux,nbf]*(B[nbf_aux,nbf])^T + 0.0*C[nbf_aux,nbf_aux]
               CALL parallel_gemm('N', 'T', tddft_env%nbf_aux, tddft_env%nbf_aux, tddft_env%nbf, &
                                  1.0_dp, this%work_aux_nao, &
                                          admm_env%A, &
                                  0.0_dp, this%work_aux_aux)

               ! Set density-like to 0.0
               CALL dbcsr_set(matrix=this%work_rho_aux_aux_symm(1)%matrix, &
                              alpha=0.0_dp)

               ! Convert aux FM to DBCSR
               CALL copy_fm_to_dbcsr(fm=this%work_aux_aux, &
                                     matrix=this%work_rho_aux_aux_symm(1)%matrix, &
                                     keep_sparsity=.TRUE.)

               ! Symmetric part
               CALL dbcsr_set(matrix=this%work_hmat_aux_symm(1)%matrix, alpha=0.0_dp)
               IF (.NOT. skip_sym) THEN
                  CALL tddft_hfx_matrix(matrix_ks=this%work_hmat_aux_symm, &
                                        rho_ao=this%work_rho_aux_aux_symm, &
                                        qs_env=qs_env)
               END IF
            ELSE
               ! Set density-like to 0.0
               CALL dbcsr_set(matrix=this%work_rho_nao_nao_symm(1)%matrix, &
                              alpha=0.0_dp)

               ! Convert FM to DBCSR
               CALL copy_fm_to_dbcsr(fm=this%work1_nao_nao, &
                                     matrix=this%work_rho_nao_nao_symm(1)%matrix, &
                                     keep_sparsity=.TRUE.)

               ! Apply HFX to symmetric part
               CALL dbcsr_set(matrix=this%work_hmat_symm(1)%matrix, alpha=0.0_dp)
               IF (.NOT. skip_sym) THEN
                  CALL tddft_hfx_matrix(matrix_ks=this%work_hmat_symm, &
                                        rho_ao=this%work_rho_nao_nao_symm, &
                                        qs_env=qs_env)
               END IF
            END IF

            !--------------------------------------
            ! Get Antisymmetric part of X vector
            !--------------------------------------
            CALL cp_fm_set_all(matrix=this%work1_nao_nao, alpha=0.0_dp)

            ! work1 = 0.5*X
            CALL cp_fm_scale_and_add( &
               alpha=0.0_dp, matrix_a=this%work1_nao_nao, &
               beta=0.5_dp, matrix_b=this%evects_nao_nao(imat))
            ! work1 -= 0.5*work2
            CALL cp_fm_scale_and_add( &
               alpha=1.0_dp, matrix_a=this%work1_nao_nao, &
               beta=-0.5_dp, matrix_b=this%work2_nao_nao)

            skip_asym = is_zero_matrix_fm(matrix=this%work1_nao_nao)

            IF (tddft_env%do_admm) THEN
               ! Project onto the auxiliary basis
               ! C[nbf_aux,nbf] = 1.0*A[nbf_aux,nbf]*B[nbf,nbf] + 0.0*C[nbf_aux,nbf]
               CALL parallel_gemm('N', 'N', tddft_env%nbf_aux, tddft_env%nbf, tddft_env%nbf, &
                                  1.0_dp, admm_env%A, &
                                          this%work1_nao_nao, &
                                  0.0_dp, this%work_aux_nao)
               ! C[nbf_aux,nbf_aux] = 1.0*A[nbf_aux,nbf]*(B[nbf_aux,nbf])^T + 0.0*C[nbf_aux,nbf_aux]
               CALL parallel_gemm('N', 'T', tddft_env%nbf_aux, tddft_env%nbf_aux, tddft_env%nbf, &
                                  1.0_dp, this%work_aux_nao, &
                                          admm_env%A, &
                                  0.0_dp, this%work_aux_aux)

               ! Set density-like to 0.0
               CALL dbcsr_set(matrix=this%work_rho_aux_aux_asymm(1)%matrix, &
                              alpha=0.0_dp)

               ! Convert aux FM to DBCSR
               CALL copy_fm_to_dbcsr(fm=this%work_aux_aux, &
                                     matrix=this%work_rho_aux_aux_asymm(1)%matrix, &
                                     keep_sparsity=.TRUE.)

               ! Apply HFX to antisymmetric part
               CALL dbcsr_set(matrix=this%work_hmat_aux_asymm(1)%matrix, alpha=0.0_dp)
               IF (.NOT. skip_asym) THEN
                  CALL tddft_hfx_matrix(matrix_ks=this%work_hmat_aux_asymm, &
                                        rho_ao=this%work_rho_aux_aux_asymm, &
                                        qs_env=qs_env)
               END IF
            ELSE
               ! Set density-like to 0.0
               CALL dbcsr_set(matrix=this%work_rho_nao_nao_asymm(1)%matrix, &
                              alpha=0.0_dp)

               ! Convert FM to DBCSR
               CALL copy_fm_to_dbcsr(fm=this%work1_nao_nao, &
                                     matrix=this%work_rho_nao_nao_asymm(1)%matrix, &
                                     keep_sparsity=.TRUE.)

               ! Apply HFX to antisymmetric part
               CALL dbcsr_set(matrix=this%work_hmat_asymm(1)%matrix, alpha=0.0_dp)
               IF (.NOT. skip_asym) THEN
                  CALL tddft_hfx_matrix(matrix_ks=this%work_hmat_asymm, &
                                        rho_ao=this%work_rho_nao_nao_asymm, &
                                        qs_env=qs_env)
               END IF
            END IF

            ! Convert DBCSR to FM
            CALL cp_fm_set_all(matrix=this%work1_nao_nao, alpha=0.0_dp)
            IF (tddft_env%do_admm) THEN
               ! Convert DBCSR to aux FM
               CALL cp_fm_set_all(matrix=this%work_aux_aux, alpha=0.0_dp)
               CALL copy_dbcsr_to_fm(matrix=this%work_hmat_aux_symm(1)%matrix, &
                                     fm=this%work_aux_aux)
               ! Project back from the auxiliary basis
               ! C[nbf,nbf_aux] = 1.0*(A[nbf_aux,nbf])^T*B[nbf_aux,nbf_aux] + 0.0*C[nbf,nbf_aux]
               CALL parallel_gemm('T', 'N', tddft_env%nbf, tddft_env%nbf_aux, tddft_env%nbf_aux, &
                                  1.0_dp, admm_env%A, &
                                          this%work_aux_aux, &
                                  0.0_dp, this%work_nao_aux)
               ! C[nbf,nbf] = 1.0*A[nbf,nbf_aux]*B[nbf_aux,nbf] + 0.0*C[nbf,nbf]
               CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf_aux, &
                                  1.0_dp, this%work_nao_aux, &
                                          admm_env%A, &
                                  0.0_dp, this%work1_nao_nao)
            ELSE
               ! Convert DBCSR to FM
               CALL copy_dbcsr_to_fm(matrix=this%work_hmat_symm(1)%matrix, &
                                     fm=this%work1_nao_nao)
            END IF

            CALL cp_fm_set_all(matrix=this%work2_nao_nao, alpha=0.0_dp)
            IF (tddft_env%do_admm) THEN
               ! Convert DBCSR to aux FM
               CALL cp_fm_set_all(matrix=this%work_aux_aux, alpha=0.0_dp)
               CALL copy_dbcsr_to_fm(matrix=this%work_hmat_aux_asymm(1)%matrix, &
                                     fm=this%work_aux_aux)
               ! Project back from the auxiliary basis
               ! C[nbf,nbf_aux] = 1.0*(A[nbf_aux,nbf])^T*B[nbf_aux,nbf_aux] + 0.0*C[nbf,nbf_aux]
               CALL parallel_gemm('T', 'N', tddft_env%nbf, tddft_env%nbf_aux, tddft_env%nbf_aux, &
                                  1.0_dp, admm_env%A, &
                                          this%work_aux_aux, &
                                  0.0_dp, this%work_nao_aux)
               ! C[nbf,nbf] = 1.0*A[nbf,nbf_aux]*B[nbf_aux,nbf] + 0.0*C[nbf,nbf]
               CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf_aux, &
                                  1.0_dp, this%work_nao_aux, &
                                          admm_env%A, &
                                  0.0_dp, this%work2_nao_nao)
            ELSE
               ! Convert DBCSR to FM
               CALL copy_dbcsr_to_fm(matrix=this%work_hmat_asymm(1)%matrix, &
                                    fm=this%work2_nao_nao)
            END IF

            CALL cp_fm_set_all(matrix=this%hmat_nao_nao(imat), alpha=0.0_dp)
            ! A = 0.0*A + 1.0*B
            CALL cp_fm_scale_and_add( &
               alpha=0.0_dp, matrix_a=this%hmat_nao_nao(imat), &
               beta=1.0_dp, matrix_b=this%work1_nao_nao)
            ! A = 1.0*A + 1.0*B
            CALL cp_fm_scale_and_add( &
               alpha=1.0_dp, matrix_a=this%hmat_nao_nao(imat), &
               beta=1.0_dp, matrix_b=this%work2_nao_nao)
            CALL cp_fm_set_all(matrix=this%work2_nao_nao, alpha=0.0_dp)

            IF (imat > 4 .OR. skip_sym) CYCLE

            !---------------
            ! Coulomb part
            !---------------

!           ! Update density on the grid
            CALL qs_rho_update_rho(rho_struct=this%rho_nao_struct, &
                                   qs_env=qs_env, &
                                   pw_env_external=pw_env, &
                                   task_list_external=task_list_orb)

            ! Get the updated grid density
            CALL qs_rho_get(this%rho_nao_struct, &
                            rho_g=rho_ia_g)

            ! Get Poisson environment
            CALL pw_env_get(pw_env, &
                            poisson_env=poisson_env)

            ! Compute Coulomb potential
            CALL pw_poisson_solve(poisson_env, &
                                  rho_ia_g(1), &
                                  pair_energy, &
                                  this%work_v_gspace)
            ! Transfer potential to real space
            CALL pw_transfer(this%work_v_gspace, this%work_v_rspace)

            ! Scale by Hartree-Fock exchange
            CALL pw_scale(this%work_v_rspace, this%spc(imat, 1))

            ! Scale by grid volume element
            CALL pw_scale(this%work_v_rspace, this%work_v_rspace%pw_grid%dvol)

            ! Integrate into AO matrix
            CALL dbcsr_set(this%hmat_coulomb(1)%matrix, 0.0_dp)
            CALL integrate_v_rspace(v_rspace=this%work_v_rspace, &
                                    hmat=this%hmat_coulomb(1), &
                                    qs_env=qs_env, &
                                    calculate_forces=.FALSE., &
                                    gapw=.FALSE., &
                                    pw_env_external=pw_env, &
                                    task_list_external=task_list_orb)

            ! Convert to full matrix and accumulate into hmat_nao_nao
            CALL copy_dbcsr_to_fm(matrix=this%hmat_coulomb(1)%matrix, &
                                  fm=this%work2_nao_nao)
            CALL cp_fm_scale_and_add(alpha=1.0_dp, &
                                     matrix_a=this%hmat_nao_nao(imat), &
                                     beta=1.0_dp, &
                                     matrix_b=this%work2_nao_nao)

         END DO

         IF (tddft_env%do_hfx_cam) THEN
            DO i_rep = 1, n_rep_hf
               DO i_thread = 1, SIZE(hfx_data, 2)
                  hfx_data(i_rep, i_thread)%potential_parameter%scale_coulomb = tddft_env%C_hf_alpha
                  hfx_data(i_rep, i_thread)%potential_parameter%scale_longrange = tddft_env%C_hf_beta
               END DO
            END DO
         ELSE
            DO i_rep = 1, n_rep_hf
               DO i_thread = 1, SIZE(hfx_data, 2)
                  hfx_data(i_rep, i_thread)%general_parameter%fraction = tddft_env%C_hf_alpha
               END DO
            END DO
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE mrsf_apply

! **************************************************************************************************
!> \brief ...
!> \param[inout] this ...
!> \param[inout] tddft_env ...
!> \param[in] gs_mos ...
!> \param[in] ivec ...
! **************************************************************************************************
   SUBROUTINE mrsf_ao_to_mo(this, tddft_env, gs_mos, ivec)
      CLASS(tddft_mrsf_type), INTENT(INOUT)              :: this
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER, INTENT(IN)                             :: gs_mos
      INTEGER, INTENT(IN)                                :: ivec

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mrsf_ao_to_mo'

      REAL(KIND=dp)                                      :: scalar, val, xlr, val2
      INTEGER                                            :: handle, i, &
         lr1, lr2, nbf, nocca, noccb, nvira, nvirb, shift

      CALL timeset(routineN, handle)

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

      ! Assign pointers from tddft_env%evects_nao_nao [nbf, nbf]
      ! agdlr => this%hmat_nao_nao(7)
      ! ado2v => this%hmat_nao_nao(1)
      ! ado1v => this%hmat_nao_nao(2)
      ! adco1 => this%hmat_nao_nao(3)
      ! adco2 => this%hmat_nao_nao(4)
      ! ao21v => this%hmat_nao_nao(5)
      ! aco12 => this%hmat_nao_nao(6)

      !-----------------------------------------------
      ! Transform back to MO basis
      !-----------------------------------------------
      ! C[1:nocca,1:nbf] = 1.0*(A[1:nbf,1:nocca])^T*B[1:nbf,1:nbf] + 0.0*C[1:nocca,1:nbf]
      CALL parallel_gemm('T', 'N', nocca, nbf, nbf, &
                         1.0_dp, gs_mos(1)%mos_occ, &    ! MOs alpha
                                 this%hmat_nao_nao(7), & ! agdlr
                         0.0_dp, this%work_noc_nao)      ! tmp

      ! C[1:nocca,1:nvirb] = 1.0*A[1:noca,1:nbf]*B[1:nbf,1:nvirb] + 0.0*C[1:nocca,1:nvirb]
      CALL parallel_gemm('N', 'N', nocca, nvirb, nbf, &
                         1.0_dp, this%work_noc_nao, &
                                 gs_mos(2)%mos_virt, &
                         0.0_dp, tddft_env%sigma_vecs(ivec))

      ! scalar = X(1:nbf,1:nbf)[lr1, lr1] or
      ! scalar = X(1:nocca,1:nvirb)[lr1, 1(lr1)]
      CALL cp_fm_get_element(tddft_env%sigma_vecs(ivec), lr1, 1, scalar)
      ! val = X(1:nbf,1:nbf)[lr2, lr2] or
      ! val = X(1:nocca,1:nvirb)[lr2, 2(lr2)]
      CALL cp_fm_get_element(tddft_env%sigma_vecs(ivec), lr2, 2, val)
      SELECT CASE (tddft_env%response_type)
      CASE (2) ! MRSF-S
         scalar = (scalar - val)*sqrt12
      CASE (3) ! MRSF-T
         scalar = (scalar + val)*sqrt12
      END SELECT
      ! X(1:nbf,1:nbf)[lr1, lr1] = scalar
      ! X(1:nocca,1:nvirb)[lr1, 1(lr1)] = scalar
      CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr1, 1, scalar)

      ! C[1:nbf,1] = 1.0*A[1:nbf,1:nbf]*B[1:nbf,lr2] + 0.0*C[1:nbf,1]
      CALL parallel_gemm('N', 'N', nbf, 1, nbf, &
                         1.0_dp, this%hmat_nao_nao(2), & ! ado1v
                                 gs_mos(2)%mos_full, &   ! MOs beta
                         0.0_dp, this%work1_nao_1, &     ! tmp
                         b_first_col=lr2)
      ! C[1:nbf,1] = 1.0*A[1:nbf,1:nbf]*B[1:nbf,lr1] + 1.0*C[1:nbf,1]
      CALL parallel_gemm('N', 'N', nbf, 1, nbf, &
                         1.0_dp, this%hmat_nao_nao(6), & ! aco12
                                 gs_mos(2)%mos_full, &   ! MOs beta
                         1.0_dp, this%work1_nao_1, &     ! tmp
                         b_first_col=lr1)
      ! C[1:noccb,1] = 1.0*(A[1:nbf,1:noccb])^T*B[1:nbf,1] + 1.0*C[1:noccb,1]
      CALL parallel_gemm('T', 'N', noccb, 1, nbf, &
                         1.0_dp, gs_mos(2)%mos_occ, &    ! MOs alpha (=beta)
                                 this%work1_nao_1, &     ! tmp
                         0.0_dp, this%work_nocb_1)       ! tmp2[1:noccb,1]
      DO i = 1, noccb
         ! scalar = work(1:nocca,1) [1:noccb, 1]
         CALL cp_fm_get_element(this%work_nocb_1, i, 1, scalar)
         CALL cp_fm_get_element(tddft_env%sigma_vecs(ivec), i, 2, val)
         val = val + scalar
         ! X(1:nbf,1:nbf)[1:noccb, lr2] = val OR
         ! X(1:nocca,1:nvirb)[1:noccb, 2(lr2)] = val
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), i, 2, val)
      END DO

      ! C[1:nbf,1] = 1.0*A[1:nbf,1:nbf]*B[1:nbf,lr1] + 0.0*C[1:nbf,1]
      CALL parallel_gemm('N', 'N', nbf, 1, nbf, &
                         1.0_dp, this%hmat_nao_nao(1), & ! ado2v
                                 gs_mos(2)%mos_full, &   ! MOs beta
                         0.0_dp, this%work1_nao_1, &     ! tmp
                         b_first_col=lr1)
      ! C[1:nbf,1] = -1.0*A[1:nbf,1:nbf]*B[1:nbf,lr2] + 1.0*C[1:nbf,1]
      CALL parallel_gemm('N', 'N', nbf, 1, nbf, &
                        -1.0_dp, this%hmat_nao_nao(6), & ! aco12
                                 gs_mos(2)%mos_full, &   ! MOs beta
                         1.0_dp, this%work1_nao_1, &     ! tmp
                         b_first_col=lr2)
      ! C[1:noccb,1] = 1.0*(A[1:nbf,1:noccb])^T*B[1:nbf,1] + 1.0*C[1:noccb,1]
      CALL parallel_gemm('T', 'N', noccb, 1, nbf, &
                         1.0_dp, gs_mos(2)%mos_occ, &    ! MOs alpha (=beta)
                                 this%work1_nao_1, &     ! tmp
                         0.0_dp, this%work_nocb_1)       ! tmp2[1:noccb,1]
      DO i = 1, noccb
         ! scalar = work(1:nocca,1) [1:noccb, 1]
         CALL cp_fm_get_element(this%work_nocb_1, i, 1, scalar)
         CALL cp_fm_get_element(tddft_env%sigma_vecs(ivec), i, 1, val)
         val = val + scalar
         ! X(1:nbf,1:nbf)[1:noccb, lr1] = val OR
         ! X(1:nocca,1:nvirb)[1:noccb, 1(lr1)] = val
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), i, 1, val)
      END DO

      ! C[1:nbf,1] = 1.0*(A[1:nbf,1:nbf])^T*B[1:nbf,lr1] + 0.0*C[1:nbf,1]
      CALL parallel_gemm('T', 'N', nbf, 1, nbf, &
                         1.0_dp, this%hmat_nao_nao(4), & ! adco2
                                 gs_mos(1)%mos_full, &   ! MOs alpha
                         0.0_dp, this%work1_nao_1, &     ! tmp
                         b_first_col=lr1)
      ! C[1:nbf,1] = 1.0*(A[1:nbf,1:nbf])^T*B[1:nbf,lr2] + 1.0*C[1:nbf,1]
      CALL parallel_gemm('T', 'N', nbf, 1, nbf, &
                         1.0_dp, this%hmat_nao_nao(5), & ! ao21v
                                 gs_mos(1)%mos_full, &   ! MOs alpha
                         1.0_dp, this%work1_nao_1, &     ! tmp
                         b_first_col=lr2)
      ! C[1:nvira,1] = 1.0*(A[1:nbf,1:vira])^T*B[1:nbf,1] + 1.0*C[1:vira,1]
      CALL parallel_gemm('T', 'N', nvira, 1, nbf, &
                         1.0_dp, gs_mos(1)%mos_virt, &   ! MOs beta (=alpha)
                                 this%work1_nao_1, &     ! tmp
                         0.0_dp, this%work_vira_1)       ! tmp2
      DO i = 1, nvira
         ! scalar = work(1:nvira,1) [1:nvira, 1]
         CALL cp_fm_get_element(this%work_vira_1, i, 1, scalar)
         CALL cp_fm_get_element(tddft_env%sigma_vecs(ivec), lr1, shift + i, val)
         val = val + scalar
         ! X(1:nbf,1:nbf)[lr1, noca+1:nbf] = val OR
         ! X(1:nocca,1:nvirb)[lr1, shift+1:nvirb] = val
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr1, shift + i, val)
      END DO

      ! C[1:nbf,1] = 1.0*(A[1:nbf,1:nbf])^T*B[1:nbf,lr2] + 0.0*C[1:nbf,1]
      CALL parallel_gemm('T', 'N', nbf, 1, nbf, &
                         1.0_dp, this%hmat_nao_nao(3), & ! adco1
                                 gs_mos(1)%mos_full, &   ! MOs alpha
                         0.0_dp, this%work1_nao_1, &     ! tmp
                         b_first_col=lr2)
      ! C[1:nbf,1] = -1.0*(A[1:nbf,1:nbf])^T*B[1:nbf,lr1] + 1.0*C[1:nbf,1]
      CALL parallel_gemm('T', 'N', nbf, 1, nbf, &
                        -1.0_dp, this%hmat_nao_nao(5), & ! ao21v
                                 gs_mos(1)%mos_full, &   ! va
                         1.0_dp, this%work1_nao_1, &     ! tmp
                         b_first_col=lr1)
      ! C[1:nvira,1] = 1.0*(A[1:nbf,1:vira])^T*B[1:nbf,1] + 1.0*C[1:vira,1]
      CALL parallel_gemm('T', 'N', nvira, 1, nbf, &
                         1.0_dp, gs_mos(1)%mos_virt, &   ! MOs beta (=alpha)
                                 this%work1_nao_1, &     ! tmp [nbf, 1]
                         0.0_dp, this%work_vira_1)       ! tmp2 [nvira, 1]

      DO i = 1, nvira
         ! scalar = work(1:nvira,1) [1:nvira, 1]
         CALL cp_fm_get_element(this%work_vira_1, i, 1, scalar)
         CALL cp_fm_get_element(tddft_env%sigma_vecs(ivec), lr2, shift + i, val)
         val = val + scalar
         ! X(1:nbf,1:nbf)[lr2, noca+1:nbf] = val OR
         ! X(1:nocca,1:nvirb)[lr2, shift+1:nvirb] = val
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr2, shift + i, val)
      END DO

      SELECT CASE (tddft_env%response_type)
      CASE (2) ! MRSF-S
         ! X(1:nbf,1:nbf)[lr2, lr2] = 0.0 or
         ! X(1:nocca,1:nvirb)[lr2, 2(lr2)] = 0.0
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr2, 2, 0.0_dp)
      CASE (3) ! MRSF-T
         ! X(1:nbf,1:nbf)[lr2, lr1] = 0.0 or
         ! X(1:nocca,1:nvirb)[lr2, 1(lr1)] = 0.0
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr2, 1, 0.0_dp)
         ! X(1:nbf,1:nbf)[lr1, lr2] = 0.0 or
         ! X(1:nocca,1:nvirb)[lr1, 2(lr2)] = 0.0
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr1, 2, 0.0_dp)
         ! X(1:nbf,1:nbf)[lr2, lr2] = 0.0 or
         ! X(1:nocca,1:nvirb)[lr2, 2(lr2)] = 0.0
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr2, 2, 0.0_dp)
      END SELECT

      CALL cp_fm_to_fm(source=tddft_env%evects_mo(ivec), &
                       destination=this%work2_noc_vir)
      CALL cp_fm_set_element(this%work2_noc_vir, lr1, 1, 0.0_dp)

      ! C[nocca,nvirb] = -1.0*A[nocca,nocca]*B[nocca,nvirb]+0.0*C[nocca,nvirb]
      CALL parallel_gemm('N', 'N', nocca, nvirb, nocca, &
                         -1.0_dp, this%focka_noc_noc, &
                                  this%work2_noc_vir, &
                          0.0_dp, this%work_noc_vir)

      ! C[nocca,nvirb] = 1.0*A[nocca,nvirb]*B[nvirb,nvirb] + 1.0*C[nocca,nvirb]
      CALL parallel_gemm('N', 'N', nocca, nvirb, nvirb, &
                         1.0_dp, this%work2_noc_vir, &
                         this%fockb_vir_vir, &
                         1.0_dp, this%work_noc_vir)

      CALL cp_fm_get_element(tddft_env%evects_mo(ivec), lr1, 1, xlr)

      SELECT CASE (tddft_env%response_type)
      CASE (2) ! MRSF-S
         DO i = 1, nvirb
            CALL cp_fm_get_element(this%fockb_vir_vir, i, 1, scalar)
            CALL cp_fm_get_element(this%work_noc_vir, lr1, i, val)
            val = val + scalar*xlr*sqrt12
            CALL cp_fm_set_element(this%work_noc_vir, lr1, i, val)

            CALL cp_fm_get_element(this%fockb_vir_vir, i, 2, scalar)
            CALL cp_fm_get_element(this%work_noc_vir, lr2, i, val)
            val = val - scalar*xlr*sqrt12
            CALL cp_fm_set_element(this%work_noc_vir, lr2, i, val)
         END DO

         DO i = 1, nocca
            CALL cp_fm_get_element(this%focka_noc_noc, i, lr1, scalar)
            CALL cp_fm_get_element(this%work_noc_vir, i, 1, val)
            val = val - scalar*xlr*sqrt12
            CALL cp_fm_set_element(this%work_noc_vir, i, 1, val)

            CALL cp_fm_get_element(this%focka_noc_noc, i, lr2, scalar)
            CALL cp_fm_get_element(this%work_noc_vir, i, 2, val)
            val = val + scalar*xlr*sqrt12
            CALL cp_fm_set_element(this%work_noc_vir, i, 2, val)
         END DO

         val2 = 0.0_dp
         DO i = 1, nocca
            CALL cp_fm_get_element(this%focka_noc_noc, lr1, i, scalar)
            CALL cp_fm_get_element(this%work2_noc_vir, i, 1, val)
            val2 = val2 - scalar*val
         END DO
         DO i = 1, nocca
            CALL cp_fm_get_element(this%focka_noc_noc, lr2, i, scalar)
            CALL cp_fm_get_element(this%work2_noc_vir, i, 2, val)
            val2 = val2 + scalar*val
         END DO
         DO i = 1, nvirb
            CALL cp_fm_get_element(this%fockb_vir_vir, 1, i, scalar)
            CALL cp_fm_get_element(this%work2_noc_vir, lr1, i, val)
            val2 = val2 + scalar*val
         END DO
         DO i = 1, nvirb
            CALL cp_fm_get_element(this%fockb_vir_vir, 2, i, scalar)
            CALL cp_fm_get_element(this%work2_noc_vir, lr2, i, val)
            val2 = val2 - scalar*val
         END DO
         val2 = val2*sqrt12

      CASE (3) ! MRSF-T
         DO i = 1, nvirb
            CALL cp_fm_get_element(this%fockb_vir_vir, i, 1, scalar)
            CALL cp_fm_get_element(this%work_noc_vir, lr1, i, val)
            val = val + scalar*xlr*sqrt12
            CALL cp_fm_set_element(this%work_noc_vir, lr1, i, val)

            CALL cp_fm_get_element(this%fockb_vir_vir, i, 2, scalar)
            CALL cp_fm_get_element(this%work_noc_vir, lr2, i, val)
            val = val + scalar*xlr*sqrt12
            CALL cp_fm_set_element(this%work_noc_vir, lr2, i, val)
         END DO

         DO i = 1, nocca
            CALL cp_fm_get_element(this%focka_noc_noc, i, lr1, scalar)
            CALL cp_fm_get_element(this%work_noc_vir, i, 1, val)
            val = val - scalar*xlr*sqrt12
            CALL cp_fm_set_element(this%work_noc_vir, i, 1, val)

            CALL cp_fm_get_element(this%focka_noc_noc, i, lr2, scalar)
            CALL cp_fm_get_element(this%work_noc_vir, i, 2, val)
            val = val - scalar*xlr*sqrt12
            CALL cp_fm_set_element(this%work_noc_vir, i, 2, val)
         END DO

         val2 = 0.0_dp
         DO i = 1, nocca
            CALL cp_fm_get_element(this%focka_noc_noc, lr1, i, scalar)
            CALL cp_fm_get_element(this%work2_noc_vir, i, 1, val)
            val2 = val2 - scalar*val
         END DO
         DO i = 1, nocca
            CALL cp_fm_get_element(this%focka_noc_noc, lr2, i, scalar)
            CALL cp_fm_get_element(this%work2_noc_vir, i, 2, val)
            val2 = val2 - scalar*val
         END DO
         DO i = 1, nvirb
            CALL cp_fm_get_element(this%fockb_vir_vir, 1, i, scalar)
            CALL cp_fm_get_element(this%work2_noc_vir, lr1, i, val)
            val2 = val2 + scalar*val
         END DO
         DO i = 1, nvirb
            CALL cp_fm_get_element(this%fockb_vir_vir, 2, i, scalar)
            CALL cp_fm_get_element(this%work2_noc_vir, lr2, i, val)
            val2 = val2 + scalar*val
         END DO
         val2 = val2*sqrt12

      END SELECT
      CALL cp_fm_get_element(this%fockb_vir_vir, 1, 1, val)
      val2 = val2 + val*xlr*0.5_dp
      CALL cp_fm_get_element(this%fockb_vir_vir, 2, 2, val)
      val2 = val2 + val*xlr*0.5_dp
      CALL cp_fm_get_element(this%focka_noc_noc, lr1, lr1, val)
      val2 = val2 - val*xlr*0.5_dp
      CALL cp_fm_get_element(this%focka_noc_noc, lr2, lr2, val)
      val2 = val2 - val*xlr*0.5_dp

      CALL cp_fm_set_element(this%work_noc_vir, lr1, 1, val2)

      ! Add corrections to sigma_mo
      ! A = 1.0 * A + 1.0 * B
      CALL cp_fm_scale_and_add( &
         alpha=1.0_dp, matrix_a=tddft_env%sigma_vecs(ivec), &
         beta=1.0_dp, matrix_b=this%work_noc_vir)

      SELECT CASE (tddft_env%response_type)
      CASE (2) ! MRSF-S
         ! X(1:nbf,1:nbf)[lr2, lr2] = 0.0 or
         ! X(1:nocca,1:nvirb)[lr2, 2(lr2)] = 0.0
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr2, 2, 0.0_dp)
      CASE (3) ! MRSF-T
         ! X(1:nbf,1:nbf)[lr2, lr1] = 0.0 or
         ! X(1:nocca,1:nvirb)[lr2, 1(lr1)] = 0.0
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr2, 1, 0.0_dp)
         ! X(1:nbf,1:nbf)[lr1, lr2] = 0.0 or
         ! X(1:nocca,1:nvirb)[lr1, 2(lr2)] = 0.0
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr1, 2, 0.0_dp)
         ! X(1:nbf,1:nbf)[lr2, lr2] = 0.0 or
         ! X(1:nocca,1:nvirb)[lr2, 2(lr2)] = 0.0
         CALL cp_fm_set_element(tddft_env%sigma_vecs(ivec), lr2, 2, 0.0_dp)
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE mrsf_ao_to_mo

! **************************************************************************************************
!> \brief ...
!> \param[inout] this ...
!> \param qs_env ...
!> \param[inout] tddft_env ...
!> \param[out] ostrength ...
!> \param[in] gs_mos ...
!> \param[in] tddfpt_control ...
!> \param[in] log_unit ...
! **************************************************************************************************
   SUBROUTINE mrsf_print(this, qs_env, tddft_env, ostrength, gs_mos, tddfpt_control, log_unit)
      CLASS(tddft_mrsf_type), INTENT(INOUT)              :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      REAL(kind=dp), DIMENSION(:), INTENT(OUT)           :: ostrength
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER, INTENT(IN)                             :: gs_mos
      TYPE(tddfpt2_control_type), POINTER, INTENT(IN)    :: tddfpt_control
      INTEGER, INTENT(IN)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mrsf_print'

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: tddfpt_section
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      CHARACTER(LEN=default_path_length)                 :: suffix
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: s2_mat, amplitudes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: indices, i_indices, a_indices
      INTEGER, DIMENSION(2)                              :: nmo_occ
      CHARACTER(len=2)                                   :: lsd_str
      CHARACTER(len=1)                                   :: mult_str
      REAL(kind=dp)                                      :: osc_strength, val, f, energ, gs_energy
      INTEGER                                            :: handle, ideriv, istate, jstate, ispin, &
         n_amplitudes, max_amplitudes, i, a, idx, lr1, lr2, nbf, nocca, noccb, nvira, nvirb, shift
      ! DEBUG variables
      REAL(kind=dp)                                      :: dipole_norm_x, dipole_norm_y, dipole_norm_z
      REAL(kind=dp)                                      :: dipole_max_x, dipole_max_y, dipole_max_z
      REAL(kind=dp)                                      :: sample_val, trden_norm, trden_max
      ! Transition density matrices
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: trden
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      REAL(kind=dp)                                      :: scalar, sqrt12, sqrt2
      INTEGER                                            :: ist, jst, i_lr1, a_lr1, i_lr2, a_lr2, i_g, a_g, i_d, a_d

      CALL timeset(routineN, handle)

      SELECT CASE (tddft_env%reference_type)
      CASE (1)
         lsd_str = 'RO'
      CASE (2)
         lsd_str = ' U'
      END SELECT

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

      ! Initialize constants for MRSF transition density calculation
      sqrt12 = 1.0_dp/SQRT(2.0_dp)
      sqrt2 =SQRT(2.0_dp)

      ! Compute (i,a) indices for special positions (MRSF active orbitals)
      i_lr1 = nocca - 1
      a_lr1 = nocca - noccb - 1

      i_lr2 = nocca
      a_lr2 = nocca - noccb

      i_g = nocca
      a_g = nocca - noccb - 1

      i_d = nocca - 1
      a_d = nocca - noccb

      CALL get_qs_env(qs_env=qs_env, &
                      matrix_s=matrix_s, &
                      energy=energy)
      gs_energy = energy%total

      DO ispin = 1, tddft_env%nSpins
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, &
                             ncol_global=nmo_occ(ispin))
      END DO

      ! *** DEBUG: Print dipole integral norms ***
      IF (log_unit > 0) THEN
         dipole_norm_x = cp_fm_frobenius_norm(tddft_env%dipole_ao_fm(1))
         dipole_norm_y = cp_fm_frobenius_norm(tddft_env%dipole_ao_fm(2))
         dipole_norm_z = cp_fm_frobenius_norm(tddft_env%dipole_ao_fm(3))
         CALL cp_fm_maxabsval(tddft_env%dipole_ao_fm(1), dipole_max_x)
         CALL cp_fm_maxabsval(tddft_env%dipole_ao_fm(2), dipole_max_y)
         CALL cp_fm_maxabsval(tddft_env%dipole_ao_fm(3), dipole_max_z)

         WRITE (log_unit, '(/,A)') "=== DEBUG: Dipole Integrals in AO Basis ==="
         WRITE (log_unit, '(A,3ES14.6)') "  Frobenius norms (x,y,z): ", dipole_norm_x, dipole_norm_y, dipole_norm_z
         WRITE (log_unit, '(A,3ES14.6)') "  Max abs values (x,y,z):  ", dipole_max_x, dipole_max_y, dipole_max_z

         ! Sample element (1,1) from each component
         CALL cp_fm_get_element(tddft_env%dipole_ao_fm(1), 1, 1, sample_val)
         WRITE (log_unit, '(A,ES14.6)') "  dipole_x(1,1): ", sample_val
         CALL cp_fm_get_element(tddft_env%dipole_ao_fm(2), 1, 1, sample_val)
         WRITE (log_unit, '(A,ES14.6)') "  dipole_y(1,1): ", sample_val
         CALL cp_fm_get_element(tddft_env%dipole_ao_fm(3), 1, 1, sample_val)
         WRITE (log_unit, '(A,ES14.6,/)') "  dipole_z(1,1): ", sample_val
      END IF

      !-----------------------------------------------------------------------------------
      ! THEORY: Transition Density Matrix and Transition Dipole Moments in MRSF-TDA
      !-----------------------------------------------------------------------------------
      ! In MRSF-TDDFT with Tamm-Dancoff Approximation (TDA), we compute the UNRELAXED
      ! transition density matrix between excited states i and j.
      !
      ! 1. UNRELAXED vs RELAXED Transition Density:
      !    - UNRELAXED: Calculated directly from TDA response amplitudes (X vector only)
      !      without orbital relaxation effects.
      !    - RELAXED: Includes orbital relaxation contributions from solving Z-vector
      !      equations. Requires additional SCF-like procedure (not implemented here).
      !
      ! 2. MRSF-TDA Transition Density Formula:
      !    For MRSF theory with pure spin states, the transition density in MO basis is:
      !        P_ij^MO = P_ij^vir/vir + P_ij^occ/occ
      !    where:
      !    - P_ij^vir/vir: virtual-virtual block contribution from X_i^T * X_j
      !    - P_ij^occ/occ: occupied-occupied block contribution from -X_i * X_j^T
      !
      ! 3. Transition Dipole Moment:
      !    The transition dipole moment between states i and j is computed as:
      !        mu_ij = Tr[P_ij^AO * mu_AO]
      !    where mu_AO is the dipole integral matrix in AO basis.
      !    For oscillator strength: f_ij = (2/3) * Delta_E_ij * |mu_ij|^2
      !
      ! 4. MO to AO Transformation:
      !    P_ij^AO = C_alpha * P_ij^MO * C_beta^T
      !    where C_alpha and C_beta are MO coefficient matrices for alpha and beta spins.
      !-----------------------------------------------------------------------------------

      ! Allocate transition density matrices [nbf, nbf] for all state pairs
      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env)
      ALLOCATE (trden(tddft_env%nStates, tddft_env%nStates))
      DO ist = 1, tddft_env%nStates
         DO jst = 1, tddft_env%nStates
            CALL cp_fm_create(matrix=trden(ist, jst), &
                              matrix_struct=this%work1_nao_nao%matrix_struct)
            CALL cp_fm_set_all(matrix=trden(ist, jst), alpha=0.0_dp)
         END DO
      END DO

      ! Compute transition densities for all state pairs
      DO ist = 1, tddft_env%nStates
         DO jst = ist, tddft_env%nStates

            ! Copy amplitudes for state ist to work_noc_vir
            CALL cp_fm_set_all(matrix=this%work_noc_vir, alpha=0.0_dp)
            CALL cp_fm_to_fm(source=tddft_env%evects_mo(ist), &
                             destination=this%work_noc_vir)

            ! Copy amplitudes for state jst to work2_noc_vir
            CALL cp_fm_set_all(matrix=this%work2_noc_vir, alpha=0.0_dp)
            CALL cp_fm_to_fm(source=tddft_env%evects_mo(jst), &
                             destination=this%work2_noc_vir)

            ! *** DEBUG: Print excitation vector norms for first transition ***
            IF (log_unit > 0 .AND. ist == 1 .AND. jst == 2) THEN
               dipole_norm_x = cp_fm_frobenius_norm(tddft_env%evects_mo(ist))
               dipole_norm_y = cp_fm_frobenius_norm(tddft_env%evects_mo(jst))
               WRITE (log_unit, '(/,A,I0,A,I0,A)') "=== DEBUG: Excitation Vectors for transition ", ist, " -> ", jst, " ==="
               WRITE (log_unit, '(A,I0,A,ES14.6)') "  Norm of evects_mo(", ist, "): ", dipole_norm_x
               WRITE (log_unit, '(A,I0,A,ES14.6)') "  Norm of evects_mo(", jst, "): ", dipole_norm_y
               ! Sample amplitude (1,1)
               CALL cp_fm_get_element(tddft_env%evects_mo(ist), 1, 1, sample_val)
               WRITE (log_unit, '(A,I0,A,ES14.6)') "  evects_mo(", ist, ")(1,1): ", sample_val
               CALL cp_fm_get_element(tddft_env%evects_mo(jst), 1, 1, sample_val)
               WRITE (log_unit, '(A,I0,A,ES14.6,/)') "  evects_mo(", jst, ")(1,1): ", sample_val
            END IF

            ! Apply multiplicity-specific MRSF corrections to amplitudes
            ! These rearrange selected active orbital amplitudes for MRSF spin-mapping
            SELECT CASE (tddft_env%response_type)
            CASE (2) ! MRSF-S (Singlet)
               ! Correction for (i_lr1, a_lr1): scale by sqrt12
               CALL cp_fm_get_element(tddft_env%evects_mo(ist), i_lr1, a_lr1, val)
               CALL cp_fm_set_element(this%work_noc_vir, i_lr1, a_lr1, val*sqrt12)

               CALL cp_fm_get_element(tddft_env%evects_mo(jst), i_lr1, a_lr1, val)
               CALL cp_fm_set_element(this%work2_noc_vir, i_lr1, a_lr1, val*sqrt12)

               ! Correction for (i_lr2, a_lr2): scale by -sqrt12 (for ist) and sqrt12 (for jst)
               CALL cp_fm_get_element(tddft_env%evects_mo(ist), i_lr1, a_lr1, val)
               CALL cp_fm_set_element(this%work_noc_vir, i_lr2, a_lr2, -val*sqrt12)

               CALL cp_fm_get_element(tddft_env%evects_mo(jst), i_lr1, a_lr1, val)
               CALL cp_fm_set_element(this%work2_noc_vir, i_lr2, a_lr2, -val*sqrt12)

            CASE (3) ! MRSF-T (Triplet)
               ! Correction for (i_lr1, a_lr1): scale by sqrt12
               CALL cp_fm_get_element(tddft_env%evects_mo(ist), i_lr1, a_lr1, val)
               CALL cp_fm_set_element(this%work_noc_vir, i_lr1, a_lr1, val*sqrt12)

               CALL cp_fm_get_element(tddft_env%evects_mo(jst), i_lr1, a_lr1, val)
               CALL cp_fm_set_element(this%work2_noc_vir, i_lr1, a_lr1, val*sqrt12)

               ! Correction for (i_lr2, a_lr2): scale by sqrt12 (same sign for triplet)
               CALL cp_fm_get_element(tddft_env%evects_mo(ist), i_lr1, a_lr1, val)
               CALL cp_fm_set_element(this%work_noc_vir, i_lr2, a_lr2, val*sqrt12)

               CALL cp_fm_get_element(tddft_env%evects_mo(jst), i_lr1, a_lr1, val)
               CALL cp_fm_set_element(this%work2_noc_vir, i_lr2, a_lr2, val*sqrt12)

               ! Set (i_g, a_g) to zero
               CALL cp_fm_set_element(this%work_noc_vir, i_g, a_g, 0.0_dp)
               CALL cp_fm_set_element(this%work2_noc_vir, i_g, a_g, 0.0_dp)

               ! Set (i_d, a_d) to zero
               CALL cp_fm_set_element(this%work_noc_vir, i_d, a_d, 0.0_dp)
               CALL cp_fm_set_element(this%work2_noc_vir, i_d, a_d, 0.0_dp)

            END SELECT

            ! Compute transition density matrix
            ! NOTE: parallel_gemm slicing parameters (x_first_row/x_first_col):
            ! - Specify starting position in the ORIGINAL matrix (before any transpose)
            ! - In Fortran: first index = row, second index = col
            ! - Dimensions (m,n,k) specify the size of data to extract from that position
            ! - Example: A[noca,nvirb] with a_first_row=i0, a_first_col=j0 and transa='T'
            !   extracts A(i0:i0+k-1, j0:j0+m-1) then transposes to A^T[m,k]

            ! Revised trden(vir/vir)
            ! C[noccb+1:nbf, noccb+1:nbf] = 1.0_dp*(A[1:nocca, 1:nvirb])^T*B[1:nocca, 1:nvirb] + 0.0_dp*C
            CALL parallel_gemm('T', 'N', nvirb, nvirb, nocca, &
                               1.0_dp, this%work2_noc_vir, &   ! xv12j [nocca, nvirb]
                                       this%work_noc_vir, &    ! xv12i [nocca, nvirb]
                               0.0_dp, trden(ist, jst), &      ! output to trden[nbf, nbf]
                               c_first_row=noccb + 1, c_first_col=noccb + 1)

            ! Correction Case 1: i in {1,2}, k in {noca-1,noca}, j NOT in {1,2}
            ! C[noccb+1:noccb+2, noccb+3:nbf] = scalar*(A[nocca-1:nocca, 1:2])^T*B[nocca-1:nocca, 3:nvirb] + 1.0_dp*C
            scalar = sqrt2 - 1.0_dp
            CALL parallel_gemm('T', 'N', 2, nvirb - 2, 2, &
                               scalar, this%work2_noc_vir, &   ! xv12j
                                       this%work_noc_vir, &    ! xv12i
                               1.0_dp, trden(ist, jst), &
                               a_first_row=nocca - 1, a_first_col=1, &
                               b_first_row=nocca - 1, b_first_col=3, &
                               c_first_row=noccb + 1, c_first_col=noccb + 3)

            ! Correction Case 2: j in {1,2}, k in {noca-1,noca}, i NOT in {1,2}
            ! C[noccb+3:nbf, noccb+1:noccb+2] = scalar*(A[nocca-1:nocca, 3:nvirb])^T*B[nocca-1:nocca, 1:2] + 1.0_dp*C
            CALL parallel_gemm('T', 'N', nvirb - 2, 2, 2, &
                               scalar, this%work2_noc_vir, &   ! xv12j
                                       this%work_noc_vir, &    ! xv12i
                               1.0_dp, trden(ist, jst), &
                               a_first_row=nocca - 1, a_first_col=3, &
                               b_first_row=nocca - 1, b_first_col=1, &
                               c_first_row=noccb + 3, c_first_col=noccb + 1)

            ! Revised trden(occ/occ)
            ! C[1:nocca, 1:nocca] = -1.0_dp*A[1:nocca, 1:nvirb]*(B[1:nocca, 1:nvirb])^T + 1.0_dp*C
            CALL parallel_gemm('N', 'T', nocca, nocca, nvirb, &
                               -1.0_dp, this%work_noc_vir, &   ! xv12i
                                        this%work2_noc_vir, &  ! xv12j
                               1.0_dp, trden(ist, jst))        ! output to trden[nbf, nbf]

            ! Correction Case 1: i in {noca-1,noca}, k in {1,2}, j NOT in {noca-1,noca}
            ! C[nocca-1:nocca, 1:nocca-2] = scalar*A[nocca-1:nocca, 1:2]*(B[1:nocca-2, 1:2])^T + 1.0_dp*C
            scalar = -(sqrt2 - 1.0_dp)
            CALL parallel_gemm('N', 'T', 2, nocca - 2, 2, &
                               scalar, this%work_noc_vir, &    ! xv12i
                                       this%work2_noc_vir, &   ! xv12j
                               1.0_dp, trden(ist, jst), &
                               a_first_row=nocca - 1, a_first_col=1, &
                               b_first_row=1, b_first_col=1, &
                               c_first_row=nocca - 1, c_first_col=1)

            ! Correction Case 2: j in {noca-1,noca}, k in {1,2}, i NOT in {noca-1,noca}
            ! C[1:nocca-2, nocca-1:nocca] = scalar*A[1:nocca-2, 1:2]*(B[nocca-1:nocca, 1:2])^T + 1.0_dp*C
            CALL parallel_gemm('N', 'T', nocca - 2, 2, 2, &
                               scalar, this%work_noc_vir, &    ! xv12i
                                       this%work2_noc_vir, &   ! xv12j
                               1.0_dp, trden(ist, jst), &
                               a_first_row=1, a_first_col=1, &
                               b_first_row=nocca - 1, b_first_col=1, &
                               c_first_row=1, c_first_col=nocca - 1)

            ! MO to AO transformation: P_ij^AO = C_alpha * P_ij^MO * C_beta^T
            ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,1:nbf]*B[1:nbf,1:nbf] + 0.0*C
            CALL parallel_gemm('N', 'N', nbf, nbf, nbf, &
                               1.0_dp, gs_mos(1)%mos_full, &
                                       trden(ist, jst), &
                               0.0_dp, this%work2_nao_nao)
            ! C[1:nbf,1:nbf] = 1.0*A[1:nbf,1:nbf]*(B[1:nbf,1:nbf])^T + 0.0*C
            CALL parallel_gemm('N', 'T', nbf, nbf, nbf, &
                               1.0_dp, this%work2_nao_nao, &
                                       gs_mos(2)%mos_full, &
                               0.0_dp, this%work1_nao_nao)

            ! *** DEBUG: Print transition density norms ***
            IF (log_unit > 0 .AND. ist == 1 .AND. jst == 2) THEN
               trden_norm = cp_fm_frobenius_norm(trden(ist, jst))
               CALL cp_fm_maxabsval(trden(ist, jst), trden_max)
               WRITE (log_unit, '(/,A)') "=== DEBUG: Transition Density in MO Basis ==="
               WRITE (log_unit, '(A,ES14.6)') "  Frobenius norm: ", trden_norm
               WRITE (log_unit, '(A,ES14.6)') "  Max abs value:  ", trden_max
               CALL cp_fm_get_element(trden(ist, jst), 1, 1, sample_val)
               WRITE (log_unit, '(A,ES14.6,/)') "  trden_MO(1,1):  ", sample_val

               trden_norm = cp_fm_frobenius_norm(this%work1_nao_nao)
               CALL cp_fm_maxabsval(this%work1_nao_nao, trden_max)
               WRITE (log_unit, '(A)') "=== DEBUG: Transition Density in AO Basis (before symmetrization) ==="
               WRITE (log_unit, '(A,ES14.6)') "  Frobenius norm: ", trden_norm
               WRITE (log_unit, '(A,ES14.6)') "  Max abs value:  ", trden_max
               CALL cp_fm_get_element(this%work1_nao_nao, 1, 1, sample_val)
               WRITE (log_unit, '(A,ES14.6,/)') "  trden_AO(1,1):  ", sample_val
            END IF

            ! Symmetrize transition density: P_sym = P + P^T
            ! work2 = P_AO
            CALL cp_fm_scale_and_add( &
               alpha=0.0_dp, matrix_a=this%work2_nao_nao, &
               beta=1.0_dp, matrix_b=this%work1_nao_nao)
            ! evects_nao_nao(1) = P_AO^T (used as tmp matrix for transpose)
            CALL cp_fm_transpose(this%work1_nao_nao, this%evects_nao_nao(1))
            ! work2 += P_AO^T => work2 = P_AO + P_AO^T
            CALL cp_fm_scale_and_add( &
               alpha=1.0_dp, matrix_a=this%work2_nao_nao, &
               beta=1.0_dp, matrix_b=this%evects_nao_nao(1))
            ! Copy symmetrized result back to work1
            CALL cp_fm_to_fm(this%work2_nao_nao, this%work1_nao_nao)

            ! *** DEBUG: Print transition density after symmetrization ***
            IF (log_unit > 0 .AND. ist == 1 .AND. jst == 2) THEN
               trden_norm = cp_fm_frobenius_norm(this%work1_nao_nao)
               CALL cp_fm_maxabsval(this%work1_nao_nao, trden_max)
               WRITE (log_unit, '(A)') "=== DEBUG: Transition Density in AO Basis (after symmetrization) ==="
               WRITE (log_unit, '(A,ES14.6)') "  Frobenius norm: ", trden_norm
               WRITE (log_unit, '(A,ES14.6)') "  Max abs value:  ", trden_max
               CALL cp_fm_get_element(this%work1_nao_nao, 1, 1, sample_val)
               WRITE (log_unit, '(A,ES14.6,/)') "  trden_AO_sym(1,1): ", sample_val
            END IF

            ! Compute transition dipole moments using AO basis integrals
            ! mu_ij = -0.5 * Tr[P_ij^AO * mu_AO]
            ! For VELOCITY form: mu_velocity = Tr[P_AO * dS/dx] / (E_j - E_i)
            DO ideriv = 1, 3
               ! Compute mu_ij^xyz = -0.5 * Tr[trden(ist,jst) * dipole_ao_fm(ideriv)]
               CALL cp_fm_trace(this%work1_nao_nao, &
                                tddft_env%dipole_ao_fm(ideriv), &
                                tddft_env%trans_dipoles_mo(ist, jst, ideriv))
               tddft_env%trans_dipoles_mo(ist, jst, ideriv) = &
                  -0.5_dp*tddft_env%trans_dipoles_mo(ist, jst, ideriv)

               ! Apply energy normalization for VELOCITY form
               IF (tddfpt_control%dipole_form == tddfpt_dipole_velocity) THEN
                  IF (ist /= jst) THEN
                     ! Energy difference for excited state transitions
                     energ = tddft_env%evals(jst) - tddft_env%evals(ist)
                     IF (ABS(energ) > 1.0E-10_dp) THEN
                        tddft_env%trans_dipoles_mo(ist, jst, ideriv) = &
                           tddft_env%trans_dipoles_mo(ist, jst, ideriv)/energ
                     END IF
                  ELSE
                     ! Diagonal elements (ist == jst) correspond to permanent dipoles
                     ! For VELOCITY form, these should be zero (no excitation energy)
                     tddft_env%trans_dipoles_mo(ist, jst, ideriv) = 0.0_dp
                  END IF
               END IF

               ! *** DEBUG: Print trace values ***
               IF (log_unit > 0 .AND. ist == 1 .AND. jst == 2) THEN
                  IF (ideriv == 1) WRITE (log_unit, '(A)') "=== DEBUG: Transition Dipole Trace Calculation ==="
                  WRITE (log_unit, '(A,I0,A,ES14.6)') "  Trace[P_AO * dipole_", ideriv, "]: ", &
                     tddft_env%trans_dipoles_mo(ist, jst, ideriv)/(-0.5_dp)
                  WRITE (log_unit, '(A,I0,A,ES14.6)') "  trans_dipole_", ideriv, " (with -0.5 factor): ", &
                     tddft_env%trans_dipoles_mo(ist, jst, ideriv)
                  IF (ideriv == 3) WRITE (log_unit, '(/)')
               END IF
            END DO
         END DO
      END DO

      ! Set diagonal elements to zero (i=j are not transition dipoles)
      ! For transition dipole moments, we only consider transitions between different states
      ! Diagonal elements (i=i) would give permanent dipole moments, not transitions
      DO ist = 1, tddft_env%nstates
         tddft_env%trans_dipoles_mo(ist, ist, 1:3) = 0.0_dp
      END DO

      ! s2_mat
      ALLOCATE (s2_mat(tddft_env%nstates), source=0.0_dp)
      SELECT CASE (tddft_env%response_type)
      CASE (2)
         s2_mat(1:tddft_env%nStates) = 0.0_dp
      CASE (3)
         s2_mat(1:tddft_env%nStates) = 2.0_dp
      END SELECT

      ! *** summary header ***
      CALL integer_to_string(3, mult_str)
      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,1X,A2,A,1X,A)') lsd_str, "-MRSF-TDDFT states of multiplicity", TRIM(mult_str)
         SELECT CASE (tddfpt_control%dipole_form)
         CASE (tddfpt_dipole_berry)
            WRITE (log_unit, '(1X,A,/)') "Transition dipoles calculated using BERRY operator formulation"
         CASE (tddfpt_dipole_length)
            WRITE (log_unit, '(1X,A,/)') "Transition dipoles calculated using LENGTH formulation"
         CASE (tddfpt_dipole_velocity)
            WRITE (log_unit, '(1X,A,/)') "Transition dipoles calculated using VELOCITY formulation"
         CASE DEFAULT
            CPABORT("Unimplemented form of the dipole operator")
         END SELECT

         WRITE (log_unit, '(T3,A,T15,A,T37,A,T69,A)') "State", "Excitation     S^2", &
            "Transition dipole (a.u.)", "Oscillator"
         WRITE (log_unit, '(T3,A,T14,A,T37,A,T49,A,T61,A,T67,A)') "number", "energy (a.u.)", &
            "x", "y", "z", "strength (a.u.)"
         WRITE (log_unit, '(T3,79("-"))')
      END IF


      ! *** summary information ***
      DO istate = 1, tddft_env%nstates
         osc_strength = 2.0_dp/3.0_dp*tddft_env%evals(istate)* &
                        accurate_dot_product(tddft_env%trans_dipoles_mo(1, istate, 1:3), &
                                             tddft_env%trans_dipoles_mo(1, istate, 1:3))
         ostrength(istate) = osc_strength
         IF (log_unit > 0) THEN
            WRITE (log_unit, '(1X,I4,1X,F22.14,T30,F4.2,T34,3(1X,ES11.4E2),T70,ES12.4E2)') &
               istate, gs_energy + tddft_env%evals(istate), s2_mat(istate), &
               tddft_env%trans_dipoles_mo(1, istate, 1:3), osc_strength
         END IF
      END DO
      DEALLOCATE (s2_mat)

      ! punch a checksum for the regs
      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T2,A,E14.6)') 'TDDFT : CheckSum  =', &
            SQRT(SUM(tddft_env%evals**2))
      END IF

      max_amplitudes = tddft_env%nocca*tddft_env%nvirb
      ALLOCATE (amplitudes(max_amplitudes))
      ALLOCATE (i_indices(max_amplitudes))
      ALLOCATE (a_indices(max_amplitudes))
      ALLOCATE (indices(max_amplitudes))

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
            WRITE (log_unit, '(T7,I8,F10.5,A)') istate, (tddft_env%evals(istate)-tddft_env%evals(1))*evolt, " eV"
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
                  i, '(adp)', tddft_env%noccb + a, '(adp)', val
            END DO
         END IF
      END DO
      IF (log_unit > 0) THEN
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      IF (log_unit > 0 .AND. .FALSE.) THEN
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
               f = 2.0_dp/3.0_dp*energ*SUM(tddft_env%trans_dipoles_mo(istate, jstate, :)**2)
               WRITE (log_unit, '(3X,I0,1X,"->",1X,I0,T11,F11.5,3X,3(1X,ES11.4E2),2X,ES12.4E2)') &
                  istate, jstate, energ*evolt, tddft_env%trans_dipoles_mo(istate, jstate, 1:3), f
            END DO
         END DO
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      DO istate = 1, tddft_env%nStates
         CALL cp_fm_set_all(matrix=this%work_noc_vir, alpha=0.0_dp)
         CALL cp_fm_to_fm(source=tddft_env%evects_mo(istate), &
                          destination=this%work_noc_vir)
         CALL parallel_gemm('N', 'T', tddft_env%nbf, tddft_env%nocca, tddft_env%nvirb, &
                            1.0_dp, gs_mos(2)%mos_virt, &
                                    this%work_noc_vir, &
                            0.0_dp, tddft_env%evects(1, istate))

      END DO
      SELECT CASE (tddft_env%response_type)
      CASE (2) ! MRSF-S
         suffix = TRIM("mrs-")
      CASE (3) ! MRSF-S
         suffix = TRIM("mrt-")
      END SELECT

      tddfpt_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TDDFPT")
      tddfpt_print_section => section_vals_get_subs_vals(tddfpt_section, "PRINT")
      CALL tddfpt_print_nto_analysis(qs_env, &
                                     tddft_env%evects, tddft_env%evals, &
                                     ostrength, &
                                     gs_mos, &
                                     matrix_s(1)%matrix, &
                                     tddfpt_print_section, suffix)

      ! Deallocate transition density matrices
      DO ist = 1, tddft_env%nStates
         DO jst = 1, tddft_env%nStates
            CALL cp_fm_release(trden(ist, jst))
         END DO
      END DO
      DEALLOCATE (trden)

      CALL timestop(handle)

   END SUBROUTINE mrsf_print

! **************************************************************************************************
!> \brief ...
!> \param[inout] this ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE mrsf_deallocate(this, qs_env)
      CLASS(tddft_mrsf_type), INTENT(INOUT)              :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mrsf_deallocate'

      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      INTEGER                                            :: handle, ivec

      CALL timeset(routineN, handle)

      !--------------------
      ! DEALLOCATE  DP
      !--------------------

      IF (ALLOCATED(this%spc)) THEN
         DEALLOCATE (this%spc)
      END IF

      !--------------------
      ! DEALLOCATE  FM
      !--------------------
      CALL deallocate_fm_matrix(this%hmat_nao_nao)
      CALL deallocate_fm_matrix(this%evects_nao_nao)
      CALL deallocate_fm_matrix(this%fock)
      CALL deallocate_fm_matrix(this%focka_noc_noc)
      CALL deallocate_fm_matrix(this%fockb_vir_vir)
      CALL deallocate_fm_matrix(this%work_nao_vir)
      CALL deallocate_fm_matrix(this%work_noc_nao)
      CALL deallocate_fm_matrix(this%work_noc_noc)
      CALL deallocate_fm_matrix(this%work_nao_noc)
      CALL deallocate_fm_matrix(this%work_vir_vir)
      CALL deallocate_fm_matrix(this%work_noc_vir)
      CALL deallocate_fm_matrix(this%work2_noc_vir)
      CALL deallocate_fm_matrix(this%work1_nao_nao)
      CALL deallocate_fm_matrix(this%work2_nao_nao)
      CALL deallocate_fm_matrix(this%work1_nao_1)
      CALL deallocate_fm_matrix(this%work2_nao_1)
      CALL deallocate_fm_matrix(this%work_vira_1)
      CALL deallocate_fm_matrix(this%work_nocb_1)
      CALL deallocate_fm_matrix(this%work_nocb_vira)
      CALL deallocate_fm_matrix(this%work_nao_nocb)
      CALL deallocate_fm_matrix(this%work_aux_aux)
      CALL deallocate_fm_matrix(this%work_aux_nao)
      CALL deallocate_fm_matrix(this%work_nao_aux)
      CALL deallocate_fm_matrix(this%dipoles_mo)

      !--------------------
      ! DEALLOCATE  DBCSR
      !--------------------

      IF (ASSOCIATED(this%hmat_coulomb)) THEN
         DO ivec = 1, SIZE(this%hmat_coulomb)
            CALL dbcsr_deallocate_matrix(this%hmat_coulomb(ivec)%matrix)
         END DO
         DEALLOCATE (this%hmat_coulomb)
      END IF

      IF (ASSOCIATED(this%work_rho_nao_nao_symm)) THEN
         DO ivec = 1, SIZE(this%work_rho_nao_nao_symm)
            CALL dbcsr_deallocate_matrix(this%work_rho_nao_nao_symm(ivec)%matrix)
         END DO
         DEALLOCATE (this%work_rho_nao_nao_symm)
      END IF
      IF (ASSOCIATED(this%work_rho_nao_nao_asymm)) THEN
         DO ivec = 1, SIZE(this%work_rho_nao_nao_asymm)
            CALL dbcsr_deallocate_matrix(this%work_rho_nao_nao_asymm(ivec)%matrix)
         END DO
         DEALLOCATE (this%work_rho_nao_nao_asymm)
      END IF
      IF (ASSOCIATED(this%work_hmat_symm)) THEN
         DO ivec = 1, SIZE(this%work_hmat_symm)
            CALL dbcsr_deallocate_matrix(this%work_hmat_symm(ivec)%matrix)
         END DO
         DEALLOCATE (this%work_hmat_symm)
      END IF
      IF (ASSOCIATED(this%work_hmat_asymm)) THEN
         DO ivec = 1, SIZE(this%work_hmat_asymm)
            CALL dbcsr_deallocate_matrix(this%work_hmat_asymm(ivec)%matrix)
         END DO
         DEALLOCATE (this%work_hmat_asymm)
      END IF
      IF (ASSOCIATED(this%work_rho_aux_aux_symm)) THEN
         DO ivec = 1, SIZE(this%work_rho_aux_aux_symm)
            CALL dbcsr_deallocate_matrix(this%work_rho_aux_aux_symm(ivec)%matrix)
         END DO
         DEALLOCATE (this%work_rho_aux_aux_symm)
      END IF
      IF (ASSOCIATED(this%work_rho_aux_aux_asymm)) THEN
         DO ivec = 1, SIZE(this%work_rho_aux_aux_asymm)
            CALL dbcsr_deallocate_matrix(this%work_rho_aux_aux_asymm(ivec)%matrix)
         END DO
         DEALLOCATE (this%work_rho_aux_aux_asymm)
      END IF
      IF (ASSOCIATED(this%work_hmat_aux_symm)) THEN
         DO ivec = 1, SIZE(this%work_hmat_aux_symm)
            CALL dbcsr_deallocate_matrix(this%work_hmat_aux_symm(ivec)%matrix)
         END DO
         DEALLOCATE (this%work_hmat_aux_symm)
      END IF
      IF (ASSOCIATED(this%work_hmat_aux_asymm)) THEN
         DO ivec = 1, SIZE(this%work_hmat_aux_asymm)
            CALL dbcsr_deallocate_matrix(this%work_hmat_aux_asymm(ivec)%matrix)
         END DO
         DEALLOCATE (this%work_hmat_aux_asymm)
      END IF

      !--------------------
      ! DEALLOCATE  PW
      !--------------------

      CALL get_qs_env(qs_env=qs_env, &
                      pw_env=pw_env)

      CALL pw_env_get(pw_env, &
                      auxbas_pw_pool=auxbas_pw_pool)

      IF (ASSOCIATED(this%work_v_gspace)) THEN
         CALL auxbas_pw_pool%give_back_pw(this%work_v_gspace)
         DEALLOCATE (this%work_v_gspace)
      END IF
      IF (ASSOCIATED(this%work_v_rspace)) THEN
         CALL auxbas_pw_pool%give_back_pw(this%work_v_rspace)
         DEALLOCATE (this%work_v_rspace)
      END IF

      CALL timestop(handle)

   END SUBROUTINE mrsf_deallocate

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION is_zero_matrix_fm(matrix)
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix

      REAL(KIND=dp)                                      :: norm

      ! Calculating the Frobenius norm
      norm = cp_fm_norm(matrix, 'F')
      is_zero_matrix_fm = (norm < threshold)
   END FUNCTION is_zero_matrix_fm

END MODULE qs_tddft_theory_mrsf
