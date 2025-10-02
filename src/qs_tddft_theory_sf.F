!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Spin-Flip Time-Dependent Density Functional Theory (SF-TDDFT) implementation
!>
!> \details This module implements the Spin-Flip TDDFT method for computing electronic excited
!>          states from a high-spin ground state reference. SF-TDDFT uses excitations that flip
!>          electron spins (alpha to beta transitions), providing access to states of lower spin
!>          multiplicity from a high-spin reference configuration.
!>
!> \par Key features
!>      - Support for both ROHF and UHF reference states
!>      - Shows the same result between SF-TDDFT and SF-TDDFPT (Sternheimer approach) in CP2K.
!>
!> \author Konstantin Komarov (constlike@gmail.com)
!> \date Aug 2025
! **************************************************************************************************
MODULE qs_tddft_theory_sf

   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: tddfpt2_control_type
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
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_trace,&
                                              cp_fm_transpose
   USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_set_all,&
                                              cp_fm_set_element,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE input_constants,                 ONLY: tddfpt_dipole_berry,&
                                              tddfpt_dipole_length,&
                                              tddfpt_dipole_velocity
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kahan_sum,                       ONLY: accurate_dot_product
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE machine,                         ONLY: m_flush
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild
   USE qs_rho_types,                    ONLY: qs_rho_create,&
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
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_theory_sf'

   PUBLIC :: tddft_sf_type

! **************************************************************************************************
!> \brief SF-TDDFT derived type implementing the abstract theory interface
!>
!> \details Contains all working arrays and procedures specific to spin-flip TDDFT calculations.
!>          The implementation handles alpha to beta excitations.
!>
!> \par Key features
!>      - Support for both ROKS/ROHF and UKS/UHF reference states
!>      - ADMM acceleration for Hartree-Fock exchange kernel
!>      - Shows the same result between SF-TDDFT and SF-TDDFPT (Sternheimer approach) in CP2K.
!>
!> \author Konstantin Komarov (constlike@gmail.com) (Luber group)
!> \date Aug 2025
! **************************************************************************************************
   TYPE, EXTENDS(tddft_theory_type) :: tddft_sf_type

      ! AO-MO transformation buffer matrices
      TYPE(cp_fm_type), POINTER                          :: work_nao_vir => NULL()   ! [nbf, nvirb]
      TYPE(cp_fm_type), POINTER                          :: work_noc_nao => NULL()   ! [nocca, nbf]
      TYPE(cp_fm_type), POINTER                          :: work_noc_noc => NULL()   ! [nocca, nocca]
      TYPE(cp_fm_type), POINTER                          :: work_nao_noc => NULL()   ! [nbf, nocca]
      TYPE(cp_fm_type), POINTER                          :: work_vir_vir => NULL()   ! [nvirb, nvirb]
      TYPE(cp_fm_type), POINTER                          :: work_noc_vir => NULL()   ! [nocca, nvirb]

      ! General AO working matrices
      TYPE(cp_fm_type), POINTER                          :: work1_nao_nao => NULL()  ! [nbf, nbf]
      TYPE(cp_fm_type), POINTER                          :: work2_nao_nao => NULL()  ! [nbf, nbf]

      ! ADMM working matrices
      TYPE(cp_fm_type), POINTER                          :: work_aux_aux => NULL()   ! [nao_aux, nao_aux]
      TYPE(cp_fm_type), POINTER                          :: work_aux_nao => NULL()   ! [nao_aux, nao]
      TYPE(cp_fm_type), POINTER                          :: work_nao_aux => NULL()   ! [nao, nao_aux]

      ! Density
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: rho_fm_nao     ! (nSpins) [nbf, nbf]

      ! Density difference matrices for electron-hole analysis
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: dens_diff_ij   ! (nStates) [nbf, nbf] hole
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: dens_diff_ab   ! (nStates) [nbf, nbf] electron

   CONTAINS
      PROCEDURE, PASS(this) :: allocate   => sf_allocate
      PROCEDURE, PASS(this) :: init_guess => sf_init_guess
      PROCEDURE, PASS(this) :: mo_to_ao   => sf_mo_to_ao
      PROCEDURE, PASS(this) :: apply      => sf_apply
      PROCEDURE, PASS(this) :: ao_to_mo   => sf_ao_to_mo
      PROCEDURE, PASS(this) :: print      => sf_print
      PROCEDURE, PASS(this) :: deallocate => sf_deallocate
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief Allocate all working arrays for SF-TDDFT calculations
!> \param[inout] this      SF theory object
!> \param[in]    qs_env    QuickStep environment
!> \param[inout] tddft_env TDDFT environment containing system dimensions
! **************************************************************************************************
   SUBROUTINE sf_allocate(this, qs_env, tddft_env)
      CLASS(tddft_sf_type), INTENT(INOUT)                :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sf_allocate'

      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_s_aux_fit
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sab_aux_fit
      TYPE(admm_type), POINTER                           :: admm_env
      INTEGER                                            :: handle, ispin

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env)

      !=== hmat_nao_nao (1) [nbf, nbf] =================
      CALL allocate_fm_matrix(matrix=this%hmat_nao_nao, &
                              nmat=1, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nbf, &
                              blacs_env=blacs_env)

      !=== evects_nao_nao (1) [nbf, nbf] ===============
      CALL allocate_fm_matrix(matrix=this%evects_nao_nao, &
                              nmat=1, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nbf, &
                              blacs_env=blacs_env)

      !=== rho_fm_nao (nSpins) [nbf, nbf] ==============
      CALL allocate_fm_matrix(matrix=this%rho_fm_nao, &
                              nmat=tddft_env%nSpins, &
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

      !=== work_nao_vir [nbf, nvirb] ===================
      CALL allocate_fm_matrix(matrix=this%work_nao_vir, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      !=== work_noc_nao [nocca, nbf] ===================
      CALL allocate_fm_matrix(matrix=this%work_noc_nao, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nbf, &
                              blacs_env=blacs_env)

      !=== work_noc_noc [nocca, nocca] =================
      CALL allocate_fm_matrix(matrix=this%work_noc_noc, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nocca, &
                              blacs_env=blacs_env)

      !=== work_nao_noc [nbf, nocca] ===================
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

      !=== work1_nao_nao [nbf, nbf] ====================
      CALL allocate_fm_matrix(matrix=this%work1_nao_nao, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nbf, &
                              blacs_env=blacs_env)

      !=== work2_nao_nao [nbf, nbf] ====================
      CALL allocate_fm_matrix(matrix=this%work2_nao_nao, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nbf, &
                              blacs_env=blacs_env)

      !=== dens_diff_ij [nbf, nbf] =====================
      CALL allocate_fm_matrix(matrix=this%dens_diff_ij, &
                              nmat=tddft_env%nStates, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nbf, &
                              blacs_env=blacs_env)

      !=== dens_diff_ab [nbf, nbf] =====================
      CALL allocate_fm_matrix(matrix=this%dens_diff_ab, &
                              nmat=tddft_env%nStates, &
                              nrow=tddft_env%nbf, &
                              ncol=tddft_env%nbf, &
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
         CALL dbcsr_create(this%work_hmat_asymm(ispin)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric, &
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

      CALL timestop(handle)

   END SUBROUTINE sf_allocate

! **************************************************************************************************
!> \brief Generate initial guess vectors for SF-TDDFT based on orbital energy differences
!> \param[inout] this      SF theory object
!> \param[in]    qs_env    QuickStep environment
!> \param[inout] tddft_env TDDFT environment
!> \param[in]    gs_mos    ground state molecular orbitals
!> \param[in]    log_unit  output unit for logging
!> \details Creates energy-ordered alpha to beta single excitations as initial trial vectors
! **************************************************************************************************
   SUBROUTINE sf_init_guess(this, qs_env, tddft_env, gs_mos, log_unit)
      CLASS(tddft_sf_type), INTENT(INOUT)                :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), POINTER, INTENT(IN)               :: gs_mos
      INTEGER, INTENT(IN)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sf_init_guess'

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_m_occ
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      INTEGER                                            :: handle, imo_occ, imo_virt, &
                                                            istate, ispin, &
                                                            nStates_selected, istate_orig

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks)

      ! Transform Fock matrices to MO basis for energy corrections
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

      ! Extract Fock matrix blocks for energy corrections
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

      ! Create energy-ordered excitation list
      nStates_selected = tddft_env%nocca*(SIZE(gs_mos(2)%evals_full) - tddft_env%noccb)

      ALLOCATE (inds(nStates_selected))
      ALLOCATE (e_virt_m_occ(nStates_selected))

      ! Compute alpha_occ tp beta_virt excitation energies
      istate = 0
      DO imo_virt = tddft_env%noccb + 1, SIZE(gs_mos(2)%evals_full)
         DO imo_occ = 1, tddft_env%nocca
            istate = (imo_virt - tddft_env%noccb - 1)*tddft_env%nocca + imo_occ
            tddft_env%e_virt_minus_occ(istate) = gs_mos(2)%evals_full(imo_virt) - gs_mos(1)%evals_full(imo_occ)
         END DO
      END DO

      e_virt_m_occ(1:nStates_selected) = tddft_env%e_virt_minus_occ(1:nStates_selected)

      ! Sort excitations by energy
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

      ! Set up initial guess vectors as single alpha to beta xcitations
      DO istate = 1, tddft_env%nVec
         tddft_env%evals(istate) = e_virt_m_occ(istate)
         istate_orig = inds(istate)
         imo_occ = MOD(istate_orig - 1, tddft_env%nocca) + 1
         imo_virt = (istate_orig - 1)/tddft_env%nocca + 1
         CALL cp_fm_set_element(tddft_env%evects_mo(istate), imo_occ, imo_virt, 1.0_dp)
         IF (istate > tddft_env%nStates) CYCLE
         IF (log_unit > 0) &
            WRITE (log_unit, '(T7,I8,T24,I8,T37,A5,T45,I8,T54,A5,T60,F14.5)') &
            istate, imo_occ, '(alp)', tddft_env%noccb + imo_virt, '(bet)', e_virt_m_occ(istate)*evolt
      END DO

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T7,A,T50,I24)') 'Number of active states:', &
            SIZE(gs_mos(1)%evals_occ)*SIZE(gs_mos(2)%evals_virt)
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      DEALLOCATE (inds, e_virt_m_occ)

      CALL timestop(handle)

   END SUBROUTINE sf_init_guess

! **************************************************************************************************
!> \brief Transform trial vector from MO to AO basis
!> \param[inout] this      SF theory object
!> \param[inout] tddft_env TDDFT environment
!> \param[in]    gs_mos    ground state molecular orbitals
!> \param[in]    ivec      index of vector to transform
! **************************************************************************************************
   SUBROUTINE sf_mo_to_ao(this, tddft_env, gs_mos, ivec)
      CLASS(tddft_sf_type), INTENT(INOUT)                :: this
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER, INTENT(IN)                             :: gs_mos
      INTEGER, INTENT(IN)                                :: ivec

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sf_mo_to_ao'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! C[nbf,nvirb] = 1.0*A[nbf,nocca]*B[nocca,nvirb] + 0.0*C[nbf,nvirb]
      CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nvirb, tddft_env%nocca, &
                         1.0_dp, gs_mos(1)%mos_occ, &
                                 tddft_env%evects_mo(ivec), &
                         0.0_dp, this%work_nao_vir)
      ! C[nbf,nbf] = 1.0*A[nbf,nvirb]*(B[nbf,nvirb])^T + 0.0*C[nbf,nbf]
      CALL parallel_gemm('N', 'T', tddft_env%nbf, tddft_env%nbf, tddft_env%nvirb, &
                         1.0_dp, this%work_nao_vir, &
                                 gs_mos(2)%mos_virt, &
                         0.0_dp, this%evects_nao_nao(1))

      CALL timestop(handle)

   END SUBROUTINE sf_mo_to_ao

! **************************************************************************************************
!> \brief Apply TDDFT operators to generate sigma vectors
!> \param[inout] this            SF theory object
!> \param[in]    qs_env          QuickStep environment
!> \param[inout] tddft_env       TDDFT environment
!> \details Applies HFX kernel by decomposing density into symmetric and antisymmetric parts
!> \note  The HFX kernel is applied to both symmetric and antisymmetric
!>        components of the density matrix. When ADMM is enabled, the
!>        expensive four-center integrals are evaluated in a smaller
!>        auxiliary basis and then transformed back to the primary basis.
! **************************************************************************************************
   SUBROUTINE sf_apply(this, qs_env, tddft_env)
      CLASS(tddft_sf_type), INTENT(INOUT)                :: this
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sf_apply'

      TYPE(admm_type), POINTER                           :: admm_env
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      !-----------------------------------------------
      ! ERI contraction: Sigma = ERI * X [DBCSR type]
      !-----------------------------------------------
      IF(tddft_env%do_hfx)THEN

         IF (tddft_env%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
         END IF

         !--------------------------------------
         ! Get Symmetric part of X vector
         !--------------------------------------
         CALL cp_fm_set_all(matrix=this%work1_nao_nao, alpha=0.0_dp)
         ! work1 = 0.5*X
         CALL cp_fm_scale_and_add( &
            alpha=0.0_dp, matrix_a=this%work1_nao_nao, &
            beta=0.5_dp, matrix_b=this%evects_nao_nao(1))
         ! work2 = X^T
         CALL cp_fm_transpose(this%evects_nao_nao(1), this%work2_nao_nao)

         ! work1 += 0.5*work2
         CALL cp_fm_scale_and_add( &
            alpha=1.0_dp, matrix_a=this%work1_nao_nao, &
            beta=0.5_dp, matrix_b=this%work2_nao_nao)

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
                                  keep_sparsity=.true.)

            ! Apply HFX-ADMM to symmetric part
            CALL dbcsr_set(matrix=this%work_hmat_aux_symm(1)%matrix, alpha=0.0_dp)
            CALL tddft_hfx_matrix(matrix_ks=this%work_hmat_aux_symm, &
                                  rho_ao=this%work_rho_aux_aux_symm, &
                                  qs_env=qs_env)
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
            CALL tddft_hfx_matrix(matrix_ks=this%work_hmat_symm, &
                                  rho_ao=this%work_rho_nao_nao_symm, &
                                  qs_env=qs_env)
         END IF

         !--------------------------------------
         ! Get Antisymmetric part of X vector
         !--------------------------------------
         CALL cp_fm_set_all(matrix=this%work1_nao_nao, alpha=0.0_dp)

         ! work1 = 0.5*X
         CALL cp_fm_scale_and_add( &
            alpha=0.0_dp, matrix_a=this%work1_nao_nao, &
            beta=0.5_dp, matrix_b=this%evects_nao_nao(1))
         ! work1 -= 0.5*work2
         CALL cp_fm_scale_and_add( &
            alpha=1.0_dp, matrix_a=this%work1_nao_nao, &
            beta=-0.5_dp, matrix_b=this%work2_nao_nao)

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
            CALL tddft_hfx_matrix(matrix_ks=this%work_hmat_aux_asymm, &
                                  rho_ao=this%work_rho_aux_aux_asymm, &
                                  qs_env=qs_env)
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
            CALL tddft_hfx_matrix(matrix_ks=this%work_hmat_asymm, &
                                  rho_ao=this%work_rho_nao_nao_asymm, &
                                  qs_env=qs_env)
         END IF

         ! Convert results back to FM from DBCSR
         CALL cp_fm_set_all(matrix=this%work1_nao_nao, alpha=0.0_dp)
         CALL cp_fm_set_all(matrix=this%work2_nao_nao, alpha=0.0_dp)

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

         CALL cp_fm_set_all(matrix=this%hmat_nao_nao(1), alpha=0.0_dp)
         ! A = 0.0*A + 1.0*B
         CALL cp_fm_scale_and_add( &
            alpha=0.0_dp, matrix_a=this%hmat_nao_nao(1), &
            beta=1.0_dp, matrix_b=this%work1_nao_nao)
         ! A = 1.0*A + 1.0*B
         CALL cp_fm_scale_and_add( &
            alpha=1.0_dp, matrix_a=this%hmat_nao_nao(1), &
            beta=1.0_dp, matrix_b=this%work2_nao_nao)
      END IF

      CALL timestop(handle)

   END SUBROUTINE sf_apply

! **************************************************************************************************
!> \brief Transform sigma vectors back to MO basis with energy corrections
!> \param[inout] this      SF theory object
!> \param[inout] tddft_env TDDFT environment
!> \param[in]    gs_mos    ground state molecular orbitals
!> \param[in]    ivec      index of vector to transform
!> \details Applies back-transformation and adds orbital energy difference corrections
! **************************************************************************************************
   SUBROUTINE sf_ao_to_mo(this, tddft_env, gs_mos, ivec)
      CLASS(tddft_sf_type), INTENT(INOUT)                :: this
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER, INTENT(IN)                             :: gs_mos
      INTEGER, INTENT(IN)                                :: ivec

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sf_ao_to_mo'

      REAL(KIND=DP)                                      :: val, sig
      INTEGER                                            :: handle, a, i, ij

      CALL timeset(routineN, handle)

      ! C[nocca,nbf] = 1.0*(A[nbf,nocca])^T * B[nbf,nbf] + 0.0*C[nocca,nbf]
      CALL parallel_gemm('T', 'N', tddft_env%nocca, tddft_env%nbf, tddft_env%nbf, &
                         1.0_dp, gs_mos(1)%mos_occ, &
                                 this%hmat_nao_nao(1), &
                         0.0_dp, this%work_noc_nao)

      ! C[nocca,nvirb] = 1.0*A[noc,nbf] * B[nbf,nvirb] + 0.0*C[nocca,nvirb]
      CALL parallel_gemm('N', 'N', tddft_env%nocca, tddft_env%nvirb, tddft_env%nbf, &
                         1.0_dp, this%work_noc_nao, &
                                 gs_mos(2)%mos_virt, &
                         0.0_dp, tddft_env%sigma_vecs(ivec))

      ! Add reference-specific corrections based on SCF type
      SELECT CASE (tddft_env%reference_type)
      CASE (1) ! ROKS/ROHF
         ! Add ROHF corrections using Fock matrix contributions
         ! Sigma = X*F_b - F_a*X
         ! C[nocca,nvirb] = -1.0*A[nocca,nocca]*B[nocca,nvirb]+0.0*C[nocca,nvirb]
         CALL parallel_gemm('N', 'N', tddft_env%nocca, tddft_env%nvirb, tddft_env%nocca, &
                           -1.0_dp, this%focka_noc_noc, &
                                    tddft_env%evects_mo(ivec), &
                            0.0_dp, this%work_noc_vir)

         ! C[nocca,nvirb] = 1.0*A[nocca,nvirb]*B[nvirb,nvirb] + 1.0*C[nocca,nvirb]
         CALL parallel_gemm('N', 'N', tddft_env%nocca, tddft_env%nvirb, tddft_env%nvirb, &
                            1.0_dp, tddft_env%evects_mo(ivec), &
                                    this%fockb_vir_vir, &
                            1.0_dp, this%work_noc_vir)

         ! Add corrections to sigma_mo
         ! A = 1.0 * A + 1.0 * B
         CALL cp_fm_scale_and_add( &
            alpha=1.0_dp, matrix_a=tddft_env%sigma_vecs(ivec), &
            beta=1.0_dp, matrix_b=this%work_noc_vir)

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

               sig = sig + tddft_env%e_virt_minus_occ(ij)*val

               CALL cp_fm_set_element(matrix=tddft_env%sigma_vecs(ivec), &
                                      irow_global=i, &
                                      icol_global=a, &
                                      alpha=sig)
            END DO
         END DO

      END SELECT

      CALL timestop(handle)

   END SUBROUTINE sf_ao_to_mo

! **************************************************************************************************
!> \brief Analyze and print SF-TDDFT excited state properties
!> \param[inout] this            SF theory object
!> \param[in]    qs_env          QuickStep environment
!> \param[inout] tddft_env       TDDFT environment
!> \param[out]   ostrength       oscillator strengths
!> \param[in]    gs_mos          ground state molecular orbitals
!> \param[in]    tddfpt_control  TDDFPT control parameters
!> \param[in]    log_unit        output unit
!> \details Computes transition properties and performs comprehensive spin contamination analysis
! **************************************************************************************************
   SUBROUTINE sf_print(this, qs_env, tddft_env, ostrength, gs_mos, tddfpt_control, log_unit)
      CLASS(tddft_sf_type), INTENT(INOUT)                :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      REAL(kind=dp), DIMENSION(:), INTENT(OUT)           :: ostrength
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER, INTENT(IN)                             :: gs_mos
      TYPE(tddfpt2_control_type), POINTER, INTENT(IN)    :: tddfpt_control
      INTEGER, INTENT(IN)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sf_print'

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: tddfpt_section
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      CHARACTER(LEN=default_path_length)                 :: suffix
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: trans_dipoles
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: trans_dipoles_exc_exc
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: s2_mat, amplitudes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: indices, i_indices, a_indices
      INTEGER, DIMENSION(2)                              :: nmo_occ
      CHARACTER(len=2)                                   :: lsd_str
      CHARACTER(len=1)                                   :: mult_str
      REAL(kind=dp)                                      :: osc_strength, val, f, energ, gs_energy, &
         n_spin_imbalance, s2_spin, s2_beta_occ, s2_cross, s2_anti, s2_offdiag
      INTEGER                                            :: handle, ideriv, istate, ispin, &
         n_amplitudes, max_amplitudes, i, a, idx, jstate, nao, nspins

      CALL timeset(routineN, handle)

      SELECT CASE (tddft_env%reference_type)
      CASE (1)
         lsd_str = 'RO'
      CASE (2)
         lsd_str = ' U'
      END SELECT

      ! Print header with method and dipole operator information
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

      ! Compute ground state density matrices for each spin
      DO ispin = 1, nspins
         CALL parallel_gemm('N', 'T', nao, nao, nmo_occ(ispin), &
                            1.0_dp, gs_mos(ispin)%mos_occ, &
                                    gs_mos(ispin)%mos_occ, &
                            0.0_dp, this%rho_fm_nao(ispin))
      END DO

      ! Convert overlap matrix to CP_FM format for spin analysis, hmat_nao_nao is tmp matrix here
      CALL cp_fm_set_all(matrix=this%hmat_nao_nao(1), alpha=0.0_dp)
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, this%hmat_nao_nao(1))

      ! Transform dipole operators to MO basis
      DO ideriv = 1, 3
         CALL parallel_gemm('T', 'N', tddft_env%nocca, tddft_env%nvirb, tddft_env%nbf, &
                            1.0_dp, tddft_env%dipole_op_mos_occ(ideriv, 1), &
                                    gs_mos(2)%mos_virt, &
                            0.0_dp, this%dipoles_mo(ideriv))
      END DO

      ! Allocate arrays for transition analysis
      ALLOCATE (trans_dipoles(tddft_env%nstates, 3), source=0.0_dp)
      ALLOCATE (s2_mat(tddft_env%nstates), source=0.0_dp)

      ! Compute excited state properties for each state
      DO istate = 1, tddft_env%nStates
         ! C[nbf,nvirb] = 1.0*A[nbf,nocca]*B[nocca,nvirb] + 0.0*C[nbf,nvirb]
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nvirb, tddft_env%nocca, &
                            1.0_dp, gs_mos(1)%mos_occ, &
                                    tddft_env%evects_mo(istate), &
                            0.0_dp, this%work_nao_vir)
         ! C[nbf,nbf] = 1.0*A[nbf,nvirb]*(B[nbf,nvirb])^T + 0.0*C[nbf,nbf]
         CALL parallel_gemm('N', 'T', tddft_env%nbf, tddft_env%nbf, tddft_env%nvirb, &
                            1.0_dp, this%work_nao_vir, &
                                    gs_mos(2)%mos_virt, &
                            0.0_dp, this%evects_nao_nao(1))

         ! Compute density difference matrices for electron-hole analysis
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
                            0.0_dp, this%work_nao_noc)
         CALL parallel_gemm('N', 'T', tddft_env%nbf, tddft_env%nbf, tddft_env%nocca, &
                            1.0_dp, this%work_nao_noc, &
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
                            0.0_dp, this%work_nao_vir)
         CALL parallel_gemm('N', 'T', tddft_env%nbf, tddft_env%nbf, tddft_env%nvirb, &
                            1.0_dp, this%work_nao_vir, &
                                    gs_mos(2)%mos_virt, &
                            0.0_dp, this%dens_diff_ab(istate))

         !--------------------------------------------------------------
         ! Compute transition dipoles for istate
         !--------------------------------------------------------------
         DO ideriv = 1, 3
            CALL cp_fm_trace(tddft_env%evects_mo(istate), &
                             this%dipoles_mo(ideriv), &
                             trans_dipoles(istate, ideriv))
         END DO

         ! Compute spin expectation value <S^2> for excited state
         n_spin_imbalance = REAL(tddft_env%nocca - tddft_env%noccb, dp)
         s2_spin = 0.25_dp*n_spin_imbalance*(n_spin_imbalance - 2.0_dp)
         s2_beta_occ = REAL(tddft_env%noccb, dp) + 1.0_dp

         !--------------------------------------------------------------
         ! Compute S^2 expectation value for excited state:
         !--------------------------------------------------------------
         ! Step 1
!        !--------------------------------------------------------------
!        ! Compute fully relaxed AO-basis density matrix for the excited state:
!        !
!        !   rho_excited = rho_ground + delta_rho
!        !               = [rho_fm_nao(1) + dens_diff_ij] (alpha)
!        !               = [rho_fm_nao(2) + dens_diff_ab] (beta)
!        !
!        ! If orbital relaxation is included, delta_rho should incorporate those corrections.
!        ! In this context, we compute the "unrelaxed" density as ground + linear response.
!        !--------------------------------------------------------------
!        ! Construct alpha excited state density: work1_nao_nao = rho_fm_nao(1) + dens_diff_ij
         CALL cp_fm_scale_and_add(0.0_dp, this%work1_nao_nao, 1.0_dp, this%rho_fm_nao(1))
         CALL cp_fm_scale_and_add(1.0_dp, this%work1_nao_nao, 1.0_dp, this%dens_diff_ij(istate))
!        ! <--- this%work1_nao_nao now contains AO density for alpha (excited state candidate)

         ! scr1 = S * (rho_alpha_excited) * S
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, this%hmat_nao_nao(1), &
                                    this%work1_nao_nao, &
                            0.0_dp, this%work2_nao_nao)
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, this%work2_nao_nao, &
                                    this%hmat_nao_nao(1), &
                            0.0_dp, this%work1_nao_nao)

!        ! Construct beta excited state density: work2_nao_nao = rho_fm_nao(2) + dens_diff_ab
         CALL cp_fm_scale_and_add(0.0_dp, this%work2_nao_nao, 1.0_dp, this%rho_fm_nao(2))
         CALL cp_fm_scale_and_add(1.0_dp, this%work2_nao_nao, 1.0_dp, this%dens_diff_ab(istate))
!        ! <--- this%work2_nao_nao now contains AO density for beta (excited state candidate)

         ! s2_cross = -tr[(rho_beta_excited) * scr1]
         CALL cp_fm_trace(this%work2_nao_nao, this%work1_nao_nao, s2_cross)
         s2_cross = -s2_cross

         ! scr1 = S * dens_diff_ij * S [ta stands for dens_diff_ij, that is, the excited-state difference (alpha block)]
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, this%hmat_nao_nao(1), &
                                    this%dens_diff_ij(istate), &
                            0.0_dp, this%work1_nao_nao)
         CALL parallel_gemm('N', 'N', tddft_env%nbf, tddft_env%nbf, tddft_env%nbf, &
                            1.0_dp, this%work1_nao_nao, &
                                    this%hmat_nao_nao(1), &
                            0.0_dp, this%work2_nao_nao)

         ! s2_anti = -tr[dens_diff_ab * scr1]
         CALL cp_fm_trace(this%dens_diff_ab(istate), this%work2_nao_nao, s2_anti)
         s2_anti = -s2_anti

         ! s2_offdiag = (1/2) * tr[symmetrized(evects_nao_nao) * Smat]
         CALL cp_fm_set_all(matrix=this%work1_nao_nao, alpha=0.0_dp)
         CALL cp_fm_transpose(matrix=this%evects_nao_nao(1), matrixt=this%work1_nao_nao)
         CALL cp_fm_scale_and_add(0.0_dp, this%work2_nao_nao, 0.5_dp, this%evects_nao_nao(1))
         CALL cp_fm_scale_and_add(1.0_dp, this%work2_nao_nao, 0.5_dp, this%work1_nao_nao)
         CALL cp_fm_trace(this%work2_nao_nao, this%hmat_nao_nao(1), s2_offdiag)
         s2_offdiag = s2_offdiag/2.0_dp

         ! Total <S^2> = s2_spin + s2_beta_occ + s2_cross - s2_anti + s2_offdiag^2
         s2_mat(istate) = s2_spin + s2_beta_occ + s2_cross - s2_anti + s2_offdiag**2
      END DO

      ! Apply spin-flip normalization factor to transition dipoles
      trans_dipoles(:, :) = SQRT(2.0_dp)*trans_dipoles(:, :)

      ! Print summary of excited states
      DO istate = 1, tddft_env%nstates
         osc_strength = 2.0_dp/3.0_dp*tddft_env%evals(istate)* &
                        accurate_dot_product(trans_dipoles(istate, 1:3), &
                                             trans_dipoles(istate, 1:3))
         ostrength(istate) = osc_strength
         IF (log_unit > 0) THEN
            WRITE (log_unit, '(1X,I4,1X,F22.14,T30,F4.2,T34,3(1X,ES11.4E2),T70,ES12.4E2)') &
               istate, gs_energy + tddft_env%evals(istate), s2_mat(istate), &
               trans_dipoles(istate, 1:3), osc_strength
         END IF
      END DO
      DEALLOCATE (trans_dipoles, s2_mat)

      ! Print checksum for regression testing
      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T2,A,E14.6)') 'TDDFT : CheckSum  =', &
            SQRT(SUM(tddft_env%evals**2))
      END IF

      ! Print detailed excitation analysis
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
                  i, '(alp)', tddft_env%noccb + a, '(bet)', val
            END DO
         END IF
      END DO
      IF (log_unit > 0) THEN
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      DEALLOCATE (amplitudes, i_indices, a_indices, indices)

      ALLOCATE (trans_dipoles_exc_exc(tddft_env%nstates, tddft_env%nstates, 3), source=0.0_dp)

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
               f = 2.0_dp/3.0_dp*energ*SUM(trans_dipoles_exc_exc(istate, jstate, :)**2)
               WRITE (log_unit, '(3X,I0,1X,"->",1X,I0,T11,F11.5,3X,3(1X,ES11.4E2),2X,ES12.4E2)') &
                  istate, jstate, energ*evolt, trans_dipoles_exc_exc(istate, jstate, 1:3), f
            END DO
         END DO
         WRITE (log_unit, "(1X,A)") &
            "-------------------------------------------------------------------------------"
      END IF

      DEALLOCATE (trans_dipoles_exc_exc)
      DO istate = 1, tddft_env%nStates
         CALL cp_fm_set_all(matrix=this%work_noc_vir, alpha=0.0_dp)
         CALL cp_fm_to_fm(source=tddft_env%evects_mo(istate), &
                          destination=this%work_noc_vir)
         CALL parallel_gemm('N', 'T', tddft_env%nbf, tddft_env%nocca, tddft_env%nvirb, &
                            1.0_dp, gs_mos(2)%mos_virt, &
                                    this%work_noc_vir, &
                            0.0_dp, tddft_env%evects(1, istate))

      END DO

      ! Perform natural transition orbital analysis
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

   END SUBROUTINE sf_print

! **************************************************************************************************
!> \brief Clean up all SF-TDDFT allocations
!> \param[inout] this   SF theory object
!> \param[in]    qs_env QuickStep environment
! **************************************************************************************************
   SUBROUTINE sf_deallocate(this, qs_env)
      CLASS(tddft_sf_type), INTENT(INOUT)                :: this
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sf_deallocate'

      INTEGER                                            :: handle, ivec

      CALL timeset(routineN, handle)

      ! Suppress unused variable warning
      CALL get_qs_env(qs_env=qs_env)

      !--------------------
      ! DEALLOCATE  FM
      !--------------------

      CALL deallocate_fm_matrix(this%hmat_nao_nao)
      CALL deallocate_fm_matrix(this%evects_nao_nao)
      CALL deallocate_fm_matrix(this%fock)
      CALL deallocate_fm_matrix(this%work_nao_vir)
      CALL deallocate_fm_matrix(this%work_noc_nao)
      CALL deallocate_fm_matrix(this%work_noc_noc)
      CALL deallocate_fm_matrix(this%work_nao_noc)
      CALL deallocate_fm_matrix(this%work_vir_vir)
      CALL deallocate_fm_matrix(this%work_noc_vir)
      CALL deallocate_fm_matrix(this%work1_nao_nao)
      CALL deallocate_fm_matrix(this%work2_nao_nao)
      CALL deallocate_fm_matrix(this%dens_diff_ij)
      CALL deallocate_fm_matrix(this%dens_diff_ab)
      CALL deallocate_fm_matrix(this%work_aux_aux)
      CALL deallocate_fm_matrix(this%work_aux_nao)
      CALL deallocate_fm_matrix(this%work_nao_aux)
      CALL deallocate_fm_matrix(this%rho_fm_nao)
      CALL deallocate_fm_matrix(this%focka_noc_noc)
      CALL deallocate_fm_matrix(this%fockb_vir_vir)
      CALL deallocate_fm_matrix(this%dipoles_mo)

      !--------------------
      ! DEALLOCATE  DBCSR
      !--------------------

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

      CALL timestop(handle)

   END SUBROUTINE sf_deallocate

END MODULE qs_tddft_theory_sf
