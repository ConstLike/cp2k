!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   TDDFT theory abstraction interface.
!> \author  Konstantin Komarov (constlike@gmail.com)
!> \date    Aug 2025
!>
!> This module defines the abstract type [tddft_theory_type] for TDDFT-based excited-state solvers.
!> The interface provides a unified contract for all theory-specific implementations (e.g. SF-TDDFT,
!> MRSF-TDDFT, etc.) and allows easy integration of new excitation models. Each derived theory
!> type must implement the deferred procedures (type-bound) for all major steps of the iterative
!> Davidson diagonalization and operator application.
!>
!> Usage:
!>   - To add a new theory (e.g. SF), create a new module (e.g. qs_tddft_theory_sf.F)
!>     that `EXTENDS(tddft_theory_type)` and implements all deferred procedures.
!>   - In the main TDDFT driver, select the theory at runtime (via input section or logic),
!>     allocate the corresponding handler, and use only the abstract interface throughout
!>     the driver code. No changes to the main loop are needed for additional theories.
!>
!> Contract (deferred procedures to implement):
!>   - allocate   : allocate any additional data required by this theory
!>   - init_guess : set up initial guess vectors in MO space
!>   - mo_to_ao   : transform trial vectors from MO to AO basis
!>   - apply      : apply all HFX/Coulomb/XC... operators in AO
!>   - ao_to_mo   : transform operator results back to MO basis
!>   - print      : print excitation analysis available for this theory
!>   - deallocate : free all theory-specific resources
!>
!> \par How to implement a NEW theory:
!>   1. Create a new module, e.g. qs_tddft_theory_NEW, with
!>      TYPE, EXTENDS(tddft_theory_type) :: tddft_NEW_type
!>   2. Implement all deferred methods for your theory.
!>   3. Add any additional fields for theory-specific working arrays, if needed.
!>   4. Select and allocate your handler in the main TDDFT driver, e.g.
!>        SELECT CASE (theory_kind)
!>          CASE('SF')
!>            ALLOCATE(tddft_SF_type :: theory)
!>          CASE('MRSF')
!>            ALLOCATE(tddft_MRSF_type :: theory)
!>          CASE('NEW')
!>            ALLOCATE(tddft_NEW_type :: theory)
!>        END SELECT
!>   5. Use only theory%init_guess(...) in all driver routines.
!>
!> [tddft_info_type] is used to store universal Davidson/infrastructure data.
!> Theory-specific arrays and parameters must be private fields of each implementation.
! **************************************************************************************************
MODULE qs_tddft_theory

   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: tddfpt2_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE qs_tddft_types,                  ONLY: tddft_info_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_theory'

   PUBLIC :: tddft_theory_type

! **************************************************************************************************
!> \brief Abstract base type for all TDDFT theory implementations.
!>        (Each derived type provides theory-specific logic for all main TDDFT steps.)
! **************************************************************************************************
   TYPE, ABSTRACT :: tddft_theory_type
      ! Dummy
      INTEGER                                            :: unused = 0

      !---------------------------------------------------------------------------
      ! Theory-specific data MUST be defined in child types (see SF, MRSF, etc.)
      !---------------------------------------------------------------------------

      !----------------------------------------------------
      ! General data cross all theories in Davidson solver:
      !----------------------------------------------------
      !> Per-spin AO Fock matrices (used in transformation to MO)
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fock

      !> Fock matrix for Alpha occupied orbitals block [nocca, nocca]
      TYPE(cp_fm_type), POINTER                          :: Focka_noc_noc => NULL()

      !> Fock matrix for Beta virtual orbitals block [nvirb, nvirb]
      TYPE(cp_fm_type), POINTER                          :: Fockb_vir_vir => NULL()

      !> Hamiltonian matrix in AO basis
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: hmat_ao_ao

      !> AO trial vectors
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: evects_ao_ao

      !> [HFX kernel] Symmetric AO density matrices [nSpins]
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: work_rho_ao_ao_symm

      !> [HFX kernel] Antisymmetric AO density matrices [nSpins]
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: work_rho_ao_ao_asymm

      !> [HFX kernel] Symmetric Hamiltonian-like matrices in AO space [nSpins]
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: work_hmat_symm

      !> [HFX kernel] Antisymmetric Hamiltonian-like matrices in AO space [nSpins]
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: work_hmat_asymm

      !> Dipoles (x,y,z) [nocca, nvirb]
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: dipoles_mo

   CONTAINS

      !> Allocate any additional theory-specific data
      PROCEDURE(allocate_iface), DEFERRED                :: allocate
      !> Set up initial guess excitation vectors (MO space)
      PROCEDURE(init_guess_iface), DEFERRED              :: init_guess
      !> MO->AO transformation for trial vectors
      PROCEDURE(mo_to_ao_iface), DEFERRED                :: mo_to_ao
      !> Apply all operators (HFX, XC, Coulomb, ...), build sigma vectors in AO space
      PROCEDURE(apply_op_iface), DEFERRED                :: apply
      !> AO->MO back transformation for resulting vectors
      PROCEDURE(ao_to_mo_iface), DEFERRED                :: ao_to_mo
      !> Print theory-specific excitation analysis
      PROCEDURE(print_iface), DEFERRED                   :: print
      !> Deallocate all theory-specific data/arrays
      PROCEDURE(deallocate_iface), DEFERRED              :: deallocate

   END TYPE tddft_theory_type

! **************************************************************************************************
!> \brief Deferred procedure interfaces for theory-specific implementation.
!>        Each derived theory type (SF, MRSF, ...) must implement these.
! **************************************************************************************************
   ABSTRACT INTERFACE

! **************************************************************************************************
!> \brief Allocate any theory-specific working arrays, buffers, etc.
!> \param[inout] this ...
!> \param qs_env ...
!> \param[inout] tddft_env ...
! **************************************************************************************************
      SUBROUTINE allocate_iface(this, qs_env, tddft_env)
         IMPORT                                          :: cp_blacs_env_type, &
                                                            qs_environment_type, &
                                                            tddft_info_type, &
                                                            tddft_theory_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(qs_environment_type), POINTER              :: qs_env
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env

      END SUBROUTINE

! **************************************************************************************************
!> \brief Set up initial guess vectors (theory-specific, e.g. SF, MRSF, ...)
!> \param[inout] this ...
!> \param qs_env ...
!> \param[inout] tddft_env ...
!> \param[in] gs_mos ...
!> \param[in] log_unit ...
! **************************************************************************************************
      SUBROUTINE init_guess_iface(this, qs_env, tddft_env, gs_mos, log_unit)
         IMPORT                                          :: qs_environment_type, &
                                                            tddfpt_ground_state_mos, &
                                                            tddft_info_type, &
                                                            tddft_theory_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(qs_environment_type), POINTER              :: qs_env
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
         TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
            POINTER, INTENT(IN)                          :: gs_mos
         INTEGER, INTENT(IN)                             :: log_unit

      END SUBROUTINE

! **************************************************************************************************
!> \brief Transform a given trial vector from MO to AO basis
!> \param[inout] this ...
!> \param[inout] tddft_env ...
!> \param[in] gs_mos ...
!> \param[in] ivec ...
! **************************************************************************************************
      SUBROUTINE mo_to_ao_iface(this, tddft_env, gs_mos, ivec)
         IMPORT                                          :: tddfpt_ground_state_mos, &
                                                            tddft_info_type, &
                                                            tddft_theory_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
         TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
            POINTER, INTENT(IN)                          :: gs_mos
         INTEGER, INTENT(IN)                             :: ivec

      END SUBROUTINE

! **************************************************************************************************
!> \brief Apply all HFX/XC/Coulomb/... operators for this theory, in AO basis
!> \param[inout] this ...
!> \param qs_env ...
!> \param[input] tddft_env ...
!> \param[in] tddfpt_control ...
! **************************************************************************************************
      SUBROUTINE apply_op_iface(this, qs_env, tddft_env, tddfpt_control)
         IMPORT                                          :: qs_environment_type, &
                                                            tddft_info_type, &
                                                            tddft_theory_type, &
                                                            tddfpt2_control_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(qs_environment_type), POINTER              :: qs_env
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
         TYPE(tddfpt2_control_type), POINTER, INTENT(IN) :: tddfpt_control

      END SUBROUTINE

! **************************************************************************************************
!> \brief Transform operator results from AO back to MO basis
!> \param[inout] this ...
!> \param[inout] tddft_env ...
!> \param[in] gs_mos ...
!> \param[in] ivec ...
! **************************************************************************************************
      SUBROUTINE ao_to_mo_iface(this, tddft_env, gs_mos, ivec)
         IMPORT                                          :: tddfpt_ground_state_mos, &
                                                            tddft_info_type, &
                                                            tddft_theory_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
         TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
            POINTER, INTENT(IN)                          :: gs_mos
         INTEGER, INTENT(IN)                             :: ivec

      END SUBROUTINE

! **************************************************************************************************
!> \brief Print theory-specific excitation analysis
!> \param[inout] this ...
!> \param qs_env ...
!> \param[inout] tddft_env ...
!> \param[out] ostrength ...
!> \param[in] gs_mos ...
!> \param[in] tddfpt_control ...
!> \param[in] log_unit ...
! **************************************************************************************************
      SUBROUTINE print_iface(this, qs_env, tddft_env, ostrength, gs_mos, tddfpt_control, log_unit)
         IMPORT                                          :: qs_environment_type, &
                                                            tddfpt2_control_type, &
                                                            tddfpt_ground_state_mos, &
                                                            tddft_info_type, &
                                                            tddft_theory_type, &
                                                            dp
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(qs_environment_type), POINTER              :: qs_env
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
         TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
            POINTER, INTENT(IN)                          :: gs_mos
         REAL(kind=dp), DIMENSION(:), INTENT(OUT)        :: ostrength
         TYPE(tddfpt2_control_type), POINTER, INTENT(IN) :: tddfpt_control
         INTEGER, INTENT(IN)                             :: log_unit

      END SUBROUTINE

! **************************************************************************************************
!> \brief Deallocate/free any theory-specific resources
!> \param this ...
!> \param qs_env ...
! **************************************************************************************************
      SUBROUTINE deallocate_iface(this, qs_env)
         IMPORT                                          :: qs_environment_type, &
                                                            tddft_theory_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(qs_environment_type), POINTER              :: qs_env

      END SUBROUTINE

   END INTERFACE

END MODULE qs_tddft_theory
