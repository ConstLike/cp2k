!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!
! **************************************************************************************************
!> \brief   TDDFT theory abstraction interface.
!> \author  Konstantin Komarov (constlike@gmail.com)
!> \date    2025-05-16
!>
!> This module defines the abstract type [tddft_theory_type] for TDDFT-based excited-state solvers.
!> The interface provides a unified contract for all theory-specific implementations (e.g. SF-TDDFT,
!> MRSF-TDDFT, etc.) and allows easy integration of new excitation models. Each derived theory
!> type must implement the deferred procedures (type-bound) for all major steps of the iterative
!> Davidson diagonalization and operator application.
!>
!> **Usage:**
!>   - To add a new theory (e.g. SF), create a new module (e.g. qs_tddft_theory_sf.F)
!>     that `EXTENDS(tddft_theory_type)` and implements all deferred procedures.
!>   - In the main TDDFT driver, select the theory at runtime (via input section or logic),
!>     allocate the corresponding handler, and use only the abstract interface throughout
!>     the driver code. No changes to the main loop are needed for additional theories.
!>
!> **Contract (deferred procedures to implement):**
!>   - initialize_guess_vectors : set up initial guess vectors in MO space
!>   - mo_to_ao_transform       : transform trial vectors from MO to AO basis
!>   - apply_operator           : apply all HFX/XC/Coulomb/... operators in AO
!>   - ao_to_mo_transform       : transform operator results back to MO basis
!>   - allocate_theory_data     : allocate any additional data required by this theory
!>   - deallocate_theory_data   : free all theory-specific resources
!>
!> [tddft_info_type] is used to store universal Davidson/infrastructure data.
!> Theory-specific arrays and parameters must be private fields of each implementation.
! **************************************************************************************************
MODULE qs_tddft_theory

   USE cp_control_types,                ONLY: tddfpt2_control_type
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_kernel_types,                 ONLY: kernel_env_type
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
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
      INTEGER :: unused = 0
      ! Theory-specific data should be defined in child types (see SF, MRSF, etc.)
   CONTAINS
      !> Set up initial guess excitation vectors (MO space)
      PROCEDURE(init_guess_iface), DEFERRED :: initialize_guess_vectors
      !> MO->AO transformation for trial vectors
      PROCEDURE(mo_to_ao_iface), DEFERRED :: mo_to_ao_transform
      !> Apply all operators (HFX, XC, Coulomb, ...), build sigma vectors in AO space
      PROCEDURE(apply_op_iface), DEFERRED :: apply_operator
      !> AO->MO back transformation for resulting vectors
      PROCEDURE(ao_to_mo_iface), DEFERRED :: ao_to_mo_transform
      !> Allocate any additional theory-specific data
      PROCEDURE(allocate_theory_iface), DEFERRED :: allocate_theory_data
      !> Deallocate all theory-specific data/arrays
      PROCEDURE(deallocate_theory_iface), DEFERRED :: deallocate_theory_data
   END TYPE tddft_theory_type

! **************************************************************************************************
!> \brief Deferred procedure interfaces for theory-specific implementation.
!>        Each derived theory type (SF, MRSF, ...) must implement these.
! **************************************************************************************************
   ABSTRACT INTERFACE
      !> Set up initial guess vectors (theory-specific, e.g. SF, MRSF, ...)
      SUBROUTINE init_guess_iface(this, tddft_env, gs_mos, tddfpt_control, log_unit)
         IMPORT                                          :: tddft_theory_type, &
                                                            tddft_info_type, &
                                                            tddfpt_ground_state_mos, &
                                                            cp_logger_type, &
                                                            tddfpt2_control_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
         TYPE(tddfpt_ground_state_mos), &
            DIMENSION(:), POINTER                        :: gs_mos
         TYPE(tddfpt2_control_type), POINTER             :: tddfpt_control
         INTEGER                                         :: log_unit
      END SUBROUTINE
      !> Transform a given trial vector from MO to AO basis
      SUBROUTINE mo_to_ao_iface(this, tddft_env, gs_mos, ivec)
         IMPORT                                          :: tddft_theory_type, &
                                                            tddft_info_type, &
                                                            tddfpt_ground_state_mos
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
         TYPE(tddfpt_ground_state_mos), &
            DIMENSION(:), POINTER                        :: gs_mos
         INTEGER, INTENT(IN)                             :: ivec
      END SUBROUTINE
      !> Apply all HFX/XC/Coulomb/... operators for this theory, in AO basis
      SUBROUTINE apply_op_iface(this, qs_env, tddft_env, kernel_env, sub_env, &
                                tddfpt_control, ivec)
         IMPORT                                          :: tddft_theory_type, &
                                                            tddft_info_type, &
                                                            qs_environment_type, &
                                                            kernel_env_type, &
                                                            tddfpt_subgroup_env_type, &
                                                            tddfpt2_control_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(qs_environment_type), POINTER              :: qs_env
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
         TYPE(kernel_env_type)                           :: kernel_env
         TYPE(tddfpt_subgroup_env_type)                  :: sub_env
         TYPE(tddfpt2_control_type), POINTER             :: tddfpt_control
         INTEGER, INTENT(IN)                             :: ivec
      END SUBROUTINE
      !> Transform operator results from AO back to MO basis
      SUBROUTINE ao_to_mo_iface(this, tddft_env, gs_mos, ivec)
         IMPORT                                          :: tddft_theory_type, &
                                                            tddft_info_type, &
                                                            tddfpt_ground_state_mos
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(tddft_info_type),    INTENT(INOUT)         :: tddft_env
         TYPE(tddfpt_ground_state_mos), &
            DIMENSION(:), POINTER                        :: gs_mos
         INTEGER, INTENT(IN)                             :: ivec
      END SUBROUTINE
      !> Allocate any theory-specific working arrays, buffers, etc.
      SUBROUTINE allocate_theory_iface(this, tddft_env)
         IMPORT :: tddft_theory_type, tddft_info_type
         CLASS(tddft_theory_type), INTENT(INOUT) :: this
         TYPE(tddft_info_type),    INTENT(INOUT) :: tddft_env
      END SUBROUTINE
      !> Deallocate/free any theory-specific resources
      SUBROUTINE deallocate_theory_iface(this, tddft_env)
         IMPORT :: tddft_theory_type, tddft_info_type
         CLASS(tddft_theory_type), INTENT(INOUT) :: this
         TYPE(tddft_info_type),    INTENT(INOUT) :: tddft_env
      END SUBROUTINE
   END INTERFACE

! **************************************************************************************************
!> \par How to implement a new theory:
!>   1. Create a new module, e.g. qs_tddft_theory_sf, with
!>      TYPE, EXTENDS(tddft_theory_type) :: tddft_sf_theory_type
!>   2. Implement all deferred methods for your theory.
!>   3. Add any additional fields for theory-specific working arrays, if needed.
!>   4. Select and allocate your handler in the main TDDFT driver, e.g.
!>        SELECT CASE (theory_kind)
!>          CASE('SF')
!>            ALLOCATE(tddft_sf_theory_type :: theory_handler)
!>          CASE('MRSF')
!>            ALLOCATE(tddft_mrsf_theory_type :: theory_handler)
!>        END SELECT
!>   5. Use only theory_handler%method(...) in all driver routines.
! **************************************************************************************************

END MODULE qs_tddft_theory
