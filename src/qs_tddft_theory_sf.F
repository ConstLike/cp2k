!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
!
MODULE qs_tddft_theory_sf

   USE qs_tddft_theory,                 ONLY: tddft_theory_type
   USE qs_tddft_types,                  ONLY: tddft_info_type
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE cp_control_types,                ONLY: tddfpt2_control_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE qs_kernel_types,                 ONLY: kernel_env_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: evolt
   USE util,                            ONLY: sort
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

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_theory_sf'

   PUBLIC :: tddft_sf_theory_type

   TYPE, EXTENDS(tddft_theory_type) :: tddft_sf_theory_type
      ! -- theory-specific fields if needed --
   CONTAINS
      PROCEDURE, PASS(this) :: allocate_theory_data     => sf_allocate_theory_data
      PROCEDURE, PASS(this) :: initialize_guess_vectors => sf_initialize_guess_vectors
      PROCEDURE, PASS(this) :: mo_to_ao_transform       => sf_mo_to_ao_transform
      PROCEDURE, PASS(this) :: apply_operator           => sf_apply_operator
      PROCEDURE, PASS(this) :: ao_to_mo_transform       => sf_ao_to_mo_transform
      PROCEDURE, PASS(this) :: print_theory_data        => sf_print_theory_data
      PROCEDURE, PASS(this) :: deallocate_theory_data   => sf_deallocate_theory_data
   END TYPE

CONTAINS

   SUBROUTINE sf_allocate_theory_data(this, tddft_env)
      CLASS(tddft_sf_theory_type), INTENT(INOUT) :: this
      TYPE(tddft_info_type),    INTENT(INOUT) :: tddft_env
      IF (.FALSE.) THEN
         print *, this%unused
         tddft_env = tddft_env
      END IF
   END SUBROUTINE sf_allocate_theory_data
   SUBROUTINE sf_initialize_guess_vectors(this, tddft_env, gs_mos, &
                                          tddfpt_control, log_unit)
      CLASS(tddft_sf_theory_type), INTENT(INOUT)         :: this
      TYPE(tddft_info_type), INTENT(INOUT)               :: tddft_env
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), POINTER                           :: gs_mos
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      INTEGER                                            :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sf_initialize_guess_vectors'

      INTEGER                                            :: handle, imo_occ, imo_virt, &
                                                            istate, &
                                                            nStates_selected, istate_orig
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: e_virt_m_occ

      CALL timeset(routineN, handle)

      IF (.FALSE.) THEN
         print *, this%unused
         NULLIFY(tddfpt_control)
      END IF

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
            WRITE (log_unit, '(T7,I8,T24,I8,T37,A5,T45,I8,T54,A5,T60,2F14.5)') &
            istate, imo_occ, '(alp)', gs_mos(2)%homo + imo_virt, '(bet)', e_virt_m_occ(istate)*evolt, &
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

   END SUBROUTINE sf_initialize_guess_vectors

   SUBROUTINE sf_mo_to_ao_transform(this, tddft_env, gs_mos, ivec)
      CLASS(tddft_sf_theory_type), INTENT(INOUT)      :: this
      TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), POINTER                        :: gs_mos
      INTEGER, INTENT(IN)                             :: ivec
      IF (.FALSE. .and. ivec == 1) THEN
         print *, this%unused
         tddft_env = tddft_env
         NULLIFY(gs_mos)
      END IF
   END SUBROUTINE sf_mo_to_ao_transform
   SUBROUTINE sf_apply_operator(this, qs_env, tddft_env, kernel_env, sub_env, &
                                tddfpt_control, ivec)
      CLASS(tddft_sf_theory_type), INTENT(INOUT)      :: this
      TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
      TYPE(qs_environment_type), POINTER              :: qs_env
      TYPE(kernel_env_type)                           :: kernel_env
      TYPE(tddfpt_subgroup_env_type)                  :: sub_env
      TYPE(tddfpt2_control_type), POINTER             :: tddfpt_control
      INTEGER, INTENT(IN)                             :: ivec
      IF (.FALSE. .AND. ivec == 1) THEN
         print *, this%unused
         tddft_env = tddft_env
         NULLIFY(qs_env)
         kernel_env = kernel_env
         sub_env = sub_env
         NULLIFY(tddfpt_control)
      END IF
   END SUBROUTINE sf_apply_operator
   SUBROUTINE sf_ao_to_mo_transform(this, tddft_env, gs_mos, ivec)
      CLASS(tddft_sf_theory_type), INTENT(INOUT)      :: this
      TYPE(tddft_info_type),    INTENT(INOUT)         :: tddft_env
      TYPE(tddfpt_ground_state_mos), &
         DIMENSION(:), POINTER                        :: gs_mos
      INTEGER, INTENT(IN)                             :: ivec
      IF (.FALSE. .and. ivec==1) THEN
         print *, this%unused
         tddft_env = tddft_env
         NULLIFY(gs_mos)
      END IF
   END SUBROUTINE sf_ao_to_mo_transform
   SUBROUTINE sf_print_theory_data(this, tddft_env)
      CLASS(tddft_sf_theory_type), INTENT(INOUT) :: this
      TYPE(tddft_info_type),    INTENT(INOUT) :: tddft_env
      IF (.FALSE.) THEN
         print *, this%unused
         tddft_env = tddft_env
      END IF
   END SUBROUTINE sf_print_theory_data
   SUBROUTINE sf_deallocate_theory_data(this, tddft_env)
      CLASS(tddft_sf_theory_type), INTENT(INOUT) :: this
      TYPE(tddft_info_type),    INTENT(INOUT) :: tddft_env
      IF (.FALSE.) THEN
         print *, this%unused
         tddft_env = tddft_env
      END IF
   END SUBROUTINE sf_deallocate_theory_data

END MODULE qs_tddft_theory_sf
