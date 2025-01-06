!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Data types and basic operations for freeze-and-thaw (fat) embedding calculations
!> \par History
!>      2022 created [Lukas Schreder]
!>      09.2024 refactored [Konstantin Komarov]
!> \author Konstantin Komarov
! **************************************************************************************************
MODULE fat_types
   USE force_env_types, ONLY: force_env_type, force_env_p_type
   USE kinds, ONLY: dp
   USE pw_pool_types, ONLY: pw_pool_type, pw_pool_p_type
   USE pw_types, ONLY: pw_r3d_rs_type, pw_c1d_gs_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: &
      fat_env_type, &
      fat_env_get, &
      fat_env_set, &
      fat_env_create, &
      fat_env_release

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fat_types'

   LOGICAL, PARAMETER, PRIVATE :: is_debug = .TRUE.

! **************************************************************************************************
!> \brief Main type for freeze-and-thaw calculations containing all necessary data
!> \par Purpose
!>      This type holds all the data needed for freeze-and-thaw embedding calculations,
!>      including reference system, fragments, energies, densities and embedding potentials
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
! **************************************************************************************************
   TYPE fat_env_type

      ! Force environments
      TYPE(force_env_type), POINTER :: force_env_ref => NULL()  !< Reference force environment
      TYPE(force_env_p_type), POINTER :: force_env_sub(:) => NULL()  !< Fragment force environments

      ! Iteration control
      INTEGER :: i_iter = 0                 !< Current iteration number
      INTEGER :: n_iter = 0                 !< Maximum number of iterations
      INTEGER :: n_frag = 0                 !< Number of fragments
      LOGICAL :: is_converged = .FALSE.     !< Convergence status
      REAL(dp) :: eps = 0.0_dp              !< Convergence threshold

      ! External calculation flags
      LOGICAL, ALLOCATABLE :: extern(:)     !< External calculation flags per fragment

      ! Energies
      REAL(dp) :: e_old = 0.0_dp            !< Previous total energy
      REAL(dp) :: e_tot = 0.0_dp            !< Current total energy
      REAL(dp) :: ekin_ref = 0.0_dp         !< Reference kinetic energy
      REAL(dp) :: exc_ref = 0.0_dp          !< Reference XC energy
      REAL(dp) :: eee = 0.0_dp              !< Electron-electron interaction energy
      REAL(dp) :: ene = 0.0_dp              !< Electron-nuclear interaction energy
      REAL(dp) :: enn = 0.0_dp              !< Nuclear-nuclear interaction energy

      ! Fragment energies
      REAL(dp), ALLOCATABLE :: e_sub(:)     !< Fragment energies
      REAL(dp), ALLOCATABLE :: e_emb(:)     !< Embedding energies
      REAL(dp), ALLOCATABLE :: ekin_sub(:)  !< Fragment kinetic energies
      REAL(dp), ALLOCATABLE :: exc_sub(:)   !< Fragment XC energies

      ! Real space (rs) densities and potentials
      TYPE(pw_r3d_rs_type), POINTER :: rho_rs(:) => NULL()  !< Real space electron densities
      TYPE(pw_r3d_rs_type), POINTER :: v_emb_rs => NULL()   !< Real space embedding potential

      ! Reciprocal (g) space (gs) densities and potentials
      TYPE(pw_c1d_gs_type), POINTER :: rho_gs(:) => NULL()  !< Reciprocal space electron densities
      TYPE(pw_c1d_gs_type), POINTER :: v_emb_gs => NULL()   !< Reciprocal space embedding potential

      ! Plane wave pools
      TYPE(pw_pool_type), POINTER :: pw_pool_ref => NULL()  !< Reference system PW pool
      TYPE(pw_pool_p_type), POINTER :: pw_pool_sub(:) => NULL()  !< Fragment system PW pools

      ! I/O
      INTEGER :: output_unit = -1           !< Output unit number

   END TYPE fat_env_type

CONTAINS

! **************************************************************************************************
!> \brief Gets information from a fat environment
!> \param fat_env The fat environment to query
!> \param is_converged Optional output for convergence status
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
! **************************************************************************************************
   SUBROUTINE fat_env_get(fat_env, is_converged)
      TYPE(fat_env_type), POINTER                      :: fat_env
      LOGICAL, OPTIONAL                                :: is_converged

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_env_get'

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      IF (PRESENT(is_converged)) is_converged = fat_env%is_converged

   END SUBROUTINE fat_env_get

! **************************************************************************************************
!> \brief Sets information in a fat environment
!> \param fat_env The fat environment to modify
!> \param is_converged Optional convergence status to set
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
! **************************************************************************************************
   SUBROUTINE fat_env_set(fat_env, is_converged)
      TYPE(fat_env_type), POINTER                      :: fat_env
      LOGICAL, OPTIONAL                                :: is_converged

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_env_set'

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      IF (PRESENT(is_converged)) fat_env%is_converged = is_converged

   END SUBROUTINE fat_env_set

! **************************************************************************************************
!> \brief Creates and initializes a new fat environment
!> \param fat_env The fat environment to create (output)
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
! **************************************************************************************************
   SUBROUTINE fat_env_create(fat_env)
      TYPE(fat_env_type), POINTER                      :: fat_env

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_env_create'

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      NULLIFY (fat_env)
      ALLOCATE (fat_env)

   END SUBROUTINE fat_env_create

! **************************************************************************************************
!> \brief Releases a fat environment and all associated resources
!> \param fat_env The fat environment to release
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
! **************************************************************************************************
   SUBROUTINE fat_env_release(fat_env)
      TYPE(fat_env_type), POINTER                      :: fat_env

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_env_release'

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      IF (ASSOCIATED(fat_env%pw_pool_sub)) DEALLOCATE (fat_env%pw_pool_sub)
      IF (ASSOCIATED(fat_env%force_env_sub)) DEALLOCATE (fat_env%force_env_sub)
      IF (ALLOCATED(fat_env%extern)) DEALLOCATE (fat_env%extern)
      IF (ALLOCATED(fat_env%e_sub)) DEALLOCATE (fat_env%e_sub)
      IF (ALLOCATED(fat_env%e_emb)) DEALLOCATE (fat_env%e_emb)
      IF (ALLOCATED(fat_env%ekin_sub)) DEALLOCATE (fat_env%ekin_sub)
      IF (ALLOCATED(fat_env%exc_sub)) DEALLOCATE (fat_env%exc_sub)

      ! Release PW objects
      IF (ASSOCIATED(fat_env%rho_rs)) DEALLOCATE (fat_env%rho_rs)
      IF (ASSOCIATED(fat_env%v_emb_rs)) DEALLOCATE (fat_env%v_emb_rs)
      IF (ASSOCIATED(fat_env%rho_gs)) DEALLOCATE (fat_env%rho_gs)
      IF (ASSOCIATED(fat_env%v_emb_gs)) DEALLOCATE (fat_env%v_emb_gs)

      DEALLOCATE (fat_env)

   END SUBROUTINE fat_env_release

END MODULE fat_types
