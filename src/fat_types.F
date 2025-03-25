MODULE fat_types
   USE force_env_types,                 ONLY: force_env_p_type,&
                                              force_env_type
   USE kinds,                           ONLY: dp
   USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                              pw_pool_type
   USE qs_local_rho_types,              ONLY: local_rho_type
   USE qs_rho0_types,                   ONLY: rho0_atom_type,&
                                              rho0_mpole_type
   USE qs_rho_atom_types,               ONLY: rho_atom_type
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

   LOGICAL, PRIVATE, PARAMETER :: is_debug = .TRUE.

   TYPE fat_env_type
      TYPE(force_env_type), POINTER   :: force_env => NULL()
      TYPE(force_env_type), POINTER   :: force_env_ref => NULL()
      TYPE(force_env_p_type), POINTER :: force_env_sub(:) => NULL()
      TYPE(pw_pool_p_type), POINTER   :: pw_pool_sub(:) => NULL()
      TYPE(pw_pool_type), POINTER     :: pw_pool_ref => NULL()
      INTEGER                         :: i_iter, n_iter, n_frag
      LOGICAL                         :: is_converged
      LOGICAL, ALLOCATABLE            :: extern(:)
      REAL(dp), ALLOCATABLE           :: e_emb(:)
      REAL(dp), ALLOCATABLE           :: e_sub(:)
      REAL(dp), ALLOCATABLE           :: ekin_sub(:)
      REAL(dp), ALLOCATABLE           :: exc_sub(:)
      REAL(dp), ALLOCATABLE           :: forces(:,:)
      REAL(dp)                        :: eps
      REAL(dp)                        :: e_old
      REAL(dp)                        :: ekin_ref
      REAL(dp)                        :: exc_ref
      REAL(dp)                        :: e_tot
      REAL(dp)                        :: ehartree
      INTEGER                         :: output_unit

      ! GAPW needs
      TYPE(rho0_mpole_type), POINTER  :: rho0_mpole => NULL()
      TYPE(rho0_atom_type), DIMENSION(:), POINTER :: rho0_atom_set => NULL()
      TYPE(rho_atom_type), DIMENSION(:), POINTER  :: rho_atom_set => NULL()
      TYPE(local_rho_type), POINTER   :: local_rho_set => NULL()
      LOGICAL                         :: use_gapw = .FALSE.
   END TYPE fat_env_type

CONTAINS

! **************************************************************************************************
! **************************************************************************************************
   SUBROUTINE fat_env_get(fat_env, is_converged)
      TYPE(fat_env_type), POINTER                        :: fat_env
      LOGICAL, OPTIONAL                                  :: is_converged

      IF (PRESENT(is_converged)) is_converged = fat_env%is_converged

   END SUBROUTINE fat_env_get

! **************************************************************************************************
! **************************************************************************************************
   SUBROUTINE fat_env_set(fat_env, is_converged)
      TYPE(fat_env_type), POINTER                        :: fat_env
      LOGICAL, OPTIONAL                                  :: is_converged

      IF (PRESENT(is_converged)) fat_env%is_converged = is_converged

   END SUBROUTINE fat_env_set

! **************************************************************************************************
! **************************************************************************************************
   SUBROUTINE fat_env_create(fat_env)
      TYPE(fat_env_type), POINTER                        :: fat_env

      NULLIFY (fat_env)
      ALLOCATE (fat_env)

   END SUBROUTINE fat_env_create

! **************************************************************************************************
! **************************************************************************************************
   SUBROUTINE fat_env_release(fat_env)
      TYPE(fat_env_type), POINTER                        :: fat_env

      IF (ASSOCIATED(fat_env%pw_pool_sub)) DEALLOCATE (fat_env%pw_pool_sub)
      IF (ASSOCIATED(fat_env%force_env_sub)) DEALLOCATE (fat_env%force_env_sub)
      IF (ALLOCATED(fat_env%extern)) DEALLOCATE (fat_env%extern)
      IF (ALLOCATED(fat_env%e_sub)) DEALLOCATE (fat_env%e_sub)
      IF (ALLOCATED(fat_env%e_emb)) DEALLOCATE (fat_env%e_emb)
      IF (ALLOCATED(fat_env%ekin_sub)) DEALLOCATE (fat_env%ekin_sub)
      IF (ALLOCATED(fat_env%exc_sub)) DEALLOCATE (fat_env%exc_sub)
      IF (ASSOCIATED(fat_env%pw_pool_sub)) DEALLOCATE (fat_env%pw_pool_sub)

      DEALLOCATE (fat_env)

   END SUBROUTINE fat_env_release

END MODULE fat_types
