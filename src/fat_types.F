module fat_types
  use force_env_types, only: force_env_type, force_env_p_type
  use kinds, only: dp
  use pw_pool_types, only: pw_pool_type, pw_pool_p_type

#include "./base/base_uses.f90"

  implicit none

  private

  public :: &
       fat_env_type, &
       fat_env_get, &
       fat_env_set, &
       fat_env_create, &
       fat_env_release

  character(len=*), private, parameter :: modulen = "fat_types"

  logical, private, parameter :: is_debug = .true.

  type fat_env_type
     type(force_env_type), pointer   :: force_env_ref
     type(force_env_p_type), pointer :: force_env_sub(:)
     type(pw_pool_p_type), pointer   :: pw_pool_sub(:) => NULL()
     type(pw_pool_type), pointer     :: pw_pool_ref => NULL()
     integer                         :: i_iter, n_iter, n_frag
     logical                         :: is_converged
     logical, allocatable            :: extern(:)
     real(kind=dp)                   :: eps, e_old, ekin_ref, exc_ref, e_tot, eee, ene, enn
     real(kind=dp), allocatable      :: e_emb(:), e_sub(:), ekin_sub(:), exc_sub(:)
     integer                         :: output_unit
  end type fat_env_type

contains

! **************************************************************************************************
! **************************************************************************************************
  subroutine fat_env_get(fat_env, is_converged)
    type(fat_env_type), pointer :: fat_env
    logical, optional :: is_converged

    if (present(is_converged)) is_converged = fat_env%is_converged

  end subroutine fat_env_get

! **************************************************************************************************
! **************************************************************************************************
  subroutine fat_env_set(fat_env, is_converged)
    type(fat_env_type), pointer :: fat_env
    logical, optional :: is_converged

    if (present(is_converged)) fat_env%is_converged = is_converged

  end subroutine fat_env_set

! **************************************************************************************************
! **************************************************************************************************
  subroutine fat_env_create(fat_env)
    type(fat_env_type), pointer :: fat_env

    nullify(fat_env)
    allocate(fat_env)

  end subroutine fat_env_create

! **************************************************************************************************
! **************************************************************************************************
  subroutine fat_env_release(fat_env)
    type(fat_env_type), pointer :: fat_env

    if (associated(fat_env%pw_pool_sub)) deallocate(fat_env%pw_pool_sub)
    if (associated(fat_env%force_env_sub)) deallocate(fat_env%force_env_sub)
    if (allocated(fat_env%extern)) deallocate(fat_env%extern)
    if (allocated(fat_env%e_sub)) deallocate(fat_env%e_sub)
    if (allocated(fat_env%e_emb)) deallocate(fat_env%e_emb)
    if (allocated(fat_env%ekin_sub)) deallocate(fat_env%ekin_sub)
    if (allocated(fat_env%exc_sub)) deallocate(fat_env%exc_sub)
    if (associated(fat_env%pw_pool_sub)) deallocate(fat_env%pw_pool_sub)

    deallocate(fat_env)

  end subroutine fat_env_release

end module fat_types
