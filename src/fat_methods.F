module fat_methods
  use atomic_kind_types, only: atomic_kind_type, get_atomic_kind, get_atomic_kind_set
  use cp_control_types, only: dft_control_type
  use cp_files, only: open_file, close_file
  use cp_log_handling, only: cp_logger_type, cp_logger_get_default_io_unit, cp_get_default_logger
  use cp_realspace_grid_cube, only: cp_pw_to_cube, cp_cube_to_pw
  use dbcsr_api, only: dbcsr_p_type
  use fat_types, only: fat_env_type
  use force_env_types, only: force_env_type
  use input_section_types, only: section_vals_type, section_vals_val_get, section_vals_get_subs_vals
  use kinds, only: dp, default_string_length
  use particle_list_types, only: particle_list_type
  use particle_types, only: particle_type
  use pw_env_types, only: pw_env_type, pw_env_get
  use pw_grid_types, only: pw_grid_type, HALFSPACE, PW_MODE_DISTRIBUTED
  use pw_methods, only: pw_axpy, pw_integral_ab, pw_transfer, pw_set, pw_zero
  use pw_poisson_methods, only: pw_poisson_solve
  use pw_poisson_types, only: pw_poisson_type
  use pw_pool_types, only: pw_pool_type
  use pw_types, only: pw_type, pw_p_type, REALDATA3D, COMPLEXDATA1D, REALSPACE, RECIPROCALSPACE
  use qs_collocate_density, only: calculate_drho_elec_dR, calculate_drho_core
  use qs_energy_types, only: qs_energy_type
  use qs_environment_types, only: qs_environment_type, get_qs_env, set_qs_env
  use qs_force, only: qs_calc_energy_force
  use qs_force_types, only: qs_force_type
  use qs_rho_types, only: qs_rho_type, qs_rho_get
  use qs_subsys_types, only: qs_subsys_type, qs_subsys_get
  use realspace_grid_cube, only: pw_to_simple_volumetric
  use string_utilities, only: integer_to_string
  use xc, only: xc_vxc_pw_create
  use xc_rho_cflags_types, only: xc_rho_cflags_type
  use xc_derivatives, only: xc_functionals_get_needs
! use physcon

#include "./base/base_uses.f90"

  implicit none

  private

  public :: &
       fat_finalize, &
       fat_init, &
       fat_is_done, &
       fat_iterate, &
       fat_results

  logical, private, parameter :: is_debug = .true.

  integer, save :: log_counter = 0

  character(len=*), private, parameter :: &
       modulen = "fat_methods", &
       fat_title_fmt = '(/, T2, "FAT| ", 4(A6, ", "), 9(A26, ", "), A26)', &
       fat_data_fmt = '(/, T2, "FAT| ", 4(I6, ", "), 9(F26.15, ", "), F26.15)'

contains

! **************************************************************************************************
!> \brief Creates a pw for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pw_create_frag(fat_env, pw, i_frag, rspace)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pw
    integer :: i_frag
    logical, optional :: rspace

    character(len=*), parameter :: routinen = "fat_pw_create_frag"

    integer :: use_data, in_space
    logical :: my_rspace

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    if (present(rspace)) then
       my_rspace = rspace
    else
       my_rspace = .true.
    end if

    if (my_rspace) then
       use_data = realdata3d
       in_space = realspace
    else
       use_data = complexdata1d
       in_space = reciprocalspace
    end if

    allocate(pw)
    call fat_env%pw_pool_sub(i_frag)%pool%create_pw(pw, use_data=use_data, in_space=in_space)
    call pw_zero(pw)

  end subroutine fat_pw_create_frag

! **************************************************************************************************
!> \brief Creates a pw for the reference system in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pw_create_ref(fat_env, pw, rspace)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), intent(out) :: pw
    logical, optional :: rspace

    character(len=*), parameter :: routinen = "fat_pw_create_ref"

    integer :: use_data, in_space
    logical :: my_rspace

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    if (present(rspace)) then
       my_rspace = rspace
    else
       my_rspace = .true.
    end if

    if (my_rspace) then
       use_data = realdata3d
       in_space = realspace
    else
       use_data = complexdata1d
       in_space = reciprocalspace
    end if

    call fat_env%pw_pool_ref%create_pw(pw, use_data=use_data, in_space=in_space)
    call pw_zero(pw)

  end subroutine fat_pw_create_ref

! **************************************************************************************************
!> \brief Releases a pw for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pw_release_frag(fat_env, pw, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pw
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_pw_release_frag"

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(pw)
    deallocate(pw)

    if (is_debug) print *, "DBG| released"

  end subroutine fat_pw_release_frag

! **************************************************************************************************
!> \brief Releases a pw for the reference system in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pw_release_ref(fat_env, pw)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), intent(inout) :: pw

    character(len=*), parameter :: routinen = "fat_pw_release_ref"

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call fat_env%pw_pool_ref%give_back_pw(pw)

    if (is_debug) print *, "DBG| released"

  end subroutine fat_pw_release_ref

! **************************************************************************************************
!> \brief Creates an array of pws for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pws_create_frag(fat_env, pws, i_frag, n, rspace)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pws(:)
    integer :: i_frag, n
    logical, optional :: rspace

    character(len=*), parameter :: routinen = "fat_pws_create_frag"

    integer :: i
    type(pw_type), pointer :: pw_tmp

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    allocate(pws(n))
    do i = 1, n
       pw_tmp = pws(i)
       call fat_pw_create_frag(fat_env, pw_tmp, i_frag, rspace)
    end do

  end subroutine fat_pws_create_frag

! **************************************************************************************************
!> \brief Creates an array of pws for the reference system in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pws_create_ref(fat_env, pws, n, rspace)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pws(:)
    integer :: n
    logical, optional :: rspace

    character(len=*), parameter :: routinen = "fat_pws_create_ref"

    integer :: i

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    allocate(pws(n))
    do i = 1, n
       call fat_pw_create_ref(fat_env, pws(i), rspace)
    end do

  end subroutine fat_pws_create_ref

! **************************************************************************************************
!> \brief Releases an array of pws for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pws_release_frag(fat_env, pws, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pws(:)
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_pws_release_frag"

    type(pw_type), pointer :: pw_tmp
    integer :: i

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    do i = 1, size(pws)
       pw_tmp => pws(i)
       call fat_pw_release_frag(fat_env, pw_tmp, i_frag)
    end do

    if (is_debug) print *, "DBG| released"

  end subroutine fat_pws_release_frag

! **************************************************************************************************
!> \brief Releases an array of pws for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pws_release_ref(fat_env, pws)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pws(:)

    character(len=*), parameter :: routinen = "fat_pws_release_ref"

    integer :: i

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    do i = 1, size(pws)
       call fat_pw_release_ref(fat_env, pws(i))
    end do
    DEALLOCATE(pws)

    if (is_debug) print *, "DBG| released"

  end subroutine fat_pws_release_ref

! **************************************************************************************************
!> \brief Gets the electron density for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_rho_frag(fat_env, rho, i_frag, rspace)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: rho(:)
    integer :: i_frag
    logical, optional :: rspace

    character(len=*), parameter :: routinen = "fat_rho_frag"

    character(len=default_string_length) :: label
    logical :: my_rspace
    type(pw_type), pointer :: rho_tmp
    type(qs_environment_type), pointer :: qs_env
    type(qs_rho_type), pointer :: qs_rho

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    if (present(rspace)) then
       my_rspace = rspace
    else
       my_rspace = .true.
    end if

    ! Initialize rho with correct data structure and grid mappings
    ! This is required even for external density to ensure proper parallel setup
    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, rho=qs_rho)
    if (my_rspace) then
       call qs_rho_get(qs_rho, rho_r=rho)
    else
       call qs_rho_get(qs_rho, rho_g=rho)
    end if

    ! We map external density into the properly initialized rho
    if (fat_env%extern(i_frag)) then
       call integer_to_string(i_frag, label)
       if (my_rspace) then
          call cp_cube_to_pw(rho(1), "extern_"//trim(label)//".cube", 1.0_dp, .false.)
       else
          call fat_pw_create_frag(fat_env, rho_tmp, i_frag)
          call cp_cube_to_pw(rho_tmp, "extern_"//trim(label)//".cube", 1.0_dp, .false.)
          call pw_transfer(rho_tmp, rho(1))
          call fat_pw_release_frag(fat_env, rho_tmp, i_frag)
       end if
    end if

!   if (is_debug) call print_pw_statistics(rho(1), "rho", i_frag)

  end subroutine fat_rho_frag


  subroutine write_pw_r_to_cube(fat_env, pw, filename)
    type(fat_env_type), intent(in) :: fat_env
    type(pw_type), intent(in) :: pw
    character(len=*), intent(in) :: filename

    integer :: unit_nr
    type(cp_logger_type), pointer :: logger
    TYPE(particle_list_type), POINTER :: particles
    type(particle_type), dimension(:), pointer :: particle_set
    type(qs_environment_type), pointer :: qs_env
    type(dft_control_type), pointer :: dft_control
    type(pw_env_type), pointer :: pw_env
    type(qs_subsys_type), pointer :: subsys
    type(pw_pool_type), pointer :: pw_pool
    integer, pointer :: stride(:)

    logger => cp_get_default_logger()
    unit_nr = cp_logger_get_default_io_unit(logger)
    log_counter = log_counter + 1

    if (log_counter>6) return

    print *, "KKK| write CUBE file number:", filename

    if (unit_nr > 0) then
      qs_env => fat_env%force_env_ref%qs_env

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, pw_env=pw_env, &
                      particle_set=particle_set, subsys=subsys)
      CALL pw_env_get(pw_env, auxbas_pw_pool=pw_pool)
      CALL qs_subsys_get(subsys, particles=particles)

      allocate(stride(3))
      stride = [1, 1, 1]

      call open_file(file_name=trim(filename), file_status="REPLACE", &
                     file_action="WRITE", unit_number=unit_nr)
      call cp_pw_to_cube(pw, unit_nr, trim(filename), particles=particles, stride=stride)
      call close_file(unit_number=unit_nr)
      deallocate(stride)
    end if
  end subroutine write_pw_r_to_cube

  subroutine write_pw_g_to_cube(fat_env, pw, filename, i_frag)
    type(fat_env_type), intent(in) :: fat_env
    type(pw_type), intent(in) :: pw
    character(len=*), intent(in) :: filename
    integer, intent(in) :: i_frag

    type(pw_type), pointer :: pw_real
    logical :: pool_exists

    if (.not. associated(fat_env%pw_pool_sub)) then
      print *, "Error: pw_pool_sub is not associated in fat_env"
      return
    end if

    if (i_frag < 1 .or. i_frag > size(fat_env%pw_pool_sub)) then
      print *, "Error: Invalid fragment index", i_frag
      return
    end if

    pool_exists = associated(fat_env%pw_pool_sub(i_frag)%pool)
    if (.not. pool_exists) then
      print *, "Error: pw_pool is not associated for fragment", i_frag
      return
    end if

    allocate(pw_real)

    ! Create a new pw_type in real space
    call fat_env%pw_pool_sub(i_frag)%pool%create_pw(pw_real, &
         use_data=realdata3d, in_space=realspace)

    if (.not. associated(pw_real)) then
      print *, "Error: pw_real was not created successfully"
      deallocate(pw_real)
      return
    end if
    call pw_zero(pw_real)

    ! Transfer data from reciprocal to real space
    call pw_transfer(pw, pw_real)

    call write_pw_r_to_cube(fat_env, pw_real, filename)

    ! Release the temporary real space pw_type
    call fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(pw_real)

    deallocate(pw_real)

  end subroutine write_pw_g_to_cube

  subroutine print_pw_statistics(pw, name, i_frag)
    type(pw_type), intent(in) :: pw
    character(len=*), intent(in) :: name
    integer, intent(in) :: i_frag

    real(dp) :: mean, std_dev, min_val, max_val, sum_val, norm_2
    integer :: count, i
    type(pw_grid_type), pointer :: grid
    logical :: print_grid = .true.

    grid => pw%pw_grid

    if (pw%in_space == REALSPACE) then
      sum_val = sum(abs(pw%cr3d))
      count = size(pw%cr3d)
      mean = sum_val / count
      std_dev = sqrt(sum((abs(pw%cr3d) - mean)**2) / count)
      min_val = minval(abs(pw%cr3d))
      max_val = maxval(abs(pw%cr3d))
      norm_2 = sqrt(sum(pw%cr3d**2))
    else if (pw%in_space == RECIPROCALSPACE) then
      sum_val = sum(abs(pw%cc))
      count = size(pw%cc)
      mean = real(sum_val) / count
      std_dev = sqrt(sum(abs(pw%cc - mean)**2) / count)
      min_val = minval(abs(pw%cc))
      max_val = maxval(abs(pw%cc))
      norm_2 = sqrt(sum(abs(pw%cc)**2))
    else
      print *, "KKK| Error: Unknown space type for PW"
      return
    end if

    ! PW Type Information
    print *, "KKK| Statistics for -> ", name, " <-(fragment ", i_frag, "):"
    print *, "KKK|  PW Type Information:"
    print *, "KKK|   In use:                ", pw%in_use
    print *, "KKK|   Space:                 ", merge("REALSPACE      ", &
      "RECIPROCALSPACE", pw%in_space == REALSPACE)
    print *, "KKK|  Data Statistics:"
    print *, "KKK|   Mean(abs):             ", mean
    print *, "KKK|   Std Dev:               ", std_dev
    print *, "KKK|   Min(abs):              ", min_val
    print *, "KKK|   Max(abs):              ", max_val
    print *, "KKK|   Sum(abs):              ", sum_val
    print *, "KKK|   Norm-2:                ", norm_2
    print *, "KKK|   Num of grid p.         ", count

    if (print_grid) then
      print *, "KKK|  Grid Information:"
      print *, "KKK|   Total Grid Points:     ", grid%ngpts
      print *, "KKK|   Grid Points in Cutoff: ", grid%ngpts_cut
      print *, "KKK|   Dimensions:            ", grid%npts
      print *, "KKK|   Local Grid Points:     ", grid%ngpts_local
      print *, "KKK|   Local Points in Cutoff:", grid%ngpts_cut_local
      print *, "KKK|   Grid Bounds:           ", grid%bounds
      print *, "KKK|   Local Grid Bounds:     ", grid%bounds_local
      print *, "KKK|   Grid Spacing:          ", grid%dr
      print *, "KKK|   Cell Matrix:"
      do i = 1, 3
        print *, "KKK|                          ", grid%dh(:,i)
      end do
      print *, "KKK|   Inverse Cell Matrix:"
      do i = 1, 3
        print *, "KKK|                          ", grid%dh_inv(:,i)
      end do
      print *, "KKK|   Orthorhombic:          ", grid%orthorhombic
      print *, "KKK|   Volume Element:        ", grid%dvol
      print *, "KKK|   Total Volume:          ", grid%vol
      print *, "KKK|   Cutoff:                ", grid%cutoff
      print *, "KKK|   Grid Span:             ", merge("HALFSPACE", &
        "FULLSPACE", grid%grid_span == HALFSPACE)
      print *, "KKK|   Has G=0:               ", grid%have_g0
      print *, "KKK|   First G≠0 index:       ", grid%first_gne0
      print *, "KKK|   Grid ID:               ", grid%id_nr
      print *, "KKK|   Reference Count:       ", grid%ref_count
      print *, "KKK|   Spherical Cutoff:      ", grid%spherical
      print *, "KKK|  Parallel Information:"
      print *, "KKK|   Mode:                  ", merge("Distributed", &
        "Local      ", grid%para%mode == PW_MODE_DISTRIBUTED)
      print *, "KKK|   Ray Distribution:      ", grid%para%ray_distribution
      print *, "KKK|   Blocked:               ", grid%para%blocked
      print *, "KKK|   Group Size:            ", grid%para%group_size
      print *, "KKK|   My Position:           ", grid%para%my_pos
      print *, "KKK|   Group Head:            ", grid%para%group_head
      print *, "KKK|   Group Head ID:         ", grid%para%group_head_id
      print *, "KKK|   Real Space Dimensions: ", grid%para%rs_dims
      print *, "KKK|   Real Space Position:   ", grid%para%rs_pos
    end if

    if (associated(grid%g)) then
      print *, "KKK|  G-vector Statistics:"
      print *, "KKK|   Mean G magnitude:      ", sum(sqrt(sum(grid%g**2, dim=1))) / size(grid%g, 2)
      print *, "KKK|   Max G magnitude:       ", sqrt(maxval(sum(grid%g**2, dim=1)))
    end if

    if (associated(grid%gsq)) then
      print *, "KKK|   Mean G^2:              ", sum(grid%gsq) / size(grid%gsq)
      print *, "KKK|   Max G^2:               ", maxval(grid%gsq)
    end if

!#if defined(__parallel)
!    print *, "KKK|  MPI Information:"
!    print *, "KKK|   MPI Mode:               ", merge("DISTRIBUTED","LOCAL      ", &
!             grid%para%mode == PW_MODE_DISTRIBUTED) 
!    print *, "KKK|   MPI Group Handle:       ", grid%para%group%get_handle()
!    print *, "KKK|   MPI Group Size:         ", grid%para%group_size
!    print *, "KKK|   MPI Position:           ", grid%para%my_pos
!    print *, "KKK|   MPI Dims:               ", grid%para%rs_dims
!    print *, "KKK|   MPI Ray Distribution:   ", grid%para%ray_distribution
!    print *, "KKK|   MPI Blocked:            ", grid%para%blocked
!    print *, "KKK|   MPI Task Rank:          ", grid%para%rs_dims
!    print *, "KKK|   Real Space Group:       ", grid%para%rs_group%get_handle()
!    print *, "KKK|   Real Space Position:    ", grid%para%rs_pos
!    print *, "KKK|   Real Space MPO:         ", grid%para%rs_mpo
!    print *, "KKK|   Has Bounds Array:       ", ALLOCATED(grid%para%bo)
!    if (ALLOCATED(grid%para%bo)) then
!      print *, "KKK|   Bounds Array Shape:     ", SHAPE(grid%para%bo)
!    endif
!    print *, "KKK|   Has YZ Ray Array:       ", ALLOCATED(grid%para%yzp)
!    print *, "KKK|   Has YZ Count Array:     ", ALLOCATED(grid%para%nyzray)
!#endif
  end subroutine print_pw_statistics
! **************************************************************************************************
!> \brief Gets the nuclear density for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_rhon_frag(fat_env, rhon, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: rhon
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_rhon_frag"

    type(qs_environment_type), pointer :: qs_env

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env=qs_env, rho_core=rhon)

  end subroutine fat_rhon_frag

! **************************************************************************************************
!> \brief Gets the gradient of the nuclear density w.r.t. a nuclear position for one fragment in the
!>        freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_drhon_frag(fat_env, drhon, i_frag, i_nuc, i_dir)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: drhon
    integer :: i_frag, i_nuc, i_dir

    character(len=*), parameter :: routinen = "fat_rhon_frag"

    type(qs_environment_type), pointer :: qs_env

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call calculate_drho_core(drhon, qs_env, i_dir, i_nuc)

  end subroutine fat_drhon_frag

! **************************************************************************************************
!> \brief Gets the gradient of the electronic density w.r.t. a nuclear position for one fragment in
!>        the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_drho_frag(fat_env, drho, i_frag, i_nuc, i_dir)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: drho
    integer :: i_frag, i_nuc, i_dir

    character(len=*), parameter :: routinen = "fat_drho_frag"

    type(dbcsr_p_type), pointer :: p_mat(:)
    type(pw_type), pointer :: drho_g
    type(qs_environment_type), pointer :: qs_env
    type(qs_rho_type), pointer :: rho

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, rho=rho)
    call qs_rho_get(rho, rho_ao=p_mat)
    call fat_pw_create_frag(fat_env, drho_g, i_frag, .false.)
    call calculate_drho_elec_dR( &
         p_mat(1)%matrix, drho=drho, drho_gspace=drho_g, qs_env=qs_env, beta=i_dir, lambda=i_nuc)
    call fat_pw_release_frag(fat_env, drho_g, i_frag)

  end subroutine fat_drho_frag

! **************************************************************************************************
!> \brief Gets the total electron density for the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_rho_ref(fat_env, pw)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pw(:)

    character(len=*), parameter :: routinen = "fat_rho_ref"

    integer :: i, j
    type(pw_type), pointer :: pw_tmp(:)

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    do i = 1, fat_env%n_frag
       call fat_rho_frag(fat_env, pw_tmp, i)
       do j = 1, size(pw_tmp)
          call pw_axpy(pw_tmp(j), pw(j), allow_noncompatible_grids=.true.)
       end do
    end do

  end subroutine fat_rho_ref

! **************************************************************************************************
!> \brief Calculate kinetic potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_ts_frag(fat_env, v, v_gga, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: v(:)
    type(pw_type), pointer :: v_gga(:)
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_v_ts_frag"

    real(dp) :: e, virial(3, 3)
    type(pw_type), pointer :: rho_g(:), rho_r(:), tau(:), v_tau(:)
    type(qs_environment_type), pointer :: qs_env
    type(section_vals_type), pointer :: input, kin_section
    type(xc_rho_cflags_type) :: needs
    logical :: use_gga = .false.
    integer :: i

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    ! Get XC section for kinetic
    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, input=input)
    kin_section => section_vals_get_subs_vals(input, "DFT%QS%OPT_EMBED%XC")
    needs = xc_functionals_get_needs( &
              section_vals_get_subs_vals(kin_section, "XC_FUNCTIONAL"), &
                                         lsd=.false., calc_potential=.true.)
    use_gga = needs%norm_drho .OR. needs%norm_drho_spin

    ! Get density
    call fat_rho_frag(fat_env, rho_r, i_frag)

    ! Initialize rho_g if needed for GGA
    if (use_gga) then
       allocate(rho_g(size(rho_r)))
       do i = 1, size(rho_r)
          call fat_env%pw_pool_sub(i_frag)%pool%create_pw(rho_g(i), &
                        use_data=complexdata1d, in_space=reciprocalspace)
          call pw_transfer(rho_r(i), rho_g(i))
       end do
    end if

    ! Calculate LDA part
    nullify(v, v_tau)
    call xc_vxc_pw_create(vxc_rho=v, vxc_tau=v_tau, rho_r=rho_r, rho_g=rho_g, &
                          tau=tau, exc=e, xc_section=kin_section, &
                          pw_pool=fat_env%pw_pool_sub(i_frag)%pool, &
                          compute_virial=.false., virial_xc=virial)
    fat_env%ekin_sub(i_frag) = e


    ! Calculate GGA part if needed
    if (use_gga) then
       call xc_vxc_pw_create(vxc_rho=v_gga, vxc_tau=v_tau, &
                             rho_r=rho_r, rho_g=rho_g, &
                             tau=tau, exc=e, xc_section=kin_section, &
                             pw_pool=fat_env%pw_pool_sub(i_frag)%pool, &
                             compute_virial=.false., virial_xc=virial)
       fat_env%ekin_sub(i_frag) = fat_env%ekin_sub(i_frag) + e

       ! Cleanup rho_g
       do i = 1, size(rho_g)
           call fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(rho_g(i))
       end do
       deallocate(rho_g)
    end if

  end subroutine fat_v_ts_frag

! **************************************************************************************************
!> \brief Calculate reference kinetic potential for the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_ts_ref(fat_env, v, v_gga)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: v(:)
    type(pw_type), pointer :: v_gga(:)

    character(len=*), parameter :: routinen = "fat_v_ts_ref"

    real(dp) :: e, virial(3, 3)
    type(pw_type), pointer :: rho_g(:), rho_r(:), tau(:), v_tau(:)
    type(qs_environment_type), pointer :: qs_env
    type(section_vals_type), pointer :: input, kin_section
    type(xc_rho_cflags_type) :: needs
    logical :: use_gga = .false.
    integer :: i

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    ! Get XC section for Kinetic
    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, input=input)
    kin_section => section_vals_get_subs_vals(input, "DFT%QS%OPT_EMBED%XC")
    needs = xc_functionals_get_needs( &
              section_vals_get_subs_vals(kin_section, "XC_FUNCTIONAL"), &
                                         lsd=.false., calc_potential=.true.)
    use_gga = needs%norm_drho .OR. needs%norm_drho_spin

    ! Create pw for total density
    call fat_pws_create_ref(fat_env, rho_r, 1)
    call fat_rho_ref(fat_env, rho_r)

    ! Initialize rho_g if needed for GGA
    if (use_gga) then
       allocate(rho_g(size(rho_r)))
       do i = 1, size(rho_r)
          call fat_env%pw_pool_ref%create_pw(rho_g(i), &
                     use_data=complexdata1d, in_space=reciprocalspace)
          call pw_transfer(rho_r(i), rho_g(i))
       end do
    end if

    ! Calculate LDA part
    nullify(v, v_tau)
    call xc_vxc_pw_create(vxc_rho=v, vxc_tau=v_tau, rho_r=rho_r, rho_g=rho_g, &
                          tau=tau, exc=e, xc_section=kin_section, &
                          pw_pool=fat_env%pw_pool_ref, &
                          compute_virial=.false., virial_xc=virial)
    fat_env%ekin_ref = e


    ! Calculate GGA part
    if (use_gga) then
        call xc_vxc_pw_create(vxc_rho=v_gga, vxc_tau=v_tau, &
                              rho_r=rho_r, rho_g=rho_g, &
                              tau=tau, exc=e, xc_section=kin_section, &
                              pw_pool=fat_env%pw_pool_ref, &
                              compute_virial=.false., virial_xc=virial)
        fat_env%ekin_ref = fat_env%ekin_ref + e

        ! Cleanup rho_g
        do i = 1, size(rho_g)
            call fat_env%pw_pool_ref%give_back_pw(rho_g(i))
        end do
        deallocate(rho_g)
    end if

    call fat_pws_release_ref(fat_env, rho_r)

  end subroutine fat_v_ts_ref

! **************************************************************************************************
!> \brief Calculate xc potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_xc_frag(fat_env, v, v_gga, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: v(:)
    type(pw_type), pointer :: v_gga(:)
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_v_xc_frag"

    real(dp) :: e, virial(3, 3)
    type(pw_env_type), pointer :: pw_env
    type(pw_type), pointer :: rho_g(:), rho_r(:), tau(:), v_tau(:)
    type(pw_pool_type), pointer :: pw_pool
    type(qs_environment_type), pointer :: qs_env
    type(section_vals_type), pointer :: input, xc_section
    type(xc_rho_cflags_type) :: needs
    logical :: use_gga = .false.
    integer :: i

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    ! Get XC section
    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, input=input)
    xc_section => section_vals_get_subs_vals(input, "DFT%XC")
    needs = xc_functionals_get_needs( &
            section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL"), &
            lsd=.false., calc_potential=.true.)
    use_gga = needs%norm_drho .OR. needs%norm_drho_spin

    call get_qs_env(qs_env, pw_env=pw_env)
    call pw_env_get(pw_env, auxbas_pw_pool=pw_pool)

    ! Get density
    call fat_rho_frag(fat_env, rho_r, i_frag)

    ! Initialize rho_g if needed for GGA
    nullify(rho_g, tau)
    if (use_gga) then
       allocate(rho_g(size(rho_r)))
       do i = 1, size(rho_r)
          call fat_env%pw_pool_sub(i_frag)%pool%create_pw(rho_g(i), &
                                use_data=complexdata1d, in_space=reciprocalspace)
          call pw_transfer(rho_r(i), rho_g(i))
       end do
    end if

    ! Calculate LDA part
    nullify(v, v_tau)
    call xc_vxc_pw_create(vxc_rho=v, vxc_tau=v_tau, rho_r=rho_r, rho_g=rho_g, &
                          tau=tau, exc=e, xc_section=xc_section, &
                          pw_pool=fat_env%pw_pool_sub(i_frag)%pool, &
                          compute_virial=.false., virial_xc=virial)
    fat_env%exc_sub(i_frag) = e


    ! Calculate GGA part if needed
    if (use_gga) then
        call xc_vxc_pw_create(vxc_rho=v_gga, vxc_tau=v_tau, &
                            rho_r=rho_r, rho_g=rho_g, &
                            tau=tau, exc=e, xc_section=xc_section, &
                            pw_pool=fat_env%pw_pool_sub(i_frag)%pool, &
                            compute_virial=.false., virial_xc=virial)
        fat_env%exc_sub(i_frag) = fat_env%exc_sub(i_frag) + e

        ! Cleanup rho_g
        do i = 1, size(rho_g)
            call fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(rho_g(i))
        end do
        deallocate(rho_g)
    end if

  end subroutine fat_v_xc_frag

! **************************************************************************************************
!> \brief Calculate reference xc potential for the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_xc_ref(fat_env, v, v_gga)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: v(:)
    type(pw_type), pointer :: v_gga(:)

    character(len=*), parameter :: routinen = "fat_v_xc_ref"

    real(dp) :: e, virial(3, 3)
    type(pw_env_type), pointer :: pw_env
    type(pw_type), pointer :: rho_g(:), rho_r(:), tau(:), v_tau(:)
    type(pw_pool_type), pointer :: pw_pool
    type(qs_environment_type), pointer :: qs_env
    type(section_vals_type), pointer :: input, xc_section
    type(xc_rho_cflags_type) :: needs
    logical :: use_gga = .false.
    integer :: i

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    ! Get XC section for XC
    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, input=input)
    xc_section => section_vals_get_subs_vals(input, "DFT%XC")
    needs = xc_functionals_get_needs( &
            section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL"), &
            lsd=.false., calc_potential=.true.)
    use_gga = needs%norm_drho .OR. needs%norm_drho_spin

    ! Get PW environment and pool
    call get_qs_env(qs_env, pw_env=pw_env)
    call pw_env_get(pw_env, auxbas_pw_pool=pw_pool)

    ! Create and initialize density arrays
    call fat_pws_create_ref(fat_env, rho_r, 1)
    call fat_rho_ref(fat_env, rho_r)

    ! Initialize rho_g if needed for GGA
    nullify(rho_g, tau)
    if (use_gga) then
       allocate(rho_g(size(rho_r)))
       do i = 1, size(rho_r)
          call fat_env%pw_pool_ref%create_pw(rho_g(i), &
                      use_data=complexdata1d, in_space=reciprocalspace)
          call pw_transfer(rho_r(i), rho_g(i))
       end do
    end if

    ! Calculate LDA part
    nullify(v, v_tau)

    ! Calculate XC potential and energy
    call xc_vxc_pw_create(vxc_rho=v, vxc_tau=v_tau, rho_r=rho_r, rho_g=rho_g, &
                          tau=tau, exc=e, xc_section=xc_section, &
                          pw_pool=pw_pool, compute_virial=.false., virial_xc=virial)

    ! Store reference XC energy
    fat_env%exc_ref = e

    ! Calculate GGA part
    if (use_gga) then
       call xc_vxc_pw_create(vxc_rho=v_gga, vxc_tau=v_tau, &
                             rho_r=rho_r, rho_g=rho_g, &
                             tau=tau, exc=e, xc_section=xc_section, &
                             pw_pool=fat_env%pw_pool_ref, &
                             compute_virial=.false., virial_xc=virial)
       fat_env%exc_ref = fat_env%exc_ref + e

       ! Cleanup rho_g
       do i = 1, size(rho_g)
           call fat_env%pw_pool_ref%give_back_pw(rho_g(i))
       end do
       deallocate(rho_g)
    end if

    ! Cleanup
    call fat_pws_release_ref(fat_env, rho_r)

  end subroutine fat_v_xc_ref

! **************************************************************************************************
!> \brief Calculate electronic potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_e_frag(fat_env, v, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: v
    integer :: i_frag

    character(len=*), parameter :: routineN = "fat_v_e_frag"

    real(dp) :: e
    type(pw_env_type), pointer :: pw_env
    type(pw_poisson_type), pointer :: poisson_env
    type(pw_type), pointer :: rho_g(:), v_g(:)
    type(qs_environment_type), pointer :: qs_env

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    nullify(rho_g)
    call fat_rho_frag(fat_env, rho_g, i_frag, .false.)

!   if (is_debug) then
!      call print_pw_statistics(rho_g(1), "rho", i_frag)
!      call write_pw_g_to_cube(fat_env, rho_g(1), "rho", i_frag)
!   end if

    allocate(v_g(1))
    call fat_env%pw_pool_sub(i_frag)%pool%create_pw(v_g(1), use_data=complexdata1d, in_space=reciprocalspace)
    call pw_zero(v_g(1))

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, pw_env=pw_env)
    call pw_env_get(pw_env, poisson_env=poisson_env)
    call pw_poisson_solve(poisson_env, rho_g(1), e, v_g(1))

    if (is_debug) &
       print *, "DBG| Electron-Electron Functional Energy (", i_frag, "): ", e
!   if (is_debug) call print_pw_statistics(v_g(1), "v_ee", i_frag)

    call pw_transfer(v_g(1), v)

    call fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(v_g(1))
    deallocate(v_g)

  end subroutine fat_v_e_frag
! **************************************************************************************************
!> \brief Calculate nuclear potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_n_frag(fat_env, v, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: v
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_v_n_frag"

    real(dp) :: e
    type(pw_env_type), pointer :: pw_env
    type(pw_poisson_type), pointer :: poisson_env
    type(pw_type), pointer :: rho_g(:), v_g(:), rho_core
    type(qs_environment_type), pointer :: qs_env

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    allocate(rho_g(1))
    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, rho_core=rho_core)
    rho_g(1) = rho_core

!   if (is_debug) call print_pw_statistics(rho_g(1), "rho_core", i_frag)
!   if (is_debug) call write_pw_g_to_cube(fat_env, rho_g(1), "rho_core", i_frag)

    allocate(v_g(1))
    call fat_env%pw_pool_sub(i_frag)%pool%create_pw(v_g(1), &
            use_data=complexdata1d, in_space=reciprocalspace)
    call pw_zero(v_g(1))

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, pw_env=pw_env)
    call pw_env_get(pw_env, poisson_env=poisson_env)
    call pw_poisson_solve(poisson_env, rho_g(1), e, v_g(1))

    if (is_debug) &
       print *, "DBG| Nuclear-Electron Functional Energy (", i_frag, "): ", e

    call pw_transfer(v_g(1), v)
    call fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(v_g(1))
    deallocate(rho_g)
    deallocate(v_g)

  end subroutine fat_v_n_frag

! **************************************************************************************************
!> \brief Calculate non-additive kinetic embedding potential for one fragment in the freeze-and-thaw
!>        environment
! **************************************************************************************************
  subroutine fat_vemb_tsna_frag(fat_env, vemb, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: vemb
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_vemb_tsna_frag"

    type(pw_type), pointer :: v(:), v_gga(:)

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call pw_zero(vemb)

    ! Reference system
    nullify(v, v_gga)
    call fat_v_ts_ref(fat_env, v, v_gga)
!   if (is_debug) call print_pw_statistics(v(1), "v_emb Kin tot", 0)
    call pw_axpy(v(1), vemb, allow_noncompatible_grids=.true.)
    call fat_pws_release_ref(fat_env, v)
    IF (ASSOCIATED(v_gga)) then
       call pw_axpy(v_gga(1), vemb, allow_noncompatible_grids=.true.)
       call fat_pws_release_ref(fat_env, v_gga)
    end if

    ! Fragment
    nullify(v, v_gga)
    call fat_v_ts_frag(fat_env, v, v_gga, i_frag)
!   if (is_debug) call print_pw_statistics(v(1), "v_emb Kin", i_frag)
    call pw_axpy(v(1), vemb, -1.0_dp, allow_noncompatible_grids=.true.)
    call fat_pws_release_frag(fat_env, v, i_frag)
    if (associated(v_gga)) then
       call pw_axpy(v_gga(1), vemb, -1.0_dp, allow_noncompatible_grids=.true.)
       call fat_pws_release_frag(fat_env, v_gga, i_frag)
    end if

  end subroutine fat_vemb_tsna_frag

! **************************************************************************************************
!> \brief Calculate non-additive xc embedding potential for one fragment in the freeze-and-thaw
!>        environment
! **************************************************************************************************
  subroutine fat_vemb_xcna_frag(fat_env, vemb, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: vemb
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_vemb_xcna_frag"

    type(pw_type), pointer :: v(:), v_gga(:)

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call pw_zero(vemb)

    ! Reference system
    nullify(v, v_gga)
    call fat_v_xc_ref(fat_env, v, v_gga)
!   if (is_debug) call print_pw_statistics(v(1), "v_emb XC tot", 0)
    call pw_axpy(v(1), vemb, allow_noncompatible_grids=.true.)
    call fat_pws_release_ref(fat_env, v)
    IF (ASSOCIATED(v_gga)) then
       call pw_axpy(v_gga(1), vemb, allow_noncompatible_grids=.true.)
       call fat_pws_release_ref(fat_env, v_gga)
    end if

    ! Fragment
    nullify(v, v_gga)
    call fat_v_xc_frag(fat_env, v, v_gga, i_frag)
!   if (is_debug) then
!     call print_pw_statistics(v(1), "v_emb XC LDA", i_frag)
!     if (associated(v_gga)) then
!       call print_pw_statistics(v_gga(1), "v_emb XC GGA", i_frag)
!     endif
!   endif
    call pw_axpy(v(1), vemb, -1.0_dp, allow_noncompatible_grids=.true.)
    call fat_pws_release_frag(fat_env, v, i_frag)
    if (associated(v_gga)) then
       call pw_axpy(v_gga(1), vemb, -1.0_dp, allow_noncompatible_grids=.true.)
       call fat_pws_release_frag(fat_env, v_gga, i_frag)
    end if

  end subroutine fat_vemb_xcna_frag

! **************************************************************************************************
!> \brief Calculate electronic embedding potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_vemb_e_frag(fat_env, vemb, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: vemb
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_vemb_e_frag"

    integer :: i
    type(pw_type), pointer :: v, v_sum
    real(kind=dp) :: integral_value

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call pw_zero(vemb)

    allocate(v_sum)
    call fat_env%pw_pool_sub(i_frag)%pool%create_pw(v_sum, use_data=realdata3d, in_space=realspace)
    call pw_zero(v_sum)

    do i = 1, fat_env%n_frag
       if (i == i_frag) cycle

       call fat_pw_create_frag(fat_env, v, i)
       call fat_v_e_frag(fat_env, v, i)

       integral_value = pw_integral_ab(v, v)
       if (is_debug) print *, "KKK| Integral of vemb_ee frag", i, ": ", integral_value
!      if (is_debug) call print_pw_statistics(v, "vemb_ee", i_frag)
!      if (is_debug) call write_pw_r_to_cube(fat_env, v, "vemb_ee")

       if (is_debug) print *, "KKK| vemb_ee allow_noncompatible_grids=True"
       call pw_axpy(v, v_sum, allow_noncompatible_grids=.true.)
       call fat_env%pw_pool_sub(i)%pool%give_back_pw(v)
       call fat_pw_release_frag(fat_env, v, i)
    end do

    call pw_transfer(v_sum, vemb)

    call fat_pw_release_frag(fat_env, v_sum, i_frag)

  end subroutine fat_vemb_e_frag
! **************************************************************************************************
!> \brief Calculate nuclear embedding potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_vemb_n_frag(fat_env, vemb, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: vemb
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_vemb_n_frag"

    integer :: i
    type(pw_type), pointer :: v, v_sum
    real(kind=dp) :: integral_value

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call pw_zero(vemb)

    allocate(v_sum)
    call fat_env%pw_pool_sub(i_frag)%pool%create_pw(v_sum, use_data=realdata3d, in_space=realspace)
    call pw_zero(v_sum)

    do i = 1, fat_env%n_frag
       if (i == i_frag) cycle

       call fat_pw_create_frag(fat_env, v, i)
       call fat_v_n_frag(fat_env, v, i)

       integral_value = pw_integral_ab(v, v)
       if (is_debug) print *, "KKK| Integral of vemb_ne frag", i, ": ", integral_value
!      if (is_debug) call print_pw_statistics(v, "vemb_ne", i_frag)
!      if (is_debug) call write_pw_r_to_cube(fat_env, v, "vemb_ne")

       if (is_debug) print *, "KKK| vemb_ne allow_noncompatible_grids=True"
       call pw_axpy(v, v_sum, allow_noncompatible_grids=.true.)
       call fat_env%pw_pool_sub(i)%pool%give_back_pw(v)
       call fat_pw_release_frag(fat_env, v, i)
    end do

    call pw_transfer(v_sum, vemb)

    call fat_pw_release_frag(fat_env, v_sum, i_frag)

  end subroutine fat_vemb_n_frag

! **************************************************************************************************
!> \brief Calculate embedding potential for one fragment in the freeze-and-thaw environment
!> \param fat_env The freeze-and-thaw environment
!> \param vemb The embedding potential (output)
!> \param i_frag The index of the current fragment
!> \param to_file Boolean flag to determine if the potential should be written to a file
! **************************************************************************************************
  subroutine fat_vemb_frag(fat_env, vemb, i_frag, to_file)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: vemb
    integer :: i_frag
    logical :: to_file

    character(len=*), parameter :: routinen = "fat_vemb_frag"

    character(len=default_string_length) :: label
    integer :: output_unit
    integer, pointer :: stride(:)
    type(particle_list_type), pointer :: particle_list
    type(pw_type), pointer :: &
      rho(:), vemb_tsna, vemb_xcna, vemb_ee, vemb_ne, weights
    type(qs_environment_type), pointer :: qs_env
    type(qs_subsys_type), pointer :: subsys

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    ! Get the electron density for the current fragment
    call fat_rho_frag(fat_env, rho, i_frag)
!   if (is_debug) call print_pw_statistics(rho(1), "rho", i_frag)

    ! Calculate and add the non-additive kinetic embedding potential
    call fat_pw_create_frag(fat_env, vemb_tsna, i_frag)
    call fat_vemb_tsna_frag(fat_env, vemb_tsna, i_frag)
!   if(is_debug) call print_pw_statistics(vemb_tsna, "v_emb Kin", i_frag)
    call pw_axpy(vemb_tsna, vemb)
    call fat_pw_release_frag(fat_env, vemb_tsna, i_frag)

    ! Calculate and add the non-additive exchange-correlation embedding potential
    call fat_pw_create_frag(fat_env, vemb_xcna, i_frag)
    call fat_vemb_xcna_frag(fat_env, vemb_xcna, i_frag)
!   if(is_debug) call print_pw_statistics(vemb_xcna, "v_emb XC", i_frag)
    call pw_axpy(vemb_xcna, vemb)
    call fat_pw_release_frag(fat_env, vemb_xcna, i_frag)

!   ! TODO: Calculate and add the non-additive exact exchange HF embedding potential
!   IF (fat_env(i_frag)%use_hfx) THEN
!      call fat_pw_create_frag(fat_env, vemb_hfx, i_frag)
!      call fat_vemb_hfx_frag(fat_env, vemb_hfx, i_frag)
!      if(is_debug) call print_pw_statistics(vemb_xcna, "v_emb XC", i_frag)
!      call pw_axpy(vemb_hfx, vemb)
!      call fat_pw_release_frag(fat_env, vemb_hfx, i_frag)
!   END IF

    ! Calculate and add the electron-electron embedding potential
    call fat_pw_create_frag(fat_env, vemb_ee, i_frag)
    call fat_vemb_e_frag(fat_env, vemb_ee, i_frag)
!   if(is_debug) call print_pw_statistics(vemb_ee, "v_emb EE", i_frag)
    call pw_axpy(vemb_ee, vemb)
    call fat_pw_release_frag(fat_env, vemb_ee, i_frag)

    ! Calculate and add the nuclear-electron embedding potential
    call fat_pw_create_frag(fat_env, vemb_ne, i_frag)
    call fat_vemb_n_frag(fat_env, vemb_ne, i_frag)
!   if(is_debug) call print_pw_statistics(vemb_ne, "v_emb NE", i_frag)
    call pw_axpy(vemb_ne, vemb)
    call fat_pw_release_frag(fat_env, vemb_ne, i_frag)

    ! If requested, write the embedding potential to files
    if (to_file) then
      IF (fat_env%force_env_ref%para_env%is_source()) THEN
        ! Convert fragment index to string for file naming
        call integer_to_string(i_frag, label)

        ! Write embedding potential to a cube file
        if (is_debug) print *, "DBG| ", modulen//"::"//routinen

        call open_file( &
             "vemb_"//trim(label)//".cube", file_status="replace",&
             & file_action="write", unit_number=output_unit)
        qs_env => fat_env%force_env_ref%qs_env
        call get_qs_env(qs_env, subsys=subsys)
        call qs_subsys_get(subsys, particles=particle_list)
        allocate(stride(3))
        stride = [1, 1, 1]
        call cp_pw_to_cube(vemb, output_unit, "vemb_"//trim(label), particle_list, stride=stride)
        deallocate(stride)

        ! Write embedding potential to a dat file (simple volumetric format)
        call open_file( &
             "vemb_"//trim(label)//".dat", file_status="REPLACE", &
             file_action="WRITE", unit_number=output_unit)
        call fat_pw_create_frag(fat_env, weights, i_frag)
        ! Set weights to the volume element of the grid
        call pw_set(weights, vemb%pw_grid%dr(1) * vemb%pw_grid%dr(2) * vemb%pw_grid%dr(3))
        ! Write the potential values along with their weights
        call pw_to_simple_volumetric(weights, output_unit, [1, 1, 1], vemb)
        call fat_pw_release_frag(fat_env, weights, i_frag)
      end if
    end if

  end subroutine fat_vemb_frag

! **************************************************************************************************
!> \brief Calculate the kinetic energy for on fragment of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_e_ts_frag(fat_env, i_frag) result(e)
    type(fat_env_type), pointer :: fat_env
    integer :: i_frag
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_e_ts_frag"

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    e = fat_env%ekin_sub(i_frag)

  end function fat_e_ts_frag

! **************************************************************************************************
!> \brief Calculate the kinetic energy for the reference system of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_e_ts_ref(fat_env) result(e)
    type(fat_env_type), pointer :: fat_env
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_e_ts_ref"

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    e = fat_env%ekin_ref

  end function fat_e_ts_ref

! **************************************************************************************************
!> \brief Calculate the xc energy for on fragment of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_e_xc_frag(fat_env, i_frag) result(e)
    type(fat_env_type), pointer :: fat_env
    integer :: i_frag
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_e_xc_frag"

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    e = fat_env%exc_sub(i_frag)

  end function fat_e_xc_frag

! **************************************************************************************************
!> \brief Calculate the xc energy for the reference system of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_e_xc_ref(fat_env) result(e)
    type(fat_env_type), pointer :: fat_env
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_e_xc_ref"

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    e = fat_env%exc_ref

  end function fat_e_xc_ref

! **************************************************************************************************
!> \brief Calculate the ee interaction energy for two fragments of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_eint_ee(fat_env, i, j) result(e)
    type(fat_env_type), pointer :: fat_env
    integer :: i, j
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_eint_ee"

    type(pw_type), pointer :: v_ee, rho, rho_tmp(:)

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call fat_pw_create_frag(fat_env, v_ee, i)
    call fat_v_e_frag(fat_env, v_ee, i)

    call fat_pw_create_frag(fat_env, rho, i)
    call fat_rho_frag(fat_env, rho_tmp, j)

    call pw_transfer(rho_tmp(1), rho)
    e = pw_integral_ab(v_ee, rho)

    call fat_pw_release_frag(fat_env, v_ee, i)
    call fat_pw_release_frag(fat_env, rho, i)

  end function fat_eint_ee

! **************************************************************************************************
!> \brief Calculate the ne interaction energy for two fragments of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_eint_ne(fat_env, i, j) result(e)
    type(fat_env_type), pointer :: fat_env
    integer :: i, j
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_eint_ne"

    type(pw_type), pointer :: v, rho, rho_tmp(:)

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call fat_pw_create_frag(fat_env, v, i)
    call fat_v_n_frag(fat_env, v, i)

    call fat_pw_create_frag(fat_env, rho, i)
    call fat_rho_frag(fat_env, rho_tmp, j)

    call pw_transfer(rho_tmp(1), rho)
    e = pw_integral_ab(v, rho)

    call fat_pw_release_frag(fat_env, v, i)
    call fat_pw_release_frag(fat_env, rho, i)

  end function fat_eint_ne

! **************************************************************************************************
!> \brief Calculate the nn interaction energy for two fragments of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_eint_nn(fat_env, i, j) result(e)
    type(fat_env_type), pointer :: fat_env
    integer :: i, j
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_eint_nn"

    type(pw_type), pointer :: v, rho, rho_tmp

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call fat_pw_create_frag(fat_env, v, i)
    call fat_v_n_frag(fat_env, v, i)

    call fat_pw_create_frag(fat_env, rho, i)
    call fat_rhon_frag(fat_env, rho_tmp, j)
!   if (is_debug) call print_pw_statistics(rho, "rho", i)
!   if (is_debug) call print_pw_statistics(rho_tmp, "rho_tmp", j)

    call pw_transfer(rho_tmp, rho)
    e = pw_integral_ab(v, rho)

    call fat_pw_release_frag(fat_env, v, i)
    call fat_pw_release_frag(fat_env, rho, i)

  end function fat_eint_nn

! **************************************************************************************************
!> \brief Calculate the total energy of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_e_tot(fat_env) result(e)
    type(fat_env_type), pointer :: fat_env
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_e_tot"

    integer :: i, j
    real(dp) :: e_tmp, e_sub_total, e_emb_total, e_ts_total, e_xc_total
    real(dp) :: e_ee_total, e_ne_total, e_nn_total

    if (is_debug) print *, "KKK| ", modulen//"::"//routinen

    e = 0.0_dp
    e_sub_total = 0.0_dp
    e_emb_total = 0.0_dp
    e_ts_total = 0.0_dp
    e_xc_total = 0.0_dp
    e_ee_total = 0.0_dp
    e_ne_total = 0.0_dp
    e_nn_total = 0.0_dp

    ! Sum of subsystem energies
    do i = 1, fat_env%n_frag
       e_sub_total = e_sub_total + fat_env%e_sub(i)
       if (is_debug) print *, "KKK| Subsystem energy (", i, "):", fat_env%e_sub(i)
    end do
    e = e + e_sub_total
    if (is_debug) print *, "KKK| Total subsystem energy:", e_sub_total

    ! Embedding potential energy
    do i = 1, fat_env%n_frag
       e_emb_total = e_emb_total + fat_env%e_emb(i)
       if (is_debug) print *, "KKK| Embedding energy (", i, "):", fat_env%e_emb(i)
    end do
    e = e - e_emb_total
    if (is_debug) print *, "KKK| Total Embedding energy:", e_emb_total

    ! Kinetic energy
    e_ts_total = fat_e_ts_ref(fat_env)
    if (is_debug) print *, "KKK| Reference Kin energy:", e_ts_total
    do i = 1, fat_env%n_frag
       e_tmp = fat_e_ts_frag(fat_env, i)
       e_ts_total = e_ts_total - e_tmp
       if (is_debug) print *, "KKK| Subsystem Kin energy (", i, "):", e_tmp
    end do
    e = e + e_ts_total
    if (is_debug) print *, "KKK| Total Kin energy:", e_ts_total

    ! Exchange-correlation energy
    e_xc_total = fat_e_xc_ref(fat_env)
    if (is_debug) print *, "KKK| Reference XC energy:", e_xc_total
    do i = 1, fat_env%n_frag
       e_tmp = fat_e_xc_frag(fat_env, i)
       e_xc_total = e_xc_total - e_tmp
       if (is_debug) print *, "KKK| Subsystem XC energy (", i, "):", e_tmp
    end do
    e = e + e_xc_total
    if (is_debug) print *, "KKK| Total XC energy:", e_xc_total

    if (fat_env%i_iter > 0) then
       ! Electron-electron interaction
       do i = 1, fat_env%n_frag
          do j = 1, fat_env%n_frag
             if (i == j) cycle
             e_tmp = fat_eint_ee(fat_env, i, j)
             e_ee_total = e_ee_total + 0.5_dp * e_tmp
             if (is_debug) print *, "KKK| EE interaction energy (", i, ",", j, "):", e_tmp
          end do
       end do
       e = e + e_ee_total
       fat_env%eee = e_ee_total
       if (is_debug) print *, "KKK| Total EE interaction energy:", e_ee_total

       ! Nuclear-electron interaction
       do i = 1, fat_env%n_frag
          do j = 1, fat_env%n_frag
             if (i == j) cycle
             e_tmp = fat_eint_ne(fat_env, i, j)
             e_ne_total = e_ne_total + e_tmp
             if (is_debug) print *, "KKK| NE interaction energy (", i, ",", j, "):", e_tmp
          end do
       end do
       e = e + e_ne_total
       fat_env%ene = e_ne_total
       if (is_debug) print *, "KKK| Total NE interaction energy:", e_ne_total

       ! Nuclear-nuclear interaction
       do i = 1, fat_env%n_frag
          do j = 1, fat_env%n_frag
             if (i == j) cycle
             e_tmp = fat_eint_nn(fat_env, i, j)
             e_nn_total = e_nn_total + 0.5_dp * e_tmp
             if (is_debug) print *, "KKK| NN interaction energy (", i, ",", j, "):", e_tmp
          end do
       end do
       e = e + e_nn_total
       fat_env%enn = e_nn_total
       if (is_debug) print *, "KKK| Total NN interaction energy:", e_nn_total
    end if

    fat_env%e_tot = e
    if (is_debug) print *, "KKK| Total energy:", e

  end function fat_e_tot
! ************************************************************************************************
!> \brief Print the intermediate results of a freeze-and-thaw iteration
! **************************************************************************************************
  subroutine fat_status(fat_env, i_frag)
    type(fat_env_type), pointer :: fat_env
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_status"

    LOGICAL :: ionode
    real(dp) :: e_tot

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    ionode = (fat_env%force_env_ref%para_env%mepos == fat_env%force_env_ref%para_env%source)

    e_tot = fat_e_tot(fat_env)

    if (ionode .and. fat_env%output_unit > 0) &
      write (fat_env%output_unit, fmt=fat_data_fmt) &
         fat_env%i_iter, &
         fat_env%n_iter, &
         i_frag, &
         fat_env%n_frag, &
         fat_env%e_sub(i_frag), &
         fat_env%e_emb(i_frag), &
         fat_env%ekin_ref, &
         fat_env%ekin_sub(i_frag), &
         fat_env%exc_ref, &
         fat_env%exc_sub(i_frag), &
         fat_env%eee, &
         fat_env%ene, &
         fat_env%enn, &
         e_tot

  end subroutine fat_status

! **************************************************************************************************
!> \brief Check the freeze-and-thaw optimization for convergence
!> \param fat_env The freeze-and-thaw environment containing all calculation data
!> \param i_frag The index of the current fragment being processed
! **************************************************************************************************
  subroutine fat_iterate_frag(fat_env, i_frag)
    type(fat_env_type), pointer :: fat_env
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_iterate_frag"

    character :: cmdmsg(100)
    character(default_string_length) :: label
    integer:: cmdstat, exitstat, input_unit, output_unit
    type(dft_control_type), pointer :: dft_control
    type(force_env_type), pointer :: force_env_tmp
    type(pw_type), pointer :: rho(:), vemb, weights !, rho_tmp
    type(qs_energy_type), pointer :: energy
    type(qs_environment_type), pointer :: qs_env_tmp
!   character(len=100) :: filename

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen, " for frag:", i_frag

    ! Get the force environment and system environment for current fragment
    force_env_tmp => fat_env%force_env_sub(i_frag)%force_env
    qs_env_tmp => force_env_tmp%qs_env

    if (is_debug) &
      print *, "KKK| ", modulen//"::"//routinen, "fat_env%i_iter:", fat_env%i_iter

    ! Get the DFT control settings and enable embedding potential
    call get_qs_env(qs_env_tmp, dft_control=dft_control)
    dft_control%apply_embed_pot = .true.

    ! Create a new embedding potential for current fragment
    call fat_pw_create_frag(fat_env, vemb, i_frag)

    ! After first iteration, calculate the embedding potential
    if (fat_env%i_iter > 0) then
       ! Calculate total embedding potential including:
       ! - Non-additive Kinetic potential
       ! - Non-additive XC potential
       ! - Electrostatic potential from other fragments
       ! (to_file = False)
       call fat_vemb_frag(fat_env, vemb, i_frag, .false.)
    end if

    ! Set the embedding potential to the qs env
    call set_qs_env(qs_env=qs_env_tmp, embed_pot=vemb)


    ! Handle external fragment calculations if configured
    if (fat_env%extern(i_frag)) then

       if (is_debug) print *, "KKK| Using extern for fragment ", i_frag

       call integer_to_string(i_frag, label)

       ! Write embedding potential to file for external program
       call open_file( &
            "vemb_"//trim(label)//".dat", file_status="REPLACE", &
            file_action="WRITE", unit_number=output_unit)

       ! Create weights for the potential based on grid spacing
       call fat_pw_create_frag(fat_env, weights, i_frag)
       call pw_set(weights, vemb%pw_grid%dr(1) * vemb%pw_grid%dr(2) * vemb%pw_grid%dr(3))

       ! Write potential with weights to output_unit file
       call pw_to_simple_volumetric(weights, output_unit, [1, 1, 1], vemb)

       ! Deallocate weights
       call fat_pw_release_frag(fat_env, weights, i_frag)

       ! Execute external calculation script
       call execute_command_line("sh extern_"//trim(label)//".sh", .true., exitstat, cmdstat, cmdmsg)
       if (cmdstat /= 0) CPABORT("External command did not execute successfully")

!      ! Get energy from DFT calculation
!      call get_qs_env(qs_env_tmp, energy=energy)
!      if (is_debug) print *, "KKK| IT IS ENERGY in fat_iterate_frag:", energy%total

       ! Get energy from external calculation
       call open_file( &
            "extern_"//trim(label)//".e", file_status="OLD", file_action="READ", unit_number=input_unit)
       read (input_unit, *) fat_env%e_sub(i_frag)

       ! Get electron density from external calculation
       call fat_rho_frag(fat_env, rho, i_frag)

       ! Calculate embedding energy for current fragment:
       ! Integrate the density * embedding potential
       fat_env%e_emb(i_frag) = pw_integral_ab(vemb, rho(1))

!      if (is_debug) then
!         ! Create temporary pw for density
!         call fat_pw_create_frag(fat_env, rho_tmp, i_frag)
!         call pw_transfer(rho(1), rho_tmp)
!         ! Write density to cube file
!         write(filename,'(A,I0,A,I0,A)') 'rho_frag', i_frag, '_iter', fat_env%i_iter, '.cube'
!         call write_pw_r_to_cube(fat_env, rho_tmp, filename)
!         ! Deallocate
!         call fat_pw_release_frag(fat_env, rho_tmp, i_frag)
!      end if

       if(is_debug) &
         print *, "KKK| Embedding energy from fragment ", i_frag, fat_env%e_emb(i_frag)

    else

      if(is_debug) &
        print *, "KKK| Using internal calculation for fragment ", i_frag

       ! Perform internal DFT calculation for current fragment:
       call qs_calc_energy_force(qs_env=qs_env_tmp, &
                                 calc_force=.true., &
                                 consistent_energies=.false., &
                                 linres=.false.)

       ! Get energy from DFT calculation
       call get_qs_env(qs_env_tmp, energy=energy)
       fat_env%e_sub(i_frag) = energy%total      ! Total energy

       ! Get embedding energy for current fragment
       fat_env%e_emb(i_frag) = energy%embed_corr

       if(is_debug) &
         print *, "KKK| Embedding energy from fragment ", i_frag, fat_env%e_emb(i_frag)

       ! Get electron density from internal calculation
       call fat_rho_frag(fat_env, rho, i_frag)

!      if (is_debug) then
!         ! Create temporary pw for density
!         call fat_pw_create_frag(fat_env, rho_tmp, i_frag)
!         call pw_transfer(rho(1), rho_tmp)
!         ! Write density to cube file
!         write(filename,'(A,I0,A,I0,A)') 'rho_frag', i_frag, '_iter', fat_env%i_iter, '.cube'
!         call write_pw_r_to_cube(fat_env, rho_tmp, filename)
!         ! Deallocate
!         call fat_pw_release_frag(fat_env, rho_tmp, i_frag)
!      end if

    end if

    if (is_debug) then
       print *, "KKK| Total subsystem (", i_frag, ") energy:", fat_env%e_sub(i_frag)
    end if

    ! Deallocate vemb potential
    call fat_pw_release_frag(fat_env, vemb, i_frag)

  end subroutine fat_iterate_frag

! **************************************************************************************************
!> \brief Performs one iteration of the freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_iterate(fat_env)
    type(fat_env_type), pointer :: fat_env

    character(len=*), parameter :: routinen = "fat_iterate"

    integer :: i_frag

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    ! Loop over all fragments to optimize their densities with embedding effects
    do i_frag = 1, fat_env%n_frag

       ! Calculate new density and energy for fragment i_frag under embedding potential from others
       call fat_iterate_frag(fat_env, i_frag)

       ! Print current optimization status
       call fat_status(fat_env, i_frag)
    end do

    ! Increment iteration counter
    fat_env%i_iter = fat_env%i_iter + 1

  end subroutine fat_iterate

! **************************************************************************************************
!> \brief Check the freeze-and-thaw optimization for convergence
! **************************************************************************************************
  function fat_is_converged(fat_env) result(is_converged)
    type(fat_env_type), pointer :: fat_env
    logical :: is_converged

    character(len=*), parameter :: routinen = "fat_is_converged"

    real(dp) :: diff

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    diff = fat_env%e_tot - fat_env%e_old
    is_converged = abs(diff) .le. fat_env%eps
    fat_env%is_converged = is_converged

    if (is_debug) then
       print *, "DBG| ", is_converged
       print *, "DBG| ", fat_env%e_tot, fat_env%e_old
       print *, "DBG| diff", diff
    end if

    fat_env%e_old = fat_env%e_tot

  end function

! **************************************************************************************************
!> \brief Check if the freeze-and-thaw optimization has converged or run into the iteration limit
! **************************************************************************************************
  function fat_is_done(fat_env) result(is_done)
    type(fat_env_type), pointer :: fat_env
    logical :: is_done

    character(len=*), parameter :: routinen = "fat_is_done"

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    is_done = fat_env%i_iter > fat_env%n_iter .or. fat_is_converged(fat_env)

  end function

! **************************************************************************************************
!> \brief Allocate the forces for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_create(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable, intent(out) :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_create"

    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, particle_set=particles)
    allocate(forces(size(particles), 3))
    forces = 0.0

  end subroutine fat_force_create

! **************************************************************************************************
!> \brief Calculate the self force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_self(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_self"

    integer :: i, i_atom, i_kind, j, k, n
    integer, allocatable :: atom_of_kind(:), kind_of(:)
    type(atomic_kind_type), pointer :: kinds(:)
    type(qs_environment_type), pointer :: qs_env
    type(qs_force_type), pointer :: qs_force(:)

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    k = 0
    do i = 1, fat_env%n_frag
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, atomic_kind_set=kinds, force=qs_force)
       call get_atomic_kind_set(kinds, natom=n, kind_of=kind_of, atom_of_kind=atom_of_kind)
       do j = 1, n
          k = k + 1
          i_kind = kind_of(j)
          i_atom = atom_of_kind(j)
          forces(k, :) = qs_force(i_kind)%total(1:3, i_atom)
          if (is_debug) print *, "DBG| k, i_kind, i_atom", k, i_kind, i_atom
       end do
       deallocate(atom_of_kind, kind_of)
    end do

  end subroutine fat_force_self

! **************************************************************************************************
!> \brief Calculate the embedding force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_emb(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_emb"

    integer :: a, i, j, k
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_type), pointer :: v, drho

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, v, i)
       call fat_vemb_frag(fat_env, v, i, .true.)
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_pw_create_frag(fat_env, drho, i)
             call fat_drho_frag(fat_env, drho, i, j, a)
             forces(k, a) = pw_integral_ab(v, drho)
             call fat_pw_release_frag(fat_env, drho, i)
          end do
       end do
       call fat_pw_release_frag(fat_env, v, i)
    end do

  end subroutine fat_force_emb

! **************************************************************************************************
!> \brief Calculate the non-additive kinetic force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_tsna(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_tsna"

    integer :: a, i, j, k
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_type), pointer :: vemb, drho

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    k = 0
    do i = 1, fat_env%n_frag  ! fragments

      ! Create embedding potential
       call fat_pw_create_frag(fat_env, vemb, i)

       ! Calculate total embedding potential (includes both LDA and GGA if present)
       CALL fat_vemb_tsna_frag(fat_env=fat_env, vemb=vemb, i_frag=i)

       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)

       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_pw_create_frag(fat_env, drho, i)
             call fat_drho_frag(fat_env, drho, i, j, a)

             ! Total force contribution (LDA + GGA if present)
             forces(k, a) = pw_integral_ab(vemb, drho)

             call fat_pw_release_frag(fat_env, drho, i)
          end do
       end do

       ! Release potential
       call fat_pw_release_frag(fat_env, vemb, i)

    end do

  end subroutine fat_force_tsna

! **************************************************************************************************
!> \brief Calculate the non-additive xc force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_xcna(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_xcna"

    integer :: a, i, j, k
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_type), pointer :: vemb, drho

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    k = 0
    do i = 1, fat_env%n_frag  ! fragments

      ! Create embedding potential
       call fat_pw_create_frag(fat_env, vemb, i)

       ! Calculate total embedding potential (includes both LDA and GGA if present)
       call fat_vemb_xcna_frag(fat_env=fat_env, vemb=vemb, i_frag=i)

       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)

       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_pw_create_frag(fat_env, drho, i)
             call fat_drho_frag(fat_env, drho, i, j, a)

             ! Total force contribution (LDA + GGA if present)
             forces(k, a) = pw_integral_ab(vemb, drho)

             call fat_pw_release_frag(fat_env, drho, i)
          end do
       end do

       ! Release potential
       call fat_pw_release_frag(fat_env, vemb, i)

    end do

  end subroutine fat_force_xcna

! **************************************************************************************************
!> \brief Calculate the ee force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_ee(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_ee"

    integer :: a, i, j, k
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_type), pointer :: v, drho

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, v, i)
       call fat_vemb_e_frag(fat_env, v, i)
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_pw_create_frag(fat_env, drho, i)
             call fat_drho_frag(fat_env, drho, i, j, a)
             forces(k, a) = pw_integral_ab(v, drho)
             call fat_pw_release_frag(fat_env, drho, i)
          end do
       end do
       call fat_pw_release_frag(fat_env, v, i)
    end do

  end subroutine fat_force_ee

! **************************************************************************************************
!> \brief Calculate the ne force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_ne(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_ee"

    integer :: a, i, j, k
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_type), pointer :: ve, vn, drho, drhon

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, ve, i)
       call fat_vemb_e_frag(fat_env, ve, i)

       call fat_pw_create_frag(fat_env, vn, i)
       call fat_vemb_n_frag(fat_env, vn, i)

       qs_env => fat_env%force_env_sub(i)%force_env%qs_env

       call get_qs_env(qs_env, particle_set=particles)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_pw_create_frag(fat_env, drho, i)
             call fat_drho_frag(fat_env, drho, i, j, a)
             forces(k, a) = pw_integral_ab(vn, drho)
             call fat_pw_release_frag(fat_env, drho, i)

             call fat_pw_create_frag(fat_env, drhon, i)
             call fat_drhon_frag(fat_env, drhon, i, j, a)
             forces(k, a) = pw_integral_ab(ve, drhon)
             call fat_pw_release_frag(fat_env, drhon, i)
          end do
       end do

       call fat_pw_release_frag(fat_env, ve, i)
       call fat_pw_release_frag(fat_env, vn, i)
    end do

  end subroutine fat_force_ne

! **************************************************************************************************
!> \brief Calculate the nn force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_nn(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_nn"

    integer :: a, i, j, k
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_type), pointer :: v, drho

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, v, i)
       call fat_vemb_n_frag(fat_env, v, i)
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_pw_create_frag(fat_env, drho, i)
             call fat_drhon_frag(fat_env, drho, i, j, a)
             forces(k, a) = pw_integral_ab(v, drho)
             call fat_pw_release_frag(fat_env, drho, i)
          end do
       end do
       call fat_pw_release_frag(fat_env, v, i)
    end do

  end subroutine fat_force_nn

! ************************************************************************************************
!> \brief Calculate the overlap force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_overlap(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_overlap"

    integer :: i, iatom, ikind, natom
    integer, allocatable :: atom_of_kind(:), kind_of(:)
    type(atomic_kind_type), pointer :: atomic_kinds(:)
    type(qs_environment_type), pointer :: qs_env
    type(qs_force_type), pointer :: force(:)

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, atomic_kind_set=atomic_kinds, force=force)
    call get_atomic_kind_set( &
         atomic_kinds, atom_of_kind=atom_of_kind, kind_of=kind_of, natom=natom)

    do iatom = 1, natom
       ikind = kind_of(iatom)
       i = atom_of_kind(iatom)
       forces(iatom, :) = force(ikind)%overlap(1:3, i)
    end do

  end subroutine fat_force_overlap

! **************************************************************************************************
!> \brief Print the forces for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_print(fat_env, forces, label)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable, intent(in) :: forces(:, :)
    character(len=*) :: label

    character(len=*), parameter :: routinen = "fat_force_print"

    integer :: i, n, output_unit
    character(len=2) :: symbol
    type(particle_type), pointer :: particles(:)
    type(atomic_kind_type), pointer :: kind
    type(qs_environment_type), pointer :: qs_env

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, particle_set=particles)
    n = size(particles)

    call open_file( &
         label//".xyz", file_status="REPLACE", file_action="WRITE", unit_number=output_unit)
    write (output_unit, '(I0)') size(particles)
    write (output_unit, *)
    do i = 1, n
       kind => particles(i)%atomic_kind
       call get_atomic_kind(kind, element_symbol=symbol)
       write (output_unit, '(A2,3(F26.15," "))') symbol, forces(i, 1), forces(i, 2), forces(i, 3)
    end do
    call close_file(output_unit)

  end subroutine fat_force_print

! **************************************************************************************************
!> \brief Calculate the forces for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force(fat_env)
    type(fat_env_type), pointer :: fat_env

    character(len=*), parameter :: routinen = "fat_force"

    integer :: i, j
    real(dp), allocatable :: &
         & forces(:, :), forces_self(:, :), forces_emb(:, :),&
         & forces_tsna(:, :), forces_xcna(:, :), forces_ee(:, :),&
         & forces_ne(:, :), forces_nn(:, :), forces_overlap(:, :)

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call fat_force_create(fat_env, forces)

    call fat_force_create(fat_env, forces_self)
    call fat_force_self(fat_env, forces_self)
    call fat_force_print(fat_env, forces_self, "force_self")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_self(i, j)
       end do
    end do
    deallocate(forces_self)

    call fat_force_create(fat_env, forces_emb)
    call fat_force_emb(fat_env, forces_emb)
    call fat_force_print(fat_env, forces_emb, "force_emb")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) - forces_emb(i, j)
       end do
    end do
    deallocate(forces_emb)

    call fat_force_create(fat_env, forces_tsna)
    call fat_force_tsna(fat_env, forces_tsna)
    call fat_force_print(fat_env, forces_tsna, "force_tsna")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_tsna(i, j)
       end do
    end do
    deallocate(forces_tsna)

    call fat_force_create(fat_env, forces_xcna)
    call fat_force_xcna(fat_env, forces_xcna)
    call fat_force_print(fat_env, forces_xcna, "force_xcna")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_xcna(i, j)
       end do
    end do
    deallocate(forces_xcna)

    call fat_force_create(fat_env, forces_ee)
    call fat_force_ee(fat_env, forces_ee)
    call fat_force_print(fat_env, forces_ee, "force_ee")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_ee(i, j)
       end do
    end do
    deallocate(forces_ee)

    call fat_force_create(fat_env, forces_ne)
    call fat_force_ne(fat_env, forces_ne)
    call fat_force_print(fat_env, forces_ne, "force_ne")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + 2*forces_ne(i, j)
       end do
    end do
    deallocate(forces_ne)

    call fat_force_create(fat_env, forces_nn)
    call fat_force_nn(fat_env, forces_nn)
    call fat_force_print(fat_env, forces_nn, "force_nn")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_nn(i, j)
       end do
    end do
    deallocate(forces_nn)

    call fat_force_create(fat_env, forces_overlap)
    call fat_force_overlap(fat_env, forces_overlap)
    call fat_force_print(fat_env, forces_overlap, "force_overlap")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_overlap(i, j)
       end do
    end do
    deallocate(forces_overlap)

    call fat_force_print(fat_env, forces, "force")

    deallocate(forces)

  end subroutine fat_force

! **************************************************************************************************
!> \brief Print the end results of a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_results(fat_env, do_forces)
    type(fat_env_type), pointer :: fat_env
    logical, optional :: do_forces

    character(len=*), parameter :: routinen = "fat_results"

    logical :: my_do_forces
    integer :: i, j
    real(dp) :: e_tmp

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    if (present(do_forces)) then
       my_do_forces = do_forces
    else
       my_do_forces = .false.
    end if

    do i = 1, fat_env%n_frag
       write (fat_env%output_unit, '(T2,"FAT| self energy fragment   ",I2,":",F26.15)') i, fat_env%e_sub(i)
    end do

    do i = 1, fat_env%n_frag
       write (fat_env%output_unit, '(T2,"FAT| emb energy fragment    ",I2,":",F26.15)') i, fat_env%e_emb(i)
    end do

    write (fat_env%output_unit, '(T2,"FAT| kin energy reference     :",F26.15)') fat_env%ekin_ref
    do i = 1, fat_env%n_frag
       write (fat_env%output_unit, '(T2,"FAT| kin energy fragment    ",I2,":",F26.15)') i, fat_env%ekin_sub(i)
    end do

    write (fat_env%output_unit, '(T2,"FAT| xc energy reference      :",F26.15)') fat_env%exc_ref
    do i = 1, fat_env%n_frag
       write (fat_env%output_unit, '(T2,"FAT| xc energy fragment     ",I2,":",F26.15)') i, fat_env%exc_sub(i)
    end do

    do i = 1, fat_env%n_frag
       do j = 1, fat_env%n_frag
          if (i == j) cycle
          e_tmp = fat_eint_ee(fat_env, i, j)
          write (fat_env%output_unit, '(T2,"FAT| ee energy fragments ",I2," ",I2,":",F26.15)') i, j, e_tmp
       end do
    end do

    do i = 1, fat_env%n_frag
       do j = 1, fat_env%n_frag
          if (i == j) cycle
          e_tmp = fat_eint_ne(fat_env, i, j)
          write (fat_env%output_unit, '(T2,"FAT| ne energy fragments ",I2," ",I2,":",F26.15)') i, j, e_tmp
       end do
    end do

    fat_env%enn = 0.0
    do i = 1, fat_env%n_frag
       do j = 1, fat_env%n_frag
          if (i == j) cycle
          e_tmp = fat_eint_nn(fat_env, i, j)
          write (fat_env%output_unit, '(T2,"FAT| nn energy fragments ",I2," ",I2,":",F26.15)') i, j, e_tmp
       end do
    end do

    write (fat_env%output_unit, '(/,T2,"ENERGY| Total FORCE_EVAL ( ",A," ) energy [a.u.]: ",T55,F26.15)') &
       "FAT", fat_env%e_tot

    if (my_do_forces) then
       call fat_force(fat_env)
    end if

  end subroutine fat_results

! **************************************************************************************************
!> \brief Cleans up a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_finalize(fat_env)
    type(fat_env_type), pointer :: fat_env

    character(len=*), parameter :: routinen = "fat_finalize"

    integer :: i
    type(dft_control_type), pointer :: dft_control
    type(qs_environment_type), pointer :: qs_env

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    do i = 1, fat_env%n_frag
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, dft_control=dft_control)
       dft_control%apply_embed_pot = .false.
    end do

  end subroutine fat_finalize

! **************************************************************************************************
!> \brief Initialize the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_init(fat_env, force_env, ref_subsys_number)
    type(fat_env_type), pointer   :: fat_env
    type(force_env_type), pointer :: force_env
    integer                       :: ref_subsys_number

    character(len=*), parameter :: routinen = "fat_init"

    integer :: handle, i, j, n_iter
    real(dp) :: eps
    type(force_env_type), pointer :: force_env_tmp
    type(pw_env_type), pointer :: pw_env_tmp
    type(qs_environment_type), pointer :: qs_env_tmp
    type(section_vals_type), pointer :: input
    TYPE(cp_logger_type), POINTER :: logger
    logical :: ionode

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    nullify(qs_env_tmp)

    ! This is a limitation of the EMBED method
    CPASSERT(ref_subsys_number == 1)

    fat_env%force_env_ref => force_env%sub_force_env(ref_subsys_number)%force_env

    ! Count fragments
    fat_env%n_frag = size(force_env%sub_force_env) - 1

    ! Get logger and initialize output unit only for I/O node
    ionode = force_env%para_env%mepos == force_env%para_env%source
    IF (ionode) then
       logger => cp_get_default_logger()
       fat_env%output_unit = cp_logger_get_default_io_unit(logger)
    else
       fat_env%output_unit = -1
    end if

    ! Allocate and setup force environments for each fragment
    allocate(fat_env%force_env_sub(fat_env%n_frag))
    j = 1
    do i = 1, size(force_env%sub_force_env)
       if (i == ref_subsys_number) cycle ! Skip reference system

       ! Store pointer to fragment's force environment
       fat_env%force_env_sub(j)%force_env => force_env%sub_force_env(i)%force_env
       j = j + 1
    end do

    ! Setup plane wave pools for density/potential calculations
    ! For reference system
    force_env_tmp => fat_env%force_env_ref
    qs_env_tmp => force_env_tmp%qs_env
    call get_qs_env(qs_env_tmp, pw_env=pw_env_tmp)
    call pw_env_get(pw_env_tmp, auxbas_pw_pool=fat_env%pw_pool_ref)

    ! For each fragment
    allocate(fat_env%pw_pool_sub(fat_env%n_frag))
    do i = 1, fat_env%n_frag
       force_env_tmp => fat_env%force_env_sub(i)%force_env
       qs_env_tmp => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env_tmp, pw_env=pw_env_tmp)
       call pw_env_get(pw_env_tmp, auxbas_pw_pool=fat_env%pw_pool_sub(i)%pool)
    end do

    ! Setup flags for external fragments, if specified
    allocate(fat_env%extern(fat_env%n_frag))
    do i = 1, fat_env%n_frag
       force_env_tmp => fat_env%force_env_sub(i)%force_env
       qs_env_tmp => force_env_tmp%qs_env
       call get_qs_env(qs_env_tmp, input=input)
       ! Check if this fragment uses external calculation
       call section_vals_val_get(input, "DFT%QS%FAT_EXTERN", l_val=fat_env%extern(i))

       if (is_debug) print *, "DBG| fat_env%extern(i)", i, fat_env%extern(i)
    end do

    ! Initial DFT calculations for all fragments (without embedding potential)
    do i = 1, fat_env%n_frag
       force_env_tmp => fat_env%force_env_sub(i)%force_env
       qs_env_tmp => force_env_tmp%qs_env
       CALL qs_calc_energy_force(qs_env=qs_env_tmp, &
                                 calc_force=.true., &
                                 consistent_energies=.false., &
                                 linres=.false.)

       if (is_debug) print *, "KKK| initial DFT calculation for fragment:", i, fat_env%n_frag
    end do
    if (is_debug) print *, "KKK| After initial DFT calculations for all fragments"

    ! Allocate arrays for storing various energy components
    allocate(fat_env%e_sub(fat_env%n_frag)) ! Subsystem energies
    fat_env%e_sub = 0.0

    allocate(fat_env%e_emb(fat_env%n_frag)) ! Embedding energies
    fat_env%e_emb = 0.0

    fat_env%ekin_ref = 0.0                     ! Reference kinetic energy
    allocate(fat_env%ekin_sub(fat_env%n_frag)) ! Subsystem kinetic energies
    fat_env%ekin_sub = 0.0

    fat_env%exc_ref = 0.0                     ! Reference XC energy
    allocate(fat_env%exc_sub(fat_env%n_frag)) ! Subsystem XC energies
    fat_env%exc_sub = 0.0

    ! Initialize interaction electrostatic energies
    fat_env%eee = 0.0 ! Electron-electron
    fat_env%ene = 0.0 ! Electron-nuclear
    fat_env%enn = 0.0 ! Nuclear-nuclear

    ! Initialize "old" energy for convergence check
    fat_env%e_old = 0.0

    ! Initialize iteration counter
    fat_env%i_iter = 0
    ! Calculate initial energies for each fragment
    do i = 1, fat_env%n_frag
       call fat_iterate_frag(fat_env, i)
    end do
    if (is_debug) print *, "KKK| After initial calculations for all fragments i_iter", i, fat_env%i_iter

    ! Get convergence parameters from input
    call get_qs_env(fat_env%force_env_ref%qs_env, input=input)
    call section_vals_val_get(input, "DFT%QS%OPT_EMBED%N_ITER", i_val=n_iter)
    fat_env%n_iter = n_iter
    call section_vals_val_get(input, "DFT%QS%OPT_EMBED%DENS_CONV_MAX", r_val=eps)
    fat_env%eps = eps

    ! Setup header for status output
    if (ionode) &
      write (fat_env%output_unit, fmt=fat_title_fmt) &
         "i_iter", &
         "n_iter", &
         "i_frag", &
         "n_frag", &
         "self", &
         "emb", &
         "kin ref", &
         "kin frag", &
         "xc ref", &
         "xc frag", &
         "ee", &
         "ne", &
         "nn", &
         "tot"

    ! Start iterative freeze-and-thaw density
    ! optimization self-consistent procedure
    call fat_iterate(fat_env)

    call timestop(handle)

  end subroutine fat_init

end module fat_methods
