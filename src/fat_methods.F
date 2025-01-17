module fat_methods
  use atomic_kind_types, only: atomic_kind_type, get_atomic_kind, get_atomic_kind_set
  use cp_control_types, only: dft_control_type
  use cp_files, only: open_file, close_file
  use cp_log_handling, only: cp_logger_type, cp_logger_get_default_io_unit, cp_get_default_logger
  use cp_realspace_grid_cube, only: cp_pw_to_cube, cp_cube_to_pw
  USE dbcsr_api, ONLY: dbcsr_type, dbcsr_p_type
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
  USE pw_types, ONLY: pw_r3d_rs_type, pw_c1d_gs_type
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
  subroutine fat_pw_create_frag(fat_env, pw_r, pw_g, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), optional, pointer :: pw_r
    type(pw_c1d_gs_type), optional, pointer :: pw_g
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_pw_create_frag"

    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    CPASSERT(PRESENT(pw_r) .OR. PRESENT(pw_g))
    CPASSERT(ASSOCIATED(fat_env%pw_pool_sub))
    CPASSERT(i_frag >= 1 .AND. i_frag <= SIZE(fat_env%pw_pool_sub))
    CPASSERT(ASSOCIATED(fat_env%pw_pool_sub(i_frag)%pool))

    IF (PRESENT(pw_r)) THEN
       ALLOCATE(pw_r)
       CALl fat_env%pw_pool_sub(i_frag)%pool%create_pw(pw_r)
       CALL pw_zero(pw_r)
    END IF

    IF (PRESENT(pw_g)) THEN
       ALLOCATE(pw_g)
       CALl fat_env%pw_pool_sub(i_frag)%pool%create_pw(pw_g)
       CALL pw_zero(pw_g)
    END IF

    call timestop(handle)

  end subroutine fat_pw_create_frag

! **************************************************************************************************
!> \brief Releases a pw for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pw_release_frag(fat_env, pw_r, pw_g, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), optional, pointer :: pw_r
    type(pw_c1d_gs_type), optional, pointer :: pw_g
    integer :: i_frag, handle

    character(len=*), parameter :: routinen = "fat_pw_release_frag"

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    CPASSERT(PRESENT(pw_r) .OR. PRESENT(pw_g))

    IF (PRESENT(pw_r)) THEN
       CPASSERT(ASSOCIATED(pw_r))
       call fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(pw_r)
       deallocate(pw_r)
    END IF

    IF (PRESENT(pw_g)) THEN
       CPASSERT(ASSOCIATED(pw_g))
       call fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(pw_g)
       deallocate(pw_g)
    END IF

    if (is_debug) print *, "DBG| released"

    call timestop(handle)

  end subroutine fat_pw_release_frag

! **************************************************************************************************
!> \brief Creates a pw for the reference system in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pw_create_ref(fat_env, pw_r, pw_g)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), optional, intent(out) :: pw_r
    type(pw_c1d_gs_type), optional, intent(out) :: pw_g

    character(len=*), parameter :: routinen = "fat_pw_create_ref"

    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    CPASSERT(PRESENT(pw_r) .OR. PRESENT(pw_g))
    CPASSERT(ASSOCIATED(fat_env%pw_pool_ref))

    IF (PRESENT(pw_r)) THEN
       call fat_env%pw_pool_ref%create_pw(pw_r)
       call pw_zero(pw_r)
    END IF

    IF (PRESENT(pw_g)) THEN
       call fat_env%pw_pool_ref%create_pw(pw_g)
       call pw_zero(pw_g)
    END IF

    call timestop(handle)

  end subroutine fat_pw_create_ref

! **************************************************************************************************
!> \brief Releases a pw for the reference system in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pw_release_ref(fat_env, pw_r, pw_g)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), optional, intent(inout) :: pw_r
    type(pw_c1d_gs_type), optional, intent(inout) :: pw_g

    character(len=*), parameter :: routinen = "fat_pw_release_ref"

    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    CPASSERT(PRESENT(pw_r) .OR. PRESENT(pw_g))

    IF (PRESENT(pw_r)) THEN
       call fat_env%pw_pool_ref%give_back_pw(pw_r)
    END IF

    IF (PRESENT(pw_g)) THEN
       call fat_env%pw_pool_ref%give_back_pw(pw_g)
    END IF

    if (is_debug) print *, "DBG| released"

    call timestop(handle)

  end subroutine fat_pw_release_ref

! **************************************************************************************************
!> \brief Creates an array of pws for the reference system in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pws_create_ref(fat_env, pws_r, pws_g, n_pws)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), optional, pointer :: pws_r(:)
    type(pw_c1d_gs_type), optional, pointer :: pws_g(:)
    integer :: n_pws

    character(len=*), parameter :: routinen = "fat_pws_create_ref"

    integer :: i, handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    CPASSERT(PRESENT(pws_r) .OR. PRESENT(pws_g))

    IF (PRESENT(pws_r)) THEN
       allocate(pws_r(n_pws))
       do i = 1, n_pws
          call fat_pw_create_ref(fat_env, pw_r=pws_r(i))
       end do
    END IF

    IF (PRESENT(pws_g)) THEN
       allocate(pws_g(n_pws))
       do i = 1, n_pws
          call fat_pw_create_ref(fat_env, pw_g=pws_g(i))
       end do
    END IF

    call timestop(handle)

  end subroutine fat_pws_create_ref

! **************************************************************************************************
!> \brief Releases an array of pws for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pws_release_ref(fat_env, pws_r, pws_g)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), optional, pointer :: pws_r(:)
    type(pw_c1d_gs_type), optional, pointer :: pws_g(:)

    character(len=*), parameter :: routinen = "fat_pws_release_ref"

    integer :: i, handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    CPASSERT(PRESENT(pws_r) .OR. PRESENT(pws_g))

    IF (PRESENT(pws_r)) THEN
       CPASSERT(ASSOCIATED(pws_r))
       do i = 1, size(pws_r)
          call fat_pw_release_ref(fat_env, pw_r=pws_r(i))
       END DO
       DEALLOCATE(pws_r)
    END IF

    IF (PRESENT(pws_g)) THEN
       CPASSERT(ASSOCIATED(pws_g))
       do i = 1, size(pws_g)
          call fat_pw_release_ref(fat_env, pw_g=pws_g(i))
       END DO
       DEALLOCATE(pws_g)
    END IF

    if (is_debug) print *, "DBG| released"

    call timestop(handle)

  end subroutine fat_pws_release_ref

! **************************************************************************************************
!> \brief Creates an array of pws for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  SUBROUTINE fat_pws_create_frag(fat_env, pws_r, pws_g, i_frag, n_pws)
    TYPE(fat_env_type), POINTER                       :: fat_env
    TYPE(pw_r3d_rs_type), OPTIONAL, POINTER           :: pws_r(:)
    TYPE(pw_c1d_gs_type), OPTIONAL, POINTER           :: pws_g(:)
    INTEGER, INTENT(IN)                               :: i_frag
    INTEGER, INTENT(IN)                               :: n_pws

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_pws_create_frag'

    INTEGER                                           :: handle, i

    CALL timeset(routineN, handle)

    IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

    CPASSERT(PRESENT(pws_r) .OR. PRESENT(pws_g))

    IF (PRESENT(pws_r)) THEN
       allocate(pws_r(n_pws))
       DO i = 1, n_pws
          call fat_env%pw_pool_sub(i_frag)%pool%create_pw(pws_r(i))
          call pw_zero(pws_r(i))
       END DO
    END IF

    IF (PRESENT(pws_g)) THEN
       allocate(pws_g(n_pws))
       DO i = 1, n_pws
          call fat_env%pw_pool_sub(i_frag)%pool%create_pw(pws_g(i))
          call pw_zero(pws_g(i))
       END DO
    END IF

    CALL timestop(handle)

   END SUBROUTINE fat_pws_create_frag

! **************************************************************************************************
!> \brief Releases an array of pws for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pws_release_frag(fat_env, pws_r, pws_g, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), optional, pointer :: pws_r(:)
    type(pw_c1d_gs_type), optional, pointer :: pws_g(:)
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_pws_release_frag"

    type(pw_r3d_rs_type), pointer :: pw_r_tmp
    type(pw_c1d_gs_type), pointer :: pw_g_tmp
    integer :: i, handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    CPASSERT(PRESENT(pws_r) .OR. PRESENT(pws_g))

    IF (PRESENT(pws_r)) THEN
       CPASSERT(ASSOCIATED(pws_r))
       do i = 1, size(pws_r)
          pw_r_tmp => pws_r(i)
          call fat_pw_release_frag(fat_env, pw_r=pw_r_tmp, i_frag=i_frag)
       end do
    END IF

    IF (PRESENT(pws_g)) THEN
       CPASSERT(ASSOCIATED(pws_g))
       do i = 1, size(pws_g)
          pw_g_tmp => pws_g(i)
          call fat_pw_release_frag(fat_env, pw_g=pw_g_tmp, i_frag=i_frag)
       end do
    END IF

    if (is_debug) print *, "DBG| released"

    call timestop(handle)

  end subroutine fat_pws_release_frag

! **************************************************************************************************
!> \brief Gets the electron density for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_get_rho_frag(fat_env, rho_r, rho_g, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), optional, pointer :: rho_r(:)
    type(pw_c1d_gs_type), optional, pointer :: rho_g(:)
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_get_rho_frag"

    character(len=default_string_length) :: label
    type(qs_environment_type), pointer :: qs_env
    type(qs_rho_type), pointer :: qs_rho
    type(pw_r3d_rs_type), pointer :: rho_r_tmp
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    CPASSERT(PRESENT(rho_r) .OR. PRESENT(rho_g))

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, rho=qs_rho)

    IF (PRESENT(rho_r)) THEN
       NULLIFY(rho_r)
       CALL qs_rho_get(qs_rho, rho_r=rho_r)
    END IF

    IF (PRESENT(rho_g)) THEN
       NULLIFY(rho_g)
       CALL qs_rho_get(qs_rho, rho_g=rho_g)
    END IF

    IF (fat_env%extern(i_frag)) THEN

       CALL integer_to_string(i_frag, label)

       IF (PRESENT(rho_r)) THEN
          CALL cp_cube_to_pw(rho_r(1), "extern_"//TRIM(label)//".cube", 1.0_dp, .FALSE.)
       END IF

       IF (PRESENT(rho_g)) THEN
          CALL fat_pw_create_frag(fat_env, pw_r=rho_r_tmp, i_frag=i_frag)
          CALL cp_cube_to_pw(rho_r_tmp, "extern_"//TRIM(label)//".cube", 1.0_dp, .FALSE.)
          CALL pw_transfer(rho_r_tmp, rho_g(1))
          CALL fat_pw_release_frag(fat_env, pw_r=rho_r_tmp, i_frag=i_frag)
       END IF

    END IF

!   if (is_debug) call print_pw_statistics(rho(1), "rho", i_frag)
!
    call timestop(handle)

  end subroutine fat_get_rho_frag

  subroutine write_pw_to_cube(fat_env, pw_r, pw_g, filename)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), optional, pointer :: pw_r
    type(pw_c1d_gs_type), optional, pointer :: pw_g
    character(len=*), intent(in) :: filename

    type(cp_logger_type), pointer :: logger
    TYPE(particle_list_type), POINTER :: particles
    type(particle_type), dimension(:), pointer :: particle_set
    type(qs_environment_type), pointer :: qs_env
    type(dft_control_type), pointer :: dft_control
    type(pw_env_type), pointer :: pw_env
    type(qs_subsys_type), pointer :: subsys
    type(pw_pool_type), pointer :: pw_pool
    type(pw_r3d_rs_type), pointer :: pw_r_tmp
    integer, pointer :: stride(:)
    integer :: unit_nr

    CPASSERT(PRESENT(pw_r) .OR. PRESENT(pw_g))

    logger => cp_get_default_logger()
    unit_nr = cp_logger_get_default_io_unit(logger)
    log_counter = log_counter + 1

    if (log_counter>6) return

    print *, "KKK| write CUBE file number:", filename

    IF (PRESENT(pw_r)) THEN
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
          call cp_pw_to_cube(pw_r, unit_nr, trim(filename), particles=particles, stride=stride)
          call close_file(unit_number=unit_nr)
          deallocate(stride)
       end if
    END IF

    IF (PRESENT(pw_g)) THEN
       if (unit_nr > 0) then
          call fat_pw_create_ref(fat_env, pw_r=pw_r_tmp)
          call pw_zero(pw_r_tmp)
          call pw_transfer(pw_g, pw_r_tmp)

          qs_env => fat_env%force_env_ref%qs_env

          CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, pw_env=pw_env, &
                          particle_set=particle_set, subsys=subsys)
          CALL pw_env_get(pw_env, auxbas_pw_pool=pw_pool)
          CALL qs_subsys_get(subsys, particles=particles)

          allocate(stride(3))
          stride = [1, 1, 1]

          call open_file(file_name=trim(filename), file_status="REPLACE", &
                         file_action="WRITE", unit_number=unit_nr)
          call cp_pw_to_cube(pw_r_tmp, unit_nr, trim(filename), particles=particles, stride=stride)
          call close_file(unit_number=unit_nr)
          deallocate(stride)

          call fat_pw_release_ref(fat_env, pw_r=pw_r_tmp)
       end if
    END IF

end subroutine write_pw_to_cube

! **************************************************************************************************
!> \brief Gets the nuclear density for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_get_rhon_frag(fat_env, rhon_g, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_c1d_gs_type), pointer :: rhon_g
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_get_rhon_frag"

    type(qs_environment_type), pointer :: qs_env
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    NULLIFY(rhon_g)

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env=qs_env, rho_core=rhon_g)

    call timestop(handle)

  end subroutine fat_get_rhon_frag

! **************************************************************************************************
!> \brief Gets the gradient of the nuclear density w.r.t. a nuclear position for one fragment in the
!>        freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_get_drhon_frag(fat_env, drhon_g, i_frag, i_nuc, i_dir)
    type(fat_env_type), pointer :: fat_env
    type(pw_c1d_gs_type), pointer :: drhon_g
    integer :: i_frag, i_nuc, i_dir

    character(len=*), parameter :: routinen = "fat_get_drhon_frag"

    type(qs_environment_type), pointer :: qs_env
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call calculate_drho_core(drhon_g, qs_env, i_dir, i_nuc)

    call timestop(handle)

  end subroutine fat_get_drhon_frag

! **************************************************************************************************
!> \brief Gets the gradient of the electronic density w.r.t. a nuclear position for one fragment in
!>        the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_get_drho_frag(fat_env, drho_r,  drho_g, i_frag, i_nuc, i_dir)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), pointer :: drho_r
    type(pw_c1d_gs_type), pointer :: drho_g
    integer :: i_frag, i_nuc, i_dir

    character(len=*), parameter :: routinen = "fat_get_drho_frag"

    type(qs_environment_type), pointer :: qs_env
    type(dbcsr_p_type), pointer :: rho_ao(:)
    type(qs_rho_type), pointer :: rho
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    CPASSERT(ASSOCIATED(drho_r) .AND. ASSOCIATED(drho_g))

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, rho=rho)
    call qs_rho_get(rho, rho_ao=rho_ao)

    call calculate_drho_elec_dR(matrix_p=rho_ao(1)%matrix, drho=drho_r, drho_gspace=drho_g, qs_env=qs_env, beta=i_dir, lambda=i_nuc)

    call timestop(handle)

  end subroutine fat_get_drho_frag

! **************************************************************************************************
!> \brief Gets the total electron density for the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_get_rho_ref(fat_env, rho_r, rho_g)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), optional :: rho_r(:)
    type(pw_c1d_gs_type), pointer, optional :: rho_g(:)

    character(len=*), parameter :: routinen = "fat_get_rho_ref"

    type(pw_r3d_rs_type), pointer :: rho_r_tmp(:)
    integer :: i, j, handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    do i = 1, fat_env%n_frag
       call fat_get_rho_frag(fat_env, rho_r=rho_r_tmp, i_frag=i)
       do j = 1, size(rho_r_tmp)
          call pw_axpy(rho_r_tmp(j), rho_r(j), allow_noncompatible_grids=.true.)
       end do
    end do

    IF (PRESENT(rho_g)) THEN
       CPASSERT(ASSOCIATED(rho_g))
       do j = 1, size(rho_r)
          call pw_transfer(rho_r(j), rho_g(j))
       end do
    END IF

    call timestop(handle)

  end subroutine fat_get_rho_ref

! **************************************************************************************************
!> \brief Calculate kinetic potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_ts_frag(fat_env, rho_r, rho_g, v_kin_pot_frag, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), pointer :: rho_r(:)
    type(pw_c1d_gs_type), pointer :: rho_g(:)
    type(pw_r3d_rs_type), pointer :: v_kin_pot_frag(:)
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_v_ts_frag"

    type(qs_environment_type), pointer :: qs_env
    type(section_vals_type), pointer :: input
    type(section_vals_type), pointer :: kin_section
    type(pw_r3d_rs_type), pointer :: tau(:)
    type(pw_r3d_rs_type), pointer :: v_tau(:)
    real(dp) :: virial(3, 3)
    real(dp) :: e_kin_frag
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    ! Get XC section for kinetic potential calculation
    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, input=input)
    kin_section => section_vals_get_subs_vals(input, "DFT%QS%OPT_EMBED%XC")

    ! Calculate fragment kinetic potential and energy
    nullify(v_kin_pot_frag, v_tau)
    call xc_vxc_pw_create(vxc_rho=v_kin_pot_frag, vxc_tau=v_tau, rho_r=rho_r, rho_g=rho_g, &
                          tau=tau, exc=e_kin_frag, xc_section=kin_section, &
                          pw_pool=fat_env%pw_pool_sub(i_frag)%pool, &
                          compute_virial=.false., virial_xc=virial)

    ! Store the fragment kinetic energy contribution
    fat_env%ekin_sub(i_frag) = e_kin_frag

    call timestop(handle)

  end subroutine fat_v_ts_frag

! **************************************************************************************************
!> \brief Calculate reference kinetic potential for the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_ts_ref(fat_env, v_kin_pot_ref)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), pointer :: v_kin_pot_ref(:)

    character(len=*), parameter :: routinen = "fat_v_ts_ref"

    type(qs_environment_type), pointer :: qs_env
    type(section_vals_type), pointer :: input
    type(section_vals_type), pointer :: kin_section
    type(pw_r3d_rs_type), pointer :: rho_r(:)
    type(pw_c1d_gs_type), pointer :: rho_g(:)
    type(pw_r3d_rs_type), pointer :: tau(:)
    type(pw_r3d_rs_type), pointer :: v_tau(:)
    real(dp) :: virial(3, 3)
    real(dp) :: e_kin_ref
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    ! Get XC section for kinetic potential calculation
    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, input=input)
    kin_section => section_vals_get_subs_vals(input, "DFT%QS%OPT_EMBED%XC")

    ! Create pws for total density
    call fat_pws_create_ref(fat_env, pws_r=rho_r, pws_g=rho_g, n_pws=1)

    ! Get density (rho_r is for LDA, rho_g is for GGA)
    call fat_get_rho_ref(fat_env, rho_r=rho_r, rho_g=rho_g)

    ! Calculate total system kinetic potential and energy
    nullify(v_kin_pot_ref, v_tau)
    call xc_vxc_pw_create(vxc_rho=v_kin_pot_ref, vxc_tau=v_tau, rho_r=rho_r, rho_g=rho_g, &
                          tau=tau, exc=e_kin_ref, xc_section=kin_section, &
                          pw_pool=fat_env%pw_pool_ref, &
                          compute_virial=.false., virial_xc=virial)

    ! Store the total kinetic energy contribution
    fat_env%ekin_ref = e_kin_ref

    ! Cleanup
    call fat_pws_release_ref(fat_env, pws_r=rho_r, pws_g=rho_g)

    call timestop(handle)

  end subroutine fat_v_ts_ref

! **************************************************************************************************
!> \brief Calculate xc potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_xc_frag(fat_env, rho_r, rho_g, v_xc_pot_frag, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), pointer :: rho_r(:)
    type(pw_c1d_gs_type), pointer :: rho_g(:)
    type(pw_r3d_rs_type), pointer :: v_xc_pot_frag(:)
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_v_xc_frag"

    type(qs_environment_type), pointer :: qs_env
    type(section_vals_type), pointer :: input
    type(section_vals_type), pointer :: xc_section
    type(pw_r3d_rs_type), pointer :: tau(:)
    type(pw_r3d_rs_type), pointer :: v_tau(:)
    real(dp) :: virial(3, 3)
    real(dp) :: e_xc_frag
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    ! Get XC section for XC potential calculation
    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, input=input)
    xc_section => section_vals_get_subs_vals(input, "DFT%XC")

    ! Calculate fragment XC potential and energy
    nullify(v_xc_pot_frag, v_tau)
    call xc_vxc_pw_create(vxc_rho=v_xc_pot_frag, vxc_tau=v_tau, rho_r=rho_r, rho_g=rho_g, &
                          tau=tau, exc=e_xc_frag, xc_section=xc_section, &
                          pw_pool=fat_env%pw_pool_sub(i_frag)%pool, &
                          compute_virial=.false., virial_xc=virial)

    ! Store the fragment XC energy contribution
    fat_env%exc_sub(i_frag) = e_xc_frag

    call timestop(handle)

  end subroutine fat_v_xc_frag

! **************************************************************************************************
!> \brief Calculate reference xc potential for the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_xc_ref(fat_env, v_xc_pot_ref)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), pointer :: v_xc_pot_ref(:)

    character(len=*), parameter :: routinen = "fat_v_xc_ref"

    type(qs_environment_type), pointer :: qs_env
    type(section_vals_type), pointer :: input
    type(section_vals_type), pointer :: xc_section
    type(pw_r3d_rs_type), pointer :: rho_r(:)
    type(pw_c1d_gs_type), pointer :: rho_g(:)
    type(pw_r3d_rs_type), pointer :: tau(:)
    type(pw_r3d_rs_type), pointer :: v_tau(:)
    real(dp) :: virial(3, 3)
    real(dp) :: e_xc_ref
    integer ::  handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    ! Get XC section for XC potential calculation
    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, input=input)
    xc_section => section_vals_get_subs_vals(input, "DFT%XC")

    ! Create pws for total density
    call fat_pws_create_ref(fat_env, pws_r=rho_r, pws_g=rho_g, n_pws=1)

    ! Get density (rho_r is for LDA, rho_g is for GGA)
    call fat_get_rho_ref(fat_env, rho_r=rho_r, rho_g=rho_g)

    ! Calculate total system XC potential and energy
    nullify(v_xc_pot_ref, v_tau)
    call xc_vxc_pw_create(vxc_rho=v_xc_pot_ref, vxc_tau=v_tau, rho_r=rho_r, rho_g=rho_g, &
                          tau=tau, exc=e_xc_ref, xc_section=xc_section, &
                          pw_pool=fat_env%pw_pool_ref, &
                          compute_virial=.false., virial_xc=virial)

    ! Store the total XC energy contribution
    fat_env%exc_ref = e_xc_ref

    ! Cleanup
    call fat_pws_release_ref(fat_env, pws_r=rho_r, pws_g=rho_g)

    call timestop(handle)

  end subroutine fat_v_xc_ref

! **************************************************************************************************
!> \brief Calculate electronic potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_e_frag(fat_env, v_elec_pot, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), pointer :: v_elec_pot
    integer :: i_frag

    character(len=*), parameter :: routineN = "fat_v_e_frag"

    type(qs_environment_type), pointer :: qs_env
    type(pw_poisson_type), pointer :: poisson_env
    type(pw_env_type), pointer :: pw_env
    type(pw_c1d_gs_type), pointer :: rho_g(:)
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    call fat_get_rho_frag(fat_env, rho_g=rho_g, i_frag=i_frag)

    ! Solve Poisson equation
    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, pw_env=pw_env)
    call pw_env_get(pw_env, poisson_env=poisson_env)
    call pw_poisson_solve(poisson_env, density=rho_g(1), vhartree=v_elec_pot)

    call timestop(handle)

  end subroutine fat_v_e_frag

! **************************************************************************************************
!> \brief Calculate nuclear potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_n_frag(fat_env, v_nuc_pot, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), pointer :: v_nuc_pot
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_v_n_frag"

    type(qs_environment_type), pointer :: qs_env
    type(pw_poisson_type), pointer :: poisson_env
    type(pw_env_type), pointer :: pw_env
    type(pw_c1d_gs_type), pointer :: rhon_g
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    call fat_get_rhon_frag(fat_env, rhon_g=rhon_g, i_frag=i_frag)

    ! Solve Poisson equation
    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, pw_env=pw_env)
    call pw_env_get(pw_env, poisson_env=poisson_env)
    call pw_poisson_solve(poisson_env, density=rhon_g, vhartree=v_nuc_pot)

    call timestop(handle)

  end subroutine fat_v_n_frag

! **************************************************************************************************
!> \brief Calculate non-additive kinetic embedding potential for one fragment in the freeze-and-thaw
!>        environment
! **************************************************************************************************
  subroutine fat_vemb_tsna_frag(fat_env, rho_r, rho_g, v_emb, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), pointer :: rho_r(:)
    type(pw_c1d_gs_type), pointer :: rho_g(:)
    type(pw_r3d_rs_type), pointer :: v_emb
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_vemb_tsna_frag"

    type(pw_r3d_rs_type), pointer :: v(:)
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    call pw_zero(v_emb)

    ! Reference system
    nullify(v)
    call fat_v_ts_ref(fat_env, v_kin_pot_ref=v)
    call pw_axpy(pw1=v(1), pw2=v_emb, allow_noncompatible_grids=.true.)
    call fat_pws_release_ref(fat_env, pws_r=v)

    ! Fragment
    nullify(v)
    call fat_v_ts_frag(fat_env, rho_r=rho_r, rho_g=rho_g, &
                       v_kin_pot_frag=v, i_frag=i_frag)
    call pw_axpy(pw1=v(1), pw2=v_emb, alpha=-1.0_dp, allow_noncompatible_grids=.true.)
    call fat_pws_release_frag(fat_env, pws_r=v, i_frag=i_frag)

    call timestop(handle)

  end subroutine fat_vemb_tsna_frag

! **************************************************************************************************
!> \brief Calculate non-additive xc embedding potential for one fragment in the freeze-and-thaw
!>        environment
! **************************************************************************************************
  subroutine fat_vemb_xcna_frag(fat_env, rho_r, rho_g, v_emb, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), pointer :: rho_r(:)
    type(pw_c1d_gs_type), pointer :: rho_g(:)
    type(pw_r3d_rs_type), pointer :: v_emb
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_vemb_xcna_frag"

    type(pw_r3d_rs_type), pointer :: v(:)
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    call pw_zero(v_emb)

    ! Reference system
    nullify(v)
    call fat_v_xc_ref(fat_env, v_xc_pot_ref=v)
    call pw_axpy(pw1=v(1), pw2=v_emb, allow_noncompatible_grids=.true.)
    call fat_pws_release_ref(fat_env, pws_r=v)

    ! Fragment
    nullify(v)
    call fat_v_xc_frag(fat_env, rho_r=rho_r, rho_g=rho_g, &
                       v_xc_pot_frag=v, i_frag=i_frag)
    call pw_axpy(pw1=v(1), pw2=v_emb, alpha=-1.0_dp, allow_noncompatible_grids=.true.)
    call fat_pws_release_frag(fat_env, pws_r=v, i_frag=i_frag)

    call timestop(handle)

  end subroutine fat_vemb_xcna_frag

! **************************************************************************************************
!> \brief Calculate electronic embedding potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_vemb_e_frag(fat_env, v_emb, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), pointer :: v_emb
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_vemb_e_frag"

    type(pw_r3d_rs_type), pointer :: v
    integer :: i, handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    call pw_zero(v_emb)

    do i = 1, fat_env%n_frag

       if (i == i_frag) cycle

       call fat_pw_create_frag(fat_env, pw_r=v, i_frag=i)

       call fat_v_e_frag(fat_env, v_elec_pot=v, i_frag=i)

       call pw_axpy(v, v_emb, allow_noncompatible_grids=.true.)

       call fat_pw_release_frag(fat_env, pw_r=v, i_frag=i)
    end do

    call timestop(handle)

  end subroutine fat_vemb_e_frag
! **************************************************************************************************
!> \brief Calculate nuclear embedding potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_vemb_n_frag(fat_env, v_emb, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), pointer :: v_emb
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_vemb_n_frag"

    type(pw_r3d_rs_type), pointer :: v
    integer :: i, handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    call pw_zero(v_emb)

    do i = 1, fat_env%n_frag

       if (i == i_frag) cycle

       call fat_pw_create_frag(fat_env, pw_r=v, i_frag=i)

       call fat_v_n_frag(fat_env, v_nuc_pot=v, i_frag=i)

       call pw_axpy(v, v_emb, allow_noncompatible_grids=.true.)

       call fat_pw_release_frag(fat_env, pw_r=v, i_frag=i)
    end do

    call timestop(handle)

  end subroutine fat_vemb_n_frag

! **************************************************************************************************
!> \brief Write embedding potential to cube and dat files
!> \param fat_env The freeze-and-thaw environment
!> \param v_emb The embedding potential to write
!> \param work_tmp Temporary workspace for potential calculations
!> \param i_frag The fragment index
!> \par History
!>      09.2024 Created [Konstantin Komarov]
! **************************************************************************************************
  SUBROUTINE fat_write_embedding_potential(fat_env, v_emb, work_tmp, i_frag)
    TYPE(fat_env_type), POINTER                       :: fat_env
    TYPE(pw_r3d_rs_type), POINTER                     :: v_emb
    TYPE(pw_r3d_rs_type), POINTER                     :: work_tmp
    INTEGER, INTENT(IN)                               :: i_frag

    CHARACTER(len=*), PARAMETER :: routineN = 'fat_write_embedding_potential'

    TYPE(qs_environment_type), POINTER                :: qs_env
    CHARACTER(len=default_string_length)              :: label
    TYPE(particle_list_type), POINTER                 :: particle_list
    TYPE(qs_subsys_type), POINTER                     :: subsys
    INTEGER, POINTER                                  :: stride(:)
    INTEGER                                           :: output_unit
    INTEGER                                           :: handle

    CALL timeset(routineN, handle)

    ! Proceed if it is the source node
    IF (.NOT. fat_env%force_env_ref%para_env%is_source()) THEN
       CALL timestop(handle)
       RETURN
    END IF

    IF (is_debug) PRINT *, "DBG| ", moduleN//"::"//routineN

    ! Convert fragment index to string for file naming
    CALL integer_to_string(i_frag, label)

    ! Write embedding potential to cube file
    CALL open_file( &
         "vemb_"//TRIM(label)//".cube", file_status="REPLACE", &
         file_action="WRITE", unit_number=output_unit)

    qs_env => fat_env%force_env_ref%qs_env
    CALL get_qs_env(qs_env, subsys=subsys)
    CALL qs_subsys_get(subsys, particles=particle_list)

    ALLOCATE(stride(3))
    stride = [1, 1, 1]
    CALL cp_pw_to_cube(v_emb, output_unit, "vemb_"//TRIM(label), particle_list, stride=stride)
    DEALLOCATE(stride)
    CALL close_file(unit_number=output_unit)

    ! Write embedding potential to dat file (simple volumetric format)
    CALL open_file( &
         "vemb_"//TRIM(label)//".dat", file_status="REPLACE", &
         file_action="WRITE", unit_number=output_unit)

    CALL pw_zero(work_tmp)

    ! Set weights (work_tmp) to the volume element of the grid
    CALL pw_set(work_tmp, &
                v_emb%pw_grid%dr(1) * &
                v_emb%pw_grid%dr(2) * &
                v_emb%pw_grid%dr(3))

    ! Write the potential values along with their weights (work_tmp)
    CALL pw_to_simple_volumetric(work_tmp, output_unit, [1, 1, 1], v_emb)

    CALL close_file(unit_number=output_unit)

    CALL timestop(handle)

  END SUBROUTINE fat_write_embedding_potential

  ! **************************************************************************************************
  !> \brief Execute external calculation script
  !> \param label Fragment label/identifier
  ! **************************************************************************************************
  SUBROUTINE fat_run_external(i_frag, gradient)
    INTEGER, INTENT(IN) :: i_frag
    LOGICAL, OPTIONAL :: gradient

    CHARACTER(len=*), PARAMETER :: routineN = 'fat_run_external'

    character(default_string_length) :: label
    CHARACTER :: cmdmsg(100)
    INTEGER :: exitstat
    INTEGER :: cmdstat
    logical :: my_grad = .FALSE.

    IF (is_debug) PRINT *, "DBG| ", moduleN//"::"//routineN

    if (PRESENT(gradient)) my_grad = gradient

    CALL integer_to_string(i_frag, label)
    if (my_grad) then
      CALL execute_command_line("sh extern_grad_"//TRIM(label)//".sh", .TRUE., exitstat, cmdstat, cmdmsg)
    else
      CALL execute_command_line("sh extern_"//TRIM(label)//".sh", .TRUE., exitstat, cmdstat, cmdmsg)
    end if

    IF (cmdstat /= 0) CPABORT("External command did not execute successfully")

  END SUBROUTINE fat_run_external

! **************************************************************************************************
!> \brief Read energy from external calculation
!> \param fat_env The freeze-and-thaw environment
!> \param i_frag Fragment index
! **************************************************************************************************
  SUBROUTINE fat_read_external_energy(fat_env, i_frag)
    TYPE(fat_env_type), POINTER :: fat_env
    INTEGER, INTENT(IN) :: i_frag

    CHARACTER(len=*), PARAMETER :: routineN = 'fat_read_external_energy'

    CHARACTER(len=default_string_length) :: label
    INTEGER :: input_unit

    IF (is_debug) PRINT *, "DBG| ", moduleN//"::"//routineN

    CALL integer_to_string(i_frag, label)
    CALL open_file("extern_"//TRIM(label)//".e", file_status="OLD", &
                  file_action="READ", unit_number=input_unit)
    READ(input_unit, *) fat_env%e_sub(i_frag)
    CALL close_file(unit_number=input_unit)

  END SUBROUTINE fat_read_external_energy

! **************************************************************************************************
!> \brief Calculate embedding potential for one fragment in the freeze-and-thaw environment
!> \param fat_env The freeze-and-thaw environment
!> \param v_emb The embedding potential (output)
!> \param i_frag The index of the current fragment
!> \param to_file Boolean flag to determine if the potential should be written to a file
! **************************************************************************************************
  subroutine fat_vemb_frag(fat_env, v_emb, i_frag, to_file)
    type(fat_env_type), pointer :: fat_env
    type(pw_r3d_rs_type), pointer :: v_emb
    integer :: i_frag
    logical :: to_file

    character(len=*), parameter :: routinen = "fat_vemb_frag"

    type(pw_r3d_rs_type), pointer :: rho_r(:)
    type(pw_c1d_gs_type), pointer :: rho_g(:)
    type(pw_r3d_rs_type), pointer :: vemb_tmp
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    call fat_get_rho_frag(fat_env, rho_r=rho_r, rho_g=rho_g, i_frag=i_frag)

    call fat_pw_create_frag(fat_env, pw_r=vemb_tmp, i_frag=i_frag)

    ! Calculate and add the non-additive kinetic embedding potential
    call fat_vemb_tsna_frag(fat_env, rho_r=rho_r, rho_g=rho_g, &
                            v_emb=vemb_tmp, i_frag=i_frag)
    call pw_axpy(vemb_tmp, v_emb)

    ! Calculate and add the non-additive exchange-correlation embedding potential
    call fat_vemb_xcna_frag(fat_env, rho_r=rho_r, rho_g=rho_g, &
                            v_emb=vemb_tmp, i_frag=i_frag)
    call pw_axpy(vemb_tmp, v_emb)

    ! Calculate and add the electron-electron embedding potential
    call fat_vemb_e_frag(fat_env, vemb_tmp, i_frag)
    call pw_axpy(vemb_tmp, v_emb)

    ! Calculate and add the nuclear-electron embedding potential
    call fat_vemb_n_frag(fat_env, vemb_tmp, i_frag)
    call pw_axpy(vemb_tmp, v_emb)

    ! If requested, write the embedding potential to files
    IF (to_file) THEN
       CALL fat_write_embedding_potential(fat_env, v_emb=v_emb, &
                                      work_tmp=vemb_tmp, i_frag=i_frag)
    END IF

    call fat_pw_release_frag(fat_env, pw_r=vemb_tmp, i_frag=i_frag)

    call timestop(handle)

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

    type(pw_r3d_rs_type), pointer :: rho_r
    type(pw_r3d_rs_type), pointer :: rho_r_tmp(:)
    type(pw_r3d_rs_type), pointer :: v_ee
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    call fat_pw_create_frag(fat_env, pw_r=v_ee, i_frag=i)
    call fat_pw_create_frag(fat_env, pw_r=rho_r, i_frag=i)

    call fat_v_e_frag(fat_env, v_ee, i)

    call fat_get_rho_frag(fat_env, rho_r=rho_r_tmp, i_frag=j)

    call pw_transfer(rho_r_tmp(1), rho_r)

    e = pw_integral_ab(v_ee, rho_r)

    call fat_pw_release_frag(fat_env, pw_r=v_ee, i_frag=i)
    call fat_pw_release_frag(fat_env, pw_r=rho_r, i_frag=i)

    call timestop(handle)

  end function fat_eint_ee

! **************************************************************************************************
!> \brief Calculate the ne interaction energy for two fragments of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_eint_ne(fat_env, i, j) result(e)
    type(fat_env_type), pointer :: fat_env
    integer :: i, j
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_eint_ne"

    type(pw_r3d_rs_type), pointer :: rho_r
    type(pw_r3d_rs_type), pointer :: rho_r_tmp(:)
    type(pw_r3d_rs_type), pointer :: v_ne
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen
    call timeset(routinen, handle)

    call fat_pw_create_frag(fat_env, pw_r=v_ne, i_frag=i)
    call fat_pw_create_frag(fat_env, pw_r=rho_r, i_frag=i)

    call fat_v_n_frag(fat_env, v_ne, i)

    call fat_get_rho_frag(fat_env, rho_r=rho_r_tmp, i_frag=j)

    call pw_transfer(rho_r_tmp(1), rho_r)

    e = pw_integral_ab(v_ne, rho_r)

    call fat_pw_release_frag(fat_env, pw_r=v_ne, i_frag=i)
    call fat_pw_release_frag(fat_env, pw_r=rho_r, i_frag=i)

    call timestop(handle)

  end function fat_eint_ne

! **************************************************************************************************
!> \brief Calculate the nn interaction energy for two fragments of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_eint_nn(fat_env, i, j) result(e)
    type(fat_env_type), pointer :: fat_env
    integer :: i, j
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_eint_nn"

    type(pw_r3d_rs_type), pointer :: v_nn
    type(pw_r3d_rs_type), pointer :: rhon_r
    type(pw_c1d_gs_type), pointer :: rhon_g_tmp
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    call fat_pw_create_frag(fat_env, pw_r=v_nn, i_frag=i)
    call fat_pw_create_frag(fat_env, pw_r=rhon_r, i_frag=i)

    call fat_v_n_frag(fat_env, v_nn, i)

    call fat_get_rhon_frag(fat_env, rhon_g=rhon_g_tmp, i_frag=j)

!   if (is_debug) call print_pw_statistics(rho, "rho", i)
!   if (is_debug) call print_pw_statistics(rhon_g_tmp, "rho_g_tmp", j)

    call pw_transfer(rhon_g_tmp, rhon_r)

    e = pw_integral_ab(v_nn, rhon_r)

    call fat_pw_release_frag(fat_env, pw_r=v_nn, i_frag=i)
    call fat_pw_release_frag(fat_env, pw_r=rhon_r, i_frag=i)

    call timestop(handle)

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

    type(dft_control_type), pointer :: dft_control
    type(force_env_type), pointer :: force_env_tmp
    type(pw_r3d_rs_type), pointer :: rho(:)
    type(pw_r3d_rs_type), pointer :: v_emb
    type(pw_r3d_rs_type), pointer :: weights
    type(qs_energy_type), pointer :: energy
    type(qs_environment_type), pointer :: qs_env_tmp
    integer :: handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen, " for frag:", i_frag

    call timeset(routinen, handle)

    ! Get the force environment and system environment for current fragment
    force_env_tmp => fat_env%force_env_sub(i_frag)%force_env
    qs_env_tmp => force_env_tmp%qs_env

    if (is_debug) &
      print *, "KKK| ", modulen//"::"//routinen, "fat_env%i_iter:", fat_env%i_iter

    ! Get the DFT control settings and enable embedding potential
    call get_qs_env(qs_env_tmp, dft_control=dft_control)
    dft_control%apply_embed_pot = .true.

!   ! Create a new embedding potential for current fragment
    call fat_pw_create_frag(fat_env, pw_r=v_emb, i_frag=i_frag)

    ! After first iteration, calculate the embedding potential
    if (fat_env%i_iter > 0) then
       ! Calculate total embedding potential including:
       ! - Non-additive Kinetic potential
       ! - Non-additive XC potential
       ! - Electrostatic potential from other fragments
       ! (to_file = False)
       call fat_vemb_frag(fat_env, v_emb, i_frag, .false.)
    end if

    ! Set the embedding potential to the qs env
    call set_qs_env(qs_env=qs_env_tmp, embed_pot=v_emb)


    ! Handle external fragment calculations if configured
    if (fat_env%extern(i_frag)) then

       if (is_debug) print *, "KKK| Using extern for fragment ", i_frag

       call fat_pw_create_frag(fat_env, pw_r=weights, i_frag=i_frag)
       CALL fat_write_embedding_potential(fat_env, v_emb=v_emb, &
                                      work_tmp=weights, i_frag=i_frag)
       call fat_pw_release_frag(fat_env, pw_r=weights, i_frag=i_frag)

       ! Execute external calculation script
       CALL fat_run_external(i_frag=i_frag)

       ! Get energy from external calculation
       CALL fat_read_external_energy(fat_env, i_frag)

       ! Get electron density from external calculation
       call fat_get_rho_frag(fat_env, rho_r=rho, i_frag=i_frag)

       ! Calculate embedding energy for current fragment:
       ! Integrate the density * embedding potential
       fat_env%e_emb(i_frag) = pw_integral_ab(v_emb, rho(1))

!      if (is_debug) then
!         ! Create temporary pw for density
!         call fat_pw_create_frag(fat_env, pw_r=rho_tmp, i_frag=i_frag)
!         call pw_transfer(rho(1), rho_tmp)
!         ! Write density to cube file
!         write(filename,'(A,I0,A,I0,A)') 'rho_frag', i_frag, '_iter', fat_env%i_iter, '.cube'
!         call write_pw_to_cube(fat_env, pw_r=rho_tmp, filename)
!         ! Deallocate
!         call fat_pw_release_frag(fat_env, pw_r=rho_tmp, i_frag=i_frag)
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

!      ! Get electron density from internal calculation
!      call fat_get_rho_frag(fat_env, rho_r=rho, i_frag=i_frag)

!      if (is_debug) then
!         ! Create temporary pw for density
!         call fat_pw_create_frag(fat_env, pw_r=rho_tmp, i_frag=i_frag)
!         call pw_transfer(rho(1), rho_tmp)
!         ! Write density to cube file
!         write(filename,'(A,I0,A,I0,A)') 'rho_frag', i_frag, '_iter', fat_env%i_iter, '.cube'
!         call write_pw_to_cube(fat_env, pw_r=rho_tmp, filename)
!         ! Deallocate
!         call fat_pw_release_frag(fat_env, pw_r=rho_tmp, i_frag=i_frag)
!      end if

    end if

    if (is_debug) then
       print *, "KKK| Total subsystem (", i_frag, ") energy:", fat_env%e_sub(i_frag)
    end if

    ! Deallocate v_emb potential
    call fat_pw_release_frag(fat_env, pw_r=v_emb, i_frag=i_frag)

    call timestop(handle)

  end subroutine fat_iterate_frag

! **************************************************************************************************
!> \brief Performs one iteration of the freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_iterate(fat_env)
    type(fat_env_type), pointer :: fat_env

    character(len=*), parameter :: routinen = "fat_iterate"

    integer :: i_frag, handle

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    ! Loop over all fragments to optimize their densities with embedding effects
    do i_frag = 1, fat_env%n_frag

       ! Calculate new density and energy for fragment i_frag under embedding potential from others
       call fat_iterate_frag(fat_env, i_frag)

       ! Print current optimization status
       call fat_status(fat_env, i_frag)
    end do

    ! Increment iteration counter
    fat_env%i_iter = fat_env%i_iter + 1

    call timestop(handle)

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
    forces = 0.0_dp

  end subroutine fat_force_create

! **************************************************************************************************
!> \brief Calculate the self force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_self(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_self"

    integer :: i, i_atom, i_kind, j, k, n, handle
    integer, allocatable :: atom_of_kind(:), kind_of(:)
    type(atomic_kind_type), pointer :: kinds(:)
    type(qs_environment_type), pointer :: qs_env
    type(qs_force_type), pointer :: qs_force(:)

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    forces = 0.0_dp

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

    call timestop(handle)

  end subroutine fat_force_self

! **************************************************************************************************
!> \brief Calculate the embedding force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_emb(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_emb"

    integer :: a, i, j, k, handle
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_r3d_rs_type), pointer :: v_r
    type(pw_r3d_rs_type), pointer :: drho_r
    type(pw_c1d_gs_type), pointer :: drho_g

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    forces = 0.0_dp

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, pw_r=v_r, i_frag=i)
       call fat_vemb_frag(fat_env, v_r, i, .true.)
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)

       call fat_pw_create_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_get_drho_frag(fat_env, drho_r=drho_r, drho_g=drho_g, i_frag=i, i_nuc=j, i_dir=a)
             forces(k, a) = pw_integral_ab(v_r, drho_r)
          end do
       end do
       call fat_pw_release_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)
       call fat_pw_release_frag(fat_env, pw_r=v_r, i_frag=i)
    end do

    call timestop(handle)

  end subroutine fat_force_emb

! **************************************************************************************************
!> \brief Calculate the non-additive kinetic force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_tsna(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_tsna"

    integer :: a, i, j, k, handle
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_r3d_rs_type), pointer :: rho_r(:)
    type(pw_c1d_gs_type), pointer :: rho_g(:)
    type(pw_r3d_rs_type), pointer :: v_emb
    type(pw_r3d_rs_type), pointer :: drho_r
    type(pw_c1d_gs_type), pointer :: drho_g

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call timeset(routinen, handle)

    forces = 0.0_dp

    k = 0
    do i = 1, fat_env%n_frag  ! fragments

      ! Create embedding potential
       call fat_pw_create_frag(fat_env, pw_r=v_emb, i_frag=i)

       call fat_get_rho_frag(fat_env, rho_r=rho_r, rho_g=rho_g, i_frag=i)

       ! Calculate total embedding potential (includes both LDA and GGA if present)
       CALL fat_vemb_tsna_frag(fat_env=fat_env, rho_r=rho_r, rho_g=rho_g, &
                               v_emb=v_emb, i_frag=i)

       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)

       call fat_pw_create_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz

             call fat_get_drho_frag(fat_env, drho_r=drho_r, drho_g=drho_g, i_frag=i, i_nuc=j, i_dir=a)

             ! Total force contribution (LDA + GGA if present)
             forces(k, a) = pw_integral_ab(v_emb, drho_r)

          end do
       end do
       call fat_pw_release_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)

       ! Release potential
       call fat_pw_release_frag(fat_env, pw_r=v_emb, i_frag=i)

    end do

    call timestop(handle)

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
    type(pw_r3d_rs_type), pointer :: v_emb
    type(pw_r3d_rs_type), pointer :: rho_r(:)
    type(pw_c1d_gs_type), pointer :: rho_g(:)
    type(pw_r3d_rs_type), pointer :: drho_r
    type(pw_c1d_gs_type), pointer :: drho_g

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    forces = 0.0_dp

    k = 0
    do i = 1, fat_env%n_frag  ! fragments

      ! Create embedding potential
       call fat_pw_create_frag(fat_env, pw_r=v_emb, i_frag=i)

       call fat_get_rho_frag(fat_env, rho_r=rho_r, rho_g=rho_g, i_frag=i)

       ! Calculate total embedding potential (includes both LDA and GGA if present)
       call fat_vemb_xcna_frag(fat_env=fat_env, rho_r=rho_r, rho_g=rho_g, &
                               v_emb=v_emb, i_frag=i)

       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)

       call fat_pw_create_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_get_drho_frag(fat_env, drho_r=drho_r, drho_g=drho_g, i_frag=i, i_nuc=j, i_dir=a)

             ! Total force contribution (LDA + GGA if present)
             forces(k, a) = pw_integral_ab(v_emb, drho_r)

          end do
       end do
       call fat_pw_release_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)

       ! Release potential
       call fat_pw_release_frag(fat_env, pw_r=v_emb, i_frag=i)

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
    type(pw_r3d_rs_type), pointer :: v_r
    type(pw_r3d_rs_type), pointer :: drho_r
    type(pw_c1d_gs_type), pointer :: drho_g

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    forces = 0.0_dp

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, pw_r=v_r, i_frag=i)
       call fat_vemb_e_frag(fat_env, v_r, i)
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)
       call fat_pw_create_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_get_drho_frag(fat_env, drho_r=drho_r, drho_g=drho_g, i_frag=i, i_nuc=j, i_dir=a)
             forces(k, a) = pw_integral_ab(v_r, drho_r)
          end do
       end do
       call fat_pw_release_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)
       call fat_pw_release_frag(fat_env, pw_r=v_r, i_frag=i)
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
    type(pw_r3d_rs_type), pointer :: ve_r
    type(pw_r3d_rs_type), pointer :: vn_r
    type(pw_r3d_rs_type), pointer :: drho_r
    TYPE(pw_r3d_rs_type), POINTER :: rho_r_tmp
    type(pw_c1d_gs_type), pointer :: drho_g
    type(pw_c1d_gs_type), pointer :: drhon_g

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    forces = 0.0_dp

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, pw_r=ve_r, i_frag=i)
       call fat_vemb_e_frag(fat_env, ve_r, i)

       call fat_pw_create_frag(fat_env, pw_r=vn_r, i_frag=i)
       call fat_vemb_n_frag(fat_env, vn_r, i)

       qs_env => fat_env%force_env_sub(i)%force_env%qs_env

       call get_qs_env(qs_env, particle_set=particles)
       call fat_pw_create_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)
       call fat_pw_create_frag(fat_env, pw_r=rho_r_tmp, pw_g=drhon_g, i_frag=i)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_get_drho_frag(fat_env, drho_r=drho_r, drho_g=drho_g, i_frag=i, i_nuc=j, i_dir=a)
             CALL pw_transfer(drho_g, rho_r_tmp)
             forces(k, a) = pw_integral_ab(vn_r, rho_r_tmp)

             call fat_get_drhon_frag(fat_env, drhon_g=drhon_g, i_frag=i, i_nuc=j, i_dir=a)
             CALL pw_transfer(drhon_g, rho_r_tmp)
             forces(k, a) = pw_integral_ab(ve_r, rho_r_tmp)
          end do
       end do
       call fat_pw_release_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)
       call fat_pw_release_frag(fat_env, pw_r=rho_r_tmp, pw_g=drhon_g, i_frag=i)

       call fat_pw_release_frag(fat_env, pw_r=ve_r, i_frag=i)
       call fat_pw_release_frag(fat_env, pw_r=vn_r, i_frag=i)
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
    type(pw_r3d_rs_type), pointer :: v_pot
    type(pw_r3d_rs_type), pointer :: rho_r_tmp
    type(pw_c1d_gs_type), pointer :: drho_g

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    forces = 0.0_dp

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, pw_r=v_pot, i_frag=i)
       call fat_vemb_n_frag(fat_env, v_pot, i)
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)
       call fat_pw_create_frag(fat_env, pw_r=rho_r_tmp, pw_g=drho_g, i_frag=i)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_get_drhon_frag(fat_env, drhon_g=drho_g, i_frag=i, i_nuc=j, i_dir=a)
             call pw_transfer(drho_g, rho_r_tmp)
             forces(k, a) = pw_integral_ab(v_pot, rho_r_tmp)
          end do
       end do
       call fat_pw_release_frag(fat_env, pw_r=rho_r_tmp, pw_g=drho_g, i_frag=i)
       call fat_pw_release_frag(fat_env, pw_r=v_pot, i_frag=i)
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

    forces = 0.0_dp

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
   SUBROUTINE fat_force_print(fat_env, forces, label)
      TYPE(fat_env_type), POINTER                        :: fat_env
      REAL(dp), ALLOCATABLE, INTENT(IN)                  :: forces(:, :)
      CHARACTER(len=*), INTENT(IN)                       :: label

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_force_print'


      INTEGER :: i, n, output_unit, log_unit
      CHARACTER(len=2) :: symbol
      CHARACTER(len=*), PARAMETER :: fmt_energy = &
        '(/," ENERGY| Total FORCE_EVAL ( FAT ) energy [a.u.]:",T61,F20.15)'
      CHARACTER(len=*), PARAMETER :: fmt_header = &
        '(/," ATOMIC FORCES in [a.u.]",/,/," # Atom   Kind   Element",T30,"X",T45,"Y",T60,"Z")'
      CHARACTER(len=*), PARAMETER :: fmt_forces = &
        '(2X,I5,2X,I5,6X,A2,T23,3F15.8)'
      CHARACTER(len=*), PARAMETER :: fmt_sum = &
        '(/," SUM OF ATOMIC FORCES",T23,4F15.8,/)'
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      TYPE(atomic_kind_type), POINTER                    :: kind
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(3)                        :: grand_total
      TYPE(cp_logger_type), POINTER                      :: logger

      IF (is_debug) PRINT *, "DBG| ", moduleN//"::"//routineN

      qs_env => fat_env%force_env_ref%qs_env
      CALL get_qs_env(qs_env, particle_set=particles)
      n = SIZE(particles)

      logger => cp_get_default_logger()
      log_unit = cp_logger_get_default_io_unit(logger)

      IF (log_unit > 0) THEN
         WRITE(log_unit, fmt_energy) fat_env%e_tot
         WRITE(log_unit, fmt_header)

         grand_total(:) = 0.0_dp
         DO i = 1, n
            kind => particles(i)%atomic_kind
            CALL get_atomic_kind(kind, element_symbol=symbol)
            WRITE(log_unit, fmt_forces) i, i, symbol, forces(i,1:3)
            grand_total(:) = grand_total(:) + forces(i,:)
         END DO

         WRITE(log_unit, fmt_sum) grand_total(:), SQRT(SUM(grand_total**2))
      END IF

      ! Write forces to XYZ file
      CALL open_file(label//".xyz", file_status="REPLACE", &
                     file_action="WRITE", unit_number=output_unit)
      WRITE (output_unit, '(I0)') n
      WRITE (output_unit, *)
      DO i = 1, n
         kind => particles(i)%atomic_kind
         CALL get_atomic_kind(kind, element_symbol=symbol)
         WRITE (output_unit, '(A2,3(F17.12," "))') &
            symbol, forces(i, 1), forces(i, 2), forces(i, 3)
      END DO
      CALL close_file(output_unit)

   END SUBROUTINE fat_force_print

! **************************************************************************************************
!> \brief Calculate the forces for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force(fat_env)
    type(fat_env_type), pointer :: fat_env

    character(len=*), parameter :: routinen = "fat_force"

    integer :: i, j
    real(kind=dp), allocatable :: forces(:, :)
    real(kind=dp), allocatable :: forces_tmp(:, :)

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    call fat_force_create(fat_env, forces)
    call fat_force_create(fat_env, forces_tmp)

    call fat_force_self(fat_env, forces_tmp)
    call fat_force_print(fat_env, forces_tmp, "force_self")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_tmp(i, j)
       end do
    end do

    call fat_force_emb(fat_env, forces_tmp)
    call fat_force_print(fat_env, forces_tmp, "force_emb")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) - forces_tmp(i, j)
       end do
    end do

    call fat_force_tsna(fat_env, forces_tmp)
    call fat_force_print(fat_env, forces_tmp, "force_tsna")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_tmp(i, j)
       end do
    end do

    call fat_force_xcna(fat_env, forces_tmp)
    call fat_force_print(fat_env, forces_tmp, "force_xcna")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_tmp(i, j)
       end do
    end do

    call fat_force_ee(fat_env, forces_tmp)
    call fat_force_print(fat_env, forces_tmp, "force_ee")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_tmp(i, j)
       end do
    end do

    call fat_force_ne(fat_env, forces_tmp)
    call fat_force_print(fat_env, forces_tmp, "force_ne")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + 2*forces_tmp(i, j)
       end do
    end do

    call fat_force_nn(fat_env, forces_tmp)
    call fat_force_print(fat_env, forces_tmp, "force_nn")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_tmp(i, j)
       end do
    end do

    call fat_force_overlap(fat_env, forces_tmp)
    call fat_force_print(fat_env, forces_tmp, "force_overlap")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_tmp(i, j)
       end do
    end do
    deallocate(forces_tmp)

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
    real(kind=dp) :: e_tmp

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

    type(dft_control_type), pointer :: dft_control
    type(qs_environment_type), pointer :: qs_env

    integer :: i

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
    real(kind=dp) :: eps
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
