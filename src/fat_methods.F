!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE fat_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_cube_to_pw,&
                                              cp_pw_to_cube,&
                                              cp_pw_to_simple_volumetric
   USE cp_subsys_types,                 ONLY: cp_subsys_type,&
                                              cp_subsys_get
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE fat_types,                       ONLY: fat_env_create,&
                                              fat_env_release,&
                                              fat_env_type
   USE force_env_types,                 ONLY: force_env_type,&
                                              force_env_get
   USE hartree_local_types,             ONLY: ecoul_1center_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: HALFSPACE,&
                                              PW_MODE_DISTRIBUTED,&
                                              pw_grid_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integral_ab,&
                                              pw_scale,&
                                              pw_set,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE scf_control_types,               ONLY: scf_control_type
   USE qs_collocate_density,            ONLY: calculate_drho_core,&
                                              calculate_drho_elec_dR,&
                                              calculate_ppl_grid
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force,                        ONLY: qs_calc_energy_force
   USE qs_force_types,                  ONLY: qs_force_type,&
                                              allocate_qs_force,&
                                              zero_qs_force
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_local_rho_types,              ONLY: local_rho_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_set,&
                                              qs_subsys_type
   USE realspace_grid_cube,             ONLY: pw_to_simple_volumetric
   USE string_utilities,                ONLY: integer_to_string
   USE xc,                              ONLY: xc_vxc_pw_create
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: fat_embedding

   LOGICAL, PRIVATE, PARAMETER :: is_debug = .FALSE.

   INTEGER, SAVE :: log_counter = 0

   CHARACTER(len=*), PRIVATE, PARAMETER :: &
      modulen = "fat_methods", &
      fat_title_fmt = '(/, T2, "FAT| ", 4(A6, ", "), 9(A26, ", "), A26)', &
      fat_data_fmt = '(/, T2, "FAT| ", 4(I6, ", "), 7(F26.15, ", "), F26.15)'

CONTAINS

! **************************************************************************************************
!> \brief Creates a pw for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_pw_create_frag(fat_env, pw_r, pw_g, rho0_s_pw, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), OPTIONAL, POINTER                   :: pw_r
      TYPE(pw_type), OPTIONAL, POINTER                   :: pw_g
      TYPE(pw_type), OPTIONAL, POINTER                   :: rho0_s_pw
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_pw_create_frag'

      INTEGER                                            :: handle

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) &
         PRINT *, "DBG| ", modulen//"::"//routinen

      CPASSERT(PRESENT(pw_r) .OR. PRESENT(pw_g))
      CPASSERT(ASSOCIATED(fat_env%pw_pool_sub))
      CPASSERT(i_frag >= 1 .AND. i_frag <= SIZE(fat_env%pw_pool_sub))
      CPASSERT(ASSOCIATED(fat_env%pw_pool_sub(i_frag)%pool))

      IF (PRESENT(pw_r)) THEN
         ALLOCATE (pw_r)
         CALL fat_env%pw_pool_sub(i_frag)%pool%create_pw(pw_r, &
                                                         use_data=realdata3d, &
                                                         in_space=realspace)
         CALL pw_zero(pw_r)
      END IF

      IF (PRESENT(pw_g)) THEN
         ALLOCATE (pw_g)
         CALL fat_env%pw_pool_sub(i_frag)%pool%create_pw(pw_g, &
                                                         use_data=complexdata1d, &
                                                         in_space=reciprocalspace)
         CALL pw_zero(pw_g)
      END IF

      IF (PRESENT(rho0_s_pw) .AND. fat_env%use_gapw) THEN
         ALLOCATE (rho0_s_pw)
         CALL fat_env%pw_pool_sub(i_frag)%pool%create_pw(rho0_s_pw, &
                                                         use_data=realdata3d, &
                                                         in_space=realspace)
         CALL pw_zero(rho0_s_pw)
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pw_create_frag

! **************************************************************************************************
!> \brief Releases a pw for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_pw_release_frag(fat_env, pw_r, pw_g, rho0_s_pw, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), OPTIONAL, POINTER                   :: pw_r
      TYPE(pw_type), OPTIONAL, POINTER                   :: pw_g
      TYPE(pw_type), OPTIONAL, POINTER                   :: rho0_s_pw
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_pw_release_frag'

      INTEGER                                            :: handle

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CALL timeset(routinen, handle)

      CPASSERT(PRESENT(pw_r) .OR. PRESENT(pw_g))

      IF (PRESENT(pw_r)) THEN
         CPASSERT(ASSOCIATED(pw_r))
         CALL fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(pw_r)
         DEALLOCATE (pw_r)
      END IF

      IF (PRESENT(pw_g)) THEN
         CPASSERT(ASSOCIATED(pw_g))
         CALL fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(pw_g)
         DEALLOCATE (pw_g)
      END IF

      IF (PRESENT(rho0_s_pw) .AND. fat_env%use_gapw) THEN
         CPASSERT(ASSOCIATED(rho0_s_pw))
         CALL fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(rho0_s_pw)
         DEALLOCATE (rho0_s_pw)
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pw_release_frag

! **************************************************************************************************
!> \brief Creates a pw for the reference system in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_pw_create_ref(fat_env, pw_r, pw_g)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), INTENT(out), OPTIONAL               :: pw_r
      TYPE(pw_type), INTENT(out), OPTIONAL               :: pw_g

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_pw_create_ref'

      INTEGER                                            :: handle

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CALL timeset(routinen, handle)

      CPASSERT(PRESENT(pw_r) .OR. PRESENT(pw_g))
      CPASSERT(ASSOCIATED(fat_env%pw_pool_ref))

      IF (PRESENT(pw_r)) THEN
         CALL fat_env%pw_pool_ref%create_pw(pw_r, &
                                            use_data=realdata3d, &
                                            in_space=realspace)
         CALL pw_zero(pw_r)
      END IF

      IF (PRESENT(pw_g)) THEN
         CALL fat_env%pw_pool_ref%create_pw(pw_g, &
                                            use_data=COMPLEXDATA1D, &
                                            in_space=RECIPROCALSPACE)
         CALL pw_zero(pw_g)
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pw_create_ref

! **************************************************************************************************
!> \brief Releases a pw for the reference system in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_pw_release_ref(fat_env, pw_r, pw_g)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), INTENT(inout), OPTIONAL             :: pw_r
      TYPE(pw_type), INTENT(inout), OPTIONAL             :: pw_g

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_pw_release_ref'

      INTEGER                                            :: handle

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CALL timeset(routinen, handle)

      CPASSERT(PRESENT(pw_r) .OR. PRESENT(pw_g))

      IF (PRESENT(pw_r)) THEN
         CALL fat_env%pw_pool_ref%give_back_pw(pw_r)
      END IF

      IF (PRESENT(pw_g)) THEN
         CALL fat_env%pw_pool_ref%give_back_pw(pw_g)
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pw_release_ref

! **************************************************************************************************
!> \brief Creates an array of pws for the reference system in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_pws_create_ref(fat_env, pws_r, pws_g, n_pws)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), OPTIONAL, POINTER                   :: pws_r(:)
      TYPE(pw_type), OPTIONAL, POINTER                   :: pws_g(:)
      INTEGER                                            :: n_pws

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_pws_create_ref'

      INTEGER                                            :: handle, i

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CALL timeset(routinen, handle)

      CPASSERT(PRESENT(pws_r) .OR. PRESENT(pws_g))

      IF (PRESENT(pws_r)) THEN
         ALLOCATE (pws_r(n_pws))
         DO i = 1, n_pws
            CALL fat_pw_create_ref(fat_env, pw_r=pws_r(i))
         END DO
      END IF

      IF (PRESENT(pws_g)) THEN
         ALLOCATE (pws_g(n_pws))
         DO i = 1, n_pws
            CALL fat_pw_create_ref(fat_env, pw_g=pws_g(i))
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pws_create_ref

! **************************************************************************************************
!> \brief Releases an array of pws for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_pws_release_ref(fat_env, pws_r, pws_g)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), OPTIONAL, POINTER                   :: pws_r(:)
      TYPE(pw_type), OPTIONAL, POINTER                   :: pws_g(:)

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_pws_release_ref'

      INTEGER                                            :: handle, i

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CALL timeset(routinen, handle)

      CPASSERT(PRESENT(pws_r) .OR. PRESENT(pws_g))

      IF (PRESENT(pws_r)) THEN
         CPASSERT(ASSOCIATED(pws_r))
         DO i = 1, SIZE(pws_r)
            CALL fat_pw_release_ref(fat_env, pw_r=pws_r(i))
         END DO
         DEALLOCATE (pws_r)
      END IF

      IF (PRESENT(pws_g)) THEN
         CPASSERT(ASSOCIATED(pws_g))
         DO i = 1, SIZE(pws_g)
            CALL fat_pw_release_ref(fat_env, pw_g=pws_g(i))
         END DO
         DEALLOCATE (pws_g)
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pws_release_ref

! **************************************************************************************************
!> \brief Creates an array of pws for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_pws_create_frag(fat_env, pws_r, pws_g, i_frag, n_pws)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), OPTIONAL, POINTER                   :: pws_r(:)
      TYPE(pw_type), OPTIONAL, POINTER                   :: pws_g(:)
      INTEGER, INTENT(IN)                                :: i_frag, n_pws

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_pws_create_frag'

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CPASSERT(PRESENT(pws_r) .OR. PRESENT(pws_g))

      IF (PRESENT(pws_r)) THEN
         ALLOCATE (pws_r(n_pws))
         DO i = 1, n_pws
            CALL fat_env%pw_pool_sub(i_frag)%pool%create_pw(pws_r(i), &
                                                            use_data=REALDATA3D, &
                                                            in_space=REALSPACE)
            CALL pw_zero(pws_r(i))
         END DO
      END IF

      IF (PRESENT(pws_g)) THEN
         ALLOCATE (pws_g(n_pws))
         DO i = 1, n_pws
            CALL fat_env%pw_pool_sub(i_frag)%pool%create_pw(pws_g(i), &
                                                            use_data=COMPLEXDATA1D, &
                                                            in_space=RECIPROCALSPACE)
            CALL pw_zero(pws_g(i))
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pws_create_frag

! **************************************************************************************************
!> \brief Releases an array of pws for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_pws_release_frag(fat_env, pws_r, pws_g, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), OPTIONAL, POINTER                   :: pws_r(:)
      TYPE(pw_type), OPTIONAL, POINTER                   :: pws_g(:)
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_pws_release_frag'

      INTEGER                                            :: handle, i
      TYPE(pw_type), POINTER                             :: pw_g_tmp
      TYPE(pw_type), POINTER                             :: pw_r_tmp

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CPASSERT(PRESENT(pws_r) .OR. PRESENT(pws_g))

      IF (PRESENT(pws_r)) THEN
         CPASSERT(ASSOCIATED(pws_r))
         DO i = 1, SIZE(pws_r)
            pw_r_tmp => pws_r(i)
            CALL fat_pw_release_frag(fat_env, pw_r=pw_r_tmp, i_frag=i_frag)
         END DO
      END IF

      IF (PRESENT(pws_g)) THEN
         CPASSERT(ASSOCIATED(pws_g))
         DO i = 1, SIZE(pws_g)
            pw_g_tmp => pws_g(i)
            CALL fat_pw_release_frag(fat_env, pw_g=pw_g_tmp, i_frag=i_frag)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pws_release_frag

! **************************************************************************************************
!> \brief Gets the electron density for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_get_rho_frag(fat_env, rho_r, rho_g, rho0_s_r, rho0_s_g, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), OPTIONAL, POINTER                   :: rho_r(:)
      TYPE(pw_type), OPTIONAL, POINTER                   :: rho_g(:)
      TYPE(pw_type), OPTIONAL, POINTER                   :: rho0_s_r
      TYPE(pw_type), OPTIONAL, POINTER                   :: rho0_s_g
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_get_rho_frag'

      INTEGER                                            :: handle
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: qs_rho
!     CHARACTER(len=default_string_length)               :: label
!     TYPE(pw_type), POINTER                             :: rho_r_tmp

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CPASSERT(PRESENT(rho_r) .OR. PRESENT(rho_g))

      qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
      CALL get_qs_env(qs_env, rho=qs_rho)

      IF (PRESENT(rho_r)) THEN
         NULLIFY (rho_r)
         CALL qs_rho_get(qs_rho, rho_r=rho_r)
!        IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) &
!           call print_pw_statistics(rho_r(1), "get_rho_r", i_frag)
      END IF

      IF (PRESENT(rho_g)) THEN
         NULLIFY (rho_g)
         CALL qs_rho_get(qs_rho, rho_g=rho_g)
!        IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) &
!           call print_pw_statistics(rho_g(1), "get_rho_g", i_frag)
      END IF

      IF (PRESENT(rho0_s_r) .AND. fat_env%use_gapw) THEN
         CALL get_qs_env(qs_env, rho0_s_rs=rho0_s_r)
      END IF

      IF (PRESENT(rho0_s_g) .AND. fat_env%use_gapw) THEN
         CALL get_qs_env(qs_env, rho0_s_gs=rho0_s_g)
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_get_rho_frag

! **************************************************************************************************
!> \brief Compute the total density (electronic + core) in real space for a fragment
!> \param fat_env The freeze-and-thaw environment
!> \param rho_tot_r The total density in real space (output)
!> \param i_frag The fragment index
! **************************************************************************************************
   SUBROUTINE fat_get_rho_tot_frag(fat_env, rho_tot_r, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: rho_tot_r
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_get_rho_tot_frag'

      INTEGER                                            :: handle, ispin
      TYPE(pw_type), POINTER                             :: rho0_s_rs
      TYPE(pw_type), POINTER                             :: rho_core_g
      TYPE(pw_type), POINTER                             :: rho_core_r
      TYPE(pw_type), POINTER                             :: rho_r(:)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CALL fat_get_rho_frag(fat_env, rho_r=rho_r, i_frag=i_frag)

      ! Sum electronic densities
      DO ispin = 1, SIZE(rho_r)
         CALL pw_axpy(rho_r(ispin), rho_tot_r)
      END DO

      ! Add core density
      IF (fat_env%use_gapw) THEN
         qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
         CALL get_qs_env(qs_env, rho0_s_rs=rho0_s_rs)
         CALL pw_axpy(rho0_s_rs, rho_tot_r)
      ELSE
         CALL fat_get_rhon_frag(fat_env, rhon_g=rho_core_g, i_frag=i_frag)
         CALL fat_pw_create_frag(fat_env, pw_r=rho_core_r, i_frag=i_frag)
         CALL pw_transfer(rho_core_g, rho_core_r)
         CALL pw_axpy(rho_core_r, rho_tot_r)
         CALL fat_pw_release_frag(fat_env, pw_r=rho_core_r, i_frag=i_frag)
      END IF

      CALL timestop(handle)
   END SUBROUTINE fat_get_rho_tot_frag

   SUBROUTINE write_pw_to_cube(fat_env, pw_r, pw_g, filename)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), OPTIONAL, POINTER                   :: pw_r
      TYPE(pw_type), OPTIONAL, POINTER                   :: pw_g
      CHARACTER(len=*), INTENT(in)                       :: filename

      INTEGER                                            :: unit_nr
      INTEGER, POINTER                                   :: stride(:)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(pw_type), POINTER                             :: pw_r_tmp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CPASSERT(PRESENT(pw_r) .OR. PRESENT(pw_g))

      logger => cp_get_default_logger()
      unit_nr = cp_logger_get_default_io_unit(logger)
      log_counter = log_counter + 1

      IF (log_counter > 6) RETURN

      PRINT *, "DBG| write CUBE file number:", filename

      IF (PRESENT(pw_r)) THEN
         IF (unit_nr > 0) THEN
            qs_env => fat_env%force_env_ref%qs_env

            CALL get_qs_env(qs_env=qs_env, &
                            dft_control=dft_control, &
                            pw_env=pw_env, &
                            particle_set=particle_set, &
                            subsys=subsys)
            CALL pw_env_get(pw_env, auxbas_pw_pool=pw_pool)
            CALL qs_subsys_get(subsys, particles=particles)

            ALLOCATE (stride(3))
            stride = [1, 1, 1]

            CALL open_file(file_name=TRIM(filename), &
                           file_status="REPLACE", &
                           file_action="WRITE", &
                           unit_number=unit_nr)
            CALL cp_pw_to_cube(pw_r, unit_nr, TRIM(filename), &
                               particles=particles, stride=stride)
            CALL close_file(unit_number=unit_nr)
            DEALLOCATE (stride)
         END IF
      END IF

      IF (PRESENT(pw_g)) THEN
         IF (unit_nr > 0) THEN
            CALL fat_pw_create_ref(fat_env, pw_r=pw_r_tmp)
            CALL pw_zero(pw_r_tmp)
            CALL pw_transfer(pw_g, pw_r_tmp)

            qs_env => fat_env%force_env_ref%qs_env

            CALL get_qs_env(qs_env=qs_env, &
                            dft_control=dft_control, &
                            pw_env=pw_env, &
                            particle_set=particle_set, &
                            subsys=subsys)
            CALL pw_env_get(pw_env, auxbas_pw_pool=pw_pool)
            CALL qs_subsys_get(subsys, particles=particles)

            ALLOCATE (stride(3))
            stride = [1, 1, 1]

            CALL open_file(file_name=TRIM(filename), &
                           file_status="REPLACE", &
                           file_action="WRITE", &
                           unit_number=unit_nr)
            CALL cp_pw_to_cube(pw_r_tmp, unit_nr, TRIM(filename), &
                               particles=particles, stride=stride)
            CALL close_file(unit_number=unit_nr)
            DEALLOCATE (stride)

            CALL fat_pw_release_ref(fat_env, pw_r=pw_r_tmp)
         END IF
      END IF

   END SUBROUTINE write_pw_to_cube

   SUBROUTINE print_pw_statistics(pw, name, i_frag)
      TYPE(pw_type), INTENT(in)                          :: pw
      CHARACTER(len=*), INTENT(in)                       :: name
      INTEGER, INTENT(in)                                :: i_frag

      INTEGER                                            :: count, i
      LOGICAL                                            :: print_grid = .TRUE.
      REAL(dp)                                           :: max_val, mean, min_val, norm_2, std_dev, &
                                                            sum_val
      TYPE(pw_grid_type), POINTER                        :: grid

      grid => pw%pw_grid

      IF (pw%in_space == REALSPACE) THEN
         sum_val = SUM(ABS(pw%cr3d))
         count = SIZE(pw%cr3d)
         mean = sum_val/count
         std_dev = SQRT(SUM((ABS(pw%cr3d) - mean)**2)/count)
         min_val = MINVAL(ABS(pw%cr3d))
         max_val = MAXVAL(ABS(pw%cr3d))
         norm_2 = SQRT(SUM(pw%cr3d**2))
      ELSE IF (pw%in_space == RECIPROCALSPACE) THEN
         sum_val = SUM(ABS(pw%cc))
         count = SIZE(pw%cc)
         mean = REAL(sum_val)/count
         std_dev = SQRT(SUM(ABS(pw%cc - mean)**2)/count)
         min_val = MINVAL(ABS(pw%cc))
         max_val = MAXVAL(ABS(pw%cc))
         norm_2 = SQRT(SUM(ABS(pw%cc)**2))
      ELSE
         PRINT *, "KKK| Error: Unknown space type for PW"
         RETURN
      END IF

      ! PW Type Information
      PRINT *, "KKK| Statistics for -> ", name, " <-(fragment ", i_frag, "):"
      PRINT *, "KKK|  PW Type Information:"
      PRINT *, "KKK|   In use:                ", pw%in_use
      PRINT *, "KKK|   Space:                 ", MERGE("REALSPACE      ", &
                                                       "RECIPROCALSPACE", &
                                                       pw%in_space == REALSPACE)
      PRINT *, "KKK|  Data Statistics:"
      PRINT *, "KKK|   Mean(abs):             ", mean
      PRINT *, "KKK|   Std Dev:               ", std_dev
      PRINT *, "KKK|   Min(abs):              ", min_val
      PRINT *, "KKK|   Max(abs):              ", max_val
      PRINT *, "KKK|   Sum(abs):              ", sum_val
      PRINT *, "KKK|   Norm-2:                ", norm_2
      PRINT *, "KKK|   Num of grid p.         ", count

      IF (print_grid) THEN
         PRINT *, "KKK|  Grid Information:"
         PRINT *, "KKK|   Total Grid Points:     ", grid%ngpts
         PRINT *, "KKK|   Grid Points in Cutoff: ", grid%ngpts_cut
         PRINT *, "KKK|   Dimensions:            ", grid%npts
         PRINT *, "KKK|   Local Grid Points:     ", grid%ngpts_local
         PRINT *, "KKK|   Local Points in Cutoff:", grid%ngpts_cut_local
         PRINT *, "KKK|   Grid Bounds:           ", grid%bounds
         PRINT *, "KKK|   Local Grid Bounds:     ", grid%bounds_local
         PRINT *, "KKK|   Grid Spacing:          ", grid%dr
         PRINT *, "KKK|   Cell Matrix:"
         DO i = 1, 3
            PRINT *, "KKK|                          ", grid%dh(:, i)
         END DO
         PRINT *, "KKK|   Inverse Cell Matrix:"
         DO i = 1, 3
            PRINT *, "KKK|                          ", grid%dh_inv(:, i)
         END DO
         PRINT *, "KKK|   Orthorhombic:          ", grid%orthorhombic
         PRINT *, "KKK|   Volume Element:        ", grid%dvol
         PRINT *, "KKK|   Total Volume:          ", grid%vol
         PRINT *, "KKK|   Cutoff:                ", grid%cutoff
         PRINT *, "KKK|   Grid Span:             ", MERGE("HALFSPACE", &
                                                          "FULLSPACE", &
                                                          grid%grid_span == HALFSPACE)
         PRINT *, "KKK|   Has G=0:               ", grid%have_g0
         PRINT *, "KKK|   First G≠0 index:       ", grid%first_gne0
         PRINT *, "KKK|   Grid ID:               ", grid%id_nr
         PRINT *, "KKK|   Reference Count:       ", grid%ref_count
         PRINT *, "KKK|   Spherical Cutoff:      ", grid%spherical
         PRINT *, "KKK|  Parallel Information:"
         PRINT *, "KKK|   Mode:                  ", MERGE("Distributed", &
                                                          "Local      ", &
                                                          grid%para%mode == PW_MODE_DISTRIBUTED)
         PRINT *, "KKK|   Ray Distribution:      ", grid%para%ray_distribution
         PRINT *, "KKK|   Blocked:               ", grid%para%blocked
         PRINT *, "KKK|   Group Size:            ", grid%para%group_size
         PRINT *, "KKK|   My Position:           ", grid%para%my_pos
         PRINT *, "KKK|   Group Head:            ", grid%para%group_head
         PRINT *, "KKK|   Group Head ID:         ", grid%para%group_head_id
         PRINT *, "KKK|   Real Space Dimensions: ", grid%para%rs_dims
         PRINT *, "KKK|   Real Space Position:   ", grid%para%rs_pos
      END IF

      IF (ASSOCIATED(grid%g)) THEN
         PRINT *, "KKK|  G-vector Statistics:"
         PRINT *, "KKK|   Mean G magnitude:      ", SUM(SQRT(SUM(grid%g**2, dim=1)))/SIZE(grid%g, 2)
         PRINT *, "KKK|   Max G magnitude:       ", SQRT(MAXVAL(SUM(grid%g**2, dim=1)))
      END IF

      IF (ASSOCIATED(grid%gsq)) THEN
         PRINT *, "KKK|   Mean G^2:              ", SUM(grid%gsq)/SIZE(grid%gsq)
         PRINT *, "KKK|   Max G^2:               ", MAXVAL(grid%gsq)
      END IF

!#if defined(__parallel)
!      PRINT *, "KKK|  MPI Information:"
!      PRINT *, "KKK|   MPI Mode:               ", MERGE("DISTRIBUTED", &
!                                                        "LOCAL      ", &
!                                                         grid%para%mode == PW_MODE_DISTRIBUTED)
!      PRINT *, "KKK|   MPI Group Handle:       ", grid%para%group%get_handle()
!      PRINT *, "KKK|   MPI Group Size:         ", grid%para%group_size
!      PRINT *, "KKK|   MPI Position:           ", grid%para%my_pos
!      PRINT *, "KKK|   MPI Dims:               ", grid%para%rs_dims
!      PRINT *, "KKK|   MPI Ray Distribution:   ", grid%para%ray_distribution
!      PRINT *, "KKK|   MPI Blocked:            ", grid%para%blocked
!      PRINT *, "KKK|   MPI Task Rank:          ", grid%para%rs_dims
!      PRINT *, "KKK|   Real Space Group:       ", grid%para%rs_group%get_handle()
!      PRINT *, "KKK|   Real Space Position:    ", grid%para%rs_pos
!      PRINT *, "KKK|   Real Space MPO:         ", grid%para%rs_mpo
!      PRINT *, "KKK|   Has Bounds Array:       ", ALLOCATED(grid%para%bo)
!      IF (ALLOCATED(grid%para%bo)) THEN
!         PRINT *, "KKK|   Bounds Array Shape:     ", SHAPE(grid%para%bo)
!      END IF
!      PRINT *, "KKK|   Has YZ Ray Array:       ", ALLOCATED(grid%para%yzp)
!      PRINT *, "KKK|   Has YZ Count Array:     ", ALLOCATED(grid%para%nyzray)
!#endif
   END SUBROUTINE print_pw_statistics

! **************************************************************************************************
!> \brief Gets the nuclear density for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_get_rhon_frag(fat_env, rhon_g, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: rhon_g
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_get_rhon_frag'

      INTEGER                                            :: handle
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      NULLIFY (rhon_g)

      qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
      CALL get_qs_env(qs_env=qs_env, rho_core=rhon_g)

      CALL timestop(handle)

   END SUBROUTINE fat_get_rhon_frag

! **************************************************************************************************
!> \brief Gets the nuclear density for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_get_ppl_frag(fat_env, v_ppl, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: v_ppl
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_get_ppl_frag'

      INTEGER                                            :: handle
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
      CALL calculate_ppl_grid(v_ppl, qs_env, do_erf_term=.TRUE.)
!     IF (is_debug) CALL print_pw_statistics(v_ppl, "v_ppl", i_frag)

      CALL timestop(handle)

   END SUBROUTINE fat_get_ppl_frag

! **************************************************************************************************
!> \brief Gets the gradient of the nuclear density w.r.t. a nuclear position for one fragment in the
!>        freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_get_drhon_frag(fat_env, drhon_g, i_frag, i_nuc, i_dir)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: drhon_g
      INTEGER                                            :: i_frag, i_nuc, i_dir

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_get_drhon_frag'

      INTEGER                                            :: handle
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
      CALL calculate_drho_core(drhon_g, qs_env, i_dir, i_nuc)

      CALL timestop(handle)

   END SUBROUTINE fat_get_drhon_frag

! **************************************************************************************************
!> \brief Gets the gradient of the electronic density w.r.t. a nuclear position for one fragment in
!>        the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_get_drho_frag(fat_env, drho_r, drho_g, i_frag, i_nuc, i_dir)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: drho_r
      TYPE(pw_type), POINTER                             :: drho_g
      INTEGER                                            :: i_frag, i_nuc, i_dir

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_get_drho_frag'

      INTEGER                                            :: handle
      TYPE(dbcsr_p_type), POINTER                        :: p_mat(:)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: rho_r

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CPASSERT(ASSOCIATED(drho_r) .AND. ASSOCIATED(drho_g))

      qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
      CALL get_qs_env(qs_env, rho=rho_r)
      CALL qs_rho_get(rho_r, rho_ao=p_mat)

      CALL calculate_drho_elec_dR(matrix_p=p_mat(1)%matrix, &
                                  drho=drho_r, &
                                  drho_gspace=drho_g, &
                                  qs_env=qs_env, &
                                  beta=i_dir, &
                                  lambda=i_nuc)

      CALL timestop(handle)

   END SUBROUTINE fat_get_drho_frag

! **************************************************************************************************
!> \brief Gets the total electron density for the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_get_rho_ref(fat_env, rho_r, rho_g)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), OPTIONAL                            :: rho_r(:)
      TYPE(pw_type), OPTIONAL, POINTER                   :: rho_g(:)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_get_rho_ref'

      INTEGER                                            :: handle, i, j
      TYPE(pw_type), POINTER                             :: rho_r_tmp(:)

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      DO i = 1, fat_env%n_frag
         CALL fat_get_rho_frag(fat_env, rho_r=rho_r_tmp, i_frag=i)
         DO j = 1, SIZE(rho_r_tmp)
            CALL pw_axpy(rho_r_tmp(j), rho_r(j), allow_noncompatible_grids=.TRUE.)
         END DO
      END DO

      IF (PRESENT(rho_g)) THEN
         CPASSERT(ASSOCIATED(rho_g))
         DO j = 1, SIZE(rho_r)
            CALL pw_transfer(rho_r(j), rho_g(j))
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_get_rho_ref

! **************************************************************************************************
!> \brief Calculate kinetic potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_v_ts_frag(fat_env, rho_r, rho_g, v_kin_pot_frag, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: rho_r(:)
      TYPE(pw_type), POINTER                             :: rho_g(:)
      TYPE(pw_type), POINTER                             :: v_kin_pot_frag(:)
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_v_ts_frag'

      INTEGER                                            :: handle
      REAL(dp)                                           :: e_kin_frag, virial(3, 3)
      TYPE(pw_type), POINTER                             :: tau(:), v_tau(:)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input, kin_section

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      ! Get XC section for kinetic potential calculation
      qs_env => fat_env%force_env_ref%qs_env
      CALL get_qs_env(qs_env, input=input)
      kin_section => section_vals_get_subs_vals(input, "DFT%QS%OPT_EMBED%XC")

      ! Calculate fragment kinetic potential and energy
      NULLIFY (v_kin_pot_frag, v_tau)
      CALL xc_vxc_pw_create(vxc_rho=v_kin_pot_frag, &
                            vxc_tau=v_tau, &
                            rho_r=rho_r, &
                            rho_g=rho_g, &
                            tau=tau, &
                            exc=e_kin_frag, &
                            xc_section=kin_section, &
                            pw_pool=fat_env%pw_pool_sub(i_frag)%pool, &
                            compute_virial=.FALSE., &
                            virial_xc=virial)

      ! Store the fragment kinetic energy contribution
      fat_env%ekin_sub(i_frag) = e_kin_frag

      CALL timestop(handle)

   END SUBROUTINE fat_v_ts_frag

! **************************************************************************************************
!> \brief Calculate reference kinetic potential for the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_v_ts_ref(fat_env, v_kin_pot_ref)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: v_kin_pot_ref(:)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_v_ts_ref'

      INTEGER                                            :: handle
      REAL(dp)                                           :: e_kin_ref, virial(3, 3)
      TYPE(pw_type), POINTER                             :: rho_g(:), rho_r(:), tau(:), v_tau(:)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input, kin_section

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      ! Get XC section for kinetic potential calculation
      qs_env => fat_env%force_env_ref%qs_env
      CALL get_qs_env(qs_env, input=input)
      kin_section => section_vals_get_subs_vals(input, "DFT%QS%OPT_EMBED%XC")

      ! Create pws for total density
      CALL fat_pws_create_ref(fat_env, pws_r=rho_r, pws_g=rho_g, n_pws=1)

      ! Get density (rho_r is for LDA, rho_g is for GGA)
      CALL fat_get_rho_ref(fat_env, rho_r=rho_r, rho_g=rho_g)

      ! Calculate total system kinetic potential and energy
      NULLIFY (v_kin_pot_ref, v_tau)
      CALL xc_vxc_pw_create(vxc_rho=v_kin_pot_ref, &
                            vxc_tau=v_tau, &
                            rho_r=rho_r, &
                            rho_g=rho_g, &
                            tau=tau, &
                            exc=e_kin_ref, &
                            xc_section=kin_section, &
                            pw_pool=fat_env%pw_pool_ref, &
                            compute_virial=.FALSE., &
                            virial_xc=virial)

      ! Store the total kinetic energy contribution
      fat_env%ekin_ref = e_kin_ref

      ! Cleanup
      CALL fat_pws_release_ref(fat_env, pws_r=rho_r, pws_g=rho_g)

      CALL timestop(handle)

   END SUBROUTINE fat_v_ts_ref

! **************************************************************************************************
!> \brief Calculate xc potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_v_xc_frag(fat_env, rho_r, rho_g, v_xc_pot_frag, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: rho_r(:)
      TYPE(pw_type), POINTER                             :: rho_g(:)
      TYPE(pw_type), POINTER                             :: v_xc_pot_frag(:)
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_v_xc_frag'

      INTEGER                                            :: handle
      REAL(dp)                                           :: e_xc_frag, virial(3, 3)
      TYPE(pw_type), POINTER                             :: tau(:), v_tau(:)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input, xc_section

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      ! Get XC section for XC potential calculation
      qs_env => fat_env%force_env_ref%qs_env
      CALL get_qs_env(qs_env, input=input)
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")

      ! Calculate fragment XC potential and energy
      NULLIFY (v_xc_pot_frag, v_tau)
      CALL xc_vxc_pw_create(vxc_rho=v_xc_pot_frag, &
                            vxc_tau=v_tau, &
                            rho_r=rho_r, &
                            rho_g=rho_g, &
                            tau=tau, &
                            exc=e_xc_frag, &
                            xc_section=xc_section, &
                            pw_pool=fat_env%pw_pool_sub(i_frag)%pool, &
                            compute_virial=.FALSE., &
                            virial_xc=virial)

      ! Store the fragment XC energy contribution
      fat_env%exc_sub(i_frag) = e_xc_frag

      CALL timestop(handle)

   END SUBROUTINE fat_v_xc_frag

! **************************************************************************************************
!> \brief Calculate reference xc potential for the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_v_xc_ref(fat_env, v_xc_pot_ref)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: v_xc_pot_ref(:)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_v_xc_ref'

      INTEGER                                            :: handle
      REAL(dp)                                           :: e_xc_ref, virial(3, 3)
      TYPE(pw_type), POINTER                             :: rho_g(:), rho_r(:), tau(:), v_tau(:)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input, xc_section

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      ! Get XC section for XC potential calculation
      qs_env => fat_env%force_env_ref%qs_env
      CALL get_qs_env(qs_env, input=input)
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")

      ! Create pws for total density
      CALL fat_pws_create_ref(fat_env, pws_r=rho_r, pws_g=rho_g, n_pws=1)

      ! Get density (rho_r is for LDA, rho_g is for GGA)
      CALL fat_get_rho_ref(fat_env, rho_r=rho_r, rho_g=rho_g)

      ! Calculate total system XC potential and energy
      NULLIFY (v_xc_pot_ref, v_tau)
      CALL xc_vxc_pw_create(vxc_rho=v_xc_pot_ref, &
                            vxc_tau=v_tau, &
                            rho_r=rho_r, &
                            rho_g=rho_g, &
                            tau=tau, &
                            exc=e_xc_ref, &
                            xc_section=xc_section, &
                            pw_pool=fat_env%pw_pool_ref, &
                            compute_virial=.FALSE., &
                            virial_xc=virial)

      ! Store the total XC energy contribution
      fat_env%exc_ref = e_xc_ref

      ! Cleanup
      CALL fat_pws_release_ref(fat_env, pws_r=rho_r, pws_g=rho_g)

      CALL timestop(handle)

   END SUBROUTINE fat_v_xc_ref

! **************************************************************************************************
!> \brief Compute the Hartree potential for a fragment using its total density (electronic + core)
!> \param fat_env The freeze-and-thaw environment
!> \param v_hartree The Hartree potential (output)
!> \param i_frag The fragment index
! **************************************************************************************************
   SUBROUTINE fat_v_hartree_frag(fat_env, v_hartree, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: v_hartree
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_v_hartree_frag'

      INTEGER                                            :: handle, ispin
      LOGICAL                                            :: gapw, gapw_xc
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_type), POINTER                             :: rho0_s_gs, rho_core, rho_g(:), &
                                                            rho_tot_g, v_hartree_g
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      ! Get QS environment and DFT control
      qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
      CALL get_qs_env(qs_env, pw_env=pw_env, dft_control=dft_control)
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc

      ! Get Poisson environment
      CALL pw_env_get(pw_env, poisson_env=poisson_env)

      ! Get electronic density
      CALL fat_get_rho_frag(fat_env, rho_g=rho_g, i_frag=i_frag)

      ! Create temporary plane waves
      CALL fat_pw_create_frag(fat_env, pw_g=rho_tot_g, i_frag=i_frag)
      CALL fat_pw_create_frag(fat_env, pw_g=v_hartree_g, i_frag=i_frag)

      ! Compute total density
      IF (gapw .OR. gapw_xc) THEN
         CALL get_qs_env(qs_env, rho0_s_gs=rho0_s_gs)
         DO ispin = 1, SIZE(rho_g)
            CALL pw_axpy(rho_g(ispin), rho_tot_g)
         END DO
         CALL pw_axpy(rho0_s_gs, rho_tot_g)
      ELSE
         CALL fat_get_rhon_frag(fat_env, rhon_g=rho_core, i_frag=i_frag)
         DO ispin = 1, SIZE(rho_g)
            CALL pw_axpy(rho_g(ispin), rho_tot_g)
         END DO
         CALL pw_axpy(rho_core, rho_tot_g)
      END IF

      ! Solve Poisson equation to get Hartree potential
      CALL pw_poisson_solve(poisson_env, density=rho_tot_g, vhartree=v_hartree_g)
      CALL pw_transfer(v_hartree_g, v_hartree)
      CALL pw_scale(v_hartree, v_hartree%pw_grid%dvol)

      ! Clean up
      CALL fat_pw_release_frag(fat_env, pw_g=rho_tot_g, i_frag=i_frag)
      CALL fat_pw_release_frag(fat_env, pw_g=v_hartree_g, i_frag=i_frag)

      CALL timestop(handle)
   END SUBROUTINE fat_v_hartree_frag

! **************************************************************************************************
!> \brief Calculate non-additive kinetic embedding potential for one fragment in the freeze-and-thaw
!>        environment
! **************************************************************************************************
   SUBROUTINE fat_vemb_tsna_frag(fat_env, rho_r, rho_g, v_emb, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: rho_r(:)
      TYPE(pw_type), POINTER                             :: rho_g(:)
      TYPE(pw_type), POINTER                             :: v_emb
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_vemb_tsna_frag'

      INTEGER                                            :: handle
      TYPE(pw_type), POINTER                             :: v(:)

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CALL pw_zero(v_emb)

      ! Reference system
      NULLIFY (v)
      CALL fat_v_ts_ref(fat_env, v_kin_pot_ref=v)
      CALL pw_axpy(pw1=v(1), pw2=v_emb, allow_noncompatible_grids=.TRUE.)
      CALL fat_pws_release_ref(fat_env, pws_r=v)

      ! Fragment
      NULLIFY (v)
      CALL fat_v_ts_frag(fat_env, rho_r=rho_r, rho_g=rho_g, &
                         v_kin_pot_frag=v, i_frag=i_frag)
      CALL pw_axpy(pw1=v(1), pw2=v_emb, alpha=-1.0_dp, allow_noncompatible_grids=.TRUE.)
      CALL fat_pws_release_frag(fat_env, pws_r=v, i_frag=i_frag)

      CALL timestop(handle)

   END SUBROUTINE fat_vemb_tsna_frag

! **************************************************************************************************
!> \brief Calculate non-additive xc embedding potential for one fragment in the freeze-and-thaw
!>        environment
! **************************************************************************************************
   SUBROUTINE fat_vemb_xcna_frag(fat_env, rho_r, rho_g, v_emb, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: rho_r(:)
      TYPE(pw_type), POINTER                             :: rho_g(:)
      TYPE(pw_type), POINTER                             :: v_emb
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_vemb_xcna_frag'

      INTEGER                                            :: handle
      TYPE(pw_type), POINTER                             :: v(:)

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CALL pw_zero(v_emb)

      ! Reference system
      NULLIFY (v)
      CALL fat_v_xc_ref(fat_env, v_xc_pot_ref=v)
      CALL pw_axpy(pw1=v(1), pw2=v_emb, allow_noncompatible_grids=.TRUE.)
      CALL fat_pws_release_ref(fat_env, pws_r=v)

      ! Fragment
      NULLIFY (v)
      CALL fat_v_xc_frag(fat_env, rho_r=rho_r, rho_g=rho_g, &
                         v_xc_pot_frag=v, i_frag=i_frag)
      CALL pw_axpy(pw1=v(1), pw2=v_emb, alpha=-1.0_dp, allow_noncompatible_grids=.TRUE.)
      CALL fat_pws_release_frag(fat_env, pws_r=v, i_frag=i_frag)

      CALL timestop(handle)

   END SUBROUTINE fat_vemb_xcna_frag

! **************************************************************************************************
!> \brief Calculate Hartree potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_vemb_hartree_frag(fat_env, v_emb, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: v_emb
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_vemb_hartree_frag'

      INTEGER                                            :: handle, j_frag
      TYPE(pw_type), POINTER                             :: v_hartree_j

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CALL pw_zero(v_emb)

      ! Sum Hartree potentials from all other fragments
      DO j_frag = 1, fat_env%n_frag

         IF (j_frag == i_frag) CYCLE

         CALL fat_pw_create_frag(fat_env, pw_r=v_hartree_j, i_frag=j_frag)

         CALL fat_v_hartree_frag(fat_env, v_hartree_j, j_frag)

         CALL pw_axpy(v_hartree_j, v_emb, allow_noncompatible_grids=.TRUE.)

         CALL fat_pw_release_frag(fat_env, pw_r=v_hartree_j, i_frag=j_frag)
      END DO

      CALL timestop(handle)

   END SUBROUTINE fat_vemb_hartree_frag

! **************************************************************************************************
!> \brief Write embedding potential to cube and dat files
!> \param fat_env The freeze-and-thaw environment
!> \param v_emb The embedding potential to write
!> \param work_tmp Temporary workspace for potential calculations
!> \param i_frag The fragment index
! **************************************************************************************************
   SUBROUTINE fat_write_embedding_potential(fat_env, v_emb, work_tmp, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: v_emb, work_tmp
      INTEGER, INTENT(IN)                                :: i_frag

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_write_embedding_potential'

      CHARACTER(len=default_string_length)               :: label
      INTEGER                                            :: handle, output_unit
      INTEGER, POINTER                                   :: stride(:)
      TYPE(cp_logger_type), POINTER                      :: logger
      CHARACTER(len=50)                                  :: filename

      CALL timeset(routineN, handle)

      IF (is_debug .and. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", moduleN//"::"//routineN
      END IF

      logger => cp_get_default_logger()
      output_unit = -1

      ALLOCATE (stride(3), source=1)
      ! Open file only on the master process
      IF (fat_env%force_env_ref%para_env%is_source()) THEN
         CALL integer_to_string(i_frag, label)
         CALL open_file(file_name="vemb-cp2k_"//TRIM(label)//".dat", &
                        file_status="REPLACE", &
                        file_action="WRITE", &
                        unit_number=output_unit)
      END IF

      CALL pw_zero(work_tmp)

      ! Set weights (work_tmp) to the volume element of the grid
      CALL pw_set(work_tmp, &
                  v_emb%pw_grid%dr(1)* &
                  v_emb%pw_grid%dr(2)* &
                  v_emb%pw_grid%dr(3))

      ! Use v_emb as the primary grid, no need for weights in the file
      CALL cp_pw_to_simple_volumetric(pw=work_tmp, &
                                      unit_nr=output_unit, &
                                      stride=stride, &
                                      pw2=v_emb)

      ! Close file on master process
      IF (fat_env%force_env_ref%para_env%is_source()) THEN
         CALL close_file(unit_number=output_unit)
      END IF

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         WRITE(filename,'(A,I0,A,I0,A)') 'embedding_dat_frag', i_frag, '.cube'
         CALL write_pw_to_cube(fat_env, pw_r=v_emb, filename=filename)
      END IF

      DEALLOCATE(stride)

      ! Synchronize all processes
      CALL fat_env%force_env_ref%para_env%sync()

      CALL timestop(handle)

   END SUBROUTINE fat_write_embedding_potential

   ! **************************************************************************************************
   !> \brief Execute external calculation script
   ! **************************************************************************************************
   SUBROUTINE fat_run_external(fat_env, i_frag, gradient)
      TYPE(fat_env_type), POINTER                        :: fat_env
      INTEGER, INTENT(IN)                                :: i_frag
      LOGICAL, OPTIONAL                                  :: gradient

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_run_external'

      CHARACTER                                          :: cmdmsg(100)
      CHARACTER(default_string_length)                   :: label
      INTEGER                                            :: cmdstat, exitstat, handle
      LOGICAL                                            :: my_grad, file_exists

      CALL timeset(routineN, handle)

      my_grad = .FALSE.
      ! Proceed if it is the source node
      IF (fat_env%force_env_ref%para_env%is_source()) THEN

         if (is_debug) PRINT *, "DBG| ", modulen//"::"//routinen

         IF (PRESENT(gradient)) my_grad = gradient

         CALL integer_to_string(i_frag, label)
         IF (my_grad) THEN
            INQUIRE(FILE="extern-grad_"//TRIM(label)//".sh", EXIST=file_exists)
            IF (.NOT. file_exists) THEN
               PRINT *, "ERROR| File extern-grad_"//TRIM(label)//".sh not found for fragment ", i_frag
               CALL cp_abort(__LOCATION__, "External gradient script file missing")
            END IF
            CALL execute_command_line("sh extern-grad_"//TRIM(label)//".sh", .TRUE., exitstat, cmdstat, cmdmsg)
         ELSE
            INQUIRE(FILE="extern_"//TRIM(label)//".sh", EXIST=file_exists)
            IF (.NOT. file_exists) THEN
               PRINT *, "ERROR| File extern_"//TRIM(label)//".sh not found for fragment ", i_frag
               CALL cp_abort(__LOCATION__, "External script file missing")
            END IF
            CALL execute_command_line("sh extern_"//TRIM(label)//".sh", .TRUE., exitstat, cmdstat, cmdmsg)
         END IF

         IF (cmdstat /= 0) CPABORT("External command did not execute successfully")
      END IF

      CALL fat_env%force_env_ref%para_env%sync()

      CALL timestop(handle)

   END SUBROUTINE fat_run_external

! **************************************************************************************************
!> \brief Read energy from external calculation
!> \param fat_env The freeze-and-thaw environment
!> \param i_frag Fragment index
! **************************************************************************************************
   SUBROUTINE fat_read_external_energy(fat_env, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      INTEGER, INTENT(IN)                                :: i_frag

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_read_external_energy'

      CHARACTER(len=default_string_length)               :: label
      INTEGER                                            :: input_unit, handle
      LOGICAL                                            :: file_exists

      CALL timeset(routineN, handle)

      IF (fat_env%force_env_ref%para_env%is_source()) THEN

         IF (is_debug) PRINT *, "DBG| ", modulen//"::"//routinen

         CALL integer_to_string(i_frag, label)

         INQUIRE(FILE="extern_"//TRIM(label)//".e", EXIST=file_exists)
         IF (.NOT. file_exists) THEN
            PRINT *, "ERROR| File extern_"//TRIM(label)//".e not found for fragment ", i_frag
            CALL cp_abort(__LOCATION__, "External energy file missing")
         END IF

         CALL open_file("extern_"//TRIM(label)//".e", &
                        file_status="OLD", &
                        file_action="READ", &
                        unit_number=input_unit)
         READ (input_unit, *) fat_env%e_sub(i_frag)
         CALL close_file(unit_number=input_unit)
      END IF

      CALL fat_env%force_env_ref%para_env%bcast(fat_env%e_sub(i_frag), &
                                                fat_env%force_env_ref%para_env%source)

      CALL timestop(handle)

   END SUBROUTINE fat_read_external_energy

! **************************************************************************************************
!> \brief Read density from external calculation
!> \param fat_env The freeze-and-thaw environment
!> \param i_frag Fragment index
! **************************************************************************************************
   SUBROUTINE fat_get_external_density(fat_env, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_get_external_density'

      INTEGER                                            :: handle
      CHARACTER(len=default_string_length)               :: label
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: qs_rho
      TYPE(pw_type), DIMENSION(:), POINTER               :: rho_r, rho_g
      LOGICAL                                            :: file_exists

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CPASSERT(fat_env%extern(i_frag))

      qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
      CALL get_qs_env(qs_env=qs_env, rho=qs_rho)
      CALL qs_rho_get(qs_rho, rho_r=rho_r, rho_g=rho_g)
!     IF (is_debug) CALL print_pw_statistics(rho_r(1), "before", i_frag)

      CALL integer_to_string(i_frag, label)
      IF (fat_env%force_env_ref%para_env%is_source()) THEN
         INQUIRE(FILE="extern_"//TRIM(label)//".cube", EXIST=file_exists)
         IF (.NOT. file_exists) THEN
            PRINT *, "ERROR| File extern_"//TRIM(label)//".cube not found for fragment ", i_frag
            CALL cp_abort(__LOCATION__, "External density file missing")
         END IF
      END IF
      CALL cp_cube_to_pw(rho_r(1), "extern_"//TRIM(label)//".cube", 1.0_dp)
!     IF (is_debug) CALL print_pw_statistics(rho_r(1), "after", i_frag)

      CALL pw_transfer(rho_r(1), rho_g(1))

      CALL qs_rho_set(qs_rho, rho_r=rho_r, rho_g=rho_g)

      CALL timestop(handle)

   END SUBROUTINE fat_get_external_density

! **************************************************************************************************
!> \brief Calculate embedding potential for one fragment in the freeze-and-thaw environment
!> \param fat_env The freeze-and-thaw environment
!> \param v_emb The embedding potential (output)
!> \param i_frag The index of the current fragment
!> \param to_file Boolean flag to determine if the potential should be written to a file
! **************************************************************************************************
   SUBROUTINE fat_vemb_frag(fat_env, v_emb, i_frag, to_file)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(pw_type), POINTER                             :: v_emb
      INTEGER                                            :: i_frag
      LOGICAL                                            :: to_file

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_vemb_frag'

      INTEGER                                            :: handle
      TYPE(pw_type), POINTER                             :: rho_g(:), rho_r(:), vemb_tmp

      CALL timeset(routineN, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CALL fat_get_rho_frag(fat_env, rho_r=rho_r, rho_g=rho_g, i_frag=i_frag)

      CALL fat_pw_create_frag(fat_env, pw_r=vemb_tmp, i_frag=i_frag)

      ! Calculate and add the non-additive kinetic potential
      CALL fat_vemb_tsna_frag(fat_env, rho_r=rho_r, rho_g=rho_g, &
                              v_emb=vemb_tmp, i_frag=i_frag)
      CALL pw_axpy(vemb_tmp, v_emb)

      ! Calculate and add the non-additive exchange-correlation potential
      CALL fat_vemb_xcna_frag(fat_env, rho_r=rho_r, rho_g=rho_g, &
                              v_emb=vemb_tmp, i_frag=i_frag)
      CALL pw_axpy(vemb_tmp, v_emb)

      ! Calculate and add the Hartree potential
      CALL fat_vemb_hartree_frag(fat_env, v_emb=vemb_tmp, i_frag=i_frag)
      CALL pw_axpy(vemb_tmp, v_emb)

      ! If requested, write the embedding potential to files
      IF (to_file) THEN
         CALL fat_write_embedding_potential(fat_env, v_emb=v_emb, &
                                            work_tmp=vemb_tmp, i_frag=i_frag)
      END IF

      CALL fat_pw_release_frag(fat_env, pw_r=vemb_tmp, i_frag=i_frag)

      CALL timestop(handle)

   END SUBROUTINE fat_vemb_frag

! **************************************************************************************************
!> \brief Calculate the kinetic energy for on fragment of a freeze-and-thaw optimization
! **************************************************************************************************
   FUNCTION fat_e_ts_frag(fat_env, i_frag) RESULT(e)
      TYPE(fat_env_type), POINTER                        :: fat_env
      INTEGER                                            :: i_frag
      REAL(dp)                                           :: e

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_e_ts_frag'

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      e = fat_env%ekin_sub(i_frag)

   END FUNCTION fat_e_ts_frag

! **************************************************************************************************
!> \brief Calculate the kinetic energy for the reference system of a freeze-and-thaw optimization
! **************************************************************************************************
   FUNCTION fat_e_ts_ref(fat_env) RESULT(e)
      TYPE(fat_env_type), POINTER                        :: fat_env
      REAL(dp)                                           :: e

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_e_ts_ref'

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      e = fat_env%ekin_ref

   END FUNCTION fat_e_ts_ref

! **************************************************************************************************
!> \brief Calculate the xc energy for on fragment of a freeze-and-thaw optimization
! **************************************************************************************************
   FUNCTION fat_e_xc_frag(fat_env, i_frag) RESULT(e)
      TYPE(fat_env_type), POINTER                        :: fat_env
      INTEGER                                            :: i_frag
      REAL(dp)                                           :: e

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_e_xc_frag'

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      e = fat_env%exc_sub(i_frag)

   END FUNCTION fat_e_xc_frag

! **************************************************************************************************
!> \brief Calculate the xc energy for the reference system of a freeze-and-thaw optimization
! **************************************************************************************************
   FUNCTION fat_e_xc_ref(fat_env) RESULT(e)
      TYPE(fat_env_type), POINTER                        :: fat_env
      REAL(dp)                                           :: e

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_e_xc_ref'

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      e = fat_env%exc_ref

   END FUNCTION fat_e_xc_ref

! **************************************************************************************************
!> \brief Compute the Hartree interaction energy between two fragments
!> \param fat_env The freeze-and-thaw environment
!> \param i First fragment index
!> \param j Second fragment index
!> \return e The interaction energy ∫ ρ_tot,i V_H,j dr
! **************************************************************************************************
   FUNCTION fat_eint_hartree(fat_env, i_frag, j_frag) RESULT(eint)
      TYPE(fat_env_type), POINTER                        :: fat_env
      INTEGER, INTENT(in)                                :: i_frag, j_frag
      REAL(dp)                                           :: eint

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_eint_hartree'

      INTEGER                                            :: handle
      TYPE(pw_type), POINTER                             :: rho_tot_i, rho_tot_j, v_hartree_i

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CALL fat_pw_create_frag(fat_env, pw_r=v_hartree_i, i_frag=i_frag)
      CALL fat_pw_create_frag(fat_env, pw_r=rho_tot_i, i_frag=i_frag)
      CALL fat_pw_create_frag(fat_env, pw_r=rho_tot_j, i_frag=j_frag)

      CALL fat_v_hartree_frag(fat_env, v_hartree_i, i_frag)

      CALL fat_get_rho_tot_frag(fat_env, rho_tot_j, j_frag)

      CALL pw_transfer(rho_tot_j, rho_tot_i)

      eint = pw_integral_ab(v_hartree_i, rho_tot_i)

      CALL fat_pw_release_frag(fat_env, pw_r=v_hartree_i, i_frag=i_frag)
      CALL fat_pw_release_frag(fat_env, pw_r=rho_tot_j, i_frag=j_frag)
      CALL fat_pw_release_frag(fat_env, pw_r=rho_tot_i, i_frag=i_frag)

      CALL timestop(handle)

   END FUNCTION fat_eint_hartree

! **************************************************************************************************
!> \brief Calculate the total energy of a freeze-and-thaw optimization
! **************************************************************************************************
   FUNCTION fat_e_tot(fat_env) RESULT(e)
      TYPE(fat_env_type), POINTER                        :: fat_env
      REAL(dp)                                           :: e

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_e_tot'

      INTEGER                                            :: i, j
      REAL(dp)                                           :: e_emb_total, e_hartree_total, &
                                                            e_sub_total, e_tmp, &
                                                            e_ts_total, e_xc_total
      LOGICAL                                            :: is_dbg, is_fine

      is_fine = fat_env%force_env_ref%para_env%is_source()
      is_dbg = is_debug .AND. is_fine

      IF (is_dbg) PRINT *, "DBG| ", modulen//"::"//routinen

      e = 0.0_dp

      ! Sum of subsystem energies
      e_sub_total = 0.0_dp
      DO i = 1, fat_env%n_frag
         e_sub_total = e_sub_total + fat_env%e_sub(i)
         IF (is_fine) PRINT *, "FAT| Subsystem energy (", i, "):", fat_env%e_sub(i)
      END DO
      e = e + e_sub_total
      IF (is_fine) PRINT *, "FAT| Total subsystem energy:", e_sub_total

      ! Embedding potential energy
      e_emb_total = 0.0_dp
      DO i = 1, fat_env%n_frag
         e_emb_total = e_emb_total + fat_env%e_emb(i)
         IF (is_fine) PRINT *, "FAT| Embedding energy (", i, "):", fat_env%e_emb(i)
         IF (i > 1 .and. is_fine) PRINT *, &
            "FAT| Embedding energy diff (", i-1, "-", i, "):", &
            fat_env%e_emb(i) - fat_env%e_emb(i-1)
      END DO
      e = e - e_emb_total
      IF (is_fine) PRINT *, "FAT| Total Embedding energy:", e_emb_total

      ! Kinetic energy
      e_ts_total = 0.0_dp
      e_ts_total = fat_e_ts_ref(fat_env)
      IF (is_fine) PRINT *, "FAT| Reference Kin energy:", e_ts_total
      DO i = 1, fat_env%n_frag
         e_tmp = fat_e_ts_frag(fat_env, i)
         e_ts_total = e_ts_total - e_tmp
         IF (is_fine) PRINT *, "FAT| Subsystem Kin energy (", i, "):", e_tmp
      END DO
      e = e + e_ts_total
      IF (is_fine) PRINT *, "FAT| Total Kin energy:", e_ts_total

      ! Exchange-correlation energy
      e_xc_total = 0.0_dp
      e_xc_total = fat_e_xc_ref(fat_env)
      IF (is_fine) PRINT *, "FAT| Reference XC energy:", e_xc_total
      DO i = 1, fat_env%n_frag
         e_tmp = fat_e_xc_frag(fat_env, i)
         e_xc_total = e_xc_total - e_tmp
         IF (is_fine) PRINT *, "FAT| Subsystem XC energy (", i, "):", e_tmp
      END DO
      e = e + e_xc_total
      IF (is_fine) PRINT *, "FAT| Total XC energy:", e_xc_total

      IF (fat_env%i_iter > 0) THEN
         ! Hartree electrostatic energy
         e_hartree_total = 0.0_dp
         DO i = 1, fat_env%n_frag
            DO j = i+1, fat_env%n_frag
               e_tmp = fat_eint_hartree(fat_env, i, j)
               e_hartree_total = e_hartree_total + e_tmp
               IF (is_fine) PRINT *, "FAT| Hartree interaction energy (", i, ",", j, "):", e_tmp
            END DO
         END DO
         e = e + e_hartree_total
         fat_env%ehartree = e_hartree_total
         IF (is_fine) PRINT *, "FAT| Total Hartree interaction energy:", e_hartree_total
      END IF

      fat_env%e_tot = e
      IF (is_fine) THEN
         PRINT *, "FAT| Total energy:", e
         PRINT *, "FAT| FAT - Frag total energies, kcal/mol:", (fat_env%e_tot-e_sub_total)*627.5094740631121212
         PRINT *, "FAT| Total Embedding energy kcal/mol:", e_emb_total*627.5094740631121212
      END IF

   END FUNCTION fat_e_tot

! ************************************************************************************************
!> \brief Print the intermediate results of a freeze-and-thaw iteration
! **************************************************************************************************
   SUBROUTINE fat_status(fat_env, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_status'

      LOGICAL                                            :: ionode
      REAL(dp)                                           :: e_tot

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      ionode = (fat_env%force_env_ref%para_env%mepos == fat_env%force_env_ref%para_env%source)

      e_tot = fat_e_tot(fat_env)

      IF (ionode .AND. fat_env%output_unit > 0) &
         WRITE (fat_env%output_unit, fmt=fat_data_fmt) &
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
         fat_env%ehartree, &
         e_tot

   END SUBROUTINE fat_status

! **************************************************************************************************
!> \brief Check the freeze-and-thaw optimization for convergence
!> \param fat_env The freeze-and-thaw environment containing all calculation data
!> \param i_frag The index of the current fragment being processed
! **************************************************************************************************
   SUBROUTINE fat_iterate_frag(fat_env, i_frag)
      TYPE(fat_env_type), POINTER                        :: fat_env
      INTEGER                                            :: i_frag

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_iterate_frag'

      INTEGER                                            :: handle
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(force_env_type), POINTER                      :: force_env_tmp
      TYPE(pw_type), POINTER                             :: rho(:), v_emb, weights
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_environment_type), POINTER                 :: qs_env_tmp
      TYPE(pw_type), POINTER                             :: rho_tmp
      CHARACTER(len=50)                                  :: filename

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen, " for frag:", i_frag
         PRINT *, "DBG| ", modulen//"::"//routinen, " fat_env%i_iter:", fat_env%i_iter
      END IF

!     ! MPI
!     CALL fat_env%force_env_ref%para_env%sync()

      ! Get the force environment and system environment for current fragment
      force_env_tmp => fat_env%force_env_sub(i_frag)%force_env
      qs_env_tmp => force_env_tmp%qs_env

      ! Get the DFT control settings and enable embedding potential
      CALL get_qs_env(qs_env_tmp, dft_control=dft_control)
      dft_control%apply_embed_pot = .TRUE.

!   ! Create a new embedding potential for current fragment
      CALL fat_pw_create_frag(fat_env, pw_r=v_emb, i_frag=i_frag)

      ! After first iteration, calculate the embedding potential
      IF (fat_env%i_iter > 0) THEN
         ! Calculate total embedding potential including:
         ! - Non-additive Kinetic potential
         ! - Non-additive XC potential
         ! - Electrostatic potential from other fragments
         CALL fat_vemb_frag(fat_env, v_emb=v_emb, i_frag=i_frag, to_file=.FALSE.)
         IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
            ! Write embedding to cube file
            WRITE(filename,'(A,I0,A,I0,A)') 'embedding_frag', i_frag, '.cube'
            CALL write_pw_to_cube(fat_env, pw_r=v_emb, filename=filename)
         END IF
      END IF

      ! Set the embedding potential to the qs env
      CALL set_qs_env(qs_env=qs_env_tmp, embed_pot=v_emb)

      ! Handle external fragment calculations if configured
      IF (fat_env%extern(i_frag)) THEN

         IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
            PRINT *, "DBG| Using extern for fragment ", i_frag
         END IF

         CALL fat_pw_create_frag(fat_env, pw_r=weights, i_frag=i_frag)
         CALL fat_write_embedding_potential(fat_env, v_emb=v_emb, &
                                            work_tmp=weights, i_frag=i_frag)
         CALL fat_pw_release_frag(fat_env, pw_r=weights, i_frag=i_frag)

         ! Execute external calculation script
         CALL fat_run_external(fat_env, i_frag=i_frag, gradient=.FALSE.)

         ! Get energy from external calculation
         CALL fat_read_external_energy(fat_env, i_frag)

         ! Get (read and set) density from external calculation
         CALL fat_get_external_density(fat_env, i_frag)
         ! Get electron density from external calculation
         ! We read the density once after the external calculation
         CALL fat_get_rho_frag(fat_env, rho_r=rho, i_frag=i_frag)

         ! Calculate embedding energy for current fragment:
         ! Integrate the density * embedding potential
         fat_env%e_emb(i_frag) = pw_integral_ab(v_emb, rho(1))

         IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
            ! Create temporary pw for density
            CALL fat_pw_create_frag(fat_env, pw_r=rho_tmp, i_frag=i_frag)
            CALL pw_transfer(rho(1), rho_tmp)
            ! Write density to cube file
            WRITE(filename,'(A,I0,A,I0,A)') 'rho_frag', i_frag, '_iter', fat_env%i_iter, '.cube'
            CALL write_pw_to_cube(fat_env, pw_r=rho_tmp, filename=filename)
            ! Deallocate
            CALL fat_pw_release_frag(fat_env, pw_r=rho_tmp, i_frag=i_frag)
         END IF
 
         IF (is_debug .and. fat_env%force_env_ref%para_env%is_source()) &
            PRINT *, "DBG| Embedding energy from fragment ", i_frag, fat_env%e_emb(i_frag)

      ELSE

         IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) &
            PRINT *, "DBG| Using internal calculation for fragment ", i_frag

         ! Perform internal DFT calculation for current fragment:
         CALL qs_calc_energy_force(qs_env=qs_env_tmp, &
                                   calc_force=.TRUE., &
                                   consistent_energies=.FALSE., &
                                   linres=.FALSE.)

         ! Get energy from DFT calculation
         CALL get_qs_env(qs_env_tmp, energy=energy)
         fat_env%e_sub(i_frag) = energy%total      ! Total energy

         ! Get embedding energy for current fragment
         fat_env%e_emb(i_frag) = energy%embed_corr

         IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) &
            PRINT *, "DBG| Embedding energy from fragment ", i_frag, fat_env%e_emb(i_frag)

         IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
            ! Get electron density from internal calculation
            CALL fat_get_rho_frag(fat_env, rho_r=rho, i_frag=i_frag)
            ! Create temporary pw for density
            CALL fat_pw_create_frag(fat_env, pw_r=rho_tmp, i_frag=i_frag)
            CALL pw_transfer(rho(1), rho_tmp)
            ! Write density to cube file
            WRITE(filename,'(A,I0,A,I0,A)') 'rho_frag', i_frag, '_iter', fat_env%i_iter, '.cube'
            CALL write_pw_to_cube(fat_env, pw_r=rho_tmp, filename=filename)
            ! Deallocate
            CALL fat_pw_release_frag(fat_env, pw_r=rho_tmp, i_frag=i_frag)
         END IF

      END IF

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) &
         PRINT *, "DBG| Total subsystem (", i_frag, ") energy:", fat_env%e_sub(i_frag)

      ! Deallocate v_emb potential
      CALL fat_pw_release_frag(fat_env, pw_r=v_emb, i_frag=i_frag)

      CALL timestop(handle)

   END SUBROUTINE fat_iterate_frag

! **************************************************************************************************
!> \brief Performs one iteration of the freeze-and-thaw optimization
! **************************************************************************************************
   SUBROUTINE fat_iterate(fat_env)
      TYPE(fat_env_type), POINTER                        :: fat_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_iterate'

      INTEGER                                            :: handle, i_frag

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      ! Loop over all fragments to optimize their densities with embedding effects
      DO i_frag = 1, fat_env%n_frag

         ! Calculate new density and energy for fragment i_frag under embedding potential from others
         CALL fat_iterate_frag(fat_env, i_frag)

         ! Print current optimization status
         CALL fat_status(fat_env, i_frag)
      END DO

      ! Increment iteration counter
      fat_env%i_iter = fat_env%i_iter + 1

      CALL timestop(handle)

   END SUBROUTINE fat_iterate

! **************************************************************************************************
!> \brief Check the freeze-and-thaw optimization for convergence
! **************************************************************************************************
   FUNCTION fat_is_converged(fat_env) RESULT(is_converged)
      TYPE(fat_env_type), POINTER                        :: fat_env
      LOGICAL                                            :: is_converged

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_is_converged'

      REAL(dp)                                           :: diff

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      diff = fat_env%e_tot - fat_env%e_old
      is_converged = ABS(diff) <= fat_env%eps
      fat_env%is_converged = is_converged

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", is_converged
         PRINT *, "DBG| ", fat_env%e_tot, fat_env%e_old
         PRINT *, "DBG| diff", diff
      END IF

      fat_env%e_old = fat_env%e_tot

   END FUNCTION

! **************************************************************************************************
!> \brief Check if the freeze-and-thaw optimization has converged or run into the iteration limit
! **************************************************************************************************
   FUNCTION fat_is_done(fat_env) RESULT(is_done)
      TYPE(fat_env_type), POINTER                        :: fat_env
      LOGICAL                                            :: is_done

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_is_done'

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      is_done = fat_env%i_iter > fat_env%n_iter .OR. fat_is_converged(fat_env)

   END FUNCTION

! **************************************************************************************************
!> \brief Allocate the forces for a freeze-and-thaw optimization
! **************************************************************************************************
   SUBROUTINE fat_force_create(fat_env, forces)
      TYPE(fat_env_type), POINTER                        :: fat_env
      REAL(dp), ALLOCATABLE, INTENT(out)                 :: forces(:, :)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_force_create'

      INTEGER                                            :: handle
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      qs_env => fat_env%force_env_ref%qs_env
      CALL get_qs_env(qs_env, particle_set=particles)
      ALLOCATE (forces(SIZE(particles), 3))
      forces = 0.0_dp

      CALL timestop(handle)

   END SUBROUTINE fat_force_create

! **************************************************************************************************
!> \brief Reads external gradients from a file and stores them in the qs_force structure
!> \param fat_env The FAT environment containing fragment and parallel information
!> \param i_frag Index of the fragment for which gradients are read
!> \param qs_force Array of force structures to store the gradients as forces
! **************************************************************************************************
   SUBROUTINE fat_get_gradients(fat_env, i_frag, qs_force)

      TYPE(fat_env_type), POINTER                        :: fat_env
      INTEGER, INTENT(IN)                                :: i_frag
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: qs_force

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_get_gradients'

      INTEGER                                            :: handle, i, i_atom, i_kind, input_unit, ios, n_atom 
      INTEGER, ALLOCATABLE                               :: atom_of_kind(:), kind_of(:)
      TYPE(atomic_kind_type), POINTER                    :: kinds(:)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(len=default_string_length)               :: label
      REAL(dp)                                           :: fx, fy, fz
      LOGICAL                                            :: file_exists

      CALL timeset(routineN, handle)

      ! Perform file reading only on the master process
      IF (fat_env%force_env_ref%para_env%is_source()) THEN

         CALL integer_to_string(i_frag, label)

         INQUIRE(FILE="extern_"//TRIM(label)//".gradients", EXIST=file_exists)
         IF (.NOT. file_exists) THEN
            PRINT *, "ERROR| File extern_"//TRIM(label)//".gradients not found for fragment ", i_frag
            CALL cp_abort(__LOCATION__, "External gradients file missing")
         END IF

         ! Retrieve fragment information from qs_env
         qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
         CALL get_qs_env(qs_env, atomic_kind_set=kinds, force=qs_force)
         CALL get_atomic_kind_set(kinds, natom=n_atom, kind_of=kind_of, atom_of_kind=atom_of_kind)

         ! Open the external gradients file
         CALL open_file(file_name="extern_"//TRIM(label)//".gradients", &
                        file_status="OLD", &
                        file_action="READ", &
                        unit_number=input_unit)

         ! Read gradients and store them as forces in qs_force
         DO i = 1, n_atom
            READ(input_unit, *, iostat=ios) fx, fy, fz
            IF (ios /= 0) THEN
               IF (is_debug) PRINT *, "ERROR| Failed to read gradients from extern_",TRIM(label),".gradients at line ", i
               CALL cp_abort(__LOCATION__, "Failed to read external gradients")
            END IF
            i_kind = kind_of(i)
            i_atom = atom_of_kind(i)
            ! Gradients from file are dE/dR, qs_force%total stores forces (+gradients)
            qs_force(i_kind)%total(1:3, i_atom) = [+fx, +fy, +fz]
         END DO

         ! Close the gradients file
         CALL close_file(unit_number=input_unit)

         ! Clean up temporary arrays
         DEALLOCATE (atom_of_kind, kind_of)
      END IF

      ! Synchronize across all MPI processes to ensure consistent qs_force data
      DO i_kind = 1, SIZE(qs_force)
         CALL fat_env%force_env_ref%para_env%bcast(qs_force(i_kind)%total, &
                                               fat_env%force_env_ref%para_env%source)
      END DO

      CALL timestop(handle)

   END SUBROUTINE fat_get_gradients

! **************************************************************************************************
!> \brief Calculate the fragments total force for a freeze-and-thaw optimization
! **************************************************************************************************
   SUBROUTINE fat_force_frag(fat_env, forces)
      TYPE(fat_env_type), POINTER                        :: fat_env
      REAL(dp), ALLOCATABLE                              :: forces(:, :)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_force_frag'

      INTEGER                                            :: handle, i_frag, i_atom, i_kind, j, k, n
      INTEGER, ALLOCATABLE                               :: atom_of_kind(:), kind_of(:)
      TYPE(atomic_kind_type), POINTER                    :: kinds(:)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_force_type), POINTER                       :: qs_force(:)

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      forces = 0.0_dp

      k = 0
      DO i_frag = 1, fat_env%n_frag
         qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
         CALL get_qs_env(qs_env, atomic_kind_set=kinds, force=qs_force)

         IF (fat_env%extern(i_frag)) THEN
            CALL fat_get_gradients(fat_env, i_frag=i_frag, qs_force=qs_force)
         END IF

         CALL get_atomic_kind_set(kinds, natom=n, kind_of=kind_of, atom_of_kind=atom_of_kind)
         DO j = 1, n
            k = k + 1
            i_kind = kind_of(j)
            i_atom = atom_of_kind(j)
            forces(k, :) = qs_force(i_kind)%total(1:3, i_atom)
         END DO
         DEALLOCATE (atom_of_kind, kind_of)
      END DO

      CALL timestop(handle)

   END SUBROUTINE fat_force_frag

! **************************************************************************************************
!> \brief Calculate the embedding force for a freeze-and-thaw optimization
! **************************************************************************************************
   SUBROUTINE fat_force_emb(fat_env, forces)
      TYPE(fat_env_type), POINTER                        :: fat_env
      REAL(dp), ALLOCATABLE                              :: forces(:, :)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_force_emb'

      INTEGER                                            :: a, handle, i, j, k
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(pw_type), POINTER                             :: drho_g, drho_r, v_r
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      forces = 0.0_dp

      k = 0
      DO i = 1, fat_env%n_frag  ! fragments

         CALL fat_pw_create_frag(fat_env, pw_r=v_r, i_frag=i)

         CALL fat_vemb_frag(fat_env, v_emb=v_r, i_frag=i, to_file=.FALSE.)

         qs_env => fat_env%force_env_sub(i)%force_env%qs_env
         CALL get_qs_env(qs_env, particle_set=particles)

         CALL fat_pw_create_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)

         DO j = 1, SIZE(particles)  ! atoms in the fragment
            k = k + 1  ! atoms in the total system
            DO a = 1, 3  ! xyz

               CALL fat_get_drho_frag(fat_env, drho_r=drho_r, drho_g=drho_g, &
                                      i_frag=i, i_nuc=j, i_dir=a)

               forces(k, a) = pw_integral_ab(v_r, drho_r)

            END DO

         END DO

         CALL fat_pw_release_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)
         CALL fat_pw_release_frag(fat_env, pw_r=v_r, i_frag=i)

      END DO

      CALL timestop(handle)

   END SUBROUTINE fat_force_emb

! **************************************************************************************************
!> \brief Calculate the non-additive kinetic force for a freeze-and-thaw optimization
! **************************************************************************************************
   SUBROUTINE fat_force_tsna(fat_env, forces)
      TYPE(fat_env_type), POINTER                        :: fat_env
      REAL(dp), ALLOCATABLE                              :: forces(:, :)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_force_tsna'

      INTEGER                                            :: a, handle, i, j, k
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(pw_type), POINTER                             :: drho_g, drho_r, rho_g(:), rho_r(:), v_emb
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      forces = 0.0_dp

      k = 0
      DO i = 1, fat_env%n_frag  ! fragments

         ! Create embedding potential
         CALL fat_pw_create_frag(fat_env, pw_r=v_emb, i_frag=i)

         CALL fat_get_rho_frag(fat_env, rho_r=rho_r, rho_g=rho_g, i_frag=i)

         ! Calculate total embedding potential (includes both LDA and GGA if present)
         CALL fat_vemb_tsna_frag(fat_env=fat_env, rho_r=rho_r, rho_g=rho_g, &
                                 v_emb=v_emb, i_frag=i)

         qs_env => fat_env%force_env_sub(i)%force_env%qs_env
         CALL get_qs_env(qs_env, particle_set=particles)

         CALL fat_pw_create_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)

         DO j = 1, SIZE(particles)  ! atoms in the fragment
            k = k + 1
            DO a = 1, 3  ! xyz

               CALL fat_get_drho_frag(fat_env, drho_r=drho_r, drho_g=drho_g, &
                                      i_frag=i, i_nuc=j, i_dir=a)

               ! Total force contribution (LDA + GGA if present)
               forces(k, a) = pw_integral_ab(v_emb, drho_r)

            END DO
         END DO
         CALL fat_pw_release_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)

         ! Release potential
         CALL fat_pw_release_frag(fat_env, pw_r=v_emb, i_frag=i)

      END DO

      CALL timestop(handle)

   END SUBROUTINE fat_force_tsna

! **************************************************************************************************
!> \brief Calculate the non-additive xc force for a freeze-and-thaw optimization
! **************************************************************************************************
   SUBROUTINE fat_force_xcna(fat_env, forces)
      TYPE(fat_env_type), POINTER                        :: fat_env
      REAL(dp), ALLOCATABLE                              :: forces(:, :)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_force_xcna'

      INTEGER                                            :: a, handle, i, j, k
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(pw_type), POINTER                             :: drho_g, drho_r, rho_g(:), rho_r(:), v_emb
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      forces = 0.0_dp

      k = 0
      DO i = 1, fat_env%n_frag  ! fragments

         ! Create embedding potential
         CALL fat_pw_create_frag(fat_env, pw_r=v_emb, i_frag=i)

         CALL fat_get_rho_frag(fat_env, rho_r=rho_r, rho_g=rho_g, i_frag=i)

         ! Calculate total embedding potential (includes both LDA and GGA if present)
         CALL fat_vemb_xcna_frag(fat_env=fat_env, rho_r=rho_r, rho_g=rho_g, &
                                 v_emb=v_emb, i_frag=i)

         qs_env => fat_env%force_env_sub(i)%force_env%qs_env
         CALL get_qs_env(qs_env, particle_set=particles)

         CALL fat_pw_create_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)

         DO j = 1, SIZE(particles)  ! atoms in the fragment
            k = k + 1
            DO a = 1, 3  ! xyz
               CALL fat_get_drho_frag(fat_env, drho_r=drho_r, drho_g=drho_g, &
                                      i_frag=i, i_nuc=j, i_dir=a)

               ! Total force contribution (LDA + GGA if present)
               forces(k, a) = pw_integral_ab(v_emb, drho_r)

            END DO
         END DO

         CALL fat_pw_release_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)

         ! Release potential
         CALL fat_pw_release_frag(fat_env, pw_r=v_emb, i_frag=i)

      END DO

      CALL timestop(handle)

   END SUBROUTINE fat_force_xcna

! **************************************************************************************************
!> \brief Calculate the Hartree force for a freeze-and-thaw optimization
! **************************************************************************************************
   SUBROUTINE fat_force_hartree(fat_env, forces)
      TYPE(fat_env_type), POINTER                        :: fat_env
      REAL(dp), ALLOCATABLE                              :: forces(:, :)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_force_hartree'

      INTEGER                                            :: a, handle, i, j, k
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(pw_type), POINTER                             :: drho_g, drho_r, drhon_g, rho_r_tmp, &
                                                            v_hartree
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CALL timeset(routinen, handle)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      forces = 0.0_dp

      k = 0
      DO i = 1, fat_env%n_frag  ! fragments
         CALL fat_pw_create_frag(fat_env, pw_r=v_hartree, i_frag=i)
         CALL fat_vemb_hartree_frag(fat_env, v_emb=v_hartree, i_frag=i)

         ! Get fragment i_frag's particles
         qs_env => fat_env%force_env_sub(i)%force_env%qs_env

         CALL get_qs_env(qs_env, particle_set=particles)
         CALL fat_pw_create_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)
         CALL fat_pw_create_frag(fat_env, pw_r=rho_r_tmp, pw_g=drhon_g, i_frag=i)

         ! Compute forces
         DO j = 1, SIZE(particles)  ! atoms in the fragment
            k = k + 1
            DO a = 1, 3  ! xyz
               ! Electronic contribution
               CALL fat_get_drho_frag(fat_env, drho_r=drho_r, drho_g=drho_g, &
                                      i_frag=i, i_nuc=j, i_dir=a)
               CALL pw_transfer(drho_g, rho_r_tmp)
               forces(k, a) = pw_integral_ab(v_hartree, rho_r_tmp)

               ! Nuclear contribution
               CALL fat_get_drhon_frag(fat_env, drhon_g=drhon_g, &
                                       i_frag=i, i_nuc=j, i_dir=a)
               CALL pw_transfer(drhon_g, rho_r_tmp)
               forces(k, a) = forces(k, a) + pw_integral_ab(v_hartree, rho_r_tmp)
            END DO
         END DO

         ! Clean up
         CALL fat_pw_release_frag(fat_env, pw_r=drho_r, pw_g=drho_g, i_frag=i)
         CALL fat_pw_release_frag(fat_env, pw_r=rho_r_tmp, pw_g=drhon_g, i_frag=i)
         CALL fat_pw_release_frag(fat_env, pw_r=v_hartree, i_frag=i)
      END DO

      CALL timestop(handle)
   END SUBROUTINE fat_force_hartree

! **************************************************************************************************
!> \brief Print the forces for a freeze-and-thaw optimization
! **************************************************************************************************
   SUBROUTINE fat_force_print(fat_env, forces, label, to_file)
      TYPE(fat_env_type), POINTER                        :: fat_env
      REAL(dp), INTENT(IN)                               :: forces(:, :)
      CHARACTER(len=*), INTENT(IN)                       :: label
      LOGICAL, OPTIONAL                                  :: to_file

      CHARACTER(len=*), PARAMETER :: &
         routineN = 'fat_force_print', &
         fmt_energy = '(/," ENERGY| Total FORCE_EVAL ( FAT ) energy [a.u.]:",T61,F20.15)', &
         fmt_forces = '(2X,I5,2X,I5,6X,A2,T23,3F15.8)', &
         fmt_sum = '(/," SUM OF ATOMIC FORCES",T23,4F15.8,/)'

      CHARACTER(len=2)                                   :: symbol
      INTEGER                                            :: i, log_unit, n, output_unit
      LOGICAL                                            :: my_to_file, is_master
      REAL(KIND=dp), DIMENSION(3)                        :: grand_total
      TYPE(atomic_kind_type), POINTER                    :: kind
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      TYPE(qs_environment_type), POINTER                 :: qs_env

      is_master = fat_env%force_env_ref%para_env%is_source()

      IF (is_debug .AND. is_master) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      my_to_file = .FALSE.
      IF (PRESENT(to_file)) my_to_file = to_file

      qs_env => fat_env%force_env_ref%qs_env
      CALL get_qs_env(qs_env, particle_set=particles)
      n = SIZE(particles)


      IF (is_master .AND. is_debug) THEN
         logger => cp_get_default_logger()
         log_unit = cp_logger_get_default_io_unit(logger)

         IF (log_unit > 0) THEN
            WRITE (log_unit, fmt_energy) fat_env%e_tot
            WRITE (log_unit, fmt='(/," ATOMIC FORCES in [a.u.] (",A,")",/,/," # Atom   Kind   Element",T30,"X",T45,"Y",T60,"Z")') &
               TRIM(label)

            grand_total(:) = 0.0_dp
            DO i = 1, n
               kind => particles(i)%atomic_kind
               CALL get_atomic_kind(kind, element_symbol=symbol)
               WRITE (log_unit, fmt_forces) i, i, symbol, forces(i, 1:3)
               grand_total(:) = grand_total(:) + forces(i, :)
            END DO

            WRITE (log_unit, fmt_sum) grand_total(:), SQRT(SUM(grand_total**2))
         END IF

         IF (my_to_file) THEN
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
         END IF
      END IF

   END SUBROUTINE fat_force_print

! **************************************************************************************************
!> \brief Store the calculated forces into force structures for all subsystems
!> \param fat_env The freeze-and-thaw environment
!> \param forces The calculated forces to store
! **************************************************************************************************
   SUBROUTINE fat_force_store(fat_env)
      TYPE(fat_env_type), POINTER                        :: fat_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_force_store'

      INTEGER                                            :: force_offset, handle, i, i_frag, iatom, &
                                                            ikind, k, natom
      INTEGER, ALLOCATABLE                               :: atom_of_kind(:), kind_of(:), natom_of_kind(:)
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(force_env_type), POINTER                      :: force_env_sub
      TYPE(particle_type), DIMENSION(:), POINTER         :: particles
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_subsys_type), POINTER                      :: subsys

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CALL timeset(routinen, handle)

      ! First update reference system
      qs_env => fat_env%force_env_ref%qs_env
      CALL get_qs_env(qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particles, &
                      force=force)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, &
                               kind_of=kind_of, &
                               natom=natom)

      NULLIFY (force, subsys)
      CALL get_qs_env(qs_env, &
                      force=force, &
                      subsys=subsys)
      IF (.NOT. ASSOCIATED(force)) THEN
         !   *** Allocate the force data structure ***
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom_of_kind=natom_of_kind)
         CALL allocate_qs_force(force, natom_of_kind)
         DEALLOCATE (natom_of_kind)
         CALL qs_subsys_set(subsys, force=force)
      END IF
      CALL zero_qs_force(force)

      ! Store forces for reference system
      DO iatom = 1, natom

         ikind = kind_of(iatom)
         i = atom_of_kind(iatom)

         ! Store in force structure
         force(ikind)%total(1:3, i) = fat_env%forces(iatom, 1:3)

         ! Store in particle_set (negative sign because forces = -dE/dR)
         particles(iatom)%f(:) = -fat_env%forces(iatom, :)
      END DO

      DEALLOCATE (atom_of_kind, kind_of)

      ! Now update forces for each subsystem
      force_offset = 0
      DO i_frag = 1, fat_env%n_frag

         force_env_sub => fat_env%force_env_sub(i_frag)%force_env
         qs_env => force_env_sub%qs_env

         ! Get subsystem structures
         CALL get_qs_env(qs_env, &
                         atomic_kind_set=atomic_kind_set, &
                         particle_set=particles, &
                         force=force)

         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  atom_of_kind=atom_of_kind, &
                                  kind_of=kind_of, &
                                  natom=natom)

         ! Store forces for this subsystem
         k = force_offset
         DO iatom = 1, natom
            k = k + 1
            ikind = kind_of(iatom)
            i = atom_of_kind(iatom)

            ! Store in force structure
            force(ikind)%total(1:3, i) = fat_env%forces(k, 1:3)

            ! Store in particle_set (negative sign because forces = -dE/dR)
            particles(iatom)%f(1:3) = -fat_env%forces(k, 1:3)
         END DO

         force_offset = force_offset + natom
         DEALLOCATE (atom_of_kind, kind_of)

      END DO

      CALL timestop(handle)

   END SUBROUTINE fat_force_store

! **************************************************************************************************
!> \brief Calculate the forces for a freeze-and-thaw optimization
! **************************************************************************************************
   SUBROUTINE fat_force(fat_env)
      TYPE(fat_env_type), POINTER                        :: fat_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_force'

      INTEGER                                            :: handle, i_frag
      REAL(kind=dp), ALLOCATABLE                         :: forces_tmp(:, :)

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      CALL timeset(routinen, handle)

      fat_env%forces = 0.0_dp
      CALL fat_force_create(fat_env, forces_tmp)

      DO i_frag = 1, fat_env%n_frag
         IF (fat_env%extern(i_frag)) THEN
            CALL fat_run_external(fat_env, i_frag, gradient=.TRUE.)
         END IF
      END DO

      CALL fat_force_frag(fat_env, forces_tmp)
      CALL fat_force_print(fat_env, forces=-forces_tmp, label="force_frag")
      fat_env%forces(:,:) = fat_env%forces(:,:) + forces_tmp(:,:)

      CALL fat_force_emb(fat_env, forces_tmp)
      CALL fat_force_print(fat_env, forces=+forces_tmp, label="force_emb")
      fat_env%forces(:, :) = fat_env%forces(:, :) - forces_tmp(:, :)

      CALL fat_force_tsna(fat_env, forces_tmp)
      CALL fat_force_print(fat_env, forces=-forces_tmp, label="force_tsna")
      fat_env%forces(:, :) = fat_env%forces(:, :) + forces_tmp(:, :)

      CALL fat_force_xcna(fat_env, forces_tmp)
      CALL fat_force_print(fat_env, forces=-forces_tmp, label="force_xcna")
      fat_env%forces(:, :) = fat_env%forces(:, :) + forces_tmp(:, :)

      CALL fat_force_hartree(fat_env, forces_tmp)
      CALL fat_force_print(fat_env, forces=-forces_tmp, label="force_hartree")
      fat_env%forces(:, :) = fat_env%forces(:, :) + forces_tmp(:, :)

      CALL fat_force_print(fat_env, forces=-fat_env%forces, label="total", to_file=.FALSE.)

      DEALLOCATE (forces_tmp)

      CALL timestop(handle)

   END SUBROUTINE fat_force

! **************************************************************************************************
!> \brief Print the end results of a freeze-and-thaw optimization
! **************************************************************************************************
   SUBROUTINE fat_results(fat_env)
      TYPE(fat_env_type), POINTER                        :: fat_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_results'

      INTEGER                                            :: i, j
      REAL(kind=dp)                                      :: e_tmp(fat_env%n_frag, fat_env%n_frag), tmp

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      IF (fat_env%force_env_ref%para_env%is_source()) THEN
         DO i = 1, fat_env%n_frag
            WRITE (fat_env%output_unit, &
                   '(T2,"FAT| Tot energy fragment   ",I2,":",F26.15)') i, fat_env%e_sub(i)
         END DO
         WRITE (fat_env%output_unit, &
                '(T2,"FAT| Tot energy of fragments [a.u.]:",F26.15)')  sum(fat_env%e_sub)

         DO i = 1, fat_env%n_frag
            WRITE (fat_env%output_unit, &
                   '(T2,"FAT| Emb energy fragment    ",I2,":",F26.15)') i, fat_env%e_emb(i)
         END DO
         WRITE (fat_env%output_unit, &
                '(T2,"FAT| Tot Emb energy [kcal/mol]:",F26.5)')  sum(fat_env%e_emb)*627.509474063112121

         WRITE (fat_env%output_unit, &
                '(T2,"FAT| Kin energy reference     :",F26.15)') fat_env%ekin_ref
         DO i = 1, fat_env%n_frag
            WRITE (fat_env%output_unit, &
                   '(T2,"FAT| Kin energy fragment    ",I2,":",F26.15)') i, fat_env%ekin_sub(i)
         END DO
         WRITE (fat_env%output_unit, &
                '(T2,"FAT| Tot Kin energy [kcal/mol]:",F26.5)') (fat_env%ekin_ref-sum(fat_env%ekin_sub))*627.509474063112121

         WRITE (fat_env%output_unit, &
                '(T2,"FAT| XC energy reference      :",F26.15)') fat_env%exc_ref
         DO i = 1, fat_env%n_frag
            WRITE (fat_env%output_unit, &
                   '(T2,"FAT| XC energy fragment     ",I2,":",F26.15)') i, fat_env%exc_sub(i)
         END DO
         WRITE (fat_env%output_unit, &
                '(T2,"FAT| Tot xc energy [kcal/mol]:",F26.5)') (fat_env%exc_ref - sum(fat_env%exc_sub))*627.509474063112121
      END IF

      DO i = 1, fat_env%n_frag
         DO j = i+1, fat_env%n_frag
            e_tmp(i, j) = fat_eint_hartree(fat_env, i, j)
            tmp = tmp + e_tmp(i, j)
         END DO
      END DO

      IF (fat_env%force_env_ref%para_env%is_source()) THEN
         DO i = 1, fat_env%n_frag
            DO j = i+1, fat_env%n_frag
               WRITE (fat_env%output_unit, &
                      '(T2,"FAT| Hartree energy fragments ",I2," ",I2,":",F26.15)') &
                      i, j, e_tmp(i, j)
            END DO
         END DO
         WRITE (fat_env%output_unit, &
                '(T2,"FAT| Tot El.static energy [kcal/mol]:",F26.5)') tmp*627.509474063112121
      END IF
      IF (fat_env%force_env_ref%para_env%is_source()) THEN
         WRITE (fat_env%output_unit, &
                '(/,T2,"FAT| Total Interaction energy [kcal/mol]: ",T55,F26.5)') &
               (fat_env%e_tot-sum(fat_env%e_sub))*627.509474063112121
         WRITE (fat_env%output_unit, &
                '(/,T2,"FAT| Total Embedding energy [kcal/mol]: ",T55,F26.5)') &
               sum(fat_env%e_emb)*627.509474063112121
      END IF

      IF (fat_env%force_env_ref%para_env%is_source()) THEN
         WRITE (fat_env%output_unit, &
                '(/,T2,"ENERGY| Total FORCE_EVAL ( ",A," ) energy [a.u.]: ",T55,F26.15)') &
               "FAT", fat_env%e_tot
      END IF

   END SUBROUTINE fat_results

! **************************************************************************************************
!> \brief Cleans up a freeze-and-thaw optimization
! **************************************************************************************************
   SUBROUTINE fat_finalize(fat_env)
      TYPE(fat_env_type), POINTER                        :: fat_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_finalize'

      INTEGER                                            :: i
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(force_env_type), POINTER                      :: force_env_tmp

      IF (is_debug .AND. fat_env%force_env_ref%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      IF (ALLOCATED(fat_env%forces)) DEALLOCATE(fat_env%forces)
      DO i = 1, fat_env%n_frag
         qs_env => fat_env%force_env_sub(i)%force_env%qs_env
         CALL get_qs_env(qs_env, dft_control=dft_control)
         dft_control%apply_embed_pot = .FALSE.

         ! Release the para_env
         force_env_tmp => fat_env%force_env_sub(i)%force_env
         IF (ASSOCIATED(force_env_tmp%para_env)) THEN
            CALL force_env_tmp%para_env%free()
         END IF
      END DO

   END SUBROUTINE fat_finalize

! **************************************************************************************************
!> \brief Initialize the freeze-and-thaw environment
! **************************************************************************************************
   SUBROUTINE fat_init(fat_env, force_env, ref_subsys_number)
      TYPE(fat_env_type), POINTER                        :: fat_env
      TYPE(force_env_type), POINTER                      :: force_env
      INTEGER                                            :: ref_subsys_number

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fat_init'

      INTEGER                                            :: handle, i, j, n_iter
      LOGICAL                                            :: ionode
      REAL(kind=dp)                                      :: eps
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(force_env_type), POINTER                      :: force_env_tmp
      TYPE(pw_env_type), POINTER                         :: pw_env_tmp
      TYPE(scf_control_type), POINTER                    :: scf_control_tmp
      TYPE(qs_environment_type), POINTER                 :: qs_env_tmp
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routinen, handle)

      IF (is_debug .AND. &
          force_env%sub_force_env(ref_subsys_number)%force_env%para_env%is_source()) THEN
         PRINT *, "DBG| ", modulen//"::"//routinen
      END IF

      NULLIFY (qs_env_tmp)

      ! This is a limitation of the EMBED method
      CPASSERT(ref_subsys_number == 1)

      fat_env%force_env => force_env
      fat_env%force_env_ref => force_env%sub_force_env(ref_subsys_number)%force_env

      ! Count fragments
      fat_env%n_frag = SIZE(force_env%sub_force_env) - 1

      ! Get logger and initialize output unit only for I/O node
      ionode = force_env%para_env%mepos == force_env%para_env%source
      IF (ionode) THEN
         logger => cp_get_default_logger()
         fat_env%output_unit = cp_logger_get_default_io_unit(logger)
      ELSE
         fat_env%output_unit = -1
      END IF

      ! Allocate and setup force environments for each fragment
      ALLOCATE (fat_env%force_env_sub(fat_env%n_frag))
      j = 1
      DO i = 1, SIZE(force_env%sub_force_env)
         IF (i == ref_subsys_number) CYCLE ! Skip reference system

         ! Store pointer to fragment's force environment
         fat_env%force_env_sub(j)%force_env => force_env%sub_force_env(i)%force_env
         j = j + 1
      END DO

      ! Setup plane wave pools for density/potential calculations
      ! For reference system
      force_env_tmp => fat_env%force_env_ref
      qs_env_tmp => force_env_tmp%qs_env
      CALL get_qs_env(qs_env_tmp, pw_env=pw_env_tmp)
      CALL pw_env_get(pw_env_tmp, auxbas_pw_pool=fat_env%pw_pool_ref)

      ! For each fragment
      ALLOCATE (fat_env%pw_pool_sub(fat_env%n_frag))
      DO i = 1, fat_env%n_frag
         force_env_tmp => fat_env%force_env_sub(i)%force_env
         qs_env_tmp => fat_env%force_env_sub(i)%force_env%qs_env
         CALL get_qs_env(qs_env_tmp, pw_env=pw_env_tmp)
         CALL pw_env_get(pw_env_tmp, auxbas_pw_pool=fat_env%pw_pool_sub(i)%pool)
      END DO

      ! Setup flags for external fragments, if specified
      ALLOCATE (fat_env%extern(fat_env%n_frag), source=.FAlSE.)
      DO i = 1, fat_env%n_frag
         force_env_tmp => fat_env%force_env_sub(i)%force_env
         qs_env_tmp => force_env_tmp%qs_env
         CALL get_qs_env(qs_env_tmp, input=input)
         ! Check if this fragment uses external calculation
         CALL section_vals_val_get(input, "DFT%QS%FAT_EXTERN", l_val=fat_env%extern(i))

         IF (is_debug .AND. ionode) &
            PRINT *, "DBG| fat_env%extern(", i, ") = ", fat_env%extern(i)
      END DO

      ! Setup the same MPI Group for all
      DO i = 1, fat_env%n_frag
         ! All fragments use the same MPI group as the reference system
         force_env_tmp => fat_env%force_env_sub(i)%force_env
         force_env_tmp%para_env => fat_env%force_env_ref%para_env
         ! Retain the reference to increase counter
         CALL force_env_tmp%para_env%retain()
      END DO

      ! Initial DFT calculations for external fragments (without embedding potential)
      do i = 1, fat_env%n_frag
         IF (.NOT. fat_env%extern(i)) cycle
         force_env_tmp => fat_env%force_env_sub(i)%force_env
         qs_env_tmp => force_env_tmp%qs_env
         CALL get_qs_env(qs_env_tmp, scf_control=scf_control_tmp)
!        scf_control_tmp%max_scf = 1
!        scf_control_tmp%ignore_convergence_failure = .TRUE.
         CALL qs_calc_energy_force(qs_env=qs_env_tmp, &
                                   calc_force=.TRUE., &
                                   consistent_energies=.FALSE., &
                                   linres=.FALSE.)
      end do

      ! Allocate arrays for storing various energy components
      ALLOCATE (fat_env%e_sub(fat_env%n_frag)) ! Subsystem energies
      fat_env%e_sub = 0.0

      ALLOCATE (fat_env%e_emb(fat_env%n_frag)) ! Embedding energies
      fat_env%e_emb = 0.0

      fat_env%ekin_ref = 0.0                     ! Reference kinetic energy
      ALLOCATE (fat_env%ekin_sub(fat_env%n_frag)) ! Subsystem kinetic energies
      fat_env%ekin_sub = 0.0

      fat_env%exc_ref = 0.0                     ! Reference XC energy
      ALLOCATE (fat_env%exc_sub(fat_env%n_frag)) ! Subsystem XC energies
      fat_env%exc_sub = 0.0

      ! Initialize interaction electrostatic energies
      fat_env%ehartree = 0.0 ! Electrostatic energy

      ! Initialize "old" energy for convergence check
      fat_env%e_old = 0.0

      ! Initialize iteration counter
      fat_env%i_iter = 0
      ! Calculate initial energies for each fragment
!     do i = 1, fat_env%n_frag
!        call fat_iterate_frag(fat_env, i)
!     end do
!     if (is_debug) print *, "KKK| After initial calculations for all fragments i_iter", i, fat_env%i_iter

      ! Get convergence parameters from input
      CALL get_qs_env(fat_env%force_env_ref%qs_env, input=input)
      CALL section_vals_val_get(input, "DFT%QS%OPT_EMBED%N_ITER", i_val=n_iter)
      fat_env%n_iter = n_iter
      CALL section_vals_val_get(input, "DFT%QS%OPT_EMBED%DENS_CONV_MAX", r_val=eps)
      fat_env%eps = eps

      ! Setup header for status output
      IF (ionode) &
         WRITE (fat_env%output_unit, fmt=fat_title_fmt) &
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
      CALL fat_iterate(fat_env)

      CALL timestop(handle)

   END SUBROUTINE fat_init

! **************************************************************************************************
!> \brief Main driver for Freeze-and-thaw frozen density embedding
!> \param force_env the force environment.
!> \param converged_embed true if the optimization has converged
!> \author Lukas Schreder
! **************************************************************************************************
subroutine fat_embedding(force_env, converged, do_forces)
      type(force_env_type), pointer                      :: force_env
      logical                                            :: converged, do_forces

      character(len=*), parameter                        :: routinen = 'fat_embedding'

      integer :: handle, iparticle
      type(fat_env_type), pointer :: fat_env
      TYPE(cp_subsys_type), POINTER :: subsys_embed
      TYPE(particle_list_type), POINTER :: particles_embed

      call timeset(routinen, handle)

      ! Get access to the particles in the embedding environment
      CALL force_env_get(force_env=force_env, subsys=subsys_embed)
      CALL cp_subsys_get(subsys=subsys_embed, particles=particles_embed)

      call fat_env_create(fat_env)
      call fat_init(fat_env, force_env, ref_subsys_number=1) ! Ref subsystem number is always 1.

      converged = .false.
      do while (.not. fat_is_done(fat_env))
         call fat_iterate(fat_env)
      end do
      converged = fat_env%is_converged

      IF (.NOT. converged) THEN
         CALL cp_warn(__LOCATION__, &
                      "FAT embedding did not converge within the maximum number of iterations.")
      END IF

      CALL fat_env%force_env_ref%para_env%sync()

      CALL fat_results(fat_env)


      ! FAT total energy - fragments total energy = embedding energy
      force_env%embed_env%pot_energy = fat_env%e_tot - sum(fat_env%e_sub)

      force_env%embed_env%energies(2:) = fat_env%e_sub

      ! Calculate and save forces
      IF (do_forces) THEN
         ! Allocate forces array if it doesn't exist yet
         IF (.NOT. ALLOCATED(fat_env%forces)) THEN
            CALL fat_force_create(fat_env, fat_env%forces)
         END IF

         ! Calculate forces and store them in fat_env%forces
         CALL fat_force(fat_env)

         CALL fat_force_store(fat_env)
         ! Save forces in the particle array of the embedding environment
         DO iparticle = 1, SIZE(particles_embed%els)
            particles_embed%els(iparticle)%f(:) = -fat_env%forces(iparticle,:)
         END DO

         IF (fat_env%force_env_ref%para_env%is_source()) THEN
            WRITE (fat_env%output_unit, &
                   '(/,T2,"ENERGY| Total FORCE_EVAL ( ",A," ) energy [a.u.]: ",T55,F26.15)') &
                  "FAT", fat_env%e_tot
         END IF

      END IF

      call fat_finalize(fat_env)
      call fat_env_release(fat_env)

      call timestop(handle)

    end subroutine fat_embedding

END MODULE fat_methods
