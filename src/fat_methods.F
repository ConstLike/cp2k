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
  USE pw_grid_types, ONLY: pw_grid_type
  USE pw_methods, ONLY: pw_zero, pw_axpy, pw_transfer, pw_integral_ab
  use pw_poisson_methods, only: pw_poisson_solve
  use pw_poisson_types, only: pw_poisson_type
  USE pw_pool_types, ONLY: pw_pool_type, pw_pool_p_type
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
!> \brief Creates a real space or g-space pw type for one fragment in the fat environment
!> \param fat_env the freeze-and-thaw environment containing calculation setup
!> \param pw_rs the pw to be created in real space grid (optional)
!> \param pw_gs the pw to be created in g-space grid (optional)
!> \param i_frag index of the current fragment
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
!> \note
!>      At least one of pw_rs or pw_gs must be present
! **************************************************************************************************
   SUBROUTINE fat_pw_create_frag(fat_env, pw_rs, pw_gs, i_frag)
      TYPE(fat_env_type), POINTER                       :: fat_env
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER           :: pw_rs
      TYPE(pw_c1d_gs_type), OPTIONAL, POINTER           :: pw_gs
      INTEGER, INTENT(IN)                               :: i_frag

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_pw_create_frag'

      INTEGER                                           :: handle

      CALL timeset(routineN, handle)

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      CPASSERT(PRESENT(pw_rs) .OR. PRESENT(pw_gs))
      CPASSERT(ASSOCIATED(fat_env%pw_pool_sub))
      CPASSERT(i_frag >= 1 .AND. i_frag <= SIZE(fat_env%pw_pool_sub))
      CPASSERT(ASSOCIATED(fat_env%pw_pool_sub(i_frag)%pool))

      IF (PRESENT(pw_rs)) THEN
         ALLOCATE(pw_rs)
         CALL fat_env%pw_pool_sub(i_frag)%pool%create_pw(pw_rs)
         CALL pw_zero(pw_rs)
      END IF

      IF (PRESENT(pw_gs)) THEN
         ALLOCATE(pw_gs)
         CALL fat_env%pw_pool_sub(i_frag)%pool%create_pw(pw_gs)
         CALL pw_zero(pw_gs)
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pw_create_frag

! **************************************************************************************************
!> \brief Releases real space or g-space pw type for one fragment in the fat environment
!> \param fat_env the freeze-and-thaw environment containing calculation setup
!> \param pw_rs the real space pw to be released from real space grid (optional)
!> \param pw_gs the g-space pw to be released from g-space grid (optional)
!> \param i_frag index of the current fragment
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
!> \note
!>      At least one of pw_rs or pw_gs must be present
! **************************************************************************************************
   SUBROUTINE fat_pw_release_frag(fat_env, pw_rs, pw_gs, i_frag)
      TYPE(fat_env_type), POINTER                       :: fat_env
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER           :: pw_rs
      TYPE(pw_c1d_gs_type), OPTIONAL, POINTER           :: pw_gs
      INTEGER, INTENT(IN)                               :: i_frag

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_pw_release_frag'

      INTEGER                                           :: handle

      CALL timeset(routineN, handle)

      CPASSERT(PRESENT(pw_rs) .OR. PRESENT(pw_gs))
      CPASSERT(ASSOCIATED(fat_env%pw_pool_sub))
      CPASSERT(i_frag >= 1 .AND. i_frag <= SIZE(fat_env%pw_pool_sub))
      CPASSERT(ASSOCIATED(fat_env%pw_pool_sub(i_frag)%pool))

      IF (PRESENT(pw_rs)) THEN
         IF (ASSOCIATED(pw_rs)) THEN
            CALL pw_rs%release()
            CALL fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(pw_rs)
            DEALLOCATE(pw_rs)
            NULLIFY(pw_rs)
         END IF
      END IF

      IF (PRESENT(pw_gs)) THEN
         IF (ASSOCIATED(pw_gs)) THEN
            CALL pw_gs%release()
            CALL fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(pw_gs)
            DEALLOCATE(pw_gs)
            NULLIFY(pw_gs)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pw_release_frag

! **************************************************************************************************
!> \brief Creates a real space or g-space pw type for the reference system
!> \param fat_env the freeze-and-thaw environment containing calculation setup
!> \param pw_rs the pw to be created in real space grid (optional)
!> \param pw_gs the pw to be created in g-space grid (optional)
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
!> \note
!>      At least one of pw_rs or pw_gs must be present
! **************************************************************************************************
   SUBROUTINE fat_pw_create_ref(fat_env, pw_rs, pw_gs)
      TYPE(fat_env_type), POINTER                       :: fat_env
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER           :: pw_rs
      TYPE(pw_c1d_gs_type), OPTIONAL, POINTER           :: pw_gs

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_pw_create_ref'

      INTEGER                                           :: handle

      CALL timeset(routineN, handle)

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      CPASSERT(PRESENT(pw_rs) .OR. PRESENT(pw_gs))
      CPASSERT(ASSOCIATED(fat_env%pw_pool_ref))

      IF (PRESENT(pw_rs)) THEN
         ALLOCATE(pw_rs)
         CALL fat_env%pw_pool_ref%create_pw(pw_rs)
         CALL pw_zero(pw_rs)
      END IF

      IF (PRESENT(pw_gs)) THEN
         ALLOCATE(pw_gs)
         CALL fat_env%pw_pool_ref%create_pw(pw_gs)
         CALL pw_zero(pw_gs)
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pw_create_ref

! **************************************************************************************************
!> \brief Releases real space or g-space pw type for the reference system
!> \param fat_env the freeze-and-thaw environment containing calculation setup
!> \param pw_rs the real space pw to be released (optional)
!> \param pw_gs the g-space pw to be released (optional)
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
!> \note
!>      At least one of pw_rs or pw_gs must be present
! **************************************************************************************************
   SUBROUTINE fat_pw_release_ref(fat_env, pw_rs, pw_gs)
      TYPE(fat_env_type), POINTER                       :: fat_env
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER           :: pw_rs
      TYPE(pw_c1d_gs_type), OPTIONAL, POINTER           :: pw_gs

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_pw_release_ref'

      INTEGER                                           :: handle

      CALL timeset(routineN, handle)

      CPASSERT(PRESENT(pw_rs) .OR. PRESENT(pw_gs))
      CPASSERT(ASSOCIATED(fat_env%pw_pool_ref))

      IF (PRESENT(pw_rs)) THEN
         IF (ASSOCIATED(pw_rs)) THEN
            CALL pw_rs%release()
            CALL fat_env%pw_pool_ref%give_back_pw(pw_rs)
            DEALLOCATE(pw_rs)
            NULLIFY(pw_rs)
         END IF
      END IF

      IF (PRESENT(pw_gs)) THEN
         IF (ASSOCIATED(pw_gs)) THEN
            CALL pw_gs%release()
            CALL fat_env%pw_pool_ref%give_back_pw(pw_gs)
            DEALLOCATE(pw_gs)
            NULLIFY(pw_gs)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pw_release_ref

! **************************************************************************************************
!> \brief Creates an array of real space or g-space pw types for a fragment in the fat environment
!> \param fat_env the freeze-and-thaw environment containing calculation setup
!> \param pws_rs array of pw's to be created in real space grid (optional)
!> \param pws_gs array of pw's to be created in g-space grid (optional)
!> \param i_frag index of the current fragment
!> \param n size of array to create
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
!> \note
!>      At least one of pws_rs or pws_gs must be present
! **************************************************************************************************
   SUBROUTINE fat_pws_create_frag(fat_env, pws_rs, pws_gs, i_frag, n)
      TYPE(fat_env_type), POINTER                       :: fat_env
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER           :: pws_rs(:)
      TYPE(pw_c1d_gs_type), OPTIONAL, POINTER           :: pws_gs(:)
      INTEGER, INTENT(IN)                               :: i_frag, n

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_pws_create_frag'

      TYPE(pw_r3d_rs_type), POINTER                     :: pw_rs_tmp
      TYPE(pw_c1d_gs_type), POINTER                     :: pw_gs_tmp
      INTEGER                                           :: handle, i

      CALL timeset(routineN, handle)

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      CPASSERT(PRESENT(pws_rs) .OR. PRESENT(pws_gs))

      IF (PRESENT(pws_rs)) THEN
         ALLOCATE(pws_rs(n))
         DO i = 1, n
            pw_rs_tmp => pws_rs(i)
            CALL fat_pw_create_frag(fat_env, pw_rs=pw_rs_tmp, i_frag=i_frag)
         END DO
      END IF

      IF (PRESENT(pws_gs)) THEN
         ALLOCATE(pws_gs(n))
         DO i = 1, n
            pw_gs_tmp => pws_gs(i)
            CALL fat_pw_create_frag(fat_env, pw_gs=pw_gs_tmp, i_frag=i_frag)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pws_create_frag

! **************************************************************************************************
!> \brief Releases an array of real space or g-space pw types for a fragment in the fat environment
!> \param fat_env the freeze-and-thaw environment containing calculation setup
!> \param pws_rs array of pw's to be released from real space grid (optional)
!> \param pws_gs array of pw's to be released from g-space grid (optional)
!> \param i_frag index of the current fragment
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
!> \note
!>      At least one of pws_rs or pws_gs must be present
! **************************************************************************************************
   SUBROUTINE fat_pws_release_frag(fat_env, pws_rs, pws_gs, i_frag)
      TYPE(fat_env_type), POINTER                       :: fat_env
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER           :: pws_rs(:)
      TYPE(pw_c1d_gs_type), OPTIONAL, POINTER           :: pws_gs(:)
      INTEGER, INTENT(IN)                               :: i_frag

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_pws_release_frag'

      TYPE(pw_r3d_rs_type), POINTER                     :: pw_rs_tmp
      TYPE(pw_c1d_gs_type), POINTER                     :: pw_gs_tmp
      INTEGER                                           :: handle, i

      CALL timeset(routineN, handle)

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      CPASSERT(PRESENT(pws_rs) .OR. PRESENT(pws_gs))

      IF (PRESENT(pws_rs)) THEN
         CPASSERT(ASSOCIATED(pws_rs))
         DO i = 1, SIZE(pws_rs)
            pw_rs_tmp => pws_rs(i)
            CALL fat_pw_release_frag(fat_env, pw_rs=pw_rs_tmp, i_frag=i_frag)
         END DO
         DEALLOCATE(pws_rs)
         NULLIFY(pws_rs)
      END IF

      IF (PRESENT(pws_gs)) THEN
         CPASSERT(ASSOCIATED(pws_gs))
         DO i = 1, SIZE(pws_gs)
            pw_gs_tmp => pws_gs(i)
            CALL fat_pw_release_frag(fat_env, pw_gs=pw_gs_tmp, i_frag=i_frag)
         END DO
         DEALLOCATE(pws_gs)
         NULLIFY(pws_gs)
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pws_release_frag

! **************************************************************************************************
!> \brief Creates an array of real space or g-space pw types for the reference system
!> \param fat_env the freeze-and-thaw environment containing calculation setup
!> \param pws_rs array of pw's to be created in real space grid (optional)
!> \param pws_gs array of pw's to be created in g-space grid (optional)
!> \param n size of array to create
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
!> \note
!>      At least one of pws_rs or pws_gs must be present
! **************************************************************************************************
   SUBROUTINE fat_pws_create_ref(fat_env, pws_rs, pws_gs, n)
      TYPE(fat_env_type), POINTER                       :: fat_env
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER           :: pws_rs(:)
      TYPE(pw_c1d_gs_type), OPTIONAL, POINTER           :: pws_gs(:)
      INTEGER, INTENT(IN)                               :: n

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_pws_create_ref'

      TYPE(pw_r3d_rs_type), POINTER                     :: pw_rs_tmp
      TYPE(pw_c1d_gs_type), POINTER                     :: pw_gs_tmp
      INTEGER                                           :: handle, i

      CALL timeset(routineN, handle)

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      CPASSERT(PRESENT(pws_rs) .OR. PRESENT(pws_gs))

      IF (PRESENT(pws_rs)) THEN
         ALLOCATE(pws_rs(n))
         DO i = 1, n
            pw_rs_tmp => pws_rs(i)
            CALL fat_pw_create_ref(fat_env, pw_rs=pw_rs_tmp)
         END DO
      END IF

      IF (PRESENT(pws_gs)) THEN
         ALLOCATE(pws_gs(n))
         DO i = 1, n
            pw_gs_tmp => pws_gs(i)
            CALL fat_pw_create_ref(fat_env, pw_gs=pw_gs_tmp)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pws_create_ref

! **************************************************************************************************
!> \brief Releases an array of real space or g-space pw types for the reference system
!> \param fat_env the freeze-and-thaw environment containing calculation setup
!> \param pws_rs array of pw's to be released from real space grid (optional)
!> \param pws_gs array of pw's to be released from g-space grid (optional)
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
!> \note
!>      At least one of pws_rs or pws_gs must be present
! **************************************************************************************************
   SUBROUTINE fat_pws_release_ref(fat_env, pws_rs, pws_gs)
      TYPE(fat_env_type), POINTER                       :: fat_env
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER           :: pws_rs(:)
      TYPE(pw_c1d_gs_type), OPTIONAL, POINTER           :: pws_gs(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_pws_release_ref'

      TYPE(pw_r3d_rs_type), POINTER                     :: pw_rs_tmp
      TYPE(pw_c1d_gs_type), POINTER                     :: pw_gs_tmp
      INTEGER                                           :: handle, i

      CALL timeset(routineN, handle)

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      CPASSERT(PRESENT(pws_rs) .OR. PRESENT(pws_gs))

      IF (PRESENT(pws_rs)) THEN
         CPASSERT(ASSOCIATED(pws_rs))
         DO i = 1, SIZE(pws_rs)
            pw_rs_tmp => pws_rs(i)
            CALL fat_pw_release_ref(fat_env, pw_rs=pw_rs_tmp)
         END DO
         DEALLOCATE(pws_rs)
         NULLIFY(pws_rs)
      END IF

      IF (PRESENT(pws_gs)) THEN
         CPASSERT(ASSOCIATED(pws_gs))
         DO i = 1, SIZE(pws_gs)
            pw_gs_tmp => pws_gs(i)
            CALL fat_pw_release_ref(fat_env, pw_gs=pw_gs_tmp)
         END DO
         DEALLOCATE(pws_gs)
         NULLIFY(pws_gs)
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_pws_release_ref

! **************************************************************************************************
!> \brief Gets electron density for one fragment in real and/or g-space
!> \param fat_env the freeze-and-thaw environment containing calculation setup
!> \param rho_rs the real space density arrays to be filled (optional)
!> \param rho_gs the g-space density arrays to be filled (optional)
!> \param i_frag index of the current fragment
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
!> \note
!>      At least one of rho_rs or rho_gs must be present
!>      For external fragments, density is read from cube file
! **************************************************************************************************
   SUBROUTINE fat_rho_frag(fat_env, rho_rs, rho_gs, i_frag)
      TYPE(fat_env_type), POINTER                       :: fat_env
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER           :: rho_rs(:)
      TYPE(pw_c1d_gs_type), OPTIONAL, POINTER           :: rho_gs(:)
      INTEGER, INTENT(IN)                               :: i_frag

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_rho_frag'
      CHARACTER(LEN=default_string_length)              :: label
      TYPE(qs_environment_type), POINTER                :: qs_env
      TYPE(qs_rho_type), POINTER                        :: qs_rho
      TYPE(pw_r3d_rs_type), POINTER                     :: rho_rs_tmp
      INTEGER                                           :: handle

      CALL timeset(routineN, handle)

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      CPASSERT(PRESENT(rho_rs) .OR. PRESENT(rho_gs))

      qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
      CALL get_qs_env(qs_env, rho=qs_rho)

      IF (PRESENT(rho_rs)) THEN
         CALL qs_rho_get(qs_rho, rho_r=rho_rs)
      END IF

      IF (PRESENT(rho_gs)) THEN
         CALL qs_rho_get(qs_rho, rho_g=rho_gs)
      END IF

      IF (fat_env%extern(i_frag)) THEN
         CALL integer_to_string(i_frag, label)
         IF (PRESENT(rho_rs)) THEN
            CALL cp_cube_to_pw(rho_rs(1), "extern_"//TRIM(label)//".cube", 1.0_dp, .FALSE.)
         END IF
         IF (PRESENT(rho_gs)) THEN
            CALL fat_pw_create_frag(fat_env, pw_rs=rho_rs_tmp, i_frag=i_frag)
            CALL cp_cube_to_pw(rho_rs_tmp, "extern_"//TRIM(label)//".cube", 1.0_dp, .FALSE.)
            CALL pw_transfer(rho_rs_tmp, rho_gs(1))
            CALL fat_pw_release_frag(fat_env, pw_rs=rho_rs_tmp, i_frag=i_frag)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE fat_rho_frag

! **************************************************************************************************
!> \brief Get the total electron density for the freeze-and-thaw environment
!> \param fat_env the freeze-and-thaw environment
!> \param rho_rs real space density array to be filled
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
! **************************************************************************************************
SUBROUTINE fat_rho_ref(fat_env, rho_rs)
   TYPE(fat_env_type), POINTER                      :: fat_env
   TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER      :: rho_rs

   CHARACTER(len=*), PARAMETER :: routineN = 'fat_rho_ref'

   TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER      :: rho_rs_tmp
   INTEGER                                          :: handle, i, j

   CALL timeset(routineN, handle)

   IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

   CPASSERT(ASSOCIATED(fat_env) .AND. ASSOCIATED(rho_rs))

   DO i = 1, SIZE(rho_rs)
      CALL pw_zero(rho_rs(i))
   END DO

   ! Sum up densities from all fragments
   DO i = 1, fat_env%n_frag
     CALL fat_rho_frag(fat_env, rho_rs=rho_rs_tmp, i_frag=i)
      DO j = 1, SIZE(rho_rs_tmp)
         CALL pw_axpy(rho_rs_tmp(j), rho_rs(j), allow_noncompatible_grids=.TRUE.)
      END DO
   END DO

   CALL timestop(handle)

END SUBROUTINE fat_rho_ref

! **************************************************************************************************
!> \brief Gets nuclear density in g-space for a fragment
!> \param fat_env the freeze-and-thaw environment
!> \param rhon_gs nuclear density in g-space
!> \param i_frag fragment index
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
! **************************************************************************************************
   SUBROUTINE fat_rhon_frag(fat_env, rhon_gs, i_frag)
      TYPE(fat_env_type), POINTER                       :: fat_env
      TYPE(pw_c1d_gs_type), POINTER                     :: rhon_gs
      INTEGER, INTENT(IN)                               :: i_frag

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_rhon_frag'

      TYPE(qs_environment_type), POINTER                :: qs_env
      INTEGER                                           :: handle

      CALL timeset(routineN, handle)

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
      CALL get_qs_env(qs_env=qs_env, rho_core=rhon_gs)

      CALL timestop(handle)

   END SUBROUTINE fat_rhon_frag

! **************************************************************************************************
!> \brief Gets derivative of electron density for fragment wrt nuclear position
!> \param fat_env the freeze-and-thaw environment containing calculation setup
!> \param drho_rs real space density derivative
!> \param drho_gs g-space density derivative
!> \param i_frag index of current fragment
!> \param i_dir derivative direction (1-3 for x,y,z)
!> \param i_nuc atom index
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
! **************************************************************************************************
   SUBROUTINE fat_drho_frag(fat_env, drho_rs, drho_gs, i_frag, i_dir, i_nuc)
      TYPE(fat_env_type), POINTER                       :: fat_env
      TYPE(pw_r3d_rs_type), POINTER                     :: drho_rs
      TYPE(pw_c1d_gs_type), POINTER                     :: drho_gs
      INTEGER, INTENT(IN)                               :: i_frag
      INTEGER, INTENT(IN)                               :: i_dir
      INTEGER, INTENT(IN)                               :: i_nuc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_drho_frag'

      TYPE(qs_environment_type), POINTER                :: qs_env
!     type(dbcsr_p_type), POINTER, DIMENSION(:)         :: rho_ao
      TYPE(qs_rho_type), POINTER                        :: rho
      INTEGER                                           :: handle

      CALL timeset(routineN, handle)

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      CPASSERT(ASSOCIATED(drho_rs) .AND. ASSOCIATED(drho_gs))
      CPASSERT(i_dir >= 1 .AND. i_dir <= 3)

      qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
      CALL get_qs_env(qs_env, rho=rho)
      print *, i_dir, i_nuc
 !    CALL qs_rho_get(rho, rho_ao=rho_ao)

!     CALL calculate_drho_elec_dR( &
!           matrix_p=rho_ao(1)%matrix, drho=drho_rs, drho_gspace=drho_gs, &
!           qs_env=qs_env, beta=i_dir, lambda=i_nuc)

      CALL timestop(handle)

   END SUBROUTINE fat_drho_frag

! **************************************************************************************************
!> \brief Gets derivative of nuclear density for fragment wrt nuclear position
!> \param fat_env the freeze-and-thaw environment
!> \param drhon_gs nuclear density derivative in g-space
!> \param i_frag fragment index
!> \param beta derivative direction (1-3 for x,y,z)
!> \param lambda atom index
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
! **************************************************************************************************
   SUBROUTINE fat_drhon_frag(fat_env, drhon_gs, i_frag, i_dir, i_nuc)
      TYPE(fat_env_type), POINTER                       :: fat_env
      TYPE(pw_c1d_gs_type), POINTER                     :: drhon_gs
      INTEGER, INTENT(IN)                               :: i_frag
      INTEGER, INTENT(IN)                               :: i_dir
      INTEGER, INTENT(IN)                               :: i_nuc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fat_drhon_frag'

      TYPE(qs_environment_type), POINTER                :: qs_env
      INTEGER                                           :: handle

      CALL timeset(routineN, handle)

      CPASSERT(i_dir >= 1 .AND. i_dir <= 3)

      qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
      CALL calculate_drho_core(drhon_gs, qs_env, i_dir, i_nuc)

      CALL timestop(handle)

   END SUBROUTINE fat_drhon_frag

! **************************************************************************************************
!> \brief Calculate kinetic potential for one fragment in the freeze-and-thaw environment
!> \param fat_env the freeze-and-thaw environment
!> \param v Kinetic potential in real space (output)
!> \param v_gga GGA contribution to kinetic potential (optional output)
!> \param i_frag index of current fragment
!> \par History
!>      09.2024 created [Konstantin Komarov]
!> \author Konstantin Komarov
! **************************************************************************************************
   SUBROUTINE fat_v_ts_frag(fat_env, v_rs, v_rs_gga,  v_rs_tau, i_frag)
      TYPE(fat_env_type), POINTER                      :: fat_env
      TYPE(pw_r3d_rs_type), POINTER                    :: v_rs(:)
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER          :: v_rs_gga(:)
      TYPE(pw_r3d_rs_type), OPTIONAL, POINTER          :: v_rs_tau(:)
      INTEGER, INTENT(IN)                              :: i_frag

      CHARACTER(len=*), PARAMETER :: routineN = 'fat_v_ts_frag'

      TYPE(qs_environment_type), POINTER               :: qs_env
      TYPE(section_vals_type), POINTER                 :: input, kin_section
      TYPE(xc_rho_cflags_type)                         :: needs
      TYPE(pw_r3d_rs_type), POINTER                    :: rho_r(:)
      TYPE(pw_c1d_gs_type), POINTER                    :: rho_g(:)
      TYPE(pw_r3d_rs_type), POINTER                    :: tau(:)
      REAL(dp)                                         :: e, virial(3, 3)
      INTEGER                                          :: handle, i
      LOGICAL                                          :: use_gga = .FALSE..

      CALL timeset(routineN, handle)

      IF (is_debug) print *, "DBG| ", moduleN//"::"//routineN

      CPASSERT(ASSOCIATED(fat_env))
      CPASSERT(i_frag >= 1 .AND. i_frag <= fat_env%n_frag)

      ! Get XC section for kinetic energy functional
      qs_env => fat_env%force_env_ref%qs_env
      CALL get_qs_env(qs_env, input=input)
      kin_section => section_vals_get_subs_vals(input, "DFT%QS%OPT_EMBED%XC")
      needs = xc_functionals_get_needs( &
         section_vals_get_subs_vals(kin_section, "XC_FUNCTIONAL"), &
         lsd=.FALSE., calc_potential=.TRUE.)
      use_gga = needs%norm_drho .OR. needs%norm_drho_spin

      ! Get real space density
      NULLIFY(rho_r, rho_g, tau)
      CALL fat_rho_frag(fat_env, rho_rs=rho_r, i_frag=i_frag)
      IF (use_gga) &
         CALL fat_pw_create_frag(fat_env, pws_gs=rho_g, i_frag=i_frag)

      ! Calculate LDA kinetic potential and energy
      NULLIFY(v_rs, v_rs_gga, v_rs_tau)
      CALL xc_vxc_pw_create( &
         vxc_rho=v_rs, vxc_tau=v_rs_tau, &
         rho_r=rho_r, rho_g=rho_g, tau=tau, &
         exc=e, xc_section=kin_section, &
         pw_pool=fat_env%pw_pool_sub(i_frag)%pool, &
         compute_virial=.FALSE., virial_xc=virial)

      fat_env%ekin_sub(i_frag) = e

      ! Calculate GGA contribution if needed
      IF (use_gga) THEN
         CPASSERT(PRESENT(v_rs_gga))

         CALL xc_vxc_pw_create( &
            vxc_rho=v_rs_gga, vxc_tau=v_rs_tau, &
            rho_r=rho_r, rho_g=rho_g, tau=tau, &
            exc=e, xc_section=kin_section, &
            pw_pool=fat_env%pw_pool_sub(i_frag)%pool, &
            compute_virial=.FALSE., virial_xc=virial)

         fat_env%ekin_sub(i_frag) = fat_env%ekin_sub(i_frag) + e

         DO i = 1, SIZE(rho_g)
            CALL fat_pw_release_frag(fat_env, rho_gs(i), i_frag)
         END DO
         DEALLOCATE(rho_g)
      END IF
      CALL fat_pw_release_frag(fat_env, rho_gs(i), i_frag)

      CALL timestop(handle)

   END SUBROUTINE fat_v_ts_frag

  subroutine fat_init(fat_env, force_env, ref_subsys_number)
    type(fat_env_type), pointer   :: fat_env
    type(force_env_type), pointer :: force_env
    integer                       :: ref_subsys_number

    character(len=*), parameter :: routinen = "fat_init"

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen
    if (ASSOCIATED(fat_env) .and. ASSOCIATED(force_env)) print *, "K", ref_subsys_number

  END SUBROUTINE

  subroutine fat_finalize(fat_env)
    type(fat_env_type), pointer :: fat_env

    character(len=*), parameter :: routinen = "fat_finalize"

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    if (ASSOCIATED(fat_env)) print *, "K"

  end subroutine fat_finalize

  function fat_is_done(fat_env) result(is_done)
    type(fat_env_type), pointer :: fat_env
    logical :: is_done

    character(len=*), parameter :: routinen = "fat_is_done"

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen

    is_done = fat_env%i_iter > fat_env%n_iter

  end function
  subroutine fat_iterate_frag(fat_env, i_frag)
    type(fat_env_type), pointer :: fat_env
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_iterate_frag"
    if (is_debug) print *, "DBG| ", modulen//"::"//routinen
    if (ASSOCIATED(fat_env)) print *, "K", i_frag
    end subroutine
  subroutine fat_results(fat_env, do_forces)
    type(fat_env_type), pointer :: fat_env
    logical, optional :: do_forces

    character(len=*), parameter :: routinen = "fat_results"
    if (is_debug) print *, "DBG| ", modulen//"::"//routinen
    if (ASSOCIATED(fat_env)) print *, "K", do_forces
    end subroutine
  subroutine fat_iterate(fat_env)
    type(fat_env_type), pointer :: fat_env

    character(len=*), parameter :: routinen = "fat_iterate"

    if (is_debug) print *, "DBG| ", modulen//"::"//routinen
    if (ASSOCIATED(fat_env)) print *, "K"
    end subroutine
end module fat_methods
