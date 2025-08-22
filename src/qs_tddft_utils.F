!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddft_utils
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              int_8
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
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_p_type,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_utils'

   PUBLIC :: allocate_fm_matrix
   PUBLIC :: deallocate_fm_matrix

   INTERFACE allocate_fm_matrix
      MODULE PROCEDURE allocate_fm_matrix_single, allocate_fm_matrix_array
   END INTERFACE

   INTERFACE deallocate_fm_matrix
      MODULE PROCEDURE deallocate_fm_matrix_single, deallocate_fm_matrix_array
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Allocates and initializes a single full matrix
!> \param[out] matrix Pointer to the matrix to allocate
!> \param[in] nrow Number of rows
!> \param[in] ncol Number of columns
!> \param[in] blacs_env BLACS environment
! **************************************************************************************************
   SUBROUTINE allocate_fm_matrix_single(matrix, nrow, ncol, blacs_env)
      TYPE(cp_fm_type), POINTER, INTENT(OUT)            :: matrix
      INTEGER, INTENT(IN)                               :: nrow, ncol
      TYPE(cp_blacs_env_type), POINTER                  :: blacs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_fm_matrix_single'

      TYPE(cp_fm_struct_type), POINTER                  :: fm_struct
      INTEGER                                           :: handle

      CALL timeset(routineN, handle)

      NULLIFY(fm_struct)
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=nrow, &
                               ncol_global=ncol, &
                               context=blacs_env)

      ALLOCATE(matrix)
      CALL cp_fm_create(matrix, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(matrix, 0.0_dp)

      CALL timestop(handle)

   END SUBROUTINE allocate_fm_matrix_single

! **************************************************************************************************
!> \brief Allocates and initializes an array of full matrices
!> \param[out] matrix Allocatable array of matrices
!> \param[in] nmat Size of the matrix array
!> \param[in] nrow Number of rows per matrix
!> \param[in] ncol Number of columns per matrix
!> \param[in] blacs_env BLACS environment
! **************************************************************************************************
   SUBROUTINE allocate_fm_matrix_array(matrix, nmat, nrow, ncol, blacs_env)
      TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE, &
         INTENT(OUT)                                    :: matrix
      INTEGER, INTENT(IN)                               :: nmat, nrow, ncol
      TYPE(cp_blacs_env_type), POINTER                  :: blacs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_fm_matrix_array'

      TYPE(cp_fm_struct_type), POINTER                  :: fm_struct
      INTEGER                                           :: handle, imat

      CALL timeset(routineN, handle)

      NULLIFY(fm_struct)
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=nrow, &
                               ncol_global=ncol, &
                               context=blacs_env)

      ALLOCATE(matrix(nmat))
      DO imat = 1, nmat
         CALL cp_fm_create(matrix(imat), fm_struct)
         CALL cp_fm_set_all(matrix(imat), 0.0_dp)
      END DO

      CALL cp_fm_struct_release(fm_struct)

      CALL timestop(handle)

   END SUBROUTINE allocate_fm_matrix_array

! **************************************************************************************************
!> \brief Deallocates a single full matrix
!> \param[inout] matrix Pointer to the matrix to deallocate
! **************************************************************************************************
   SUBROUTINE deallocate_fm_matrix_single(matrix)
      TYPE(cp_fm_type), POINTER, INTENT(INOUT)          :: matrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_fm_matrix_single'

      INTEGER                                           :: handle

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(matrix)) THEN
         CALL cp_fm_release(matrix)
         DEALLOCATE(matrix)
      END IF

      CALL timestop(handle)

   END SUBROUTINE deallocate_fm_matrix_single

! **************************************************************************************************
!> \brief Deallocates an array of full matrices
!> \param[inout] matrix Allocatable array of matrices
! **************************************************************************************************
   SUBROUTINE deallocate_fm_matrix_array(matrix)
      TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: matrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_fm_matrix_array'

      INTEGER                                           :: handle, imat

      CALL timeset(routineN, handle)

      IF (ALLOCATED(matrix)) THEN
         DO imat = 1, SIZE(matrix)
            CALL cp_fm_release(matrix(imat))
         END DO
         DEALLOCATE(matrix)
      END IF

      CALL timestop(handle)

   END SUBROUTINE deallocate_fm_matrix_array

END MODULE qs_tddft_utils
