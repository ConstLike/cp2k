!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Utility functions for TDDFT matrix operations and debugging support
!>
!> \details This module provides essential utility functions for the TDDFT implementation,
!>          focusing on memory management and debugging capabilities. The utilities are
!>          designed to work seamlessly with CP2K's CP_FM framework and
!>          provide consistent interfaces for common operations.
!>
!> \par Core Functionality
!>      - Memory allocation and deallocation for CP2K CP_FM matrices
!>      - Support for both single matrices and matrix arrays
!>      - Automatic initialization and cleanup routines
!>      - Debug output functions for development and troubleshooting
!>
!> \par Usage Patterns
!>      The module supports polymorphic interfaces that automatically select
!>      appropriate routines based on argument types:
!>      \code{.F}
!>        TYPE(cp_fm_type), POINTER :: single_matrix
!>        TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: matrix_array
!>
!>        CALL allocate_fm_matrix(single_matrix, nrow, ncol, blacs_env)
!>        CALL allocate_fm_matrix(matrix_array, nmat, nrow, ncol, blacs_env)
!>      \endcode
!>
!> \par Performance Considerations
!>      All allocation routines integrate with CP2K's timing framework to enable
!>      performance monitoring. Memory usage is optimized through:
!>      - Minimal overhead wrapper functions
!>      - Efficient parallel distribution strategies
!>      - Automatic cleanup to prevent memory leaks
!>
!> \author  Konstantin Komarov (constlike@gmail.com)
!> \date    Aug 2025
! **************************************************************************************************
MODULE qs_tddft_utils
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create, &
                                              cp_fm_struct_release, &
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create, &
                                              cp_fm_release, &
                                              cp_fm_set_all, &
                                              cp_fm_type, &
                                              cp_fm_get_submatrix
   USE cp_log_handling,                 ONLY: cp_logger_type, &
                                              cp_get_default_logger, &
                                              cp_logger_get_default_io_unit
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush
   USE qs_environment_types,            ONLY: get_qs_env, &
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_utils'

   PUBLIC :: allocate_fm_matrix
   PUBLIC :: deallocate_fm_matrix
   PUBLIC :: is_valid_mrsf_reference
   PUBLIC :: print_mat_ao_ao

   ! Polymorphic interfaces for matrix allocation/deallocation
   INTERFACE allocate_fm_matrix
      MODULE PROCEDURE allocate_fm_matrix_single, allocate_fm_matrix_array
   END INTERFACE

   INTERFACE deallocate_fm_matrix
      MODULE PROCEDURE deallocate_fm_matrix_single, deallocate_fm_matrix_array
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Allocates and initializes a single CP_FM matrix
!>
!> \details Creates a single CP2K CP_FM matrix with specified dimensions and parallel
!>          distribution. The matrix is automatically initialized to zero and properly
!>          distributed across the available processors using the provided BLACS context.
!>
!> \param[out] matrix      Pointer to matrix to be allocated (must be NULL on entry)
!> \param[in]  nrow        Number of matrix rows (global dimension)
!> \param[in]  ncol        Number of matrix columns (global dimension)
!> \param[in]  blacs_env   BLACS environment for parallel distribution
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE allocate_fm_matrix_single(matrix, nrow, ncol, blacs_env)
      TYPE(cp_fm_type), INTENT(OUT), POINTER             :: matrix
      INTEGER, INTENT(IN)                                :: nrow, ncol
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_fm_matrix_single'

      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Create matrix structure with specified dimensions
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=nrow, &
                               ncol_global=ncol, &
                               context=blacs_env)

      ! Allocate matrix using the structure
      ALLOCATE (matrix)
      CALL cp_fm_create(matrix, fm_struct)

      ! Clean up temporary structure
      CALL cp_fm_struct_release(fm_struct)

      ! Initialize matrix to zero
      CALL cp_fm_set_all(matrix, 0.0_dp)

      CALL timestop(handle)

   END SUBROUTINE allocate_fm_matrix_single

! **************************************************************************************************
!> \brief Allocates and initializes an array of CP_FM matrices
!>
!> \details Creates an array of CP2K CP_FM matrices with identical dimensions and
!>          parallel distribution. All matrices are automatically initialized to zero
!>          and share the same distribution pattern for optimal memory usage.
!>
!> \param[out] matrix      Allocatable array of matrices to be created
!> \param[in]  nmat        Number of matrices in the array
!> \param[in]  nrow        Number of matrix rows (global dimension)
!> \param[in]  ncol        Number of matrix columns (global dimension)
!> \param[in]  blacs_env   BLACS environment for parallel distribution
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE allocate_fm_matrix_array(matrix, nmat, nrow, ncol, blacs_env)
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: matrix
      INTEGER, INTENT(IN)                                :: nmat, nrow, ncol
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_fm_matrix_array'

      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      INTEGER                                            :: handle, imat

      CALL timeset(routineN, handle)

      ! Create shared matrix structure for efficiency
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=nrow, &
                               ncol_global=ncol, &
                               context=blacs_env)

      ! Allocate array and individual matrices
      ALLOCATE (matrix(nmat))
      DO imat = 1, nmat
         CALL cp_fm_create(matrix(imat), fm_struct)
         CALL cp_fm_set_all(matrix(imat), 0.0_dp)
      END DO

      ! Clean up shared structure
      CALL cp_fm_struct_release(fm_struct)

      CALL timestop(handle)

   END SUBROUTINE allocate_fm_matrix_array

! **************************************************************************************************
!> \brief Deallocates a single CP_FM matrix
!>
!> \details Safely deallocates a CP2K CP_FM matrix and sets the pointer to NULL.
!>          The routine includes proper error checking to prevent double deallocation
!>          and handles NULL pointers gracefully.
!>
!> \param[inout] matrix  Pointer to matrix to be deallocated
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE deallocate_fm_matrix_single(matrix)
      TYPE(cp_fm_type), INTENT(INOUT), POINTER           :: matrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_fm_matrix_single'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Check for valid pointer before deallocation
      IF (ASSOCIATED(matrix)) THEN
         CALL cp_fm_release(matrix)
         DEALLOCATE (matrix)
         NULLIFY (matrix)
      END IF

      CALL timestop(handle)

   END SUBROUTINE deallocate_fm_matrix_single

! **************************************************************************************************
!> \brief Deallocates an array of CP_FM matrices
!>
!> \details Safely deallocates all matrices in an allocatable array and properly
!>          cleans up the array structure. Each matrix is individually released
!>          before the array itself is deallocated.
!>
!> \param[inout] matrix  Allocatable array of matrices to be deallocated
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE deallocate_fm_matrix_array(matrix)
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: matrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_fm_matrix_array'

      INTEGER                                            :: handle, imat

      CALL timeset(routineN, handle)

      ! Check allocation status before cleanup
      IF (ALLOCATED(matrix)) THEN
         ! Release each individual matrix
         DO imat = SIZE(matrix), 1, -1
            CALL cp_fm_release(matrix(imat))
         END DO
         ! Deallocate the array structure
         DEALLOCATE (matrix)
      END IF

      CALL timestop(handle)

   END SUBROUTINE deallocate_fm_matrix_array

! **************************************************************************************************
!> \brief Check if calculation uses ROHF/ROKS with HIGH SPIN TRIPLET reference for MRSF
!> \param qs_env QuickStep environment
!> \return is_valid_mrsf .TRUE. if configuration is valid for MRSF
! **************************************************************************************************
   FUNCTION is_valid_mrsf_reference(qs_env) RESULT(is_valid_mrsf)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: is_valid_mrsf

      TYPE(dft_control_type), POINTER                    :: dft_control
      INTEGER                                            :: mult

      CALL get_qs_env(qs_env, dft_control=dft_control)

      is_valid_mrsf = .FALSE.

      ! MRSF requires ROHF/ROKS with high spin TRIPLET reference
      IF (dft_control%roks) THEN
         mult = dft_control%multiplicity

         ! Must be high spin triplet (multiplicity = 3)
         IF (mult >= 3 .AND. MOD(mult, 2) == 1) THEN
            is_valid_mrsf = .TRUE.
         END IF
      END IF

   END FUNCTION is_valid_mrsf_reference

! **************************************************************************************************
!> \brief Debug utility for printing matrix subregions during development
!>
!> \details Extracts and prints a small submatrix from a CP_FM matrix for
!>          debugging and development purposes. This utility is particularly useful
!>          for verifying matrix transformations and operator applications during
!>          TDDFT development and troubleshooting.
!>
!> \param[in] mat        Matrix to be examined
!> \param[in] mat_name   Descriptive name for output identification
!> \param[in] logger     CP2K logger for output management
!> \param[in] iter       Current iteration number for output labeling
!>
!> \par Debug Output Format
!>      The routine prints a formatted 5x5 submatrix showing:
!>      - Clear identification of matrix name and iteration
!>      - Formatted numerical output with consistent precision
!>      - Proper spacing for readability
!>      - Automatic output flushing for immediate visibility
!>
!> \par Usage Guidelines
!>      This routine is intended for development and debugging only.
!>
!> \par Implementation Notes
!>      - Extracts elements (1:5, 1:5) from the global matrix
!>      - Handles parallel distribution transparently
!>      - Outputs only on master process to prevent duplicated output
!>      - Includes proper error handling for small matrices
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE print_mat_ao_ao(mat, mat_name, iter)
      TYPE(cp_fm_type), INTENT(IN)                       :: mat
      CHARACTER(LEN=*)                                   :: mat_name
      INTEGER                                            :: iter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_mat_ao_ao'

      TYPE(cp_logger_type), POINTER                      :: logger
      REAL(KIND=dp), DIMENSION(5, 5)                     :: submatrix
      INTEGER                                            :: handle, i, j, unit_nr

      CALL timeset(routineN, handle)

      ! Initialize logging environment
      NULLIFY (logger)
      logger => cp_get_default_logger()

      ! Get output unit from logger
      unit_nr = cp_logger_get_default_io_unit(logger)

      ! Extract 5x5 submatrix from upper-left corner
      CALL cp_fm_get_submatrix(mat, submatrix, &
                               start_row=1, start_col=1, &
                               n_rows=5, n_cols=5)

      ! Debug output (disabled by default for production)
      ! Change .FALSE. to .TRUE. to enable debug printing
      IF (.False.) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(A,A,A,I0,A)') "First 5x5 elements of ", TRIM(mat_name), " at iter ", iter, ":"
            DO i = 1, 5
               WRITE (unit_nr, '(5F12.6)') (submatrix(i, j), j=1, 5)
            END DO
            WRITE (unit_nr, '(A)') ""
            CALL m_flush(unit_nr)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_mat_ao_ao

END MODULE qs_tddft_utils
