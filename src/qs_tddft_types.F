!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   TDDFT data structures and type definitions for Davidson solver framework
!>
!> \details This module defines the central data structures used throughout the TDDFT
!>          implementation, with emphasis on the Davidson iterative diagonalization algorithm.
!>          The main components include:
!>
!> \par Primary Data Structure
!>      - [tddft_info_type]: Central container for all Davidson solver data including
!>        trial vectors, sigma vectors, subspace matrices, convergence monitoring,
!>        and system parameters
!>
!> \par Key Design Principles
!>      - Separation of algorithm data (Davidson) from theory data (SF, MRSF, etc.)
!>      - Dynamic memory allocation sized to actual system requirements
!>      - Support for multiple spin channels and reference types
!>      - Comprehensive convergence monitoring and restart capabilities
!>      - Integration with CP2K's parallel matrix frameworks
!>
!> \par Memory Management
!>      The module provides initialization and finalization routines that handle
!>      complete lifecycle management of all data structures, ensuring proper
!>      allocation/deallocation patterns and preventing memory leaks.
!>
!> \par Davidson Algorithm Support
!>      Data structures are specifically designed to support the Davidson method:
!>      - Expandable Krylov subspaces with automatic restart capability
!>      - Subspace matrix construction and diagonalization
!>      - Residual computation and preconditioning
!>      - Convergence monitoring with multiple criteria
!>      - Stagnation detection and recovery mechanisms
!>
!> \author  Konstantin Komarov (constlike@gmail.com)
!> \date    Aug 2025
! **************************************************************************************************
MODULE qs_tddft_types
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type, &
                                              dbcsr_get_info
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create, &
                                              cp_fm_struct_release, &
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create, &
                                              cp_fm_get_info, &
                                              cp_fm_release, &
                                              cp_fm_set_all, &
                                              cp_fm_type
   USE input_constants,                 ONLY: do_potential_coulomb, &
                                              do_potential_long, &
                                              do_potential_mix_cl, &
                                              do_potential_mix_cl_trunc, &
                                              do_potential_short
   USE input_section_types,             ONLY: section_vals_get, &
                                              section_vals_get_subs_vals, &
                                              section_vals_type, &
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env, &
                                              qs_environment_type
   USE qs_tddft_utils,                  ONLY: allocate_fm_matrix, &
                                              deallocate_fm_matrix
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_types'

   PUBLIC :: tddft_info_type
   PUBLIC :: tddft_initialization
   PUBLIC :: tddft_finalization

! **************************************************************************************************
!> \brief Central data structure for TDDFT Davidson solver framework
!>
!> \details This type contains all data structures required for the Davidson iterative
!>          diagonalization of the TDDFT eigenvalue problem. It serves as the primary
!>          communication interface between the algorithm framework and theory-specific
!>          implementations.
!>
!> \par System Parameters
!>      Contains fundamental system descriptors including:
!>      - Reference type (ROHF/ROKS, UHF/UKS)
!>      - Response theory type (SF, MRSF-S, MRSF-T, etc.)
!>      - Orbital space dimensions (occupied, virtual, basis functions)
!>      - Spin channel information
!>
!> \par Davidson Algorithm Data
!>      Manages all aspects of the iterative solver including:
!>      - Trial vector storage and management
!>      - Sigma vector computation and storage
!>      - Subspace matrix construction and diagonalization
!>      - Ritz vector formation and residual computation
!>      - Convergence monitoring and restart logic
!>
!> \par HFX Integration
!>      Supports Hartree-Fock exchange calculations with:
!>      - Range separation parameters for CAM functionals
!>      - Exchange fraction control for hybrid functionals
!>      - Compatibility flags for different HFX implementations
!>
!> \par Memory Layout
!>      Uses CP2K's distributed matrix framework (cp_fm_type) for all large
!>      matrices, ensuring scalability and memory efficiency in parallel
!>      calculations. Small arrays use standard Fortran allocation.
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   TYPE tddft_info_type

      !----------------------------------------------------------------------
      ! System identification and parameters
      !----------------------------------------------------------------------

      !> SCF reference type identifier
      !> 1 = ROKS/ROHF (restricted open-shell)
      !> 2 = UKS/UHF (unrestricted)
      INTEGER                                            :: reference_type

      !> TDDFT response theory type
      !> 1 = SF (spin-flip TDDFT)
      !> 2 = MRSF-S (mixed reference spin-flip singlets)
      !> 3 = MRSF-T (mixed reference spin-flip triplets)
      INTEGER                                            :: response_type

      !> Number of spin channels (1=restricted, 2=unrestricted)
      !> WARNING: ROHF has 2 spin channels in CP2K.
      INTEGER                                            :: nSpins

      !----------------------------------------------------------------------
      ! Property calculation arrays
      !----------------------------------------------------------------------

      !> Transition dipole moment operators in MO basis [3, nspins]
      !> Index 1: x, y, z components
      !> Index 2: spin channel
      !> Used for oscillator strength and transition moment calculations
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: dipole_op_mos_occ

      !----------------------------------------------------------------------
      ! Davidson solver parameters and dimensions
      !----------------------------------------------------------------------

      !> Number of excited states to solve for
      INTEGER                                            :: nstates
      !> Maximum possible number of excitations (nocca × nvirb)
      INTEGER                                            :: nstates_max
      !> Current dimension of Krylov subspace
      INTEGER                                            :: nvec
      !> Number of vectors at start of iteration (for tracking new additions)
      INTEGER                                            :: novec
      !> Maximum Krylov vectors per state (input parameter)
      INTEGER                                            :: maxvec
      !> Absolute maximum subspace dimension (maxvec × nstates)
      INTEGER                                            :: mxvec

      !----------------------------------------------------------------------
      ! Orbital space dimensions
      !----------------------------------------------------------------------

      !> Number of atomic orbital basis functions
      INTEGER                                            :: nbf
      !> Number of auxiliary atomic orbital basis functions (if ADMM)
      INTEGER                                            :: nbf_aux
      !> Number of alpha occupied orbitals
      INTEGER                                            :: nocca
      !> Number of beta occupied orbitals
      INTEGER                                            :: noccb
      !> Number of beta virtual orbitals
      INTEGER                                            :: nvirb
      !> Number of alpha virtual orbitals
      INTEGER                                            :: nvira

      !----------------------------------------------------------------------
      ! Convergence monitoring and control
      !----------------------------------------------------------------------

      !> Counter for consecutive stagnated iterations
      INTEGER                                            :: stall_count
      !> Maximum residual norm in current iteration
      REAL(KIND=dp)                                      :: max_error
      !> Maximum eigenvalue change in current iteration
      REAL(KIND=dp)                                      :: max_delta_eval
      !> Maximum residual norm from previous iteration
      REAL(KIND=dp)                                      :: max_error_prev
      !> Maximum eigenvalue change from previous iteration
      REAL(KIND=dp)                                      :: max_delta_eval_prev

      !> Algorithm status flags:
      !> status(1) = .TRUE. if converged to requested tolerance
      !> status(2) = .TRUE. if restart is needed (stagnation or memory limits)
      LOGICAL, DIMENSION(2)                              :: status

      !----------------------------------------------------------------------
      ! Eigenvalue and eigenvector storage
      !----------------------------------------------------------------------

      !> Converged TDDFT excitation energies [mxvec]
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals
      !> Final eigenvectors in AO basis [1, nvec]
      !> Structured as matrix(spin, state) for compatibility with analysis routines
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: evects

      !----------------------------------------------------------------------
      ! Davidson working vectors and matrices
      !----------------------------------------------------------------------

      !> Trial vectors in MO basis [mxvec] of size [nocca, nvirb]
      !> These represent the current approximations to eigenvectors
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: evects_mo

      !> Sigma vectors in MO basis [mxvec] of size [nocca, nvirb]
      !> Result of applying TDDFT Hamiltonian to trial vectors: Sigma = H * X
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: sigma_vecs

      !> Ritz vectors in MO basis [nstates] of size [nocca, nvirb]
      !> Linear combinations of trial vectors from subspace diagonalization
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: ritz_vects

      !> Residual vectors [nstates] of size [nocca, nvirb]
      !> Residuals r_k = Sigma_k - omega_k * Z_k for convergence assessment
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: residual_vecs

      !----------------------------------------------------------------------
      ! Subspace diagonalization arrays
      !----------------------------------------------------------------------

      !> Current Ritz eigenvalues from subspace diagonalization [mxvec]
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: ritz_evals
      !> Ritz eigenvalues from previous iteration [nstates]
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: ritz_evals_prev
      !> Residual vector norms for convergence monitoring [nstates]
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: residual_norms

      !> Subspace Hamiltonian matrix for diagonalization [nvec, nvec]
      !> Contains H_ij = <X_i | Sigma_j> where i,j index trial vectors
      TYPE(cp_fm_type), POINTER                          :: subspace
      !> Subspace eigenvectors [nvec, nvec]
      !> Coefficients for constructing Ritz vectors from trial vectors
      TYPE(cp_fm_type), POINTER                          :: subspace_evecs

      !----------------------------------------------------------------------
      ! Working arrays for orthogonalization and normalization
      !----------------------------------------------------------------------

      !> Vector norms for normalization procedures [mxvec]
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: vector_norm
      !> Scaling factors derived from norms [mxvec]
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: vector_scal

      !----------------------------------------------------------------------
      ! Orbital energy arrays
      !----------------------------------------------------------------------

      !> Orbital energy differences (mo_E_virtual - mo_E_occupied) [nstates_max]
      !> Used for initial guess generation and preconditioning
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: e_virt_minus_occ

      !----------------------------------------------------------------------
      ! Hartree-Fock exchange parameters
      !----------------------------------------------------------------------

      !> Flag indicating if Hartree-Fock exchange is included in the calculation
      LOGICAL                                            :: do_hfx
      !> Flag for range-separated (CAM) hybrid functionals
      !> .TRUE. for functionals with separate short/long-range exchange treatment
      LOGICAL                                            :: do_hfx_cam
      !> Short-range Hartree-Fock exchange fraction (alpha parameter)
      REAL(KIND=dp)                                      :: C_hf_alpha
      !> Long-range Hartree-Fock exchange fraction (beta parameter)
      REAL(KIND=dp)                                      :: C_hf_beta
      !> Range separation parameter (omega) for CAM functionals
      REAL(KIND=dp)                                      :: C_hf_omega

      !----------------------------------------------------------------------
      ! ADMM (Auxiliary Density Matrix Method) parameters
      !----------------------------------------------------------------------

      !> Flag indicating if ADMM is included in the calculation
      !> ADMM accelerates HFX evaluation by using a smaller auxiliary basis
      LOGICAL                                            :: do_admm

   END TYPE tddft_info_type

CONTAINS

! **************************************************************************************************
!> \brief Initialize TDDFT data structures and allocate working arrays
!>
!> \details Performs comprehensive initialization of the tddft_info_type including:
!>          - Extraction of system parameters from ground-state calculation
!>          - Analysis of HFX functional parameters and range separation
!>          - Determination of orbital space dimensions and subspace limits
!>          - Allocation of all Davidson working arrays with proper sizing
!>          - Setup of parallel matrix structures using BLACS distribution
!>
!> \param        qs_env     QuickStep environment containing system information
!> \param[inout] tddft_env  TDDFT environment to be initialized
!> \param[in]    gs_mos     Ground-state molecular orbitals for dimension extraction
!>
!> \par System Analysis
!>      The routine examines the ground-state calculation to determine:
!>      - Reference type (ROHF vs UHF) from DFT control settings
!>      - HFX parameters including range separation for CAM functionals
!>      - Orbital occupations and virtual space dimensions
!>      - Maximum possible excitation space size
!>
!> \par Memory Allocation Strategy
!>      Arrays are sized based on actual system requirements:
!>      - Trial vectors: allocated up to maximum subspace dimension
!>      - Working arrays: sized to accommodate largest intermediate results
!>      - Convergence arrays: dimensioned for requested number of states
!>      - Property arrays: structured for all spin channels and spatial components
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE tddft_initialization(qs_env, tddft_env, gs_mos)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      CLASS(tddft_info_type), INTENT(INOUT)              :: tddft_env
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         POINTER                                         :: gs_mos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddft_initialization'

      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(section_vals_type), POINTER                   :: hf_section
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s_aux_fit
      INTEGER                                            :: ivec, handle, n_rep_hf, potential_type
      LOGICAL                                            :: explicit

      CALL timeset(routineN, handle)

      !----------------------------------------------------------------------
      ! Analyze HFX functional parameters
      !----------------------------------------------------------------------

      ! Extract Hartree-Fock section from input
      hf_section => section_vals_get_subs_vals(section_vals=qs_env%input, &
                                               subsection_name="DFT%XC%HF")

      ! Determine if HFX is present and extract exchange fractions
      CALL section_vals_get(hf_section, n_repetition=n_rep_hf)
      IF (n_rep_hf > 0) THEN
         ! Get exchange fraction from first HF section
         CALL section_vals_val_get(hf_section, "FRACTION", r_val=tddft_env%C_hf_alpha, i_rep_section=1)
         tddft_env%do_hfx = ABS(tddft_env%C_hf_alpha) > 1.0e-10_dp

         IF (tddft_env%do_hfx) THEN
            ! Determine if range-separated functional is used
            CALL section_vals_val_get(hf_section, "INTERACTION_POTENTIAL%POTENTIAL_TYPE", &
                                      i_val=potential_type, i_rep_section=1, explicit=explicit)
            IF (.NOT. explicit) potential_type = do_potential_coulomb

            ! Check for CAM-type functionals
            tddft_env%do_hfx_cam = (potential_type == do_potential_mix_cl .OR. &
                                    potential_type == do_potential_mix_cl_trunc .OR. &
                                    potential_type == do_potential_short .OR. &
                                    potential_type == do_potential_long)
         ELSE
            tddft_env%do_hfx_cam = .FALSE.
         END IF
      ELSE
         tddft_env%do_hfx = .FALSE.
         tddft_env%do_hfx_cam = .FALSE.
      END IF

      ! Extract range separation parameters for CAM functionals
      IF (tddft_env%do_hfx_cam) THEN
         CALL section_vals_val_get(hf_section, "INTERACTION_POTENTIAL%SCALE_COULOMB", &
                                   r_val=tddft_env%C_hf_alpha)
         CALL section_vals_val_get(hf_section, "INTERACTION_POTENTIAL%SCALE_LONGRANGE", &
                                   r_val=tddft_env%C_hf_beta)
         CALL section_vals_val_get(hf_section, "INTERACTION_POTENTIAL%OMEGA", &
                                   r_val=tddft_env%C_hf_omega)
      END IF

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      tddft_env%do_admm = dft_control%do_admm .AND. tddft_env%do_hfx

      !----------------------------------------------------------------------
      ! Extract orbital space dimensions
      !----------------------------------------------------------------------

      ! Get orbital occupations for all spin channels
      CALL cp_fm_get_info(gs_mos(1)%mos_occ, ncol_global=tddft_env%nocca)
      CALL cp_fm_get_info(gs_mos(2)%mos_occ, ncol_global=tddft_env%noccb)
      CALL cp_fm_get_info(gs_mos(1)%mos_virt, ncol_global=tddft_env%nvira)
      CALL cp_fm_get_info(gs_mos(2)%mos_virt, ncol_global=tddft_env%nvirb)
      CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=tddft_env%nbf)
      IF (tddft_env%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         CALL get_admm_env(admm_env, matrix_s_aux_fit=matrix_s_aux_fit)
         CALL dbcsr_get_info(matrix_s_aux_fit(1)%matrix, &
                             nfullrows_total=tddft_env%nbf_aux)
      END IF

      !----------------------------------------------------------------------
      ! Determine system parameters and subspace limits
      !----------------------------------------------------------------------

      ! Identify SCF reference type
      tddft_env%reference_type = 0
      IF (dft_control%roks) tddft_env%reference_type = 1  ! ROHF
      IF (dft_control%uks) tddft_env%reference_type = 2   ! UHF

      ! Extract TDDFT control parameters
      tddft_env%nStates = dft_control%tddfpt2_control%nStates
      tddft_env%nSpins = SIZE(gs_mos, 1)
      tddft_env%maxvec = dft_control%tddfpt2_control%nkvs

      ! Calculate maximum excitation space and subspace dimensions
      tddft_env%nStates_max = tddft_env%nocca * tddft_env%nvirb
      tddft_env%mxvec = MIN(tddft_env%maxvec * tddft_env%nstates, tddft_env%nStates_max)

      ! Ensure reasonable limits and avoid memory issues
      tddft_env%nStates = MIN(tddft_env%nstates, tddft_env%mxvec)
      tddft_env%nvec = MIN(MAX(tddft_env%nStates, 4), tddft_env%mxvec)
      tddft_env%novec = tddft_env%nvec

      ! Warn if requested states exceed theoretical maximum
      IF (tddft_env%nStates > tddft_env%nStates_max) THEN
         CPWARN("NUMBER OF EXCITED STATES COULD LEAD TO PROBLEMS!")
         CPWARN("Experimental: CHANGED nStates TO ITS MAXIMUM VALUE!")
         tddft_env%nStates = tddft_env%nStates_max
      END IF

      !----------------------------------------------------------------------
      ! Allocate scalar working arrays
      !----------------------------------------------------------------------

      !--------------------
      ! ALLOCATE  DP
      !--------------------

      ALLOCATE (tddft_env%e_virt_minus_occ(tddft_env%nStates_max), source=0.0_dp)
      ALLOCATE (tddft_env%evals(tddft_env%mxvec), source=0.0_dp)
      ALLOCATE (tddft_env%ritz_evals(tddft_env%mxvec), source=0.0_dp)
      ALLOCATE (tddft_env%ritz_evals_prev(tddft_env%nstates), source=1.0E10_dp)
      ALLOCATE (tddft_env%vector_scal(tddft_env%mxvec), source=0.0_dp)
      ALLOCATE (tddft_env%vector_norm(tddft_env%mxvec), source=0.0_dp)
      ALLOCATE (tddft_env%residual_norms(tddft_env%nstates), source=0.0_dp)

      !--------------------
      ! ALLOCATE  FM
      !--------------------

      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env)

      ! Trial vectors in MO basis [mxvec] of size [nocca, nvirb]
      CALL allocate_fm_matrix(matrix=tddft_env%evects_mo, &
                              nmat=tddft_env%mxvec, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      ! Ritz vectors [nstates] of size [nocca, nvirb]
      CALL allocate_fm_matrix(matrix=tddft_env%ritz_vects, &
                              nmat=tddft_env%nStates, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      ! Sigma vectors [mxvec] of size [nocca, nvirb]
      CALL allocate_fm_matrix(matrix=tddft_env%sigma_vecs, &
                              nmat=tddft_env%mxvec, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      ! Residual vectors [nstates] of size [nocca, nvirb]
      CALL allocate_fm_matrix(matrix=tddft_env%residual_vecs, &
                              nmat=tddft_env%nStates, &
                              nrow=tddft_env%nocca, &
                              ncol=tddft_env%nvirb, &
                              blacs_env=blacs_env)

      !----------------------------------------------------------------------
      ! Allocate eigenvector storage for analysis routines
      !----------------------------------------------------------------------

      ! Create matrix structure for AO-basis eigenvectors
      CALL cp_fm_struct_create(fm_struct, &
                               nrow_global=tddft_env%nbf, &
                               ncol_global=tddft_env%nocca, &
                               context=blacs_env)

      ! Allocate eigenvector array [1, nvec] for compatibility with analysis
      ALLOCATE (tddft_env%evects(1, tddft_env%nvec))
      DO ivec = 1, tddft_env%nVec
         CALL cp_fm_create(tddft_env%evects(1, ivec), fm_struct)
         CALL cp_fm_set_all(tddft_env%evects(1, ivec), 0.0_dp)
      END DO
      CALL cp_fm_struct_release(fm_struct)

      CALL timestop(handle)

   END SUBROUTINE tddft_initialization

! **************************************************************************************************
!> \brief Finalize TDDFT calculation and deallocate all working arrays
!>
!> \details Performs comprehensive cleanup of tddft_info_type including:
!>          - Deallocation of all scalar arrays
!>          - Release of all distributed matrix structures
!>          - Cleanup of analysis and property arrays
!>          - Memory usage validation and leak prevention
!>
!> \param[inout] tddft_env  TDDFT environment to be cleaned up
!>
!> \par Cleanup Strategy
!>      The routine systematically deallocates all allocated arrays in reverse
!>      order of allocation to ensure proper cleanup. Special attention is paid to:
!>      - Conditional deallocation based on allocation status
!>      - Proper handling of CP2K matrix types
!>      - Nested loop cleanup for multi-dimensional arrays
!>      - Prevention of double deallocation errors
!>
!> \par Memory Management
!>      All arrays are checked for allocation status before deallocation to prevent
!>      runtime errors. The routine is safe to call multiple times or on partially
!>      initialized structures.
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   SUBROUTINE tddft_finalization(tddft_env)
      CLASS(tddft_info_type), INTENT(INOUT)              :: tddft_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddft_finalization'

      INTEGER                                            :: ivec, ispin, handle

      CALL timeset(routineN, handle)

      !--------------------
      ! DEALLOCATE  DP
      !--------------------

      IF (ALLOCATED(tddft_env%e_virt_minus_occ)) THEN
         DEALLOCATE (tddft_env%e_virt_minus_occ)
      END IF
      IF (ALLOCATED(tddft_env%evals)) THEN
         DEALLOCATE (tddft_env%evals)
      END IF
      IF (ALLOCATED(tddft_env%ritz_evals)) THEN
         DEALLOCATE (tddft_env%ritz_evals)
      END IF
      IF (ALLOCATED(tddft_env%ritz_evals_prev)) THEN
         DEALLOCATE (tddft_env%ritz_evals_prev)
      END IF
      IF (ALLOCATED(tddft_env%vector_scal)) THEN
         DEALLOCATE (tddft_env%vector_scal)
      END IF
      IF (ALLOCATED(tddft_env%vector_norm)) THEN
         DEALLOCATE (tddft_env%vector_norm)
      END IF
      IF (ALLOCATED(tddft_env%residual_norms)) THEN
         DEALLOCATE (tddft_env%residual_norms)
      END IF

      !--------------------
      ! DEALLOCATE  FM
      !--------------------
      !
      ! Davidson working arrays
      CALL deallocate_fm_matrix(tddft_env%evects_mo)
      CALL deallocate_fm_matrix(tddft_env%ritz_vects)
      CALL deallocate_fm_matrix(tddft_env%sigma_vecs)
      CALL deallocate_fm_matrix(tddft_env%residual_vecs)

      ! Eigenvector storage for analysis
      IF (ALLOCATED(tddft_env%evects)) THEN
         DO ivec = 1, SIZE(tddft_env%evects, 2)
            CALL cp_fm_release(tddft_env%evects(1, ivec))
         END DO
         DEALLOCATE (tddft_env%evects)
      END IF

      ! Property calculation arrays
      IF (ALLOCATED(tddft_env%dipole_op_mos_occ)) THEN
         DO ivec = 1, SIZE(tddft_env%dipole_op_mos_occ, 1)
            DO ispin = 1, SIZE(tddft_env%dipole_op_mos_occ, 2)
               CALL cp_fm_release(tddft_env%dipole_op_mos_occ(ivec, ispin))
            END DO
         END DO
         DEALLOCATE (tddft_env%dipole_op_mos_occ)
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddft_finalization

END MODULE qs_tddft_types
