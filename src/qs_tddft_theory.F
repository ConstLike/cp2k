!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   TDDFT theory abstraction interface for excitation-state calculations
!>
!> \details This module defines the abstract base type [tddft_theory_type] for TDDFT-based
!>          excited-state solvers. The interface provides a unified contract for all theory-specific
!>          implementations and enables seamless integration of different excitation models into
!>          the Davidson iterative solver framework.
!>
!> \par Supported Theories
!>      Current implementations include:
!>      - Spin-flip TDDFT (SF): Single-reference excited states with spin contamination
!>      - Mixed Reference Spin-Flip (MRSF): Multi-reference states with pure spin separation
!>
!> \par Design Philosophy
!>      The abstraction separates algorithm-specific logic (Davidson iteration, convergence
!>      monitoring, subspace management) from theory-specific computations (operator applications,
!>      transformations, analysis). This design facilitates:
!>      - Addition of new theories without modifying core Davidson solver
!>      - Consistent interface for all excitation models
!>      - Modular testing and development of individual theories
!>      - Runtime theory selection based on input parameters
!>
!> \par Implementation Contract
!>      Each theory must implement the following deferred procedures:
!>      - allocate    : Theory-specific memory allocation and initialization
!>      - init_guess  : Construction of initial trial vectors in MO space
!>      - mo_to_ao    : Transformation from molecular orbital to atomic orbital basis
!>      - apply       : Application of TDDFT operators (HFX, Coulomb, XC kernel)
!>      - ao_to_mo    : Back-transformation from AO to MO basis
!>      - print       : Analysis and output of computed excited states
!>      - deallocate  : Memory cleanup and resource deallocation
!>
!> \par Usage Pattern
!>      The polymorphic interface enables theory-agnostic driver code:
!>      \code{.F}
!>        CLASS(tddft_theory_type), ALLOCATABLE :: theory
!>        SELECT CASE (theory_name)
!>          CASE('SF')
!>            ALLOCATE(tddft_sf_type :: theory)
!>          CASE('MRSF')
!>            ALLOCATE(tddft_mrsf_type :: theory)
!>        END SELECT
!>        CALL theory%init_guess(...)
!>        CALL theory%apply(...)
!>      \endcode
!>
!> \par Extension Guidelines
!>      To implement a new theory:
!>      1. Create a derived type extending tddft_theory_type
!>      2. Add theory-specific data members for working arrays and parameters
!>      3. Implement all deferred procedures according to the interface specification
!>      4. Add theory selection logic to the main TDDFT driver
!>      5. Provide appropriate documentation and test cases
!>
!> \author  Konstantin Komarov (constlike@gmail.com)
!> \date    Aug 2025
! **************************************************************************************************
MODULE qs_tddft_theory

   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: tddfpt2_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE qs_tddft_types,                  ONLY: tddft_info_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddft_theory'

   PUBLIC :: tddft_theory_type

! **************************************************************************************************
!> \brief Abstract base type for all TDDFT theory implementations
!>
!> \details This abstract type defines the common interface that all TDDFT theories must
!>          implement. It contains data members shared across different theories and declares
!>          deferred procedures that encapsulate theory-specific operations.
!>
!> \par Common Data Members
!>      The base type includes matrices and arrays that are used by multiple theories:
!>      - Fock matrices for orbital transformations
!>      - Hamiltonian matrices for operator applications
!>      - Working arrays for AO-MO transformations
!>      - Dipole moment operators for property calculations
!>
!> \par Theory-Specific Extensions
!>      Derived types should add their own data members for theory-specific parameters,
!>      working arrays, and intermediate results. Common additions include:
!>      - Spin-coupling parameters (MRSF theories)
!>      - Exchange correlation kernel matrices
!>      - Density difference matrices
!>      - Specialized transformation buffers
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
   TYPE, ABSTRACT :: tddft_theory_type

      ! Placeholder to prevent empty type
      INTEGER                                            :: unused = 0

      !---------------------------------------------------------------------------
      ! Common data members shared across theories
      ! Theory-specific data MUST be defined in child types (see SF, MRSF, etc.)
      !---------------------------------------------------------------------------

      !> Per-spin Fock matrices in AO basis [nspins]
      !> Used for orbital energy corrections and MO transformations
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fock

      !> Alpha occupied orbitals Fock matrix block [nocca, nocca]
      !> Contains occupied-occupied block of Fock matrix for energy corrections
      TYPE(cp_fm_type), POINTER                          :: Focka_noc_noc => NULL()

      !> Beta virtual orbitals Fock matrix block [nvirb, nvirb]
      !> Contains virtual-virtual block of Fock matrix for energy corrections
      TYPE(cp_fm_type), POINTER                          :: Fockb_vir_vir => NULL()

      !> Hamiltonian matrices in AO basis [nmat]
      !> Store results of operator applications before MO transformation
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: hmat_ao_ao

      !> Trial vectors in AO basis [nmat]
      !> Temporary storage for AO representation of excitation vectors
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: evects_ao_ao

      !----------------------------------------------------------------------
      ! HFX kernel working arrays
      !----------------------------------------------------------------------

      !> Symmetric AO density matrices for HFX calculations [nSpins]
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: work_rho_ao_ao_symm

      !> Antisymmetric AO density matrices for HFX calculations [nSpins]
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: work_rho_ao_ao_asymm

      !> Symmetric Hamiltonian matrices from HFX kernel [nSpins]
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: work_hmat_symm

      !> Antisymmetric Hamiltonian matrices from HFX kernel [nSpins]
      TYPE(dbcsr_p_type), POINTER, DIMENSION(:)          :: work_hmat_asymm

      !----------------------------------------------------------------------
      ! Property calculation arrays
      !----------------------------------------------------------------------

      !> Transition dipole moment operators in MO basis [3, nSpins]
      !> Contains x, y, z components for oscillator strength calculations
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: dipoles_mo

   CONTAINS

      !----------------------------------------------------------------------
      ! Deferred procedures - must be implemented by all derived theories
      !----------------------------------------------------------------------

      !> Allocate theory-specific working arrays and initialize parameters
      PROCEDURE(allocate_iface), DEFERRED                :: allocate
      !> Generate initial guess vectors based on orbital energy differences
      PROCEDURE(init_guess_iface), DEFERRED              :: init_guess
      !> Transform trial vectors from molecular orbital to atomic orbital basis
      PROCEDURE(mo_to_ao_iface), DEFERRED                :: mo_to_ao
      !> Apply TDDFT operators (HFX, Coulomb, XC kernel) to generate sigma vectors
      PROCEDURE(apply_op_iface), DEFERRED                :: apply
      !> Transform sigma vectors from atomic orbital back to molecular orbital basis
      PROCEDURE(ao_to_mo_iface), DEFERRED                :: ao_to_mo
      !> Analyze and print excited state properties and transition moments
      PROCEDURE(print_iface), DEFERRED                   :: print
      !> Deallocate theory-specific arrays and clean up resources
      PROCEDURE(deallocate_iface), DEFERRED              :: deallocate

   END TYPE tddft_theory_type

! **************************************************************************************************
!> \brief Interface specifications for theory-specific procedure implementations
!>
!> \details These abstract interfaces define the calling conventions and parameter lists
!>          that all theory implementations must follow. Each procedure serves a specific
!>          role in the Davidson iteration cycle and must maintain compatibility with
!>          the algorithm framework.
! **************************************************************************************************
   ABSTRACT INTERFACE

! **************************************************************************************************
!> \brief Allocate theory-specific working arrays and initialize parameters
!>
!> \details This procedure is called once at the beginning of each theory calculation
!>          to set up all required working arrays, parameters, and data structures.
!>          Common tasks include:
!>          - Allocating transformation buffers sized for the specific system
!>          - Setting theory-specific parameters from input sections
!>          - Initializing sparse matrix structures for HFX kernels
!>          - Creating parallel work distribution objects
!>
!> \param[inout] this        Theory object instance to be initialized
!> \param        qs_env      QuickStep environment containing system information
!> \param[inout] tddft_env   TDDFT environment with Davidson solver data structures
!>
!> \par Implementation Requirements
!>      - Must allocate all arrays needed for subsequent operations
!>      - Should validate system size compatibility with theory assumptions
!>      - Must initialize parallel data structures if using distributed arrays
!>      - Should set default values for all theory-specific parameters
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
      SUBROUTINE allocate_iface(this, qs_env, tddft_env)
         IMPORT                                          :: cp_blacs_env_type, &
                                                            qs_environment_type, &
                                                            tddft_info_type, &
                                                            tddft_theory_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(qs_environment_type), POINTER              :: qs_env
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env

      END SUBROUTINE

! **************************************************************************************************
!> \brief Generate initial guess vectors for Davidson iteration
!>
!> \details Constructs initial trial vectors based on ground-state orbital properties.
!>          Common approaches include:
!>          - Single excitation vectors (i to a transitions)
!>          - Vectors based on orbital energy differences
!>          - Configuration interaction singles (CIS) solutions
!>          - Previous calculation results
!>
!> \param[inout] this       Theory object instance
!> \param        qs_env     QuickStep environment
!> \param[inout] tddft_env  TDDFT environment to store initial vectors
!> \param[in]    gs_mos     Ground-state molecular orbitals for all spins
!> \param[in]    log_unit   Output unit for logging guess vector information
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
      SUBROUTINE init_guess_iface(this, qs_env, tddft_env, gs_mos, log_unit)
         IMPORT                                          :: qs_environment_type, &
                                                            tddfpt_ground_state_mos, &
                                                            tddft_info_type, &
                                                            tddft_theory_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(qs_environment_type), POINTER              :: qs_env
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
         TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
            POINTER, INTENT(IN)                          :: gs_mos
         INTEGER, INTENT(IN)                             :: log_unit

      END SUBROUTINE

! **************************************************************************************************
!> \brief Transform trial vector from molecular orbital to atomic orbital basis
!>
!> \details Performs the transformation X_AO = C_occ * X_MO * C_virt^T where C_occ and
!>          C_virt are occupied and virtual molecular orbital coefficients. The AO
!>          representation is required for applying two-electron integrals and exchange-
!>          correlation kernels.
!>
!> \param[inout] this       Theory object instance
!> \param[inout] tddft_env  TDDFT environment containing MO vectors and AO storage
!> \param[in]    gs_mos     Ground-state molecular orbitals for transformation
!> \param[in]    ivec       Index of trial vector to transform
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
      SUBROUTINE mo_to_ao_iface(this, tddft_env, gs_mos, ivec)
         IMPORT                                          :: tddfpt_ground_state_mos, &
                                                            tddft_info_type, &
                                                            tddft_theory_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
         TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
            POINTER, INTENT(IN)                          :: gs_mos
         INTEGER, INTENT(IN)                             :: ivec

      END SUBROUTINE

! **************************************************************************************************
!> \brief Apply TDDFT operators to generate sigma vectors in AO basis
!>
!> \details This is the core computational step where TDDFT operators are applied to
!>          trial vectors to generate sigma vectors: Sigma = A * X. The operators include:
!>          - Hartree-Fock exchange integrals
!>          - Coulomb electron repulsion integrals
!>          - Exchange-correlation kernel contributions
!>          - Orbital energy difference corrections
!>
!> \param[inout] this            Theory object instance
!> \param        qs_env          QuickStep environment containing integral data
!> \param[inout] tddft_env       TDDFT environment with trial vectors in AO basis
!> \param[in]    tddfpt_control  Control parameters for operator selection and accuracy
!>
!> \par Operator Applications
!>      The procedure typically performs:
!>      1. Decomposition of AO density into symmetric/antisymmetric parts
!>      2. Application of HFX integrals using existing CP2K infrastructure
!>      3. Coulomb potential evaluation via Poisson solver
!>      4. XC kernel evaluation on real-space grids
!>      5. Assembly of total sigma vector in AO basis
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
      SUBROUTINE apply_op_iface(this, qs_env, tddft_env, tddfpt_control)
         IMPORT                                          :: qs_environment_type, &
                                                            tddft_info_type, &
                                                            tddft_theory_type, &
                                                            tddfpt2_control_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(qs_environment_type), POINTER              :: qs_env
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
         TYPE(tddfpt2_control_type), POINTER, INTENT(IN) :: tddfpt_control

      END SUBROUTINE

! **************************************************************************************************
!> \brief Transform sigma vectors from atomic orbital back to molecular orbital basis
!>
!> \details Performs the back-transformation Sigma_MO = C_occ^T * Sigma_AO * C_virt to
!>          return sigma vectors to MO representation for subspace construction.
!>          May also apply orbital energy corrections and theory-specific modifications.
!>
!> \param[inout] this       Theory object instance
!> \param[inout] tddft_env  TDDFT environment containing AO sigma vectors
!> \param[in]    gs_mos     Ground-state molecular orbitals for transformation
!> \param[in]    ivec       Index of sigma vector to transform
!>
!> \par Theory-Specific Operations
!>      Beyond simple back-transformation, this procedure often includes:
!>      - Addition of diagonal orbital energy differences
!>      - Application of Fock matrix corrections for ROHF references
!>      - Enforcement of spin symmetry constraints
!>      - Scaling factors for different excitation types
!>
!> \par Output Requirements
!>      - Must store result in tddft_env%sigma_vecs(ivec)[nocca,nvirtb]
!>      - Result must be properly normalized and orthogonal to ground state
!>      - Should preserve all relevant quantum numbers (spin, symmetry)
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
      SUBROUTINE ao_to_mo_iface(this, tddft_env, gs_mos, ivec)
         IMPORT                                          :: tddfpt_ground_state_mos, &
                                                            tddft_info_type, &
                                                            tddft_theory_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
         TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
            POINTER, INTENT(IN)                          :: gs_mos
         INTEGER, INTENT(IN)                             :: ivec

      END SUBROUTINE

! **************************************************************************************************
!> \brief Analyze and print excited state properties and transition moments
!>
!> \details Performs post-convergence analysis of computed excited states including:
!>          - Excitation energy and oscillator strength calculations
!>          - Transition dipole moment evaluations
!>          - Natural transition orbital analysis
!>          - Spin expectation values and multiplicity assignments
!>          - Dominant configuration analysis
!>
!> \param[inout] this            Theory object instance
!> \param        qs_env          QuickStep environment
!> \param[inout] tddft_env       TDDFT environment containing converged states
!> \param[out]   ostrength       Computed oscillator strengths for all states
!> \param[in]    gs_mos         Ground-state molecular orbitals
!> \param[in]    tddfpt_control  Control parameters for output format and level
!> \param[in]    log_unit       Output unit for results printing
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
      SUBROUTINE print_iface(this, qs_env, tddft_env, ostrength, gs_mos, tddfpt_control, log_unit)
         IMPORT                                          :: qs_environment_type, &
                                                            tddfpt2_control_type, &
                                                            tddfpt_ground_state_mos, &
                                                            tddft_info_type, &
                                                            tddft_theory_type, &
                                                            dp
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(qs_environment_type), POINTER              :: qs_env
         TYPE(tddft_info_type), INTENT(INOUT)            :: tddft_env
         TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
            POINTER, INTENT(IN)                          :: gs_mos
         REAL(kind=dp), DIMENSION(:), INTENT(OUT)        :: ostrength
         TYPE(tddfpt2_control_type), POINTER, INTENT(IN) :: tddfpt_control
         INTEGER, INTENT(IN)                             :: log_unit

      END SUBROUTINE

! **************************************************************************************************
!> \brief Deallocate theory-specific arrays and clean up resources
!>
!> \details Performs complete cleanup of all theory-specific data structures allocated
!>          during the calculation. This includes full matrices, sparse matrices,
!>          plane wave grids, and any other dynamically allocated resources.
!>
!> \param[inout] this    Theory object instance to be cleaned up
!> \param        qs_env  QuickStep environment (for parallel cleanup coordination)
!>
!> \par Cleanup Responsibilities
!>      Must deallocate:
!>      - All full matrix arrays (cp_fm_type)
!>      - All sparse matrix arrays (dbcsr_type)
!>      - Plane wave grid structures
!>      - Parallel working arrays and communication contexts
!>      - Theory-specific parameter arrays
!>
!> \par Memory Management
!>      Proper cleanup is essential for:
!>      - Preventing memory leaks in long calculations
!>      - Allowing multiple theory runs in sequence
!>      - Clean program termination
!>      - Accurate memory usage reporting
!>
!> \author  Konstantin Komarov
!> \date    Aug 2025
! **************************************************************************************************
      SUBROUTINE deallocate_iface(this, qs_env)
         IMPORT                                          :: qs_environment_type, &
                                                            tddft_theory_type
         CLASS(tddft_theory_type), INTENT(INOUT)         :: this
         TYPE(qs_environment_type), POINTER              :: qs_env

      END SUBROUTINE

   END INTERFACE

END MODULE qs_tddft_theory
