!-----------------------------------
!
! This module initially was written by Janhunen Jani Salomon
! A number of changes listed below was made by DS (me)
!
module PETSc_Solver
  use ParallelOperationValues
  use CurrentProblemValues
  use BlockAndItsBoundaries
  use petsc
  implicit none
  private !declarations are private by default

#include <petsc/finclude/petsc.h>

  KSP :: ksp     ! solver object
  PC ::  pc      ! pre-conditioner object
  
  Mat :: Amat                ! Matrix type
  Vec :: bvec                ! Right hand side
  Vec :: xvec                ! Solution
  PetscInt :: ntoteq         ! total number of equations
  PetscInt :: nloceq         ! local number of equations
  
  integer :: parComm
  
  public SolverInitialization, Solve, InsertData, FetchData, SolverDestroy, MyPetscInitialize
  
contains
  
!-----------------------------------------
!
  subroutine Solve
  implicit none
  PetscErrorCode :: ierr  

! sets KSP options from the options database
    call KSPSetOperators(ksp, Amat, Amat, ierr)
    call KSPSetFromOptions(ksp,ierr)
! solves linear system
    call KSPSolve(ksp, bvec, xvec, ierr)

  end subroutine Solve
  
!-----------------------------------------
!
  subroutine MyPetscInitialize
    implicit none  
    character(20) :: petsc_config='petsc.rc'   ! PETSc database file
    PetscErrorCode :: ierr
    ! Initializes the petsc database and mpi
    call PetscInitialize(petsc_config, ierr)
  end subroutine MyPetscInitialize
  
!-----------------------------------------
!
  subroutine SolverInitialization(comm, MatVecsCreated)

    implicit none

    INCLUDE 'mpif.h'
    
    integer, intent(IN) ::  comm  !??? not used anywhere???
    logical, intent(IN) :: MatVecsCreated

    logical :: bottom = .false.
    logical :: left =   .false.
    logical :: right =  .false.
    logical :: top =    .false. ! these designate whether the block is adjacent to respective boundaries 

    PetscInt :: one, four, five, six, nine

!    character(20) :: petsc_config='petsc.rc'   ! PETSc database file
    PetscErrorCode :: ierr

    integer jbegin, jend, ibegin, iend, solved_row
    PetscInt :: irow_global
    PetscInt ::   jcolumn_global(1:9)
    PetscScalar :: value_at_jcol(1:9)
    integer i, j

    INTEGER nio, position_flag

    REAL(8), ALLOCATABLE :: eps_ishifted_j(:,:)
    REAL(8), ALLOCATABLE :: eps_i_jshifted(:,:)

    REAL(8), ALLOCATABLE :: X11_ishifted_j(:,:), X12_ishifted_j(:,:)
    REAL(8), ALLOCATABLE :: X21_i_jshifted(:,:), X22_i_jshifted(:,:)

!    REAL(8) eps_xx, eps_xy, eps_yx, eps_yy
!    REAL(8) e_Exy, e_Wxy, e_Nyx, e_Syx
!    REAL(8) e_Wxx, e_Exx, e_Nyy, e_Syy
!    REAL(8) c_NS, c_EW  

    INTEGER ALLOC_ERR

!    integer            :: m, n, nx, ny, i, j, k, ix, jy
!    PetscInt :: nrows, ncols, one=1, five=5, temp
!    PetscInt :: irow(5), jcol(5), tmpcol(5)
!    PetscScalar,parameter  :: wv(5)=(/ 0.25, 0.25, -1.0, 0.25, 0.25 /)
!    PetscScalar :: v(5)
!    PetscViewer :: viewer
!    PetscScalar :: tol
!    INTEGER local_x_max,local_x_min,local_y_min,local_y_max
!    INTEGER i_local,j_local

    one =  1
    four = 4
    five=  5
    six =  6
    nine = 9
    
! Initializes the petsc database and mpi
!!    call PetscInitialize(petsc_config, ierr) moved to CurrentProblemValues

! Remember the communicator for future reference
    parComm=PETSC_COMM_WORLD

! Creates a matrix where the type is determined from either a call to MatSetType() 
! or from the options database with a call to MatSetFromOptions()
! The default matrix type is AIJ
    IF (.not.MatVecsCreated) THEN
       call MatCreate(parComm, Amat, ierr)

       ntoteq = N_to_solve_total
       nloceq = block_N_of_nodes_to_solve   !######## ??????????????????????????? NEW ##########

! Set MPI distribution for A matrix [by S.J.]

! Sets the local and global sizes and checks to determine compatibility
! below we use PETSC_DECIDE instead of number of local rows/columns
!####??????    call MatSetSizes(Amat, PETSC_DECIDE, PETSC_DECIDE, ntoteq, ntoteq, ierr)
       call MatSetSizes(Amat, nloceq, nloceq, ntoteq, ntoteq, ierr)    !######## ??????????????????????????? NEW ##########

! Builds matrix object for a particular matrix type
! MATAIJ = "aij" - a matrix type to be used for sparce matrices
       call MatSetType(Amat, MATAIJ, ierr)

! Creates a matrix where the type is determined from the options database
       call MatSetFromOptions(Amat, ierr)
          five=5    !?????  why???

! For good matrix assembly performance the user should preallocate the matrix storage
! the second argument is the number of nonzeros per row (same for all rows)
       call MatSeqAIJSetPreallocation(Amat, nine, PETSC_NULL_INTEGER, ierr)
       five=5    !????? why???

! Preallocates memory for a sparce parallel matrix in AIJ format (the default parallel petsc format)
! the second argument is the number of nonzeros per row in DIAGONAL portion of local submatrix
! the fourth argument is the number of nonzeros per row in the OFF-DIAGONAL portion of local submatrix
       call MatMPIAIJSetPreallocation(Amat, nine, PETSC_NULL_INTEGER, nine, PETSC_NULL_INTEGER, ierr)
       one=1    !????  why???
    ELSE
    END IF

! ny - the number of columns
! jcol - global indices of columns
! v - a logically two-dimensional array of values

!----------------------------------------
! Initialization of matrix coefficients written by DS
! globally indexed matrix elements are initialized in parallel

    jbegin = indx_y_min+1
    jend   = indx_y_max-1
    ibegin = indx_x_min+1
    iend   = indx_x_max-1

    IF (Rank_of_process_left.LT.0)  THEN 
        ibegin = indx_x_min
        left = .true.
    END IF    
    IF (Rank_of_process_right.LT.0) THEN
        iend   = indx_x_max
        right = .true.
    END IF    
    IF (Rank_of_process_below.LT.0) THEN
       jbegin = indx_y_min
       bottom = .true.
    END IF 
    IF (Rank_of_process_above.LT.0) THEN
        jend   = indx_y_max
        top = .true.
    END IF   

    solved_row = iend - ibegin + 1

    ALLOCATE(eps_ishifted_j(ibegin:iend+1, jbegin:jend), STAT=ALLOC_ERR)   ! eps_ishifted_j(i,j) is between nodes {i-1,j} and {i,j}
    ALLOCATE(eps_i_jshifted(ibegin:iend, jbegin:jend+1), STAT=ALLOC_ERR)   ! eps_i_jshifted(i,j) is between nodes {i,j-1} and {i,j}

    eps_ishifted_j = 1.0_8
    eps_i_jshifted = 1.0_8

    DO j = jbegin, jend
       DO i = ibegin, iend+1
          CALL SET_EPS_ISHIFTED(i, j, eps_ishifted_j(i,j))
       END DO
    END DO

    DO j = jbegin, jend+1
       DO i = ibegin, iend
          CALL SET_EPS_JSHIFTED(i, j, eps_i_jshifted(i,j))
       END DO
    END DO

    ALLOCATE(X11_ishifted_j(ibegin:iend+1, jbegin:jend), STAT=ALLOC_ERR)   ! eps_ishifted_j(i,j) is between nodes {i-1,j} and {i,j}
    ALLOCATE(X12_ishifted_j(ibegin:iend+1, jbegin:jend), STAT=ALLOC_ERR)   ! eps_ishifted_j(i,j) is between nodes {i-1,j} and {i,j}

    ALLOCATE(X21_i_jshifted(ibegin:iend, jbegin:jend+1), STAT=ALLOC_ERR)   ! eps_i_jshifted(i,j) is between nodes {i,j-1} and {i,j}
    ALLOCATE(X22_i_jshifted(ibegin:iend, jbegin:jend+1), STAT=ALLOC_ERR)   ! eps_i_jshifted(i,j) is between nodes {i,j-1} and {i,j}

    X11_ishifted_j = 0.0_8
    X12_ishifted_j = 0.0_8
    X21_i_jshifted = 0.0_8
    X22_i_jshifted = 0.0_8

    DO j = jbegin, jend
       DO i = ibegin, iend + 1
          CALL SET_X11_X12_ISHIFTED(i, j, X11_ishifted_j(i,j), X12_ishifted_j(i,j))
       END DO
    END DO

    DO j = jbegin, jend + 1
       DO i = ibegin, iend
          CALL SET_X21_X22_JSHIFTED(i, j, X21_i_jshifted(i,j), X22_i_jshifted(i,j))
       END DO
    END DO

    irow_global = global_offset ! first global node number in the block will be (global_offset + 1)

    j = indx_y_min !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    IF (bottom) THEN
! boundary object along bottom border
       DO i = ibegin, iend !nodes on the border y = 0
          irow_global = irow_global + 1
          jcolumn_global(1) = nodenum(i,j)        ! C
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
       END DO
    END IF

    DO j = indx_y_min + 1, indx_y_max - 1 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                       
       i = indx_x_min

       IF (left) THEN
! boundary object along left border
          irow_global = irow_global + 1
          IF (.NOT.block_has_symmetry_plane_X_left) THEN
! Dirichlet (given potential) boundary
             jcolumn_global(1) = nodenum(i,j)        ! C
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
          ELSE
! the left border is a symmetry plane

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(i, j, nio, position_flag)

! note that we are at the left edge of the domain and the symmetry is applied here,
! so the boundary object must be symmetric relative to x=0 as well
! therefore we are either at the bottom surface, or inside, or at the top surface of the inner object

             SELECT CASE (position_flag)
                CASE (9)
! metal
                   jcolumn_global(1) = nodenum(i,j)        ! C
                   value_at_jcol(1) = 1.0_8
                   call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
                CASE DEFAULT 
! inside dielectric or plasma
! assume *_ishifted_j(i,j) = *_ishifted_j(i+1,j)

                   jcolumn_global(1) = nodenum(i,j-1)      ! S
                   jcolumn_global(2) = nodenum(i,j)        ! C
                   jcolumn_global(3) = nodenum(i+1,j)      ! E
                   jcolumn_global(4) = nodenum(i,j+1)      ! N
          
                   jcolumn_global(5) = nodenum(i+1,j+1)    ! NE, for center node on symmetry line
                   jcolumn_global(6) = nodenum(i+1,j-1)    ! SE
! base stencil
                   value_at_jcol(1) =   eps_i_jshifted(i,j)                                                                               ! S
                   value_at_jcol(2) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i+1,j) + eps_ishifted_j(i+1,j))      ! C
                   value_at_jcol(3) =   eps_ishifted_j(i+1,j) + eps_ishifted_j(i+1,j)                                                     ! E
                   value_at_jcol(4) =   eps_i_jshifted(i,j+1)                                                                             ! N
! implicit susceptibility contribution
                   value_at_jcol(1) = value_at_jcol(1) + X22_i_jshifted(i,j) - X12_ishifted_j(i+1,j) + X12_ishifted_j(i,j)                             ! S
                   value_at_jcol(2) = value_at_jcol(2) - X11_ishifted_j(i,j) - X11_ishifted_j(i+1,j) - X22_i_jshifted(i,j) - X22_i_jshifted(i,j+1)     ! C

                   value_at_jcol(3) = value_at_jcol(3) + X11_ishifted_j(i+1,j) + X21_i_jshifted(i,j+1) - X21_i_jshifted(i,j) + &                       ! E
                                                       & X11_ishifted_j(i,j)   - X21_i_jshifted(i,j+1) + X21_i_jshifted(i,j)

                   value_at_jcol(4) = value_at_jcol(4) + X22_i_jshifted(i,j+1) + X12_ishifted_j(i+1,j) - X12_ishifted_j(i,j)                           ! N

                   value_at_jcol(5) = X12_ishifted_j(i+1,j) + X21_i_jshifted(i,j+1) - &                                                                ! NE
                                    & X12_ishifted_j(i,j)   - X21_i_jshifted(i,j+1)

                   value_at_jcol(6) = -X12_ishifted_j(i+1,j) - X21_i_jshifted(i,j) + &                                                                 ! SE
                                    &  X12_ishifted_j(i,j)   + X21_i_jshifted(i,j)

                   call MatSetValues(Amat, one, irow_global, six, jcolumn_global(1:6), value_at_jcol(1:6), INSERT_VALUES, ierr) 
             END SELECT

          END IF   !### IF (.NOT.block_has_symmetry_plane_X_left) THEN
       END IF      !### IF (left) THEN

       DO i = indx_x_min + 1, indx_x_max - 1

          irow_global = irow_global + 1

          IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN
             jcolumn_global(1) = nodenum(i,j)        ! C
             value_at_jcol(1) = 1.0_8
             call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
          ELSE

! check whether the point is inside or at the surface of any inner object
             CALL FIND_INNER_OBJECT_CONTAINING_POINT(i, j, nio, position_flag)

             SELECT CASE (position_flag)
                CASE (9)
! metal
                   jcolumn_global(1) = nodenum(i,j)      ! C
                   value_at_jcol(1) = 1.0_8
                   call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr)
                CASE DEFAULT 
! inside dielectric, inside plasma, or in a corner of a dielectric object

                   jcolumn_global(1) = nodenum(i  ,j-1)      ! S
                   jcolumn_global(2) = nodenum(i-1,j  )      ! W
                   jcolumn_global(3) = nodenum(i  ,j  )      ! C
                   jcolumn_global(4) = nodenum(i+1,j  )      ! E
                   jcolumn_global(5) = nodenum(i  ,j+1)      ! N
          
                   jcolumn_global(6) = nodenum(i+1,j+1)      ! NE
                   jcolumn_global(7) = nodenum(i+1,j-1)      ! SE
                   jcolumn_global(8) = nodenum(i-1,j-1)      ! SW
                   jcolumn_global(9) = nodenum(i-1,j+1)      ! SE
! base stencil
                   value_at_jcol(1) =   eps_i_jshifted(i,j)                                                                                  ! S
                   value_at_jcol(2) =   eps_ishifted_j(i,j)                                                                                  ! W
                   value_at_jcol(3) = -(eps_i_jshifted(i,j) + eps_i_jshifted(i,j+1) + eps_ishifted_j(i,j) + eps_ishifted_j(i+1,j))           ! C
                   value_at_jcol(4) =   eps_ishifted_j(i+1,j)                                                                                ! E
                   value_at_jcol(5) =   eps_i_jshifted(i,j+1)                                                                                ! N
! implicit susceptibility contribution
                   value_at_jcol(1) = value_at_jcol(1) + X22_i_jshifted(i,j) - X12_ishifted_j(i+1,j) + X12_ishifted_j(i,j)                             ! S
                   value_at_jcol(2) = value_at_jcol(2) + X11_ishifted_j(i,j) - X21_i_jshifted(i,j+1) + X21_i_jshifted(i,j)                             ! W
                   value_at_jcol(3) = value_at_jcol(3) - X11_ishifted_j(i,j) - X11_ishifted_j(i+1,j) - X22_i_jshifted(i,j) - X22_i_jshifted(i,j+1)     ! C
                   value_at_jcol(4) = value_at_jcol(4) + X11_ishifted_j(i+1,j) + X21_i_jshifted(i,j+1) - X21_i_jshifted(i,j)                           ! E
                   value_at_jcol(5) = value_at_jcol(5) + X22_i_jshifted(i,j+1) + X12_ishifted_j(i+1,j) - X12_ishifted_j(i,j)                           ! N

                   value_at_jcol(6) =  X12_ishifted_j(i+1,j) + X21_i_jshifted(i,j+1)                                                            ! NE
                   value_at_jcol(7) = -X12_ishifted_j(i+1,j) - X21_i_jshifted(i,j)                                                              ! SE
                   value_at_jcol(8) =  X12_ishifted_j(i,j)   + X21_i_jshifted(i,j)                                                              ! SW
                   value_at_jcol(9) = -X12_ishifted_j(i,j)   - X21_i_jshifted(i,j+1)                                                            ! SE

                   call MatSetValues(Amat, one, irow_global, nine, jcolumn_global(1:9), value_at_jcol(1:9), INSERT_VALUES, ierr) 
             END SELECT
          END IF   !### IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) THEN

       END DO   !###   DO i = indx_x_min + 1, indx_x_max - 1

       i = indx_x_max

       IF (iend.EQ.indx_x_max) THEN
! boundary object along right border
          irow_global = irow_global + 1
          jcolumn_global(1) = nodenum(i,j)     ! C
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
       END IF

    END DO   !###   DO j = indx_y_min + 1, indx_y_max - 1

    j = indx_y_max !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.

    IF (jend.EQ.indx_y_max) THEN
! boundary object along top border
       DO i = ibegin, iend
          irow_global = irow_global + 1
          jcolumn_global(1) = nodenum(i,j)     ! C
          value_at_jcol(1) = 1.0_8
          call MatSetValues(Amat, one, irow_global, one, jcolumn_global(1:1), value_at_jcol(1:1), INSERT_VALUES, ierr) 
       END DO
    END IF

! end of initialization of matrix coefficients written by DS ---------------

    DEALLOCATE(eps_ishifted_j, STAT=ALLOC_ERR)   ! eps_ishifted_j(i,j) is between nodes {i-1,j} and {i,j}
    DEALLOCATE(eps_i_jshifted, STAT=ALLOC_ERR)   ! eps_i_jshifted(i,j) is between nodes {i,j-1} and {i,j}

    DEALLOCATE(X11_ishifted_j, STAT=ALLOC_ERR)
    DEALLOCATE(X12_ishifted_j, STAT=ALLOC_ERR)

    DEALLOCATE(X21_i_jshifted, STAT=ALLOC_ERR)
    DEALLOCATE(X22_i_jshifted, STAT=ALLOC_ERR)

! call MatSetValues(Amat, one, ix, ny, jcol, v, INSERT_VALUES,ierr) 
! inserts or adds a block of values into a matrix
! MatAssemblyBegin() and MatAssemblyEnd() must be called after all calls to matSetValues have been completed.
! Amat - the matrix
! one - the number of rows
! ix - global indices of rows
! ny - the number of columns
! jcol - global indices of columns
! v - a logically two-dimensional array of values
! INSERT_VALUES replaces existing entries with new values

! Commit changes to the matrix [S.J.]

! Begins assembling the matrix
    call MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY, ierr)
! Completes assembling the matrix
    call MatAssemblyEnd(Amat, MAT_FINAL_ASSEMBLY, ierr)
    
! Save matrix for testing ....
!    if(.true.) then
!       call  PetscViewerBinaryOpen(parComm, 'Amat.bin', FILE_MODE_WRITE, viewer, ierr)
!       call  MatView(Amat, viewer, ierr)
!       call  PetscViewerDestroy(viewer,ierr)        
!    end if
    
    IF (.not.MatVecsCreated) THEN
! Create solver [S.J.]

! KSP is abstract petsc object that manages all Krylov methods. 
! It manages the linear solves (even when krylov accelerators are not used like in direct solvers

! Creates the default KSP context
! here ksp is location where to put KSP context
       call KSPCreate(parComm, ksp, ierr)

! Sets the matrix associated with the linear system and a (possibly) different one associated with the preconditioner
! argument #3 is the matrix to be used in constructing the preconditioner, 
! usually the same as argument #2 (the matrix that defines the linear system)
       call KSPSetOperators(ksp, Amat, Amat, ierr)

! . . Set method (this set-up uses the LU-decomposition as the solver) [S.J.]

! returns a pointer to the preconditioner context set with KSPSetPc
       call KSPGetPc(ksp, pc, ierr)

! sets the preconditioner to be used to calculate the application of the preconditioner on a vector
       call KSPSetUp(ksp, ierr)

! Create B and X vectors [S.J.]

! note that ### bvec ### is the vector of the right hand side and ### xvec ### is the solution vector

! get vector(s) compatible with the matrix, i.e. with the same parallel layout
       call MatCreateVecs(Amat, PETSC_NULL_VEC, bvec, ierr)

! configures the vector from the options database
       call VecSetFromOptions(bvec, ierr)

! creates a new vector of the same type as an existing vector
       call VecDuplicate(bvec, xvec, ierr)
    ELSE
       CALL KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)        
    END IF        
    
    ! Done with setting up vectors and matrices!
    return
  end subroutine SolverInitialization
  
!-------------------------------------
!
! this subroutine inserts dvec which in the actual call is the real(8) array rhs
! into vector bvec
  subroutine InsertData(nnodes,dvec)

    PetscInt, intent(IN) :: nnodes
    PetscScalar, intent(IN) :: dvec(:)

    PetscInt, save :: n, val
    PetscInt, allocatable, save :: jvec(:)

    integer ibegin, jbegin, iend, jend

    integer :: i, j !, k
!    PetscViewer :: viewer
!    character(20) :: filename='bvector.dat'
    PetscErrorCode :: ierr
    
    if(.not.allocated(jvec)) then

       allocate(jvec(1 : nnodes))

! calculation of jvec modified by DS

       jbegin = indx_y_min+1
       jend   = indx_y_max-1
       ibegin = indx_x_min+1
       iend   = indx_x_max-1

       IF (Rank_of_process_left.LT.0)  ibegin = indx_x_min
       IF (Rank_of_process_right.LT.0) iend   = indx_x_max
       IF (Rank_of_process_below.LT.0) jbegin = indx_y_min
       IF (Rank_of_process_above.LT.0) jend   = indx_y_max

       n=0
       do j = jbegin, jend  !indx_y_min+1, indx_y_max-1
          do i = ibegin, iend  !indx_x_min+1, indx_x_max-1
!          if(i<1.or.i>=global_maximal_i.or.j<1.or.j>=global_maximal_j) cycle
             n = n+1
             jvec(n)=global_offset + n  !(N_grid_block_x*N_grid_block_y)*Rank_of_process+N_grid_block_x*(j-indx_y_min-1)+i-indx_x_min-1 
          end do
       end do

    end if

    val = n 
    if(.not.(nnodes.eq.n)) print *,'Not conformant!'
    if(size(dvec)<n) print *,'dvec too small!', size(dvec), n

! inserts or adds values into certain locations of a vector
! bvec - vector to insert in
! val - number of elements to add
! jvec - indices where to add
! dvec - array of values
! INSERT_VALUES - replaces existing entries with new values
    call VecSetValues(bvec, val, jvec, dvec, INSERT_VALUES, ierr)
    
! begins assembling the vector, should be called after completing all calls to VecSetValues
    call VecAssemblyBegin(bvec, ierr)

! completes assembling the vector
    call VecAssemblyEnd(bvec, ierr)

!    call PetscViewerCreate(parComm,viewer,ierr)
!    call PetscViewerASCIIOpen(parComm,filename,viewer,ierr)
!!    call PetscViewerSetType(viewer,PETSCVIEWERMATLAB,ierr)
!    call VecView(bvec,viewer,ierr)
!    call PetscViewerDestroy(parComm,viewer,ierr)
    return
  end subroutine InsertData

!-------------------------------------
!
  subroutine FetchData(nnodes, dvec)

    PetscInt, intent(IN) :: nnodes
    PetscScalar, intent(OUT) :: dvec(:)

    integer ibegin, jbegin, iend, jend

    integer :: n, i, j   !, k, val
    PetscInt, allocatable, save :: jvec(:)
!    PetscViewer :: viewer
!    character(20) :: filename='soln.mat'
    PetscErrorCode :: ierr
    
    if(.not.allocated(jvec)) then
       allocate(jvec(1 : nnodes))

! calculation of jvec modified by DS

       jbegin = indx_y_min+1
       jend   = indx_y_max-1
       ibegin = indx_x_min+1
       iend   = indx_x_max-1

       IF (Rank_of_process_left.LT.0)  ibegin = indx_x_min
       IF (Rank_of_process_right.LT.0) iend   = indx_x_max
       IF (Rank_of_process_below.LT.0) jbegin = indx_y_min
       IF (Rank_of_process_above.LT.0) jend   = indx_y_max

       n=0
       do j = jbegin, jend  !indx_y_min+1, indx_y_max-1
          do i = ibegin, iend  !indx_x_min+1, indx_x_max-1
!            if(i<1.or.i>=global_maximal_i.or.j<1.or.j>=global_maximal_j) cycle
             n = n+1
             jvec(n)=global_offset + n  !(N_grid_block_x*N_grid_block_y)*Rank_of_process+N_grid_block_x*(j-indx_y_min-1)+i-indx_x_min-1
          end do
       end do
    end if

! gets values from certain locations of a vector
    call VecGetValues(xvec, nnodes, jvec, dvec, ierr)
    
! begins assembling the vector
    call VecAssemblyBegin(xvec, ierr)
! completes assembling the vector
    call VecAssemblyEnd(xvec, ierr)

!  Get the local part of the solution
!    call PetscViewerCreate(parComm,viewer,ierr)
!    call PetscViewerASCIIOpen(parComm,filename,viewer,ierr)
!!    call PetscViewerSetType(viewer,PETSCVIEWERMATLAB,ierr)
!    call VecView(bvec,viewer,ierr)
!    call PetscViewerDestroy(parComm,viewer,ierr)
    return
  end subroutine FetchData
  

!-------------------------------------
!
  subroutine SolverDestroy
    PetscErrorCode :: ierr
    call VecDestroy(xvec, ierr)
    call VecDestroy(bvec, ierr)
    call MatDestroy(Amat,ierr)
    call KSPDestroy(ksp,ierr)
    return
  end subroutine SolverDestroy
  
end module PETSc_Solver

!-------------------------------
! 7---6---5
! |       |
! 8   0   4
! |       |
! 1---2---3
!
SUBROUTINE FIND_INNER_OBJECT_CONTAINING_POINT(i, j, nio, position_flag)

  USE CurrentProblemValues

  IMPLICIT NONE

  INTEGER, INTENT(IN) ::  i, j  ! x,y indices of the point of interest
  INTEGER, INTENT(OUT) :: nio   ! number of inner object containing the point
                                ! it is meaningful only if position_flag>=0
  INTEGER, INTENT(OUT) :: position_flag   ! 9 for a metal object
                                          ! 0 inside a dielectric object
                                          ! 1 to 8 at the surface of the dielectric object
                                          ! -1 if there is no object containing the given point
  INTEGER n

  nio = -1
  position_flag = -1

  IF (N_of_inner_objects.EQ.0) RETURN

  DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
       IF (i.LT.whole_object(n)%ileft) CYCLE
       IF (i.GT.whole_object(n)%iright) CYCLE
       IF (j.LT.whole_object(n)%jbottom) CYCLE
       IF (j.GT.whole_object(n)%jtop) CYCLE

! since we are here, the point is either at the surface or inside inner material object n

       nio = n  ! save the number of the inner object

       IF (whole_object(n)%object_type.EQ.METAL_WALL) THEN
! the material object is metal, stop scanning inner objects, dielectric properties are omitted (relevant if metal is inside dielectric)
          position_flag = 9
          RETURN
       END IF

       position_flag = 0
       IF (i.EQ.whole_object(n)%ileft) THEN
! point on the left side of a dielectric object
          IF (j.EQ.whole_object(n)%jbottom) THEN
! left bottom corner
             position_flag = 1
          ELSE IF (j.EQ.whole_object(n)%jtop) THEN
! left top corner
             position_flag = 7
          ELSE
! left side, not a corner
             position_flag = 8
          END IF
       ELSE IF (i.EQ.whole_object(n)%iright) THEN
! point on the right side of a dielectric object
          IF (j.EQ.whole_object(n)%jbottom) THEN
! right bottom corner
             position_flag = 3
          ELSE IF (j.EQ.whole_object(n)%jtop) THEN
! right top corner
             position_flag = 5
          ELSE
! right side, not a corner
             position_flag = 4
          END IF
       ELSE IF (j.EQ.whole_object(n)%jbottom) THEN
! point on the bottom side of a dielectric object
          position_flag = 2
       ELSE IF (j.EQ.whole_object(n)%jtop) THEN
! point on the top side of a dielectric object
          position_flag = 6
       END IF
! note that we do not stop scanning inner material objects here since there may be metal objects inside the dielectric
    END DO

    RETURN

END SUBROUTINE FIND_INNER_OBJECT_CONTAINING_POINT

!-------------------------------
! 7---6---5
! |       |
! 8   0   4
! |       |
! 1---2---3
!
SUBROUTINE CHECK_IF_INNER_OBJECT_CONTAINS_POINT(myobject, i, j, position_flag)

  USE CurrentProblemValues

  IMPLICIT NONE

  TYPE(boundary_object), INTENT(IN) :: myobject
  INTEGER, INTENT(IN) :: i, j  ! x,y indices of the point of interest
  INTEGER, INTENT(OUT) :: position_flag   ! 9 for a metal object
                                          ! 0 inside a dielectric object
                                          ! 1 to 8 at the surface of the dielectric object
                                          ! -1 if there is no object containing the given point
  INTEGER n

  position_flag = -1

!  IF (N_of_inner_objects.EQ.0) RETURN

  IF (i.LT.myobject%ileft) RETURN
  IF (i.GT.myobject%iright) RETURN
  IF (j.LT.myobject%jbottom) RETURN
  IF (j.GT.myobject%jtop) RETURN

! since we are here, the point is either at the surface or inside inner material object n

  IF (myobject%object_type.EQ.METAL_WALL) THEN
! the material object is metal, stop scanning inner objects, dielectric properties are omitted (relevant if metal is inside dielectric)
     position_flag = 9
     RETURN
  END IF

  position_flag = 0
  IF (i.EQ.myobject%ileft) THEN
! point on the left side of a dielectric object
     IF (j.EQ.myobject%jbottom) THEN
! left bottom corner
        position_flag = 1
     ELSE IF (j.EQ.myobject%jtop) THEN
! left top corner
        position_flag = 7
     ELSE
! left side, not a corner
        position_flag = 8
     END IF
  ELSE IF (i.EQ.myobject%iright) THEN
! point on the right side of a dielectric object
     IF (j.EQ.myobject%jbottom) THEN
! right bottom corner
        position_flag = 3
     ELSE IF (j.EQ.myobject%jtop) THEN
! right top corner
        position_flag = 5
     ELSE
! right side, not a corner
        position_flag = 4
     END IF
  ELSE IF (j.EQ.myobject%jbottom) THEN
! point on the bottom side of a dielectric object
     position_flag = 2
  ELSE IF (j.EQ.myobject%jtop) THEN
! point on the top side of a dielectric object
     position_flag = 6
  END IF
! note that we do not stop scanning inner material objects here since there may be metal objects inside the dielectric

  RETURN

END SUBROUTINE CHECK_IF_INNER_OBJECT_CONTAINS_POINT

!--------------------------------------------------
!
REAL(8) FUNCTION Get_Surface_Charge_Inner_Object(i,j,position_flag, myobject)

  USE CurrentProblemValues
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
 
  INTEGER i, j            ! x,y indices of the point of interest
  INTEGER position_flag   ! 9 for a metal object
                          ! 0 inside a dielectric object
                          ! 1 to 8 at the surface of the dielectric object
                          ! -1 if there is no object containing the given point
  TYPE(boundary_object) myobject

  INTEGER pos

  INTEGER i_left_top, i_right_top, i_right_bottom

  Get_Surface_Charge_Inner_Object = 0.0_8 

  IF (N_of_inner_objects.EQ.0) RETURN
!  IF (nio.LT.N_of_boundary_objects+1) RETURN
!  IF (nio.GT.N_of_boundary_and_inner_objects) RETURN
  IF (position_flag.LT.1) RETURN
  IF (position_flag.GT.8) RETURN

  pos=-1

! ilt   -  ---- irt
!  |             |
!  |             |
!  1 ilbb ---- irb
!
  i_left_top        =                  myobject%jtop   - myobject%jbottom + 1
  i_right_top       = i_left_top     + myobject%iright - myobject%ileft
  i_right_bottom    = i_right_top    + myobject%jtop   - myobject%jbottom
!  i_left_bottom_bis = i_right_bottom + myobject%iright - myobject%ileft - 1

  SELECT CASE (position_flag)
     CASE (1)
! left bottom corner
        pos=1
     CASE (8)
! left side
        pos = j - myobject%jbottom  + 1
     CASE (7)
! left top corner
        pos = i_left_top
     CASE (6)
! top side
        pos = i_left_top + i - myobject%ileft
     CASE (5) 
! right top corner
        pos = i_right_top
     CASE (4)
! right side
        pos = i_right_top + myobject%jtop - j
     CASE (3)
! right bottom corner
        pos = i_right_bottom
     CASE (2)
  ! bottom side
        pos = i_right_bottom + myobject%iright - i
  END SELECT

  IF ((pos.GE.1).AND.(pos.LE.myobject%N_boundary_nodes)) THEN
     Get_Surface_Charge_Inner_Object = myobject%surface_charge(pos)
     RETURN
  ELSE
     print '("Error in Get_Surface_Charge_Inner_Object")'
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

END FUNCTION Get_Surface_Charge_Inner_Object

!-----------------------------------------------
!
SUBROUTINE SET_EPS_ISHIFTED(i, j, eps)  ! here point {i,j} is between nodes {i-1,j} and {i,j}

  USE CurrentProblemValues

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
 
  INTEGER, INTENT(IN) :: i, j
  REAL(8), INTENT(OUT) :: eps

  INTEGER count  ! counts dielectric objects
  LOGICAL segment_is_inside_dielectric_object
  LOGICAL segment_is_on_surface_of_metal_object
  INTEGER n1

! find all inner objects owning segment {i-1,j}-{i,j}
! assume that only two dielectric objects may own a common segment
! allow a metal object to be added on top of that

  eps = 0.0_8
  count = 0
  segment_is_inside_dielectric_object = .FALSE.

! first, look through the dielectric objects only
  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.EQ.METAL_WALL) CYCLE

! find an inner object containing point {i-1,j}

     IF (i-1.LT.whole_object(n1)%ileft) CYCLE
     IF (i-1.GT.whole_object(n1)%iright) CYCLE
     IF (j.LT.whole_object(n1)%jbottom) CYCLE
     IF (j.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i-1,j} is either at the surface or inside inner object n1

! find whether this object contains point {i,j}
 
     IF (i.GT.whole_object(n1)%iright) CYCLE

! since we are here, point {i,j} is either at the surface or inside inner object n1
! and segment {i-1,j}-{i,j} is owned by the object
     
! find whether segment {i-1,j}-{i,j} is at the surface or inside
! it only can be on the surface if this is the top or the bottom edge of the object

     IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is on the surface

        count = count + 1
        eps = eps + whole_object(n1)%eps_diel

     ELSE   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is inside
        segment_is_inside_dielectric_object = .TRUE.
        eps = whole_object(n1)%eps_diel
        count = 1
        EXIT

     END IF   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! at this point
! if count==0 then the segment is either in/at metal object or in vacuum
! if count==1 then the segment either inside a dielectric object or on the surface of ONE (not two) dielectric object
! if count==2 then the segment is at the interface between two dielectrics
! if count==3 this is an error

! second, look through the metal objects only
! if segment belongs to a metal object 
! return zero if it is inside
! if it is on the surface
! return (#1) epsilon of adjacent dielectric object or (#2) 1 if the metal object is in vacuum or (#3) 0 if the metal object is attached to another metal object

  segment_is_on_surface_of_metal_object = .FALSE.

  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.NE.METAL_WALL) CYCLE

! find an inner object containing point {i-1,j}

     IF (i-1.LT.whole_object(n1)%ileft) CYCLE
     IF (i-1.GT.whole_object(n1)%iright) CYCLE
     IF (j.LT.whole_object(n1)%jbottom) CYCLE
     IF (j.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i-1,j} is either at the surface or inside inner object n1

! find whether this object contains point {i,j}
 
     IF (i.GT.whole_object(n1)%iright) CYCLE

! since we are here, point {i,j} is either at the surface or inside inner object n1
! and segment {i-1,j}-{i,j} is owned by the object
 
! find whether segment {i-1,j}-{i,j} is at the surface or inside
! it only can be on the surface if this is the top or the bottom edge of the object
     
     IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is on the surface

        IF (segment_is_on_surface_of_metal_object) THEN
! #3
! another metal object already has this segment on its border, so the segment is between two metal objects
! return zero epsilon
           eps = 0.0_8
           RETURN
        END IF
! remember that a metal object was found
        segment_is_on_surface_of_metal_object = .TRUE.

     ELSE   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is inside

        eps = 0.0_8
        RETURN

     END IF   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! note that cases when segment is inside a metal object or in between two attached metal objects are already processed above, followed by RETURN calls

  IF (segment_is_on_surface_of_metal_object) THEN
     IF (count.EQ.0) THEN
! segment on metal surface facing vacuum
        eps = 1.0_8
     ELSE
! segment on metal surface facing dielectric
        eps = eps / count
     END IF
     RETURN
  END IF

  IF (count.EQ.0) THEN
! segment is in vacuum
     eps = 1.0_8
     RETURN
  ELSE IF (count.EQ.1) THEN
! segment is inside a dielectric object
     IF (segment_is_inside_dielectric_object) RETURN
! segment is on the surface of a dielectric object, facing vacuum
     eps = (1.0_8 + eps) / 2.0_8
     RETURN
  ELSE IF (count.EQ.2) THEN
     eps = eps / DBLE(count)
     RETURN
  ELSE
     PRINT '("Error-3 in SET_EPS_ISHIFTED for i/j ",2x,i4,2x,i4)', i, j
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

END SUBROUTINE SET_EPS_ISHIFTED


!-----------------------------------------------
!
SUBROUTINE SET_EPS_JSHIFTED(i, j, eps)  ! here point {i,j} is between nodes {i,j-1} and {i,j}

  USE CurrentProblemValues

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
 
  INTEGER, INTENT(IN) :: i, j
  REAL(8), INTENT(OUT) :: eps

  INTEGER count  ! counts dielectric objects
  LOGICAL segment_is_inside_dielectric_object
  LOGICAL segment_is_on_surface_of_metal_object
  INTEGER n1

! find all inner objects owning segment {i,j-1}-{i,j}
! assume that only two dielectric objects may own a common segment
! allow a metal object to be added on top of that

  eps = 0.0_8
  count = 0
  segment_is_inside_dielectric_object = .FALSE.

! first, look through the dielectric objects only
  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.EQ.METAL_WALL) CYCLE

! find an inner object containing point {i,j-1}

     IF (i.LT.whole_object(n1)%ileft) CYCLE
     IF (i.GT.whole_object(n1)%iright) CYCLE
     IF (j-1.LT.whole_object(n1)%jbottom) CYCLE
     IF (j-1.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i,j-1} is either at the surface or inside inner object n1

! find whether this object contains point {i,j}
 
     IF (j.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i,j} is either at the surface or inside inner object n1
! and segment {i,j-1}-{i,j} is owned by the object
     
! find whether segment {i,j-1}-{i,j} is at the surface or inside
! it only can be on the surface if this is the left or the right edge of the object

     IF ((i.EQ.whole_object(n1)%ileft).OR.(i.EQ.whole_object(n1)%iright)) THEN
! segment is on the surface

        count = count + 1
        eps = eps + whole_object(n1)%eps_diel

     ELSE   !### IF ((i.EQ.whole_object(n1)%ileft).OR.(i.EQ.whole_object(n1)%iright)) THEN
! segment is inside
        segment_is_inside_dielectric_object = .TRUE.
        eps = whole_object(n1)%eps_diel
        count = 1
        EXIT

     END IF   !### ((i.EQ.whole_object(n1)%ileft).OR.(i.EQ.whole_object(n1)%iright)) THEN

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! at this point
! if count==0 then the segment is either in/on metal object or in vacuum
! if count==1 then the segment is either inside a dielectric object or on the surface of ONE (not two) dielectric object
! if count==2 then the segment is at the interface between two dielectrics
! if count>=3 this is an error

! second, look through the metal objects only
! if segment belongs to a metal object 
! return zero if it is inside
! if it is on the surface
! return (#1) epsilon of adjacent dielectric object or (#2) 1 if the metal object is in vacuum or (#3) 0 if the metal object is attached to another metal object

  segment_is_on_surface_of_metal_object = .FALSE.

  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

     IF (whole_object(n1)%object_type.NE.METAL_WALL) CYCLE

! find an inner object containing point {i,j-1}

     IF (i.LT.whole_object(n1)%ileft) CYCLE
     IF (i.GT.whole_object(n1)%iright) CYCLE
     IF (j-1.LT.whole_object(n1)%jbottom) CYCLE
     IF (j-1.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i-1,j} is either at the surface or inside inner object n1

! find whether this object contains point {i,j}
 
     IF (j.GT.whole_object(n1)%jtop) CYCLE

! since we are here, point {i,j} is either at the surface or inside inner object n1
! and segment {i,j-1}-{i,j} is owned by the object
 
! find whether segment {i-1,j}-{i,j} is at the surface or inside
! it only can be on the surface if this is the top or the bottom edge of the object
     
     IF ((i.EQ.whole_object(n1)%ileft).OR.(i.EQ.whole_object(n1)%iright)) THEN
! segment is on the surface

        IF (segment_is_on_surface_of_metal_object) THEN
! #3
! another metal object already has this segment on its border, so the segment is between two metal objects
! return zero epsilon
           eps = 0.0_8
           RETURN
        END IF
! remember that a metal object was found
        segment_is_on_surface_of_metal_object = .TRUE.

     ELSE   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN
! segment is inside

        eps = 0.0_8
        RETURN

     END IF   !### IF ((j.EQ.whole_object(n1)%jtop).OR.(j.EQ.whole_object(n1)%jbottom)) THEN

  END DO   !###  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! note that cases when segment is inside a metal object or in between two attached metal objects are already processed above, followed by RETURN calls

  IF (segment_is_on_surface_of_metal_object) THEN
     IF (count.EQ.0) THEN
! segment on metal surface facing vacuum
        eps = 1.0_8
     ELSE
! segment on metal surface facing dielectric
        eps = eps / count
     END IF
     RETURN
  END IF

  IF (count.EQ.0) THEN
! segment is in vacuum
     eps = 1.0_8
     RETURN
  ELSE IF (count.EQ.1) THEN
! segment is inside a dielectric object
     IF (segment_is_inside_dielectric_object) RETURN
! segment is on the surface of a dielectric object, facing vacuum
     eps = (1.0_8 + eps) / 2.0_8
     RETURN
  ELSE IF (count.EQ.2) THEN
     eps = eps / DBLE(count)
     RETURN
  ELSE
     PRINT '("Error-3 in SET_EPS_JSHIFTED for i/j ",2x,i4,2x,i4)', i, j
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

END SUBROUTINE SET_EPS_JSHIFTED

!-----------------------------------------------
!
SUBROUTINE SET_X11_X12_ISHIFTED(i, j, x_11, x_12)  ! here point {i,j} is between nodes {i-1,j} and {i,j}

  USE CurrentProblemValues
  USE BlockAndItsBoundaries, ONLY : indx_x_max, block_has_symmetry_plane_X_left

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
 
  INTEGER, INTENT(IN) :: i, j
  REAL(8), INTENT(OUT) :: x_11, x_12

  INTEGER count_surfs  ! counts surfaces of inner objects
  INTEGER n1
  REAL(8) myx

  x_11 = 0.0_8
  x_12 = 0.0_8

  IF (i.EQ.0) THEN
! point half-a-cell left of the left simulation domain edge
     IF (.NOT.block_has_symmetry_plane_X_left) RETURN
  END IF

  IF (i.EQ.indx_x_max+1) THEN
! point half-a-cell right of the right simulation domain edge
     RETURN
  END IF

  count_surfs = 0

  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! find an inner object containing shifted point {i,j}

     IF (j.LT.whole_object(n1)%jbottom) CYCLE
     IF (j.GT.whole_object(n1)%jtop) CYCLE

     myx = DBLE(i) - 0.5_8

     IF ( (myx.GT.DBLE(whole_object(n1)%ileft)) .AND. (myx.LT.DBLE(whole_object(n1)%iright)) ) THEN

        IF (j.EQ.whole_object(n1)%jbottom) THEN
! point on the bottom surface
           count_surfs = count_surfs + 1
        ELSE IF (j.EQ.whole_object(n1)%jtop) THEN
! point on the top surface
           count_surfs = count_surfs + 1
        ELSE
! point is inside, return zeros
           RETURN
        END IF

     END IF

     IF (count_surfs.GE.2) THEN
! point is between two objects, return zeros
        RETURN 
     END IF

  END DO

! since we are here, point is not inside an inner object (though it may be at its surface)
! return average over neighbor nodes surrounding the shifted point

  IF ((i.EQ.0).AND.(block_has_symmetry_plane_X_left)) THEN
! assume symmetry x_11(-0.5) = x_11(0.5)
     x_11 = 0.5_8 * (X11(i,j) + X11(i+1,j))
     x_12 = 0.5_8 * (X12(i,j) + X12(i+1,j))
     RETURN
  END IF

  x_11 = 0.5_8 * (X11(i-1,j) + X11(i,j))
  x_12 = 0.5_8 * (X12(i-1,j) + X12(i,j))

  RETURN

END SUBROUTINE SET_X11_X12_ISHIFTED


!-----------------------------------------------
!
SUBROUTINE SET_X21_X22_JSHIFTED(i, j, x_21, x_22)  ! here point {i,j} is between nodes {i,j-1} and {i,j}

  USE CurrentProblemValues
  USE BlockAndItsBoundaries, ONLY : indx_y_max

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
 
  INTEGER, INTENT(IN) :: i, j
  REAL(8), INTENT(OUT) :: x_21, x_22

  INTEGER count_surfs  ! counts surfaces of inner objects
  INTEGER n1
  REAL(8) myy

  x_21 = 0.0_8
  x_22 = 0.0_8

  IF (j.EQ.0) THEN
! point half-a-cell under the bottom simulation domain edge
     RETURN
  END IF

  IF (j.EQ.indx_y_max+1) THEN
! point half-a-cell above the top simulation domain edge
     RETURN
  END IF

  count_surfs = 0
  
  DO n1 = N_of_boundary_objects+1, N_of_boundary_and_inner_objects

! find an inner object containing shifted point {i,j}

     IF (i.LT.whole_object(n1)%ileft) CYCLE
     IF (i.GT.whole_object(n1)%iright) CYCLE

     myy = DBLE(j) - 0.5_8

     IF ( (myy.GT.DBLE(whole_object(n1)%jbottom)) .AND. (myy.LT.DBLE(whole_object(n1)%jtop)) ) THEN

        IF (i.EQ.whole_object(n1)%ileft) THEN
! point on the left surface
           count_surfs = count_surfs + 1
        ELSE IF (i.EQ.whole_object(n1)%iright) THEN
! point on the right surface
           count_surfs = count_surfs + 1
        ELSE
! point is inside, return zeros
           RETURN
        END IF

     END IF

     IF (count_surfs.GE.2) THEN
! point is between two objects, return zeros
        RETURN 
     END IF

  END DO

! since we are here, point is not inside an inner object (though it may be at its surface)
! return average over neighbor nodes surrounding the shifted point

  x_21 = 0.5_8 * (X21(i,j-1) + X21(i,j))
  x_22 = 0.5_8 * (X22(i,j-1) + X22(i,j))

  RETURN

END SUBROUTINE SET_X21_X22_JSHIFTED
