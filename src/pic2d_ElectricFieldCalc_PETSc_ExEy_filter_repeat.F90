
!--------------------------------------
!
! procedures InsertData, Solve, and FetchData used below were written by Janhunen Jani Salomon
!
!
SUBROUTINE SOLVE_POTENTIAL_WITH_PETSC

  USE PETSc_Solver
  USE ParallelOperationValues
  USE CurrentProblemValues
  USE BlockAndItsBoundaries
  USE Diagnostics, ONLY : Save_probes_data_T_cntr, N_of_probes_block, Probe_params_block_list, probe_F_block

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  REAL(8), ALLOCATABLE :: queue(:)     ! array of phi-values aligned left-to-right--bottom-to-top
  REAL(8), ALLOCATABLE :: rhsvalue(:)  ! array of right-hand-side-values corresponding to the queue
  INTEGER ALLOC_ERR

  REAL(8) factor_rho

  INTEGER jbegin, jend, ibegin, iend
  INTEGER irow_global
  INTEGER nn

  INTEGER i, j, n, m, i_start, i_end, j_start, j_end, nobj

  INTEGER n1, n3
  REAL(8), ALLOCATABLE :: rbufer(:)

  INTEGER npb
   
  INTEGER i_repeat

  ALLOCATE(   queue(1:block_N_of_nodes_to_solve), STAT = ALLOC_ERR)
  ALLOCATE(rhsvalue(1:block_N_of_nodes_to_solve), STAT = ALLOC_ERR)
 
  n=0
  rhsvalue = 0.0_8
  queue = 0.0_8
  phi = 0.0_8
 
!  factor_rho = -0.25_8 / DBLE(N_of_particles_cell)
  factor_rho = -1.0_8 / DBLE(N_of_particles_cell)

  jbegin = indx_y_min+1
  jend   = indx_y_max-1
  ibegin = indx_x_min+1
  iend   = indx_x_max-1

  IF (Rank_of_process_left.LT.0)  ibegin = indx_x_min
  IF (Rank_of_process_right.LT.0) iend   = indx_x_max
  IF (Rank_of_process_below.LT.0) jbegin = indx_y_min
  IF (Rank_of_process_above.LT.0) jend   = indx_y_max

  irow_global = global_offset

  nn=0

!    j = indx_y_min !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
  IF (jbegin.EQ.indx_y_min) THEN
! boundary object along bottom border
     DO i = ibegin, iend
        nn = nn + 1
        DO n = 1, N_of_local_object_parts_below
           m = index_of_local_object_part_below(n)
           i_start = local_object_part(m)%istart
           i_end   = local_object_part(m)%iend
           nobj   = local_object_part(m)%object_number
           IF ((i.GE.i_start).AND.(i.LE.i_end)) THEN
              IF (whole_object(nobj)%object_type.EQ.METAL_WALL) THEN
                 rhsvalue(nn) = whole_object(nobj)%phi
              ELSE IF (whole_object(nobj)%object_type.EQ.VACUUM_GAP) THEN
                 rhsvalue(nn) = whole_object(nobj)%phi_profile(i)
              END IF
              EXIT
           END IF
        END DO
     END DO
  END IF

  DO j = indx_y_min+1, indx_y_max-1 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!     i = indx_x_min

     IF (ibegin.EQ.indx_x_min) THEN
! boundary object along left border
        nn = nn + 1
        IF (.NOT.block_has_symmetry_plane_X_left) THEN
           DO n = 1, N_of_local_object_parts_left
              m = index_of_local_object_part_left(n)
              j_start = local_object_part(m)%jstart
              j_end   = local_object_part(m)%jend
              nobj   = local_object_part(m)%object_number
              IF ((j.GE.j_start).AND.(j.LE.j_end)) THEN
                 IF (whole_object(nobj)%object_type.EQ.METAL_WALL) THEN
                    rhsvalue(nn) = whole_object(nobj)%phi
                 ELSE IF (whole_object(nobj)%object_type.EQ.VACUUM_GAP) THEN
                    rhsvalue(nn) = whole_object(nobj)%phi_profile(j)
                 END IF
                 EXIT
              END IF
           END DO
        ELSE !IF (whole_object(nobj)%object_type.EQ.SYMMETRY_PLANE) THEN
           rhsvalue(nn) = factor_rho * (rho_i(indx_x_min,j) - rho_e(indx_x_min,j))
        END IF
     END IF

     DO i = indx_x_min+1, indx_x_max-1
        nn = nn + 1
        rhsvalue(nn) = factor_rho * (rho_i(i,j) - rho_e(i,j)) !factor_rho is negative

        IF ((i.EQ.i_given_F_double_period_sys).AND.(j.EQ.j_given_F_double_period_sys)) rhsvalue(nn) = given_F_double_period_sys

     END DO

!     i = indx_x_max

     IF (iend.EQ.indx_x_max) THEN
! boundary object along right border
        nn = nn + 1
        DO n = 1, N_of_local_object_parts_right
           m = index_of_local_object_part_right(n)
           j_start = local_object_part(m)%jstart
           j_end   = local_object_part(m)%jend
           nobj   = local_object_part(m)%object_number
           IF ((j.GE.j_start).AND.(j.LE.j_end)) THEN
              IF (whole_object(nobj)%object_type.EQ.METAL_WALL) THEN
                 rhsvalue(nn) = whole_object(nobj)%phi
              ELSE IF (whole_object(nobj)%object_type.EQ.VACUUM_GAP) THEN
                 rhsvalue(nn) = whole_object(nobj)%phi_profile(j)
              END IF
              EXIT
           END IF
        END DO
     END IF

  END DO   !### DO j = indx_y_min+1, indx_y_max-1

!    j = indx_y_max !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
  IF (jend.EQ.indx_y_max) THEN
! boundary object along top border
     DO i = ibegin, iend
        nn = nn + 1
        DO n = 1, N_of_local_object_parts_above
           m = index_of_local_object_part_above(n)
           i_start = local_object_part(m)%istart
           i_end   = local_object_part(m)%iend
           nobj   = local_object_part(m)%object_number
           IF ((i.GE.i_start).AND.(i.LE.i_end)) THEN
              IF (whole_object(nobj)%object_type.EQ.METAL_WALL) THEN
                 rhsvalue(nn) = whole_object(nobj)%phi
              ELSE IF (whole_object(nobj)%object_type.EQ.VACUUM_GAP) THEN
                 rhsvalue(nn) = whole_object(nobj)%phi_profile(i)
              END IF
              EXIT
           END IF
        END DO
     END DO
  END IF

! fool proof
  IF (nn.NE.block_N_of_nodes_to_solve) THEN
     PRINT '("proc ",i4," :: Error-1 in SOLVE_POTENTIAL_WITH_PETSC ",2(2x,i9))', Rank_of_process, nn, block_N_of_nodes_to_solve
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

! find metal inner objects 
  nn = 0
  DO j = jbegin, jend
     DO i = ibegin, iend
        nn = nn + 1
        DO nobj = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF (whole_object(nobj)%object_type.NE.METAL_WALL) CYCLE
           IF (i.LT.whole_object(nobj)%ileft) CYCLE
           IF (i.GT.whole_object(nobj)%iright) CYCLE
           IF (j.LT.whole_object(nobj)%jbottom) CYCLE
           IF (j.GT.whole_object(nobj)%jtop) CYCLE
           rhsvalue(nn) = whole_object(nobj)%phi
        END DO
     END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!print '("proc ",i4," before InsertData")', Rank_of_process

  call InsertData(nn, rhsvalue)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!print '("proc ",i4," before Solve")', Rank_of_process
!stop

  call Solve

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!print '("proc ",i4," before FetchData")', Rank_of_process
!stop

  call FetchData(nn, queue)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!print '("proc ",i4," after FetchData")', Rank_of_process
!stop

  nn = 0
  DO j = jbegin, jend ! "nodes to solve"
     DO i = ibegin, iend
        nn = nn+1
        phi(i,j) = queue(nn)
     END DO
  END DO

! exchange potential values in overlap nodes with neighbor processes (blocks)
  n1 = indx_y_max - indx_y_min + 1
  n3 = indx_x_max - indx_x_min + 1

  IF (WHITE) THEN  
! "white processes"

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:n1), STAT=ALLOC_ERR)

     IF (Rank_of_process_right.GE.0) THEN
! ## 1 ## send right potential in column left of the right edge
        rbufer(1:n1) = phi(indx_x_max-1, indx_y_min:indx_y_max)
        CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (Rank_of_process_left.GE.0) THEN
! ## 2 ## send left potential in column right of the left edge
        rbufer(1:n1) = phi(indx_x_min+1, indx_y_min:indx_y_max)
        CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (Rank_of_process_left.GE.0) THEN
! ## 3 ## receive from left potential along the left edge
        CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process_left, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_min, indx_y_min:indx_y_max) = rbufer(1:n1)   
     END IF

     IF (Rank_of_process_right.GE.0) THEN
! ## 4 ## receive from right potential along the right edge
        CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process_right, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_max, indx_y_min:indx_y_max) = rbufer(1:n1)
     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

     IF (Rank_of_process_above.GE.0) THEN
! ## 5 ## send up potential in the row below the top edge
        rbufer(1:n3) = phi(indx_x_min:indx_x_max, indx_y_max-1)
        CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_above, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (Rank_of_process_below.GE.0) THEN
! ## 6 ## send down potential in the row above the bottom edge
        rbufer(1:n3) = phi(indx_x_min:indx_x_max, indx_y_min+1)
        CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_below, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (Rank_of_process_below.GE.0) THEN
! ## 7 ## receive from below potential in the bottom edge
        CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_below, Rank_of_process_below, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_min:indx_x_max, indx_y_min) = rbufer(1:n3)
     END IF

     IF (Rank_of_process_above.GE.0) THEN
! ## 8 ## receive from above potential in the top edge
        CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_above, Rank_of_process_above, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_min:indx_x_max, indx_y_max) = rbufer(1:n3)
     END IF

  ELSE
! "black" processes

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:n1), STAT=ALLOC_ERR)

     IF (Rank_of_process_left.GE.0) THEN
! ## 1 ## receive from left potential along the left edge
        CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process_left, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_min, indx_y_min:indx_y_max) = rbufer(1:n1)   
     END IF

     IF (Rank_of_process_right.GE.0) THEN
! ## 2 ## receive from right potential along the right edge
        CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process_right, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_max, indx_y_min:indx_y_max) = rbufer(1:n1)
     END IF

     IF (Rank_of_process_right.GE.0) THEN
! ## 3 ## send right potential in column left of the right edge
        rbufer(1:n1) = phi(indx_x_max-1, indx_y_min:indx_y_max)
        CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (Rank_of_process_left.GE.0) THEN
! ## 4 ## send left potential in column right of the left edge
        rbufer(1:n1) = phi(indx_x_min+1, indx_y_min:indx_y_max)
        CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

     IF (Rank_of_process_below.GE.0) THEN
! ## 5 ## receive from below potential in the bottom edge
        CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_below, Rank_of_process_below, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_min:indx_x_max, indx_y_min) = rbufer(1:n3)
     END IF

     IF (Rank_of_process_above.GE.0) THEN
! ## 6 ## receive from above potential in the top edge
        CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_above, Rank_of_process_above, MPI_COMM_WORLD, stattus, ierr)
        phi(indx_x_min:indx_x_max, indx_y_max) = rbufer(1:n3)
     END IF

     IF (Rank_of_process_above.GE.0) THEN
! ## 7 ## send up potential in the row below the top edge
        rbufer(1:n3) = phi(indx_x_min:indx_x_max, indx_y_max-1)
        CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_above, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

     IF (Rank_of_process_below.GE.0) THEN
! ## 8 ## send down potential in the row above the bottom edge
        rbufer(1:n3) = phi(indx_x_min:indx_x_max, indx_y_min+1)
        CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_below, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END IF

  END IF !black/white selection

! filtering:

  i_repeat = 1
  DO WHILE(i_repeat.LE.n_repeat)

     i_repeat = i_repeat + 1
     ALLOCATE(phi_avg(indx_x_min:indx_x_max, indx_y_min:indx_y_max), STAT= ALLOC_ERR)
!    2D binomial filter
     phi_avg(indx_x_min, indx_y_min:indx_y_max) = phi(indx_x_min, indx_y_min:indx_y_max)
     phi_avg(indx_x_max, indx_y_min:indx_y_max) = phi(indx_x_max, indx_y_min:indx_y_max)
     phi_avg(indx_x_min:indx_x_max, indx_y_min) = phi(indx_x_min:indx_x_max, indx_y_min)
     phi_avg(indx_x_min:indx_x_max, indx_y_max) = phi(indx_x_min:indx_x_max, indx_y_max)

     DO i = indx_x_min, indx_x_max 
        DO j = indx_y_min + 1, indx_y_max - 1
           phi_avg(i, j) = 0.25_8 * ( phi(i, j-1) + 2.0_8 * phi(i, j) + phi(i, j+1) )
        END DO
     END DO

     DO j = indx_y_min + 1, indx_y_max - 1
        DO i = indx_x_min + 1, indx_x_max - 1
           phi(i, j) = 0.25_8 * ( phi_avg(i-1, j) + 2.0_8 * phi_avg(i, j) + phi_avg(i+1, j) )
        END DO
     END DO

     IF (block_has_symmetry_plane_X_left) THEN
        DO j = indx_y_min + 1, indx_y_max - 1
           phi(indx_x_min, j) = 0.5_8 * (phi_avg(indx_x_min, j) + phi_avg(indx_x_min + 1, j))
        END DO
     END IF  

     DEALLOCATE(phi_avg, STAT = ALLOC_ERR)

     IF (WHITE) THEN
!    "white processes"

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n1), STAT=ALLOC_ERR)

        IF (Rank_of_process_right.GE.0) THEN
!    ## 1 ## send right potential in column left of the right edge
           rbufer(1:n1) = phi(indx_x_max-1, indx_y_min:indx_y_max)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process, MPI_COMM_WORLD, request, ierr)
        END IF

        IF (Rank_of_process_left.GE.0) THEN
!    ## 2 ## send left potential in column right of the left edge
           rbufer(1:n1) = phi(indx_x_min+1, indx_y_min:indx_y_max)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process, MPI_COMM_WORLD, request, ierr)
        END IF

        IF (Rank_of_process_left.GE.0) THEN
!    ## 3 ## receive from left potential along the left edge
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process_left, MPI_COMM_WORLD, stattus, ierr)
           phi(indx_x_min, indx_y_min:indx_y_max) = rbufer(1:n1)
        END IF

        IF (Rank_of_process_right.GE.0) THEN
!    ## 4 ## receive from right potential along the right edge
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process_right, MPI_COMM_WORLD, stattus, ierr)
           phi(indx_x_max, indx_y_min:indx_y_max) = rbufer(1:n1)
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

        IF (Rank_of_process_above.GE.0) THEN
!    ## 5 ## send up potential in the row below the top edge
           rbufer(1:n3) = phi(indx_x_min:indx_x_max, indx_y_max-1)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_above, Rank_of_process, MPI_COMM_WORLD, request, ierr)
        END IF

        IF (Rank_of_process_below.GE.0) THEN
!    ## 6 ## send down potential in the row above the bottom edge
           rbufer(1:n3) = phi(indx_x_min:indx_x_max, indx_y_min+1)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_below, Rank_of_process, MPI_COMM_WORLD, request, ierr)
        END IF

        IF (Rank_of_process_below.GE.0) THEN
!    ## 7 ## receive from below potential in the bottom edge
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_below, Rank_of_process_below, MPI_COMM_WORLD, stattus, ierr)
           phi(indx_x_min:indx_x_max, indx_y_min) = rbufer(1:n3)
        END IF

        IF (Rank_of_process_above.GE.0) THEN
!    ## 8 ## receive from above potential in the top edge
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_above, Rank_of_process_above, MPI_COMM_WORLD, stattus, ierr)
           phi(indx_x_min:indx_x_max, indx_y_max) = rbufer(1:n3)
        END IF

     ELSE
!    "black" processes  
        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n1), STAT=ALLOC_ERR)

        IF (Rank_of_process_left.GE.0) THEN
!    ## 1 ## receive from left potential along the left edge
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process_left, MPI_COMM_WORLD, stattus, ierr)
           phi(indx_x_min, indx_y_min:indx_y_max) = rbufer(1:n1)
        END IF

        IF (Rank_of_process_right.GE.0) THEN
!    ## 2 ## receive from right potential along the right edge
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process_right, MPI_COMM_WORLD, stattus, ierr)
           phi(indx_x_max, indx_y_min:indx_y_max) = rbufer(1:n1)
        END IF

        IF (Rank_of_process_right.GE.0) THEN
!    ## 3 ## send right potential in column left of the right edge
           rbufer(1:n1) = phi(indx_x_max-1, indx_y_min:indx_y_max)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_right, Rank_of_process, MPI_COMM_WORLD, request, ierr)
        END IF

        IF (Rank_of_process_left.GE.0) THEN
!    ## 4 ## send left potential in column right of the left edge
           rbufer(1:n1) = phi(indx_x_min+1, indx_y_min:indx_y_max)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_of_process_left, Rank_of_process, MPI_COMM_WORLD, request, ierr)
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

        IF (Rank_of_process_below.GE.0) THEN
!    ## 5 ## receive from below potential in the bottom edge
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_below, Rank_of_process_below, MPI_COMM_WORLD, stattus, ierr)
           phi(indx_x_min:indx_x_max, indx_y_min) = rbufer(1:n3)
        END IF

        IF (Rank_of_process_above.GE.0) THEN
!    ## 6 ## receive from above potential in the top edge
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_above, Rank_of_process_above, MPI_COMM_WORLD, stattus, ierr)
           phi(indx_x_min:indx_x_max, indx_y_max) = rbufer(1:n3)
        END IF

        IF (Rank_of_process_above.GE.0) THEN
!    ## 7 ## send up potential in the row below the top edge
           rbufer(1:n3) = phi(indx_x_min:indx_x_max, indx_y_max-1)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_above, Rank_of_process, MPI_COMM_WORLD, request, ierr)
        END IF

        IF (Rank_of_process_below.GE.0) THEN
!    ## 8 ## send down potential in the row above the bottom edge
           rbufer(1:n3) = phi(indx_x_min:indx_x_max, indx_y_min+1)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_of_process_below, Rank_of_process, MPI_COMM_WORLD, request, ierr)
        END IF
     
     END IF  !black/white selecton  

  END DO ! do-while filtering iterations

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)

! ################ diagnostics, electrostatic potential #################
  IF (1+T_cntr.EQ.Save_probes_data_T_cntr) THEN
     DO npb = 1, N_of_probes_block
        i = Probe_params_block_list(1,npb)
        j = Probe_params_block_list(2,npb)
        probe_F_block(npb) = phi(i,j)
     END DO
  END IF

  DEALLOCATE(rhsvalue, STAT=ALLOC_ERR)
  DEALLOCATE(queue,   stat=alloc_err)
 
END SUBROUTINE SOLVE_POTENTIAL_WITH_PETSC


!-------------------------------------------------------
!
SUBROUTINE CALCULATE_ELECTRIC_FIELD

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE BlockAndItsBoundaries
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  REAL(8) factor1_E_from_F, f1

  INTEGER i, j, k, pos, n1, n3, pos1, pos2
  INTEGER bufsize 

  REAL(8), ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR

  factor1_E_from_F = F_scale_V / (E_scale_Vm * delta_x_m)
  f1 = factor1_E_from_F

  IF (cluster_rank_key.EQ.0) THEN ! process is the field master
     c_phi = 0.0_8
! values of the potential in the blocks are already set, including those on edges           
     DO j = indx_y_min, indx_y_max 
        DO i = indx_x_min,  indx_x_max 
           c_phi(i, j) = phi(i, j)
        END DO
     END DO

! receive the values of potential from field calculators, row by row
     DO k = 2, cluster_N_blocks

        bufsize = (field_calculator(k)%indx_x_max - field_calculator(k)%indx_x_min + 1) * &
                & (field_calculator(k)%indx_y_max - field_calculator(k)%indx_y_min + 1)
        ALLOCATE(rbufer(bufsize), STAT = ALLOC_ERR)
        CALL MPI_RECV(rbufer, bufsize, MPI_DOUBLE_PRECISION, field_calculator(k)%rank, field_calculator(k)%rank, MPI_COMM_WORLD, stattus, ierr)
        
        pos = 0
        DO j = field_calculator(k)%indx_y_min, field_calculator(k)%indx_y_max
           DO i = field_calculator(k)%indx_x_min, field_calculator(k)%indx_x_max
              pos = pos + 1
              c_phi(i, j) = rbufer(pos)
           END DO
        END DO

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)

     END DO ! loop over field calculators

! Exchange boundary values of c_phi with neighbor masters. 
! All values at (c_indx_x_min : c_indx_x_max, c_indx_y_min : c_indx_y_max) are already set.
! All clusters need values of 
! Ex(c_indx_x_min - 1 : c_indx_x_max, c_indx_y_min : c_indx_y_max) and
! Ey(c_indx_x_min : c_indx_x_max, c_indx_y_min - 1 : c_indx_y_max).
! Therefore, c_phi is needed on "ghost edges" at 
! (c_indx_x_min - 1, c_indx_y_min : c_indx_y_max),
! (c_indx_x_max + 1, c_indx_y_min : c_indx_y_max),
! (c_indx_x_min : c_indx_x_max, c_indx_y_min - 1),
! (c_indx_x_min : c_indx_x_max, c_indx_y_max + 1).
! Either receive extra values from adjacent clusters, or assign at the outer boundary.

     n1 = c_indx_y_max - c_indx_y_min + 3 ! extra rows and columns, in fact, run from min to max only
     n3 = c_indx_x_max - c_indx_x_min + 3 ! extra corner nodes are not used 

! master processes exchange the values at "ghost" nodes     
     IF (WHITE_CLUSTER) THEN  
! "white processes"
        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n1), STAT=ALLOC_ERR)

        IF (Rank_horizontal_right.GE.0) THEN
! ## 1 ## send right c_phi in the column left of the right edge
           rbufer(1:n1) = c_phi(c_indx_x_max - 2, c_indx_y_min - 1 : c_indx_y_max + 1)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 2 ## send left c_phi in the column to the right of the left edge
           rbufer(1:n1) = c_phi(c_indx_x_min + 2, c_indx_y_min - 1 : c_indx_y_max + 1)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 3 ## receive from left c_phi in the extra column on the left
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           c_phi(c_indx_x_min - 1, c_indx_y_min - 1 : c_indx_y_max + 1) = rbufer(1:n1)   
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 4 ## receive from right c_phi in the extra column on the right
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           c_phi(c_indx_x_max + 1, c_indx_y_min - 1 : c_indx_y_max + 1) = rbufer(1:n1)
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

        IF (Rank_horizontal_above.GE.0) THEN
! ## 5 ## send up c_phi in the row below the top edge
           rbufer(1:n3) = c_phi(c_indx_x_min - 1 : c_indx_x_max + 1, c_indx_y_max - 2)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 6 ## send down c_phi in the row above the bottom edge
           rbufer(1:n3) = c_phi(c_indx_x_min - 1 : c_indx_x_max + 1, c_indx_y_min + 2)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 7 ## receive from below c_phi in the extra row below
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           c_phi(c_indx_x_min - 1 : c_indx_x_max + 1, c_indx_y_min - 1) = rbufer(1:n3)
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 8 ## receive from above c_phi in the extra row above
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           c_phi(c_indx_x_min - 1 : c_indx_x_max + 1, c_indx_y_max + 1) = rbufer(1:n3)
        END IF

     ELSE
! "black" processes
        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n1), STAT=ALLOC_ERR)

        IF (Rank_horizontal_left.GE.0) THEN
! ## 1 ## receive from left c_phi in the extra column on the left
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           c_phi(c_indx_x_min - 1, c_indx_y_min - 1 : c_indx_y_max + 1) = rbufer(1:n1)   
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 2 ## receive from right c_phi in the extra column on the right:
           CALL MPI_RECV(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           c_phi(c_indx_x_max + 1, c_indx_y_min - 1 : c_indx_y_max + 1) = rbufer(1:n1)
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 3 ## send right c_phi in the column to the left of the right edge
           rbufer(1:n1) = c_phi(c_indx_x_max - 2, c_indx_y_min - 1 : c_indx_y_max + 1)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 4 ## send left c_phi in the column to the right of the left edge
           rbufer(1:n1)= c_phi(c_indx_x_min + 2, c_indx_y_min - 1 : c_indx_y_max + 1)
           CALL MPI_SEND(rbufer, n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:n3), STAT=ALLOC_ERR)

        IF (Rank_horizontal_below.GE.0) THEN
! ## 5 ## receive from below c_phi in the extra row below
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           c_phi(c_indx_x_min - 1 : c_indx_x_max + 1, c_indx_y_min - 1) = rbufer(1:n3)
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 6 ## receive from above c_phi in the extra row above
           CALL MPI_RECV(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           c_phi(c_indx_x_min - 1 : c_indx_x_max + 1, c_indx_y_max + 1) = rbufer(1:n3)
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 7 ## send up c_phi in the row below the top edge
           rbufer(1:n3) = c_phi(c_indx_x_min - 1 : c_indx_x_max + 1, c_indx_y_max - 2)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 8 ## send down c_phi in the row above the bottom edge
           rbufer(1:n3) = c_phi(c_indx_x_min - 1 : c_indx_x_max + 1, c_indx_y_min + 2)
           CALL MPI_SEND(rbufer, n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

     END IF !black/white selection  for "ghost" values transfer

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)

! ready to calculate Ex and Ey on respective grids
! First, go over the nodes where Ex and Ey are defined on any cluster:
     DO j = c_indx_y_min, c_indx_y_max - 1 
        DO i = c_indx_x_min, c_indx_x_max 
           EY(i, j) = f1 * (c_phi(i, j) - c_phi(i, j+1))
        END DO
     END DO   

     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max - 1
           EX(i, j) = f1 * (c_phi(i, j) - c_phi(i+1, j))
        END DO
     END DO
! Here the values of Ex, Ey at the ghost nodes located outside the domain can be set according to adopted extrapolation method.
! These values are used for interpolating the field to particle positions
     i = c_indx_x_min - 1
     IF (Rank_horizontal_left.LT.0) THEN
!     IF (Rank_of_master_left.LT.0) THEN
        DO j = c_indx_y_min, c_indx_y_max
           EX(i, j) = EX(c_indx_x_min, j)
        END DO
     ELSE !use the values of c_phi(c_indx_x_min - 1, j) taken from adjacent cluster:
        DO j = c_indx_y_min, c_indx_y_max
           EX(i, j) = f1 * (c_phi(i, j) - c_phi(i+1, j))     
        END DO   
     END IF   
   
     i = c_indx_x_max 
     IF (Rank_horizontal_right.LT.0) THEN
!     IF (Rank_of_master_right.LT.0) THEN
        DO j = c_indx_y_min, c_indx_y_max
           EX(i, j) = EX(c_indx_x_max - 1, j)
        END DO
     ELSE !use the values of c_phi at i = c_indx_x_max + 1 taken from adjacent cluster:
        DO j = c_indx_y_min, c_indx_y_max
           EX(i, j) = f1 * (c_phi(i, j) - c_phi(i+1, j))
        END DO     
     END IF
    
     j = c_indx_y_max 
     IF (Rank_horizontal_above.LT.0) THEN
!     IF (Rank_of_master_above.LT.0) THEN
        DO i = c_indx_x_min, c_indx_x_max
           EY(i, j) = EY(i, c_indx_y_max - 1)
         END DO
      ELSE !use the values of c_phi at j = c_indx_y_max + 1 taken from adjacent cluster:
         DO i = c_indx_x_min, c_indx_x_max
           EY(i, j) = f1 * (c_phi(i, j) - c_phi(i, j+1))
         END DO     
      END IF   

      j = c_indx_y_min - 1
      IF (Rank_horizontal_below.LT.0) THEN
!      IF (Rank_of_master_below.LT.0) THEN
        DO i = c_indx_x_min, c_indx_x_max
           EY(i, j) = EY(i, c_indx_y_min)
        END DO
      ELSE !use the values of c_phi at j = c_indx_y_min - 1 transferred from adjacent cluster:
         DO i = c_indx_x_min, c_indx_x_max
           EY(i, j) = f1 * (c_phi(i, j) - c_phi(i, j+1))
         END DO  
      END IF

     IF (symmetry_plane_X_left) THEN
! enforce boundary condition EX = 0 at the symmetry plane
        EX(c_indx_x_min - 1, c_indx_y_min:c_indx_y_max) = - EX(c_indx_x_min, c_indx_y_min:c_indx_y_max)
     END IF

! ready to send complete field array to all members of the cluster

  ELSE !the block is NOT the field calculator master in its cluster

! send values of potential to the block's field master, row by row
     bufsize = (indx_x_max - indx_x_min + 1) * (indx_y_max - indx_y_min + 1)
     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)
     ALLOCATE(rbufer(1:bufsize), STAT = ALLOC_ERR)
     pos = 0
     DO j = indx_y_min, indx_y_max
        DO i = indx_x_min, indx_x_max
           pos = pos + 1
           rbufer(pos) = phi(i, j)
        END DO  
     END DO
     CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, field_master, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     DEALLOCATE(rbufer, STAT=ALLOC_ERR)
! the block is ready to receive complete field array from the master of the cluster

  END IF ! field master selection

! master of the cluster distributes complete Ex, Ey arrays on respective grids:
  call mpi_barrier(mpi_comm_world, ierr)
  bufsize = (c_indx_x_max - c_indx_x_min + 2) * (c_indx_y_max - c_indx_y_min + 1)
  CALL MPI_BCAST(EX, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr) 

  call mpi_barrier(mpi_comm_world, ierr)
  bufsize = (c_indx_x_max - c_indx_x_min + 1) * (c_indx_y_max - c_indx_y_min + 2)
  CALL MPI_BCAST(EY, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr)

! each member (master and non-master) in each cluster accumulates fields for ions in the whole cluster domain 
! (in this case there is no need for communications at all)
  DO j = c_indx_y_min, c_indx_y_max
     DO i = c_indx_x_min - 1, c_indx_x_max
        acc_EX(i,j) = acc_EX(i,j) + EX(i,j)
     END DO
  END DO
  DO i = c_indx_x_min, c_indx_x_max
     DO j = c_indx_y_min - 1, c_indx_y_max
        acc_EY(i,j) = acc_EY(i,j) + EY(i,j)
     END DO
  END DO
  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)
  
END SUBROUTINE CALCULATE_ELECTRIC_FIELD

!------------------------------------------------------
!
SUBROUTINE CLEAR_ACCUMULATED_FIELDS

  USE CurrentProblemValues
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INTEGER i, j

  DO j = c_indx_y_min, c_indx_y_max
     DO i = c_indx_x_min - 1, c_indx_x_max
        acc_EX(i,j) = 0.0_8
     END DO
  END DO
  DO i = c_indx_x_min, c_indx_x_max
     DO j = c_indx_y_min - 1, c_indx_y_max
        acc_EY(i,j) = 0.0_8
     END DO
  END DO
 
END SUBROUTINE CLEAR_ACCUMULATED_FIELDS
