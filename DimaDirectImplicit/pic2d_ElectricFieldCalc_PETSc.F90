
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
  USE ClusterAndItsBoundaries

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

  INTEGER pos, k
  INTEGER bufsize

  ALLOCATE(   queue(1:block_N_of_nodes_to_solve), STAT = ALLOC_ERR)
  ALLOCATE(rhsvalue(1:block_N_of_nodes_to_solve), STAT = ALLOC_ERR)
 
  n=0
  rhsvalue = 0.0_8
  queue = 0.0_8

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
           rhsvalue(nn) = factor_rho * rho(indx_x_min,j)
        END IF
     END IF

     DO i = indx_x_min+1, indx_x_max-1
        nn = nn + 1
        rhsvalue(nn) = factor_rho * rho(i,j)

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

  call InsertData(nn, rhsvalue)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  call Solve

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  call FetchData(nn, queue)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

!call save_2d_array('rhs',ibegin,iend,jbegin,jend,rhsvalue)
!call save_2d_array('que',ibegin,iend,jbegin,jend,queue)

  DEALLOCATE(rhsvalue, STAT=ALLOC_ERR)

  IF (cluster_rank_key.NE.0) THEN
! non-master processes send their pieces to masters and quit
     CALL MPI_SEND(queue, block_N_of_nodes_to_solve, MPI_DOUBLE_PRECISION, field_master, Rank_of_process, MPI_COMM_WORLD, request, ierr)
!print '(2(2x,i4),2x,i6)', field_master, Rank_of_process, block_N_of_nodes_to_solve 
     DEALLOCATE(queue, STAT = ALLOC_ERR)
     RETURN
  END IF

! the rest is performed by masters only

! save own contribution
  pos = 0
  DO j = jbegin, jend
     DO i = ibegin, iend
        pos = pos + 1
        c_phi_ext(i,j) = queue(pos)
     END DO
  END DO
  DEALLOCATE(queue, STAT = ALLOC_ERR)

! receive pieces from field calculators
  DO k = 2, cluster_N_blocks

     bufsize = (field_calculator(k)%iend - field_calculator(k)%ibegin + 1) * &
             & (field_calculator(k)%jend - field_calculator(k)%jbegin + 1)

!print '(2(2x,i4),2x,i6,5(2x,i4))', Rank_of_process, field_calculator(k)%rank, bufsize, k, field_calculator(k)%ibegin, field_calculator(k)%iend, field_calculator(k)%jbegin, field_calculator(k)%jend

     ALLOCATE(rbufer(bufsize), STAT = ALLOC_ERR)
     CALL MPI_RECV(rbufer, bufsize, MPI_DOUBLE_PRECISION, field_calculator(k)%rank, field_calculator(k)%rank, MPI_COMM_WORLD, stattus, ierr)
     pos = 0
     DO j = field_calculator(k)%jbegin, field_calculator(k)%jend
        DO i = field_calculator(k)%ibegin, field_calculator(k)%iend
           pos = pos + 1
           c_phi_ext(i, j) = rbufer(pos)
        END DO
     END DO    
     DEALLOCATE(rbufer, STAT = ALLOC_ERR)
  END DO ! loop over field calculators

!call save_phi

  n1 = c_indx_y_max - c_indx_y_min + 3
  n3 = c_indx_x_max - c_indx_x_min + 3

! master processes exchange the values at "ghost" nodes     
  IF (WHITE_CLUSTER) THEN  
! "white processes"
     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     bufsize = n1+n1
     ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

     IF (Rank_horizontal_right.GE.0) THEN
! ## 1 ## send right 2 columns left of the non-extended right edge
        rbufer(1:n1)         = c_phi_ext(c_indx_x_max-2, c_indx_y_min-1:c_indx_y_max+1)
        rbufer(n1+1:bufsize) = c_phi_ext(c_indx_x_max-1, c_indx_y_min-1:c_indx_y_max+1)
        CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_left.GE.0) THEN
! ## 2 ## send left 2 columns right of the non-extended left edge
        rbufer(1:n1)         = c_phi_ext(c_indx_x_min + 1, c_indx_y_min-1:c_indx_y_max+1)
        rbufer(n1+1:bufsize) = c_phi_ext(c_indx_x_min + 2, c_indx_y_min-1:c_indx_y_max+1)
        CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_left.GE.0) THEN
! ## 3 ## receive from left 2 left-most columns
        CALL MPI_RECV(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
        c_phi_ext(c_indx_x_min-1, c_indx_y_min-1:c_indx_y_max+1) = rbufer(1:n1)
        c_phi_ext(c_indx_x_min,   c_indx_y_min-1:c_indx_y_max+1) = rbufer(n1+1:bufsize)
     END IF

     IF (Rank_horizontal_right.GE.0) THEN
! ## 4 ## receive from right 2 right-most columns
        CALL MPI_RECV(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
        c_phi_ext(c_indx_x_max,   c_indx_y_min-1:c_indx_y_max+1) = rbufer(1:n1)
        c_phi_ext(c_indx_x_max+1, c_indx_y_min-1:c_indx_y_max+1) = rbufer(n1+1:bufsize)
     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     bufsize = n3+n3
     ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

     IF (Rank_horizontal_above.GE.0) THEN
! ## 5 ## send up 2 rows below the non-extended top edge
        rbufer(1:n3)         = c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_max-2)
        rbufer(n3+1:bufsize) = c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_max-1)
        CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_below.GE.0) THEN
! ## 6 ## send down 2 rows above the non-extended bottom edge
        rbufer(1:n3)         = c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_min+1)
        rbufer(n3+1:bufsize) = c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_min+2)
        CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_below.GE.0) THEN
! ## 7 ## receive from below 2 bottom-most rows
        CALL MPI_RECV(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
        c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_min-1) = rbufer(1:n3)
        c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_min)   = rbufer(n3+1:bufsize)
     END IF

     IF (Rank_horizontal_above.GE.0) THEN
! ## 8 ## receive from above 2 top-most rows
        CALL MPI_RECV(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
        c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_max)   = rbufer(1:n3)
        c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_max+1) = rbufer(n3+1:bufsize)
     END IF

  ELSE
! "black" processes
     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     bufsize = n1+n1
     ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

     IF (Rank_horizontal_left.GE.0) THEN
! ## 1 ## receive from left 2 left-most columns
        CALL MPI_RECV(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
        c_phi_ext(c_indx_x_min-1, c_indx_y_min-1:c_indx_y_max+1) = rbufer(1:n1)
        c_phi_ext(c_indx_x_min,   c_indx_y_min-1:c_indx_y_max+1) = rbufer(n1+1:bufsize)
     END IF
     
     IF (Rank_horizontal_right.GE.0) THEN
! ## 2 ## receive from right 2 right-most columns
        CALL MPI_RECV(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
        c_phi_ext(c_indx_x_max,   c_indx_y_min-1:c_indx_y_max+1) = rbufer(1:n1)
        c_phi_ext(c_indx_x_max+1, c_indx_y_min-1:c_indx_y_max+1) = rbufer(n1+1:bufsize)
     END IF

     IF (Rank_horizontal_right.GE.0) THEN
! ## 3 ## send right 2 columns left of the non-extended right edge
        rbufer(1:n1)         = c_phi_ext(c_indx_x_max-2, c_indx_y_min-1:c_indx_y_max+1)
        rbufer(n1+1:bufsize) = c_phi_ext(c_indx_x_max-1, c_indx_y_min-1:c_indx_y_max+1)
        CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_left.GE.0) THEN
! ## 4 ## send left 2 columns right of the non-extended left edge
        rbufer(1:n1)         = c_phi_ext(c_indx_x_min + 1, c_indx_y_min-1:c_indx_y_max+1)
        rbufer(n1+1:bufsize) = c_phi_ext(c_indx_x_min + 2, c_indx_y_min-1:c_indx_y_max+1)
        CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     bufsize = n3+n3
     ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

     IF (Rank_horizontal_below.GE.0) THEN
! ## 5 ## receive from below 2 bottom-most rows
        CALL MPI_RECV(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
        c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_min-1) = rbufer(1:n3)
        c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_min)   = rbufer(n3+1:bufsize)
     END IF

     IF (Rank_horizontal_above.GE.0) THEN
! ## 6 ## receive from above 2 top-most rows
        CALL MPI_RECV(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
        c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_max)   = rbufer(1:n3)
        c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_max+1) = rbufer(n3+1:bufsize)
     END IF

     IF (Rank_horizontal_above.GE.0) THEN
! ## 7 ## send up 2 rows below the non-extended top edge
        rbufer(1:n3)         = c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_max-2)
        rbufer(n3+1:bufsize) = c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_max-1)
        CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_below.GE.0) THEN
! ## 8 ## send down 2 rows above the non-extended bottom edge
        rbufer(1:n3)         = c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_min+1)
        rbufer(n3+1:bufsize) = c_phi_ext(c_indx_x_min-1:c_indx_x_max+1, c_indx_y_min+2)
        CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

  END IF !black/white selection  for "ghost" values transfer

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)

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

  INTEGER npc, npa

  factor1_E_from_F = F_scale_V / (E_scale_Vm * delta_x_m)
  f1 = factor1_E_from_F

  IF (cluster_rank_key.EQ.0) THEN 
! process is the field master

! EX

     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min-1, c_indx_x_max
           EX(i, j) = f1 * (c_phi_ext(i, j) - c_phi_ext(i+1, j))
        END DO
     END DO

     IF (Rank_horizontal_left.LT.0) THEN
        DO j = c_indx_y_min, c_indx_y_max
           EX(c_indx_x_min-1, j) = EX(c_indx_x_min, j)
        END DO
     END IF

     IF (Rank_horizontal_right.LT.0) THEN
        DO j = c_indx_y_min, c_indx_y_max
           EX(c_indx_x_max, j) = EX(c_indx_x_max-1, j)
        END DO
     END IF

! EY

     DO j = c_indx_y_min-1, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max 
           EY(i, j) = f1 * (c_phi_ext(i, j) - c_phi_ext(i, j+1))
        END DO
     END DO

     IF (Rank_horizontal_above.LT.0) THEN
        DO i = c_indx_x_min, c_indx_x_max
           EY(i, c_indx_y_max) = EY(i, c_indx_y_max-1)
        END DO
     END IF

     IF (Rank_horizontal_below.LT.0) THEN
        DO i = c_indx_x_min, c_indx_x_max
           EY(i, c_indx_y_min-1) = EY(i, c_indx_y_min)
        END DO
     END IF

     IF (symmetry_plane_X_left) THEN
! enforce boundary condition EX = 0 at the symmetry plane
        EX(c_indx_x_min-1, c_indx_y_min:c_indx_y_max) = - EX(c_indx_x_min, c_indx_y_min:c_indx_y_max)
     END IF

! ready to send complete field array to all members of the cluster

  END IF ! field master selection   !###   IF (cluster_rank_key.EQ.0) THEN

! master of the cluster distributes complete Ex, Ey arrays on respective grids:
  call mpi_barrier(mpi_comm_world, ierr)
  bufsize = (c_indx_x_max - c_indx_x_min + 2) * (c_indx_y_max - c_indx_y_min + 1)
  CALL MPI_BCAST(EX, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr) 

  call mpi_barrier(mpi_comm_world, ierr)
  bufsize = (c_indx_x_max - c_indx_x_min + 1) * (c_indx_y_max - c_indx_y_min + 2)
  CALL MPI_BCAST(EY, bufsize, MPI_DOUBLE_PRECISION, 0, COMM_CLUSTER, ierr)

END SUBROUTINE CALCULATE_ELECTRIC_FIELD
