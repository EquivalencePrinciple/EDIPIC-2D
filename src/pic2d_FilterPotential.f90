!-----------------------------------------
!
SUBROUTINE FILTER_POTENTIAL

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER i, j
  REAL(8) temp

  INTEGER n1, n3
  REAL(8), ALLOCATABLE :: rbufer(:)
  INTEGER ALLOC_ERR

  INTEGER pos, k
  INTEGER bufsize

! masters only
  IF (cluster_rank_key.NE.0) RETURN

! filter along x
  DO j = MAX(c_indx_y_min,1), MIN(c_indx_y_max,global_maximal_j-1)
! sweep in positive x-direction
     temp = c_phi_ext(c_indx_x_min,j)
     DO i = c_indx_x_min, c_indx_x_max-1
        c_phi_ext(i,j) = c_phi_ext(i,j) + c_phi_ext(i+1,j)
     END DO
! sweep in negative x-direction
     DO i =  c_indx_x_max-1, c_indx_x_min+1, -1
        c_phi_ext(i,j) = c_phi_ext(i,j) + c_phi_ext(i-1,j)
     END DO
     c_phi_ext(c_indx_x_min,j) = temp
  END DO

! filter along y
  DO i = c_indx_x_min+1, c_indx_x_max-1
! sweep in positive y-direction
     temp = c_phi_ext(i,c_indx_y_min)
     DO j = c_indx_y_min, c_indx_y_max-1
        c_phi_ext(i,j) = c_phi_ext(i,j) + c_phi_ext(i,j+1)
     END DO
! sweep in negative x-direction
     DO j = c_indx_y_max-1, c_indx_y_min+1, -1
        c_phi_ext(i,j) = c_phi_ext(i,j) + c_phi_ext(i,j-1)
     END DO
     c_phi_ext(i,c_indx_y_min) = temp
  END DO

! multiply inner points by by 1/16
  DO j = c_indx_y_min+1, c_indx_y_max-1
     DO i = c_indx_x_min+1, c_indx_x_max-1
        c_phi_ext(i,j) = c_phi_ext(i,j) * 0.0625_8
     END DO
  END DO

! synchronize overlaps

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

END SUBROUTINE FILTER_POTENTIAL
