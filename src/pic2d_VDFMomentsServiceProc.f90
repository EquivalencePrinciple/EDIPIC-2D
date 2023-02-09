
!--------------------------------------------------------
!
SUBROUTINE SYNCHRONIZE_MOMENTS_IN_OVERLAP_NODES

  USE ParallelOperationValues
  USE ClusterAndItsBoundaries
  USE Snapshots

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER n1  ! number of nodes in the y-direction
  INTEGER n3  ! number of nodes in the x-direction

  INTEGER ii, jj
  INTEGER i, j

  INTEGER pos1, pos2

  INTEGER bufsize
  INTEGER ALLOC_ERR
  REAL, ALLOCATABLE :: rbufer(:)

  n1 = c_indx_y_max_ext - c_indx_y_min_ext + 1
  n3 = c_indx_x_max_ext - c_indx_x_min_ext + 1

  IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! special case of self-connected X-periodic cluster
     DO j = c_indx_y_min_ext, c_indx_y_max_ext

        DO ii = 1, 1+ext_overlap
           cs_N(c_indx_x_min+ii, j) = cs_N(c_indx_x_min+ii, j) + cs_N(c_indx_x_max+ii-1, j)
           cs_N(c_indx_x_max-ii, j) = cs_N(c_indx_x_max-ii, j) + cs_N(c_indx_x_min-ii+1, j)
        END DO

!        cs_N(c_indx_x_min+1, j) = cs_N(c_indx_x_min+1, j) + cs_N(c_indx_x_max, j) 
!        cs_N(c_indx_x_max-1, j) = cs_N(c_indx_x_max-1, j) + cs_N(c_indx_x_min, j)

!### maybe a bug below, should've waited till vertical exchange is completed
!        cs_N(c_indx_x_min, j) = cs_N(c_indx_x_max-1, j)
!        cs_N(c_indx_x_max, j) = cs_N(c_indx_x_min+1, j)

     END DO

     DO j = c_indx_y_min_ext, c_indx_y_max_ext
        DO ii = 1, 1+ext_overlap
           cs_VX(c_indx_x_min+ii, j) = cs_VX(c_indx_x_min+ii, j) + cs_VX(c_indx_x_max+ii-1, j)
           cs_VX(c_indx_x_max-ii, j) = cs_VX(c_indx_x_max-ii, j) + cs_VX(c_indx_x_min-ii+1, j)
        END DO
     END DO

     DO j = c_indx_y_min_ext, c_indx_y_max_ext
        DO ii = 1, 1+ext_overlap
           cs_VY(c_indx_x_min+ii, j) = cs_VY(c_indx_x_min+ii, j) + cs_VY(c_indx_x_max+ii-1, j)
           cs_VY(c_indx_x_max-ii, j) = cs_VY(c_indx_x_max-ii, j) + cs_VY(c_indx_x_min-ii+1, j)
        END DO
     END DO

     DO j = c_indx_y_min_ext, c_indx_y_max_ext
        DO ii = 1, 1+ext_overlap
           cs_VZ(c_indx_x_min+ii, j) = cs_VZ(c_indx_x_min+ii, j) + cs_VZ(c_indx_x_max+ii-1, j)
           cs_VZ(c_indx_x_max-ii, j) = cs_VZ(c_indx_x_max-ii, j) + cs_VZ(c_indx_x_min-ii+1, j)
        END DO
     END DO

     DO j = c_indx_y_min_ext, c_indx_y_max_ext
        DO ii = 1, 1+ext_overlap
           cs_WX(c_indx_x_min+ii, j) = cs_WX(c_indx_x_min+ii, j) + cs_WX(c_indx_x_max+ii-1, j)
           cs_WX(c_indx_x_max-ii, j) = cs_WX(c_indx_x_max-ii, j) + cs_WX(c_indx_x_min-ii+1, j)
        END DO
     END DO

     DO j = c_indx_y_min_ext, c_indx_y_max_ext
        DO ii = 1, 1+ext_overlap
           cs_WY(c_indx_x_min+ii, j) = cs_WY(c_indx_x_min+ii, j) + cs_WY(c_indx_x_max+ii-1, j)
           cs_WY(c_indx_x_max-ii, j) = cs_WY(c_indx_x_max-ii, j) + cs_WY(c_indx_x_min-ii+1, j)
        END DO
     END DO

     DO j = c_indx_y_min_ext, c_indx_y_max_ext
        DO ii = 1, 1+ext_overlap
           cs_WZ(c_indx_x_min+ii, j) = cs_WZ(c_indx_x_min+ii, j) + cs_WZ(c_indx_x_max+ii-1, j)
           cs_WZ(c_indx_x_max-ii, j) = cs_WZ(c_indx_x_max-ii, j) + cs_WZ(c_indx_x_min-ii+1, j)
        END DO
     END DO

     DO j = c_indx_y_min_ext, c_indx_y_max_ext
        DO ii = 1, 1+ext_overlap
           cs_VXVY(c_indx_x_min+ii, j) = cs_VXVY(c_indx_x_min+ii, j) + cs_VXVY(c_indx_x_max+ii-1, j)
           cs_VXVY(c_indx_x_max-ii, j) = cs_VXVY(c_indx_x_max-ii, j) + cs_VXVY(c_indx_x_min-ii+1, j)
        END DO
     END DO

     DO j = c_indx_y_min_ext, c_indx_y_max_ext
        DO ii = 1, 1+ext_overlap
           cs_VXVZ(c_indx_x_min+ii, j) = cs_VXVZ(c_indx_x_min+ii, j) + cs_VXVZ(c_indx_x_max+ii-1, j)
           cs_VXVZ(c_indx_x_max-ii, j) = cs_VXVZ(c_indx_x_max-ii, j) + cs_VXVZ(c_indx_x_min-ii+1, j)
        END DO
     END DO

     DO j = c_indx_y_min_ext, c_indx_y_max_ext
        DO ii = 1, 1+ext_overlap
           cs_VYVZ(c_indx_x_min+ii, j) = cs_VYVZ(c_indx_x_min+ii, j) + cs_VYVZ(c_indx_x_max+ii-1, j)
           cs_VYVZ(c_indx_x_max-ii, j) = cs_VYVZ(c_indx_x_max-ii, j) + cs_VYVZ(c_indx_x_min-ii+1, j)
        END DO
     END DO

     DO j = c_indx_y_min_ext, c_indx_y_max_ext
        DO ii = 1, 1+ext_overlap
           cs_QX(c_indx_x_min+ii, j) = cs_QX(c_indx_x_min+ii, j) + cs_QX(c_indx_x_max+ii-1, j)
           cs_QX(c_indx_x_max-ii, j) = cs_QX(c_indx_x_max-ii, j) + cs_QX(c_indx_x_min-ii+1, j)
        END DO
     END DO

     DO j = c_indx_y_min_ext, c_indx_y_max_ext
        DO ii = 1, 1+ext_overlap
           cs_QY(c_indx_x_min+ii, j) = cs_QY(c_indx_x_min+ii, j) + cs_QY(c_indx_x_max+ii-1, j)
           cs_QY(c_indx_x_max-ii, j) = cs_QY(c_indx_x_max-ii, j) + cs_QY(c_indx_x_min-ii+1, j)
        END DO
     END DO

     DO j = c_indx_y_min_ext, c_indx_y_max_ext
        DO ii = 1, 1+ext_overlap
           cs_QZ(c_indx_x_min+ii, j) = cs_QZ(c_indx_x_min+ii, j) + cs_QZ(c_indx_x_max+ii-1, j)
           cs_QZ(c_indx_x_max-ii, j) = cs_QZ(c_indx_x_max-ii, j) + cs_QZ(c_indx_x_min-ii+1, j)
        END DO
     END DO

  END IF

! include neighbor contributions in nodes which are 1:1+ext_overlap lines away from the cluster boundary

  IF (WHITE_CLUSTER) THEN  
! "white processes"

     bufsize = 13 * n1 * (ext_overlap + 1)

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

     IF (Rank_horizontal_right.GE.0) THEN
! ## 1 ## send right moments in the right edge
        pos2=0
        DO ii = 0, ext_overlap
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_N(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VX(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VY(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VZ(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_WX(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_WY(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_WZ(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VXVY(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_QX(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_QY(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_QZ(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
        END DO

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_left.GE.0) THEN
! ## 2 ## send left moments in the left edge
        pos2=0
        DO ii = ext_overlap, 0, -1
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_N(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VX(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VY(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VZ(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_WX(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_WY(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_WZ(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_QX(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_QY(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_QZ(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
        END DO

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_left.GE.0) THEN
! ## 3 ## receive from left moments in the vertical line next to the left edge
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)

        pos2=0
        DO ii = 1, 1+ext_overlap
           pos1=pos2+1
           pos2=pos2+n1
           cs_N(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_N(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_VX(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VX(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_WX(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_WX(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_WY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_WY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_WZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_WZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_VXVY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VXVY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VXVZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VXVZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VYVZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VYVZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_QX(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_QX(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_QY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_QY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_QZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_QZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
        END DO
     END IF

     IF (Rank_horizontal_right.GE.0) THEN
! ## 4 ## receive from right moments in the vertical line next to the right edge
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
        pos2=0
        DO ii = ext_overlap+1, 1, -1
           pos1=pos2+1
           pos2=pos2+n1
           cs_N(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_N(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_VX(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VX(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_WX(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_WX(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_WY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_WY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_WZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_WZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_VXVY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VXVY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VXVZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VXVZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VYVZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VYVZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_QX(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_QX(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_QY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_QY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_QZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_QZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
        END DO
     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     bufsize = 13 * n3 * (ext_overlap + 1)
     ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

     IF (Rank_horizontal_above.GE.0) THEN
! ## 5 ## send up moments in the top edge
        pos2 = 0
        DO jj = 0, ext_overlap
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_N(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_WX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_WY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_WZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_QX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_QY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_QZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
        END DO

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_below.GE.0) THEN
! ## 6 ## send down moments in the bottom edge
        pos2 = 0
        DO jj = ext_overlap, 0, -1
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_N(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_WX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_WY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_WZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_QX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_QY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_QZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
        END DO

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_below.GE.0) THEN
! ## 7 ## receive from below moments in the vertical line above the bottom line
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)

        pos2 = 0
        DO jj = 1, 1+ext_overlap
           pos1=pos2+1
           pos2=pos2+n3
           cs_N(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_N(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_VX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_VX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_VY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_VZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_WX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_WX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_WY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_WY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_WZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_WZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_VXVY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_VXVY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VXVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_VXVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VYVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_VYVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_QX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_QX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_QY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_QY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_QZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_QZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
        END DO
     END IF

     IF (Rank_horizontal_above.GE.0) THEN
! ## 8 ## receive from above moments in the vertical line under the top line
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)

        pos2 = 0
        DO jj = 1+ext_overlap, 1, -1
           pos1=pos2+1
           pos2=pos2+n3
           cs_N(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_N(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_VX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_VX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_VY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_VZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_WX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_WX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_WY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_WY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_WZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_WZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_VXVY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_VXVY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VXVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_VXVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VYVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_VYVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_QX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_QX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_QY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_QY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_QZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_QZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
        END DO
     END IF

  ELSE
! "black" processes

     bufsize = 13 * n1 * (ext_overlap + 1)

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

     IF (Rank_horizontal_left.GE.0) THEN
! ## 1 ## receive from left moments in the vertical line next to the left edge
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)

        pos2=0
        DO ii = 1, 1+ext_overlap
           pos1=pos2+1
           pos2=pos2+n1
           cs_N(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_N(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_VX(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VX(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_WX(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_WX(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_WY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_WY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_WZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_WZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_VXVY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VXVY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VXVZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VXVZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VYVZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VYVZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_QX(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_QX(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_QY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_QY(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_QZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_QZ(c_indx_x_min+ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
        END DO
     END IF

     IF (Rank_horizontal_right.GE.0) THEN
! ## 2 ## receive from right moments in the vertical line next to the right edge
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
        pos2=0
        DO ii = ext_overlap+1, 1, -1
           pos1=pos2+1
           pos2=pos2+n1
           cs_N(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_N(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_VX(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VX(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_WX(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_WX(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_WY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_WY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_WZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_WZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_VXVY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VXVY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VXVZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VXVZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_VYVZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_VYVZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n1
           cs_QX(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_QX(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_QY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_QY(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n1
           cs_QZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) = cs_QZ(c_indx_x_max-ii, c_indx_y_min_ext:c_indx_y_max_ext) + rbufer(pos1:pos2)
        END DO
      END IF

     IF (Rank_horizontal_right.GE.0) THEN
! ## 3 ## send right moments in the right edge
        pos2=0
        DO ii = 0, ext_overlap
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_N(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VX(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VY(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VZ(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_WX(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_WY(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_WZ(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VXVY(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_QX(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_QY(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_QZ(c_indx_x_max+ii, c_indx_y_min_ext:c_indx_y_max_ext)
        END DO

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_left.GE.0) THEN
! ## 4 ## send left moments in the left edge
        pos2=0
        DO ii = ext_overlap, 0, -1
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_N(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VX(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VY(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VZ(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_WX(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_WY(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_WZ(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)

           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_QX(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_QY(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
           pos1=pos2+1
           pos2=pos2+n1
           rbufer(pos1:pos2) = cs_QZ(c_indx_x_min-ii, c_indx_y_min_ext:c_indx_y_max_ext)
        END DO

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     bufsize = 13 * n3 * (ext_overlap + 1)
     ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

     IF (Rank_horizontal_below.GE.0) THEN
! ## 5 ## receive from below moments in the vertical line above the bottom line
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)

        pos2 = 0
        DO jj = 1, 1+ext_overlap
           pos1=pos2+1
           pos2=pos2+n3
           cs_N(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_N(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_VX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_VX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_VY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_VZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_WX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_WX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_WY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_WY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_WZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_WZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_VXVY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_VXVY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VXVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_VXVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VYVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_VYVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_QX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_QX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_QY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_QY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_QZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) = cs_QZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min+jj) + rbufer(pos1:pos2)
        END DO
     END IF

     IF (Rank_horizontal_above.GE.0) THEN
! ## 6 ## receive from above moments in the vertical line under the top line
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)

        pos2 = 0
        DO jj = 1+ext_overlap, 1, -1
           pos1=pos2+1
           pos2=pos2+n3
           cs_N(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_N(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_VX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_VX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_VY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_VZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_WX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_WX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_WY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_WY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_WZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_WZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_VXVY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_VXVY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VXVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_VXVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_VYVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_VYVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)

           pos1=pos2+1
           pos2=pos2+n3
           cs_QX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_QX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_QY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_QY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
           pos1=pos2+1
           pos2=pos2+n3
           cs_QZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) = cs_QZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max-jj) + rbufer(pos1:pos2)
        END DO
     END IF

     IF (Rank_horizontal_above.GE.0) THEN
! ## 7 ## send up moments in the top edge
        pos2 = 0
        DO jj = 0, ext_overlap
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_N(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_WX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_WY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_WZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_QX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_QY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_QZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_max+jj)
        END DO

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
      END IF

     IF (Rank_horizontal_below.GE.0) THEN
! ## 8 ## send down moments in the bottom edge
        pos2 = 0
        DO jj = ext_overlap, 0, -1
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_N(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_WX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_WY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_WZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)

           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_QX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_QY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
           pos1=pos2+1
           pos2=pos2+n3
           rbufer(pos1:pos2) = cs_QZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min-jj)
        END DO

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

  END IF   !###   IF (WHITE_CLUSTER) THEN

! set values at cluster boundaries equal to those of the inner nodes in neighbors if they exist (to get continuous coverage for the moments)
! note that while arrays cs_N, cs_VX/Y/Z, etc have extended index limits, we only set values at the ordinary [not extended] cluster boundary

  n1 = c_indx_y_max - c_indx_y_min + 1
  n3 = c_indx_x_max - c_indx_x_min + 1

  IF (WHITE_CLUSTER) THEN  
! "white processes"

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     bufsize = 13 * n1
     ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

     IF (Rank_horizontal_right.GE.0) THEN
! ## 1 ## send right moments in the right edge MINUS ONE
        pos1=1
        pos2=n1
        rbufer(pos1:pos2) = cs_N(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VX(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WX(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QX(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_left.GE.0) THEN
! ## 2 ## send left moments in the left edge PLUS ONE
        pos1=1
        pos2=n1
        rbufer(pos1:pos2) = cs_N(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VX(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WX(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QX(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_left.GE.0) THEN
! ## 3 ## receive from left moments in the vertical line AT the left edge
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)

        pos2 = 0

        pos1=pos2+1
        pos2=pos2+n1
        cs_N(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_VX(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VY(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VZ(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_WX(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_WY(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_WZ(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_VXVY(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VXVZ(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VYVZ(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_QX(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_QY(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_QZ(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
     END IF
     
     IF (Rank_horizontal_right.GE.0) THEN
! ## 4 ## receive from right moments in the vertical line AT the right edge
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)

        pos2 = 0

        pos1=pos2+1
        pos2=pos2+n1
        cs_N(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_VX(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VY(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VZ(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_WX(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_WY(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_WZ(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_VXVY(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VXVZ(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VYVZ(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_QX(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_QY(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_QZ(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
    END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     bufsize = 13 * n3
     ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

     IF (Rank_horizontal_above.GE.0) THEN
! ## 5 ## send up moments in the top edge MINUS ONE
        pos1=1
        pos2=n3
        rbufer(pos1:pos2) = cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_below.GE.0) THEN
! ## 6 ## send down moments in the bottom edge PLUS ONE
        pos1=1
        pos2=n3
        rbufer(pos1:pos2) = cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_below.GE.0) THEN
! ## 7 ## receive from below moments in the horizontal line AT the bottom
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)

        pos2 = 0

        pos1=pos2+1
        pos2=pos2+n3
        cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
     END IF

     IF (Rank_horizontal_above.GE.0) THEN
! ## 8 ## receive from above moments in the horizontal line AT the top
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)

        pos2 = 0

        pos1=pos2+1
        pos2=pos2+n3
        cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
     END IF
     
  ELSE
! "black" processes

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     bufsize = 13 * n1
     ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

     IF (Rank_horizontal_left.GE.0) THEN
! ## 1 ## receive from left moments in the vertical line AT the left edge
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)

        pos2 = 0

        pos1=pos2+1
        pos2=pos2+n1
        cs_N(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_VX(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VY(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VZ(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_WX(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_WY(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_WZ(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_VXVY(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VXVZ(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VYVZ(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_QX(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_QY(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_QZ(c_indx_x_min, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
     END IF

     IF (Rank_horizontal_right.GE.0) THEN
! ## 2 ## receive from right moments in the vertical line AT the right edge
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)

        pos2 = 0

        pos1=pos2+1
        pos2=pos2+n1
        cs_N(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_VX(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VY(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VZ(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_WX(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_WY(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_WZ(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_VXVY(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VXVZ(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_VYVZ(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n1
        cs_QX(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_QY(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n1
        cs_QZ(c_indx_x_max, c_indx_y_min:c_indx_y_max) = rbufer(pos1:pos2)
     END IF

     IF (Rank_horizontal_right.GE.0) THEN
! ## 3 ## send right moments in the right edge MINUS ONE
        pos1=1
        pos2=n1
        rbufer(pos1:pos2) = cs_N(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VX(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WX(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QX(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QY(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_left.GE.0) THEN
! ## 4 ## send left moments in the left edge PLUS ONE
        pos1=1
        pos2=n1
        rbufer(pos1:pos2) = cs_N(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VX(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WX(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QX(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QY(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
        pos1=pos2+1
        pos2=pos2+n1
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     bufsize = 13 * n3
     ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)

     IF (Rank_horizontal_below.GE.0) THEN
! ## 5 ## receive from below moments in the horizontal line AT the bottom
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)

        pos2 = 0

        pos1=pos2+1
        pos2=pos2+n3
        cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_min) = rbufer(pos1:pos2)
     END IF
     
     IF (Rank_horizontal_above.GE.0) THEN
! ## 6 ## receive from above moments in the horizontal line AT the top
        CALL MPI_RECV(rbufer, bufsize, MPI_REAL, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)

        pos2 = 0

        pos1=pos2+1
        pos2=pos2+n3
        cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)

        pos1=pos2+1
        pos2=pos2+n3
        cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
        pos1=pos2+1
        pos2=pos2+n3
        cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_max) = rbufer(pos1:pos2)
     END IF
     
     IF (Rank_horizontal_above.GE.0) THEN
! ## 7 ## send up moments in the top edge MINUS ONE
        pos1=1
        pos2=n3
        rbufer(pos1:pos2) = cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF

     IF (Rank_horizontal_below.GE.0) THEN
! ## 8 ## send down moments in the bottom edge PLUS ONE
        pos1=1
        pos2=n3
        rbufer(pos1:pos2) = cs_N(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VX(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WX(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_WZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VXVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_VYVZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QX(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QY(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
        pos1=pos2+1
        pos2=pos2+n3
        rbufer(pos1:pos2) = cs_QZ(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)

        CALL MPI_SEND(rbufer, bufsize, MPI_REAL, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
     END IF
     
  END IF

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer,  STAT=ALLOC_ERR)

END SUBROUTINE SYNCHRONIZE_MOMENTS_IN_OVERLAP_NODES

!----------------------------------------------------
!
SUBROUTINE ADJUST_DENSITY_AT_WALL_BOUNDARIES

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE ClusterAndItsBoundaries
!  USE BlockAndItsBoundaries, ONLY : indx_x_min, indx_x_max, indx_y_min, indx_y_max
  USE Snapshots

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER n1  ! number of nodes in the y-direction
  INTEGER n2  !
  INTEGER n3  ! number of nodes in the x-direction

  REAL, ALLOCATABLE :: rbufer_n(:)
  REAL, ALLOCATABLE :: rbufer_vx(:)
  REAL, ALLOCATABLE :: rbufer_vy(:)
  REAL, ALLOCATABLE :: rbufer_vz(:)
  REAL, ALLOCATABLE :: rbufer_vx2(:)
  REAL, ALLOCATABLE :: rbufer_vy2(:)
  REAL, ALLOCATABLE :: rbufer_vz2(:)
  INTEGER ALLOC_ERR
  INTEGER bufsize

  INTEGER i, j, k
  INTEGER pos_i_j, pos_ip1_j, pos_i_jp1, pos_ip1_jp1
  REAL ax_ip1, ax_i, ay_jp1, ay_j
  REAL vij, vip1j, vijp1, vip1jp1
  REAL v_ij, v_ip1j, v_ijp1, v_ip1jp1

  INTEGER pos1, pos2

  REAL, ALLOCATABLE :: rbufer(:)

    IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! special case of self-connected X-periodic cluster
     IF (Rank_of_master_above.LT.0) THEN
        DO i = c_indx_x_min, c_indx_x_max
           cs_N(i, c_indx_y_max) = 2.0 * cs_N(i, c_indx_y_max)
        END DO
     END IF
  
     IF (Rank_of_master_below.LT.0) THEN
        DO i = c_indx_x_min, c_indx_x_max
           cs_N(i, c_indx_y_min) = 2.0 * cs_N(i, c_indx_y_min)
        END DO
     END IF

  ELSE
 
     IF (Rank_of_master_left.LT.0) THEN
        DO j = c_indx_y_min+1, c_indx_y_max-1
           cs_N(c_indx_x_min, j) = 2.0 * cs_N(c_indx_x_min, j)
        END DO
     END IF

     IF (Rank_of_master_right.LT.0) THEN
        DO j = c_indx_y_min+1, c_indx_y_max-1
           cs_N(c_indx_x_max, j) = 2.0 * cs_N(c_indx_x_max, j)
        END DO
     END IF
  
     IF (Rank_of_master_above.LT.0) THEN
        DO i = c_indx_x_min+1, c_indx_x_max-1
           cs_N(i, c_indx_y_max) = 2.0 * cs_N(i, c_indx_y_max)
        END DO
     END IF
  
     IF (Rank_of_master_below.LT.0) THEN
        DO i = c_indx_x_min+1, c_indx_x_max-1
           cs_N(i, c_indx_y_min) = 2.0 * cs_N(i, c_indx_y_min)
        END DO
     END IF

     SELECT CASE (c_left_bottom_corner_type)
        CASE (SURROUNDED_BY_WALL)
           cs_N(c_indx_x_min, c_indx_y_min) = 4.0 * cs_N(c_indx_x_min, c_indx_y_min) 
        CASE (FLAT_WALL_LEFT)
           cs_N(c_indx_x_min, c_indx_y_min) = 2.0 * cs_N(c_indx_x_min, c_indx_y_min) 
        CASE (FLAT_WALL_BELOW)
           cs_N(c_indx_x_min, c_indx_y_min) = 2.0 * cs_N(c_indx_x_min, c_indx_y_min) 
        CASE (EMPTY_CORNER_WALL_LEFT)
           cs_N(c_indx_x_min, c_indx_y_min+1) = 0.66666666666 * cs_N(c_indx_x_min, c_indx_y_min+1)    ! 2*2/3=4/3=1/(3/4)
        CASE (EMPTY_CORNER_WALL_BELOW)
           cs_N(c_indx_x_min+1, c_indx_y_min) = 0.66666666666 * cs_N(c_indx_x_min+1, c_indx_y_min)    ! 2*2/3=4/3=1/(3/4)
     END SELECT
     
     SELECT CASE (c_left_top_corner_type)
        CASE (SURROUNDED_BY_WALL)
           cs_N(c_indx_x_min, c_indx_y_max) = 4.0 * cs_N(c_indx_x_min, c_indx_y_max) 
        CASE (FLAT_WALL_LEFT)
           cs_N(c_indx_x_min, c_indx_y_max) = 2.0 * cs_N(c_indx_x_min, c_indx_y_max) 
        CASE (FLAT_WALL_ABOVE)
           cs_N(c_indx_x_min, c_indx_y_max) = 2.0 * cs_N(c_indx_x_min, c_indx_y_max) 
        CASE (EMPTY_CORNER_WALL_LEFT)
           cs_N(c_indx_x_min, c_indx_y_max-1) = 0.66666666666 * cs_N(c_indx_x_min, c_indx_y_max-1)    ! 2*2/3=4/3=1/(3/4)
        CASE (EMPTY_CORNER_WALL_ABOVE)
           cs_N(c_indx_x_min+1, c_indx_y_max) = 0.66666666666 * cs_N(c_indx_x_min+1, c_indx_y_max)    ! 2*2/3=4/3=1/(3/4)
     END SELECT

     SELECT CASE (c_right_bottom_corner_type)
        CASE (SURROUNDED_BY_WALL)
           cs_N(c_indx_x_max, c_indx_y_min) = 4.0 * cs_N(c_indx_x_max, c_indx_y_min) 
        CASE (FLAT_WALL_RIGHT)
           cs_N(c_indx_x_max, c_indx_y_min) = 2.0 * cs_N(c_indx_x_max, c_indx_y_min) 
        CASE (FLAT_WALL_BELOW)
           cs_N(c_indx_x_max, c_indx_y_min) = 2.0 * cs_N(c_indx_x_max, c_indx_y_min) 
        CASE (EMPTY_CORNER_WALL_RIGHT)
           cs_N(c_indx_x_max, c_indx_y_min+1) = 0.66666666666 * cs_N(c_indx_x_max, c_indx_y_min+1)    ! 2*2/3=4/3=1/(3/4)
        CASE (EMPTY_CORNER_WALL_BELOW)
           cs_N(c_indx_x_max-1, c_indx_y_min) = 0.66666666666 * cs_N(c_indx_x_max-1, c_indx_y_min)    ! 2*2/3=4/3=1/(3/4)
     END SELECT
     
     SELECT CASE (c_right_top_corner_type)
        CASE (SURROUNDED_BY_WALL)
           cs_N(c_indx_x_max, c_indx_y_max) = 4.0 * cs_N(c_indx_x_max, c_indx_y_max) 
        CASE (FLAT_WALL_RIGHT)
           cs_N(c_indx_x_max, c_indx_y_max) = 2.0 * cs_N(c_indx_x_max, c_indx_y_max) 
        CASE (FLAT_WALL_ABOVE)
           cs_N(c_indx_x_max, c_indx_y_max) = 2.0 * cs_N(c_indx_x_max, c_indx_y_max) 
        CASE (EMPTY_CORNER_WALL_RIGHT)
           cs_N(c_indx_x_max, c_indx_y_max-1) = 0.66666666666 * cs_N(c_indx_x_max, c_indx_y_max-1)    ! 2*2/3=4/3=1/(3/4)
        CASE (EMPTY_CORNER_WALL_ABOVE)
           cs_N(c_indx_x_max-1, c_indx_y_max) = 0.66666666666 * cs_N(c_indx_x_max-1, c_indx_y_max)    ! 2*2/3=4/3=1/(3/4)
     END SELECT

  END IF

END SUBROUTINE ADJUST_DENSITY_AT_WALL_BOUNDARIES
