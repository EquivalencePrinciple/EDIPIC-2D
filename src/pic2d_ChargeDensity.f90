
!------------------------------
!
SUBROUTINE GATHER_STREAMING_CHARGE_DENSITY

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE IonParticles
  USE ClusterAndItsBoundaries
  USE BlockAndItsBoundaries, ONLY : indx_x_min, indx_x_max, indx_y_min, indx_y_max

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request


  INTEGER n1  ! number of nodes in the y-direction
  INTEGER n2  !
  INTEGER n3  ! number of nodes in the x-direction

  REAL(8), ALLOCATABLE :: rbufer(:)
  REAL(8), ALLOCATABLE :: rbufer2(:)
  INTEGER ALLOC_ERR
  INTEGER bufsize

  REAL(8), ALLOCATABLE :: c_rho(:,:)   ! these arrays are used by cluster masters only
  REAL(8), ALLOCATABLE :: c_X11(:,:)   ! they are not necessary ooutside this procedure (AFAIK now)
  REAL(8), ALLOCATABLE :: c_X12(:,:)
  REAL(8), ALLOCATABLE :: c_X21(:,:)
  REAL(8), ALLOCATABLE :: c_X22(:,:)

  INTEGER i, j, k, s
  INTEGER pos_i_j, pos_ip1_j, pos_i_jp1, pos_ip1_jp1
  REAL(8) ax_ip1, ax_i, ay_jp1, ay_j
  REAL(8) vij, vip1j, vijp1

  INTEGER pos

  REAL(8) alfa_x, alfa_y, alfa_z
  REAL(8) alfa_x2, alfa_y2, alfa_z2
  REAL(8) theta2, invtheta
  REAL(8) K11, K12, K21, K22
  REAL(8) A11, A12, A21, A22

  INTEGER shift
  INTEGER shift_X11, shift_X12, shift_X21, shift_X22

  INTEGER nio, position_flag

! functions
  REAL(8) Get_Surface_Charge_Inner_Object
  REAL(8) Bx, By, Bz

  n1 = c_indx_y_max - c_indx_y_min + 1
! was  n2 = -c_indx_x_min + 1 - c_indx_y_min * n1
  n3 = c_indx_x_max - c_indx_x_min + 1
  n2 = -c_indx_x_min + 1 - c_indx_y_min * n3

  bufsize = n1 * n3
  ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer2(1:bufsize), STAT=ALLOC_ERR)

  rbufer = 0.0_8

  DO k = 1, N_electrons
     
     i = INT(electron(k)%X)
     j = INT(electron(k)%Y)

     IF (electron(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
     IF (electron(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

     if ((i.lt.c_indx_x_min).or.(i.gt.(c_indx_x_max-1)).or.(j.lt.c_indx_y_min).or.(j.gt.(c_indx_y_max-1))) then
        print '("Process ",i4," : Error-1 in GATHER_ELECTRON_CHARGE_DENSITY : index out of bounds")', Rank_of_process
        print '("Process ",i4," : k/N_electrons : ",i8,2x,i8)', Rank_of_process, k, N_electrons
        print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag
        print '("Process ",i4," : minx/maxx/miny/maxy : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     end if

!     pos = i - c_indx_x_min + 1 + (j - c_indx_y_min) * (c_indx_x_max - c_indx_x_min + 1)

     pos_i_j     = i + j * n3 + n2
     pos_ip1_j   = pos_i_j + 1
     pos_i_jp1   = pos_i_j + n3
     pos_ip1_jp1 = pos_i_jp1 + 1

     ax_ip1 = electron(k)%X - DBLE(i)
     ax_i   = 1.0_8 - ax_ip1

     ay_jp1 = electron(k)%Y - DBLE(j)
     ay_j = 1.0_8 - ay_jp1

     vij   = ax_i   * ay_j
     vip1j = ax_ip1 * ay_j
     vijp1 = ax_i   * ay_jp1

     if ((pos_i_j.gt.bufsize)) then
        write (*,*) "OOPS"     
        print '(2x,8(2x,i10))', Rank_of_process, bufsize, pos_i_j, i, j, k, n1, n2
     end if

     rbufer(pos_i_j)     = rbufer(pos_i_j)     + vij                         !ax_i   * ay_j
     rbufer(pos_ip1_j)   = rbufer(pos_ip1_j)   + vip1j                       !ax_ip1 * ay_j
     rbufer(pos_i_jp1)   = rbufer(pos_i_jp1)   + vijp1                       !ax_i   * ay_jp1
     rbufer(pos_ip1_jp1) = rbufer(pos_ip1_jp1) + 1.0_8 - vij - vip1j - vijp1 !ax_ip1 * ay_jp1

  END DO

! collect densities from all processes in a cluster
  rbufer2 = 0.0_8
  CALL MPI_REDUCE(rbufer, rbufer2, bufsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_CLUSTER, ierr)

  IF (cluster_rank_key.EQ.0) THEN

     ALLOCATE (c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE (c_X11(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE (c_X12(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE (c_X21(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
     ALLOCATE (c_X22(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     pos = 0
     DO j = c_indx_y_min, c_indx_y_max
        DO i = c_indx_x_min, c_indx_x_max
           pos = pos+1

           c_rho(i,j) = -rbufer2(pos)   ! minus because it's electrons

           alfa_x = -0.5_8 * Bx(DBLE(i), DBLE(j))
           alfa_y = -0.5_8 * By(DBLE(i), DBLE(j))
           alfa_z = -0.5_8 * Bz(DBLE(i), DBLE(j))

           alfa_x2 = alfa_x**2
           alfa_y2 = alfa_y**2
           alfa_z2 = alfa_z**2

           theta2 = alfa_x2 + alfa_y2 + alfa_z2
           invtheta = 1.0_8 / (1.0_8 + theta2)
    
! matrix K, same as R in Gibbons & Hewett; A_inverse = (I + K)/2 = (I + R)/2
           K11 =  (1.0_8 - theta2 + 2.0_8 * alfa_x2) * invtheta
           K12 =  2.0_8 * (alfa_x * alfa_y + alfa_z) * invtheta

           K21 =  2.0_8 * (alfa_x * alfa_y - alfa_z) * invtheta
           K22 =  (1.0_8 - theta2 + 2.0_8 * alfa_y2) * invtheta

           A11 = 0.5_8 * (K11 + 1.0_8)
           A12 = 0.5_8 * K12

           A21 = 0.5_8 * K21
           A22 = 0.5_8 * (K22 + 1.0_8)

           c_X11(i,j) = A11 * rbufer2(pos)
           c_X12(i,j) = A12 * rbufer2(pos)
           c_X21(i,j) = A21 * rbufer2(pos)
           c_X22(i,j) = A22 * rbufer2(pos)
        END DO
     END DO

  END IF

! ions
  
  DO s = 1, N_spec

     rbufer = 0.0_8
     
     DO k = 1, N_ions(s)
     
        i = INT(ion(s)%part(k)%X)
        j = INT(ion(s)%part(k)%Y)
        ! check if the particle within the cluster to which the node belongs:
        IF (ion(s)%part(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
        IF (ion(s)%part(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

!     pos = i - c_indx_x_min + 1 + (j - c_indx_y_min) * (c_indx_x_max - c_indx_x_min + 1)

        pos_i_j     = i + j * n3 + n2
        pos_ip1_j   = pos_i_j + 1
        pos_i_jp1   = pos_i_j + n3
        pos_ip1_jp1 = pos_i_jp1 + 1

        ax_ip1 = ion(s)%part(k)%X - DBLE(i)
        ax_i   = 1.0_8 - ax_ip1

        ay_jp1 = ion(s)%part(k)%Y - DBLE(j)
        ay_j = 1.0_8 - ay_jp1

        vij   = ax_i   * ay_j
        vip1j = ax_ip1 * ay_j
        vijp1 = ax_i   * ay_jp1

        if ((pos_i_j.gt.bufsize)) then
           write (*,*) "OOPS"     
           print '(2x,8(2x,i10))', Rank_of_process, bufsize, pos_i_j, i, j, k, n1, n2
        end if
        
        rbufer(pos_i_j)     = rbufer(pos_i_j)     + vij                         !ax_i   * ay_j
        rbufer(pos_ip1_j)   = rbufer(pos_ip1_j)   + vip1j                       !ax_ip1 * ay_j
        rbufer(pos_i_jp1)   = rbufer(pos_i_jp1)   + vijp1                       !ax_i   * ay_jp1
        rbufer(pos_ip1_jp1) = rbufer(pos_ip1_jp1) + 1.0_8 - vij - vip1j - vijp1 !ax_ip1 * ay_jp1

     END DO

! collect densities from all processes in a cluster
     rbufer2 = 0.0_8
     CALL MPI_REDUCE(rbufer, rbufer2, bufsize, MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM_CLUSTER, ierr)

     IF (cluster_rank_key.EQ.0) THEN

        pos = 0
        DO j = c_indx_y_min, c_indx_y_max
           DO i = c_indx_x_min, c_indx_x_max
              pos = pos+1

              c_rho(i,j) = c_rho(i,j) + Qs(s) * rbufer2(pos)
        
              alfa_x = QM2s(s) * Bx(DBLE(i), DBLE(j))
              alfa_y = QM2s(s) * By(DBLE(i), DBLE(j))
              alfa_z = QM2s(s) * Bz(DBLE(i), DBLE(j))

              alfa_x2 = alfa_x**2
              alfa_y2 = alfa_y**2
              alfa_z2 = alfa_z**2

              theta2 = alfa_x2 + alfa_y2 + alfa_z2
              invtheta = 1.0_8 / (1.0_8 + theta2)

!    matrix K, same as R in Gibbons & Hewett; A_inverse = (I + K)/2 = (I + R)/2
              K11 =  (1.0_8 - theta2 + 2.0_8 * alfa_x2) * invtheta
              K12 =  2.0_8 * (alfa_x * alfa_y + alfa_z) * invtheta

              K21 =  2.0_8 * (alfa_x * alfa_y - alfa_z) * invtheta
              K22 =  (1.0_8 - theta2 + 2.0_8 * alfa_y2) * invtheta

              A11 = 0.5_8 * (K11 + 1.0_8)
              A12 = 0.5_8 * K12
           
              A21 = 0.5_8 * K21
              A22 = 0.5_8 * (K22 + 1.0_8)

              c_X11(i,j) = c_X11(i,j) + A11 * Qs(s)**2 * rbufer2(pos)
              c_X12(i,j) = c_X12(i,j) + A12 * Qs(s)**2 * rbufer2(pos)
              c_X21(i,j) = c_X21(i,j) + A21 * Qs(s)**2 * rbufer2(pos)
              c_X22(i,j) = c_X22(i,j) + A22 * Qs(s)**2 * rbufer2(pos)
           END DO
        END DO

     END IF   !###   IF (cluster_rank_key.EQ.0) THEN
           
  END DO   !###   DO s = 1, N_spec

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

! now cluster masters exchange information about densities in overlapping nodes
  IF (cluster_rank_key.EQ.0) THEN

     c_X11 = c_X11 * Chi_factor
     c_X12 = c_X12 * Chi_factor
     c_X21 = c_X21 * Chi_factor
     c_X22 = c_X22 * Chi_factor

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! special case of self-connected X-periodic cluster
        DO j = c_indx_y_min, c_indx_y_max
           c_rho(c_indx_x_min+1, j) = c_rho(c_indx_x_min+1, j) + c_rho(c_indx_x_max, j) 
           c_rho(c_indx_x_max-1, j) = c_rho(c_indx_x_max-1, j) + c_rho(c_indx_x_min, j)
           c_rho(c_indx_x_min, j) = c_rho(c_indx_x_max-1, j)
           c_rho(c_indx_x_max, j) = c_rho(c_indx_x_min+1, j)

           c_X11(c_indx_x_min+1, j) = c_X11(c_indx_x_min+1, j) + c_X11(c_indx_x_max, j) 
           c_X11(c_indx_x_max-1, j) = c_X11(c_indx_x_max-1, j) + c_X11(c_indx_x_min, j)
           c_X11(c_indx_x_min, j) = c_X11(c_indx_x_max-1, j)
           c_X11(c_indx_x_max, j) = c_X11(c_indx_x_min+1, j)
           
           c_X12(c_indx_x_min+1, j) = c_X12(c_indx_x_min+1, j) + c_X12(c_indx_x_max, j) 
           c_X12(c_indx_x_max-1, j) = c_X12(c_indx_x_max-1, j) + c_X12(c_indx_x_min, j)
           c_X12(c_indx_x_min, j) = c_X12(c_indx_x_max-1, j)
           c_X12(c_indx_x_max, j) = c_X12(c_indx_x_min+1, j)

           c_X21(c_indx_x_min+1, j) = c_X21(c_indx_x_min+1, j) + c_X21(c_indx_x_max, j) 
           c_X21(c_indx_x_max-1, j) = c_X21(c_indx_x_max-1, j) + c_X21(c_indx_x_min, j)
           c_X21(c_indx_x_min, j) = c_X21(c_indx_x_max-1, j)
           c_X21(c_indx_x_max, j) = c_X21(c_indx_x_min+1, j)

           c_X22(c_indx_x_min+1, j) = c_X22(c_indx_x_min+1, j) + c_X22(c_indx_x_max, j) 
           c_X22(c_indx_x_max-1, j) = c_X22(c_indx_x_max-1, j) + c_X22(c_indx_x_min, j)
           c_X22(c_indx_x_min, j) = c_X22(c_indx_x_max-1, j)
           c_X22(c_indx_x_max, j) = c_X22(c_indx_x_min+1, j)
        END DO
     END IF

! exchange 1, with summation:
     IF (WHITE_CLUSTER) THEN  
! "white processes"

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:5*n1), STAT=ALLOC_ERR)

        shift = -c_indx_y_min+1
        shift_X11 = shift + n1
        shift_X12 = shift_X11 + n1
        shift_X21 = shift_X12 + n1
        shift_X22 = shift_X21 + n1

        IF (Rank_horizontal_right.GE.0) THEN
! ## 1 ## send right values in the right edge
           rbufer(     1:  n1) = c_rho(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           rbufer(  n1+1:2*n1) = c_X11(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           rbufer(2*n1+1:3*n1) = c_X12(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           rbufer(3*n1+1:4*n1) = c_X21(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           rbufer(4*n1+1:5*n1) = c_X22(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 2 ## send left values in the left edge
           rbufer(     1:  n1) = c_rho(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           rbufer(  n1+1:2*n1) = c_X11(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           rbufer(2*n1+1:3*n1) = c_X12(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           rbufer(3*n1+1:4*n1) = c_X21(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           rbufer(4*n1+1:5*n1) = c_X22(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 3 ## receive from left values in the vertical line next to the left edge
           CALL MPI_RECV(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho(c_indx_x_min+1, j) = c_rho(c_indx_x_min+1, j) + rbufer(j+shift)
              c_X11(c_indx_x_min+1, j) = c_X11(c_indx_x_min+1, j) + rbufer(j+shift_X11)
              c_X12(c_indx_x_min+1, j) = c_X12(c_indx_x_min+1, j) + rbufer(j+shift_X12)
              c_X21(c_indx_x_min+1, j) = c_X21(c_indx_x_min+1, j) + rbufer(j+shift_X21)
              c_X22(c_indx_x_min+1, j) = c_X22(c_indx_x_min+1, j) + rbufer(j+shift_X22)
           END DO
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 4 ## receive from right values in the vertical line next to the right edge
           CALL MPI_RECV(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho(c_indx_x_max-1, j) = c_rho(c_indx_x_max-1, j) + rbufer(j+shift)
              c_X11(c_indx_x_max-1, j) = c_X11(c_indx_x_max-1, j) + rbufer(j+shift_X11)
              c_X12(c_indx_x_max-1, j) = c_X12(c_indx_x_max-1, j) + rbufer(j+shift_X12)
              c_X21(c_indx_x_max-1, j) = c_X21(c_indx_x_max-1, j) + rbufer(j+shift_X21)
              c_X22(c_indx_x_max-1, j) = c_X22(c_indx_x_max-1, j) + rbufer(j+shift_X22)
           END DO           
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:5*n3), STAT=ALLOC_ERR)

        shift = -c_indx_x_min+1
        shift_X11 = shift + n3
        shift_X12 = shift_X11 + n3
        shift_X21 = shift_X12 + n3
        shift_X22 = shift_X21 + n3

        IF (Rank_horizontal_above.GE.0) THEN
! ## 5 ## send up values in the top edge
           rbufer(     1:  n3) = c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           rbufer(  n3+1:2*n3) = c_X11(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           rbufer(2*n3+1:3*n3) = c_X12(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           rbufer(3*n3+1:4*n3) = c_X21(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           rbufer(4*n3+1:5*n3) = c_X22(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           CALL MPI_SEND(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 6 ## send down values in the bottom edge
           rbufer(     1:  n3) = c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           rbufer(  n3+1:2*n3) = c_X11(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           rbufer(2*n3+1:3*n3) = c_X12(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           rbufer(3*n3+1:4*n3) = c_X21(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           rbufer(4*n3+1:5*n3) = c_X22(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           CALL MPI_SEND(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 7 ## receive from below values in the horizontal line above the bottom line
           CALL MPI_RECV(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_min+1) = c_rho(i, c_indx_y_min+1) + rbufer(i+shift)
              c_X11(i, c_indx_y_min+1) = c_X11(i, c_indx_y_min+1) + rbufer(i+shift_X11)
              c_X12(i, c_indx_y_min+1) = c_X12(i, c_indx_y_min+1) + rbufer(i+shift_X12)
              c_X21(i, c_indx_y_min+1) = c_X21(i, c_indx_y_min+1) + rbufer(i+shift_X21)
              c_X22(i, c_indx_y_min+1) = c_X22(i, c_indx_y_min+1) + rbufer(i+shift_X22)
           END DO
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 8 ## receive from above values in the horizontal line under the top line
           CALL MPI_RECV(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_max-1) = c_rho(i, c_indx_y_max-1) + rbufer(i+shift)
              c_X11(i, c_indx_y_max-1) = c_X11(i, c_indx_y_max-1) + rbufer(i+shift_X11)
              c_X12(i, c_indx_y_max-1) = c_X12(i, c_indx_y_max-1) + rbufer(i+shift_X12)
              c_X21(i, c_indx_y_max-1) = c_X21(i, c_indx_y_max-1) + rbufer(i+shift_X21)
              c_X22(i, c_indx_y_max-1) = c_X22(i, c_indx_y_max-1) + rbufer(i+shift_X22)
           END DO
        END IF

     ELSE
! "black" processes

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:5*n1), STAT=ALLOC_ERR)

        shift = -c_indx_y_min+1
        shift_X11 = shift + n1
        shift_X12 = shift_X11 + n1
        shift_X21 = shift_X12 + n1
        shift_X22 = shift_X21 + n1

        IF (Rank_horizontal_left.GE.0) THEN
! ## 1 ## receive from left values in the vertical line next to the left edge
           CALL MPI_RECV(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho(c_indx_x_min+1, j) = c_rho(c_indx_x_min+1, j) + rbufer(j+shift)
              c_X11(c_indx_x_min+1, j) = c_X11(c_indx_x_min+1, j) + rbufer(j+shift_X11)
              c_X12(c_indx_x_min+1, j) = c_X12(c_indx_x_min+1, j) + rbufer(j+shift_X12)
              c_X21(c_indx_x_min+1, j) = c_X21(c_indx_x_min+1, j) + rbufer(j+shift_X21)
              c_X22(c_indx_x_min+1, j) = c_X22(c_indx_x_min+1, j) + rbufer(j+shift_X22)
           END DO
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 2 ## receive from right values in the vertical line next to the right edge
           CALL MPI_RECV(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho(c_indx_x_max-1, j) = c_rho(c_indx_x_max-1, j) + rbufer(j+shift)
              c_X11(c_indx_x_max-1, j) = c_X11(c_indx_x_max-1, j) + rbufer(j+shift_X11)
              c_X12(c_indx_x_max-1, j) = c_X12(c_indx_x_max-1, j) + rbufer(j+shift_X12)
              c_X21(c_indx_x_max-1, j) = c_X21(c_indx_x_max-1, j) + rbufer(j+shift_X21)
              c_X22(c_indx_x_max-1, j) = c_X22(c_indx_x_max-1, j) + rbufer(j+shift_X22)
           END DO           
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 3 ## send right values in the right edge
           rbufer(     1:  n1) = c_rho(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           rbufer(  n1+1:2*n1) = c_X11(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           rbufer(2*n1+1:3*n1) = c_X12(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           rbufer(3*n1+1:4*n1) = c_X21(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           rbufer(4*n1+1:5*n1) = c_X22(c_indx_x_max, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 4 ## send left values in the left edge
           rbufer(     1:  n1) = c_rho(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           rbufer(  n1+1:2*n1) = c_X11(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           rbufer(2*n1+1:3*n1) = c_X12(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           rbufer(3*n1+1:4*n1) = c_X21(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           rbufer(4*n1+1:5*n1) = c_X22(c_indx_x_min, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:5*n3), STAT=ALLOC_ERR)

        shift = -c_indx_x_min+1
        shift_X11 = shift + n3
        shift_X12 = shift_X11 + n3
        shift_X21 = shift_X12 + n3
        shift_X22 = shift_X21 + n3

        IF (Rank_horizontal_below.GE.0) THEN
! ## 5 ## receive from below values in the horizontal line above the bottom line
           CALL MPI_RECV(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_min+1) = c_rho(i, c_indx_y_min+1) + rbufer(i+shift)
              c_X11(i, c_indx_y_min+1) = c_X11(i, c_indx_y_min+1) + rbufer(i+shift_X11)
              c_X12(i, c_indx_y_min+1) = c_X12(i, c_indx_y_min+1) + rbufer(i+shift_X12)
              c_X21(i, c_indx_y_min+1) = c_X21(i, c_indx_y_min+1) + rbufer(i+shift_X21)
              c_X22(i, c_indx_y_min+1) = c_X22(i, c_indx_y_min+1) + rbufer(i+shift_X22)
           END DO
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 6 ## receive from above values in the horizontal line under the top line
           CALL MPI_RECV(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_max-1) = c_rho(i, c_indx_y_max-1) + rbufer(i+shift)
              c_X11(i, c_indx_y_max-1) = c_X11(i, c_indx_y_max-1) + rbufer(i+shift_X11)
              c_X12(i, c_indx_y_max-1) = c_X12(i, c_indx_y_max-1) + rbufer(i+shift_X12)
              c_X21(i, c_indx_y_max-1) = c_X21(i, c_indx_y_max-1) + rbufer(i+shift_X21)
              c_X22(i, c_indx_y_max-1) = c_X22(i, c_indx_y_max-1) + rbufer(i+shift_X22)
           END DO
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 7 ## send up values in the top edge
           rbufer(     1:  n3) = c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           rbufer(  n3+1:2*n3) = c_X11(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           rbufer(2*n3+1:3*n3) = c_X12(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           rbufer(3*n3+1:4*n3) = c_X21(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           rbufer(4*n3+1:5*n3) = c_X22(c_indx_x_min:c_indx_x_max, c_indx_y_max)
           CALL MPI_SEND(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 8 ## send down values in the bottom edge
           rbufer(     1:  n3) = c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           rbufer(  n3+1:2*n3) = c_X11(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           rbufer(2*n3+1:3*n3) = c_X12(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           rbufer(3*n3+1:4*n3) = c_X21(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           rbufer(4*n3+1:5*n3) = c_X22(c_indx_x_min:c_indx_x_max, c_indx_y_min)
           CALL MPI_SEND(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

     END IF !black/white selection
 
     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)

! set values at the edges of cluster domains
! in most cases we do not need the charge desnity at the edges, only the X** arrays
! but for simplicity we update/synchronize everything 

     IF (WHITE_CLUSTER) THEN
! "white processes"

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:5*n1), STAT=ALLOC_ERR)

        shift = -c_indx_y_min+1
        shift_X11 = shift + n1
        shift_X12 = shift_X11 + n1
        shift_X21 = shift_X12 + n1
        shift_X22 = shift_X21 + n1

        IF (Rank_horizontal_right.GE.0) THEN
! ## 09 ## send right values to the left of the right edge
           rbufer(     1:  n1) = c_rho(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           rbufer(  n1+1:2*n1) = c_X11(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           rbufer(2*n1+1:3*n1) = c_X12(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           rbufer(3*n1+1:4*n1) = c_X21(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           rbufer(4*n1+1:5*n1) = c_X22(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
            CALL MPI_SEND(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr)
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 10 ## send left values to the right of the left edge
           rbufer(     1:  n1) = c_rho(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           rbufer(  n1+1:2*n1) = c_X11(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           rbufer(2*n1+1:3*n1) = c_X12(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           rbufer(3*n1+1:4*n1) = c_X21(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           rbufer(4*n1+1:5*n1) = c_X22(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr)
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 11 ## receive from left values for the left edge 
           CALL MPI_RECV(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho(c_indx_x_min, j) = rbufer(j+shift)
              c_X11(c_indx_x_min, j) = rbufer(j+shift_X11)
              c_X12(c_indx_x_min, j) = rbufer(j+shift_X12)
              c_X21(c_indx_x_min, j) = rbufer(j+shift_X21)
              c_X22(c_indx_x_min, j) = rbufer(j+shift_X22)
           END DO
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 12 ## receive from right values for the right edge
           CALL MPI_RECV(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho(c_indx_x_max, j) = rbufer(j+shift)
              c_X11(c_indx_x_max, j) = rbufer(j+shift_X11)
              c_X12(c_indx_x_max, j) = rbufer(j+shift_X12)
              c_X21(c_indx_x_max, j) = rbufer(j+shift_X21)
              c_X22(c_indx_x_max, j) = rbufer(j+shift_X22)
           END DO
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:5*n3), STAT=ALLOC_ERR)

        shift = -c_indx_x_min+1
        shift_X11 = shift + n3
        shift_X12 = shift_X11 + n3
        shift_X21 = shift_X12 + n3
        shift_X22 = shift_X21 + n3

        IF (Rank_horizontal_above.GE.0) THEN
! ## 13 ## send up values below the top edge
           rbufer(     1:  n3) = c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           rbufer(  n3+1:2*n3) = c_X11(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           rbufer(2*n3+1:3*n3) = c_X12(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           rbufer(3*n3+1:4*n3) = c_X21(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           rbufer(4*n3+1:5*n3) = c_X22(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           CALL MPI_SEND(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 14 ## send down values above the bottom edge
           rbufer(     1:  n3) = c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           rbufer(  n3+1:2*n3) = c_X11(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           rbufer(2*n3+1:3*n3) = c_X12(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           rbufer(3*n3+1:4*n3) = c_X21(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           rbufer(4*n3+1:5*n3) = c_X22(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           CALL MPI_SEND(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 15 ## receive from below values for the bottom edge
           CALL MPI_RECV(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_min) = rbufer(i+shift)
              c_X11(i, c_indx_y_min) = rbufer(i+shift_X11)
              c_X12(i, c_indx_y_min) = rbufer(i+shift_X12)
              c_X21(i, c_indx_y_min) = rbufer(i+shift_X21)
              c_X22(i, c_indx_y_min) = rbufer(i+shift_X22)
           END DO
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 16 ## receive from above values for the top edge
           CALL MPI_RECV(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_max) = rbufer(i+shift)
              c_X11(i, c_indx_y_max) = rbufer(i+shift_X11)
              c_X12(i, c_indx_y_max) = rbufer(i+shift_X12)
              c_X21(i, c_indx_y_max) = rbufer(i+shift_X21)
              c_X22(i, c_indx_y_max) = rbufer(i+shift_X22)
           END DO
        END IF   

     ELSE ! "black" cluster

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:5*n1), STAT=ALLOC_ERR)

        shift = -c_indx_y_min+1
        shift_X11 = shift + n1
        shift_X12 = shift_X11 + n1
        shift_X21 = shift_X12 + n1
        shift_X22 = shift_X21 + n1

        IF (Rank_horizontal_left.GE.0) THEN
! ## 09 ## receive from left values for the left edge
           CALL MPI_RECV(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal_left, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho(c_indx_x_min, j) = rbufer(j+shift)
              c_X11(c_indx_x_min, j) = rbufer(j+shift_X11)
              c_X12(c_indx_x_min, j) = rbufer(j+shift_X12)
              c_X21(c_indx_x_min, j) = rbufer(j+shift_X21)
              c_X22(c_indx_x_min, j) = rbufer(j+shift_X22)
           END DO
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 10 ## receive from right values for the right edge
           CALL MPI_RECV(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal_right, COMM_HORIZONTAL, stattus, ierr)
           DO j = c_indx_y_min, c_indx_y_max
              c_rho(c_indx_x_max, j) = rbufer(j+shift)
              c_X11(c_indx_x_max, j) = rbufer(j+shift_X11)
              c_X12(c_indx_x_max, j) = rbufer(j+shift_X12)
              c_X21(c_indx_x_max, j) = rbufer(j+shift_X21)
              c_X22(c_indx_x_max, j) = rbufer(j+shift_X22)
           END DO
        END IF

        IF (Rank_horizontal_right.GE.0) THEN
! ## 11 ## send right values to the left of the right edge
           rbufer(     1:  n1) = c_rho(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           rbufer(  n1+1:2*n1) = c_X11(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           rbufer(2*n1+1:3*n1) = c_X12(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           rbufer(3*n1+1:4*n1) = c_X21(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           rbufer(4*n1+1:5*n1) = c_X22(c_indx_x_max-1, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_right, Rank_horizontal, COMM_HORIZONTAL, request, ierr)
        END IF

        IF (Rank_horizontal_left.GE.0) THEN
! ## 12 ## send left values to the right of the left edge
           rbufer(     1:  n1) = c_rho(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           rbufer(  n1+1:2*n1) = c_X11(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           rbufer(2*n1+1:3*n1) = c_X12(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           rbufer(3*n1+1:4*n1) = c_X21(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           rbufer(4*n1+1:5*n1) = c_X22(c_indx_x_min+1, c_indx_y_min:c_indx_y_max)
           CALL MPI_SEND(rbufer, 5*n1, MPI_DOUBLE_PRECISION, Rank_horizontal_left, Rank_horizontal, COMM_HORIZONTAL, request, ierr)
        END IF

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:5*n3), STAT=ALLOC_ERR)

        shift = -c_indx_x_min+1
        shift_X11 = shift + n3
        shift_X12 = shift_X11 + n3
        shift_X21 = shift_X12 + n3
        shift_X22 = shift_X21 + n3

        IF (Rank_horizontal_below.GE.0) THEN
! ## 13 ## receive from below values for the bottom edge
           CALL MPI_RECV(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal_below, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_min) = rbufer(i+shift)
              c_X11(i, c_indx_y_min) = rbufer(i+shift_X11)
              c_X12(i, c_indx_y_min) = rbufer(i+shift_X12)
              c_X21(i, c_indx_y_min) = rbufer(i+shift_X21)
              c_X22(i, c_indx_y_min) = rbufer(i+shift_X22)
           END DO
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 14 ## receive from above values for the top edge
           CALL MPI_RECV(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal_above, COMM_HORIZONTAL, stattus, ierr)
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_max) = rbufer(i+shift)
              c_X11(i, c_indx_y_max) = rbufer(i+shift_X11)
              c_X12(i, c_indx_y_max) = rbufer(i+shift_X12)
              c_X21(i, c_indx_y_max) = rbufer(i+shift_X21)
              c_X22(i, c_indx_y_max) = rbufer(i+shift_X22)
           END DO
        END IF

        IF (Rank_horizontal_above.GE.0) THEN
! ## 15 ## send up values below the top edge
           rbufer(     1:  n3) = c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           rbufer(  n3+1:2*n3) = c_X11(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           rbufer(2*n3+1:3*n3) = c_X12(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           rbufer(3*n3+1:4*n3) = c_X21(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           rbufer(4*n3+1:5*n3) = c_X22(c_indx_x_min:c_indx_x_max, c_indx_y_max-1)
           CALL MPI_SEND(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_above, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

        IF (Rank_horizontal_below.GE.0) THEN
! ## 16 ## send down values above the bottom edge
           rbufer(     1:  n3) = c_rho(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           rbufer(  n3+1:2*n3) = c_X11(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           rbufer(2*n3+1:3*n3) = c_X12(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           rbufer(3*n3+1:4*n3) = c_X21(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           rbufer(4*n3+1:5*n3) = c_X22(c_indx_x_min:c_indx_x_max, c_indx_y_min+1)
           CALL MPI_SEND(rbufer, 5*n3, MPI_DOUBLE_PRECISION, Rank_horizontal_below, Rank_horizontal, COMM_HORIZONTAL, request, ierr) 
        END IF

     END IF ! for black/white selection for setting charge density at the edge 

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)

     CALL MPI_BARRIER(COMM_HORIZONTAL, ierr)

! adjust values c_rho at the boundaries with material walls

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! special case of self-connected X-periodic cluster
        IF (Rank_of_master_above.LT.0) THEN
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_max) = 2.0_8 * c_rho(i, c_indx_y_max)
           END DO
        END IF
  
        IF (Rank_of_master_below.LT.0) THEN
           DO i = c_indx_x_min, c_indx_x_max
              c_rho(i, c_indx_y_min) = 2.0_8 * c_rho(i, c_indx_y_min)
           END DO
        END IF

     ELSE !self-periodicity check
 
        IF (Rank_of_master_left.LT.0) THEN !includes symmetry line, if present
           DO j = c_indx_y_min+1, c_indx_y_max-1
              c_rho(c_indx_x_min, j) = 2.0_8 * c_rho(c_indx_x_min, j)
           END DO
        END IF

        IF (Rank_of_master_right.LT.0) THEN
           DO j = c_indx_y_min+1, c_indx_y_max-1
              c_rho(c_indx_x_max, j) = 2.0_8 * c_rho(c_indx_x_max, j)
           END DO
        END IF
  
        IF (Rank_of_master_above.LT.0) THEN
           DO i = c_indx_x_min+1, c_indx_x_max-1
              c_rho(i, c_indx_y_max) = 2.0_8 * c_rho(i, c_indx_y_max)
           END DO
        END IF
  
        IF (Rank_of_master_below.LT.0) THEN
           DO i = c_indx_x_min+1, c_indx_x_max-1
              c_rho(i, c_indx_y_min) = 2.0_8 * c_rho(i, c_indx_y_min)
           END DO
        END IF

        SELECT CASE (c_left_bottom_corner_type)
           CASE (SURROUNDED_BY_WALL)
              c_rho(c_indx_x_min, c_indx_y_min) = 4.0_8 * c_rho(c_indx_x_min, c_indx_y_min) 
           CASE (EMPTY_CORNER_WALL_LEFT)
              c_rho(c_indx_x_min, c_indx_y_min+1) = 0.66666666666666_8 * c_rho(c_indx_x_min, c_indx_y_min+1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_BELOW)
              c_rho(c_indx_x_min+1, c_indx_y_min) = 0.66666666666666_8 * c_rho(c_indx_x_min+1, c_indx_y_min)    ! 2*2/3=4/3=1/(3/4)
        END SELECT
     
        SELECT CASE (c_left_top_corner_type)
           CASE (SURROUNDED_BY_WALL)
              c_rho(c_indx_x_min, c_indx_y_max) = 4.0_8 * c_rho(c_indx_x_min, c_indx_y_max) 
           CASE (EMPTY_CORNER_WALL_LEFT)
              c_rho(c_indx_x_min, c_indx_y_max-1) = 0.66666666666666_8 * c_rho(c_indx_x_min, c_indx_y_max-1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_ABOVE)
              c_rho(c_indx_x_min+1, c_indx_y_max) = 0.66666666666666_8 * c_rho(c_indx_x_min+1, c_indx_y_max)    ! 2*2/3=4/3=1/(3/4)
        END SELECT

        SELECT CASE (c_right_bottom_corner_type)
           CASE (SURROUNDED_BY_WALL)
              c_rho(c_indx_x_max, c_indx_y_min) = 4.0_8 * c_rho(c_indx_x_max, c_indx_y_min) 
           CASE (EMPTY_CORNER_WALL_RIGHT)
              c_rho(c_indx_x_max, c_indx_y_min+1) = 0.66666666666666_8 * c_rho(c_indx_x_max, c_indx_y_min+1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_BELOW)
              c_rho(c_indx_x_max-1, c_indx_y_min) = 0.66666666666666_8 * c_rho(c_indx_x_max-1, c_indx_y_min)    ! 2*2/3=4/3=1/(3/4)
        END SELECT
     
        SELECT CASE (c_right_top_corner_type)
           CASE (SURROUNDED_BY_WALL)
              c_rho(c_indx_x_max, c_indx_y_max) = 4.0_8 * c_rho(c_indx_x_max, c_indx_y_max) 
           CASE (EMPTY_CORNER_WALL_RIGHT)
              c_rho(c_indx_x_max, c_indx_y_max-1) = 0.66666666666666_8 * c_rho(c_indx_x_max, c_indx_y_max-1)    ! 2*2/3=4/3=1/(3/4)
           CASE (EMPTY_CORNER_WALL_ABOVE)
              c_rho(c_indx_x_max-1, c_indx_y_max) = 0.66666666666666_8 * c_rho(c_indx_x_max-1, c_indx_y_max)    ! 2*2/3=4/3=1/(3/4)
        END SELECT

     END IF  !###   IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN

! prepare and send accumulated charge density to field calculators
     DO k = 2, cluster_N_blocks

        shift = (field_calculator(k)%indx_x_max - field_calculator(k)%indx_x_min + 1) * &
              & (field_calculator(k)%indx_y_max - field_calculator(k)%indx_y_min + 1)
        bufsize = 5*shift

        shift_X11 = shift
        shift_X12 = shift_X11 + shift
        shift_X21 = shift_X12 + shift
        shift_X22 = shift_X21 + shift

        IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
        ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)
        pos=0
        DO j = field_calculator(k)%indx_y_min, field_calculator(k)%indx_y_max
           DO i = field_calculator(k)%indx_x_min, field_calculator(k)%indx_x_max
              pos = pos + 1
              rbufer(pos) = c_rho(i,j)
              DO nio = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
                 CALL CHECK_IF_INNER_OBJECT_CONTAINS_POINT(whole_object(nio), i, j, position_flag)
                 rbufer(pos) = rbufer(pos) + Get_Surface_Charge_Inner_Object(i, j, position_flag, whole_object(nio))   !### was -
              END DO
              rbufer(pos+shift_X11) = c_X11(i,j)
              rbufer(pos+shift_X12) = c_X12(i,j)
              rbufer(pos+shift_X21) = c_X21(i,j)
              rbufer(pos+shift_X22) = c_X22(i,j)
           END DO
        END DO

        CALL MPI_SEND(rbufer, bufsize, MPI_DOUBLE_PRECISION, field_calculator(k)%rank, Rank_of_process, MPI_COMM_WORLD, request, ierr) 
     END DO

! cluster master is a field calculator too, prepare its own charge density
     DO j = indx_y_min, indx_y_max
        DO i = indx_x_min, indx_x_max
           rho(i, j) = c_rho(i, j)
           DO nio = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
              CALL CHECK_IF_INNER_OBJECT_CONTAINS_POINT(whole_object(nio), i, j, position_flag)
              rho(i, j) = rho(i, j) + Get_Surface_Charge_Inner_Object(i, j, position_flag, whole_object(nio))   !### was -
           END DO
           X11(i,j) = c_X11(i,j)
           X12(i,j) = c_X12(i,j)
           X21(i,j) = c_X21(i,j)
           X22(i,j) = c_X22(i,j)
        END DO
     END DO

     IF (ALLOCATED(c_rho)) DEALLOCATE(c_rho, STAT=ALLOC_ERR)
     IF (ALLOCATED(c_X11)) DEALLOCATE(c_X11, STAT=ALLOC_ERR)
     IF (ALLOCATED(c_X12)) DEALLOCATE(c_X12, STAT=ALLOC_ERR)
     IF (ALLOCATED(c_X21)) DEALLOCATE(c_X21, STAT=ALLOC_ERR)
     IF (ALLOCATED(c_X22)) DEALLOCATE(c_X22, STAT=ALLOC_ERR)

  ELSE !field master selection    IF (cluster_rank_key.EQ.0) THEN ! nonzero rank within cluster:
       
     shift = (indx_x_max - indx_x_min + 1) * (indx_y_max - indx_y_min + 1)  
     bufsize = 5 * shift

     shift_X11 = shift
     shift_X12 = shift_X11 + shift
     shift_X21 = shift_X12 + shift
     shift_X22 = shift_X21 + shift

     IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
     ALLOCATE(rbufer(1:bufsize), STAT=ALLOC_ERR)
        
     CALL MPI_RECV(rbufer, bufsize, MPI_DOUBLE_PRECISION, field_master, field_master, MPI_COMM_WORLD, stattus, ierr)
        
     pos = 0
     DO j = indx_y_min, indx_y_max
        DO i = indx_x_min, indx_x_max
           pos = pos + 1
           rho(i,j) = rbufer(pos)
           X11(i,j) = rbufer(pos+shift_X11)
           X12(i,j) = rbufer(pos+shift_X12)
           X21(i,j) = rbufer(pos+shift_X21)
           X22(i,j) = rbufer(pos+shift_X22)
        END DO
     END DO
     
  END IF
  
! cleanup
  IF (ALLOCATED(rbufer))  DEALLOCATE(rbufer,  STAT=ALLOC_ERR)
  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT=ALLOC_ERR)

END SUBROUTINE GATHER_STREAMING_CHARGE_DENSITY
