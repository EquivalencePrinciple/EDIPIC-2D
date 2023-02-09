!-------------------------------------------
!
SUBROUTINE ADVANCE_ELECTRONS_FINAL_PLUS

  USE ParallelOperationValues, ONLY : cluster_rank_key
  USE CurrentProblemValues, ONLY : T_cntr
  USE Snapshots
  USE ClusterAndItsBoundaries, ONLY : c_indx_x_min, c_indx_x_max, c_indx_y_min, c_indx_y_max, &
                                    & c_indx_x_min_ext, c_indx_x_max_ext, c_indx_y_min_ext, c_indx_y_max_ext
  USE Diagnostics, ONLY : Save_probes_data_T_cntr

  IMPLICIT NONE

  LOGICAL collect_electron_moments_2d_now
  INTEGER ALLOC_ERR

  collect_electron_moments_2d_now = .FALSE.

  IF ((current_snap.GE.1).AND.(current_snap.LE.N_of_all_snaps)) THEN
     IF (T_cntr.EQ.Tcntr_snapshot(current_snap)) THEN

        IF ( save_data(4).OR. &
           & save_data(5).OR. &
           & save_data(6).OR. &
           & save_data(7).OR. &
           & save_data(8).OR. &
           & save_data(9).OR. &
           & save_data(10).OR. &
           & save_data(11).OR. &
           & save_data(12).OR. &
           & save_data(13).OR. &
           & save_data(14).OR. &
           & save_data(15).OR. &
           & save_data(16).OR. &
           & save_data(17).OR. &
           & save_data(18).OR. &
           & save_data(19).OR. &
           & save_data(20).OR. &
           & save_data(21).OR. &
           & save_data(22) ) collect_electron_moments_2d_now = .TRUE.
     END IF
  END IF

!print *, "ADVANCE_ELECTRONS_FINAL_PLUS ", collect_electron_moments_2d_now

  IF (collect_electron_moments_2d_now) THEN

     IF (cluster_rank_key.EQ.0) THEN

! arrays required by procedures collecting moments of electron and ion vdfs

        ALLOCATE(cs_N(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min_ext:c_indx_y_max_ext), STAT=ALLOC_ERR)

        ALLOCATE(cs_VX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min_ext:c_indx_y_max_ext), STAT=ALLOC_ERR)
        ALLOCATE(cs_VY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min_ext:c_indx_y_max_ext), STAT=ALLOC_ERR)
        ALLOCATE(cs_VZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min_ext:c_indx_y_max_ext), STAT=ALLOC_ERR)

        ALLOCATE(cs_WX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min_ext:c_indx_y_max_ext), STAT=ALLOC_ERR)
        ALLOCATE(cs_WY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min_ext:c_indx_y_max_ext), STAT=ALLOC_ERR)
        ALLOCATE(cs_WZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min_ext:c_indx_y_max_ext), STAT=ALLOC_ERR)

        ALLOCATE(cs_VXVY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min_ext:c_indx_y_max_ext), STAT=ALLOC_ERR)
        ALLOCATE(cs_VXVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min_ext:c_indx_y_max_ext), STAT=ALLOC_ERR)
        ALLOCATE(cs_VYVZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min_ext:c_indx_y_max_ext), STAT=ALLOC_ERR)

        ALLOCATE(cs_QX(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min_ext:c_indx_y_max_ext), STAT=ALLOC_ERR)
        ALLOCATE(cs_QY(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min_ext:c_indx_y_max_ext), STAT=ALLOC_ERR)
        ALLOCATE(cs_QZ(c_indx_x_min_ext:c_indx_x_max_ext, c_indx_y_min_ext:c_indx_y_max_ext), STAT=ALLOC_ERR)

! output arrays
        ALLOCATE(cs_Ne(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

        ALLOCATE(cs_VXe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_VYe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_VZe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

        ALLOCATE(cs_WXe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_WYe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_WZe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

        ALLOCATE(cs_QXe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_QYe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)
        ALLOCATE(cs_QZe(c_indx_x_min:c_indx_x_max, c_indx_y_min:c_indx_y_max), STAT=ALLOC_ERR)

     ELSE   !###   IF (cluster_rank_key.EQ.0) THEN

! the moments of the distribution functions are calculated using MPI_REDUCE which takes sum of values from all particle calculators
! the results are stored in master processes in arrays cs_N, cs_VX, etc
! the particle calculators, in general, don't need these arrays at all, 
! so allocating these arrays in non-master processes is just a waste of memory
! however, the compiler reports an error if the code is compiled with -C flag (check everything)
! this can be avoided if at least some minimal size arrays are allocated in the non-master processes
        ALLOCATE(cs_N(1,1), STAT=ALLOC_ERR)

        ALLOCATE(cs_VX(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_VY(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_VZ(1,1), STAT=ALLOC_ERR)

        ALLOCATE(cs_WX(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_WY(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_WZ(1,1), STAT=ALLOC_ERR)

        ALLOCATE(cs_VXVY(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_VXVZ(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_VYVZ(1,1), STAT=ALLOC_ERR)

        ALLOCATE(cs_QX(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_QY(1,1), STAT=ALLOC_ERR)
        ALLOCATE(cs_QZ(1,1), STAT=ALLOC_ERR)

     END IF   !###   IF (cluster_rank_key.EQ.0) THEN

! clear arrays
     cs_N = 0.0

     cs_VX = 0.0
     cs_VY = 0.0
     cs_VZ = 0.0

     cs_WX = 0.0
     cs_WY = 0.0
     cs_WZ = 0.0

     cs_VXVY = 0.0
     cs_VXVZ = 0.0
     cs_VYVZ = 0.0

     cs_QX = 0.0
     cs_QY = 0.0
     cs_QZ = 0.0

     CALL ADVANCE_ELECTRONS_FINAL_WITH_MOMENTS_2D

! cleanup

     DEALLOCATE(cs_N, STAT = ALLOC_ERR)

     DEALLOCATE(cs_VX, STAT = ALLOC_ERR)
     DEALLOCATE(cs_VY, STAT = ALLOC_ERR)
     DEALLOCATE(cs_VZ, STAT = ALLOC_ERR)

     DEALLOCATE(cs_WX, STAT = ALLOC_ERR)
     DEALLOCATE(cs_WY, STAT = ALLOC_ERR)
     DEALLOCATE(cs_WZ, STAT = ALLOC_ERR)

     DEALLOCATE(cs_VXVY, STAT = ALLOC_ERR)
     DEALLOCATE(cs_VXVZ, STAT = ALLOC_ERR)
     DEALLOCATE(cs_VYVZ, STAT = ALLOC_ERR)

     DEALLOCATE(cs_QX, STAT = ALLOC_ERR)
     DEALLOCATE(cs_QY, STAT = ALLOC_ERR)
     DEALLOCATE(cs_QZ, STAT = ALLOC_ERR)

  ELSE   !###   IF (collect_electron_moments_2d_now) THEN

     IF (T_cntr.EQ.Save_probes_data_T_cntr) THEN
        CALL ADVANCE_ELECTRONS_FINAL_WITH_MOMENTS_PROBES
     ELSE
        CALL ADVANCE_ELECTRONS_FINAL
     END IF

  END IF   !###   IF (collect_electron_moments_2d_now) THEN

END SUBROUTINE ADVANCE_ELECTRONS_FINAL_PLUS

!---------------------------------------
!streaming density not adjusted at the bounding surfaces
!adjusted at the symmetry line only
SUBROUTINE ADVANCE_ELECTRONS_FINAL

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE ClusterAndItsBoundaries
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER k
  INTEGER i, j, ix, jx, iy, jy
  REAL(8) X_Ex, Y_Ey
  REAL(8) ax_ip1, ax_i, ax_jp1, ax_j !meaning of these coefficients has changed
  REAL(8) ay_ip1, ay_i, ay_jp1, ay_j
  REAL(8) E_X, E_Y
  REAL(8) xarg, yarg
  REAL(8) alfa_x, alfa_y, alfa_z
  REAL(8) alfa_x2, alfa_y2, alfa_z2
  REAL(8) theta2, invtheta
  REAL(8) K11, K12, K21, K22, K31, K32
  REAL(8) A11, A12, A21, A22, A31, A32
  REAL(8) dVx, dVy, dVz

  INTEGER n
  LOGICAL collision_with_inner_object_occurred

! functions
  REAL(8) Bx, By, Bz, Ez

!print *, "ADVANCE_ELECTRONS_FINAL ", Rank_of_process

! clear counters of particles to be sent to neighbor processes
  N_e_to_send_left  = 0
  N_e_to_send_right = 0
  N_e_to_send_above = 0
  N_e_to_send_below = 0

!###??? do not do this ???
! clear counters of particles that hit boundary objects
!###  whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count = 0
!###  e_colls_with_bo(1:N_of_boundary_and_inner_objects)%N_of_saved_parts = 0

! clear counters of particles emitted by boundary objects in order to account for the secondary electron emission
!  whole_object(1:N_of_boundary_and_inner_objects)%electron_emit_count = 0  !### ?????

!print *, "enter", Rank_of_process, N_electrons, max_N_electrons

  k=0
  DO WHILE (k.LT.N_electrons)

     k = k + 1

     if ( (electron(k)%X.lt.c_X_area_min).or. &
        & (electron(k)%X.gt.c_X_area_max).or. &
        & (electron(k)%Y.lt.c_Y_area_min).or. &
        & (electron(k)%Y.gt.c_Y_area_max) ) then
        print '("Process ",i4," : Error-1 in ADVANCE_ELECTRONS : particle out of bounds xmin/xmax/ymin/ymax : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
        print '("Process ",i4," : k/N_electrons : ",i8,2x,i8)', Rank_of_process, k, N_electrons
        print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     end if

! interpolate electric field

     i = FLOOR(electron(k)%X)
     j = FLOOR(electron(k)%Y)

     IF (electron(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
     IF (electron(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

     ix = FLOOR(electron(k)%X - 0.5_8)   ! Ex x-indices limits: c_indx_x_min - 1 : c_indx_x_max
     jx = j

     iy = i
     jy = FLOOR(electron(k)%Y - 0.5_8)   ! Ey y-indices limits: c_indx_y_min - 1 : c_indx_y_max

     X_Ex = DBLE(ix) + 0.5_8             ! x-coordinate of Ex node ix
     Y_Ey = DBLE(jy) + 0.5_8             ! y-coordinate of Ey node jy

     ax_ip1 = electron(k)%X - X_Ex
     ax_i   = 1.0_8 - ax_ip1

     ax_jp1 = electron(k)%Y - DBLE(jx)
     ax_j   = 1.0_8 - ax_jp1

     ay_ip1 = electron(k)%X - DBLE(iy)
     ay_i   = 1.0_8 - ay_ip1

     ay_jp1 = electron(k)%Y - Y_Ey
     ay_j   = 1.0_8 - ay_jp1

! interpolate Ex and Ey:
     E_X = EX(ix, jx) * ax_i * ax_j + EX(ix+1, jx) * ax_ip1 * ax_j + EX(ix, jx+1) * ax_i * ax_jp1 + EX(ix+1, jx+1) * ax_ip1 * ax_jp1
     E_Y = EY(iy, jy) * ay_i * ay_j + EY(iy+1, jy) * ay_ip1 * ay_j + EY(iy, jy+1) * ay_i * ay_jp1 + EY(iy+1, jy+1) * ay_ip1 * ay_jp1

! calculate magnetic field factors

! shift the electrons back to the previous step to get the mag. field
     Xarg = electron(k)%X - electron(k)%VX !! coordinate at t^n, restored, only used to evaluate mag. field fron analyt. expressions
     Yarg = electron(k)%Y - electron(k)%VY !

     alfa_x = -0.5_8 * Bx(Xarg, Yarg)
     alfa_y = -0.5_8 * By(Xarg, Yarg)
     alfa_z = -0.5_8 * Bz(Xarg, Yarg)

     alfa_x2 = alfa_x**2
     alfa_y2 = alfa_y**2
     alfa_z2 = alfa_z**2

     theta2 = alfa_x2 + alfa_y2 + alfa_z2
     invtheta = 1.0_8 / (1.0_8 + theta2)
    
!    matrix K, same as R in Gibbons & Hewett; A_inverse = (K + I)/2 = (R + I)/2
     K11 =  (1.0_8 - theta2 + 2.0_8 * alfa_x2) * invtheta
     K12 =  2.0_8 * (alfa_x * alfa_y + alfa_z) * invtheta
!     K13 =  2.0_8 * (alfa_x * alfa_z - alfa_y) * invtheta

     K21 =  2.0_8 * (alfa_x * alfa_y - alfa_z) * invtheta
     K22 =  (1.0_8 - theta2 + 2.0_8 * alfa_y2) * invtheta
!     K23 =  2.0_8 * (alfa_y * alfa_z + alfa_x) * invtheta

     K31 =  2.0_8 * (alfa_x * alfa_z + alfa_y) * invtheta
     K32 =  2.0_8 * (alfa_y * alfa_z - alfa_x) * invtheta
!     K33 =  (1.0_8 - theta2 + 2.0_8 * alfa_z2) * invtheta

     A11 = 0.5_8 * (K11 + 1.0_8)
     A12 = 0.5_8 * K12
!     A13 = 0.5_8 * K13

     A21 = 0.5_8 * K21
     A22 = 0.5_8 * (K22 + 1.0_8)
!     A23 = 0.5_8 * K23

     A31 = 0.5_8 * K31
     A32 = 0.5_8 * K32
!     A33 = 0.5_8 * (K33 + 1.0_8)

! Gibbons & Hewett Eq.(2.7)             
     dVx = -0.5_8 * (A11 * E_X + A12 * E_Y)
     dVy = -0.5_8 * (A21 * E_X + A22 * E_Y)
     dVz = -0.5_8 * (A31 * E_X + A32 * E_Y)
        
     electron(k)%VX = electron(k)%VX + dVx
     electron(k)%VY = electron(k)%VY + dVy
     electron(k)%VZ = electron(k)%VZ + dVz

     electron(k)%X = electron(k)%X + dVx
     electron(k)%Y = electron(k)%Y + dVy

     ! a particle crossed symmetry plane, reflect it (both predictive and final push)
     IF (symmetry_plane_X_left) THEN
        IF (electron(k)%X.LT.c_X_area_min) THEN
           electron(k)%X = MAX(c_X_area_min, c_X_area_min + c_X_area_min - electron(k)%X)
           electron(k)%VX = - electron(k)%VX
!!***06/12/22           electron(k)%VZ = -electron(k)%VZ ! enforces symmetry wih mag. field?
!           electron(k)%AX = -electron(k)%AX !not sure, but may not matter as AX = 0 on the symmetry line
!###          electron(k)%VZ = -electron(k)%VZ   !###??? do we not have to do this when BY is on ???
        END IF
     END IF

! check whether a collision with an inner object occurred
     collision_with_inner_object_occurred = .FALSE.
     DO n = N_of_boundary_objects + 1, N_of_boundary_and_inner_objects
        IF (electron(k)%X.LE.whole_object(n)%Xmin) CYCLE
        IF (electron(k)%X.GE.whole_object(n)%Xmax) CYCLE
        IF (electron(k)%Y.LE.whole_object(n)%Ymin) CYCLE
        IF (electron(k)%Y.GE.whole_object(n)%Ymax) CYCLE
! collision detected
!        write (*, *) "TRY 2"
        CALL TRY_ELECTRON_COLL_WITH_INNER_OBJECT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
        CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
        collision_with_inner_object_occurred = .TRUE.
        EXIT
     END DO

     IF (collision_with_inner_object_occurred) CYCLE

! most probable situation when the particle remains inside the area
     IF ( (electron(k)%X.GE.c_X_area_min) .AND. &
        & (electron(k)%X.LE.c_X_area_max) .AND. &
        & (electron(k)%Y.GE.c_Y_area_min) .AND. &
        & (electron(k)%Y.LE.c_Y_area_max) ) CYCLE

! since we are here, a particle did not collide with an inner object but crossed an area boundary
! note, in a periodic system the particle still may collide with an inner object after transfer to the other domain
! therefore it is still necessary to check add list for collisions after process receives particles from neighbors

     CALL PROCESS_ELECTRON_LEAVING_CLUSTER(k)

  END DO

END SUBROUTINE ADVANCE_ELECTRONS_FINAL

!---------------------------------------
!streaming density not adjusted at the bounding surfaces
!adjusted at the symmetry line only
SUBROUTINE ADVANCE_ELECTRONS_FINAL_WITH_MOMENTS_2D

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE ClusterAndItsBoundaries
  
!------------------------------------------>>>
  USE Snapshots
  USE Diagnostics
!------------------------------------------<<<

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER k
  INTEGER i, j, ix, jx, iy, jy
  REAL(8) X_Ex, Y_Ey
  REAL(8) ax_ip1, ax_i, ax_jp1, ax_j !meaning of these coefficients has changed
  REAL(8) ay_ip1, ay_i, ay_jp1, ay_j
  REAL(8) E_X, E_Y
  REAL(8) xarg, yarg
  REAL(8) alfa_x, alfa_y, alfa_z
  REAL(8) alfa_x2, alfa_y2, alfa_z2
  REAL(8) theta2, invtheta
  REAL(8) K11, K12, K21, K22, K31, K32
  REAL(8) A11, A12, A21, A22, A31, A32
  REAL(8) dVx, dVy, dVz

  INTEGER n
  LOGICAL collision_with_inner_object_occurred

!------------------------------------------>>>
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

  REAL, ALLOCATABLE :: rbufer_vxvy(:)
  REAL, ALLOCATABLE :: rbufer_vxvz(:)
  REAL, ALLOCATABLE :: rbufer_vyvz(:)

  REAL, ALLOCATABLE :: rbufer_vx3(:)
  REAL, ALLOCATABLE :: rbufer_vy3(:)
  REAL, ALLOCATABLE :: rbufer_vz3(:)

  INTEGER ALLOC_ERR
  INTEGER bufsize

  REAL(8) Xhalf, Yhalf

  INTEGER pos_i_j, pos_ip1_j, pos_i_jp1, pos_ip1_jp1
  REAL vij, vip1j, vijp1, vip1jp1

  REAL rvx, rvy, rvz
  REAL rvx2, rvy2, rvz2
  REAL rvxvy, rvxvz, rvyvz
  REAL rvx3, rvy3, rvz3

  INTEGER npc, npa

  REAL inv_N, rtemp
!------------------------------------------<<<

! functions
  REAL(8) Bx, By, Bz, Ez

!print *, "ADVANCE_ELECTRONS_FINAL_WITH_MOMENTS_2D ", Rank_of_process

!------------------------------------------->>>
  n1 = c_indx_y_max_ext - c_indx_y_min_ext + 1
  n3 = c_indx_x_max_ext - c_indx_x_min_ext + 1
  n2 = -c_indx_x_min_ext + 1 - c_indx_y_min_ext * n3

  bufsize = n1 * n3
  ALLOCATE(rbufer_n(1:bufsize), STAT=ALLOC_ERR)

  ALLOCATE(rbufer_vx(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vy(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vz(1:bufsize), STAT=ALLOC_ERR)

  ALLOCATE(rbufer_vx2(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vy2(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vz2(1:bufsize), STAT=ALLOC_ERR)

  ALLOCATE(rbufer_vxvy(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vxvz(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vyvz(1:bufsize), STAT=ALLOC_ERR)

  ALLOCATE(rbufer_vx3(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vy3(1:bufsize), STAT=ALLOC_ERR)
  ALLOCATE(rbufer_vz3(1:bufsize), STAT=ALLOC_ERR)

  rbufer_n = 0.0

  rbufer_vx = 0.0
  rbufer_vy = 0.0
  rbufer_vz = 0.0

  rbufer_vx2 = 0.0
  rbufer_vy2 = 0.0
  rbufer_vz2 = 0.0

  rbufer_vxvy = 0.0
  rbufer_vxvz = 0.0
  rbufer_vyvz = 0.0

  rbufer_vx3 = 0.0
  rbufer_vy3 = 0.0
  rbufer_vz3 = 0.0
!-------------------------------------------<<<

! clear counters of particles to be sent to neighbor processes
  N_e_to_send_left  = 0
  N_e_to_send_right = 0
  N_e_to_send_above = 0
  N_e_to_send_below = 0

!###??? do not do this ???
! clear counters of particles that hit boundary objects
!###  whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count = 0
!###  e_colls_with_bo(1:N_of_boundary_and_inner_objects)%N_of_saved_parts = 0

! clear counters of particles emitted by boundary objects in order to account for the secondary electron emission
!  whole_object(1:N_of_boundary_and_inner_objects)%electron_emit_count = 0  !### ?????

!print *, "enter", Rank_of_process, N_electrons, max_N_electrons

  k=0
  DO WHILE (k.LT.N_electrons)

     k = k + 1

     if ( (electron(k)%X.lt.c_X_area_min).or. &
        & (electron(k)%X.gt.c_X_area_max).or. &
        & (electron(k)%Y.lt.c_Y_area_min).or. &
        & (electron(k)%Y.gt.c_Y_area_max) ) then
        print '("Process ",i4," : Error-1 in ADVANCE_ELECTRONS : particle out of bounds xmin/xmax/ymin/ymax : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
        print '("Process ",i4," : k/N_electrons : ",i8,2x,i8)', Rank_of_process, k, N_electrons
        print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     end if

! interpolate electric field

     i = FLOOR(electron(k)%X)
     j = FLOOR(electron(k)%Y)

     IF (electron(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
     IF (electron(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

     ix = FLOOR(electron(k)%X - 0.5_8)   ! Ex x-indices limits: c_indx_x_min - 1 : c_indx_x_max
     jx = j

     iy = i
     jy = FLOOR(electron(k)%Y - 0.5_8)   ! Ey y-indices limits: c_indx_y_min - 1 : c_indx_y_max

     X_Ex = DBLE(ix) + 0.5_8             ! x-coordinate of Ex node ix
     Y_Ey = DBLE(jy) + 0.5_8             ! y-coordinate of Ey node jy

     ax_ip1 = electron(k)%X - X_Ex
     ax_i   = 1.0_8 - ax_ip1

     ax_jp1 = electron(k)%Y - DBLE(jx)
     ax_j   = 1.0_8 - ax_jp1

     ay_ip1 = electron(k)%X - DBLE(iy)
     ay_i   = 1.0_8 - ay_ip1

     ay_jp1 = electron(k)%Y - Y_Ey
     ay_j   = 1.0_8 - ay_jp1

! interpolate Ex and Ey:
     E_X = EX(ix, jx) * ax_i * ax_j + EX(ix+1, jx) * ax_ip1 * ax_j + EX(ix, jx+1) * ax_i * ax_jp1 + EX(ix+1, jx+1) * ax_ip1 * ax_jp1
     E_Y = EY(iy, jy) * ay_i * ay_j + EY(iy+1, jy) * ay_ip1 * ay_j + EY(iy, jy+1) * ay_i * ay_jp1 + EY(iy+1, jy+1) * ay_ip1 * ay_jp1

! calculate magnetic field factors

! shift the electrons back to the previous step to get the mag. field
     Xarg = electron(k)%X - electron(k)%VX !! coordinate at t^n, restored, only used to evaluate mag. field fron analyt. expressions
     Yarg = electron(k)%Y - electron(k)%VY !

     alfa_x = -0.5_8 * Bx(Xarg, Yarg)
     alfa_y = -0.5_8 * By(Xarg, Yarg)
     alfa_z = -0.5_8 * Bz(Xarg, Yarg)

     alfa_x2 = alfa_x**2
     alfa_y2 = alfa_y**2
     alfa_z2 = alfa_z**2

     theta2 = alfa_x2 + alfa_y2 + alfa_z2
     invtheta = 1.0_8 / (1.0_8 + theta2)
    
!    matrix K, same as R in Gibbons & Hewett; A_inverse = (K + I)/2 = (R + I)/2
     K11 =  (1.0_8 - theta2 + 2.0_8 * alfa_x2) * invtheta
     K12 =  2.0_8 * (alfa_x * alfa_y + alfa_z) * invtheta
!     K13 =  2.0_8 * (alfa_x * alfa_z - alfa_y) * invtheta

     K21 =  2.0_8 * (alfa_x * alfa_y - alfa_z) * invtheta
     K22 =  (1.0_8 - theta2 + 2.0_8 * alfa_y2) * invtheta
!     K23 =  2.0_8 * (alfa_y * alfa_z + alfa_x) * invtheta

     K31 =  2.0_8 * (alfa_x * alfa_z + alfa_y) * invtheta
     K32 =  2.0_8 * (alfa_y * alfa_z - alfa_x) * invtheta
!     K33 =  (1.0_8 - theta2 + 2.0_8 * alfa_z2) * invtheta

     A11 = 0.5_8 * (K11 + 1.0_8)
     A12 = 0.5_8 * K12
!     A13 = 0.5_8 * K13

     A21 = 0.5_8 * K21
     A22 = 0.5_8 * (K22 + 1.0_8)
!     A23 = 0.5_8 * K23

     A31 = 0.5_8 * K31
     A32 = 0.5_8 * K32
!     A33 = 0.5_8 * (K33 + 1.0_8)

! Gibbons & Hewett Eq.(2.7)             
     dVx = -0.5_8 * (A11 * E_X + A12 * E_Y)
     dVy = -0.5_8 * (A21 * E_X + A22 * E_Y)
     dVz = -0.5_8 * (A31 * E_X + A32 * E_Y)
        
     electron(k)%VX = electron(k)%VX + dVx   ! velocities at t^{n+1/2}
     electron(k)%VY = electron(k)%VY + dVy
     electron(k)%VZ = electron(k)%VZ + dVz

     electron(k)%X = electron(k)%X + dVx     ! coordinates at t^{n+1}
     electron(k)%Y = electron(k)%Y + dVy

! coordinates at t^n+1/2
     Xhalf = 0.5_8 * (Xarg + electron(k)%X)
     Yhalf = 0.5_8 * (Yarg + electron(k)%Y)

     ! a particle crossed symmetry plane, reflect it (both predictive and final push)
     IF (symmetry_plane_X_left) THEN
        IF (electron(k)%X.LT.c_X_area_min) THEN
           electron(k)%X = MAX(c_X_area_min, c_X_area_min + c_X_area_min - electron(k)%X)
           electron(k)%VX = - electron(k)%VX
!!***06/12/22           electron(k)%VZ = -electron(k)%VZ ! enforces symmetry wih mag. field?
!           electron(k)%AX = -electron(k)%AX !not sure, but may not matter as AX = 0 on the symmetry line
!###          electron(k)%VZ = -electron(k)%VZ   !###??? do we not have to do this when BY is on ???
        END IF
        IF (Xhalf.LT.c_X_area_min) THEN
           Xhalf = MAX(c_X_area_min, c_X_area_min + c_X_area_min - Xhalf)
        END IF
     END IF

! check whether a collision with an inner object occurred
     collision_with_inner_object_occurred = .FALSE.
     DO n = N_of_boundary_objects + 1, N_of_boundary_and_inner_objects
        IF (electron(k)%X.LE.whole_object(n)%Xmin) CYCLE
        IF (electron(k)%X.GE.whole_object(n)%Xmax) CYCLE
        IF (electron(k)%Y.LE.whole_object(n)%Ymin) CYCLE
        IF (electron(k)%Y.GE.whole_object(n)%Ymax) CYCLE
! collision detected
!        write (*, *) "TRY 2"
        CALL TRY_ELECTRON_COLL_WITH_INNER_OBJECT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
        CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
        collision_with_inner_object_occurred = .TRUE.
        EXIT
     END DO

     IF (collision_with_inner_object_occurred) CYCLE

     IF ( (Xhalf.GE.c_X_area_min_ext) .AND. &
        & (Xhalf.LE.c_X_area_max_ext) .AND. &
        & (Yhalf.GE.c_Y_area_min_ext) .AND. &
        & (Yhalf.LE.c_Y_area_max_ext) ) THEN
! account for the contribution of the particle to the EVDF moments using velocities and coordinates at t^{n+1/2}

        i = FLOOR(Xhalf)
        j = FLOOR(Yhalf)
        
        pos_i_j     = i + j * n3 + n2
        pos_ip1_j   = pos_i_j + 1
        pos_i_jp1   = pos_i_j + n3
        pos_ip1_jp1 = pos_i_jp1 + 1

        ax_ip1 = Xhalf - DBLE(i)
        ax_i   = 1.0_8 - ax_ip1

        ay_jp1 = Yhalf - DBLE(j)
        ay_j = 1.0_8 - ay_jp1

        vij   = REAL(ax_i   * ay_j)
        vip1j = REAL(ax_ip1 * ay_j)
        vijp1 = REAL(ax_i   * ay_jp1)
        vip1jp1 = 1.0 - vij - vip1j - vijp1

        rbufer_n(pos_i_j)     = rbufer_n(pos_i_j)     + vij     !ax_i   * ay_j
        rbufer_n(pos_ip1_j)   = rbufer_n(pos_ip1_j)   + vip1j   !ax_ip1 * ay_j
        rbufer_n(pos_i_jp1)   = rbufer_n(pos_i_jp1)   + vijp1   !ax_i   * ay_jp1
        rbufer_n(pos_ip1_jp1) = rbufer_n(pos_ip1_jp1) + vip1jp1 !ax_ip1 * ay_jp1

        rvx = REAL(electron(k)%VX)
        rvy = REAL(electron(k)%VY)
        rvz = REAL(electron(k)%VZ)

        rvx2 = rvx * rvx
        rvy2 = rvy * rvy
        rvz2 = rvz * rvz

        rvxvy = rvx * rvy
        rvxvz = rvx * rvz
        rvyvz = rvy * rvz

        rvx3 = (rvx2 + rvy2 + rvz2) * rvx
        rvy3 = (rvx2 + rvy2 + rvz2) * rvy
        rvz3 = (rvx2 + rvy2 + rvz2) * rvz

        rbufer_vx(pos_i_j)     = rbufer_vx(pos_i_j)     + rvx * vij
        rbufer_vx(pos_ip1_j)   = rbufer_vx(pos_ip1_j)   + rvx * vip1j
        rbufer_vx(pos_i_jp1)   = rbufer_vx(pos_i_jp1)   + rvx * vijp1
        rbufer_vx(pos_ip1_jp1) = rbufer_vx(pos_ip1_jp1) + rvx * vip1jp1

        rbufer_vy(pos_i_j)     = rbufer_vy(pos_i_j)     + rvy * vij
        rbufer_vy(pos_ip1_j)   = rbufer_vy(pos_ip1_j)   + rvy * vip1j
        rbufer_vy(pos_i_jp1)   = rbufer_vy(pos_i_jp1)   + rvy * vijp1
        rbufer_vy(pos_ip1_jp1) = rbufer_vy(pos_ip1_jp1) + rvy * vip1jp1

        rbufer_vz(pos_i_j)     = rbufer_vz(pos_i_j)     + rvz * vij
        rbufer_vz(pos_ip1_j)   = rbufer_vz(pos_ip1_j)   + rvz * vip1j
        rbufer_vz(pos_i_jp1)   = rbufer_vz(pos_i_jp1)   + rvz * vijp1
        rbufer_vz(pos_ip1_jp1) = rbufer_vz(pos_ip1_jp1) + rvz * vip1jp1

        rbufer_vx2(pos_i_j)     = rbufer_vx2(pos_i_j)     + rvx2 * vij
        rbufer_vx2(pos_ip1_j)   = rbufer_vx2(pos_ip1_j)   + rvx2 * vip1j
        rbufer_vx2(pos_i_jp1)   = rbufer_vx2(pos_i_jp1)   + rvx2 * vijp1
        rbufer_vx2(pos_ip1_jp1) = rbufer_vx2(pos_ip1_jp1) + rvx2 * vip1jp1

        rbufer_vy2(pos_i_j)     = rbufer_vy2(pos_i_j)     + rvy2 * vij
        rbufer_vy2(pos_ip1_j)   = rbufer_vy2(pos_ip1_j)   + rvy2 * vip1j
        rbufer_vy2(pos_i_jp1)   = rbufer_vy2(pos_i_jp1)   + rvy2 * vijp1
        rbufer_vy2(pos_ip1_jp1) = rbufer_vy2(pos_ip1_jp1) + rvy2 * vip1jp1

        rbufer_vz2(pos_i_j)     = rbufer_vz2(pos_i_j)     + rvz2 * vij
        rbufer_vz2(pos_ip1_j)   = rbufer_vz2(pos_ip1_j)   + rvz2 * vip1j
        rbufer_vz2(pos_i_jp1)   = rbufer_vz2(pos_i_jp1)   + rvz2 * vijp1
        rbufer_vz2(pos_ip1_jp1) = rbufer_vz2(pos_ip1_jp1) + rvz2 * vip1jp1

        rbufer_vxvy(pos_i_j)     = rbufer_vxvy(pos_i_j)     + rvxvy * vij
        rbufer_vxvy(pos_ip1_j)   = rbufer_vxvy(pos_ip1_j)   + rvxvy * vip1j
        rbufer_vxvy(pos_i_jp1)   = rbufer_vxvy(pos_i_jp1)   + rvxvy * vijp1
        rbufer_vxvy(pos_ip1_jp1) = rbufer_vxvy(pos_ip1_jp1) + rvxvy * vip1jp1

        rbufer_vxvz(pos_i_j)     = rbufer_vxvz(pos_i_j)     + rvxvz * vij
        rbufer_vxvz(pos_ip1_j)   = rbufer_vxvz(pos_ip1_j)   + rvxvz * vip1j
        rbufer_vxvz(pos_i_jp1)   = rbufer_vxvz(pos_i_jp1)   + rvxvz * vijp1
        rbufer_vxvz(pos_ip1_jp1) = rbufer_vxvz(pos_ip1_jp1) + rvxvz * vip1jp1

        rbufer_vyvz(pos_i_j)     = rbufer_vyvz(pos_i_j)     + rvyvz * vij
        rbufer_vyvz(pos_ip1_j)   = rbufer_vyvz(pos_ip1_j)   + rvyvz * vip1j
        rbufer_vyvz(pos_i_jp1)   = rbufer_vyvz(pos_i_jp1)   + rvyvz * vijp1
        rbufer_vyvz(pos_ip1_jp1) = rbufer_vyvz(pos_ip1_jp1) + rvyvz * vip1jp1

        rbufer_vx3(pos_i_j)     = rbufer_vx3(pos_i_j)     + rvx3 * vij
        rbufer_vx3(pos_ip1_j)   = rbufer_vx3(pos_ip1_j)   + rvx3 * vip1j
        rbufer_vx3(pos_i_jp1)   = rbufer_vx3(pos_i_jp1)   + rvx3 * vijp1
        rbufer_vx3(pos_ip1_jp1) = rbufer_vx3(pos_ip1_jp1) + rvx3 * vip1jp1

        rbufer_vy3(pos_i_j)     = rbufer_vy3(pos_i_j)     + rvy3 * vij
        rbufer_vy3(pos_ip1_j)   = rbufer_vy3(pos_ip1_j)   + rvy3 * vip1j
        rbufer_vy3(pos_i_jp1)   = rbufer_vy3(pos_i_jp1)   + rvy3 * vijp1
        rbufer_vy3(pos_ip1_jp1) = rbufer_vy3(pos_ip1_jp1) + rvy3 * vip1jp1

        rbufer_vz3(pos_i_j)     = rbufer_vz3(pos_i_j)     + rvz3 * vij
        rbufer_vz3(pos_ip1_j)   = rbufer_vz3(pos_ip1_j)   + rvz3 * vip1j
        rbufer_vz3(pos_i_jp1)   = rbufer_vz3(pos_i_jp1)   + rvz3 * vijp1
        rbufer_vz3(pos_ip1_jp1) = rbufer_vz3(pos_ip1_jp1) + rvz3 * vip1jp1

     END IF   !###   IF ( (Xhalf.GE.c_X_area_min_ext).AND.(Xhalf.LE.c_X_area_max_ext).AND.(Yhalf.GE.c_Y_area_min_ext).AND.(Yhalf.LE.c_Y_area_max_ext) ) THEN

! most probable situation when the particle remains inside the area
     IF ( (electron(k)%X.GE.c_X_area_min) .AND. &
        & (electron(k)%X.LE.c_X_area_max) .AND. &
        & (electron(k)%Y.GE.c_Y_area_min) .AND. &
        & (electron(k)%Y.LE.c_Y_area_max) ) CYCLE

! since we are here, a particle did not collide with an inner object but crossed an area boundary
! note, in a periodic system the particle still may collide with an inner object after transfer to the other domain
! therefore it is still necessary to check add list for collisions after process receives particles from neighbors

     CALL PROCESS_ELECTRON_LEAVING_CLUSTER(k)

  END DO

!------------------------------------------->>>
! collect moments from all processes in a cluster
!
! here size of arrays cs_* is (c_indx_x_max_ext:c_indx_x_min_ext, c_indx_y_max_ext:c_indx_y_min_ext)
!
  CALL MPI_REDUCE(rbufer_n, cs_N, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_REDUCE(rbufer_vx, cs_VX, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vy, cs_VY, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vz, cs_VZ, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  
  CALL MPI_REDUCE(rbufer_vx2, cs_WX, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vy2, cs_WY, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vz2, cs_WZ, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  
  CALL MPI_REDUCE(rbufer_vxvy, cs_VXVY, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vxvz, cs_VXVZ, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vyvz, cs_VYVZ, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
    
  CALL MPI_REDUCE(rbufer_vx3, cs_QX, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vy3, cs_QY, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)
  CALL MPI_REDUCE(rbufer_vz3, cs_QZ, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  DEALLOCATE(rbufer_n, STAT=ALLOC_ERR)

  DEALLOCATE(rbufer_vx, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vy, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vz, STAT=ALLOC_ERR)

  DEALLOCATE(rbufer_vx2, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vy2, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vz2, STAT=ALLOC_ERR)

  DEALLOCATE(rbufer_vxvy, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vxvz, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vyvz, STAT=ALLOC_ERR)

  DEALLOCATE(rbufer_vx3, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vy3, STAT=ALLOC_ERR)
  DEALLOCATE(rbufer_vz3, STAT=ALLOC_ERR)

! particle calculators which are not masters are done
  IF (cluster_rank_key.NE.0) RETURN

  IF (T_cntr.EQ.Save_probes_data_T_cntr) THEN
! save electron VDF moments in probes
!
!### probes relevant to cluster domain are within the extended area now ###
!
! do this before synchronization because otherwise nodes in overlapped areas will not be processed correctly
! contributions from multiple clusters will be assembled in DO_PROBE_DIAGNOSTICS* procedures
     DO npc = 1, N_of_probes_cluster_ext
        npa = List_of_probes_cluster_ext(npc)    ! npa = n-umber of p-robe, a-ll   npc = n-umber of p-robe in c-luster
        i = Probe_position(1, npa)
        j = Probe_position(2, npa)

        probe_Ne_cluster(npc) = cs_N(i,j)

        probe_VXe_cluster(npc) = cs_VX(i,j)
        probe_VYe_cluster(npc) = cs_VY(i,j)
        probe_VZe_cluster(npc) = cs_VZ(i,j)

        probe_WXe_cluster(npc) = cs_WX(i,j)
        probe_WYe_cluster(npc) = cs_WY(i,j)
        probe_WZe_cluster(npc) = cs_WZ(i,j)

        probe_VXVYe_cluster(npc) = cs_VXVY(i,j)
        probe_VXVZe_cluster(npc) = cs_VXVZ(i,j)
        probe_VYVZe_cluster(npc) = cs_VYVZ(i,j)

        probe_QXe_cluster(npc) = cs_QX(i,j)
        probe_QYe_cluster(npc) = cs_QY(i,j)
        probe_QZe_cluster(npc) = cs_QZ(i,j)
     END DO
  END IF

  CALL SYNCHRONIZE_MOMENTS_IN_OVERLAP_NODES   !### NEW :: increased overlap width

! calculate average dimensional velocities, energies, and heat flows
  DO j = c_indx_y_min, c_indx_y_max
     DO i = c_indx_x_min, c_indx_x_max
        IF (cs_N(i,j).GT.1.0e-9) THEN    ! note this is small but not zero

           inv_N = 1.0 / cs_N(i,j)

           cs_VX(i, j) = cs_VX(i, j) * inv_N
           cs_VY(i, j) = cs_VY(i, j) * inv_N
           cs_VZ(i, j) = cs_VZ(i, j) * inv_N

           cs_WX(i, j) = cs_WX(i, j) * inv_N
           cs_WY(i, j) = cs_WY(i, j) * inv_N
           cs_WZ(i, j) = cs_WZ(i, j) * inv_N

           cs_VXVY(i, j) = 2.0 * cs_VXVY(i, j) * inv_N
           cs_VXVZ(i, j) = 2.0 * cs_VXVZ(i, j) * inv_N
           cs_VYVZ(i, j) = 2.0 * cs_VYVZ(i, j) * inv_N

           rvx = cs_VX(i, j)
           rvy = cs_VY(i, j)
           rvz = cs_VZ(i, j)

           rvx2 = cs_WX(i, j)
           rvy2 = cs_WY(i, j)
           rvz2 = cs_WZ(i, j)

!           rtemp = 2.0 * (rvx * rvx + rvy * rvy + rvz * rvz) - rvx2 - rvy2 - rvz2
           rtemp = rvx * rvx + rvy * rvy + rvz * rvz
           rtemp = rtemp + rtemp - rvx2 - rvy2 - rvz2

           cs_QX(i,j) = cs_QX(i,j) * inv_N + &
                      & (rtemp - rvx2 - rvx2) * rvx - &
                      & cs_VXVY(i,j) * rvy - & 
                      & cs_VXVZ(i,j) * rvz

           cs_QY(i,j) = cs_QY(i,j) * inv_N - &
                      & cs_VXVY(i,j) * rvx + &
                      & (rtemp - rvy2 - rvy2) * rvy - &
                      & cs_VYVZ(i,j) * rvz

           cs_QZ(i,j) = cs_QZ(i,j) * inv_N - &
                      & cs_VXVZ(i,j) * rvx - &
                      & cs_VYVZ(i,j) * rvy + &
                      & (rtemp - rvz2 - rvz2) * rvz

        ELSE
           cs_N(i,j) = 0.0
           cs_VX(i, j) = 0.0
           cs_VY(i, j) = 0.0
           cs_VZ(i, j) = 0.0
           cs_WX(i, j) = 0.0
           cs_WY(i, j) = 0.0
           cs_WZ(i, j) = 0.0
           cs_VXVY(i, j) = 0.0
           cs_VXVZ(i, j) = 0.0
           cs_VYVZ(i, j) = 0.0
           cs_QX(i, j) = 0.0
           cs_QY(i, j) = 0.0
           cs_QZ(i, j) = 0.0
        END IF

     END DO
  END DO

! adjust moments at the boundaries with material walls

  CALL ADJUST_DENSITY_AT_WALL_BOUNDARIES

! transfer moments to output arrays
  DO j = c_indx_y_min, c_indx_y_max
     
     cs_Ne(c_indx_x_min:c_indx_x_max, j) = cs_N(c_indx_x_min:c_indx_x_max, j)

     cs_VXe(c_indx_x_min:c_indx_x_max, j) = cs_VX(c_indx_x_min:c_indx_x_max, j)
     cs_VYe(c_indx_x_min:c_indx_x_max, j) = cs_VY(c_indx_x_min:c_indx_x_max, j)
     cs_VZe(c_indx_x_min:c_indx_x_max, j) = cs_VZ(c_indx_x_min:c_indx_x_max, j)

     cs_WXe(c_indx_x_min:c_indx_x_max, j) = cs_WX(c_indx_x_min:c_indx_x_max, j)
     cs_WYe(c_indx_x_min:c_indx_x_max, j) = cs_WY(c_indx_x_min:c_indx_x_max, j)
     cs_WZe(c_indx_x_min:c_indx_x_max, j) = cs_WZ(c_indx_x_min:c_indx_x_max, j)

     cs_QXe(c_indx_x_min:c_indx_x_max, j) = cs_QX(c_indx_x_min:c_indx_x_max, j)
     cs_QYe(c_indx_x_min:c_indx_x_max, j) = cs_QY(c_indx_x_min:c_indx_x_max, j)
     cs_QZe(c_indx_x_min:c_indx_x_max, j) = cs_QZ(c_indx_x_min:c_indx_x_max, j)

  END DO
  
END SUBROUTINE ADVANCE_ELECTRONS_FINAL_WITH_MOMENTS_2D

!---------------------------------------
!streaming density not adjusted at the bounding surfaces
!adjusted at the symmetry line only
SUBROUTINE ADVANCE_ELECTRONS_FINAL_WITH_MOMENTS_PROBES

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE ClusterAndItsBoundaries
  
  USE Diagnostics

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

!------------------------------------------>>>
  INTEGER bufsize
  REAL, ALLOCATABLE :: rbufer_local(:)
  REAL, ALLOCATABLE :: rbufer_local_2(:)
  INTEGER ALLOC_ERR
!------------------------------------------<<<

  INTEGER k
  INTEGER i, j, ix, jx, iy, jy
  REAL(8) X_Ex, Y_Ey
  REAL(8) ax_ip1, ax_i, ax_jp1, ax_j !meaning of these coefficients has changed
  REAL(8) ay_ip1, ay_i, ay_jp1, ay_j
  REAL(8) E_X, E_Y
  REAL(8) xarg, yarg
  REAL(8) alfa_x, alfa_y, alfa_z
  REAL(8) alfa_x2, alfa_y2, alfa_z2
  REAL(8) theta2, invtheta
  REAL(8) K11, K12, K21, K22, K31, K32
  REAL(8) A11, A12, A21, A22, A31, A32
  REAL(8) dVx, dVy, dVz

!------------------------------------------>>>
  REAL(8) Xhalf, Yhalf
  LOGICAL no_probe_in_cell_corner
  INTEGER npc, npa, pos
  REAL weight
  REAL rvx, rvy, rvz
  REAL rvx2, rvy2, rvz2
  REAL rvxvy, rvxvz, rvyvz
  REAL rvx3, rvy3, rvz3
  REAL inv_N, rtemp
!------------------------------------------<<<

  INTEGER n
  LOGICAL collision_with_inner_object_occurred

! functions
  REAL(8) Bx, By, Bz, Ez

!print *, "ADVANCE_ELECTRONS_FINAL_WITH_MOMENTS_PROBES ", Rank_of_process

! clear counters of particles to be sent to neighbor processes
  N_e_to_send_left  = 0
  N_e_to_send_right = 0
  N_e_to_send_above = 0
  N_e_to_send_below = 0

!###??? do not do this ???
! clear counters of particles that hit boundary objects
!###  whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count = 0
!###  e_colls_with_bo(1:N_of_boundary_and_inner_objects)%N_of_saved_parts = 0

! clear counters of particles emitted by boundary objects in order to account for the secondary electron emission
!  whole_object(1:N_of_boundary_and_inner_objects)%electron_emit_count = 0  !### ?????

!------------------------------------------->>>
  IF (N_of_probes_cluster_ext.GT.0) THEN
     bufsize = 13 * N_of_probes_cluster_ext  ! {N,VX,VY,VZ,WX,WY,WZ,VXVY,VXVZ,VYVZ,QX,QY,QZ}
     ALLOCATE(rbufer_local(bufsize), STAT=ALLOC_ERR)
     ALLOCATE(rbufer_local_2(bufsize), STAT=ALLOC_ERR)
     rbufer_local   = 0.0
     rbufer_local_2 = 0.0
  END IF
!-------------------------------------------<<<

!print *, "enter", Rank_of_process, N_electrons, max_N_electrons

  k=0
  DO WHILE (k.LT.N_electrons)

     k = k + 1

     if ( (electron(k)%X.lt.c_X_area_min).or. &
        & (electron(k)%X.gt.c_X_area_max).or. &
        & (electron(k)%Y.lt.c_Y_area_min).or. &
        & (electron(k)%Y.gt.c_Y_area_max) ) then
        print '("Process ",i4," : Error-1 in ADVANCE_ELECTRONS : particle out of bounds xmin/xmax/ymin/ymax : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
        print '("Process ",i4," : k/N_electrons : ",i8,2x,i8)', Rank_of_process, k, N_electrons
        print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     end if

! interpolate electric field

     i = FLOOR(electron(k)%X)
     j = FLOOR(electron(k)%Y)

     IF (electron(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
     IF (electron(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

     ix = FLOOR(electron(k)%X - 0.5_8)   ! Ex x-indices limits: c_indx_x_min - 1 : c_indx_x_max
     jx = j

     iy = i
     jy = FLOOR(electron(k)%Y - 0.5_8)   ! Ey y-indices limits: c_indx_y_min - 1 : c_indx_y_max

     X_Ex = DBLE(ix) + 0.5_8             ! x-coordinate of Ex node ix
     Y_Ey = DBLE(jy) + 0.5_8             ! y-coordinate of Ey node jy

     ax_ip1 = electron(k)%X - X_Ex
     ax_i   = 1.0_8 - ax_ip1

     ax_jp1 = electron(k)%Y - DBLE(jx)
     ax_j   = 1.0_8 - ax_jp1

     ay_ip1 = electron(k)%X - DBLE(iy)
     ay_i   = 1.0_8 - ay_ip1

     ay_jp1 = electron(k)%Y - Y_Ey
     ay_j   = 1.0_8 - ay_jp1

! interpolate Ex and Ey:
     E_X = EX(ix, jx) * ax_i * ax_j + EX(ix+1, jx) * ax_ip1 * ax_j + EX(ix, jx+1) * ax_i * ax_jp1 + EX(ix+1, jx+1) * ax_ip1 * ax_jp1
     E_Y = EY(iy, jy) * ay_i * ay_j + EY(iy+1, jy) * ay_ip1 * ay_j + EY(iy, jy+1) * ay_i * ay_jp1 + EY(iy+1, jy+1) * ay_ip1 * ay_jp1

! calculate magnetic field factors

! shift the electrons back to the previous step to get the mag. field
     Xarg = electron(k)%X - electron(k)%VX !! coordinate at t^n, restored, only used to evaluate mag. field fron analyt. expressions
     Yarg = electron(k)%Y - electron(k)%VY !

     alfa_x = -0.5_8 * Bx(Xarg, Yarg)
     alfa_y = -0.5_8 * By(Xarg, Yarg)
     alfa_z = -0.5_8 * Bz(Xarg, Yarg)

     alfa_x2 = alfa_x**2
     alfa_y2 = alfa_y**2
     alfa_z2 = alfa_z**2

     theta2 = alfa_x2 + alfa_y2 + alfa_z2
     invtheta = 1.0_8 / (1.0_8 + theta2)
    
!    matrix K, same as R in Gibbons & Hewett; A_inverse = (K + I)/2 = (R + I)/2
     K11 =  (1.0_8 - theta2 + 2.0_8 * alfa_x2) * invtheta
     K12 =  2.0_8 * (alfa_x * alfa_y + alfa_z) * invtheta
!     K13 =  2.0_8 * (alfa_x * alfa_z - alfa_y) * invtheta

     K21 =  2.0_8 * (alfa_x * alfa_y - alfa_z) * invtheta
     K22 =  (1.0_8 - theta2 + 2.0_8 * alfa_y2) * invtheta
!     K23 =  2.0_8 * (alfa_y * alfa_z + alfa_x) * invtheta

     K31 =  2.0_8 * (alfa_x * alfa_z + alfa_y) * invtheta
     K32 =  2.0_8 * (alfa_y * alfa_z - alfa_x) * invtheta
!     K33 =  (1.0_8 - theta2 + 2.0_8 * alfa_z2) * invtheta

     A11 = 0.5_8 * (K11 + 1.0_8)
     A12 = 0.5_8 * K12
!     A13 = 0.5_8 * K13

     A21 = 0.5_8 * K21
     A22 = 0.5_8 * (K22 + 1.0_8)
!     A23 = 0.5_8 * K23

     A31 = 0.5_8 * K31
     A32 = 0.5_8 * K32
!     A33 = 0.5_8 * (K33 + 1.0_8)

! Gibbons & Hewett Eq.(2.7)             
     dVx = -0.5_8 * (A11 * E_X + A12 * E_Y)
     dVy = -0.5_8 * (A21 * E_X + A22 * E_Y)
     dVz = -0.5_8 * (A31 * E_X + A32 * E_Y)
        
     electron(k)%VX = electron(k)%VX + dVx   ! velocities at t^{n+1/2}
     electron(k)%VY = electron(k)%VY + dVy
     electron(k)%VZ = electron(k)%VZ + dVz

     electron(k)%X = electron(k)%X + dVx     ! coordinates at t^{n+1}
     electron(k)%Y = electron(k)%Y + dVy

! coordinates at t^n+1/2
     Xhalf = 0.5_8 * (Xarg + electron(k)%X)
     Yhalf = 0.5_8 * (Yarg + electron(k)%Y)

     ! a particle crossed symmetry plane, reflect it (both predictive and final push)
     IF (symmetry_plane_X_left) THEN
        IF (electron(k)%X.LT.c_X_area_min) THEN
           electron(k)%X = MAX(c_X_area_min, c_X_area_min + c_X_area_min - electron(k)%X)
           electron(k)%VX = - electron(k)%VX
!!***06/12/22           electron(k)%VZ = -electron(k)%VZ ! enforces symmetry wih mag. field?
!           electron(k)%AX = -electron(k)%AX !not sure, but may not matter as AX = 0 on the symmetry line
!###          electron(k)%VZ = -electron(k)%VZ   !###??? do we not have to do this when BY is on ???
        END IF
        IF (Xhalf.LT.c_X_area_min) THEN
           Xhalf = MAX(c_X_area_min, c_X_area_min + c_X_area_min - Xhalf)
        END IF
     END IF

! check whether a collision with an inner object occurred
     collision_with_inner_object_occurred = .FALSE.
     DO n = N_of_boundary_objects + 1, N_of_boundary_and_inner_objects
        IF (electron(k)%X.LE.whole_object(n)%Xmin) CYCLE
        IF (electron(k)%X.GE.whole_object(n)%Xmax) CYCLE
        IF (electron(k)%Y.LE.whole_object(n)%Ymin) CYCLE
        IF (electron(k)%Y.GE.whole_object(n)%Ymax) CYCLE
! collision detected
!        write (*, *) "TRY 2"
        CALL TRY_ELECTRON_COLL_WITH_INNER_OBJECT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
        CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
        collision_with_inner_object_occurred = .TRUE.
        EXIT
     END DO

     IF (collision_with_inner_object_occurred) CYCLE

! account for the contribution of the particle to the EVDF moments using velocities and coordinates at t^{n+1/2}

     i = FLOOR(Xhalf)
     j = FLOOR(Yhalf)
        
!--------------------------------->>>
! collect particle contribution to the electron moments in all close probes
     DO npc = 1, N_of_probes_cluster_ext

        npa = List_of_probes_cluster_ext(npc)
        no_probe_in_cell_corner = .TRUE.

        IF (i.EQ.Probe_position(1,npa)) THEN
           IF (j.EQ.Probe_position(2,npa)) THEN
! the probe is in the left bottom corner of the cell containing the particle
              no_probe_in_cell_corner = .FALSE.
              ax_ip1 = Xhalf - DBLE(i)
              ax_i   = 1.0_8 - ax_ip1
              ay_jp1 = Yhalf - DBLE(j)
              ay_j = 1.0_8 - ay_jp1
              weight = REAL(ax_i * ay_j)   ! wij
           ELSE IF ((j+1).EQ.Probe_position(2,npa)) THEN
! the probe is in the left top corner of the cell containing the particle
              no_probe_in_cell_corner = .FALSE.
              ax_ip1 = Xhalf - DBLE(i)
              ax_i   = 1.0_8 - ax_ip1
              ay_jp1 = Yhalf - DBLE(j)
              ay_j = 1.0_8 - ay_jp1
              weight = REAL(ax_i * ay_jp1)   ! vijp1
           END IF
        ELSE IF ((i+1).EQ.Probe_position(1,npa)) THEN
           IF (j.EQ.Probe_position(2,npa)) THEN
! the probe is in the right bottom corner of the cell containing the particle
              no_probe_in_cell_corner = .FALSE.
              ax_ip1 = Xhalf - DBLE(i)
              ax_i   = 1.0_8 - ax_ip1
              ay_jp1 = Yhalf - DBLE(j)
              ay_j = 1.0_8 - ay_jp1
              weight = REAL(ax_ip1 * ay_j)   ! vip1j
           ELSE IF ((j+1).EQ.Probe_position(2,npa)) THEN
! the probe is in the right top corner of the cell containing the particle
              no_probe_in_cell_corner = .FALSE.
              ax_ip1 = Xhalf - DBLE(i)
              ax_i   = 1.0_8 - ax_ip1
              ay_jp1 = Yhalf - DBLE(j)
              ay_j = 1.0_8 - ay_jp1
              weight = REAL(ax_ip1 * ay_jp1)   ! vip1jp1
           END IF
        END IF

        IF (no_probe_in_cell_corner) CYCLE

        pos = (npc-1) * 13

        rvx = REAL(electron(k)%VX)
        rvy = REAL(electron(k)%VY)
        rvz = REAL(electron(k)%VZ)

        rvx2 = rvx * rvx
        rvy2 = rvy * rvy
        rvz2 = rvz * rvz

        rvxvy = rvx * rvy
        rvxvz = rvx * rvz
        rvyvz = rvy * rvz

        rvx3 = (rvx2 + rvy2 + rvz2) * rvx
        rvy3 = (rvx2 + rvy2 + rvz2) * rvy
        rvz3 = (rvx2 + rvy2 + rvz2) * rvz

        rbufer_local(pos+1) = rbufer_local(pos+1) + weight

        rbufer_local(pos+2) = rbufer_local(pos+2) + weight * rvx
        rbufer_local(pos+3) = rbufer_local(pos+3) + weight * rvy
        rbufer_local(pos+4) = rbufer_local(pos+4) + weight * rvz

        rbufer_local(pos+5) = rbufer_local(pos+5) + weight * rvx2
        rbufer_local(pos+6) = rbufer_local(pos+6) + weight * rvy2
        rbufer_local(pos+7) = rbufer_local(pos+7) + weight * rvz2

        rbufer_local(pos+8)  = rbufer_local(pos+8)  + weight * rvxvy
        rbufer_local(pos+9)  = rbufer_local(pos+9)  + weight * rvxvz
        rbufer_local(pos+10) = rbufer_local(pos+10) + weight * rvyvz

        rbufer_local(pos+11) = rbufer_local(pos+11) + weight * rvx3
        rbufer_local(pos+12) = rbufer_local(pos+12) + weight * rvy3
        rbufer_local(pos+13) = rbufer_local(pos+13) + weight * rvz3

     END DO   !###   DO npc = 1, N_of_probes_cluster
!--------------------------------->>>

! most probable situation when the particle remains inside the area
     IF ( (electron(k)%X.GE.c_X_area_min) .AND. &
        & (electron(k)%X.LE.c_X_area_max) .AND. &
        & (electron(k)%Y.GE.c_Y_area_min) .AND. &
        & (electron(k)%Y.LE.c_Y_area_max) ) CYCLE

! since we are here, a particle did not collide with an inner object but crossed an area boundary
! note, in a periodic system the particle still may collide with an inner object after transfer to the other domain
! therefore it is still necessary to check add list for collisions after process receives particles from neighbors

     CALL PROCESS_ELECTRON_LEAVING_CLUSTER(k)

  END DO  !### DO WHILE (k.LT.N_electrons)
!------------------------------------------->>>

  IF (N_of_probes_cluster_ext.LE.0) RETURN

! collect moments from all processes in a cluster
  CALL MPI_REDUCE(rbufer_local, rbufer_local_2, bufsize, MPI_REAL, MPI_SUM, 0, COMM_CLUSTER, ierr)

  IF (cluster_rank_key.EQ.0) THEN
     pos=0
     DO npc = 1, N_of_probes_cluster_ext
        probe_Ne_cluster(npc) = rbufer_local_2(pos+1)

        probe_VXe_cluster(npc) = rbufer_local_2(pos+2)
        probe_VYe_cluster(npc) = rbufer_local_2(pos+3)
        probe_VZe_cluster(npc) = rbufer_local_2(pos+4)

        probe_WXe_cluster(npc) = rbufer_local_2(pos+5)
        probe_WYe_cluster(npc) = rbufer_local_2(pos+6)
        probe_WZe_cluster(npc) = rbufer_local_2(pos+7)

        probe_VXVYe_cluster(npc) = rbufer_local_2(pos+8)
        probe_VXVZe_cluster(npc) = rbufer_local_2(pos+9)
        probe_VYVZe_cluster(npc) = rbufer_local_2(pos+10)
        
        probe_QXe_cluster(npc) = rbufer_local_2(pos+11)
        probe_QYe_cluster(npc) = rbufer_local_2(pos+12)
        probe_QZe_cluster(npc) = rbufer_local_2(pos+13)
        
        pos = pos+13
     END DO
  END IF

! cleanup
  DEALLOCATE(rbufer_local, STAT = ALLOC_ERR)
  DEALLOCATE(rbufer_local_2, STAT = ALLOC_ERR)

END SUBROUTINE ADVANCE_ELECTRONS_FINAL_WITH_MOMENTS_PROBES








