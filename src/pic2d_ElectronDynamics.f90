!---------------------------------------
!streaming density not adjusted at the bounding surfaces
!adjusted at the symmetry line only
SUBROUTINE ADVANCE_ELECTRONS_PRELIMINARY

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE ClusterAndItsBoundaries
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER k
  REAL(8) E_Z 
  REAL(8) alfa_x, alfa_y, alfa_z
  REAL(8) alfa_x2, alfa_y2, alfa_z2
  REAL(8) theta2, invtheta
  REAL(8) K11, K12, K13, K21, K22, K23, K31, K32, K33
  REAL(8) A11, A12, A13, A21, A22, A23, A31, A32, A33
  REAL(8) VX_temp, VY_temp, VZ_temp
  REAL(8) vec1, vec2, vec3

  INTEGER n
  LOGICAL collision_with_inner_object_occurred

! functions
  REAL(8) Bx, By, Bz, Ez

! clear counters of particles to be sent to neighbor processes
  N_e_to_send_left  = 0
  N_e_to_send_right = 0
  N_e_to_send_above = 0
  N_e_to_send_below = 0

! clear counters of particles that hit boundary objects
  whole_object(1:N_of_boundary_and_inner_objects)%electron_hit_count = 0
  e_colls_with_bo(1:N_of_boundary_and_inner_objects)%N_of_saved_parts = 0

! clear counters of particles emitted by boundary objects in order to account for the secondary electron emission
  whole_object(1:N_of_boundary_and_inner_objects)%electron_emit_count = 0  !### ?????

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

     E_Z = Ez(electron(k)%X, electron(k)%Y)

! calculate magnetic field factors

     alfa_x = -0.5_8 * Bx(electron(k)%X, electron(k)%Y)
     alfa_y = -0.5_8 * By(electron(k)%X, electron(k)%Y)
     alfa_z = -0.5_8 * Bz(electron(k)%X, electron(k)%Y)

     alfa_x2 = alfa_x**2
     alfa_y2 = alfa_y**2
     alfa_z2 = alfa_z**2

     theta2 = alfa_x2 + alfa_y2 + alfa_z2
     invtheta = 1.0_8 / (1.0_8 + theta2)
    
! matrix K, same as R in Gibbons & Hewett; A_inverse = (I + K)/2 = (I + R)/2
     K11 =  (1.0_8 - theta2 + 2.0_8 * alfa_x2) * invtheta
     K12 =  2.0_8 * (alfa_x * alfa_y + alfa_z) * invtheta
     K13 =  2.0_8 * (alfa_x * alfa_z - alfa_y) * invtheta

     K21 =  2.0_8 * (alfa_x * alfa_y - alfa_z) * invtheta
     K22 =  (1.0_8 - theta2 + 2.0_8 * alfa_y2) * invtheta
     K23 =  2.0_8 * (alfa_y * alfa_z + alfa_x) * invtheta

     K31 =  2.0_8 * (alfa_x * alfa_z + alfa_y) * invtheta
     K32 =  2.0_8 * (alfa_y * alfa_z - alfa_x) * invtheta
     K33 =  (1.0_8 - theta2 + 2.0_8 * alfa_z2) * invtheta

     A11 = 0.5_8 * (K11 + 1.0_8)
     A12 = 0.5_8 * K12
     A13 = 0.5_8 * K13

     A21 = 0.5_8 * K21
     A22 = 0.5_8 * (K22 + 1.0_8)
     A23 = 0.5_8 * K23

     A31 = 0.5_8 * K31
     A32 = 0.5_8 * K32
     A33 = 0.5_8 * (K33 + 1.0_8)

! predictive move
! Gibbons & Hewett Eq.(2.6):             
     vec1 =  0.5_8 * electron(k)%AX
     vec2 =  0.5_8 * electron(k)%AY
     vec3 = -E_Z

     VX_temp = electron(k)%VX
     VY_temp = electron(k)%VY
     VZ_temp = electron(k)%VZ

     electron(k)%VX = K11 * VX_temp + K12 * VY_temp + K13 * VZ_temp + A11 * vec1 + A12 * vec2 + A13 * vec3
     electron(k)%VY = K21 * VX_temp + K22 * VY_temp + K23 * VZ_temp + A21 * vec1 + A22 * vec2 + A23 * vec3
     electron(k)%VZ = K31 * VX_temp + K32 * VY_temp + K33 * VZ_temp + A31 * vec1 + A32 * vec2 + A33 * vec3
        
     electron(k)%X = electron(k)%X + electron(k)%VX
     electron(k)%Y = electron(k)%Y + electron(k)%VY

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

END SUBROUTINE ADVANCE_ELECTRONS_PRELIMINARY

!---------------------------------------
!
SUBROUTINE UPDATE_ELECTRON_ACCELERATIONS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles
  USE ClusterAndItsBoundaries
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER k
  INTEGER i, j, ix, jx, iy, jy
  REAL(8) ax_ip1, ax_i, ax_jp1, ax_j !meaning of these coefficients has changed
  REAL(8) ay_ip1, ay_i, ay_jp1, ay_j
  REAL(8) X_Ex, Y_Ey
  REAL(8) E_X, E_Y

  DO k = 1, N_electrons

     if ( (electron(k)%X.lt.c_X_area_min).or. &
        & (electron(k)%X.gt.c_X_area_max).or. &
        & (electron(k)%Y.lt.c_Y_area_min).or. &
        & (electron(k)%Y.gt.c_Y_area_max) ) then
        print '("Process ",i4," : Error-1 in UPDATE_ELECTRON_ACCELERATIONS : particle out of bounds xmin/xmax/ymin/ymax : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
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

     electron(k)%AX = 0.5_8 * (electron(k)%AX - E_X)
     electron(k)%AY = 0.5_8 * (electron(k)%AY - E_Y)

  END DO

END SUBROUTINE UPDATE_ELECTRON_ACCELERATIONS

!-----------------------------------------------------------------------------------------------
!
SUBROUTINE PROCESS_ELECTRON_LEAVING_CLUSTER(k)

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE ElectronParticles

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER, INTENT(INOUT) :: k

  IF (electron(k)%X.LT.c_X_area_min) THEN

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the electron by the period length
        electron(k)%X = electron(k)%X + L_period_X
        IF (electron(k)%Y.LT.c_Y_area_min) THEN
! particle is somewhere near the left bottom cell of the area
           IF (Rank_of_master_below.LT.0) THEN
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           ELSE
              CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)                       
           END IF
           CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
        ELSE IF (electron(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
           IF (Rank_of_master_above.LT.0) THEN
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           ELSE
              CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
           END IF
           CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
        END IF
        RETURN
     END IF

     IF ( (electron(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
        & (electron(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

        IF (Rank_of_master_left.GE.0) THEN
! left neighbor cluster exists
           CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
        ELSE
! left neighbor cluster does not exist
           CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX,  electron(k)%VY, electron(k)%VZ, electron(k)%tag)   ! left
        END IF

     ELSE IF (electron(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the left bottom cell of the area

        SELECT CASE (c_left_bottom_corner_type)
           CASE (HAS_TWO_NEIGHBORS)
              IF ((c_X_area_min-electron(k)%X).LT.(c_Y_area_min-electron(k)%Y)) THEN
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
              END IF

           CASE (FLAT_WALL_BELOW)
              CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)

           CASE (FLAT_WALL_LEFT)
              IF (electron(k)%Y.GE.c_Y_area_min) THEN                 
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)                       
              END IF

           CASE (SURROUNDED_BY_WALL)
              IF ((c_X_area_min-electron(k)%X).LT.(c_Y_area_min-electron(k)%Y)) THEN
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              END IF

           CASE (EMPTY_CORNER_WALL_LEFT)
              CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
               
           CASE (EMPTY_CORNER_WALL_BELOW)
              CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)

        END SELECT

     ELSE IF (electron(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top left cell of the area

        SELECT CASE (c_left_top_corner_type)
           CASE (HAS_TWO_NEIGHBORS)
              IF ((c_X_area_min-electron(k)%X).LT.(electron(k)%Y-c_Y_area_max)) THEN
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
              END IF

           CASE (FLAT_WALL_ABOVE)
              CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)

           CASE (FLAT_WALL_LEFT)
              IF (electron(k)%Y.LE.c_Y_area_max) THEN                 
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
              END IF

           CASE (SURROUNDED_BY_WALL)
              IF ((c_X_area_min-electron(k)%X).LT.(electron(k)%Y-c_Y_area_max)) THEN
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              END IF

           CASE (EMPTY_CORNER_WALL_LEFT)
              CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)

           CASE (EMPTY_CORNER_WALL_ABOVE)
              CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)

        END SELECT

     ELSE
! ERROR, we shouldn't be here
        PRINT '("ERROR-1 in PROCESS_ELECTRON_LEAVING_CLUSTER: we should not be here")'
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     END IF
     CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
     RETURN
  END IF

  IF (electron(k)%X.GT.c_X_area_max) THEN

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the electron by the period length
        electron(k)%X = electron(k)%X - L_period_X
        IF (electron(k)%Y.LT.c_Y_area_min) THEN
! particle is somewhere near the left bottom cell of the area
           IF (Rank_of_master_below.LT.0) THEN
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           ELSE
              CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)                       
           END IF
           CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
        ELSE IF (electron(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
           IF (Rank_of_master_above.LT.0) THEN
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           ELSE
              CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
           END IF
           CALL REMOVE_ELECTRON(k)  ! this subroutine does  N_electrons = N_electrons - 1 and k = k-1
        END IF
        RETURN
     END IF

     IF ( (electron(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
        & (electron(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

        IF (Rank_of_master_right.GE.0) THEN
! right neighbor cluster exists
           CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
        ELSE
! right neighbor cluster does not exist
           CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
        END IF

     ELSE IF (electron(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the right bottom cell of the area

        SELECT CASE (c_right_bottom_corner_type)
           CASE (HAS_TWO_NEIGHBORS)
              IF ((electron(k)%X-c_X_area_max).LT.(c_Y_area_min-electron(k)%Y)) THEN
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
              END IF

           CASE (FLAT_WALL_BELOW)
              CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
                
           CASE (FLAT_WALL_RIGHT)
              IF (electron(k)%Y.GE.c_Y_area_min) THEN
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)                       
              END IF

           CASE (SURROUNDED_BY_WALL)
              IF ((electron(k)%X-c_X_area_max).LT.(c_Y_area_min-electron(k)%Y)) THEN
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              END IF

           CASE (EMPTY_CORNER_WALL_RIGHT)
              CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)

           CASE (EMPTY_CORNER_WALL_BELOW)
              CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)

        END SELECT

     ELSE IF (electron(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top right cell of the area

        SELECT CASE (c_right_top_corner_type)
           CASE (HAS_TWO_NEIGHBORS)
              IF ((electron(k)%X-c_X_area_max).LT.(electron(k)%Y-c_Y_area_max)) THEN
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
              END IF

           CASE (FLAT_WALL_ABOVE)
              CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)

           CASE (FLAT_WALL_RIGHT)
              IF (electron(k)%Y.LE.c_Y_area_max) THEN
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)                    
              END IF

           CASE (SURROUNDED_BY_WALL)
              IF ((electron(k)%X-c_X_area_max).LT.(electron(k)%Y-c_Y_area_max)) THEN
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              ELSE
                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
              END IF

           CASE (EMPTY_CORNER_WALL_RIGHT)
              CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)

           CASE (EMPTY_CORNER_WALL_ABOVE)
              CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)

        END SELECT

     ELSE
! ERROR, we shouldn't be here
        PRINT '("ERROR-2 in PROCESS_ELECTRON_LEAVING_CLUSTER: we should not be here")'
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     END IF
     CALL REMOVE_ELECTRON(k)  !       this subroutine does  N_electrons = N_electrons - 1
     RETURN
  END IF

! since we are here, c_X_area_min <= electron(k)%Y <= c_X_area_max

  IF (electron(k)%Y.GT.c_Y_area_max) THEN

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
        IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
           CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
        ELSE
           CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
        END IF
        CALL REMOVE_ELECTRON(k)  !       this subroutine does  N_electrons = N_electrons - 1
        RETURN
     END IF

     IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
        CALL ADD_ELECTRON_TO_SEND_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
     ELSE
! neighbor cluster above does not exist
        IF ((electron(k)%X.GE.(c_X_area_min+1.0_8)).AND.(electron(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
           CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
        ELSE IF (electron(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left top corner
           IF (c_left_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
              CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
           ELSE
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           END IF
        ELSE IF (electron(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right top corner
           IF (c_right_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
              CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX,  electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
           ELSE
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           END IF
        END IF
     END IF
     CALL REMOVE_ELECTRON(k)  !       this subroutine does  N_electrons = N_electrons - 1
     RETURN
  END IF

  IF (electron(k)%Y.LT.c_Y_area_min) THEN

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
        IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster above exists
           CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
        ELSE
           CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
        END IF
        CALL REMOVE_ELECTRON(k)  !       this subroutine does  N_electrons = N_electrons - 1
        RETURN
     END IF

     IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster below exists, remove particle and prepare to send it to the neighbor below
        CALL ADD_ELECTRON_TO_SEND_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
     ELSE
! neighbor cluster below does not exist
        IF ((electron(k)%X.GE.(c_X_area_min+1.0_8)).AND.(electron(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
           CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
        ELSE IF (electron(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left bottom corner
           IF (c_left_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
              CALL ADD_ELECTRON_TO_SEND_LEFT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
           ELSE
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           END IF
        ELSE IF (electron(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right bottom corner
           IF (c_right_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
              CALL ADD_ELECTRON_TO_SEND_RIGHT(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%AX, electron(k)%AY, electron(k)%tag)
           ELSE
              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron(k)%X, electron(k)%Y, electron(k)%VX, electron(k)%VY, electron(k)%VZ, electron(k)%tag)
           END IF
        END IF
     END IF
     CALL REMOVE_ELECTRON(k)  !       this subroutine does  N_electrons = N_electrons - 1
     RETURN
  END IF

END SUBROUTINE PROCESS_ELECTRON_LEAVING_CLUSTER
! 
!----------------------------------------
!
SUBROUTINE REMOVE_ELECTRON(k)

  USE ParallelOperationValues
  USE ElectronParticles

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER, INTENT(INOUT) :: k

  IF ((k.LT.1).OR.(k.GT.N_electrons)) THEN
     PRINT '("Process ",i6," : ERROR in REMOVE_ELECTRON : index k invalid")', Rank_of_process
     PRINT '("Process ",i6," : k= ", i7," N_electrons= ",i7)', Rank_of_process, k, N_electrons
     PRINT '("Process ",i6," : PROGRAM TERMINATED")', Rank_of_process
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

  IF (k.LT.N_electrons) THEN

     electron(k)%X   = electron(N_electrons)%X
     electron(k)%Y   = electron(N_electrons)%Y
     electron(k)%VX  = electron(N_electrons)%VX
     electron(k)%VY  = electron(N_electrons)%VY
     electron(k)%VZ  = electron(N_electrons)%VZ

     electron(k)%AX  = electron(N_electrons)%AX
     electron(k)%AY  = electron(N_electrons)%AY

     electron(k)%tag = electron(N_electrons)%tag

  END IF

  N_electrons = N_electrons - 1
  k = k-1                  ! to ensure that the new k-particle is processed

END SUBROUTINE REMOVE_ELECTRON

!------------------------------------------
!
SUBROUTINE ADD_ELECTRON_TO_SEND_LEFT(x, y, vx, vy, vz, ax, ay, tag)

  USE ElectronParticles, ONLY : N_e_to_send_left, max_N_e_to_send_left, electron_to_send_left
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_X_left, L_period_x

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

  REAL(8) x, y, vx, vy, vz, ax, ay
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     real(8) ax
     real(8) ay
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_e_to_send_left = N_e_to_send_left + 1

  IF (N_e_to_send_left.GT.max_N_e_to_send_left) THEN
! increase the size of the list array
     current_N = max_N_e_to_send_left
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = electron_to_send_left(k)%X
        bufer(k)%Y   = electron_to_send_left(k)%Y
        bufer(k)%VX  = electron_to_send_left(k)%VX
        bufer(k)%VY  = electron_to_send_left(k)%VY
        bufer(k)%VZ  = electron_to_send_left(k)%VZ

        bufer(k)%AX  = electron_to_send_left(k)%AX
        bufer(k)%AY  = electron_to_send_left(k)%AY

        bufer(k)%tag = electron_to_send_left(k)%tag
     END DO
     IF (ALLOCATED(electron_to_send_left)) DEALLOCATE(electron_to_send_left, STAT=DEALLOC_ERR)
     max_N_e_to_send_left = max_N_e_to_send_left + MAX(50, max_N_e_to_send_left/10)
     ALLOCATE(electron_to_send_left(1:max_N_e_to_send_left), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        electron_to_send_left(k)%X   = bufer(k)%X
        electron_to_send_left(k)%Y   = bufer(k)%Y
        electron_to_send_left(k)%VX  = bufer(k)%VX
        electron_to_send_left(k)%VY  = bufer(k)%VY
        electron_to_send_left(k)%VZ  = bufer(k)%VZ

        electron_to_send_left(k)%AX  = bufer(k)%AX
        electron_to_send_left(k)%AY  = bufer(k)%AY

        electron_to_send_left(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  IF (periodic_boundary_X_left) THEN
     electron_to_send_left(N_e_to_send_left)%X = x + L_period_X
  ELSE
     electron_to_send_left(N_e_to_send_left)%X = x
  END IF
  electron_to_send_left(N_e_to_send_left)%Y   = y
  electron_to_send_left(N_e_to_send_left)%VX  = vx
  electron_to_send_left(N_e_to_send_left)%VY  = vy
  electron_to_send_left(N_e_to_send_left)%VZ  = vz

  electron_to_send_left(N_e_to_send_left)%AX  = ax
  electron_to_send_left(N_e_to_send_left)%AY  = ay

  electron_to_send_left(N_e_to_send_left)%tag = tag

END SUBROUTINE ADD_ELECTRON_TO_SEND_LEFT

!------------------------------------------
!
SUBROUTINE ADD_ELECTRON_TO_SEND_RIGHT(x, y, vx, vy, vz, ax, ay, tag)

  USE ElectronParticles, ONLY : N_e_to_send_right, max_N_e_to_send_right, electron_to_send_right
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_X_right, L_period_x

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

  REAL(8) x, y, vx, vy, vz, ax, ay
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     real(8) ax
     real(8) ay
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_e_to_send_right = N_e_to_send_right + 1

  IF (N_e_to_send_right.GT.max_N_e_to_send_right) THEN
! increase the size of the list array
     current_N = max_N_e_to_send_right
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = electron_to_send_right(k)%X
        bufer(k)%Y   = electron_to_send_right(k)%Y
        bufer(k)%VX  = electron_to_send_right(k)%VX
        bufer(k)%VY  = electron_to_send_right(k)%VY
        bufer(k)%VZ  = electron_to_send_right(k)%VZ

        bufer(k)%AX  = electron_to_send_right(k)%AX
        bufer(k)%AY  = electron_to_send_right(k)%AY

        bufer(k)%tag = electron_to_send_right(k)%tag
     END DO
     IF (ALLOCATED(electron_to_send_right)) DEALLOCATE(electron_to_send_right, STAT=DEALLOC_ERR)
     max_N_e_to_send_right = max_N_e_to_send_right + MAX(50, max_N_e_to_send_right/10)
     ALLOCATE(electron_to_send_right(1:max_N_e_to_send_right), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        electron_to_send_right(k)%X   = bufer(k)%X
        electron_to_send_right(k)%Y   = bufer(k)%Y
        electron_to_send_right(k)%VX  = bufer(k)%VX
        electron_to_send_right(k)%VY  = bufer(k)%VY
        electron_to_send_right(k)%VZ  = bufer(k)%VZ

        electron_to_send_right(k)%AX  = bufer(k)%AX
        electron_to_send_right(k)%AY  = bufer(k)%AY

        electron_to_send_right(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  IF (periodic_boundary_X_right) THEN
     electron_to_send_right(N_e_to_send_right)%X = x - L_period_x
  ELSE
     electron_to_send_right(N_e_to_send_right)%X = x
  END IF
  electron_to_send_right(N_e_to_send_right)%Y   = y
  electron_to_send_right(N_e_to_send_right)%VX  = vx
  electron_to_send_right(N_e_to_send_right)%VY  = vy
  electron_to_send_right(N_e_to_send_right)%VZ  = vz

  electron_to_send_right(N_e_to_send_right)%AX  = ax
  electron_to_send_right(N_e_to_send_right)%AY  = ay

  electron_to_send_right(N_e_to_send_right)%tag = tag

END SUBROUTINE ADD_ELECTRON_TO_SEND_RIGHT

!------------------------------------------
!
SUBROUTINE ADD_ELECTRON_TO_SEND_ABOVE(x, y, vx, vy, vz, ax, ay, tag)

  USE ElectronParticles, ONLY : N_e_to_send_above, max_N_e_to_send_above, electron_to_send_above
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_Y_above, L_period_y

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

  REAL(8) x, y, vx, vy, vz, ax, ay
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     real(8) ax
     real(8) ay
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_e_to_send_above = N_e_to_send_above + 1

  IF (N_e_to_send_above.GT.max_N_e_to_send_above) THEN
! increase the size of the list array
     current_N = max_N_e_to_send_above
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = electron_to_send_above(k)%X
        bufer(k)%Y   = electron_to_send_above(k)%Y
        bufer(k)%VX  = electron_to_send_above(k)%VX
        bufer(k)%VY  = electron_to_send_above(k)%VY
        bufer(k)%VZ  = electron_to_send_above(k)%VZ

        bufer(k)%AX  = electron_to_send_above(k)%AX
        bufer(k)%AY  = electron_to_send_above(k)%AY

        bufer(k)%tag = electron_to_send_above(k)%tag
     END DO
     IF (ALLOCATED(electron_to_send_above)) DEALLOCATE(electron_to_send_above, STAT=DEALLOC_ERR)
     max_N_e_to_send_above = max_N_e_to_send_above + MAX(50, max_N_e_to_send_above/10)
     ALLOCATE(electron_to_send_above(1:max_N_e_to_send_above), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        electron_to_send_above(k)%X   = bufer(k)%X
        electron_to_send_above(k)%Y   = bufer(k)%Y
        electron_to_send_above(k)%VX  = bufer(k)%VX
        electron_to_send_above(k)%VY  = bufer(k)%VY
        electron_to_send_above(k)%VZ  = bufer(k)%VZ

        electron_to_send_above(k)%AX  = bufer(k)%AX
        electron_to_send_above(k)%AY  = bufer(k)%AY

        electron_to_send_above(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  electron_to_send_above(N_e_to_send_above)%X   = x
  IF (periodic_boundary_Y_above) THEN
     electron_to_send_above(N_e_to_send_above)%Y = y - L_period_y
  ELSE
     electron_to_send_above(N_e_to_send_above)%Y = y
  END IF
  electron_to_send_above(N_e_to_send_above)%VX  = vx
  electron_to_send_above(N_e_to_send_above)%VY  = vy
  electron_to_send_above(N_e_to_send_above)%VZ  = vz

  electron_to_send_above(N_e_to_send_above)%AX  = ax
  electron_to_send_above(N_e_to_send_above)%AY  = ay

  electron_to_send_above(N_e_to_send_above)%tag = tag

END SUBROUTINE ADD_ELECTRON_TO_SEND_ABOVE

!------------------------------------------
!
SUBROUTINE ADD_ELECTRON_TO_SEND_BELOW(x, y, vx, vy, vz, ax, ay, tag)

  USE ElectronParticles, ONLY : N_e_to_send_below, max_N_e_to_send_below, electron_to_send_below
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_Y_below, L_period_y

  IMPLICIT NONE

  REAL(8) x, y, vx, vy, vz, ax, ay
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     real(8) ax
     real(8) ay
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_e_to_send_below = N_e_to_send_below + 1

  IF (N_e_to_send_below.GT.max_N_e_to_send_below) THEN
! increase the size of the list array
     current_N = max_N_e_to_send_below
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = electron_to_send_below(k)%X
        bufer(k)%Y   = electron_to_send_below(k)%Y
        bufer(k)%VX  = electron_to_send_below(k)%VX
        bufer(k)%VY  = electron_to_send_below(k)%VY
        bufer(k)%VZ  = electron_to_send_below(k)%VZ

        bufer(k)%AX  = electron_to_send_below(k)%AX
        bufer(k)%AY  = electron_to_send_below(k)%AY

        bufer(k)%tag = electron_to_send_below(k)%tag
     END DO
     IF (ALLOCATED(electron_to_send_below)) DEALLOCATE(electron_to_send_below, STAT=DEALLOC_ERR)
     max_N_e_to_send_below = max_N_e_to_send_below + MAX(50, max_N_e_to_send_below/10)
     ALLOCATE(electron_to_send_below(1:max_N_e_to_send_below), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        electron_to_send_below(k)%X   = bufer(k)%X
        electron_to_send_below(k)%Y   = bufer(k)%Y
        electron_to_send_below(k)%VX  = bufer(k)%VX
        electron_to_send_below(k)%VY  = bufer(k)%VY
        electron_to_send_below(k)%VZ  = bufer(k)%VZ

        electron_to_send_below(k)%AX  = bufer(k)%AX
        electron_to_send_below(k)%AY  = bufer(k)%AY

        electron_to_send_below(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  electron_to_send_below(N_e_to_send_below)%X   = x
  IF (periodic_boundary_Y_below) THEN
     electron_to_send_below(N_e_to_send_below)%Y = y + L_period_y
  ELSE
     electron_to_send_below(N_e_to_send_below)%Y = y
  END IF
  electron_to_send_below(N_e_to_send_below)%VX  = vx
  electron_to_send_below(N_e_to_send_below)%VY  = vy
  electron_to_send_below(N_e_to_send_below)%VZ  = vz

  electron_to_send_below(N_e_to_send_below)%AX  = ax
  electron_to_send_below(N_e_to_send_below)%AY  = ay

  electron_to_send_below(N_e_to_send_below)%tag = tag

END SUBROUTINE ADD_ELECTRON_TO_SEND_BELOW

!----------------------------------------
!
SUBROUTINE ADD_ELECTRON_TO_ADD_LIST(x, y, vx, vy, vz, ax, ay, tag)

  USE ElectronParticles, ONLY : N_e_to_add, max_N_e_to_add, electron_to_add
!  USE ParallelOperationValues

  IMPLICIT NONE

  REAL(8) x, y, vx, vy, vz, ax, ay
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     real(8) ax
     real(8) ay
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_e_to_add = N_e_to_add + 1

  IF (N_e_to_add.GT.max_N_e_to_add) THEN
! increase the size of the list array
     current_N = max_N_e_to_add
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = electron_to_add(k)%X
        bufer(k)%Y   = electron_to_add(k)%Y
        bufer(k)%VX  = electron_to_add(k)%VX
        bufer(k)%VY  = electron_to_add(k)%VY
        bufer(k)%VZ  = electron_to_add(k)%VZ

        bufer(k)%AX  = electron_to_add(k)%AX
        bufer(k)%AY  = electron_to_add(k)%AY

        bufer(k)%tag = electron_to_add(k)%tag
     END DO
     IF (ALLOCATED(electron_to_add)) DEALLOCATE(electron_to_add, STAT=DEALLOC_ERR)
     max_N_e_to_add = max_N_e_to_add + MAX(50, max_N_e_to_add/10)
     ALLOCATE(electron_to_add(1:max_N_e_to_add), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        electron_to_add(k)%X   = bufer(k)%X
        electron_to_add(k)%Y   = bufer(k)%Y
        electron_to_add(k)%VX  = bufer(k)%VX
        electron_to_add(k)%VY  = bufer(k)%VY
        electron_to_add(k)%VZ  = bufer(k)%VZ

        electron_to_add(k)%AX  = bufer(k)%AX
        electron_to_add(k)%AY  = bufer(k)%AY

        electron_to_add(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  electron_to_add(N_e_to_add)%X   = x
  electron_to_add(N_e_to_add)%Y   = y
  electron_to_add(N_e_to_add)%VX  = vx
  electron_to_add(N_e_to_add)%VY  = vy
  electron_to_add(N_e_to_add)%VZ  = vz

  electron_to_add(N_e_to_add)%AX  = ax
  electron_to_add(N_e_to_add)%AY  = ay

  electron_to_add(N_e_to_add)%tag = tag

END SUBROUTINE ADD_ELECTRON_TO_ADD_LIST

!----------------------------------------
!
SUBROUTINE REMOVE_ELECTRON_FROM_ADD_LIST(k)

  USE ParallelOperationValues
  USE ElectronParticles

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER, INTENT(INOUT) :: k

  IF ((k.LT.1).OR.(k.GT.N_e_to_add)) THEN
     PRINT '("Process ",i6," : ERROR in REMOVE_ELECTRON_FROM_ADD_LIST : index k invalid")', Rank_of_process
     PRINT '("Process ",i6," : k= ", i7," N_e_to_add= ",i7)', Rank_of_process, k, N_e_to_add
     PRINT '("Process ",i6," : PROGRAM TERMINATED")', Rank_of_process
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

  IF (k.LT.N_e_to_add) THEN

     electron_to_add(k)%X   = electron_to_add(N_e_to_add)%X
     electron_to_add(k)%Y   = electron_to_add(N_e_to_add)%Y
     electron_to_add(k)%VX  = electron_to_add(N_e_to_add)%VX
     electron_to_add(k)%VY  = electron_to_add(N_e_to_add)%VY
     electron_to_add(k)%VZ  = electron_to_add(N_e_to_add)%VZ

     electron_to_add(k)%AX  = electron_to_add(N_e_to_add)%AX
     electron_to_add(k)%AY  = electron_to_add(N_e_to_add)%AY

     electron_to_add(k)%tag = electron_to_add(N_e_to_add)%tag

  END IF

  N_e_to_add = N_e_to_add - 1
  k = k-1                  ! to ensure that the new k-particle is processed

!print '("Process ",i4," called REMOVE_ELECTRON_FROM_ADD_LIST(k), T_cntr= ",i7," k= ",i4)', Rank_of_process, T_cntr, k

END SUBROUTINE REMOVE_ELECTRON_FROM_ADD_LIST

!----------------------------------------
! This subroutine is called after ADVANCE_ELECTRONS, before exchange of electrons takes place.
! It removes particles which do not belong to the cluster domain.
! Such particles may appear after emission from surface of inner objects.
! It is expected that the emission is never directed into boundary objects aligned along the main simulation domain's boundary.
! The algorithm below, which places alien particles into proper SEND* lists is the same as in ADVANCE_ELECTRONS.
! However, collision with a boundary object at this stage is not expected and is considered an error.
!
SUBROUTINE FIND_ALIENS_IN_ELECTRON_ADD_LIST

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE ElectronParticles, ONLY : N_e_to_add, electron_to_add

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER k, n 

  IF (N_of_inner_objects.EQ.0) RETURN

! find, process, and exclude electrons which collided with inner objects
  k=0
  DO WHILE (k.LT.N_e_to_add)
     k = k+1

     IF (symmetry_plane_X_left) THEN
        IF (electron_to_add(k)%X.LT.c_X_area_min) THEN
! since we are here, a particle to be added is beyond the symmetry plane at X=0
! this is a rare but possible situation, for example:
! a particle crosses the symmetry plane and collides with an inner object at the same time
! such a particle is reflected before trying for collision, but the crossing point for the reflected particle may still be at X<0
! [because for collisions with inner objects we do "ray tracing" to find the point of collision]
! if the collision is followed by emission of a secondary electron, then the emitted electron will have X<0
! so we simply move this particle rightward, symmetrically relative to plane X=0
           electron_to_add(k)%X = MAX(c_X_area_min, c_X_area_min + c_X_area_min - electron_to_add(k)%X)
! unlike in ADVANCE_ELECTRONS, here we do not change velocity because it is random
!           electron(k)%VX = -electron(k)%VX
!!###          electron(k)%VZ = -electron(k)%VZ   !###??? do we not have to do this when BY is on ??? 
        END IF
     END IF


! most probable situation when the particle is inside the area
     IF ( (electron_to_add(k)%X.GE.c_X_area_min) .AND. &
        & (electron_to_add(k)%X.LE.c_X_area_max) .AND. &
        & (electron_to_add(k)%Y.GE.c_Y_area_min) .AND. &
        & (electron_to_add(k)%Y.LE.c_Y_area_max) ) CYCLE

! since we are here, particle k is outside the domain of this cluster

     IF (electron_to_add(k)%X.LT.c_X_area_min) THEN

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the electron by the period length
           electron_to_add(k)%X = electron_to_add(k)%X + L_period_X
           IF (electron_to_add(k)%Y.LT.c_Y_area_min) THEN
! particle is below the bottom side of the area
              IF (Rank_of_master_below.LT.0) THEN
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-1 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)                       
              END IF
              CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
           ELSE IF (electron_to_add(k)%Y.GT.c_Y_area_max) THEN
! particle is above the top side of the area
              IF (Rank_of_master_above.LT.0) THEN
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-2 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
              END IF
              CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
           END IF
           CYCLE
        END IF

        IF ( (electron_to_add(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
           & (electron_to_add(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

           IF (Rank_of_master_left.GE.0) THEN
! left neighbor cluster exists
              CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
           ELSE
! left neighbor cluster does not exist
!#              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX,  electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)   ! left
! error
                 print '("error-3 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF

        ELSE IF (electron_to_add(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the left bottom cell of the area

           SELECT CASE (c_left_bottom_corner_type)
              CASE (HAS_TWO_NEIGHBORS)
                 IF ((c_X_area_min-electron_to_add(k)%X).LT.(c_Y_area_min-electron_to_add(k)%Y)) THEN
                    CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
                 END IF

              CASE (FLAT_WALL_BELOW)
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)

              CASE (FLAT_WALL_LEFT)
                 IF (electron_to_add(k)%Y.GE.c_Y_area_min) THEN                 
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                    print '("error-4 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)                       
                 END IF

              CASE (SURROUNDED_BY_WALL)
!                 IF ((c_X_area_min-electron_to_add(k)%X).LT.(c_Y_area_min-electron_to_add(k)%Y)) THEN
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 ELSE
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 END IF
! error
                 print '("error-5 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

              CASE (EMPTY_CORNER_WALL_LEFT)
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)

              CASE (EMPTY_CORNER_WALL_BELOW)
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)

           END SELECT

        ELSE IF (electron_to_add(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top left cell of the area

           SELECT CASE (c_left_top_corner_type)
              CASE (HAS_TWO_NEIGHBORS)
                 IF ((c_X_area_min-electron_to_add(k)%X).LT.(electron_to_add(k)%Y-c_Y_area_max)) THEN
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
                 END IF

              CASE (FLAT_WALL_ABOVE)
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)

              CASE (FLAT_WALL_LEFT)
                 IF (electron_to_add(k)%Y.LE.c_Y_area_max) THEN                 
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                    print '("error-6 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
                 END IF

              CASE (SURROUNDED_BY_WALL)
!                 IF ((c_X_area_min-electron_to_add(k)%X).LT.(electron_to_add(k)%Y-c_Y_area_max)) THEN
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 ELSE
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 END IF
! error
                 print '("error-7 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

              CASE (EMPTY_CORNER_WALL_LEFT)
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)

              CASE (EMPTY_CORNER_WALL_ABOVE)
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)

           END SELECT
        END IF
        CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
        CYCLE
     END IF         !### IF (electron_to_add(k)%X.LT.c_X_area_min) THEN

     IF (electron_to_add(k)%X.GT.c_X_area_max) THEN

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the electron by the period length
           electron_to_add(k)%X = electron_to_add(k)%X - L_period_X
           IF (electron_to_add(k)%Y.LT.c_Y_area_min) THEN
! particle is below the bottom side of the area
              IF (Rank_of_master_below.LT.0) THEN
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-8 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)                       
              END IF
              CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
           ELSE IF (electron_to_add(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
              IF (Rank_of_master_above.LT.0) THEN
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-9 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              ELSE
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
              END IF
              CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
           END IF
           CYCLE
        END IF

        IF ( (electron_to_add(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
           & (electron_to_add(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

           IF (Rank_of_master_right.GE.0) THEN
! right neighbor cluster exists
              CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
           ELSE
! right neighbor cluster does not exist
!              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
              print '("error-10 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF

        ELSE IF (electron_to_add(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the right bottom cell of the area

           SELECT CASE (c_right_bottom_corner_type)
              CASE (HAS_TWO_NEIGHBORS)
                 IF ((electron_to_add(k)%X-c_X_area_max).LT.(c_Y_area_min-electron_to_add(k)%Y)) THEN
                    CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
                 END IF

              CASE (FLAT_WALL_BELOW)
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
                
              CASE (FLAT_WALL_RIGHT)
                 IF (electron_to_add(k)%Y.GE.c_Y_area_min) THEN
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                    print '("error-11 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)                       
                 END IF

              CASE (SURROUNDED_BY_WALL)
!                 IF ((electron_to_add(k)%X-c_X_area_max).LT.(c_Y_area_min-electron_to_add(k)%Y)) THEN
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 ELSE
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 END IF
! error
                    print '("error-12 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

              CASE (EMPTY_CORNER_WALL_RIGHT)
                 CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)

              CASE (EMPTY_CORNER_WALL_BELOW)
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)

           END SELECT

        ELSE IF (electron_to_add(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top right cell of the area

           SELECT CASE (c_right_top_corner_type)
              CASE (HAS_TWO_NEIGHBORS)
                 IF ((electron_to_add(k)%X-c_X_area_max).LT.(electron_to_add(k)%Y-c_Y_area_max)) THEN
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
                 END IF

              CASE (FLAT_WALL_ABOVE)
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)

              CASE (FLAT_WALL_RIGHT)
                 IF (electron_to_add(k)%Y.LE.c_Y_area_max) THEN
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                    print '("error-13 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)                    
                 END IF

              CASE (SURROUNDED_BY_WALL)
!                 IF ((electron_to_add(k)%X-c_X_area_max).LT.(electron_to_add(k)%Y-c_Y_area_max)) THEN
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 ELSE
!                    CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
!                 END IF
! error
                    print '("error-14 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

              CASE (EMPTY_CORNER_WALL_RIGHT)
                 CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)

              CASE (EMPTY_CORNER_WALL_ABOVE)
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)

           END SELECT

        END IF
        CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
        CYCLE
     END IF         !### IF (electron_to_add(k)%X.GT.c_X_area_max) THEN

! since we are here, c_X_area_min <= electron_to_add(k)%Y <= c_X_area_max

     IF (electron_to_add(k)%Y.GT.c_Y_area_max) THEN

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
           IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
              CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
           ELSE
!              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
              print '("error-15 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF
           CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
           CYCLE
        END IF

        IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
           CALL ADD_ELECTRON_TO_SEND_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
        ELSE
! neighbor cluster above does not exist
           IF ((electron_to_add(k)%X.GE.(c_X_area_min+1.0_8)).AND.(electron_to_add(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
!              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
              print '("error-16 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           ELSE IF (electron_to_add(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left top corner
              IF (c_left_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
              ELSE
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-17 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
           ELSE IF (electron_to_add(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right top corner
              IF (c_right_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX,  electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
              ELSE
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_ABOVE(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-18 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
           END IF
        END IF
        CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
        CYCLE
     END IF     !### IF (electron_to_add(k)%Y.GT.c_Y_area_max) THEN

     IF (electron_to_add(k)%Y.LT.c_Y_area_min) THEN

        IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
           IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster above exists
              CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
           ELSE
!              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
              print '("error-19 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF
           CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
           CYCLE
        END IF

        IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster below exists, remove particle and prepare to send it to the neighbor below
           CALL ADD_ELECTRON_TO_SEND_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
        ELSE
! neighbor cluster below does not exist
           IF ((electron_to_add(k)%X.GE.(c_X_area_min+1.0_8)).AND.(electron_to_add(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
!              CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
              print '("error-20 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           ELSE IF (electron_to_add(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left bottom corner
              IF (c_left_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
                 CALL ADD_ELECTRON_TO_SEND_LEFT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
              ELSE
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-21 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
           ELSE IF (electron_to_add(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right bottom corner
              IF (c_right_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
                 CALL ADD_ELECTRON_TO_SEND_RIGHT(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%AX, electron_to_add(k)%AY, electron_to_add(k)%tag)
              ELSE
!                 CALL PROCESS_ELECTRON_COLL_WITH_BOUNDARY_BELOW(electron_to_add(k)%X, electron_to_add(k)%Y, electron_to_add(k)%VX, electron_to_add(k)%VY, electron_to_add(k)%VZ, electron_to_add(k)%tag)
! error
                 print '("error-22 in FIND_ALIENS_IN_ELECTRON_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
           END IF
        END IF
        CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
        CYCLE
     END IF     !### IF (electron_to_add(k)%Y.LT.c_Y_area_min) THEN

  END DO   !### DO WHILE (k.LT.N_e_to_add)

END SUBROUTINE FIND_ALIENS_IN_ELECTRON_ADD_LIST

!----------------------------------------
!
SUBROUTINE FIND_INNER_OBJECT_COLL_IN_ELECTRON_ADD_LIST

!  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ElectronParticles, ONLY : N_e_to_add, electron_to_add

  IMPLICIT NONE

  INTEGER k, n 
  REAL(8) X_move, Y_move

  IF (N_of_inner_objects.EQ.0) RETURN

! find, process, and exclude electrons which collided with inner objects
  k=0
  DO WHILE (k.LT.N_e_to_add)
     k = k+1
     DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
        IF (electron_to_add(k)%X.LE.whole_object(n)%Xmin) CYCLE
        IF (electron_to_add(k)%X.GE.whole_object(n)%Xmax) CYCLE
        IF (electron_to_add(k)%Y.LE.whole_object(n)%Ymin) CYCLE
        IF (electron_to_add(k)%Y.GE.whole_object(n)%Ymax) CYCLE
! collision detected
        CALL TRY_ELECTRON_COLL_WITH_INNER_OBJECT( electron_to_add(k)%X, &
                                                & electron_to_add(k)%Y, &
                                                & electron_to_add(k)%VX, &
                                                & electron_to_add(k)%VY, &
                                                & electron_to_add(k)%VZ, &
                                                & electron_to_add(k)%tag ) 
        CALL REMOVE_ELECTRON_FROM_ADD_LIST(k)  ! this subroutine does  N_e_to_add = N_e_to_add - 1 and k = k-1
        EXIT
     END DO
  END DO

END SUBROUTINE FIND_INNER_OBJECT_COLL_IN_ELECTRON_ADD_LIST

!----------------------------------------
!
SUBROUTINE PROCESS_ADDED_ELECTRONS

  USE ParallelOperationValues
  USE CurrentProblemValues

  USE ElectronParticles, ONLY : N_e_to_add, N_electrons, max_N_electrons, electron, electron_to_add

  USE rng_wrapper  !###???

  IMPLICIT NONE

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     real(8) ax
     real(8) ay
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR
  INTEGER k, current_N

  INTEGER random_j
  INTEGER temptag
  REAL(8) tempX, tempY, tempVX, tempVY, tempVZ, tempAX, tempAY
  
  IF (N_e_to_add.GT.(max_N_electrons-N_electrons)) THEN
! increase the size of the main electron array
     current_N = max_N_electrons
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = electron(k)%X
        bufer(k)%Y   = electron(k)%Y
        bufer(k)%VX  = electron(k)%VX
        bufer(k)%VY  = electron(k)%VY
        bufer(k)%VZ  = electron(k)%VZ

        bufer(k)%AX  = electron(k)%AX
        bufer(k)%AY  = electron(k)%AY

        bufer(k)%tag = electron(k)%tag
     END DO
     DEALLOCATE(electron, STAT=ALLOC_ERR)
     max_N_electrons = max_N_electrons + MAX(N_e_to_add-(max_N_electrons-N_electrons), max_N_electrons/10)
     ALLOCATE(electron(1:max_N_electrons), STAT=ALLOC_ERR)
     DO k = 1, current_N
        electron(k)%X   = bufer(k)%X
        electron(k)%Y   = bufer(k)%Y
        electron(k)%VX  = bufer(k)%VX
        electron(k)%VY  = bufer(k)%VY
        electron(k)%VZ  = bufer(k)%VZ

        electron(k)%AX  = bufer(k)%AX
        electron(k)%AY  = bufer(k)%AY

        electron(k)%tag = bufer(k)%tag
     END DO
     DEALLOCATE(bufer, STAT=ALLOC_ERR)
  END IF
  
  DO k = 1, N_e_to_add
     electron(k+N_electrons)%X   = electron_to_add(k)%X
     electron(k+N_electrons)%Y   = electron_to_add(k)%Y
     electron(k+N_electrons)%VX  = electron_to_add(k)%VX
     electron(k+N_electrons)%VY  = electron_to_add(k)%VY
     electron(k+N_electrons)%VZ  = electron_to_add(k)%VZ

     electron(k+N_electrons)%AX  = electron_to_add(k)%AX
     electron(k+N_electrons)%AY  = electron_to_add(k)%AY

     electron(k+N_electrons)%tag = electron_to_add(k)%tag
  END DO

! update electron counter
  N_electrons = N_electrons + N_e_to_add

!???????????? shuffle the added and the available particles ?????????????
  DO k = N_electrons - N_e_to_add + 1, N_electrons
!     random_j = MAX(1, MIN( N_electrons - N_e_to_add, INT(well_random_number() * (N_electrons-N_e_to_add))) )
     random_j = MAX(1, MIN( N_electrons, INT(well_random_number() * N_electrons)))

     IF (random_j.EQ.k) CYCLE

     tempX   = electron(random_j)%X
     tempY   = electron(random_j)%Y
     tempVX  = electron(random_j)%VX
     tempVY  = electron(random_j)%VY
     tempVZ  = electron(random_j)%VZ

     tempAX  = electron(random_j)%AX
     tempAY  = electron(random_j)%AY

     temptag = electron(random_j)%tag

     electron(random_j)%X   = electron(k)%X
     electron(random_j)%Y   = electron(k)%Y
     electron(random_j)%VX  = electron(k)%VX
     electron(random_j)%VY  = electron(k)%VY
     electron(random_j)%VZ  = electron(k)%VZ

     electron(random_j)%AX  = electron(k)%AX
     electron(random_j)%AY  = electron(k)%AY

     electron(random_j)%tag = electron(k)%tag

     electron(k)%X   = tempX
     electron(k)%Y   = tempY
     electron(k)%VX  = tempVX
     electron(k)%VY  = tempVY
     electron(k)%VZ  = tempVZ

     electron(k)%AX  = tempAX
     electron(k)%AY  = tempAY

     electron(k)%tag = temptag
  END DO

! clear counter of particles to be added to this process
  N_e_to_add = 0  

END SUBROUTINE PROCESS_ADDED_ELECTRONS
