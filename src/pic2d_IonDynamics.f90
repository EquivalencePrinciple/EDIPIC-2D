!---------------------------------------
!
SUBROUTINE ADVANCE_IONS_PRELIMINARY

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE IonParticles
  USE ClusterAndItsBoundaries
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER s, k
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
  N_ions_to_send_left  = 0
  N_ions_to_send_right = 0
  N_ions_to_send_above = 0
  N_ions_to_send_below = 0

! clear counters of particles that hit the boundary objects
  DO k = 1, N_of_boundary_and_inner_objects
     whole_object(k)%ion_hit_count(1:N_spec) = 0
     ion_colls_with_bo(k)%N_of_saved_parts = 0
  END DO

! cycle over ion species
  DO s = 1, N_spec

! cycle over particles of the ion species
     k=0
     DO WHILE (k.LT.N_ions(s))

        k = k + 1

        if ( (ion(s)%part(k)%X.lt.c_X_area_min).or. &
           & (ion(s)%part(k)%X.gt.c_X_area_max).or. &
           & (ion(s)%part(k)%Y.lt.c_Y_area_min).or. &
           & (ion(s)%part(k)%Y.gt.c_Y_area_max) ) then
           print '("Process ",i4," : Error-1 in ADVANCE_IONS : particle out of bounds xmin/xmax/ymin/ymax : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
           print '("Process ",i4," : s/k/N_ions(s) : ",i3,2x,i8,2x,i8)', Rank_of_process, s, k, N_ions(s)
           print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        end if

        IF (ions_sense_Ez) THEN
           E_Z = Ez(ion(s)%part(k)%X, ion(s)%part(k)%Y)
        ELSE
           E_Z = 0.0_8
        END IF
        
        IF (ions_sense_magnetic_field) THEN
! magnetic field accounted for
! calculate magnetic field factors   !##### MODIFY FOR IONS ########

           alfa_x = QM2s(s) * Bx(ion(s)%part(k)%X, ion(s)%part(k)%Y)
           alfa_y = QM2s(s) * By(ion(s)%part(k)%X, ion(s)%part(k)%Y)
           alfa_z = QM2s(s) * Bz(ion(s)%part(k)%X, ion(s)%part(k)%Y)

           alfa_x2 = alfa_x**2
           alfa_y2 = alfa_y**2
           alfa_z2 = alfa_z**2

           theta2 = alfa_x2 + alfa_y2 + alfa_z2
           invtheta = 1.0_8 / (1.0_8 + theta2)

!    matrix K, same as R in Gibbons & Hewett; A_inverse = (I + K)/2 = (I + R)/2
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
           vec1 = 0.5_8 * ion(s)%part(k)%AX
           vec2 = 0.5_8 * ion(s)%part(k)%AY
           vec3 = E_Z * QMs(s)                  ! QMs(s) = Qs(s) / Ms(s)

           VX_temp = ion(s)%part(k)%VX
           VY_temp = ion(s)%part(k)%VY
           VZ_temp = ion(s)%part(k)%VZ

           ion(s)%part(k)%VX = K11 * VX_temp + K12 * VY_temp + K13 * VZ_temp + A11 * vec1 + A12 * vec2 + A13 * vec3
           ion(s)%part(k)%VY = K21 * VX_temp + K22 * VY_temp + K23 * VZ_temp + A21 * vec1 + A22 * vec2 + A23 * vec3
           ion(s)%part(k)%VZ = K31 * VX_temp + K32 * VY_temp + K33 * VZ_temp + A31 * vec1 + A32 * vec2 + A33 * vec3

        ELSE   !### IF (ions_sense_mangetic_field) THEN
! magnetic field effects omitted
! K=R=I, A_inverse = I

! predictive move
! Gibbons & Hewett Eq.(2.6):             
           ion(s)%part(k)%VX = ion(s)%part(k)%VX + 0.5_8 * ion(s)%part(k)%AX
           ion(s)%part(k)%VY = ion(s)%part(k)%VY + 0.5_8 * ion(s)%part(k)%AY
           ion(s)%part(k)%VZ = ion(s)%part(k)%VZ + E_Z * QMs(s)              ! QMs(s) = Qs(s) / Ms(s)

        END IF   !### IF (ions_sense_mangetic_field) THEN

! coordinate advance

        ion(s)%part(k)%X = ion(s)%part(k)%X + ion(s)%part(k)%VX
        ion(s)%part(k)%Y = ion(s)%part(k)%Y + ion(s)%part(k)%VY

! a particle crossed symmetry plane, reflect it
        IF (symmetry_plane_X_left) THEN
           IF (ion(s)%part(k)%X.LT.c_X_area_min) THEN
              ion(s)%part(k)%X = MAX(c_X_area_min, c_X_area_min + c_X_area_min - ion(s)%part(k)%X)
              ion(s)%part(k)%VX = -ion(s)%part(k)%VX
           END IF
        END IF

! check whether a collision with an inner object occurred
        collision_with_inner_object_occurred = .FALSE.
        DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF (ion(s)%part(k)%X.LE.whole_object(n)%Xmin) CYCLE
           IF (ion(s)%part(k)%X.GE.whole_object(n)%Xmax) CYCLE
           IF (ion(s)%part(k)%Y.LE.whole_object(n)%Ymin) CYCLE
           IF (ion(s)%part(k)%Y.GE.whole_object(n)%Ymax) CYCLE
! collision detected
           CALL TRY_ION_COLL_WITH_INNER_OBJECT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
           CALL REMOVE_ION(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           collision_with_inner_object_occurred = .TRUE.
           EXIT
        END DO

        IF (collision_with_inner_object_occurred) CYCLE

! most probable situation when the particle remains inside the area
        IF ( (ion(s)%part(k)%X.GE.c_X_area_min) .AND. &
           & (ion(s)%part(k)%X.LE.c_X_area_max) .AND. &
           & (ion(s)%part(k)%Y.GE.c_Y_area_min) .AND. &
           & (ion(s)%part(k)%Y.LE.c_Y_area_max) ) CYCLE

! since we are here, a particle did not collide with an inner object but crossed an area boundary
! note, in a periodic system the particle still may collide with an inner object after transfer to the other domain
! therefore it is still necessary to check add list for collisions after process receives particles from neighbors

        CALL PROCESS_ION_LEAVING_CLUSTER(s, k)

     END DO  ! end of cycle over particles of ion species

  END DO ! end of cycle over ion species

END SUBROUTINE ADVANCE_IONS_PRELIMINARY

!---------------------------------------
!
SUBROUTINE UPDATE_ION_ACCELERATIONS

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE IonParticles
  USE ClusterAndItsBoundaries
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER s, k
  INTEGER i, j, ix, jx, iy, jy
  REAL(8) ax_ip1, ax_i, ax_jp1, ax_j !meaning of these coefficients has changed
  REAL(8) ay_ip1, ay_i, ay_jp1, ay_j
  REAL(8) X_Ex, Y_Ey
  REAL(8) E_X, E_Y

! cycle over ion species
  DO s = 1, N_spec
! cycle over particles of the ion species

     DO k = 1, N_ions(s)

        if ( (ion(s)%part(k)%X.lt.c_X_area_min).or. &
           & (ion(s)%part(k)%X.gt.c_X_area_max).or. &
           & (ion(s)%part(k)%Y.lt.c_Y_area_min).or. &
           & (ion(s)%part(k)%Y.gt.c_Y_area_max) ) then
           print '("Process ",i4," : Error-1 in ADVANCE_IONS : particle out of bounds xmin/xmax/ymin/ymax : ",4(2x,e14.7))', Rank_of_process, c_X_area_min, c_X_area_max, c_Y_area_min, c_Y_area_max
           print '("Process ",i4," : s/k/N_ions(s) : ",i3,2x,i8,2x,i8)', Rank_of_process, s, k, N_ions(s)
           print '("Process ",i4," : x/y/vx/vy/vz/tag : ",5(2x,e14.7),2x,i4)', Rank_of_process, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag
           CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
        end if

! interpolate electric field

        i = FLOOR(ion(s)%part(k)%X)
        j = FLOOR(ion(s)%part(k)%Y)

        IF (ion(s)%part(k)%X.EQ.c_X_area_max) i = c_indx_x_max-1
        IF (ion(s)%part(k)%Y.EQ.c_Y_area_max) j = c_indx_y_max-1

        ix = FLOOR(ion(s)%part(k)%X - 0.5_8)   ! Ex x-indices limits: c_indx_x_min - 1 : c_indx_x_max
        jx = j
        
        iy = i
        jy = FLOOR(ion(s)%part(k)%Y - 0.5_8)   ! Ey y-indices limits: c_indx_y_min - 1 : c_indx_y_max

        X_Ex = DBLE(ix) + 0.5_8             ! x-coordinate of Ex node ix
        Y_Ey = DBLE(jy) + 0.5_8             ! y-coordinate of Ey node jy

        ax_ip1 = ion(s)%part(k)%X - X_Ex
        ax_i   = 1.0_8 - ax_ip1

        ax_jp1 = ion(s)%part(k)%Y - DBLE(jx)
        ax_j   = 1.0_8 - ax_jp1

        ay_ip1 = ion(s)%part(k)%X - DBLE(iy)
        ay_i   = 1.0_8 - ay_ip1

        ay_jp1 = ion(s)%part(k)%Y - Y_Ey
        ay_j   = 1.0_8 - ay_jp1

! interpolate Ex and Ey:
        E_X = EX(ix, jx) * ax_i * ax_j + EX(ix+1, jx) * ax_ip1 * ax_j + EX(ix, jx+1) * ax_i * ax_jp1 + EX(ix+1, jx+1) * ax_ip1 * ax_jp1
        E_Y = EY(iy, jy) * ay_i * ay_j + EY(iy+1, jy) * ay_ip1 * ay_j + EY(iy, jy+1) * ay_i * ay_jp1 + EY(iy+1, jy+1) * ay_ip1 * ay_jp1

        ion(s)%part(k)%AX = 0.5_8 * ion(s)%part(k)%AX + QM2s(s) * E_X   ! QM2s(s) = 0.5_8 * Qs(s) / Ms(s)
        ion(s)%part(k)%AY = 0.5_8 * ion(s)%part(k)%AY + QM2s(s) * E_Y

     END DO
  END DO

END SUBROUTINE UPDATE_ION_ACCELERATIONS

!-----------------------------------------------------------------------------------------------
!
SUBROUTINE PROCESS_ION_LEAVING_CLUSTER(s, k)

  USE ParallelOperationValues
  USE IonParticles
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER, INTENT(IN) :: s
  INTEGER, INTENT(INOUT) :: k

  IF (ion(s)%part(k)%X.LT.c_X_area_min) THEN

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the ion by the period length
        ion(s)%part(k)%X = ion(s)%part(k)%X + L_period_X
        IF (ion(s)%part(k)%Y.LT.c_Y_area_min) THEN
! particle is somewhere near the left bottom cell of the area
           IF (Rank_of_master_below.LT.0) THEN
              CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
           ELSE
              CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)                       
           END IF
           CALL REMOVE_ION(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
        ELSE IF (ion(s)%part(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
           IF (Rank_of_master_above.LT.0) THEN
              CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
           ELSE
              CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
           END IF
           CALL REMOVE_ION(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
        END IF
        RETURN
     END IF

     IF ( (ion(s)%part(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
       & (ion(s)%part(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

        IF (Rank_of_master_left.GE.0) THEN
! left neighbor cluster exists
           CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
        ELSE
! left neighbor cluster does not exist
           CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX,  ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)   ! left
        END IF

     ELSE IF (ion(s)%part(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the left bottom cell of the area

        SELECT CASE (c_left_bottom_corner_type)
           CASE (HAS_TWO_NEIGHBORS)
              IF ((c_X_area_min-ion(s)%part(k)%X).LT.(c_Y_area_min-ion(s)%part(k)%Y)) THEN
                 CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
              ELSE
                 CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
              END IF

           CASE (FLAT_WALL_BELOW)
              CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)

           CASE (FLAT_WALL_LEFT)
              IF (ion(s)%part(k)%Y.GE.c_Y_area_min) THEN                 
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE
                 CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)                       
              END IF

           CASE (SURROUNDED_BY_WALL)
              IF ((c_X_area_min-ion(s)%part(k)%X).LT.(c_Y_area_min-ion(s)%part(k)%Y)) THEN
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              END IF

           CASE (EMPTY_CORNER_WALL_LEFT)
              CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
               
           CASE (EMPTY_CORNER_WALL_BELOW)
              CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)

        END SELECT

     ELSE IF (ion(s)%part(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top left cell of the area

        SELECT CASE (c_left_top_corner_type)
           CASE (HAS_TWO_NEIGHBORS)
              IF ((c_X_area_min-ion(s)%part(k)%X).LT.(ion(s)%part(k)%Y-c_Y_area_max)) THEN
                 CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
              ELSE
                 CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
              END IF

           CASE (FLAT_WALL_ABOVE)
              CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)

           CASE (FLAT_WALL_LEFT)
              IF (ion(s)%part(k)%Y.LE.c_Y_area_max) THEN                 
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE
                 CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
              END IF

           CASE (SURROUNDED_BY_WALL)
              IF ((c_X_area_min-ion(s)%part(k)%X).LT.(ion(s)%part(k)%Y-c_Y_area_max)) THEN
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              END IF

           CASE (EMPTY_CORNER_WALL_LEFT)
              CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)

           CASE (EMPTY_CORNER_WALL_ABOVE)
              CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)

        END SELECT

     ELSE
! ERROR, we shouldn't be here
        PRINT '("ERROR-1 in ADVANCE_IONS: we should not be here")'
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     END IF

     CALL REMOVE_ION(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
     RETURN

  END IF

  IF (ion(s)%part(k)%X.GT.c_X_area_max) THEN

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the ion by the period length
        ion(s)%part(k)%X = ion(s)%part(k)%X - L_period_X
        IF (ion(s)%part(k)%Y.LT.c_Y_area_min) THEN
! particle is somewhere near the left bottom cell of the area
           IF (Rank_of_master_below.LT.0) THEN
              CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
           ELSE
              CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)                       
           END IF
           CALL REMOVE_ION(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
        ELSE IF (ion(s)%part(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
           IF (Rank_of_master_above.LT.0) THEN
              CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
           ELSE
              CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
           END IF
           CALL REMOVE_ION(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
        END IF
        RETURN
     END IF
           
     IF ( (ion(s)%part(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
        & (ion(s)%part(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

        IF (Rank_of_master_right.GE.0) THEN
! right neighbor cluster exists
           CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
        ELSE
! right neighbor cluster does not exist
           CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
        END IF

     ELSE IF (ion(s)%part(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the right bottom cell of the area

        SELECT CASE (c_right_bottom_corner_type)
           CASE (HAS_TWO_NEIGHBORS)
              IF ((ion(s)%part(k)%X-c_X_area_max).LT.(c_Y_area_min-ion(s)%part(k)%Y)) THEN
                 CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
              ELSE
                 CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
              END IF

           CASE (FLAT_WALL_BELOW)
              CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
                
           CASE (FLAT_WALL_RIGHT)
              IF (ion(s)%part(k)%Y.GE.c_Y_area_min) THEN
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE
                 CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)                       
              END IF

           CASE (SURROUNDED_BY_WALL)
              IF ((ion(s)%part(k)%X-c_X_area_max).LT.(c_Y_area_min-ion(s)%part(k)%Y)) THEN
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              END IF

           CASE (EMPTY_CORNER_WALL_RIGHT)
              CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)

           CASE (EMPTY_CORNER_WALL_BELOW)
              CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)

        END SELECT

     ELSE IF (ion(s)%part(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top right cell of the area

        SELECT CASE (c_right_top_corner_type)
           CASE (HAS_TWO_NEIGHBORS)
              IF ((ion(s)%part(k)%X-c_X_area_max).LT.(ion(s)%part(k)%Y-c_Y_area_max)) THEN
                 CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
              ELSE
                 CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
              END IF

           CASE (FLAT_WALL_ABOVE)
              CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)

           CASE (FLAT_WALL_RIGHT)
              IF (ion(s)%part(k)%Y.LE.c_Y_area_max) THEN
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE
                 CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)                    
              END IF

           CASE (SURROUNDED_BY_WALL)
              IF ((ion(s)%part(k)%X-c_X_area_max).LT.(ion(s)%part(k)%Y-c_Y_area_max)) THEN
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              ELSE
                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
              END IF

           CASE (EMPTY_CORNER_WALL_RIGHT)
              CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)

           CASE (EMPTY_CORNER_WALL_ABOVE)
              CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)

        END SELECT

     ELSE
! ERROR, we shouldn't be here
        PRINT '("ERROR-2 in ADVANCE_IONS: we should not be here")'
        CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
     END IF

     CALL REMOVE_ION(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
     RETURN

  END IF

! since we are here, c_X_area_min <= ion(s)&part(k)%Y <= c_X_area_max

  IF (ion(s)%part(k)%Y.GT.c_Y_area_max) THEN

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
        IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
           CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
        ELSE
           CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
        END IF
        CALL REMOVE_ION(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
        RETURN
     END IF

     IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
        CALL ADD_ION_TO_SEND_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
     ELSE
! neighbor cluster above does not exist
        IF ((ion(s)%part(k)%X.GE.(c_X_area_min+1.0_8)).AND.(ion(s)%part(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
           CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
        ELSE IF (ion(s)%part(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left top corner
           IF (c_left_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
              CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
           ELSE
              CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
           END IF
        ELSE IF (ion(s)%part(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right top corner
           IF (c_right_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
              CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX,  ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
           ELSE
              CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
           END IF
        END IF
     END IF
     CALL REMOVE_ION(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
     RETURN
  END IF

  IF (ion(s)%part(k)%Y.LT.c_Y_area_min) THEN

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
        IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster above exists
           CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
        ELSE
           CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
        END IF
        CALL REMOVE_ION(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
        RETURN
     END IF

     IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster below exists, remove particle and prepare to send it to the neighbor below
        CALL ADD_ION_TO_SEND_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
     ELSE
! neighbor cluster below does not exist
        IF ((ion(s)%part(k)%X.GE.(c_X_area_min+1.0_8)).AND.(ion(s)%part(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
           CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
        ELSE IF (ion(s)%part(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left bottom corner
           IF (c_left_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
              CALL ADD_ION_TO_SEND_LEFT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
           ELSE
              CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
           END IF
        ELSE IF (ion(s)%part(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right bottom corner
           IF (c_right_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
              CALL ADD_ION_TO_SEND_RIGHT(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%AX, ion(s)%part(k)%AY, ion(s)%part(k)%tag)
           ELSE
              CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion(s)%part(k)%X, ion(s)%part(k)%Y, ion(s)%part(k)%VX, ion(s)%part(k)%VY, ion(s)%part(k)%VZ, ion(s)%part(k)%tag)
           END IF
        END IF
     END IF
     CALL REMOVE_ION(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
     RETURN
  END IF

END SUBROUTINE PROCESS_ION_LEAVING_CLUSTER

!----------------------------------------
!
SUBROUTINE REMOVE_ION(s, k)

  USE ParallelOperationValues
  USE IonParticles

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER, INTENT(IN) :: s
  INTEGER, INTENT(INOUT) :: k

  IF ((s.LT.1).OR.(s.GT.N_spec)) THEN
     PRINT '("Process ",i6," : ERROR-1 in REMOVE_ION : index s invalid")', Rank_of_process
     PRINT '("Process ",i6," : s = ",i3," k= ", i7," N_ions(s)= ",i7)', Rank_of_process, s, k, N_ions(s)
     PRINT '("Process ",i6," : PROGRAM TERMINATED")', Rank_of_process
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

  IF ((k.LT.1).OR.(k.GT.N_ions(s))) THEN
     PRINT '("Process ",i6," : ERROR-2 in REMOVE_ION : index k invalid")', Rank_of_process
     PRINT '("Process ",i6," : s = ",i3," k= ", i7," N_ions(s)= ",i7)', Rank_of_process, s, k, N_ions(s)
     PRINT '("Process ",i6," : PROGRAM TERMINATED")', Rank_of_process
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

  IF (k.LT.N_ions(s)) THEN

     ion(s)%part(k)%X   = ion(s)%part(N_ions(s))%X
     ion(s)%part(k)%Y   = ion(s)%part(N_ions(s))%Y
     ion(s)%part(k)%VX  = ion(s)%part(N_ions(s))%VX
     ion(s)%part(k)%VY  = ion(s)%part(N_ions(s))%VY
     ion(s)%part(k)%VZ  = ion(s)%part(N_ions(s))%VZ
     ion(s)%part(k)%AX  = ion(s)%part(N_ions(s))%AX
     ion(s)%part(k)%AY  = ion(s)%part(N_ions(s))%AY
     ion(s)%part(k)%tag = ion(s)%part(N_ions(s))%tag

  END IF

  N_ions(s) = N_ions(s) - 1
  k = k - 1                  ! to ensure that the new k-particle is processed

END SUBROUTINE REMOVE_ION

!------------------------------------------
!
SUBROUTINE ADD_ION_TO_SEND_LEFT(s, x, y, vx, vy, vz, ax, ay, tag)

  USE IonParticles , ONLY : N_ions_to_send_left, max_N_ions_to_send_left, ion_to_send_left
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_X_left, L_period_x

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz, ax, ay
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     real(8) AX
     real(8) AY
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_ions_to_send_left(s) = N_ions_to_send_left(s) + 1

  IF (N_ions_to_send_left(s).GT.max_N_ions_to_send_left(s)) THEN
! increase the size of the list array
     current_N = max_N_ions_to_send_left(s)
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = ion_to_send_left(s)%part(k)%X
        bufer(k)%Y   = ion_to_send_left(s)%part(k)%Y
        bufer(k)%VX  = ion_to_send_left(s)%part(k)%VX
        bufer(k)%VY  = ion_to_send_left(s)%part(k)%VY
        bufer(k)%VZ  = ion_to_send_left(s)%part(k)%VZ
        bufer(k)%AX  = ion_to_send_left(s)%part(k)%AX
        bufer(k)%AY  = ion_to_send_left(s)%part(k)%AY
        bufer(k)%tag = ion_to_send_left(s)%part(k)%tag
     END DO
     IF (ALLOCATED(ion_to_send_left(s)%part)) THEN 
        DEALLOCATE(ion_to_send_left(s)%part, STAT=DEALLOC_ERR)
        !NULLIFY(ion_to_send_left(s)%part)
     END IF
     max_N_ions_to_send_left(s) = max_N_ions_to_send_left(s) + MAX(50, max_N_ions_to_send_left(s)/10)
     ALLOCATE(ion_to_send_left(s)%part(1:max_N_ions_to_send_left(s)), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        ion_to_send_left(s)%part(k)%X   = bufer(k)%X
        ion_to_send_left(s)%part(k)%Y   = bufer(k)%Y
        ion_to_send_left(s)%part(k)%VX  = bufer(k)%VX
        ion_to_send_left(s)%part(k)%VY  = bufer(k)%VY
        ion_to_send_left(s)%part(k)%VZ  = bufer(k)%VZ
        ion_to_send_left(s)%part(k)%AX  = bufer(k)%AX
        ion_to_send_left(s)%part(k)%AY  = bufer(k)%AY
        ion_to_send_left(s)%part(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  IF (periodic_boundary_X_left) THEN
     ion_to_send_left(s)%part(N_ions_to_send_left(s))%X = x + L_period_X
  ELSE
     ion_to_send_left(s)%part(N_ions_to_send_left(s))%X = x
  END IF
  ion_to_send_left(s)%part(N_ions_to_send_left(s))%Y   = y
  ion_to_send_left(s)%part(N_ions_to_send_left(s))%VX  = vx
  ion_to_send_left(s)%part(N_ions_to_send_left(s))%VY  = vy
  ion_to_send_left(s)%part(N_ions_to_send_left(s))%VZ  = vz
  ion_to_send_left(s)%part(N_ions_to_send_left(s))%AX  = ax
  ion_to_send_left(s)%part(N_ions_to_send_left(s))%AY  = ay
  ion_to_send_left(s)%part(N_ions_to_send_left(s))%tag = tag

END SUBROUTINE ADD_ION_TO_SEND_LEFT

!------------------------------------------
!
SUBROUTINE ADD_ION_TO_SEND_RIGHT(s, x, y, vx, vy, vz, ax, ay, tag)

  USE IonParticles, ONLY : N_ions_to_send_right, max_N_ions_to_send_right, ion_to_send_right
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_X_right, L_period_x

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz, ax, ay
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     real(8) AX
     real(8) AY
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_ions_to_send_right(s) = N_ions_to_send_right(s) + 1

  IF (N_ions_to_send_right(s).GT.max_N_ions_to_send_right(s)) THEN     !######## allocate somewhere ion_to_send_right first ######
! increase the size of the list array
     current_N = max_N_ions_to_send_right(s)
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = ion_to_send_right(s)%part(k)%X
        bufer(k)%Y   = ion_to_send_right(s)%part(k)%Y
        bufer(k)%VX  = ion_to_send_right(s)%part(k)%VX
        bufer(k)%VY  = ion_to_send_right(s)%part(k)%VY
        bufer(k)%VZ  = ion_to_send_right(s)%part(k)%VZ
        bufer(k)%AX  = ion_to_send_right(s)%part(k)%AX
        bufer(k)%AY  = ion_to_send_right(s)%part(k)%AY
        bufer(k)%tag = ion_to_send_right(s)%part(k)%tag
     END DO
     IF (ALLOCATED(ion_to_send_right(s)%part)) THEN 
        DEALLOCATE(ion_to_send_right(s)%part, STAT=DEALLOC_ERR)
        !NULLIFY(ion_to_send_right(s)%part)
     END IF
     max_N_ions_to_send_right(s) = max_N_ions_to_send_right(s) + MAX(50, max_N_ions_to_send_right(s)/10)
     ALLOCATE(ion_to_send_right(s)%part(1:max_N_ions_to_send_right(s)), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        ion_to_send_right(s)%part(k)%X   = bufer(k)%X
        ion_to_send_right(s)%part(k)%Y   = bufer(k)%Y
        ion_to_send_right(s)%part(k)%VX  = bufer(k)%VX
        ion_to_send_right(s)%part(k)%VY  = bufer(k)%VY
        ion_to_send_right(s)%part(k)%VZ  = bufer(k)%VZ
        ion_to_send_right(s)%part(k)%AX  = bufer(k)%AX
        ion_to_send_right(s)%part(k)%AY  = bufer(k)%AY
        ion_to_send_right(s)%part(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  IF (periodic_boundary_X_right) THEN
     ion_to_send_right(s)%part(N_ions_to_send_right(s))%X = x - L_period_X
  ELSE
     ion_to_send_right(s)%part(N_ions_to_send_right(s))%X = x
  END IF
  ion_to_send_right(s)%part(N_ions_to_send_right(s))%Y   = y
  ion_to_send_right(s)%part(N_ions_to_send_right(s))%VX  = vx
  ion_to_send_right(s)%part(N_ions_to_send_right(s))%VY  = vy
  ion_to_send_right(s)%part(N_ions_to_send_right(s))%VZ  = vz
  ion_to_send_right(s)%part(N_ions_to_send_right(s))%AX  = ax
  ion_to_send_right(s)%part(N_ions_to_send_right(s))%AY  = ay
  ion_to_send_right(s)%part(N_ions_to_send_right(s))%tag = tag

END SUBROUTINE ADD_ION_TO_SEND_RIGHT

!------------------------------------------
!
SUBROUTINE ADD_ION_TO_SEND_ABOVE(s, x, y, vx, vy, vz, ax, ay, tag)

  USE IonParticles, ONLY : N_ions_to_send_above, max_N_ions_to_send_above, ion_to_send_above
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_Y_above, L_period_y

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz, ax, ay
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     real(8) AX
     real(8) AY
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_ions_to_send_above(s) = N_ions_to_send_above(s) + 1

  IF (N_ions_to_send_above(s).GT.max_N_ions_to_send_above(s)) THEN
! increase the size of the list array
     current_N = max_N_ions_to_send_above(s)
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = ion_to_send_above(s)%part(k)%X
        bufer(k)%Y   = ion_to_send_above(s)%part(k)%Y
        bufer(k)%VX  = ion_to_send_above(s)%part(k)%VX
        bufer(k)%VY  = ion_to_send_above(s)%part(k)%VY
        bufer(k)%VZ  = ion_to_send_above(s)%part(k)%VZ
        bufer(k)%AX  = ion_to_send_above(s)%part(k)%AX
        bufer(k)%AY  = ion_to_send_above(s)%part(k)%AY
        bufer(k)%tag = ion_to_send_above(s)%part(k)%tag
     END DO
     IF (ALLOCATED(ion_to_send_above(s)%part)) THEN 
        DEALLOCATE(ion_to_send_above(s)%part, STAT=DEALLOC_ERR)
        !NULLIFY(ion_to_send_above(s)%part)
     END IF
     max_N_ions_to_send_above(s) = max_N_ions_to_send_above(s) + MAX(50, max_N_ions_to_send_above(s)/10)
     ALLOCATE(ion_to_send_above(s)%part(1:max_N_ions_to_send_above(s)), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        ion_to_send_above(s)%part(k)%X   = bufer(k)%X
        ion_to_send_above(s)%part(k)%Y   = bufer(k)%Y
        ion_to_send_above(s)%part(k)%VX  = bufer(k)%VX
        ion_to_send_above(s)%part(k)%VY  = bufer(k)%VY
        ion_to_send_above(s)%part(k)%VZ  = bufer(k)%VZ
        ion_to_send_above(s)%part(k)%AX  = bufer(k)%AX
        ion_to_send_above(s)%part(k)%AY  = bufer(k)%AY
        ion_to_send_above(s)%part(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  ion_to_send_above(s)%part(N_ions_to_send_above(s))%X   = x
  IF (periodic_boundary_Y_above) THEN
     ion_to_send_above(s)%part(N_ions_to_send_above(s))%Y = y - L_period_y
  ELSE
     ion_to_send_above(s)%part(N_ions_to_send_above(s))%Y = y
  END IF
  ion_to_send_above(s)%part(N_ions_to_send_above(s))%VX  = vx
  ion_to_send_above(s)%part(N_ions_to_send_above(s))%VY  = vy
  ion_to_send_above(s)%part(N_ions_to_send_above(s))%VZ  = vz
  ion_to_send_above(s)%part(N_ions_to_send_above(s))%AX  = ax
  ion_to_send_above(s)%part(N_ions_to_send_above(s))%AY  = ay
  ion_to_send_above(s)%part(N_ions_to_send_above(s))%tag = tag

END SUBROUTINE ADD_ION_TO_SEND_ABOVE

!------------------------------------------
!
SUBROUTINE ADD_ION_TO_SEND_BELOW(s, x, y, vx, vy, vz, ax, ay, tag)

  USE IonParticles, ONLY : N_ions_to_send_below, max_N_ions_to_send_below, ion_to_send_below
  USE ClusterAndItsBoundaries, ONLY : periodic_boundary_Y_below, L_period_y

!use CurrentProblemValues, only : T_cntr
!USE ParallelOperationValues, only : Rank_of_process

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz, ax, ay
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     real(8) AX
     real(8) AY
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_ions_to_send_below(s) = N_ions_to_send_below(s) + 1

  IF (N_ions_to_send_below(s).GT.max_N_ions_to_send_below(s)) THEN
! increase the size of the list array
     current_N = max_N_ions_to_send_below(s)
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = ion_to_send_below(s)%part(k)%X
        bufer(k)%Y   = ion_to_send_below(s)%part(k)%Y
        bufer(k)%VX  = ion_to_send_below(s)%part(k)%VX
        bufer(k)%VY  = ion_to_send_below(s)%part(k)%VY
        bufer(k)%VZ  = ion_to_send_below(s)%part(k)%VZ
        bufer(k)%AX  = ion_to_send_below(s)%part(k)%AX
        bufer(k)%AY  = ion_to_send_below(s)%part(k)%AY
        bufer(k)%tag = ion_to_send_below(s)%part(k)%tag
     END DO
     IF (ALLOCATED(ion_to_send_below(s)%part)) THEN 
        DEALLOCATE(ion_to_send_below(s)%part, STAT=DEALLOC_ERR)
        !NULLIFY(ion_to_send_below(s)%part)
     END IF
     max_N_ions_to_send_below(s) = max_N_ions_to_send_below(s) + MAX(50, max_N_ions_to_send_below(s)/10)
     ALLOCATE(ion_to_send_below(s)%part(1:max_N_ions_to_send_below(s)), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        ion_to_send_below(s)%part(k)%X   = bufer(k)%X
        ion_to_send_below(s)%part(k)%Y   = bufer(k)%Y
        ion_to_send_below(s)%part(k)%VX  = bufer(k)%VX
        ion_to_send_below(s)%part(k)%VY  = bufer(k)%VY
        ion_to_send_below(s)%part(k)%VZ  = bufer(k)%VZ
        ion_to_send_below(s)%part(k)%AX  = bufer(k)%AX
        ion_to_send_below(s)%part(k)%AY  = bufer(k)%AY
        ion_to_send_below(s)%part(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  ion_to_send_below(s)%part(N_ions_to_send_below(s))%X   = x
  IF (periodic_boundary_Y_below) THEN
     ion_to_send_below(s)%part(N_ions_to_send_below(s))%Y = y + L_period_y
  ELSE
     ion_to_send_below(s)%part(N_ions_to_send_below(s))%Y = y
  END IF
  ion_to_send_below(s)%part(N_ions_to_send_below(s))%VX  = vx
  ion_to_send_below(s)%part(N_ions_to_send_below(s))%VY  = vy
  ion_to_send_below(s)%part(N_ions_to_send_below(s))%VZ  = vz
  ion_to_send_below(s)%part(N_ions_to_send_below(s))%AX  = ax
  ion_to_send_below(s)%part(N_ions_to_send_below(s))%AY  = ay
  ion_to_send_below(s)%part(N_ions_to_send_below(s))%tag = tag

END SUBROUTINE ADD_ION_TO_SEND_BELOW

!----------------------------------------
!
SUBROUTINE ADD_ION_TO_ADD_LIST(s, x, y, vx, vy, vz, ax, ay, tag)

  USE IonParticles, ONLY : N_ions_to_add, max_N_ions_to_add, ion_to_add

  IMPLICIT NONE

  INTEGER s
  REAL(8) x, y, vx, vy, vz, ax, ay
  INTEGER tag

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     real(8) AX
     real(8) AY
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER k, current_N
  
  N_ions_to_add(s) = N_ions_to_add(s) + 1

  IF (N_ions_to_add(s).GT.max_N_ions_to_add(s)) THEN
! increase the size of the list array
     current_N = max_N_ions_to_add(s)
     ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
     DO k = 1, current_N
        bufer(k)%X   = ion_to_add(s)%part(k)%X
        bufer(k)%Y   = ion_to_add(s)%part(k)%Y
        bufer(k)%VX  = ion_to_add(s)%part(k)%VX
        bufer(k)%VY  = ion_to_add(s)%part(k)%VY
        bufer(k)%VZ  = ion_to_add(s)%part(k)%VZ
        bufer(k)%AX  = ion_to_add(s)%part(k)%AX
        bufer(k)%AY  = ion_to_add(s)%part(k)%AY
        bufer(k)%tag = ion_to_add(s)%part(k)%tag
     END DO
     IF (ALLOCATED(ion_to_add(s)%part)) DEALLOCATE(ion_to_add(s)%part, STAT=DEALLOC_ERR)
     max_N_ions_to_add(s) = max_N_ions_to_add(s) + MAX(50, max_N_ions_to_add(s)/10)
     ALLOCATE(ion_to_add(s)%part(1:max_N_ions_to_add(s)), STAT=DEALLOC_ERR)
     DO k = 1, current_N
        ion_to_add(s)%part(k)%X   = bufer(k)%X
        ion_to_add(s)%part(k)%Y   = bufer(k)%Y
        ion_to_add(s)%part(k)%VX  = bufer(k)%VX
        ion_to_add(s)%part(k)%VY  = bufer(k)%VY
        ion_to_add(s)%part(k)%VZ  = bufer(k)%VZ
        ion_to_add(s)%part(k)%AX  = bufer(k)%AX
        ion_to_add(s)%part(k)%AY  = bufer(k)%AY
        ion_to_add(s)%part(k)%tag = bufer(k)%tag
     END DO
     IF (ALLOCATED(bufer)) DEALLOCATE(bufer, STAT=DEALLOC_ERR)
  END IF
  
  ion_to_add(s)%part(N_ions_to_add(s))%X   = x
  ion_to_add(s)%part(N_ions_to_add(s))%Y   = y
  ion_to_add(s)%part(N_ions_to_add(s))%VX  = vx
  ion_to_add(s)%part(N_ions_to_add(s))%VY  = vy
  ion_to_add(s)%part(N_ions_to_add(s))%VZ  = vz
  ion_to_add(s)%part(N_ions_to_add(s))%AX  = ax
  ion_to_add(s)%part(N_ions_to_add(s))%AY  = ay
  ion_to_add(s)%part(N_ions_to_add(s))%tag = tag

END SUBROUTINE ADD_ION_TO_ADD_LIST

!----------------------------------------
!
SUBROUTINE REMOVE_ION_FROM_ADD_LIST(s, k)

  USE ParallelOperationValues
  USE IonParticles

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER, INTENT(IN)    :: s
  INTEGER, INTENT(INOUT) :: k

  IF ((k.LT.1).OR.(k.GT.N_ions_to_add(s))) THEN
     PRINT '("Process ",i6," : ERROR in REMOVE_ION_FROM_ADD_LIST : index k invalid")', Rank_of_process
     PRINT '("Process ",i6," : k= ", i7," N_ions_to_add(",i2,")= ",i7)', Rank_of_process, k, s, N_ions_to_add(s)
     PRINT '("Process ",i6," : PROGRAM TERMINATED")', Rank_of_process
     CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
  END IF

  IF (k.LT.N_ions_to_add(s)) THEN

     ion_to_add(s)%part(k)%X   = ion_to_add(s)%part(N_ions_to_add(s))%X
     ion_to_add(s)%part(k)%Y   = ion_to_add(s)%part(N_ions_to_add(s))%Y
     ion_to_add(s)%part(k)%VX  = ion_to_add(s)%part(N_ions_to_add(s))%VX
     ion_to_add(s)%part(k)%VY  = ion_to_add(s)%part(N_ions_to_add(s))%VY
     ion_to_add(s)%part(k)%VZ  = ion_to_add(s)%part(N_ions_to_add(s))%VZ
     ion_to_add(s)%part(k)%AX  = ion_to_add(s)%part(N_ions_to_add(s))%AX
     ion_to_add(s)%part(k)%AY  = ion_to_add(s)%part(N_ions_to_add(s))%AY
     ion_to_add(s)%part(k)%tag = ion_to_add(s)%part(N_ions_to_add(s))%tag

  END IF

  N_ions_to_add(s) = N_ions_to_add(s) - 1
  k = k-1                  ! to ensure that the new k-particle is processed

!print '("Process ",i4," called REMOVE_ION_FROM_ADD_LIST(k), T_cntr= ",i7," k= ",i4)', Rank_of_process, T_cntr, k

END SUBROUTINE REMOVE_ION_FROM_ADD_LIST

!----------------------------------------
! This subroutine is called after ADVANCE_IONS, before exchange of ions takes place.
! It removes particles which do not belong to the cluster domain.
! Such particles may appear after emission from surface of inner objects.
! It is expected that the emission is never directed into boundary objects aligned along the main simulation domain's boundary.
! The algorithm below, which places alien particles into proper SEND* lists is the same as in ADVANCE_IONS.
! However, collision with a boundary object at this stage is not expected and is considered an error.
!
SUBROUTINE FIND_ALIENS_IN_ION_ADD_LIST

  USE ParallelOperationValues
  USE CurrentProblemValues
  USE ClusterAndItsBoundaries
  USE IonParticles, ONLY : N_ions_to_add, ion_to_add, N_spec

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER ierr

  INTEGER s, k, n 

  IF (N_of_inner_objects.EQ.0) RETURN

! cycle over ion species
  DO s = 1, N_spec
! cycle over particles of the ion species
     k=0
     DO WHILE (k.LT.N_ions_to_add(s))
        k = k + 1

        IF (symmetry_plane_X_left) THEN
           IF (ion_to_add(s)%part(k)%X.LT.c_X_area_min) THEN
! a particle crossed symmetry plane... this is an ion!!! YIKES !!!!
! since we are here, a particle to be added is beyond the symmetry plane at X=0
! the reason may be injection of ions from surfaces of inner material objects, which is not implemented
! so, presently, this should not happen, and we report it as an error
! note that for electrons this is possible (due to SEE), and such particles are moved symmetrically relative to plane x=0
!              ion_to_add(s)%part(k)%X = MAX(c_X_area_min, c_X_area_min + c_X_area_min - ion_to_add(s)%part(k)%X)
!              ion(s)%part(k)%VX = -ion(s)%part(k)%VX
              PRINT '("Proc ",i4," Error-00 in FIND_ALIENS_IN_ION_ADD_LIST, particle ",i8," of species ",i2," is beyond symmetry plane ",5(2x,e12.5),2x,i2)', Rank_of_process, k, s, &
                   & ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag
              CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
           END IF
        END IF

! most probable situation when the particle is inside the area
        IF ( (ion_to_add(s)%part(k)%X.GE.c_X_area_min) .AND. &
           & (ion_to_add(s)%part(k)%X.LE.c_X_area_max) .AND. &
           & (ion_to_add(s)%part(k)%Y.GE.c_Y_area_min) .AND. &
           & (ion_to_add(s)%part(k)%Y.LE.c_Y_area_max) ) CYCLE

! since we are here, particle s,k is outside the domain of this cluster

        IF (ion_to_add(s)%part(k)%X.LT.c_X_area_min) THEN

           IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the ion by the period length
              ion_to_add(s)%part(k)%X = ion_to_add(s)%part(k)%X + L_period_X
              IF (ion_to_add(s)%part(k)%Y.LT.c_Y_area_min) THEN
! particle is below the bottom side of the area
                 IF (Rank_of_master_below.LT.0) THEN
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-1 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)                       
                 END IF
                 CALL REMOVE_ION_FROM_ADD_LIST(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              ELSE IF (ion_to_add(s)%part(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
                 IF (Rank_of_master_above.LT.0) THEN
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-2 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                 END IF
                 CALL REMOVE_ION_FROM_ADD_LIST(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              END IF
              CYCLE
           END IF

           IF ( (ion_to_add(s)%part(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
              & (ion_to_add(s)%part(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

              IF (Rank_of_master_left.GE.0) THEN
! left neighbor cluster exists
                 CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
              ELSE
! left neighbor cluster does not exist
!                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX,  ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)   ! left
! error
                 print '("error-3 in FIND_ALIENS_IN_ION_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF

           ELSE IF (ion_to_add(s)%part(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the left bottom cell of the area

              SELECT CASE (c_left_bottom_corner_type)
                 CASE (HAS_TWO_NEIGHBORS)
                    IF ((c_X_area_min-ion_to_add(s)%part(k)%X).LT.(c_Y_area_min-ion_to_add(s)%part(k)%Y)) THEN
                       CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                    END IF

                 CASE (FLAT_WALL_BELOW)
                    CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)

                 CASE (FLAT_WALL_LEFT)
                    IF (ion_to_add(s)%part(k)%Y.GE.c_Y_area_min) THEN                 
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                       print '("error-4 in FIND_ALIENS_IN_ION_ADD_LIST")'
                       CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                    ELSE
                       CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)                       
                    END IF

                 CASE (SURROUNDED_BY_WALL)
!                    IF ((c_X_area_min-ion_to_add(s)%part(k)%X).LT.(c_Y_area_min-ion_to_add(s)%part(k)%Y)) THEN
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    ELSE
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    END IF
! error
                    print '("error-5 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

                 CASE (EMPTY_CORNER_WALL_LEFT)
                    CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
               
                 CASE (EMPTY_CORNER_WALL_BELOW)
                    CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)

              END SELECT

           ELSE IF (ion_to_add(s)%part(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top left cell of the area

              SELECT CASE (c_left_top_corner_type)
                 CASE (HAS_TWO_NEIGHBORS)
                    IF ((c_X_area_min-ion_to_add(s)%part(k)%X).LT.(ion_to_add(s)%part(k)%Y-c_Y_area_max)) THEN
                       CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                    END IF

                 CASE (FLAT_WALL_ABOVE)
                    CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)

                 CASE (FLAT_WALL_LEFT)
                    IF (ion_to_add(s)%part(k)%Y.LE.c_Y_area_max) THEN                 
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                       print '("error-6 in FIND_ALIENS_IN_ION_ADD_LIST")'
                       CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                    ELSE
                       CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                    END IF

                 CASE (SURROUNDED_BY_WALL)
!                    IF ((c_X_area_min-ion_to_add(s)%part(k)%X).LT.(ion_to_add(s)%part(k)%Y-c_Y_area_max)) THEN
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    ELSE
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    END IF
! error
                    print '("error-7 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

                 CASE (EMPTY_CORNER_WALL_LEFT)
                    CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)

                 CASE (EMPTY_CORNER_WALL_ABOVE)
                    CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)

              END SELECT
           END IF
           CALL REMOVE_ION_FROM_ADD_LIST(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           CYCLE
        END IF   !### IF (ion_to_add(s)%part(k)%X.LT.c_X_area_min) THEN

        IF (ion_to_add(s)%part(k)%X.GT.c_X_area_max) THEN

           IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
! if the cluster is periodically [X] connected on itself we avoid exchange in the X-direction
! shift the ion by the period length
              ion_to_add(s)%part(k)%X = ion_to_add(s)%part(k)%X - L_period_X
              IF (ion_to_add(s)%part(k)%Y.LT.c_Y_area_min) THEN
! particle is somewhere near the left bottom cell of the area
                 IF (Rank_of_master_below.LT.0) THEN
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-8 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)                       
                 END IF
                 CALL REMOVE_ION_FROM_ADD_LIST(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              ELSE IF (ion_to_add(s)%part(k)%Y.GT.c_Y_area_max) THEN
! particle is somewhere near the top left cell of the area
                 IF (Rank_of_master_above.LT.0) THEN
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-9 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 ELSE
                    CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                 END IF
                 CALL REMOVE_ION_FROM_ADD_LIST(s, k)  ! this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              END IF
              CYCLE
           END IF
           
           IF ( (ion_to_add(s)%part(k)%Y.GE.(c_Y_area_min+1.0_8)).AND. &
              & (ion_to_add(s)%part(k)%Y.LE.(c_Y_area_max-1.0_8)) ) THEN
! most probable situation when the particle stays at least one cell away from the Y-boundaries of the area

              IF (Rank_of_master_right.GE.0) THEN
! right neighbor cluster exists
                 CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
              ELSE
! right neighbor cluster does not exist
!                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                 print '("error-10 in FIND_ALIENS_IN_ION_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF

           ELSE IF (ion_to_add(s)%part(k)%Y.LT.(c_Y_area_min+1.0_8)) THEN
! particle is somewhere near the right bottom cell of the area

              SELECT CASE (c_right_bottom_corner_type)
                 CASE (HAS_TWO_NEIGHBORS)
                    IF ((ion_to_add(s)%part(k)%X-c_X_area_max).LT.(c_Y_area_min-ion_to_add(s)%part(k)%Y)) THEN
                       CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                    END IF

                 CASE (FLAT_WALL_BELOW)
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                
                 CASE (FLAT_WALL_RIGHT)
                    IF (ion_to_add(s)%part(k)%Y.GE.c_Y_area_min) THEN
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                       print '("error-11 in FIND_ALIENS_IN_ION_ADD_LIST")'
                       CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

                    ELSE
                       CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)                       
                    END IF

                 CASE (SURROUNDED_BY_WALL)
!                    IF ((ion_to_add(s)%part(k)%X-c_X_area_max).LT.(c_Y_area_min-ion_to_add(s)%part(k)%Y)) THEN
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    ELSE
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    END IF
! error
                       print '("error-12 in FIND_ALIENS_IN_ION_ADD_LIST")'
                       CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

                 CASE (EMPTY_CORNER_WALL_RIGHT)
                    CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)

                 CASE (EMPTY_CORNER_WALL_BELOW)
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)

              END SELECT

           ELSE IF (ion_to_add(s)%part(k)%Y.GT.(c_Y_area_max-1.0_8)) THEN
! particle is somewhere near the top right cell of the area

              SELECT CASE (c_right_top_corner_type)
                 CASE (HAS_TWO_NEIGHBORS)
                    IF ((ion_to_add(s)%part(k)%X-c_X_area_max).LT.(ion_to_add(s)%part(k)%Y-c_Y_area_max)) THEN
                       CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                    ELSE
                       CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                    END IF

                 CASE (FLAT_WALL_ABOVE)
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)

                 CASE (FLAT_WALL_RIGHT)
                    IF (ion_to_add(s)%part(k)%Y.LE.c_Y_area_max) THEN
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                       print '("error-13 in FIND_ALIENS_IN_ION_ADD_LIST")'
                       CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                    ELSE
                       CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)                    
                    END IF

                 CASE (SURROUNDED_BY_WALL)
!                    IF ((ion_to_add(s)%part(k)%X-c_X_area_max).LT.(ion_to_add(s)%part(k)%Y-c_Y_area_max)) THEN
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    ELSE
!                       CALL PROCESS_ION_COLL_WITH_BOUNDARY_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
!                    END IF
! error
                       print '("error-14 in FIND_ALIENS_IN_ION_ADD_LIST")'
                       CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

                 CASE (EMPTY_CORNER_WALL_RIGHT)
                    CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)

                 CASE (EMPTY_CORNER_WALL_ABOVE)
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)

              END SELECT

           END IF
           CALL REMOVE_ION_FROM_ADD_LIST(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           CYCLE
        END IF        !###  IF (ion_to_add(s)%part(k)%X.GT.c_X_area_max) THEN

! since we are here, c_X_area_min <= ion_to_add(s)&part(k)%Y <= c_X_area_max

        IF (ion_to_add(s)%part(k)%Y.GT.c_Y_area_max) THEN

           IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
              IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
                 CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
              ELSE
!                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                 print '("error-15 in FIND_ALIENS_IN_ION_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
              CALL REMOVE_ION_FROM_ADD_LIST(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              CYCLE
           END IF

           IF (Rank_of_master_above.GE.0) THEN
! neighbor cluster above exists
              CALL ADD_ION_TO_SEND_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
           ELSE
! neighbor cluster above does not exist
              IF ((ion_to_add(s)%part(k)%X.GE.(c_X_area_min+1.0_8)).AND.(ion_to_add(s)%part(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
!                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                 print '("error-16 in FIND_ALIENS_IN_ION_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              ELSE IF (ion_to_add(s)%part(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left top corner
                 IF (c_left_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
                    CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                 ELSE
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-17 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 END IF
              ELSE IF (ion_to_add(s)%part(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right top corner
                 IF (c_right_top_corner_type.EQ.EMPTY_CORNER_WALL_ABOVE) THEN
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX,  ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                 ELSE
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_ABOVE(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-18 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 END IF
              END IF
           END IF
           CALL REMOVE_ION_FROM_ADD_LIST(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           CYCLE
        END IF    !### IF (ion_to_add(s)%part(k)%Y.GT.c_Y_area_max) THEN

        IF (ion_to_add(s)%part(k)%Y.LT.c_Y_area_min) THEN

           IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN
              IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster above exists
                 CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
              ELSE
!                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                 print '("error-19 in FIND_ALIENS_IN_ION_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              END IF
              CALL REMOVE_ION(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
              CYCLE
           END IF

           IF (Rank_of_master_below.GE.0) THEN
! neighbor cluster below exists, remove particle and prepare to send it to the neighbor below
              CALL ADD_ION_TO_SEND_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
           ELSE
! neighbor cluster below does not exist
              IF ((ion_to_add(s)%part(k)%X.GE.(c_X_area_min+1.0_8)).AND.(ion_to_add(s)%part(k)%X.LE.(c_X_area_max-1.0_8))) THEN
! most probable situation when the particle stays at least one cell away from the X-boundaries of the area
!                 CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                 print '("error-20 in FIND_ALIENS_IN_ION_ADD_LIST")'
                 CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
              ELSE IF (ion_to_add(s)%part(k)%X.LT.(c_X_area_min+1.0_8)) THEN
! particle near the left bottom corner
                 IF (c_left_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
                    CALL ADD_ION_TO_SEND_LEFT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                 ELSE
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-21 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                 END IF
              ELSE IF (ion_to_add(s)%part(k)%X.GT.(c_X_area_max-1.0_8)) THEN
! particle near the right bottom corner
                 IF (c_right_bottom_corner_type.EQ.EMPTY_CORNER_WALL_BELOW) THEN
                    CALL ADD_ION_TO_SEND_RIGHT(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%AX, ion_to_add(s)%part(k)%AY, ion_to_add(s)%part(k)%tag)
                 ELSE
!                    CALL PROCESS_ION_COLL_WITH_BOUNDARY_BELOW(s, ion_to_add(s)%part(k)%X, ion_to_add(s)%part(k)%Y, ion_to_add(s)%part(k)%VX, ion_to_add(s)%part(k)%VY, ion_to_add(s)%part(k)%VZ, ion_to_add(s)%part(k)%tag)
! error
                    print '("error-22 in FIND_ALIENS_IN_ION_ADD_LIST")'
                    CALL MPI_ABORT(MPI_COMM_WORLD, ierr)
                END IF
              END IF
           END IF
           CALL REMOVE_ION_FROM_ADD_LIST(s, k)  !       this subroutine does  N_ions(s) = N_ions(s) - 1 and k = k-1
           CYCLE
        END IF   !### IF (ion_to_add(s)%part(k)%Y.LT.c_Y_area_min) THEN

     END DO   !### DO WHILE (k.LT.N_ions_to_add(s))
  END DO    !### DO s = 1, N_spec

END SUBROUTINE FIND_ALIENS_IN_ION_ADD_LIST

!----------------------------------------
!
SUBROUTINE FIND_INNER_OBJECT_COLL_IN_ION_ADD_LIST

!  USE ParallelOperationValues
  USE CurrentProblemValues
  USE IonParticles, ONLY : N_ions_to_add, ion_to_add, N_spec

  IMPLICIT NONE

  INTEGER s, k, n 
  REAL(8) X_move, Y_move

  IF (N_of_inner_objects.EQ.0) RETURN

  DO s = 1, N_spec
  
! find, process, and exclude ions which collided with inner objects
     k=0
     DO WHILE (k.LT.N_ions_to_add(s))
        k = k+1
        DO n = N_of_boundary_objects+1, N_of_boundary_and_inner_objects
           IF (ion_to_add(s)%part(k)%X.LE.whole_object(n)%Xmin) CYCLE
           IF (ion_to_add(s)%part(k)%X.GE.whole_object(n)%Xmax) CYCLE
           IF (ion_to_add(s)%part(k)%Y.LE.whole_object(n)%Ymin) CYCLE
           IF (ion_to_add(s)%part(k)%Y.GE.whole_object(n)%Ymax) CYCLE
! collision detected
           X_move = 2.0_8 * ion_to_add(s)%part(k)%VX * N_subcycles
           Y_move = 2.0_8 * ion_to_add(s)%part(k)%VY * N_subcycles
           CALL TRY_ION_COLL_WITH_INNER_OBJECT( s, &
                                              & ion_to_add(s)%part(k)%X, &
                                              & ion_to_add(s)%part(k)%Y, &
                                              & ion_to_add(s)%part(k)%VX, &
                                              & ion_to_add(s)%part(k)%VY, &
                                              & ion_to_add(s)%part(k)%VZ, &
                                              & X_move,                   &
                                              & Y_move,                   &
                                              & ion_to_add(s)%part(k)%tag )  !, &
!                                                  & whole_object(n) )
           CALL REMOVE_ION_FROM_ADD_LIST(s, k)  ! this subroutine does  N_ions_to_add(s) = N_ions_to_add(s) - 1 and k = k-1
           EXIT
        END DO
     END DO

  END DO

END SUBROUTINE FIND_INNER_OBJECT_COLL_IN_ION_ADD_LIST

!----------------------------------------
!
SUBROUTINE PROCESS_ADDED_IONS

  USE ParallelOperationValues
  USE CurrentProblemValues

  USE IonParticles, ONLY : N_ions_to_add, N_ions, max_N_ions, ion, ion_to_add, N_spec

  IMPLICIT NONE

  TYPE particle
     real(8) X
     real(8) Y
     real(8) VX
     real(8) VY
     real(8) VZ
     real(8) AX
     real(8) AY
     integer tag
  END TYPE particle

  TYPE(particle), ALLOCATABLE :: bufer(:)

  INTEGER ALLOC_ERR, DEALLOC_ERR
  INTEGER s, k, current_N

  DO s = 1, N_spec
  
     IF (N_ions_to_add(s).GT.(max_N_ions(s)-N_ions(s))) THEN
! increase the size of the main ion array
        current_N = max_N_ions(s)
        ALLOCATE(bufer(1:current_N), STAT=ALLOC_ERR)
        DO k = 1, current_N
           bufer(k)%X   = ion(s)%part(k)%X
           bufer(k)%Y   = ion(s)%part(k)%Y
           bufer(k)%VX  = ion(s)%part(k)%VX
           bufer(k)%VY  = ion(s)%part(k)%VY
           bufer(k)%VZ  = ion(s)%part(k)%VZ
           bufer(k)%AX  = ion(s)%part(k)%AX
           bufer(k)%AY  = ion(s)%part(k)%AY
           bufer(k)%tag = ion(s)%part(k)%tag
        END DO
        DEALLOCATE(ion(s)%part, STAT=DEALLOC_ERR)
        !NULLIFY(ion(s)%part)
        max_N_ions(s) = max_N_ions(s) + MAX(N_ions_to_add(s)-(max_N_ions(s)-N_ions(s)), max_N_ions(s)/10)
        ALLOCATE(ion(s)%part(1:max_N_ions(s)), STAT=ALLOC_ERR)
        DO k = 1, current_N
           ion(s)%part(k)%X   = bufer(k)%X
           ion(s)%part(k)%Y   = bufer(k)%Y
           ion(s)%part(k)%VX  = bufer(k)%VX
           ion(s)%part(k)%VY  = bufer(k)%VY
           ion(s)%part(k)%VZ  = bufer(k)%VZ
           ion(s)%part(k)%AX  = bufer(k)%AX
           ion(s)%part(k)%AY  = bufer(k)%AY
           ion(s)%part(k)%tag = bufer(k)%tag
        END DO
        DEALLOCATE(bufer, STAT=DEALLOC_ERR)
     END IF
  
     DO k = 1, N_ions_to_add(s)
        ion(s)%part(k+N_ions(s))%X   = ion_to_add(s)%part(k)%X
        ion(s)%part(k+N_ions(s))%Y   = ion_to_add(s)%part(k)%Y
        ion(s)%part(k+N_ions(s))%VX  = ion_to_add(s)%part(k)%VX
        ion(s)%part(k+N_ions(s))%VY  = ion_to_add(s)%part(k)%VY
        ion(s)%part(k+N_ions(s))%VZ  = ion_to_add(s)%part(k)%VZ
        ion(s)%part(k+N_ions(s))%AX  = ion_to_add(s)%part(k)%AX
        ion(s)%part(k+N_ions(s))%AY  = ion_to_add(s)%part(k)%AY
        ion(s)%part(k+N_ions(s))%tag = ion_to_add(s)%part(k)%tag
     END DO

! update electron counter
     N_ions(s) = N_ions(s) + N_ions_to_add(s)

  END DO

! clear counter of particles to be added to this process
  N_ions_to_add = 0  

END SUBROUTINE PROCESS_ADDED_IONS

