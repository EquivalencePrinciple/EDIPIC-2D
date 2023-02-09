!===========================================
PROGRAM MainProg

  USE CurrentProblemValues
  USE ParallelOperationValues
  USE LoadBalancing
  USE ClusterAndItsBoundaries
  USE Checkpoints
  USE PETSc_Solver

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  logical final_step_flag
  INTEGER ierr
  REAL(8) start, finish
  INTEGER n

  CHARACTER(54) rmandmkdir_command    ! rm -rfv checkdir_TTTTTTTT ; mkdir -v checkdir_TTTTTTTT
                                      ! ----x----I----x----I----x----I----x----I----x----I----

  REAL(8) t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20

  INTERFACE
     FUNCTION convert_int_to_txt_string(int_number, length_of_string)
       CHARACTER*(length_of_string) convert_int_to_txt_string
       INTEGER int_number
       INTEGER length_of_string
     END FUNCTION convert_int_to_txt_string
  END INTERFACE

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, Rank_of_process, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_of_processes, ierr)

  CALL SET_PHYSICAL_CONSTANTS

  CALL PrepareMaxwellDistribIntegral

  Start_T_cntr = 0
  T_cntr_global_load_balance = Start_T_cntr
  T_cntr_cluster_load_balance = Start_T_cntr

  CALL INITIATE_PARAMETERS
!print *, "did INITIATE_PARAMETERS"

  CALL INITIATE_ELECTRON_NEUTRAL_COLLISIONS

  CALL INITIATE_ION_NEUTRAL_COLLISIONS

!print *, "did INITIATE_ELECTRON_NEUTRAL_COLLISIONS"
!CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 
!CALL MPI_ABORT(MPI_COMM_WORLD, ierr)

  CALL INITIATE_PROBE_DIAGNOSTICS

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL INITIATE_WALL_DIAGNOSTICS

  CALL INITIATE_GENERAL_DIAGNOSTICS

  CALL INITIATE_en_COLL_DIAGNOSTICS

  CALL INITIATE_in_COLL_DIAGNOSTICS

  CALL INITIATE_SNAPSHOTS

  CALL ADJUST_T_CNTR_SAVE_CHECKPOINT

  IF ((.NOT.use_mpiio_checkpoint).AND.(Rank_of_process.EQ.0).AND.(T_cntr_save_checkpoint.GE.Start_T_cntr)) THEN
! create a directory for a future checkpoint to make sure that the directory exists before the checkpoint files are saved
     rmandmkdir_command = 'rm -rfv checkdir_TTTTTTTT ; mkdir -v checkdir_TTTTTTTT'
                         ! ----x----I----x----I----x----I----x----I----x----I----
     rmandmkdir_command(18:25) = convert_int_to_txt_string(T_cntr_save_checkpoint, 8)
     rmandmkdir_command(47:54) = convert_int_to_txt_string(T_cntr_save_checkpoint, 8)
     CALL SYSTEM(rmandmkdir_command)
  END IF

  start = MPI_WTIME()

  DO T_cntr = Start_T_cntr, Max_T_cntr

     if (Rank_of_process.eq.0) print *, T_cntr

     IF (T_cntr.EQ.T_cntr_save_checkpoint) THEN
        IF (use_mpiio_checkpoint) THEN
           CALL SAVE_CHECKPOINT_MPIIO_2
        ELSE
           CALL SAVE_CHECKPOINT_POSIX
        END IF
        T_cntr_save_checkpoint = T_cntr_save_checkpoint + dT_save_checkpoint
        CALL ADJUST_T_CNTR_SAVE_CHECKPOINT
        IF ((.NOT.use_mpiio_checkpoint).AND.(Rank_of_process.EQ.0)) THEN
 !create a directory for a future checkpoint to make sure that the directory exists before the checkpoint files are saved
           rmandmkdir_command = 'rm -rfv checkdir_TTTTTTTT ; mkdir -v checkdir_TTTTTTTT'
                               ! ----x----I----x----I----x----I----x----I----x----I----
           rmandmkdir_command(18:25) = convert_int_to_txt_string(T_cntr_save_checkpoint, 8)
           rmandmkdir_command(47:54) = convert_int_to_txt_string(T_cntr_save_checkpoint, 8)
           CALL SYSTEM(rmandmkdir_command)
        END IF
     END IF
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     call report_total_number_of_particles

     IF (T_cntr.EQ.T_cntr_global_load_balance) THEN
        CALL GLOBAL_LOAD_BALANCE  ! includes calls to SET_COMMUNICATIONS 
                                  !                   DISTRIBUTE_CLUSTER_PARAMETERS
        T_cntr_global_load_balance = T_cntr_global_load_balance + dT_global_load_balance
     END IF
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     IF (T_cntr.EQ.T_cntr_cluster_load_balance) THEN
        CALL BALANCE_LOAD_WITHIN_CLUSTER
        T_cntr_cluster_load_balance = T_cntr_cluster_load_balance + dT_cluster_load_balance
     END IF
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!--------------- PRE-PUSH -----------------------------------------

     CALL ADVANCE_ELECTRONS_PRELIMINARY                         !   velocity:  n-1/2 --->"streaming", coordinate :: n ---> "streaming"
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL ADVANCE_IONS_PRELIMINARY                              !   velocity:  n-1/2 --->"streaming", coordinate :: n ---> "streaming"
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL FIND_ALIENS_IN_ION_ADD_LIST
     CALL FIND_ALIENS_IN_ELECTRON_ADD_LIST

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN  ! send/receive BOTH electrons and ions crossing the borders
        CALL EXCHANGE_PARTICLES_WITH_ABOVE_BELOW_NEIGHBOURS            ! need only one X-pass for self-connected X-periodic clusters
     ELSE                                                              !
        CALL EXCHANGE_PARTICLES_WITH_LEFT_RIGHT_NEIGHBOURS             ! in general, three passes X-Y-X are needed
        CALL EXCHANGE_PARTICLES_WITH_ABOVE_BELOW_NEIGHBOURS            !
        CALL EXCHANGE_PARTICLES_WITH_LEFT_RIGHT_NEIGHBOURS             !
     END IF
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL FIND_INNER_OBJECT_COLL_IN_ELECTRON_ADD_LIST
     CALL FIND_INNER_OBJECT_COLL_IN_ION_ADD_LIST

     CALL PROCESS_ADDED_ELECTRONS   
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL PROCESS_ADDED_IONS
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL GATHER_SURFACE_CHARGE_DENSITY      ! does c_local_object_part(n)%surface_charge = 0 .0 for non-master processes
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL GATHER_SURFACE_CHARGE_DENSITY_INNER_OBJECTS   ! whole_object%surface_charge_variation=0 is done here
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!------------- FINAL PUSH ------------------------------------

     CALL GATHER_STREAMING_CHARGE_DENSITY ! here surface charge density on inner dielectric objects is subtracted from the electron volume charge density
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL UPDATE_WALL_POTENTIALS(T_cntr) !if waveform is specified
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

if (Rank_of_process.eq.0) print '("before SolverInitialization ...")'

     CALL SolverInitialization(MPI_COMM_WORLD, MatVecsCreated) !fills the matrix
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL SOLVE_POTENTIAL_WITH_PETSC       
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     DO n = 1, N_filter
        CALL FILTER_POTENTIAL
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
     END DO

     CALL CALCULATE_ELECTRIC_FIELD
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!call save_phi
!call save_EX
!call save_EY
!call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
!if (T_cntr.eq.5) stop

     CALL ADVANCE_ELECTRONS_FINAL_PLUS                    !   velocity: "streaming" ---> n+1/2, coordinate :: "streaming" ---> n+1
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL SAVE_ELECTRONS_COLLIDED_WITH_BOUNDARY_OBJECTS   ! diagnostics
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL ADVANCE_IONS_FINAL_PLUS                         !   velocity: "streaming" ---> n+1/2, coordinate :: "streaming" ---> n+1
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL SAVE_IONS_COLLIDED_WITH_BOUNDARY_OBJECTS        ! diagnostics
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL SAVE_PROBES_DATA
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL CREATE_SNAPSHOT
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL FIND_ALIENS_IN_ION_ADD_LIST
     CALL FIND_ALIENS_IN_ELECTRON_ADD_LIST

     IF (periodic_boundary_X_left.AND.periodic_boundary_X_right) THEN  ! send/receive BOTH electrons and ions crossing the borders
        CALL EXCHANGE_PARTICLES_WITH_ABOVE_BELOW_NEIGHBOURS            ! need only one X-pass for self-connected X-periodic clusters
     ELSE                                                              !
        CALL EXCHANGE_PARTICLES_WITH_LEFT_RIGHT_NEIGHBOURS             ! in general, three passes X-Y-X are needed
        CALL EXCHANGE_PARTICLES_WITH_ABOVE_BELOW_NEIGHBOURS            !
        CALL EXCHANGE_PARTICLES_WITH_LEFT_RIGHT_NEIGHBOURS             !
     END IF
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL FIND_INNER_OBJECT_COLL_IN_ELECTRON_ADD_LIST
     CALL FIND_INNER_OBJECT_COLL_IN_ION_ADD_LIST

     CALL PROCESS_ADDED_ELECTRONS   
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

!if (Rank_of_process.eq.0) print '("A##############")'

     CALL COLLECT_PARTICLE_BOUNDARY_HITS
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!if (Rank_of_process.eq.0) print '("B##############")'

     CALL PERFORM_ELECTRON_NEUTRAL_COLLISIONS
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL SAVE_en_COLLISIONS
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL SAVE_en_COLLISIONS_2D               ! if snapshot is created at T_cntr, corresponding ionization rate is saved at T_cntr-1
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL PERFORM_RESONANT_CHARGE_EXCHANGE    ! the only ion-neutral collision kind for now, does not produce new ions
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL SAVE_in_COLLISIONS
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL PROCESS_ADDED_IONS                  ! add the new ions to the main array
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL UPDATE_ELECTRON_ACCELERATIONS        !### new, do it here because after final push some particles may change clusters
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!if (Rank_of_process.eq.0) print '("C##############")'

     CALL UPDATE_ION_ACCELERATIONS             !### new, do it here because after final push some particles may change clusters
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL PERFORM_ELECTRON_EMISSION_SETUP
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL PERFORM_ELECTRON_EMISSION_SETUP_INNER_OBJECTS
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL PROCESS_ADDED_ELECTRONS             ! add the new electrons to the main array
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL SAVE_BOUNDARY_PARTICLE_HITS_EMISSIONS           ! diagnostics
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL GATHER_SURFACE_CHARGE_DENSITY      ! does c_local_object_part(n)%surface_charge = 0 .0 for non-master processes
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     CALL GATHER_SURFACE_CHARGE_DENSITY_INNER_OBJECTS   ! whole_object%surface_charge_variation=0 is done here
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) 

!if (Rank_of_process.eq.0) print '("D##############")'

  END DO ! main loop

  finish = MPI_WTIME()

  PRINT '(2x,"**** Process ",i3" : Simulation time is  : ", f12.3," sec")', Rank_of_process, finish - start

  CALL FINISH_SNAPSHOTS

  CALL SolverDestroy
  CALL PetscFinalize(ierr)
  CALL MPI_FINALIZE(ierr)

END PROGRAM MainProg

