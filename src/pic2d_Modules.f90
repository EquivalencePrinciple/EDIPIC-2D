!------------------------------------
!
MODULE ExternalFields

  REAL(8) Bx_ext
  REAL(8) By_ext
  REAL(8) Bz_ext

  REAL(8) Ez_ext

! parameters magnetic field as in the paper of Boeuf and Garrigues
  REAL(8) y_Bmax
  REAL(8) Bz_0
  REAL(8) Bz_max
  REAL(8) Bz_Lsys
  REAL(8) half_over_sigma2_1
  REAL(8) half_over_sigma2_2
  REAL(8) a1, a2
  REAL(8) b1, b2

! parameters of set of wires with JZ current to get magnetic field with Bx, By
  INTEGER N_JZ_wires
  REAL(8), ALLOCATABLE :: JZwire_X(:)    ! x-coordinate of the wire
  REAL(8), ALLOCATABLE :: JZwire_Y(:)    ! y-coordinate of the wire
  REAL(8), ALLOCATABLE :: JZwire_JZ(:)    ! JZ current in the wire

END MODULE ExternalFields

!------------------------------------
!
MODULE ParallelOperationValues

  INTEGER, PARAMETER :: SHIFT1 = 0 !1000000
  INTEGER, PARAMETER :: SHIFT2 = 0 !2000000
  INTEGER, PARAMETER :: SHIFT3 = 0 !2000000

  INTEGER Rank_of_process
  INTEGER N_of_processes

  INTEGER Rank_of_process_left
  INTEGER Rank_of_process_right
  INTEGER Rank_of_process_above
  INTEGER Rank_of_process_below

  LOGICAL WHITE           ! for field solver, block-related

  LOGICAL WHITE_CLUSTER   ! for particle mover, cluster-related
  
  INTEGER field_master
  INTEGER particle_master
  INTEGER cluster_rank_key

  INTEGER Rank_of_master_left
  INTEGER Rank_of_master_right
  INTEGER Rank_of_master_above
  INTEGER Rank_of_master_below

  INTEGER COMM_BOUNDARY
  INTEGER Rank_boundary
  INTEGER N_processes_boundary

  INTEGER COMM_CLUSTER
  INTEGER Rank_cluster
  INTEGER N_processes_cluster

  INTEGER COMM_HORIZONTAL
  INTEGER Rank_horizontal
  INTEGER N_processes_horizontal

  INTEGER N_processes_cluster_left
  INTEGER N_processes_cluster_right
  INTEGER N_processes_cluster_above
  INTEGER N_processes_cluster_below

  INTEGER Rank_horizontal_left
  INTEGER Rank_horizontal_right
  INTEGER Rank_horizontal_above
  INTEGER Rank_horizontal_below

  INTEGER Rank_of_bottom_left_cluster_master

END MODULE  ParallelOperationValues

!------------------------------------
!
MODULE CurrentProblemValues

  INTEGER, PARAMETER :: VACUUM_GAP = 0
  INTEGER, PARAMETER :: METAL_WALL = 1
  INTEGER, PARAMETER :: PERIODIC_PIPELINE_X = 2
  INTEGER, PARAMETER :: PERIODIC_PIPELINE_Y = 3
  INTEGER, PARAMETER :: DIELECTRIC = 4
  INTEGER, PARAMETER :: SYMMETRY_PLANE = 5

  INTEGER, PARAMETER :: PERIODICITY_NONE    = 0
  INTEGER, PARAMETER :: PERIODICITY_X       = 1
  INTEGER, PARAMETER :: PERIODICITY_X_PETSC = 10
  INTEGER, PARAMETER :: PERIODICITY_X_Y     = 2

  REAL(8), PARAMETER :: e_Cl     = 1.602189d-19      ! Charge of single electron [Cl]
  REAL(8), PARAMETER :: m_e_kg   = 9.109534d-31      ! Mass of single electron [kg]
  REAL(8), PARAMETER :: true_eps_0_Fm = 8.854188d-12      ! The dielectric constant [F/m]
  REAL(8), PARAMETER :: mu_0_Hm  = 1.25663706d-6     ! vacuum permeability [H/m] or [m kg s^-2 A^-2]
  REAL(8), PARAMETER :: amu_kg   = 1.660565d-27      ! atomic mass unit [kg]
  REAL(8), PARAMETER :: kB_JK    = 1.38064852d-23    ! Boltzmann constant [J/K]

  REAL(8), PARAMETER :: pi = 3.141592653589793_8

  REAL(8) eps_0_Fm

  INTEGER N_filter

  INTEGER i_given_F_double_period_sys    ! in a system which is periodic in both X and Y directions, if there is no metal objects with given potential
  INTEGER j_given_F_double_period_sys    ! it is necessary to specify a point with some given potential, otherwise the field solver converges only if no dielectric objects are inside (pure plasma)
  REAL(8) given_F_double_period_sys      ! so here we specify the node (i,j) and the potential in the node

  INTEGER periodicity_flag   ! shows the presence of periodicity, is used to switch between periodic and non-periodic field solvers

  REAL(8) T_e_eV        ! scale electron temperature [eV]
  REAL(8) N_plasma_m3   ! scale electron density [m^-3]

!  INTEGER N_of_cells_debye   ! number of cells per scale electron Debye length
  REAL N_of_cells_debye   ! number of cells per scale electron Debye length

  INTEGER N_max_vel          ! maximal expected velocity [units of scale thermal electron velocity]
  
  INTEGER N_blocks_x    ! number of blocks (processes) along the X (horizontal) direction (numbering starts at 1)
  INTEGER N_blocks_y    ! number of blocks (processes) along the Y (vertical) direction (numbering starts at 1)

  INTEGER N_grid_block_x  ! number of cells along the X-direction in a block
  INTEGER N_grid_block_y  ! number of cells along the Y-direction in a block

  INTEGER N_of_particles_cell  ! number of macroparticles per cell for the scale density

  INTEGER cluster_N_blocks_x  ! number of blocks along the X-direction in a cluster
  INTEGER cluster_N_blocks_y  ! number of blocks along the Y-direction in a cluster
  INTEGER cluster_N_blocks    ! total number of processes (blocks) in a cluster

  INTEGER global_maximal_i    ! maximal values of indices in the X- and Y-directions
  INTEGER global_maximal_j    ! note that the origin has i/j=0/0

  INTEGER N_clusters_x        ! for rectangular domains, number of clusters along X
  INTEGER N_clusters_y        ! for rectangular domains, number of clusters along Y

  REAL(8) init_Te_eV   ! initial electron temperature [eV]
  REAL(8) init_Ne_m3   ! initial electron density [m^-3]

  REAL(8) v_Te_ms
  REAL(8) W_plasma_s1
  REAL(8) L_debye_m
  REAL(8) delta_x_m
  REAL(8) delta_t_s

  REAL(8) E_scale_Vm
  REAL(8) B_scale_T
  REAL(8) F_scale_V
  REAL(8) V_scale_ms

  REAL(8) N_scale_part_m3
  REAL(8) current_factor_Am2
  REAL(8) energy_factor_eV

  REAL(8) heat_flow_factor_Wm2
  REAL(8) temperature_factor_eV
   
  REAL(8) Chi_factor               !for "implicit permittivity"

  INTEGER T_cntr
  INTEGER Start_T_cntr
  INTEGER Max_T_cntr
  INTEGER N_subcycles              ! number of electron cycles per ion cycle (must be odd)

  LOGICAL MatVecsCreated           ! to control initial creation of the matrix and vectors

  TYPE segment_of_object
     INTEGER istart
     INTEGER jstart
     INTEGER iend
     INTEGER jend
!     REAL(8), ALLOCATABLE :: surface_charge(:)

     LOGICAL, ALLOCATABLE :: cell_is_covered(:)

! surface_n_flux_e(:)
! surface_n_flux_i(:,:)
! surface_energy_flux_e(:)
! surface_energy_flux_i(:)

  END TYPE segment_of_object

  TYPE boundary_object
     INTEGER object_id_number
     INTEGER object_type
     INTEGER electron_hit_count
     INTEGER electron_emit_count
     INTEGER, ALLOCATABLE :: ion_hit_count(:)
     REAL(8) total_charge
     REAL(8) phi                  ! electrostatic potential, used with metal walls and to calculate external field with dielectric walls
     REAL(8) phi_const            ! constant part of the electrostatic potential
     REAL(8) phi_var              ! time varying part of the electrostatic potential
     REAL(8) omega                ! frequency of time varying part
     REAL(8) phase                ! phase of time varying part

! waveform defines periodic non-harmonic variation of potential, the shape is defined by a user via data file
     LOGICAL use_waveform
     INTEGER N_wf_points                    ! number of waveform data points, must be no less than 2
     REAL,    ALLOCATABLE :: wf_phi(:)      ! array of potential values of waveform data points
     INTEGER, ALLOCATABLE :: wf_T_cntr(:)   ! array of times (in units of timesteps) of waveform data points

     REAL    N_electron_constant_emit      ! constant number of electron macroparticles to be injected each time step (for example due to emission from a thermocathode)
     INTEGER model_constant_emit
     REAL(8) factor_convert_vinj_normal_constant_emit    ! factor to be used to convert values provided by Get*Velocity procedures to desired temperature
     REAL(8) factor_convert_vinj_parallel_constant_emit  ! factor to be used to convert values provided by Get*Velocity procedures to desired temperature
     REAL(8) v_ebeam_constant_emit                       ! velocity of the electron beam

     REAL(8) eps_diel                 ! relative dielectric permittivity, used with dielectric walls only

     INTEGER number_of_segments
     INTEGER L                    ! total length of the boundary object
     INTEGER n_connected_to_start ! number of the boundary object connected to the start point of this boundary object
     INTEGER n_connected_to_end   ! number of the boundary object connected to the start point of this boundary object
     TYPE(segment_of_object), ALLOCATABLE :: segment(:)
     REAL(8), ALLOCATABLE :: phi_profile(:)   ! ######## must be redone ??? #######, just put it here to compile now, must be in the segment

     INTEGER ileft               ! for a rectangular inner object we specify only left, right, bottom, and top index limits
     INTEGER iright
     INTEGER jbottom
     INTEGER jtop
     INTEGER N_boundary_nodes    ! number of nodes 

     REAL(8) Xmin
     REAL(8) Xmax
     REAL(8) Ymin
     REAL(8) Ymax

! are -1 by default
! for inner objects placed across periodic boundaries
     INTEGER object_copy_periodic_X_right  ! number of a copy inner object placed across opposite X periodic boundary (on the right)
     INTEGER object_copy_periodic_Y_above  ! number of a copy inner object placed across opposite Y periodic boundary (above)

     LOGICAL object_does_NOT_cross_symmetry_plane_X  ! default is .TRUE.

     REAL(8), ALLOCATABLE :: surface_charge(:)
     REAL(8), ALLOCATABLE :: surface_charge_variation(:)

     CHARACTER(6) material

! variables below define electron-induced emission of secondary electrons
! based on the 1D EDIPIC
     
! we have three possible processes when an electron hits the wall, producing a secondary electron:
! 1 - elastic backscattering
! 2 - inelastic backscattering
! 3 - true secondary emission

     LOGICAL SEE_enabled              ! switches on/off secondary electron emission

     INTEGER Emitted_model(1:3)          ! for each possible process determines the way of processing
     
     INTEGER Elast_refl_type             ! type of reflection for elastically reflected electron: specular / random

! ELASTIC, MODEL 1
     REAL(8) setD_see_elastic            ! constant coefficient (ratio of emitted and incident electron fluxess)
     REAL(8) minE_see_elastic            ! LOWER energy boundary, [dimensionless] 
     REAL(8) maxE_see_elastic            ! UPPER energy boundary, [dimensionless]

! ELASTIC, MODEL 2
     REAL(8) E_elast_0                   ! threshold energy, [dimensionless]
     REAL(8) E_elast_max                 ! maximum emission energy, [dimensionless]
     REAL(8) maxD_elast                  ! maximum emission yield
     REAL(8) dE_elast                    ! rate of decaying (like half-width) at energies > E_elast_max_eV, [dimensionless]
     REAL(8) Frac_elast_highenergy       ! fraction of total SEE yield at high energies, >=0, <<1
     
! INELASTIC, MODEL 1
     REAL(8) setD_see_inelastic          ! constant coefficient (ratio of emitted and incident electron fluxes)
     REAL(8) minE_see_inelastic          ! LOWER energy boundary, [dimensionless] 
     REAL(8) maxE_see_inelastic          ! UPPER energy boundary, [dimensionless]

! INELASTIC, MODEL 2
     REAL(8) Frac_inelastic              ! fraction of total SEE yield, >=0, <<1

! TRUE SECONDARY, MODEL 1
     REAL(8) setD_see_true               ! the coefficient (ratio of emitted to incident electrons)
     REAL(8) minE_see_true               ! LOWER energy boundary, [dimensionless] 
     REAL(8) maxE_see_true               ! UPPER energy boundary, [dimensionless]

! TRUE SECONDARY, MODEL 2 / CLASSIC SEE COEFFICIENT
     REAL(8) E_see_0                     ! threshold energy, [dimensionless]
     REAL(8) E_see_max                   ! maximal emission energy, [dimensionless]
     REAL(8) maxD_see_classic            ! maximal emission coefficient (normal to the surface)
     REAL(8) k_smooth                    ! smoothness factor (from 0 for very rough to 2 for polished surfaces)

     REAL(8) T_see_true_eV               ! Temperature of injected true secondary electrons, [eV]
     REAL(8) factor_convert_seetrue_vinj ! factor to be used to convert values provided by Get*Velocity procedures to desired temperature
     REAL(8) lowest_energy_for_see       ! threshold energy below which no emission occurs

! variables below define emission of electrons caused by ion impact on the material surface
! based on the 1D EDIPIC
! here ii_ee stands for ion-induced-electron-emission

     LOGICAL reflects_all_ions        ! switches on/off reflection of all ions (presently specular)
     LOGICAL ion_induced_EE_enabled   ! switches on/off electron emission caused by ions hitting the wall

     REAL(8), ALLOCATABLE :: setD_ii_ee_true(:)                    ! the coefficient (ratio of emitted to incident electrons)
     REAL(8), ALLOCATABLE :: minE_ii_ee_true(:)                    ! LOWER energy boundary, [dimensionless] 
     REAL(8), ALLOCATABLE :: maxE_ii_ee_true(:)                    ! UPPER energy boundary, [dimensionless]
     REAL(8), ALLOCATABLE :: T_ii_ee_true_eV(:)                    ! Temperature of injected true secondary electrons, [eV]
     REAL(8), ALLOCATABLE :: factor_convert_ii_ee_true_vinj(:)     ! factor to be used to convert values provided by Get*Velocity procedures to desired temperature

  END TYPE boundary_object

  INTEGER N_of_boundary_objects
  INTEGER N_of_inner_objects
! inner objects can only be of type METAL_WALL or DIELECTRIC
! inner objects may overlap
! inner objects may touch the boundary of the domain
  INTEGER N_of_boundary_and_inner_objects  ! = N_of_boundary_objects + N_of_inner_objects

  TYPE(boundary_object), ALLOCATABLE :: whole_object(:)

  TYPE collided_particle
     INTEGER token
     REAL coll_coord
     REAL VX
     REAL VY
     REAL VZ
  END TYPE collided_particle

  TYPE boundary_object_statistics
     LOGICAL must_be_saved
     INTEGER max_N_of_saved_parts
     INTEGER N_of_saved_parts
     TYPE(collided_particle), ALLOCATABLE :: part(:)
  END TYPE boundary_object_statistics

  TYPE(boundary_object_statistics), ALLOCATABLE :: ion_colls_with_bo(:)
  TYPE(boundary_object_statistics), ALLOCATABLE :: e_colls_with_bo(:)

  REAL(8), ALLOCATABLE :: EX(:,:)      ! these arrays cover the whole cluster
  REAL(8), ALLOCATABLE :: EY(:,:)      !

!  REAL(8), ALLOCATABLE :: phi(:,:)     ! these arrays cover a single field-calculating block
  REAL(8), ALLOCATABLE :: rho(:,:)     !
  REAL(8), ALLOCATABLE :: X11(:,:)     !
  REAL(8), ALLOCATABLE :: X12(:,:)     !
  REAL(8), ALLOCATABLE :: X21(:,:)     !
  REAL(8), ALLOCATABLE :: X22(:,:)     !

END MODULE CurrentProblemValues

!-----------------------------   ????????????????????
!
MODULE ArraysOfGridValues

END MODULE ArraysOfGridValues

!-----------------------------
!
MODULE ElectronParticles

  TYPE electron_particle
     REAL(8) X
     REAL(8) Y
     REAL(8) VX
     REAL(8) VY
     REAL(8) VZ
     REAL(8) AX
     REAL(8) AY
     INTEGER tag
  END TYPE electron_particle

  INTEGER     N_electrons ! number of electron macroparticles
  INTEGER max_N_electrons ! current size (number of particles) of array of electrons macroparticles

  INTEGER     N_e_to_add  ! number of electron macroparticles to be added
  INTEGER max_N_e_to_add  ! current size (number of particles) of array electron_to_add

  TYPE(electron_particle), ALLOCATABLE :: electron(:)
  TYPE(electron_particle), ALLOCATABLE :: electron_to_add(:)

! electron counters
  INTEGER N_e_to_send_left
  INTEGER N_e_to_send_right
  INTEGER N_e_to_send_above
  INTEGER N_e_to_send_below

! sizes of arrays of sent particles
  INTEGER max_N_e_to_send_left
  INTEGER max_N_e_to_send_right
  INTEGER max_N_e_to_send_above
  INTEGER max_N_e_to_send_below
  
! electron arrays of particles to be sent to neighbors
  TYPE(electron_particle), ALLOCATABLE :: electron_to_send_left(:)
  TYPE(electron_particle), ALLOCATABLE :: electron_to_send_right(:)
  TYPE(electron_particle), ALLOCATABLE :: electron_to_send_above(:)
  TYPE(electron_particle), ALLOCATABLE :: electron_to_send_below(:)

END MODULE ElectronParticles

!-----------------------------
!
MODULE IonParticles

  LOGICAL ions_sense_magnetic_field
  LOGICAL ions_sense_EZ

!  INTEGER, PARAMETER :: N_spec = 1         ! number of ion species
  INTEGER N_spec        ! number of ion species

  TYPE particle
     REAL(8) X
     REAL(8) Y
     REAL(8) VX
     REAL(8) VY
     REAL(8) VZ
     REAL(8) AX
     REAL(8) AY
     INTEGER tag
  END TYPE particle

  TYPE particle_array
!     TYPE(particle), POINTER :: part(:)
     TYPE(particle), ALLOCATABLE :: part(:)
  END TYPE particle_array

  INTEGER, ALLOCATABLE :: Qs(:) !(1:N_spec)               ! q/e, [relative] charge of each ion species (e positive)
                                                          ! note, in the present version the ion charge can be either positive or negative 
                                                          ! but the absolute value should not exceed 3e
  REAL(8), ALLOCATABLE :: M_i_amu(:) !(1:N_spec)          ! ion mass, atomic mass unit

  REAL(8), ALLOCATABLE :: init_Ti_eV(:)                   ! initial ion temperature [eV]
  REAL(8), ALLOCATABLE :: init_NiNe(:)                    ! initial relative density

  ! NEW ADDITION : arrays to hold initial directed ion velocities [in units of ion thermal speed for specific species]
  REAL(8), ALLOCATABLE :: init_Uxi(:)  
  REAL(8), ALLOCATABLE :: init_Uyi(:)
  REAL(8), ALLOCATABLE :: init_Uzi(:)
  !--------------------------    

  REAL(8), ALLOCATABLE :: Ms(:) !(1:N_spec)               ! Mi/me
  REAL(8), ALLOCATABLE :: QMs(:) !(1:N_spec)              ! (q/e)(me/Mi)
  REAL(8), ALLOCATABLE :: QM2s(:) !(1:N_spec)             ! (q/e)(me/Mi)/2
!  REAL(8), ALLOCATABLE :: QM2sNsub(:) !(1:N_spec)         ! (q/e)(me/Mi)(N_subcycles/2)

  INTEGER, ALLOCATABLE ::      N_ions(:) !(1:N_spec)       ! number of ion macroparticles for each ion species
  INTEGER, ALLOCATABLE ::  max_N_ions(:) !(1:N_spec)       ! current size (number of particles) of array of macroparticles for each ion species

  INTEGER, ALLOCATABLE ::      N_ions_to_add(:) !(1:N_spec)  ! number of ion macroparticles to be added for each ion species
  INTEGER, ALLOCATABLE ::  max_N_ions_to_add(:) !(1:N_spec)  ! current size (number of particles) of array ion_to_add

  TYPE(particle_array), ALLOCATABLE :: ion(:)          ! Array of pointers on arrays of ions of different species
  TYPE(particle_array), ALLOCATABLE :: ion_to_add(:)   ! Array of pointers on arrays of ions of different species to be added

! ion species counters for particles sent to neighbors
  INTEGER, ALLOCATABLE ::  N_ions_to_send_left(:) !(1:N_spec)
  INTEGER, ALLOCATABLE ::  N_ions_to_send_right(:) !(1:N_spec)
  INTEGER, ALLOCATABLE ::  N_ions_to_send_above(:) !(1:N_spec)
  INTEGER, ALLOCATABLE ::  N_ions_to_send_below(:) !(1:N_spec)

! sizes of arrays of sent particles
  INTEGER, ALLOCATABLE ::  max_N_ions_to_send_left(:) !(1:N_spec)
  INTEGER, ALLOCATABLE ::  max_N_ions_to_send_right(:) !(1:N_spec)
  INTEGER, ALLOCATABLE ::  max_N_ions_to_send_above(:) !(1:N_spec)
  INTEGER, ALLOCATABLE ::  max_N_ions_to_send_below(:) !(1:N_spec)
  
! ion species arrays of particles to be sent to neighbors
  TYPE(particle_array), ALLOCATABLE :: ion_to_send_left(:)
  TYPE(particle_array), ALLOCATABLE :: ion_to_send_right(:)
  TYPE(particle_array), ALLOCATABLE :: ion_to_send_above(:)
  TYPE(particle_array), ALLOCATABLE :: ion_to_send_below(:)

END MODULE IonParticles

!----------------------------
!
MODULE ClusterAndItsBoundaries

  INTEGER, PARAMETER :: HAS_TWO_NEIGHBORS       = 10
  INTEGER, PARAMETER :: SURROUNDED_BY_WALL      = 20
  INTEGER, PARAMETER :: FLAT_WALL_LEFT          = 31
  INTEGER, PARAMETER :: FLAT_WALL_RIGHT         = 32
  INTEGER, PARAMETER :: FLAT_WALL_ABOVE         = 33
  INTEGER, PARAMETER :: FLAT_WALL_BELOW         = 34
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_LEFT  = 41
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_RIGHT = 42
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_BELOW = 43
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_ABOVE = 44

  INTEGER c_row
  INTEGER c_column

  REAL(8) c_X_area_min
  REAL(8) c_X_area_max
  REAL(8) c_Y_area_min
  REAL(8) c_Y_area_max

  REAL(8) c_X_area_min_ext
  REAL(8) c_X_area_max_ext
  REAL(8) c_Y_area_min_ext
  REAL(8) c_Y_area_max_ext

  INTEGER c_indx_x_min
  INTEGER c_indx_x_max
  INTEGER c_indx_y_min
  INTEGER c_indx_y_max

  INTEGER, PARAMETER :: ext_overlap = 1   ! natural overlap is 1 cell
                                          ! extended overlap is 1 + ext_overlap cells
                                          
  INTEGER c_indx_x_min_ext
  INTEGER c_indx_x_max_ext
  INTEGER c_indx_y_min_ext
  INTEGER c_indx_y_max_ext

  INTEGER c_left_bottom_corner_type
  INTEGER c_left_top_corner_type
  INTEGER c_right_bottom_corner_type
  INTEGER c_right_top_corner_type

  LOGICAL periodic_boundary_X_left
  LOGICAL periodic_boundary_X_right
  LOGICAL periodic_boundary_Y_below
  LOGICAL periodic_boundary_Y_above

  INTEGER i_period_x   ! integer values are easier to transmit to particle calculators
  INTEGER i_period_y   ! 
  REAL(8) L_period_X
  REAL(8) L_period_Y

  TYPE object_link
     INTEGER object_number
     INTEGER segment_number
     INTEGER istart
     INTEGER jstart
     INTEGER iend
     INTEGER jend
!     TYPE(segment_of_object) segment     
     REAL(8), ALLOCATABLE :: surface_charge(:) 

  END TYPE object_link

! if the object is bigger than the cluster area these flags show how to connect with neighbor clusters to get their contributions to the surface charge

  LOGICAL connect_left      ! if a cluster has connect_left=TRUE,  it has a neighbor on the left  with connect_right = TRUE
  LOGICAL connect_right     ! if a cluster has connect_right=TRUE, it has a neighbor on the right with connect_left  = TRUE
  LOGICAL connect_below     ! if a cluster has connect_below=TRUE, it has a neighbor below with connect_above = TRUE
  LOGICAL connect_above     ! if a cluster has connect_above=TRUE, it has a neighbor above with connect_below = TRUE

  LOGICAL symmetry_plane_X_left

  INTEGER n_left(1:2)       ! contain indices of objects in array c_local_object_part
  INTEGER n_right(1:2)      ! which endpoints will participate in surface charge exchange
  INTEGER n_below(1:2)      ! with left/right/below/above neighbor cluster
  INTEGER n_above(1:2)

  INTEGER, PARAMETER :: c_max_N_of_local_object_parts = 20
  INTEGER c_N_of_local_object_parts
  TYPE(object_link) c_local_object_part(1:c_max_N_of_local_object_parts)

  INTEGER c_N_of_local_object_parts_left
  INTEGER c_N_of_local_object_parts_right
  INTEGER c_N_of_local_object_parts_above
  INTEGER c_N_of_local_object_parts_below

  INTEGER c_index_of_local_object_part_left( 1:c_max_N_of_local_object_parts)
  INTEGER c_index_of_local_object_part_right(1:c_max_N_of_local_object_parts)
  INTEGER c_index_of_local_object_part_above(1:c_max_N_of_local_object_parts)
  INTEGER c_index_of_local_object_part_below(1:c_max_N_of_local_object_parts)

  TYPE field_calc_proc
     INTEGER rank

     INTEGER indx_x_min
     INTEGER indx_x_max
     INTEGER indx_y_min
     INTEGER indx_y_max

     INTEGER ibegin
     INTEGER iend
     INTEGER jbegin
     INTEGER jend
  END TYPE field_calc_proc

  TYPE(field_calc_proc), ALLOCATABLE :: field_calculator(:)

  REAL(8), ALLOCATABLE :: c_phi_ext(:,:)   ! potential in cluster domain + additional layer of nodes on each side

END MODULE ClusterAndItsBoundaries

!-----------------------------------
!
MODULE LoadBalancing

  INTEGER T_cntr_global_load_balance
  INTEGER T_cntr_cluster_load_balance
  INTEGER dT_global_load_balance            ! time step interval for global load balancing, must be an integer number of N_subcycles
  INTEGER dT_cluster_load_balance           ! time step interval for balancing load within a cluster, must be an integer number of dT_cluster_load_balance

  TYPE cluster_params
     INTEGER particle_master
     INTEGER N_processes
     INTEGER N_particles

     INTEGER N_processes_balanced
     REAL    avg_load_balanced
  END TYPE cluster_params

  TYPE(cluster_params), ALLOCATABLE :: cluster(:)

END MODULE LoadBalancing


!----------------------------
!
MODULE BlockAndItsBoundaries

!  INTEGER, PARAMETER :: VACUUM_GAP = 0
!  INTEGER, PARAMETER :: METAL_WALL = 1

  INTEGER, PARAMETER :: HAS_TWO_NEIGHBORS       = 10
  INTEGER, PARAMETER :: SURROUNDED_BY_WALL      = 20
  INTEGER, PARAMETER :: FLAT_WALL_LEFT          = 31
  INTEGER, PARAMETER :: FLAT_WALL_RIGHT         = 32
  INTEGER, PARAMETER :: FLAT_WALL_ABOVE         = 33
  INTEGER, PARAMETER :: FLAT_WALL_BELOW         = 34
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_LEFT  = 41
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_RIGHT = 42
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_BELOW = 43
  INTEGER, PARAMETER :: EMPTY_CORNER_WALL_ABOVE = 44

  INTEGER block_row        ! from 1 to N_blocks_y, number of row in the matrix of processes
  INTEGER block_column     ! from 1 to N_blocks_x, number of column in the matrix of processes

  INTEGER, ALLOCATABLE :: offset_of_block(:) !used to place corner nodes of 9-point stencil
  INTEGER, ALLOCATABLE :: nodenum(:,:)

  REAL(8) X_area_min
  REAL(8) X_area_max
  REAL(8) Y_area_min
  REAL(8) Y_area_max

  INTEGER indx_x_min
  INTEGER indx_x_max
  INTEGER indx_y_min
  INTEGER indx_y_max

  INTEGER left_bottom_corner_type
  INTEGER left_top_corner_type
  INTEGER right_bottom_corner_type
  INTEGER right_top_corner_type

  TYPE object_link
     INTEGER object_number
     INTEGER istart
     INTEGER jstart
     INTEGER iend
     INTEGER jend
!     TYPE(segment_of_object) segment     
!     REAL(8), ALLOCATABLE :: surface_charge(:) 
  END TYPE object_link

  INTEGER, PARAMETER :: max_N_of_local_object_parts = 20
  INTEGER N_of_local_object_parts
  TYPE(object_link) local_object_part(1:max_N_of_local_object_parts)

  LOGICAL block_has_symmetry_plane_X_left

  INTEGER N_of_local_object_parts_left
  INTEGER N_of_local_object_parts_right
  INTEGER N_of_local_object_parts_above
  INTEGER N_of_local_object_parts_below

  INTEGER index_of_local_object_part_left(1:max_N_of_local_object_parts)
  INTEGER index_of_local_object_part_right(1:max_N_of_local_object_parts)
  INTEGER index_of_local_object_part_above(1:max_N_of_local_object_parts)
  INTEGER index_of_local_object_part_below(1:max_N_of_local_object_parts)

  LOGICAL block_periodic_boundary_X_left
  LOGICAL block_periodic_boundary_X_right
  LOGICAL block_periodic_boundary_Y_below
  LOGICAL block_periodic_boundary_Y_above

  INTEGER N_to_solve_total
  INTEGER block_N_of_nodes_to_solve
  INTEGER global_offset

  INTEGER process_left_bottom_right_inner_node
  INTEGER process_left_solved_nodes_row_length

  INTEGER process_right_bottom_left_inner_node
  INTEGER process_right_solved_nodes_row_length

  INTEGER process_above_left_bottom_inner_node
  INTEGER process_below_left_top_inner_node
  
END MODULE BlockAndItsBoundaries


!----------------------------
!
MODULE Diagnostics

! diagnostic control
  INTEGER, PARAMETER :: Max_N_of_probes = 100

  INTEGER Save_probes_data_step        ! WriteOut_step   ! Time interval (in steps) for writing into the file
  INTEGER Save_probes_data_T_cntr      ! WriteStart_step ! Time (in steps) when writing starts

  INTEGER Save_probes_data_T_cntr_rff  ! the value of Save_probes_data_T_cntr read from file (rff) init_probes.dat
                                       ! we need to save this value for consistent initialization of snapshots 
                                       ! when checkpoints are used to continue simulation

  INTEGER TextOut_skip    ! Periods of writing to be skipped between text outputs, 0 = write each, 1 = skip one, etc.

  INTEGER text_output_counter    ! this counter is used to skip extra diagnostics periods between two text outputs

  INTEGER N_of_saved_records     ! number of records saved in each time dependence data file, is used to trim protocol data files
                                 ! when checkpoint is used to initialize the system

  INTEGER N_of_probes               ! Total number of probes
  INTEGER N_of_probes_cluster       ! number of probes in a cluster within the ordinary overlapping region
  INTEGER N_of_probes_cluster_ext   ! number of probes in a cluster within the extended overlapping region

  INTEGER, ALLOCATABLE :: Probe_position(:,:)           ! 1:2,1:N_of_probes_all ; (1,n)=x_n, (2,n)=y_n
  INTEGER, ALLOCATABLE :: List_of_probes_cluster(:)     ! 1:N_of_probes_cluster
  INTEGER, ALLOCATABLE :: List_of_probes_cluster_ext(:) ! 1:N_of_probes_cluster_ext

  REAL, ALLOCATABLE :: probe_Ne_cluster(:)        ! these arrays keep diagnostics values obtained in different subroutines

  REAL, ALLOCATABLE :: probe_VXe_cluster(:)       ! till they are saved in the main probe diagnostics routine
  REAL, ALLOCATABLE :: probe_VYe_cluster(:)       !
  REAL, ALLOCATABLE :: probe_VZe_cluster(:)       !

  REAL, ALLOCATABLE :: probe_WXe_cluster(:)       !
  REAL, ALLOCATABLE :: probe_WYe_cluster(:)       !
  REAL, ALLOCATABLE :: probe_WZe_cluster(:)       !

  REAL, ALLOCATABLE :: probe_VXVYe_cluster(:)     !
  REAL, ALLOCATABLE :: probe_VXVZe_cluster(:)     !
  REAL, ALLOCATABLE :: probe_VYVZe_cluster(:)     !

  REAL, ALLOCATABLE :: probe_QXe_cluster(:)       !
  REAL, ALLOCATABLE :: probe_QYe_cluster(:)       !
  REAL, ALLOCATABLE :: probe_QZe_cluster(:)       !

  REAL, ALLOCATABLE :: probe_Ni_cluster(:,:)      !

  REAL, ALLOCATABLE :: probe_VXi_cluster(:,:)     !
  REAL, ALLOCATABLE :: probe_VYi_cluster(:,:)     !
  REAL, ALLOCATABLE :: probe_VZi_cluster(:,:)     !

  REAL, ALLOCATABLE :: probe_VXVYi_cluster(:,:)   !
  REAL, ALLOCATABLE :: probe_VXVZi_cluster(:,:)   !
  REAL, ALLOCATABLE :: probe_VYVZi_cluster(:,:)   !

  REAL, ALLOCATABLE :: probe_WXi_cluster(:,:)     !
  REAL, ALLOCATABLE :: probe_WYi_cluster(:,:)     !
  REAL, ALLOCATABLE :: probe_WZi_cluster(:,:)     !

  REAL, ALLOCATABLE :: probe_QXi_cluster(:,:)     !
  REAL, ALLOCATABLE :: probe_QYi_cluster(:,:)     !
  REAL, ALLOCATABLE :: probe_QZi_cluster(:,:)     !

  TYPE diagnostic_process
     INTEGER rank
     INTEGER N_probes
     INTEGER N_probes_ext
     INTEGER, ALLOCATABLE :: probe_number(:)
     INTEGER, ALLOCATABLE :: probe_number_ext(:)
  END TYPE diagnostic_process

  TYPE(diagnostic_process), ALLOCATABLE :: diag_cluster(:)  ! this array wil be used by process with rank 0 which writes to the data files

END MODULE Diagnostics

!--------------------------
!
MODULE Snapshots

  INTEGER N_of_all_snaps                        ! number of all snapshots  
  INTEGER, ALLOCATABLE ::     Tcntr_snapshot(:)     ! timesteps when the snapshot files are written
  INTEGER, ALLOCATABLE :: save_evdf_snapshot(:)     ! flags controlling how evdf is saved
  INTEGER, ALLOCATABLE ::   save_pp_snapshot(:)     ! flags controlling saving of phase planes

  LOGICAL, ALLOCATABLE :: save_ionization_rates_2d(:)   ! turns on/off accumulation and saving of ioniization rates

  LOGICAL, ALLOCATABLE :: save_ions_collided_with_bo(:)    ! turns on/off saving of ions collided with boundary objects
                                                           ! (this must be confirmed by individual request for each boundary object)

  LOGICAL, ALLOCATABLE :: save_e_collided_with_bo(:)    ! turns on/off saving of electrons collided with boundary objects
                                                        ! (this must be confirmed by individual request for each boundary object)

  INTEGER current_snap                          ! index of current snapshot (which must be created)

! below, prefix cs stands for c-luster s-napshot
 
!  REAL, ALLOCATABLE :: cs_phi(:,:)
!  REAL, ALLOCATABLE :: cs_EX(:,:)
!  REAL, ALLOCATABLE :: cs_EY(:,:)
!  REAL, ALLOCATABLE :: cs_JXsum(:,:)
!  REAL, ALLOCATABLE :: cs_JYsum(:,:)
!  REAL, ALLOCATABLE :: cs_JZsum(:,:)

! arrays for calculation of VDF moments for both electrons and ions
! these arrays are made global to simplify synchronization of overlapping nodes
! these arrays are larger than the ordinary cluster domain to account for possible particles with midpoints outside the cluster domain
! these arrays are transferred to output arrays and discarded immediately after calculation and synchronization
  REAL, ALLOCATABLE :: cs_N(:,:)
  REAL, ALLOCATABLE :: cs_VX(:,:)
  REAL, ALLOCATABLE :: cs_VY(:,:)
  REAL, ALLOCATABLE :: cs_VZ(:,:)
  REAL, ALLOCATABLE :: cs_WX(:,:)
  REAL, ALLOCATABLE :: cs_WY(:,:)
  REAL, ALLOCATABLE :: cs_WZ(:,:)
  REAL, ALLOCATABLE :: cs_VXVY(:,:)
  REAL, ALLOCATABLE :: cs_VXVZ(:,:)
  REAL, ALLOCATABLE :: cs_VYVZ(:,:)
  REAL, ALLOCATABLE :: cs_QX(:,:)
  REAL, ALLOCATABLE :: cs_QY(:,:)
  REAL, ALLOCATABLE :: cs_QZ(:,:)

! electron output arrays
  REAL, ALLOCATABLE :: cs_Ne(:,:)
  REAL, ALLOCATABLE :: cs_VXe(:,:)
  REAL, ALLOCATABLE :: cs_VYe(:,:)
  REAL, ALLOCATABLE :: cs_VZe(:,:)
  REAL, ALLOCATABLE :: cs_WXe(:,:)
  REAL, ALLOCATABLE :: cs_WYe(:,:)
  REAL, ALLOCATABLE :: cs_WZe(:,:)
  REAL, ALLOCATABLE :: cs_QXe(:,:)    ! arrays required to calculate heat flow vector
  REAL, ALLOCATABLE :: cs_QYe(:,:)
  REAL, ALLOCATABLE :: cs_QZe(:,:)

! ion output arrays
  REAL, ALLOCATABLE :: cs_Ni(:,:,:)
  REAL, ALLOCATABLE :: cs_VXi(:,:,:)
  REAL, ALLOCATABLE :: cs_VYi(:,:,:)
  REAL, ALLOCATABLE :: cs_VZi(:,:,:)
  REAL, ALLOCATABLE :: cs_WXi(:,:,:)
  REAL, ALLOCATABLE :: cs_WYi(:,:,:)
  REAL, ALLOCATABLE :: cs_WZi(:,:,:)
  REAL, ALLOCATABLE :: cs_QXi(:,:,:)    ! arrays required to calculate heat flow vector
  REAL, ALLOCATABLE :: cs_QYi(:,:,:)
  REAL, ALLOCATABLE :: cs_QZi(:,:,:)

! flags for turning output of various parameters on/off

  LOGICAL save_data(1:38)

! variables for calculation of 1D and 2D velocity distribution functions

  INTEGER, PARAMETER :: NOANYVDF = 0
  INTEGER, PARAMETER :: ONLY1D = 1
  INTEGER, PARAMETER :: ONLY2D = 2
  INTEGER, PARAMETER :: BOTH1DAND2D = 3

  INTEGER N_vdfbox_x       ! number of spatial boxes in the X-direction in a cluster
  INTEGER N_vdfbox_y       !           same as above in the Y-direction
  INTEGER N_vdfbox_all     ! total number of spatial boxes in a cluster

  INTEGER N_max_vel_e      ! maximal electron velocity in units of v_Te_ms
  INTEGER N_vbins_e        ! number of velocity bins per v_Te_ms for electrons

  INTEGER N_max_vel_i      ! maximal electron velocity in units of v_Te_ms * sqrt(me/Ms)
  INTEGER N_vbins_i        ! number of velocity bins per v_Te_ms * sqrt(me/Ms) for ions

  INTEGER indx_v_min_e     ! limits of velocity bin index for electrons
  INTEGER indx_v_max_e     !

  INTEGER indx_v_min_i     ! limits of velocity bin index for ions
  INTEGER indx_v_max_i     !

  INTEGER, ALLOCATABLE :: evxdf(:,:)
  INTEGER, ALLOCATABLE :: evydf(:,:)
  INTEGER, ALLOCATABLE :: evzdf(:,:)

  INTEGER, ALLOCATABLE :: evxvydf(:,:,:)

  INTEGER, ALLOCATABLE :: isvxdf(:,:,:)
  INTEGER, ALLOCATABLE :: isvydf(:,:,:)
  INTEGER, ALLOCATABLE :: isvzdf(:,:,:)

! variables for saving phase planes

  INTEGER, PARAMETER :: NOANYPP = 0
  INTEGER, PARAMETER :: ONLYelectronPP = 1
  INTEGER, PARAMETER :: ONLYionPP = 2
  INTEGER, PARAMETER :: BOTHelectronANDionPP = 3

  INTEGER N_pp_boxes

  TYPE index_limits
     INTEGER imin
     INTEGER jmin
     INTEGER imax
     INTEGER jmax
  END TYPE index_limits

  TYPE(index_limits), ALLOCATABLE :: pp_box(:)

  TYPE specific_coll_diag
     REAL, ALLOCATABLE :: counter_local(:,:)
  END TYPE specific_coll_diag

  TYPE coll_diag
!     INTEGER N_of_activated_colproc
     TYPE(specific_coll_diag), ALLOCATABLE :: activated_collision(:)
  END TYPE coll_diag

  TYPE(coll_diag), ALLOCATABLE :: diagnostics_neutral(:)

END MODULE Snapshots

!*******************************************************************************************
! This module is used by two subroutines, calculating the arbitrary velocity
! according to the isotropic maxwell distribution
MODULE MaxwellVelocity
  REAL(8), PARAMETER :: U_max = 3.0_8
  INTEGER, PARAMETER :: R_max     = 9000 !300
  INTEGER, PARAMETER :: R_max_inj = 4500 !150
  REAL(8) v(0:R_max)
  REAL(8) v_inj(0:R_max_inj)
END MODULE MaxwellVelocity

!--------------------------------------------------------------
!
!MODULE SetupValues

!  USE IonParticles, ONLY : N_spec

!  LOGICAL ht_use_e_emission_from_cathode_zerogradf
!  LOGICAL ht_use_e_emission_from_cathode
!  LOGICAL ht_emission_constant
!  LOGICAL ht_grid_requested
!  LOGICAL ht_soft_grid_requested
!  LOGICAL ht_injection_inside
!  LOGICAL ht_use_ionization_source

!  INTEGER N_macro_constant_injection
!  REAL(8) injection_y
!  INTEGER grid_j
!  REAL(8) F_grid

!  INTEGER total_cathode_N_e_to_inject   ! this variable is used to calculate total value of negatively charged particles
                                        ! that left the system without being balanced by an escaping ion
                                        ! this value becomes negative if the number of escaping ions exceeds 
                                        ! that of escaping electrons

!  INTEGER total_cathode_N_e_injected    ! actual number of electron macroparticles injected, used for diagnostics only

!  INTEGER j_ion_source_1   ! full source boundaries
!  INTEGER j_ion_source_2   !

!  INTEGER c_j_ion_source_1   ! source boundaries within a cluster 
!  INTEGER c_j_ion_source_2   !

!  INTEGER, ALLOCATABLE :: N_to_ionize_total(:)      !#########(1:N_spec)
!  INTEGER, ALLOCATABLE :: N_to_ionize_cluster(:)    !#########(1:N_spec)

!  INTEGER, PARAMETER :: c_R_max = 1000              ! maximal value of integral of ionization rate distribution within the range allocated to a cluster
!  REAL(8) yi(0:c_R_max)

!  REAL(8) factor_convert_vinj                       ! velocity conversion factors, injected electrons
!  REAL(8) factor_convert_vion_e                     ! ionization electrons
!  REAL(8), ALLOCATABLE ::  factor_convert_vion_i(:) !########(1:N_spec)  ! ionization ions

!END MODULE SetupValues

!--------------------------------------------------------------
!
MODULE Checkpoints

  LOGICAL use_mpiio_checkpoint

  INTEGER use_checkpoint             ! 0/1/2 = don't use / use to continue older run / use to start a new run
  INTEGER T_cntr_to_continue

  INTEGER dT_save_checkpoint
  INTEGER T_cntr_save_checkpoint

  INTEGER Save_probes_data_T_cntr_check
  INTEGER N_of_saved_records_check
  INTEGER current_snap_check

END MODULE Checkpoints

!-------------------------
!
MODULE MCCollisions

  LOGICAL en_collisions_turned_off
  LOGICAL no_ionization_collisions
  LOGICAL no_rcx_collisions

  INTEGER N_neutral_spec

  TYPE collision_type
     LOGICAL activated
     INTEGER type
     INTEGER ion_species_produced
     REAL(8) threshold_energy_eV
     INTEGER N_crsect_points
     REAL(8), ALLOCATABLE :: energy_eV(:)
     REAL(8), ALLOCATABLE :: crsect_m2(:)
  END TYPE collision_type

  TYPE neutral_type
     CHARACTER(6) name
     REAL(8) M_amu
     REAL(8) N_m3
     REAL(8) T_K
     REAL(8) Ux         ! neutral directed velocity [units of neutral thermal speed]
     REAL(8) Uy
     REAL(8) Uz
     INTEGER N_en_colproc
     INTEGER N_of_energy_segments
     TYPE(collision_type), ALLOCATABLE :: en_colproc(:)
     REAL(8), ALLOCATABLE :: energy_segment_boundary_value(:)
     REAL(8), ALLOCATABLE :: energy_segment_step(:)

     logical rcx_on
     integer rcx_ion_species_index
     real(8) sigma_rcx_m2_1eV
     real(8) alpha_rcx
  END TYPE neutral_type

  TYPE(neutral_type), ALLOCATABLE :: neutral(:)

  TYPE brief_collision_type
     INTEGER id_number
     INTEGER type
     INTEGER ion_species_produced
     REAL(8) threshold_energy_eV
     REAL(8) ion_velocity_factor
  END TYPE brief_collision_type

  TYPE selected_collision_probability
     INTEGER N_of_energy_segments
     INTEGER N_of_energy_values
     INTEGER N_of_activated_colproc
     REAL(8) max_colliding_fraction
     TYPE(brief_collision_type), ALLOCATABLE :: colproc_info(:)
     INTEGER, ALLOCATABLE :: counter(:)
     REAL(8), ALLOCATABLE :: energy_segment_boundary_value(:)
     REAL(8), ALLOCATABLE :: energy_segment_step(:)
     INTEGER, ALLOCATABLE :: energy_segment_boundary_index(:)
     REAL(8), ALLOCATABLE :: prob_colproc_energy(:,:)
     REAL(8), ALLOCATABLE :: energy_eV(:)
  END TYPE selected_collision_probability

  TYPE(selected_collision_probability), ALLOCATABLE :: collision_e_neutral(:)

  TYPE rcx_collision_data
     logical rcx_on
     integer neutral_species_index
     real(8) probab_thermal
!     real(8) factor_eV
     real(8) vfactor !to scale the Maxwellian distribution
     INTEGER counter  ! to be used for diagnostics
  END TYPE rcx_collision_data

  TYPE(rcx_collision_data), allocatable :: collision_rcx(:)

! LINKED LIST, which store the numbers of particles participating in collisions

  TYPE binary_tree
     INTEGER number

     INTEGER neutral
     INTEGER indx_coll

     TYPE (binary_tree), POINTER :: Larger
     TYPE (binary_tree), POINTER :: Smaller     
  END TYPE binary_tree

  TYPE(binary_tree), POINTER :: Collided_particle

END MODULE MCCollisions
