PROGRAM Redepo
  ! Redeposition Code
!******************************************************************************
! Calculates the pits and grooves erosion profile as a function of time
! including redeposition of sputtered material
!
! This is a fortran version of the code originally developed in Matlab
! by Paddy Simha and John Kelland
!******************************************************************************
 use iomesh_mod
 use sputter_mod
 use ifport
 USE, INTRINSIC :: IEEE_ARITHMETIC ! For IS_NAN, IS_INF, and Z_NAN
 USE, INTRINSIC :: ieee_features
 implicit none

 !-----------------------------------------------------------------------------
 !  Control Parameters
 !-----------------------------------------------------------------------------
 real*8, parameter :: moly_mass = 1.593E-22  ! Mass of a molybdenum (or grid material) atom in grams
 integer, parameter :: nionsample = 9E2  ! Unitless, amount of CEX ions used for analysis
 ! VERY IMPORTANT ^^^^^^ PSB CHANGED WAS 2.5E5 -----------------------------------------------------
 integer, parameter :: degree_of_precision_element_integration=2  ! degree of precision for triangle quadrature for redepo emission site
 real*8, parameter :: erosionStepSize = 25  ! Time step size in hrs
 real*8, parameter :: maxSimulationTime = 8200  ! Time to run the erosion iterations in hours
 real*8, parameter :: current_correction_factor=0.125  ! navg
 real*8 :: particle_current_factor = 1  ! leave this at 1, do not edit
 integer, parameter :: maxReflections_redep = 10  ! only 6 needed at max for deep pits and shallow grooves
 integer, parameter :: maxReflections_ions = 1000
 integer, parameter :: ion_data_source = 1  !  0 = CEX3D raw output file, 1 = processed ion input file, 2 = uniform flux (test case)
 integer, parameter :: mesh_data_source = 0  !  0 = initial flat face mesh, 1 = restart file, 2 = v-groove test case
 type (ion_type), dimension(nionsample) :: ionsample
 character(LEN = 50) :: ionfilename
 integer, parameter :: nerosionSteps = ceiling(maxSimulationTime/erosionStepSize) + 1
 real*8, dimension(maxnpt) :: delta_x_old, delta_y_old, delta_z_old
 INTEGER, PARAMETER :: DP = KIND(1.0d0) !PSB added


 integer :: erosionStep
 real*8 :: erosionTime
 real*8, dimension(nerosionSteps) :: max_pit_depth
 real*8, dimension(nerosionSteps) :: time_stamp
 integer :: time_count
 real*8, dimension(nerosionSteps) :: max_swing_depth
 real*8, dimension(nerosionSteps) :: groove_centre_depth
 real*8, dimension(nerosionSteps) :: pit_groove_ratio
 real*8, dimension(nerosionSteps) :: pitRate
 real*8, dimension(nerosionSteps) :: grooveCentreRate
 real*8, dimension(maxnpt,nerosionSteps) :: x_data_cell
 real*8, dimension(maxnpt,nerosionSteps) :: y_data_cell
 real*8, dimension(maxnpt,nerosionSteps) :: z_data_cell
 integer :: i, j, k, m, z , inode, icell, iion, irefl
 integer, dimension(maxnpt) :: node_elements
 integer, dimension(6) :: element_list
 integer, dimension(3, maxnpt) :: node_element_listing
 real*8, dimension(3, maxntri) :: xcoor
 real*8, dimension(maxntri) :: t, u, v
 logical, dimension(maxntri) :: intersect
 integer :: tri, nelem
 integer :: num_hits
 integer, dimension(1) :: impact_tri
 integer :: impact_cell
 integer :: nlostboys    ! number of ions that did not impact based on triangle_ray_intersection
 integer :: nbottomboys    ! number of ions that impacted the bottom of the domain based on triangle_ray_intersection
 real*8 :: inc_angle
 real*8 :: sum
 real*8 :: yield, mass_loss
 real*8, dimension(3) :: reflectedray
 real*8, parameter :: e = 1.60217662E-19 ! Elementary charge in C
 real*8, parameter :: massCarbon = 1.994307e-23 ! Carbon grid atom mass in g
 real*8, parameter :: densityCarbon = 0.0022 ! Carbon grid density in g/mm3 (pyrolytic graphite)
 real*8, dimension(3) :: new_coords  !  updated coords for nodes after erosion step
 real*8 :: startTime, stopTime
 real*8, dimension(3, maxntri) :: vert1
 real*8, dimension(3, maxntri) :: vert2
 real*8, dimension(3, maxntri) :: vert3
 character(:), allocatable :: lebfilename ! PSB
 integer :: nray  !PSB
 real*8 :: accel_thickness !PSB
 CHARACTER (LEN=256) :: command_string !PSB
CHARACTER (LEN=7) :: step_str !PSB
  CHARACTER (LEN=50) :: plot_subfolder  !PSB
 CHARACTER (LEN=50) :: plot_name_str !PSB
 INTEGER,DIMENSION(nionsample) :: ion_reflection_counts !PSB -> added to track number of reflections per ion
 REAL*8 :: curr_x, curr_y, curr_z !PSB
 REAL*8 :: total_mass_lost_this_step ! PSB This variable will hold the sum for the current step
 CHARACTER (LEN=250) :: python_script_path !PSB 
 CHARACTER (LEN=250) :: output_folder!PSB 
 CHARACTER (LEN=250) :: file_identifier_str !PSB 
 INTEGER, PARAMETER :: ION_IMPACT_UNIT = 114 ! PSB Choose a new, unused file unit number 
 CHARACTER (LEN=256) :: ion_impact_filename ! PSB Name for the detailed ion impact file
 REAL(DP), DIMENSION(3) :: current_coords !PSB
 LOGICAL :: is_nan_val, is_out_of_range! PSB
 INTEGER :: write_idx, i_clear !PSB
 REAL(DP), PARAMETER :: MIN_COORD_VAL = -1.0_DP  !PSB
 REAL(DP), PARAMETER :: MAX_COORD_VAL = 10000.0_DP !PSB
 INTEGER :: max_til_idx
 ! PSB for updateGeometry
LOGICAL :: delta_old_exist = .FALSE.  ! PSB for updateGeometry To track if delta_x_old/y_old/z_old are initialized
REAL*8, DIMENSION(maxnpt) :: delta_x_new, delta_y_new, delta_z_new ! PSB for updateGeometry
REAL*8 :: meanArea ! PSB for updateGeometry For adaptive refinement
REAL*8 :: erosion_time      ! Current total erosion time
    REAL*8, DIMENSION(maxnpt) :: node_mass_loss ! Mass loss per node (proxy for recession)

    ! Local Variables
    INTEGER :: new_npt                      ! Temporary for new number of points after refinement
    
    REAL*8, DIMENSION(3) :: current_node_normal ! Normal vector for current node
    REAL*8 :: node_recession               ! Recession distance for current node
    REAL*8 :: pre_sputter_x, pre_sputter_y, pre_sputter_z ! Original node coordinates before this step's movement
    REAL*8 :: current_x, current_y, current_z ! Current node coordinates (for clarity in calculations)
    REAL*8 :: coord_diff_x, coord_diff_y, coord_diff_z ! Differences for distance calculation
    REAL*8 :: distance_sq, min_dist_sq      ! Squared distances for performance
    INTEGER :: min_dist_idx                ! Index of closest node
    INTEGER :: rebuildIter                 ! Loop for mesh rebuilds
    INTEGER :: row_count                   ! Loop for triangle rows
    INTEGER, DIMENSION(3) :: row_indices   ! Indices of vertices in a triangle
    REAL*8 :: x_val_1, x_val_2, x_val_3     ! X-coordinates of triangle vertices
    REAL*8 :: y_val_1, y_val_2, y_val_3     ! Y-coordinates of triangle vertices
    LOGICAL :: is_on_hypotenuse             ! Flag if triangle is on hypotenuse boundary
    INTEGER, ALLOCATABLE, DIMENSION(:) :: flagged_elements ! Indices of triangles to remove
    INTEGER :: value_count                  ! Counter for flagged_elements
    INTEGER :: new_count                    ! Counter for new_til
    INTEGER, DIMENSION(3, maxntri) :: new_til ! Temporary storage for new triangulation
    REAL*8, DIMENSION(maxnpt) :: temp_x, temp_y, temp_z ! For filtering NaNs
    INTEGER :: temp_count                  ! Counter for temp arrays
    REAL*8 :: x_orig_hole_bound, y_orig_hypot_bound ! Original boundary coordinates for checks
    REAL*8, DIMENSION(3) :: current_centroid_coords ! Centroid coordinates of a cell
    REAL*8, DIMENSION(3,3) :: pointVec      ! Stores 3D coordinates of triangle vertices
    REAL*8, DIMENSION(3,3) :: joinVectors   ! Vectors forming sides of a triangle
    REAL*8, DIMENSION(3) :: angles          ! Angles of a triangle
    REAL*8 :: maxAngle                      ! Maximum angle in a triangle
    INTEGER :: maxInd                       ! Index of vertex with max angle
    REAL*8, DIMENSION(3) :: newPoint        ! New point coordinates for refinement
    REAL*8, DIMENSION(maxnpt) :: distances  ! Distances to current node
    INTEGER, DIMENSION(maxnpt) :: distSortIndices ! Sorted indices by distance
    REAL*8 :: sortedDistances_2nd           ! Second smallest distance for coarsening
    INTEGER :: coarsenCounter              ! Counter for coarsening operations
    INTEGER :: maxCoarsenCounter           ! Max coarsening operations allowed
    REAL*16 :: dx_last_x                     ! For calculating a characteristic length
    REAL*8 :: length_x_grid                 ! A measure of the x-dimension of the grid
    REAL*8 :: current_node_x, current_node_y, current_node_z ! Current node coordinates from vcl3d
    INTEGER :: ierr                        ! Error status from dtris2
    INTEGER, DIMENSION(size(vcl, dim = 2)) :: ind_dtris2 ! Indices for dtris2
    INTEGER, DIMENSION(3, size(til, dim = 2)) :: tnbr ! Neighbors for dtris2
    INTEGER, DIMENSION(size(vcl, dim = 2)) :: stack_dtris2 ! Stack for dtris2
    REAL*8 :: time_step_per_erosion_step ! Conversion factor for node movement
    LOGICAL :: is_flagged
    INTEGER :: erosion_step
    REAL*8, DIMENSION(3) :: distanceVec
    REAL*8 :: maxDist
    INTEGER :: maxDistInd
    REAL*8 :: COARSEN_FACTOR 
INTEGER :: num_points_to_target_for_coarsening
INTEGER :: ierr_coarsening_result
REAL*8 :: min_dist_1_value 
REAL*8 :: min_dist_2_value , current_distance , COARSENING_DISTANCE_THRESHOLD_FACTOR
INTEGER :: min_dist_1_idx_found 
  REAL*16, PARAMETER :: TOL_PROF = 1.0E-8   ! Tolerance for floating point comparisons
    REAL*16, PARAMETER :: TAN_30_DEG = 0.577350269190 ! TAN(30 degrees)
    INTEGER, PARAMETER :: MAX_GRID_POINTS = 10000 ! Maximum possible points
REAL*16 :: current_dist_sq
REAL*16 , PARAMETER :: num_NAN = 0
INTEGER, ALLOCATABLE :: point_map(:)
INTEGER :: k_index
REAL(8), ALLOCATABLE :: vcl3d_temp(:,:), vcl_temp(:,:)
TYPE(node_type), ALLOCATABLE :: node_temp(:)
    integer :: new_ntri
      integer :: ind(size(vcl, dim = 2))
      integer, dimension(size(vcl, dim = 2)) :: stack
min_dist_1_value = HUGE(0.0D0) ! Stores the smallest *non-zero* distance 
min_dist_2_value = HUGE(0.0D0) ! Stores the second smallest *non-zero* distance 
min_dist_1_idx_found = 0 ! Stores the index of the point for min_dist_1_value
COARSENING_DISTANCE_THRESHOLD_FACTOR = 0.20 

!-----------------------------------------------------------------------------
!  Read in CEX ion data
!-----------------------------------------------------------------------------
select case (ion_data_source)
  case (0)     
    !  Read in CEX ion data from original CEX3D output file
  case (1)     
    !  Read in CEX ion data from processed ion input file
    ionfilename = 'CEX_ion_data'
    call read_ion_file(ionfilename, domain, simulation_time, npandgions)
  case (2)     
    !  Create a uniform flux of ions for testing
  case default
    print*, "Improper input for ion_data_source (must be 0, 1, or 2)"
end select

!-----------------------------------------------------------------------------
!  Setup initial grid
!-----------------------------------------------------------------------------
select case (mesh_data_source)
    case (0)    
      !  Generate mesh for initial flat accel grid face
      call mesh_flat_surface(domain, npt, vcl, vcl3d, maxntri, nfacept, accel_thickness)
      call triangulate_surface(domain, npt, vcl, vcl3d, ntri, til, nfacetri, cell)

      !---------------------------------------------------------------------------
      ! Save the vertex coordinate list and mesh triangle index list for plotting
      !---------------------------------------------------------------------------
      open(13, FILE = 'VCL') 
        do i = 1, npt
          write(13,*) vcl(1,i), vcl(2,i)
        end do
      close(13)
      open(12, FILE = 'Triangles')
        do i = 1, ntri
          write(12,*) til(1,i), til(2,i), til(3,i)
        end do
      close(12)

      call mesh_boundaries(domain, npt, vcl, vcl3d, ntri, til, cell, nfacept)
   
      ! Define vertices of each cell
      do icell = 1, ntri
        vert1(:, icell) = vcl3d(:, til(1, icell))  
        vert2(:, icell) = vcl3d(:, til(2, icell))
        vert3(:, icell) = vcl3d(:, til(3, icell))
      end do
    startTime = dclock()
    do icell = 1, ntri
      ! Define the triangle incidence list for cell
      cell(icell)%til(:) = til(:, icell)
      ! Define vertices of cell
      cell(icell)%vert1(:) = vcl3d(:, cell(icell)%til(1))
      cell(icell)%vert2(:) = vcl3d(:, cell(icell)%til(2))
      cell(icell)%vert3(:) = vcl3d(:, cell(icell)%til(3))
      ! Calculate other cell parameters
      call calc_cell_params(cell(icell), icell)
    end do
    stopTime = dclock()
    !print *, "cell: elapsed time:", stopTime - startTime
    startTime = dclock()
    do inode = 1, nfacept
      ! Define coordinates of node from vertex coordinate list
      node(inode)%coords(:) = vcl3d(:, inode)
      ! Calculate other node parameters
      call calc_node_params(node(inode), inode, domain, nfacetri, cell)
    end do
    stopTime = dclock()
    !print *, "node: elapsed time:", stopTime - startTime

case (1)    !  Read in mesh from a restart file
case (2)    !  Generate a mesh for a v-groove test case
case default
  !print*, "Improper input for mesh_data_source (must be 0, 1, or 2)"
end select

print *, " - - - - Post Boundaries - - - - "
print *, "Number of points (npt): ", npt
print *, "Maximum allowed points (maxnpt): ", maxnpt
print *, "Number of triangles (ntri): ", ntri
print *, "Number of face points (nfacept): ", nfacept
print *, "Maximum allowed triangles (maxntri): ", maxntri

 max_til_idx = 0
      DO i = 1, ntri
          max_til_idx = MAX(max_til_idx, til(1,i), til(2,i), til(3,i))
      END DO
      PRINT *, "Max index in TIL after boundaries: ", max_til_idx
      PRINT *, "Current npt (should match max_til_idx) after boundaires: ", npt
 

!-----------------------------------------------------------------------------
!  Start simulation of erosion profiles
!-----------------------------------------------------------------------------
!
!  Load array of Lebedev points used in sphere quadrature
!
lebfilename = '/LebedevRays.txt'
call read_lebedev_file(lebfilename, ray, nray)
!
!  Initialize Stuff
!
time_count = 1
max_pit_depth(time_count) = 0
max_swing_depth(time_count) = 0
groove_centre_depth(time_count) = 0
time_stamp(time_count) = 0
time_count = time_count + 1

delta_x_old = vcl(1, :)
delta_y_old = vcl(1, :)
delta_z_old = vcl(1, :)

!---------------------------------------------------------------------------
!  Update VCL3d with the vertex coordinate list and mesh triangle index list for plotting
!---------------------------------------------------------------------------
open(113, FILE = 'VCL3D')
  ! header write(113,*) "time, x, y, z"
  do inode = 1,npt
    write(113,'("         0",3(",",ES13.6))') vcl3d(1,inode), vcl3d(2,inode), vcl3d(3,inode)
  end do
close(113)

!---------------------------------------------------------------------------
! PSB: Export 3D mesh data for Python visualization (directly uses VCL3D and Triangles files)
!---------------------------------------------------------------------------
! No need to write separate files like mesh_points_3d.txt or mesh_triangles_3d.txt.
! The Python script now reads vertex coordinates from VCL3D and triangle indices from Triangles.
plot_subfolder = 'mesh_plots'
CALL EXECUTE_COMMAND_LINE("mkdir -p " // TRIM(plot_subfolder), wait=.true., exitstat=j)
IF (j /= 0) THEN
    PRINT *, "Warning: Could not create directory ", TRIM(plot_subfolder), ". Error code: ", j
    PRINT *, "Plot will be saved in the current directory instead."
    plot_subfolder = '.'
END IF
plot_name_str = "post_boundaries" ! Identifier for the plot filename ! Command to call the Python script, passing the output folder and plot identifier
command_string = "python plot_mesh.py " // TRIM(plot_subfolder) // " " // TRIM(plot_name_str)
CALL EXECUTE_COMMAND_LINE(TRIM(command_string), wait=.true.)

! ---------------------------------------------------------------------------
! EROSION ITERATION LOOP
! ---------------------------------------------------------------------------
!print *, "Starting erosion iterations..!"

EROSION: do erosionStep = 1, 1  ! PSB: Iterate beyond 1:5
    print *, "------------------------------------------------------------------------"   ! PSB
    print *, "                          Erosion Step:", erosionStep                     ! PSB
    print *, "Number of Points: " , npt
    erosionTime = erosionStep * erosionStepSize

    !-----------------------------------------------------------------------------
    !  Select random subset of ions out of npandgions total that hit the pits and grooves pattern
    !-----------------------------------------------------------------------------
    do iion = 1, npandgions
      ionindex(iion) = iion
    end do
    call randperm(ionindex,npandgions,nionsample)
    do iion = 1, nionsample
      ionsample(iion) = ion(ionindex(iion))
    end do

    !-----------------------------------------------------------------------------
    !  Sputtering section
    !-----------------------------------------------------------------------------
    ! PSB: Iterate through each triangular cell on the face of the accelerator grid and set mass_loss of each cell to 0.0
    do icell = 1, nfacetri
      cell(icell)%mass_loss = 0.0d0
    end do
    !PSB: Reset number of ions that were lost or hit the bottom face
    nlostboys = 0
    nbottomboys = 0
    ion_reflection_counts = 0 !PSB added 
    
    DOION: do iion = 1, nionsample !PSb: Iterate through each of the randomly selected ions
        CALL trace_ion_path(iion, 0, ion(iion)%origin, ion(iion)%dir, &
                          nlostboys, nbottomboys, ion_reflection_counts, &
                          vert1, vert2, vert3, t, u, v, xcoor, intersect, &
                          cell, node, til, maxReflections_ions)
    end do DOION
    close(112)
    print *, "Ion path traced"

    !print*, "No. ions that didn't hit anything: ", nlostboys
    !print*, "No. ions that hit the bottom: ", nbottomboys
    !print*, "No. of ions total (including ions that hit bottom or nothin): ", npt

    !-----------------------------------------------------------------------------
    !  Geometry update section PSB
    !-----------------------------------------------------------------------------
    do inode = 1, nfacept
      node(inode)%vel = node(inode)%mass_loss/simulation_time/(node(inode)%area * densityCarbon) * node(inode)%normal  ! vel is positive inward
      !If node is below the thruster grid thickness but is rising upwards, ignore it.
        ! IF (vcl3d(3,i) <= domain%zmin .AND. node_recession <= 0.0d0) THEN
        !     ! If the node is at or below the original Z-min (bottom of the grid) AND it's trying to move downwards , then ignore the recession for now.
        !     ! If node_recession is positive (sputtering), it moves along -normal, inward of material.
        !     ! This prevents the grid from "growing" into the thruster body from below.
        !     node_recession = 0.0d0 
        !     !node_mass_loss(i) = 0.0d0
        !     print *, "Recession at punch through: " , vcl3d(3,i) , "    " , (domain%zmax - accel_thickness)
        ! END IF
        call ABintegrate(node(inode)%coords, node(inode)%vel, node(inode)%vel_old, erosionStepSize * 3600, node(inode)%firstStep)
        vcl3d(:,inode) = node(inode)%coords
        vcl(:,inode) = node(inode)%coords(1:2)
        node(inode)%vel_old = node(inode)%vel
    end do
    nfacept = npt

    ! Filter out NaN values from node coordinates and update VCL 
            temp_count = 0
            DO i = 1, npt
                IF ((vcl3d(1,i) == vcl3d(1,i)) .AND. &
                    (vcl3d(2,i) == vcl3d(2,i)) .AND. &
                    (vcl3d(3,i) == vcl3d(3,i))) THEN
                    temp_count = temp_count + 1
                    temp_x(temp_count) = vcl3d(1,i)
                    temp_y(temp_count) = vcl3d(2,i)
                    temp_z(temp_count) = vcl3d(3,i)
                ELSE
                    PRINT *, "NaN value at ", i, ":", vcl3d(1,i), vcl3d(2,i), vcl3d(3,i)
                END IF
            END DO
            npt = temp_count
            DO i = 1, npt
                vcl3d(1,i) = temp_x(i)
                vcl3d(2,i) = temp_y(i)
                vcl3d(3,i) = temp_z(i)
                node(i)%coords(1) = temp_x(i)
                node(i)%coords(2) = temp_y(i)
                node(i)%coords(3) = temp_z(i)
            END DO
            DO i = 1, npt
                vcl(1,i) = vcl3d(1,i)
                vcl(2,i) = vcl3d(2,i)
            END DO

    ! PSB wrote this section to remove duplicate points or 0,0 not in matlab code
            CALL clean_vcl3d(vcl3d, npt, maxnpt)
            
            ! Debug file update
                OPEN(UNIT=21, FILE='DEBUGPREMESHONEA', ACTION='WRITE') 
                        WRITE(21, '(A)') '# ErosionTime, X, Y, Z'
                            DO inode = 1, npt 
                                curr_x = node(inode)%coords(1) 
                                curr_y = node(inode)%coords(2) 
                                curr_z = node(inode)%coords(3)
                                WRITE(21, '(F10.2, "," , ES25.18, ",", ES25.18, ",", ES25.18)')erosionTime, curr_x, curr_y, curr_z
                            END DO 
                CLOSE(21) 
                OPEN(UNIT=36, FILE='DEBUGPREMESHONEB', ACTION='WRITE') 
                    WRITE(36, '(A)') 'triangles'
                    do i = 1, ntri
                    write(36,*) til(1,i), til(2,i), til(3,i)
                    end do
                CLOSE(36) 

    ! Remesh once
            print *, " - - - - PRE REMESH 1 ----------"
            print *, "Number of points (npt): ", npt
            print *, "Number of triangles (ntri): ", ntri
            print *, "Number of face points (nfacept): ", nfacept

            call dtris2 (npt, maxnpt, vcl, ind, ntri, til, tnbr, stack, ierr)
            IF (ierr /= 0) THEN
                print *, "ERROR IN DTRIS2 - ERROR NUMBER " , ierr
            END IF
            !call triangulate_surface(domain, npt, vcl, vcl3d, ntri, til, nfacetri, cell) ! Line 1256 MATLAB
            nfacept = npt 

            print *, " - - - - POST REMESH 1 ----------"
            print *, "Number of points (npt): ", npt
            print *, "Number of triangles (ntri): ", ntri
            print *, "Number of face points (nfacept): ", nfacept

    ! FLAGGED ELEMENTS CLEAN UP MESHING 
            ! Clean up the meshing: Assess the triangulation for improper elements
            ALLOCATE(flagged_elements(maxntri)) ! Allocate dynamically or use a fixed large size
            flagged_elements = 0 ! Initialize to 0 or some indicator that it's not flagged
            value_count = 0
            ! Re-evaluate `nfacetri` based on the new `ntri` from `dtris2` +  will store valid triangles in `new_til` and then copy back to `til`.
            new_count = 0
            DO i = 1, ntri
                row_indices = til(:,i)
                ! Check if the indices are within bounds of `vcl3d`. If `npt` has decreased, some `til` entries might point to old, removed indices. This check is crucial.
                IF (row_indices(1) > npt .OR. row_indices(2) > npt .OR. row_indices(3) > npt) THEN
                    print *, "skip because npt = " , npt , " and we have " , row_indices(1), row_indices(2), row_indices(3)
                    CONTINUE ! Skip this triangle, its vertices are no longer valid 
                END IF
                x_val_1 = vcl3d(1, row_indices(1))
                x_val_2 = vcl3d(1, row_indices(2))
                x_val_3 = vcl3d(1, row_indices(3))
                y_val_1 = vcl3d(2, row_indices(1))
                y_val_2 = vcl3d(2, row_indices(2))
                y_val_3 = vcl3d(2, row_indices(3))
                ! Check if an element is composed entirely of nodes along the hypotenuse
                is_on_hypotenuse = (ABS(y_val_1 - (x_val_1 * TAN_30_DEG)) < TOL_PROF) .AND. (ABS(y_val_2 - (x_val_2 * TAN_30_DEG)) < TOL_PROF) .AND. (ABS(y_val_3 - (x_val_3 * TAN_30_DEG)) < TOL_PROF)
                IF (is_on_hypotenuse) THEN
                    value_count = value_count + 1
                    flagged_elements(value_count) = i ! Store the index of the flagged triangle
                END IF

                ! ! Check sum of angles is not more than 180
                !         ! Extract node coordinates for the current triangle
                !         pointVec(1,:) = vcl3d(:, row_indices(1))
                !         pointVec(2,:) = vcl3d(:, row_indices(2))
                !         pointVec(3,:) = vcl3d(:, row_indices(3))
                !         ! Centroid check: `sqrt(centroid_new(1)^2+centroid_new(2)^2)<0.7207*domain_x_max`
                !         current_centroid_coords(1) = (pointVec(1,1) + pointVec(2,1) + pointVec(3,1)) / 3.0d0
                !         current_centroid_coords(2) = (pointVec(1,2) + pointVec(2,2) + pointVec(3,2)) / 3.0d0
                !         current_centroid_coords(3) = (pointVec(1,3) + pointVec(2,3) + pointVec(3,3)) / 3.0d0
                !         ! Originally 0.7207
                !         IF (SQRT(current_centroid_coords(1)**2 + current_centroid_coords(2)**2) < (0.7207 * domain%xmax)) THEN ! MATCH WITH LINE 418 IN IOMESH - IMPORTANT
                !         ! The earea outside the 0.55 boudnary is where the more important sputtering phenomena occurs, otherwise, too close tothe grid hole to oberve any sputtering effects
                !         CONTINUE ! Skip if centroid is too close to the origin
                !         END IF
                !         ! Calculate normalized edge vectors
                !         joinVectors(1,:) = pointVec(2,:) - pointVec(1,:)
                !         joinVectors(2,:) = pointVec(3,:) - pointVec(2,:)
                !         joinVectors(3,:) = pointVec(1,:) - pointVec(3,:)

                !         joinVectors(1,:) = joinVectors(1,:) / norm2(joinVectors(1,:))
                !         joinVectors(2,:) = joinVectors(2,:) / norm2(joinVectors(2,:))
                !         joinVectors(3,:) = joinVectors(3,:) / norm2(joinVectors(3,:))

                !         ! Calculate internal angles using dot product and arccos
                !         angles(1) = ACOS(-DOT_PRODUCT(joinVectors(3,:), joinVectors(1,:)))
                !         angles(2) = ACOS(-DOT_PRODUCT(joinVectors(1,:), joinVectors(2,:)))
                !         angles(3) = ACOS(-DOT_PRODUCT(joinVectors(2,:), joinVectors(3,:)))
                        
                !         IF (abs((angles(1) + angles(2) + angles(3)) - PI) > 1e-8) THEN
                !             print *, 'NOTE: Triangle angles dont sum to pi for index ' , i , 'since we have ' , angles(1) , angles(2) , angles(3)
                !             value_count = value_count + 1
                !             flagged_elements(value_count) = i ! Store the index of the flagged triangle
                !         END IF
                            
            END DO ! Store all proper elements into new_til
            new_count = 0
            DO i = 1, ntri ! Check if the current triangle index 'i' is in the list of flagged elements
                    is_flagged = .FALSE.
                    DO j = 1, value_count
                        IF (flagged_elements(j) == i) THEN
                            is_flagged = .TRUE.
                            print *, "Flagged triangle found"
                            EXIT
                        END IF
                    END DO
                    IF (.NOT. is_flagged) THEN
                        new_count = new_count + 1
                        new_til(:, new_count) = til(:, i)
                    END IF
            END DO
            DEALLOCATE(flagged_elements) ! Replace the original triangulation with its modified counterpart
            ntri = new_count
            DO i = 1, ntri
                til(:, i) = new_til(:, i)
            END DO
            nfacetri = ntri ! Update nfacetri to reflect the new number of surface triangles

    ! UPDATE cell parameters with new data
            ! Now, update cell parameters and node parameters for the new mesh
            ! Recalculate cell parameters
            DO i = 1, ntri
                cell(i)%til(:) = til(:, i)
                cell(i)%vert1(:) = vcl3d(:, cell(i)%til(1))
                cell(i)%vert2(:) = vcl3d(:, cell(i)%til(2))
                cell(i)%vert3(:) = vcl3d(:, cell(i)%til(3))
                CALL calc_cell_params(cell(i), i)
            END DO
            ! Recalculate node parameters for surface nodes
            DO i = 1, nfacept ! Only surface nodes have these properties
                node(i)%coords(:) = vcl3d(:, i)
                CALL calc_node_params(node(i), i, domain, nfacetri, cell)
            END DO

    !-----------------------------------------------------------------------------
    !  Adaptive Refinement 
    !-----------------------------------------------------------------------------
    ! Calculate mean area for triangle elements 
            IF (nionsample > 0 ) THEN  ! always
                    ! Equivalent to MATLAB's findElements=find(element_centroids(:,1)>0.8); and meanArea=mean(element_areas(findElements));
                    sum = 0.0d0
                    temp_count = 0
                    DO i = 1, nfacetri ! Assuming nfacetri covers the relevant surface elements
                        IF (cell(i)%centroid(1) > 0.8d0) THEN
                            sum = sum + cell(i)%area
                            temp_count = temp_count + 1
                        END IF
                    END DO
                    IF (temp_count > 0) THEN
                        meanArea = sum / DBLE(temp_count)
                    ELSE
                        meanArea = 0.0d0 ! Avoid division by zero
                    END IF
            END IF

    ! Coarsen by obtuse angles, number of points does not change, a vertex of an obtuse triangle is moved to the midpoint of opposite edge
            IF (MOD(erosion_step, 5) == 1 .OR. erosion_time <= 250.0d0) THEN 
                    NFACEITER: DO i = 1, nfacetri ! Iterate over surface triangles
                        row_indices = cell(i)%til(:)
                        ! Boundary check: `max(element_pick(2:4))>length(thruster_grid_x)` ! Corresponds to `max(row_indices) > npt`
                        IF (MAXVAL(row_indices) > npt) THEN
                            print *, "ERROR: Invalid vertex detected, should not happen: " , MAXVAL(row_indices), " and npt = " , npt
                            CONTINUE ! Skip if any vertex index is out of bounds (might be from old invalid element)
                        END IF
                        ! Extract node coordinates for the current triangle
                        pointVec(1,:) = vcl3d(:, row_indices(1))
                        pointVec(2,:) = vcl3d(:, row_indices(2))
                        pointVec(3,:) = vcl3d(:, row_indices(3))
                        ! Centroid check: `sqrt(centroid_new(1)^2+centroid_new(2)^2)<0.7207*domain_x_max`
                        current_centroid_coords(1) = (pointVec(1,1) + pointVec(2,1) + pointVec(3,1)) / 3.0d0
                        current_centroid_coords(2) = (pointVec(1,2) + pointVec(2,2) + pointVec(3,2)) / 3.0d0
                        current_centroid_coords(3) = (pointVec(1,3) + pointVec(2,3) + pointVec(3,3)) / 3.0d0
                        ! Originally 0.7207
                        IF (SQRT(current_centroid_coords(1)**2 + current_centroid_coords(2)**2) < (0.7207 * domain%xmax)) THEN ! MATCH WITH LINE 418 IN IOMESH - IMPORTANT
                        ! The earea outside the 0.55 boudnary is where the more important sputtering phenomena occurs, otherwise, too close tothe grid hole to oberve any sputtering effects
                        CONTINUE ! Skip if centroid is too close to the origin
                        END IF
                        ! Calculate normalized edge vectors
                        joinVectors(1,:) = pointVec(2,:) - pointVec(1,:)
                        joinVectors(2,:) = pointVec(3,:) - pointVec(2,:)
                        joinVectors(3,:) = pointVec(1,:) - pointVec(3,:)

                        joinVectors(1,:) = joinVectors(1,:) / norm2(joinVectors(1,:))
                        joinVectors(2,:) = joinVectors(2,:) / norm2(joinVectors(2,:))
                        joinVectors(3,:) = joinVectors(3,:) / norm2(joinVectors(3,:))

                        ! Calculate internal angles using dot product and arccos
                        angles(1) = ACOS(-DOT_PRODUCT(joinVectors(3,:), joinVectors(1,:)))
                        angles(2) = ACOS(-DOT_PRODUCT(joinVectors(1,:), joinVectors(2,:)))
                        angles(3) = ACOS(-DOT_PRODUCT(joinVectors(2,:), joinVectors(3,:)))

                        ! Find the maximum angle and its corresponding index
                        maxAngle = angles(1)
                        maxInd = 1
                        IF (angles(2) > maxAngle) THEN
                            maxAngle = angles(2)
                            maxInd = 2
                        END IF
                        IF (angles(3) > maxAngle) THEN
                            maxAngle = angles(3)
                            maxInd = 3
                        END IF

                        ! If the maximum angle is obtuse (> 3*pi/4 radians = 135 degrees)
                        MAXANGLELOOP: IF (maxAngle > (3.0d0 * PI / 4.0d0)) THEN
                            ! Coarsen by moving the vertex opposite the obtuse angle to the midpoint of the opposite edge.
                            SELECT CASE (maxInd)
                                CASE (1) ! Vertex 1 has the obtuse angle, move it to midpoint of edge 2-3
                                    newPoint(1) = 0.5d0 * (pointVec(2,1) + pointVec(3,1))
                                    newPoint(2) = 0.5d0 * (pointVec(2,2) + pointVec(3,2))
                                    newPoint(3) = 0.5d0 * (pointVec(2,3) + pointVec(3,3))
                                CASE (2) ! Vertex 2 has the obtuse angle, move it to midpoint of edge 1-3
                                    newPoint(1) = 0.5d0 * (pointVec(1,1) + pointVec(3,1))
                                    newPoint(2) = 0.5d0 * (pointVec(1,2) + pointVec(3,2))
                                    newPoint(3) = 0.5d0 * (pointVec(1,3) + pointVec(3,3))
                                CASE (3) ! Vertex 3 has the obtuse angle, move it to midpoint of edge 1-2
                                    newPoint(1) = 0.5d0 * (pointVec(1,1) + pointVec(2,1))
                                    newPoint(2) = 0.5d0 * (pointVec(1,2) + pointVec(2,2))
                                    newPoint(3) = 0.5d0 * (pointVec(1,3) + pointVec(2,3))
                            END SELECT

                            ! Update the coordinates of the moved node
                            vcl3d(:, row_indices(maxInd)) = newPoint(:)
                            node(row_indices(maxInd))%coords(:) = newPoint(:)
                            print *, "NOTE: Obtuse angle triangle with new vertices" , vcl3d(:, row_indices(maxInd))
                        END IF MAXANGLELOOP
                    END DO NFACEITER
            END IF

    ! ADD REFINE BY AREAS

    ! Refine based on lengths: add more points
            LENGTHREFINE: IF (MOD(erosion_step, 5) == 1 .OR. erosion_time <= 250.0d0) THEN !
                ! Refine based on areas Line 1358 Matlab  This section adds new nodes (centroids of large elements).
                    NFACETRILOOP: DO i = 1, nfacetri
                        row_indices = cell(i)%til(:)
                        ! Boundary check
                        IF (MAXVAL(row_indices) > npt) THEN
                            CONTINUE
                        END IF
                        ! Extract node coordinates
                        pointVec(1,:) = vcl3d(:, row_indices(1))
                        pointVec(2,:) = vcl3d(:, row_indices(2))
                        pointVec(3,:) = vcl3d(:, row_indices(3))

                        IF (cell(i)%area > (1.25d0 * meanArea)) THEN
                            ! Centroid check
                            current_centroid_coords(1) = (pointVec(1,1) + pointVec(2,1) + pointVec(3,1)) / 3.0d0
                            current_centroid_coords(2) = (pointVec(1,2) + pointVec(2,2) + pointVec(3,2)) / 3.0d0
                            current_centroid_coords(3) = (pointVec(1,3) + pointVec(2,3) + pointVec(3,3)) / 3.0d0

                            IF (SQRT(current_centroid_coords(1)**2 + current_centroid_coords(2)**2) < (0.8182d0 * domain%xmax)) THEN
                                CONTINUE
                            END IF

                            ! Add new point if space is available
                            ! Do -10 to account for the fact that mesh boundaries adds points based on boundary locations (~ 6 pts)
                            IF (npt < (maxnpt * 0.6)) THEN
                                npt = npt + 1
                                vcl3d(1, npt) = current_centroid_coords(1)
                                vcl3d(2, npt) = current_centroid_coords(2)
                                vcl3d(3, npt) = current_centroid_coords(3)
                                node(npt)%coords(:) = current_centroid_coords(:)
                                ! Initialize other node properties if necessary (e.g., normal, on_bndry)
                                ! These will be properly calculated after remeshing.
                                print *, "NOTE: Length Refinement pt added at index:", npt , " with coords ",  current_centroid_coords(:)
                                node(npt)%on_bndry = 0
                                node(npt)%on_corner = 0
                                node(npt)%area = 0.0d0
                                node(npt)%mass_loss = 0.0d0
                                node(npt)%vel(:) = 0.0d0
                                node(npt)%vel_old(:) = 0.0d0
                                node(npt)%firstStep = 1
                                ! Clear in_cells and cell_angle
                                node(npt)%in_cells = 0
                                node(npt)%cell_angle = 0.0d0
                            ELSE
                                PRINT *, "Warning: Max number of points (maxnpt) reached during area refinement."
                                EXIT ! Exit the loop if no more points can be added
                            END IF
                        END IF
                    END DO NFACETRILOOP
            END IF LENGTHREFINE
    
    ! Filter out NaN values from node coordinates and update VCL 
            temp_count = 0
            DO i = 1, npt
                    IF ((vcl3d(1,i) == vcl3d(1,i)) .AND. &
                        (vcl3d(2,i) == vcl3d(2,i)) .AND. &
                        (vcl3d(3,i) == vcl3d(3,i))) THEN
                            IF ((vcl3d(1,i) /= num_NAN) .AND. (vcl3d(2,i) /= num_NAN) .AND. (vcl3d(3,i) /= num_NAN)) THEN
                                temp_count = temp_count + 1
                                temp_x(temp_count) = vcl3d(1,i)
                                temp_y(temp_count) = vcl3d(2,i)
                                temp_z(temp_count) = vcl3d(3,i)
                            END IF
                    ELSE
                        PRINT *, "NaN value at ", i, ":", vcl3d(1,i), vcl3d(2,i), vcl3d(3,i)
                    END IF
            END DO
            npt = temp_count
            DO i = 1, npt
                vcl3d(1,i) = temp_x(i)
                vcl3d(2,i) = temp_y(i)
                vcl3d(3,i) = temp_z(i)
                node(i)%coords(1) = temp_x(i)
                node(i)%coords(2) = temp_y(i)
                node(i)%coords(3) = temp_z(i)
            END DO
            DO i = 1, npt
                vcl(1,i) = vcl3d(1,i)
                vcl(2,i) = vcl3d(2,i)
            END DO

    OPEN(UNIT=33, FILE='DEBUGONEA', ACTION='WRITE') 
        WRITE(33, '(A)') '# ErosionTime, X, Y, Z'
            DO inode = 1, npt 
            curr_x = node(inode)%coords(1) 
            curr_y = node(inode)%coords(2) 
            curr_z = node(inode)%coords(3)
            WRITE(33, '(F10.2, "," , ES25.18, ",", ES25.18, ",", ES25.18)')erosionTime, curr_x, curr_y, curr_z
            END DO 
    CLOSE(33) 

    OPEN(UNIT=34, FILE='DEBUGONEB', ACTION='WRITE') 
        WRITE(34, '(A)') '# ErosionTime, X, Y, Z'
            WRITE(34, '(A)') 'triangles'
            do i = 1, ntri
              write(34,*) til(1,i), til(2,i), til(3,i)
            end do
    CLOSE(34) 

    !-----------------------------------------------------------------------------
    !  Coarsening Algorithm
    !-----------------------------------------------------------------------------
            coarsenCounter = 0
            maxCoarsenCounter = FLOOR(DBLE(npt) / 5.0d0)
            NPTLOOPS: DO i = 1, npt
                    IF ((vcl3d(1,i) /= vcl3d(1,i)) .OR. (.NOT. vcl3d(2,i) == vcl3d(2,i)) .OR. (.NOT. vcl3d(3,i) == vcl3d(3,i))) THEN
                        CONTINUE 
                    END IF
                    current_node_x = vcl3d(1,i)
                    current_node_y = vcl3d(2,i)
                    current_node_z = vcl3d(3,i)
                    IF (SQRT(current_node_x**2 + current_node_y**2) < (0.8182d0 * domain%xmax)) THEN
                        CONTINUE
                    END IF
                    IF (ABS(current_node_x - domain%xmax) < TOL_PROF .AND. &
                        ABS(current_node_y - (TAN_30_DEG * domain%xmax)) < TOL_PROF) THEN
                        CONTINUE
                    END IF
                    IF (ABS(current_node_x - domain%xmax) < TOL_PROF .AND. &
                        ABS(current_node_y - 0.0d0) < TOL_PROF) THEN
                        CONTINUE
                    END IF
                    min_dist_sq = HUGE(0.0d0) 
                    sortedDistances_2nd = HUGE(0.0d0)
                    min_dist_idx = 0
                    
                    DO j = 1, npt
                        IF (j == i) THEN
                            CYCLE
                        END IF
                        coord_diff_x = vcl3d(1,j) - current_node_x
                        coord_diff_y = vcl3d(2,j) - current_node_y
                        coord_diff_z = vcl3d(3,j) - current_node_z
                        distances(j) = SQRT(coord_diff_x**2 + coord_diff_y**2 + coord_diff_z**2)
                        ! print *, 'distance of point ' , i , ' to ' , j , ' is ' , distances(j) 
                    END DO 

                    ! SORTING
                    min_dist_sq = HUGE(0.0d0) 
                    sortedDistances_2nd = HUGE(0.0d0)
                    min_dist_idx = 0
                    DO j = 1, npt
                       IF (j == i) THEN
                           CYCLE ! Skip self-distance
                       END IF
                       coord_diff_x = vcl3d(1,j) - current_node_x
                       coord_diff_y = vcl3d(2,j) - current_node_y
                       coord_diff_z = vcl3d(3,j) - current_node_z
                       current_dist_sq = coord_diff_x**2 + coord_diff_y**2 + coord_diff_z**2

                       IF (current_dist_sq < min_dist_sq) THEN
                           sortedDistances_2nd = min_dist_sq ! The previous min becomes the second min
                           min_dist_sq = current_dist_sq
                           min_dist_idx = j
                       ELSE IF (current_dist_sq < sortedDistances_2nd) THEN
                           sortedDistances_2nd = current_dist_sq
                       END IF
                   END DO ! Now, SQRT(min_dist_sq) is the closest distance, and min_dist_idx is its index. ! SQRT(sortedDistances_2nd) is the second closest distance.
                   print *, "Closest pt to ", i , " is ", min_dist_idx, " with dist " , SQRT(min_dist_sq)
                   !print *, "2nd closest other distance is ", SQRT(sortedDistances_2nd)

                   COARSEN_FACTOR = 0.08
                   print *, "Value at which coarsening occurs is hardcoded and is " , COARSEN_FACTOR 
                   ! If the distance to the second closest point (closest *other* point) is below the threshold
                   print *, SQRT(sortedDistances_2nd)
                   IF (SQRT(sortedDistances_2nd) < COARSEN_FACTOR) THEN
                       ! Calculate the midpoint of the current point (i) and its closest neighbor (min_dist_idx)
                       newPoint(1) = 0.5d0 * (vcl3d(1,i) + vcl3d(1,min_dist_idx))
                       newPoint(2) = 0.5d0 * (vcl3d(2,i) + vcl3d(2,min_dist_idx))
                       newPoint(3) = 0.5d0 * (vcl3d(3,i) + vcl3d(3,min_dist_idx))
                       ! Mark the current point (i) as 0,0,0
                       vcl3d(1,i) = 0
                       vcl3d(2,i) = 0
                       vcl3d(3,i) = 0
                       ! Update node structure for point 'i' as well
                       node(i)%coords(1) = num_NAN
                       node(i)%coords(2) = num_NAN
                       node(i)%coords(3) = num_NAN
                       ! Boundary protection checks for the point *to be moved* (min_dist_idx)
                       ! If min_dist_idx is a protected boundary point, we skip moving it.
                        IF (ABS(vcl3d(1,min_dist_idx) - domain%xmax) < TOL_PROF .AND. &
                            ABS(vcl3d(2,min_dist_idx) - 0.0d0) < TOL_PROF) THEN
                            ! This point (min_dist_idx) is a protected "groove centre", do not move it.
                            ! No action needed here, as the ELSE block below contains the move operation.
                        ELSE IF (ABS(vcl3d(1,min_dist_idx) - domain%xmax) < TOL_PROF .AND. &
                                    ABS(vcl3d(2,min_dist_idx) - (TAN_30_DEG * domain%xmax)) < TOL_PROF) THEN
                            ! This point (min_dist_idx) is a protected "pit centre", do not move it.
                            ! No action needed here either.
                        ELSE
                            ! If min_dist_idx is NOT a protected boundary point, update its coordinates to the midpoint.
                            vcl3d(1,min_dist_idx) = newPoint(1)
                            vcl3d(2,min_dist_idx) = newPoint(2)
                            vcl3d(3,min_dist_idx) = newPoint(3)
                            print *, 'Coarsened New points with coordinates ' , vcl3d(1,min_dist_idx), vcl3d(2,min_dist_idx), vcl3d(3,min_dist_idx)
                            ! Update node structure for min_dist_idx as well
                            node(min_dist_idx)%coords(1) = newPoint(1)
                            node(min_dist_idx)%coords(2) = newPoint(2)
                            node(min_dist_idx)%coords(3) = newPoint(3)
                        END IF
                       coarsenCounter = coarsenCounter + 1
                   END IF

                   IF (coarsenCounter >= maxCoarsenCounter) THEN
                        EXIT 
                    END IF
            end do NPTLOOPS
            print *, 'Coarsen Counter: ' , coarsenCounter
            print *, "New npt: " , npt

    ! Filter out NaN values from node coordinates and update VCL 
            temp_count = 0
            DO i = 1, npt
                    IF ((vcl3d(1,i) == vcl3d(1,i)) .AND. &
                        (vcl3d(2,i) == vcl3d(2,i)) .AND. &
                        (vcl3d(3,i) == vcl3d(3,i))) THEN
                            IF ((vcl3d(1,i) /= num_NAN) .AND. (vcl3d(2,i) /= num_NAN) .AND. (vcl3d(3,i) /= num_NAN)) THEN
                                temp_count = temp_count + 1
                                temp_x(temp_count) = vcl3d(1,i)
                                temp_y(temp_count) = vcl3d(2,i)
                                temp_z(temp_count) = vcl3d(3,i)
                            END IF
                    ELSE
                        PRINT *, "NaN value at ", i, ":", vcl3d(1,i), vcl3d(2,i), vcl3d(3,i)
                    END IF
            END DO
            npt = temp_count
            DO i = 1, npt
                vcl3d(1,i) = temp_x(i)
                vcl3d(2,i) = temp_y(i)
                vcl3d(3,i) = temp_z(i)
                node(i)%coords(1) = temp_x(i)
                node(i)%coords(2) = temp_y(i)
                node(i)%coords(3) = temp_z(i)
            END DO
            DO i = 1, npt
                vcl(1,i) = vcl3d(1,i)
                vcl(2,i) = vcl3d(2,i)
            END DO
            ! nfacept = npt 

    !-----------------------------------------------------------------------------
    !  Remesh Code
    !-----------------------------------------------------------------------------
    REMESHLOOP: DO k = 1, 1 ! Line 1510 matlab
        print *, " - - - - PRE REMESH LOOP #1 ITERATION ", k , " ----------"
        print *, "Number of points (npt): ", npt
        print *, "Number of triangles (ntri): ", ntri
        print *, "Number of face points (nfacept): ", nfacept
        ! Generate the Delaunay triangulation
        do i = 1, npt
            ind(i) = i
        end do
        call dtris2 (npt, maxnpt, vcl, ind, ntri, til, tnbr, stack, ierr)
        IF (ierr /= 0) THEN
            print *, "ERROR IN DTRIS2 - ERROR NUMBER " , ierr
        END IF
        !call triangulate_surface(domain, npt, vcl, vcl3d, ntri, til, nfacetri, cell)
        !call mesh_boundaries(domain, npt, vcl, vcl3d, ntri, til, cell, nfacept) !ASK DR POLK
        !nfacept = npt
        print *, " - - - - POST REMESH #1 ITERATION ", k , " ----------"
        print *, "Number of points (npt): ", npt
        print *, "Number of triangles (ntri): ", ntri
        print *, "Number of face points (nfacept): ", nfacept
        print *, "Max triangles: " , maxntri

        print *, "Clean up flagged elements"
        ! FLAGGED ELEMENTS CLEAN UP MESHING 
                ! Clean up the meshing: Assess the triangulation for improper elements
                ALLOCATE(flagged_elements(maxntri)) ! Allocate dynamically or use a fixed large size
                flagged_elements = 0 ! Initialize to 0 or some indicator that it's not flagged
                value_count = 0
                ! Re-evaluate `nfacetri` based on the new `ntri` from `dtris2` +  will store valid triangles in `new_til` and then copy back to `til`.
                new_count = 0
                NTRILOOP: DO i = 1, ntri
                    row_indices = til(:,i)
                    ! Check if the indices are within bounds of `vcl3d`. If `npt` has decreased, some `til` entries might point to old, removed indices. This check is crucial.
                    IF (row_indices(1) > npt .OR. row_indices(2) > npt .OR. row_indices(3) > npt) THEN
                        PRINT "(A, I5, A, I5, A, I5, A, I5)", "skip because npt = ", npt, " and we have ", row_indices(1), row_indices(2), row_indices(3)
                        CONTINUE ! Skip this triangle, its vertices are no longer valid 
                    END IF
                    x_val_1 = vcl3d(1, row_indices(1))
                    x_val_2 = vcl3d(1, row_indices(2))
                    x_val_3 = vcl3d(1, row_indices(3))
                    y_val_1 = vcl3d(2, row_indices(1))
                    y_val_2 = vcl3d(2, row_indices(2))
                    y_val_3 = vcl3d(2, row_indices(3))
                    ! Check if an element is composed entirely of nodes along the hypotenuse
                    is_on_hypotenuse = (ABS(y_val_1 - (x_val_1 * TAN_30_DEG)) < TOL_PROF) .AND. (ABS(y_val_2 - (x_val_2 * TAN_30_DEG)) < TOL_PROF) .AND. (ABS(y_val_3 - (x_val_3 * TAN_30_DEG)) < TOL_PROF)
                    IF (is_on_hypotenuse) THEN
                        value_count = value_count + 1
                        print *, "This triangle lies along a hypotenuse"
                        flagged_elements(value_count) = i ! Store the index of the flagged triangle
                    END IF
                    ! Check sum of angles is not more than 180
                            ! ! Extract node coordinates for the current triangle
                            ! pointVec(1,:) = vcl3d(:, row_indices(1))
                            ! pointVec(2,:) = vcl3d(:, row_indices(2))
                            ! pointVec(3,:) = vcl3d(:, row_indices(3))
                            ! ! Centroid check: `sqrt(centroid_new(1)^2+centroid_new(2)^2)<0.7207*domain_x_max`
                            ! current_centroid_coords(1) = (pointVec(1,1) + pointVec(2,1) + pointVec(3,1)) / 3.0d0
                            ! current_centroid_coords(2) = (pointVec(1,2) + pointVec(2,2) + pointVec(3,2)) / 3.0d0
                            ! current_centroid_coords(3) = (pointVec(1,3) + pointVec(2,3) + pointVec(3,3)) / 3.0d0
                            ! ! Originally 0.7207
                            ! IF (SQRT(current_centroid_coords(1)**2 + current_centroid_coords(2)**2) < (0.7207 * domain%xmax)) THEN ! MATCH WITH LINE 418 IN IOMESH - IMPORTANT
                            ! ! The earea outside the 0.55 boudnary is where the more important sputtering phenomena occurs, otherwise, too close tothe grid hole to oberve any sputtering effects
                            ! CONTINUE ! Skip if centroid is too close to the origin
                            ! END IF
                            ! ! Calculate normalized edge vectors
                            ! joinVectors(1,:) = pointVec(2,:) - pointVec(1,:)
                            ! joinVectors(2,:) = pointVec(3,:) - pointVec(2,:)
                            ! joinVectors(3,:) = pointVec(1,:) - pointVec(3,:)

                            ! joinVectors(1,:) = joinVectors(1,:) / norm2(joinVectors(1,:))
                            ! joinVectors(2,:) = joinVectors(2,:) / norm2(joinVectors(2,:))
                            ! joinVectors(3,:) = joinVectors(3,:) / norm2(joinVectors(3,:))

                            ! ! Calculate internal angles using dot product and arccos
                            ! angles(1) = ACOS(-DOT_PRODUCT(joinVectors(3,:), joinVectors(1,:)))
                            ! angles(2) = ACOS(-DOT_PRODUCT(joinVectors(1,:), joinVectors(2,:)))
                            ! angles(3) = ACOS(-DOT_PRODUCT(joinVectors(2,:), joinVectors(3,:)))
                            
                            ! IF (abs((angles(1) + angles(2) + angles(3)) - PI) > 1e-8) THEN
                            !     print *, 'NOTE: Triangle angles dont sum to pi for index ' , i , 'since we have ' , angles(1) , angles(2) , angles(3)
                            !     value_count = value_count + 1
                            !     flagged_elements(value_count) = i ! Store the index of the flagged triangle
                            ! END IF            
                END DO NTRILOOP! Store all proper elements into new_til
                
                new_count = 0
                DO i = 1, ntri ! Check if the current triangle index 'i' is in the list of flagged elements
                        is_flagged = .FALSE.
                        DO j = 1, value_count
                            IF (flagged_elements(j) == i) THEN
                                is_flagged = .TRUE.
                                print *, "Flagged triangle found"
                                EXIT
                            END IF
                        END DO
                        IF (.NOT. is_flagged) THEN
                            new_count = new_count + 1
                            new_til(:, new_count) = til(:, i)
                        END IF
                    END DO
                DEALLOCATE(flagged_elements) ! Replace the original triangulation with its modified counterpart
                ntri = new_count
                DO i = 1, ntri
                    til(:, i) = new_til(:, i)
                END DO
                nfacetri = ntri ! Update nfacetri to reflect the new number of surface triangles

        print *, "Filter out NaN values from node coordinates and update VCL "
        ! Filter out NaN values from node coordinates and update VCL 
                temp_count = 0
                DO i = 1, npt
                    IF ((vcl3d(1,i) == vcl3d(1,i)) .AND. &
                        (vcl3d(2,i) == vcl3d(2,i)) .AND. &
                        (vcl3d(3,i) == vcl3d(3,i))) THEN
                            IF ((vcl3d(1,i) /= num_NAN) .AND. (vcl3d(2,i) /= num_NAN) .AND. (vcl3d(3,i) /= num_NAN)) THEN
                                temp_count = temp_count + 1
                                temp_x(temp_count) = vcl3d(1,i)
                                temp_y(temp_count) = vcl3d(2,i)
                                temp_z(temp_count) = vcl3d(3,i)
                            END IF
                    ELSE
                        PRINT *, "NaN value at ", i, ":", vcl3d(1,i), vcl3d(2,i), vcl3d(3,i)
                    END IF
                END DO
                npt = temp_count
                DO i = 1, npt
                    vcl3d(1,i) = temp_x(i)
                    vcl3d(2,i) = temp_y(i)
                    vcl3d(3,i) = temp_z(i)
                    node(i)%coords(1) = temp_x(i)
                    node(i)%coords(2) = temp_y(i)
                    node(i)%coords(3) = temp_z(i)
                END DO
                DO i = 1, npt
                    vcl(1,i) = vcl3d(1,i)
                    vcl(2,i) = vcl3d(2,i)
                END DO
                ! nfacept = npt 

        print *, "Extract thruster grid points from triangulation after removing bad triangles (1573 Matlab)"
        ! Extract thruster grid points from triangulation after removing bad triangles (1573 Matlab)
        print *, "a"    
                ALLOCATE(point_map(maxnpt)) ! MAKES A BIG DIFFERENCE
                point_map = 0
                k_index = 0
                print *, "b"
                DO i = 1, ntri
                     print *, "c"
                    DO j = 1, 3
                        ! Add a boundary check for safety, as npt might have decreased
                        print *, "here again"
                        IF (til(j, i) <= npt) THEN
                            IF (point_map(til(j, i)) == 0) THEN
                                k_index = k_index + 1
                                point_map(til(j, i)) = k_index
                                print *, "We have " , point_map(til(j, i)), " = " , k_index
                            END IF
                        END IF
                    END DO
                END DO
                new_npt = k_index
                print *, "npt is " , npt , "and new npt is " , new_npt
                IF (new_npt < npt+1) THEN
                    ALLOCATE(vcl3d_temp(3, maxnpt), vcl_temp(2, maxnpt))
                    ALLOCATE(node_temp(maxnpt))
                    k_index = 0
                    DO i = 1, npt
                        IF (point_map(i) > 0) THEN
                            k_index = k_index + 1
                            vcl3d_temp(1, k_index) = vcl3d(1, i)
                            vcl3d_temp(2, k_index) = vcl3d(2, i)
                            vcl3d_temp(3, k_index) = vcl3d(3, i)
                            vcl_temp(1, k_index) = vcl(1, i)
                            vcl_temp(2, k_index) = vcl(2, i)
                            node_temp(k_index) = node(i)
                        END IF
                    END DO
                    npt = new_npt
                    !nfacept = npt 
                    
                    ! Loop through and copy data from the temporary arrays to the main arrays.
                    DO i = 1, npt
                        print *, "npt is " , npt , "and new npt is " , new_npt
                        print *,"i = " , i , "(" , vcl3d_temp(1, i), vcl3d_temp(3, i) , vcl3d_temp(2, i) , ")"
                        vcl3d(1, i) = vcl3d_temp(1, i)
                        vcl3d(2, i) = vcl3d_temp(2, i)
                        vcl3d(3, i) = vcl3d_temp(3, i)
                        vcl(1, i) = vcl_temp(1, i)
                        vcl(2, i) = vcl_temp(2, i)
                        node(i) = node_temp(i) ! ERROR HERE
                    END DO
                    DEALLOCATE(vcl3d_temp, vcl_temp, node_temp)
                END IF
                DEALLOCATE(point_map)
                
                ! Debug file update
                        OPEN(UNIT=53, FILE='DEBUGPOSTMESHONEB', ACTION='WRITE') 
                                    do i = 1, ntri
                                    write(53,*) til(1,i), til(2,i), til(3,i)
                                    end do
                        CLOSE(53) 
                        OPEN(UNIT=39, FILE='DEBUGPOSTMESHONEA', ACTION='WRITE') 
                            DO inode = 1, npt 
                                curr_x = node(inode)%coords(1) 
                                curr_y = node(inode)%coords(2) 
                                curr_z = node(inode)%coords(3)
                                WRITE(39, '(F10.2, "," , ES25.18, ",", ES25.18, ",", ES25.18)')erosionTime, curr_x, curr_y, curr_z
                            END DO 
                        CLOSE(39) 
        
        ! Debug prints
                print *, " - - - -  Post Flagged Refinement " , k , " - - - - "
                print *, "Number of points (npt): ", npt
                print *, "Number of triangles (ntri): ", ntri
                print *, "Number of face points (nfacept): ", nfacept
                max_til_idx = 0
                DO i = 1, ntri
                    max_til_idx = MAX(max_til_idx, til(1,i), til(2,i), til(3,i))
                END DO
                PRINT *, "Max index in TIL after boundaries: ", max_til_idx
                PRINT *, "Current npt (should match max_til_idx) after boundaires:  ", npt

        ! REMESH
                print *, " - - - - REMESH #2 LOOP ", k , " ----------"
                call triangulate_surface(domain, npt, vcl, vcl3d, ntri, til, nfacetri, cell)
                !call mesh_boundaries(domain, npt, vcl, vcl3d, ntri, til, cell, nfacept) !ASK DR POLK
                nfacept = npt
                print *, " - - - -  Post Remesh # 2 for " , k , " - - - - "
                print *, "Number of points (npt): ", npt
                print *, "Number of triangles (ntri): ", ntri
                print *, "Number of face points (nfacept): ", nfacept

        ! Debug file update
                OPEN(UNIT=28, FILE='DEBUGPOSTMESHTWOA', ACTION='WRITE') 
                            DO inode = 1, npt 
                                curr_x = node(inode)%coords(1) 
                                curr_y = node(inode)%coords(2) 
                                curr_z = node(inode)%coords(3)
                                WRITE(28, '(F10.2, "," , ES25.18, ",", ES25.18, ",", ES25.18)')erosionTime, curr_x, curr_y, curr_z
                            END DO 
                CLOSE(28) 
                OPEN(UNIT=35, FILE='DEBUGPOSTMESHTWOB', ACTION='WRITE') 
                    do i = 1, ntri
                    write(35,*) til(1,i), til(2,i), til(3,i)
                    end do
                CLOSE(35) 

        print *, "FLAGGED ELEMENTS CLEAN UP MESHING "
        ! FLAGGED ELEMENTS CLEAN UP MESHING 
                ! Clean up the meshing: Assess the triangulation for improper elements
                ALLOCATE(flagged_elements(maxntri)) ! Allocate dynamically or use a fixed large size
                flagged_elements = 0 ! Initialize to 0 or some indicator that it's not flagged
                value_count = 0
                ! Re-evaluate `nfacetri` based on the new `ntri` from `dtris2` +  will store valid triangles in `new_til` and then copy back to `til`.
                new_count = 0
                NTRILOOPTWO: DO i = 1, ntri
                    row_indices = til(:,i)
                    ! Check if the indices are within bounds of `vcl3d`. If `npt` has decreased, some `til` entries might point to old, removed indices. This check is crucial.
                    IF (row_indices(1) > npt .OR. row_indices(2) > npt .OR. row_indices(3) > npt) THEN
                        PRINT *, "skip because npt = ", npt, " and we have ", row_indices(1), row_indices(2), row_indices(3)
                        CONTINUE ! Skip this triangle, its vertices are no longer valid 
                    END IF
                    x_val_1 = vcl3d(1, row_indices(1))
                    x_val_2 = vcl3d(1, row_indices(2))
                    x_val_3 = vcl3d(1, row_indices(3))
                    y_val_1 = vcl3d(2, row_indices(1))
                    y_val_2 = vcl3d(2, row_indices(2))
                    y_val_3 = vcl3d(2, row_indices(3))
                    ! Check if an element is composed entirely of nodes along the hypotenuse
                    is_on_hypotenuse = (ABS(y_val_1 - (x_val_1 * TAN_30_DEG)) < TOL_PROF) .AND. (ABS(y_val_2 - (x_val_2 * TAN_30_DEG)) < TOL_PROF) .AND. (ABS(y_val_3 - (x_val_3 * TAN_30_DEG)) < TOL_PROF)
                    IF (is_on_hypotenuse) THEN
                        value_count = value_count + 1
                        flagged_elements(value_count) = i ! Store the index of the flagged triangle
                    END IF

                    ! ! Check sum of angles is not more than 180
                    !         ! Extract node coordinates for the current triangle
                    !         pointVec(1,:) = vcl3d(:, row_indices(1))
                    !         pointVec(2,:) = vcl3d(:, row_indices(2))
                    !         pointVec(3,:) = vcl3d(:, row_indices(3))
                    !         ! Centroid check: `sqrt(centroid_new(1)^2+centroid_new(2)^2)<0.7207*domain_x_max`
                    !         current_centroid_coords(1) = (pointVec(1,1) + pointVec(2,1) + pointVec(3,1)) / 3.0d0
                    !         current_centroid_coords(2) = (pointVec(1,2) + pointVec(2,2) + pointVec(3,2)) / 3.0d0
                    !         current_centroid_coords(3) = (pointVec(1,3) + pointVec(2,3) + pointVec(3,3)) / 3.0d0
                    !         ! Originally 0.7207
                    !         IF (SQRT(current_centroid_coords(1)**2 + current_centroid_coords(2)**2) < (0.7207 * domain%xmax)) THEN ! MATCH WITH LINE 418 IN IOMESH - IMPORTANT
                    !         ! The earea outside the 0.55 boudnary is where the more important sputtering phenomena occurs, otherwise, too close tothe grid hole to oberve any sputtering effects
                    !         CONTINUE ! Skip if centroid is too close to the origin
                    !         END IF
                    !         ! Calculate normalized edge vectors
                    !         joinVectors(1,:) = pointVec(2,:) - pointVec(1,:)
                    !         joinVectors(2,:) = pointVec(3,:) - pointVec(2,:)
                    !         joinVectors(3,:) = pointVec(1,:) - pointVec(3,:)

                    !         joinVectors(1,:) = joinVectors(1,:) / norm2(joinVectors(1,:))
                    !         joinVectors(2,:) = joinVectors(2,:) / norm2(joinVectors(2,:))
                    !         joinVectors(3,:) = joinVectors(3,:) / norm2(joinVectors(3,:))

                    !         ! Calculate internal angles using dot product and arccos
                    !         angles(1) = ACOS(-DOT_PRODUCT(joinVectors(3,:), joinVectors(1,:)))
                    !         angles(2) = ACOS(-DOT_PRODUCT(joinVectors(1,:), joinVectors(2,:)))
                    !         angles(3) = ACOS(-DOT_PRODUCT(joinVectors(2,:), joinVectors(3,:)))
                            
                    !         IF (abs((angles(1) + angles(2) + angles(3)) - PI) > 1e-8) THEN
                    !             print *, 'NOTE: Triangle angles dont sum to pi for index ' , i , 'since we have ' , angles(1) , angles(2) , angles(3)
                    !             value_count = value_count + 1
                    !             flagged_elements(value_count) = i ! Store the index of the flagged triangle
                    !         END IF            
                END DO NTRILOOPTWO! Store all proper elements into new_til

                new_count = 0
                DO i = 1, ntri ! Check if the current triangle index 'i' is in the list of flagged elements
                        is_flagged = .FALSE.
                        DO j = 1, value_count
                            IF (flagged_elements(j) == i) THEN
                                is_flagged = .TRUE.
                                print *, "Flagged triangle found"
                                EXIT
                            END IF
                        END DO
                        IF (.NOT. is_flagged) THEN
                            new_count = new_count + 1
                            new_til(:, new_count) = til(:, i)
                        END IF
                    END DO
                DEALLOCATE(flagged_elements) ! Replace the original triangulation with its modified counterpart
                ntri = new_count
                DO i = 1, ntri
                    til(:, i) = new_til(:, i)
                END DO
                nfacetri = ntri ! Update nfacetri to reflect the new number of surface triangles

        print *, "Filter out NaN values from node coordinates and update VCL "
        ! Filter out NaN values from node coordinates and update VCL 
                temp_count = 0
                DO i = 1, npt
                    IF ((vcl3d(1,i) == vcl3d(1,i)) .AND. &
                        (vcl3d(2,i) == vcl3d(2,i)) .AND. &
                        (vcl3d(3,i) == vcl3d(3,i))) THEN
                        temp_count = temp_count + 1
                        temp_x(temp_count) = vcl3d(1,i)
                        temp_y(temp_count) = vcl3d(2,i)
                        temp_z(temp_count) = vcl3d(3,i)
                    ELSE
                        PRINT *, "NaN value at ", i, ":", vcl3d(1,i), vcl3d(2,i), vcl3d(3,i)
                    END IF
                END DO
                npt = temp_count
                DO i = 1, npt
                    vcl3d(1,i) = temp_x(i)
                    vcl3d(2,i) = temp_y(i)
                    vcl3d(3,i) = temp_z(i)
                    node(i)%coords(1) = temp_x(i)
                    node(i)%coords(2) = temp_y(i)
                    node(i)%coords(3) = temp_z(i)
                END DO
                DO i = 1, npt
                    vcl(1,i) = vcl3d(1,i)
                    vcl(2,i) = vcl3d(2,i)
                END DO
                ! nfacept = npt 
        
        print *, "UPDATE Files with new data"
        ! UPDATE Files with new data
                ! Now, update cell parameters and node parameters for the new mesh
                ! Recalculate cell parameters
                temp_count = 0
                DO i = 1, npt
                    !print *, (vcl3d(1,i) == vcl3d(1,i)) , (vcl3d(2,i) == vcl3d(2,i)) ,  (vcl3d(3,i) == vcl3d(3,i)), ' 0 = ? = ' , vcl3d(1,i), vcl3d(2,i) , vcl3d(3,i) , ((vcl3d(1,i) /= num_NAN) .AND. (vcl3d(2,i) /= num_NAN) .AND. (vcl3d(3,i) /= num_NAN))
                    IF ((vcl3d(1,i) == vcl3d(1,i)) .AND. (vcl3d(2,i) == vcl3d(2,i)) .AND. (vcl3d(3,i) == vcl3d(3,i))) THEN
                        IF ((vcl3d(1,i) /= num_NAN) .AND. (vcl3d(2,i) /= num_NAN) .AND. (vcl3d(3,i) /= num_NAN)) THEN
                            temp_count = temp_count + 1
                            temp_x(temp_count) = vcl3d(1,i)
                            temp_y(temp_count) = vcl3d(2,i)
                            temp_z(temp_count) = vcl3d(3,i)
                            !print *, 'VCL3d is ' , temp_x(temp_count), temp_y(temp_count), temp_z(temp_count)
                        END IF 
                    END IF 
                END DO
                npt = temp_count
                DO i = 1, npt
                    vcl3d(1,i) = temp_x(i)
                    vcl3d(2,i) = temp_y(i)
                    vcl3d(3,i) = temp_z(i)
                    node(i)%coords(1) = temp_x(i)
                    node(i)%coords(2) = temp_y(i)
                    node(i)%coords(3) = temp_z(i)
                END DO
                DO i = 1, npt
                    vcl(1,i) = vcl3d(1,i)
                    vcl(2,i) = vcl3d(2,i)
                END DO
                do icell = 1, ntri
                    ! Define the triangle incidence list for cell
                    cell(icell)%til(:) = til(:,icell)
                    ! Define vertices of cell
                    cell(icell)%vert1(:) = vcl3d(:,cell(icell)%til(1))
                    cell(icell)%vert2(:) = vcl3d(:,cell(icell)%til(2))
                    cell(icell)%vert3(:) = vcl3d(:,cell(icell)%til(3))
                    ! Calculate other cell parameters
                    call calc_cell_params(cell(icell), icell)
                end do
                ! Recalculate node parameters for surface nodes

                ! Debug file update
                OPEN(UNIT=54, FILE='DEBUGPOSTLOOPA', ACTION='WRITE') 
                            DO inode = 1, npt 
                                curr_x = node(inode)%coords(1) 
                                curr_y = node(inode)%coords(2) 
                                curr_z = node(inode)%coords(3)
                                WRITE(54, '(F10.2, "," , ES25.18, ",", ES25.18, ",", ES25.18)')erosionTime, curr_x, curr_y, curr_z
                            END DO 
                CLOSE(54) 
                OPEN(UNIT=56, FILE='DEBUGPOSTLOOPB', ACTION='WRITE') 
                    do i = 1, ntri
                    write(56,*) til(1,i), til(2,i), til(3,i)
                    end do
                CLOSE(56) 

                DO i = 1, nfacept ! Only surface nodes have these properties
                    print *, "nfacept is " , nfacept
                    print *, "size of node is " , SIZE(node)
                    node(i)%coords(:) = vcl3d(:, i)
                    CALL calc_node_params(node(i), i, domain, nfacetri, cell)
                END DO
                print *, "after"
        
        ! Debug prints
            print *, " - - - -  End loop " , k , " - - - - "
            print *, "Number of points (npt): ", npt
            print *, "Number of triangles (ntri): ", ntri
            print *, "Number of face points (nfacept): ", nfacept
            max_til_idx = 0
            DO i = 1, ntri
                max_til_idx = MAX(max_til_idx, til(1,i), til(2,i), til(3,i))
            END DO
            PRINT *, "Max index in TIL after boundaries: ", max_til_idx
            PRINT *, "Current npt (should match max_til_idx) after boundaires:  ", npt

    END DO REMESHLOOP

    !---------------------------------------------------------------------------
    ! Save the vertex coordinate list and mesh triangle index list for plotting
    !---------------------------------------------------------------------------
    !     open(13, FILE = 'VCL') 
    !     do i = 1, npt
    !       write(13,*) i, vcl(1,i), vcl(2,i)
    !     end do
    !     close(13)

    ! open(12, FILE = 'Triangles')
    !         do i = 1, ntri
    !           write(12,*) til(1,i), til(2,i), til(3,i)
    !         end do
    !         close(12)
   
    !---------------------------------------------------------------------------
    ! PSB Save 3D mesh data for Python plotting (now reads VCL3D and Triangles directly)
    !---------------------------------------------------------------------------
    ! No need to open/write mesh_points_3d.txt and mesh_triangles_3d.txt anymore
    ! Since the Python script will read VCL3D and Triangles
    plot_subfolder = 'mesh_plots'
    OPEN(UNIT=13, FILE='VCL3D', ACTION='WRITE') 
        WRITE(13, '(A)') '# ErosionTime, X, Y, Z'
        DO inode = 1, nfacept 
          curr_x = node(inode)%coords(1) 
          curr_y = node(inode)%coords(2) 
          curr_z = node(inode)%coords(3)
          WRITE(13, '(F10.2, "," , ES25.18, ",", ES25.18, ",", ES25.18)')erosionTime, curr_x, curr_y, curr_z
        END DO 
    CLOSE(13) 
    open(12, FILE = 'Triangles')
            do i = 1, ntri
              write(12,*) til(1,i), til(2,i), til(3,i)
            end do
            close(12)
    CALL EXECUTE_COMMAND_LINE("mkdir -p " // TRIM(plot_subfolder), wait=.true., exitstat=j)
    IF (j /= 0) THEN
        PRINT *, "Warning: Could not create directory ", TRIM(plot_subfolder), ". Error code: ", j
        PRINT *, "Plot will be saved in the current directory instead."
        plot_subfolder = '.'
    END IF
    WRITE(step_str, '(I0)') erosionStep ! Convert erosionStep to string
    plot_name_str = "remesh_step_" // TRIM(step_str) ! Identifier for the plot filename
    ! Command to call the Python script, passing the output folder and plot identifier
    command_string = "python plot_mesh.py " // TRIM(plot_subfolder) // " " // plot_name_str
    CALL EXECUTE_COMMAND_LINE(TRIM(command_string), wait=.true.)
  
END DO EROSION

contains
subroutine randperm(x, n, k)
 !******************************************************************************
 ! Returns a random permutation of k out of n objects
 ! Based on Green, Behavior Research Methods and Instrum., 1977, Vol 9 (6), 559
 !
 !******************************************************************************
 implicit none
 integer, intent(inout), dimension(*) :: x
 integer, intent(in) :: n  ! total number of objects
 integer, intent(in) :: k  ! number of objects in random permutation
 integer :: i, j
 real*8 :: r, t

 do j = 1, k
   call random_seed()
   call random_number(r)
   i = int((float(n - j + 1)) * r) + j
   t = x(j)
   x(j) = x(i)
   x(i) = t
 end do
 return
end subroutine randperm

 ! >>> PSB WROTE THIS SUBROUTINE: trace_ion_path <<<
!******************************************************************************
! RECURSIVE SUBROUTINE: trace_ion_path
! This subroutine simulates the path of a single ion through the simulation
! domain, handling reflections off symmetry boundaries and calculating
! sputtering events. It uses recursion to manage the reflection chain.
!
! Arguments:
!   iion_idx          (INTENT(IN))    : The original index of the ion in the sampled bundle (1 to nionsample).
!   current_refl_count(INTENT(IN))    : The number of reflections this ion has already undergone in this path trace.
!   ion_current_origin(INTENT(INOUT)) : The current 3D position of the ion (updates on reflection).
!   ion_current_dir   (INTENT(INOUT)) : The current 3D direction vector of the ion (updates on reflection).
!   nlostboys         (INTENT(INOUT)) : Counter for ions that leave the domain without hitting anything.
!   nbottomboys       (INTENT(INOUT)) : Counter for ions that hit the bottom boundary.
!   ion_refl_counts_arr(INTENT(INOUT)): Array to store total reflections for each sampled ion.
!   vert1, vert2, vert3 (INTENT(IN))  : 3D coordinates of triangle vertices (defines the mesh).
!   t, u, v, xcoor, intersect (INTENT(OUT)): Outputs from triangle_ray_intersection.
!                                         Note: these are temporary and get overwritten in recursive calls.
!   cell              (INTENT(INOUT)) : Array of cell_type structures (contains mass_loss accumulators, normals, etc.).
!                                         Mass loss is accumulated here.
!   node              (INTENT(INOUT)) : Array of node_type structures (contains mass_loss accumulators).
!                                         Mass loss is accumulated here.
!   til               (INTENT(IN))    : Triangle Incidence List (node indices for each triangle).
!   max_reflections   (INTENT(IN))    : Maximum allowed reflections for an ion.
!
! External Functions/Subroutines Used:
!   triangle_ray_intersection (from iomesh_mod or sputter_mod)
!   count (Fortran intrinsic)
!   findloc (Fortran intrinsic)
!   dot (Fortran intrinsic)
!   angle (from sputter_mod)
!   sputter_yield_carbon (from sputter_mod)
!******************************************************************************
RECURSIVE SUBROUTINE trace_ion_path(iion_idx, current_refl_count, &
                                   ion_current_origin, ion_current_dir, &
                                   nlostboys, nbottomboys, ion_refl_counts_arr, &
                                   vert1, vert2, vert3, t, u, v, xcoor, intersect, &
                                   cell, node, til, max_reflections)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: iion_idx
    INTEGER, INTENT(IN) :: current_refl_count
    REAL*8, DIMENSION(3), INTENT(INOUT) :: ion_current_origin ! Position can update on reflection
    REAL*8, DIMENSION(3), INTENT(INOUT) :: ion_current_dir    ! Direction updates on reflection
    INTEGER, INTENT(INOUT) :: nlostboys
    INTEGER, INTENT(INOUT) :: nbottomboys
    INTEGER, DIMENSION(*) , INTENT(INOUT) :: ion_refl_counts_arr ! Array to update for this ion_idx

    ! Mesh data (passed from main program)
    REAL*8, DIMENSION(3, maxntri), INTENT(IN) :: vert1, vert2, vert3
    INTEGER, DIMENSION(3, maxntri), INTENT(IN) :: til
    TYPE(cell_type), DIMENSION(maxntri), INTENT(INOUT) :: cell ! Mass loss accumulated
    TYPE(node_type), DIMENSION(maxnpt), INTENT(INOUT) :: node   ! Mass loss accumulated

    ! Outputs from triangle_ray_intersection (local to this call, or passed through)
    REAL*8, DIMENSION(maxntri), INTENT(OUT) :: t, u, v
    REAL*8, DIMENSION(3, maxntri), INTENT(OUT) :: xcoor
    LOGICAL, DIMENSION(maxntri), INTENT(OUT) :: intersect

    ! Parameters
    INTEGER, INTENT(IN) :: max_reflections
    ! Accessing global parameters from host scope (Redepo program)
    ! These are parameters, so they are constant and can be accessed directly.
    ! If they were variables, they would need to be passed as arguments or put in a module.
    REAL*8, PARAMETER :: e = 1.60217662E-19
    REAL*8, PARAMETER :: massCarbon = 1.994307e-23

    ! Local variables
    INTEGER :: num_hits
    INTEGER, DIMENSION(1) :: impact_tri_arr
    INTEGER :: impact_cell
    REAL*8 :: inc_angle
    REAL*8 :: yield, mass_loss
    REAL*8, DIMENSION(3) :: reflectedray

    ! For files, psb
    REAL*8, DIMENSION(3) :: node1_old_coords, node2_old_coords, node3_old_coords 
    REAL*8, DIMENSION(3) :: node1_new_coords, node2_new_coords, node3_new_coords ! <--- NEW DECLARATION 
    INTEGER :: n1_idx, n2_idx, n3_idx

    ! Debug print
    ! print *, "trace_ion_path was called"

    ! Base Case 1: Exceeded maximum reflections
    IF (current_refl_count >= max_reflections) THEN
        ! This ion's path is terminated as it has reflected too many times, treat as a "lost boy" for practical purposes 
        WRITE(112,"(2I8,6E15.6)") iion_idx, current_refl_count, ion_current_origin, ion_current_dir
        ! Debug print
        print *, "Max reflection count achieved for " , iion_idx
        RETURN ! Exit the recursive call
    END IF

    ! Perform ray-triangle intersection for the current ion state
    CALL triangle_ray_intersection(iion_idx, ion_current_origin, ion_current_dir, &
                                  vert1, vert2, vert3, intersect, t, u, v, xcoor, &
                                  borderIn = 'inclusive')

    num_hits = COUNT(intersect)
    ! Debug print
    ! print *, "numhits for " , iion_idx ," is ", num_hits
    ! Find the closest intersected triangle.
    ! Assuming triangle_ray_intersection or a subsequent sort ensures impact_tri_arr(1) is the closest.
    ! A more robust check would be to find the minimum 't' among all true 'intersect' entries.
    IF (num_hits > 0) THEN
        impact_tri_arr = FINDLOC(intersect, VALUE = .TRUE.)
        impact_cell = impact_tri_arr(1) ! Take the first one found (assumed closest)
    ELSE
        ! No hits in the current step. This is a base case.
        ! Ion did not strike any surface in domain (small fraction of total at edges or with vz near zero)
        WRITE(112,"(A,2I8,6E15.6)") "NO_HIT", iion_idx, current_refl_count, ion_current_origin, ion_current_dir   
        PRINT*, "Lost, No hit "
        nlostboys = nlostboys + 1
        RETURN ! Exit the recursive call
    END IF

    ! Base Case / Recursive Step Logic based on impact type
    IF (cell(impact_cell)%on_bndry == -1) THEN
        ! Base Case 2: Ion left through top of domain (should not happen)
        WRITE(112,"(A,2I8,6E15.6)") "TOP_EXIT", iion_idx, current_refl_count, ion_current_origin, ion_current_dir, impact_cell  
        PRINT*, "ion number ", iion_idx, "left through top of domain (should not happen)"
        RETURN ! Exit the recursive call
    ELSE IF (cell(impact_cell)%on_bndry == -2) THEN
        ! Base Case 3: Ion hit bottom of domain (should only happen after punch-through)
        WRITE(112,"(A,2I8,6E15.6)") "PUNCHTHRU", iion_idx, current_refl_count, ion_current_origin, ion_current_dir, impact_cell  
        nbottomboys = nbottomboys + 1
        PRINT*, "Hit bottom."
        RETURN ! Exit the recursive call
    ELSE IF (cell(impact_cell)%on_bndry == 0) THEN
        ! Base Case 4: Ion hit sputtering surface (calculate sputtered mass and share to nodes)

        ! --- Capture old node coordinates BEFORE mass loss is applied --- 
      n1_idx = til(1, impact_cell) 
      n2_idx = til(2, impact_cell) 
      n3_idx = til(3, impact_cell) 
      node1_old_coords = node(n1_idx)%coords 
      node2_old_coords = node(n2_idx)%coords 
      node3_old_coords = node(n3_idx)%coords 
      ! -----------------------------------------------------------------

        inc_angle = angle(ion_current_dir(:), -1.0d0 * cell(impact_cell)%normal(:))
        yield = sputter_yield_carbon(ion(iion_idx)%KE, inc_angle) ! Use original ion's KE
        mass_loss = yield * ion(iion_idx)%qCEXion/e * massCarbon ! Use original ion's macro-charge

        ! Accumulate mass loss for the impacted cell
        cell(impact_cell)%mass_loss = cell(impact_cell)%mass_loss + mass_loss

        ! Distribute mass loss to the three nodes of the impacted triangle using barycentric coordinates
        ! Note: u and v are from the *current* intersection (impact_cell)
        node(til(2,impact_cell))%mass_loss = node(til(2,impact_cell))%mass_loss + u(impact_cell) * mass_loss
        node(til(3,impact_cell))%mass_loss = node(til(3,impact_cell))%mass_loss + v(impact_cell) * mass_loss
        node(til(1,impact_cell))%mass_loss = node(til(1,impact_cell))%mass_loss + (1 - u(impact_cell) - v(impact_cell)) * mass_loss
        sum = node(til(2,impact_cell))%mass_loss + node(til(3,impact_cell))%mass_loss + node(til(1,impact_cell))%mass_loss

        !PRINT*, "sp"
        !PRINT*, "Sputtered | Ion Index:", iion_idx, "Reflections:", current_refl_count, &
        !"| Node", til(1,impact_cell), "Loss:", (1 - u(impact_cell) - v(impact_cell)) * mass_loss, &
        !"| Node", til(2,impact_cell), "Loss:", u(impact_cell) * mass_loss, &
        !"| Node", til(3,impact_cell), "Loss:", v(impact_cell) * mass_loss
        WRITE(112,"(2I8,6E15.6)") iion_idx, current_refl_count, ion_current_origin, ion_current_dir, impact_cell
        
        ! --- Write detailed impact data to file --- 
        ! OPEN(UNIT=ION_IMPACT_UNIT, FILE=ion_impact_filename, STATUS='OLD', ACTION='WRITE', POSITION='APPEND') 
        ! WRITE(ION_IMPACT_UNIT, '(I0, ",           ", I0, ",        ", I0, ",                ", 3(I0, ","), ' // & 
        !   '3("            ", ES26.18, ","), 3("            ", ES26.18, ","), 3("            ", ES26.18, ","), ' // & 
        !   '" ", ES26.18, ",", " ", ES26.18, ",", " ", ES26.18)') & 
        !   iion_idx, ionindex(iion_idx), impact_cell, &
        !   n1_idx, n2_idx, n3_idx, & 
        !   node1_old_coords(1), node1_old_coords(2), node1_old_coords(3), & 
        !   node2_old_coords(1), node2_old_coords(2), node2_old_coords(3), & 
        !   node3_old_coords(1), node3_old_coords(2), node3_old_coords(3), & 
        !   xcoor(1, impact_cell), xcoor(2, impact_cell), xcoor(3, impact_cell)
        !   CLOSE(ION_IMPACT_UNIT) 


        ! ------------------------------------------

        RETURN ! Exit the recursive call (ion's path is resolved)

    ELSE !PSB changed from ELSE IF (cell(impact_cell)%on_bndry > 0) THEN
        !print *, "-------------------------------- REFLECTION ACHIEVED ----------------------------------------"
        ! Recursive Step: Ion hit symmetry boundary (reflect)
        ! Increment reflection count for this specific ion
        ion_refl_counts_arr(iion_idx) = ion_refl_counts_arr(iion_idx) + 1

        ! Calculate reflected ray direction
        reflectedray = ion_current_dir - 2 * DOT_PRODUCT(ion_current_dir, cell(impact_cell)%normal) * cell(impact_cell)%normal

        ! Update ion's current direction. Origin remains at the intersection point.
        ! For the next recursive call, the origin should be slightly offset from the surface
        ! to prevent immediate re-intersection with the same triangle due to floating point precision.
        ! IE: The idea of stepping away from the wall so that the ray stays just clear of the surface it just reflected from and doesn't keep bouncing
        ! A small epsilon along the new direction is common practice.
        ion_current_origin = xcoor(:, impact_cell) + reflectedray * 1.0d-6 ! Move slightly off surface
        ion_current_dir = reflectedray

        ! Recursive Call: Continue tracing the ion's path from the new origin and direction
        CALL trace_ion_path(iion_idx, current_refl_count + 1, &
                            ion_current_origin, ion_current_dir, &
                            nlostboys, nbottomboys, ion_refl_counts_arr, &
                            vert1, vert2, vert3, t, u, v, xcoor, intersect, &
                            cell, node, til, max_reflections)
        RETURN ! Exit this recursive call (the result comes from the deeper call)
    END IF
  END SUBROUTINE trace_ion_path

    ! This subroutine should be placed within your redepo file, not in a module.
    SUBROUTINE clean_vcl3d(redepo_vcl, redepo_npt, maxnptval)
    ! Cleans up the redepo_vcl array in-place by removing duplicate points and points at (0,0).
    !
    ! Arguments:
    ! redepo_vcl: An input/output array of real numbers representing 3D coordinates.
    ! redepo_npt:   An input/output integer representing the number of points.
    ! maxnptval:    The maximum number of points allowed in the array.

    IMPLICIT NONE

    REAL*8, PARAMETER :: TOLERANCE = 1.0E-9
    INTEGER, INTENT(IN) :: maxnptval
    REAL*8, INTENT(INOUT), DIMENSION(3, maxnptval) :: redepo_vcl
    INTEGER, INTENT(INOUT) :: redepo_npt

    INTEGER :: m, j
    INTEGER :: write_idx
    LOGICAL :: is_duplicate, is_zero

    print *, "clean_vcl3d called"

    ! Use a new index to write valid points to, effectively compacting the array.
    write_idx = 0

    OUTER_LOOP: DO m = 1, redepo_npt
        ! Check for (0,0,0) point
        is_zero = (ABS(redepo_vcl(1, m)) < TOLERANCE) .AND. &
                  (ABS(redepo_vcl(2, m)) < TOLERANCE) .AND. &
        IF (is_zero) THEN
            CYCLE OUTER_LOOP
        END IF

        ! Check for duplicates against points already written to the beginning of the array.
        is_duplicate = .FALSE.
        DO j = 1, write_idx
            IF ((ABS(redepo_vcl(1, m) - redepo_vcl(1, j)) < TOLERANCE) .AND. &
                (ABS(redepo_vcl(2, m) - redepo_vcl(2, j)) < TOLERANCE)) THEN
                is_duplicate = .TRUE.
                EXIT
            END IF
        END DO

        ! If the point is not a duplicate, copy it to the next available position.
        IF (.NOT. is_duplicate) THEN
            write_idx = write_idx + 1
            ! Only copy if we need to; if m == write_idx, the point is already in place.
            IF (m /= write_idx) THEN
                redepo_vcl(:, write_idx) = redepo_vcl(:, m)
            END IF
        END IF
    END DO OUTER_LOOP

    ! The number of valid points is the final value of write_idx.
    redepo_npt = write_idx
    print *, "new nodes is " , redepo_npt

END SUBROUTINE clean_vcl3d
  end program Redepo


  