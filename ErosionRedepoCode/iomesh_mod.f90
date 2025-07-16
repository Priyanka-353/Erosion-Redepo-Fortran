module iomesh_mod
  use sputter_mod
  implicit none

  ! private
  ! public geom_type, ion_type, cell_type

  type geom_type
    integer :: type
    real*8 :: hole_diam
    real*8 :: xmin
    real*8 :: xmax
    real*8 :: ymin
    real*8 :: ymax
    real*8 :: zmin
    real*8 :: zmax
  end type

  type ion_type
    real*8, dimension(3) :: birth_loc
    real*8, dimension(3) :: death_loc
    real*8, dimension(3) :: origin
    real*8, dimension(3) :: vel
    real*8, dimension(3) :: dir
    real*8 :: KE
    real*8 :: qCEXion
    real*8 :: xydisp
    integer :: gridHit
    integer :: source
    integer :: index
  end type

  type cell_type
    real*8, dimension(3) :: til
    real*8, dimension(3) :: vert1
    real*8, dimension(3) :: vert2
    real*8, dimension(3) :: vert3
    real*8, dimension(3) :: edge1
    real*8, dimension(3) :: edge2
    real*8, dimension(3) :: edge3
    real*8, dimension(3) :: edgelength
    real*8, dimension(3) :: vertangle
    real*8 :: area
    real*8 :: mass_loss
    real*8, dimension(3) :: normal
    real*8, dimension(3) :: centroid
    integer :: on_bndry   ! -2 = bottom face, -1 = top face,  0 = accel face, 1 = front major wall,
        !  2 = right minor wall, 3 = hypotenuse wall or rear major wall, 4 = left major wall
  end type

  type node_type
    real*8, dimension(3) :: coords
    real*8, dimension(3) :: normal
    real*8 :: area  ! area associated with node (1/3 of area of cells containing node)
    integer, dimension(12) :: in_cells
    integer :: on_bndry   !  0 = none, 1 = hole, 2 = major wall, 3 = right minor wall,
        !  4 = hypotenuse or upper major wall, 5 = left minor wall
    integer :: on_corner  !  0 = none, 1 = lower left (major left minor intersection),
        !  2 = lower right (major minor intersection), 3 = upper right (minor-hypotenuse intersection
        !  or upper major right minor intersection), 4 = upper left (upper major left minor intersection)
    real*8, dimension(12) :: cell_angle
    real*8 :: mass_loss  !  mass loss assigned the node in a given time step
    real*8, dimension(3) :: vel  !  current node recession velocity
    real*8, dimension(3) :: vel_old   ! node recession velocity at previous time step
    integer :: firstStep   !  flag for first step in geometry update (use Euler integration, otherwise Adams-Bashforth)
  end type

  type ray_type
    real*8, dimension(3) :: coords
    real*8 :: weight
    real*8 :: polar_angle, azim_angle
  end type

  type (geom_type) :: domain
  type (ion_type), dimension(:), allocatable :: ion
  integer, dimension(:), allocatable :: ionindex

  real*8 :: simulation_time
  integer :: npandgions

  integer, parameter :: maxnpt = 300 ! 300 is maximum number of vertices in mesh--change if necessary
  integer, parameter :: maxntri = 600 ! 600 is maximum number of triangles in mesh--change if necessary
  integer, parameter :: maxnray = 2030 ! max number of Lebedev rays
  type (node_type), dimension(maxnpt) :: node
  type (cell_type), dimension(maxntri) :: cell
  type (ray_type), dimension(maxnray) :: ray
  real*8, dimension(2, maxnpt) :: vcl
  real*8, dimension(3, maxnpt) :: vcl3d
  ! real*8, dimension(3, maxnpt) :: v_old
  integer, dimension(3, maxntri) :: til
  integer :: npt, nfacept
  integer :: ntri, nfacetri

  contains
  !*******************************************************************************
    subroutine read_lebedev_file(filename, ray, nray)
    !-----------------------------------------------------------------------------
    !  Read in angles and weights for Lebedev rays
    !
    !  History:
    !
    !    Jay Polk
    !    Mar 17, 2022
    !
    !  Parameters:
    !
    !    Input, character (len = 50), filename, the name of the Lebedev ray datafile
    !    Output, type (ray_type) ray, a data structure containing geometry info for the domain
    !    Output, integer, nray, the number of Lebedev rays read out of file
    !-----------------------------------------------------------------------------

      character(len = 50), intent(in) :: filename
      type (ray_type), intent(inout) :: ray(:) !PSB
      integer, intent(out) :: nray
      integer :: i
      integer :: IOstatus

      !---------------------------------------------------------------------------
      !  Read data from specified file
      !---------------------------------------------------------------------------
      open (unit = 10, file = filename,form = 'UNFORMATTED')
      rewind(10)
      i = 1
      IOstatus = 0
      do
        read(10, IOSTAT=IOstatus) ray(i)%azim_angle, ray(i)%polar_angle, ray(i)%weight 
        if (IOstatus == 0) then ! Read until end of file !PSB
          if (ray(i)%polar_angle < 90d0) then  ! ignore entries for polar angles > 90; inside surface
            ray(i)%azim_angle = ray(i)%azim_angle*PI/180d0
            ray(i)%polar_angle = ray(i)%polar_angle*PI/180d0
            ray(i)%weight = ray(i)%weight * 2d0
            ray(i)%coords(1) = cos(ray(i)%azim_angle) * sin(ray(i)%polar_angle) !PSB
            ray(i)%coords(2) = sin(ray(i)%azim_angle) * sin(ray(i)%polar_angle) !PSB
            ray(i)%coords(3) = cos(ray(i)%polar_angle) !PSB
            i = i + 1
          end if

        else if (IOstatus > 0) then
          print*, "Error reading Lebedev ray points"
          exit ! PSB

        else ! End of file reached
          nray = i - 1 !PSB
          exit
        end if
      end do
      print *, "Read lebedev file"
      close(10)
    end subroutine read_lebedev_file

  !*******************************************************************************
    subroutine read_ion_file(filename, domain, simulation_time, npandgions)
    !-----------------------------------------------------------------------------
    !  Read in and/or define domain geometry, sim time, and number of ions, then
    !  populate ion data structure with data for each CEX ion hitting the accel
    !  grid.
    !
    !  History:
    !    Jay Polk
    !    May 6, 2021
    !
    !  Parameters:
    !    Input, character (len = 50), filename, the name of the ion datafile
    !    Output, type (geom_type) domain, a data structure containing geometry info for the domain
    !    Output (as global), type (ion_type) ion, a data structure containing CEX ion data from CEX3D. This
    !       is allocated dynamically in this subroutine after reading the number of ions, so it is created
    !       at the top of the module and accessed here globally.  In other subroutines it is passed as a
    !       parameter.
    !    Output, real (kind = 8), simulation_time, the CEX3D simulation time
    !    Output, integer, npandgions, the number of CEX ions striking the downstream face of the accel grid
    !-----------------------------------------------------------------------------

      character(len = *), intent(in) :: filename !PSB
      type (geom_type), intent(out) :: domain
      real*8, intent(out) :: simulation_time
      integer, intent(out) :: npandgions
      character(LEN = 50) :: actualname
      real*8 :: accel_thickness
      integer :: i

      !---------------------------------------------------------------------------
      !  Read data from specified file
      !---------------------------------------------------------------------------
      open (unit = 10, file = filename,form = 'UNFORMATTED')
      rewind(10)
      read(10) actualname
      print*, "CEX ion data extracted from ", actualname
      read(10) domain%zmax         ! accel_z
      read(10) domain%hole_diam    ! accel_hole_diameter
      read(10) accel_thickness
      domain%zmin = domain%zmax - accel_thickness
      read(10) domain%xmax
      read(10) simulation_time
      read(10) npandgions
      domain%xmin = 0.0d0
      domain%ymax = domain%xmax * tan(30.0d0 * pi/180.0d0)
      domain%ymin = 0.0d0

      print*, "zmax = ", domain%zmax
      print*, "hole_diam = ", domain%hole_diam
      print*, "accel_thickness = ", accel_thickness
      print*, "xmax = ", domain%xmax
      print*, "Number of CEX ions hitting in pits and grooves = ",npandgions
      print*, "Simulation_time = ", simulation_time

      !---------------------------------------------------------------------------
      !  Allocate array of ion structures--this should have visibility in module
      !  and any unit which uses this module
      !---------------------------------------------------------------------------
      allocate(ion(npandgions))
      allocate(ionindex(npandgions))

      !---------------------------------------------------------------------------
      !  Read data for CEX ions into the array of ion data structure
      !---------------------------------------------------------------------------
      i = 1
      !open(16, FILE = 'birthness')
      !open(17, FILE = 'deathness')
      !open(18, FILE = 'velness')
      !open(19, FILE = 'misc1')
      !open(20, FILE = 'misc2')
      do while (i <= npandgions)
        read(10) ion(i)%index, ion(i)%birth_loc(1), ion(i)%birth_loc(2), ion(i)%birth_loc(3), ion(i)%gridHit,  &
        ion(i)%death_loc(1), ion(i)%death_loc(2), ion(i)%death_loc(3), ion(i)%vel(1), ion(i)%vel(2), ion(i)%vel(3), &
        ion(i)%KE, ion(i)%qCEXion, ion(i)%source, ion(i)%xydisp, ion(i)%dir(1), ion(i)%dir(2), ion(i)%dir(3), &
        ion(i)%origin(1), ion(i)%origin(2), ion(i)%origin(3)
        if (ion(i)%death_loc(1) > 1.11d0) then
          print*, "death_loc > 1.11!!"
          !write(16,*) ion(i)%birth_loc(1), ion(i)%birth_loc(2), ion(i)%birth_loc(3)
          !write(19,*) ion(i)%gridHit, ion(i)%KE, ion(i)%qCEXion
          !write(20,*) ion(i)%source, ion(i)%xydisp
          !write(17,*) ion(i)%death_loc(1), ion(i)%death_loc(2), ion(i)%death_loc(3)
          !write(18,*) ion(i)%vel(1), ion(i)%vel(2), ion(i)%vel(3)
        end if
        i = i + 1
      end do
      !close(16)
      !close(17)
      !close(18)
      !close(19)
      !close(20)

      !PSB added the following code
      i = 1
      print*, "Ion #", i
      print*, "  Death Location (mm):       ", ion(i)%death_loc(1), ion(i)%death_loc(2), ion(i)%death_loc(3)
      print*, "  Velocity (mm/us):          ", ion(i)%vel(1), ion(i)%vel(2), ion(i)%vel(3)
      print*, "  Kinetic Energy (eV):       ", ion(i)%KE
      print*, "  Macro-particle Charge:     ", ion(i)%qCEXion
      print*, "  Source Region ID:          ", ion(i)%source
      print*, "  Radial Disp. (xydisp, mm): ", ion(i)%xydisp
      print*, "  Direction Vector:          ", ion(i)%dir(1), ion(i)%dir(2), ion(i)%dir(3)
      print*, "  Origin Location (mm):      ", ion(i)%origin(1), ion(i)%origin(2), ion(i)%origin(3)
      print*, "------------------------------------------------------------"
      i = npandgions
      print*, "Ion #", i
      print*, "  Death Location (mm):       ", ion(i)%death_loc(1), ion(i)%death_loc(2), ion(i)%death_loc(3)
      print*, "  Velocity (mm/us):          ", ion(i)%vel(1), ion(i)%vel(2), ion(i)%vel(3)
      print*, "  Kinetic Energy (eV):       ", ion(i)%KE
      print*, "  Macro-particle Charge:     ", ion(i)%qCEXion
      print*, "  Source Region ID:          ", ion(i)%source
      print*, "  Radial Disp. (xydisp, mm): ", ion(i)%xydisp
      print*, "  Direction Vector:          ", ion(i)%dir(1), ion(i)%dir(2), ion(i)%dir(3)
      print*, "  Origin Location (mm):      ", ion(i)%origin(1), ion(i)%origin(2), ion(i)%origin(3)
      print*, "------------------------------------------------------------"

      !PSB print*, ion(i)%death_loc(1), ion(i)%death_loc(2), ion(i)%death_loc(3), ion(i)%vel(1), ion(i)%vel(2), ion(i)%vel(3), &
      !PSB ion(i)%KE, ion(i)%qCEXion, ion(i)%source, ion(i)%xydisp, ion(i)%dir(1), ion(i)%dir(2), ion(i)%dir(3), &
      !PSB ion(i)%origin(1), ion(i)%origin(2), ion(i)%origin(3)
      !PSB i = npandgions
      !PSB print*, ion(i)%death_loc(1), ion(i)%death_loc(2), ion(i)%death_loc(3), ion(i)%vel(1), ion(i)%vel(2), ion(i)%vel(3), &
      !PSB ion(i)%KE, ion(i)%qCEXion, ion(i)%source, ion(i)%xydisp, ion(i)%dir(1), ion(i)%dir(2), ion(i)%dir(3), &
      !PSB ion(i)%origin(1), ion(i)%origin(2), ion(i)%origin(3)
      close(10)
    end subroutine read_ion_file

  !*******************************************************************************
    ! PSB ADDED ACCEL_THICKNESS
    subroutine mesh_flat_surface(domain, npt, vcl, vcl3d, maxntri, nfacept, accel_thickness)
    !-----------------------------------------------------------------------------
    !  Uses the geometry data in domain structure to generate points on downstream
    !  surface of the accel grid.  Assumes that the surface is flat (beginning of
    !  life geometry).
    !
    !  History:
    !
    !    Jay Polk
    !    May 9, 2021
    !    Nov 28, 2021--pulled final triangulation out into separate subroutine
    !
    !  Parameters:
    !
    !    Input, type (geom_type) domain, a data structure containing geometry info for the domain
    !    Output, real (kind = 8) npt, the number of points in the vcl
    !    Output, real (kind = 8) vcl, a (2,maxpt) array with the vertex coordinate list for the grid face mesh.
    !       This is the point array that is passed to dtris2 for delaunay triangulation.
    !    Output, real (kind = 8) vcl3d, a (3,maxpt) array with the vertex coordinate list for the grid face mesh
    !       (including the z-coordinate)
    !    Output, real (kind = 8) maxntri, the max number of triangular cells generated by dtris2
    !    Output, real (kind = 8) nfacept, the number of points in the vcl on the accel grid face
    !-----------------------------------------------------------------------------
      type (geom_type), intent(in) :: domain
      integer, intent(out) :: npt
      real*8, dimension(:,:), intent(out) :: vcl
      real*8, dimension(:,:), intent(out) :: vcl3d
      integer, intent(in) :: maxntri
      integer, intent(out) :: nfacept

      integer :: i, j, k
      integer :: x1numpnts
      integer :: x2numpnts = 8
      integer :: ierr
      integer :: numcircpts
      real*8, dimension(:), ALLOCATABLE :: x, y
      real*8 :: dx1, dx2
      real*8 :: rad, dtheta
      integer :: ind(size(vcl, dim = 2))
      integer :: numtri
      integer, dimension(3, maxntri) :: tilist   ! Local version of til for mesh refinement
                                                 ! using dtris2; not til for final mesh
      integer, dimension(3, size(til, dim = 2)) :: tnbr
      integer, dimension(size(vcl, dim = 2)) :: stack
      real*8, dimension(2) :: centroid


      ! PSB Added the following variables
      INTEGER, PARAMETER :: REMESH_ITERS = 3 
      INTEGER :: remeshIters 
      REAL*8, PARAMETER :: SMALL_TOL = 1.0E-5_8 ! Tolerance for floating point comparisons 
      REAL*8 :: x_val_1, x_val_2, x_val_3 
      REAL*8 :: y_val_1, y_val_2, y_val_3 
      REAL*8 :: z_val_1, z_val_2, z_val_3 
      REAL*8, DIMENSION(3) :: v12, v13, crossProd, normal_vector 
      REAL*8 :: bottom_val, xCentre, yCentre 
      LOGICAL :: isOnWall_major, isOnWall_minor, isOnWall_hypot, isInHole
      ! Dynamic arrays for elements and flagged_elements (Fortran equivalent of Matlab cell arrays/dynamic lists) 
      ! We'll use a dynamic array for triangles to filter them 
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: TRI_filtered ! To store valid triangles after filtering 
      INTEGER :: current_ntri_filtered ! Current number of triangles in TRI_filtered 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: flagged_elements ! To store indices of triangles to be removed 
      INTEGER :: num_flagged_elements
      INTEGER :: node_idx1, node_idx2, node_idx3
      LOGICAL :: is_flagged
      INTEGER :: row_count
      real*8 :: accel_thickness

      !---------------------------------------------------------------------------
      ! Define points along major and minor walls to build initial grid
      !---------------------------------------------------------------------------
      dx2 = (domain%xmax - domain%hole_diam/2.0d0)/(x2numpnts - 0.25d0)
      x1numpnts = NINT((domain%hole_diam/2.0d0)/dx2)

      allocate(x(x1numpnts + x2numpnts))
      allocate(y(x1numpnts + x2numpnts))

      dx1 = (domain%hole_diam/2.0D+00)/(x1numpnts - 1)
      do i = 1,x1numpnts
        x(i) = (i - 1) * dx1
        y(i) = x(i) * tan(30.0d0 * pi/180.0d0)
      end do

      do i = 1,x2numpnts
        x(i + x1numpnts) = (domain%hole_diam/2.0d0) + 0.75d0 * dx2 + (i - 1) * dx2
        y(i + x1numpnts) = x(i + x1numpnts) * tan(30 * pi/180)
      end do
      !---------------------------------------------------------------------------
      ! Define initial mesh on webbing
      !---------------------------------------------------------------------------
      npt = 0
      do i = 1, (x1numpnts + x2numpnts)
        do j = 1, i
          rad = sqrt(x(i)**2 + y(j)**2)
          if (rad > (1.04d0 * 0.5d0 * domain%hole_diam)) then
            npt = npt + 1
            vcl3d(1,npt) = x(i)
            vcl3d(2,npt) = y(j)
            vcl3d(3,npt) = domain%zmax
          end if
        end do
      end do
      !---------------------------------------------------------------------------
      ! Define mesh points on edge of hole
      !---------------------------------------------------------------------------
      numcircpts = domain%hole_diam * pi/12.0d0/(0.5d0 * dx2)
      dtheta = 30.0d0 * pi/180.0d0/numcircpts
      do i = 1, numcircpts + 1
        npt = npt + 1
        vcl3d(1,npt) = 0.5d0 * domain%hole_diam * cos((i - 1) * dtheta)
        vcl3d(2,npt) = 0.5d0 * domain%hole_diam * sin((i - 1) * dtheta)
        vcl3d(3,npt) = domain%zmax
      end do
      !---------------------------------------------------------------------------
      ! Add mesh points at midpoints of edges along major wall and hypotenuse
      !---------------------------------------------------------------------------
      do i = 1, x2numpnts-1
        npt = npt + 1
        vcl3d(1,npt) = x(i + x1numpnts) + dx2/2.0d0
        vcl3d(2,npt) = vcl3d(1,npt) * tan(30.0d0 * pi/180.0d0)
        vcl3d(3,npt) = domain%zmax
      end do

      do i = 1, x2numpnts-1
        npt = npt + 1
        vcl3d(1,npt) = x(i + x1numpnts) + dx2/2.0d0
        vcl3d(2,npt) = 0.0d0
        vcl3d(3,npt) = domain%zmax
      end do
      !---------------------------------------------------------------------------
      ! Create initial vertex coordinate list and index
      !---------------------------------------------------------------------------
  REMESH_LOOP: DO remeshIters = 1, REMESH_ITERS !PSB - Resmesh iterations
      do i = 1, npt
        vcl(1,i) = vcl3d(1,i)
        vcl(2,i) = vcl3d(2,i)
        ind(i) = i
      end do
      !---------------------------------------------------------------------------
      ! Generate initial Delaunay triangulation of mesh points
      !---------------------------------------------------------------------------
      call dtris2 (npt, maxnpt, vcl, ind, numtri, tilist, tnbr, stack, ierr)

      !---------------------------------------------------------------------------
      ! PSB CHANGE: Mesh cleaning, identifies and flags undesirable 
      ! triangles for removal to ensure the mesh accurately represents
      !---------------------------------------------------------------------------

      ! flagged_elements will temporarily store the row indices of triangles  that we identify as "bad" during the cleaning process.
      ALLOCATE(flagged_elements(numtri))
      num_flagged_elements = 0 ! Initialize counter for flagged triangles
      ! Calculate the Z-coordinate representing the "bottom" of the grid to identify "punch-through" triangles that have incorrectly extended below the intended lower surface of the grid.
      bottom_val = domain%zmax - accel_thickness 

      TRIANGLE_CLEANING_LOOP: DO row_count = 1, numtri
          ! Get node indices for the current triangle
          node_idx1 = tilist(1, row_count)
          node_idx2 = tilist(2, row_count)
          node_idx3 = tilist(3, row_count)
          ! Filter 1: Triangles with zero node indices
          ! Why: A node index of zero usually indicates an uninitialized or invalid entry.
          IF (node_idx1 == 0 .OR. node_idx2 == 0 .OR. node_idx3 == 0) THEN
                    num_flagged_elements = num_flagged_elements + 1
                    flagged_elements(num_flagged_elements) = row_count
                    print *, "A triangle had zero node indices"
                    CYCLE TRIANGLE_CLEANING_LOOP ! Skip to the next triangle
          END IF
          ! Filter 2: Triangles where all nodes are on the hypotenuse line
          ! Why: In this specific grid geometry, there's a defined "hypotenuse" boundary at 30 degrees.
          ! Get 3D coordinates of the current triangle's nodes      
          x_val_1 = vcl3d(1, node_idx1)
          y_val_1 = vcl3d(2, node_idx1)
          z_val_1 = vcl3d(3, node_idx1)
          x_val_2 = vcl3d(1, node_idx2)
          y_val_2 = vcl3d(2, node_idx2)
          z_val_2 = vcl3d(3, node_idx2)
          x_val_3 = vcl3d(1, node_idx3)
          y_val_3 = vcl3d(2, node_idx3)
          z_val_3 = vcl3d(3, node_idx3)
          isOnWall_hypot = (ABS(y_val_1 - (x_val_1 * TAN(30.0D0 * PI/180.0D0))) < SMALL_TOL) .AND. &
                            (ABS(y_val_2 - (x_val_2 * TAN(30.0D0 * PI/180.0D0))) < SMALL_TOL) .AND. &
                            (ABS(y_val_3 - (x_val_3 * TAN(30.0D0 * PI/180.0D0))) < SMALL_TOL)

          IF (isOnWall_hypot) THEN
              num_flagged_elements = num_flagged_elements + 1
              flagged_elements(num_flagged_elements) = row_count
              CYCLE TRIANGLE_CLEANING_LOOP
          END IF
      END DO TRIANGLE_CLEANING_LOOP

      ! Build the new_TRI (TRI_filtered) with only the good elements PSB
             IF (num_flagged_elements > 0) THEN
            ! Allocate the new filtered array with just enough space for the good triangles.
            ALLOCATE(TRI_filtered(3, numtri - num_flagged_elements))
            current_ntri_filtered = 0
            NEW_TRI_BUILD_LOOP: DO row_count = 1, numtri
                is_flagged = .FALSE.
                ! Check if the current triangle was flagged for removal.
                DO k = 1, num_flagged_elements
                    IF (row_count == flagged_elements(k)) THEN
                        is_flagged = .TRUE.
                        EXIT ! Found it, no need to check further
                    END IF
                END DO
                IF (.NOT. is_flagged) THEN
                    ! If not flagged, copy this triangle's node indices to the new filtered list.
                    current_ntri_filtered = current_ntri_filtered + 1
                    TRI_filtered(1, current_ntri_filtered) = tilist(1, row_count)
                    TRI_filtered(2, current_ntri_filtered) = tilist(2, row_count)
                    TRI_filtered(3, current_ntri_filtered) = tilist(3, row_count)
                END IF
            END DO NEW_TRI_BUILD_LOOP
            DEALLOCATE(TRI_filtered) ! PSB 2
          END IF
          DEALLOCATE(flagged_elements)

      !---------------------------------------------------------------------------
      ! Refine mesh around pits and grooves region (add points at centroids of
      ! triangular elements from initial triangulation)
      !---------------------------------------------------------------------------
      !Refine mesh after we know the current mesh is clean
      MESHREFINE: do i = 1,numtri
        ! CHECK INVALID NODES
          node_idx1 = tilist(1, row_count)
          node_idx2 = tilist(2, row_count)
          node_idx3 = tilist(3, row_count)
          
          IF (node_idx1<1 .OR. &
              node_idx2<1 .OR. &
              node_idx3<1 ) THEN
          print *, "Invalid node index in triangle" ,i,  "with indices " , node_idx1, node_idx2, node_idx3
          CYCLE MESHREFINE
          END IF

        centroid(1) = (vcl(1,tilist(1,i)) + vcl(1,tilist(2,i)) + vcl(1,tilist(3,i)))/3.0d0
        centroid(2) = (vcl(2,tilist(1,i)) + vcl(2,tilist(2,i)) + vcl(2,tilist(3,i)))/3.0d0
        ! PSB added to check if the new point will exceed max array size

        IF (npt < SIZE(vcl, DIM = 2)) then
          if (sqrt(centroid(1)**2 + centroid(2)**2) > 0.8182d0 * domain%xmax) then
              npt = npt + 1
              vcl3d(1,npt) = centroid(1)
              vcl3d(2,npt) = centroid(2)
              vcl3d(3,npt) = domain%zmax
              vcl(1,npt) = vcl3d(1,npt)
              vcl(2,npt) = vcl3d(2,npt)
              ind(npt) = npt
          end if
        ELSE  !PSB
          print *, "Warning: Max number of points reached during refinement. Cannot add more points." !PSB
          EXIT REMESH_LOOP
        END IF
      end do MESHREFINE
    end do REMESH_LOOP! End of remesh iters 
    nfacept = npt     ! Define number of points in accel face mesh
    
    end subroutine mesh_flat_surface

  !*******************************************************************************
    subroutine triangulate_surface(domain, npt, vcl, vcl3d, ntri, til, nfacetri, cell)
    !-----------------------------------------------------------------------------
    !  Uses the geometry data in domain structure to generate a mesh on downstream
    !  surface of the accel grid.  Assumes that the surface is flat (beginning of
    !  life geometry).
    !
    !  History:
    !
    !    Jay Polk
    !    May 9, 2021
    !
    !  Parameters:
    !
    !    Input, type (geom_type) domain, a data structure containing geometry info for the domain
    !    Input, real (kind = 8) npt, the number of points in the vcl
    !    Input, real (kind = 8) vcl, a (2,maxpt) array with the vertex coordinate list for the grid face mesh.
    !       This is the point array that is passed to dtris2 for delaunay triangulation.
    !    Input, real (kind = 8) vcl3d, a (3,maxpt) array with the vertex coordinate list for the grid face mesh
    !       (including the z-coordinate)
    !    Output, real (kind = 8) ntri, the number of triangular cells generated by dtris2
    !    Output, integer til, the triangle incidence list generated by dtris2; elements
    !       are indices of vcl; vertices of triangles are in counterclockwise order.
    !    Output, real (kind = 8) nfacetri, the number of triangular cells on the accel grid face
    !    Input/Output, type (cell_type) cell, data structure containing cell info; this subroutine only populates
    !       the on_bndry component = 0 (all cells generated here are on accel grid face)
    !-----------------------------------------------------------------------------
      type (geom_type), intent(in) :: domain
      integer, intent(in) :: npt
      real*8, dimension(:,:), intent(in) :: vcl
      real*8, dimension(:,:), intent(in) :: vcl3d
      integer, intent(inout) :: ntri
      integer, dimension(:,:), intent(inout) :: til
      integer, intent(out) :: nfacetri
      type (cell_type), dimension(:), intent(inout) :: cell

      integer :: i, k
      integer :: new_ntri
      integer :: ierr
      integer :: ind(size(vcl, dim = 2))
      integer, dimension(3, size(til, dim = 2)) :: tnbr, new_til
      integer, dimension(size(vcl, dim = 2)) :: stack
      logical :: isOnWall_major, isOnWall_minor, isOnWall_hypot, isInHole

      print *, "Number of Points (npt):", npt !PSB
      print *, "Coordinates of Point 226 (vcl3d(:,226)):"  !PSB
      print *, "  x:", vcl3d(1,226) !PSB
      print *, "  y:", vcl3d(2,226) !PSB
      print *, "  z:", vcl3d(3,226) !PSB

      !print*, "triangulation routine: npt, vcl3d(226)", npt,vcl3d(:,226) PSB
      !---------------------------------------------------------------------------
      ! Generate the final Delaunay triangulation
      !---------------------------------------------------------------------------
      do i = 1, npt
        ind(i) = i
      end do
      call dtris2 (npt, npt, vcl, ind, ntri, til, tnbr, stack, ierr)
      !---------------------------------------------------------------------------
      ! Clean up the mesh (remove bad triangles in hole or coincident with boundaries)
      !---------------------------------------------------------------------------
      new_ntri = 0
      do i = 1,ntri
        cell(i)%centroid(1) = (vcl(1,til(1,i)) + vcl(1,til(2,i)) + vcl(1,til(3,i)))/3.0d0
        cell(i)%centroid(2) = (vcl(2,til(1,i)) + vcl(2,til(2,i)) + vcl(2,til(3,i)))/3.0d0

        isOnWall_major = (abs(vcl(1,til(1,i))) <= 1d-4) .and. (abs(vcl(1,til(2,i))) <= 1d-4) &
        .and. (abs(vcl(1,til(3,i))) <= 1d-4)
        isOnWall_minor = (abs(vcl(1,til(1,i)) - domain%xmax) <= 1e-4) .and. (abs(vcl(1,til(2,i)) - domain%xmax)  &
        < 1d-4) .and. (abs(vcl(1,til(3,i)) - domain%xmax) <= 1d-4)
        isOnWall_hypot = (abs(vcl(1,til(1,i)) * tan(30.0d0) - vcl(2,til(1,i))) <= 1d-4) &
        .and. (abs(vcl(1,til(2,i)) * tan(30.0d0)-vcl(2,til(2,i))) <= 1d-4) &
        .and. (abs(vcl(1,til(3,i)) * tan(30.0d0)-vcl(2,til(3,i))) <= 1d-4)
        isInHole = sqrt(cell(i)%centroid(1)**2 + cell(i)%centroid(2)**2) <= (0.5d0 * domain%hole_diam)

        if (isOnWall_major .or. isOnWall_minor .or. isOnWall_hypot .or. isInHole) then
          ! Don't include in new version of til(:, :)
        else
          new_ntri = new_ntri + 1
          new_til(1,new_ntri) = til(1,i)
          new_til(2,new_ntri) = til(2,i)
          new_til(3,new_ntri) = til(3,i)
          cell(new_ntri)%on_bndry = 0         !  cells are located on grid face, not boundaries
        end if
      end do

      do i = 1, new_ntri
        til(:,i) = new_til(:,i)
      end do
      ntri = new_ntri
      nfacetri = ntri       ! Define number of triangles in accel face mesh
      ! print*, "triangulation routine: ntri, til(100)", ntri,til(:,100) !PSB
      print *, "Number of Triangles (ntri):", ntri !PSB
      print *, "Vertices of Triangle 100 (til(:,100)):" !PSB
      print *, "  Vertex 1 index:", til(1,100) !PSB
      print *, "  Vertex 2 index:", til(2,100) !PSB
      print *, "  Vertex 3 index:", til(3,100) !PSB

    end subroutine triangulate_surface

  !*******************************************************************************
    subroutine mesh_boundaries(domain, npt, vcl, vcl3d, ntri, til, cell, nfacept)
    !-----------------------------------------------------------------------------
    !  Uses the geometry data in domain structure to generate a mesh on the top
    !  and the symmetry boundaries.  The triangle indices are defined in counter-
    !  clockwise order looking from inside the domain, so that cell normal vectors
    !  point inward.
    !
    !  History:
    !
    !    Jay Polk
    !    May 9, 2021
    !
    !  Parameters:
    !
    !    Input, type (geom_type) domain, a data structure containing geometry info for the domain
    !    Input/Output, integer npt, the number of points in the vcl. Modified to include boundary points.
    !    Input/Output, real (kind = 8) vcl, a (2,maxpt) array with the vertex coordinate list for the grid face mesh.
    !       Points on the domain boundaries are added to the vcl for the accel grid face.
    !    Input/Output, real (kind = 8) vcl3d, a (3,maxpt) array with the vertex coordinate list for the grid face mesh
    !       (including the z-coodinate)
    !    Input/Output, integer ntri, the number of triangular cells. Modified to include boundary cells.
    !    Input/Output, integer til, the triangle incidence list; elements are indices of vcl;
    !       vertices of triangles are in counterclockwise order.
    !    Input/Output, type (cell_type) cell, data structure containing cell info; this subroutine only populates
    !       the on_bndry component
    !-----------------------------------------------------------------------------
      type (geom_type), intent(in) :: domain
      integer, intent(inout) :: npt
      real*8, dimension(:,:), intent(inout) :: vcl
      real*8, dimension(:,:), intent(inout) :: vcl3d
      integer, intent(inout) :: ntri
      integer, dimension(:,:), intent(inout) :: til
      type (cell_type), dimension(:), intent(inout) :: cell
      integer, intent(inout) :: nfacept
      integer :: i, k
      integer, dimension(size(vcl, dim = 2)) :: ind
      !---------------------------------------------------------------------------
      ! Add vertices and triangles for side walls and top
      !---------------------------------------------------------------------------
      ! Define vertices on front major wall (tri and rect domains)
      npt = npt + 1                 ! Point 1
      vcl3d(1,npt) = domain%xmin
      vcl3d(2,npt) = domain%ymin
      vcl3d(3,npt) = domain%zmin
      npt = npt + 1                 ! Point 2
      vcl3d(1,npt) = domain%xmax
      vcl3d(2,npt) = domain%ymin
      vcl3d(3,npt) = domain%zmin
      npt = npt + 1                 ! Point 3
      vcl3d(1,npt) = domain%xmax
      vcl3d(2,npt) = domain%ymin
      vcl3d(3,npt) = domain%zmax
      npt = npt + 1                 ! Point 4
      vcl3d(1,npt) = domain%xmin
      vcl3d(2,npt) = domain%ymin
      vcl3d(3,npt) = domain%zmax
      ! Define vertices on corner of hypotenuse wall or back major wall and right minor wall (tri and rect domains)
      npt = npt + 1                 ! Point 5
      vcl3d(1,npt) = domain%xmax
      vcl3d(2,npt) = domain%ymax
      vcl3d(3,npt) = domain%zmin
      npt = npt + 1                 ! Point 6
      vcl3d(1,npt) = domain%xmax
      vcl3d(2,npt) = domain%ymax
      vcl3d(3,npt) = domain%zmax
      ! Define vertices on corner of back major wall and left minor wall (rect domains)
      if (domain%type == 'rect') then
        npt = npt + 1                 ! Point 7
        vcl3d(1,npt) = domain%xmin
        vcl3d(2,npt) = domain%ymax
        vcl3d(3,npt) = domain%zmax
        npt = npt + 1                 ! Point 8
        vcl3d(1,npt) = domain%xmin
        vcl3d(2,npt) = domain%ymax
        vcl3d(3,npt) = domain%zmin
      end if
      do i = nfacept + 1, npt
        vcl(1,i) = vcl3d(1,i)
        vcl(2,i) = vcl3d(2,i)
        ind(i) = i
      end do

      ! Define triangles on front major wall (tri and rect domains)
      ntri = ntri + 1                 ! Triangle 1 on Jay's GoodNotes diagram
      til(1,ntri) = nfacept + 1
      til(2,ntri) = nfacept + 4
      til(3,ntri) = nfacept + 2
      cell(ntri)%on_bndry = 1
      ntri = ntri + 1                 ! Triangle 2 on Jay's GoodNotes diagram
      til(1,ntri) = nfacept + 2
      til(2,ntri) = nfacept + 4
      til(3,ntri) = nfacept + 3
      cell(ntri)%on_bndry = 1
      ! Define triangles on right minor wall (tri and rect domains)
      ntri = ntri + 1                 ! Triangle 3 on Jay's GoodNotes diagram
      til(1,ntri) = nfacept + 2
      til(2,ntri) = nfacept + 6
      til(3,ntri) = nfacept + 5
      cell(ntri)%on_bndry = 2
      ntri = ntri + 1                 ! Triangle 4 on Jay's GoodNotes diagram
      til(1,ntri) = nfacept + 2
      til(2,ntri) = nfacept + 3
      til(3,ntri) = nfacept + 6
      cell(ntri)%on_bndry = 2
      if (domain%type == 'tri') then
        ! Define triangles on hypotenuse wall (tri domains)
        ntri = ntri + 1                 ! Triangle 5 on Jay's GoodNotes diagram
        til(1,ntri) = nfacept + 4
        til(2,ntri) = nfacept + 1
        til(3,ntri) = nfacept + 5
        cell(ntri)%on_bndry = 3
        ntri = ntri + 1                 ! Triangle 6 on Jay's GoodNotes diagram
        til(1,ntri) = nfacept + 4
        til(2,ntri) = nfacept + 5
        til(3,ntri) = nfacept + 6
        cell(ntri)%on_bndry = 3
        ! Define triangle on top wall (tri domains)
        ntri = ntri + 1                 ! Triangle 7 on Jay's GoodNotes diagram
        til(1,ntri) = nfacept + 4
        til(2,ntri) = nfacept + 6
        til(3,ntri) = nfacept + 3
        cell(ntri)%on_bndry = -1
        ! Define triangles on bottom of domain; used in test for punch-through (tri and rect domains)
        ntri = ntri + 1                 ! Triangle -1 on Jay's GoodNotes diagram
        til(1,ntri) = nfacept + 1
        til(2,ntri) = nfacept + 2
        til(3,ntri) = nfacept + 5
        cell(ntri)%on_bndry = -2
      else      ! rect domain
        ! Define triangles on back major wall (rect domains)
        ntri = ntri + 1                 ! Triangle 5 on Jay's GoodNotes diagram
        til(1,ntri) = nfacept + 7
        til(2,ntri) = nfacept + 8
        til(3,ntri) = nfacept + 5
        cell(ntri)%on_bndry = 3
        ntri = ntri + 1                 ! Triangle 6 on Jay's GoodNotes diagram
        til(1,ntri) = nfacept + 7
        til(2,ntri) = nfacept + 5
        til(3,ntri) = nfacept + 6
        cell(ntri)%on_bndry = 3
        ! Define triangles on other minor wall (rect domains)
        ntri = ntri + 1                 ! Triangle 7 on Jay's GoodNotes diagram
        til(1,ntri) = nfacept + 1
        til(2,ntri) = nfacept + 7
        til(3,ntri) = nfacept + 4
        cell(ntri)%on_bndry = 4
        ntri = ntri + 1                 ! Triangle 8 on Jay's GoodNotes diagram
        til(1,ntri) = nfacept + 1
        til(2,ntri) = nfacept + 8
        til(3,ntri) = nfacept + 7
        cell(ntri)%on_bndry = 4
        ! Define triangles on top face (rect domains)
        ntri = ntri + 1                 ! Triangle 9 on Jay's GoodNotes diagram
        til(1,ntri) = nfacept + 4
        til(2,ntri) = nfacept + 6
        til(3,ntri) = nfacept + 3
        cell(ntri)%on_bndry = -1
        ntri = ntri + 1                 ! Triangle 10 on Jay's GoodNotes diagram
        til(1,ntri) = nfacept + 4
        til(2,ntri) = nfacept + 7
        til(3,ntri) = nfacept + 6
        cell(ntri)%on_bndry = -1
        ! Define triangles on bottom of domain; used in test for punch-through (tri and rect domains)
        ntri = ntri + 1                 ! Triangle -1 on Jay's GoodNotes diagram
        til(1,ntri) = nfacept + 1
        til(2,ntri) = nfacept + 2
        til(3,ntri) = nfacept + 5
        cell(ntri)%on_bndry = -2
        ntri = ntri + 1                 ! Triangle -2 on Jay's GoodNotes diagram
        til(1,ntri) = nfacept + 1
        til(2,ntri) = nfacept + 5
        til(3,ntri) = nfacept + 6
        cell(ntri)%on_bndry = -2
      end if
    end subroutine mesh_boundaries
  !*******************************************************************************
    subroutine calc_cell_params(acell, cellnum)
    !-----------------------------------------------------------------------------
    !  Calculates cell parameters to populate parts of the cell data structure
    !
    !  History:
    !
    !    Jay Polk
    !    Dec 10, 2021
    !
    !    This version operates on individual cells--loop over all cells is now in
    !    main program, making this subroutine easier to test.
    !
    !  Parameters:
    !
    !    Input/output, type (cell_type) acell, data structure for a single cell;
    !       this subroutine populates many of the components
    !    Input integer numcell, the array index of the individual cell passed in
    !-----------------------------------------------------------------------------
      type (cell_type), intent(inout) :: acell
      integer, intent(in) :: cellnum

      real*8, dimension(3) :: normalvec
      real*8 :: norm
      real*8 :: sum
      !---------------------------------------------------------------------------
      ! Calculate cell geometric parameters (populate cell data structure)
      !---------------------------------------------------------------------------

      ! Find vectors for two edges sharing vert1 and for third edge
      acell%edge1(:) = acell%vert2(:) - acell%vert1(:)
      acell%edge2(:) = acell%vert3(:) - acell%vert1(:)
      acell%edge3(:) = acell%vert3(:) - acell%vert2(:)

      ! Calculate edge lengths
      acell%edgelength(1) = norm2(acell%edge1(:))
      acell%edgelength(2) = norm2(acell%edge2(:))
      acell%edgelength(3) = norm2(acell%edge3(:))

      ! Calculate angles associated with vertices
      acell%vertangle(1) = angle(acell%edge1(:), acell%edge2(:))
      acell%vertangle(2) = angle(acell%edge3(:), -1.0d0 * acell%edge1(:))
      acell%vertangle(3) = angle(-1.0d0 * acell%edge2(:), -1.0d0 * acell%edge3(:))
      sum = acell%vertangle(1) + acell%vertangle(2) + acell%vertangle(3)
      if (abs(sum - pi) > 1e-8) then
        print*, "Error in cell angle calculation--3 angles do not sum to pi"
        print*, "cellnum, angles, sum of angles = ", cellnum, acell%vertangle(1), acell%vertangle(2), acell%vertangle(3), sum
      end if

      !  Define cell-centered normal and area
      normalvec = cross(acell%edge1(:),acell%edge2(:))
      norm = norm2(normalvec(:))
      acell%area = norm/2
      acell%normal(:) = normalvec(:)/norm

      !  Define cell centroid coords
      acell%centroid(1) = (vcl3d(1,acell%til(1)) + vcl3d(1,acell%til(2)) + vcl3d(1,acell%til(3)))/3.0d0
      acell%centroid(2) = (vcl3d(2,acell%til(1)) + vcl3d(2,acell%til(2)) + vcl3d(2,acell%til(3)))/3.0d0
      acell%centroid(3) = (vcl3d(3,acell%til(1)) + vcl3d(3,acell%til(2)) + vcl3d(3,acell%til(3)))/3.0d0
    end subroutine calc_cell_params

  !*******************************************************************************
    subroutine calc_node_params(node, nodenum, domain, nfacetri, cell)
    !-----------------------------------------------------------------------------
    !  Calculates node parameters to populate parts of the node data structure
    !
    !  History:
    !
    !    Jay Polk
    !    May 10, 2021
    !
    !  Parameters:
    !
    !    Input/Output, type (node_type) node, data structure containing node info;
    !       this subroutine populates many of the components
    !    Input, integer nodenum, the array index of the node passed in
    !    Input, type (geom_type) domain, a data structure containing geometry info for the domain
    !    Input, integer nfacetri, the number of triangles (cells) on the grid face
    !    Input, type (cell_type) cell, data structure array containing cell info
    !-----------------------------------------------------------------------------
      type (node_type), intent(inout) :: node
      integer, intent(in) :: nodenum
      type (geom_type), intent(in) :: domain
      integer, intent(in) :: nfacetri
      type (cell_type), dimension(:), intent(in) :: cell

      real*8, dimension(3) :: normalvec
      real*8 :: normal(3)
      integer :: i, j, k, member

      !---------------------------------------------------------------------------
      ! Determine if node is on edges or corners of domain
      !---------------------------------------------------------------------------
      node%on_bndry = 0      ! initialize to zero (assume internal node)
      node%on_corner = 0
      if (abs(node%coords(2) - domain%ymin) <= 1d-4) then
        ! On lower major edge
        if (abs(node%coords(1) - domain%xmin) <= 1d-4) then
          ! On lower left corner
          node%on_corner = 1
        else if (abs(node%coords(1) - domain%xmax) <= 1d-4) then
          ! On lower right corner
          node%on_corner = 2
        else
          ! Not on a corner
          node%on_bndry = 2
        end if
      else if (abs(node%coords(1) - domain%xmax) <= 1d-4) then
        ! On right minor edge
        if (abs(node%coords(2) - domain%ymax) <= 1d-4) then
          ! On upper right corner
          node%on_corner = 3
        else
          ! Not on a corner
          node%on_bndry = 3
        end if
      else if ((abs(node%coords(1) * tan(30.0d0) - node%coords(2)) <= 1d-4) .and. (domain%type == 'tri')) then
        ! On hypotenuse
        if (abs((node%coords(1)**2 + node%coords(2)**2) - (domain%hole_diam/2)**2) <= 1d-4) then
          ! Corner of hypotenuse and aperture wall
          node%on_corner = 4
        else
          ! Not on corner
          node%on_bndry = 4
        end if
      else if ((abs(node%coords(2) - domain%ymax) <= 1d-4) .and. (domain%type == 'rect')) then
        ! On upper major edge (rect)
        if (abs(node%coords(1) - domain%xmin) <= 1d-4) then
          ! On upper left corner (rect)
          node%on_corner = 4
        else
          ! Not on a corner
          node%on_bndry = 4
        end if
      else if ((abs(node%coords(1) - domain%xmin) <= 1d-4) .and. (domain%type == 'rect')) then
        ! On left minor edge (rect)
        node%on_bndry = 1
      else if (abs((node%coords(1)**2 + node%coords(2)**2) - (domain%hole_diam/2)**2) <= 1d-4) then
        ! On accel aperture edge
        node%on_bndry = 1
      else
        ! Internal node; already set on_bndry and on_corner to zero initially
      end if

      !---------------------------------------------------------------------------
      ! Identify which cells contain node and store the associated angles and areas
      !---------------------------------------------------------------------------
      node%in_cells = 0          ! Set in_cells to all zeroes initially
      node%area = 0              ! Set area sum to zero initially
      member = 1
      do j = 1, nfacetri
        do k = 1, 3
          if (cell(j)%til(k) == i) then ! Triangle incidence list for cell(j) contains node i at k
            node%in_cells(member) = j                        ! Add cell j to in_cells list
            node%cell_angle(member) = cell(j)%vertangle(k)   ! Add angle to cell_angle list
            node%area = node%area + 1.0d0/3.0d0 * cell(j)%area
            member = member + 1
            exit
          end if
        end do
      end do
      !---------------------------------------------------------------------------
      ! Calculate the node normal vector as the angle-weighted average of surrounding cells
      !---------------------------------------------------------------------------
      if (node%on_corner > 0) then
        node%normal(:) = (/ 0.0d0, 0.0d0, 1.0d0 /)
      else
        member = 1
        normal(:) = 0.0d0
        do while (node%in_cells(member) /= 0)
          normal(:) = normal(:) + node%cell_angle(member) * cell(node%in_cells(member))%normal(:)
          member = member + 1
        end do
        node%normal(:) = normal(:)/norm2(normal(:))
      end if
      !---------------------------------------------------------------------------
      ! Initialize variables used in erosion/geom update calculations
      !---------------------------------------------------------------------------
      ! node%vel = (/ 0.0d0, 0.0d0, 0.0d0 /)
      ! node%vel_old = (/ 0.0d0, 0.0d0, 0.0d0 /)
      node%mass_loss = 0.0d0
      ! node%firstStep = 1
    end subroutine calc_node_params

  !*******************************************************************************
    subroutine calc_cell_params_old(ntri, til, vcl3d, cell)
    !-----------------------------------------------------------------------------
    !  Calculates cell parameters to populate parts of the cell data structure
    !
    !  History:
    !
    !    Jay Polk
    !    May 9, 2021
    !
    !  Parameters:
    !
    !    Input, integer ntri, the number of triangular cells.
    !    Input, integer til, the triangle incidence list; elements are indices of vcl;
    !       vertices of triangles are in counterclockwise order.
    !    Input, real (kind = 8) vcl3d, a (3,maxpt) array with the vertex coordinate list for the grid face mesh
    !       (including the z-coodinate)
    !    Output, type (cell_type) cell, data structure containing cell info; this subroutine populates many of
    !       the components
    !-----------------------------------------------------------------------------
      integer, intent(in) :: ntri
      integer, dimension(:,:), intent(in) :: til
      real*8, dimension(:,:), intent(in) :: vcl3d
      type (cell_type), dimension(:), intent(out) :: cell
      real*8, dimension(3) :: normalvec
      real*8 :: norm
      real*8 :: sum
      integer :: i
      !---------------------------------------------------------------------------
      ! Calculate cell geometric parameters (populate cell data structure)
      !---------------------------------------------------------------------------
      do i = 1,ntri
        ! Define the triangle incidence list for each cell
        cell(i)%til(:) = til(:,i)

        ! Define vertices of each cell
        cell(i)%vert1(:) = vcl3d(:,cell(i)%til(1))
        cell(i)%vert2(:) = vcl3d(:,cell(i)%til(2))
        cell(i)%vert3(:) = vcl3d(:,cell(i)%til(3))

        ! Find vectors for two edges sharing vert1 and for third edge
        cell(i)%edge1(:) = vcl3d(:,cell(i)%til(2)) - vcl3d(:,cell(i)%til(1))
        cell(i)%edge2(:) = vcl3d(:,cell(i)%til(3)) - vcl3d(:,cell(i)%til(1))
        cell(i)%edge3(:) = vcl3d(:,cell(i)%til(3)) - vcl3d(:,cell(i)%til(2))

        ! Calculate edge lengths
        cell(i)%edgelength(1) = norm2(cell(i)%edge1(:))
        cell(i)%edgelength(2) = norm2(cell(i)%edge2(:))
        cell(i)%edgelength(3) = norm2(cell(i)%edge3(:))

        ! Calculate angles associated with vertices
        cell(i)%vertangle(1) = angle(cell(i)%edge1(:), cell(i)%edge2(:))
        cell(i)%vertangle(2) = angle(cell(i)%edge3(:), -1.0d0 * cell(i)%edge1(:))
        cell(i)%vertangle(3) = angle(-1.0d0 * cell(i)%edge2(:), -1.0d0 * cell(i)%edge3(:))
        sum = cell(i)%vertangle(1) + cell(i)%vertangle(2) + cell(i)%vertangle(3)
        if (abs(sum - pi) > 1e-8) then
          print*, "Error in cell angle calculation--3 angles do not sum to pi"
          print*, "i, cell, sum of angles = ", i, cell(i)%vertangle(1), cell(i)%vertangle(2), cell(i)%vertangle(3), sum
        end if

        !  Define cell-centered normal and area
        normalvec = cross(cell(i)%edge1(:),cell(i)%edge2(:))
        norm = norm2(normalvec(:))
        cell(i)%area = norm/2
        cell(i)%normal(:) = normalvec(:)/norm

        !  Define cell centroid coords
        cell(i)%centroid(1) = (vcl3d(1,cell(i)%til(1)) + vcl3d(1,cell(i)%til(2)) + vcl3d(1,cell(i)%til(3)))/3.0d0
        cell(i)%centroid(2) = (vcl3d(2,cell(i)%til(1)) + vcl3d(2,cell(i)%til(2)) + vcl3d(2,cell(i)%til(3)))/3.0d0
        cell(i)%centroid(3) = (vcl3d(3,cell(i)%til(1)) + vcl3d(3,cell(i)%til(2)) + vcl3d(3,cell(i)%til(3)))/3.0d0
      end do

      open(16, FILE = 'NormVecs')
      do i = 1,ntri
        WRITE(16,*) cell(i)%normal(1), cell(i)%normal(2), cell(i)%normal(3)
      end do
      close(16)
    end subroutine calc_cell_params_old

  !*******************************************************************************
    subroutine calc_node_params_old(nfacept, vcl3d, domain, nfacetri, cell, node)
    !-----------------------------------------------------------------------------
    !  Calculates node parameters to populate parts of the node data structure
    !
    !  History:
    !
    !    Jay Polk
    !    May 10, 2021
    !
    !  Parameters:
    !
    !    Input, integer nfacept, the number of points (nodes) on the grid face
    !    Input, real (kind = 8) vcl3d, a (3,maxpt) array with the vertex coordinate list for the grid face mesh
    !       (including the z-coodinate)
    !    Input, type (geom_type) domain, a data structure containing geometry info for the domain
    !    Input, integer nfacetri, the number of triangles (cells) on the grid face
    !    Input, type (cell_type) cell, data structure containing cell info
    !    Input/Output, type (node_type) node, data structure containing node info;
    !       this subroutine populates many of the components
    !-----------------------------------------------------------------------------
      integer, intent(in) :: nfacept
      real*8, dimension(:,:), intent(in) :: vcl3d
      type (geom_type), intent(in) :: domain
      integer, intent(in) :: nfacetri
      type (cell_type), dimension(:), intent(in) :: cell
      type (node_type), dimension(:), intent(inout) :: node
      real*8, dimension(3) :: normalvec
      real*8 :: normal(3)
      integer :: i, j, k, member

      do i = 1,nfacept
        !---------------------------------------------------------------------------
        ! Define coordinates of node from vertex coordinate list
        !---------------------------------------------------------------------------
        node(i)%coords(:) = vcl3d(:,i)
        !---------------------------------------------------------------------------
        ! Determine if node is on edges or corners of domain
        !---------------------------------------------------------------------------
        node(i)%on_bndry = 0      ! initialize to zero (assume internal node)
        node(i)%on_corner = 0
        if (abs(node(i)%coords(2) - domain%ymin) <= 1d-4) then
          ! On lower major edge
          if (abs(node(i)%coords(1) - domain%xmin) <= 1d-4) then
            ! On lower left corner
            node(i)%on_corner = 1
          else if (abs(node(i)%coords(1) - domain%xmax) <= 1d-4) then
            ! On lower right corner
            node(i)%on_corner = 2
          else
            ! Not on a corner
            node(i)%on_bndry = 2
          end if
        else if (abs(node(i)%coords(1) - domain%xmax) <= 1d-4) then
          ! On right minor edge
          if (abs(node(i)%coords(2) - domain%ymax) <= 1d-4) then
            ! On upper right corner
            node(i)%on_corner = 3
          else
            ! Not on a corner
            node(i)%on_bndry = 3
          end if
        else if ((abs(node(i)%coords(1) * tan(30.0d0) - node(i)%coords(2)) <= 1d-4) .and. (domain%type == 'tri')) then
          ! On hypotenuse
          if (abs((node(i)%coords(1)**2 + node(i)%coords(2)**2) - (domain%hole_diam/2)**2) <= 1d-4) then
            ! Corner of hypotenuse and aperture wall
            node(i)%on_corner = 4
          else
            ! Not on corner
            node(i)%on_bndry = 4
          end if
        else if ((abs(node(i)%coords(2) - domain%ymax) <= 1d-4) .and. (domain%type == 'rect')) then
          ! On upper major edge (rect)
          if (abs(node(i)%coords(1) - domain%xmin) <= 1d-4) then
            ! On upper left corner (rect)
            node(i)%on_corner = 4
          else
            ! Not on a corner
            node(i)%on_bndry = 4
          end if
        else if ((abs(node(i)%coords(1) - domain%xmin) <= 1d-4) .and. (domain%type == 'rect')) then
          ! On left minor edge (rect)
          node(i)%on_bndry = 1
        else if (abs((node(i)%coords(1)**2 + node(i)%coords(2)**2) - (domain%hole_diam/2)**2) <= 1d-4) then
          ! On accel aperture edge
          node(i)%on_bndry = 1
        else
          ! Internal node; already set on_bndry and on_corner to zero initially
        end if
      end do

      !---------------------------------------------------------------------------
      ! Identify which cells contain node and store the associated angles and areas
      !---------------------------------------------------------------------------
      do i = 1, nfacept
        node(i)%in_cells = 0          ! Set in_cells to all zeroes initially
        node(i)%area = 0              ! Set area sum to zero initially
        member = 1
        do j = 1, nfacetri
          do k = 1, 3
            if (cell(j)%til(k) == i) then ! Triangle incidence list for cell(j) contains node i at k
              node(i)%in_cells(member) = j                        ! Add cell j to in_cells list
              node(i)%cell_angle(member) = cell(j)%vertangle(k)   ! Add angle to cell_angle list
              node(i)%area = node(i)%area + 1.0d0/3.0d0 * cell(j)%area
              member = member + 1
              exit
            end if
          end do
        end do
        !---------------------------------------------------------------------------
        ! Correct area if node is on boundary edge or corner (symmetry surfaces)
        !---------------------------------------------------------------------------
        select case (node(i)%on_corner)
          case (1)
            if (domain%type == 'rect') then
              ! lower left corner
              ! node(i)%area = 4.0d0 * node(i)%area
            else
              ! lower left corner of aperture wall edge (tri)
              ! node(i)%area = 2.0d0 * node(i)%area
            end if
          case (2:3)   ! node is on lower right or upper right corner
            ! node(i)%area = 4.0d0 * node(i)%area
          case (4)
            if (domain%type == 'rect') then
              ! upper right corner
              ! node(i)%area = 4.0d0 * node(i)%area
            else
              ! upper right corner of aperture wall edge (tri)
              ! node(i)%area = 2.0d0 * node(i)%area
            end if
          case default
            ! internal node; area is correct
          end select
          select case (node(i)%on_bndry)
            case (1)
              if (domain%type == 'rect') then
                ! left boundary of rect domain
                ! node(i)%area = 2.0d0 * node(i)%area
              end if
            case (2:4)   ! node is on lower, right, or upper (hypotenuse) boundary
              ! node(i)%area = 2.0d0 * node(i)%area
            case default
              ! internal node; area is correct
          end select
      end do
      do i = 1, nfacept    !Loop for normal vector calc and initializing geom update params
        !---------------------------------------------------------------------------
        ! Calculate the node normal vector as the angle-weighted average of surrounding cells
        !---------------------------------------------------------------------------
        if (node(i)%on_corner > 0) then
          node(i)%normal(:) = (/ 0.0d0, 0.0d0, 1.0d0 /)
          cycle
        else
          member = 1
          normal(:) = 0.0d0
          do while (node(i)%in_cells(member) /= 0)
            normal(:) = normal(:) + node(i)%cell_angle(member) * cell(node(i)%in_cells(member))%normal(:)
            member = member + 1
          end do
          node(i)%normal(:) = normal(:)/norm2(normal(:))
        end if
        !---------------------------------------------------------------------------
        ! Initialize variables used in erosion/geom update calculations
        !---------------------------------------------------------------------------
        ! node(i)%vel = (/ 0.0d0, 0.0d0, 0.0d0 /)
        ! node(i)%vel_old = (/ 0.0d0, 0.0d0, 0.0d0 /)
        node(i)%mass_loss = 0.0d0
        ! node(i)%firstStep = 1
      end do
    end subroutine calc_node_params_old

  !*******************************************************************************
    subroutine ABintegrate(coords, vel, vel_old, stepSize, firstStep)
    !-----------------------------------------------------------------------------
    !  Adams-Bashforth integrator: updates the position of a point based on calculated
    !  velocity over a given time step.
    !
    !  History:
    !
    !    Jay Polk
    !    November 24, 2021
    !
    !  Parameters:
    !
    !    Input/Output, real (kind = 8) coords, a dim (3) array with the current coordinates
    !    Input/Output, real (kind = 8) vel, a dim (3) array with the current velocity
    !    Input/Output, real (kind = 8) vel_old, a dim (3) array with the velocity from the previous step
    !    Input, real (kind = 8) stepSize, integration time step
    !    Input/Output, integer firstStep, flag indicating if this is the first step of the integration
    !     (algorithm uses Euler integration for first step because vel_old is not available)
    !-----------------------------------------------------------------------------
      real*8, dimension(3), intent(inout) :: coords
      real*8, dimension(3), intent(inout) :: vel
      real*8, dimension(3), intent(inout) :: vel_old
      real*8, intent(in) :: stepSize
      integer, intent(inout) :: firstStep

      if (firstStep == 1) then ! use Euler integration
        coords = coords - stepSize * vel
        firstStep = 0
      else ! use Adams-Bashforth 2-step integration
        coords = coords - 3.0d0/2.0d0 * stepSize * vel  + &
          1.0d0/2.0d0 * stepSize * vel_old
      end if
    end subroutine ABintegrate
  !*****************************************************************************
    function angle (A, B) result (C)
    !---------------------------------------------------------------------------
    !  angle returns the angle between two vectors
    !
    !  History:
    !
    !    Jay Polk
    !    May 10, 2021
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, a dim (3) vector
    !    Input, real ( kind = 8 ) B, a dim (3) vector
    !    Output, real ( kind = 8 ) C, the angle between A and B in radians
    !---------------------------------------------------------------------------
      implicit none

      real*8, dimension(3), intent(in) :: A
      real*8, dimension(3), intent(in) :: B
      real*8 :: C
      real*8 :: cosC

      cosC = dot(A,B)/(norm2(A)*norm2(B))
      C = acos(cosC)

    end function angle

end module iomesh_mod
