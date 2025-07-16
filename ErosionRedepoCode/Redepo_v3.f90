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
 implicit none

 !-----------------------------------------------------------------------------
 !  Control Parameters
 !-----------------------------------------------------------------------------
 real*8, parameter :: moly_mass = 1.593E-22  ! Mass of a molybdenum (or grid material) atom in grams
 integer, parameter :: nionsample = 1000  ! Unitless, amount of CEX ions used for analysis
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
 integer :: i, j, inode, icell, iion, irefl
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
  CHARACTER (LEN=20) :: step_str  !PSB
  CHARACTER (LEN=50) :: plot_subfolder  !PSB
 CHARACTER (LEN=20) :: plot_name_str !PSB
 INTEGER,DIMENSION(nionsample) :: ion_reflection_counts !PSB -> added to track number of reflections per ion
 print *, "Variables have been declared" !PSB
 !-----------------------------------------------------------------------------
 !  Read in CEX ion data
 !-----------------------------------------------------------------------------
 select case (ion_data_source)
   case (0)     !  Read in CEX ion data from original CEX3D output file
   case (1)     !  Read in CEX ion data from processed ion input file
     ionfilename = 'CEX_ion_data'
     call read_ion_file(ionfilename, domain, simulation_time, npandgions)
   case (2)     !  Create a uniform flux of ions for testing
   case default
     print*, "Improper input for ion_data_source (must be 0, 1, or 2)"
 end select

 !-----------------------------------------------------------------------------
 !  Setup initial grid
 !-----------------------------------------------------------------------------
 select case (mesh_data_source)
   case (0)    !  Generate mesh for initial flat accel grid face
     call mesh_flat_surface(domain, npt, vcl, vcl3d, maxntri, nfacept, accel_thickness)
     call triangulate_surface(domain, npt, vcl, vcl3d, ntri, til, nfacetri, cell)
     PRINT *, "Number of points (npt): ", npt 
     PRINT *, "Number of triangles (ntri): ", ntri 
     PRINT *, "Number of face points (nfacept): ", nfacept 
     PRINT *, "Maximum allowed triangles (maxntri): ", maxntri
     !---------------------------------------------------------------------------
     ! Save the vertex coordinate list and mesh triangle index list for plotting
     !---------------------------------------------------------------------------
     open(13, FILE = 'VCL') 
     do i = 1,npt
       write(13,*) i, vcl(1,i), vcl(2,i)
     end do
     close(13)
     open(12, FILE = 'Triangles')
     do i = 1,ntri
       WRITE(12,*) til(1,i), til(2,i), til(3,i)
     end do
     close(12)


     !--------------------------------------------------------------------------- 
     ! PSB Save 3D mesh data for Python plotting (new files) 
     !--------------------------------------------------------------------------- 
     OPEN(UNIT=14, FILE='mesh_points_3d.txt', STATUS='REPLACE', ACTION='WRITE') 
     WRITE(14, '(A)') '# X Y Z' ! Header for columns 
      DO i = 1, npt 
        WRITE(14, '(3(ES15.6))') vcl3d(1,i), vcl3d(2,i), vcl3d(3,i) 
      END DO 
    CLOSE(14) 
    OPEN(UNIT=15, FILE='mesh_triangles_3d.txt', STATUS='REPLACE', ACTION='WRITE') 
    WRITE(15, '(A)') '# Node1 Node2 Node3 (1-indexed)' ! Header for columns 
      DO i = 1, ntri ! Write 1-indexed triangle node IDs 
        WRITE(15, '(3(I8))') til(1,i), til(2,i), til(3,i) 
      END DO 
    CLOSE(15) 
    ! Call Python script to generate plot 
    plot_subfolder = 'mesh_plots' ! Define your desired subfolder name 
    ! Create the subfolder if it doesn't exist 
    ! Use platform-independent way if possible, but 'mkdir' is common on Linux/macOS 
    ! For Windows, 'mkdir' also works. 
    CALL EXECUTE_COMMAND_LINE("mkdir -p " // TRIM(plot_subfolder), wait=.true., exitstat=j) 
    IF (j /= 0) THEN 
      PRINT *, "Warning: Could not create directory ", TRIM(plot_subfolder), ". Error code: ", j 
      PRINT *, "Plot will be saved in the current directory instead." 
      plot_subfolder = '.' ! Fallback to current directory 
    END IF ! Convert the integer 'erosionStep' to a string 
    ! Construct the command string to pass parameters to Python 
    plot_name_str = "pre_boundaries"
    command_string = "python plot_mesh.py " // TRIM(plot_subfolder) // " " // TRIM(plot_name_str) 
    CALL EXECUTE_COMMAND_LINE(TRIM(command_string), wait=.true.)



    call mesh_boundaries(domain, npt, vcl, vcl3d, ntri, til, cell, nfacept)
     do icell = 1, ntri
       vert1(:, icell) = vcl3d(:,til(1,icell))     ! Define vertices of each cell
       vert2(:, icell) = vcl3d(:,til(2,icell))
       vert3(:, icell) = vcl3d(:,til(3,icell))
     end do
     startTime = dclock()
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
     stopTime = dclock()
     print *, "cell: elapsed time:", stopTime - startTime
     startTime = dclock()
     do inode = 1,nfacept
       ! Define coordinates of node from vertex coordinate list
       node(inode)%coords(:) = vcl3d(:,inode)
       ! Calculate other node parameters
       call calc_node_params(node(inode), inode, domain, nfacetri, cell)
     end do
     stopTime = dclock()
     print *, "node: elapsed time:", stopTime - startTime
   case (1)    !  Read in mesh from  a restart file
   case (2)    !  Generate a mesh for a v-groove test case
   case default
     print*, "Improper input for mesh_data_source (must be 0, 1, or 2)"
 end select

 
 

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
 time_count=1
 max_pit_depth(time_count) = 0
 max_swing_depth(time_count) = 0
 groove_centre_depth(time_count) = 0
 time_stamp(time_count) = 0
 time_count = time_count + 1

 delta_x_old = vcl(1,:)
 delta_y_old = vcl(1,:)
 delta_z_old = vcl(1,:)
 !
 !  Erosion iteration loop
 !
 print*, "Starting erosion iterations..!"
 open(113, FILE = 'VCL3D')
 write(113,*) "time, x, y, z"
 do inode = 1,nfacept
   write(113,'("         0",3(",",ES13.6))') vcl3d(1,inode), vcl3d(2,inode), vcl3d(3,inode)
 end do

EROSION: do erosionStep = 1, 2 ! PSB: IMPORTANT CHANGED 5 TO 2 
  print *, "----------------------------------------------------------------------" !PSB
  print *, "                          Erosion Step:", erosionStep !PSB
   erosionTime = erosionStep * erosionStepSize
   !  Store existing surface geometry in arrays(timestep)--save history of surface geometry

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

   print *, "Number of ions chosen to be sampled is nionsample which is", nionsample
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
   open(112, FILE = 'LostBoys')
   DOION: do iion = 1, nionsample !PSb: Iterate through each of the randomly selected ions
    ! PSB: This subroutine checks whether a ray (or line) intersects a triangle in 3D space.
    ! outputs a boolean array indicating whether the ray intersects with each triangle, 
    ! the distance from the ray origin to the intersection point, 
    ! barycentric coordinates of the intersection point,
    ! cartesian coordinates of the intersection point, if requested
     call triangle_ray_intersection (ion(iion)%origin, ion(iion)%dir, vert1, vert2, vert3, intersect, t, u, v, xcoor,borderIn = 'inclusive')
     
     ! PSB: Calculates reflection attempts or interactions
     ! Track an ions reflective path through the domain (Reflections off symmetry boundaries until it hits the sputtering surface, leaves domain, or has reflected max times)
     ! Track how many times each ion reflected
     do irefl = 1, maxReflections_ions
       num_hits = count(intersect) !PSB: num_hits is the total number of triangles that the current ion's ray intersected
       impact_tri = findloc(intersect, value = .true.)  ! PSB: returns indices of all triangles the ray intersected
       
       for each value of impact_tri {
       impact_cell = impact_tri(i) ! PSB: start with the first impact cell
       if (num_hits == 0) then
         ! ion did not strike any surface in domain (small fraction of total at edges or with vz near zero)
         WRITE(112,"(2I8,6E15.6)") iion, irefl, ion(iion)%origin, ion(iion)%dir
         nlostboys = nlostboys + 1
         !PSB: Example
         ! If iion = 12345 (ion number), ierefl = 1 (first attempt at tracing this ion and there was no hit)
         ! ion(iion)%origin = ion starting position like [0.1234..., 0.987654..., 2.500...], similar with direction in x, y, and z
         ! Lost boys would look like this:
         ! 12345 1 1.234567E-001 9.876543E-001 2.500000E+000 1.234568E-002 -9.876543E-002 -9.900000E-001
         ! PSB: I SEE MANY IONS HAVE A NAN DIRECTION - CHECK WHY THAT IS
         exit
       else if (cell(impact_cell)%on_bndry == -1) then
         ! ion left through top of domain (should not happen)
         print*, "ion number ", iion, "left through top of domain (should not happen)"
         exit
       else if (cell(impact_cell)%on_bndry == -2) then
         ! ion hit bottom of domain (should only happen after punch-through)
         ! print*, "ion number ", i, "hit bottom of domain (should only happen after punch-through)"
         nbottomboys = nbottomboys + 1
         exit
       else if (cell(impact_cell)%on_bndry > 0) then
         ! ion hit symmetry boundary (reflect)
         ! PSB: Vector reflection formula
         ! reflectedray is a temporary 3 element vector to store the newly calculated reflected direction
         reflectedray = -2 * dot(ion(iion)%dir,cell(impact_cell)%normal) + ion(iion)%dir
         ! update ion directon to reflected ray
         ion(iion)%dir = reflectedray
         ion_reflection_counts(iion) = ion_reflection_counts(iion) + 1 !psb added
         call triangle_ray_intersection ! PSB added because if an ion is reflected, even if we know its direction, we need to know what triangle it will hit
          do again for impact_tri = findloc(intersect, value = .true.)
       else
         ! ion hit surface (calculate sputtered mass and share to nodes)
         inc_angle = angle(ion(iion)%dir(:), -1.0d0 * cell(impact_cell)%normal(:)) !PSB: incidence angle of ion with respect to surface normal
         yield = sputter_yield_carbon(ion(iion)%KE, inc_angle) ! PSB: The function returns the sputter yield in atoms per ion, which represents how many atoms of carbon are ejected from the surface when one Xenon ion strikes it.
         mass_loss = yield * ion(iion)%qCEXion/e * massCarbon ! PSB: Calculate the mass loss due to this single ion impact
         print*, "yield (atoms/ion) = ", yield
         cell(impact_cell)%mass_loss = cell(impact_cell)%mass_loss + mass_loss !PSB: Add the mass loss form the current ion impact to total mass_loss accumulated for the specific cell where the triangle was hit
         
         ! PSB: Use barycentric interpolation: distribute mass_loss from triangle to its nodes
         ! u + v + w = 1 where u is for second vertex, v, for third, etc. are weights for eriosion
         node(til(2,impact_cell))%mass_loss = node(til(2,impact_cell))%mass_loss + u(impact_cell) * mass_loss 
         node(til(3,impact_cell))%mass_loss = node(til(3,impact_cell))%mass_loss + v(impact_cell) * mass_loss
         node(til(1,impact_cell))%mass_loss = node(til(1,impact_cell))%mass_loss + (1 - u(impact_cell) - v(impact_cell)) * mass_loss
         sum = node(til(2,impact_cell))%mass_loss + node(til(3,impact_cell))%mass_loss + &
           node(til(1,impact_cell))%mass_loss !PSB sum should be equal to the mass losses distributed to the three nodes of impact_cell, ie, cell(impact_cell)%mass_los
       end if
       print*, "Ion Number:", iion, "  Triangle Cells the ion hit:", impact_tri, "  Total cell loss:" , cell(impact_cell)%mass_loss
       }

     end do
  end do DOION

   close(112)

   print*, "No. ions that didn't hit anything: ", nlostboys
   print*, "No. ions that hit the bottom: ", nbottomboys
   print*, "No. of ions total: ", npt

    OPEN(UNIT=114, FILE='ion_reflection_counts.txt', STATUS='REPLACE') 
    WRITE(114, '(A)') '# Sampled_Ion_ID Original_Ion_ID Reflection_Count' 
    DO iion = 1, nionsample ! WRITE(114, '(3I12)') iion, ionindex(iion), ion_reflection_counts(iion) 
    END DO
    CLOSE(114)
    PRINT *, "Ion reflection counts saved to ion_reflection_counts.txt"


   !-----------------------------------------------------------------------------
   !  Redeposition section
   !-----------------------------------------------------------------------------
   !-----------------------------------------------------------------------------
   !  Geometry update section
   !-----------------------------------------------------------------------------

  !  do inode = 1, nfacept
  !    node(inode)%vel = node(inode)%mass_loss/simulation_time/(node(inode)%area * densityCarbon) * node(inode)%normal  ! vel is positive inward
  !    call ABintegrate(node(inode)%coords, node(inode)%vel, node(inode)%vel_old, erosionStepSize * 3600, node(inode)%firstStep)
  !    vcl3d(:,inode) = node(inode)%coords
  !    vcl(:,inode) = node(inode)%coords(1:2)
  !    node(inode)%vel_old = node(inode)%vel
  !  end do
  !  print *, "Erosion Time         :", erosionTime !PSB
  !  print *, "Node 226 Mass Loss   :", node(226)%mass_loss !PSB
   
  !   print *, "Node 226 Coordinates:"
  !   print '(A, F10.4)', "  x: ", node(226)%coords(1)! PSB
  !   print '(A, F10.4)', "  y: ", node(226)%coords(2)! PSB
  !   print '(A, F10.4)', "  z: ", node(226)%coords(3)! PSB


  !  npt = nfacept
  !  call triangulate_surface(domain, npt, vcl, vcl3d, ntri, til, nfacetri, cell)
  !  call mesh_boundaries(domain, npt, vcl, vcl3d, ntri, til, cell, nfacept)
  !  do icell = 1, ntri
  !    cell(icell)%til(:) = til(:,icell)! Define the triangle incidence list for cell
  !    cell(icell)%vert1(:) = vcl3d(:,cell(icell)%til(1)) ! Define vertices of cell
  !    cell(icell)%vert2(:) = vcl3d(:,cell(icell)%til(2))
  !    cell(icell)%vert3(:) = vcl3d(:,cell(icell)%til(3))
  !    call calc_cell_params(cell(icell), icell)! Calculate other cell parameters
  !  end do
  !  do inode = 1,nfacept

  !    node(inode)%coords(:) = vcl3d(:,inode)! Define coordinates of node from vertex coordinate list

  !    call calc_node_params(node(inode), inode, domain, nfacetri, cell)! Calculate other node parameters
  !  end do
  !  do inode = 1,nfacept
  !    write(113,'(F10.2,3(",",ES13.6))') erosionTime, vcl3d(1,inode), vcl3d(2,inode), vcl3d(3,inode)
  !  end do





   ! ------------------------------ PYTHON PLOTTING -------------------------------

  !  OPEN(UNIT=14, FILE='mesh_points_3d.txt', STATUS='REPLACE', ACTION='WRITE') 
  !  WRITE(14, '(A)') '# X Y Z (Current Erosion Step)' 
  !  DO i = 1, npt 
  !   WRITE(14, '(3(ES15.6))') vcl3d(1,i), vcl3d(2,i), vcl3d(3,i) 
  ! END DO 
  ! CLOSE(14) 
  ! OPEN(UNIT=15, FILE='mesh_triangles_3d.txt', STATUS='REPLACE', ACTION='WRITE') 
  ! WRITE(15, '(A)') '# Node1 Node2 Node3 (1-indexed)' 
  ! DO i = 1, ntri 
  !   WRITE(15, '(3(I8))') til(1,i), til(2,i), til(3,i) 
  ! END DO 
  ! CLOSE(15)
  !   plot_subfolder = 'mesh_plots' 
  !   CALL EXECUTE_COMMAND_LINE("mkdir -p " // TRIM(plot_subfolder), wait=.true., exitstat=j) 
  !   IF (j /= 0) THEN 
  !     PRINT *, "Warning: Could not create directory ", TRIM(plot_subfolder), ". Error code: ", j 
  !     PRINT *, "Plot will be saved in the current directory instead." 
  !     plot_subfolder = '.' ! Fallback to current directory 
  !   END IF ! Convert the integer 'erosionStep' to a string 
  !   WRITE(step_str, '(I0)') erosionStep ! I0 format specifier for no leading/trailing spaces
  !   plot_name_str = "erosion_step_" // TRIM(step_str)
  !   command_string = "python plot_mesh.py " // TRIM(plot_subfolder) // " " // TRIM(plot_name_str) 
  !   CALL EXECUTE_COMMAND_LINE(TRIM(command_string), wait=.true.)

 end do EROSION
 !-----------------------------------------------------------------------------
 !  Deallocate dynamic arrays
 !-----------------------------------------------------------------------------
 deallocate(ion)
 close(113)



end program Redepo



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