module sputter_mod
  implicit none
! sputter
  real*8, parameter :: pi = 3.14159265358979d0

  contains
  !*****************************************************************************
    function sputter_yield_carbon (energy, inc_angle) result(yield)

    !---------------------------------------------------------------------------
    !
    !  sputter_yield_carbon returns the yield for Xe on C in atoms/ion
    !
    !  Discussion:
    !
    !    This routine calculates the total yield in atoms/ion for a Xe
    !    ion impacting C at normal incidence using the Eckstein model.
    !    A multiplication factor to correct for the variation in yield
    !    with incidence angle based on the Wei model is then applied. aL
    !    and epsilon_energy calculated from these parameters:
    !    real*8, parameter :: a_0 = 0.5292d0     ! Bohr radius (angstroms)
    !    real*8, parameter :: eps_0 = 1.4184e-40 ! permittivity of free space (C^2/eV/Angstrom)
    !    real*8, parameter :: e = 1.602e-19      ! electron charge (C)
    !    real*8, parameter :: M_Xe = 131.29d0    ! Xe mass (AMU)
    !    real*8, parameter :: M_C = 12.0107d0    ! C mass (AMU)
    !    real*8, parameter :: Z_Xe = 54.0d0      ! Xe atomic number
    !    real*8, parameter :: Z_C = 6.0d0        ! C atomic number
    !    aL = (((9.0d0 * pi**2)/128.0d0 )**(1.0d0/3.0d0)) * a_0 * (Z_Xe**(2.0d0/3.0d0) &
    !    + Z_C**(2.0d0/3.0d0))**(-1.0d0/2.0d0)
    !    epsilon_energy = aL/(Z_Xe * Z_C) * (4.0d0 * pi * eps_0/(e**2)) * (M_C/(M_C + M_Xe))
    !
    !  History:
    !
    !    Jay Polk
    !    March 15, 2021
    !    Based on matlab subroutine in Paddy Simha's redeposition code
    !
    !  Reference:
    !
    !    Based on fits to the Eckstein yield model and Wei angular yield
    !    model from John Yim's paper:
    !    Yim, J.T., "A survey of xenon ion sputter yield data and fits
    !    relevant to electric propulsion spacecraft integration,"
    !    IEPC-2017-060, 2017.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) energy, ion energy in eV
    !
    !    Input, real ( kind = 8 ) inc_angle, incidence angle in radians
    !
    !    Output, real ( kind = 8 ) yield, yield in atoms/ion
    !---------------------------------------------------------------------------
      implicit none

      real*8, intent(in) :: energy
      real*8, intent(in) :: inc_angle
      real*8 :: yield


      real*8, parameter :: E_th = 21.0d0        ! Threshold energy (see Yim)
      real*8, parameter :: Q = 4.0d0            ! Fit parameter for Eckstein model
      real*8, parameter :: lambda = 0.8d0       ! Fit parameter for Eckstein model
      real*8, parameter :: mu = 1.8d0           ! Fit parameter for Eckstein model
      real*8, parameter :: beta_alpha = 0.88d0  ! beta/alpha in Wei model
      real*8, parameter :: a_alpha = 2.20d0   ! a/alpha in Wei model
      real*8, parameter :: pi = 4.0d0 * atan(1.0d0)
      real*8, parameter :: aL = 0.111716093898588  ! Lindhard screening length
      real*8, parameter :: epsilon_energy = 2.007116283414568d-006  ! epsilon/energy
      real*8 epsilon                            ! Reduced energy
      real*8 w                                  ! Parameter in Eckstein model
      real*8 s_n                                ! reduced nuclear stopping power
      real*8 P                                  ! Correction factor for non-normal incidence

      epsilon = epsilon_energy * energy
      w = epsilon + 0.1728d0 * sqrt(epsilon) + 0.008d0 * epsilon**(0.1504d0)
      s_n = (0.5d0 * log(1.0d0 + 1.2288d0 * epsilon))/w
      yield = Q * s_n * (((energy/E_th) - 1.0d0)**mu)/((lambda/w) + ((energy/E_th) - 1.0d0)**mu)

      P = (1.0d0/sqrt(1.0d0 + (beta_alpha**2) * (tan(inc_angle))**2)) * exp(0.5d0 * &
      (a_alpha**2) * (1.0d0 - (1.0d0/(1.0d0 + (beta_alpha**2) * (tan(inc_angle))**2))))
      yield = P * yield

    end function sputter_yield_carbon
  !*****************************************************************************

    function crossarray (A, B, mask) result (C)
    !---------------------------------------------------------------------------
    !  crossarray returns the cross product for an array of vectors
    !
    !  History:
    !
    !    Jay Polk
    !    March 28, 2021
    !    Renamed from cross to crossaray on May 8, 2021 to distinguish from new cross subroutine
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, an (3,N) array of N vectors
    !    Input, real ( kind = 8 ) B, an (3,N) array of N vectors
    !    Optional Input, logical mask, a vector of N values that can be used to mask the calculations
    !    Output, real ( kind = 8 ) C, a (3,N) array of N cross products of A and B vectors
    !---------------------------------------------------------------------------
      implicit none

      real*8, dimension(:,:), intent(in) :: A
      real*8, dimension(:,:), intent(in) :: B
      logical, dimension (:), optional, intent(in) :: mask
      real*8, dimension(size(A, dim = 1), size(A, dim = 2)) :: C
      real*8 :: N, M
      integer :: i, j

      M = size(A, dim = 1)
      N = size(A, dim = 2)
      ! Check that size(A) == size(B) and that dim 1 is 3
      if (.not.(present(mask))) then
        do i = 1,N
          C(1,i) = A(2,i) * B(3,i) - A(3,i) * B(2,i)
          C(2,i) = A(3,i) * B(1,i) - A(1,i) * B(3,i)
          C(3,i) = A(1,i) * B(2,i) - A(2,i) * B(1,i)
        end do
      else
        where (mask)
          C(1,:) = A(2,:) * B(3,:) - A(3,:) * B(2,:)
          C(2,:) = A(3,:) * B(1,:) - A(1,:) * B(3,:)
          C(3,:) = A(1,:) * B(2,:) - A(2,:) * B(1,:)
        end where
      end if
    end function crossarray
  !*****************************************************************************

    function cross (A, B) result (C)
    !---------------------------------------------------------------------------
    !  cross returns the cross product of two vectors
    !
    !  History:
    !
    !    Jay Polk
    !    May 8, 2021
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, a dim (3) vector
    !    Input, real ( kind = 8 ) B, a dim (3) vector
    !    Output, real ( kind = 8 ) C, a dim (3) vector cross product of A and B vectors
    !---------------------------------------------------------------------------
      implicit none

      real*8, dimension(3), intent(in) :: A
      real*8, dimension(3), intent(in) :: B
      real*8, dimension(3) :: C

      integer :: i, j

      ! Check that size(A) == size(B) == 3

      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)

    end function cross
  !*****************************************************************************

    function dotarray (A, B, mask) result (C)
    !---------------------------------------------------------------------------
    !  dotarray returns the dot product for an array of vectors
    !
    !  History:
    !
    !    Jay Polk
    !    March 30, 2021
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, an (3,N) array of N vectors
    !    Input, real ( kind = 8 ) B, an (3,N) array of N vectors
    !    Optional Input, logical mask, a vector of N values that can be used to mask the calculations
    !    Output, real ( kind = 8 ) C, a vector of N dot products of A and B vectors
    !---------------------------------------------------------------------------
      implicit none

      real*8, dimension(:,:), intent(in) :: A
      real*8, dimension(:,:), intent(in) :: B
      logical, dimension (:), optional, intent(in) :: mask
      real*8, dimension(size(A, dim = 2)) :: C
      real*8 :: N, M
      integer :: i, j

      M = size(A, dim = 1)
      N = size(A, dim = 2)
      ! Check that size(A) == size(B) and that dim 1 is 3

      C(:) = 0.0d0
      if (.not.(present(mask))) then
        C(:) = A(1,:) * B(1,:) + A(2,:) * B(2,:) + A(3,:) * B(3,:)
      else
        where (mask)
          C(:) = A(1,:) * B(1,:) + A(2,:) * B(2,:) + A(3,:) * B(3,:)
        end where
      end if
    end function dotarray
  !*****************************************************************************
  !*****************************************************************************

    function dot (A, B) result (C)
    !---------------------------------------------------------------------------
    !  dot returns the dot product of two vectors
    !
    !  History:
    !
    !    Jay Polk
    !    May 15, 2021
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, a vector
    !    Input, real ( kind = 8 ) B, a second vector
    !    Output, real ( kind = 8 ) C, the dot product of A and B
    !---------------------------------------------------------------------------
      implicit none

      real*8, dimension(:), intent(in) :: A
      real*8, dimension(:), intent(in) :: B
      real*8 :: C
      integer :: i 

      C = 0.0d0
      do i = 1, size(A)
        C = C + A(i) * B(i)
      end do

    end function dot
  !*****************************************************************************

    subroutine triangle_ray_intersection (ion_idx_in, origpt, dirvec, vert1, vert2, vert3, intersect, t, u, v, xcoor, pTypeIn, lTypeIn, borderIn, epsIn, fullReturnIn)
    !---------------------------------------------------------------------------
    !TRIANGLERAYINTERSECTION Ray/triangle intersection.
    !    INTERSECT = TriangleRayIntersection(ORIG, DIR, VERT1, VERT2, VERT3)
    !      calculates ray/triangle intersections using the algorithm proposed
    !      BY moduleller and Trumbore (1997), implemented as highly vectorized
    !      MATLAB code. The ray starts at ORIG and points toward DIR. The
    !      triangle is defined by vertix points: VERT1, VERT2, VERT3. All input
    !      arrays are in Nx3 or 1x3 format, where N is number of triangles or
    !      rays.
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) orig, ion energy in eV
    !
    !
    !   [INTERSECT, T, U, V, XCOOR] = TriangleRayIntersection(...)
    !     Returns:
    !     * Intersect - boolean array of length N (num triangles) informing which line and
    !                 triangle pair intersect
    !     * t   - distance from the ray origin to the intersection point in
    !             units of |dir|. Provided only for line/triangle pair that
    !             intersect unless 'fullReturn' parameter is true.
    !     * u,v - barycentric coordinates of the intersection point
    !     * xcoor - cartesian coordinates of the intersection point
    !
    !   TriangleRayIntersection(...,'param','value','param','value'...) allows
    !    additional param/value pairs to be used. Allowed parameters:
    !    * planeType - 'one sided' or 'two sided' (default) - how to treat
    !        triangles. In 'one sided' version only intersections in single
    !        direction are counted and intersections with back facing
    !           tringles are ignored
    !    * lineType - 'ray' (default), 'line' or 'segment' - how to treat rays:
    !        - 'line' means infinite (on both sides) line;
    !        - 'ray' means infinite (on one side) ray comming out of origin;
    !        - 'segment' means line segment bounded on both sides
    !    * border - controls border handling:
    !        - 'normal'(default) border - triangle is exactly as defined.
    !           Intersections with border points can be easily lost due to
    !           rounding errors.
    !        - 'inclusive' border - triangle is marginally larger.
    !           Intersections with border points are always captured but can
    !           lead to double counting when working with surfaces.
    !        - 'exclusive' border - triangle is marginally smaller.
    !           Intersections with border points are not captured and can
    !           lead to under-counting when working with surfaces.
    !    * epsilon - (default = 1e-5) controls border size
    !    * fullReturn - (default = false) controls returned variables t, u, v,
    !        and xcoor. By default in order to save time, not all t, u & v are
    !        calculated, only t, u & v for intersections can be expected.
    !        fullReturn set to true will force the calculation of them all.
    !
    ! ALGORITHM:
    !  Function solves
    !        |t|
    !    M * |u| = (o-v0)
    !        |v|
    !  for [t; u; v] where M = [-d, v1-v0, v2-v0]. u,v are barycentric coordinates
    !  and t - the distance from the ray origin in |d| units
    !  ray/triangle intersect if u>=0, v>=0 and u+v<=1
    !
    ! NOTE:
    !  The algorithm is able to solve several types of problems:
    !  * many faces / single ray  intersection
    !  * one  face  / many   rays intersection
    !  * one  face  / one    ray  intersection
    !  * many faces / many   rays intersection
    !  In order to allow that to happen all imput arrays are expected in Nx3
    !  format, where N is number of vertices or rays. In most cases number of
    !  vertices is different than number of rays, so one of the imputs will
    !  have to be cloned to have the right size. Use "repmat(A,size(B,1),1)".
    !
    ! Based on:
    !  *"Fast, minimum storage ray-triangle intersection". Tomas M?ller and
    !    Ben Trumbore. Journal of Graphics Tools, 2(1):21--28, 1997.
    !    http://www.graphics.cornell.edu/pubs/1997/MT97.pdf
    !  * http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/raytri/
    !  * http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/raytri/raytri.c
    !
    ! Author:
    !    Jarek Tuszynski (jaroslaw.w.tuszynski@leidos.com)
    !
    ! License: BSD license (http://en.wikipedia.org/wiki/BSD_licenses)
    !---------------------------------------------------------------------------
    implicit none
    real*8, dimension(3), intent(in) :: origpt
    real*8, dimension(3), intent(in) :: dirvec
    real*8, dimension(:,:), intent(in) :: vert1
    real*8, dimension(:,:), intent(in) :: vert2
    real*8, dimension(:,:), intent(in) :: vert3
    logical, dimension(size(vert1, dim = 2)), intent(out) :: intersect
    real*8, dimension(size(vert1, dim = 2)), intent(out) :: t
    real*8, dimension(size(vert1, dim = 2)), intent(out) :: u
    real*8, dimension(size(vert1, dim = 2)), intent(out) :: v
    real*8, dimension(3, size(vert1, dim = 2)), intent(out) :: xcoor
    character(len = *), optional, intent(in) :: pTypeIn
    character(len = *), optional, intent(in) :: lTypeIn
    character(len = *), optional, intent(in) :: borderIn
    logical, optional, intent(in) :: fullReturnIn
    real*8, optional, intent(in) :: epsIn

    character(len = 10) :: planeType
    character(len = 10) :: lineType
    character(len = 10) :: border
    logical :: fullReturn
    real*8 :: eps
    real*8, dimension(3, size(vert1, dim = 2)) :: orig
    real*8, dimension(3, size(vert1, dim = 2)) :: dir
    real*8, dimension(3, size(vert1, dim = 2)) :: edge1
    real*8, dimension(3, size(vert1, dim = 2)) :: edge2
    real*8, dimension(3, size(vert1, dim = 2)) :: tvec
    real*8, dimension(3, size(vert1, dim = 2)) :: pvec
    real*8, dimension(3, size(vert1, dim = 2)) :: qvec
    real*8, dimension(size(vert1, dim = 2)) :: det
    logical, dimension(size(vert1, dim = 2)) :: angleOK
    logical, dimension(size(vert1, dim = 2)) :: ok
    integer :: N
    integer :: i, j
    real*8, parameter :: d_qnan = transfer((/ Z'00000000', Z'7FF80000' /),1.0_8)
    real*8 :: zero

    ! PSB added input ion_idx_in
    integer, parameter :: INTERSECT_LOG_UNIT = 20 ! PSB Choose an unused unit number
    logical, save :: first_intersect_write = .true. ! PSB Flag for writing header
    integer, dimension(size(vert1, dim=2)) :: intersect_int
    integer, intent(in) :: ion_idx_in ! New Ion index to be saved

    ! Debug print
    if (present(epsIn)) then
      eps = epsIn
    else
      eps = 1.0d-5
    end if
    if (present(pTypeIn)) then
      planeType  = pTypeIn
    else
      planeType = 'one sided'
    end if
    if (present(lTypeIn)) then
      lineType   = lTypeIn
    else
      lineType = 'ray'
    end if
    If (present(borderIn)) then
      border = borderIn
    else
      border = 'normal'
    end if
    if (present(fullReturnIn)) then
      fullReturn = fullReturnIn
    else
      fullReturn = .false.
    end if

    !! Test array dimensions here...

    !! Replicate orig and dir to make arrays the same size as the triangle vertex arrays
    N = size(vert1, dim = 2)
    orig(:,:) = spread(origpt(:), 2, N)
    dir(:,:) = spread(dirvec(:), 2, N)


    !! Set up border parameter
    select case (border)
      case ('normal')
        zero=0.0d0
      case ('inclusive')
        zero=eps
      case ('exclusive')
        zero=-eps
      case default
        ! Do error handling here
        ! error('Border parameter must be either "normal", "inclusive" or "exclusive"')
      end select

    !! initialize default output  (logical intersect should be false by default; others seem unnecessary)
    intersect = .false.
    t = huge(1.0d0)
    u = t
    v = t

    !! Find faces parallel to the ray

    do j = 1, 3
      edge1(j,:) = vert2(j,:) - vert1(j,:)         ! find vectors for two edges sharing vert1
      edge2(j,:) = vert3(j,:) - vert1(j,:)
      tvec(j,:)  = orig(j,:) - vert1(j,:)          ! vector from vert1 to ray origin
    end do

    pvec  = crossarray(dir, edge2)  ! begin calculating determinant - also used to calculate U parameter

    det = dotarray(edge1, pvec)          ! determinant of the matrix M = dotarray(edge1,pvec)
    select case (planeType)
      case ('two sided')            ! treats triangles as two sided
        angleOK(:) = (abs(det(:)) > eps) ! if determinant is near zero then ray lies in the plane of the triangle
      case ('one sided')            ! treats triangles as one sided
        angleOK(:) = (det(:) > eps)
      case default
        ! Do error handling here
        ! error('Triangle parameter must be either "one sided" or "two sided"')
    end select

    if (count(angleOK == .true.) == 0) return  ! if all parallel then no intersections

    !! Different behavior depending on one or two sided triangles
    where (.not.angleOK)              ! change to avoid division by zero
      det = d_qnan
    end where

    u(:) = dotarray(tvec, pvec)/det(:)     ! 1st barycentric coordinate

    if (fullReturn) then
      ! calculate all variables for all line/triangle pairs
      qvec = crossarray(tvec, edge1)    ! prepare to test V parameter

      v(:) = dotarray(dir, qvec)/det(:)        ! 2nd barycentric coordinate

      t = dotarray(edge2, qvec)/det(:)    ! 'position on the line' coordinate

      ! test if line/plane intersection is within the triangle
      ok(:)   = (angleOK(:) .and. (u(:) >= (-zero)) .and. (v(:) >= (-zero)) .and. ((u(:) + v(:)) <= (1.0 + zero)))
    else
      ! limit some calculations only to line/triangle pairs where it makes
      ! a difference. It is tempting to try to push this concept of
      ! limiting the number of calculations to only the necessary to "u"
      ! and "t" but that produces slower code

      v = d_qnan
      t = v

      ok(:) = (angleOK(:) .and. (u(:) >= (-zero)) .and. (u(:) <= (1.0d0 + zero)))  ! mask
      ! if all line/plane intersections are outside the triangle then no intersections
      if (.not.any(ok)) return

      qvec = crossarray(tvec, edge1, ok) ! prepare to test V parameter

      v(:) = dotarray(dir,qvec,ok)/det(:)  ! 2nd barycentric coordinate

      t(:) = dotarray(edge2, qvec, ok)/det(:)

      ! test if line/plane intersection is within the triangle
      ok(:) = (ok(:) .and. (v(:) >= -zero) .and. ((u(:) + v(:)) <= (1.0d0 + zero)))
    end if

    !! Test where along the line the line/plane intersection occurs
    select case (lineType)
      case ('line')      ! infinite line
        intersect(:) = ok(:)
      case ('ray')       ! ray is bound on one side
        intersect(:) = (ok(:) .and. (t >= (-zero))); ! intersection on the correct side of the origin
      case ('segment')   ! segment is bound on two sides
        intersect(:) = (ok(:) .and. (t >= (-zero)) .and. (t <= (1.0d0 + zero))); ! intersection between origin and destination
      case default
        ! error('lineType parameter must be either "line", "ray" or "segment"');
        ! Do error handling here
    end select

    !! calculate intersection coordinates if requested
    where (ok)
      xcoor(1,:) = vert1(1,:) + edge1(1,:) * u(:) + edge2(1,:) * v(:)
      xcoor(2,:) = vert1(2,:) + edge1(2,:) * u(:) + edge2(2,:) * v(:)
      xcoor(3,:) = vert1(3,:) + edge1(3,:) * u(:) + edge2(3,:) * v(:)
    end where

    ! PSB ------- Added for file writing ----------------
    ! --- Write intersect array to file ---
    OPEN(UNIT=INTERSECT_LOG_UNIT, FILE='triangle_intersection_log.txt', &
         STATUS='OLD', ACTION='WRITE', POSITION='APPEND', IOSTAT=i)
    IF (i /= 0) THEN ! Check if file exists and can be opened
        ! If file doesn't exist (e.g., first run or deleted), create it.
        OPEN(UNIT=INTERSECT_LOG_UNIT, FILE='triangle_intersection_log.txt', &
             STATUS='REPLACE', ACTION='WRITE', IOSTAT=i)
        IF (i == 0) THEN
            WRITE(INTERSECT_LOG_UNIT, '(A)') '# ion_idx, intersect(1), intersect(2), ...'
            first_intersect_write = .false. ! Reset flag after header write
        ELSE
            PRINT *, "Error: Could not open or create triangle_intersection_log.txt (IOSTAT:", i, ")"
            RETURN ! Exit subroutine if cannot write log
        END IF
    ELSE IF (first_intersect_write) THEN
        first_intersect_write = .false. ! Only write header once per program execution
    END IF
    ! Method 1: Convert to integer array (more explicit)
    intersect_int = merge(1, 0, intersect) ! 1 if .true., 0 if .false.
    WRITE(INTERSECT_LOG_UNIT, '(I8, A, *(I1, : , ","), I1)') ion_idx_in, ", ", intersect_int(1:N-1), intersect_int(N)
    ! Format: Ion Index, ", ", followed by each boolean (0 or 1)
    ! This assumes N (number of triangles) isn't excessively large, otherwise, the line could be very long.
    CLOSE(INTERSECT_LOG_UNIT)

    end subroutine triangle_ray_intersection

end module sputter_mod
