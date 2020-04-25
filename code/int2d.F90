module int2d

    !various subroutines involved in the interpolation and extrapolation in 2D planes
    
    use pwl_interp_2d_scattered
    
    implicit none
    
contains
    
    
    !These subroutines linearly interpolate 2D data using Delauney Triangulation
    
    subroutine getplaneeq(node_xy,zd,cross,dim_num,node_num)

    !get the coefficients of the equation for a plane based on 3 points

implicit none
  real :: vector(3,2) !2 lots of x,y,z vectors of 3 points on a plane
  real, intent(out) :: cross(3) !coefficients for plane equation
  integer, intent(in) :: dim_num,node_num
  real ( kind = 8 ), intent(in), dimension (dim_num,node_num) :: node_xy !input x,y coordinates with data
  real ( kind = 8 ), intent(in) :: zd(node_num) !input z data matching the coordinates
  real temp(4)

    !get 2 vectors on the plane
    !vector 1
          vector(1,1) = node_xy(1,2) - node_xy(1,1) !x
          vector(2,1) = node_xy(2,2) - node_xy(2,1) !y
          vector(3,1) = zd(2) - zd(1)               !z
          !vector 2
          vector(1,2) = node_xy(1,3) - node_xy(1,1)
          vector(2,2) = node_xy(2,3) - node_xy(2,1)
          vector(3,2) = zd(3) - zd(1)                           !height is the same
      
          !cross product
          temp(1) = vector(2,1) * vector(3,2) - vector(3,1) * vector(2,2)
          temp(2) = vector(3,1) * vector(1,2) - vector(1,1) * vector(3,2)
          temp(3) = vector(1,1) * vector(2,2) - vector(2,1) * vector(1,2)
        
            !substitute in first point to get 'd' coefficient
          temp(4)  = temp(1)*node_xy(1,1) + temp(2)*node_xy(2,1) + temp(3)*zd(1)
          
          !restructure the equation to be in terms of z
          cross(1) = -temp(1)/temp(3)
          cross(2) = -temp(2)/temp(3)
          cross(3) = temp(4)/temp(3)
          

          
    end subroutine
    
    
    
    
!This subroutine performs interpolation of a number of layers, including any pre-processing of the delauney triangulation. An all-in-one solution.
!This subroutine also does other checks, for example to see if the points are colinear, in which case triangulation will fail.
subroutine interp_layers(dim_num,node_num,nxew,nyew,node_xy_in,zd_in,xyi,zi,nlayer,plocation,lpdepths)

implicit none


  integer ( kind = 4 ), intent(in) :: dim_num
  integer ( kind = 4 ), intent(in) :: nxew,nyew !dimensions of the site
  integer ( kind = 4 ), intent(in) :: node_num !number of input coordinates to interpolate from
  real ( kind = 8 ),intent(in) :: node_xy_in(:,:) !input x,y coordinates with data
  real ( kind = 8 ) :: node_xy(dim_num,node_num) !input x,y coordinates with data
  real (kind = 8),intent(in) :: zd_in(:,:)
  real ( kind = 8 ) :: zd(nlayer-1,node_num) ! input z data matching the coordinates
  integer, intent(in) :: nlayer !number of layers to interpolate
  real(8) tempval(nlayer-1)
  
  !extra variables in case points are colinear and need to be duplicated with an offset
  real(8), allocatable :: node_xy2(:,:),zd2(:,:)
  integer node_num2

  real ( kind = 8 ), intent(in) :: xyi(dim_num,nxew*nyew) !coordinates of points to interpolate
  real ( kind = 8 ), intent(out) :: zi(nlayer-1,nxew,nyew) !z data of points to interpolate
  
  integer i,j,k
  
  !triangulation variables
  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ), allocatable :: element_neighbor(:,:)
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), allocatable :: triangle(:,:)
  
  real rise(node_num-1),run(node_num-1),slopes(node_num-1)
  real expandfactor !factor by which to expand the hull
  
  !point interpolation at pile location stuff
  integer, optional, intent(in) :: plocation(:,:)
  real(8), optional, intent(out) :: lpdepths(:,:)
  !logical, intent(in) :: getsurface
  
  !convex hull stuff
  REAL :: Hull(2,0:node_num+1)          !x,y coords of the convex hull
  REAL :: zHull(nlayer-1,0:node_num+1)   !layer height information of the convex hull
  real temp(2)
  real xcomp(2),ycomp(2),vlen(2) !x,y components and length of 2 adjacent vectors from a point on the convex hull
  integer nHull
  
  !copy the date, since it'll be manipulated but we don't want to tamper with the original
  node_xy = node_xy_in
  zd = zd_in
  
  
    !Check to see if the points are colinear
    !If they are, then create a duplicate set of points and put an offset between them
    !at 90 degrees to the line, forming a series of rectangles to triangulate.
  if (node_num > 2) then
    do i = 1,node_num-1
          rise(i) = node_xy(2,1) - node_xy(2,i+1)
          run(i) = node_xy(1,1) - node_xy(1,i+1)
    end do
    if (rise(1) == 0) then !calculate slopes, try to ensure no divide by zero
        slopes = rise/run
    else
        slopes = run/rise
    end if
    if(all(abs(slopes-slopes(1)) <= 0.001)) then !points are indeed colinear, or close to it
        node_num2 = node_num*2  !solve the problem by duplicating the points and providing offsets at 90 degrees to the slope
        allocate(node_xy2(2,node_num2),zd2(nlayer-1,node_num2))
        expandfactor = min(abs(run(1)),abs(rise(1)))
        if(expandfactor == 0) expandfactor = 1
        node_xy2(1,:node_num) = node_xy(1,:) - 2*nxew*rise(1)/expandfactor  !offset one way
        node_xy2(2,:node_num) = node_xy(2,:) + 2*nyew*run(1)/expandfactor
        node_xy2(1,node_num+1:) = node_xy(1,:) + 2*nxew*rise(1)/expandfactor !and the other way
        node_xy2(2,node_num+1:) = node_xy(2,:) - 2*nyew*run(1)/expandfactor
        zd2(:,:node_num) = zd
        zd2(:,node_num+1:) = zd
    else if(node_num > 3) then !if points are not colinear, and there are more than 3 (since 3 points are reserved to define a plane) then do regular interpolation
        !if(node_num == 3) return
        
        !Extrapolate by creating and expanding a convex hull. Seeds sorted coordinates
        
        !sort in increasing order for x coordinate
        do k = node_num,2,-1	 
			do j = node_num-1,node_num-k+1,-1
				if (node_xy(1,j+1) < node_xy(1,j)) then 
					temp = node_xy(:,j)
					node_xy(:,j) = node_xy(:,j+1)
					node_xy(:,j+1) = temp
                    tempval = zd(:,j)
                    zd(:,j) = zd(:,j+1)
                    zd(:,j+1) = tempval
				end if
			end do
        end do   
        
        !secondary sort in increasing order for y coordinate
        do k = node_num,2,-1	 
			do j = node_num-1,node_num-k+1,-1
				if (node_xy(1,j+1) <= node_xy(1,j) .and. node_xy(2,j+1) < node_xy(2,j)) then 
					temp = node_xy(:,j)
					node_xy(:,j) = node_xy(:,j+1)
					node_xy(:,j+1) = temp
                    tempval = zd(:,j)
                    zd(:,j) = zd(:,j+1)
                    zd(:,j+1) = tempval
				end if
			end do
		end do   
        
        call ConvHull(node_num,node_xy,zd,nlayer,nHull,Hull,zHull)
        
        node_num2 = node_num + nHull    
        allocate(node_xy2(2,node_num2),zd2(nlayer-1,node_num2))
        node_xy2(:,:node_num) = node_xy
        zd2(:,:node_num) = zd
        
        !find tangents on the convex hull by finding the vectors of the points either side of the point of interest
        !and normalising it by the vector's length. The average slope of the two resulting slopes is the tangent.
        do i = 1,nHull
            !get the vector components
            xcomp(1) = Hull(1,i+1) - Hull(1,i)
            xcomp(2) = Hull(1,i) - Hull(1,i-1)
            ycomp(1) = Hull(2,i+1) - Hull(2,i)
            ycomp(2) = Hull(2,i) - Hull(2,i-1)
            
            !normalise by the vector lengths
            vlen = sqrt(xcomp**2+ycomp**2)
            xcomp = xcomp/vlen
            ycomp = ycomp/vlen
            
            !get tangent vector
            xcomp(1) = xcomp(1)+xcomp(2)
            ycomp(1) = ycomp(1)+ycomp(2)
            
            !Get new extrapolated points by expanding the convex hull by tangential vectors and keeping the heights constant
            expandfactor = min(abs(xcomp(1)),abs(ycomp(1)))
            if(expandfactor == 0) expandfactor = 1
            node_xy2(1,node_num+i) = Hull(1,i) + 2*nxew*ycomp(1)/expandfactor
            node_xy2(2,node_num+i) = Hull(2,i) - 2*nyew*xcomp(1)/expandfactor
            zd2(:,node_num+i) = zHull(:,i)      !copy heights directly

  
        end do
        
        !node_num2 = node_num            !otherwise use original data
        !allocate(node_xy2(2,node_num2),zd2(node_num2))
        !node_xy2 = node_xy
        !zd2 = zd
    else !there should only be 3 points here. Generate plane.
        node_num2 = node_num            !otherwise use original data
        allocate(node_xy2(2,node_num2),zd2(nlayer-1,node_num2))
        node_xy2 = node_xy
        zd2 = zd
    end if
  else
    node_num2 = node_num            !otherwise use original data
    allocate(node_xy2(2,node_num2),zd2(nlayer-1,node_num2))
    node_xy2 = node_xy
    zd2 = zd
      
  end if
  

    
    allocate(element_neighbor(3,2*node_num2),triangle(element_order,2*node_num2))
    
    !open(20457,file='points.txt')
    !do i =1,node_num2
    !    write(20457,*), node_xy2(1,i),node_xy2(2,i),zd2(1,i)
    !    
    !end do
    !close(20457)
  
  !  Set up the Delaunay triangulation if there are sufficient points (otherwise not needed).
  if (node_num2 > 3) then

    call r8tris2 ( node_num2, node_xy2, element_num, triangle, element_neighbor )
  end if

if (.not. present(plocation) .and. .not. present(lpdepths)) then
	!Interpolate each layer boundary, other than the top and bottom, which is equal to 1 and nzew respectively.
	 do i = 1,nlayer - 1
		 call interp_surface(dim_num,node_num2,nxew*nyew,node_xy2,zd2(i,:),xyi,zi(i,:,:), element_num, triangle, element_neighbor)
	 end do
 
 else !just export the layer depths at each pile location
 	 do i = 1,nlayer - 1

		 	call interp_points(dim_num,node_num2,size(plocation,2),node_xy2,zd2(i,:),dble(transpose(plocation)),lpdepths(i,:), element_num, triangle, element_neighbor)

	 end do
 
 
 end if
 
 
             if (.false.) then !save the various boundary components to disk for checking
                open(999,file='bfldsi.txt')
                do i = 1,nyew
                    write(999,'(10000000(F8.4,X))')  zi(1,:,i)
                end do
                close(999)
                open(999,file='si_heights.txt')
                write(999,'(10000000(F12.4,X))') zd2(1,:)
                write(999,'(10000000(F12.4,X))') node_xy2(1,:)
                write(999,'(10000000(F12.4,X))') node_xy2(2,:)
                close(999)
                stop
             end if
             

 
 deallocate(triangle,element_neighbor)


end subroutine

   
    
    
!interpolate a single plane 
subroutine interp_surface(dim_num,node_num,ni,node_xy,zd,xyi,zi, element_num, triangle, element_neighbor)

!*****************************************************************************80
!
!! TEST02 tests PWL_INTERP_2D_SCATTERED_VALUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), intent(in) :: dim_num
  integer ( kind = 4 ), intent(in) :: ni !number of points to interpolate
  integer ( kind = 4 ), intent(in) :: node_num !number of input coordinates to interpolate from
  real ( kind = 8 ), intent(in) :: node_xy(:,:) !input x,y coordinates with data
  real ( kind = 8 ), intent(in) :: zd(:) ! zd(node_num) !input z data matching the coordinates
  
  !extra variables in case points are colinear and need to be duplicated with an offset
  real(8), allocatable :: node_xy2(:,:),zd2(:)

  real ( kind = 8 ), intent(in) :: xyi(:,:) !xyi(dim_num,ni) !coordinates of points to interpolate
  real ( kind = 8 ), intent(out) :: zi(:,:) !z data of points to interpolate
  !real(8) :: zi2(ni)
  
  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ), intent(in) :: element_neighbor(:,:) ! element_neighbor(3,2*node_num)
  integer ( kind = 4 ), intent(in) :: element_num
  integer ( kind = 4 ), intent(in) :: triangle(:,:) !triangle(element_order,2*node_num)
  
  integer x,y,i,counter !loop counters
  real point3(2)
  integer perp_vector(2) !vector at 90 degrees to the plan view representation of the vector connecting 2 points in the 2 borehole case (helps create a third point to define a plane
  integer vector1(3,2),vector2(3,2) !x,y,z vectors of 3 points on a plane
  real cross(3) !coefficients for plane equation:  z = -x(a/c) -y(b/c) + d/c
  
  !Copies of the input data. Ensure at least 3 points are taken when building the plane
  real(8) z2(3),xy2(dim_num,3) 

  

  

!
!  Some special cases for when there are 3 or less points to interpolate (since triangulation won't work in these cases)
!
  if (node_num == 1) then
      zi = zd(1)
  else if (node_num == 2) then
      if (zd(2) == zd(1)) then !both points have the same depth; set boundary at this constant depth
          zi = zd(1) 
      else  !get equation of the plan from the two points, plus a third point which is perpenicular to the [point1, point2] vector, but at the same height as point1.
          z2(:node_num) = zd !nint(zd)    !fit points into bigger array
          xy2(:,:node_num) = node_xy    !fit points into bigger array
          !calculate third point from perpendicular vector
          perp_vector(2) = node_xy(1,2) - node_xy(1,1)      !get perpependicular vector         
          perp_vector(1) = -(node_xy(2,2) - node_xy(2,1))
          xy2(1,3) = xy2(1,1) + perp_vector(1)      !get new point by translating first point by vector
          xy2(2,3) = xy2(2,1) + perp_vector(2)
          z2(3) = z2(1) !height of new point is same as old point
          

          call getplaneeq(xy2,z2,cross,dim_num,3)
          counter = 0
          do y = 1,size(zi,2)
              do x=1,size(zi,1)    
                    counter = counter + 1
                    zi(x,y) = cross(1)*xyi(1,counter) + cross(2)*xyi(2,counter) + cross(3)  !interpolate the layer based on the plane equation
              end do
          end do
      end if
  else if (node_num == 3) then
     if (zd(2) == zd(1) .and. zd(3) == zd(2)) then !all points have the same depth; set boundary at this constant depth
          zi = zd(1) 
     else
          call getplaneeq(node_xy,zd,cross,dim_num,3)
          counter = 0
          do y = 1,size(zi,2)
              do x=1,size(zi,1)    
                    counter = counter + 1
                    zi(x,y) = cross(1)*xyi(1,counter) + cross(2)*xyi(2,counter) + cross(3)  !interpolate the layer based on the plane equation
              end do
          end do
     end if
  
  else !There are more than 3 points to work with: do piecewise linear interpolation
  

!
!  Evaluate the interpolant using Delaunay trigangulation
!

  call pwl_interp_2d_scattered_value ( node_num, node_xy, zd, element_num, triangle(:,:element_num), element_neighbor(:,:element_num), ni, xyi, zi )
  !call planint(zi,160,160,element_num,triangle,node_num,sngl(node_xy),sngl(zd)) !,zi2)
  
  

  
  end if

end subroutine




   
!interpolate individual points on the plane
subroutine interp_points(dim_num,node_num,ni,node_xy,zd,xyi,zi, element_num, triangle, element_neighbor)

!*****************************************************************************80
!
!! TEST02 tests PWL_INTERP_2D_SCATTERED_VALUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), intent(in) :: dim_num
  integer ( kind = 4 ), intent(in) :: ni !number of points to interpolate
  integer ( kind = 4 ), intent(in) :: node_num !number of input coordinates to interpolate from
  real ( kind = 8 ), intent(in) :: node_xy(:,:) !input x,y coordinates with data
  real ( kind = 8 ), intent(in) :: zd(:) ! zd(node_num) !input z data matching the coordinates
  
  !extra variables in case points are colinear and need to be duplicated with an offset
  real(8), allocatable :: node_xy2(:,:),zd2(:)


  real ( kind = 8 ), intent(in) :: xyi(:,:) !xyi(dim_num,ni) !coordinates of points to interpolate
  real ( kind = 8 ), intent(out) :: zi(:) !zi(ni) !z data of points to interpolate
  !real(8) :: zi2(ni)
  
  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ), intent(in) :: element_neighbor(:,:) ! element_neighbor(3,2*node_num)
  integer ( kind = 4 ), intent(in) :: element_num
  integer ( kind = 4 ), intent(in) :: triangle(:,:) !triangle(element_order,2*node_num)
  
  integer x,y,i,counter !loop counters
  real point3(2)
  integer perp_vector(2) !vector at 90 degrees to the plan view representation of the vector connecting 2 points in the 2 borehole case (helps create a third point to define a plane
  integer vector1(3,2),vector2(3,2) !x,y,z vectors of 3 points on a plane
  real cross(3) !coefficients for plane equation:  z = -x(a/c) -y(b/c) + d/c
  
  !Copies of the input data. Ensure at least 3 points are taken when building the plane
  real(8) z2(3),xy2(dim_num,3) 

  

!  Some special cases for when there are 3 or less points to interpolate (since triangulation won't work in these cases)
!
  if (node_num == 1) then
      zi = zd(1)
  else if (node_num == 2) then
      if (zd(2) == zd(1)) then !both points have the same depth; set boundary at this constant depth
          zi = zd(1) 
      else  !get equation of the plan from the two points, plus a third point which is perpenicular to the [point1, point2] vector, but at the same height as point1.
          z2(:node_num) = zd    !fit points into bigger array
          xy2(:,:node_num) = node_xy    !fit points into bigger array
          !calculate third point from perpendicular vector
          perp_vector(2) = node_xy(1,2) - node_xy(1,1)      !get perpependicular vector         
          perp_vector(1) = -(node_xy(2,2) - node_xy(2,1))
          xy2(1,3) = xy2(1,1) + perp_vector(1)      !get new point by translating first point by vector
          xy2(2,3) = xy2(2,1) + perp_vector(2)
          z2(3) = z2(1) !height of new point is same as old point

          call getplaneeq(xy2,z2,cross,dim_num,3)
          do x=1,ni
            zi(x) = cross(1)*xyi(1,x) + cross(2)*xyi(2,x) + cross(3)  !interpolate the layer based on the plane equation
          end do
      end if
  else if (node_num == 3) then
     if (zd(2) == zd(1) .and. zd(3) == zd(2)) then !all points have the same depth; set boundary at this constant depth
          zi = zd(1) 
     else
          call getplaneeq(node_xy,zd,cross,dim_num,3)
          do x=1,ni
            zi(x) = cross(1)*xyi(1,x) + cross(2)*xyi(2,x) + cross(3)  !interpolate the layer based on the plane equation
          end do
     end if
  
  else !There are more than 3 points to work with: do piecewise linear interpolation
  

!
!  Evaluate the interpolant using Delaunay trigangulation
!

  call pwl_interp_2d_scattered_value ( node_num, node_xy, zd, element_num, triangle(:,:element_num), element_neighbor(:,:element_num), ni, xyi, zi )
  !call planint(zi,160,160,element_num,triangle,node_num,sngl(node_xy),sngl(zd)) !,zi2)
  
  

  
  end if

end subroutine




!subroutines and functions for creating a convex hull

FUNCTION Cross(v1,v2,v3)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(2)    !< input vector 1
REAL,INTENT(IN) :: v2(2)    !< input vector 2
REAL(8),INTENT(IN) :: v3(2)    !< input vector 3
!-----------------------------------------------
! OUTPUT VARIABLES
REAL            :: Cross    !< cross product
!-----------------------------------------------
! LOCAL VARIABLES
!===============================================
Cross=(v2(1)-v1(1))*(v3(2)-v1(2))-(v2(2)-v1(2))*(v3(1)-v1(1))
END FUNCTION Cross

SUBROUTINE ConvHull(nPoints,Points,zpoints,nlayer,nHull,Hull,zHull)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: nPoints
REAL(8),INTENT(IN)     :: Points(2,0:nPoints-1)
real(8), intent(in)    :: zpoints(nlayer-1,0:nPoints-1)
integer, intent(in) :: nlayer
!------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT) :: nHull
! NOTE: allocate Hull always one point greater than Points, because we save the first value twice
REAL,INTENT(OUT)    :: Hull(2,0:nPoints+1)
REAL,INTENT(OUT)    :: zHull(nlayer-1,0:nPoints+1)
!------------------------------------------------
! LOCAL VARIABLES
REAL                :: Lower(2,0:nPoints-1)
REAL                :: Upper(2,0:nPoints-1)
real                :: zLower(nlayer-1,0:nPoints-1),zUpper(nlayer-1,0:nPoints-1) !z coordinates
INTEGER             :: i,iLower,iUpper
!================================================
IF(nPoints.LE.1)THEN
  Hull  = Points
  nHull = nPoints
ELSE
  iLower = 0
  Lower  = -HUGE(1.)
  DO i=0,nPoints-1
    DO WHILE(iLower.GE.2)
      if(Cross(Lower(:,iLower-2),Lower(:,iLower-1),Points(:,i))>=0.) exit !originally > 0. >= 0 ensures that points on a straight line are included
      Lower(:,iLower) = -HUGE(1.)
      iLower          = iLower - 1
    END DO
    Lower(:,iLower) = Points(:,i)
    zLower(:,iLower) = zpoints(:,i)
    iLower = iLower + 1
  END DO

  iUpper = 0
  Upper  = HUGE(1.)
  DO i=nPoints-1,0,-1
    DO WHILE(iUpper.GE.2)
      if(Cross(Upper(:,iUpper-2),Upper(:,iUpper-1),Points(:,i))>=0.) exit !originally > 0. >= 0 ensures that points on a straight line are included
      Upper(:,iUpper) = HUGE(1.)
      iUpper          = iUpper - 1
    END DO
    Upper(:,iUpper) = Points(:,i)
    zUpper(:,iUpper) = zpoints(:,i)
    iUpper = iUpper + 1
  END DO

  iLower = iLower-1
  iUpper = iUpper-1
  nHull  = iLower+iUpper !+1 (eliminate this plus 1. We don't want to close the polygon, we only want a list of points defining it).
  
  ! NOTE: Initialize Hull with zeros
  Hull   = 0.

  ! NOTE: save values in Hull
  Hull(:,0     :iLower       -1) = Lower(:,0:iLower-1)  !x,y coords
  Hull(:,iLower:iLower+iUpper-1) = Upper(:,0:iUpper-1)
  
  zHull(:,0     :iLower       -1) = zLower(:,0:iLower-1) !layer height information
  zHull(:,iLower:iLower+iUpper-1) = zUpper(:,0:iUpper-1)

  ! NOTE: save first value twice, and the 2nd value twice (to help with tangent calculations later
  Hull(:,       iLower+iUpper  ) = Hull (:,0         )
  Hull(:,       iLower+iUpper +1 ) = Hull (:,1         )
  
  zHull(:,       iLower+iUpper  ) = zHull (:,0         )
  zHull(:,       iLower+iUpper +1 ) = zHull (:,1         )
  
  
  
    ! NOTE: save values in Hull
  Hull(:,1     :iLower       ) = Lower(:,0:iLower-1)  !x,y coords
  Hull(:,iLower+1:iLower+iUpper) = Upper(:,0:iUpper-1)
  
  zHull(:,1     :iLower       ) = zLower(:,0:iLower-1) !layer height information
  zHull(:,iLower+1:iLower+iUpper) = zUpper(:,0:iUpper-1)

  ! NOTE: save first value twice, and the 2nd value twice (to help with tangent calculations later
  Hull(:,       0  ) = Hull (:,iLower+iUpper       )
  Hull(:,       iLower+iUpper +1 ) = Hull (:,1         )
  
  zHull(:,       0  ) = zHull (:,iLower+iUpper         )
  zHull(:,       iLower+iUpper +1 ) = zHull (:,1         )
  
  
END IF

END SUBROUTINE


subroutine planint(bfld,nxew,nyew,num_triangles,triangle,node_num,node_xy,zd) !,bfldo)

!This a custom subroutine to interpolate a triangulated plane. It works by looping through the triangles
!and loops from the top of the triangle downwards, from the left to the right. 
!It builds the left and right sides by sorting the points first such that points 1 and 3 are at the top and bottom respectively,
!while tending to be at the left and right respectively should point 2 be at the same level as either.

implicit none

integer, intent(in) :: triangle(3,2*node_num)
real(4),intent(in) :: node_xy(2,node_num)
real(4),intent(in) :: zd(node_num)
integer, intent(in) :: num_triangles
real(8), intent(out) :: bfld(nxew,nyew)
!real(8), intent(in) :: bfldo(nxew,nyew)
real bfld2(nxew,nyew)
integer, intent(in) :: nxew,nyew,node_num

integer left(nyew),right(nyew),tempside(nyew)
integer trinode(2,3), triz(3)
integer tempvec(2)
real tempval
real cross(3) !parameters for plane equation
integer position

integer j,k,tri,x,y,i,p
integer ystart,yend

real start, finish


        
    do tri = 1,num_triangles
        
            !extract the 3 points of the current triangle
        !call CPU_TIME(start)
        !do p=1,10000

            do j = 1,3 
                triz(j) = zd(triangle(j,tri))
                do i=1,2
                 trinode(i,j) = nint(node_xy(i,triangle(j,tri))) !node_xy should theoretically already be integers from earlier in the program
               end do
            end do
            
        ! --- order the points in increasing y coordinates then x coordinates
        do k = 3,2,-1	 
			do j = 2,4-k,-1
				if (trinode(2,j+1) < trinode(2,j)) then 
					tempvec = trinode(:,j)
					trinode(:,j) = trinode(:,j+1)
					trinode(:,j+1) = tempvec
                    tempval = triz(j)
                    triz(j) = triz(j+1)
                    triz(j+1) = tempval
				end if
			end do
        end do   
        
        !secondary sort in increasing order for y coordinate
        do k = 3,2,-1	 
			do j = 2,4-k,-1
				if (trinode(2,j+1) <= trinode(2,j) .and. trinode(1,j+1) < trinode(1,j)) then 
					tempvec = trinode(:,j)
					trinode(:,j) = trinode(:,j+1)
					trinode(:,j+1) = tempvec
                    tempval = triz(j)
                    triz(j) = triz(j+1)
                    triz(j+1) = tempval
				end if
			end do
        end do   
        
        
        
        position = (trinode(1,3)-trinode(1,1))*(trinode(2,2)-trinode(2,1)) - (trinode(2,3)-trinode(2,1))*(trinode(1,2)-trinode(1,1))
        

        
        !linearly interpolate x coords
        ystart = max(1,trinode(2,1))
        yend = min(nyew,trinode(2,3))
        do i = ystart,yend
            right(i) = trinode(1,1) +  (i-trinode(2,1)) * (trinode(1,3)-trinode(1,1))/(trinode(2,3)-trinode(2,1))
        end do
        
        
        if( trinode(2,2) > trinode(2,1)) then
            ystart = max(1,trinode(2,1))
            yend = min(nyew,trinode(2,2))
            do i = ystart,yend
                left(i) = trinode(1,1) +  (i-trinode(2,1)) * (trinode(1,2)-trinode(1,1))/(trinode(2,2)-trinode(2,1))
            end do
        end if
            
        if( trinode(2,2) < trinode(2,3)) then
            ystart = max(1,trinode(2,2))
            yend = min(nyew,trinode(2,3))
            do i = ystart,yend
                left(i) = trinode(1,2) +  (i-trinode(2,2)) * (trinode(1,3)-trinode(1,2))/(trinode(2,3)-trinode(2,2))
            end do
        end if
            
        !The above code assumes point 2 lies on one side of the others. Flip the left and right vectors if this isn't the case.
        if (position<0) then
            tempside = right
            right = left
            left = tempside
        end if
        
        
        ystart = max(1,minval(trinode(2,:)))
        yend = min(nyew,maxval(trinode(2,:)))
        
        do i = ystart,yend
            if(right(i) < 1) right(i) = 1
            if(right(i) > nxew) right(i) = nxew
            if(left(i) < 1) left(i) = 1
            if(left(i) > nxew) left(i) = nxew
        end do
        
        
        !end do
        !call CPU_TIME(finish)
        !write(*,*) (finish-start)/10000
        
        
        !call CPU_TIME(start)
        !do p=1,10000
        

        !get equation for plane
        call getplaneeq(dble(trinode),dble(triz),cross,2,3)
        do x = ystart,yend
          do y=left(x),right(x)
            bfld2(x,y) = cross(1)*y + cross(2)*x + cross(3)  !interpolate the layer based on the plane equation
            if(bfld2(x,y) > 100) then
                write(*,*)
            !if(abs(bfld2(x,y)-bfldo(x,y)) > 0.000001) then
            !    write(*,*)
            end if
          end do
        end do
        
        !end do
        !call CPU_TIME(finish)
        !write(*,*) (finish-start)/10000

        
        
    end do
            
    bfld = bfld2
            



end subroutine



!---------------------------------------------------------------------------------------------------------------------------

! The following subroutines are to do with pre-processing the delauney triangulation and interpolating individual points as opposed to a grid

!---------------------------------------------------------------------------------------------------------------------------


 
    
!This subroutine calculates the delauney triangulation
!Delauney triangulation is not an inexpensive procedure, and it only needs to be done once for a set of boreholes.
!This subroutine also does other checks, for example to see if the points are colinear, in which case triangulation will fail.
subroutine interp_layers_prep(dim_num,node_num,nxew,nyew,node_xy_in,node_xy_out,zd_out,element_neighbor_out,triangle_out,element_num,extranum)

!Note, this version is intended to pre-process the delauney triangulation for later use.

implicit none


  integer ( kind = 4 ), intent(in) :: dim_num
  integer ( kind = 4 ), intent(in) :: nxew,nyew !dimensions of the site
  integer ( kind = 4 ), intent(in) :: node_num !number of input coordinates to interpolate from
  integer, intent(in) :: node_xy_in(:,:) !node_xy_in(node_num,dim_num)
  real, intent(out) :: node_xy_out(:,:) !node_xy_out(node_num,dim_num)
  real ( kind = 8 ) :: node_xy(dim_num,node_num) !input x,y coordinates with data
  real ( kind = 8 )  :: node_xy_copy(dim_num,node_num) !input x,y coordinates with data copy
  integer, intent(out) :: zd_out(:) !zd_out(node_num) !input z data matching the coordinates
  real ( kind = 8 ) :: zd(node_num) !input z data matching the coordinates
  real(8) tempval
  integer, intent(out) :: extranum !number of extrapolated points
  
  !extra variables in case points are colinear and need to be duplicated with an offset
  real(8), allocatable :: node_xy2(:,:),zd2(:)
  integer node_num2

  
  integer i,j,k
  
  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ), allocatable :: element_neighbor(:,:)
  integer ( kind = 4 ) element_num
  integer ( kind = 4 ), allocatable :: triangle(:,:)
  
  
  integer, intent(out) :: element_neighbor_out(:,:)
  integer, intent(out) :: triangle_out(:,:)
  
  real rise(node_num-1),run(node_num-1),slopes(node_num-1)
  real expandfactor !factor by which to expand the hull
  
  !logical, intent(in) :: getsurface
  
  !convex hull stuff
  REAL :: Hull(2,0:node_num+1)          !x,y coords of the convex hull
  REAL :: zHull(0:node_num+1)   !layer height information of the convex hull
  real temp(2)
  real xcomp(2),ycomp(2),vlen(2) !x,y components and length of 2 adjacent vectors from a point on the convex hull
  integer nHull
  
  !transpose input coordinates and store as the required double precision
  node_xy = transpose(node_xy_in)
  extranum = 0
  
  
    !Check to see if the points are colinear
    !If they are, then create a duplicate set of points and put an offset between them
    !at 90 degrees to the line, forming a series of rectangles to triangulate.
  if (node_num > 2) then
    do i = 1,node_num-1
          rise(i) = node_xy(2,1) - node_xy(2,i+1)
          run(i) = node_xy(1,1) - node_xy(1,i+1)
    end do
    if (rise(1) == 0) then !calculate slopes, try to ensure no divide by zero
        slopes = rise/run
    else
        slopes = run/rise
    end if
    if(all(abs(slopes-slopes(1)) <= 0.001)) then !points are indeed colinear, or close to it
        node_num2 = node_num*2  !solve the problem by duplicating the points and providing offsets at 90 degrees to the slope
        allocate(node_xy2(2,node_num2),zd2(node_num2))
        expandfactor = min(abs(run(1)),abs(rise(1)))
        if(expandfactor == 0) expandfactor = 1
        node_xy2(1,:node_num) = node_xy(1,:) - 2*nxew*rise(1)/expandfactor  !offset one way
        node_xy2(2,:node_num) = node_xy(2,:) + 2*nyew*run(1)/expandfactor
        node_xy2(1,node_num+1:) = node_xy(1,:) + 2*nxew*rise(1)/expandfactor !and the other way
        node_xy2(2,node_num+1:) = node_xy(2,:) - 2*nyew*run(1)/expandfactor
    else if(node_num > 3) then !if points are not colinear, and there are more than 3 (since 3 points are reserved to define a plane) then do regular interpolation
        !Extrapolate by creating and expanding a convex hull. Seeds sorted coordinates
        
        zd = [ (j, j = 1,node_num) ] !create list of array indices corresponding to original coordinates
        
        
        node_xy_copy = node_xy
        
        !sort in increasing order for x coordinate
        do k = node_num,2,-1	 
			do j = node_num-1,node_num-k+1,-1
				if (node_xy(1,j+1) < node_xy(1,j)) then 
					temp = node_xy(:,j)
					node_xy(:,j) = node_xy(:,j+1)
					node_xy(:,j+1) = temp
                    tempval = zd(j)
                    zd(j) = zd(j+1)
                    zd(j+1) = tempval
				end if
			end do
        end do   
        
        !secondary sort in increasing order for y coordinate
        do k = node_num,2,-1	 
			do j = node_num-1,node_num-k+1,-1
				if (node_xy(1,j+1) <= node_xy(1,j) .and. node_xy(2,j+1) < node_xy(2,j)) then 
					temp = node_xy(:,j)
					node_xy(:,j) = node_xy(:,j+1)
					node_xy(:,j+1) = temp
                    tempval = zd(j)
                    zd(j) = zd(j+1)
                    zd(j+1) = tempval
				end if
			end do
        end do   
        
        !zd = [ (j = 1,node_num) ] !create list of array indices corresponding to sorted coordinates
        
        call ConvHull1z(node_num,node_xy,zd,nHull,Hull,zHull)
        
        node_num2 = node_num + nHull    
        allocate(node_xy2(2,node_num2),zd2(node_num2))
        node_xy2(:,:node_num) = node_xy_copy !use original order of data now that the convex hull is finished
        zd2(:node_num) = zd
        
        !find tangents on the convex hull by finding the vectors of the points either side of the point of interest
        !and normalising it by the vector's length. The average slope of the two resulting slopes is the tangent.
        do i = 1,nHull
            !get the vector components
            xcomp(1) = Hull(1,i+1) - Hull(1,i)
            xcomp(2) = Hull(1,i) - Hull(1,i-1)
            ycomp(1) = Hull(2,i+1) - Hull(2,i)
            ycomp(2) = Hull(2,i) - Hull(2,i-1)
            
            !normalise by the vector lengths
            vlen = sqrt(xcomp**2+ycomp**2)
            xcomp = xcomp/vlen
            ycomp = ycomp/vlen
            
            !get tangent vector
            xcomp(1) = xcomp(1)+xcomp(2)
            ycomp(1) = ycomp(1)+ycomp(2)
            
            !Get new extrapolated points by expanding the convex hull by tangential vectors and keeping the heights constant
            expandfactor = min(abs(xcomp(1)),abs(ycomp(1)))
            if(expandfactor == 0) expandfactor = 1
            node_xy2(1,node_num+i) = Hull(1,i) + 2*nxew*ycomp(1)/expandfactor
            node_xy2(2,node_num+i) = Hull(2,i) - 2*nyew*xcomp(1)/expandfactor
            zd2(node_num+i) = zHull(i)      !copy heights directly

  
        end do
        
        
        !node_num2 = node_num            !otherwise use original data
        !allocate(node_xy2(2,node_num2),zd2(node_num2))
        !node_xy2 = node_xy
        !zd2 = zd

    else !there should only be 3 points here. Generate plane.
        node_num2 = node_num            !otherwise use original data
        allocate(node_xy2(2,node_num2),zd2(node_num2))
        node_xy2 = node_xy
        zd2 = zd
    end if
  else
      
    node_num2 = node_num            !otherwise use original data
    allocate(node_xy2(2,node_num2),zd2(node_num2))
    node_xy2 = node_xy
    zd2 = zd
      
  end if
  

    
    
    

  
  !  Set up the Delaunay triangulation if there are sufficient points (otherwise not needed).
  if (node_num2 > 3) then
    allocate(element_neighbor(3,2*node_num2),triangle(element_order,2*node_num2))
    call r8tris2 ( node_num2, node_xy2, element_num, triangle, element_neighbor )
    
    
      !--- export final parameters ---
      extranum = node_num2-node_num
      node_xy_out(:extranum,:) = transpose(node_xy2(:,node_num+1:))  !save x/y coords of the convex hull
      zd_out(:extranum) = nint(zd2(node_num+1:)) !save the z index of the original data for the convex hull
  
      !store the triangulation 
      element_neighbor_out(:,:2*node_num2) = element_neighbor
      triangle_out(:,:2*node_num2) = triangle
    
     deallocate(triangle,element_neighbor)
    
  end if
  



!if (.not. present(plocation) .and. .not. present(lpdepths)) then
!	!Interpolate each layer boundary, other than the top and bottom, which is equal to 1 and nzew respectively.
!	 do i = 1,nlayer - 1
!		 call interp_surface(dim_num,node_num2,nxew*nyew,node_xy2,zd2(i,:),xyi,zi(i,:,:), element_num, triangle, element_neighbor)
!	 end do
! 
! else !just export the layer depths at each pile location
! 	 do i = 1,nlayer - 1
!
!		 	call interp_surface(dim_num,node_num2,size(plocation,1),node_xy2,zd2(i,:),dble(transpose(plocation)),lpdepths(i,:), element_num, triangle, element_neighbor)
!
!	 end do
! 
! 
! end if
 


 



end subroutine




!This subroutine performs interpolation of a number of layers.
!Delauney triangulation is not an inexpensive procedure, and it only needs to be done once for a set of boreholes.
!This subroutine also does other checks, for example to see if the points are colinear, in which case triangulation will fail.
subroutine interp_layers_proc(dim_num,node_num,nxew,nyew,node_xy,zd,xyi,zi,nlayer,node_xy_out,zd_out,element_neighbor,triangle,element_num,extranum,plocation,lpdepths)

!Note, this version is intended to do the interpolation given that interp_layers_prep has already been called

implicit none


  integer ( kind = 4 ), intent(in) :: dim_num
  integer ( kind = 4 ), intent(in) :: nxew,nyew !dimensions of the site
  integer ( kind = 4 ), intent(in) :: node_num !number of input coordinates to interpolate from
  real ( kind = 8 )  :: node_xy(:,:) ! (dim_num,node_num) !input x,y coordinates with data
  real ( kind = 8 ), intent(in) :: zd(:,:) !(nlayer-1,node_num) !input z data matching the coordinates
  integer, intent(in) :: nlayer !number of layers to interpolate
  
  !extra variables in case points are colinear and need to be duplicated with an offset
  real(8), allocatable :: node_xy2(:,:),zd2(:,:)
  integer node_num2

  real ( kind = 8 ), intent(in) :: xyi(:,:) ! (dim_num,nxew*nyew) !coordinates of points to interpolate
  real ( kind = 8 ), intent(out) :: zi(:,:,:) ! (nlayer-1,nxew,nyew) !z data of points to interpolate
  
  real, intent(in) ::  node_xy_out(:,:) !node_xy_out(node_num,dim_num)
  integer, intent(in) :: zd_out(:) !zd_out(node_num) !input z data matching the coordinates
  integer, intent(in) :: extranum !number of extrapolated points
  
  integer i,j,k
  
  integer ( kind = 4 ), parameter :: element_order = 3
  integer ( kind = 4 ), intent(in) :: element_neighbor(:,:)
  integer ( kind = 4 ), intent(in) :: element_num
  integer ( kind = 4 ), intent(in)  :: triangle(:,:)
  
  real rise(node_num-1),run(node_num-1),slopes(node_num-1)
  real expandfactor !factor by which to expand the hull
  
  
  !point interpolation at pile location stuff
  real(8), optional, intent(in) :: plocation(:,:)
  real(8), optional, intent(out) :: lpdepths(:,:)
  !logical, intent(in) :: getsurface
  
  !convex hull stuff
  REAL :: Hull(2,0:node_num+1)          !x,y coords of the convex hull
  REAL :: zHull(nlayer-1,0:node_num+1)   !layer height information of the convex hull
  real temp(2)
  real xcomp(2),ycomp(2),vlen(2) !x,y components and length of 2 adjacent vectors from a point on the convex hull
  integer nHull
  
  
  if(extranum > 0) then !if there are extrapolated data points, add them in
    node_num2 = node_num + extranum   
    allocate(node_xy2(2,node_num2),zd2(nlayer-1,node_num2))
    node_xy2(:,:node_num) = node_xy 
    zd2(:,:node_num) = zd
        
    node_xy2(:,node_num+1:) = transpose(node_xy_out(:extranum,:))
        
    zd2(:,node_num+1:) = zd(:,zd_out(:extranum))
      
      
  else
    node_num2 = node_num            !otherwise use original data
    allocate(node_xy2(2,node_num2),zd2(nlayer-1,node_num2))
    node_xy2 = node_xy
    zd2 = zd
      
  end if
  


if (.not. present(plocation) .and. .not. present(lpdepths)) then
	!Interpolate each layer boundary, other than the top and bottom, which is equal to 1 and nzew respectively.
	 do i = 1,nlayer - 1
		 call interp_surface(dim_num,node_num2,nxew*nyew,node_xy2,zd2(i,:),xyi,zi(i,:,:), element_num, triangle, element_neighbor)
	 end do
 
 else !just export the layer depths at each pile location
 	 do i = 1,nlayer - 1
		 	call interp_points(dim_num,node_num2,size(plocation,2),node_xy2,zd2(i,:),plocation,lpdepths(i,:), element_num, triangle, element_neighbor)

	 end do
 
 
 end if



end subroutine




!same as the above ConvHull subroutine, except works for a single set of Z values (1 layer) instead of multiple layers.
!Intended to work by sorting the Z points as array indices, to extract the corresponding depth values at a later time
SUBROUTINE ConvHull1z(nPoints,Points,zpoints,nHull,Hull,zHull)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: nPoints
REAL(8),INTENT(IN)     :: Points(2,0:nPoints-1)
real(8), intent(in)    :: zpoints(0:nPoints-1)
!------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT) :: nHull
! NOTE: allocate Hull always one point greater than Points, because we save the first value twice
REAL,INTENT(OUT)    :: Hull(2,0:nPoints+1)
REAL,INTENT(OUT)    :: zHull(0:nPoints+1)
!------------------------------------------------
! LOCAL VARIABLES
REAL                :: Lower(2,0:nPoints-1)
REAL                :: Upper(2,0:nPoints-1)
real                :: zLower(0:nPoints-1),zUpper(0:nPoints-1) !z coordinates
INTEGER             :: i,iLower,iUpper
!================================================
IF(nPoints.LE.1)THEN
  Hull  = Points
  nHull = nPoints
ELSE
  iLower = 0
  Lower  = -HUGE(1.)
  DO i=0,nPoints-1
    DO WHILE(iLower.GE.2)
      if(Cross(Lower(:,iLower-2),Lower(:,iLower-1),Points(:,i))>=0.) exit !originally > 0. >= 0 ensures that points on a straight line are included
      Lower(:,iLower) = -HUGE(1.)
      iLower          = iLower - 1
    END DO
    Lower(:,iLower) = Points(:,i)
    zLower(iLower) = zpoints(i)
    iLower = iLower + 1
  END DO

  iUpper = 0
  Upper  = HUGE(1.)
  DO i=nPoints-1,0,-1
    DO WHILE(iUpper.GE.2)
      if(Cross(Upper(:,iUpper-2),Upper(:,iUpper-1),Points(:,i))>=0.) exit !originally > 0. >= 0 ensures that points on a straight line are included
      Upper(:,iUpper) = HUGE(1.)
      iUpper          = iUpper - 1
    END DO
    Upper(:,iUpper) = Points(:,i)
    zUpper(iUpper) = zpoints(i)
    iUpper = iUpper + 1
  END DO

  iLower = iLower-1
  iUpper = iUpper-1
  nHull  = iLower+iUpper !+1 (eliminate this plus 1. We don't want to close the polygon, we only want a list of points defining it).
  
  ! NOTE: Initialize Hull with zeros
  Hull   = 0.

  ! NOTE: save values in Hull
  Hull(:,0     :iLower       -1) = Lower(:,0:iLower-1)  !x,y coords
  Hull(:,iLower:iLower+iUpper-1) = Upper(:,0:iUpper-1)
  
  zHull(0     :iLower       -1) = zLower(0:iLower-1) !layer height information
  zHull(iLower:iLower+iUpper-1) = zUpper(0:iUpper-1)

  ! NOTE: save first value twice, and the 2nd value twice (to help with tangent calculations later
  Hull(:,       iLower+iUpper  ) = Hull (:,0         )
  Hull(:,       iLower+iUpper +1 ) = Hull (:,1         )
  
  zHull(iLower+iUpper  ) = zHull (0         )
  zHull(iLower+iUpper +1 ) = zHull (1         )
  
  
  
    ! NOTE: save values in Hull
  Hull(:,1     :iLower       ) = Lower(:,0:iLower-1)  !x,y coords
  Hull(:,iLower+1:iLower+iUpper) = Upper(:,0:iUpper-1)
  
  zHull(1     :iLower       ) = zLower(0:iLower-1) !layer height information
  zHull(iLower+1:iLower+iUpper) = zUpper(0:iUpper-1)

  ! NOTE: save first value twice, and the 2nd value twice (to help with tangent calculations later
  Hull(:,       0  ) = Hull (:,iLower+iUpper       )
  Hull(:,       iLower+iUpper +1 ) = Hull (:,1         )
  
  zHull(       0  ) = zHull (iLower+iUpper         )
  zHull(       iLower+iUpper +1 ) = zHull (1         )
  
  
END IF

END SUBROUTINE





    end module

    
    
    
    

