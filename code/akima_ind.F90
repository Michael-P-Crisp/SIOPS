

subroutine getdesind(xvals, tvals, n_in, xlen, ylen, tlen, n_out, yval)
!get pile design as array index

!Linearly interpolate data to get pile design or true pile settlements
!This has been adapted from the akima subroutine, so it has the same interface
!Linear interpolation has 2 benefits: there's no pre-processing involved for finding 
! curve coefficients (well, slopes could be pre-calculated but that's pretty minimal), 
! and furthermore the actual interpolation itself is a simpler mathematical operation.
!data must be decreasing

!The downside compared to the akima subroutine is that you need a much higher density of
!data to maintain the same degree of accuracy. This could potentially take up a non-trivial
!amount of RAM (hundreds of megabytes depending on settings), and takes more binary search
!iterations to find the location of the interpolant location, although this latter point 
!shouldn't have much impact due to its excellent log(n) efficiency.

!	If dim xvals > 1, it is assumed the xvalues vary and yvals are constant
!	If dim yvals > 1, it is assumed the yvalues vary and xvals are constant
!	If dim xvals == yvals, unique curves are assumed.
!	If dim tvals > 1, the points to be interpolated are assumed to be unique to
!		the curves and are interpolated accordingly. Otherwise all points are 
!		interpolated for all curves.

! Michael Crisp February 2018

use extrafuncs

  implicit none

  integer ( kind = 4 ), intent(in) :: n_in !number of data points of original curve
  integer ( kind = 4 ), intent(in) :: n_out !number of data points to interpolate 
  integer ( kind = 4 ), intent(in) :: xlen,ylen, tlen !number of data sets to use

  integer ( kind = 4 ) i,j
  !integer ( kind = 4 ) ibcbeg
  !integer ( kind = 4 ) ibcend
  real ( kind = 4 ), intent(in) :: xvals(n_in,xlen)
  real ( kind = 4 ), intent(in) :: tvals(n_out,tlen)
  !real ( kind = 4 ) ybcbeg
  !real ( kind = 4 ) ybcend
  !real ( kind = 4 ) ypp(n_in) !knot information
  !real ( kind = 8 ) yppval
  !real ( kind = 8 ) ypval
  integer ( kind = 4 ), intent(out) :: yval(n_out,max(xlen,ylen))
  logical check(xlen) !have a dedicated vector for checking order of x values
  
  
  real nan
  
    !real, dimension(n), intent(out) :: dydx  ! derivative of y w.r.t. x
    !real, dimension(n, npt), intent(out) :: dydxpt, dydypt  ! derivative of y w.r.t. xpt and ypt
  
  !real delta_x 
  
!    delta_x: half-width of the smoothing interval added in the valley of absolute-value function
!             this allows the derivatives with respect to the data points (dydxpt, dydypt)
!             to also be C1 continuous.
!             set to parameter to 0 to get the original Akima function (but only if
!             you don't need dydxpt, dydypt)


        ! evaluate polynomial (and derivative) 		-cut it out to increase speed. Place in lower code blocks and tweak to match code.
!         dx = (x(i) - xpt(j))
!         y(i) = p0(j) + p1(j)*dx + p2(j)*dx**2 + p3(j)*dx**3
!         dydx(i) = p1(j) + 2.0_dp*p2(j)*dx + 3.0_dp*p3(j)*dx**2
! 
!         do k = 1, npt
!             dydxpt(i, k) = dp0dxpt(j, k) + dp1dxpt(j, k)*dx + dp2dxpt(j, k)*dx**2 + dp3dxpt(j, k)*dx**3
!             if (k == j) then
!                 dydxpt(i, k) = dydxpt(i, k) - dydx(i)
!             end if
!             dydypt(i, k) = dp0dypt(j, k) + dp1dypt(j, k)*dx + dp2dypt(j, k)*dx**2 + dp3dypt(j, k)*dx**3
!         end do


!create nan 



!  Check input data is suitable (commented out for speed purposes, since these issues shouldn't ever happen in theory)

	!Need at least 2 points
  !if ( n_in <= 1 ) then
  !  write ( *, '(a)' ) ' Error: The number of points in curve must be at least 2.'
  !  write ( *, '(a,i8)' ) '  The input value of N = ', n_in
  !  read(*,*)
  !  stop 
  !end if
  !
  !!X values must be decreasing.
  !do i = 1, n_in - 1
  !  check = xvals(i+1,:) >= xvals(i,:)
  !  if ( any(check) ) then
  !    write ( *, '(a)' ) 'Error: Curve x-values must be strictly increasing, but'
  !    write ( *, '(a,i8,a,g14.6)' ) '  X(',  i,') = ', xvals(i,:)
  !    write ( *, '(a,i8,a,g14.6)' ) '  X(',i+1,') = ', xvals(i+1,:)
  !    read(*,*)
  !    stop 1
  !  end if
  !end do

  

	if (xlen == ylen) then
	
		if (tlen == 1)  then	
			do j = 1,xlen
				call r8vec_bracket_findright ( xvals(:,j), tvals(:,1), yval(:,j))						!get intervals of points in curve
				!yval(:,j) = yvals(left,j) + dt*(b(left) + dt*(c(left) + dt*d(left)))/error					!interpolate values
			end do
			
		else if (tlen == xlen)  then
		
			do j = 1,xlen
				call r8vec_bracket_findright ( xvals(:,j), tvals(:,j), yval(:,j))								
				!yval(:,j) = yvals(left,j) + dt*(b(left) + dt*(c(left) + dt*d(left)))/error
			end do
			
		else 
			write(*,*) "Error, there is mismatch between size of curve points and"
			write(*,*) "points to be interpolated. Check 2nd dimensions."
            read(*,*)
			stop
		end if
		
	else if (xlen == 1 .and. ylen > 1) then
	
	
		if (tlen == 1)  then
	
			call r8vec_bracket_findright ( xvals(:,1), tvals(:,1), yval(:,j)) !This is contant, can be outside the loop
			
		else if (tlen == ylen)  then
		
			do j = 1,ylen

			
				call r8vec_bracket_findright ( xvals(:,1), tvals(:,j), yval(:,j))

			end do
			
		else 
			write(*,*) "Error, there is mismatch between size of curve points and"
			write(*,*) "points to be interpolated. Check 2nd dimensions."
            read(*,*)
			stop
		end if
				
		
	else if (ylen == 1 .and. xlen > 1) then
		
		if (tlen == 1)  then
		
			do j = 1,xlen
									
				call r8vec_bracket_findright ( xvals(:,j), tvals(:,1), yval(:,j))
														
			end do
		
		else if (tlen == xlen)  then
		
			do j = 1,xlen

				call r8vec_bracket_findright ( xvals(:,j), tvals(:,j), yval(:,j))

			end do
			
		else 
			write(*,*) "Error, there is mismatch between size of curve points and"
			write(*,*) "points to be interpolated. Check 2nd dimensions."
            read(*,*)
			stop
		end if
		
		
	else
		write(*,*) "Error: Unmatching sets of x and y curve inputs. Check 2nd dimensions."
        read(*,*)
		stop
	end if
		

end subroutine
