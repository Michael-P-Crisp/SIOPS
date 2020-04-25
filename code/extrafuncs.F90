module extrafuncs

contains



subroutine r8vec_bracket ( x, xval, left, right, dt,h )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!

! Modified to work on a vector of points.
! Also uses binary search instead of linear search.

  implicit none


  integer ( kind = 4 ) i,j,n
  integer ( kind = 4 ), intent(out) ::  left(:)
  integer ( kind = 4 ), intent(out) ::  right(:)
  real ( kind = 4 ) x(:)
  real ( kind = 4 ) xval(:)
  integer mid
  

  real ( kind = 4 ), intent(out) ::  dt(:)
  real ( kind = 4 ), intent(out) ::  h(:)
  
  
  n = size(x)
  
  !Interval spans full bounds by default
  right = n-1
  left = 2
  
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!



	do j = 1,size(xval)
	
	!Check if it's in the outer bounds
	  if (xval(j) < x(2)) then
	     right(j) = 2
	     left(j) = 1
	     cycle
	  else if(xval(j) >= x(n-1)) then
	  	 right(j) = n
	  	 left(j) = n-1
	  	 cycle

	  end if
	  
	!Otherwise, use binary search to find interval

    do while (right(j) > left(j) + 1)

      mid = ( left(j) + right(j) + 1 ) / 2

      if ( x(mid) <= xval(j) ) then
        left(j) = mid
      else
        right(j) = mid
      end if
	
	end do
   end do
  
  dt = xval - x(left)
  h = x(right) - x(left)
  
  !old linear search
!   	  do i = 2, n - 1
! 
! 		if ( xval(j) < x(i) ) then
! 		  right(j) = i
! 		  exit
! 		end if
! 
! 	   end do

  return
end subroutine



subroutine spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

!*****************************************************************************80
!
!! SPLINE_CUBIC_SET computes the second derivatives of a piecewise cubic spline.
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to
!    determine the second derivative data, passing in the data to be
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output,
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) )
!             + C(IVAL) * ( T - T(IVAL) )^2
!             + D(IVAL) * ( T - T(IVAL) )^3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )^2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))^2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      =
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL)
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!
!    Boundary conditions must be applied at the first and last knots.
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points; N must be
!    at least 2.
!
!    Input, real ( kind = 8 ) T(N), the points where data is specified.
!    The values should be distinct, and increasing.
!
!    Input, real ( kind = 8 ) Y(N), the data values to be interpolated.
!
!    Input, integer ( kind = 4 ) IBCBEG, the left boundary condition flag:
!    0: the spline should be a quadratic over the first interval;
!    1: the first derivative at the left endpoint should be YBCBEG;
!    2: the second derivative at the left endpoint should be YBCBEG;
!    3: Not-a-knot: the third derivative is continuous at T(2).
!
!    Input, real ( kind = 8 ) YBCBEG, the left boundary value, if needed.
!
!    Input, integer ( kind = 4 ) IBCEND, the right boundary condition flag:
!    0: the spline should be a quadratic over the last interval;
!    1: the first derivative at the right endpoint should be YBCEND;
!    2: the second derivative at the right endpoint should be YBCEND;
!    3: Not-a-knot: the third derivative is continuous at T(N-1).
!
!    Input, real ( kind = 8 ) YBCEND, the right boundary value, if needed.
!
!    Output, real ( kind = 8 ) YPP(N), the second derivatives of
!    the cubic spline.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) a1(n)
  real ( kind = 4 ) a2(n)
  real ( kind = 4 ) a3(n)
  real ( kind = 4 ) a4(n)
  real ( kind = 4 ) a5(n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) info
  real ( kind = 4 ) t(n)
  real ( kind = 4 ) y(n)
  real ( kind = 4 ) ybcbeg
  real ( kind = 4 ) ybcend
  real ( kind = 4 ) ypp(n)


!
!  Zero out the matrix.
!
  a1 = 0.0D+00
  a2 = 0.0D+00
  a3 = 0.0D+00
  a4 = 0.0D+00
  a5 = 0.0D+00
!
!  Set the first equation.
!
  if ( ibcbeg == 0 ) then
    b(1) = 0.0D+00
    a3(1) =  1.0D+00
    a4(1) = -1.0D+00
  else if ( ibcbeg == 1 ) then
    b(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
    a3(1) = ( t(2) - t(1) ) / 3.0D+00
    a4(1) = ( t(2) - t(1) ) / 6.0D+00
  else if ( ibcbeg == 2 ) then
    b(1) = ybcbeg
    a3(1) = 1.0D+00
    a4(1) = 0.0D+00
  else if ( ibcbeg == 3 ) then
    b(1) = 0.0D+00
    a3(1) = - ( t(3) - t(2) )
    a4(1) =   ( t(3)        - t(1) )
    a5(1) = - (        t(2) - t(1) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCBEG must be 0, 1, 2 or 3.'
    write ( *, '(a,i8)' ) '  The input value is IBCBEG = ', ibcbeg
    read(*,*)
    stop 1
  end if
!
!  Set the intermediate equations.
!
  do i = 2, n - 1
    b(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) &
         - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
    a2(i) = ( t(i+1) - t(i)   ) / 6.0D+00
    a3(i) = ( t(i+1) - t(i-1) ) / 3.0D+00
    a4(i) = ( t(i)   - t(i-1) ) / 6.0D+00
  end do
!
!  Set the last equation.
!
  if ( ibcend == 0 ) then
    b(n) = 0.0D+00
    a2(n) = -1.0D+00
    a3(n) = 1.0D+00
  else if ( ibcend == 1 ) then
    b(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
    a2(n) = ( t(n) - t(n-1) ) / 6.0D+00
    a3(n) = ( t(n) - t(n-1) ) / 3.0D+00
  else if ( ibcend == 2 ) then
    b(n) = ybcend
    a2(n) = 0.0D+00
    a3(n) = 1.0D+00
  else if ( ibcend == 3 ) then
    b(n) = 0.0D+00
    a1(n) = - ( t(n) - t(n-1) )
    a2(n) =   ( t(n)          - t(n-2) )
    a3(n) = - (        t(n-1) - t(n-2) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCEND must be 0, 1, 2 or 3.'
    write ( *, '(a,i8)' ) '  The input value is IBCEND = ', ibcend
    read(*,*)
    stop 1
  end if
!
!  Special case:
!    N = 2, IBCBEG = IBCEND = 0.
!
  if ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then

    ypp(1) = 0.0D+00
    ypp(2) = 0.0D+00
!
!  Solve the linear system.
!
  else

    call penta ( n, a1, a2, a3, a4, a5, b, ypp )

  end if

  return
end subroutine


subroutine penta ( n, a1, a2, a3, a4, a5, b, x )

!*****************************************************************************80
!
!! PENTA solves a pentadiagonal system of linear equations.
!
!  Discussion:
!
!    The matrix A is pentadiagonal.  It is entirely zero, except for
!    the main diagaonal, and the two immediate sub- and super-diagonals.
!
!    The entries of Row I are stored as:
!
!      A(I,I-2) -> A1(I)
!      A(I,I-1) -> A2(I)
!      A(I,I)   -> A3(I)
!      A(I,I+1) -> A4(I)
!      A(I,I-2) -> A5(I)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Cheney, Kincaid,
!    Numerical Mathematics and Computing,
!    1985, pages 233-236.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), A3(N), A4(N), A5(N), the nonzero
!    elements of the matrix.  Note that the data in A2, A3 and A4
!    is overwritten by this routine during the solution process.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 4 ) a1(n)
  real ( kind = 4 ) a2(n)
  real ( kind = 4 ) a3(n)
  real ( kind = 4 ) a4(n)
  real ( kind = 4 ) a5(n)
  real ( kind = 4 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 4 ) x(n)
  real ( kind = 4 ) xmult

  do i = 2, n - 1
    xmult = a2(i) / a3(i-1)
    a3(i) = a3(i) - xmult * a4(i-1)
    a4(i) = a4(i) - xmult * a5(i-1)
    b(i) = b(i) - xmult * b(i-1)
    xmult = a1(i+1) / a3(i-1)
    a2(i+1) = a2(i+1) - xmult * a4(i-1)
    a3(i+1) = a3(i+1) - xmult * a5(i-1)
    b(i+1) = b(i+1) - xmult * b(i-1)
  end do

  xmult = a2(n) / a3(n-1)
  a3(n) = a3(n) - xmult * a4(n-1)
  x(n) = ( b(n) - xmult * b(n-1) ) / a3(n)
  x(n-1) = ( b(n-1) - a4(n-1) * x(n) ) / a3(n-1)
  do i = n - 2, 1, -1
    x(i) = ( b(i) - a4(i) * x(i+1) - a5(i) * x(i+2) ) / a3(i)
  end do

  return
end subroutine


subroutine spline_cubic_val_vec(y, ypp,left,right,dt,h,yval)


!*****************************************************************************80
!
!! SPLINE_CUBIC_VAL evaluates a piecewise cubic spline at a point.
!
!  Discussion:
!
!    SPLINE_CUBIC_SET must have already been called to define the
!    values of YPP.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A
!             + B * ( T - T(IVAL) )
!             + C * ( T - T(IVAL) )^2
!             + D * ( T - T(IVAL) )^3
!
!    Here:
!      A = Y(IVAL)
!      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C = YPP(IVAL) / 2
!      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) T(N), the knot values.
!
!    Input, real ( kind = 8 ) Y(N), the data values at the knots.
!
!    Input, real ( kind = 8 ) YPP(N), the second derivatives of the
!    spline at the knots.
!
!    Input, real ( kind = 8 ) TVAL, a point, typically between T(1) and
!    T(N), at which the spline is to be evalulated.  If TVAL lies outside
!    this range, extrapolation is used.
!
!    Output, real ( kind = 8 ) YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.

! NOTE: This has been modified so as to not calculate any derivatives, which are not needed.
! The calculation has also been modified to be vectorized

!
  implicit none

  real ( kind = 4 ) y(:)
  real ( kind = 4 ) ypp(:)
  real ( kind = 4 ) dt(:)
  real ( kind = 4 ) h(:)
  integer ( kind = 4 ) left(:)
  integer ( kind = 4 ) right(:)

  !real ( kind = 8 ) yppval
  !real ( kind = 8 ) ypval
  real ( kind = 4 ) yval(:)
  
  yval = y(left) &
       + dt * ( ( y(right) - y(left) ) / h &
              - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h &
       + dt * ( 0.5D+00 * ypp(left) &
       + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0D+00 * h ) ) ) )
       

  !ypval = ( y(right) - y(left) ) / h &
  !     - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h &
  !     + dt * ( ypp(left) &
  !     + dt * ( 0.5D+00 * ( ypp(right) - ypp(left) ) / h ) )

  !yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h


  return
  
end subroutine


!---------- akima subroutines --------------


subroutine abs_smooth(x, delta_x, y)
    ! absolute value function with small quadratic in the valley
    ! so that it is C1 continuous

    implicit none
    integer, parameter :: dp = kind(0.d0)

    ! in
    real(4), intent(in) :: x, delta_x

    ! out
    real(4), intent(out) :: y


    if (x >= delta_x) then
        y = x
    elseif (x <= -delta_x) then
        y = -x
    else
        y = x**2/(2.0_dp*delta_x) + delta_x/2.0_dp
    end if

end subroutine abs_smooth



subroutine setupa(n, xpt, ypt, p1, p2, p3) !, delta_x)
    ! setup the Akima spline function

    implicit none

    !integer, parameter :: dp = kind(0.d0)
    real(4), parameter :: eps = 1d-30

    ! in
    integer, intent(in) :: n
    real(4), dimension(n), intent(in) :: xpt, ypt  ! given points
    !real(4), intent(in) :: delta_x
    !f2py real(dp), intent(in) :: delta_x = 0.1

    ! out
    real(4), dimension(n-1), intent(out) :: p1(:), p2(:), p3(:)  ! spline coefficients

    ! local
    integer :: i
    real(4), dimension(-1:n+1) :: m
    real(4), dimension(n) :: t
    real(4) :: w1, w2 !,m1, m2, m3, m4
    real(4) ::  dx, t1,t2 !dx(n-1) !,t1(n-1), t2(n-1)

    ! compute segment slopes
    do i = 1, n-1
        m(i) = (ypt(i+1) - ypt(i)) / (xpt(i+1) - xpt(i))
    end do
    
    ! compute segment slopes (array syntax for speed?) - not sure if working
!     dx = xpt(2:) - xpt(:n-1)
!     m(1:n-1) = (ypt(2:) - ypt(:n-1)) / (dx)
    
    

    ! estimation for end points
    m(0) = 2.0*m(1) - m(2)
    m(-1) = 2.0*m(0) - m(1)
    m(n) = 2.0*m(n-1) - m(n-2)
    m(n+1) = 2.0*m(n) - m(n-1)

    ! slope at points
    do i = 1, n
        !m1 = m(i-2)
        !m2 = m(i-1)
        !m3 = m(i)
        !m4 = m(i+1)
        w1 = abs(m(i+1) - m(i))
        w2 = abs(m(i-1) - m(i-2))
        !call abs_smooth(m4 - m3, delta_x, w1)	!swap commenting with above 2 lines if derivatives wanted
        !call abs_smooth(m2 - m1, delta_x, w2)
        if ( w1 < eps .and. w2 < eps ) then
            t(i) = 0.5*(m(i-1) + m(i))  ! special case to avoid divide by zero
        else
            t(i) = (w1*m(i-1) + w2*m(i)) / (w1 + w2)
        end if
    end do

    ! polynomial coefficients
    do i = 1, n-1
        dx = xpt(i+1) - xpt(i)
        t1 = t(i)
        t2 = t(i+1)
        p1(i) = t1
        p2(i) = (3.0*m(i) - 2.0*t1 - t2)/dx
        p3(i) = (t1 + t2 - 2.0*m(i))/dx**2
    end do


	! polynomial coefficients (array syntax. dx calculated above) - not sure is working
	!p2 = (3.0*m(1:n-1) - 2.0*t(:n-1) - t(2:))/dx
	!p3 = (t(:n-1) + t(2:) - 2.0 * m(1:n-1))/dx**2


end subroutine setupa


!same as top of file, but doesn't calculate h or right bracket
subroutine r8vec_bracket_noh ( x, xval, left, dt,error)

  implicit none


  integer ( kind = 4 ) i,j,n
  integer ( kind = 4 ), intent(out) ::  left(:)
  integer right
  real ( kind = 4 ) x(:)
  real ( kind = 4 ) xval(:)
  integer mid
  

  real ( kind = 4 ), intent(out) ::  dt(:)
  logical, intent(out) ::  error(:)
 
  
  n = size(x)
  
  !Interval spans full bounds by default
  left = 1
  error = .false. 
  
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!



	do j = 1,size(xval)
	
		!Check if it's in the outer bounds
		  if (xval(j) < x(1) .or. xval(j) > x(n)) then
			 error(j) = .true.
			 cycle
		  end if
	  
		!Otherwise, use binary search to find interval

		right = n
		do while (right > left(j) + 1)

			  mid = ( left(j) + right + 1 ) / 2

			  if ( x(mid) <= xval(j) ) then
				left(j) = mid
			  else
				right = mid
			  end if
	
		end do
   end do
  
  dt = xval - x(left)

  return
end subroutine




!same as top of file, just finds the point to the right
!The idea being that no interpolation is performed later,
!Rather, the pile depth is rounded up to the nearest element
!The settlement is therefore known exactly

!NOTE: ASSUMES X VALUES ARE DECREASING
subroutine r8vec_bracket_findright ( x, xval, right)

  implicit none


  integer ( kind = 4 ) i,j,n
  integer ( kind = 4 ), intent(out) ::  right(:)
  integer left
  real ( kind = 4 ) x(:)
  real ( kind = 4 ) xval(:)
  integer mid

 
  
  n = size(x)
  

  
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!



	do j = 1,size(xval)
	
		!Check if it's in the outer bounds
		  if (xval(j) > x(1)) then
               right(j) = -100
               cycle
          else if(xval(j) < x(n)) then
			 !error(j) = .true.
              right(j) = -200
			 cycle
		  end if
	  
		!Otherwise, use binary search to find interval

		right(j) = n
        left = 1
		do while (right(j) > left + 1)

			  mid = ( left + right(j) + 1 ) / 2

			  if ( x(mid) > xval(j) ) then
				left = mid
			  else
				right(j) = mid
			  end if
	
		end do
   end do
  

  return
end subroutine

end module
