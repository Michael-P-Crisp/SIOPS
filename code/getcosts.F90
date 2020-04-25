

subroutine getcosts(diffset,costs,x1,x2,y1,y2,niter,ninv)
         
!Calculate failure costs based on max differential settlement
!Assumes cost is given by a bounded linear function, where:
! The cost is constant (the minimum) below the bound
! The cost is constant (the maximum) above the bound
! The cost is linearly interpolated between the bounds

! Michael Crisp February 2018

use extrafuncs

  implicit none

	integer, intent(in) :: niter,ninv !no. of realisations, investigations
	
	real, intent(in) :: x1,x2 !lower and upper bound differential settlements corresponding to costs
	real, intent(in) :: y1,y2 !lower and upper bound costs

	real slope !slope of linear function
	
	real :: diffset(niter,ninv) !differential settlement input
	
	integer :: iter,inv !loop counters
	
	!output variables
	real, intent(out) :: costs(niter,ninv)
	
	!calculate slope here since it's constant
	slope = (y2-y1)/(x2-x1)
  	
  	do inv = 1,ninv
		do iter = 1,niter
            if(diffset(iter,inv) < 0) then
                costs(iter,inv) = -1
			else if(diffset(iter,inv) <= x1) then
				costs(iter,inv) = y1
			else if(diffset(iter,inv) >= x2) then
				costs(iter,inv) = y2
			else
				costs(iter,inv) = (diffset(iter,inv) - x1) * slope + y1
			end if 
		end do
	end do


end subroutine
