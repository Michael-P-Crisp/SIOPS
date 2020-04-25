subroutine getheights(bfld,indices,plocation,xdim,ydim,npl,pdim,power,radius,nlayer,height)

!get soil layer depths at individual piles. Specifically, the 'effective' layer depths,
!calculated as an inverse-distance weighted average, with distances to an arbitrary power.
!This reflects how soil further away from a pile has a smaller impact on its performance.
!A contant height for the layers at each pile is required for both the 2D axisymmetric analysis
!as well as the Mylonakis and Gazetas method for pile settlement.
    
implicit none

integer, intent(in) :: xdim,ydim,npl,power
real, intent(in) :: bfld(nlayer-1,xdim,ydim)
integer, intent(in) :: indices(2,xdim,ydim)
integer, intent(in) :: plocation(npl,2),radius(2)

integer, intent(in) :: pdim(2)
integer, intent(in) :: nlayer 

integer :: xlow,xhigh,ylow,yhigh
real :: dist(xdim,ydim)
integer :: p,i,j
real sumweights

real, intent(out) :: height(nlayer-1,npl)

real start,finish




!Loop through the piles
do p = 1,npl
    !Get the upper and lower bounds of soil indices in each dimension. I.e. a square of height information around the pile.
    
    call CPU_TIME(start)
    !do j = 1,1000
    
	    xlow = plocation(p,1) - radius(1)
	    ylow = plocation(p,2) - radius(2)
	    xhigh= plocation(p,1) + radius(1) + pdim(1) - 1 
	    yhigh= plocation(p,2) + radius(2) + pdim(2) - 1 

	    dist = sqrt((indices(1,:,:)-plocation(p,1)-float(pdim(1))/4)**2 + (indices(2,:,:)-plocation(p,2)-float(pdim(2))/4)**2)  !Calculate the distances of each soil element from the pile using pythagoras
	    dist = dist**(-power) !Multiply them by an arbitrary power (negative for inverse).
	    dist(plocation(p,1):plocation(p,1)+pdim(1)-1,plocation(p,2):plocation(p,2)+pdim(2)-1) = 0 !The elements within the pile radius must be zero, as they are replaced by the pile itself.
        sumweights = sum(dist(xlow:xhigh,ylow:yhigh)) !get the sum of the weights as part of the weighted average calculation
    
    !end do
    call cpu_time(finish)
    !write(*,*) 'prep: ',(finish - start)/1000
    
    call CPU_TIME(start)
    !do j = 1,1000
    
        do i = 1,nlayer-1 !get the effective heights for each layer.        
	        height(i,p) = sum(bfld(i,xlow:xhigh,ylow:yhigh) * dist(xlow:xhigh,ylow:yhigh)) / sumweights
        end do
    
    !end do
    call cpu_time(finish)
    !write(*,*) 'height: ',(finish - start)/1000
    
    !read(*,*)

end do

end subroutine