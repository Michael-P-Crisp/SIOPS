subroutine getbhdepths(nxe,nye,nze,bfld,nlayer,n_inv,n_bh,inv_bh,inv_coords,inv_depths,bound_depths)

!*****************************************************************************80
! Get SI performance (reduced young's modulus of each layer + effective layer depth for each pile)
!
  implicit none
  
  !all of this is site investigation stuff
  integer,		intent(in) :: nxe,nye,nze						!number of soil elements in x,y,z dimensions
  real, 		intent(in) :: bfld(nlayer+1,nxe,nye) 			!n x 2D array of n+1 layer boundaries. Bottom boundary for layer i is "bfld(i+i,:,:) + 1".
  integer,		intent(in) :: nlayer							!number of soil layers
  integer, 		intent(in) :: n_inv 							!total number of investigations
  integer, 		intent(in) :: n_bh								!max number of boreholes
  integer, 		intent(in) :: inv_bh(n_inv)						!number of boreholes in each investigation
  integer,		intent(in) :: inv_coords(n_inv,n_bh,2)			!x,y coordinates of boreholes for each investigation
  integer, 		intent(in) :: inv_depths(n_inv,3) 				!sampling depth start, end, interval


! -- output variables --
  real,		intent(out) :: bound_depths(n_inv,nlayer+1,n_bh)			!layer boundary information as encountered by each borehole and investigation    bound_depths(n_inv,nlayer+1,n_inv,n_bh)
  
! -- local variables --

  integer		inv,bh,j,k,i,d,layer                        !loop counters
  integer		depths(nze),scount							!vector containing depths of samples, no. samples
  logical 		inside(nze) 								!samples which are inside a given layer
  integer       tscount                                     !total sample count in current investigation
  integer		stepback									!number of elements to step backwards for discrete sampling (depth interval/2)
  integer       ctrue,ctrue2                                !count of true values
  


! #########################		GET LAYER DEPTHS AT BOREHOLES ACCORDING TO INVESTIGATION		#########################	
  
  
  !set upper and lower boundary depths to top and bottom of field
  bound_depths(:,1,:) = 0
  bound_depths(:,nlayer+1,:) = nze 

  do inv = 1,n_inv
  
	    !depths = -101
	    !scount = 0
	    !do d = inv_depths(inv,1),inv_depths(inv,2),inv_depths(inv,3)
		   ! scount = scount + 1
		   ! depths(scount) = d
	    !end do

	
	    !get stepback (0 for continuous)

	
	    !if more than one layer, collect intermediate boundaries
        !height of current layer is based on the top of the layer below?
	    !do layer = nlayer,2,-1       
		   ! do bh = 1,inv_bh(inv)	!loop through boreholes
			  !  inside(:scount) = depths(:scount) > bfld(layer,inv_coords(inv,bh,1),inv_coords(inv,bh,2)) &
					!     .and. depths(:scount) <=  bfld(layer+1,inv_coords(inv,bh,1),inv_coords(inv,bh,2)) !get all depths in current layer
			  !  ctrue = count(inside(:scount)) !count samples in current layer
			  !
			  !  !if not penetrated, set the boundary to the upper boundary of the lower layer (or some kind of nan?)
			  !  if (ctrue == 0) then
				 !   bound_depths(inv,layer,bh) = bound_depths(inv,layer+1,bh)
				 !   cycle
			  !  else !take first sample in new boundary as layer boundary if continuous, or halfway between present and previous depth of first encountered sample within layer
				 !   bound_depths(inv,layer,bh) = max(minval(pack(depths(:scount),inside(:scount))) - stepback - 1,0) 		!The -1 accounts for fortran/python conversion
     !           end if
     !           
     !
		   ! end do
     !   end do
        
            do bh = 1,inv_bh(inv)	!loop through boreholes
			    bound_depths(inv,:,bh) = bfld(:,inv_coords(inv,bh,1),inv_coords(inv,bh,2)) !get direct depths from boundary
            end do
            
            !If the current test is discrete round the depths up to the nearest sample, and subtract half the increment
            if(inv_depths(inv,3) > 1) then
                bound_depths(inv,:,:) = inv_depths(inv,3)*ceiling(bound_depths(inv,:,:)/inv_depths(inv,3))-inv_depths(inv,3)/2
            end if
    end do
	
	

end subroutine


