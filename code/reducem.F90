module reducem
    
    contains
    
    subroutine reducesi(nze,nlayer,n_bh,n_bh_max,inv_coords,inv_reduction,evals,percentile,s_dev,npiles,plocation,mindist,power,goodvals,nsamples)

! This subroutine reduces site investigation data to a given representative value per layer (optionally specific to a each footing)

! Note that all units are in terms of soil elements. If you are working in metres, please convert before hand.

! Subroutine returns the reduced value of each investigation, according to either:
!             1 - standard arithmetic average
!             2 - geometric average
!             3 - harmonic average
!             4 - minimum
!             5 - first quartile
!			  6 - standard deviation below mean (geometric/lognormal)

! Note that the 1Q reduction technique is relatively slow as it requires first sorting the samples in 
! ascending order. Partial sorting up to the desired percentile is used to increase speed.

! Inverse distance weighted reductions of arbitrary power are provided for the arithmetic, geometric and harmonic averages,
! as well as 'borehole' minimums, which take the minimum of the average of each borehole. Theoretically, these distance weighting
! and minimum techniques could be applied to other reduction methods such as the first quartile and standard deviation below the mean,
! however that's probably overkill in terms of complexity as well as overly-conservative. The reason why these methods have only been
! applied to the averages is that the average should *theoretically* be a reasonable representation of the effective modulus at a given
! location, and doing anything on top of that builds in a degree of redundancy. The other methods already have some redundancy built-in,
! an so do not require any more.

! See "QUANTIFYING THE RISK OF GEOTECHNICAL SITE INVESTIGATIONS" - 2006 PhD Thesis by J. Goldsworthy for details.
!
! Note: The unlisted reduction methods need to be tweaked so that they return the -101 value if there are zero samples in the current layer

! RETURNS:
! An array of reduced values for each soil layer from each investigation (evals)

! Michael Crisp February 2018



  implicit none

! -- input variables --
  integer,		intent(in) :: nze								!maximum no. samples per borehole (field depth)
  integer,		intent(in) :: nlayer							!number of soil layers
  integer, 		intent(in) :: n_bh_max,n_bh                     !max number of boreholes, number of boreholes in current investigation
  integer,		intent(in) :: inv_coords (:,:) !(n_bh_max,2)            !x,y coordinates of boreholes for each investigation
  integer, 		intent(in) :: inv_reduction						!reduction method for each test	
  real,			intent(in) :: percentile						!percentile to use in reduction method (proportion)
  real,			intent(in) :: s_dev								!standard deviations below mean for reduction method (geometric)
  integer,		intent(in) :: npiles							!number of piles (relevant for inverse-distance methods)
  integer,		intent(in) :: plocation(:,:) !(npiles,2)				!x,y coordinates for each pile (relevant for inverse-distance methods)
  integer,		intent(in) :: mindist							!minimum distance between footing and investigation for inverse distance methods
  real,			intent(in) :: power								!power of inverse-distance weighting
  real,         intent(in) :: goodvals(:,:,:) ! (nlayer,n_bh_max,nze)         !store valid values of SI samples in 1D array
  integer,		intent(in) :: nsamples(:,:) !(nlayer,n_bh_max)				!Number of viable samples in each layer in each borehole

! -- output variables --
  real,			intent(out) :: evals(:,:) ! (nlayer,npiles)				!reduced value for each investigation
  
! -- local variables --

  integer		inv,bh,j,k,i,d,layer,pile                              !loop counters
  real          nan
  
  !various arrays for storing samples, and masks of valid values
  

  real          lgoodvals(n_bh*nze)                  		!vector version of the borehole data, sometimes the log of the values
  
  ! variables related to removing invalid values
  real          gmn,gsd                                     !geometric mean and standard deviation of SI samples
  integer 		cumsum										!cumulative sum
  
  !reduction variables
  real          Q                                           !approximate index of 1st quartile
  real          temp                                        !temp storage in sorting algorithm
  real 			dist(n_bh)									!distance between foundation and borehole
  real			tempred(n_bh)								!reduced values per borehole
  real			weightsum(n_bh)								!sum of weights for a given borehole
  
!functions
  integer iseed

	
	!set reduced values to negative
	!anything still negative at the end means that there were no viable samples in a particular layer
    nan = 0.0 !generate undefined value
  	nan = 0.0/nan
	evals = nan
	
	!write(*,*) goodvals(layer,bh,:0)
	
        
!----- Global reduction methods ------
        
        !arithmetic average
        if (inv_reduction == 1) then
        	do layer = 1,nlayer
        		cumsum = 0
				do bh = 1,n_bh
					lgoodvals(cumsum + 1:cumsum + nsamples(layer,bh)) = goodvals(layer,bh,:nsamples(layer,bh)) 
					cumsum = cumsum + nsamples(layer,bh)
				end do
				if(cumsum > 0) then	!if there are no viable samples in the entire layer, return error value
					evals(layer,:) = sum(lgoodvals(:cumsum))/cumsum
				else
					evals(layer,:) = nan
				end if
            end do
        
        !geometric average
        else if (inv_reduction == 2) then
        	do layer = 1,nlayer
        		cumsum = 0
				do bh = 1,n_bh
					lgoodvals(cumsum + 1:cumsum + nsamples(layer,bh)) = goodvals(layer,bh,:nsamples(layer,bh))
					cumsum = cumsum + nsamples(layer,bh)
				end do
				if (cumsum > 0) then
            		evals(layer,:) = exp(sum(log(lgoodvals(:cumsum))/cumsum))
            	else
            		evals(layer,:) = nan
            	end if
            end do
        
        !harmonic average
        else if (inv_reduction == 3) then
        	do layer = 1,nlayer
        		cumsum = 0
				do bh = 1,n_bh
					lgoodvals(cumsum + 1:cumsum + nsamples(layer,bh)) = goodvals(layer,bh,:nsamples(layer,bh))
					cumsum = cumsum + nsamples(layer,bh)
				end do
				if (cumsum > 0) then
	            	evals(layer,:) = cumsum/(sum(1/lgoodvals(:cumsum)))
	            else
	            	evals(layer,:) = nan
	            end if
            end do
        
        !minimum
        else if (inv_reduction == 4) then
        	do layer = 1,nlayer
        		cumsum = 0
				do bh = 1,n_bh
					lgoodvals(cumsum + 1:cumsum + nsamples(layer,bh)) = goodvals(layer,bh,:nsamples(layer,bh))
					cumsum = cumsum + nsamples(layer,bh)
				end do
				if (cumsum > 0) then
					evals(layer,:) = minval(lgoodvals(:cumsum))
				else
					evals(layer,:) = nan
				end if
        	end do
            
        !percentile
        else if (inv_reduction == 5) then
        
        	do layer = 1,nlayer
        		cumsum = 0
				do bh = 1,n_bh
					lgoodvals(cumsum + 1:cumsum + nsamples(layer,bh)) = goodvals(layer,bh,:nsamples(layer,bh))
					cumsum = cumsum + nsamples(layer,bh)
                end do
                
                !goodvals(1,1,:14)
                
                if (cumsum == 0) then
                	evals(layer,:) = nan
				else
			
					Q = percentile*(cumsum-1)+1 !index of 1st quartile

					!need to sort samples
					!Loop end is normally 2, but this is faster (only relevant bit is sorted)
					!E.g for 25th percentile, only need to sort 25% of the data.
					do k = cumsum,max(2,cumsum-ceiling(Q)),-1	 
						do j = cumsum-1,cumsum-k+1,-1
							if (lgoodvals(j+1) < lgoodvals(j)) then 
								temp = lgoodvals(j)
								lgoodvals(j) = lgoodvals(j+1)
								lgoodvals(j+1) = temp
							end if
						end do
					end do   
				  
					!If percentile falls between two values, linearly interpolate the two bounds by the corresponding proportion
					evals(layer,:) = (lgoodvals(ceiling(Q)) - lgoodvals(floor(Q)))*(percentile*(cumsum-1)-floor(percentile*(cumsum-1))) + lgoodvals(floor(Q))
				end if
            end do
            
        !geometric standard deviations below geometric mean
        else if (inv_reduction == 6) then
        
        	do layer = 1,nlayer
        		cumsum = 0
				do bh = 1,n_bh
					lgoodvals(cumsum + 1:cumsum + nsamples(layer,bh)) = goodvals(layer,bh,:nsamples(layer,bh))
					cumsum = cumsum + nsamples(layer,bh)
				end do
				
				if (cumsum == 0) then
				
					evals(layer,:) = nan
        		else
					!take log of values as first step to calculate geometric mean
					lgoodvals(:cumsum) = log(lgoodvals(:cumsum))
		
					!take mean of values for 2nd step
					gmn = sum(lgoodvals(:cumsum))/cumsum
		
					!calculate standard deviation
					gsd = 0
					do i = 1,cumsum
						gsd = gsd + (lgoodvals(i)-gmn)**2
					end do
		
					!Actual geometric mean and standard deviation
					gsd = exp(sqrt(gsd/cumsum))
					gmn = exp(gmn)
			
					!value gsd standard deviations below mean
					evals(layer,:) = gmn/(gsd**s_dev)
				end if
			end do
			
			
!-------- inverse distance methods ------
			
		!inverse distance arithmetic
		else if (inv_reduction == 7) then
		
			do pile = 1,npiles
				do layer = 1,nlayer
					cumsum = 0
					do bh = 1,n_bh
						!only look at boreholes where there are valid samples
						if(nsamples(layer,bh) > 0) then
							!get inverse distances between foundation and borehole, ensuring that distance is greater than a specified minimum (typically zero or footing radius)
							dist(bh) = sqrt(real((inv_coords(bh,1)-plocation(pile,1))**2 + (inv_coords(bh,2)-plocation(pile,2))**2))	!get inverse distance
							if(dist(bh) < mindist) then
								dist(bh) = (1/mindist)**power		!convert to weight
							else
								dist(bh) = (1/dist(bh))**power		!convert to weight
							end if
				
							!weight borehole values by their distance
							lgoodvals(cumsum + 1:cumsum + nsamples(layer,bh)) = goodvals(layer,bh,:nsamples(layer,bh)) * dist(bh)
							weightsum(bh) = nsamples(layer,bh) * dist(bh)	!calculate sum of weights for borehole
							cumsum = cumsum + nsamples(layer,bh)			
						else
							!no data at current borehole, set to zero weight. Do not add to cumsum.
							weightsum(bh) = 0.0
						end if
					end do
					!perform weighted average by summing weighted values and dividing by sum of weights
					evals(layer,pile) = sum(lgoodvals(:cumsum))/sum(weightsum)
				end do
			end do
			
			
		!inverse distance geometric
		else if (inv_reduction == 8) then
		
			do pile = 1,npiles
				do layer = 1,nlayer
					cumsum = 0
					do bh = 1,n_bh
						!only look at boreholes where there are valid samples
						if(nsamples(layer,bh) > 0) then
							!get inverse distances between foundation and borehole, ensuring that distance is greater than a specified minimum (typically zero or footing radius)
							dist(bh) = sqrt(real((inv_coords(bh,1)-plocation(pile,1))**2 + (inv_coords(bh,2)-plocation(pile,2))**2))	!get inverse distance
							if(dist(bh) < mindist) then
								dist(bh) = (1/mindist)**power		!convert to weight
							else
								dist(bh) = (1/dist(bh))**power		!convert to weight
							end if
				
							!weight borehole values by their distance
							lgoodvals(cumsum + 1:cumsum + nsamples(layer,bh)) = log(goodvals(layer,bh,:nsamples(layer,bh))) * dist(bh)
							weightsum(bh) = nsamples(layer,bh) * dist(bh)	!calculate sum of weights for borehole
							cumsum = cumsum + nsamples(layer,bh)
						else
							!no data at current borehole, set to zero weight. Do not add to cumsum.
							weightsum(bh) = 0.0
						end if
					end do
					!perform weighted average by summing weighted values and dividing by sum of weights
					evals(layer,pile) = exp(sum(lgoodvals(:cumsum)/sum(weightsum)))
				end do
			end do
				
			
		!inverse distance harmonic
		else if (inv_reduction == 9) then
		
			do pile = 1,npiles
				do layer = 1,nlayer
					cumsum = 0
					do bh = 1,n_bh
						!only look at boreholes where there are valid samples
						if(nsamples(layer,bh) > 0) then
							!get inverse distances between foundation and borehole, ensuring that distance is greater than a specified minimum (typically zero or footing radius)
							dist(bh) = sqrt(real((inv_coords(bh,1)-plocation(pile,1))**2 + (inv_coords(bh,2)-plocation(pile,2))**2))	!get inverse distance
							if(dist(bh) < mindist) then
								dist(bh) = (1/mindist)**power		!convert to weight
							else
								dist(bh) = (1/dist(bh))**power		!convert to weight
							end if
				
							!weight borehole values by their distance
							lgoodvals(cumsum + 1:cumsum + nsamples(layer,bh)) = dist(bh)/goodvals(layer,bh,:nsamples(layer,bh))
							weightsum(bh) = nsamples(layer,bh) * dist(bh)	!calculate sum of weights for borehole
							cumsum = cumsum + nsamples(layer,bh)
						else
							!no data at current borehole, set to zero weight. Do not add to cumsum.
							weightsum(bh) = 0.0
						end if
					end do
					!perform weighted average by summing weighted values and dividing by sum of weights
					evals(layer,pile) = sum(weightsum)/sum(lgoodvals(:cumsum))
				end do
			end do
			
			
!-------- borehole minimum distance methods ------

        !arithmetic minimum
        else if (inv_reduction == 10) then
        	do layer = 1,nlayer
        		cumsum = 0
				do bh = 1,n_bh
					if(nsamples(layer,bh) > 0) then														!only take average if samples exist
						cumsum = cumsum + 1
						tempred(cumsum) = sum(goodvals(layer,bh,:nsamples(layer,bh)))/nsamples(layer,bh)
					end if
				end do
				evals(layer,:) = minval(tempred(:cumsum))										!take worst case of each borehole
            end do
        
        !geometric minimum
        else if (inv_reduction == 11) then
        	do layer = 1,nlayer
        		cumsum = 0
				do bh = 1,n_bh
					if(nsamples(layer,bh) > 0) then														!only take average if samples exist
						cumsum = cumsum + 1
						tempred(cumsum) = exp(sum(log(goodvals(layer,bh,:nsamples(layer,bh)))/nsamples(layer,bh)))
					end if
				end do
				evals(layer,:) = minval(tempred(:cumsum))										!take worst case of each borehole
            end do
        
        !harmonic minimum
        else if (inv_reduction == 12) then
        	do layer = 1,nlayer
        		cumsum = 0
				do bh = 1,n_bh
					if(nsamples(layer,bh) > 0) then														!only take average if samples exist
						cumsum = cumsum + 1
						tempred(cumsum) = nsamples(layer,bh)/(sum(1/goodvals(layer,bh,:nsamples(layer,bh))))
					end if
				end do
				evals(layer,:) = minval(tempred(:cumsum))										!take worst case of each borehole
            end do
            
        else
            write(*,*) "Error: Invalid reduction method chosen for investigation",inv_reduction
            write(*,*) "Valid values are:"
            write(*,*) "1 - standard arithmetic average"
            write(*,*) "2 - geometric average"
            write(*,*) "3 - harmonic average"
            write(*,*) "4 - minimum"
            write(*,*) "5 - percentile"
            write(*,*) "6 - geometric SD"
            stop
        end if

    end subroutine
    
end module
