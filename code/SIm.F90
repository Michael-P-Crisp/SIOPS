

    !multiple layer implementation of site investigations. Needs both the layer's Young's modulus and an array of boundary information.
    
subroutine do_SI(goodvals,nsamples,nxe,nye,nze,emn,bfld,nlayer,n_inv,n_bh,inv_bh,inv_coords,inv_depths,num_tests,num_errors,inv_test,add_errors,use_CI,test_errors,kseed,CI,rand_realisation,lmean_ln,lsd_ln,efld2D,C01D)
                                               !^revert to efld

use si
use LAS1D

! This subroutine does site investigations of a 3D soil field.
! Options include adding test errors (bias, random, transformation) (zero and negative values are discarded)
! It can also remove values that don't conform to a geometric confidence interval of a user-specified level (recommend 99%)

! Note that all units are in terms of soil elements. If you are working in metres, please convert before hand.

! Subroutine returns the reduced value of each investigation, according to either:
!             1 - standard arithmetic average
!             2 - geometric average
!             3 - harmonic average
!             4 - minimum
!             5 - first quartile
!			  6 - standard deviation below mean (geometric/lognormal)

! Each investigation is entirely independant, and can be specified by its:
!	Number of boreholes
!	Borehole coordinates
!	Sampling depths (given as start,stop and interval)
!	Test type
!	Reduction method (including S.D. below mean and percentile)
!	Confidence interval to truncate (if specified)
!   

! Random errors are added based on unit-mean lognormal random values. These values are generated based on
! a random seed value. Using the same seed each time should produce the same values, all else being equal.
! Set the seed to zero to use the system clock as an initial seed (although not recommended).
! Errors represent 'driller' error applied on a per-borehole basis, measurement error applied on per-sample basis, 
! and parameter transformation error, applied on a global basis.

! Note that the 1Q reduction technique is relatively slow as it requires first sorting the samples in 
! ascending order. Partial sorting up to the desired percentile is used to increase speed.

! Note due to the complexity and shear amount of inputs, no input safety checks have been implemented.
! Essentially, just make sure that all values are >= 0 (besides kseed and percentile) and it should be fine.

! See "QUANTIFYING THE RISK OF GEOTECHNICAL SITE INVESTIGATIONS" - 2006 PhD Thesis by J. Goldsworthy for details.

! RETURNS:
! evals - An array of reduced values for each soil layer from each investigation. The viable samples are stacked
! at the front of the nze index.
! nsamples - The number of viable samples from each borehole in each investigation. Used to extract evals.
! bound_depths - Height information for each layer as encountered by testing

! Michael Crisp February 2018

!NOTE: THIS CODE ORIGINALLY WORKED WITH MULTI-LAYERED PROFILES THAT HAD VARIABLE SOIL WITHIN EACH LAYER
!       IT HAS BEEN MODIFIED SO THAT THE AVERAGE LAYER PROPERTY IS GIVEN DIRECTLY AS OPPOSED TO AN EXTRACTED
!       COLUMN OF VALUES. SEARCH FOR "revert" IF YOU WANT TO UNDO THIS.


  implicit none

! -- input variables --
  integer,		intent(in) :: nxe,nye,nze						!number of soil elements in x,y,z dimensions
  !real, 		intent(in) :: efld(nxe,nye,nze) 				!3D array of soil E values <- revert
  real,         intent(in) :: emn(nlayer)                       !the mean E value in each layer. Assumes soil is uniform.
  real, 		intent(in) :: bfld(nlayer+1,nxe,nye) 			!n x 2D array of n+1 layer boundaries. Bottom boundary for layer i is "bfld(i+i,:,:) + 1".
  integer,		intent(in) :: nlayer							!number of soil layers
  integer, 		intent(in) :: n_inv 							!total number of investigations
  integer, 		intent(in) :: n_bh								!max number of boreholes
  integer, 		intent(in) :: inv_bh(n_inv)						!number of boreholes in each investigation
  integer,		intent(in) :: inv_coords(n_inv,n_bh,2)			!x,y coordinates of boreholes for each investigation
  integer, 		intent(in) :: inv_depths(n_bh,n_inv,3) 				!sampling depth start, end, interval
  integer,		intent(in) :: num_tests							!number of tests
  integer, 		intent(in) :: num_errors						!number of error sources for tests
  integer,	    intent(in) :: inv_test(n_bh,n_inv)					!test type for each investigation
  logical,		intent(in) :: add_errors						!whether to add test errors or not
  logical,		intent(in) :: use_CI							!whether to truncate based on confidence interval
  real,			intent(in) :: CI    							!confidence interval to be applied to each investigation
  real, 		intent(in) :: test_errors(num_tests,num_errors)	!matrix of test error stats (num tests x num errors)
  integer,      intent(in) :: kseed                             !initial random seed
  
  real, intent(in) :: efld2D(nlayer,nxe,nye)     !2D random field of soil properties representing Young's modulus
  real, intent(in) :: C01D(nze) ! correlation data for generating 1D random fields
  real, intent(in) :: lmean_ln(nlayer),lsd_ln(nlayer) !lognormal statistics for each layer for generating random values      <- revert; comment this out
  integer, intent(in) :: rand_realisation ! 3 = generate 1D random field for each borehole.  

! -- output variables --
  real,			intent(out) :: goodvals(n_inv,nlayer,n_bh,nze)                   !store valid values of SI samples in 1D array
  integer,		intent(out) :: nsamples(n_inv,nlayer,n_bh)						!keep track of how many viable samples are in each layer in each borehole
  !integer,		intent(out) :: bound_depths(n_inv,nlayer+1,n_bh)			!layer boundary information as encountered by each borehole and investigation    bound_depths(n_inv,nlayer+1,n_inv,n_bh)
  
! -- local variables --

  integer		inv,bh,j,k,i,d,layer                        !loop counters
  integer		depths(nze,n_bh),scount(n_bh)							!vector containing depths of samples, no. samples
  logical 		inside(nze) 								!samples which are inside a given layer
  integer       tscount                                     !total sample count in current investigation
  integer		stepback									!number of elements to step backwards for discrete sampling (depth interval/2)
  
  !various arrays for storing samples, and masks of valid values
 
  
  !real			sidata(n_bh,nze)							!store SI samples (up to possible maximum)
  real			sidata_rand(nlayer,n_bh,nze)                !store fake random SI samples (up to possible maximum)
  
  real          lgoodvals(n_bh*nze)                  		!vector version of the above info, sometimes the log of the values
  logical       gvmask(nze)				                    !mask invalid values of the above info
  

  ! variables related to test errors 
  integer       seed                                        !random seed
  real          rand_vals(n_bh*nze), rand_vals_small(nze)  !random values
  real          bhmeans(n_bh)                               !mean of borehole values
  real          Ltest_errors(num_tests,num_errors)          !lognormal test errors
  real          layersum(num_tests,nlayer)                                  !track the sum of values in each layer
  integer       layersumcount(num_tests,nlayer)                             !track the number of samples in each layer
  
  ! variables related to removing invalid values
  integer       ctrue,ctrue2                                !count of true values
  real          gmn,gsd                                     !geometric mean and standard deviation of SI samples
  real          lower,upper                                 !lower and upper bounds of geometric confidence interval
  integer 		cumsum										!cumulative sum
  real          nan
  
  
!functions
  integer iseed
  
  
	
	!initialise random number generator
	seed = iseed(max(0,kseed)) !ensure kseed is at least zero
	
	!set reduced values to negative
	!anything still negative at the end means that there were no viable samples in a particular layer
    nan = 0.0 !generate undefined value
  	nan = 0.0/nan
	goodvals = nan
    
    
    !convert normal COVs to lognormal parameters
    Ltest_errors = sqrt(log(1.0 + test_errors ** 2))
	
	
	do inv = 1,n_inv
	
! #########################		SITE INVESTIGATION FOR SOIL PROPERTIES		#########################	
        
        
        !----- do investigation and associated calcs -------

        ! build vector of sample depths
		tscount = 0 !total number of samples in current investigation
        scount = 0 !number of samples in current borehole
        do bh = 1,inv_bh(inv)
		    do d = inv_depths(bh,inv,1),inv_depths(bh,inv,2),inv_depths(bh,inv,3) 
			    scount(bh) = scount(bh) + 1
			    depths(scount(bh),bh) = d
            end do
            tscount = tscount + scount(bh)
        end do
		
		!write(*,*) 'h1'

	
		!Do investigation; extract all relevant samples and repack into a 1D vector of values
		!It'll extract the full column for now. Discrete samples at the relevant depths will be extracted later
        !NOTE: THIS HAS BEEN COMMENTED OUT <- revert
		!do bh = 1,inv_bh(inv)	!loop through boreholes
		!	!sidata(bh,:inv_depths(inv,2)) = efld(inv_coords(inv,bh,1),inv_coords(inv,bh,2),:inv_depths(inv,2))
		!	sidata(bh,:scount) = efld(inv_coords(inv,bh,1),inv_coords(inv,bh,2),depths(:scount))
		!end do

		
		!vectorised equivalent, although need to tweak the ordering of the reshape for it to be correct.
		!Apparently this is slower than the above anyway
        !sidata(:inv_bh(inv),:scount) = reshape(efld(inv_coords(inv,:inv_bh(inv),1),inv_coords(inv,:inv_bh(inv),2),depths(:scount)),[inv_bh(inv),scount])
        
    	!write(*,*) 'h2'
        
        
        ! Get sample values
        if (rand_realisation <= 2) then ! Take from layer mean
            do layer = 1,nlayer
                sidata_rand(layer,:,:) = emn(layer)
            end do
        else if (rand_realisation == 3) then ! Take constant value from 2D random field
            do layer = 1,nlayer
                do bh = 1,inv_bh(inv)
                    sidata_rand(layer,bh,:scount(bh)) = exp(efld2D(layer,inv_coords(inv,bh,1),inv_coords(inv,bh,2)))
                end do
            end do
        else		! Take from 2D random field as above, but add a 1D random field
            do layer = 1,nlayer
                do bh = 1,inv_bh(inv)
                    call las1g( rand_vals_small(:inv_depths(bh,inv,2)), inv_depths(bh,inv,2), C01D) !generate zero-mean, unit variance random noise
                    sidata_rand(layer,bh,:scount(bh)) = exp(efld2D(layer,inv_coords(inv,bh,1),inv_coords(inv,bh,2)) + rand_vals_small(depths(:scount(bh),bh)) * lsd_ln(layer))
                end do
            end do
        end if
        
        
        !-----add random errors to samples-----
        
        if (add_errors) then
        	
        	
        
            !borehole average:
            !bhmeans(:inv_bh(inv)) = sum(sidata(:inv_bh(inv),:scount),1)/scount                         !get borehole means                                          
        
            call vnorm(rand_vals,inv_bh(inv))                                                           !generate zero-mean, unit-variance normal values
            layersum = 0
            layersumcount = 0
            do bh = 1,inv_bh(inv)
                rand_vals(bh) = exp(rand_vals(bh) * Ltest_errors(inv_test(bh,inv),1))      !transform to unit-mean lognormal values of desired SD
                do layer = 1,nlayer
            	
            		!get logical array of valid values in current layer, for masking
            		gvmask(:scount(bh)) = depths(:scount(bh),bh) > bfld(layer,inv_coords(inv,bh,1),inv_coords(inv,bh,2)) &
            			 .and. depths(:scount(bh),bh) <=  bfld(layer+1,inv_coords(inv,bh,1),inv_coords(inv,bh,2))
					nsamples(inv,layer,bh) = count(gvmask(:scount(bh)))													!count valid values for current layer and borehole
                    goodvals(inv,layer,bh,:nsamples(inv,layer,bh)) = sidata_rand(layer,bh,:nsamples(inv,layer,bh)) !pack(sidata(bh,:scount),gvmask(:scount))			!collect <- revert
                    
                    !apply bias to boreholes on a per-layer basis
                    goodvals(inv,layer,bh,:nsamples(inv,layer,bh)) = goodvals(inv,layer,bh,:nsamples(inv,layer,bh)) + ( rand_vals(bh) - 1 ) * sum(goodvals(inv,layer,bh,:nsamples(inv,layer,bh)))/nsamples(inv,layer,bh)
				
                    !generate random per-sample errors
                    call vnorm(rand_vals_small,nsamples(inv,layer,bh))
                    rand_vals_small(:nsamples(inv,layer,bh)) = exp(rand_vals_small(:nsamples(inv,layer,bh)) * Ltest_errors(inv_test(bh,inv),2))
                    goodvals(inv,layer,bh,:nsamples(inv,layer,bh)) = goodvals(inv,layer,bh,:nsamples(inv,layer,bh)) * rand_vals_small(:nsamples(inv,layer,bh)) !apply random errors
                    
                    !track the layer statistics
                    layersum(inv_test(bh,inv),layer) = layersum(inv_test(bh,inv),layer) + sum(goodvals(inv,layer,bh,:nsamples(inv,layer,bh)))
                    layersumcount(inv_test(bh,inv),layer) = layersumcount(inv_test(bh,inv),layer) + nsamples(inv,layer,bh)
                end do
            end do
        
 
            do layer = 1,nlayer
                call vnorm(rand_vals,1)
                rand_vals(2:num_tests+1) = exp(rand_vals(1) * Ltest_errors(:,3))
                do bh = 1,inv_bh(inv)
                    ! apply global transformation errors errors
                    goodvals(inv,layer,bh,:nsamples(inv,layer,bh))  = goodvals(inv,layer,bh,:nsamples(inv,layer,bh)) + (rand_vals(inv_test(bh,inv)+1) - 1) * layersum(inv_test(bh,inv),layer)/layersumcount(inv_test(bh,inv),layer)
                    
                    !--Remove values less than or equal to zero--
                                
                    !gvmask(:nsamples(inv,layer,bh)) = goodvals(inv,layer,bh,:nsamples(inv,layer,bh)) > 0
					!ctrue = count(gvmask(:nsamples(inv,layer,bh)))																!count number of valid values
					!goodvals(inv,layer,bh,:ctrue) = pack(goodvals(inv,layer,bh,:nsamples(inv,layer,bh)),gvmask(:nsamples(inv,layer,bh)))	!pack good values back into the array
	
                    !this code is faster than the above version
                    ctrue = 0
                    do i = 1,nsamples(inv,layer,bh)
                        if(goodvals(inv,layer,bh,i) > 0) then
                            ctrue = ctrue + 1
                            goodvals(inv,layer,bh,ctrue) = goodvals(inv,layer,bh,i)
                        end if
                    end do
                    
                    nsamples(inv,layer,bh) = ctrue	
                end do
            end do
            
    
            !write(*,*) 'h4'
            
            
            ! Sort out borehole data within each layer, and remove invalid values (i.e values less than zero, which aren't physically possible)
            !This also takes out samples at discrete depths from the continuous borehole sample column
            !do layer = 1,nlayer
            !	do bh = 1,inv_bh(inv)
            		!get logical array of valid values in current layer, for masking
            !		gvmask(:scount) = depths(:scount) > bfld(layer,inv_coords(inv,bh,1),inv_coords(inv,bh,2)) &
            !			 .and. depths(:scount) <=  bfld(layer+1,inv_coords(inv,bh,1),inv_coords(inv,bh,2)) &
            !			 .and. sidata(bh,:scount) > 0
			!		nsamples(inv,layer,bh) = count(gvmask(:scount))													!count valid values for current layer and borehole
			!		goodvals(inv,layer,bh,:nsamples(inv,layer,bh)) = pack(sidata(bh,:scount),gvmask(:scount))			!collect
                    !write(*,*) goodvals(1,1,1,:14)
                    !write(*,*) goodvals(1,2,1,:14)
			!	end do
			!end do
			
			!write(*,*) 'h5'
        
        else	!otherwise just use the full original data
        	do layer = 1,nlayer
            	do bh = 1,inv_bh(inv)
            		!get logical array of valid values in current layer, for masking
            		gvmask(:scount(bh)) = depths(:scount(bh),bh) > bfld(layer,inv_coords(inv,bh,1),inv_coords(inv,bh,2)) &
            			 .and. depths(:scount(bh),bh) <=  bfld(layer+1,inv_coords(inv,bh,1),inv_coords(inv,bh,2))
					nsamples(inv,layer,bh) = count(gvmask(:scount(bh)))													!count valid values for current layer and borehole
					goodvals(inv,layer,bh,:nsamples(inv,layer,bh)) = sidata_rand(layer,bh,:nsamples(inv,layer,bh)) !pack(sidata(bh,:scount),gvmask(:scount))			!collect <- revert
				end do
			end do
        end if
        
        
   
        
        !---- remove invalid values from data ----
        
        if (use_CI) then ! remove values above and below 99% lognormal CI within given layer
        
			
			do layer = 1,nlayer
			
				!pack all values from a particular layer to a single vector
				!this do block is reused a few times, but IMO its not worth putting it into a subroutine due to its small size
				!and the fact that there is some (minor) overhead in calling a subroutine
				cumsum = 0
				do bh = 1,inv_bh(inv)
					lgoodvals(cumsum + 1:cumsum + nsamples(inv,layer,bh)) = goodvals(inv,layer,bh,:nsamples(inv,layer,bh))
					cumsum = cumsum + nsamples(inv,layer,bh)
                end do
                !nsamples(inv,:,bh)
				
				!write(*,*) 'h6'

				!take log of values as first step to calculate geometric mean
				lgoodvals(:cumsum) = log(lgoodvals(:cumsum))
		
				!take mean of values for 2nd step
				gmn = sum(lgoodvals(:cumsum))/cumsum
		
				!calculate standard deviation
				gsd = 0
				do i = 1,cumsum
					gsd = gsd + (lgoodvals(i)-gmn)**2
				end do
				
				!write(*,*) 'h7'
		
				!Actual geometric mean and standard deviation
				gsd = exp(sqrt(gsd/cumsum))
				gmn = exp(gmn)
		
				!lower and upper confidence interval bounds
				lower = gmn / (gsd * CI)
				upper = gmn * (gsd * CI)
		
				!remove values outside the CI for current layer
				do bh = 1,inv_bh(inv)
					!get mask vector for invalid values
					!gvmask(:nsamples(inv,layer,bh)) = goodvals(inv,layer,bh,:nsamples(inv,layer,bh)) > lower .and. goodvals(inv,layer,bh,:nsamples(inv,layer,bh)) < upper
					!ctrue = count(gvmask(:nsamples(inv,layer,bh)))																!count number of valid values
					!goodvals(inv,layer,bh,:ctrue) = pack(goodvals(inv,layer,bh,:nsamples(inv,layer,bh)),gvmask(:nsamples(inv,layer,bh)))	!pack good values back into the array
					
                    !This code is a slightly faster version of the above
                    ctrue = 0
                    do i = 1,nsamples(inv,layer,bh)
                        if(goodvals(inv,layer,bh,i) > lower .and. goodvals(inv,layer,bh,i) < upper) then
                            ctrue = ctrue + 1
                            goodvals(inv,layer,bh,ctrue) = goodvals(inv,layer,bh,i)
                        end if
                    end do
                
                    nsamples(inv,layer,bh) = ctrue																				!update number of good values
                    !write(*,*) goodvals(1,1,1,:14)
                    !write(*,*) goodvals(1,2,1,:14) goodvals(inv,layer,bh,:)
                end do
            end do			
        end if
		
       
	
		!write(*,*) 'h8'
		
! #########################		GET LAYER BOUNDARIES ACCORDING TO INVESTIGATION		#########################	
        
        !set upper and lower boundary depths to top and bottom of field
!        bound_depths(inv,1,:) = 0
!        bound_depths(inv,nlayer+1,:) = nze !   <---nope!!!!!!
        
        !get stepback (0 for continuous)
!        stepback = inv_depths(inv,3)/2
        
        !if more than one layer, collect intermediate boundaries
!        do layer = nlayer,2,-1       
!        	do bh = 1,inv_bh(inv)	!loop through boreholes
!        		inside(:scount) = depths(:scount) > bfld(layer,inv_coords(inv,bh,1),inv_coords(inv,bh,2)) &
!            			 .and. depths(:scount) <=  bfld(layer+1,inv_coords(inv,bh,1),inv_coords(inv,bh,2)) !get all depths in current layer
!        		ctrue = count(inside(:scount)) !count samples in current layer
        		
        		!if not penetrated, set the boundary to the upper boundary of the lower layer (or some kind of nan?)
!        		if (ctrue == 0) then
!        			bound_depths(inv,layer,bh) = bound_depths(inv,layer+1,bh)
!        			cycle
!        		else !take first sample in new boundary as layer boundary if continuous, or halfway between present and previous depth of first encountered sample within layer
!        			bound_depths(inv,nlayer,bh) = max(minval(pack(depths(:scount),inside(:scount))) - stepback - 1,0) 		!The -1 accounts for fortran/python conversion
!        		end if

!			end do
!		end do
        
        
    
        								

		
	end do
	



end subroutine


