
    
module SI

contains
    
subroutine s_inv(nxe,nye,nze,efld,n_inv,n_bh,inv_bh,inv_coords,inv_depths,num_tests,num_errors,inv_test,add_errors,use_CI,test_errors,inv_reduction,evals,kseed,CI,percentile,s_dev)

! This subroutine does site investigations of a 3D soil field.
! Options include adding test errors (bias, random, transformation) (zero and negative values are discarded)
! It can also remove values that don't conform to a geometric confidence interval of a user-specified level (recommend 99%)

! Note that all units are in terms of soil elements. If you are working in metres, please convert before hand.

! Subroutine returns the reduced value of each investigation, according to either:
!             1 - standard arithmetic average
!             2 - geometric average
!             3 - harmonic average
!             4 - first quartile
!			  5 - standard deviation below mean (geometric/lognormal)
!             6 - minimum
! Generally speaking, the reduction methods are in order of increasing conservatism, although this isn't strict.

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
! Set the seed to zero to use the system clock as an initial seed.

! Note that the 1Q reduction technique is relatively slow as it requires first sorting the samples in 
! ascending order. Partial sorting up to the desired percentile is used to increase speed, and in some cases
! only values below the mean are used, which helps for larger sets.

! Note due to the complexity and shear amount of inputs, no input safety checks have been implemented.
! Essentially, just make sure that all values are >= 0 (besides kseed and percentile) and it should be fine.

! See "QUANTIFYING THE RISK OF GEOTECHNICAL SITE INVESTIGATIONS" - 2006 PhD Thesis by J. Goldsworthy for details.

! Michael Crisp February 2018



  implicit none

! -- input variables --
  integer,		intent(in) :: nxe,nye,nze						!number of soil elements in x,y,z dimensions
  real, 		intent(in) :: efld(:,:,:) !efld(nxe,nye,nze) 				!3D array of soil E values
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
  real, 		intent(in) :: test_errors(num_tests,num_errors)	!matrix of test error stats (num tests x num errors)
  integer, 		intent(in) :: inv_reduction(n_inv)				!reduction method for each test
  integer,      intent(in) :: kseed                             !initial random seed
  real,			intent(in) :: CI    							!confidence interval factor
  real,			intent(in) :: percentile(n_inv)					!percentile to use in reduction method (proportion)
  real,			intent(in) :: s_dev(n_inv)						!standard deviations below mean for reduction method (geometric)


! -- output variables --
  real,			intent(out) :: evals(:)						!reduced value for each investigation
  
! -- local variables --

  integer		inv,bh,j,k,i,d                              !loop counters
  integer		depths(nze,n_bh),scount(n_bh)					!vector containing depths of samples, no. samples
  integer       tscount                                     !total sample count in current investigation
  integer       counter

  
  !various arrays for storing samples, and masks of valid values
  
  !real test(10000)
  
  real			sidata(n_bh,nze)							!store SI samples (up to possible maximum)
  real          goodvals(n_bh*nze)                          !store valid values of SI samples in 1D array
  real          goodvals2(n_bh*nze)                         !copy
  real          lgoodvals(n_bh*nze)                         !log values of the above info
  logical       gvmask(n_bh*nze)                            !mask invalid values of the above info

  ! variables related to test errors 
  integer       seed                                        !random seed
  real          rand_vals(n_bh*nze)                         !random values
  real          bhmeans(n_bh)                               !mean of borehole values
  real          Ltest_errors(num_tests,num_errors)          !lognormal test errors
  integer       trans_count(num_tests)
  real          trans_sum(num_tests)
  
  ! variables related to removing invalid values
  integer       ctrue,ctrue2                                !count of true values
  real          gmn,gsd                                     !geometric mean and standard deviation of SI samples
  real          lower,upper                                 !lower and upper bounds of geometric confidence interval
  
  

  
!functions
  integer iseed
  
  
	
	!initialise random number generator
	!seed = iseed(max(0,kseed)) !ensure kseed is at least zero - not needed since it's already initialized
	
	
	do inv = 1,n_inv
        
        
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

        !tscount = inv_bh(inv)*scount !number of samples in current investigation
        
 !       if (inv_bh(inv) == 1) then
 !           open(555,file='sidata1.txt', status="old", position="append", action="write")
 !       else if (inv_bh(inv) == 2) then
 !           open(555,file='sidata2.txt', status="old", position="append", action="write")
 !       else if (inv_bh(inv) == 3) then
 !            open(555,file='sidata3.txt', status="old", position="append", action="write")
 !       else
 !           open(555,file='sidata4.txt', status="old", position="append", action="write")
 !       end if
	
		!Do investigation; extract all relevant samples and repack into a 1D vector of values
		do bh = 1,inv_bh(inv)	!loop through boreholes
			sidata(bh,:scount(bh)) = efld(inv_coords(inv,bh,1),inv_coords(inv,bh,2),depths(:scount(bh),bh))
         !   write(555,'(1000000F8.4,X)') sidata(bh,:scount)
        end do
        
        !close(555)
        
        !inv_coords(:,:,1)
		
		!vectorised equivalent, although need to tweak the ordering of the reshape for it to be correct.
		!Apparently this is slower than the above anyway
        !sidata(:inv_bh(inv),:scount) = reshape(efld(inv_coords(inv,:inv_bh(inv),1),inv_coords(inv,:inv_bh(inv),2),depths(:scount)),[inv_bh(inv),scount])
        
        !-----add random errors to samples-----
        
        if (add_errors) then
        	
        	!convert normal COVs to lognormal parameters
  			Ltest_errors = sqrt(log(1.0 + test_errors ** 2))
        
            !borehole average:
            !bhmeans(:inv_bh(inv)) = sum(sidata(:inv_bh(inv),:scount),1)/scount                         !get borehole means                                          
        
            call vnorm(rand_vals,inv_bh(inv))                                                           !generate zero-mean, unit-variance normal values
            do bh = 1,inv_bh(inv)                                                                       !apply bias to boreholes
                !sidata(bh,:scount) = sidata(bh,:scount) + bhmeans(bh) * ( rand_vals(bh) - 1 )
                rand_vals(bh) = exp(rand_vals(bh) * Ltest_errors(inv_test(bh,inv),1))      !transform to unit-mean lognormal values of desired SD
                sidata(bh,:scount(bh)) = sidata(bh,:scount(bh)) + ( rand_vals(bh) - 1 ) * sum(sidata(bh,:scount(bh)))/scount(bh)
            end do

        

        
        	!apply per-sample random error
            trans_sum = 0
            trans_count = 0
            do bh = 1,inv_bh(inv)
                !random values
                call vnorm(rand_vals,scount(bh))
                rand_vals(:scount(bh)) = exp(rand_vals(:scount(bh)) * Ltest_errors(inv_test(bh,inv),2))
                sidata(bh,:scount(bh)) = sidata(bh,:scount(bh)) * rand_vals(:scount(bh))
                
                !keep track of sum and size to calculate global averages later for global transformation error
                trans_sum(inv_test(bh,inv)) = trans_sum(inv_test(bh,inv)) + sum(sidata(bh,:scount(bh)))
                trans_count(inv_test(bh,inv)) = trans_count(inv_test(bh,inv)) + scount(bh)
            end do
        	!goodvals(:tscount) = pack(sidata(:inv_bh(inv),:scount),.true.) * rand_vals(:tscount)
        	
    
            
            !apply transformation error. Make sure the proportional change is consistent across all tests
            call vnorm(rand_vals,1)
            rand_vals(2:num_tests+1) = exp(rand_vals(1) * Ltest_errors(:,3))
            
            counter = 1
            do bh = 1,inv_bh(inv)
                goodvals(counter:counter+scount(bh)-1) = sidata(bh,:scount(bh)) + (rand_vals(inv_test(bh,inv)+1) - 1)  * trans_sum(inv_test(bh,inv))/trans_count(inv_test(bh,inv))
                counter = counter + scount(bh)
            end do

            ! remove values less than zero
   !         gvmask(:tscount) = goodvals(:tscount) > 0													!mask invalid values
			!ctrue = count(gvmask(:tscount))																!count invalid values
			!goodvals(:ctrue) = pack(goodvals(:tscount),gvmask(:tscount))								!collect
            
            !This code is a slightly faster version of the above, but untested (only about 5% difference overall)
            ctrue = 0
            do i = 1,tscount
                if(goodvals(i) > 0) then
                    ctrue = ctrue + 1
                    goodvals(ctrue) = goodvals(i)
                end if
            end do
        
        else	!otherwise just use the full original data
            
            counter = 0
            do bh = 1,inv_bh(inv)
                do j = 1,scount(bh)
                    counter = counter + 1
                    goodvals(counter) = sidata(bh,j)
                end do
            end do
        
        	!goodvals(:tscount) = pack(sidata(:inv_bh(inv),:scount),.true.) 								!flatten array
        	ctrue = tscount
        
        end if
        
        
        !---- remove invalid values from data ----
        
        if (use_CI) then
        
			! remove values above and below 99% lognormal CI
		
			lgoodvals(:ctrue) = log(goodvals(:ctrue))
		
			gmn = sum(lgoodvals(:ctrue))/ctrue
		
			gsd = 0
			do i = 1,ctrue
				gsd = gsd + (lgoodvals(i)-gmn)**2
			end do
		
			!geometric mean and standard deviation
			gsd = exp(sqrt(gsd/ctrue))
			gmn = exp(gmn)
		
			lower = gmn / (gsd * CI)
			upper = gmn * (gsd * CI)
		
			!remove values outside the CI
			!gvmask(:ctrue) = goodvals(:ctrue) > lower .and. goodvals(:ctrue) < upper
			!ctrue2 = count(gvmask(:ctrue))
			!goodvals(:ctrue2) = pack(goodvals(:ctrue),gvmask(:ctrue))
            
            !This code is a slightly faster version of the above
            ctrue2 = 0
            do i = 1,ctrue
                if(goodvals(i) > lower .and. goodvals(i) < upper) then
                    ctrue2 = ctrue2 + 1
                    goodvals(ctrue2) = goodvals(i)
                end if
            end do
			
			ctrue = ctrue2
			
        end if
        
        !apply nominated reduction method
        
        call reduce_single(evals(inv),inv_reduction(inv),ctrue,goodvals,gvmask,lgoodvals,percentile(inv),s_dev(inv))
        
    end do


    end subroutine
    
    
    !do reduction for single layer soils
    subroutine reduce_single(evals,inv_reduction,ctrue,goodvals,gvmask,lgoodvals,percentile,s_dev)
    
    implicit none
    
    integer, intent(in) :: inv_reduction,ctrue !choice of reduction method, number of samples
    real, intent(in) :: percentile,s_dev !nominated percentile and standard dviation for the latter methods
    real :: goodvals(:)
    !this variables are defined as inputs simply so that they're not allocated in memory upon each call
    real :: lgoodvals(:)
    logical :: gvmask(:)
    
    
    
    !local variables
      !1st quartile variables
      real          Q                                           !approximate index of 1st quartile
      real          temp                                        !temp storage in sorting algorithm
      real          gmn,gsd                                     !geometric mean and standard deviation of SI samples
      integer ctrue2
      integer k,j,i !loop counters
    
    
    real, intent(out) :: evals !reduced value
    
    
            !----- reduce remaining good values ------
        
        if (inv_reduction == 1) then
            evals = sum(goodvals(:ctrue))/ctrue
        else if (inv_reduction == 2) then
            evals = exp(sum(log(goodvals(:ctrue)))/ctrue)
        else if (inv_reduction == 3) then
            evals = 1/(sum(1/goodvals(:ctrue))/ctrue)

        else if (inv_reduction == 4) then
            
            Q = percentile*(ctrue-1)+1 !index of 1st quartile
            
! if the desired percentile is less than 50%, and if the vector is large enough,
! employ the optimisation of extracting and sorting values less than the mean.
! As samples should be lognormally distributed, the median should always be less
! than the mean. If the sample count is less than ~120, the overhead in this extraction
! results in a slowdown. It's also not recommended to use a percentile close to 50%,
! as the reliability cannot be guarenteed when test errors are applied.

            if(percentile < 0.5 .and. ctrue > 120) then
				!Extract values smaller than the median
				!gvmask(:ctrue) = goodvals(:ctrue) < sum(goodvals(:ctrue))/ctrue
				!ctrue2 = count(gvmask(:ctrue))
				!lgoodvals(:ctrue2) = pack(goodvals(:ctrue),gvmask(:ctrue))
                
                ctrue2 = 0
                gmn = sum(goodvals(:ctrue))/ctrue
                do i = 1,ctrue
                    if(goodvals(i) < gmn) then
                        ctrue2 = ctrue2 + 1
                        lgoodvals(ctrue2) = goodvals(i)
                    end if
                end do
			
			else
				ctrue2 = ctrue
                lgoodvals = goodvals
			end if
            
            !need to sort samples
            do k = ctrue2,max(2,ctrue2-ceiling(Q)),-1	!Loop end is normally 2, but this is faster (only relevant bit is sorted)  max(2,ctrue-ceiling(Q))
				do j = ctrue2-1,ctrue2-k+1,-1
					if (lgoodvals(j+1) < lgoodvals(j)) then 
						temp = lgoodvals(j)
						lgoodvals(j) = lgoodvals(j+1)
						lgoodvals(j+1) = temp
					end if
				end do
            end do   
                  
            evals = (lgoodvals(ceiling(Q)) - lgoodvals(floor(Q)))*(percentile*(ctrue-1)-floor(percentile*(ctrue-1))) + lgoodvals(floor(Q))
            
            
        else if (inv_reduction == 5) then
        
        	lgoodvals(:ctrue) = log(goodvals(:ctrue))
		
			gmn = sum(lgoodvals(:ctrue))/ctrue
		
			gsd = 0
			do i = 1,ctrue
				gsd = gsd + (lgoodvals(i)-gmn)**2
			end do
		
			!geometric mean and standard deviation
			gsd = exp(sqrt(gsd/ctrue))
			gmn = exp(gmn)
			
			evals = gmn/(gsd**s_dev)
            
        else if (inv_reduction == 6) then
            evals = minval(goodvals(:ctrue))
            
        else
            write(*,*) "Error: Invalid reduction method chosen for investigation" ,inv_reduction
            write(*,*) "Valid values are:"
            write(*,*) "1 - standard arithmetic average"
            write(*,*) "2 - geometric average"
            write(*,*) "3 - harmonic average"
            write(*,*) "4 - percentile"
            write(*,*) "5 - geometric SD"
            read(*,*)
            stop
        end if

		

    
    end subroutine
    
    !a version of the vnorm subroutine within a module for a safe interface. Original version in vnorm.F90 kept for compatibility with old code
      subroutine vnorm_module( g, n )
      real, intent(out) :: g(:)
      real randu
      data twopi/6.2831853071795864769/
      data two/2.0/

      nn = n/2
      nn = 2*nn
      do i = 1, nn, 2
         a = twopi*randu(0)
         r = sqrt(-two*alog(randu(0)))
         g(i)   = r*cos(a)
         g(i+1) = r*sin(a)
      end do 

      if( n .eq. nn ) return
!                                     for n odd, set the last value
      a = twopi*randu(0)
      r = sqrt(-two*alog(randu(0)))
      g(n) = r*cos(a)

      return
      end

end module
