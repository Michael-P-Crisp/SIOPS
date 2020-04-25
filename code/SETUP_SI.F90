module setup_SI
    
implicit none

contains 

subroutine deterministic_SI(soffset,swidth,nbhmax,in_tests,in_depths,in_reductions,bhnums,testerrors,sampfreq,in_sdev,in_percentile, &		!input variables
		inv_coords,inv_depths,inv_nbh,inv_test,ninv,inv_reduction,percentile,s_dev,si_performance)	 !output variables

	!define the site investigation information for a deterministic analysis. 
	!Hardcoded to do a combination of some tests for each number of boreholes and reduction method
	!Note that prime numbers of boreholes greater than 3 are not considered, and are removed.
	!Boreholes are evenly spread out over the site investigation area in a regular grid layout.


		!input
		integer, intent(in) :: soffset(2)		!x,y offset of borehole from corner of soil (elements)
		integer, intent(in) :: swidth(2)		!x,y dimensions of site investigation area
		!integer, intent(in) :: bhdepth,ntest	!No. boreholes, borehole depth (elements), No. tests
		real, allocatable, intent(in) :: testerrors(:,:) !Transformation error, bias error, measurement error for each test
		integer, allocatable, intent(in) :: sampfreq(:) !Sampling frequency for each test (elements)
		real, allocatable :: si_performance(:)
		integer, intent(in) :: in_tests(:),in_depths(:) !input list of tests, depths and reduction methods to investigate in analysis
		character(2), intent(in) :: in_reductions(:)		!input list of reduction methods
		integer, intent(in) :: bhnums(:) !A list of the number of boreholes in each investigation,
		real, intent(in) :: in_sdev,in_percentile
		
		!output
		integer, allocatable, intent(out) :: inv_coords(:,:,:)			!x,y coordinates of boreholes for each investigation
		integer, allocatable, intent(out) :: inv_depths(:,:,:) 
		integer, allocatable, intent(out) :: inv_nbh(:)						!number of boreholes in each investigation
		integer, allocatable, intent(out) 	:: inv_test(:,:)					!test type for each investigation
        integer, allocatable, intent(out) 	:: inv_reduction(:)					!reduction method for each investigation
        real, allocatable, intent(out) 	:: percentile(:)					!percentile for use in reduction method
        real, allocatable, intent(out) 	:: s_dev(:)					!geometric standard deviation below geometric mean in reduction method
		integer, intent(out) 	:: ninv					!total number of investigations
        integer, intent(out) :: nbhmax !max number of boreholes
		
		!local variables
		integer i,x,y,counter,bh,test,inv,j,red,depth			!loop counters
        integer :: n_red    !technically 6 reduction methods, although let's not use 'minimum' due to it being excessively conservative. So 5 at most.
		integer :: n_depths !number of depths
		integer, allocatable :: bhnums_noprime(:) !A list of the number of boreholes in each investigation without prime numbers.
		logical, allocatable :: bhnums_griddable(:) !true for boreholes numbers that can be placed on a grid
		integer :: nbhact !number of boreholes minus the prime numbers, 
		integer numprime !number of prime numbers in borehole list
		integer factors(2) !spacing in the x,y directions between boreholes, as a factor of the site investigation area dimensions
		integer inn
		integer nbh,n_test
		character(2), parameter :: rednames(5) = (/ 'SA','GA','HA','1Q','SD'/) !reduction method names; hard-coded in SI module
		
	
			! ------------get site investigation coordinates and sampling depths ------------
		
		nbh = size(bhnums)
		n_test = size(in_tests)
		n_red = size(in_reductions)
		n_depths = size(in_depths)
		
		allocate(bhnums_griddable(nbh))

		!---remove prime numbers from the list (since they can't be regularly spaced on a grid)---
		! find the griddable numbers 
		bhnums_griddable =.false. !Assume it's not griddable for now
		do i = 1,nbh
			if (bhnums(i) <= 3 .or. bhnums(i) == 5) then !we'll keep 1, 2, 3 and 5 boreholes as a special case
				bhnums_griddable(i) = .true.
				cycle
			end if 
			do j = 2,bhnums(i)-1
				if (mod(bhnums(i),j) == 0) then	!As soon as a single factor is found, stop the loop
					bhnums_griddable(i) = .true.
					exit
				end if
			end do
        end do
		
		!now remove the prime numbers
		nbhact = count(bhnums_griddable)
		allocate(bhnums_noprime(nbhact))
		bhnums_noprime = pack(bhnums,bhnums_griddable)
        nbhmax = maxval(bhnums_noprime)
		
		! ----------------------------------------------
		
		
		!total number of valid investigations
		ninv = nbhact*n_test*n_red*n_depths
		
		!allocate site investigation arrays
		allocate(inv_coords(ninv,nbhmax,2))		!x,y coordinates for each borehole in each investigation
		allocate(inv_depths(nbhmax,ninv,3)) 			!sample depth information for each investigation
		allocate(inv_nbh(ninv))					!number of boreholes in each investigation
		allocate(inv_test(nbhmax,ninv))				!test type for each investigation
        allocate(inv_reduction(ninv))			!reduction for each investigation
        allocate(percentile(ninv))
        allocate(s_dev(ninv))
		allocate(si_performance(ninv))
        
        !set percentile and standard deviation values in those reduction methods. Constant across all realisations for now.
        percentile = in_percentile
        s_dev = in_sdev
        
		
		inv = 0 !keep track of current investigation
        do red = 1,n_red
            do test= 1,n_test
            	do depth = 1,n_depths
                    
					!do the rest of the boreholes
					do i=1,nbhact
						inv = inv+1
						inv_depths(1,inv,:) = [1,in_depths(depth),sampfreq(in_tests(test))]	!sampling depth information
                        do j = 1,3
                            inv_depths(:,inv,j) = inv_depths(1,inv,j)
                        end do
						inv_nbh(inv) = bhnums_noprime(i)								!number of boreholes
						inv_test(:,inv) = in_tests(test)							!test for investigation
						!get reduction method number from test name
						do j = 1,size(rednames)
							if( in_reductions(red) == rednames(j) ) then
								inv_reduction(inv) = j
								exit
							end if
						end do
						!x,y coordinates for each borehole in each investigation
						if (bhnums_noprime(i) == 1) then !if a single borehole, place in the centre
							inv_coords(inv,1,1) = swidth(1)/2 + soffset(1) - 1
							inv_coords(inv,1,2) = swidth(2)/2 + soffset(2) - 1
						else if(bhnums_noprime(i) == 2) then !if 2 boreholes, place at opposite corners
							inv_coords(inv,1,1) = soffset(1)
							inv_coords(inv,1,2) = soffset(2)
							inv_coords(inv,2,1) = swidth(1) + soffset(1) - 1
							inv_coords(inv,2,2) = swidth(2) + soffset(2)  - 1     
						else if(bhnums_noprime(i) == 3) then !if 3 boreholes, place on 3 of the corners
							inv_coords(inv,1,1) = soffset(1) 
							inv_coords(inv,1,2) = soffset(2) 
							inv_coords(inv,2,1) = swidth(1) + soffset(1) - 1
							inv_coords(inv,2,2) = swidth(2) + soffset(2)  - 1 
							inv_coords(inv,3,1) = soffset(1)
							inv_coords(inv,3,2) = swidth(2) + soffset(2) - 1
                            !inv_coords(inv,3,1) = swidth(1)/2 + soffset(1) - 1 !actually, place the third in the centre
							!inv_coords(inv,3,2) = swidth(2)/2 + soffset(2) - 1
                        else if(bhnums_noprime(i) == 5) then !if 5 boreholes, place 4 at the corners and 1 centrally
                            inv_coords(inv,1,1) = soffset(1)
							inv_coords(inv,1,2) = soffset(2)
							inv_coords(inv,2,1) = swidth(1) + soffset(1) - 1
							inv_coords(inv,2,2) = swidth(2) + soffset(2) - 1
							inv_coords(inv,3,1) = soffset(1)
							inv_coords(inv,3,2) = swidth(2) + soffset(2) - 1
                            inv_coords(inv,4,1) = swidth(1) + soffset(1) - 1
                            inv_coords(inv,4,2) = soffset(2) 
                            inv_coords(inv,5,1) = swidth(1)/2 + soffset(1) - 1
							inv_coords(inv,5,2) = swidth(2)/2 + soffset(2) - 1
						else !spread as a grid over the full area
							!need to find largest small factor to make the grid most square-like
							do j = 1,floor(sqrt(real(bhnums_noprime(i)))) !largest possible small factor must be the square root
								if ( mod(bhnums_noprime(i),j) == 0) then
									factors(1) = bhnums_noprime(i)/j
									factors(2) = bhnums_noprime(i)/factors(1)
								end if
							end do
							counter = 0
							do y=1,factors(2)    !get coordinates based on spacing
								do x=1,factors(1)
									counter = counter + 1
									inv_coords(inv,counter,1) = (x-1) * (swidth(1)-1)/(factors(1)-1) + soffset(1) 
									inv_coords(inv,counter,2) = (y-1) * (swidth(2)-1)/(factors(2)-1) + soffset(2) 
								end do
                            end do
						end if
					
					end do
                end do
            end do
        end do
        
		
		
end subroutine
        
subroutine EA_SI(soffset,swidth,nbhmax,in_tests,in_depths,in_reductions,bhnums,testerrors,sampfreq,in_sdev,in_percentile, &		!input variables
		inv_coords,inv_depths,inv_nbh,inv_test,ninv,inv_reduction,percentile,s_dev,si_performance)	 !output variables

	!define the site investigation information for input to the evolutionary algorithm. 
	!Hardcoded to do a combination of some tests for each number of boreholes and reduction method.
	!Recommend keeping everything constant besides the number of boreholes. However, you may try combinations of other parameters for a sensitivity analysis.
	!Borehole locations are to be randomly assigned later in the evolutionary algorithm setup.

		!input
		integer, intent(in) :: soffset(2)		!x,y offset of borehole from corner of soil (elements)
		integer, intent(in) :: swidth(2)		!x,y dimensions of site investigation area
		!integer, intent(in) :: bhdepth,ntest	!No. boreholes, borehole depth (elements), No. tests
		real, allocatable, intent(in) :: testerrors(:,:) !Transformation error, bias error, measurement error for each test
		integer, allocatable, intent(in) :: sampfreq(:) !Sampling frequency for each test (elements)
		real, allocatable :: si_performance(:)
		integer, intent(in) :: in_tests(:),in_depths(:) !input list of tests, depths and reduction methods to investigate in analysis
		character(2), intent(in) :: in_reductions(:)		!input list of reduction methods
		integer, intent(in) :: bhnums(:) !A list of the number of boreholes in each investigation,
		real, intent(in) :: in_sdev,in_percentile
		
		!output
		integer, allocatable, intent(out) :: inv_coords(:,:,:)			!x,y coordinates of boreholes for each investigation
		integer, allocatable, intent(out) :: inv_depths(:,:,:) 
		integer, allocatable, intent(out) :: inv_nbh(:)						!number of boreholes in each investigation
		integer, allocatable, intent(out) 	:: inv_test(:,:)					!test type for each investigation
        integer, allocatable, intent(out) 	:: inv_reduction(:)					!reduction method for each investigation
        real, allocatable, intent(out) 	:: percentile(:)					!percentile for use in reduction method
        real, allocatable, intent(out) 	:: s_dev(:)					!geometric standard deviation below geometric mean in reduction method
		integer, intent(out) 	:: ninv					!total number of investigations
        integer, intent(out) :: nbhmax !max number of boreholes
		
		!local variables
		integer i,x,y,counter,bh,test,inv,j,red,depth			!loop counters
        integer :: n_red    !technically 6 reduction methods, although let's not use 'minimum' due to it being excessively conservative. So 5 at most.
		integer :: n_depths !number of depths
		integer nbh,n_test
		character(2), parameter :: rednames(5) = (/ 'SA','GA','HA','1Q','SD' /) !reduction method names; hard-coded in SI module
		
	
			! ------------get site investigation coordinates and sampling depths ------------
		
		nbh = size(bhnums)
		n_test = size(in_tests)
		n_red = size(in_reductions)
		n_depths = size(in_depths)
        nbhmax = maxval(bhnums)
		
		! ----------------------------------------------
		
		!total number of valid investigations
		ninv = nbh*n_test*n_red*n_depths
		
		!allocate site investigation arrays
		allocate(inv_coords(ninv,nbh,2))		!x,y coordinates for each borehole in each investigation
		allocate(inv_depths(nbh,ninv,3)) 			!sample depth information for each investigation
		allocate(inv_nbh(ninv))					!number of boreholes in each investigation
		allocate(inv_test(nbh,ninv))				!test type for each investigation
        allocate(inv_reduction(ninv))			!reduction for each investigation
        allocate(percentile(ninv))
        allocate(s_dev(ninv))
		allocate(si_performance(ninv))
        
        !set percentile and standard deviation values in those reduction methods. Constant across all realisations for now.
        percentile = in_percentile
        s_dev = in_sdev
		
		
		inv = 0 !keep track of current investigation
        do red = 1,n_red
            do test= 1,n_test
            	do depth = 1,n_depths
                    
					!do the rest of the boreholes
					do i=1,nbh
						inv = inv+1
						inv_depths(1,inv,:) = [1,in_depths(depth),sampfreq(in_tests(test))]	!sampling depth information
                        do j = 1,3
                            inv_depths(:,inv,j) = inv_depths(1,inv,j)
                        end do
						inv_nbh(inv) = bhnums(i)								!number of boreholes
						inv_test(:,inv) = in_tests(test)							!test for investigation
						!get reduction method number from test name
						do j = 1,size(rednames)
							if( in_reductions(red) == rednames(j) ) then
								inv_reduction(inv) = j
								exit
							end if
						end do
					end do
                end do
            end do
        end do
        
		
		
        end subroutine
        
        
        
        
subroutine getgrid_SI(soffset,swidth,bhnum,inv_coords,startpop)	 !output variables


	!Boreholes are evenly spread out over the site investigation area in a regular grid layout.
    !This is done for the first population member for the genetic algorithm.
    !If the boreholes can't fit on a grid, and don't fit into one of the special cases, 
        !restort to default random locations done in the calling subroutine
    !Note that if 

    !difference between this and the first subroutine is that this only applies to one borehole


		!input
		integer, intent(in) :: soffset(2)		!x,y offset of borehole from corner of soil (elements)
		integer, intent(in) :: swidth(2)		!x,y dimensions of site investigation area
		!integer, intent(in) :: bhdepth,ntest	!No. boreholes, borehole depth (elements), No. tests
		integer, intent(in) :: bhnum !A list of the number of boreholes in each investigation,
		
		!output
		integer, intent(out) :: inv_coords(:,:)			!x,y coordinates of boreholes for each investigation
        integer, intent(out) :: startpop !start at 2nd population member if a regular grid can be achieved, otherwise start at 1.
		
		!local variables
		integer i,x,y,counter,bh,j		!loop counters
		logical :: bhnums_griddable !true for boreholes numbers that can be placed on a grid
		integer numprime !number of prime numbers in borehole list
		integer factors(2) !spacing in the x,y directions between boreholes, as a factor of the site investigation area dimensions
		integer inn
		
	

		! find the griddable numbers 
		bhnums_griddable =.false. !Assume it's not griddable for now
			if (bhnum <= 3 .or. bhnum == 5) then !we'll keep 1, 2, 3 and 5 boreholes as a special case
				bhnums_griddable = .true.
			end if 
			do j = 2,bhnum-1
				if (mod(bhnum,j) == 0) then	!As soon as a single factor is found, stop the loop
					bhnums_griddable = .true.
					exit
				end if
            end do

		
		!assume grid can be achieved for now
            startpop = 2
		

			!x,y coordinates for each borehole in each investigation
			if (bhnum == 1) then !if a single borehole, place in the centre
				inv_coords(1,1) = swidth(1)/2 + soffset(1) - 1
				inv_coords(1,2) = swidth(2)/2 + soffset(2) - 1
			else if(bhnum == 2) then !if 2 boreholes, place at opposite corners
				inv_coords(1,1) = soffset(1)
				inv_coords(1,2) = soffset(2)
				inv_coords(2,1) = swidth(1) + soffset(1) - 1
				inv_coords(2,2) = swidth(2) + soffset(2)  - 1     
			else if(bhnum == 3) then !if 3 boreholes, place on 3 of the corners
				inv_coords(1,1) = soffset(1)
				inv_coords(1,2) = soffset(2)
				inv_coords(2,1) = swidth(1) + soffset(1) - 1
				inv_coords(2,2) = swidth(2) + soffset(2)  - 1 
				inv_coords(3,1) = soffset(1)
				inv_coords(3,2) = swidth(2) + soffset(2) - 1
            else if(bhnum == 5) then !if 5 boreholes, place 4 at the corners and 1 centrally
                inv_coords(1,1) = soffset(1)
				inv_coords(1,2) = soffset(2)
				inv_coords(2,1) = swidth(1) + soffset(1) - 1
				inv_coords(2,2) = swidth(2) + soffset(2)  - 1 
				inv_coords(3,1) = soffset(1)
				inv_coords(3,2) = swidth(2) + soffset(2) - 1
                inv_coords(4,1) = swidth(1) + soffset(1) - 1
                inv_coords(4,2) = soffset(2) 
                inv_coords(5,1) = swidth(1)/2 + soffset(1) - 1
				inv_coords(5,2) = swidth(2)/2 + soffset(2) - 1
            else !spread as a grid over the full area
                if(bhnums_griddable) then

				    !need to find largest small factor to make the grid most square-like
				    do j = 1,floor(sqrt(real(bhnum))) !largest possible small factor must be the square root
					    if ( mod(bhnum,j) == 0) then
						    factors(1) = bhnum/j
						    factors(2) = bhnum/factors(1)
					    end if
				    end do
				    counter = 0
				    do y=1,factors(2)    !get coordinates based on spacing
					    do x=1,factors(1)
						    counter = counter + 1
						    inv_coords(counter,1) = (x-1) * (swidth(1)-1)/(factors(1)-1) + soffset(1) 
                            inv_coords(counter,2) = (y-1) * (swidth(2)-1)/(factors(2)-1) + soffset(2)
					    end do
                    end do
                else !boreholes aren't gridable. Resort to random.
                    startpop = 1
                end if
			end if
	
		
		
end subroutine


subroutine bh_at_piles(plocation,prad,bhnum,inv_coords,startpop)	 !output variables


	!if the current number of boreholes in the GA is equal to a multiple of the number of piles, place the boreholes
    !at the pile locations. This may create some boreholes occupying the same physical space, but that will be fixed at a later time


		!input
        integer, intent(in) :: plocation(:,:)		!Pile x,y coordinates in elements
        integer, intent(in) :: prad(2)          !pile x,y dimensions in elements
		!integer, intent(in) :: bhdepth,ntest	!No. boreholes, borehole depth (elements), No. tests
		integer, intent(in) :: bhnum !A list of the number of boreholes in each investigation,
		
		!output
		integer, intent(out) :: inv_coords(:,:)			!x,y coordinates of boreholes for each investigation
        integer, intent(out) :: startpop !start at 2nd population member if a regular grid can be achieved, otherwise start at 1.
		
		!local variables
		integer i,x,y,counter,counter2,bh,p		!loop counters
		integer bhfactor !ratio of boreholes to piles
		integer inn
		
	
        bhfactor = bhnum/size(plocation,1)

		
		!assume grid can be achieved for now
            startpop = 2
            counter2 = 0
		    do p = 1,size(plocation,1)
                counter = 0
                loop: do !attempt to spread the borehole locations evenly over the elements within the pile interior
				    do y=plocation(p,2),plocation(p,2) + prad(2) - 1    !get coordinates based on spacing
					    do x=plocation(p,1),plocation(p,1) + prad(1) - 1
						    counter = counter + 1
                            if (counter > bhfactor) then !if we've reached the bh-per-pile limit, exit
                                exit loop
                            end if
                            counter2 = counter2 + 1
						    inv_coords(counter2,1) = x
                            inv_coords(counter2,2) = y
					    end do
                    end do
                end do loop
            end do
	
		
		
end subroutine


!return coordinates for a single borehole on every single element over the soil. Can be used to generate a 'heat map' of site investigation performance
!This subroutine is looped in the main program, and will create investigations for heat maps for all combinations of borehole depth, test and reduction method
subroutine full_spread(soffset,swidth,sstep,in_tests,in_depths,in_reductions,testerrors,sampfreq,in_sdev,in_percentile, &		!input variables
		inv_coords,inv_depths,inv_nbh,inv_test,ninv,inv_reduction,percentile,s_dev,si_performance,iter,invcount,fcost, pcost, probfail, avediff, diffgeo)	
        

		!input
		integer, intent(in) :: soffset(2)		!x,y offset of borehole from corner of soil (elements)
		integer, intent(in) :: swidth(2)		!x,y dimensions of site investigation area
        integer, intent(in) :: sstep(2)            !borehole step size in each dimension for the heat map mode
		!integer, intent(in) :: bhdepth,ntest	!No. boreholes, borehole depth (elements), No. tests
		real, allocatable, intent(in) :: testerrors(:,:) !Transformation error, bias error, measurement error for each test
		integer, allocatable, intent(in) :: sampfreq(:) !Sampling frequency for each test (elements)
		real, allocatable :: si_performance(:)
		integer, intent(in) :: in_tests(:),in_depths(:) !input list of tests, depths and reduction methods to investigate in analysis
		character(2), intent(in) :: in_reductions(:)		!input list of reduction methods
		real, intent(in) :: in_sdev,in_percentile
        integer, intent(in) :: iter
		
		!output
		integer, allocatable, intent(out) :: inv_coords(:,:,:)			!x,y coordinates of boreholes for each investigation
		integer, allocatable, intent(out) :: inv_depths(:,:,:) 
		integer, allocatable, intent(out) :: inv_nbh(:)						!number of boreholes in each investigation
		integer, allocatable, intent(out) 	:: inv_test(:,:)					!test type for each investigation
        integer, allocatable, intent(out) 	:: inv_reduction(:)					!reduction method for each investigation
        real, allocatable, intent(out) 	:: percentile(:)					!percentile for use in reduction method
        real, allocatable, intent(out) 	:: s_dev(:)					!geometric standard deviation below geometric mean in reduction method
		integer, intent(out) 	:: ninv					!total number of investigations
        integer, allocatable ,intent(out) :: invcount(:)
        real, allocatable :: fcost(:), pcost(:), probfail(:), avediff(:), diffgeo(:) !failure cost, pile cost, prob. failure, ave. diff. set, geometric statistic
		
		!local variables
		integer i,x,y,counter,bh,test,inv,j,red,depth			!loop counters
        integer :: n_red    !technically 6 reduction methods, although let's not use 'minimum' due to it being excessively conservative. So 5 at most.
		integer :: n_depths !number of depths
		integer nbh,n_test
		character(2), parameter :: rednames(6) = (/ 'SA','GA','HA','1Q','SD','MN' /) !reduction method names; hard-coded in SI module
        integer i_red,i_test,i_depth !indices for the current parameters
	
	
			! ------------get site investigation coordinates and sampling depths ------------
		
        !get some parameter sizes
		n_test = size(in_tests)
		n_red = size(in_reductions)
		n_depths = size(in_depths)
        
        
        
        !get current parameters based on iteration
        inv = 0
        firstloop: do red = 1,n_red
            do test= 1,n_test
            	do depth = 1,n_depths
                    inv = inv + 1
                    if(inv == iter) then
                        i_test = test
                        i_red = red
                        i_depth = depth
                        exit firstloop
                    end if
                end do
            end do
        end do firstloop
		
		! ----------------------------------------------
		
		!total number of valid investigations
		ninv = ((swidth(1)-soffset(1))/sstep(1)+1) * ((swidth(2)-soffset(2))/sstep(2)+1)
		
		!allocate site investigation arrays
		allocate(inv_coords(ninv,1,2))		!x,y coordinates for each borehole in each investigation
		allocate(inv_depths(1,ninv,3)) 			!sample depth information for each investigation
		allocate(inv_nbh(ninv))					!number of boreholes in each investigation
		allocate(inv_test(1,ninv))				!test type for each investigation
        allocate(inv_reduction(ninv))			!reduction for each investigation
        allocate(percentile(ninv))
        allocate(s_dev(ninv))
        allocate(invcount(ninv))
         !only allocate the variables if they haven't been allocated already. Assumes that the number of boreholes doesn't change, which is currently the case. 
        !the above variables seem to be automatically de-allocated since they're specified as intent(out)
        if(.not. allocated(si_performance)) then
            allocate(fcost(ninv), pcost(ninv), probfail(ninv), avediff(ninv), diffgeo(ninv))
            allocate(si_performance(ninv))
        end if
        
        !These parameters are constant across all realisations for now.
        percentile = in_percentile
        s_dev = in_sdev
		inv_nbh = 1
        inv_test = in_tests(i_test)
        !get reduction method number from test name
		do j = 1,size(rednames)
			if( in_reductions(i_red) == rednames(j) ) then
				inv_reduction = j
				exit
			end if
		end do
		
		inv = 0 !keep track of current investigation
		!Set borehole coordinates and test types
		do i=soffset(1),swidth(1),sstep(1) !x coordinates
            do j = soffset(2),swidth(2),sstep(2) !y coordinates
                inv= inv+1
                inv_coords(inv,1,1) = i
                inv_coords(inv,1,2) = j
                inv_depths(1,inv,:) = [1,in_depths(i_depth),sampfreq(in_tests(i_test))]	!sampling depth information
                do x = 1,3
                    inv_depths(:,inv,x) = inv_depths(1,inv,x)
                end do
            end do
        end do


    end subroutine
    
    
        
end module