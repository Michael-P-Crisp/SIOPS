module setup_EA
    
    use GA
    use variables
    use setup_SI
    
    implicit none
    
    
    
    contains
    

    subroutine EAprep(soilseeds, & !soil generation variables
                      inv_bh_in,inv_depths_in,inv_test_in,add_errors,use_CI,test_errors,inv_reduction_in,conf_int,percentile_in,s_dev_in,soffset,swidth,sstep &					!site investigation variables
                      ,npdepths,detdisp,plocation,pdepths,ck_set,nrep_mc,rel_loads,load_con,num_loads,preps,buildingweight,difftol,failurevals,costvals,si_performance,pilecost,testcost,testnames,costmetric,abstol &
    ,maxEA,deterministic,max_cons,popsize,perc_error,mut_rate,frac_keep,finaloutput,mut_lims,imut,inv_no,EAseed,siea_mode,num_elite,stopmode,unistart,opt_local, & !EA variables
    femvech,femvecv,prad,usepie,npl2,goodpiles,goodcases,sdist,sumweights,extents,indices,CKheight)  !multiple layer site investigations
    

    
    ! --- soil generation variables ----
    
      integer, intent(in) :: soilseeds(:)
      
      
! ---- evolutionary algorithm variables

        integer, intent(in) :: maxEA !max number of EA generations
        integer, intent(in) :: deterministic !SI mode. 1 = 'deterministic' analysis of investigation performance, 2 = 'heat map', 3 = optimise investigations with EA
        integer, intent(in) :: max_cons     !max number of consequtively equal EA performance before stopping
        integer, intent(in) :: popsize     !population size for evolutionary algorithm (in this context, the number of identical investigations to run)
        real, intent(in) :: perc_error     !percentage error for EA stop critera (out of 100)
        real, intent(in) :: mut_rate       !mutation rate for EA
        real,intent(in) :: mut_lims(2)      ! min and max mutations if dynamic mutation is specified
        integer,intent(in) :: imut !mutation mode; 1 = constant, 2 = adapt based on fitness, 3 = adapt based on proximity in parameter space
        
        real, intent(in) :: frac_keep      !fraction of population to keep after each EA iteration
        integer, intent(in) :: finaloutput !output mode for EA. 1 = final only, 2 = save best of each EA generation, 3 = save full final population, 4 = save everything
        
        integer, intent(in) :: inv_no !the number of the current investigation being optimised
        integer, intent(in) :: EAseed !random seed for evolutionary algorithm (set to 0 for clock-based, or higher for consistent numbers each run)
        integer, intent(in) :: siea_mode !manner of controlling the borehole locations in the EA. 1 Means start across full soil. 2 Means Start within building footprint. 3 Means ALWAYS keep within building footprint.
	
        integer, intent(in) :: num_elite !number of elite population members that are carried over into the next generation, unmutated
        logical, intent(in) :: stopmode !.true. means wait for EA to have converged (remain unchanged for a number of realisations), .false. means wait for a number of iterations past the current GLOBAL optimum before stopping (in case fitness starts increasing).
        logical unistart !true means have random initial population according to uniform distribution, otherwise attempt to have normal offsets from a regular grid
        
! ---- site investigation variables -----


	integer,  intent(in) :: inv_bh_in						!number of boreholes in each investigation
	integer,  intent(in) :: inv_depths_in(:) 				!sampling depth start, end, interval
	integer,  intent(in) :: inv_test_in					!test type for each investigation
	logical,  intent(in) :: add_errors						!whether to add test errors or not
	real,  intent(in) :: test_errors(:,:)	!matrix of test error stats (num tests x num errors)
	integer,  intent(in) :: inv_reduction_in				!reduction method for each test
	real :: evals(nrep_mc,popsize)						!reduced value for each investigation
	real,  intent(in) :: percentile_in,s_dev_in !percentile and standard deviation below mean for reduction methods
    real, intent(in) :: testcost(:)     !cost per metre for each test
    character(4), intent(in) :: testnames(:) !names for each test
    integer, intent(in) :: soffset(2)		!x,y offset of borehole from corner of soil (elements)
    integer, intent(in) :: swidth(2)		!x,y dimensions of site investigation area
    integer, intent(in) :: sstep(2)            !borehole step size in each dimension for the heat map mode
	
	integer :: nbh,bhdepth,ntest	!No. boreholes, borehole depth (elements), No. tests

	
	real :: conf_int !confidence interval for truncating SI sample values 
	logical :: use_ci 
    
!--- multi-layer site investigation variables
    
    integer, intent(in) :: preps(2)		!No. piles in x,y directions - repetitions
    
    real, intent(in) :: femvech(:),femvecv(:) !size of FEM elements in the horizontal and vertical directions
    integer, intent(in) :: prad(2)			!pile width in x,y directions
    logical, intent(in) :: usepie
    real, intent(in) :: sdist(:,:,:),sumweights(:) !distance of elements from each pile, sum of distance weights
    integer, intent(in) :: extents(:,:,:)!extent of soil around pile
    integer :: goodpiles(preps(1)*preps(2))        !vector of piles that have an associated applied load
    integer :: goodcases                               !number of good piles in goodpiles 
    integer :: npl2 !good number of piles under certain conditions
    real CKheight(:,:,:)  !(nrep_MC,nlayer-1,preps(1)*preps(2)) !effective layer depths at each pile associated with the full, original soil
    integer :: indices(:,:,:) !(2,nxew,nyew)
    

! --- pile performance variables
	integer, intent(in) :: npdepths							!number of pile depths to analyse
	integer, intent(in) :: pdepths(:)			!pile depths to analyse (elements)
	real, intent(in) :: detdisp(:)				!deterministic settlement for each pile depth (1 kN applied load. Soil stiffness 1 MPa)
	!integer, intent(in) :: preps(2)		!No. piles in x,y directions - repetitions
	integer, intent(in) :: plocation(:,:)		!Pile x,y coordinates
	
	real, intent(in) :: rel_loads(:)		!set of relative pile loads based on tributary area of building
	integer, intent(in) :: load_con(:)			!connectivity vector saying which pile has which load
	
	real, intent(in) :: buildingweight,difftol
	
	real, intent(in) :: ck_set(:,:,:) !true ck pile settlement for each realisation, pile and depth
    
    integer, intent(in) :: nrep_mc !number of MC realisations
    integer, intent(in) :: num_loads !number of load cases
    
    
    real, intent(in) :: failurevals(2), costvals(2) !2 x,y coordinates defining the line to interpolate costs from differential settlement
    real, intent(in) :: pilecost !cost for pile, per metre depth
    
    real, intent(out) :: si_performance !performance of the current optimal investigation used as the objective function
    real :: opt_stats(5),opt_stats2(5) !a set of all performance metrics for the optimal investigation
    integer, intent(in) :: costmetric !Performance metric to use. -1 = failure cost, 0 = probability of failure, 1 or more = average differntial settlement to the power of the costmetric value (higher values more heavily penalise excessive values)
    real, intent(in) :: abstol !absolute settlement tolerance  (mm). Positive values are used directly. Negative values make it a function of differential settlement * minimum pile spacing (i.e., automatically determined).
        
    
    
    !------local variables-----
    
    integer :: inv_bh(popsize)						!number of boreholes in each investigation
	integer	:: inv_coords(popsize,inv_bh_in,2)			!x,y coordinates of boreholes for each investigation
    integer :: mininv_coords(inv_bh_in,2)			!x,y coordinates associated with the minimum cost
    real :: mincost,mincost2     !minimum cost from the GA
	integer :: inv_depths(inv_bh_in,popsize,3) 				!sampling depth start, end, interval
	integer	:: inv_test(inv_bh_in,popsize)					!test type for each investigation
	integer :: inv_reduction(popsize)				!reduction method for each test
	real :: percentile(popsize),s_dev(popsize) !percentile and standard deviation below mean for reduction method
    real randu !random number function
    character(1000) str2
    character(100) str3
    real smin(2),sside(2) !offset and range for site investigation
    real lo, hi !lower and upper limits for the investigations in the x and y directions
    real temp_arr2(inv_bh_in,2) !array of random values for initialising population
    real temp_arr3(2)
    integer samevals !keep track of identical members of the population
    integer :: numvar !number of variables per population member (no. BH x 2)
    integer :: startpop !start at 2nd population member if a regular grid can be achieved, otherwise start at 1.
    integer iga,iga2 !generation counter
    real propbadmc !proportion of bad monte carlo realisations for best population member
    
    real :: fulltime,fulltime2 !time for GA to complete 
    
    real si_bound !proportion of the side to bound inwards. i.e 0.8 means the middle 80% 
    
    !Time and seed related things
    INTEGER::how_big   ! Used in the RANDOM_SEED subroutine
    INTEGER :: vals(8)       ! Contains values from the time/date call
    INTEGER,ALLOCATABLE,DIMENSION(:)::seeds      ! Vector w/ vals for RANDOM_SEED brtn
    integer status
    
    !Variables related to the 2nd-phase local optimisation
    logical, intent(in) :: opt_local !whether to do a quick GA analysis on the final solution to see if a better one is nearby
    real normsd !standard deviation (elements) for the random offsets to the final solution
    
    !loop counters
    integer i,x,y,j,k
    
    
    numvar = 2*inv_bh_in
    normsd = sqrt(2.0)*(swidth(1)+swidth(2))/8 !make the normal standard deviation = sqrt(2)*average investigation area length/4
    
    !_______________________________________________________
! Initialize random number generator
! Some machines may require more care in calling the random number generator
! This program sets a seed randomly from computer clock

  CALL RANDOM_SEED(SIZE=how_big) ! Finds the size of array expected by subroutine
  ALLOCATE(seeds(how_big),STAT=status)
  IF(status/=0) THEN
    WRITE(*,*) '  An error occurred allocating the array ‘seeds’ in the main program.'
  END IF

  if(EAseed > 0) then           
      seeds = [(j*EAseed,j=1,size(seeds))]   !set fixed seed based on input
  else  
      CALL RANDOM_SEED            ! Initializes the seed
      CALL DATE_AND_TIME(VALUES=vals)   !These values depend on the current time
      IF(vals(8)==0) THEN               ! We only want a non-zero value
        vals(8)=vals(5)+vals(6)+vals(7) ! Substitute in the case of zero (HH+MM+SS)
      END IF
      CALL RANDOM_SEED(GET=seeds) ! Gets the seed
      seeds=seeds*vals(8)         ! Adjusts seed so it is different each time
    
  end if
  CALL RANDOM_SEED(PUT=seeds) ! Seeds the random number generator
  DEALLOCATE(seeds)

    
    
    
    !Create a set of site investigation options, mostly copied from the single input case
    inv_bh = inv_bh_in
    do i = 1,3
        inv_depths(:,:,i) = inv_depths_in(i) !all boreholes in all investigations have same depth
    end do
    inv_test = inv_test_in
    inv_reduction = inv_reduction_in
    percentile = percentile_in
    s_dev = s_dev_in
    
    inv_coords = -1 !set invalid default
    
    !Get lower bound of investigation for the initial conditions, as well as the investigation widths
    if(siea_mode == 1) then !based on the soil width
        smin = 1
        sside(1) = nxew-1
        sside(2) = nyew-1
    else
        smin = soffset+1    !based on the building width
        sside = smin+(swidth-1)
    end if
    
    !same as above, except this is for after the EA starts
    if (siea_mode < 3) then
        lo = 1
        hi = min(nxew,nyew)
    else
        lo = minval(soffset)
        hi = minval(soffset)+maxval(swidth)
    end if
    
    
    
    if( mod(inv_bh_in,goodcases) == 0 .and. .true.) then
        !if the number of boreholes is a multiple of the number of piles, concentrate them at the piles if specified.
        call bh_at_piles(plocation(goodpiles(:goodcases),:),prad,inv_bh_in,inv_coords(1,:,:),startpop)
    else
        !Set the first population member as having regular spacing over the building footprint
        call getgrid_SI(soffset,swidth,inv_bh_in,inv_coords(1,:,:),startpop)
    end if
    
    !randomly generate borehole coordinates for the remaining population  
    if(unistart .or. popsize == 1 .or. startpop == 1) then
        
        do i = startpop, popsize
    	    call RANDOM_NUMBER(temp_arr2) !get random values for current population (between 0 and 1)
    	    temp_arr2(:,1) = (sside(1)-smin(1))*temp_arr2(:,1)+smin(1) !scale them to the limits
		    temp_arr2(:,2) = (sside(2)-smin(1))*temp_arr2(:,2)+smin(2) !scale them to the limits
		    inv_coords(i,:,:) = nint(temp_arr2)	  !fill them in
            
        end do
        
        
    else !create an initial population based on random offsets from a regular grid, if one exists
        do i = 2,popsize !add random offsets to all population members except the first one, which shall remain unchanged in case it's optimal
            
            call vnorm( temp_arr2, numvar ) !get random values with normal distribution N(0,1)
            inv_coords(i,:,:) = inv_coords(1,:,:) + nint(temp_arr2*normsd)
        end do


        !check that the new coordinates are within the limits
        do i = 1,2
            where(inv_coords(:,:,i) < smin(i)) inv_coords(:,:,i) = smin(i)
            where(inv_coords(:,:,i) > sside(i)) inv_coords(:,:,i) = sside(i)
        end do
        

        
    end if
    
    !save initial population if debugging
    if(finaloutput == 5) then
        write(str2,'(A,I0,A)')'si_results/Initial_population-',inv_no,'.txt'
        open(2934,file=str2)
        write(2934,'(100000(I4,X))') (j,j=1,inv_bh(1)),(j,j=1,inv_bh(1))
        do i = 1,popsize !add random offsets to all population members except the first one, which shall remain unchanged in case it's optimal
            write(2934,'(100000(I4,X))') inv_coords(i,:,1),inv_coords(i,:,2)
        end do
        close(2934)
    end if
    

    !Make sure all borehole locations are valid with respect to intra- and inter-investigation similarity
    call force_valid_SI(inv_coords,inv_bh_in,popsize,smin,sside)
    

    !Output evolution
    if(finaloutput == 2 .or. finaloutput == 4) then
        if (singletrue) then
            write(str2,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'si_results/EvolutionEA-coords_Inv-',inv_no,'_piles-',goodcases,'_BArea-',nint(buildingarea),'_BH-',inv_bh(1),'_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy,'.txt'
        else
            write(str2,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'si_results/EvolutionEA-coords_Inv-',inv_no,'_piles-',goodcases,'_BArea-',nint(buildingarea),'_BH-',inv_bh(1),'_nlayers-',nlayer,'_ratio2l-',nint(lmean_ave(1)/minval(lmean_ave(:2))),'+',nint(lmean_ave(2)/minval(lmean_ave(:2))),'_depth2l-',nint(dz*ldepths(1)),'_bSD-',nint(bsd*dz),'.txt'
        end if
        open(503,file=str2) !,BUFFERED='no')
        write(503,'(A14,X,A10,X,1000000(I8))') 'Fitness','%_Bad_Reps',(j,j=1,inv_bh(1)),(j,j=1,inv_bh(1))
    end if

	!----perform evolutionary algorithm----
	
	write(*,'(A6,2X,A8,3X,A8,2X,A8,X,A)') 'EA Gen', 'Min Cost', 'Mutation', 'Time (s)', 'Bad Reps (%)'
  	write(*,'(A)') '---------------------------------------'

    call my_cga(soilseeds, & !soil generation variables
                      popsize,inv_bh_in,inv_bh,inv_depths,inv_test,add_errors,use_CI,test_errors,inv_reduction,conf_int,percentile,s_dev,inv_coords,mininv_coords,mincost &					!site investigation variables
                      ,npdepths,detdisp,plocation,pdepths,ck_set,nrep_mc,rel_loads,load_con,num_loads,preps,buildingweight,difftol,failurevals,costvals,si_performance,pilecost,testcost,testnames,costmetric,abstol &
    ,maxEA,deterministic,max_cons,popsize,perc_error,mut_rate,frac_keep,finaloutput,numvar,lo,hi,mut_lims,imut,inv_no,num_elite,stopmode,iga,propbadmc,fulltime, & !EA variables
    femvech,femvecv,prad,usepie,npl2,goodpiles,goodcases,sdist,sumweights,extents,indices,CKheight,sstep,opt_stats) !multiple layer site investigations
    
    
    !Create slight variations of the fittest population member, and do a quick optimisation of those.
    !Assumes that the final borehole location set from the above  optimisation is in the near vicinity of the true global optimisation
    !The uniform mutation system used above isn't very good for finding the exact global optima.
    !This works by applying random normal offsets to ALL coordinate copies.
    
    iga2 = 0
    fulltime2= 0 
    if(opt_local) then
        inv_coords(1,:,:) = mininv_coords !Copies the fittest set 
        do i = 2,popsize !add random offsets to all population members except the first one, which shall remain unchanged in case it's optimal
            
            call vnorm( temp_arr2, numvar ) !get random values with normal distribution N(0,1)
            inv_coords(i,:,:) = mininv_coords + nint(temp_arr2*normsd)
            
        end do

        !check that the new coordinates are within the limits
        do i = 1,2
            where(inv_coords(:,:,i) < smin(i)) inv_coords(:,:,i) = smin(i)
            where(inv_coords(:,:,i) > sside(i)) inv_coords(:,:,i) = sside(i)
        end do
        
        
        !Make sure all borehole locations are valid with respect to intra- and inter-investigation similarity
        call force_valid_SI(inv_coords,inv_bh_in,popsize,smin,sside)
        
        write(*,*) "Starting phase 2"
        
        call my_cga(soilseeds, & !soil generation variables
                      popsize,inv_bh_in,inv_bh,inv_depths,inv_test,add_errors,use_CI,test_errors,inv_reduction,conf_int,percentile,s_dev,inv_coords,mininv_coords,mincost2 &					!site investigation variables
                      ,npdepths,detdisp,plocation,pdepths,ck_set,nrep_mc,rel_loads,load_con,num_loads,preps,buildingweight,difftol,failurevals,costvals,si_performance,pilecost,testcost,testnames,costmetric,abstol &
    ,maxEA,deterministic,max_cons/2,popsize,perc_error,mut_rate,frac_keep,finaloutput,numvar,lo,hi,mut_lims,imut,inv_no,num_elite,stopmode,iga2,propbadmc,fulltime2, & !EA variables
    femvech,femvecv,prad,usepie,npl2,goodpiles,goodcases,sdist,sumweights,extents,indices,CKheight,sstep,opt_stats2) !multiple layer site investigations
        
        !ensure that the true optimal case is saved. Shouldn't really need this, especially if elitism is used, but it's good to be sure
        if(mincost2 < mincost) then
            mincost = mincost2
            opt_stats = opt_stats2
        end if
        
    end if
    
    
      !--output best case to file--
  !Site investigation attributes and final cost
  write(500,'(I6,2X,A4,2X,A3,3X,F9.2,2X,F10.6,X,I8,X,F8.1,X,E15.8,X,E15.8,X,E15.8,X,E15.8,X,E15.8)') inv_bh_in,testnames(inv_test_in),rednames(inv_reduction_in),dz*inv_depths_in(2),propbadmc,iga+iga2,fulltime+fulltime2,opt_stats
  !x coords
  write(501,'(1000000(F8.2))') (dz*mininv_coords(j,1),j=1,inv_bh_in)
  !y coords
  write(502,'(1000000(F8.2))') (dz*mininv_coords(j,2),j=1,inv_bh_in)
  
  !make sure all files are written to
  flush(500) 
  flush(501)
  flush(502)
    
    if(finaloutput == 2 .or. finaloutput == 4) then
        close(503)
    end if
    
    end subroutine
    
    
    
    !This subroutine checks that all investigations are valid. In other words, that no two investigations
    !in a population are identical (all boreholes in the same location), or that any two boreholes within a
    !single investigation are at the same location; this latter point will cause a crash in the multi-layer
    !mode due to the Delauney triangulation used.
    !If any such cases are detected, then the offending borehole(s) are replaced by random values.
    
    !Note that the first check does not account for order. In other words, two investigations may have all boreholes
    !in the same location, but the rank/order of those boreholes is different. This is unlikely to occur, at least
    !in the initial population (and probably within the GA itself).
    ! It's not a critical problem if two investigations are the same, especially if the occurances are small compared
    ! to the population size.
    
    subroutine force_valid_SI(inv_coords,inv_bh_in,popsize,smin,sside)

    real,intent(in) :: smin(2),sside(2)
    integer, intent(in) :: inv_bh_in, popsize !number of boreholes and population size
    integer, intent(inout) :: inv_coords(:,:,:)
    
    real temp_arr2(inv_bh_in,2) !array of random values for initialising population
    real temp_arr3(2)
    integer samevals
    integer i,j,k
    
        !Loop through the inputs to ensure that no two members of the population are identical
    samevals = 1
    do while (samevals > 0) !Keep looping until all members are unique (because there's no guarentee that )
    	samevals = 0		!This number should remain 0 if all population members are unique
		do i = 1,popsize-1	!Compare each member with every other member
			do j = i+1,popsize 
				if(all(inv_coords(i,:,:)==inv_coords(j,:,:))) then
                    samevals = samevals + 1
                    call RANDOM_NUMBER(temp_arr2) !get replacement random values
					temp_arr2(:,1) = (sside(1))*temp_arr2(:,1)+smin(1) !scale them to the limits
					temp_arr2(:,2) = (sside(2))*temp_arr2(:,2)+smin(2) !scale them to the limits
					inv_coords(j,:,:) = nint(temp_arr2)	  !fill them in
				end if
			end do
		end do
    end do    
    
    if(inv_bh_in > 1) then !Make sure that no boreholes WITHIN an investigation overlap - this will cause a problem for the multi-layer triangulation
        samevals = 1
        do while (samevals > 0) !Keep looping until all members are unique (because there's no guarentee that )
    	    samevals = 0		!This number should remain 0 if all population members are unique
		    do i = 1,popsize	!Compare each member with every other member
			    do j = 1,inv_bh_in-1
                    do k = j+1,inv_bh_in
				        if(inv_coords(i,j,1) == inv_coords(i,k,1) .and. inv_coords(i,j,2) == inv_coords(i,k,2) ) then
                            samevals = samevals + 1
					        call RANDOM_NUMBER(temp_arr3) !get replacement random values
					        temp_arr3 = sside*temp_arr3+smin !scale them to the limits
					        inv_coords(i,k,:) = nint(temp_arr3)	  !fill them in
                        end if
                    end do
			    end do
		    end do
        end do    
    end if
    
    
    end subroutine
                      
    end module
