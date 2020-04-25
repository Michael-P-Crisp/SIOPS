module GA 
    
use variables
    
    contains
    
    
    !Notes: 
    !Currently there is no check for coincident boreholes within an investigation, however this is unlikely to occur. Rare instances are probably ok. 
    
!=======================================================
subroutine my_cga(soilseeds, & !soil generation variables
                      ninv,nbh,inv_bh,inv_depths,inv_test,add_errors,use_CI,test_errors,inv_reduction,conf_int,percentile,s_dev,inv_coords,mininv_coords,mincost &					!site investigation variables
                      ,npdepths,detdisp,plocation,pdepths,ck_set,nrep_mc,rel_loads,load_con,num_loads,preps,buildingweight,difftol,failurevals,costvals,si_performance,pilecost,testcost,testnames,costmetric,abstol &
    ,maxit,deterministic,max_same,popsize,tol,mutrate_init,selection,finaloutput,npar,lo,hi,mut_lims,imut,datafolder,inv_no,num_elite,stopmode,iga,propbadmc,fulltime, & !EA variables
    femvech,femvecv,prad,usepie,npl2,goodpiles,goodcases,sdist,sumweights,extents,indices,CKheight,sstep,opt_stats) !multiple layer site investigations
! Main genetic algorithm program for Haupt and Haupt
! 2003 - uses High Performance Fortran
! Purpose: Optimizes variables based on a cost
! function using a genetic algorithm. Based on pseudocode in
! Haupt and Haupt, 1998
!
! Date        Programmer     Description of Changes 
! =====       ============   ==========================
! 3 July 2003 Jaymon Knight  Code based on pseudocode
! 19 Nov 2003 Sue Haupt      Revised for 2nd ed of Haupt and Haupt
! June 2017   Nikola Mirkov  Adapt to "ordinary" fortran from HPF version
! May 2019    Michael Crisp  Customised GA for site investigation problem
!                            Properly implemented the stop criterion as well
!                               as adding a new one.
!                            Added an adaptive mutation rate (use with caution)
!                            Implemented Elitism.
!                            Fixed a bug so less than half the population can be
!                               kept for mating.
!                            Added various outputs related to the program.
!
!
  !USE funct
  USE qsort_c_module
  USE process_SI
  USE process_SI_multi
  USE variables

  IMPLICIT NONE

  ! Define parameters
  ! Define GA parameters
  ! Use these to tune the code to your problem

  INTEGER, intent(in)::maxit     ! Maximum number of iterations
  INTEGER, intent(in)::max_same  ! Maximum# of consecutively equal vals
  INTEGER, intent(in)::popsize  ! Size of population
  INTEGER, intent(in)::npar !2*nbh       ! Number of parameters
  REAL, intent(in)::tol         ! Percent error for stop criteria
  REAL, intent(in) ::mutrate_init     ! Initial Mutation rate
  real mutrate					!mutation rate
  REAL, intent(in)::selection   ! Fraction of population to keep
  Real, intent(in)::lo !1       ! Minimum parameter value
  REAL, intent(in)::hi !min(nxew,nyew) !10.       ! Maximum parameter value
  integer, intent(in) :: inv_no !the number of the current investigation being optimised
  integer, intent(in) :: num_elite !number of elite population members that are carried over into the next generation, unmutated
  logical, intent(in) :: stopmode !.true. means wait for EA to have converged (remain unchanged for a number of realisations), .false. means wait for a number of iterations past the current GLOBAL optimum before stopping (in case fitness starts increasing).
  real :: gentime !time to do an entire generation
  real, intent(out) :: fulltime !time for GA to complete 
  
  integer soffset(2),swidth(2) !placeholder variables
  integer, intent(in) :: sstep(2)            !borehole step size in each dimension for the heat map mode

  ! Define variables
  INTEGER::status    ! Error flag
  INTEGER::keep      ! Number kept from each generation
  INTEGER::nM         ! Number of matings
  INTEGER::nmut      ! Total number of mutations
  INTEGER::iga       ! Generation counter
  INTEGER::i,j       ! Indices
  INTEGER::same      ! Counter for consecutively equal values
  INTEGER::gsame     ! Counter for number of iterations since current best global value reached
  INTEGER::bad_sort  ! Counts number of bad sorts from hpf grade_up - maybe not needed now?
  REAL::minc         ! Minimum cost
  REAL::temp         ! Temporary variable
  REAL::xy           ! Mix from ma and pa
  
  !dynamic mutation rate variables (copied this from an old implementation of the GA PIKAIA implementation, not sure how effective it is).
    real, intent(in) ::           mut_lims(2)   !min and max mutation rate
    integer, intent(in) ::        imut  !mutation mode; 1=constant, 2=based on fitness, 3=based on proximity in parameter space
    real              rdif                                  !difference metric between best and median population members
    real,parameter ::   rdiflo=0.05, rdifhi=0.25, delta=1.5   !max and min difference metric thresholds, and the rate of mutation change once threshold is reached
  
  integer,allocatable :: invcount(:) !number of valid values for goodinv

  integer ij,iy,iz,iloop,k !loop counters
  real cr !timing variables
  integer allfinish,allstart,start,finish !timing variables
  integer parsum !sum of differences in parameters between any two population sizes
  character(1000), intent(in) :: datafolder !directory of data information
  character(1000) str2
  integer samevals !variable to ensure there are no identical population members

  ! Define matrix variables
  
  INTEGER,ALLOCATABLE,DIMENSION(:)::ind2       ! For sorting mutated population
  INTEGER,ALLOCATABLE,DIMENSION(:)::ma,pa      ! Parents (indices)


  INTEGER,ALLOCATABLE,DIMENSION(:)::xp         ! Crossover point
  INTEGER,ALLOCATABLE,DIMENSION(:)::ix         ! Index of mate #1
  INTEGER,ALLOCATABLE,DIMENSION(:)::mrow,mcol  ! Used for sorting mutations

  REAL,ALLOCATABLE,DIMENSION(:)::odds,odds_tmp ! Involved in pairing
  REAL,ALLOCATABLE,DIMENSION(:)::pick1,pick2   ! Mates one and two
  REAL,ALLOCATABLE,DIMENSION(:)::temp_arr_1    ! Temporary 1-d array
  REAL,ALLOCATABLE,DIMENSION(:)::r             ! Mixing parameter
  
  REAL::par(popsize,npar),par2(popsize,npar)   ! Matrix of population values
  REAL::cost(popsize)         ! Cost function evaluated
  INTEGER ::ind(popsize)        ! Sorted indices from cost function
  real temp_arr2(npar)
  
  
  
  real mincost !minimum cost for the entire process
  real mincost2 !check for no improvement past current optimum
  real propbadmc !proportion of bad monte carlo realisations for best population member
  integer mininv_coords(nbh,2)			!x,y coordinates associated with the minimum cost
  
  !-----------------------   VARIABLES TO PASS THROUGH TO THE SI PERFORMANCE SUBROUTINE ------
  
  integer, intent(in) :: deterministic !SI mode. 1 = 'deterministic' analysis of investigation performance, 2 = 'heat map', 3 = optimise investigations with EA
  integer, intent(in) :: finaloutput !output mode for EA. 1 = final only, 2 = save best of each EA generation, 3 = save full final population, 4 = save everything
    
  
    integer, intent(in) :: inv_bh(popsize)						!number of boreholes in each investigation
	integer	:: inv_coords(popsize,nbh,2)			!x,y coordinates of boreholes for each investigation
	integer, intent(in) :: inv_depths(nbh,popsize,3) 				!sampling depth start, end, interval
	integer, intent(in)	:: inv_test(nbh,popsize)					!test type for each investigation
	integer, intent(in) :: inv_reduction(popsize)				!reduction method for each test
	real, intent(in) :: percentile(popsize),s_dev(popsize) !percentile and standard deviation below mean for reduction method
    integer, intent(in) :: nbh,ninv
  
      ! --- soil generation variables ----
    
  
      integer, intent(in) :: soilseeds(:)

      
! ---- site investigation variables -----



	logical,  intent(in) :: add_errors						!whether to add test errors or not
	real,  intent(in) :: test_errors(:,:)	!matrix of test error stats (num tests x num errors)
	real :: evals(nrep_mc,popsize)						!reduced value for each investigation
    real, intent(in) :: testcost(:)     !cost per metre for each test
    character(4), intent(in) :: testnames(:) !names for each test
	
	integer :: bhdepth,ntest	!No. boreholes, borehole depth (elements), No. tests

	
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
    real CKheight(nrep_MC,nlayer-1,preps(1)*preps(2)) !effective layer depths at each pile associated with the full, original soil
    integer :: indices(2,nxew,nyew)
    
    

! --- pile performance variables
	integer, intent(in) :: npdepths							!number of pile depths to analyse
	integer, intent(in) :: pdepths(:)			!pile depths to analyse (elements)
	real, intent(in) :: detdisp(:)				!deterministic settlement for each pile depth (1 kN applied load. Soil stiffness 1 MPa)

	integer, intent(in) :: plocation(:,:)		!Pile x,y coordinates
	
	real, intent(in) :: rel_loads(:)		!set of relative pile loads based on tributary area of building
	integer, intent(in) :: load_con(:)			!connectivity vector saying which pile has which load
	
	real, intent(in) :: buildingweight,difftol
	
	real, intent(in) :: ck_set(:,:,:) !true ck pile settlement for each realisation, pile and depth
    
    integer, intent(in) :: nrep_mc !number of MC realisations
    integer, intent(in) :: num_loads !number of load cases
    
    
    real, intent(in) :: failurevals(2), costvals(2) !2 x,y coordinates defining the line to interpolate costs from differential settlement
    real, intent(in) :: pilecost !cost for pile, per metre depth
    
    real, intent(out) :: si_performance, opt_stats(5) !performance of the current optimal investigation used as the objective function, and a set of all performance metrics for the optimal investigation
    real :: fcost(popsize), pcost(popsize), probfail(popsize), avediff(popsize), diffgeo(popsize) !failure cost, pile cost, prob. failure, ave. diff. set, geometric statistic
    integer, intent(in) :: costmetric !Which objective function to use for optimisation; -1 = cost, 0 = probability of failure, 1 or more = average differential settlement weighted to the power of costmetric
    real, intent(in) :: abstol !absolute settlement tolerance  (mm). Positive values are used directly. Negative values make it a function of differential settlement * minimum pile spacing (i.e., automatically determined).
        
    logical check
    real temp_arr3(2)
  !-----------------------
  

  ! Calculate variables
  mutrate = mutrate_init !set the mutation rate to the initial value 
  keep=FLOOR(selection*popsize)                ! Number to keep from each generation
  nM=CEILING(REAL(popsize-keep)/2.)             ! Number of matings
  

  !Allocate arrays (block 1)
  ALLOCATE(odds(keep+1),odds_tmp(keep+1), &
           ma(nM),pa(nM),pick1(nM),pick2(nM), &
           r(nM),xp(nM),ix(CEILING(REAL(popsize-keep)/2.)),STAT=status)    !ix(CEILING(REAL(keep)/2.)),STAT=status)
  IF(status/=0) THEN
    WRITE(*,*) 'Error allocating arrays in main program.'
    STOP
  END IF
  
  allocate(invcount(popsize))

!_______________________________________________________
! Create the initial population, evaluate costs, sort
  
  !Already created!

  !CALL RANDOM_NUMBER(par) ! Fills par matrix w/ random numbers 
  !par=(hi-lo)*par+lo      ! between hi & lo !Normalizes values

!_______________________________________________________
! Start generations

  iga=0
  minc=0.
  same=0
  gsame=0
  bad_sort=0
  fulltime = 0
  mincost = huge(mincost) !make the mincost the largest possible value to start off
  mincost2 = mincost

  !OPEN(UNIT=10,FILE='data.dat',STATUS='REPLACE',ACTION='WRITE',IOSTAT=status)
  !IF(status/=0) THEN
  !  WRITE(*,*)"Error opening file 'out.dat'."
  !END IF
  
  
  !pack the investigation coordinates into the parameter array
  do iz = 1,popsize
      iloop = 0
      do iy = 1,nbh
          do ij = 1,2
              iloop = iloop + 1 
              par(iz,iloop)  = inv_coords(iz,iy,ij)
          end do
      end do
  end do
  
  !timing information
  CALL system_clock(count_rate=cr)
  



  
  !----------- PERFORM EVOLUTIONARY ALGORITHM GENERATIONS ----------- 
  
  DO WHILE(iga<maxit)
      iga=iga+1          ! Increment counter
      

      
    
      
    !get start time
    CALL system_clock(allstart)
    


    !CALL ff(par,cost)  ! <<<< Calculates cost using function ff <<<<
    
 !repack the parameter value array into the inv_coords shape    
    !count potentially use a single reshape statement instead, but this guarentees the order is correct and consistent both ways.
      do iz = 1,popsize
          iloop = 0
          do iy = 1,nbh
              do ij = 1,2
                  iloop = iloop + 1 
                  inv_coords(iz,iy,ij) = nint(par(iz,iloop))
              end do
          end do
      end do
      
    
    ! ----- GET FITNESS OF THE POPULATION! -----
    !CALL system_clock(start)
    if (singletrue) then
        call get_si_perf(soilseeds, & !soil generation variables
                          ninv,nbh,1,1,1,inv_bh,inv_coords,inv_depths,inv_test,add_errors,use_CI,test_errors,inv_reduction,conf_int,percentile,s_dev,soffset,swidth,sstep  &					!site investigation variables
                          ,npdepths,detdisp,plocation,pdepths,ck_set,nrep_mc,rel_loads,load_con,num_loads,preps,buildingweight,difftol,failurevals,costvals,cost,pilecost,testcost,testnames,costmetric,inv_no, &
                          deterministic,finaloutput,abstol,datafolder,invcount,goodpiles,goodcases,fcost, pcost, probfail, avediff, diffgeo)
    else
        call get_si_perf_multi(soilseeds, & !soil generation variables
                    ninv,nbh,1,1,1,inv_bh,inv_coords,inv_depths,inv_test,add_errors,use_CI,test_errors,inv_reduction,conf_int,percentile,s_dev,soffset,swidth,sstep  &					!site investigation variables
                    ,npdepths,detdisp,plocation,pdepths,ck_set,nrep_mc,rel_loads,load_con,num_loads,preps,buildingweight,difftol,failurevals,costvals,cost,pilecost,testcost,testnames,costmetric,inv_no, &
                    deterministic,abstol,femvech,femvecv,prad,datafolder,usepie,npl2,goodpiles,goodcases,sdist,sumweights,extents,indices,CKheight,finaloutput,invcount,fcost, pcost, probfail, avediff, diffgeo)
                                
        
    end if
     !CALL system_clock(finish)
    !--------------------------------------------
    

    ! Arange (1,popsize) and store in ind
    do i=1,popsize
      ind(i) = i
    enddo


    ! Quicksort array and return sorted array and new index order relative to old one
    ! Min cost in element 1,
    ! Cost in order stored in ind

    
    call QsortCIndx(cost,ind)

    ! WRITE(*,*) minc,cost(1),iga
    
    !save the 
    if (cost(1)<mincost) then 
    	mincost = cost(1) !save new lowest cost
    	mininv_coords = inv_coords(ind(1),:,:)
        propbadmc = 100*real(nrep_MC-invcount(ind(1)))/nrep_MC
        
        opt_stats(1) = fcost(ind(1))
        opt_stats(2) = pcost(ind(1))
        opt_stats(3) = probfail(ind(1))
        opt_stats(4) = avediff(ind(1))
        opt_stats(5) = diffgeo(ind(1)) !diffgeo(ind(1))
    end if
        
    if ((mincost2-cost(1))/cost(1)<=tol/100.) then !we've found a better (or equal) solution!
        gsame = gsame + 1
    else
        gsame = 0
        mincost2 = cost(1)
    end if
    
    IF(ABS((cost(1)-minc)/cost(1))<=tol/100.) THEN
      same=same+1
    ELSE
      same=0
    END IF
    
    minc=cost(1) !save current cost for next generation calculation
    par=par(ind,:) ! Puts par in the order stored in ind.

    
    
    
    !Some outputs
    if(finaloutput == 2 .or. finaloutput >= 4) then
     write(503,'(E14.6,X,F10.6,X,1000000(F8.2))') cost(1),100*real(nrep_MC-invcount(ind(1)))/nrep_MC,(dz*inv_coords(ind(1),j,1),j=1,inv_bh(ind(1))),(dz*inv_coords(ind(1),j,2),j=1,inv_bh(ind(1)))
     flush(503)
    end if
    
    
    if((same == max_same .and. stopmode) .or. (gsame == max_same .and. .not. stopmode)) then
        write(*,*) 'stopping criteria reached'
        exit
    end if


    
    
    !_______________________________________________________
    ! Update mutation rate if neccessary
    
!     dynamical adjustment of mutation rate;
!        imut=2 : adjustment based on fitness differential between best and median individuals
!        imut=3 : adjustment based on metric distance between best and median individuals
    
    if( imut > 1 ) then

      if(imut.eq.2.) then
!     Adjustment based on fitness differential 
         rdif=abs(cost(1)-cost(popsize/2))/(cost(1)+cost(popsize/2))
      else if(imut.eq.3.)then
!     Adjustment based on normalized metric distance
         rdif=0.
         do i=1,npar
            rdif=rdif+( par(1,i)-par(popsize/2,i) )**2
         end do
         rdif=sqrt( rdif ) / real(npar)
      endif

      if(rdif.le.rdiflo)then
         mutrate=min(mut_lims(2),mutrate*delta)
      else if(rdif.ge.rdifhi)then
         mutrate=max(mut_lims(1),mutrate/delta)
      endif
      
    end if
    nmut=CEILING((popsize-1-num_elite)*npar*mutrate)       ! Number of mutations
 
    !_______________________________________________________
    ! Pair chromosomes and produce offspring
    odds(1)=0. !first spot is zero

    ! Fills remainder of odds matrix w/ values keep-1
    DO i=1,keep
      odds(i+1)=keep+1-i
    END DO

    !# HPF function
    !# odds(2:keep+1)=SUM_PREFIX(odds(2:keep+1))

    ! Instead of SUM_PREFIX function
    do i=2,keep+1
      odds_tmp(i) = sum(odds(2:i))
    enddo
    odds(2:keep+1)=odds_tmp(2:keep+1) ! weights chromosomes based upon position in the list

    temp=odds(keep+1)
    odds(2:keep+1)=odds(2:keep+1)/temp

    ! Probability distribution function
    CALL RANDOM_NUMBER(pick1) !mate #1
    CALL RANDOM_NUMBER(pick2) !mate #2

    ! ma and pa contain the indices of the chromosomes that will mate
    DO i=1,nM
      DO j=2,keep+1
        IF(pick1(i)<=odds(j) .AND. pick1(i)>odds(j-1)) THEN
          ma(i)=j-1
        END IF
        IF(pick2(i)<=odds(j) .AND. pick2(i)>odds(j-1)) THEN
          pa(i)=j-1
        END IF
      END DO
    END DO
    

    !_______________________________________________________
    ! Performs mating using single point crossover

    i=0

    DO i=1,nM !CEILING(REAL(keep)/2.)
      ix(i)=2*i-1
    END DO

    !Allocate temporary array (block 2) (Subroutine requires a real argument)
    ALLOCATE(temp_arr_1(nM),STAT=status)
    IF(status/=0) THEN
      WRITE(*,*) 'Error allocating the arrays of allocation block 2 of the main program.'
      STOP
    END IF

    CALL RANDOM_NUMBER(temp_arr_1)

    xp=CEILING(temp_arr_1*REAL(npar))

    DEALLOCATE(temp_arr_1)

    CALL RANDOM_NUMBER(r)

    par2=par
    

    DO i=1,nM
      xy=par2(ma(i),xp(i))-par2(pa(i),xp(i))             ! mix from ma & pa
      par2(keep+ix(i),:)=par2(ma(i),:)                   ! first offspring variable
      par2(keep+ix(i)+1,:)=par2(pa(i),:)                 ! second offspring variable
      par2(keep+ix(i),xp(i))=par2(ma(i),xp(i))-r(i)*xy   ! first offspring variable
      par2(keep+ix(i)+1,xp(i))=par2(pa(i),xp(i))+r(i)*xy ! second offspring variable

      ! Perform crossover when last variable not selected
      IF(xp(i)<npar) THEN                                
        DO j=1,xp(i)
          par2(keep+ix(i),j)=par2(keep+ix(i),j)
          par2(keep+ix(i)+1,j)=par2(keep+ix(i)+1,j)
        END DO

        DO j=xp(i)+1,npar
          par2(keep+ix(i),j)=par2(keep+ix(i)+1,j)
          par2(keep+ix(i)+1,j)=par2(keep+ix(i),j)
        END DO

      END IF
    END DO

    par=par2
    

    !_______________________________________________________
    ! Mutate the population

    ! Allocate arrays (block 3)
    ALLOCATE(temp_arr_1(nmut),mrow(nmut),mcol(nmut),ind2(nmut),STAT=status)

    IF(status/=0) THEN
      WRITE(*,*) 'Error allocating the arrays of allocation block 3 of the main program.'
      STOP
    END IF

    !get members of population to mutate
    CALL RANDOM_NUMBER(temp_arr_1)
    mrow=CEILING(temp_arr_1*(popsize-1-num_elite))+1+num_elite
    
    call QsortCInt(mrow)

    CALL RANDOM_NUMBER(temp_arr_1)
    mcol=CEILING(temp_arr_1*npar)

    CALL RANDOM_NUMBER(temp_arr_1)
    temp_arr_1=nint((hi-lo)*temp_arr_1+lo) ! Normalizes values between hi & lo


    DO i=1,nmut
      par(mrow(i),mcol(i))=temp_arr_1(i)
    END DO

    

    !IF(MINVAL(cost)/=cost(1) .or. minval(mrow) == 1) THEN
    !  bad_sort=bad_sort+1
    !  IF(bad_sort<=1) THEN
    !    !WRITE(10,*)cost
    !  END IF
    !END IF
    
    DEALLOCATE(mrow,mcol,temp_arr_1,ind2)
    
    
        !Loop through the inputs to ensure that no two members of the population are identical
    !samevals = 1
    !do while (samevals > 0) !Keep looping until all members are unique (because there's no guarentee that )
    !	samevals = 0	
        do i = 1,popsize-1
            do j = i+1,popsize 
                if(all(nint(par(i,:))==nint(par(j,:)))) then
    !                samevals = samevals + 1
                    call RANDOM_NUMBER(temp_arr2)
                    temp_arr2 = nint((hi-lo)*temp_arr2+lo)
                    par(j,:) = temp_arr2
                end if
            end do
        end do
    !end do
        
    if(nbh > 1) then !Make sure that no boreholes WITHIN an investigation overlap - this will cause a problem for the multi-layer triangulation
        samevals = 1
        do while (samevals > 0) !Keep looping until all members are unique (because there's no guarentee that )
    	    samevals = 0		!This number should remain 0 if all population members are unique
		    do i = 1,popsize	!Compare each member with every other member
			    do j = 1,npar-2,2
                    do k = j+2,npar,2
                        if( nint(par(i,j)) == nint(par(i,k)) .and. nint(par(i,j+1)) == nint(par(i,k+1)) ) then
                            if (i==1) then
                                write(*,*)
                            end if
                            !write(*,*) 'hereo'
                            samevals = samevals + 1
					        call RANDOM_NUMBER(temp_arr3) !get replacement random values
					        !temp_arr3 = sside*temp_arr3+smin !scale them to the limits
                            temp_arr3 = nint((hi-lo)*temp_arr3+lo)
					        par(i,k:k+1) = nint(temp_arr3)	  !fill them in
                        end if
                    end do
			    end do
		    end do
        end do    
    end if
    
    
    
    
    !make sure the parameters are still within a reasonable margin, by moving them to the nearest valid point
    par = nint(par)
    par2 = par
    where(par < lo) par = lo
    where(par > hi) par = hi
    
    if(.not. all(nint(par) == nint(par))) then
        write(*,*)
    end if

    

    

    !Timing procedures
    CALL system_clock(allfinish)
    gentime = real(allfinish-allstart)/cr
    fulltime = fulltime + gentime
    write(*,'(I6,X,E11.5,X,F8.6,X,F8.1,X,I0)') iga,mincost,mutrate,gentime,nint(100*real(nrep_MC-invcount(ind(1)))/nrep_MC)

  END DO

! GA FINISHED!
!_______________________________________________________

  DEALLOCATE(odds,pick1,pick2,ma,pa,r,xp,ix)
  
  !CLOSE(10)
  
  !--output best case to file--
  !Site investigation attributes and final cost
  !write(500,'(I6,2X,A4,2X,A3,3X,F9.2,2X,F10.6,X,I8,X,E14.7)') inv_bh(ind(1)),testnames(inv_test(ind(1))),rednames(inv_reduction(ind(1))),dz*inv_depths(ind(1),2),propbadmc,iga,mincost
  !x coords
  !write(501,'(1000000(F8.2))') (dz*mininv_coords(j,1),j=1,inv_bh(ind(1)))
  !y coords
  !write(502,'(1000000(F8.2))') (dz*mininv_coords(j,2),j=1,inv_bh(ind(1)))
  


  deallocate(invcount)
  WRITE(*,"(I4,' iterations were required to obtain ',I4,' consecutive values of ',G18.8)") iga,max_same,mincost

    END SUBROUTINE

End module