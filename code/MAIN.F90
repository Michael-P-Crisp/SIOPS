!c  *********************************************************************
!c  *                                                                   *
!c  *                            program soilgen                        *
!c  *                                                                   *
!c  *********************************************************************
!c  Single Precision Version 1.0
!c  Written by Michael P Crisp, adapted from work by Gordon A. Fenton and 
!   D. Vaughan Griffiths

!c  Latest Update: September, 2018
!c
!c  PURPOSE Generate multi-layer soil profiles
!c
!c  DESCRIPTION: This minimalistic program uses the GAF library involving the local 
!   average subdivision (LAS) method to generate lognormally distributed random fields
!   representing virtual soils. See the book "Risk Assessment in Geotechnical Engineering
!   (2008)" for details on the theory involved beyond that stated in these subroutines.
!
!   The program also determines pile settlement within the generated virtual soils, as well as
!   site investigations, and piles designed from tha tsite investigation information.
!   Note that the software must be run in "CK" mode first (only once per soil configuration!)
!    to generate the true settlement values,then subsequently run in "SI" mode to conduct site 
!   investigations, and compare the results to the true settlement values.
!
!   This program has been tested on intel fortran 12.1 under optimisation level /O1 (/O2 and above
!   yielded some incorrect results, possibly a compiler bug for this ~6 year old compiler).
!   It has also been tested under gfortran 6.10 with optimisation -O3 (no bugs detected - good!).
!   Note that intel fortran seems to be faster, even with the minimal /O1 optimisation, however
!   they were not tested on the same machine, so this could be a hardware matter.
!
!   WARNING: This program has defined invalid values through the generation of NaN values by division
!   by zero, and tracked them through the non-standard (?) isnan function. While this works on both
!   compilers, if this ever stops working in the future, the NaN functionally can be replaced by the
!   intrinsic IEEE modules.

!   Note that there are a lot of variables in this program, leading to convoluted subroutine interfaces.
!   One option to simplify things is by using shared module variables for commonly-used variables. This 
!   would clean up interfaces and variable declarations. However I'd advise against it - "if it ain't broke,
!   don't fix it." Not to mention variable declarations make it easier to see what's being used.
!
!   See "Framework for the Optimisation of Site Investigations for Pile Designs in Complex Multi-Layered Soil"
!   Crisp et al. (2018) for detailed description on the functionality of this software. (Note that multi-layer
!   functionality has been stripped out to simplify the program, and because the required FEA to get pile
!   settlement is not practical without a supercomputer. Instead, the PIE method is a sufficient substitute.
!   (see process_CK.F90 for details on pile settlement).
!   The paper "INFLUENCE OF SITE INVESTIGATION BOREHOLE PATTERN AND AREA ON PILE FOUNDATION PERFORMANCE" 
!   Crisp et al. (2018) is also a good reference as a highly succinct document of the overall process.
!
!   In the unlikely event that this program 'crashes', one of two things will happen. Firstly, the program might
!   stop and display an error message (hitting enter will close the program). Alternatively, the program may
!   simply quit. In this case, a text file should be generated in the program's directory explaining the error.
!   Note that while reasonsable steps have been taken to ensure the validity of this software, it is not guarenteed
!   to work correctly under all circumstances, as it may contain unknown bugs. It is suggested that results should
!   be checked to see that they are reasonable and logical given the program inputs.
!c
!c
!c  REVISION HISTORY:
!	v2: Program is now better suited to purely deterministic Monte Carlo analysis as opposed to Evolutionary algorithms
!		Inputs are more intuitive.
!
!	v2.1: Fixed an issue where only a single soil was being saved to disk instead of all of them.
!			Make the system more automatic in knowing when it saves the files to disk. Now, if the soil SOF and COV isn't
!				changed, then it can re-use the saved soils, saving time.
!			Instead of generating soils with a particular mean stiffness, the program now always generates unit soil 
!				stiffness, and instead scales the applied load instead. This gives exactly the same results, but is faster.
!
!	v3.0: Implemented the genetic algorithm designed to optimised borehole locations. It could potentially be adapted to
!			optimise other variables as well, such as test type, borehole depth etc, although the code would need to be
!			modified to put the population values into the proper arrays for fitness assessment, and that they're scaled to
!			the appropriate maximum and minimum values.
!		Added a 'PP' mode (pre-process). This generates the deterministic pile settlement curve as a function of pile length,
!			as well as the associated soil weights for each length, both of which are needed by the pseudo-incremental energy 
!			method (PIE method) uses to determine true pile performance in a variable soil. Previously these were computed by 
!			a seperate program. This information is generated by linear-elastic FEM
!		Added an 'automatic' mode which should theoretically be able to perform all components ('PP','CK','SI) in a row based
!			on which files are available in the data folder. 
!		All data information besides inputs are now in the data folder. The data folder directory can be specified as an input
!			which is useful for running several versions of the program simultaneously with different inputs without
!			requiring the (often quite large) data to be stored in several locations.
!		Added a subroutine to export soil as a CSV.
    
!   v4.0: Added a multi-layer analysis system, however this currently not hooked up to the genetic algorithm. While it is possible
!           to do this, the ML mode is an order of magnitude slower than the SL mode, and so it is less feasible to run.
!           This mode can have undulating layers, however the soil properties within each layer are constant, due to a limitation
!           of the settlement method used. Research has shown that 2 highly constrasting layers could be considered a worst case,
!           as additional layers beyond this do not significantly add to the uncertainty.
    
!   v4.1 Added an option to generate a very large soil and take random subsets for each Monte Carlo realisation for the single layer
!           mode, as opposed to generating a new soil each time, or reading in from the disk. This is notably faster, and recommended.
!       Modified the pile design process to pre-process all settlement curves to increase their resolution by a factor of 10. This is
!           intended to produce 0.1m pile length increments. Instead of designing piles by interpolating curves, the design is 
!           undertaken by finding the discrete locations that bracket the true design, taking the bigger one (essentially rounding up
!           to the nearest 0.1m). This design is given as an array index, which later gives the true pile settlement directly, as 
!           opposed to further interpolation. This is notably faster.
!       Added options for an initial population to be generated using normal offsets from a regular grid.
!       Added an optional 2nd phase evolutionary algorithm, which recreates the population using normal offsets from the fittest
!           member of phase 1. Theoretically, once the intitial phase has slowed down, this means that a global optimum is nearby,
!           and this offset-based repopulation should produce a good chance of finding it quickly.
!c-------------------------------------------------------------------------
program main

	use readinm
	use variables
    use setup_SI
    use process_CK
    use process_SI
    use process_SI_multi
    use setup_EA
    use writesoils
    use checkfiles
    use weights
    use detcheck

	implicit none


	!----------------variables for virtual soil generation (no need to change)------------------
    !note: moved to variables module to be accessed directly from subroutines, as opposed to arguments
    character(100) soildsc !soil description string for saving files
    

	
	!----------------variables for pile settlement information------------------
	integer :: npdepths							!number of pile depths to analyse
	integer,ALLOCATABLE :: pdepths(:)			!pile depths to analyse (elements)
	real, allocatable :: detdisp(:)				!deterministic settlement for each pile depth (1 kN applied load. Soil stiffness 1 MPa)
	real(8) :: cg_tol 					!max no. iterations for FEM subroutine
	integer :: cg_limit 				!tolerance for FEM
	integer :: preps(2)		!No. piles in x,y directions - repetitions
	integer :: prad(2)			!pile width in x,y directions
	integer, allocatable :: plocation(:,:)		!Pile x,y coordinates
	
	integer :: femrad,femdepth !radius and depth of soil around pile (elements)
    integer, allocatable :: startstress(:,:) !boundary for PIE field in terms of number of elements away from the pile in each dimension
	real(8) :: pieval			!effective soil elastic modulus for PIE method
	real, allocatable :: soilweight(:,:,:,:) !Soil weights required for the PIE method of getting CK settlement
	integer :: nels									!total number of soil elements in PIE method
	real :: disp !pile settlement from FEM
    
    real :: costvals(2),failurevals(2) 
    real :: buildingweight !total weight of building
    
    integer, allocatable :: load_con(:) !connectivity array specifying which rel_load is applied to each pile index
	real, allocatable :: rel_loads(:) !relative proportion of pile loads for various load cases
	integer num_loads !number of relative pile load cases
    
    real, allocatable :: femvech(:),femvecv(:) !size of FEM elements in the horizontal and vertical directions

    
    real :: pilecost !cost for pile, per metre depth
    real :: difftol !differential settlement design tolerance
    logical :: usepie !whether to use FEM approximation such as the PIE method for single layers and Mylonakis and Gazetas method for multi-layer soils
    logical :: regmesh !use regular or variable FEM mesh (recommend latter - false)
    real :: abstol !absolute settlement tolerance  (mm). Positive values are used directly. Negative values make it a function of differential settlement * minimum pile spacing (i.e., automatically determined).
        
	
	
	!---------------- site investigation variables--------------------
	
	integer,  allocatable :: inv_nbh(:)						!number of boreholes in each investigation
	integer,  allocatable	:: inv_coords(:,:,:)			!x,y coordinates of boreholes for each investigation
	integer,  allocatable :: inv_depths(:,:,:) 				!sampling depth start, end, interval
	integer,  allocatable	:: inv_test(:,:)				!test type for each investigation
	logical	:: add_errors						!whether to add test errors or not
	integer,  allocatable :: inv_reduction(:)				!reduction method for each test
	real, allocatable :: testerrors(:,:) !Transformation error, bias error, measurement error for each test
	real, allocatable :: percentile(:),s_dev(:) !percentile and standard deviation below mean for reduction methods
	integer, allocatable :: sampfreq(:) !Sampling frequency for each test (elements)
	
	integer :: soffset(2)		!x,y offset of borehole from corner of soil (elements)
	integer :: swidth(2)		!x,y dimensions of site investigation area
    integer :: sstep(2)            !borehole step size in each dimension for the heat map mode
	integer :: nbhmax,bhdepth,ntest	!No. boreholes, borehole depth (elements), No. tests
	integer :: ninv				!number of investigations
    
    real, allocatable :: testcost(:) !cost per metre depth for each test
    character(4), allocatable :: testnames(:) !names for each test
    
    integer, allocatable :: in_tests(:),in_depths(:) !input list of tests, depths and reduction methods to investigate in analysis
	character(2), allocatable :: in_reductions(:)		!input list of reduction methods
	integer, allocatable :: bhnums(:) !A list of the number of boreholes in each investigation,
	real :: in_sdev,in_percentile !standard deviations below the mean in the Geometric "SD" reduction method. Percentile to use for the percentile reduction method. Recommend 1 and 0.25 currently.
	
	real, allocatable :: si_performance(:)
	
	real :: conf_int  !confidence interval for truncating site investigation sample values (2.576 corresponds to 99%)
	
	logical :: use_ci
	integer read_si !If true, will read in the configuration of individual investigations using the 'input/si.txt'. Otherwise will generate a set of all combinations of investigation inputs
	                !0 = don't read, 1 = read attributes and borehole locations, 2 = also read borehole depths and tests on a per-borehole basis

    !effective pile layer 
    real, allocatable :: sdist(:,:,:),sumweights(:) !distance of elements from each pile, sum of distance weights
    integer, allocatable :: extents(:,:,:), indices(:,:,:) !extent of soil around pile, x/y coords of layer boundary indices
    
    integer, allocatable :: goodpiles(:)        !vector of piles that have an associated applied load
    integer :: goodcases                               !number of good piles in goodpiles 
    integer :: npl2 !good number of piles under certain conditions

	real, allocatable :: CKheight(:,:,:) !effective layer depths at each pile associated with the full, original soil
	

	
	
	!---------------- Evolutionary algorithm variables ------------
	
	integer nrep_MC !number of monte carlo realisations 
	integer maxEA !max number of evolutionary algorithm iterations
    integer :: max_cons     !max number of consequtively equal EA performance before stopping
    integer :: popsize     !population size for evolutionary algorithm
    real :: perc_error     !percentage error for EA stop critera (out of 100)
    real :: mut_rate       !mutation rate for EA
    real :: frac_keep      !fraction of population to keep after each EA iteration
    
    real :: mut_lims(2)      ! min and max mutations if dynamic mutation is specified
    integer :: imut !mutation mode; 1 = constant, 2 = adapt based on fitness, 3 = adapt based on proximity in parameter space
	
	integer, allocatable :: soilseeds(:) !unique random seed for each Monte Carlo realisation
    integer :: costmetric !Performance metric to use. -1 = failure cost, 0 = probability of failure, 1 or more = average differntial settlement to the power of the costmetric value (higher values more heavily penalise excessive values)
    
    integer :: deterministic !SI mode. 1 = 'deterministic' analysis of investigation performance, 2 = 'heat map', 3 = optimise investigations with EA
    integer :: finaloutput !output mode for EA. 1 = final only, 2 = save best of each EA generation, 3 = save full final population, 4 = save everything
    
    integer,allocatable :: invcount(:) !number of valid values for goodinv
    integer :: siea_mode !manner of controlling the borehole locations in the EA. 1 Means start across full soil. 2 Means Start within building footprint. 3 Means ALWAYS keep within building footprint.
	
    integer :: num_elite !number of elite population members that are carried over into the next generation, unmutated
    logical :: stopmode !.true. means wait for EA to have converged (remain unchanged for a number of realisations), .false. means wait for a number of iterations past the current GLOBAL optimum before stopping (in case fitness starts increasing).
    
    logical :: unistart !true means have random initial population according to uniform distribution, otherwise attempt to have normal offsets from a regular grid
    logical :: opt_local !whether to do a quick GA analysis on the final solution to see if a better one is nearby
    
    !some site investigation performance metrics
    real, allocatable :: fcost(:), pcost(:), probfail(:), avediff(:), diffgeo(:) !failure cost, pile cost, prob. failure, ave. diff. set, geometric statistic
	
	!---------------- general variables--------------------
	integer EAseed !random seed for evolutionary algorithm (set to 0 for clock-based, or higher for consistent numbers each run)
	integer i,j,ix,iy,iz,outpx,outpy,outpz,L,istat,jr,x,y,z,pd,count !loop counters
	character(1000) str2 !string for output
    character(100) str3
	character(100) rfname !name for file output
	
	real, allocatable :: ck_pset(:,:,:) !ck pile settlement for each realisation, for the piles in the 2D grid, per pile depth
	character(100) :: runmode 	!whether to run the CK or SI modes
	logical :: thismode(3) !whether to do the PP, CK or SI analysis respectively
	
	!external functions
	real randu  !function for random number generation
	integer iseed !seed initialisation function
    character(1000) :: datafolder !the directory the data is stored in
    
    real, allocatable :: evals2d(:,:)
    real pdisp
    integer pstart
    
    integer num_args !number of arguments passed to the program (overwrites input file)

    real dummy

	!---------------------------- READ IN INFORMATION AND PERFORM INITIAL SET UP ------------------------------- -check:arg_temp_created

	write(*,*) 'Starting SIOPS Software'

	istat=99 !file number for problem output
	open(istat,file='potential_problems.txt')
	
	num_args = command_argument_count()
    
    !check whether the 4 key input files exist. If they don't, give the use the choice to generate examples
	call inputs_exist()
	
    !get evolutionary algorithm input
	call evol_input(nrep_MC,maxEA,runmode,EAseed,soilseeds,costmetric,deterministic,max_cons,popsize,perc_error,mut_rate,frac_keep,finaloutput,mut_lims,imut,datafolder,siea_mode,num_elite,stopmode,unistart,opt_local)

	!get soil input variables
	call readin_soil(nrep_MC)
		
    !get pile and building inputs
    call readin_pile(femrad,femdepth,prad,npdepths,pdepths,cg_tol,cg_limit,preps,plocation,nels,costvals,failurevals, & !nbh,bhdepth
		buildingweight,load_con,rel_loads,num_loads,dz,detdisp,femvech,femvecv,pilecost,difftol,usepie,regmesh,abstol,emean,datafolder,goodpiles,goodcases)
    
    !get site investigation inputs
    call readin_si(soffset,swidth,sstep,ntest,testerrors,sampfreq, & !nbh,bhdepth
		testcost,testnames,use_CI,conf_int,add_errors,in_tests,in_depths,in_reductions,bhnums,read_si,in_sdev,in_percentile,save_effE)     
        
    !do deterministic analysis to check average pile designs, if requested
    if (deterministic == 0) then
        call det_check(datafolder,pdepths,npdepths,num_loads,abstol,difftol,buildingweight,preps,rel_loads,load_con,in_sdev,in_percentile,prad,femrad,goodpiles,goodcases,plocation)
        stop
    end if
	
	!calculate the size of the correlation arrays and allocate them. Build the autocorrelation matrices for virtual soil generation.
	call soilsizes(istat)
    
    

	!if runmode == 'al' (all), then figure out which components need to be done based on the 
	!existence and size of files in the data folder, compared to the inputs
	call checkdata(runmode,datafolder,npdepths,preps(1)*preps(2),nrep_MC,thismode,startstress,prad,soilweight)
	
	!Export a range of virtual soils as CSVs if soil_reps(1) is greater than 1.
	!X is varying first, then Y, then Z the slowest.

	call exportsoil(soilseeds,nrep_MC,datafolder)

	
		!------------------------------ GET FEA PILE INFORMATION  --------------------------------
    
    if((runmode == 'PP' .or. thismode(1)) .and. singletrue) then
        
        !allocate(evals2d(nxew,nzew))
        !evals2d = 1
        !pstart = size(femvech)/2+1
        !do i=1,npdepths
        !    call prep_fem2d(evals2d,0.3,femvech,femvecv,prad(1),pdisp,pdepths(i),dz)
        !    write(*,*) pdisp
        !end do
        
        
        
        write(*,*) 'Generating pile settlement curve, and'
        write(*,*) 'associated soil weights for PIE method.'
        !The CK and SI modes both depend on the determininistic pile settlement curve generated with this subroutine.
        !Additionally, the CK settlement curves are calculated using the PIE method, and need the weights from this
        !subroutine in order to produce the weighted average soil at each pile length and the associated settlement.
        call getweights(dz,femrad,femdepth,prad,npdepths,pdepths,cg_tol,cg_limit,startstress,datafolder)
        if(runmode == 'PP') stop
    end if
    
    
    	!------------------------------ TRUE PILE PERFORMANCE  --------------------------------
	
	
	if((runmode == 'CK' .or. thismode(2))) then ! .and. singletrue) then	
		write(*,*) "Doing Complete Knowledge soil analysis."
		!Get CK pile performance by getting settlement for each pile and pile length
		
		
		call get_ck_set(soilseeds & !soil generation variables
                      ,preps,npdepths,detdisp,plocation,pdepths,cg_tol,cg_limit,prad,femrad,femdepth,pieval,soilweight,nels,ck_pset,nrep_mc,femvech,femvecv,usepie,regmesh,startstress,datafolder)

		if(runmode == 'CK') stop
    end if
    
    
    	!------------------------------ SITE INVESTIGATION PERFORMANCE --------------------------------
	!		This entire section is a part of the 'objective function' of the evolutionary algorithm
    
    
    if(runmode == 'SI' .or. thismode(3) .or. (.not. singletrue)) then
    	
		
		!get complete knowledge pile performance
        if(singletrue .or. .not. (singletrue .or. usepie))  call read_ck(nrep_MC,npdepths,preps,ck_pset,datafolder,detdisp,pdepths,prad)
        if((.not. singletrue)) call prepmultiSI(npl2,goodpiles,goodcases,sdist,sumweights,extents,indices,CKheight,preps,load_con,usepie,femrad,soilseeds,prad,nrep_MC,plocation)
        
        !Nodisksave: 0 = keep in ram (much slower, but can save a fair amount of disk space), 1=write and read soil, 2=read soil only
        !if(deterministic) nodisksave = 0 !hardcode it to be zero in the deterministic case for now, as this is safer if you're not paying attention to inputs.
        
        !This subroutine writes the soils to disk so that they may be loaded in at a later time.
		!This greatly increases the speed of the process in the Evolutionary algorithm, where the 
		!same soils are used over and over
        !NOTE: This option has been made redundant by the superset option which is much faster, and can be used in parallel processes.

		!call storesoils(soilseeds,nrep_MC,datafolder)

            !Get a very big soil from which to extract a subset (for superset mode)
        !Note that superset won't work properly with anisotropy or with FEM specified for CK analysis.
	    if (singletrue .and. superset) then
	    	write(*,*) 'Generating soil superset'
            kseed = randu(soilseeds(1)) * 1234567890 !ensure random numbers are consistent across MC realisations
            call RF3D(soil_dummy_var)
            !open(3579, file='supersoil.dat', access="stream")
            !!!call outputstats(nxe,nye,nze,dz,dz,dz,efld(:nxe,:nye,:nze),15.,1.)
            !!
            !write(3579) efld(:nxe,:nye,:nze)
            !!
            !close(3579)
            !stop
        end if
        
        
        !save soil description to a string for use in input and output 
        if (singletrue) then 
			write(soildsc,'(A,I0,A,I0,A,I0)') '_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy
		else
			write(soildsc,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '_nlayers-',nlayer,'_ratio2l-',nint(lmean_ave(1)/minval(lmean_ave(:2))),'+',nint(lmean_ave(2)/minval(lmean_ave(:2))),'_depth2l-',nint(dz*ldepths(1)),'_bSD-',nint(bsd*dz)
		end if
        
        
        
        
        
        !deterministic analysis: Only do the Monte Carlo simulation once, on a fixed set of site investigation inputs.
    	!Note, this does not save the soils to disk. It generates them on the fly. 
    	!Recommend 2000 MC realisations when calculating probability of failure, and 8000 MC realisations when calculating failure cost.
        if(deterministic == 1) then
            
            if (read_si > 1) then
				!read in the full set of data for individual investigations to analyse. 
				!This input can be typed out manually, but is usually the output/product of the evolutionary subroutine for the purpose of a deterministic costing analysis
				call input_SI(nbhmax,sampfreq,inv_coords,inv_depths,inv_nbh,inv_test,ninv,inv_reduction,percentile,s_dev,si_performance,dz,testnames,deterministic,in_sdev,in_percentile,soildsc,goodcases, read_si)
			else
				!Generate a set of investigations to analyse: all combinations of the specified parameters.
				call deterministic_SI(soffset,swidth,nbhmax,in_tests,in_depths,in_reductions,bhnums,testerrors,sampfreq,in_sdev,in_percentile, &		!input variables
					inv_coords,inv_depths,inv_nbh,inv_test,ninv,inv_reduction,percentile,s_dev,si_performance)	 !output
            end if
            
            write(*,'(A,X,I0,X,A)') "Performing fixed analysis on",ninv,"fixed investigation options."

			allocate(invcount(ninv),fcost(ninv), pcost(ninv), probfail(ninv), avediff(ninv), diffgeo(ninv))
            
            if (singletrue) then

                !Get site investigation performance in variable, single layer soils
    	        call get_si_perf(soilseeds & !soil generation variables
                          ,ninv,nbhmax,size(in_tests),size(in_depths),size(in_reductions),inv_nbh,inv_coords,inv_depths,inv_test,add_errors,use_CI,testerrors,inv_reduction,conf_int,percentile,s_dev,soffset,swidth,sstep &					!site investigation variables
                          ,npdepths,detdisp,plocation,pdepths,ck_pset,nrep_mc,rel_loads,load_con,num_loads,preps,buildingweight,difftol,failurevals,costvals,si_performance,pilecost,testcost,testnames,costmetric,1,deterministic,finaloutput,abstol,datafolder,invcount,goodpiles,goodcases,fcost, pcost, probfail, avediff, diffgeo)
            
            else

                !Get site investigation performance in multiple layer soils with uniform soil properties
                call get_si_perf_multi(soilseeds & !soil generation variables
                          ,ninv,nbhmax,size(in_tests),size(in_depths),size(in_reductions),inv_nbh,inv_coords,inv_depths,inv_test,add_errors,use_CI,testerrors,inv_reduction,conf_int,percentile,s_dev,soffset,swidth,sstep &					!site investigation variables
                          ,npdepths,detdisp,plocation,pdepths,ck_pset,nrep_mc,rel_loads,load_con,num_loads,preps,buildingweight,difftol,failurevals,costvals,si_performance,pilecost,testcost,testnames,costmetric,finaloutput,deterministic,abstol, &
                    femvech,femvecv,prad,datafolder,usepie,npl2,goodpiles,goodcases,sdist,sumweights,extents,indices,CKheight,1,invcount,fcost, pcost, probfail, avediff, diffgeo)  !multiple layer site investigations
            end if
            
            
        else if(deterministic == 2) then
            
            write(*,'(A,X,I0,X,A)') "Performing heat map analysis on",size(in_tests) * size(in_reductions) * size(in_depths),"investigation options."

            
            do i = 1,size(in_tests) * size(in_reductions) * size(in_depths)
                
                
                
                 call full_spread(soffset,swidth,sstep,in_tests,in_depths,in_reductions,testerrors,sampfreq,in_sdev,in_percentile, &		!input variables
		                inv_coords,inv_depths,inv_nbh,inv_test,ninv,inv_reduction,percentile,s_dev,si_performance,i,invcount,fcost, pcost, probfail, avediff, diffgeo)	
                 
                 write(*,'(A,I0,A,F5.1,A)') 'Starting investigation ',i,': '//testnames(inv_test(1,i))//' to a depth of ',dz*inv_depths(1,i,2),' m, with the '//rednames(inv_reduction(i))//' method.'
            
                 
     !            deterministic_SI(soffset,swidth,nbhmax,in_tests,in_depths,in_reductions,bhnums,testerrors,sampfreq,in_sdev,in_percentile, &		!input variables
					!inv_coords,inv_depths,inv_nbh,inv_test,ninv,inv_reduction,percentile,s_dev,si_performance)	 !output
            
                if (singletrue) then
            
                    !Get site investigation performance in variable, single layer soils
    	            call get_si_perf(soilseeds & !soil generation variables
                              ,ninv,1,size(in_tests),size(in_depths),size(in_reductions),inv_nbh,inv_coords,inv_depths,inv_test,add_errors,use_CI,testerrors,inv_reduction,conf_int,percentile,s_dev,soffset,swidth,sstep &					!site investigation variables
                              ,npdepths,detdisp,plocation,pdepths,ck_pset,nrep_mc,rel_loads,load_con,num_loads,preps,buildingweight,difftol,failurevals,costvals,si_performance,pilecost,testcost,testnames,costmetric,i,deterministic,finaloutput,abstol,datafolder,invcount,goodpiles,goodcases,fcost, pcost, probfail, avediff, diffgeo)
            
                else
            	
            
                    !Get site investigation performance in multiple layer soils with uniform soil properties
                    call get_si_perf_multi(soilseeds & !soil generation variables
                              ,ninv,1,size(in_tests),size(in_depths),size(in_reductions),inv_nbh,inv_coords,inv_depths,inv_test,add_errors,use_CI,testerrors,inv_reduction,conf_int,percentile,s_dev,soffset,swidth,sstep &					!site investigation variables
                              ,npdepths,detdisp,plocation,pdepths,ck_pset,nrep_mc,rel_loads,load_con,num_loads,preps,buildingweight,difftol,failurevals,costvals,si_performance,pilecost,testcost,testnames,costmetric,i,deterministic,abstol, &
                        femvech,femvecv,prad,datafolder,usepie,npl2,goodpiles,goodcases,sdist,sumweights,extents,indices,CKheight,finaloutput,invcount,fcost, pcost, probfail, avediff, diffgeo)  !multiple layer site investigations
                end if
                
            end do
            
            
        
        else if(deterministic == 3) then
        	
			
			!Set up inputs for each investigation.
			call EA_SI(soffset,swidth,nbhmax,in_tests,in_depths,in_reductions,bhnums,testerrors,sampfreq,in_sdev,in_percentile, &		!input variables
				inv_coords,inv_depths,inv_nbh,inv_test,ninv,inv_reduction,percentile,s_dev,si_performance)	 !output variables
	
            write(*,'(A,X,I0,X,A)') "Performing Evolutionary Algorithm to optimise",ninv,"site investigations (borehole locations)."
            
            	! ---- output prep ---
			if(costmetric == -4) str3 = 'Expected_failure_cost'
			if(costmetric == -3) str3 = 'Expected_total_cost'
            if(costmetric == -2) str3 = 'Average_differential_settlement'
            if(costmetric == -1) str3 = 'Prob._Failure'
            if(costmetric >= 0 ) write(str3,'(I0,A)') costmetric,'_geometric_SD'

	
			write(str2,'(A,I0,A,I0,A,I0,A,A)') 'si_results/FinalEA-stats_piles-',goodcases,'_BArea-',nint(buildingarea),'_BH-',inv_nbh(1),trim(soildsc),'.txt'
	
			open(500,file=str2)!,buffered='no')
			write(500,'(I0,X,A)') popsize,	'!total number of investigations. Using metric: '//str3
			write(500,'(I0,X,A)') nbhmax, '!max no. boreholes'
			write(500,'(I0,X,A)') size(in_tests), '!No. Tests'
			write(500,'(I0,X,A)') size(in_depths), '!No. BH depths'
			write(500,'(I0,X,A)') size(in_reductions), '!No reduction methods'
			write(500,'(A6,2X,A4,2X,A4,2X,A9,2X,A10,2X,A8,X,A8,X,A15,X,A15,X,A15,X,A15,X,A15)') 'No._BH','Test','Red.','Depth_(m)','%_Bad_Reps','No._Gens','Time_(s)','Fail._cost_($)','Pile_cost_($)','Prob._Fail._(%)','Ave._Diff._Set.','Geo._statistic'
	
			!output all borehole coordinates in metres from the corner of the virtual soil
			!Each column corresponds to a borehole. Each row corresponds to an investigation/population member.
			write(str2,'(A,I0,A,I0,A,I0,A,A)') 'si_results/FinalEA-Xcoords_piles-',goodcases,'_BArea-',nint(buildingarea),'_BH-',inv_nbh(1),trim(soildsc),'.txt'
			open(501,file=str2) !,buffered='no')
			write(str2,'(A,I0,A,I0,A,I0,A,A)') 'si_results/FinalEA-Ycoords_piles-',goodcases,'_BArea-',nint(buildingarea),'_BH-',inv_nbh(1),trim(soildsc),'.txt'
			open(502,file=str2) !,buffered='no')

			write(501,'(1000000(I8))') (j,j=1,nbhmax)
			write(502,'(1000000(I8))') (j,j=1,nbhmax)
      
            !----call genetic algorithm----
            do i = 1,ninv !Loop through each investigation set to optimise them individually
            	write(*,'(A,I0,A,I0,A,F5.1,A)') 'Starting investigation ',i,': ',inv_nbh(i),' boreholes with a '//testnames(inv_test(1,i))//' to a depth of ',dz*inv_depths(1,i,2),' m, with the '//rednames(inv_reduction(i))//' method.'
            
                !Prepare current investigation for the Evolutionary Algorithm here. 
                call EAprep(soilseeds, & !soil generation variables
                      inv_nbh(i),inv_depths(1,i,:),inv_test(1,i),add_errors,use_CI,testerrors,inv_reduction(i),conf_int,percentile(i),s_dev(i),soffset,swidth,sstep &					!site investigation variables
                      ,npdepths,detdisp,plocation,pdepths,ck_pset,nrep_mc,rel_loads,load_con,num_loads,preps,buildingweight,difftol,failurevals,costvals,si_performance(i),pilecost,testcost,testnames,costmetric,abstol &
                    ,maxEA,deterministic,max_cons,popsize,perc_error,mut_rate,frac_keep,finaloutput,mut_lims,imut,datafolder,i,EAseed,siea_mode,num_elite,stopmode,unistart,opt_local, & !EA variables
                    femvech,femvecv,prad,usepie,npl2,goodpiles,goodcases,sdist,sumweights,extents,indices,CKheight) !multi layer site investigation variables
            
            end do
            
            close(500)
            close(501)
            close(502)
        
        end if
        
    end if          
	

	close(istat)

end program
