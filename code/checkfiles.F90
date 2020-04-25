module checkfiles
    
    use variables
    
    implicit none
    
    contains
    
    !check to see that the neccessary data has been previously generated. If it hasn't, then automatically generate what's missing.
    !This subroutine also reads in the soil weights

	subroutine checkdata(runmode,datafolder,npdepths,npl,nrep_MC,thismode,startstress,prad,soilweight)
	
	!input variables
	integer, intent(in) :: npdepths,npl,nrep_MC !number of pile depths, piles and monte carlo realisations
	character(1000), intent(in) :: datafolder !location of data directory
	character(2), intent(in) :: runmode
    real, allocatable, intent(out) :: soilweight(:,:,:,:)
    integer, intent(in) :: prad(2)
    integer, allocatable, intent(out) :: startstress(:,:)
	
	!output
	logical, intent(out) :: thismode(3) !whether to do the PP, CK or SI analysis respectively
	
	!local variables
	logical exists !check if current file exists
	integer readstatus !check read status
	integer sumfiles !tracker for number of files
	real placeholder,p2 !placeholder variables
	character(1000) str2
	integer i,j !loop counter
	
	
	thismode = .false. !set default
    
        

	
	if(singletrue) then !do the file checks to see what processess need to be done
	
	
	!------Check deterministic settlement curve and soil weights for each pile depth (SI mode)------
	
		!look at deterministic settlement curve
        
        if(runmode == 'AL' .or. runmode == 'CK') then
		
            write(str2,'(A,A,I0,A,F4.2,A)') trim(datafolder),'settlement_prad-',prad(1),'_esize-',dz,'.txt'
		    inquire(file=str2,exist=exists) !check if the file exists 
		    if(.not. exists) then
			    if (runmode == 'AL') then
                    write(*,*) "Warning, deterministic pile settlement curve doesn't exist"
                    write(*,*) "Generating new curve."
			        thismode(1) = .true.
                else
                    write(*,*) "Error, deterministic pile settlement curve doesn't exist"
                    write(*,*) "Please run program with 'AL' or 'PP' stage setting."
                    read(*,*)
                    stop
                end if
			
		    else
				
			    open(667, file= str2,status='old')
			    read(667,*) 
			    do i = 1,npdepths
				    read(667,*,iostat=readstatus) p2,placeholder
				    if(readstatus < 0) then	
                        if  (runmode == 'AL') then
					        write(*,*) "Warning, number of pile depths doesn't match current input."
					        write(*,*) "Generating deterministic pile settlement curve"
					        thismode(1) = .true.
                        else
                            write(*,*) "Error, number of pile depths doesn't match current input."
                            write(*,*) "Please run program with 'AL' or 'PP' stage setting."
                            read(*,*)
                            stop
                        end if
					    exit
				    end if
			    end do
			    close(667)
			
		    end if
		

		    !look at soil weights for each pile length
	
		    !Check to see how many files exist

            write(str2,'(A,A,I0,A,F4.2,A)') trim(datafolder),'soilweights_prad-',prad(1),'_esize-',dz,'.dat'
		    inquire(file=str2,exist=exists) !check if the file exists 

	
		    !If all the files aren't here, generate a new set
		    if(.not. exists) then
	            if(runmode == 'AL') then
                    write(*,*) "Warning, PIE soil weights haven't been generated for"
                    write(*,*) "current pile radius and element size"
			        write(*,*) "Generating new set of soil weights"
			        thismode(1) = .true.
                else
                    write(*,*) "Error, PIE soil weights haven't been generated for"
                    write(*,*) "current pile radius and element size"
                    write(*,*) "Please run program with 'AL' or 'PP' stage setting."
                    read(*,*)
                    stop
                end if
			
            else !if all the files are here, check to see that the files themselves are complete
            
            
                allocate(startstress(3,npdepths))
                write(str2,'(A,A,I0,A,F4.2,A)') trim(datafolder),'model_bounds_prad-',prad(1),'_esize-',dz,'.txt'
	            open(667,file=str2,status='old')
                read(667,*) 
                do j = 1,npdepths
	                read(667,*) startstress(:,j)
                end do
                close(667)
    
                allocate(soilweight(npdepths,maxval(startstress(1,:))*2+prad(1),maxval(startstress(2,:))*2+prad(2),maxval(startstress(3,:))))
	
                write(str2,'(A,A,I0,A,F4.2,A)') trim(datafolder),'soilweights_prad-',prad(1),'_esize-',dz,'.dat'
			    open(666,file=str2,access='stream')
			    read(666,iostat=readstatus) soilweight								!read in data
			    !readstatus < 0 means end of file has been reached, in which case the current inputs aren't compatible with previous input
			    !regenerate data
			    if(readstatus < 0) then	
                    if(runmode == 'AL') then
				        write(*,*) "Warning, PIE soil weight field size doesn't match current input."
				        write(*,*) "Generating new set of pile settlement curves"
				        thismode(1) = .true.
                    else
                        write(*,*) "Error, PIE soil weight field size doesn't match current input."
                        write(*,*) "Please run program with 'AL' or 'PP' stage setting."
                        read(*,*)
                        stop
                    end if
			    end if
			    close(666)
		    end if
        end if
        
        if(runmode == 'AL' .or. runmode == 'SI') then
	
    
	
	    !------Check pile settlement files for current soil input (CK mode)------
	
		    !Check to see how many files exist
		    sumfiles = 0
		    do i = 1,npl
                write(str2,'(A,A,I0,A,I0,A,F4.2,A,I0,A,I0,A,I0,A)') trim(datafolder),'ck_pile-',i,'_prad-',prad(1),'_esize-',dz,'_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy,'.txt'
			    !write(str2,'(A,A,I0,A,I0,A,I0,A,I0,A)') trim(datafolder),'ck_pile-',i,'_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy,'.txt'
			    inquire(file=str2,exist=exists) !check if the file exists 
			    if(exists) sumfiles = sumfiles + 1
		    end do
	
		    !If all the files aren't here, generate a new set
		    if(sumfiles < npl) then
	
			    if(sumfiles == 0) then
				    write(*,*) "Warning, true pile settlement functions don't exist for current inputs"
			    else
				    write(*,*) "Warning, not all piles have associated settlement data."
                end if
		
                if(runmode == 'AL') then
			        write(*,*) "Generating new set of pile settlement curves"
			        thismode(2) = .true.
                else
                    write(*,*) "Run program with 'CK' or 'PP' stage setting."
                    read(*,*)
                    stop
                end if
		
            else !if all the files are here, check to see that the files themselves are complete
                write(str2,'(A,A,I0,A,I0,A,F4.2,A,I0,A,I0,A,I0,A)') trim(datafolder),'ck_pile-',1,'_prad-',prad(1),'_esize-',dz,'_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy,'.txt'
			    !write(str2,'(A,A,I0,A,I0,A,I0,A,I0,A)') trim(datafolder),'ck_pile-',1,'_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy,'.txt'
			    open(666,file=str2) 
			    do i = 1,nrep_MC
				    read(666,'(100000(E14.7,X))',iostat=readstatus) (placeholder,j=1,npdepths)
				    !readstatus < 0 means end of file has been reached, in which case the current inputs aren't compatible with previous input
				    !regenerate data
				    if(readstatus < 0) then	
                        if(runmode == 'AL') then
					        write(*,*) "Warning, not all pile depths or MC realisations have associated settlement."
					        write(*,*) "Generating new set of pile settlement curves"
					        thismode(2) = .true.
                        else
                            write(*,*) "Error, not all pile depths or MC realisations have associated settlement."
                            write(*,*) "Please run program with 'AL' or 'CK' stage setting."
                            read(*,*)
                            stop
                        end if
                        exit
				    end if
			    end do
			    close(666)
	
            end if
            
        end if
    
	
    end if
    
    if (runmode == 'AL') then
        
        	!------ Do site investigation component --------
	
		thismode(3) = .true. 
        
    end if
	
	
    end subroutine
    
    
    !Creates a dummy set of inputs in the main directory
    subroutine save_templates()
    
    integer unit
    
    unit=1951
    open(unit,file='si_input.txt')
    
write(unit,*) "!Site investigation information"
write(unit,*) "41,41		!initial borehole x,y offset from corner of soil in elements (ideally coincides with building plan). This is ignored with the EA."
write(unit,*) "80,80		!dimensions of site investigation area in elements. This is ignored with the EA."
write(unit,*) "2,2		!borehole step size in each dimension for the heat map mode"
write(unit,*) "6		!No. test types (leave at 6)"
write(unit,*) ".true.		!use confidence interval to truncate unrealistic sample values (recommended, especially if used with 'add_errors' option."
write(unit,*) "2.576		!confidence interval z score for the above truncation, recommend 99% (2.576)"
write(unit,*) ".true.		!'add_errors' - add random errors to tested samples (highly recommended for realistic sampling)"
write(unit,*) "!Below are test stats: Transformation error, bias error, measurement error, sampling frequency (elements), cost per metre."
write(unit,*) "0 , 0 , 0,         3,  120, dct			!discrete testing, no errors"
write(unit,*) "0 , 0 , 0,         1,  77,  cts			!continuous testing, no errors"
write(unit,*) ".20 , .40 , .25,   3,  156, SPT			!SPT"
write(unit,*) ".15 , .20 , .15,   1,  77,  CPT			!CPT"
write(unit,*) ".20, .20, 10e-30,  3,  330, TT			!Triaxial test"
write(unit,*) ".15 , .15, .10,    3,  120, DMT			!DMT"
write(unit,*) "1,0.25		!standard deviations below the mean in the Geometric 'SD' reduction method. Percentile to use for the percentile reduction method. Recommend 1 and 0.25 currently."
write(unit,*) "!note: the following information is only relevant is read_si is set to false, in which case the investigations to assess are built using all combinations of the below parameters."
write(unit,*) "0		!read_si parameter. 0 = Generate investigations from the below variables. 1 = read investigation configurations from 'input/si.txt' data. 2 = also read tests and depths on a per-borehole basis."
write(unit,*) "30		!Number of borehole cases to investigate	Note: These 4 lines correspond to the 4 lines below them."
write(unit,*) "4		!Number of tests to investigate"
write(unit,*) "5		!Number of reduction methods to investigate"
write(unit,*) "1		!Number of borehole depth cases to analyse"
write(unit,*) "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,24,25,30,36,40,45,49,55,64,81,100"
write(unit,*) "3,4,5,6	!Test numbers corresponding to the above tests"
write(unit,*) "SA,GA,HA,1Q,SD	!Reduction methods"
write(unit,*) "80	!Borehole depths (elements)"

    
    close(unit)
    
    
    open(unit,file='soil_input.txt')
    
write(unit,*) ".false.						!.true. for variable, single layer analysis, .false. for multiple layer analysis using uniform soil properties within each layer"
write(unit,*) 
write(unit,*) "192 192 96	  	  			!number of x,y,z elements in original generated soil (must fit a*2^b)"
write(unit,*) "80  80  40					!x,y,z dimensions of site (m)"
write(unit,*) ".true.						!store soils in memory. Should use one very large soil in single layer mode."
write(unit,*) "5						!The upscale factor for single layer soils. Recommend 5."
write(unit,*) "dlavx3	dlavx2					!correlation function for 3D and 2D LAS"
write(unit,*) "0.5						!size of elements (cube length)"
write(unit,*) "2048,	6					!limit for stage 0 matrix size, where this value is k1*k2*k3, max no. subdivisions"
write(unit,*) "-1 1						!Lower and upper range of MC realisations to export soils for. 1st value -ve doesn't export. 2nd value -ve exits after exporting."
write(unit,*) "!Single layer soil options"
write(unit,*) "l						!soil distribution; n = normal, l = lognormal, b = beta"
write(unit,*) "15.0  0.8   					!Soil stiffness: Mean and COV (as a proportion of the mean). (Single layer only)"
write(unit,*) "30.						!horizontal SOF (single layer only)"
write(unit,*) "!Multiple layer analysis options"
write(unit,*) "100						!boundary SOF"
write(unit,*) "8						!boundary SD (elements)"
write(unit,*) "2						!number of layers"
write(unit,*) "20,40,30,40,50,60,70				!Layer depth (elements)	(only needed for the n-1'th top layers)"
write(unit,*) "5,50,40,4,1,10					!Young's modulus of each layer"
write(unit,*) "1, 10						!standard deviation of Young's modulus for each layer"
write(unit,*) ".false.						!ignore depths and properties of each layer, and instead read in explicit borehole depths from other file"
write(unit,*) ".true.						!enforce fixed layer depths at boreholes despite added randomness, if the above option is set to true"
write(unit,*) ""
write(unit,*) ""


    close(unit)
    
    
    open(unit,file='EA_input.txt')
    
write(unit,*) "1		!Program run mode 0 = pile design test, 1 = fixed mode, 2 = heat map mode, 3 = evolutionary mode"
write(unit,*) "8000		!number of Monte Carlo realisations"
write(unit,*) "-4		!Performance metric to use. -4 failure cost, -3 total cost, -2 ave. diff. Settlement, -1 prob. Failure, >=0 geometric statistic"
write(unit,*) "100,100		!random seeds for virtual soil and genetic algorithm. Set to 0 to be random (clock based) or higher for consistent random numbers. VERY strongly recommend soilseed be kept as a positive integer (e.g 100)."
write(unit,*) "AL  		!single layer processing option. PP = stage 1 preprocessing. CK = stage 2 preprocessing. SI = site investigation analysis. AL = automatic."
write(unit,*) "4		!output mode for heat maps and GA. 1 outputs minimal, essential information. 4 outputs all information. See manual for other values."
write(unit,*) "'.\data\'	!directory of the data folder"
write(unit,*) "!evolutionary variables"
write(unit,*) "200		!max number of evolutionary algorithm iterations"
write(unit,*) "5		!max number of consecutively equal values"
write(unit,*) ".false.		!Stopping mode: .true. means wait for it to have converged (remain unchanged for a number of realisations), .false. means wait for a number of iterations past the current GLOBAL optimum before stopping (in case fitness starts increasing)."
write(unit,*) "1000		!Population size"
write(unit,*) "1.0		!Percentage error for stop criteria (out of 100)"
write(unit,*) "0.1,0.001,0.2	!Mutation rates: Initial, minimum, maximum"
write(unit,*) "1		!mutation mode; 1 = constant, 2 = adapt based on fitness, 3 = adapt based on proximity in parameter space"
write(unit,*) "0.5		!Fraction of population to keep as parents"
write(unit,*) "1		!Number of population members to keep UNCHANGED from previous population (this concept is called elitism)."
write(unit,*) "1		!manner of controlling the borehole locations in the EA. 1 Means start across full soil. 2 Means Start within building footprint. 3 Means ALWAYS keep within building footprint."
write(unit,*) ".false.		!generate initial population using a uniform distribution, otherwise use normally-distributed offsets from a regular grid, iff possible"
write(unit,*) ".true.		!conduct a 2nd phase of genetic algorithm optimisation using normal offsets from the fittest investigation of phase 1"

    
    close(unit)
    
    
    open(unit,file='pile_input.txt')
    
 write(unit,*) "!Pile information"
write(unit,*) "1, 1		!width of the pile in each dimension (elements)"
write(unit,*) "0.002		!Differential settlement design tolerance (mm/mm). Becomes an absolute settlement design tolerance as a multiple of pile spacing."
write(unit,*) "-1		!If positive, this value becomes the absolute design tolerance (mm)."
write(unit,*) "20		!number of pile depths to assess"
write(unit,*) "0.0001, 5000	!convergence tolerance and max number of iterations for FEM subroutine"
write(unit,*) "1		!how to treat layer boundaries at pile: (1) a single point in SI and CK, (2) inverse weighted in CK, (3) inverse weighted in SI and CK"
write(unit,*) ".false.		!If true give the coordinates directly in a row of x coords and row of y coords, and each column is a pile. Otherwise, spread evenly over a grid."
write(unit,*) "5		!Specify number of piles when giving direct coordinates"
write(unit,*) "1		!List of pile x coordinates (these lists are ignored if the above value is false)"
write(unit,*) "1		!List of pile y coordinates"
write(unit,*) "1		!Index of relative load for each pile"
write(unit,*) "5,5		!Specify a grid of X * Y piles (e.g. 5*5 = 25 total)"
write(unit,*) "41,41		!initial pile x,y offset from corner of soil in elements (must be >= soil radius)"
write(unit,*) "80,80		!foundation plan dimensions in the x,y directions (elements)"
write(unit,*) "1 1 1 1 1	!index of relative load for each pile"
write(unit,*) "1 1 1 1 1	!i.e with current settings, 1=0.25, 2=0.5, 3=1.0 based on tributary area. These values correspond to corner, edge and internal piles respectively assuming a rectangular building that is symmetrical in the x and y directions."
write(unit,*) "1 1 1 1 1"
write(unit,*) "1 1 1 1 1"
write(unit,*) "1 1 1 1 1"
write(unit,*) "0.25 0.5 1.0	!relative pile loads for each pile load cases"
write(unit,*) "1600		!Plan view area of a single floor in the building (i.e. 1600 = 40 m by 40 m, 400 = 20m by 20m. Should be less than or equal to the foundation area ideally). This is used in the building weight and structural cost calculation."
write(unit,*) "3		!Number of floors"
write(unit,*) "8		!uniformly distributed load over 1 square metre. Recommend 8 kPa (5 dead load + 3 live load). Used in building weight calculation. Weight = unit applied load * plan area * No. Floors"
write(unit,*) "0.003,0.009	!lower and upper bound for linear differential settlement -> failure cost function"
write(unit,*) "47500000	!building cost"
write(unit,*) "200		!pile cost per metre"


    
    close(unit)
    
    
    open(unit,file='readme.txt')
    
    write(unit,*) "This program requires 4 files in an 'input' folder in the program's directory:"
    write(unit,*) "EA_input.txt   - overall control and genetic algorithm parameters"
    write(unit,*) "pile_input.txt - variabiles related to the structure and foundation"
    write(unit,*) "soil_input.txt - control over the soil properties and choice of single layer variable vs multi-layer constant"
    write(unit,*) "si_input.txt   - site investigation parameters as well as test inaccuracy data"
    write(unit,*)
    write(unit,*) "Furthermore, there should be a 'si_results' folder for results to be saved in"
    write(unit,*) "and a 'data' folder, the location of which must be given in the EA_input.txt file"
    write(unit,*) 
    write(unit,*) "Note: There is a SIOPS manual avaiable for download."
    
    close(unit)
    
    end subroutine

end module
