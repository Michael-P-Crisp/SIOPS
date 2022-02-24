!c  *********************************************************************
!c  *                                                                   *
!c  *                         subroutine reduce readinm                 *
!c  *                                                                   *
!c  *********************************************************************
!c  Single Precision Version 1.0
!c  Written by Michael P Crisp
!c  Latest Update: Jun 2, 2016
!c
!c  PURPOSE Reads in various inputs regarding investigation and design
!c
!c-------------------------------------------------------------------------

module readinm

use variables
use checkfiles

implicit none

!note that the mac/linux "/" will need to be the "\" on windows
character(1) :: slash = "/"
logical exists
integer num_args !number of arguments passed to the program (overwrites input file)
integer :: cur_arg=0 !current argument
character(100) argval !storage for argument
integer numpiles

private :: num_args,cur_arg,argval,numpiles

contains


	!this function will return the argument
	real function getarg_real() 
		cur_arg = cur_arg + 1
		call get_command_argument(cur_arg,argval)
		read(argval,*) getarg_real
	end function
	
	
	!whether the 4 key input files exist. If they don't, give the use the choice to generate examples
	subroutine inputs_exist()
    
        logical exists(4) !whether the 4 key input files exist
        integer choice,istat
    
        !note that the "/" will need to be the "\" on windows
        inquire(file='input'//slash//'EA_input.txt',exist=exists(1)) !check if the file exists
        inquire(file='input'//slash//'soil_input.txt',exist=exists(2)) !check if the file exists
        inquire(file='input'//slash//'si_input.txt',exist=exists(3)) !check if the file exists
        inquire(file='input'//slash//'pile_input.txt',exist=exists(4)) !check if the file exists 
		if(.not. all(exists)) then
            write(*,*) "Error, at least one of the following key input files" 
            write(*,*) "do not exist in the 'input' folder:"
            write(*,*) '    EA_input.txt'
            write(*,*) '    soil_input.txt'
            write(*,*) '    si_input.txt'
            write(*,*) '    pile_input.txt'
            write(*,*)
            write(*,*) "Would you like to generate example files in the program directory?"
            write(*,*) "If so press '1' then 'enter'."
            write(*,*) "Otherwise, press 'enter' key to exit."
            read(*,'(I1)',iostat=istat) choice
            !if a number was entered, and the number is '1', then save the output.
            if(istat == 0 .and. choice == 1) call save_templates()
            stop
        end if
    
    
    
    
    end subroutine
	
	


	!read in LAS soil generation data

	subroutine readin_soil(nrep_MC)



		integer  inn !soil property input  
        integer, intent(in) :: nrep_MC
        
        real, allocatable :: lsd(:) !standard deviation of E for each layer

		!integer, intent(out) :: ndpvec(:)
		integer :: xdim,ydim,zdim !dimensions of working area in m (to be converted into elements)
		integer j,i
		real pvr !temporary lognormal variable
        integer :: superscale !the multiplier from which to upsize the original soil in each dimension. E.g 5 = 5*5*5 = 125 times the volume.
        real :: stiff_mult !stiffness multiplier when using command line arguments

		integer istat
        
        

		inn= 500

		open(inn, file='input'//slash//'soil_input.txt',status='old')
		
		
		!check to see if any arguments were passed to the program
		num_args = command_argument_count()
		cur_arg = 0
		

        !True for variable, single layer analysis, false for uniform, multi-layer analysis
        read(inn,*) singletrue
        read(inn,*)
        
        !Note that superset won't work properly with anisotropy or with FEM specified for CK analysis.

		! read in simulation options
		read (inn,*) nxe, nye, nze
		read (inn,*) xdim, ydim, zdim
        read (inn,*) superset
        read (inn,*) superscale
		read (inn,*) varfnc,bvarfnc
		read (inn,*) dz
		read (inn,*) MXK,MXM
        read (inn,*) soil_reps(1),soil_reps(2)
        read (inn,*) 

		!read in single layer options

		read (inn,*) distribution
		if(distribution == 'n' .or. distribution == 'l') then
			read (inn,*) sdata(1,1),sdata(1,2)
			if(num_args >= 2) then !if an argument is present, then use that instead
				stiff_mult = getarg_real()
				sdata(1,1) = sdata(1,1) * stiff_mult
				sdata(1,2) = getarg_real()
			end if
			emean = sdata(1,1)
			sdata(1,1) = 1	!Scale the 
			sdata(1,2) = sdata(1,1) * sdata(1,2) !Get the standard deviation from the COV.
		else if(distribution == 'b') then
			read (inn,*) sdata(1,1),sdata(1,2),sdata(1,3),sdata(1,4)
		end if
		read (inn,*) soilth(1), soilth(2)
        nodisksave = -1 !deprecated; changing this value shouldn't affect anything anymore
        
        
        !Read in multiple layer options
        read(inn,*)
        read(inn,*) bsof
        read(inn,*) bsd
        read(inn,*) nlayer
        if(num_args >= 3) nlayer = nint(getarg_real())
        allocate(lmean_ave(nlayer),lsd(nlayer),ldepths(nlayer))
        read(inn,*) (ldepths(i),i=1,nlayer-1)
        read(inn,*) (lmean_ave(i),i=1,nlayer)
        read(inn,*) (lsd(i),i=1,nlayer)
        read(inn,*) rand_realisations
        if(num_args >= 2) then
        	lmean_ave = lmean_ave * stiff_mult
        end if
        read(inn,*) uselcoords
        read(inn,*) enforcelcoords
        
        !overwrite some more variables if arguments are present
        if(num_args >= 5) then
        	soilth(1) = getarg_real()
            soilth(2) = getarg_real() !soilth(1)
! 			nlayer = getarg_real()
! 			do i = 1,nlayer
! 				ldepths(i) = getarg_real()
! 			end do
! 			do i = 1,nlayer
! 				lmean_ave(i) = getarg_real()
! 			end do
        end if
        
        ! The soil is considered isotropic if the SOF in the horizontal and vertical directions is within a 1 mm tolerance (conservative tolerance)
        is_isotropic = abs(soilth(1) - soilth(2)) < 0.001
        anisotropy = nint(soilth(1)/soilth(2))	!for label purposes, round anisotropy to the nearest integer
        
        
        if(num_args >= 6) bsd = nint(getarg_real())
        
        if(num_args >= 7) ldepths(1) = nint(getarg_real())
        
        if(num_args >= 9) then
            deallocate(lmean_ave)
            allocate(lmean_ave(max(nlayer,3)))
            lmean_ave(1) = getarg_real()
            lmean_ave(2) = getarg_real()
            lmean_ave(3) = lmean_ave(1)
        end if
        

		!convert some input variables into element-compatible units:
		nxew=xdim/dz
		nyew=ydim/dz
		nzew=zdim/dz
        
        ! determine whether each layer's Young's properties are to be represented by a 2D random field for the multi-layer case
        var_props = any(lsd > 0)

		
		!calculate lognormal parameters for single layer soil
		if (distribution == 'l') then
			pvr = log(1.0 + sdata(1,2) ** 2 / sdata(1,1) ** 2)
			sdata(1,3) = log(sdata(1,1)) - 0.5 * pvr
			sdata(1,4) = sqrt(pvr)
        end if 
        
        if (singletrue) then
            if (superset > 1) then
                nxe = nxe * superscale
                nye = nye * superscale
                nze = nze * superscale
                !allocate full soil field
                if (is_isotropic) then
                    allocate(efld(nxe,nye,5*nze/4)) !LAS needs a bit of extra room for storage 
                else
                    allocate(efld(nxe,nye,nze))
                end if
            else
                allocate(efld(nxew,nyew,nzew))
            end if
        end if
        
        
        if (singletrue) then
        	if (is_isotropic) then
            	write(*,*) 'Running in Single-layer mode (isotropic).'
            else
            	write(*,*) 'Running in Single-layer mode (anisotropic).'
            end if
        else
            write(*,*) 'Running in Multi-layer mode.'
        end if



		!---------print some error messages from incorrect input --------


		if(nxew > nxe .or. nyew > nye .or. nzew > nze) then
			write(istat,*) "error: working field is bigger than generated field"
			write(istat,*) "working field dimensions (x,y,z) are:" ,nxew,nyew,nzew
			write(istat,*) "generated field dimensions are:" ,nxe,nye,nze/anisotropy
			write(istat,*) "decrease working site size, element size,"
			write(istat,*) "or increase generated dimensions."
            read(*,*)
			stop
		end if
		
		if(distribution /= 'n' .and. distribution /= 'l' .and. distribution /= 'b') then
			write(istat,*) "error: Incorrect field distribution chosen."
			write(istat,*) "Must choose 'n' for normal, 'l' for lognormal, or 'b' for beta"
			write(istat,*) "You have chosen:" ,distribution
            read(*,*)
			stop
		end if


		if( varfnc .ne. 'dlavx3' .and. varfnc .ne. 'dlsep3' .and. varfnc .ne.  'dlspx3' .and. varfnc .ne. 'dlafs3' .and. varfnc .ne. 'dlsfr3' ) then

			write(istat,*)'*** ERROR: unknown covariance function:',varfnc
			write(istat,*)'choose one of the following: dlavx3, dlsep3, dlspx3, dlafs3, dlsfr3'
			write(istat,*)'dlavx3 is recommended'
            read(*,*)
			stop
		end if


		close(inn)
        

        ! If specified, read in multi-layer soil profiles from a file.
        if(uselcoords) then
            deallocate(lmean_ave,lsd)
            open(678,file='input'//slash//'layer_data.txt')
            read(678,*) nbhmulti
            read(678,*) nlayer
            allocate(lmean_ave(nlayer),lsd(nlayer))
            read(678,*) (lmean_ave(j), j=1,nlayer)
            read(678,*) (lsd(j), j=1,nlayer)
            allocate(bhxymulti(2,nbhmulti),bhdepthsmulti(nlayer-1,nbhmulti))    
            read(678,*) (bhxymulti(1,j), j=1,nbhmulti)   !read in x coordinates
            read(678,*) (bhxymulti(2,j), j=1,nbhmulti)   !read in y coordinates
            read(678,*)
            do i = 1,nlayer-1
                read(678,*) (bhdepthsmulti(i,j), j=1,nbhmulti)
            end do
        end if
        
        if (.not. singletrue) then
            allocate(bfld(nlayer+1,nxew,nyew))
        end if
		
        
        ! if vertical random fields are going to be added to the multi-layer soil for the site investigation component, reduce the overall
        ! variability by 2/3. Note that if both the 2D random fields used for stiffness values and the 1D random fields are normally distributed,
        ! a reduction of 1/sqrt(2) would be used to achieve the original variability. This slightly lower value of 2/3 is arbitrary, but helps
        ! to account for the higher resulting variability when two lognormal distributions are combined.
        if (rand_realisations == 3) then
    	   lsd = 2 * lsd / 3
        else if (rand_realisations == 0) then! if zero is specified, force variability to zero (this probably isn't necessary, but can't hurt)
        	lsd = 0 
        end if 
   
   		!conver the normal statistics into lognormal for multiple layers
        allocate(lmeans(nlayer,nrep_MC),lmean_ln(nlayer),lsd_ln(nlayer))
        do j = 1,nlayer
            pvr = log(1.0 + lsd(j) ** 2 / lmean_ave(j) ** 2)
            lmean_ln(j) = log(lmean_ave(j)) - 0.5 * pvr
            lsd_ln(j) = sqrt(pvr)
        end do
        
        
        ! This stores the mean stiffness for each layer in each Monte Carlo realisation.
        do j = 1,nlayer
            if(rand_realisations == 1 .and. .not. all(lsd <= 0)) then 
            !generate lognormally-distributed values as the uniform E for each soil layer and realisation
                call vnorm(lmeans(j,:), nrep_MC)
                lmeans(j,:) = exp(lmean_ln(j) + lmeans(j,:) * lsd_ln(j))
            else  !if specified or the standard deviations are zero, just use the mean directly
                lmeans(j,:) = lmean_ave(j)
            end if
        end do
            
   


	!-------------------------------------


	end subroutine
	
	
		!read in pile data. Calculate pile coordinates
	
		subroutine readin_pile(femrad,femdepth,prad,npdepths,pdepths,cg_tol,cg_limit,preps,plocation,nels,costvals,failurevals, & !nbh,bhdepth
		buildingweight,load_con,rel_loads,num_loads,dz,detdisp,femvech,femvecv,pilecost,difftol,usepie,regmesh,abstol,emean,goodpiles,goodcases,istat)
		
		
		integer, intent(out) :: femrad,femdepth !radius and depth of soil around pile (elements)
		integer, intent(out) :: prad(2)			!radius of pile (elements)
		integer, intent(out) :: npdepths			!number of pile depths to test
		integer, allocatable, intent(out) :: pdepths(:) !pile depths (elements)
		real(8), intent(out) :: cg_tol 		!max no. iterations for FEM subroutine
		integer, intent(out) :: cg_limit 		!tolerance for FEM
		integer :: poffset(2)		!x,y offset of pile from corner of soil (elements)
		integer, intent(out) :: preps(2)		!No. piles in x,y directions - repetitions
		
		integer, intent(out) :: nels									!total number of soil elements in PIE method
		real, intent(out) :: difftol !differential pile settlement limit
		integer, allocatable, intent(out) :: load_con(:) !connectivity array specifying which rel_load is applied to each pile index
		real, allocatable, intent(out) :: rel_loads(:) !relative proportion of pile loads for various load cases
		integer, intent(out) :: num_loads !number of pile load options
		real, intent(out) :: failurevals(2) !lower and upper bounds for differential settlement failure limits to be linearly interpolated. Recommend 0.003 and 0.009.
		real, intent(out) :: costvals(2) !Costs assocated with the lower and upper bound of the above limits
        real, intent(out) :: pilecost !cost for pile, per metre depth
        logical, intent(out) :: usepie !whether to use FEM approximation such as the PIE method for single layers and Mylonakis and Gazetas method for multi-layer soils
        logical, intent(out) :: regmesh !use regular or variable FEM mesh (recommend latter - false)
        real, intent(out) :: buildingweight
		real(8), intent(in) :: dz
        real, intent(out) :: abstol !absolute settlement tolerance  (mm). Positive values are used directly. Negative values make it a function of differential settlement * minimum pile spacing (i.e., automatically determined).
        real, intent(in) :: emean !mean soil stiffness

        
        integer, allocatable, intent(out) :: goodpiles(:)        !vector of piles that have an associated applied load
    	integer, intent(out) :: goodcases   
        
        !temporary fem input information
        integer, allocatable :: setnumh(:),setnumv(:)
        real, allocatable :: seth(:),setv(:)
        integer femh,femv
		
		
		!some variables based on the input
		integer, allocatable, intent(out) :: plocation(:,:)	!x,y coordinates of each pile, to reflect 2D pile grid
		real, allocatable, intent(out) :: femvech(:),femvecv(:)
		
		!local variables
		integer i,j,x,y,count,bh,test,inv			!loop counters
		integer inn
		real pspacing(2) !pile spacing in the x,y directions based on number of piles and building size
		integer foundationwidth(2) !foundation area dimensions width in x,y directions (elements)
		real applied_load !uniformly distributed load over 1 square metre. Recommend 8 kPa (5 dead load + 3 live load)
		real power_coeffs(2) !a and b power equation coefficients for failure cost = area * a * number of floors ** b.
        logical :: givecoords !specify coordinates of piles directly instead of a regular grid
        integer npl !number of piles
        integer :: numfloors	!number of floors in building
		
		real, allocatable, intent(out) :: detdisp(:)
        
        logical exists !check if input files exist
	
		real :: temp									!placeholder to do with input formatting
        
        integer num_args !number of command line input arguments
        character(1000) str2

        integer, intent(in) :: istat ! debug file

        !If number of arguments is greater than zero, assume that multiple instances of the program are running in parallel.
        !Don't save the pile locations, as it's duplicate information and may cause output conflicts.
        num_args = command_argument_count()
	
		inn= 500


		open(inn, file='input'//slash//'pile_input.txt',status='old')



		! read in simulation options
		read (inn,*) 
		!read (inn,*) femrad
		!read (inn,*) femdepth
        !read (inn,*) femh
        !read (inn,*) femv
        !allocate(setnumh(femh),seth(femh),setnumv(femv),setv(femv))
        !read (inn,*) (setnumh(i), i=1,femh)
        !read (inn,*) (setnumv(i), i=1,femv)
        !read (inn,*) (seth(i), i=1,femh)
        !read (inn,*) (setv(i), i=1,femv)
        
        !hard-coded example for a 0.5m wide pile in a 40x40x40m block of soil with 0.5m cubic soil elements. Max pile depth 20m.
        !FEM mesh elements get larger further away from the pile
        femh = 13		!number of unique FEM elements in the horizontal direction
        femv = 4		!number of unique FEM elements in the vertical direction
        allocate(setnumh(femh),seth(femh),setnumv(femv),setv(femv))
        setnumh = [1,1,1,1,1,1,3,1,1,1,1,1,1]	!repetition of FEM element sizes in the horizontal direction
        setnumv = [50,5,4,2]			!repetition of FEM element sizes in the vertical direction
        seth = [6.5,4.5,3.,2.,1.5,1.,0.5,1.,1.5,2.,3.,4.5,6.5]		!series of FEM element sizes in the horizontal direction
        setv = [0.5,1.,1.5,2.]				!series of FEM element sizes in the vertical direction
        
		read (inn, *) pile_foundation
        read (inn,*) prad(1),prad(2)
        read (inn,*) difftol
        read (inn,*) abstol
		read (inn,*) npdepths
        if (pile_foundation) then ! pile foundation
            femrad = nint(npdepths/dz)
            femdepth = femrad*2
            npdepths = npdepths + 1
            allocate(pdepths(npdepths))
            pdepths = [(nint(i/dz),i=0,npdepths-1)]
        else ! pad footing
            femrad = maxval(prad) * 5 ! make FEM boundary 5x the distance from the pad of the maximum footing width
            femdepth = maxval(prad) * 4
            npdepths = ceiling(real(maxval(prad))/2)
            allocate(pdepths(npdepths))  ! number of footing widths to assess
            pdepths = [(i, i=1,maxval(prad),2)]
        end if

        if ((.not. singletrue) .and. (.not. pile_foundation)) then
			write(istat,*) "error: pad footings aren't supported in multi-layer mode"
			stop
        end if
        
							
		!read (inn,*) (pdepths(i),i=1,npdepths)		!The footing depth values to analyse
		read (inn,*) cg_tol,cg_limit
        !read (inn,*) usepie
        usepie = .true.
        !read (inn,*) regmesh
        regmesh = .true.
        read (inn, *) multitype
        read (inn,*) givecoords
		read (inn,*) npl

        if ( givecoords) then !if coords are given directly, read them in 
            numpiles = npl
            preps(1) = numpiles
            preps(2) = 1
            allocate(plocation(npl,2),goodpiles(npl),load_con(npl))
            read (inn,*) plocation(:,1)
            read (inn,*) plocation(:,2)
            read(inn,*) (load_con(x),x= 1, preps(1))
            read (inn,*)
        else
            read (inn,*)
            read (inn,*)
            read (inn,*)
            read (inn,*) preps(1),preps(2)
            npl = preps(1)*preps(2)
            numpiles = npl
            allocate(plocation(npl,2),goodpiles(npl),load_con(npl))
        end if
        read (inn,*) poffset(1),poffset(2)
		read (inn,*) foundationwidth(1),foundationwidth(2)
        if ( givecoords) then !if coords are given directly, read them in 
            do y= 1,preps(2)
		        read(inn,*) 
            end do
        else
		    do y= 1,preps(2)
		        read(inn,*) (load_con(x),x= (y-1)*preps(1) + 1, y*preps(1))
            end do
        end if
		if(num_args >= numpiles + 9) then
			i=0
			do y= 1,preps(2)
				do x = 1,preps(1)
					i = i+1
					load_con(i) = getarg_real()			
				end do
			end do
		end if
		num_loads = maxval(load_con)
		allocate(rel_loads(num_loads))
		read (inn,*) rel_loads
		read (inn,*) buildingarea
		if(num_args > numpiles + 10) buildingarea = getarg_real()	
		read (inn,*) numfloors
		read (inn,*) applied_load
		read (inn,*) failurevals
		!read (inn,*) power_coeffs
        read (inn,*) costvals(2)
        read (inn,*) pilecost
		
		!Calculate failure cost limits
		costvals(1) = 0 !lower limit is always 0 (zero cost for zero failure)
        !---note that the building cost is now given directly as opposed to calculating it from the area and number of floors---
		!costvals(2) = buildingarea * power_coeffs(1) * numfloors ** power_coeffs(2) !upper limit follows this trend based on the number of floors, and is linearly proportional to area
		
		!Calculate building weight
        if (singletrue) then !do single layer analysis prep
		    buildingweight = applied_load * buildingarea * numfloors 
		    buildingweight = buildingweight / emean !Instead of changing mean soil stiffness, we are changing the applied load instead. Bigger mean = less settlement = 'smaller load'
		else
            buildingweight = applied_load * buildingarea * numfloors !If multi-layered mode, use actual building weight directly
        end if
            
		allocate(detdisp(npdepths))
        
        !---build vectors of FEM mesh information when variable-mesh FEM is used instead of the PIE method for CK settlement---
        !NOTE: The variable-mesh mode has been commented out for the publicly-released version as it's too time consuming to use.
        allocate(femvech(int(sum(setnumh))),femvecv(int(sum(setnumv))))
 		count = 0
         do i=1,femh
             do j = 1,setnumh(i)
                 count = count+1
                 femvech(count) = seth(i)
 			end do
 		end do
 		count = 0
         do i=1,femv
             do j = 1,setnumv(i)
                 count = count+1
                 femvecv(count) = setv(i)
 			end do
 		end do
		
		!create a vector of x,y pile locations, looping in x direction fastest
		!if statements ensure that piles are within building plan (negate pile width)
		if (.not. givecoords) then
		    count = 0
            do x = 1,2  !get pile spacing in both directions
                if (preps(x) > 1) then
                    pspacing(x) = real(foundationwidth(x)-prad(x))/(preps(x)-1)          
                else
                    pspacing(x) = 1 !if there's only one row/column, make it equal to itself
		        end if
            end do
		    do y=1,preps(2)
			    do x=1,preps(1)
				    count = count + 1
				    plocation(count,1) = poffset(1) + nint((x-1) * pspacing(1))
				    plocation(count,2) = poffset(2) + nint((y-1) * pspacing(2))
			    end do
            end do
        end if
        
        
        !process which piles have an applied load (load_con relates pile index to relative load, with 0 having no load associated)
		goodcases = 0
		do i = 1,preps(1)*preps(2)
			if(load_con(i) > 0) then
				goodcases = goodcases + 1
				goodpiles(goodcases) = i
			end if
		end do
        
        
		
		nels = (femrad*2 + prad(1)) * (femrad*2 + prad(2)) *femdepth
        
                !output pile locations in metres. This is pretty much a constant, so would only have to be output once, i.e at the start of the EA
        if(num_args == 0) then
        	write(str2,'(A,I0,A,I0,A)') 'si_results/pile_locations-',goodcases,'_BArea-',nint(buildingarea),'.txt'
            open(505,file=str2)
            write(505,'(A8,X,A8)') 'X','Y'
            do i = 1,goodcases
                write(505,'(F8.2,X,F8.2)') dz*plocation(goodpiles(i),:)
            end do
            close(505)
        end if



	end subroutine
	
        !read in site investigation data. Process later.
	
		subroutine readin_si(soffset,swidth,sstep,ntest,testerrors,sampfreq, & !nbh,bhdepth
		testcost,testnames,use_CI,conf_int,add_errors,in_tests,in_depths,in_reductions,bhnums,read_si,in_sdev,in_percentile,save_effE)
		
		integer, intent(out) :: soffset(2)		!x,y offset of borehole from corner of soil (elements)
		integer, intent(out) :: swidth(2)		!x,y dimensions of site investigation area
        integer, intent(out) :: sstep(2)            !borehole step size in each dimension for the heat map mode
		!integer, intent(out) :: nbh,bhdepth,ntest	!No. boreholes, borehole depth (elements), No. tests
		integer, intent(out) :: ntest !number of test types
		real, allocatable, intent(out) :: testerrors(:,:) !Transformation error, bias error, measurement error for each test
		integer, allocatable, intent(out) :: sampfreq(:) !Sampling frequency for each test (elements)
        real, allocatable, intent(out) :: testcost(:) !cost per metre depth for each test
        character(4), allocatable, intent(out) :: testnames(:) !names for each test
        logical, intent(out) :: use_CI !use a confidence interval to truncate unrealistic sample values. Recommended.
		real, intent(out) :: conf_int !confidence interval z score for the above truncation, recommend 99% (2.576) 
        logical, intent(out) :: add_errors !add random errors to tested samples (highly recommended for realistic sampling)
		integer, intent(out) :: read_si !If true, will read in the configuration of individual investigations using the 'input/si.txt'. Otherwise will generate a set of all combinations of investigation inputs
		real, intent(out) :: in_sdev,in_percentile !standard deviations below the mean in the Geometric "SD" reduction method. Percentile to use for the percentile reduction method. Recommend 1 and 0.25 currently.
		
		integer, allocatable, intent(out) :: in_tests(:),in_depths(:) !input list of tests, depths and reduction methods to investigate in analysis
		character(2), allocatable, intent(out) :: in_reductions(:)		!input list of reduction methods
		integer, allocatable, intent(out) :: bhnums(:) !A list of the number of boreholes in each investigation,	
        integer, intent(out) :: save_effE !output effective soil modulus of each investigation in each MC realisation. Not needed for your project. 0 = No. 1 = save 2= load.
		
		!local variables
		integer i,j,x,y,count,bh,test,inv			!loop counters
		integer inn
		integer nbh,n_depths, n_red,n_test !sizes of the lists of borehole number, depth, reduction and test inputs to analyse
		integer temparg
		
	
		real :: temp									!placeholder to do with input formatting
	
		inn= 500

		open(inn, file='input'//slash//'si_input.txt',status='old')

		read (inn,*)
		read (inn,*) soffset(1),soffset(2)
		read (inn,*) swidth(1),swidth(2)
        read (inn,*) sstep(1),sstep(2)    
        !sstep = [1,1]
		if(num_args >= numpiles + 11) then
			swidth(1) = getarg_real()
			swidth(2) = swidth(1) !assume it's square for now
		end if
		read (inn,*) ntest
        read (inn,*) use_CI
        read (inn,*) conf_int
        read (inn,*) add_errors
		read (inn,*)
		allocate(testerrors(ntest,3),sampfreq(ntest),testcost(ntest),testnames(ntest))
		do i=1,ntest
			read (inn,*) testerrors(i,1),testerrors(i,2),testerrors(i,3),sampfreq(i),testcost(i),testnames(i)
		end do
		read(inn,*)  in_sdev,in_percentile
        !read (inn,*) save_effE 
        save_effE  = .false.
		
		read(inn,*)
		read(inn,*) read_si !If true, will read in the configuration of individual investigations using the 'input/si.txt'. Otherwise will generate a set of all combinations of investigation inputs
		
		!Read in the number of values for each parameter
		read(inn,*) nbh
		read(inn,*) n_test
		read(inn,*) n_red
		read(inn,*) n_depths
		
		!specify some internal things
		allocate(bhnums(nbh),in_tests(n_test),in_depths(n_depths),in_reductions(n_red))
        read (inn,*) (bhnums(i),i=1,nbh)		!Read in list of borehole numbers
		read (inn,*) (in_tests(i),i=1,n_test)
		if(num_args >= numpiles + 12) in_tests(1) = getarg_real()
		read (inn,*) (in_reductions(i),i=1,n_red)
		if(num_args >= numpiles + 13) then
			temparg = getarg_real()
			in_reductions(1) = rednames(temparg)
		end if
		read (inn,*) (in_depths(i),i=1,n_depths)
		if(num_args >= numpiles + 14) in_depths(1) = getarg_real()
		if(num_args >= numpiles + 15) bhnums(1) = getarg_real()
		
        !Check for validity of reduction method
        do i = 1,n_test
            if(in_tests(i) > ntest) then
                write(*,*) 'Error, invalid test number chosen:',in_tests(i)
                write(*,*) 'Please choose a value between 1 and',ntest
                write(*,*) 'Press enter to continue.'
                read(*,*)
                stop
            end if
        end do
        
        !Check for validity of test type
        do i = 1,n_depths
            if(in_depths(i) > nzew) then
                write(*,*) 'Error, borehole is too deep:',in_depths(i),'elements'
                write(*,*) 'Please choose a value between 1 and',nzew
                write(*,*) 'Press enter to continue.'
                read(*,*)
                stop
            end if
        end do
        
        do i = 1,n_red
            if(.not. any(in_reductions(i) == rednames)) then
                write(*,*) 'Error, invalid reduction method chosen:',in_reductions(i)
                write(*,*) 'Please choose among:',rednames
                write(*,*) 'Press enter to continue.'
                read(*,*)
                stop
            end if
        end do

		close(inn)


	end subroutine
	
	
	
	!get input related to usage of the evolutionary algorithm
	subroutine evol_input(nrep_MC,maxEA,runmode,EAseed,soilseeds,costmetric,deterministic,max_cons,popsize,perc_error,mut_rate,frac_keep,finaloutput,mut_lims,imut,siea_mode,num_elite,stopmode,unistart,local_opt)

	
	character(100) :: str
	character(100), intent(out) :: runmode
	integer inn,i
	integer, intent(out) :: nrep_MC		!number of monte carlo realistions
	integer, intent(out) :: maxEA		!max number of evolutionary algorithm iterations
    integer, intent(out) :: max_cons     !max number of consequtively equal EA performance before stopping
    integer, intent(out) :: popsize     !population size for evolutionary algorithm
    real, intent(out) :: perc_error     !percentage error for EA stop critera (out of 100)
    real, intent(out) :: mut_rate,mut_lims(2)       !mutation rate for EA, as well as the min and max mutations if dynamic mutation is specified
    integer, intent(out) :: imut !mutation mode; 1 = constant, 2 = adapt based on fitness, 3 = adapt based on proximity in parameter space
    real, intent(out) :: frac_keep      !fraction of population to keep after each EA iteration
	integer, allocatable, intent(out) :: soilseeds(:) !vector of soil seeds, one for each MC realisation
    integer, intent(out) :: costmetric !Performance metric to use. -1 = failure cost, 0 = probability of failure, 1 or more = average differntial settlement to the power of the costmetric value (higher values more heavily penalise excessive values)
    integer, intent(out) :: deterministic !whether a site investigation is assessed for performance once (true) or optimised in an evolutionary algorithm (false)
	integer, intent(out) :: finaloutput !output mode for EA. 1 = final only, 2 = save best of each EA generation, 3 = save full final population, 4 = save everything
    integer soilseed,EAseed !random seeds for virtual soil and genetic algorithm. Set to 0 to be random (clock based) or higher for consistent random numbers. VERY strongly recommend soilseed be kept as a positive integer (e.g 100).
    !character(1000),intent(out) :: datafolder
    integer, intent(out) :: siea_mode !manner of controlling the borehole locations in the EA. 1 Means start across full soil. 2 Means Start within building footprint. 3 Means ALWAYS keep within building footprint.
	integer, intent(out) :: num_elite !number of elite population members that are carried over into the next generation, unmutated
    logical, intent(out) :: stopmode !.true. means wait for EA to have converged (remain unchanged for a number of realisations), .false. means wait for a number of iterations past the current GLOBAL optimum before stopping (in case fitness starts increasing).
    logical, intent(out) :: unistart !true means have random initial population according to uniform distribution, otherwise attempt to have normal offsets from a regular grid
    logical, intent(out) :: local_opt !whether to do a quick GA analysis on the final solution to see if a better one is nearby
    
    !external functions
	real randu
    integer iseed

	
	
		inn= 500
		open(inn, file='input'//slash//'EA_input.txt',status='old')
        
		read(inn,*) deterministic
        read(inn,*) nrep_MC
        read(inn,*) costmetric
        read(inn,*) soilseed,EAseed
        read(inn,*) runmode
        read(inn,*) finaloutput
        read(inn,*) datafolder
        read(inn,*)
		read(inn,*) maxEA
        read(inn,*) max_cons
        read(inn,*) stopmode
        read(inn,*) popsize
        read(inn,*) perc_error
        read(inn,*) mut_rate,mut_lims(1),mut_lims(2)
        read(inn,*) imut
        read(inn,*) frac_keep
        read(inn,*) num_elite   
        read(inn,*) siea_mode
        read(inn,*) unistart
        read(inn,*) local_opt
       
		
		allocate(soilseeds(nrep_MC))
		
		
		!initialise random number generator
		soilseed = iseed(max(0,soilseed)) ! !ensure kseed is at least zero 
        
		
		!save vector of unique random seeds to be used for each soil
		do i=1,nrep_MC
			soilseeds(i) = randu(0) * 1234567890	!zero does not reset the seed, rather it continues from the previous value
		end do
		soilseed = randu(soilseed) * 1234567890 !reset the seed regardless of number of realisations
		
		close(inn)
        
        !--- check inputs for validity ---
        
        if (costmetric < -4) then
            write(*,*) 'Warning, invalid fitness cost type chosen. Choose:'
            write(*,*) '-4 for expected failure cost'
            write(*,*) '-3 for expected total cost'
            write(*,*) '-2 for average differential settlement'
            write(*,*) '-1 for probability of failure'
            write(*,*) 'n>=0 for n geometric standard deviations above the geometric mean'
            write(*,*) 'You have chosen:',costmetric
            write(*,*) 'Hit any key to exit.'
            read(*,*)
            stop
        end if
        
 
        
        if (imut < 1 .or. imut > 3) then
            write(*,*) 'Warning, invalid mutation mode. Choose:'
            write(*,*) '1 for constant'
            write(*,*) '2 for dynamic based on fitness'
            write(*,*) '3 for dynamic based on proximity in parameter space'
            write(*,*) 'You have chosen:',imut
            write(*,*) 'Hit any key to exit.'
            read(*,*)
            stop
        end if
        
        if(all(['PP','SI','CK','AL'] .ne. runmode)) then
        	write(*,*) 'Warning, invalid single layer processing: '//runmode//'. Choose:'
        	write(*,*) "'PP' to Pre-Process the foundation."
        	write(*,*) "'CK' to get true pile settlements using Complete Knowledge of the soil."
        	write(*,*) "'SI' to do Site Investigation analysis."
        	write(*,*) "'AL' (all) to automatically do all modes required up to and including SI."
            write(*,*) "Press any key to exit"
            read(*,*)
            stop
        end if
        
        if(siea_mode<1 .or. siea_mode>3) then
            write(*,*) 'Error, invalid SI EA control mode:',siea_mode
            write(*,*) 'Enter a value between 1 and 3, inclusive.'
            write(*,*) "Press any key to exit"
            read(*,*)
            stop
        end if
        
        !this will crash the program if violated
        if(mod(CEILING((popsize-frac_keep*popsize)),2) /= 0 ) then
            write(*,*) 'Error, number of population members to keep must be an equal number.'
            write(*,*) 'Suggest increasing the population size of changing the proportion of'
            write(*,*) '    population members to keep.'
            write(*,*) 'Value currently:', CEILING((popsize-frac_keep*popsize))
            read(*,*)
            stop
        end if
        
        if (deterministic < 0 .or. deterministic > 3) then
            write(*,*) 'Error, invalid program run mode. Choose:'
            write(*,*) '0 for approximate determistic pile design'
            write(*,*) '1 for fixed site investigation analysis'
            write(*,*) '2 for heat map site investigation analysis'
            write(*,*) '3 for evolutinoary site investigation analysis'
            write(*,*) 'You have chosen:',deterministic
            write(*,*) 'Hit any key to exit.'
            read(*,*)
            stop
        end if
        
        
		
	end subroutine
	
	
	
	subroutine read_ck(nrep_MC,npdepths,preps,ck_set,detdisp,pdepths,prad)
	
	
	real, allocatable, intent(out) :: ck_set(:,:,:)
    real, allocatable :: pilecurves(:,:)
	integer counter,x,y,pd,i
	character(1000) str2
    real, intent(out) :: detdisp(:)
    integer, intent(in) :: prad(2)
    
    integer, parameter :: curvefactor=10
	real pdepthssi((npdepths-1)*curvefactor+1)
    integer, intent(in) :: pdepths(:)			!pile depths to analyse (elements)
	
	
	integer, intent(in) :: nrep_MC,npdepths
	integer, intent(in) :: preps(2)
    real temp
    
    integer,parameter :: scalefactor=10
	
	logical exists
    
    write(str2,'(A,A,I0,A,F4.2,A)') trim(datafolder),'settlement_prad-',prad(1),'_esize-',dz,'.txt'
    open(667, file= str2,status='old')
	read(667,*) 
	do i = 1,npdepths
		read(667,*) temp,detdisp(i)
	end do
	close(667)
    
	
	allocate(ck_set((npdepths-1)*curvefactor+1,nrep_MC,preps(1)*preps(2)))
    allocate(pilecurves(npdepths,nrep_MC))
    
    
    !build higher resolution settlement curve to help ensure that the inverse of the function matches it closely (boundary condition interpretation means this isn't neccessarily the case)
	pdepthssi(1) = 0 !first point
	counter = 1
	do i=1,npdepths-1 !loop through pairs of original x points
		do x = 1,curvefactor !linearly interpolate new x points between pair
			counter = counter + 1
			pdepthssi(counter) = x*real(pdepths(i+1) - pdepths(i))/curvefactor + pdepths(i)	
		end do
    end do

    
	
        
	    do x = 1,preps(1)*preps(2)
		    if(singletrue) then !single layer
                write(str2,'(A,A,I0,A,I0,A,F4.2,A,I0,A,I0,A,I0,A)') trim(datafolder),'ck_pile-',x,'_prad-',prad(1),'_esize-',dz,'_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy,'.txt'
            else
                write(str2,'(A,A,I0,A,I0,A,F4.2,A,I0,A,I0,A,I0,A)') trim(datafolder),'ck_pile-',x,'_prad-',prad(1),'_esize-',dz,'_nlayers-',nlayer,'_ratio2l-',nint(lmean_ave(2)/lmean_ave(1)),'_depth2l-',nint(dz*ldepths(1)),'.txt'
		    end if
			open(500,file=str2,status='old')
			do i = 1,nrep_MC
				read(500,'(100000(E14.7,X))') (pilecurves(pd,i),pd=1,npdepths)
			end do
			close(500)
            call akima(real(pdepths), pilecurves, pdepthssi, npdepths, 1, nrep_MC, 1, size(pdepthssi), ck_set(:,:,x)) !interpolate high-res curve ck_set(:,1,1) pilecurves(:,1)

        end do

    
    !do x = 1,preps(1)*preps(2)
    !    write(str2,'(A,I0,A)') 'results/ck_test',x,'.txt'
    !    open(500,file=str2)
!		do i = 1,nrep_MC
!			write(500,'(10000(E14.7,X))') (pilecurve(pd,i,x),pd=1,npdepths)
!		end do
!		close(500)
 !   end do

    
    !scale CK pile settlement according to current soil stiffness
    !ck_set = ck_set !/ sdata(1)  !Don't scale the curves anymore. Scale the applied load instead.
    
        

	
	deallocate(pilecurves)
	
    end subroutine
    
    
    subroutine input_SI(nbh,sampfreq,inv_coords,inv_depths,inv_nbh,inv_test,ninv,inv_reduction,percentile,s_dev,si_performance,dz,testnames,deterministic,in_sdev,in_percentile,soildsc,goodcases,read_si)	 !output variables

	!define the site investigation information. Hardcoded to do a combination of some tests for each number of boreholes and reduction method


        !input
		integer, allocatable, intent(in) :: sampfreq(:) !Sampling frequency for each test (elements)
		real, allocatable :: si_performance(:)
        real(8), intent(in) :: dz
        character(4), intent(in) :: testnames(:)
        integer, intent(in) :: deterministic !SI mode. 1 = 'deterministic' analysis of investigation performance, 2 = 'heat map', 3 = optimise investigations with EA
        real, intent(in) :: in_sdev,in_percentile
        character(100), intent(in) :: soildsc
        integer, intent(in) :: goodcases
        integer, intent(in) :: read_si !1 = read borehole attributes and locations, 2 = also read tests and depths on a per-borehole basis

		
		!output
		integer, allocatable, intent(out) :: inv_coords(:,:,:)			!x,y coordinates of boreholes for each investigation
		integer, allocatable, intent(out) :: inv_depths(:,:,:) 
		integer, allocatable, intent(out) :: inv_nbh(:)						!number of boreholes in each investigation
		integer, allocatable, intent(out) 	:: inv_test(:,:)					!test type for each investigation
        integer, allocatable, intent(out) 	:: inv_reduction(:)					!reduction method for each investigation
        real, allocatable, intent(out) 	:: percentile(:)					!percentile for use in reduction method
        real, allocatable, intent(out) 	:: s_dev(:)					!geometric standard deviation below geometric mean in reduction method
		integer, intent(out) :: ninv					!total number of investigations
        integer, intent(out) :: nbh
		
		!local variables
		integer i,x,y,count,bh,test,inv,j,red			!loop counters
        integer, parameter :: n_red = 5     !technically 6 reduction methods, although let's not use 'minimum' due to it being excessively conservative
	    character(3) :: temptest
        character(2) :: tempred
        real :: tempdepth
        real, allocatable :: inv_ctemp(:,:,:)
        character(1000) str2
        
        
        !will first try to open "si.txt". If this doesn't exist, use a longer name based on the soil description
        inquire(file='input'//slash//'si.txt',exist=exists) !check if the file exists 
		if(.not. exists) then
            write(str2,'(A,I0,A,I0,A,A)') 'input/Population-stats_piles-',goodcases,'_BArea-',nint(buildingarea),trim(soildsc),'.txt'
            inquire(file=str2,exist=exists) 
            if(.not. exists) then
                write(*,*) 'Error, file does not exist:'
                write(*,*) 'input'//slash//'si.txt'
                write(*,*) 'Press enter to continue.'
                read(*,*)
                stop
            else
                open(666,file=str2,status='old')
            end if
        else
            open(666,file='input'//slash//'si.txt',status='old')
        end if
        
        !check that bh_depths and bh_depths files exist
        if(read_si == 3) then
            inquire(file='input'//slash//'bh_depths.txt',exist=exists) !check if the file exists 
		    if(.not. exists) then
                    write(*,*) 'Error, file does not exist:'
                    write(*,*) 'input'//slash//'bh_depths.txt'
                    write(*,*) 'Press enter to continue.'
                    read(*,*)
                    stop
            end if
        
            inquire(file='input'//slash//'bh_tests.txt',exist=exists) !check if the file exists 
		    if(.not. exists) then
                    write(*,*) 'Error, file does not exist:'
                    write(*,*) 'input'//slash//'bh_tests.txt'
                    write(*,*) 'Press enter to continue.'
                    read(*,*)
                    stop
            end if
            
            !open files
            open(728,file='input'//slash//'bh_depths.txt',status='old')
            open(729,file='input'//slash//'bh_tests.txt',status='old')
            
            !first line should be blank
            read(728,*)
            read(729,*)
        
        end if

        
        !read in number of investigations and max number of boreholes
        read(666,*) ninv
        read(666,*) nbh
        read(666,*)
        read(666,*)
        read(666,*)
        read(666,*)
        
        allocate(inv_coords(ninv,nbh,2))		!x,y coordinates for each borehole in each investigation
        allocate(inv_ctemp(ninv,nbh,2))
		allocate(inv_depths(nbh,ninv,3)) 			!sample depth information for each investigation
		allocate(inv_nbh(ninv))					!number of boreholes in each investigation
		allocate(inv_test(nbh,ninv))				!test type for each investigation
        allocate(inv_reduction(ninv))			!reduction for each investigation
        allocate(percentile(ninv))
        allocate(s_dev(ninv))
		allocate(si_performance(ninv))
        
        
        !set some defaults
        inv_depths(:,:,1) = 1
        s_dev = in_sdev
        percentile = in_percentile
        
        !loop through investigations and get site investigation parameters; number of investigations, test type, reduction method, and depth
        
        do inv = 1,ninv

                read(666,*) inv_nbh(inv), temptest, tempred, tempdepth
                
                if(read_si == 2) then !process borehole depths and tests on a per-investigation basis
                
                    !convert depth from metres to elements
                    inv_depths(:,inv,2) = nint(tempdepth/dz)

                    !get test number from test name
                    do i = 1,size(testnames)
                        if( temptest == testnames(i) ) then
                            inv_test(:,inv) = i
                            exit
                        end if
                    end do
                    
                    !get test sampling frequency
                    inv_depths(:,inv,3) = sampfreq(inv_test(1,inv))
                
                else    !read in borehole depths and tests on a per-borehle basis
                    
                    read(728,*) inv_depths(:inv_nbh(inv),inv,2) !get depths
                    read(729,*) inv_test(:inv_nbh(inv),inv) !get tests
                    inv_depths(:inv_nbh(inv),inv,3) = sampfreq(inv_test(:inv_nbh(inv),inv)) !get sampling frequency based on tests
                    
                end if
                
                !get reduction method number from test name
                do i = 1,size(rednames)
                    if( tempred == rednames(i) ) then
                        inv_reduction(inv) = i
                        exit
                    end if
                end do
                
                
                

        end do
        
        
        close(666)
        if(read_si == 3) then
            close(728)
            close(729)
        end if
        
        !if a single, deterministic investigation is to be conducted, then read in the borehole locations
        !otherwise, the locations will be randomly generated at a later time in the evolutionary algorithm
        
        if (deterministic == 1) then
            !will first try to open "si_X/Ycoords.txt". If this doesn't exist, use a longer name based on the soil description
            !X coords
            inquire(file='input'//slash//'si_Xcoords.txt',exist=exists) !check if the file exists 
		    if(.not. exists) then
                write(str2,'(A,I0,A,I0,A,A)') 'input/Population-Xcoords_piles-',goodcases,'_BArea-',nint(buildingarea),trim(soildsc),'.txt'
                inquire(file=str2,exist=exists)
                if(.not. exists) then
                    write(*,*) 'Error, file does not exist:'
                    write(*,*) 'input'//slash//'si_Xcoords.txt'
                    write(*,*) 'Press enter to continue.'
                    read(*,*)
                    stop
                else
                    open(666,file=str2,status='old')
                end if
            else
                open(666,file='input'//slash//'si_Xcoords.txt',status='old')
            end if
            !Y coords
            inquire(file='input'//slash//'si_Ycoords.txt',exist=exists) !check if the file exists 
		    if(.not. exists) then
                write(str2,'(A,I0,A,I0,A,A)') 'input/Population-Ycoords_piles-',goodcases,'_BArea-',nint(buildingarea),trim(soildsc),'.txt'
                inquire(file=str2,exist=exists)
                if(.not. exists) then
                    write(*,*) 'Error, file does not exist:'
                    write(*,*) 'input'//slash//'si_Ycoords.txt'
                    write(*,*) 'Press enter to continue.'
                    read(*,*)
                    stop
                else
                    open(667,file=str2,status='old')
                end if
            else
                open(667,file='input'//slash//'si_Ycoords.txt',status='old')
            end if
        
            !Read in data from the files
            read(666,*)
            read(667,*)
            
            do inv = 1,ninv
		        read(666,*) (inv_ctemp(inv,j,1),j=1,inv_nbh(inv))
                read(667,*) (inv_ctemp(inv,j,2),j=1,inv_nbh(inv))
            end do
            
            
            !convert coordinates to unit of elements from meters
            inv_coords = nint(inv_ctemp/dz)
        
            close(666)
            close(667)
        
        end if

        
        
        
    end subroutine
            
	

end module
