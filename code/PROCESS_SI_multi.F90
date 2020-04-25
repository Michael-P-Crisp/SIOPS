
   
module process_SI_multi

use getperfmulti
use variables
use soilgen
use fem_prep
use ESETT
use getdiff
use SETUP_SI
use si_stats
use output_results

implicit none

contains
 
!Get site investigation performance information, done by the following steps:
!   1. Conduct site investigations across all realisations of random soils
!   2. Design piles according to those investigations, 
!   3. Get true performance in terms of differential settlement (requires program be run in 'CK' mode to get true pile behaviour first, unless FEM approximation is used).
!   4. Calculate failure costs based on differential settlement, along with pile construction and site investigation costs (add them together)
!   5. Site investigation performance is the average cost across realisations

!Alternatively, the investigation performance can be given in terms of the probability of failure, which is calculated from steps 1-3 above,
!and then counting the number of times the differntial settlement exceeds the minimum failure threshold (i.e failurevals(1) ), divided by the
!total number of realisations. Generally speaking, this metric would require fewer monte realisations (i.e program can be faster), as there are
!no extreme outliers which can skew the results, which would otherwise decrease stability. 
!On the other hand, the cost metric allows you to penalise excessively strong foundations and overly-thorough investigations to show what investigation
!is best overall in terms of minimal cost. 
!The probability of failure metric can maximize accuracy for a given investigation in terms of borehole depth and location, but across investigations,
!it will favour more conservative approaches. Still, this could form a Pareto front, where the knee of the curve could potentially be defined as optimal.

!WARNING: There could easily be the case where some Monte Carlo realisations are invalid. This occurs considering that the pile has a finite range of possible
! lengths, i.e 0 m - 20 m. If the soil conditions are such that the pile is "designed" to be less than the minimum or greater than the maximum, then the whole 
! group of piles cannot be used. The code has been written to keep track of and account for these invalid realisations (through NaNs - see main.F90), however if 
! the number of invalid realisations becomes excessively large for any reason, then the results will not be accurate.
    
    subroutine get_si_perf_multi(soilseeds, & !soil generation variables
                      ninv,nbh,in_tests,in_depths,in_reductions,inv_bh,inv_coords,inv_depths,inv_test,add_errors,use_CI,test_errors,inv_reduction,conf_int,percentile,s_dev,soffset,swidth,sstep &					!site investigation variables
                      ,npdepths,detdisp,plocation,pdepths,ck_set,nrep_mc,rel_loads,load_con,num_loads,preps,buildingweight,difftol,failurevals,costvals,si_performance,pilecost,testcost,testnames,costmetric,EA_generation,deterministic,abstol, &
        femvech,femvecv,prad,datafolder,usepie,npl2,goodpiles,goodcases,sdist,sumweights,extents,indices,CKheight,finaloutput,invcount,fcost, pcost, probfail, avediff, diffgeo)
                      
                      
! ----- soil generation variables (don't touch) ----
      integer iseed
      integer, intent(in) :: soilseeds(:)
      character(1000),intent(in) :: datafolder !a string representing the directory the data is stored in
      character(200) :: soildsc !soil description string containing soil attributes, to help distinguish between different cases, particularly single layer vs multiple layers



! ---- site investigation variables -----

	integer,  intent(in) :: inv_bh(:)						!number of boreholes in each investigation
	integer,  intent(in)	:: inv_coords(:,:,:)			!x,y coordinates of boreholes for each investigation
	integer,  intent(in) :: inv_depths(:,:,:) 				!sampling depth start, end, interval
	integer,  intent(in)	:: inv_test(:,:)					!test type for each investigation
	logical,  intent(in)	:: add_errors						!whether to add test errors or not
	real,  intent(in) :: test_errors(:,:)	!matrix of test error stats (num tests x num errors)
	integer,  intent(in) :: inv_reduction(:)				!reduction method for each test
	!real :: evals(nrep_mc,ninv)						!reduced value for each investigation
	real,  intent(in) :: percentile(:),s_dev(:) !percentile and standard deviation below mean for reduction methods
    real, intent(in) :: testcost(:)     !cost per metre for each test
    character(4), intent(in) :: testnames(:) !names for each test
	
	integer :: nbh,bhdepth,ntest	!No. boreholes, borehole depth (elements), No. tests
	integer :: ninv				!number of investigations
    integer,intent(in) :: in_tests,in_depths,in_reductions !number of tests, depths, and reductions investigated
	
	real :: conf_int !confidence interval for truncating SI sample values 
	logical :: use_ci 
    
    logical :: goodinv(nrep_MC,ninv) !array for valid site investigation performance information
    integer,allocatable :: goodloc(:,:) !array for valid site investigation performance locations
    integer :: invcount(ninv) !number of valid values for goodinv
    
    
    
    logical  allones !check whether every investigation is a single borehole
    
    integer, intent(in) :: soffset(2)		!x,y offset of borehole from corner of soil (elements)
	integer, intent(in) :: swidth(2)		!x,y dimensions of site investigation area
    integer, intent(in) :: sstep(2)            !borehole step size in each dimension for the heat map mode
    
    !triangulation pre-processing variables
    
    real :: node_xy_out(ninv,nbh,2)
    integer :: element_neighbor(ninv,3,4*nbh),triangle(ninv,3,4*nbh)
    integer :: zd_out(ninv,nbh)
    integer :: element_num(ninv) !number of triangles in interpolated surface
    integer:: extranum(ninv) !number of extrapolated points
    real(8) :: xyi(2,nxew*nyew) !coordinates of points to interpolate
    

! --- pile performance variables
    integer, intent(in) :: prad(2)			!pile width in x,y directions
	integer, intent(in) :: npdepths							!number of pile depths to analyse
	integer, intent(in) :: pdepths(:)			!pile depths to analyse (elements)
	real, intent(in) :: detdisp(:)				!deterministic settlement for each pile depth (1 kN applied load. Soil stiffness 1 MPa)
	integer, intent(in) :: preps(2)		!No. piles in x,y directions - repetitions
	integer, intent(in) :: plocation(:,:)		!Pile x,y coordinates
    real(8) :: plocation_dble(2,npl2*product(prad))
	
	real, intent(in) :: rel_loads(:)		!set of relative pile loads based on tributary area of building
	integer, intent(in) :: load_con(:)			!connectivity vector saying which pile has which load
    integer :: goodpiles(preps(1)*preps(2))        !vector of piles that have an associated applied load
    integer :: goodcases                               !number of good piles in goodpiles 
    integer :: goodcases_design               !number of good pile cases to design
    integer :: npl2 !good number of piles under certain conditions
	
	real, intent(in) :: buildingweight,difftol
    real :: settol
    real :: dist((preps(1)*preps(2))**2) !distances between piles
    integer, intent(in) :: indices(2,nxew,nyew)
	
	real, allocatable :: sides(:,:,:)          !sides(num_loads,nrep_mc,ninv)   !SI pile designs
    real, allocatable :: avesides(:,:)   ! average SI pile designs
    !real :: ckdes(preps(1)*preps(2),nrep_mc)   !CK pile designs preps(1)*preps(2),nrep_mc
	real, intent(in) :: ck_set(:,:,:) !true ck pile settlement for each realisation, pile and depth
    
    integer, intent(in) :: nrep_mc !number of MC realisations
    integer, intent(in) :: num_loads !number of load cases
    
    !values to get the actual total length of piles based on sum of pile cases
    integer num_loadSI
    integer load_conSI(num_loads)
    
    !expanded relative load information
    integer load_con2(preps(1)*preps(2))
    real rel_loads2(preps(1)*preps(2))
    
    real rel_loadsCK(preps(1)*preps(2)) !relative loads for CK piles
    integer load_conCK(preps(1)*preps(2))
    
    real, intent(in) :: failurevals(2), costvals(2) !2 x,y coordinates defining the line to interpolate costs from differential settlement
    real, intent(in) :: pilecost !cost for pile, per metre depth
    real, intent(in) :: abstol !absolute settlement tolerance  (mm). Positive values are used directly. Negative values make it a function of differential settlement * minimum pile spacing (i.e., automatically determined).
        
    real, intent(in) :: femvech(:),femvecv(:) !size of FEM elements in the horizontal and vertical directions
    
    

    real icost(ninv) !average costs associate with investigation
    real, intent(out) :: fcost(:), pcost(:), probfail(:), avediff(:), diffgeo(:) !failure cost, pile cost, prob. failure, ave. diff. set, geometric statistic
    real :: diffgeo2(ninv)


	real :: costs(nrep_MC,ninv) !costs for each investigation and realisation
	real :: diffset(nrep_MC,ninv) !max differential settlement for each investigation and realisation
    
    logical, intent(in) :: usepie !true to use FEM approximation such as the PIE method or M&K method, otherwise use FEM directly
    
    real CKheight(nrep_MC,nlayer-1,preps(1)*preps(2)) !effective layer depths at each pile associated with the full, original soil
    
    !effective pile layer stuff
    real, intent(in) :: sdist(:,:,:),sumweights(:) !distance of elements from each pile, sum of distance weights
    integer, intent(in) :: extents(:,:,:)!extent of soil around pile
    
    !other variables
    
    real, intent(out) :: si_performance(:) !site investigation performances
    character(1000) str2
    integer, intent(in) :: costmetric !Performance metric to use. -1 = failure cost, 0 = probability of failure, 1 or more = average differntial settlement to the power of the costmetric value (higher values more heavily penalise excessive values)
    integer, intent(in) :: EA_generation !current generation of the EA algorithm
    integer, intent(in) :: finaloutput !output mode for EA. 1 = final only, 2 = save best of each EA generation, 3 = save full final population, 4 = save everything
    
    integer pd,iter,x,y,jr,counter,i,j,pile,p2 !loop counters
    real randu !random generator
    logical, parameter :: debug =.false.
    character(2), parameter :: rednames(6) = (/ 'SA','GA','HA','1Q','SD','MN' /) !reduction method names; hard-coded in SI module
    integer, intent(in) :: deterministic !SI mode. 1 = 'deterministic' analysis of investigation performance, 2 = 'heat map', 3 = optimise investigations with EA
    integer tempEAit !temporary EA iteration count
    
	!want a higher res curve to get the inverse of the settlement function. 4 should be enough.
	integer, parameter :: curvefactor=10
	real detsetsi((npdepths-1)*curvefactor+1)
	real pdepthssi((npdepths-1)*curvefactor+1)
	
    
    !time variables
    real :: start, finish,allstart,allfinish
    
    
    !multiple layer site investigation variables
    
    !integer, parameter :: nlayer=2
    !real bfld(nlayer+1,nxew,nyew)
    !real        :: goodvals(ninv,nlayer,nbh,nze)                   !store valid values of SI samples in 1D array
    !integer	 :: nsamples(ninv,nlayer,nbh)						!keep track of how many viable samples are in each layer in each borehole
    !integer	 :: bound_depths(ninv,nlayer+1,nbh)			!layer boundary information as encountered by each borehole and investigation    bound_depths(n_inv,nlayer+1,n_inv,n_bh)
    integer, parameter :: npl = 2,invpower=2,mindist=1
    real, parameter ::  power=2
    real, allocatable :: layer_at_piles(:,:,:,:)
    real :: evals(nrep_MC,ninv,nlayer) !evals(ninv,size(inv_reduction),nlayer)
    real, allocatable :: si_pset(:,:,:,:)    !pile settlement curves for piles from soil model
    real, allocatable :: si_tempset(:,:,:) !subset of pile curves with increased curve resolution
    
    
    real testfield(160,160),ofield(160,160)
    
    logical ave_horiz_mode !
    logical, parameter :: use_avedepths = .false. !true if a horizontal boundary for each layer is used, where the depth is the average of all boreholes
    

    call cpu_time(allstart)
    allones = all(inv_bh == 1)
    
    ave_horiz_mode = use_avedepths .or. allones
    
    
    !allones = .false.
    
    if (ave_horiz_mode .and. usepie) then !If there is only a single borehole, the layers for each pile are constant. Only store once.
    	allocate(layer_at_piles(nrep_MC,ninv,nlayer-1,1))
    	allocate(sides(size(rel_loads),nrep_mc,ninv)) 
        goodcases_design = size(rel_loads)
    else
    	allocate(layer_at_piles(nrep_MC,ninv,nlayer-1,npl2)) !otherwise store for every pile
    	allocate(sides(preps(1)*preps(2),nrep_mc,ninv))
        goodcases_design = goodcases
    end if
    
    !get indices for interpolation
    x = 0
    do i=1,nyew
		do j = 1,nxew
			x = x + 1
			xyi(1,x) = j
			xyi(2,x) = i
		end do
    end do
    
        counter = 0
        do x = 1,npl2
            do i = 1,prad(1)
                do j = 1,prad(2)
                    counter =  counter +1
                    plocation_dble(1,counter) = plocation(goodpiles(x),1) + i - 1 !plocation_dble(2,:)
                    plocation_dble(2,counter) = plocation(goodpiles(x),2) + j - 1
                end do
            end do
        end do

    !############################# DO MONTE CARLO ANALYSIS SITE INVESTIGATIONS ###############################

    !inv_coords(:,:,2)

    !if (save_effE == 2 .and. deterministic == 1) then !Read in site investigation results (full pile settlement curves from soil model) from file. Only makes sense in deterministic analysis.
    !        allocate(si_pset(npdepths,nrep_MC,ninv,preps(1)*preps(2)))
    !        write(str2,'(A,A,I0,A,I0,A,I0,A,I0,A)') trim(datafolder),'si-piles_ninv-',ninv,'_nlayers-',nlayer,'_ratio2l-',nint(lmean_ave(2)/lmean_ave(1)),'_depth2l-',nint(dz*ldepths(1)),'.dat'
    !        open(500,file=str2)
    !        read(500) si_pset
    !        close(500)
    !    
    !else
    
        !Do multiple layer site investigations
 
        !get timing information
		!call cpu_time(start)
        
        bfld(1,:,:) = 1
        bfld(nlayer+1,:,:) = nzew
        
        
        !pre-process triangulation
        if (nbh > 2) then
            do i = 1,ninv
                call interp_layers_prep(2,inv_bh(i),nxew,nyew,inv_coords(i,:inv_bh(i),:),node_xy_out(i,:inv_bh(i),:),zd_out(i,:inv_bh(i)),element_neighbor(i,:,:4*inv_bh(i)),triangle(i,:,:4*inv_bh(i)),element_num(i),extranum(i))   
            end do
        else
            extranum = 0
        end if
        


		!---loop through Monte Carlo realisations---
        !if specified, the CK soils have been pre-processed along with the CK layer depths at pile locations. No need to generate them again.
    if (superset) then
		do iter = 1,nrep_MC    
			!call cpu_time(start)
            !write(*,*) iter
            
            bfld(2:nlayer,:,:) = fullbfld(iter,:,:,:)
            
            
			!---Generate virtual soil---
			kseed = randu(soilseeds(iter)) * 1234567890
            !call soil_layers()                             !get layer boundaries
            !efld = lmeans(nlayer)                     !Fill with oldest, deepest soil
            !do j=nlayer-1,1,-1  !Progressively work forwards through time, adding newer layer
            !    do x=1,nxew     !fill in x direction
            !        !write(*,'(1000000000I4)') nint(bfld(n+1,x,:))
            !        do y=1,nyew !fill in y direction
            !            efld(x,y,:nint(bfld(j+1,x,y))) = lmeans(j)                              
            !        end do
            !    end do
            !end do
            
            !Get layer heights at each pile for true, full, original soil (later used to determine actual differential settlement)
            !This only needs to be done once per MC simulation, so theoretically it could be pre-processed earlier in the program, but it's pretty quick so it might as well be here.
            !call getheights(bfld(2:nlayer,:,:),indices,plocation(goodpiles(:goodcases),:),nxew,nyew,goodcases,prad,invpower,radius,nlayer,tempheight) !
            !CKheight(iter,:,goodpiles(:goodcases)) = tempheight
            
!             do j = 1,goodcases
!                 do i = 1,nlayer-1 !get the effective heights for each layer.        
! 	                CKheight(iter,i,goodpiles(j)) = dz * sum(bfld(i+1,extents(j,1,1):extents(j,1,2),extents(j,2,1):extents(j,2,2)) * sdist(j,extents(j,1,1):extents(j,1,2),extents(j,2,1):extents(j,2,2))) / sumweights(j)
!                 end do
!             end do

            !bfld(2,:,:) = transpose(bfld(2,:,:))
            
            !conduct site investigations
            call si_multi_stuff(nxew,nyew,nzew,lmeans(:,iter),bfld,nlayer,ninv,nbh,inv_bh,inv_coords,inv_depths,size(test_errors,1),inv_test,add_errors,use_CI,conf_int,test_errors,size(test_errors,2),kseed, &
						plocation(goodpiles(:goodcases),:),prad,npl2,power,mindist,percentile(1),s_dev(1),size(inv_reduction),inv_reduction,invpower, &
						evals(iter,:,:),layer_at_piles(iter,:,:,:),sdist,sumweights,extents, &
                        node_xy_out,zd_out,element_neighbor,triangle,element_num,extranum,allones,use_avedepths,xyi,multitype,str2,plocation_dble,iter,finaloutput,soil_reps,rand_realisations,lmean_ln,lsd_ln) 
            
            !if (iter < 4) then
            !    write(*,*) ckheight(iter,1,goodpiles(:goodcases))
            !    write(*,*) layer_at_piles(iter,1,1,:)
            !    write(*,*) layer_at_piles(iter,94,1,:)
            !
            !    write(*,*)
            !end if
            !
            !bfld(2,:,:) = transpose(bfld(2,:,:))
            
            !write(*,*) 'done'
            !read(*,*)
                
			!call cpu_time(finish)
			!write(*,*) iter,(finish - start) !<- Don't write the progress to screen if using the EA
        end do
    else
        do iter = 1,nrep_MC    
			!call cpu_time(start)
            !write(*,*) iter

			!---Generate virtual soil---
			kseed = randu(soilseeds(iter)) * 1234567890
            
            call soil_layers(xyi) !get layer boundary
            
 
            
            !Get layer heights at each pile for true, full, original soil (later used to determine actual differential settlement)
            if(multitype == 1) then !just use the layer boundaries within the pile
                do j = 1,goodcases
                    do i = 1,nlayer-1        
	                    CKheight(iter,i,goodpiles(j)) = dz * sum(bfld(i+1,plocation(goodpiles(j),1):plocation(goodpiles(j),1)+prad(1)-1,plocation(goodpiles(j),2):plocation(goodpiles(j),2)+prad(2)-1))/product(prad)
                    end do
                end do

            else !get the effective heights for each layer based on inverse distance weighted average of surrounding soil
                 do j = 1,goodcases
                     do i = 1,nlayer-1     
 	                    CKheight(iter,i,goodpiles(j)) = dz * sum(bfld(i+1,extents(j,1,1):extents(j,1,2),extents(j,2,1):extents(j,2,2)) * sdist(j,extents(j,1,1):extents(j,1,2),extents(j,2,1):extents(j,2,2))) / sumweights(j)
                     end do
                 end do
            end if

            
            !conduct site investigations
            call si_multi_stuff(nxew,nyew,nzew,lmeans(:,iter),bfld,nlayer,ninv,nbh,inv_bh,inv_coords,inv_depths,size(test_errors,1),inv_test,add_errors,use_CI,conf_int,test_errors,size(test_errors,2),kseed, &
						plocation(goodpiles(:goodcases),:),prad,npl2,power,mindist,percentile(1),s_dev(1),size(inv_reduction),inv_reduction,invpower, &
						evals(iter,:,:),layer_at_piles(iter,:,:,:),sdist,sumweights,extents, &
                        node_xy_out,zd_out,element_neighbor,triangle,element_num,extranum,allones,use_avedepths,xyi,multitype,str2,plocation_dble,iter,finaloutput,soil_reps,rand_realisations,lmean_ln,lsd_ln) 

        end do

    end if
		!evals = evals*sdata(1) !scale up effective soil stiffness to actual values. Commented because now scaling weight instead of soil mean.
		!get and output timing information
		!call cpu_time(finish)
		!write(*,*) 'SI: ',(finish - start)
    !end if

    
            !############################# PROCESS PILE CONFIGURATION ###############################
            
            !expand the relative loads for the good piles to work with the design subroutine
            do i = 1,goodcases
                rel_loads2(i) = rel_loads(load_con(goodpiles(i)))
                load_con2(i) = i
            end do
            
    
            !calculate distances between each pile here instead of in the heavily nested loop below
	        counter = 0
	        do pile = 1,goodcases-1
		        do p2 = pile+1,goodcases
			        counter = counter + 1
			        dist(counter) = sqrt(real((plocation(goodpiles(pile),1)-plocation(goodpiles(p2),1))**2 + (plocation(goodpiles(pile),2)-plocation(goodpiles(p2),2))**2))
		        end do
            end do
            
    
            !get absolute settlement design tolerance based on minimum pile spacing and differential settlement limit, or directly use abstol value if a positive value is given.
            if (abstol > 0) then
                settol = abstol
            else
                settol = ceiling(minval(dist(:counter))*dz*difftol*1000) !round settlement tolerance up to nearest millimetre
            end if 

            
            !get CK loads and connectivity
	        do i=1,goodcases
		        rel_loadsCK(i) = rel_loads(load_con(goodpiles(i)))
                load_conCK(i) = i
            end do

	        
	        !build higher resolution settlement curve to help ensure that the inverse of the function matches it closely (boundary condition interpretation means this isn't neccessarily the case)
	        pdepthssi(1) = 0 !first point
	        counter = 1
	        do i=1,npdepths-1 !loop through pairs of original x points
				do j = 1,curvefactor !linearly interpolate new x points between pair
					counter = counter + 1
					pdepthssi(counter) = j*real(pdepths(i+1) - pdepths(i))/curvefactor + pdepths(i)	
				end do
            end do
	        !pdepthsssi(npdepths*curvefactor) = pdepths(npdepths)
            
	        call akima(real(pdepths), detdisp, pdepthssi, npdepths, 1, 1, 1, size(pdepthssi), detsetsi)

   
            !############################# GET RAW SITE INVESTIGATION PERFORMANCE INFORMATION ###############################
	
            !design foundations accordint to site investigation information
            if (usepie) then !apply Mylonakis and Gazetas settlement method to get pile design and differential settlement
            
                layer_at_piles = layer_at_piles * dz !convert layer depths to metres
                !call cpu_time(start)
                call multides_1D(prad,nrep_MC,ninv,goodcases,nlayer,layer_at_piles,evals,sides(:goodcases_design,:,:),rel_loads2(:goodcases),load_con2(:goodcases), settol, buildingweight,dz,nzew,plocation(goodpiles(:goodcases),:)*sngl(dz)*1000,lmeans,CKheight(:,:,goodpiles(:goodcases)),diffset,rel_loads,ave_horiz_mode,load_con(goodpiles(:goodcases)))
                                                       !multides_1D(pdepths,npdepths,prad,nrep_MC,ninv,npl,nlayer,layer_at_piles,evals,sides,rel_loads,load_con, destol, buildingweight,dz,nzew,plocation,sdata,CKheights,diffset)
                !call cpu_time(finish)
                !write(*,*) 'design', finish-start 
            else
                
	            if (save_effE == 1) then  !Calculate full pile settlement curves
                    if(.not. allocated(si_pset)) allocate(si_pset(npdepths,nrep_MC,ninv,preps(1)*preps(2)))
                    call multisets(pdepths,npdepths,femvech,femvecv,prad,nrep_MC,ninv,preps(1)*preps(2),layer_at_piles,evals,si_pset)
                    !call prep_fem2d(evals,0.3,femvech,femvecv,pradin,disp,pd,dz) !replace this with a wrapper for the investigations based on young's modulus and layer dept
                end if
                if (save_effE == 1 .or. save_effE == 2) then !design pile from settlement curves (curves have either been created fresh, or loaded in)
                    if(.not. allocated(si_pset)) allocate(si_pset(npdepths,nrep_MC,ninv,preps(1)*preps(2)))
                    allocate(si_tempset((npdepths-1)*curvefactor+1,nrep_MC,preps(1)*preps(2)))
                    !tweak this to work for curves, not young's modulus
                    do i = 1,ninv
                        call akima(real(pdepths), si_pset(:,:,i,:), pdepthssi, npdepths, 1, nrep_MC*preps(1)*preps(2), 1, size(pdepthssi), si_tempset)
                        call despile(si_tempset(:,:,goodpiles(:goodcases)), real(pdepthssi*dz), rel_loads2(:goodcases), load_con2(:goodcases), settol, buildingweight, detsetsi,sides(:goodcases,:,i), nrep_MC, goodcases,size(pdepthssi),1,size(pdepthssi),goodcases,goodcases)                                                                                                                                                                 
                    end do
                    deallocate(si_tempset)
                else !design pile iteratively using binary-search (bisection method) to get the design. This is usually faster than generating the full settlement curve, and is MUCH faster for long piles.
                    call multides(pdepths,npdepths,femvech,femvecv,prad,nrep_MC,ninv,goodcases,layer_at_piles,evals,sides,rel_loads2(:goodcases),load_con2(:goodcases), settol, buildingweight, ave_horiz_mode)
                end if
                
                !calculate differential settlement of the above designs based on the true pile performance information
	            call getdiff_cts(ck_set(:,:,goodpiles(:goodcases)), sides, real(pdepthssi*dz), rel_loads2(:goodcases), load_con2(:goodcases), real(plocation(goodpiles(:goodcases),:)*dz*1000),buildingweight,diffset, nrep_MC, goodcases,size(pdepthssi),ninv,goodcases)		!get get true max differential settlement
            
            end if
            
        
            
            !Interpolate failure costs from differential settlement
            call getcosts(diffset,costs,failurevals(1),failurevals(2),costvals(1),costvals(2),nrep_mc,ninv) 
            
    
	        !Get true pile designs
            !This subroutine appears to be causing an error, so it's left out. Fortunately, it's not needed in the analysis.
	        !call despile(ck_set(npdepths:1:-1,:,goodpiles(:goodcases)), real(pdepths(npdepths:1:-1)*dz), rel_loadsCK(:goodcases), load_conCK(:goodcases), settol, buildingweight, detdisp,ckdes(:goodcases,:), nrep_MC, goodcases,npdepths,1,npdepths,goodcases,goodcases)	!design foundations from site investigations
    
            !############################# CALCULATE OVERALL PERFORMANCE OF EACH INVESTIGATION ###############################
    
            !get locations of good information
            goodinv = costs >= 0
            allocate(goodloc(nrep_MC,ninv))
            invcount = count(goodinv,1)
            do j=1,ninv
                goodloc(:invcount(j),j) = pack([(i,i=1,nrep_MC)],goodinv(:,j))
            end do
    
                
                !open(505,file='sides_spx.dat',access='stream')
                !write(505) sides(1,:,:)
                !close(505)
                !

            !get total length of each pile case based on number of piles for each case 
            allocate(avesides(nrep_mc,ninv))
                if(ave_horiz_mode) then
                num_loadSI = 0
                do i = 1,num_loads
                    if(count(load_con==i) > 0) then
                        num_loadSI = num_loadSI + 1
                        sides(i,:,:) = sides(i,:,:) * count(load_con==i)
                        load_conSI(num_loadSI) = i
                    end if
                end do
                avesides = sum(sides(load_conSI(:num_loadSI),:,:),1)
            else
                num_loadSI = npl2
                avesides = sum(sides(:num_loadSI,:,:),1)   
            end if
            
            !avesides(:invcount(i),:) = sum(sides(load_conSI(:num_loadSI),goodloc(:invcount(i),i),:),1) !save total cumulative pile length in first index
            
            
                      
            !add average costs together for each investigation
            
			!--------Get a variety of performance metrics for each investigation, typically some form of average relating to differential settlement --------- 
		
			call calc_perf(ninv,nrep_MC,costmetric,pilecost,testcost,inv_test,inv_depths,inv_bh,invcount,goodloc,failurevals,avesides,costs,diffset,fcost,pcost,icost,probfail,avediff,diffgeo,diffgeo2,si_performance)
            




    ! ---------------- output results to files --------------------
            
                !print average and standard deviation of differential settlement and pile lengths for current input
    !for predominantly debugging purposes
    if (finaloutput == 7) then
        open(505,file='diffset.txt')
        write(505,'(A)') '! rows are investigations, columns are Monte Carlo realisations'
        do i = 1,ninv
             !tempvec(:invcount(i)) = diffset(goodloc(:invcount(i),i),i)
             !tempmean = sum(tempvec(:invcount(i)),1)/invcount(i)
             !tempvec(:invcount(i)) = tempvec(:invcount(i)) - tempmean
             !tempsd = sum(tempvec(:invcount(i))**2)
             !tempsd = sqrt(tempsd/invcount(i))
             
             !write(505,*) tempmean,tempsd
             !do j = 1,invcount(i)
             !    write(505,*) diffset(goodloc(j,i),i)
             !end do
            write(505,'(100000(E8.3,X))') diffset(:,i) 
        end do
        close(505)
        
        
        do j = 1, npl2
            write(str2,'(A,I0,A)') 'sides_pile-',j,'.txt'
            open(505,file=str2)
            write(505,'(A)') '! rows are investigations, columns are Monte Carlo realisations'
            do i = 1,ninv
            
                 !tempvec(:invcount(i)) = avesides(goodloc(:invcount(i),i),i)
                 !tempmean = sum(tempvec(:invcount(i)),1)/invcount(i)
                 !tempvec(:invcount(i)) = tempvec(:invcount(i)) - tempmean
                 !tempsd = sum(tempvec(:invcount(i))**2)
                 !tempsd = sqrt(tempsd/invcount(i))
                 !write(505,*) tempmean,tempsd
                 !do j = 1,invcount(i)
                 !    write(505,*) avesides(goodloc(j,i),i)
                 !end do
                write(505,'(100000(E8.3,X))') sides(j,:,i) 
            end do
            close(505)
        end do
        


    end if
    
        
       ! do i = 1,160
       !     counter = 0
       !     do j=(i-1)*160+1,(i-1)*160+160
       !         counter = counter + 1
			    !testfield(i,counter) = layer_at_piles(1,j,1,1)*2
       !     end do
       ! end do
        
  
        
        !ofield = fullbfld(1,1,:,:)
        
                !do j = 1,goodcases   
	               !     write(*,*) CKheight(1,1,goodpiles(j)),sum(ofield(plocation(goodpiles(j),1):plocation(goodpiles(j),1)+prad(1)-1,plocation(goodpiles(j),2):plocation(goodpiles(j),2)+prad(2)-1))*0.5/product(prad),sum(testfield(plocation(goodpiles(j),1):plocation(goodpiles(j),1)+prad(1)-1,plocation(goodpiles(j),2):plocation(goodpiles(j),2)+prad(2)-1))*0.5/product(prad)
                !end do
                !write(*,*)
                !do j = 1,goodcases   
	               !     write(*,*) CKheight(1,1,goodpiles(j))*2,sum(ofield(plocation(goodpiles(j),1):plocation(goodpiles(j),1)+prad(1)-1,plocation(goodpiles(j),2):plocation(goodpiles(j),2)+prad(2)-1))/product(prad),sum(testfield(plocation(goodpiles(j),1):plocation(goodpiles(j),1)+prad(1)-1,plocation(goodpiles(j),2):plocation(goodpiles(j),2)+prad(2)-1))/product(prad)
                !end do
                !
    
        !open(505,file='lap_spx.dat',access='stream')
        !write(505) layer_at_piles(:,:,1,1)          ! layer_at_piles(1,:,1,1)
        !close(505)
    !!    
    !    open(505,file='ckheight_spx.dat',access='stream')
    !    write(505) CKheight
    !    close(505)
    !!
    !
    !!
    !!!sides(size(rel_loads),nrep_mc,ninv)
    !!
    !open(505,file='diffset_3L.dat',access='stream')
    !write(505) diffset
    !close(505)
    !

    
    

    
    write(soildsc,'(A,I0,A,I0,A,I0,A,I0,A,I0)') '_nlayers-',nlayer,'_ratio2l-',nint(lmean_ave(1)/minval(lmean_ave(:2))),'+',nint(lmean_ave(2)/minval(lmean_ave(:2))),'_depth2l-',nint(dz*ldepths(1)),'_bSD-',nint(bsd*dz)
    call save_output(soildsc,finaloutput,deterministic,goodcases,goodloc,invcount,nrep_MC,buildingarea,ninv,nbh,in_tests,in_depths,in_reductions,dz,inv_coords,inv_bh,inv_depths,inv_reduction,inv_test,swidth,soffset,sstep,testnames,rednames, &
			EA_generation,costmetric,plocation,avesides,diffset,si_performance,costs,fcost,pcost,icost,probfail,avediff,diffgeo,diffgeo2)
    deallocate(avesides)
    
    !if it's running the GA, make all investigations that have more than half invalid realisations have a bad score
    if(deterministic == 3) then
        where(invcount < ceiling(real(nrep_MC)/5)) si_performance = huge(si_performance(1))
    else
    
    !This should never really happen, so I'm treating it as a critical error.
            if(any(invcount == 0)) then
                write(*,*) 'Error: no valid realisations.'
                write(*,*) 'Try alternate soil strength, or'
                write(*,*) 'increase realisation count.'
                write(*,*) "Press 'enter' to continue."
                read(*,*)
                stop
            end if
            
    end if
            
    call cpu_time(allfinish)
	!write(*,*) EA_generation,(allfinish - allstart)
    !do i = 1, ninv
    !
    !    write(*,*) 'max, min, mean, sd', maxval(diffset(goodloc(:invcount(i),i),i)),minval(diffset(goodloc(:invcount(i),i),i)),avediff(i),sqrt(sum((diffset(goodloc(:invcount(i),i),i) - avediff(i))**2)/invcount(i))
    !end do 
    !read(*,*)
    !read(*,*)
    
        end subroutine
        
        
        !get pile settlement curve using 2D axisymmetric FEM for site investigation (deprecated)
        subroutine multisets(pdepths,npdepths,femvech,femvecv,prad,nrep_MC,ninv,npl,layer_at_piles,evals,si_pset)
        
            real efld2d(nxew,nzew) !Virtual soil 2D representation; Uniform properties within each layer, horizontal layer boundary
            integer, intent(in) :: npdepths							!number of pile depths to analyse
	        integer, intent(in) :: pdepths(:)			!pile depths to analyse (elements)
            real, intent(in) :: femvech(:),femvecv(:) !size of FEM elements in the horizontal and vertical directions
            integer, intent(in) :: prad(2)			!pile width in x,y directions
            integer i,counter,pile,iter,inv,pd                   !loop counters
            real(4) disp                        !pile settlement
            integer,intent(in) ::     nrep_MC,ninv,npl
            
            
            real, intent(in) :: layer_at_piles(nrep_MC,ninv,nlayer-1,npl)
            integer :: lap(ninv,nlayer-1,npl) !integer version of the above array
            real, intent(in) :: evals(nrep_MC,ninv,nlayer) !evals(ninv,size(inv_reduction),nlayer)
            
            real, intent(out) :: si_pset(npdepths,nrep_MC,ninv,npl)
        
            do iter = 1,nrep_MC                                     !loop through realisations
                
                !Loop counter for 10% increments
                if((iter-int(10*iter/nrep_MC)*nrep_MC/10 == 0)) write(*,'(I0,A,X)',advance='no') iter*100/nrep_MC,'%'
                
                !Finally round the heights to the nearest element, making sure the resulting values are within the field
                lap = nint(layer_at_piles(iter,:,:,:))
                where(lap < 1) lap = 1
                where(lap > nzew) lap = nzew
                
                do inv = 1,ninv                                     !loop through investigations
                    do pile = 1,npl                                 !loop through piles 
        
                        !construct 2D soil profile
                        efld2d = evals(iter,inv,nlayer)                     !Fill with oldest, deepest soil
                        do i=nlayer-1,1,-1
                            efld2d(:,:lap(inv,i,pile)) = evals(iter,inv,i)      !Progressively work forwards through time, adding newer layers
                        end do
                        !lap(inv,:,pile)
        
                        !Construct pile settlement curve 
                        counter = 0
                        do pd = 1,npdepths
                            counter = counter + 1
					        call prep_fem2d(efld2d,0.3,femvech,femvecv,prad(1),disp,pdepths(pd),dz)
                            si_pset(counter,iter,inv,pile) = disp
                        end do
                    end do
                end do
            end do
            
        
        end subroutine
            
        
             !Design pile using using 2D axisymmetric FEM for site investigation (deprecated)
            !Find the design using binary search (i.e. bisection method) - a root finding method that's generally faster than a linear search, especially with long piles.
            !It's an iterative bracketing method that halves the search space each iteration until consequtive designs above and below the tolerance is found (then the bigger one is kept).
        subroutine multides(pdepths,npdepths,femvech,femvecv,prad,nrep_MC,ninv,npl,layer_at_piles,evals,sides,rel_loads,load_con, destol, buildingweight,allones)
        
            real efld2d(nxew,nzew) !Virtual soil 2D representation; Uniform properties within each layer, horizontal layer boundary
            integer, intent(in) :: npdepths							!number of pile depths to analyse
	        integer, intent(in) :: pdepths(:)			!pile depths to analyse (elements)
            real, intent(in) :: femvech(:),femvecv(:) !size of FEM elements in the horizontal and vertical directions
            integer, intent(in) :: prad(2)			!pile width in x,y directions
            integer i,counter,pile,iter,inv,pd                   !loop counters
            real(4) disp                        !pile settlement
            integer,intent(in) ::     nrep_MC,ninv,npl
            real, intent(in) :: rel_loads(:)		!set of relative pile loads based on tributary area of building
	        integer, intent(in) :: load_con(:)			!connectivity vector saying which pile has which load
            real, intent(in) :: destol !design tolerance (mm)
            real :: temp_destol !scaled design tolerance 
            real, intent(in) :: buildingweight
            real :: sum_loads !sum of pile loads
            
            
            real, intent(in) :: layer_at_piles(nrep_MC,ninv,nlayer-1,npl) !layer depth information
            integer :: lap(ninv,nlayer-1,npl) !integer version of the above array
            real, intent(in) :: evals(nrep_MC,ninv,nlayer) !evals(ninv,size(inv_reduction),nlayer)
            
            
            real, intent(out) :: sides(npl,nrep_mc,ninv)   !SI pile designs
            
            logical, intent(in) :: allones !whether all investigations have a single borehole
            real start,finish
            
            real nan                        !nan
            real leftdisp,rightdisp,cdisp   !displacements associated with the left bracket, right bracket and centre
            integer left,right,mid          !pile lengths (elements) associated with the left/right brackets and centre
            
            integer npltemp
            
            nan = 0.0 !generate undefined value
  	        nan = 0.0/nan
            
            if(allones) then
                npltemp = 1
            else
                npltemp = npl
            end if
            
            
            !get sum of pile loads
	        sum_loads = 0.0
	        do pile = 1,npl
                if(load_con(pile) == 0) cycle  !skip piles with no associated load
		        sum_loads = sum_loads + rel_loads(load_con(pile))
            end do
            
            
        
            do iter = 1,nrep_MC                                     !loop through realisations
            !do concurrent (iter=1:nrep_MC)
      
                !Loop counter for 10% increments
                if((iter-int(10*iter/nrep_MC)*nrep_MC/10 == 0)) write(*,'(I0,A,X)',advance='no') iter*100/nrep_MC,'%'
                
                !Finally round the heights to the nearest element, making sure the resulting values are within the field
                lap = nint(layer_at_piles(iter,:,:,:))

                
                call cpu_time(start)
                
                invloop: do inv = 1,ninv                                     !loop through investigations
                    
                    do pile = 1,npltemp !                    !loop through piles 
                        
                        !scale the design tolerance (mm) according to the applied load
                        temp_destol = destol * sum_loads/ (buildingweight * rel_loads(pile))
        
                        !construct 2D soil profile
                        efld2d = evals(iter,inv,nlayer)                     !Fill with oldest, deepest soil
                        do i=nlayer-1,1,-1
                            efld2d(:,:lap(inv,i,pile)) = evals(iter,inv,i)      !Progressively work forwards through time, adding newer layers
                        end do

                    ! ---- perform binary search! ------
	
		                !Check if it's in the outer bounds
                        left = pdepths(1)
                        right = pdepths(npdepths)
                        call prep_fem2d(efld2d,0.3,femvech,femvecv,prad(1),leftdisp,left,dz)
                        call prep_fem2d(efld2d,0.3,femvech,femvecv,prad(1),rightdisp,right,dz)          
                        if (leftdisp < temp_destol) then 
                            sides(:,iter,inv) = -100 
                            cycle invloop
                        else if(rightdisp > temp_destol) then
                            sides(:,iter,inv) = -200
 			                cycle invloop
                        end if
	  
		                !Otherwise, use binary search to find interval containing true design

		                do while (right > left + 1) !exit after the brackets are consequtive elements.

			                    mid = ( left + right + 1 ) / 2
                                call prep_fem2d(efld2d,0.3,femvech,femvecv,prad(1),cdisp,mid,dz)

			                    if ( cdisp >= temp_destol ) then
				                    left = mid
			                    else
				                    right = mid
			                    end if
	
                        end do
 
                        !binary search would have converged by now assuming a valid design exists.
                        sides(pile,iter,inv) = right !make sure it uses the LARGER pile length (greater than the design tol)

                    

                    end do
                end do invloop
                
                call cpu_time(finish)
                write(*,*) iter, finish-start
                
            end do
            
            if(allones) then
                do pile = 2,npl
                    sides(pile,:,:) = sides(1,:,:)
                end do
            end if

            

        end subroutine
        
        

    end module
