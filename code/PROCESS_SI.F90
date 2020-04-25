
   
module process_SI

use SI
use variables
use getdiff
use setup_SI
use si_stats
use output_results

implicit none

contains
 
!Get site investigation performance information, done by the following steps:
!   1. Conduct site investigations across all realisations of random soils
!   2. Design piles according to those investigations, 
!   3. Get true performance in terms of differential settlement (requires program be run in 'CK' mode to get true pile behaviour first).
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
    
    subroutine get_si_perf(soilseeds,ninv,nbh,in_tests,in_depths,in_reductions,inv_bh,inv_coords,inv_depths,inv_test,add_errors,use_CI,test_errors,inv_reduction,conf_int,percentile,s_dev,soffset,swidth,sstep &					!site investigation variables
                      ,npdepths,detdisp,plocation,pdepths,ck_set,nrep_mc,rel_loads,load_con,num_loads,preps,buildingweight,difftol,failurevals,costvals,si_performance,pilecost,testcost,testnames,costmetric,EA_generation &
                      ,deterministic,finaloutput,abstol,datafolder,invcount,goodpiles,goodcases,fcost, pcost, probfail, avediff, diffgeo)
                      
                      
! ----- soil generation variables (don't touch) ----

      integer iseed
      integer, intent(in) :: soilseeds(:)
      real, allocatable :: efldave(:,:,:)
      integer xpos,ypos,zpos !offsets of random subset for the superset mode 
      character(200) :: soildsc !soil description string containing soil attributes, to help distinguish between different cases, particularly single layer vs multiple layers


! ---- site investigation variables -----

	integer,  intent(in) :: inv_bh(:)						!number of boreholes in each investigation
	integer,  intent(in)	:: inv_coords(:,:,:)			!x,y coordinates of boreholes for each investigation
	integer,  intent(in) :: inv_depths(:,:,:) 				!sampling depth start, end, interval
	integer,  intent(in)	:: inv_test(:,:)					!test type for each investigation
	logical,  intent(in)	:: add_errors						!whether to add test errors or not
	real,  intent(in) :: test_errors(:,:)	!matrix of test error stats (num tests x num errors)
	integer,  intent(in) :: inv_reduction(:)				!reduction method for each test
	real :: evals(nrep_mc,ninv)						!reduced value for each investigation
	real,  intent(in) :: percentile(:),s_dev(:) !percentile and standard deviation below mean for reduction methods
    real, intent(in) :: testcost(:)     !cost per metre for each test
    character(4), intent(in) :: testnames(:) !names for each test
	
	integer :: nbh,bhdepth,ntest	!No. boreholes, borehole depth (elements), No. tests
	integer :: ninv				!number of investigations
	integer,intent(in) :: in_tests,in_depths,in_reductions !number of tests, depths, and reductions investigated
	
	real :: conf_int !confidence interval for truncating SI sample values 
	logical :: use_ci 
    
    logical :: goodinv(nrep_MC,ninv) !array for valid site investigation performance information
    integer :: goodloc(nrep_MC,ninv) !array for valid site investigation performance locations
    integer :: invcount(ninv) !number of valid values for goodinv

    

    
    integer, intent(in) :: soffset(2)		!x,y offset of borehole from corner of soil (elements)
	integer, intent(in) :: swidth(2)		!x,y dimensions of site investigation area
    integer, intent(in) :: sstep(2)            !borehole step size in each dimension for the heat map mode
    

! --- pile performance variables
	integer, intent(in) :: npdepths							!number of pile depths to analyse
	integer, intent(in) :: pdepths(:)			!pile depths to analyse (elements)
    
	real, intent(in) :: detdisp(:)				!deterministic settlement for each pile depth (1 kN applied load. Soil stiffness 1 MPa)
	integer, intent(in) :: preps(2)		!No. piles in x,y directions - repetitions
	integer, intent(in) :: plocation(:,:)		!Pile x,y coordinates
	
	real, intent(in) :: rel_loads(:)		!set of relative pile loads based on tributary area of building
	integer, intent(in) :: load_con(:)			!connectivity vector saying which pile has which load
    integer,intent(in) :: goodpiles(:)        !vector of piles that have an associated applied load
    integer,intent(in) :: goodcases                               !number of good piles in goodpiles 
	
	real, intent(in) :: buildingweight,difftol
    real :: settol
    real :: dist((preps(1)*preps(2))**2) !distances between piles
	
	!real :: sides(num_loads,nrep_mc,ninv)   !SI pile designs   stored in terms of m
    integer :: sides(nrep_mc,ninv,num_loads)   !SI pile designs stored a an array index
    real :: avesides(nrep_mc,ninv)   ! average SI pile designs
    real :: ckdes(preps(1)*preps(2),nrep_mc)   !CK pile designs
	real, intent(in) :: ck_set(:,:,:) !true ck pile settlement for each realisation, pile and depth
    
    integer, intent(in) :: nrep_mc !number of MC realisations
    integer, intent(in) :: num_loads !number of load cases
    
    integer num_loadSI
    integer load_conSI(num_loads)
    
    real rel_loadsCK(preps(1)*preps(2)) !relative loads for CK piles
    integer load_conCK(preps(1)*preps(2))
    
    real, intent(in) :: failurevals(2), costvals(2) !2 x,y coordinates defining the line to interpolate costs from differential settlement
    real, intent(in) :: pilecost !cost for pile, per metre depth
    real, intent(in) :: abstol !absolute settlement tolerance  (mm). Positive values are used directly. Negative values make it a function of differential settlement * minimum pile spacing (i.e., automatically determined).
        

    real icost(ninv) !average costs associate with investigation
    real, intent(out) :: fcost(:), pcost(:), probfail(:), avediff(:), diffgeo(:) !failure cost, pile cost, prob. failure, ave. diff. set, geometric statistic
    !real avediff2(ninv) !average differential settlement squared
    real :: diffgeo2(ninv)
    !real diffgeo0(ninv)

	real :: costs(nrep_MC,ninv) !costs for each investigation and realisation
	real :: diffset(nrep_MC,ninv) !max differential settlement for each investigation and realisation
    !real :: evalsmean(ninv)
    
    !other variables
    
    real, intent(out) :: si_performance(:) !site investigation performances
    character(1000) str2
    integer, intent(in) :: costmetric !Performance metric to use. -1 = failure cost, 0 = probability of failure, 1 or more = average differntial settlement to the power of the costmetric value (higher values more heavily penalise excessive values)
    integer, intent(in) :: EA_generation !current generation of the EA algorithm
    integer, intent(in) :: finaloutput !output mode for EA. 1 = final only, 2 = save best of each EA generation, 3 = save full final population, 4 = save everything
    character(1000),intent(in) :: datafolder !a string representing the directory the data is stored in
    
    integer pd,iter,x,y,jr,counter,i,j,pile,p2 !loop counters
    real randu !random generator
    logical, parameter :: debug =.false.
    character(2), parameter :: rednames(6) = (/ 'SA','GA','HA','1Q','SD','MN' /) !reduction method names; hard-coded in SI module
    integer, intent(in) :: deterministic !SI mode. 1 = 'deterministic' analysis of investigation performance, 2 = 'heat map', 3 = optimise investigations with EA
    integer tempEAit !temporary EA iteration count
    
	!want a higher res curve to get the inverse of the settlement function. 4 should be enough.
	integer, parameter :: curvefactor=10
	real detsetsi((npdepths-1)*curvefactor+1)
	integer pdepthssi((npdepths-1)*curvefactor+1)
    !real pdepthssi((npdepths-1)*curvefactor+1)
    integer :: pdepths_ind(npdepths)
	
    
    !time variables
    real(8) :: start, finish,allstart,allfinish,rate
    integer cm,cr,cend,c1,c2
    real s1,f1 !timing variables


    !############################# DO MONTE CARLO ANALYSIS SITE INVESTIGATIONS ###############################

    !inv_coords(:,:,2)
    
    call cpu_time(s1)
    
    !if (save_effE == 2 .and. deterministic == 1) then !Read in site investigation results (effective young's modulus) from file. Only makes sense in deterministic analysis.
    !    write(str2,'(A,I0,A,I0,A,I0,A)') 'si_results/effective_E_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy,'.txt'
    !    open(507,file=str2)
    !    do i = 1,nrep_mc
    !        read(507,'(I6,X,100000F8.3)')j,evals(i,:)
    !    end do
    !    close(507)
    !    
    !else

        if(superset) then
            
            
            call cpu_time(s1)
            do iter = 1,nrep_MC    
                
                write(*,*) iter
                		
               !generate random offsets within superset soil
               kseed = randu(soilseeds(iter)) * 1234567890 !ensure random numbers are consistent across MC realisations
		       xpos=NINT(randu(0)*(nxe-nxew))
		       ypos=NINT(randu(0)*(nye-nyew))
		       zpos=NINT(randu(0)*(nze-nzew))

			    !--conduct site investigations--
                
                call s_inv(nxew,nyew,nzew,efld(xpos+1:xpos+nxew,ypos+1:ypos+nyew,zpos+1:zpos+nzew),ninv,nbh,inv_bh,inv_coords,inv_depths,size(test_errors,1),size(test_errors,2),inv_test,add_errors,use_CI,test_errors,inv_reduction,evals(iter,:),kseed,conf_int,percentile,s_dev)


            end do
            call cpu_time(f1)
            !write(*,*) 'si',f1-s1
            !read(*,*)
            !write(*,*) 'done'
            
        
	    !If getting site investigation performance once, generate virtual soil
	    else  !generate soil for each MC reaisation
	    
	        allocate(efldave(nxew,nyew,nzew))
			!open(500,file= trim(datafolder)//'soilave.dat',access='stream',status='old',action='read')
			!read(500) efldave
			!close(500)
	    
		    do iter = 1,nrep_MC    
                !write(*,*) iter
			    !CALL system_clock(c1)
			    ! --- generate soil ---
                kseed = randu(soilseeds(iter)) * 1234567890 !ensure random numbers are consistent across MC realisations
                ! You must guarentee that the soils generated here (efld) are IDENTICAL to those produced in the PROCESS_CK subroutine.
                
                
                call sim3dw()
							      
				!efld = efld - efldave !subtract average field to ensure the mean everywhere is zero
                
			    !conduct site investigations
			    call s_inv(nxew,nyew,nzew,efld,ninv,nbh,inv_bh,inv_coords,inv_depths,size(test_errors,1),size(test_errors,2),inv_test,add_errors,use_CI,test_errors,inv_reduction,evals(iter,:),kseed,conf_int,percentile,s_dev)
			    !CALL system_clock(c2)
			    !write(*,*) iter,dble(c2 - c1)/rate
            end do
        	deallocate(efldave)
        
        end if	
      !      !otherwise, read it in 
      !
      !      !get timing information
      !      !CALL system_clock(count_rate=cr)
      !      !CALL system_clock(cm)
      !      !rate = REAL(cr)
      !      
      !      open(507,file=trim(datafolder)//'soil.dat',access='stream',status='old',action='read')  !shared !share='DENYNONE' <- neither of these seem to work on gfortran (allow multiple processes to read same file simultaneously)
      !      
		    !do iter = 1,nrep_MC    
      !
			   ! !---read in virtual soil with unit mean stiffness---
      !          
			   ! read(507) efld
      !
			   ! !--conduct site investigations--
      !          
      !          kseed = randu(soilseeds(iter)) * 1234567890
      !          call s_inv(nxew,nyew,nzew,efld,ninv,nbh,inv_bh,inv_coords,inv_depths,size(test_errors,1),size(test_errors,2),inv_test,add_errors,use_CI,test_errors,inv_reduction,evals(iter,:),kseed,conf_int,percentile,s_dev)
      !
      !
      !      end do
      !      close(507)
		    !!evals = evals*sdata(1,1) !scale up effective soil stiffness to actual values. Commented because now scaling weight instead of soil mean.
		    !!get and output timing information
      !      !CALL system_clock(cend)
		    !!write(*,*) 'Site investigation:',real(cend-cm)/rate
      !  end if
    !end if

        
        !write(str2,'(A,I0,A,I0,A,I0,A)') 'si_results/effective_E_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy,'.txt'
        !open(507,file=str2)
        !do i = 1,nrep_mc
        !    write(507,'(100000(F9.4,X))'),evals(i,:)
        !end do
        !close(507)
        !stop
        !

            !call cpu_time(s1)
    
            !############################# PROCESS PILE CONFIGURATION ###############################
                
    
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
	        !pdepthssi(1) = 0 !first point for metres version
            pdepthssi(1) = 1
            pdepths_ind(1) = 1
	        counter = 1
	        do i=1,npdepths-1 !loop through pairs of original x points
                pdepths_ind(i+1) = i*curvefactor+1
				do j = 1,curvefactor !linearly interpolate new x points between pair
					counter = counter + 1
					!pdepthssi(counter) = j*real(pdepths(i+1) - pdepths(i))/curvefactor + pdepths(i)	 <- real version in terms of metres
                    pdepthssi(counter) = counter     !<- array index
				end do
            end do
	        !pdepthsssi(npdepths*curvefactor) = pdepths(npdepths)
            
	        call akima(real(pdepths_ind), detdisp, real(pdepthssi), npdepths, 1, 1, 1, size(pdepthssi), detsetsi)
            !call akima(real(pdepths), detdisp, pdepthssi, npdepths, 1, 1, 1, size(pdepthssi), detsetsi) <- metres version
            
    !    call cpu_time(f1)
    !write(*,*) 'other1',f1-s1
    !call cpu_time(s1)
   
            !############################# GET RAW SITE INVESTIGATION PERFORMANCE INFORMATION ###############################
	
            !design foundations accordint to site investigation information
	        !call despile(evals, real(pdepthssi*dz), rel_loads, load_con, settol, buildingweight, detsetsi,sides, nrep_MC, preps(1)*preps(2),size(pdepthssi),ninv,1,num_loads,1)	!<- metres version
            call despile_discrete(evals, rel_loads, load_con, settol, buildingweight, detsetsi,sides, nrep_MC, preps(1)*preps(2),size(pdepthssi),ninv,1,num_loads,1)
            
            !calculate differential settlement of the above designs based on the true pile performance information
	        call getdiff_discrete(ck_set(:,:,goodpiles(:goodcases)), sides, real(pdepthssi*dz), rel_loads, load_con(goodpiles(:goodcases)), real(plocation(goodpiles(:goodcases),:)*dz*1000),buildingweight,diffset, nrep_MC, goodcases,size(pdepthssi),ninv,num_loads)		!get get true max differential settlement

            !Interpolate failure costs from differential settlement
            call getcosts(diffset,costs,failurevals(1),failurevals(2),costvals(1),costvals(2),nrep_mc,ninv) 

    
	        !Get true pile designs
            !This subroutine appears to be causing an error, so it's left out. Fortunately, it's not needed in the analysis.
	        !call despile(ck_set(npdepths:1:-1,:,goodpiles(:goodcases)), real(pdepths(npdepths:1:-1)*dz), rel_loadsCK(:goodcases), load_conCK(:goodcases), settol, buildingweight, detdisp,ckdes(:goodcases,:), nrep_MC, goodcases,npdepths,1,npdepths,goodcases,goodcases)	!design foundations from site investigations
    
            !############################# CALCULATE OVERALL PERFORMANCE OF EACH INVESTIGATION ###############################
    
            !get locations of good information
            goodinv = costs >= 0
            invcount = count(goodinv,1)
            do j=1,ninv
                goodloc(:invcount(j),j) = pack([(i,i=1,nrep_MC)],goodinv(:,j))
            end do
    
    !    call cpu_time(f1)
    !write(*,*) 'other4',f1-s1
    !call cpu_time(s1)
        

            !get total length of each pile case based on number of piles for each case            
            num_loadSI = 0
            do i = 1,num_loads
                if(count(load_con==i) > 0) then
                    num_loadSI = num_loadSI + 1
                    sides(:,:,i) = sides(:,:,i) * count(load_con==i)
                    load_conSI(num_loadSI) = i
                end if
            end do
            
            !avesides(:invcount(i),:) = sum(sides(load_conSI(:num_loadSI),goodloc(:invcount(i),i),:),1) !save total cumulative pile length in first index
            avesides = sum(sides(:,:,load_conSI(:num_loadSI)),3)
            avesides = avesides / curvefactor
            
    !        call cpu_time(f1)
    !write(*,*) 'other5',f1-s1 sides(1,1,3)
    !call cpu_time(s1)
         
    !--------Get a variety of performance metrics for each investigation, typically some form of average relating to differential settlement --------- 

	
	call calc_perf(ninv,nrep_MC,costmetric,pilecost,testcost,inv_test,inv_depths,inv_bh,invcount,goodloc,failurevals,avesides,costs,diffset,fcost,pcost,icost,probfail,avediff,diffgeo,diffgeo2,si_performance)
    



	! ---------------- output results to files --------------------
	write(soildsc,'(A,I0,A,I0,A,I0)') '_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy
    call save_output(soildsc,finaloutput,deterministic,goodcases,goodloc,invcount,nrep_MC,buildingarea,ninv,nbh,in_tests,in_depths,in_reductions,dz,inv_coords,inv_bh,inv_depths,inv_reduction,inv_test,swidth,soffset,sstep,testnames,rednames, &
			EA_generation,costmetric,plocation,avesides,diffset,si_performance,costs,fcost,pcost,icost,probfail,avediff,diffgeo,diffgeo2)

    !print average and standard deviation of differential settlement and pile lengths for current input
    !for predominantly debugging purposes
    if (finaloutput == 5) then
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
        
        
        do j = 1, num_loadSI
            write(str2,'(A,I0,A)') 'sides_pileload-',j,'.txt'
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
                write(505,'(100000(E8.3,X))') real(sides(:,i,load_conSI(j)))/curvefactor
            end do
            close(505)
        end do
        
        

    end if
    
    
        !if it's running the GA, make all investigations that have more than half invalid realisations have a bad score
    if(deterministic == 3) then
        where(invcount < ceiling(real(nrep_MC)/3)) si_performance = huge(si_performance(1))
    end if
    
    
    !This should never really happen, so I'm treating it as a critical error.
            if(any(invcount == 0)) then
                write(*,*) 'Error: no valid realisations.'
                write(*,*) 'Try alternate soil strength, or'
                write(*,*) 'increase realisation count.'
                write(*,*) "Press 'enter' to continue."
                read(*,*)
                stop
            end if
            
    call cpu_time(allfinish)
    
    !call cpu_time(s1)
    !write(*,*) 'write',s1-f1
     !call cpu_time(f1)
    !write(*,*) 'si',f1-s1
    

    end subroutine

    end module
