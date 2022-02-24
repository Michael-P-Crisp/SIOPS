module output_results

!This module contains the majority of code for saving final investigation performance results and related information to text files



implicit none

contains


subroutine save_output(soildsc,finaloutput,deterministic,goodcases,goodloc,invcount,nrep_MC,buildingarea,ninv,nbh,in_tests,in_depths,in_reductions,dz,inv_coords,inv_bh,inv_depths,inv_reduction,inv_test,swidth,soffset,sstep,testnames,rednames, &
			EA_generation,costmetric,plocation,avesides,diffset,si_performance,costs,fcost,pcost,icost,probfail,avediff,diffgeo,diffgeo2)


!input
character(200),intent(in) :: soildsc !soil description string containing soil attributes, to help distinguish between different cases, particularly single layer vs multiple layers

real, intent(in) :: fcost(ninv),pcost(ninv),icost(ninv) !average costs associate with failure, pile construction, and investigation
real, intent(in) :: probfail(ninv) !probability of failure
real, intent(in) :: avediff(ninv) !average differential settlement
!real avediff2(ninv) !average differential settlement squared
real, intent(in) :: diffgeo(ninv) !Alternate metric for investigation quality. A geometric standard deviation above the geometric mean
real, intent(in) :: diffgeo2(ninv)
!real diffgeo0(ninv)
real, intent(out) :: si_performance(:) !site investigation performances

real, intent(in) :: buildingarea !building area
integer, intent (in) :: ninv, nbh, in_tests, in_depths, in_reductions !number of investigations, boreholes, tests, depths, reductions
integer,  intent(in) :: inv_bh(:)						!number of boreholes in each investigation
integer,  intent(in) :: inv_coords(:,:,:)			!x,y coordinates of boreholes for each investigation
integer,  intent(in) :: inv_depths(:,:,:) 				!sampling depth start, end, interval
integer,  intent(in) :: inv_test(:,:)					!test type for each investigation
integer,  intent(in) :: inv_reduction(:)				!reduction method for each test
character(4), intent(in) :: testnames(:) !names for each test
character(2), intent(in) :: rednames(:)  !reduction method names; hard-coded in SI module
integer, intent(in) :: nrep_MC !number of Monte Carlo realisations
real(8), intent(in) :: dz !element size

integer, intent(in) :: soffset(2)		!x,y offset of borehole from corner of soil (elements)
integer, intent(in) :: swidth(2)		!x,y dimensions of site investigation area
integer, intent(in) :: sstep(2)            !borehole step size in each dimension for the heat map mode

real, intent(in) :: costs(nrep_MC,ninv) !costs for each investigation and realisation
real, intent(in) :: diffset(nrep_MC,ninv) !max differential settlement for each investigation and realisation
integer, intent(in) :: deterministic !SI mode. 1 = 'deterministic' analysis of investigation performance, 2 = 'heat map', 3 = optimise investigations with EA

integer, intent(in) :: goodloc(nrep_MC,ninv) !array for valid site investigation performance locations
integer, intent(in) :: invcount(ninv) !number of valid values for goodinv

integer, intent(in) :: plocation(:,:)		!Pile x,y coordinates
real, intent(in) :: avesides(nrep_mc,ninv)   ! average SI pile designs

integer, intent(in) :: costmetric !Performance metric to use. -1 = failure cost, 0 = probability of failure, 1 or more = average differntial settlement to the power of the costmetric value (higher values more heavily penalise excessive values)
integer, intent(in) :: EA_generation !current generation of the EA algorithm
integer, intent(in) :: finaloutput !output mode for EA. 1 = final only, 2 = save best of each EA generation, 3 = save full final population, 4 = save everything

integer,intent(in) :: goodcases                               !number of good piles in goodpiles 

!local variables
integer heat_extents(2) !the scope of the site investigation field when using the heat map mode
integer i,j !loop counter
integer tempEAit !temporary EA iteration count
character(1000) str2,str3,str4

!local variables
real tempvec(nrep_MC) !temporary array
real tempmean !temp mean
real tempsd !temp standard deviation



    if (deterministic < 3 .or. finaloutput == 3 .or. finaloutput >= 4) then
    
    	!Output site investigation information from a deterministic run
    	if (deterministic == 1) then

    	  	write(str2,'(A,I0,A,I0,A,A)') 'si_results/Population-stats_piles-',goodcases,'_BArea-',nint(buildingarea),trim(soildsc),'.txt'
    		write(str3,'(A,I0,A,I0,A,A)') 'si_results/Population-Xcoords_piles-',goodcases,'_BArea-',nint(buildingarea),trim(soildsc),'.txt'
			write(str4,'(A,I0,A,I0,A,A)') 'si_results/Population-Ycoords_piles-',goodcases,'_BArea-',nint(buildingarea),trim(soildsc),'.txt'

			
  
    	
			open(505,file=str2)
			write(505,'(I0,X,A)') ninv,	'!total number of investigations'
			write(505,'(I0,X,A)') nbh, '!max no. boreholes'
			write(505,'(I0,X,A)') in_tests, '!No. Tests'
			write(505,'(I0,X,A)') in_depths, '!No. BH depths'
			write(505,'(I0,X,A)') in_reductions, '!No reduction methods'
			write(505,'(A6,X,A4,X,A4,X,A9,X,A14,X,A14,X,A11,X,A10,X,A15,XA14,X,A16)') 'No._BH','Test','Red.','Depth_(m)','Fail._cost_($)','Pile_cost_($)','SI_cost_($)','P._Failure','Ave._Diff._Set.','Geo._statistic','No._bad_reps_(%)'
			do i = 1,ninv
				write(505,'(I6,X,A4,X,A4,X,F9.2,X,F14.2,X,F14.2,X,F11.2,X,F10.6,X,E15.6,X,F14.11,X,F16.5)') inv_bh(i),trim(testnames(inv_test(1,i))),rednames(inv_reduction(i)),dz*inv_depths(1,i,2),fcost(i),pcost(i),icost(i),probfail(i),avediff(i),diffgeo(i),100*real(nrep_MC-invcount(i))/nrep_MC
			end do
			close(505)
			
			!output all borehole coordinates in metres from the corner of the virtual soil
			!Each column corresponds to a borehole. Each row corresponds to an investigation/population member.
			open(505,file=str3)
			open(506,file=str4)
			write(505,'(1000000(I8))') (j,j=1,nbh)
			write(506,'(1000000(I8))') (j,j=1,nbh)
			do i = 1,ninv
				write(505,'(1000000(F8.2))') (dz*inv_coords(i,j,1),j=1,inv_bh(i))
				write(506,'(1000000(F8.2))') (dz*inv_coords(i,j,2),j=1,inv_bh(i))
			end do
			close(505)
			close(506)
            
        else if (deterministic == 2) then !export heat map to file
        
        	!calculate dimensions of the heatmap
         	heat_extents(1) = ((swidth(1)-soffset(1))/sstep(2)+1)
            heat_extents(2) = ((swidth(2)-soffset(2))/sstep(2)+1)

			
    		!if finalouput <= 2 then save a single heatmap under the default name
            if(finaloutput <= 2) then
            	write(str2,'(A,I0,A,I0,A,A,A,A,A,I0,A,A)') 'si_results/Heatmap_piles-',goodcases,'_BArea-',nint(buildingarea),'_test-',trim(testnames(inv_test(1,1))),'_reduction-',rednames(inv_reduction(1)),'_depth-',nint(inv_depths(1,1,2)*dz),trim(soildsc),'.txt'
            	call export_heatmap(heat_extents,soffset,swidth,sstep,dz,str2,si_performance)
            
            !if finaloutput > 2, then total cost, probability of failure, and differentia settlement (potentially weighted) will be exported as individual maps
            else
				if (costmetric == -1 .or.  finaloutput > 2) then	
					write(str2,'(A,I0,A,I0,A,A,A,A,A,I0,A,A)') 'si_results/Heatmap-totalcost_piles-',goodcases,'_BArea-',nint(buildingarea),'_test-',trim(testnames(inv_test(1,1))),'_reduction-',rednames(inv_reduction(1)),'_depth-',nint(inv_depths(1,1,2)*dz),trim(soildsc),'.txt'
					si_performance = fcost + pcost + icost !store the total cost, not including the site investigation cost
					call export_heatmap(heat_extents,soffset,swidth,sstep,dz,str2,si_performance)
				end if
				if (costmetric == 0 .or.  finaloutput > 2) then 
					write(str2,'(A,I0,A,I0,A,A,A,A,A,I0,A,A)') 'si_results/Heatmap-probfail_piles-',goodcases,'_BArea-',nint(buildingarea),'_test-',trim(testnames(inv_test(1,1))),'_reduction-',rednames(inv_reduction(1)),'_depth-',nint(inv_depths(1,1,2)*dz),trim(soildsc),'.txt'
					call export_heatmap(heat_extents,soffset,swidth,sstep,dz,str2,probfail)			
				end if
				if (costmetric > 0 .or.  finaloutput > 2) then !average (potentially weighted) differential settlement
					write(str2,'(A,I0,A,I0,A,A,A,A,A,I0,A,A)') 'si_results/Heatmap-avediff_piles-',goodcases,'_BArea-',nint(buildingarea),'_test-',trim(testnames(inv_test(1,1))),'_reduction-',rednames(inv_reduction(1)),'_depth-',nint(inv_depths(1,1,2)*dz),trim(soildsc),'.txt'
					call export_heatmap(heat_extents,soffset,swidth,sstep,dz,str2,avediff)		
				end if
			
				!if finaloutput > 3, also save the failure cost and pile cost as individual heatmaps. 
				if (finaloutput > 3) then	
					write(str2,'(A,I0,A,I0,A,A,A,A,A,I0,A,A)') 'si_results/Heatmap-failcost_piles-',goodcases,'_BArea-',nint(buildingarea),'_test-',trim(testnames(inv_test(1,1))),'_reduction-',rednames(inv_reduction(1)),'_depth-',nint(inv_depths(1,1,2)*dz),trim(soildsc),'.txt'
					call export_heatmap(heat_extents,soffset,swidth,sstep,dz,str2,fcost)
					write(str2,'(A,I0,A,I0,A,A,A,A,A,I0,A,A)') 'si_results/Heatmap-pilecost_piles-',goodcases,'_BArea-',nint(buildingarea),'_test-',trim(testnames(inv_test(1,1))),'_reduction-',rednames(inv_reduction(1)),'_depth-',nint(inv_depths(1,1,2)*dz),trim(soildsc),'.txt'
					call export_heatmap(heat_extents,soffset,swidth,sstep,dz,str2,pcost)
                    write(str2,'(A,I0,A,I0,A,A,A,A,A,I0,A,A)') 'si_results/Heatmap-geostat1_piles-',goodcases,'_BArea-',nint(buildingarea),'_test-',trim(testnames(inv_test(1,1))),'_reduction-',rednames(inv_reduction(1)),'_depth-',nint(inv_depths(1,1,2)*dz),trim(soildsc),'.txt'
					call export_heatmap(heat_extents,soffset,swidth,sstep,dz,str2,diffgeo)
                    write(str2,'(A,I0,A,I0,A,A,A,A,A,I0,A,A)') 'si_results/Heatmap-geostat0_piles-',goodcases,'_BArea-',nint(buildingarea),'_test-',trim(testnames(inv_test(1,1))),'_reduction-',rednames(inv_reduction(1)),'_depth-',nint(inv_depths(1,1,2)*dz),trim(soildsc),'.txt'
					call export_heatmap(heat_extents,soffset,swidth,sstep,dz,str2,diffgeo2)
                    !write(str2,'(A,I0,A,I0,A,A,A,A,A,I0,A,I0,A,I0,A,I0,A)') 'si_results/Heatmap-geostat2_piles-',goodcases,'_BArea-',nint(buildingarea),'_test-',testnames(inv_test(1)),'_reduction-',rednames(inv_reduction(1)),'_depth-',nint(inv_depths(1,2)*dz),'_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy,'.txt'
					!call export_heatmap(heat_extents,soffset,swidth,sstep,dz,str2,diffgeo2)
                    !write(str2,'(A,I0,A,I0,A,A,A,A,A,I0,A,I0,A,I0,A,I0,A)') 'si_results/Heatmap-geostat0_piles-',goodcases,'_BArea-',nint(buildingarea),'_test-',testnames(inv_test(1)),'_reduction-',rednames(inv_reduction(1)),'_depth-',nint(inv_depths(1,2)*dz),'_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy,'.txt'
					!call export_heatmap(heat_extents,soffset,swidth,sstep,dz,str2,diffgeo0)
					!write(str2,'(A,I0,A,I0,A,A,A,A,A,I0,A,I0,A,I0,A,I0,A)') 'si_results/Heatmap-avediff2_piles-',goodcases,'_BArea-',nint(buildingarea),'_test-',testnames(inv_test(1)),'_reduction-',rednames(inv_reduction(1)),'_depth-',nint(inv_depths(1,2)*dz),'_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy,'.txt'
                    !call export_heatmap(heat_extents,soffset,swidth,sstep,dz,str2,avediff2)
				end if
			end if

			!if finaloutput > 1, then save a heatmap for the number of reps		
			if(finaloutput > 1) then
				write(str2,'(A,I0,A,I0,A,A,A,A,A,I0,A,A)') 'si_results/Heatmap-badreps_piles-',goodcases,'_BArea-',nint(buildingarea),'_test-',trim(testnames(inv_test(1,1))),'_reduction-',rednames(inv_reduction(1)),'_depth-',nint(inv_depths(1,1,2)*dz),trim(soildsc),'.txt'
				si_performance = 100*real(nrep_MC-invcount)/nrep_MC
				call export_heatmap(heat_extents,soffset,swidth,sstep,dz,str2,si_performance)
            end if
		
	    else !if specified, output the final population for each investigation run after the EA has finished

    		tempEAit = EA_generation
    		write(str2,'(A,I0,A,I0,A,I0,A,I0,A,A)') 'si_results/Population-stats_Inv-',tempEAit,'_BH-',inv_bh(1),'_piles-',goodcases,'_BArea-',nint(buildingarea),trim(soildsc),'.txt'
        
            !output all borehole coordinates in metres from the corner of the virtual soil
			!Each row is a member of the population. First column is the fitness of the investigation, followed by all the x coordinates of the boreholes, followed by all the y coordinates.
			open(505,file=str2)
			write(505,'(A14,X,A10,X,1000000(I8))') 'Fitness','%_Bad_Reps',(j,j=1,inv_bh(1)),(j,j=1,inv_bh(1))
			do i = 1,ninv
				write(505,'(E14.7,X,F10.6,X,1000000(F8.2))') si_performance(i),100*real(nrep_MC-invcount(i))/nrep_MC,(dz*inv_coords(i,j,1),j=1,inv_bh(i)),(dz*inv_coords(i,j,2),j=1,inv_bh(i))
			end do
			close(505)
		end if
        
    end if
        

        
        !if (save_effE == 1) then !0 = No. 1 = save 2= load.
        !    write(str2,'(A,I0,A,I0,A,I0,A)') 'si_results/effective_E_',trim(soildsc),'.txt'
        !    open(505,file=str2)
        !    do i = 1,nrep_mc
        !        write(505,'(I6,X,100000F8.3)')i,evals(i,:)
        !    end do
        !    close(505)

        !    write(str2,'(A,A,I0,A,A)') trim(datafolder),'si-piles_ninv-',ninv,trim(soildsc),'.dat'
        !    open(505,file=str2,access='stream')
        !    write(505) si_pset
        !    close(505)
        !end if
        
        




end subroutine





    !save heatmap to text files. This subroutine is called from PROCESS_SI and PROCESS_SI_MULTI where the file is named
    subroutine export_heatmap(heat_extents,soffset,swidth,sstep,dz,str2,si_performance)
    
    	integer, intent(in) :: heat_extents(2),soffset(2),swidth(2),sstep(2)
    	real, intent(in) :: si_performance(:)
    	real(8), intent(in) :: dz !element size
    	character(1000), intent(in) :: str2 !name of file
    	
    	integer i,j,unit
    	
    	unit = 505
    
    	open(unit,file=str2)
					
		write(unit,'(A)') '   From,      To, Resolution'
		write(unit,'(F8.2,X,F8.2,X,A,I0)') soffset(1)*dz,swidth(1)*dz,'1/',sstep(1)
		write(unit,'(F8.2,X,F8.2,X,A,I0)') soffset(2)*dz,swidth(2)*dz,'1/',sstep(2)
	
		do i = 1,heat_extents(1)
			write(unit,'(100000(E15.6,X))') (si_performance(j),j=(i-1)*heat_extents(1)+1,(i-1)*heat_extents(1)+heat_extents(2))
		end do
			
		close(unit)
    
    end subroutine


!open and prep files that store the final borehole coordinates and investigation statistics from the GA
!these files are open for the majority of the program runtime when in GA mode.
subroutine finalstatfiles()






end subroutine




end module
