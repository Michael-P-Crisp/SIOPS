!c  *********************************************************************
!c  *                                                                   *
!c  *                            program soilgen                        *
!c  *                                                                   *
!c  *********************************************************************
!c  Single Precision Version 1.0
!c  Written by Michael P Crisp

!c  Latest Update: March, 2017
!c
!c  PURPOSE Generate multi-layer soil profiles
!c
!	This subroutine produces and outputs the final soil, consisting of an 
!	arbitrary number of layers, each with a uniform set of soil properties.
!
!	This module should be able to be placed within a larger program that makes use
!	of spatially-variable random soils, including (but not limited to):
!		Settlement modelling,
!		Site investigation optimisation
!		Slope stability analysis
!		Reliability-based design
!		Groundwater modelling
!
!	Note it uses pre-calculated correlation matrices generated via the las3i subroutine.
!c
!c  REVISION HISTORY:
!c-------------------------------------------------------------------------
module soilgen

    use variables
    use int2d
    use sim2sd

	implicit none

	contains

	subroutine soil_layers(xyi)	
		

		integer i,j,ix,iy,iz,outpx,outpy,outpz,L,istat !loop counters
		real randu  !function for random number generation

		logical reset !flag indicating whether the soil correlation matrix needs to be (re)calculated - only calculated for initialization 
		real(8) :: tbfld(nlayer-1,nxew,nyew),tbfld2(nlayer-1,nxew,nyew)  !temp arrays for boundary information
		integer x,y,z !loop counters
        real bsdata(4)
        !integer, optional, intent(in) :: nrep

        real(8) fielddepthsmulti(nlayer-1,nbhmulti)      !layer information from random field at borehole locations (to cancel out discrepencies)
        real(8), intent(in) :: xyi(:,:)
        
        logical, parameter :: debug = .false.
        !character(100) str2


		!------------------------------ GENERATE MULTI-LAYER SOIL PROFILE --------------------------------
		reset=.false. !don't regenerate soil correlation matrices
        
        bsdata= [0.0d+0,bsd,0.0d+0,0.0d+0]
        
        !Set top and bottom to be top and bottom of the field
        bfld(1,:,:) = 1
        bfld(nlayer+1,:,:) = nzew
        

		

		!increment random number generator
		kseed = randu(kseed)* 1234567890    
        
        !get 2D random fields for soil layers
        if ( bsd == 0) then
            bfld(2:nlayer,:,:) = 0.0
        else
            do i=1,nlayer-1
		
				!call LAS for zero mean, normally-distributed 2D random field for random noise of layer boundary (just uses top slice atm)
				call sim2d(bfld(i+1,:,:),nxe,nye,nxew,nyew,ceiling(real(5*nye)/4.0),dx,dy,kseed,MXM,MXK, &
				bC0,bCT,bCC,bCS,bCI,bAT,bAC,bAS,bAI, &
				bM,bk1,bk2,bkk,bsdata,'n') !hardcode distribution to be normal with zero mean for now
                
                !if(present(nrep)) then
                !    
                !    write(str2,'(A,I0,A)') 'bfld',nrep,'.txt'
                !    open(1351,file=str2)
                !
                !    do j = 1,nxew
                !        write(1351,'(10000(F6.2,X)))') bfld(2,j,:)
                !    end do
                !
                !    close(1351)
                !    
                !end if 


            end do
        end if
            
            !Create layers by merging random noise with mean layer geology
            if (uselcoords) then !Interpolate boundary from coordinates
                !write(*,*) bhxymulti
                call interp_layers(2,nbhmulti,nxew,nyew,bhxymulti,bhdepthsmulti,xyi,tbfld,nlayer) !,empty1,empty2,.true.)
                !write(*,*) bhxymulti
                if(enforcelcoords) then
                    do j = 1,nbhmulti
                        fielddepthsmulti(1:nlayer-1,j) = bfld(2:nlayer,bhxymulti(1,j),bhxymulti(2,j)) !get locations of random field at borehole locations
                    end do
                    call interp_layers(2,nbhmulti,nxew,nyew,bhxymulti,fielddepthsmulti,xyi,tbfld2,nlayer) !generate interpolated layer based on layer discrepencies
                    bfld(2:nlayer,:,:) = bfld(2:nlayer,:,:) + tbfld - tbfld2        !Create final field by adding random noise to interpolated layer, minus interpolated discrepencies at boreholes       
                else
                    bfld(2:nlayer,:,:) = bfld(2:nlayer,:,:) + tbfld                 !create field by adding random noise to interpolated layer.
                end if
            else !otherwise just use the given depths
                do i = 1,nlayer-1
                    bfld(i+1,:,:) = bfld(i+1,:,:) + ldepths(i)
                end do
            end if
            
            if (debug) then !save the various boundary components to disk for checking
                open(999,file='bfldrandom.txt')
                open(998,file='bfldmean.txt')
                open(997,file='bfldanti.txt')
                do y = 1,nyew
                    write(999,'(10000000(F8.4,X))')  bfld(2,:,y) !- tbfld(1,:,y) + tbfld2(1,:,y)
                    write(998,'(10000000(F8.4,X))') tbfld(1,:,y)
                    write(997,'(10000000(F8.4,X))') tbfld2(1,:,y)
                end do
                close(999)
                close(998)
                close(997)
            end if
            
            
            
            
            !Mimic processess of erosion and deposition by ensuring that if a newer layer is lower than an older layer,
            !the new layer eats into the older layer.
            do i = nlayer-1,1,-1
                where(bfld(i,:,:) > bfld(i+1,:,:)) bfld(i+1,:,:) = bfld(i,:,:)  
            end do
            
            !Ensure boundaries are always within the soil limits
            do y = 1,nyew
                do x = 1,nxew
                    do i = nlayer-1,2,-1
                        if(bfld(i,x,y) > nzew) bfld(i,x,y) = nzew
                        if(bfld(i,x,y) < 1) bfld(i,x,y) = 1
                        if(bfld(i,x,y) > bfld(i+1,x,y)) bfld(i+1,x,y) = bfld(i,x,y)
                    end do
                end do
            end do
            
                !open(999,file='bfld1.txt')
                !open(998,file='bfld2.txt')
                !do y = 1,nyew
                !    write(999,'(10000000(F8.4,X))')  bfld(2,:,y)
                !    write(998,'(10000000(F8.4,X))')  bfld(3,:,y)
                !end do
                !close(999)
                !close(998)
                
                



	end subroutine

end module
