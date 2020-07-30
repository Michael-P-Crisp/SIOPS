
subroutine despile(siset, pdepths, rel_loads, load_con, settol, buildingforce, detset,sides, niter, npile,npdepths,ninv,nsi,num_loads,np)
         
!Design pile foundations given either a series of linear elastic pile settlement curves, or
!a young's modulus value and a single deterministic pile curve (useful for single layer site
!investigation analysis)

!Can input multiple different pile load cases, designs will be given for each case

!for the siset data, it needs to be scaled according to the per-pile load. However, because it
!can sometimes be a large matrix requiring a lot of RAM, the code was designed such that instead
!of simply making a scaling a copy of the data, it scales the input data and continually scales it
!back to the original value before proceeding.


! Michael Crisp February 2018

use extrafuncs

  implicit none

	integer, intent(in) :: niter,npile,npdepths,ninv !no. of realisations, ck piles, pile depths, investigations
	integer, intent(in) :: nsi !either number of pile depths or 1, which implies that youngs moduli values are given
	integer, intent(in) :: num_loads !number of relative pile loads

	real :: siset(nsi,niter,ninv,np)
	real, intent(in) :: pdepths(npdepths)
	real, intent(in) :: rel_loads(num_loads)		!set of relative pile loads based on tributary area of building
	integer, intent(in) :: load_con(npile)			!connectivity vector saying which pile has which load
	real, intent(in) :: detset(npdepths)			!deterministic pile settlement (for when soil stiffness is provided for SI)
	
	real, intent(in) :: settol !pile design tolerance (mm)
	real, intent(in) :: buildingforce !Total weight of building
	integer, intent(in) :: np !number of piles associated with distinct Young's modulus data (either =1 or =npile)
	
	
	real sum_loads !sum of pile loads
	real nan !not a number
	
	logical :: maskarr(npile*npile) !mask array. Declaring it here so a temporary array isn't repeatedly created
	logical :: maskarr2(npile)
	
	integer :: i,pile,inv,iter,p2,count !loop counters
	integer :: sisize !number of points to design from site investigation data
	integer :: pile_con(npile)	!pile connectivity vector pointing to subset of site investigation data

	
	!output variables
	real :: sides(num_loads,niter,ninv)         ! <----------------- OUTPUT!!! ----------------
	
	!get sum of pile loads
	sum_loads = 0.0
	do i = 1,npile
        if(load_con(i) == 0) cycle  !skip piles with no associated load
		sum_loads = sum_loads + rel_loads(load_con(i))
	end do
	
	nan = 0.0 !generate undefined value
  	nan = 0.0/nan
  	
  	!sisize = no. realisation * no. investigations
  	sisize = size(siset(1,:,:,1))
  	
  	!which site investigation data to refer to for each pile case design
  	if (np == 1) then				!only 1 pile dataset provided, scale all values from this
  		pile_con = 1
  	else if (np == num_loads) then	!unique Young's modulus given for each pile. Use appropriate subset.
  		pile_con = load_con
  	else 
  		write(*,*) "Error, no. input piles must equal 1 or no. of pile loads"
  		write(*,*) "No. piles is:", np
  		write(*,*) "No. piles loads is:",num_loads
        read(*,*)
  		stop
  	end if
  
! ---- get pile designs ----

	if(nsi == 1) then !assume young's moduli was given, calculate using fast method (e.g for SI piles)
		
		siset = sum_loads * siset * settol /  buildingforce
		do i = 1,num_loads
			siset(1,:,:,pile_con(i)) = siset(1,:,:,pile_con(i)) / rel_loads(i) !scale siset by the ratio of current to previous relative load
			call akima(detset(npdepths:1:-1), pdepths(npdepths:1:-1), siset(1,:,:,pile_con(i)), npdepths, 1, 1, 1, sisize, sides(i,:,:))
			siset(1,:,:,pile_con(i)) = siset(1,:,:,pile_con(i)) * rel_loads(i)
		end do
		siset = buildingforce * siset / (sum_loads * settol) !return siset to original values

		
	else if (nsi == npdepths) then !interpolate individual curves	(e.g for CK piles)
		siset = siset * buildingforce / sum_loads
		do i = 1,num_loads
			siset(:,:,:,pile_con(i)) = siset(:,:,:,pile_con(i)) * rel_loads(i)
			call akima(siset(npdepths:1:-1,:,:,pile_con(i)), pdepths(npdepths:1:-1), [settol], npdepths, sisize, 1, 1, 1, sides(i,:,:))	
			siset(:,:,:,pile_con(i)) = siset(:,:,:,pile_con(i)) / rel_loads(i) 
		end do
		siset = sum_loads * siset / buildingforce !return siset to original values
        
        
	else
	  	write(*,*) "Pile depth input points must be equal to number of pile depths, or"
  		write(*,*) "equal to 1 (implying a Young's modulus is given to scale deterministic settlement)"
  		write(*,*) "You have given:", nsi
        read(*,*)
  		stop
	end if
	
	
	


    end subroutine
    
    
    
    
    
    
    
    
!Same as above, but find pile depth in terms of array index, rounded up to the nearest one
subroutine despile_discrete(siset, rel_loads, load_con, settol, buildingforce, detset,sides, niter, npile,npdepths,ninv,nsi,num_loads,np)
         
!Design pile foundations given either a series of linear elastic pile settlement curves, or
!a young's modulus value and a single deterministic pile curve (useful for single layer site
!investigation analysis)

!Can input multiple different pile load cases, designs will be given for each case

!for the siset data, it needs to be scaled according to the per-pile load. However, because it
!can sometimes be a large matrix requiring a lot of RAM, the code was designed such that instead
!of simply making a scaling a copy of the data, it scales the input data and continually scales it
!back to the original value before proceeding.


! Michael Crisp February 2018

use extrafuncs

  implicit none

	integer, intent(in) :: niter,npile,npdepths,ninv !no. of realisations, ck piles, pile depths, investigations
	integer, intent(in) :: nsi !either number of pile depths or 1, which implies that youngs moduli values are given
	integer, intent(in) :: num_loads !number of relative pile loads

	real :: siset(nsi,niter,ninv,np)
	!integer, intent(in) :: pdepths(npdepths)
	real, intent(in) :: rel_loads(num_loads)		!set of relative pile loads based on tributary area of building
	integer, intent(in) :: load_con(npile)			!connectivity vector saying which pile has which load
	real, intent(in) :: detset(npdepths)			!deterministic pile settlement (for when soil stiffness is provided for SI)
	
	real, intent(in) :: settol !pile design tolerance (mm)
	real, intent(in) :: buildingforce !Total weight of building
	integer, intent(in) :: np !number of piles associated with distinct Young's modulus data (either =1 or =npile)
	
	
	real sum_loads !sum of pile loads
	real nan !not a number
	
	logical :: maskarr(npile*npile) !mask array. Declaring it here so a temporary array isn't repeatedly created
	logical :: maskarr2(npile)
	
	integer :: i,j,pile,inv,iter,p2,count !loop counters
	integer :: sisize !number of points to design from site investigation data
	integer :: pile_con(npile)	!pile connectivity vector pointing to subset of site investigation data
    
    integer :: actnumloads !actual number of pile loads, as opposed to the full potential set
    integer :: loadcases(num_loads)		!actual set of pile load cases

	
	!output variables
	integer :: sides(niter,ninv,num_loads)         ! <----------------- OUTPUT!!! ----------------
	
	!get sum of pile loads
	sum_loads = 0.0
	do i = 1,npile
        if(load_con(i) == 0) cycle  !skip piles with no associated load
		sum_loads = sum_loads + rel_loads(load_con(i))
	end do
	
	nan = 0.0 !generate undefined value
  	nan = 0.0/nan
  	
  	!sisize = no. realisation * no. investigations
  	sisize = size(siset(1,:,:,1))
  	
  	!which site investigation data to refer to for each pile case design
  	if (np == 1) then				!only 1 pile dataset provided, scale all values from this
  		pile_con = 1
        actnumloads = 0
        !search for load cases in current set
        do i = 1,num_loads !loop through cases
            do j = 1,npile !search through vector
                if (load_con(j) == i) then 
                    actnumloads = actnumloads + 1
                    loadcases(actnumloads) = i
                    exit
                end if
            end do
        end do
  	else if (np == num_loads) then	!unique Young's modulus given for each pile. Use appropriate subset.
  		pile_con = load_con
  	else 
  		write(*,*) "Error, no. input piles must equal 1 or no. of pile loads"
  		write(*,*) "No. piles is:", np
  		write(*,*) "No. piles loads is:",num_loads
        read(*,*)
  		stop
    end if
    
    sides = 0
  
! ---- get pile designs ----

	if(nsi == 1) then !assume young's moduli was given, calculate using fast method (e.g for SI piles)
		
		siset = sum_loads * siset * settol /  buildingforce
		do i = 1,actnumloads
			siset(1,:,:,pile_con(i)) = siset(1,:,:,pile_con(i)) / rel_loads(loadcases(i)) !scale siset by the ratio of current to previous relative load
			call getdesind(detset, siset(:,:,:,pile_con(i)), npdepths, 1, 1, 1, sisize, sides(:,:,loadcases(i)))
			siset(1,:,:,pile_con(i)) = siset(1,:,:,pile_con(i)) * rel_loads(loadcases(i))
		end do
		siset = buildingforce * siset  / (sum_loads * settol) !return siset to original values

		
	else if (nsi == npdepths) then !interpolate individual curves	(e.g for CK piles)
		siset = siset * buildingforce / sum_loads
		do i = 1,num_loads
			siset(:,:,:,pile_con(i)) = siset(:,:,:,pile_con(i)) * rel_loads(i)
			call getdesind(siset(:,:,:,pile_con(i)), [settol], npdepths, sisize, 1, 1, 1, sides(:,:,i))	
			siset(:,:,:,pile_con(i)) = siset(:,:,:,pile_con(i)) / rel_loads(i) 
		end do
		siset = sum_loads * siset / buildingforce !return siset to original values
        
        
	else
	  	write(*,*) "Pile depth input points must be equal to number of pile depths, or"
  		write(*,*) "equal to 1 (implying a Young's modulus is given to scale deterministic settlement)"
  		write(*,*) "You have given:", nsi
        read(*,*)
  		stop
	end if
	
	
	


end subroutine

