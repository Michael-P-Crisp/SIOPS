
module getdiff

contains

subroutine getdiff_cts(ckset, sides, pdepths, rel_loads, load_con, plocation,buildingforce,diffset, niter, npile,npdepths,ninv,num_loads)
         
!calculate maximum differential settlement between piles that support a building,
!based on linear elastic settlement curves with depth that represent true pile deformation
!and pile lengths designed from a site investigation

!supports different loads applied to each pile

! Michael Crisp February 2018

use extrafuncs

  implicit none

	integer, intent(in) :: niter,npile,npdepths,ninv !no. of realisations, ck piles, pile depths, investigations
	integer, intent(in) :: num_loads !number of relative pile loads

	real, intent(in) :: ckset(:,:,:) !(npdepths,niter,npile) !ck pile settlement curve
	real, intent(in) :: pdepths(:) !(npdepths)
	real, intent(in) :: rel_loads(num_loads)		!set of relative pile loads based on tributary area of building
	integer, intent(in) :: load_con(npile)			!connectivity vector saying which pile has which load
	real, intent(in) :: plocation(npile,2)			!x,y locations for each pile
	real, intent(in) :: sides(num_loads,niter,ninv)
	
	real, intent(in) :: buildingforce !Total weight of building
	
	real cktemp(size(ckset,1))
	
	real sum_loads !sum of pile loads
	real trueset(ninv,npile)
	real dist(npile*npile) !distances between piles
	real difftemp(npile*npile) !distances between piles
	real nan !not a number
	real rel_loads2(0:num_loads) !want to store a 1 in the zero'th element
	
	!logical :: maskarr(npile*npile) !mask array. Declaring it here so a temporary array isn't repeatedly created
	logical :: maskarr2(npile)
	
	integer :: i,pile,inv,iter,p2,count !loop counters
	
	!output variables
	real, intent(out) :: diffset(niter,ninv)
	
	
	!get sum of pile loads
	sum_loads = 0.0
	do i = 1,npile
		sum_loads = sum_loads + rel_loads(load_con(i))
	end do
	
	nan = 0.0 !generate undefined value
  	nan = 0.0/nan
  	
  	!make copy of relative load vector with 1 at the start
  	rel_loads2(0) = 1.0
  	rel_loads2(1:) = rel_loads
	
	
! ---- Get differential settlements ----
	
	!calculate distances between each pile here instead of in the heavily nested loop below
	count = 0
	do pile = 1,npile-1
		do p2 = pile+1,npile
			count = count + 1
			dist(count) = sqrt((plocation(pile,1)-plocation(p2,1))**2 + (plocation(pile,2)-plocation(p2,2))**2)
		end do
	end do
	
	!calculate the max differential settlements
	do iter = 1,niter
		!get true settlement of all individual piles for current realisation and investigation
		do pile = 1,npile
			cktemp = ckset(:,iter,pile) * buildingforce * rel_loads(load_con(pile)) / sum_loads
			call akima(pdepths, cktemp,sides(load_con(pile),iter,:),npdepths,1,1,1,ninv,trueset(:,pile)) !sides(load_con(1),iter,:) trueset(:,1) ckset(:,iter,1)
		end do
		!calculate differential settlement between all combinations of piles
		do inv = 1,ninv
			!maskarr2 = isnan(trueset(inv,:))	!if any of the piles don't have settlement values, discard the realisation for that investigation
			maskarr2 = trueset(inv,:) < 0
            if( any( maskarr2 )) then
				diffset(iter,inv) = -100
				cycle
			end if
			count = 0
			do pile = 1,npile-1
				do p2 = pile+1,npile
					count = count + 1
					difftemp(count) = abs(trueset(inv,pile) - trueset(inv,p2))/dist(count)
				end do
			end do
			!maskarr(:count) = difftemp(:count) /= nan		!Take max differential settlement, ignoring nan values
			!diffset(iter,inv) = maxval(difftemp(:count),mask=maskarr(:count))
			diffset(iter,inv) = maxval(difftemp(:count))
		end do
	end do
	
	


    end subroutine

    
    
    
    

subroutine getdiff_discrete(ckset, sides, pdepths, rel_loads, load_con, plocation,buildingforce,diffset, niter, npile,npdepths,ninv,num_loads)
         
!calculate maximum differential settlement between piles that support a building,
!based on linear elastic settlement curves with depth that represent true pile deformation
!and pile lengths designed from a site investigation

!supports different loads applied to each pile

! Michael Crisp February 2018

use extrafuncs

  implicit none

	integer, intent(in) :: niter,npile,npdepths,ninv !no. of realisations, ck piles, pile depths, investigations
	integer, intent(in) :: num_loads !number of relative pile loads

	real, intent(in) :: ckset(:,:,:) !ckset(npdepths,niter,npile) !ck pile settlement curve
	real, intent(in) :: pdepths(:) !pdepths(npdepths)
	real, intent(in) :: rel_loads(:) !rel_loads(num_loads)		!set of relative pile loads based on tributary area of building
	integer, intent(in) :: load_con(:) !load_con(npile)			!connectivity vector saying which pile has which load
	real, intent(in) :: plocation(:,:) !plocation(npile,2)			!x,y locations for each pile
	integer, intent(in) :: sides(:,:,:) !sides(niter,ninv,num_loads)
	
	real, intent(in) :: buildingforce !Total weight of building
	
	real cktemp(npdepths)
    !real diff1,diff2,diff3,diff4
    !integer cr,s1,f1
	
	real sum_loads !sum of pile loads
	real trueset(npile)
	real dist(npile*npile) !distances between piles
	real difftemp(npile*npile) !distances between piles
	real nan !not a number
	real rel_loads2(0:num_loads) !want to store a 1 in the zero'th element
	
	!logical :: maskarr(npile*npile) !mask array. Declaring it here so a temporary array isn't repeatedly created
	logical :: maskarr2(npile)
	
	integer :: i,pile,inv,iter,p2,count !loop counters
	
	!output variables
	real, intent(out) :: diffset(niter,ninv)
	
	
	!get sum of pile loads
	sum_loads = 0.0
	do i = 1,npile
		sum_loads = sum_loads + rel_loads(load_con(i))
	end do
	
	nan = 0.0 !generate undefined value
  	nan = 0.0/nan
  	
  	!make copy of relative load vector with 1 at the start
  	rel_loads2(0) = 1.0
  	rel_loads2(1:) = rel_loads
	
	
! ---- Get differential settlements ----
	
	!calculate distances between each pile here instead of in the heavily nested loop below
	count = 0
	do pile = 1,npile-1
		do p2 = pile+1,npile
			count = count + 1
			dist(count) = sqrt((plocation(pile,1)-plocation(p2,1))**2 + (plocation(pile,2)-plocation(p2,2))**2)
		end do
    end do
    
    
    !CALL system_clock(count_rate=cr)
	
	!calculate the max differential settlements
	do iter = 1,niter
		!get true settlement of all individual piles for current realisation and investigation

		!calculate differential settlement between all combinations of piles
		do inv = 1,ninv
            

            !CALL system_clock(s1)
            
			maskarr2 = sides(iter,inv,load_con) < 0	!if any of the piles don't have settlement values, discard the realisation for that investigation
			if( any( maskarr2 )) then
				diffset(iter,inv) = -1
				cycle
            end if
            
            !CALL system_clock(f1)
            !diff1=diff1+real(f1-s1)/cr
            !CALL system_clock(s1)
            
            !get true settlement, as design index is stored in sides array.
            do pile = 1,npile
                trueset(pile) = ckset(sides(iter,inv,load_con(pile)),iter,pile) * buildingforce * rel_loads(load_con(pile)) / sum_loads !ckset(:,1,2)
            end do
            
        !CALL system_clock(f1)
        !diff2=diff2+real(f1-s1)/cr
        !CALL system_clock(s1)

			count = 0
			do pile = 1,npile-1
				do p2 = pile+1,npile
					count = count + 1
					difftemp(count) = (trueset(pile) - trueset(p2))/dist(count)
				end do
            end do
            
            
            !CALL system_clock(f1)
            !diff3=diff3+real(f1-s1)/cr
            !CALL system_clock(s1)
            
            difftemp(:count) = abs(difftemp(:count))
            
			!maskarr(:count) = difftemp(:count) /= nan		!Take max differential settlement, ignoring nan values
			!diffset(iter,inv) = maxval(difftemp(:count),mask=maskarr(:count))
			diffset(iter,inv) = maxval(difftemp(:count))
            
            !
            !CALL system_clock(f1)
            !diff4=diff4+real(f1-s1)/cr

		end do
    end do
	
    !write(*,*)'diff1',diff1
    !write(*,*)'diff2',diff2
    !write(*,*)'diff3',diff3
    !write(*,*)'diff4',diff4
    !read(*,*)
	   !


end subroutine

end module

    