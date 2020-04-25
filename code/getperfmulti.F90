module getperfmulti
    
    use int2d
    use reducem
    
    contains
    
subroutine si_multi_stuff(nxe,nye,nze,emn,bfld,nlayer,n_inv,n_bh,inv_bh,inv_coords,inv_depths,num_tests,inv_test,add_errors,use_CI,CI,test_errors,num_errors,kseed, &
						plocation,pdim,npl,power,mindist,percentile,s_dev,n_reds,inv_reduction,invpower,evals,layer_at_piles,sdist,sumweights,extents, &
                        node_xy_out,zd_out,element_neighbor,triangle,element_num,extranum,allones,use_avedepths,xyi,multitype,str2,plocation_dble,iter,finaloutput,soil_reps,rand_realisations,lmean_ln,lsd_ln)

!*****************************************************************************
! Get SI performance (reduced young's modulus of each layer + effective layer depth for each pile)
!

  implicit none
  
  integer, intent(in) :: pdim(2) !pile radius
  integer, intent(in) :: iter !current monte carlo realisation
  integer, intent(in) :: finaloutput !integer relating to the quantity of information to output
  integer, intent(in) :: soil_reps(2) !range of monte carlo realisations for which to output virtual soils
  
  !all of this is site investigation stuff 
  integer,		intent(in) :: nxe,nye,nze						!number of soil elements in x,y,z dimensions
  !real, 		intent(in) :: efld(nxe,nye,nze) 				!3D array of soil E values
  real,			intent(in) :: emn(:)							!mean soil stiffness in each layer
  real, 		intent(in) :: bfld(:,:,:)  !bfld(nlayer+1,nxe,nye)			!n x 2D array of n+1 layer boundaries. Bottom boundary for layer i is "bfld(i+i,:,:) + 1".
  integer,		intent(in) :: nlayer							!number of soil layers
  integer, 		intent(in) :: n_inv 							!total number of investigations
  integer, 		intent(in) :: n_bh								!max number of boreholes
  integer, 		intent(in) :: inv_bh(n_inv)						!number of boreholes in each investigation
  integer,		intent(in) :: inv_coords(n_inv,n_bh,2)			!x,y coordinates of boreholes for each investigation
  integer, 		intent(in) :: inv_depths(n_bh,n_inv,3) 				!sampling depth start, end, interval
  integer,		intent(in) :: num_tests							!number of tests
  integer, 		intent(in) :: num_errors						!number of error sources for tests
  integer,	    intent(in) :: inv_test(n_bh,n_inv)					!test type for each investigation
  logical,		intent(in) :: add_errors						!whether to add test errors or not
  logical,		intent(in) :: use_CI							!whether to truncate based on confidence interval
  real,			intent(in) :: CI    							!confidence interval to be applied to each investigation
  real, 		intent(in) :: test_errors(num_tests,num_errors)	!matrix of test error stats (num tests x num errors)
  integer,      intent(in) :: kseed                             !initial random seed
  
  real, intent(in) :: lmean_ln(:),lsd_ln(:) !lognormal statistics for each layer for generating random values      <- revert; comment this out
  logical, intent(in) :: rand_realisations !randomize true soil uniform values (true) or apply white noise to soil samples (false)
  
  real(8) :: inv_coords_real(2,n_bh)			!transposed, single precision equivallent of the inv_coords array for subroutine compatibility
  real :: sdist(npl,nxe,nye),sumweights(npl) !distance of elements from each pile, sum of distance weights
  integer :: extents(npl,2,2)!extent of soil around pile
  integer :: indices(2,nxe,nye)
  logical, intent(in) :: allones !true if all investigations use a single borehole
  logical, intent(in) :: use_avedepths !true if a horizontal boundary for each layer is used, where the depth is the average of all boreholes
  
  
! -- output variables --
  real			:: goodvals(n_inv,nlayer,n_bh,nze)                   !store valid values of SI samples in 1D array
  integer		:: nsamples(n_inv,nlayer,n_bh)						!keep track of how many viable samples are in each layer in each borehole
  real		:: bound_depths(n_inv,nlayer-1,n_bh)			!layer boundary information as encountered by each borehole and investigation    bound_depths(n_inv,nlayer+1,n_inv,n_bh)
  
  real(8)		:: bound_depths_dble(nlayer-1,n_bh)			!layer boundary information as encountered by each borehole and investigation    bound_depths(n_inv,nlayer+1,n_inv,n_bh)
  
  real, intent(out) :: evals(:,:) !evals(n_inv,nlayer)
  real, intent(out) :: layer_at_piles(:,:,:) !n_inv,nlayer-1,npl
  real(8) :: lap(nlayer-1,npl*product(pdim)) !temp version of the above
  
! -- local variables --

  integer		inv,bh,j,k,i,d,layer,red,n,x,y                        !loop counters
  integer		depths(nze),scount							!vector containing depths of samples, no. samples
  logical 		inside(nze) 								!samples which are inside a given layer
  integer       tscount                                     !total sample count in current investigation
  integer		stepback									!number of elements to step backwards for discrete sampling (depth interval/2)
  real			empty1(1,1),empty2(1,1)
  
  !2D interpolation
  real(8), intent(in)  :: xyi(:,:) !coordinates of points to interpolate
  real(8) :: bfldt(nlayer-1,nxe,nye)
  real :: bfldt_single(nlayer-1,nxe,nye)
  
  !pre-processed triangulation variables
  
    real :: node_xy_out(:,:,:)
    integer :: element_neighbor(:,:,:),triangle(:,:,:)
    integer :: zd_out(:,:)
    integer :: element_num(:) !number of triangles in interpolated surface
    integer:: extranum(:) !number of extrapolated points
  
  
  !height stuff

  integer, intent(in) :: plocation(:,:) !plocation(npl,2)
  real(8), intent(in) :: plocation_dble(:,:) !2,npl*product(pdim))
  integer, intent(in) :: npl
  integer,intent(in) :: invpower !inverse distance power for effective layer depths
  integer, intent(in) :: multitype !whether the layer depths at piles are treated as (1) a single point in SI and CK, (2) inverse weighted in CK, (3) inverse weighted in SI and CK
  integer parea !number of elements contained within pile
  

  logical :: truedepths = .false.
  
  !the depth at which to place layers which are not encountered by the borehole. Typically just below the base of the borehole, unless this is outside the soil limit
  integer belowdepth 
  logical isbelow(n_bh) !whether all depths for a layer are below the boreholes
  
  
  !reduction stuff
  integer,		intent(in) :: mindist							!minimum distance between footing and investigation for inverse distance methods
  real,			intent(in) :: power								!power of inverse-distance weighting
  integer, 		intent(in) :: inv_reduction(n_inv) !n_reds)						!reduction method for each test	
  real,			intent(in) :: percentile						!percentile to use in reduction method (proportion)
  real,			intent(in) :: s_dev								!standard deviations below mean for reduction method (geometric)
  integer,		intent(in) :: n_reds
  real			 :: redvals(nlayer,npl)				!reduced value for each investigation
  
  !timing stuff
  real start,finish
  
  character(1000) str2 !this variable is only an input so that it can be shared withthe calling subroutine, to avoid allocating the memory 1000's of times
  

  
  !write(*,*) 'here0'
  ! #########################		GET STIFFNESS PROPERTIES OF EACH LAYER		#########################	
  
  !Do site investigations and collect samples from soil
  
   if(add_errors .or. .not. rand_realisations) then !if dealing with errors, do the full site investigation proceedure

		call do_SI(goodvals,nsamples,nxe,nye,nze,emn,bfld,nlayer,n_inv,n_bh,inv_bh,inv_coords,inv_depths,num_tests,num_errors,inv_test,add_errors,use_CI,test_errors,kseed,CI,rand_realisations,lmean_ln,lsd_ln)

  

		!write(*,*) 'here1'
	
	!    goodvals(1,1,1,:14)

	  !do red = 1,n_reds

		do inv = 1,n_inv
	
		!do the reduction method independently on each investigation
			call reducesi(nze,nlayer,inv_bh(inv),n_bh,inv_coords(inv,:,:),inv_reduction(inv),redvals,percentile,s_dev,npl,plocation,mindist,power,goodvals(inv,:,:,:),nsamples(inv,:,:))

			!write(*,*) 'here1a'
	
			evals(inv,:) = redvals(:,1) !don't include any inverse-distance stuff for now. If you want this, you'll need to restructure the evals array to have an index for all piles
			!save the reduction method for each investigation
            
            ! some python code for ensuring that layers are properly extended by replacing the properties of fully-eroded layers with logical values
            ! this shouldn't be needed since if such a layer is encountered by the site investigation, then the layer thickness should be zero, and it's ignored
            ! if this is to be implemented in the future, it'll need to be converted to fortran syntax.
            !i=0
            !while np.isnan(evals[inv,0]) #ensure that the first layer has soil properties, even if deeper layers are the first sampled
            !    i += 1
            !    evals[inv,0] = evals[inv,i]
            !
            !for i in range(1,nlayer):   #if layers aren't encountered, use the properties of previous layers
            !    if np.isnan(evals[inv,i]): evals[inv,i] = evals[inv,i-1]
            
		end do

	
	  !end do
  
   else
  !Otherwise, the young's modulus for each layer is already known, and can be saved directly (much faster)
        
  		do j = 1,nlayer	
			evals(:,j) = emn(j)
        end do
  
  end if



	
	! #########################		GET LAYER DEPTHS AT BOREHOLES		#########################	
	
	!This only needs to be done for the discrete and continuous test cases, as the heights do not change with testing inaccuracies
	
	!get layer depths at each borehole depending on where layer is encountered in each investigation
  

	!call getbhdepths(nxe,nye,nze,bfld,nlayer,n_inv,n_bh,inv_bh,inv_coords,inv_depths,bound_depths)
  

  
     if(allones) then       !If all the investigations have a single borehole (i.e. with the heatmap mode, and on other occasions), layer depth is constant. Set depth and exit subroutine.
         do inv=1,n_inv
            layer_at_piles(inv,:,1) = bfld(2:nlayer,inv_coords(inv,1,1),inv_coords(inv,1,2))  
            do i = 1,nlayer-1
                !if a layer is below the bottom of the borehole, put it at the bottom of the borehole
                if(layer_at_piles(inv,i,1) > inv_depths(1,inv,2)) then
                    layer_at_piles(inv,i,1) = nze
                !If the current test is discrete round the depths up to the nearest sample, and subtract half the increment 
                else if(inv_depths(1,inv,3) > 1) then
                    layer_at_piles(inv,i,1) = inv_depths(1,inv,3)*ceiling(layer_at_piles(inv,i,1)/inv_depths(1,inv,3))-inv_depths(1,inv,3)/2
                end if
            end do
         end do
        return
     else           !otherwise, get layer depths at boreholes
         
        !the depth at which to place layers which are not encountered by the borehole. Typically just below the base of the borehole, unless this is outside the soil limit
        
    
          do inv = 1,n_inv
    
                do i = 1,nlayer-1
                    isbelow = .false.
                    do bh = 1,inv_bh(inv)	!loop through boreholes
                        belowdepth = min(inv_depths(bh,inv,2) + 1,nze)
			            bound_depths(inv,i,bh) = bfld(i+1,inv_coords(inv,bh,1),inv_coords(inv,bh,2)) !get direct depths from boundary
                        !if a layer is below the bottom of the borehole, put it just below the borehole
                        if(bound_depths(inv,i,bh) > inv_depths(bh,inv,2)) then
                            bound_depths(inv,i,bh) = belowdepth
                            isbelow(bh) = .true.
                        !If the current test is discrete round the depths up to the nearest sample, and subtract half the increment
                        else if(inv_depths(bh,inv,3) > 1) then
                            bound_depths(inv,i,bh) = inv_depths(bh,inv,3)*ceiling(bound_depths(inv,i,bh)/inv_depths(bh,inv,3))-inv_depths(bh,inv,3)/2
                        end if
                    end do
                    !if the layer is below the depth of all boreholes, then set all depths to the bottom of the soil, effectively eliminating it 
                    if(all(isbelow(:inv_bh(inv)))) then
                        bound_depths(inv,i,:) = nze 
                    end if
                end do
          end do
     end if

    
    
    ! #########################		GET EFFECTIVE LAYER DEPTHS AT PILES		#########################
     
     
     if ( use_avedepths ) then !use average horizontal boundaries by taking the average of layer depths at boreholes
        
         do inv=1,n_inv
            layer_at_piles(inv,:,1) = sum(bound_depths(inv,:,:inv_bh(inv)),2)/inv_bh(inv)
         end do
         return
         
     else !otherwise perform interpolation/extrapolation
     
     
	
        if (multitype == 3) then !take effective depths for soil model as inverse-distance weighted averages of interpolated boundary

        

	 
	            do inv = 1,n_inv
			            inv_coords_real = transpose(inv_coords(inv,:,:))
			            bound_depths_dble = bound_depths(inv,:,:)
            
                        !2D interpolation subroutine    
                        call interp_layers_proc(2,inv_bh(inv),nxe,nye,inv_coords_real(:,:inv_bh(inv)),bound_depths_dble(:,:inv_bh(inv)),xyi,bfldt,nlayer,node_xy_out(inv,:inv_bh(inv),:),zd_out(inv,:inv_bh(inv)),element_neighbor(inv,:,:),triangle(inv,:,:),element_num(inv),extranum(inv)) !,empty1,empty2,.true.)
                        bfldt_single = bfldt
                    
                        !-------- exporting the site investigation soil model boundaries ------
                        if(finaloutput == 5 .and. (iter >= abs(soil_reps(1)) .and. iter <= abs(soil_reps(2)))) then
                            do n = 2,nlayer
                                write(str2,'(A,I0,A,I0,A)') 'sibound_iter-',iter,'_layer-',n,'_inv-',inv,'.txt'
	                            open(200,file=str2)
                                do y=1,size(bfldt_single,3) !fill in y direction
                                    write(200,'(1000000000I4)') nint(bfldt_single(n-1,:,y)) !output current row to file
                                end do
                                close(200)
                            end do
                        end if
                        !

                    
                        !Ensure boundaries are always within the soil limits
                        !Mimic processess of erosion and deposition by ensuring that if a newer layer is lower than an older layer,
                        !the new layer eats into the older layer.
                        !This is only really makes a difference if the piles are located outside the of the site investigation, where extrapolation occurs
                    
                        ! Edit: this probably doesn't make any difference considering it's done to the final layer depths,
                        ! certainly not enough to warrent the extra computational time
                        !do y = 1,nye
                        !    do x = 1,nxe
                        !        do i = nlayer-1,1,-1
                        !            if(bfldt_single(i,x,y) > nze) bfldt_single(i,x,y) = nze
                        !            if(bfldt_single(i,x,y) < 1) bfldt_single(i,x,y) = 1
                        !        end do
                        !        do i = nlayer-2,1,-1
                        !            if(bfldt_single(i,x,y) > bfldt_single(i+1,x,y)) bfldt_single(i+1,x,y) = bfldt_single(i,x,y)
                        !        end do
                        !    end do
                        !end do bfldt_single(41,41) bfldt_single(120,120)
                    
                        !--- a check to make sure the interpolated  boundary is the same as the original at the borehole locations ----
                        !if (inv == 1) then
                        !    do bh=1,inv_bh(inv)
                        !        write(*,*) inv_coords(inv,bh,:), bfldt_single(1,inv_coords(inv,bh,1),inv_coords(inv,bh,2)),bound_depths(inv,1,bh)
                        !    
                        !    end do
                        !end if
                        !write(*,*)

            
                        !save height for the relevant investigations for each pile, for each layer
                        do d = 1,npl
                            !write(*,*) extents(d,1,1),extents(d,1,2),extents(d,2,1),extents(d,2,2)
                            do i = 1,nlayer-1 !get the effective heights for each layer.        
	                            layer_at_piles(inv,i,d) = sum(bfldt_single(i,extents(d,1,1):extents(d,1,2),extents(d,2,1):extents(d,2,2)) * sdist(d,extents(d,1,1):extents(d,1,2),extents(d,2,1):extents(d,2,2))) / sumweights(d)
                            
                            end do
                        
                            !Ensure that the erosion rules are followed
                            !Ensure layer depths at pile are within soil limits
                            !This can still be an issue due despite applying this operation to the soil layers due to finite math
                            do i = nlayer-2,1,-1
                                if(layer_at_piles(inv,i,d) > layer_at_piles(inv,i+1,d)) layer_at_piles(inv,i+1,d) = layer_at_piles(inv,i,d) 
                                if(layer_at_piles(inv,i,d) > nze) layer_at_piles(inv,i,d) = nze
                                if(layer_at_piles(inv,i,d) < 1) layer_at_piles(inv,i,d) = 1
                            end do
                            if(layer_at_piles(inv,nlayer-1,d) > nze) layer_at_piles(inv,nlayer-1,d) = nze
                            if(layer_at_piles(inv,nlayer-1,d) < 1) layer_at_piles(inv,nlayer-1,d) = 1
                        
                        end do
                    
                end do
      
            
  
            
        else !simply take the soil depth at each pile as according to the interpolated boundary (much faster)
        
            !call cpu_time(start)


                    parea = product(pdim)

                    do inv = 1,n_inv
            
                        inv_coords_real = transpose(inv_coords(inv,:,:))
                        bound_depths_dble = bound_depths(inv,:,:)
                    
                        call interp_layers_proc(2,inv_bh(inv),nxe,nye,inv_coords_real(:,:inv_bh(inv)),bound_depths_dble(:,:inv_bh(inv)),xyi,bfldt,nlayer,node_xy_out(inv,:inv_bh(inv),:),zd_out(inv,:inv_bh(inv)),element_neighbor(inv,:,:),triangle(inv,:,:),element_num(inv),extranum(inv),plocation_dble,lap)
                    
                    
                    
                        !Ensure that the erosion rules are followed
                        !Ensure layer depths at pile are within soil limits
                        do n = 1,npl
                            do i = nlayer-2,1,-1
                                layer_at_piles(inv,i,n:n) = sum(lap(i,parea*(n-1)+1:parea*n))/real(parea) !repack the layer heights into the array, taking the average of samples within the layer
                                if(layer_at_piles(inv,i,n) > layer_at_piles(inv,i+1,n)) layer_at_piles(inv,i+1,n) = layer_at_piles(inv,i,n) 
                                if(layer_at_piles(inv,i,n) > nze) layer_at_piles(inv,i,n) = nze
                                if(layer_at_piles(inv,i,n) < 1) layer_at_piles(inv,i,n) = 1
                            
                            
                            end do
                            if(layer_at_piles(inv,nlayer-1,n) > nze) layer_at_piles(inv,nlayer-1,n) = nze
                            if(layer_at_piles(inv,nlayer-1,n) < 1) layer_at_piles(inv,nlayer-1,n) = 1
                            layer_at_piles(inv,nlayer-1,n:n) = sum(lap(nlayer-1,parea*(n-1)+1:parea*n))/real(parea)
                        end do
                    
                        !if (inv == 1) then
                        !    write(*,*) layer_at_piles(inv,1,:)
                        !    write(*,*) (bfld(2,inv_coords(inv,i,1),inv_coords(inv,i,2)), i=1,inv_bh(inv)) 
                        !    write(*,*)
                        !    write(*,*) (inv_coords(inv,i,:),i=1,inv_bh(inv))
                        !    do j = 1,npl
                        !        write(*,*) nint(plocation_dble(:,j))
                        !    end do
                        !    write(*,*)
                        !    write(*,*)
                        !end if

                

                    end do
            

        
            !call cpu_time(finish)
            !write(*,*) (finish-start)/1000

        
        
        end if
        
    end if
        
        
	

	

        end subroutine
  

    end module

