module detcheck
    
    
    use SI
    use fem_prep
    use ESETT
    use soilgen
    use variables
    
implicit none



    contains
   
    
    
       subroutine det_check(datafolder,pdepths,npdepths,num_loads,abstol,difftol,buildingweight,preps,rel_loads,load_con,in_sdev,in_percentile,prad,femrad,goodpiles,goodcases,plocation) 
       
       
       character(100), intent(in) :: datafolder !directory of data
       integer, intent(in) :: npdepths !number of pile depths in single layer mode
       integer, intent(in) :: num_loads
       real, intent(in) :: abstol,difftol,buildingweight
       integer, intent(in) :: preps(2)
       real, intent(in) :: rel_loads(:)
       integer, intent(in) :: load_con(:)
       real, intent(in) :: in_sdev,in_percentile
       integer, intent(in) :: prad(2) !pile radius in elements
       integer, intent(in) :: femrad !mesh dimensions 
       integer, intent(in) :: goodcases !number of valid piles
       integer, intent(in) :: goodpiles(:) !index of good piles
       integer, intent(in) :: plocation(:,:) !pile locations in elements
       integer, intent(in) :: pdepths(npdepths)			!pile depths to analyse (elements)
       
       
       	integer, parameter :: curvefactor=10 !factor by which to increase settlement curve resolution
        integer, parameter :: numsamples = 10000 !number of random samples to generate
        real, parameter :: vval = 0.3 !poisson's ratio set to 0.3
        character(2), parameter :: rednames(5) = (/ 'SA','GA','HA','1Q','SD' /) !reduction method names; hard-coded in SI module
        
	    real detsetsi((npdepths-1)*curvefactor+1)
	    integer pdepthssi((npdepths-1)*curvefactor+1)
        !real pdepthssi((npdepths-1)*curvefactor+1)
        integer :: pdepths_ind(npdepths)
        integer :: sides(5,1,num_loads)
        real :: detdisp(npdepths)
        
        real :: samples(numsamples)
        real    lgoodvals(numsamples)                         !log values of the above info
        logical gvmask(numsamples)                            !mask invalid values of the above info
        
        
        
        !local variables
        integer i,j,k,counter !loop counters
        real temp
        real eval(1,1,1,1) !reduced effective Young's modulus
        character(1000) :: str2
        logical exists 
        real efld2d(femrad,nzew) !soil field for the 2D FEA mesh
        real :: femvech(femrad), femvecv(nzew)
        real :: dist((preps(1)*preps(2))**2) !distances between piles
        real :: settol
        real :: realsides(5,maxval(load_con))
        integer :: load_conSI(num_loads)
        real(8) :: bsdtemp
        
        !--multi layer stuff--

        real(8) :: xyi(2,nxew*nyew) !coordinates of points to interpolate
        real multisides(goodcases)
        real lap(nlayer-1,goodcases)
        
                    !expanded relative load information
        integer load_con2(preps(1)*preps(2))
        real rel_loads2(preps(1)*preps(2))
        
        
        
        
        
        
        open(2547,file='deterministic_report.txt')
        write(2547,'(A)') 'This file contains approximate average pile designs for each load case and all reduction methods for the single layer soil input.'
        write(2547,'(A)') 'It also shows the average pile design for all piles for the multi-layer soil input.'
        write(2547,'(A)') 'These results indicate whether the current settings produce feasible pile designs prior to using SIOPS.'
        write(2547,'(A)') 'Note: Negative designs are invalid. For example:'
        write(2547,'(A)') '-100 means the required pile length is shorter than allowed (0 m). Increase soil stiffness or reduce pile load.'
        write(2547,'(A)') '-200 means the required pile length is longer than allowed. Decrease soil stiffness or reduce pile load.'
        write(2547,*)
        write(2547,*)
        write(2547,*) 'Single layer pile designs (m):'
        write(2547,*) '--------------------------------'
        write(2547,'(A,X,F6.1)') 'Maximum pile length:',pdepths(npdepths)*dz
        
        
      !#################### SINGLE LAYER ASSESSMENT #########################
        !Get the average pile lengths for each load case and reduction method

    
            !------get fake virtual soil by generating white noise -----
            !generate white noise (zero mean, unit variance)
            call vnorm( samples, numsamples )
            
            !transform into desired statistics
            if(distribution == 'n') then
				samples = sdata(1,1) + samples * sdata(1,2)	
		    else if(distribution == 'l') then
				samples = exp(sdata(1,3) + samples * sdata(1,4))
		    else if(distribution == 'b') then
				samples = sdata(1,1) + 0.5*(sdata(1,2)-sdata(1,1))*(1.0 + tanh((sdata(1,3)+sdata(1,4)*samples)/6.2831853))
		    else 
				write(*,*) 'Incorrect distribution selected. n = normal, l = lognormal, b = bounded.'
            end if
            if(singletrue) then
                allocate(bfld(nlayer+1,nxew,nyew))
            else
                samples = samples * emean
            end if

            
            !-----Get pile settlement curve values----
            write(str2,'(A,A,I0,A,F4.2,A)') trim(datafolder),'settlement_prad-',prad(1),'_esize-',dz,'.txt'
            inquire(file=str2,exist=exists) !check if the file exists
            if (exists) then    
                open(667, file= str2,status='old')
	            read(667,*) 
	            do i = 1,npdepths
		            read(667,*) temp,detdisp(i)
	            end do
	            close(667)
            else
                efld2d = 1
                femvech = dz
                femvecv = dz
                do i = 1,npdepths
                    call prep_fem2d(efld2d,vval,femvech,femvecv,(prad(1)+prad(2))/2,detdisp(i),pdepths(i),dz)
                end do
                
            end if
            
            !calculate distances between each pile here instead of in the heavily nested loop below
	        counter = 0
	        do i = 1,goodcases-1
		        do j = i+1,goodcases
			        counter = counter + 1
			        dist(counter) = sqrt(real((plocation(goodpiles(i),1)-plocation(goodpiles(j),1))**2 + (plocation(goodpiles(i),2)-plocation(goodpiles(j),2))**2))
		        end do
            end do
        
            !get absolute settlement design tolerance based on minimum pile spacing and differential settlement limit, or directly use abstol value if a positive value is given.
            if (abstol > 0) then
                settol = abstol
            else
                settol = ceiling(minval(dist(:counter))*dz*difftol*1000) !round settlement tolerance up to nearest millimetre
            end if 
    
    
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
    
           !loop through the reduction methods
            do i = 1,5
                call reduce_single(eval(1,1,1,1),i,numsamples,samples,gvmask,lgoodvals,in_percentile,in_sdev)
                call despile_discrete(eval, rel_loads, load_con, settol, buildingweight, detsetsi,sides(i:i,:,:), 1, preps(1)*preps(2),size(pdepthssi),1,1,num_loads,1)
                counter = 0
                do j = 1,num_loads
                    if(count(load_con==j) > 0) then
                        counter = counter + 1
                        realsides(i,counter) = sides(i,1,j)
                        load_conSI(counter) = j
                    end if
                end do
                
            end do
            
            where(realsides < 0) realsides = realsides * curvefactor
        
            write(2547,'(A15,10000(X,F6.3))') 'Pile load case:',rel_loads(load_conSI(:counter))
            write(2547,*)
            write(2547,'(A15,10000(X,F6.1))') 'Pile design: '//rednames(1),realsides(1,:counter)/curvefactor
            do i = 2,5
                write(2547,'(A15,10000(X,F6.1))') rednames(i),realsides(i,:counter)/curvefactor
            end do
            
            write(2547,*)
            if(minval(realsides) < 0 .or. maxval(realsides) > pdepths(npdepths)) then
                write(2547,*) 'Inputs INVALID'
            else
                write(2547,*) 'Inputs VALID'
            end if
            
            
            
            !#################### MULTIPLE LAYER ASSESSMENT #########################
            
            write(2547,*)
            write(2547,*)
            write(2547,*) 'Multiple layer pile designs (m):'
            write(2547,*) '--------------------------------'
            write(2547,'(A,X,F6.1,X,A)') 'Maximum pile length:',nzew*dz,'(soil depth)'
            
            !expand the relative loads for the good piles to work with the design subroutine
            do i = 1,goodcases
                rel_loads2(i) = rel_loads(load_con(goodpiles(i)))
                load_con2(i) = i
            end do
            
            !get indexing for layer interpolation
            counter = 0
            do i=1,nyew
                do j=1,nxew
                    counter = counter+1
                    xyi(1,counter) = j
                    xyi(2,counter) = i
                end do
            end do
            
            !generate layer boundaries
            bsdtemp = bsd
            bsd = 0
            call soil_layers(xyi)	
            bsd = bsdtemp
            
            !get layer boundaries at pile locations. Just use the average height within the pile here.
            do j = 1,goodcases
                do i = 1,nlayer-1    
	                lap(i,j) = dz * sum(bfld(i+1,plocation(goodpiles(j),1):plocation(goodpiles(j),1)+prad(1)-1,plocation(goodpiles(j),2):plocation(goodpiles(j),2)+prad(2)-1))/product(prad)
                end do
            end do
            
            write(*,*) lmean_ave
            
            !design piles
            if(singletrue) then !acccount for scaled building weight for single layer mode
                call multides_1D_det(prad,goodcases,nlayer,lap,lmean_ave,multisides,rel_loads2(:goodcases),load_con2(:goodcases), settol, buildingweight*emean,dz,nzew)
            else
                call multides_1D_det(prad,goodcases,nlayer,lap,lmean_ave,multisides,rel_loads2(:goodcases),load_con2(:goodcases), settol, buildingweight,dz,nzew)
                call multisetcurve(prad,goodcases,nlayer,lap,lmean_ave,rel_loads2(:goodcases),load_con2(:goodcases), settol, buildingweight,dz,nzew)
            end if
            
            !output results
            write(2547,'(A15,10000(X,F6.3))') 'Pile load case:',rel_loads2(:goodcases)
            
            
            write(2547,'(A15,10000(X,F6.1))') 'X coords:',plocation(goodpiles(:goodcases),1)*dz
            write(2547,'(A15,10000(X,F6.1))') 'Y coords:',plocation(goodpiles(:goodcases),2)*dz
            write(2547,*)
            write(2547,'(A15,10000(X,F6.1))') 'Pile design:',multisides
            
            write(2547,*)
            if(minval(multisides) < 0 .or. maxval(multisides) > nzew) then
                write(2547,*) 'Inputs INVALID'
            else
                write(2547,*) 'Inputs VALID'
            end if

            
            
            close(2547)
            if(singletrue) deallocate(bfld)
            
            
    
    end subroutine
    
    
end module