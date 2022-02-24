module writesoils

    
    use variables
    use soilgen
    
contains


!This subroutine writes the soils to disk so that they may be loaded in at a later time (deprecated)
!This greatly increases the speed of the process in the Evolutionary algorithm, where the 
!same soils are used over and over

subroutine storesoils(soilseeds,nrep_MC)

implicit none


      !real(4) :: efldave(nxew,nyew,nzew) !virtual soil averaged across the MC realisations
      integer, intent(in) :: soilseeds(nrep_MC)

      
      integer, intent(in) :: nrep_MC
      
      !local variables
      character(1000) str2 !string for output
      integer tempsof, tempcov, tempnrep
      integer cov
      integer kseed
      integer i
      logical exists
      real randu
      logical generate
      real(4) :: sdata_ck(4) !CK soil parameters, same as SI but with unit mean stiffness

      real(4) :: pvr !lognormal variable

      
      cov = nint(100*sdata(1,2)/sdata(1,1))
      generate = .false.
      
        if(nodisksave >= nrep_MC .and. singletrue .and. superset==1) then !Read the soils into disk later. Make sure the soils have the correct statistics.
        !f the files don't exist, or are the wrong statistics, then regenerate them. 
            write(str2,'(A,A)') trim(datafolder),'soilstats.txt'
            inquire(file=str2,exist=exists) !check if the file exists 
        	if(exists) then
				open(500,file=str2,status='old')
				read(500,*) tempsof, tempcov, tempnrep
				close(500)
				if(tempsof /= nint(soilth(1)) .or. tempcov /= cov .or. tempnrep < nrep_MC) then
					generate = .true.
					write(*,'(A,X,I0,X,A,X,I0,X,A,X,I0,A,X,I0,A)') 'Warning, the soils',tempnrep,'saved to disk have SOF',tempsof,'m and COV',tempcov,'%'
					write(*,'(A,X,I0,X,A,X,I0,A,X,I0,A,X,I0)') 'The current inputs are SOF',nint(soilth(1)),'m, COV',cov,'% and No. MC realisations',nrep_MC
					write(*,*) 'Generating new soils with current inputs, and saving to disk, over-writing old soils.'     	
				end if 
			else
				generate = .true.
				write(*,*) 'Warning, no soils are saved to disk.'
				write(*,*) 'Generating new soils with current inputs, and saving to disk.'
			end if
			
		else 
			write(*,*) "Generating soils within MC realisations. No disk space required."			
		end if
			
        if (generate) then !write soils to file. Required disk space = 4*x*y*z*(No. MC realisations)*(element size)**-3. X,y,z are the dimensions of site in m.
        
        
            !convert and use unit mean (1 MPa) soil stiffness for CK and scale at a later time.
            ! NOTE: Now the soil is generated as unit mean, and the load is scaled instead
			!sdata_ck(1) = 1.0
			!sdata_ck(2) = sdata(2)/sdata(1)
			!pvr = log(1.0 + sdata_ck(2) ** 2 / sdata_ck(1) ** 2)
			!sdata_ck(3) = log(sdata_ck(1)) - 0.5 * pvr
			!sdata_ck(4) = sqrt(pvr)
			
			
			!get the average virtual soil (see avesoil subroutine for details)
			!call avesoil(soilseeds,nrep_MC,datafolder,efldave)
        
            write(*,'(A,X,I0,X,A,X,F6.2,A)') 'Saving', nrep_MC,'virtual soils to disk, needing', real(nrep_MC)*4*nxew*nyew*nzew/1000000000, 'GB of free space.'
            !Generate the virtual soils once and save them to file, then load them as needed in the EA
            write(str2,'(A,A)') trim(datafolder),'soil.dat'
            open(500,file=str2,access='stream') !,buffered='yes')
			do i = 1,nrep_MC     !Loop through Monte Carlo realisations
			
				!Progress indicator
				if(mod(i,nrep_MC/10)==0) write(*,'(I0,A,X)',advance='no') i*100/nrep_MC,'%'
			
                kseed = randu(soilseeds(i)) * 1234567890    !ensure random numbers are consistent across MC realisations
				call RF3D(soil_dummy_var)
				!write(str2,'(A,I0,A)') 'soildata/soil',i,'.dat'  
				
				!save soil to file in a compact and fast binary format (this format might be machine or compiler-specific..). Subtract to get
				!efld = efld-efldave
				write(500) efld     			!Subtract the average, for reasons given above.
				!write(*,*) minval(efld),maxval(efld),sum(efld)/size(efld)
				!write(*,'(X,I0)',advance='no') i
            end do
            close(500)
            write(*,*)
            write(str2,'(A,A)') trim(datafolder),'soilstats.txt'
            open(500,file=str2) !Save the mean and SOF of the soil to a file to ensure it's the same as the one loaded in.
            write(500,*) nint(soilth), cov, nrep_MC
            close(500)
        end if


end subroutine



subroutine avesoil(soilseeds,nrep_MC,efldave)

	!calculate the average soil across Monte Carlo realisations

	!There is a very, very small small correlation in the soils across Monte Carlo realisations such that the average
	!is not exactly 1.0, as it theoretically should be. The average is itself a correlated random field with weak and strong regions.
	!Although this is negligible in most applications (1% COV in the average when a COV of 80% is specified), it can mess around with
	!the genetic algorithm a bit, as the borehole gravitate towards the average weak zones.
	
	!To combat this, the average field is calculated first, and then subtracted from each individual random field. 
	!This should bring the average of each soil across MC realisations to a constant value of 1.0 everywhere.
	
	implicit none
	

      real(8) :: efldavedb(nxew,nyew,nzew) !virtual soil averaged across the MC realisations, double precision
      
      integer, intent(in) :: soilseeds(nrep_MC)
      
      integer, intent(in) :: nrep_MC
      
      real(4), intent(out) :: efldave(nxew,nyew,nzew) !virtual soil averaged across the MC realisations
      
      !local variables
      character(1000) str2 !string for output
      integer kseed
      real newmin
      integer i,x,y,z
      logical exists
      real randu
      real(4) :: sdata_ck(4) !CK soil parameters, same as SI but with unit mean stiffness

      real(4) :: pvr !lognormal variable

	
	
	write(str2,'(A,A,I0,A,I0,A,I0,A,I0,A)') trim(datafolder),'soilave_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_nMC-',nrep_MC,'_anis-',anisotropy,'.dat' 
	inquire(file=str2,exist=exists) !check if the file exists 
	if(exists) then
	
		!file already exists, read it in
		open(500,file=str2,access='stream')
		read(500) efldave
		close(500)
		
	else
		write(*,*) 'Calculating averaged virtual soil.'
		!generate the average soil file
		efldavedb = 0
		do i = 1,nrep_MC     !Loop through Monte Carlo realisations
		
			!Progress indicator
			if(mod(i,nrep_MC/10)==0) write(*,'(I0,A,X)',advance='no') i*100/nrep_MC,'%'
			
			kseed = randu(soilseeds(i)) * 1234567890    !ensure random numbers are consistent across MC realisations
			call RF3D(soil_dummy_var)
				 
			efldavedb = efldavedb + efld !Do the sum part of the average.

		end do
		write(*,*)
		efldave = efldavedb/nrep_MC !Divide to get the actual average.
		newmin = maxval(efldave)
		!Give the array an average of zero to turn it into bias, 
        !or at least as close to zero as possible while ensuring all values remain positive (hence the maxval)
		efldave = efldave - newmin !sum(efldave)/size(efldave) 
		
		write(*,*) 'Mean increased by',1-newmin

	
		!Save the average soil to file for later

		open(500,file=str2,access='stream')
		write(500) efldave !save as single precision
		close(500)
        
	
	end if
	
end subroutine


subroutine exportsoil(soilseeds,nrep_MC)
						 

	!This subroutine exports soils to .
	
	implicit none
	
	  integer NGS

      integer, intent(in) :: soilseeds(nrep_MC)

      !character(1000),intent(in) :: datafolder !a string representing the directory the data is stored in
      
      integer, intent(in) :: nrep_MC

      
      !local variables
      character(1000) str2 !string for output
      integer i,x,y,z,n,j,xpos,ypos,zpos
      logical exists
      real randu
      real(4) :: sdata_ck(4) !CK soil parameters, same as SI but with unit mean stiffness
      integer b,l,counter
      real(8), allocatable :: xyi(:,:)
      real efld2D(nlayer,nxew,nyew) !Young's modulus in each layer as a 2D random field
      real sdata_temp(4)

      real(4) :: pvr !lognormal variable

      !only save soils if the first number is greater than 0
	  if (soil_reps(1) > 0) then
	  
	  	write(*,'(A,I0,A,I0)') 'Exporting soils to CSVs. Doing range ',soil_reps(1),' - ',abs(soil_reps(2))
        
        
        if(.not. singletrue) then
                          !get index values
            counter=0
            allocate(xyi(2,nxew*nyew))
            allocate(efld(nxew,nyew,nzew))
        
            do i=1,nyew
                do j=1,nxew
                    counter = counter+1
                    xyi(1,counter) = j
                    xyi(2,counter) = i
                end do
            end do
        end if
        
        !generate big single layer here if specified
        if(singletrue .and. superset > 1) call RF3D(soil_dummy_var)
        
        
		do i = soil_reps(1),abs(soil_reps(2))     !Loop through Monte Carlo realisations

			kseed = randu(soilseeds(i)) * 1234567890    !ensure random numbers are consistent across MC realisations

            
            if (singletrue) then !output single layer

                if (superset > 1) then !generate subsets from big soils
                   xpos=NINT(randu(0)*(nxe-nxew))
		           ypos=NINT(randu(0)*(nye-nyew))
		           zpos=NINT(randu(0)*(nze-nzew))
                   
                   call savetofile(efld(xpos+1:xpos+nxew,ypos+1:ypos+nyew,zpos+1:zpos+nzew),str2,i)
                else !otherwise generate individual soils
                   call RF3D(soil_dummy_var)       !call single layer soil
                   !efld = efld * emean
                   call savetofile(efld,str2,i)
                end if
            else !save multiple layer soils, both the 3D one (need to construct manually) and boundary layers
                !Young's modulus in each layer optionally represented by a 2D random field if standard deviation is greater than zero. Generate fields.
                do n = 1,nlayer
                    sdata_temp(3) = lmean_ln(n)
                    sdata_temp(4) = lsd_ln(n)
                    ! generate logarithm of 2D random field for Young's modulus (assumes that it's lognormally-distributed)
                    call sim2d(efld2D(n,:,:),nxe,nye,nxew,nyew,5*nye/4,dz,dz,kseed,MXM,MXK, &
                    bC0(2,:),bCT(2,:,:),bCC(2,:,:,:),bCS(2,:,:,:),bCI(2,:,:),bAT(2,:,:,:),bAC(2,:,:,:,:),bAS(2,:,:,:,:),bAI(2,:,:,:), bM, bk1, bk2, bkk,sdata_temp,distribution)
                end do

                call soil_layers(xyi)   !get soil layer boundaries
                do n=nlayer,1,-1  !Progressively work forwards through time, adding newer layer
                    write(str2,'(A,I0,A,I0,A)') 'bound',i,'_layer',n,'.txt'
	                open(200,file=str2)
                    do y=1,nyew !fill in y direction
                        do x=1,nxew     !fill in x direction
                            if (.false.) then
                                efld(x,y,:nint(bfld(n+1,x,y))) = exp(efld2D(n,x,y))      
                            else
                                call las1g( efld(x,y,:nint(bfld(n+1,x,y))), nint(bfld(n+1,x,y)), C01D) !generate zero-mean, unit variance random noise
                                efld(x,y,:nint(bfld(n+1,x,y))) = exp(efld2D(n,x,y) + efld(x,y,:nint(bfld(n+1,x,y))) * lsd_ln(n))
                            end if
                        end do
                        if (n < nlayer) write(200,'(1000000000I4)') nint(bfld(n+1,:,y)) !output current row to file
                    end do
                    close(200)
                end do
                
                call savetofile(efld,str2,i)
                
            end if 
					 

            
				 

        end do
        
        if( .not. singletrue) deallocate(efld)
        
        if(soil_reps(2) < 0) stop

	
	  end if



end subroutine


subroutine savetofile(efld,str2,i)
!save a single realisation of the 3D soil to file

implicit none

real, intent(in) :: efld(:,:,:)
integer, intent(in) :: i    !realisation
character(1000) :: str2
integer x,y,z


	!Save the file in a long row of fixed-width values. X is varying first, then Y, then Z the slowest.
	write(str2,'(A,I0,A)') 'soil',i,'.txt'
	open(200,file=str2)
    write(200,'(I0,A)') nxew,' X dimension (varying fastest)'
    write(200,'(I0,A)') nyew,' Y dimension'
    write(200,'(I0,A)') nzew,' Z dimension (varying slowest)'
	!write(200,'(1000000000(E14.7,A))') ((((efld(x,y,z),' '),x=1,nxew),y=1,nyew),z=1,nzew)
    do z = 1,nzew
        do y = 1,nyew
            do x=1,nxew
                write(200,*) efld(x,y,z)
            end do
        end do 
    end do
            
	close(200)


end subroutine


subroutine savetofile_binary(efld,str2,i)
!save a single realisation of the 3D soil to file in a binary format (faster and much less storage space required)

implicit none

real, intent(in) :: efld(:,:,:)
integer, intent(in) :: i    !realisation
character(1000) :: str2


	!Save the file in a long row of fixed-width values. X is varying first, then Y, then Z the slowest.
	write(str2,'(A,I0,A)') 'soil',i,'.dat'
	open(200,file=str2,access='stream')
	!write(200,'(1000000000(E14.7,A))') ((((efld(x,y,z),' '),x=1,nxew),y=1,nyew),z=1,nzew)

    write(200) efld
            
	close(200)


end subroutine

end module