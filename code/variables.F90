module variables
    
    use LAS1D
    use sim3de 
    use piecewise_CMD
    
    implicit none
    
        !correlation matrix information for soil layers
		real, allocatable :: C0(:),CC(:,:,:), CE(:,:,:), CS(:,:,:), CI(:,:)
		real, allocatable :: AC(:,:,:,:), AE(:,:,:,:), AS(:,:,:,:), AI(:,:,:)
		real, allocatable :: ATC(:,:,:), ATS(:,:,:), ATI(:,:)
		real, allocatable :: CTC(:,:), CTS(:,:), CTI(:)
		integer :: M, k1, k2, k3, kk
        
        !same as above, but for the 2D layer boundary LAS profiles
		real, allocatable :: bC0(:)
    	real, allocatable :: bCT(:,:), bCC(:,:,:), bCS(:,:,:), bCI(:,:)
    	real, allocatable :: bAT(:,:,:), bAC(:,:,:,:), bAS(:,:,:,:), bAI(:,:,:)
		integer :: bM, bk1, bk2, bkk
        
        !same as above, but for the 1D CMD random field generator
		real, allocatable :: C01D(:)
        
        !same as above, but for the piecewise 3D CMD random field generator
        real(4), allocatable :: R0x(:), R0y(:), R0z(:)
        
        ! Pointer to desired 3D random field subroutine
        procedure(no_arg_sub), pointer :: RF3D => NULL()
        real soil_dummy_var ! a dummy argument is needed when calling the soil subroutine though a procedural pointer due to a bug in intel compiler version 16
        
        interface
            subroutine no_arg_sub(dummy)
            real, intent(in) :: dummy
            end subroutine
        end interface
        
        
        
        character(2), parameter :: rednames(6) = (/ 'SA','GA','HA','1Q','SD','MN' /) !reduction method names; hard-coded in SI module
        
        !whether to generate a very large soil in RAM, and take random subsets, as opposed to generating multiple independant soils
        !This will also pre-process and store the layer boundaries in RAM, and pre-process the CK layer depths at piles.
            !There is no performance penalty for setting this to false in the multi-layer mode unless using the genetic algorithm.
        logical :: superset   
        

        real, allocatable :: lmeans(:,:),ldepths(:) !mean and average boundary of each layer
        real, allocatable :: lmean_ave(:)
        real, allocatable :: lmean_ln(:),lsd_ln(:) !lognormal statistics for each layer for generating random values
        real(8) bsof !boundary scale of fluctuation (m)
        real(8) bsd  !boundary standard deviation (elements)
        integer nlayer !number of layers
        logical uselcoords !use layer coordinates from another file
        logical enforcelcoords !ensure that the layer boundaries at each input borehole are at their correct location after having random noise added
        
        integer :: nxe,nye,nze,nxew,nyew,nzew !dimensions of original aN**b single-layer LAS profile, dimensions of final Working field
        integer :: zroom !extra space needed in LAS for storing previous stages of the field
        integer :: MXM,MXK !max number of subdivisions, max dimension of stage-0 profile
        
        integer :: anisotropy !degree of anisotropy (hardcoded to be constant for all layers for now) reduces vertical SOF by this factor
        character(1) distribution !distribution for soil. 'n' = normal, 'l' = lognormal (recommended), 'b' = beta
	    logical :: singletrue !true if using a single, variable layer. False if using multi-layered soils with uniform properties
        
        !first two slots store the normal/expected mean and standard deviation, where the last two store the lognormal parameters
	    !for the beta distribution, the 4 slots correspond to the 4 beta distribution parameters
	    real(4) :: sdata(1,4) !soil statistics 
        
        real(8) :: dx,dy,dz !size of elements (m)
		character(6) :: varfnc,bvarfnc !correlation function to use for soil
        
	    real(8) :: soilth(2) !horizontal and vertical SOF
        
        real(4), allocatable :: efld(:,:,:) !array for storing the output of a single layer LAS field by the sim3de algorithm
        real(4), allocatable :: bfld(:,:,:) !layer boundary depths for each layer
        real(4), allocatable :: fullbfld(:,:,:,:) !storing all intermediate layer boundaries for all monte carlo realisations
        
        integer :: save_effE !output effective soil modulus of each investigation in each MC realisation. Not needed for your project. 0 = No. 1 = save 2= load.
        real :: emean !soil mean stiffness
    
        integer :: soil_reps(2) !range of Monte Carlo realisations from which to export the virtual soils as a CSV
        integer :: nodisksave !Nodisksave: 0 = keep in ram (much slower, but can save a fair amount of disk space), 1=write and read soil, 2=read soil only
        integer kseed !random seed
        
        !Data for when user wants to set layer boundaries based on existing borehole information from real sites
        integer :: nbhmulti                            !number of input boreholes with layer information
        real(8), allocatable :: bhxymulti(:,:)                  !Coordinates of the above boreholes
        real(8), allocatable :: bhdepthsmulti(:,:)         !layer information from the above boreholes
        
        
        
        
        
        !pile-related variables 
        
        real :: buildingarea !building area (m**2)
        integer multitype !whether the layer depths at piles are treated as (1) a single point in SI and CK, (2) inverse weighted in CK, (3) inverse weighted in SI and CK
		
        
        !True: randomize each uniform layer stiffness in the true soil as opposed to having a single constant value across all Monte Carlo realisations
        !False: apply the randomness as a white noise field for the borehole samples
        logical :: rand_realisations = .true. 
        
        
        
    
    
    contains

    

    !Get the sizes of the arrays for the soil correlation matrix, allocate them and calculate the required correlations.
    !This prepares the 3D random field generators in the single layer mode. Alternatively, it prepares the 2D and 1D random
    !   field generators for the multi layer mode.
    !It also uses different 3D random field generators depending on whether the soil is isotropic.
	subroutine soilsizes(istat)

		implicit none

		integer, intent(in) :: istat
		integer j1,j2,j3,u1,u2,u3,ndiv,msize,j,i
        
        real(4), allocatable :: soiltemp(:,:,:)
		
        real testrand(nzew)
        integer input
        real dummy
		
        integer nlayer
        
        logical is_isotropic ! true if the soil is isotropic (vertical SOF = horizontal SOF)

        zroom = 5*nze/4


		1  format(a,a,a)
		2  format(a,e13.4)
3          format(a,i4,a,i4,a,i4,a,i4,a)
           
        ! The soil is considered isotropic if the SOF in the horizontal and vertical directions is within a 1 mm tolerance (conservative tolerance)
        is_isotropic = abs(soilth(1) - soilth(2)) < 0.001
           
           
      
           
        if (singletrue) then
            
            if(is_isotropic) then
                !Use LAS if it's isotropic

                u1 = nxe
                u2 = nye
                u3 = nze

                do ndiv = 0, MXM
                    msize = u1*u2*u3
                    if( msize .le. MXK ) exit
                    j1 = u1/2
                    j2 = u2/2
                    j3 = u3/2
                    if( 2*j1 .ne. u1 .or. 2*j2 .ne. u2 .or. 2*j3 .ne. u3 ) exit
                    u1 = j1
                    u2 = j2
                    u3 = j3
                end do

                if(  msize .gt. MXK ) then
                    write(istat,1)'Error: unable to determine an acceptable combination of k1, k2, k3 and m'
                    write(istat,1)'       such that k1*2**m = N1, k2*2**m = N2 and k3*2**m = N3 for 3D LAS field.'
                    write(istat,3)'       k1 = ',u1,', k2 = ',u2,', k3 = ',u3,', m = ',ndiv
                    write(istat,3)'       (k1*k2*k3 must be less than ',MXK,' and m must be less than ',MXM,')'
                    write(istat,1)'       Try changing N1, N2, and/or N3.'
                    read(*,*)
                    stop
                end if

                allocate(C0(msize*(msize + 1)/2),CC(28,8,ndiv), CE(28,12,ndiv), CS(28,6,ndiv), CI(28,ndiv) &
                ,AC(8,7,8,ndiv), AE(12,7,12,ndiv), AS(18,7,6,ndiv), AI(27,7,ndiv) &
                ,ATC(4,7,4), ATS(6,7,4), ATI(9,7) &
                ,CTC(28,4), CTS(28,4), CTI(28))
                
                M = ndiv
                k1=u1
                k2=u2
                k3=u3
                kk=msize

                !get soil correlation arrays
                call sim3d_init(nxe,nye,nze,dz,dz,dz,soilth(1),soilth(1),soilth(1),varfnc,MXM,MXK, &
                C0(:),CC(:,:,:),CE(:,:,:),CS(:,:,:),CI(:,:),AC(:,:,:,:),AE(:,:,:,:),AS(:,:,:,:),AI(:,:,:), &
                ATC(:,:,:), ATS(:,:,:), ATI(:,:),CTC(:,:),CTS(:,:),CTI(:),M,k1, k2, k3, kk)

                if(superset) then
                    RF3D => sim3dw_nosubset ! don't take a random subset, as this is done later
                else
                    RF3D => sim3dw  ! return a field of the exact working size
                end if

            
            else ! do anisotropic soils
                !use the piecewise covariance matrix decomposition method
                
                if(superset) then !If storing a single big soil in memory...
                    allocate(R0x(nxe*nxe), R0y(nye*nye), R0z(nze*nze))  ! allocate a larger soil
                    RF3D => cmd_pw_loop     ! this loop implementation is typically faster for large soils, so use it here
                    ! get correlation arrays
                    call piecewise_init(nxe,nye,nze,dz,soilth(1),soilth(2),bvarfnc,R0x,R0y,R0z)
                else
                    allocate(R0x(nxew*nxew), R0y(nyew*nyew), R0z(nzew*nzew))
                    RF3D => cmd_pw_matmul   
                    ! get correlation arrays
                    call piecewise_init(nxew,nyew,nzew,dz,soilth(1),soilth(2),bvarfnc,R0x,R0y,R0z)
                end if

            end if

        else
            
        !prepare random field generator for the 2D layer boundary LAS profiles

		
            u1 = nxe
            u2 = nye

            do ndiv = 0, MXM
                msize = u1*u2
                if( msize .le. MXK ) exit
                j1 = u1/2
                j2 = u2/2
                if( 2*j1 .ne. u1 .or. 2*j2 .ne. u2) exit
                u1 = j1
                u2 = j2
            end do


            if(  msize .gt. MXK ) then
                write(istat,1)'Error: unable to determine an acceptable combination of k1, k2 and m'
                write(istat,1)'       such that k1*2**m = N1 and k2*2**m = N2 for 2D LAS field.'
                write(istat,3)'       k1 = ',u1,', k2 = ',u2,', m = ',ndiv
                write(istat,3)'       (k1*k2 must be less than ',MXK,' and m must be less than ',MXM,')'
                write(istat,1)'       Try changing N1 and/or N2.'
                stop
            end if
            
            allocate(bC0((msize*(msize + 1))/2),bCT(6,2), bCC(6,4,ndiv), bCS(6,4,ndiv), bCI(6,ndiv) &
            ,bAT(3,3,2), bAC(4,3,4,ndiv), bAS(6,3,4,ndiv), bAI(9,3,ndiv))
            
            bM = ndiv
            bk1=u1
            bk2=u2
            bkk=msize

                call sim2sd_init(nxe,nye,dz,dz,bsof,bsof,bvarfnc,MXM,MXK, &
                bC0(:),bCT(:,:),bCC(:,:,:),bCS(:,:,:),bCI(:,:),bAT(:,:,:),bAC(:,:,:,:),bAS(:,:,:,:),bAI(:,:,:), &
                bM,bk1,bk2,bkk) !boundary random noise generation


        !Prepare 1D random field generator, which will be used in the multi-layer hybrid mode

            !setup 1D random field generator; covariance matrix decomposition method
            allocate(C01D((nxe*(nxe + 1))/2))
            
            call sim1sd_init(nzew,dz,soilth(1),bvarfnc,C01D)
            
            !test it (seems to work)
            call las1g( testrand, nzew ,C01D)

        end if
        
        
    
        


    end subroutine
    
    
    !wrapper for the sim3d subroutine that generates virtual soils
    subroutine sim3dw(dummy)
    
    real, intent(in) :: dummy
    call sim3d(efld,nxe,nye,nze,zroom,nxew,nyew,nzew,dx,dy,dz,kseed,MXM,MXK, &
                      C0,CC,CE,CS,CI,AC,AE,AS,AI,ATC, ATS, ATI,CTC,CTS,CTI,M,k1,k2,k3,kk,sdata,distribution,anisotropy)
    
    end subroutine
    
    
    !wrapper for the sim3d subroutine that generates virtual soils
    subroutine sim3dw_nosubset(dummy)
    
    real, intent(in) :: dummy
    call sim3d_nosubset(efld,nxe,nye,nze,zroom,nxew,nyew,nzew,dx,dy,dz,kseed,MXM,MXK, &
                      C0,CC,CE,CS,CI,AC,AE,AS,AI,ATC, ATS, ATI,CTC,CTS,CTI,M,k1,k2,k3,kk,sdata,distribution,anisotropy)
    
    end subroutine
    
    
    !wrapper for the sim3d_init subroutine that initialises the random soil generator
    subroutine sim3d_initw()
    
    call sim3d_init(nxe,nye,nze,dx,dy,dz,soilth(1),soilth(1),soilth(1),varfnc,MXM,MXK, &
                      C0,CC,CE,CS,CI,AC,AE,AS,AI,ATC, ATS, ATI,CTC,CTS,CTI,M,k1,k2,k3,kk)
    
    end subroutine
    
            
    
    subroutine cmd_pw_loop(dummy)
    
    real, intent(in) :: dummy
		
    call cmd_piecewise_loop(efld, nxe,nye,nze ,R0x, R0y, R0z,sdata, distribution)
    
    end subroutine
        
    
    
    subroutine cmd_pw_matmul(dummy)
    
    real, intent(in) :: dummy
    
    call cmd_piecewise_matmul(efld, nxew,nyew,nzew ,R0x, R0y, R0z,sdata, distribution)
    
    
    end subroutine
    
    
    
	

    end module

    
