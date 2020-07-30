!c  *********************************************************************
!c  *                                                                   *
!c  *                         subroutine ESETT                          *
!c  *                                                                   *
!c  *********************************************************************
!c  Single Precision Version 1.0
!c  Written by Michael P Crisp
!c  Latest Update: Jun 4, 2018
!c
!c  PURPOSE  Get elastic pile settlement in a 1D multi-layered soil with uniform soil properties.
!c
!c  DESCRIPTION: Pile is designed using the Mylonakis and Gazetas method: 
!   "Settlement and additional internal forces of grouped piles in layered soil", 1998.
    
!   This method treats a pile in a multi-layered soil as individual pile components.
!   The stiffness at the base of each component is taken as the overall stiffness of the
!   component directly below it. The method works from the bottom-most component upwards.

!   Some assumptions are that the pile tip is a rigid circular disk on the surface of a homogeneous elastic stratum.
!   The pile is also perfectly bonded to the surrounding soil.
!   The pile can be treated either as either entirely floating, with the below soil extending to
!       infinity, or can terminate at rigid bedrock (the latter is applied here, with bedrock
!       occuring at an arbitrary depth below the pile, specified by the depth of the soil).
!   The method is based on Randolph and Wroth (1998), in that the soil around a pile shaft is represented
!       by distributed springs (Winkler assumption).

!   A modificiation to the method is how base stiffness of the bottom component is calculated.
!   Traditionally, this is equal to the material properties of the layer the pile tip is resting
!   on or in. Here, the base stiffness is calculated as a weighted harmonic average of the soil
!   below the tip. The weights are determined by integrating the area of exponentially-decaying function
!   over the soil below the tip. This reflects the general assumption that stresses due to the pile are
!   also exponentially decaying with distance. A harmonic average is used as it is highly representative
!   of soils that are constant in the horizontal direction (see: Three-dimensional probabilistic 
!   foundation settlement, Fenton and Griffiths 2005).
!   This base stiffness modification has 2 advantages.
!   1. The pile is affected soil that is deeper in the ground, as is the case in reality.
!   2. It results in a continuous settlement function. Without it, there is a discontinunity across layers
!       as the base stiffness suddenly changes.
!
!   While it is difficult to determine what the best decay rate really is, preliminary testing has 
!   shown that a half-life of 3 m is reasonable for a range of pile widths.
!
!c
!c  PARAMETERS:
!c  
!c  Ap, Bp, Cp, des area, width, circumference and length of pile respectively
!c
!c  DL, Pall is the design load and allowable load (kN)
!c
!c  qclip is an intermediate parameter in the skin resistance calculation
!c
!c  qca is the clipped average
!c
!c  upper, lower are the upper and lower bounds for the qclip values in terms
!c  of array index
!c
!c  start, finish are the start and end depths for the base qc to be calculated
!c  
!c  Pb, Ps are the base and skin resistance of the pile
!c
!c  kc nad psi are the penetrometer bearing capacity factor and a 
!c  soil/construction adjustment parameter
!c
!c  qsi, qsimax are the limit unit skin friction and its maximum allowable value
!c
!c  a is 1.5xpile width, used for the qca calculation
!c
!c  REVISION HISTORY:
!c-------------------------------------------------------------------------
      
      module ESETT
      
      implicit none
      
    !these are local variables for the esett1D subroutine. 
    !I put them here in the hope that having them constantly in RAM as opposed to reinitialising them a billion times will speed things up a bit.
      real, allocatable :: effe(:),weights(:) !variables for weighted harmonic average of base stiffness
      real, allocatable :: thicknesses(:) !thickness of each layer according to the pile
      real thickness !thickness of current layer
      real Gave ! average shear modulus over pile
      real, parameter :: dr = log(2.0)/3 !exponential decay rate for a half-life of 3 m.
      integer penelayer !no. layers penetrated
      real Kb, omega, mu, delta, K, prevKb,theta !variables used in the settlement method
      integer i, j, n !loop counters
      real Eps !very small tolerance for layer thickness
      
      private :: effe,weights,thickness,dr,penelayer,Kb,omega,mu,delta,K,prevKb,theta,i,j,n,Eps
      
    contains
    
    
    
    !A single instance of the settlement calculation    
    function esett1D(E,Bp,DL,Ep,nlayer,des,aveloc,Ap,pi,G,v) result(settle) !get the settlement for a particular pile depth

    real, intent(in) :: pi, E(:) !E(nlayer)
    real, intent(in) ::  G(:) !G(nlayer)
    real, intent(in) :: Bp, Ep, Ap !pile diameter, stiffness, area
    
    integer, intent(in) ::  nlayer ! number of layers
    real,intent(in) :: DL !design load
    
    real :: settle ! settlement of final design
    real, intent(in) :: des
    
    real, intent(in) :: aveloc(:) !aveloc(nlayer+1) !location of layer boundaries in terms of depth
    
    
    
    real, intent(in) :: v !poisson's ratio
    
      
    
     do i = 1,nlayer !find out which layer the pile gets in to
		if(des >= aveloc(i) .and. des < aveloc(i+1)) then
		    penelayer = i
		    exit
		end if
     end do
     
     !do base stiffness at bottom of pile
    do n = nlayer,penelayer+1,-1      !get weights from below the layer the pile is founded in
        weights(n) = -exp(-dr*(aveloc(n+1)-des))/dr + exp(-dr*(aveloc(n)-des))/dr !area = upper integral - lower integral
        effe(n) = weights(n)/E(n) !weighted young's modulus for current component (inverted for harmonic average calculation)
    end do
    weights(penelayer) = -exp(-dr*(aveloc(n+1)-des))/dr + 1/dr !get weight from layer pile is founded in 
    effe(penelayer) = weights(penelayer)/E(penelayer)
    effe(1) = sum(weights(penelayer:))/sum(effe(penelayer:)) !store weighted harmonic average in first element
            
	Kb = ((Bp*effe(1))/(1-v**2)) !* (1+0.65*Bp/(aveloc(nlayer+1)-des)) !this (1 + 0.65...) bit accounts for rock at a certain depth
     
			  
	theta = 2*pi/log(5*des*(1-v)/Bp)
	do i = penelayer,1,-1
		!if(aveloc(i+1) - aveloc(i) <= Eps) cycle !if the current layer is a non-existant lense, skip
		thickness = min(des,aveloc(i+1)) - aveloc(i)
        !if(thickness < 1 .and. i == 1) then
        !    write(*,*)
        !end if
		mu =sqrt((theta*G(i))/(Ep*Ap))

		omega = Kb/(Ep*Ap*mu)
				
		K=Ep*Ap*mu*(omega+tanh(thickness*mu))/(1+omega*tanh(thickness*mu))
        
        Kb = K !update base stiffness of next (upper) segment as the total stiffness of the current segment
				
	end do
			  
	settle = DL/K

    end function


        !A single instance of the settlement calculation    
    ! This version averages the shear modulus along the pile shaft rather than assessing the pile in each layer
    function esett1D_ave_friction(E,Bp,DL,Ep,nlayer,des,aveloc,Ap,pi,G,v,thicknesses) result(settle) !get the settlement for a particular pile depth

        real, intent(in) :: pi, E(:) !E(nlayer)
        real, intent(in) ::  G(:) !G(nlayer)
        real, intent(in) :: Bp, Ep, Ap !pile diameter, stiffness, area
        
        integer, intent(in) ::  nlayer ! number of layers
        real,intent(in) :: DL !design load
        
        real :: settle ! settlement of final design
        real, intent(in) :: des
        
        real, intent(in) :: aveloc(:) !aveloc(nlayer+1) !location of layer boundaries in terms of depth
        
        real thicknesses(:)
        
        
        real, intent(in) :: v !poisson's ratio
        
          
        
         do i = 1,nlayer !find out which layer the pile gets in to
            if(des >= aveloc(i) .and. des < aveloc(i+1)) then
                penelayer = i
                exit
            end if
        end do
                  
        theta = 2*pi/log(5*des*(1-v)/Bp)

        ! Cacluate base stiffness
        do n = nlayer,penelayer+1,-1      !get weights from below the layer the pile is founded in
            weights(n) = -exp(-dr*(aveloc(n+1)-des))/dr + exp(-dr*(aveloc(n)-des))/dr !area = upper integral - lower integral
            effe(n) = weights(n)/E(n) !weighted young's modulus for current component (inverted for harmonic average calculation)
        end do
        weights(penelayer) = -exp(-dr*(aveloc(n+1)-des))/dr + 1/dr !get weight from layer pile is founded in 
        effe(penelayer) = weights(penelayer)/E(penelayer)
        effe(1) = sum(weights(penelayer:))/sum(effe(penelayer:)) !store weighted harmonic average in first element
        
        Kb = ((Bp*effe(1))/(1-v**2)) * (1+0.65*Bp/(aveloc(nlayer+1)-des)) !this (1 + 0.65...) bit accounts for rock at a certain depth
        
        !get the thickness of the pile in each layer prior to the one the pile is based in
        thicknesses(:penelayer-1) = aveloc(2:penelayer) - aveloc(:penelayer-1)
        thicknesses(penelayer) = min(des,aveloc(penelayer+1)) - aveloc(penelayer) ! do base layer

        
        ! take the weighted arithmetic average of the shear modulus values along the shaft
        Gave = sum(thicknesses(:penelayer) * G(:penelayer))/des
        !Gave = product(G(:penelayer) ** thicknesses(:penelayer))**(1/des)
        

        mu =sqrt((theta*Gave)/(Ep*Ap))
    
        omega = Kb/(Ep*Ap*mu)
                    
        K=Ep*Ap*mu*(omega+tanh(des*mu))/(1+omega*tanh(des*mu))

                  
        settle = DL/K
                  
    end function
    
    
    
            !A single instance of the settlement calculation    
    ! Works on a rigid pile
    function esett1D_rigid(E,Bp,DL,Ep,nlayer,des,aveloc,Ap,pi,G,v,thicknesses) result(settle) !get the settlement for a particular pile depth

        real, intent(in) :: pi, E(:) !E(nlayer)
        real, intent(in) ::  G(:) !G(nlayer)
        real, intent(in) :: Bp, Ep, Ap !pile diameter, stiffness, area
        
        integer, intent(in) ::  nlayer ! number of layers
        real,intent(in) :: DL !design load
        
        real :: settle ! settlement of final design
        real, intent(in) :: des
        
        real, intent(in) :: aveloc(:) !aveloc(nlayer+1) !location of layer boundaries in terms of depth
        
        real thicknesses(:)
        
        
        real, intent(in) :: v !poisson's ratio
        
          
        
         do i = 1,nlayer !find out which layer the pile gets in to
            if(des >= aveloc(i) .and. des < aveloc(i+1)) then
                penelayer = i
                exit
            end if
        end do
                  
        

        ! Cacluate base stiffness
        do n = nlayer,penelayer+1,-1      !get weights from below the layer the pile is founded in
            weights(n) = -exp(-dr*(aveloc(n+1)-des))/dr + exp(-dr*(aveloc(n)-des))/dr !area = upper integral - lower integral
            effe(n) = weights(n)/G(n) !weighted young's modulus for current component (inverted for harmonic average calculation)
        end do
        weights(penelayer) = -exp(-dr*(aveloc(n+1)-des))/dr + 1/dr !get weight from layer pile is founded in 
        effe(penelayer) = weights(penelayer)/G(penelayer)
        effe(1) = sum(weights(penelayer:))/sum(effe(penelayer:)) !store weighted harmonic average in first element
        
        
        !get the thickness of the pile in each layer prior to the one the pile is based in
        thicknesses(:penelayer-1) = aveloc(2:penelayer) - aveloc(:penelayer-1)
        thicknesses(penelayer) = min(des,aveloc(penelayer+1)) - aveloc(penelayer) ! do base layer

        ! do stiffness for pile segment in bottom-most layer
        theta = 2*pi/log(5*max(Bp,thicknesses(penelayer))*(1-v)/Bp)
        K = G(penelayer) * Bp * ((2.0/(1-v)) * effe(1) / G(penelayer) + theta * thicknesses(penelayer) / Bp) !* (1 + 0.65 * Bp / (aveloc(nlayer+1) - des))

        ! loop upwards through the layer segments, using the total stiffness of the segment immediately below as the base stiffness
        do i = penelayer-1,1,-1
            theta = 2*pi/log(5*max(Bp,thicknesses(i))*(1-v)/Bp)
            K =  K + G(i) * theta * thicknesses(i)
        end do

                  
        settle = DL/K
                  
    end function
    
    
    
        !Closed form 2 layer solution (doesn't appear to work properly)
    function esett1D_2L(E,Bp,DL,Ep,nlayer,des,aveloc,Ap,pi,G,v,thicknesses) result(settle) !get the settlement for a particular pile depth

        real, intent(in) :: pi, E(:) !E(nlayer)
        real, intent(in) ::  G(:) !G(nlayer)
        real, intent(in) :: Bp, Ep, Ap !pile diameter, stiffness, area
        
        integer, intent(in) ::  nlayer ! number of layers
        real,intent(in) :: DL !design load
        
        real :: settle ! settlement of final design
        real, intent(in) :: des
        
        real, intent(in) :: aveloc(:) !aveloc(nlayer+1) !location of layer boundaries in terms of depth
        
        real thicknesses(:)
        real mu1,mu2,muave
        
        
        real, intent(in) :: v !poisson's ratio
        
          
        
         do i = 1,nlayer !find out which layer the pile gets in to
            if(des >= aveloc(i) .and. des < aveloc(i+1)) then
                penelayer = i
                exit
            end if
        end do
                  
        theta = 2*pi/log(5*des*(1-v)/Bp)

        ! Cacluate base stiffness
        do n = nlayer,penelayer+1,-1      !get weights from below the layer the pile is founded in
            weights(n) = -exp(-dr*(aveloc(n+1)-des))/dr + exp(-dr*(aveloc(n)-des))/dr !area = upper integral - lower integral
            effe(n) = weights(n)/E(n) !weighted young's modulus for current component (inverted for harmonic average calculation)
        end do
        weights(penelayer) = -exp(-dr*(aveloc(n+1)-des))/dr + 1/dr !get weight from layer pile is founded in 
        effe(penelayer) = weights(penelayer)/E(penelayer)
        effe(1) = sum(weights(penelayer:))/sum(effe(penelayer:)) !store weighted harmonic average in first element
        
        Kb = ((Bp*effe(1))/(1-v**2)) * (1+0.65*Bp/(aveloc(nlayer+1)-des)) !this (1 + 0.65...) bit accounts for rock at a certain depth
        
        !get the thickness of the pile in each layer prior to the one the pile is based in
        thicknesses = 0
        thicknesses(:penelayer-1) = aveloc(2:penelayer) - aveloc(:penelayer-1)
        thicknesses(penelayer) = min(des,aveloc(penelayer+1)) - aveloc(penelayer) ! do base layer

        
        ! take the weighted arithmetic average of the shear modulus values along the shaft
        Gave = sum(thicknesses(:penelayer) * G(:penelayer))/des
        !Gave = product(G(:penelayer) ** thicknesses(:penelayer))**(1/des)
        

        mu1 =sqrt((theta*G(1))/(Ep*Ap))
        mu2 =sqrt((theta*G(2))/(Ep*Ap))
        muave = sqrt((theta*Gave)/(Ep*Ap))
    
        omega = Kb/(Ep*Ap*muave)
                    
        if (penelayer == 1) then
            !omega = Kb/(Ep*Ap*mu1)
            K=Ep*Ap*mu1*(omega+tanh(des*mu1))/(1+omega*tanh(des*mu1))
        
        else
            !omega = Kb/(Ep*Ap*mu2)
            !K=Ep*Ap*mu1* (mu1*tanh(thicknesses(1)*mu1) + mu1*omega*tanh(thicknesses(1)*mu1)*tanh(thicknesses(2)*mu2) + mu2*omega + mu2*tanh(thicknesses(2)*mu2) ) / &
            !         (mu1 + mu1+omega*tanh(thicknesses(2)*mu2) + mu2*omega*tanh(thicknesses(1)*mu1) + mu2*tanh(thicknesses(1)*mu1)*tanh(thicknesses(2)*mu2))
            K=Ep*Ap*mu1* (mu1*tanh(thicknesses(1)*mu1) + mu1*omega*tanh(thicknesses(1)*mu1)*tanh(thicknesses(2)*mu2) + mu2*omega + mu2*tanh(thicknesses(2)*mu2) ) / &
                     (mu1 + mu1+omega*tanh(thicknesses(2)*mu2) + mu2*omega*tanh(thicknesses(1)*mu1) + mu2*tanh(thicknesses(1)*mu1)*tanh(thicknesses(2)*mu2))
        end if
                  
        settle = DL/K
                  
    end function
    
    
    

    
      
    !Get the pile designs from each investigation and MC realisation.
    !Works by iteratively improving the solution estimate through the bisection method.
    !Design is undertaken according to settlement, determined using the Mylonakis and Gazetas method.
    
    !This subroutine also finds the true pile settlement in the original, actual soil, as well as differential settlement.
      
            subroutine multides_1D(prad,nrep_MC,ninv,npl,nlayer,layer_at_piles,evals,sides,rel_loads,load_con, destol, buildingweight,dz,nzew,plocation,sdata,CKheights,diffset,rel_loads_app,allones,load_con_app,CKprops)
        

            integer, intent(in) :: prad(2)			!pile width in x,y directions
            integer, intent(in) :: nlayer
            integer counter,pile,iter,inv,pd,p2                   !loop counters
            real(4) disp(nrep_MC,ninv)                        !pile settlement
            integer,intent(in) ::     nrep_MC,ninv,npl
            real, intent(in) :: rel_loads(:)		!srelative pile loads based on tributary area of building
            real, intent(in) :: rel_loads_app(:)    !reduced set of relative applied loads
            integer, intent(in) :: load_con_app(:)  !connectivity vector for applied loads
            real :: rel_loads_temp(npl) !temporary applied loads
            integer :: load_con_temp(npl)
	        integer, intent(in) :: load_con(:)			!connectivity vector saying which pile has which load
            real, intent(in) :: destol !design tolerance (mm)
            real :: temp_destol !scaled design tolerance 
            real, intent(in) :: buildingweight
            real :: sum_loads !sum of pile loads
            real, parameter :: pi = 3.14159
            real(8), intent(in) :: dz !length of a soil element
            integer, intent(in) :: nzew !number of soil elements in vertical direction
            real Bp,Ap,Ep,DL !pile diameter, area, stiffness and applied load
            real v,G(nlayer),Eave(nlayer) !soil poisson's ratio, shear modulus of each layer
            real G_ck(nlayer),sdata_tmp(nlayer) !shear modulus and current young's modulus for complete knowledge profile
            logical, intent(in) :: allones !true if all investigations have a single borehole (special case that can be highly optimised)
            
            
            real, intent(in) :: CKprops(:,:,:) !effective young's modulus at each pile in each layer according to the CK soil
            real, intent(in) :: layer_at_piles(:,:,:,:) ! layer_at_piles(nrep_MC,ninv,nlayer-1,npl) !layer depth information
            real, intent(in) :: evals(nrep_MC,ninv,nlayer) !evals(ninv,size(inv_reduction),nlayer)
            real, intent(in) :: plocation(npl,2)
            integer :: pileindex(npl)
            
            
            !complete knowledge stuff
            real, intent(in) :: sdata(nlayer,nrep_MC)
            real, intent(in) ::  CKheights(:,:,:) !CKheights(nrep_MC,nlayer-1,npl)
            real :: dist((npl)**2) !distances between piles
            real :: trueset(npl) !true, actual CK settlement for each pile
            real :: truesetfull(npl,nrep_MC,ninv)
            real, intent(out) :: diffset(nrep_MC,ninv) !maximum differential settlment
            real difftemp(npl**2) !differential settlements of piles
            logical :: maskarr2(npl)
            
            real, intent(out) :: sides(:,:,:) !sides(npl,nrep_mc,ninv)   !SI pile designs
            
            real nan                        !nan
            real leftdisp,rightdisp,cdisp,olddisp   !displacements associated with the left bracket, right bracket and centre
            real left,right,mid          !pile lengths (elements) associated with the left/right brackets and centre
            real aveloc(nlayer+1)           !layer locations of current case, including top and bottom of soil
            real desres                 !design resolution. E.g. 0.1 rounds up to the nearest 0.1 metres.

            real mscale !scale factor for the regular falsi method
            real mcheck !make sure the mscale is feasible
            !real testsets(300)
            
            real(8) start,finish ,diff1,diff2,diff3,diff4,diff0,rate
            integer(8) allstart,allfinish
            
            real left2,right2,oldleft,oldright,hardright
            real leftdisp2,rightdisp2
            integer counter2
            integer side
            real value
            integer npltemp !temporary number of piles
            logical sameloads !check to see if all piles have the same load
            
            logical, parameter :: bisection = .false.
            
            
            diff1=0
            diff2=0
            diff3=0
            diff4=0
            diff0=0
            
            !write(*,*) 'here'
            !read(*,*)
            
            
            nan = 0.0 !generate undefined value
  	        nan = 0.0/nan
            
            
            !some parameters 
            Bp = (prad(1)+prad(2))/2*dz
            Ap = pi*Bp**2/4 !get area of pile based on diameter
            !Ep = 10**20 !huge(Ep) !30000
            Ep = 1
            do i = 1,20 !setting it directly doesn't seem to work for some reason
                Ep = Ep * 10
            end do
            v=0.3
            desres = 0.1
            Eps = epsilon(v)*1000 !this tolerance is less than a milimeter, but should be enough for the finite/imperfect math
            sameloads = all(rel_loads == rel_loads(1)) !check whether all piles have the same load
            
            !allocate local variables in the settlement subroutine
            allocate(effe(nlayer),weights(nlayer))
      
            CALL system_clock(count_rate=rate)
            
            
            
            !get sum of pile loads
	        sum_loads = 0.0
	        do pile = 1,npl
                if(load_con(pile) == 0) cycle  !skip piles with no associated load
		        sum_loads = sum_loads + rel_loads(load_con(pile))
            end do
            
            !calculate distances between each pile here instead of in the heavily nested loop below
	        counter = 0
	        do pile = 1,npl-1
		        do p2 = pile+1,npl
			        counter = counter + 1
			        dist(counter) = sqrt((plocation(pile,1)-plocation(p2,1))**2 + (plocation(pile,2)-plocation(p2,2))**2)
		        end do
            end do
            
            
            !an optimisation for single borehole case
            if(allones) then
                npltemp = size(rel_loads_app) !only do as many piles as there are pile load cases
                rel_loads_temp(:npltemp) = rel_loads_app
                load_con_temp = load_con_app
                pileindex = 1 !in this case, pile layer depth information is only stored in the first index
            else
                npltemp = npl
                rel_loads_temp = rel_loads
                load_con_temp = load_con
                pileindex = [(pd,pd=1,npl)]
            end if
            


            mcheck = 100 * epsilon(mcheck)
            mscale = 1
            
            !Check if it's in the outer bounds
            oldleft = Bp             !pdepths(1)*dz        !This method doesn't seem to be valid for pile lengths shorter than its diameter, so set diameter as the lower bound.
            !Set the upper bound at bedrock depth minus a very, very small amount. Otherwise there will be a divide-by-zero in the esett subroutine if the pile is founded exactly on bedrock.
            hardright = nzew*sngl(dz)*(1-epsilon(nan)*10) !pdepths(npdepths)*dz       !use bedrock as upper bound instead of the longest pile
  
            oldright = hardright
            aveloc = [0.0,(sum(layer_at_piles(:,:,j,:))/size(layer_at_piles(:,:,j,:)),j=1,nlayer-1),sngl(nzew*dz)] 
        !     Eave = sum(sum(evals,1),1)/(nrep_MC*ninv)
        !     G = Eave/(2*(1+v)) !convert young's modulus to shear modulus
        !     DL = buildingweight*sum(rel_loads)/(sum_loads*size(rel_loads))
        !     
        !     do while (oldright > oldleft + desres) !exit after the brackets are consequtive elements.
        !               
        !                         
 			    !mid = (oldleft + oldright) / 2
        !         !call prep_fem2d(efld2d,0.3,femvech,femvecv,prad(1),cdisp,mid,dz)
        !         cdisp = esett1D(Eave,Bp,DL,Ep,nlayer,mid,aveloc,Ap,pi,G,v)
        !               
 			    !if ( cdisp >= destol ) then
 				   ! oldleft = mid
 			    !else
 				   ! oldright = mid
 			    !end if
 	      !               
        !     end do
        !    
        !    !oldleft = Bp
        !    !oldright = hardright
            
        
            do iter = 1,nrep_MC                                     !loop through realisations
                
                
                !write(*,*) iter

                
                invloop: do inv = 1,ninv                                     !loop through investigations
                    

                    

                    
                    !call cpu_time(start)
                    
                    do pile = 1,npltemp                                 !loop through piles 
                        
                        !scale the design tolerance (mm) according to the applied load
                        
                        !CALL system_clock(allstart)
                        
                        DL = buildingweight*rel_loads_temp(pile)/sum_loads
                        

                        
                        !aveloc = [0.0,(layer_at_piles(iter,inv,j,pile),j=1,nlayer-1),max(maxval(layer_at_piles(iter,inv,:,pile)),sngl(nzew*dz))]
                        aveloc(2:nlayer) = layer_at_piles(iter,inv,:,pileindex(pile)) !aveloc(2)
                        G = evals(iter,inv,:)/(2*(1+v)) !convert young's modulus to shear modulus
                        
        
                        
               !         CALL system_clock(allfinish)
			            !diff0 = diff0 + dble(allfinish-allstart)/rate
               !         CALL system_clock(allstart)
                        
                        
                        
                        !value = Muller(Bp,real(Bp+nzew*dz)/2,real(nzew*dz),evals(iter,inv,:),Bp,DL,Ep,nlayer,aveloc,Ap,pi,G,v,destol)

                    ! ---- perform bisection method to design pile ------
	
           
                        !call cpu_time(start)
                        
		
                        
                        leftdisp = esett1D(evals(iter,inv,:),Bp,DL,Ep,nlayer,oldleft,aveloc,Ap,pi,G,v)
                        
                        
   
                        
                        !call cpu_time(finish)
                        !diff0 = diff0 +  finish-start
                        !call cpu_time(start)
                        
                        
                        
                        rightdisp = esett1D(evals(iter,inv,:),Bp,DL,Ep,nlayer,oldright,aveloc,Ap,pi,G,v)
   
             
                        
               !         CALL system_clock(allfinish)
			            !diff2 = diff2 + dble(allfinish-allstart)/rate
               !         CALL system_clock(allstart)
                        
                        
                        !call cpu_time(finish)
                        !diff1 = diff1 +  finish-start
                        !call cpu_time(start)
                        
                        !Check to see if the design lies within the current bounds.
                        ! if (leftdisp < destol .or. rightdisp > destol) then
                        ! 
                        !     !If it doesn't, move the left and right brackets outward as needed, unless there is no more room
                        !     do while (leftdisp < destol .and. oldleft - desres >= bp)
                        !         oldleft = oldleft - desres
                        !         leftdisp = esett1d(evals(iter,inv,:),bp,dl,ep,nlayer,oldleft,aveloc,ap,pi,g,v)
                        !     end do
                        ! 
                        !     do while (rightdisp > destol .and. oldright + desres <= hardright)
                        !         oldright = oldright + desres
                        !         rightdisp = esett1d(evals(iter,inv,:),bp,dl,ep,nlayer,oldright,aveloc,ap,pi,g,v)
                        !     end do
                        !     
                        !     !If the design still isn't in the root, then mark it invalid
                        !     if (leftdisp < destol .or. rightdisp > destol) then
                        !         sides(:,iter,inv) = nan !if a single pile is invalid, no point doing the rest. Skip to next
                        !         diffset(iter,inv) = nan
 			                    !exit invloop
                        !     end if
                        ! end if
!                         
                              !If the design still isn't in the root, then mark it invalid
                             if (leftdisp < destol .or. rightdisp > destol) then
                                 sides(:,iter,inv) = -1 !if a single pile is invalid, no point doing the rest. Skip to next
                                 diffset(iter,inv) = -1
 			                    cycle invloop
                             end if
                        
                        !call cpu_time(finish)
                        !diff2 = diff2 +  finish-start
                        !call cpu_time(start)
                        
               !         CALL system_clock(allfinish)
			            !diff2 = diff2 + dble(allfinish-allstart)/rate
               !         CALL system_clock(allstart)
               !
               !         
               !         CALL system_clock(allfinish)
			            !diff2 = diff2 + dble(allfinish-allstart)/rate
               !         CALL system_clock(allstart)
                        
                        ! call cpu_time(finish)
                        !diff3 = diff3 +  finish-start
                        !call cpu_time(start)
                      !  
                      !      leftdisp2 = leftdisp
                      !      rightdisp2 = rightdisp
                      !  
                      !  !if(bisection) then
                      !      leftdisp2 = leftdisp
                      !      rightdisp2 = rightdisp
                      !      left2=left
                      !      right2=right
	                     !  !
		                    !!Otherwise, use bisection method to find interval containing true design
                      !  
                      !      counter2 = 0
                      !  
                      !
                      		right = oldright
                            left = oldleft
                      
		                    do while (right > left + desres) !exit after the brackets are consequtive elements.
                      
                               !counter2 = counter2 + 1
                               
			                       mid = (left + right) / 2
                                   !call prep_fem2d(efld2d,0.3,femvech,femvecv,prad(1),cdisp,mid,dz)
                                   cdisp = esett1D(evals(iter,inv,:),Bp,DL,Ep,nlayer,mid,aveloc,Ap,pi,G,v)
                      
			                       if ( cdisp >= destol ) then
				                       left = mid
                                       !leftdisp = cdisp
			                       else
				                       right = mid
                                       !rightdisp = cdisp
			                       end if
	                     
                            end do
               !             
               !             
               !          CALL system_clock(allfinish)
			            !diff3 = diff3 + dble(allfinish-allstart)/rate
               !         CALL system_clock(allstart)
                            
                      !      
                      !      !binary search would have converged by now assuming a valid design exists.
                      !      sides(pile,iter,inv) = right !make sure it uses the LARGER pile length (greater than the design tol)
                        
                       !else
                        
                        !Otherwise, use bisection method to find interval containing true design
                            !counter = 0
                            !side = 0
                            !rightdisp = rightdisp - destol
                            !leftdisp = leftdisp - destol
                            !olddisp = rightdisp
                            !!mid = (left * rightdisp - right * leftdisp)/(rightdisp - leftdisp)
                            !
                            !mid = (left+right)/2
                            !call esett1D(evals(iter,inv,:),Bp,DL,Ep,nlayer,mid,cdisp,aveloc,Ap,pi,G,v)
                            !cdisp = cdisp - destol
                            !
                            !do while (cdisp > 0)
                            !    
                            !    counter = counter + 1
                            !    
                            !    oldleft = mid
                            !    mid = max((left * cdisp - mid * leftdisp) / (cdisp - leftdisp),oldleft+2)
                            !    left = oldleft
                            !    
                            !    call esett1D(evals(iter,inv,:),Bp,DL,Ep,nlayer,mid,cdisp,aveloc,Ap,pi,G,v)
                            !    cdisp = cdisp - destol
                            !    right = mid !save the righthand point in case it's over the root
                            !    
                            !end do
                            
                       !      rightdisp = rightdisp - destol
                       !      leftdisp = leftdisp - destol
                       !
                       !      !counter = 0
                       !      side = 0
                       !      right = oldright
                       !      left = oldleft
                       !      
 		                    !do while (right > left + desres) !exit after the brackets are consequtive elements.
                       !          
                       !          !counter = counter + 1
                       !
 			                   !     !mid = (left + right) / 2
                       !              
                       !              mid = (left * rightdisp - right * leftdisp)/(rightdisp - leftdisp)
                       !              
                       !              !if(isnan(mid)) then
                       !              !    write(*,*) 'uh oh'
                       !              !end if
                       !              
                       !              !call prep_fem2d(efld2d,0.3,femvech,femvecv,prad(1),cdisp,mid,dz)
                       !              cdisp = esett1D(evals(iter,inv,:),Bp,DL,Ep,nlayer,mid,aveloc,Ap,pi,G,v) - destol
                       !              
 			                   !     if ( cdisp > 0 ) then
                       !                  if(side == 1) then
                       !                      mscale = 1 - cdisp/leftdisp
                       !                      if (.not. mscale >= mcheck)  mscale = 0.5 !note the .not. is neccessary to avoid NaNs when the exact value is found
                       !                      rightdisp = rightdisp * mscale
                       !                  end if
                       !                  left = mid
                       !                  leftdisp = cdisp
                       !                  side = 1
                       !              else if (cdisp < 0) then
                       !                  if(side == -1) then
                       !                      mscale = 1 - cdisp/rightdisp
                       !                      if (.not. mscale >= mcheck) mscale = 0.5
                       !                      leftdisp = leftdisp * mscale
                       !                  end if
 				                  !      right = mid
                       !                  rightdisp = cdisp
                       !                  side = -1
                       !              else !exact solution found
                       !                  right = mid
                       !                  exit
                       !              end if        
                       !              
                       !
                       !              
                       !              !if(rightdisp + destol <= mcheck .or. leftdisp + destol <= mcheck .or. mscale <= mcheck) then
                       !              !    write(*,*) 'found'
                       !              !    read(*,*)
                       !              !end if
                       !                  
 	                     !
                       !      end do
                            
                            !write(*,*) counter
                            
                            !binary search would have converged by now assuming a valid design exists.
                            sides(pile,iter,inv) = right !make sure it uses the LARGER pile length (greater than the design tol)
                        
                        !end if
                        

                    !call cpu_time(finish)
                    !diff4 = diff4 +  finish-start
                    !call cpu_time(start)
                    
           !         CALL system_clock(allfinish)
			        !diff4 = diff4 + dble(allfinish-allstart)/rate
           !         CALL system_clock(allstart)
                        
                        

                    end do
                    
                    
                    !if(inv_bh(inv) == 1 .and. sameloads) sides(:,iter,inv) = sides(1,iter,inv)
                    
              
                    
                    !call cpu_time(finish)
                    !diff1 = diff1 +  finish-start
                    !call cpu_time(start)
                    !

                ! ---- Get differential settlements for current set of piles ----
                    
                    !if(any(sides(:npl,iter,inv) <= 0.5 ) .or. any(sides(:npl,iter,inv) >= 39.9999 )) then
                    !    write(*,*) 'damn'
                    !    read(*,*)
                    !end if
                    
                    
                    !get true pile settlement for individual piles
                    do pile = 1,npl
                        !sdata_tmp = sdata(:,iter)
                        sdata_tmp = CKprops(iter,:,pile)
                        G_ck = sdata_tmp/(2*(1+v)) !convert young's modulus to shear modulus
                        !aveloc = [0.0,(CKheights(iter,j,pile),j=1,nlayer-1),max(maxval(CKheights(iter,:,pile)),sngl(nzew*dz))]
                        aveloc(2:nlayer) = CKheights(iter,:,pile) 
                        DL = buildingweight*rel_loads(pile)/sum_loads
                        trueset(pile) = esett1D(sdata_tmp,Bp,DL,Ep,nlayer,sides(load_con_temp(pile),iter,inv),aveloc,Ap,pi,G_ck,v)
                        truesetfull(pile,iter,inv) = trueset(pile)
                    end do
                    
                    
                    

                    

                    !write(*,*) sides(:,iter,inv)
                    !write(*,*) CKheights(iter,:,:)
                    !write(*,*) trueset(:)
                    !write(*,*) sdata_tmp
                    !write(*,*)
                    
                    !call cpu_time(finish)
                    !diff2 = diff2 +  finish-start
                    !call cpu_time(start)
                    
 
		            !calculate differential settlement between all combinations of piles
			        counter = 0
			        do pile = 1,npl-1
				        do p2 = pile+1,npl
					        counter = counter + 1
					        difftemp(counter) = abs(trueset(pile) - trueset(p2))/dist(counter)
				        end do
                    end do
                    
                    
                    !call cpu_time(finish)
                    !diff3 = diff3 +  finish-start
                    !call cpu_time(start)
                    
			        !maskarr(:count) = difftemp(:count) /= nan		!Take max differential settlement, ignoring nan values
			        !diffset(iter,inv) = maxval(difftemp(:count),mask=maskarr(:count))
			        diffset(iter,inv) = maxval(difftemp(:counter)) !save maximum differential settlement of pile group              
                    
                    !
                    !call cpu_time(finish)
                    !diff4 = diff4 +  finish-start
                    !
                    
                end do invloop
            end do
            
            !write(*,*) diff0,diff1,diff2,diff3,diff4
            !read(*,*)
            
            !write(*,*) 'here2'
            !read(*,*)
            
            if(.false.) then
            open(2034,file='tempdisp.txt')
            open(2044,file='tempE1.txt')
            open(2045,file='tempE2.txt')
            DL = 1
            do iter = 1,nrep_MC
                do inv = 1,ninv
                    
                    aveloc = [0.0,(layer_at_piles(iter,inv,j,pile),j=1,nlayer-1),max(maxval(layer_at_piles(iter,inv,:,pile)),sngl(nzew*dz))]
                    G = evals(iter,inv,:)/(2*(1+v)) !convert young's modulus to shear modulus
                    
                    if(any(isnan(evals(iter,inv,:)))) then
                        disp(iter,inv) = nan
                    else
                        disp(iter,inv) = esett1D(evals(iter,inv,:),Bp,DL,Ep,nlayer,1.0,aveloc,Ap,pi,G,v)
                    end if
                end do
                write(2034,'(100000(F15.5,X))') disp(iter,:)
                write(2044,'(100000(F15.5,X))') evals(iter,:,1)
                write(2045,'(100000(F15.5,X))') evals(iter,:,2)
            end do
            close(2034)
            close(2044)
            close(2045)
            end if
            
                !open(505,file='pset_spx.dat',access='stream')
                !write(505) truesetfull
                !close(505)


            ! open(unit=2821,file='trueset.txt')
            ! do i = 1,nrep_mc
            !     write(2821,'(10000(F5.1,X))')  truesetfull(:,i,1)
            ! end do
            ! close(2821)
            ! stop
            
            deallocate(effe,weights)
            

            end subroutine
            
            
            
            
            !This is a streamlined version of the above subroutine that gets pile designs for a single deterministic instance
            subroutine multides_1D_det(prad,npl,nlayer,layer_at_piles,evals,sides,rel_loads,load_con, destol, buildingweight,dz,nzew)
        

            integer, intent(in) :: prad(2)			!pile width in x,y directions
            integer, intent(in) :: nlayer
            integer counter,pile,iter,inv,pd,p2                   !loop counters
    
            integer,intent(in) ::     npl
            real, intent(in) :: rel_loads(:)		!srelative pile loads based on tributary area of building
            real :: rel_loads_temp(npl) !temporary applied loads
            integer :: load_con_temp(npl)
	        integer, intent(in) :: load_con(:)			!connectivity vector saying which pile has which load
            real, intent(in) :: destol !design tolerance (mm)
            real :: temp_destol !scaled design tolerance 
            real, intent(in) :: buildingweight
            real :: sum_loads !sum of pile loads
            real, parameter :: pi = 3.14159
            real(8), intent(in) :: dz !length of a soil element
            integer, intent(in) :: nzew !number of soil elements in vertical direction
            real Bp,Ap,Ep,DL !pile diameter, area, stiffness and applied load
            real v,G(nlayer),Eave(nlayer) !soil poisson's ratio, shear modulus of each layer
            
            real, intent(in) :: layer_at_piles(:,:) ! layer_at_piles(nrep_MC,ninv,nlayer-1,npl) !layer depth information
            real, intent(in) :: evals(nlayer) !evals(ninv,size(inv_reduction),nlayer)
            
            
            
            real, intent(out) :: sides(:) !sides(npl,nrep_mc,ninv)   !SI pile designs
            
            real nan                        !nan
            real leftdisp,rightdisp,cdisp,olddisp   !displacements associated with the left bracket, right bracket and centre
            real left,right,mid          !pile lengths (elements) associated with the left/right brackets and centre
            real aveloc(nlayer+1)           !layer locations of current case, including top and bottom of soil
            real desres                 !design resolution. E.g. 0.1 rounds up to the nearest 0.1 metres.


            !real testsets(300)
            
            integer(8) allstart,allfinish
            
            real left2,right2,oldleft,oldright,hardright
            integer counter2
            integer side
            real value
            integer npltemp !temporary number of piles
            logical sameloads !check to see if all piles have the same load
            
            logical, parameter :: bisection = .false.
            

            
            nan = 0.0 !generate undefined value
  	        nan = 0.0/nan
            
            
            !some parameters 
            Bp = (prad(1)+prad(2))/2*dz
            Ap = pi*Bp**2/4 !get area of pile based on diameter
            !Ep = 10**20 !huge(Ep) !30000
            Ep = 1
            do i = 1,20 !setting it directly doesn't seem to work for some reason
                Ep = Ep * 10
            end do
            v=0.3
            desres = 0.1
            Eps = epsilon(v)*1000 !this tolerance is less than a milimeter, but should be enough for the finite/imperfect math
            sameloads = all(rel_loads == rel_loads(1)) !check whether all piles have the same load
            
            !allocate local variables in the settlement subroutine
            allocate(effe(nlayer),weights(nlayer))
      
            
            
            
            !get sum of pile loads
	        sum_loads = 0.0
	        do pile = 1,npl
                if(load_con(pile) == 0) cycle  !skip piles with no associated load
		        sum_loads = sum_loads + rel_loads(load_con(pile))
            end do

            
            !Check if it's in the outer bounds
            oldleft = Bp             !pdepths(1)*dz        !This method doesn't seem to be valid for pile lengths shorter than its diameter, so set diameter as the lower bound.
            !Set the upper bound at bedrock depth minus a very, very small amount. Otherwise there will be a divide-by-zero in the esett subroutine if the pile is founded exactly on bedrock.
            hardright = nzew*sngl(dz)*(1-epsilon(nan)*10) !pdepths(npdepths)*dz       !use bedrock as upper bound instead of the longest pile
  
            oldright = hardright
            aveloc = [0.0,(sum(layer_at_piles(j,:))/size(layer_at_piles(j,:)),j=1,nlayer-1),sngl(nzew*dz)] 

        
                
                    do pile = 1,npl                                !loop through piles 
                        
                        !scale the design tolerance (mm) according to the applied load
             
                        DL = buildingweight*rel_loads(pile)/sum_loads
                        

                        
                        !aveloc = [0.0,(layer_at_piles(iter,inv,j,pile),j=1,nlayer-1),max(maxval(layer_at_piles(iter,inv,:,pile)),sngl(nzew*dz))]
                        aveloc(2:nlayer) = layer_at_piles(:,pile) !aveloc(2)
                        G = evals/(2*(1+v)) !convert young's modulus to shear modulus


                    ! ---- perform bisection method to design pile ------

                        leftdisp = esett1D(evals,Bp,DL,Ep,nlayer,oldleft,aveloc,Ap,pi,G,v)

                        rightdisp = esett1D(evals,Bp,DL,Ep,nlayer,oldright,aveloc,Ap,pi,G,v)
                    
                              !If the design still isn't in the root, then mark it invalid
                             if (leftdisp < destol) then 
                                 sides(pile) = -100
                                 cycle 
                             else if(rightdisp > destol) then
                                 sides(pile) = -200
 			                    cycle 
                             end if
                        

                      		right = oldright
                            left = oldleft
                      
		                    do while (right > left + desres) !exit after the brackets are consequtive elements.
                      
                               !counter2 = counter2 + 1
                               
			                       mid = (left + right) / 2
                                   !call prep_fem2d(efld2d,0.3,femvech,femvecv,prad(1),cdisp,mid,dz)
                                   cdisp = esett1D(evals,Bp,DL,Ep,nlayer,mid,aveloc,Ap,pi,G,v)
                      
			                       if ( cdisp >= destol ) then
				                       left = mid
                                       !leftdisp = cdisp
			                       else
				                       right = mid
                                       !rightdisp = cdisp
			                       end if
	                     
                            end do
          
                            !binary search would have converged by now assuming a valid design exists.
                            sides(pile) = right !make sure it uses the LARGER pile length (greater than the design tol)


                    end do
  
                    

    

            
            deallocate(effe,weights)
            

            end subroutine
            
            
            
            
            
            
                !Get settlement curves using the esett function. This is more for debugging than anything else
            subroutine multisetcurve(disp,disp3,pdepths,prad,npl,nlayer,layer_at_piles,evals,rel_loads,load_con, destol, buildingweight,dz,nzew)
            
            implicit none
        
            
            integer, parameter :: curveres = 100!number of points on ththe settlement curve
            real :: pdepths(:),disp(:),disp3(:) !pile lengths and associated displacements
            real :: disp2(size(disp))

            integer, intent(in) :: prad(2)			!pile width in x,y directions
            integer, intent(in) :: nlayer
            integer counter,pile,iter,inv,pd,p2                   !loop counters
    
            integer,intent(in) ::     npl
            real, intent(in) :: rel_loads(:)		!srelative dpile loads based on tributary area of building
            real :: rel_loads_temp(npl) !temporary applied loads
            integer :: load_con_temp(npl)
	        integer, intent(in) :: load_con(:)			!connectivity vector saying which pile has which load
            real, intent(in) :: destol !design tolerance (mm)
            real :: temp_destol !scaled design tolerance 
            real, intent(in) :: buildingweight
            real :: sum_loads !sum of pile loads
            real, parameter :: pi = 3.14159
            real(8), intent(in) :: dz !length of a soil element
            integer, intent(in) :: nzew !number of soil elements in vertical direction
            real Bp,Ap,Ep,DL !pile diameter, area, stiffness and applied load
            real v,G(nlayer),Eave(nlayer) !soil poisson's ratio, shear modulus of each layer
            
            real, intent(in) :: layer_at_piles(:,:) ! layer_at_piles(nrep_MC,ninv,nlayer-1,npl) !layer depth information
            real, intent(in) :: evals(nlayer) !evals(ninv,size(inv_reduction),nlayer)
            

            
            real nan                        !nan
            real leftdisp,rightdisp,cdisp,olddisp   !displacements associated with the left bracket, right bracket and centre
            real left,right,mid          !pile lengths (elements) associated with the left/right brackets and centre
            real aveloc(nlayer+1)           !layer locations of current case, including top and bottom of soil
            real desres                 !design resolution. E.g. 0.1 rounds up to the nearest 0.1 metres.


            !real testsets(300)
            
            integer(8) allstart,allfinish
            real start, finish
            
            
            integer counter2
            integer side
            real value
            integer npltemp !temporary number of piles
            logical sameloads !check to see if all piles have the same load
            
            logical, parameter :: bisection = .false.
            

            
            nan = 0.0 !generate undefined value
  	        nan = 0.0/nan
            
            
            !some parameters 
            Bp = (prad(1)+prad(2))/2*dz
            Ap = pi*Bp**2/4 !get area of pile based on diameter
            !Ep = 10**20 !huge(Ep) !30000
            Ep = 1
            do i = 1,20 !setting it directly doesn't seem to work for some reason
                Ep = Ep * 10
            end do
            Ep = 30000
            v=0.3
            desres = 0.1
            Eps = epsilon(v)*1000 !this tolerance is less than a milimeter, but should be enough for the finite/imperfect math
            sameloads = all(rel_loads == rel_loads(1)) !check whether all piles have the same load
            
            !allocate local variables in the settlement subroutine
            allocate(effe(nlayer),weights(nlayer),thicknesses(nlayer))
      
            
            
            
            !get sum of pile loads
	        sum_loads = 0.0
	        do pile = 1,npl
                if(load_con(pile) == 0) cycle  !skip piles with no associated load
		        sum_loads = sum_loads + rel_loads(load_con(pile))
            end do

            
            aveloc = [0.0,(sum(layer_at_piles(j,:))/size(layer_at_piles(j,:)),j=1,nlayer-1),sngl(nzew*dz)] 
            
        
                
                    do pile = 1,1 !npl                                !loop through piles  (first pile for now)
                        
                        !scale the design tolerance (mm) according to the applied load
             
                        DL = buildingweight*rel_loads(pile)/sum_loads
                        

                        
                        !aveloc = [0.0,(layer_at_piles(iter,inv,j,pile),j=1,nlayer-1),max(maxval(layer_at_piles(iter,inv,:,pile)),sngl(nzew*dz))]
                        aveloc(2:nlayer) = layer_at_piles(:,pile) !aveloc(2)
                        G = evals/(2*(1+v)) !convert young's modulus to shear modulus


		                    do pd = 1,size(pdepths)

                                   disp(pd) = esett1D(evals,Bp,DL,Ep,nlayer,pdepths(pd),aveloc,Ap,pi,G,v)
                                   !disp(pd) = esett1D_ave_friction(evals,Bp,DL,Ep,nlayer,pdepths(pd),aveloc,Ap,pi,G,v,thicknesses)
                                   disp3(pd) = esett1D_rigid(evals,Bp,DL,Ep,nlayer,pdepths(pd),aveloc,Ap,pi,G,v,thicknesses)
                                   !disp(pd) = esett1D_2L(evals,Bp,DL,Ep,nlayer,pdepths(pd),aveloc,Ap,pi,G,v,thicknesses)

                            end do
                    end do
   
  
            deallocate(effe,weights)
            

            end subroutine
     
	 
	 end module
		
