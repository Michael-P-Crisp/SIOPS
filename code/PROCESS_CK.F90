
    module process_CK

    use fem_prep
    use writesoils
    use variables
    use soilgen
    use sim2sd
    use weights


    implicit none

    contains

    !Get true pile performance (CK - complete knowledge), using the full information of the random field. Two options are provided:
    !   1. The Pseudo-incremental energy (PIE) method, an EXTREMELY fast approximation of linear-elastic finite element analysis, with only minor errors.
    !   2. A linear-elastic Finite Element Analysis (FEA) subroutine. Uses a variable-element mesh to significantly speed up the process, where LAS
    !       elements are geometrically averaged.

    ! PIE is 2 orders of magnitude faster than FEM, although it is 'relatively' inaccurate. PIE is really the only practical option for analysis if a
    ! super-computer is not available. Note that it only works for single-layer soil profiles with a constant mean.
    ! See "Effective Young's modulus of a spatially variable sol mass under a footing", Ching et al. (2018) for details.
    ! See description in 'fem_3d.F90' file for FEM details.

    ! A check for notable errors produced by PIE inaccuracy, as well as a very, very crude compensation mechanism is provided. See below for details.


    subroutine get_ck_set(soilseeds & !soil generation variables
    ,preps,npdepths,detdisp,plocation,pdepths,cg_tol,cg_limit,prad,femrad,femdepth,pieval,soilweight,nels,ck_pset,nrep_mc,femvech,femvecv,usepie,regmesh,startstress,datafolder)


    ! ----- soil generation variables (don't touch) ----

    integer NGS

    !virtual soil
    real(4) :: efldave(nxew,nyew,nzew)
    !averaged virtual soil across MC realisations
    integer, intent(in) :: soilseeds(:)

    real(4) :: sdata_ck(4) !CK soil parameters, same as SI but with unit mean stiffness
    real(4) :: pvr !lognormal variable
    integer iseed
    real randu

    integer xpos,ypos,zpos !offsets of random subset for the superset mode

    integer kseed !random seed

    ! ---- PIE and FEM variables -----

    integer, intent(in) :: npdepths							!number of pile depths to analyse
    integer, intent(in) :: pdepths(:)			!pile depths to analyse (elements)
    real, intent(out) :: detdisp(:)				!deterministic settlement for each pile depth (1 kN applied load. Soil stiffness 1 MPa)
    real(8), intent(in) :: cg_tol 					!max no. iterations for FEM subroutine
    integer, intent(in) :: cg_limit 				!tolerance for FEM
    integer, intent(in) :: preps(2)		!No. piles in x,y directions - repetitions
    integer, intent(in) :: prad(2)			!pile width in x,y directions
    integer, intent(in) :: plocation(:,:)		!Pile x,y coordinates
    logical, intent(in) :: regmesh !whether to use a regular mesh at resolution of the LAS field, or the specified variable element sizes

    real, intent(in) :: femvech(:), femvecv(:)
    integer, intent(in) :: startstress(:,:)
    integer centerpie(2) !the centre of the truncated PIE weight array

    integer, intent(in) :: femrad,femdepth !radius and depth of soil around pile (elements)

    real, intent(in) :: soilweight(:,:,:,:) !Soil weights required for the PIE method of getting CK settlement
    integer, intent(in) :: nels									!total number of soil elements in PIE method


    real, allocatable, intent(out) :: ck_pset(:,:,:) !true ck pile settlement for each realisation, pile and depth



    ! ---- finite element analysis variables -----
    integer, intent(in) :: nrep_mc

    logical,intent(in) :: usepie !whether to use FEM approximation such as the PIE method for single layers and Mylonakis and Gazetas method for multi-layer soils

    real :: disp !pile settlement from FEM
    real(8) :: pieval,pieval2			!effective soil elastic modulus for PIE method
    real :: sweightsum(npdepths) !,sweightsum2(npdepths) !the sum of soil weights for each pile length
    !real :: soilweight2(21,81,81,80)
    character(1000) str2
    real temp
    logical dekcheck(npdepths-1)

    !other variables
    integer pd,iter,x,y,z,jr,counter,i,j,k !loop counters
    real start,finish !time variables
    character(1000),intent(in) :: datafolder !a string reprebsenting the directory the data is stored in
    logical exists
    real(8), allocatable :: xyi(:,:)

    !multi layer variables
    integer :: indices(2,nxew,nyew) !x,y indices of layer boundary
    real height(nlayer-1,preps(1)*preps(2))
    real efld2d(nxew,nzew) !Virtual soil 2D representation; Uniform properties within each layer, horizontal layer boundary
    integer femrads(2)

    real piedisp1(21),piedisp2(21)


    !read in pile settlement curve values
    write(str2,'(A,A,I0,A,F4.2,A)') trim(datafolder),'settlement_prad-',prad(1),'_esize-',dz,'.txt'
    open(667, file= str2,status='old')
    read(667,*)
    do i = 1,npdepths
        read(667,*) temp,detdisp(i)
    end do
    close(667)



    if(singletrue) then
        centerpie(1) = maxval(startstress(1,:)) + 1
        centerpie(2) = maxval(startstress(2,:)) + 1

        !get sum of soil weights for weighted average calculation
        do pd = 1,npdepths
            sweightsum(pd) = sum(soilweight(pd,centerpie(1)-startstress(1,pd):centerpie(1)+startstress(1,pd)+prad(1)-1,centerpie(2)-startstress(2,pd):centerpie(2)+startstress(2,pd)+prad(2)-1,:startstress(3,pd)))
        end do
    end if


    allocate(ck_pset(npdepths,nrep_MC,preps(1)*preps(2)))

    !convert and use unit mean (1 MPa) soil stiffness for CK and scale at a later time.
    ! NOTE: Now the soil is generated as unit mean, and the load is scaled instead
    !sdata_ck(1) = 1.0
    !sdata_ck(2) = sdata(2)/sdata(1)
    !pvr = log(1.0 + sdata_ck(2) ** 2 / sdata_ck(1) ** 2)
    !sdata_ck(3) = log(sdata_ck(1)) - 0.5 * pvr
    !sdata_ck(4) = sqrt(pvr)


    !get the average virtual soil (see avesoil subroutine for details and explanation)
    !if (usepie <= 2) then
    !    call avesoil(soilseeds, nrep_MC,datafolder,efldave)
    !end if

    !get indices for layer interpolation
    ! do i=1,nxew
    !     indices(1,i,:) = i
    ! end do
    !do j = 1,nyew
    !	indices(2,:,j) = j
    ! end do
    !get 2D boundary element indices
    do i=1,nxew
        do j = 1,nyew
            indices(1,i,j) = i
            indices(2,i,j) = j
        end do
    end do
    femrads = femrad







    !Get a very big soil from which to extract a subset (for superset mode)
    !Note that superset won't work properly with anisotropy or with FEM specified for CK analysis.
    if (singletrue .and. superset) then
        write(*,*) 'Generating soil superset'
        kseed = randu(soilseeds(1)) * 1234567890 !ensure random numbers are consistent across MC realisations
        call RF3D(soil_dummy_var)!Get soil
        write(*,*)
        efld = log(efld) !log it here so it doesn't have to be done later; part of the PIE weighted geometric average calculation
    end if
    if(.not. singletrue) then
        !get index values
        counter=0
        allocate(xyi(2,nxew*nyew))

        do i=1,nyew
            do j=1,nxew
                counter = counter+1
                xyi(1,counter) = j
                xyi(2,counter) = i
            end do
        end do
    end if
    !open(19851,file='bigtest.txt')
    !do x=1,nxe
    !    do y = 1,nye
    !        do z = 1,nze
    !            write(19851,*) efld(x,y,z)
    !        end do
    !    end do
    !end do
    !close(19851)


    !output the average random field from the superset across Monte Carlo realisations to check that there
    !is no or negligible bias in soil properties.
    if(.false.) then
        efldave = 0
        do iter = 1,nrep_MC
            if((iter-int(10*iter/nrep_MC)*nrep_MC/10 == 0)) write(*,'(I0,A,X)',advance='no') iter*100/nrep_MC,'%'
            kseed = randu(soilseeds(iter)) * 1234567890 !ensure random numbers are consistent across MC realisations

            !generate random offsets
            xpos=NINT(randu(0)*(nxe-nxew))
            ypos=NINT(randu(0)*(nye-nyew))
            zpos=NINT(randu(0)*(nze-nzew))


            efldave = efldave + efld(xpos:xpos+nxew,ypos:ypos+nyew,zpos+1:zpos+nzew)

        end do

        open(4151,file='outputave.txt')
        do x=1,nxew
            do y = 1,nyew
                do z = 1,nzew
                    write(4151,*) efldave(x,y,z)
                end do
            end do
        end do

        close(4151)
        stop
    end if


    !get sum of soil weights for weighted average calculation


    write(*,*) 'Determining true pile settlements'
    call cpu_time(start)
    do iter = 1,nrep_MC


        !write(*,*) iter

        !Progress indicator
        if((iter-int(10*iter/nrep_MC)*nrep_MC/10 == 0)) write(*,'(I0,A,X)',advance='no') iter*100/nrep_MC,'%'


        if (singletrue .and. superset) then
            kseed = randu(soilseeds(iter)) * 1234567890 !ensure random numbers are consistent across MC realisations

            !generate random offsets
            xpos=NINT(randu(0)*(nxe-nxew))
            ypos=NINT(randu(0)*(nye-nyew))
            zpos=NINT(randu(0)*(nze-nzew))

            !open(4151,file='output.txt') startstress(1,:)

            do x=1,preps(1)*preps(2)
                !do concurrent (x=1:preps(1)*preps(2))
                do pd = 1,npdepths
                    !get effective young's modulus
                    !pieval = exp(sum(efld(xpos+plocation(x,1)-femrad:xpos+plocation(x,1)+femrad+prad(1)-1,ypos+plocation(x,2)-femrad:ypos+plocation(x,2)+femrad+prad(2)-1,zpos+1:zpos+femdepth) * soilweight(pd,:,:,:))/sweightsum(pd))
                    pieval = exp(sum(efld(xpos+plocation(x,1)-startstress(1,pd):xpos+plocation(x,1)+startstress(1,pd)+prad(1)-1,ypos+plocation(x,2)-startstress(2,pd):ypos+plocation(x,2)+startstress(2,pd)+prad(2)-1,zpos+1:zpos+startstress(3,pd)) * &
                        soilweight(pd,centerpie(1)-startstress(1,pd):centerpie(1)+startstress(1,pd)+prad(1)-1,centerpie(2)-startstress(2,pd):centerpie(2)+startstress(2,pd)+prad(2)-1,:startstress(3,pd)))/sweightsum(pd))
                    !Get pile settlement by scaling the deterministic value
                    ck_pset(pd,iter,x) = detdisp(pd)/pieval


                    !write(4151,*) ck_pset(:,iter,x)

                end do

                !close(4151)
                !stop
            end do



            !ck_pset(:,1,1)
        else





            !---- get pile settlement with depth, for each pile in the group ---
            if (singletrue) then !do single layer processing if in single layer mode
                if(usepie) then !use pie method to get settlement
                    ! --- generate soil ---
                    kseed = randu(soilseeds(iter)) * 1234567890 !ensure random numbers are consistent across MC realisations
                    call RF3D(soil_dummy_var)
                    !efld = efld - efldave !subtract average field for reasons given above


                    efld = log(efld) !required to take log to calculate weighted geometric average. Faster to do the whole thing at once.

                    do x=1,preps(1)*preps(2)
                        !do concurrent (x=1:preps(1)*preps(2))
                        do pd = 1,npdepths
                            !get effective young's modulus
                            !pieval = exp(sum(efld(plocation(x,1)-femrad:plocation(x,1)+femrad+prad(1)-1,plocation(x,2)-femrad:plocation(x,2)+femrad+prad(2)-1,:femdepth) * soilweight(pd,:,:,:))/sweightsum(pd))
                            pieval = exp(sum(efld(plocation(x,1)-startstress(1,pd):plocation(x,1)+startstress(1,pd)+prad(1)-1,plocation(x,2)-startstress(2,pd):plocation(x,2)+startstress(2,pd)+prad(2)-1,:startstress(3,pd)) * &
                                soilweight(pd,centerpie(1)-startstress(1,pd):centerpie(1)+startstress(1,pd)+prad(1)-1,centerpie(2)-startstress(2,pd):centerpie(2)+startstress(2,pd)+prad(2)-1,:startstress(3,pd)))/sweightsum(pd))
                            !Get pile settlement by scaling the deterministic value
                            ck_pset(pd,iter,x) = detdisp(pd)/pieval

                        end do
                    end do
                else         !use 3D finite element analysis to get settlement

                    ! --- generate soil ---
                    kseed = randu(soilseeds(iter)) * 1234567890 !ensure random numbers are consistent across MC realisations
                    call RF3D(soil_dummy_var)
                    efld = efld - efldave !subtract average field for reasons given above

                    do x=1,preps(1)*preps(2)
                        do pd = 1,npdepths
                            call prep_fem3D(efld(plocation(x,1)-femrad:plocation(x,1)+femrad+prad(1)-1,plocation(x,2)-femrad:plocation(x,2)+femrad+prad(2)-1,:femdepth),dz,femrad*2+prad(1),femdepth,prad,pdepths(pd),cg_tol,cg_limit,femvech,femvecv,disp,regmesh) !assumes elements are cubes
                            ck_pset(pd,iter,x) = disp
                        end do
                    end do
                end if
            else !multiple layers
                if (usepie) then
                    return !no pre-processing needed in multi-layer mode if FEM isn't used
                else
                    kseed = randu(soilseeds(iter)) * 1234567890

                    !Get layer boundary depths
                    call soil_layers(xyi)


                    !Get the equivallent uniform layer depths at each pile based on inverse distance weighting
                    call getheights(bfld,indices,plocation,nxew,nyew,preps(1)*preps(2),prad,2,femrads,nlayer,height)

                    !do x=1,size(height,1)
                    !    write(*,'(10000F8.3)') (height(x,pd),pd=1,size(height,2))
                    !end do
                    !write(*,*) bfld(2,:5,:5)
                    !read(*,*)

                    do x=1,preps(1)*preps(2) !Use 2D linear elastic analysis (assumes properties are constant in the horizontal direction)
                        !Build soil layer for current pile
                        efld2d = lmeans(nlayer,iter)                     !Fill with oldest, deepest soil
                        do i=nlayer-1,1,-1
                            efld2d(:,:nint(height(i,x))) = lmeans(i,iter)      !Progressively work forwards through time, adding newer layers
                        end do
                        !Get settlement for each pile depth
                        do pd = 1,npdepths
                            call prep_fem2d(efld2d,0.3,femvech,femvecv,prad(1),disp,pdepths(pd),dz)
                            ck_pset(pd,iter,x) = disp
                        end do
                    end do

                end if
            end if
        end if




    end do
    write(*,*) !put a newline symbol

    call cpu_time(finish)
    write(*,*) iter, finish-start

    !check that pile settlement is decreasing as pile length increases. This MUST be the case given the linear-elastic model assumptions, and
    !because the akima interpolation REQUIRES this to be the case. However, if there's a section of the curve that's particularly flat (horizontal)
    !then it's possible that the above assumption is not strictly true as a result of minor random errors.
    !Therefore, if pile settlement happens to increase with depth, it is assumed that the curve is fairly flat, and the offending point
    !is replaced with linear interpolation of the adjacent depths. If the offending point is one of the end point, then simply add or subtract
    !the smallest possible number to it.
    !Actually, forget the above compensation; it's simplest to just swap the values in the event that there are multiple offending bad values in a row
    if (npdepths >= 2) then
        do x = 1,preps(1)*preps(2)
            do i = 1,nrep_MC
                dekcheck = ck_pset(2:,i,x) >= ck_pset(:npdepths-1,i,x) !check if any settlement increases
                if(any(dekcheck)) then  !if so, then sort the whole vector in decreasing order
                    do k = npdepths,2,-1
                        do j = npdepths-1,npdepths-k+1,-1
                            if ( ck_pset(j+1,i,x) > ck_pset(j,i,x)) then
                                temp = ck_pset(j,i,x)
                                ck_pset(j,i,x) = ck_pset(j+1,i,x)
                                ck_pset(j+1,i,x) = temp
                            end if
                        end do
                    end do
                end if
            end do
        end do
    end if

    !if (ck_pset(1,i,x) <= ck_pset(2,i,x) ) then
    !    !ck_pset(1,i,x) = ck_pset(1,i,x) + (ck_pset(2,i,x) - ck_pset(1,i,x))/100000              !ck_pset(2,i,x) + tiny(0.0)
    !    temp = ck_pset(1,i,x)
    !    ck_pset(1,i,x) = ck_pset(2,i,x)
    !    ck_pset(2,i,x) = temp
    !end if
    !do pd = 2,npdepths-1
    !    if (ck_pset(pd,i,x) <= ck_pset(pd+1,i,x) ) then ! ck_pset(:,i,x)
    !        !ck_pset(pd,i,x) = ck_pset(pd-1,i,x) + (pdepths(pd-1)-pdepths(pd))*(ck_pset(pd-1,i,x)-ck_pset(pd+1,i,x))/(pdepths(pd-1)-pdepths(pd+1))
    !        temp = ck_pset(pd,i,x)
    !        ck_pset(pd,i,x) = ck_pset(pd+1,i,x)
    !        ck_pset(pd+1,i,x) = temp
    !    end if
    !end do
    !if (ck_pset(npdepths-1,i,x) <= ck_pset(npdepths,i,x) ) then
    !    !ck_pset(npdepths,i,x) = ck_pset(npdepths,i,x) - (ck_pset(npdepths,i,x)-ck_pset(npdepths-1,i,x))/100000        !ck_pset(npdepths-1,i,x) - tiny(0.0)
    !    temp = ck_pset(npdepths-1,i,x)
    !    ck_pset(npdepths-1,i,x) = ck_pset(npdepths,i,x)
    !    ck_pset(npdepths,i,x) = temp
    !end if


    !restore original soil properties for SI analysis


    !output results to files
    write(*,*) 'Saving settlement to disk'
    do x = 1,preps(1)*preps(2)
        if (singletrue) then
            write(str2,'(A,A,I0,A,I0,A,F4.2,A,I0,A,I0,A,I0,A)') trim(datafolder),'ck_pile-',x,'_prad-',prad(1),'_esize-',dz,'_sof-',nint(soilth(1)),'_cov-',nint(100*sdata(1,2)/sdata(1,1)),'_anis-',anisotropy,'.txt'
        else
            write(str2,'(A,A,I0,A,I0,A,F4.2,A,I0,A,I0,A,I0,A)') trim(datafolder),'ck_pile-',x,'_prad-',prad(1),'_esize-',dz,'_nlayers-',nlayer,'_ratio2l-',nint(lmean_ave(2)/lmean_ave(1)),'_depth2l-',nint(dz*ldepths(1)),'.txt'
        end if
        open(500,file=str2)
        do i = 1,nrep_MC
            write(500,'(100000(E14.7,X))') (ck_pset(pd,i,x),pd=1,npdepths)
            !'(1000000(E11.6,X))'
        end do
        close(500)
    end do


    end subroutine





    !Get the true depths of the CK layers at the pile boundaries.
    !If specified, generate all layer boundaries and save into a big array
    !This subroutine must be called before process_si_multi since it sets up some arrays and values
    subroutine prepmultiSI(npl2,goodpiles,goodcases,sdist,sumweights,extents,indices,CKheight,preps,load_con,usepie,radius,soilseeds,prad,nrep_MC,plocation)


    integer, intent(in) :: preps(2) !number of piles in each dimension
    integer, intent(in) :: load_con(preps(1)*preps(2))			!connectivity vector saying which pile has which load
    logical, intent(in) :: usepie !true to use FEM approximation such as the PIE method or M&K method, otherwise use FEM directly
    integer, intent(in) :: radius
    integer, intent(in) :: soilseeds(nrep_MC),prad(2),nrep_MC,plocation(preps(1)*preps(2),2)

    !effective pile layer stuff
    real, allocatable, intent(out) :: sdist(:,:,:),sumweights(:) !distance of elements from each pile, sum of distance weights
    integer, allocatable, intent(out) :: extents(:,:,:), indices(:,:,:) !extent of soil around pile, x/y coords of layer boundary indices

    integer, intent(in) :: goodpiles(:)        !vector of piles that have an associated applied load
    integer, intent(in) :: goodcases                               !number of good piles in goodpiles
    integer, intent(out) :: npl2 !good number of piles under certain conditions

    real, allocatable, intent(out) :: CKheight(:,:,:) !effective layer depths at each pile associated with the full, original soil

    real, allocatable :: soilweight_2D(:,:) ! 2D cross section of axisymmetric soil weights around pile
    real, allocatable :: sweights(:,:) ! soil weights around pile in plan view, according to axisymmetric FEA


    !local variables
    integer :: i,j,iter,counter
    real :: randu
    real, parameter ::  power=1.5
    real start,ender
    integer tempradius !instead of using the max pile length for the radius used for the weighted average calculation, use 5*pile diameter
    real(8) :: xyi(2,nxew*nyew)

    real bfldave(nxew,nyew)
    real efld2D(nxew,nyew)
    real sdata_temp(4)



    !tempradius = 5*(prad(1)+prad(2))/2
    tempradius = nint(5.0/dz)  ! use 5 metres

    allocate(indices(2,nxew,nyew))
    allocate(CKheight(nrep_MC,nlayer-1,preps(1)*preps(2)))
    allocate(CKprops(nrep_MC,nlayer,preps(1)*preps(2)))
    !allocate(goodpiles(preps(1)*preps(2)))



    if (usepie) then
        npl2 = goodcases
    else
        npl2 = preps(1)*preps(2)
    end if

    allocate(sdist(npl2,nxew,nyew), sumweights(npl2),extents(npl2,2,2))

    !get index values
    counter=0

    do i=1,nyew
        do j=1,nxew
            counter = counter+1
            xyi(1,counter) = j
            xyi(2,counter) = i
            indices(1,j,i) = j
            indices(2,j,i) = i
        end do
    end do


    !call cpu_time(allstart)

    !---generate soil weights around pile in a 2D axisymmetric manner. This code isn't finished! so it's commented out ---
    !allocate(soilweight_2D(tempradius,nzew),sweights(tempradius*2+prad(1),tempradius*2+prad(2)))
    !call getweights_2d(tempradius,nzew,(prad(1)+prad(2))/2,nzew/2,dz,soilweight_2D)
    !call make_plan_weights(soilweight_2D,sweights)


    do j = 1,npl2
        !Get the upper and lower bounds of soil indices in each dimension. I.e. a square of height information around the pile.
        extents(j,1,1) = plocation(goodpiles(j),1) - tempradius !xlow
        extents(j,2,1) = plocation(goodpiles(j),2) - tempradius !ylow
        extents(j,1,2) = plocation(goodpiles(j),1) + tempradius + prad(1) - 1 !xhigh
        extents(j,2,2) = plocation(goodpiles(j),2) + tempradius + prad(2) - 1  !yhigh

        sdist(j,:,:) = sqrt((indices(1,:,:)-plocation(goodpiles(j),1)-float(prad(1))/4)**2 + (indices(2,:,:)-plocation(goodpiles(j),2)-float(prad(2))/4)**2)  !Calculate the distances of each soil element from the pile using pythagoras
        sdist(j,:,:) = (sdist(j,:,:)*dz)**(-power) !Multiply them by an arbitrary power (negative for inverse).
        sdist(j,plocation(goodpiles(j),1):plocation(goodpiles(j),1)+prad(1)-1,plocation(goodpiles(j),2):plocation(goodpiles(j),2)+prad(2)-1) = 0 !elements within the pile radius must be near-zero as they are likely replaced by the pile itself
        sdist(j,plocation(goodpiles(j),1):plocation(goodpiles(j),1)+prad(1)-1,plocation(goodpiles(j),2):plocation(goodpiles(j),2)+prad(2)-1) = maxval(sdist)/(10 * prad(1)*prad(2))  !Set the weighting within the pile to be 10% of the maximum outside the pile.
        sumweights(j) = sum(sdist(j,extents(j,1,1):extents(j,1,2),extents(j,2,1):extents(j,2,2))) !get the sum of the weights as part of the weighted average calculation
    end do




    !if all the soil layers are to be stored in-memory, allocate the array.
    !This also pre-processes the CK layer depths at pile locations
    if (superset) then
        allocate(fullbfld(nrep_MC,nlayer-1,nxew,nyew))


        if (.not. var_props) then
            do i = 1,nlayer !if the soil properties aren't presented by a variable random field, then set the CK properties at the deterministic value
                CKprops(:,i,:) = lmean_ln(i)
            end do
        end if

            write(*,*) 'Pre-prcessing soil layers'

            !call cpu_time(start)

                    !     	bfldave = 0

        !loop through Monte Carlo realisations
        do iter = 1,nrep_MC
            !call cpu_time(start)
            if((iter-int(10*iter/nrep_MC)*nrep_MC/10 == 0)) write(*,'(I0,A,X)',advance='no') iter*100/nrep_MC,'%'


            !---Generate virtual soil---
            kseed = randu(soilseeds(iter)) * 1234567890

            !Emulate variable soil properties through a 2D random field (hybrid soil model approach)
            if (var_props) then
                do i = 1,nlayer
                    sdata_temp(3) = lmean_ln(i)
                    sdata_temp(4) = lsd_ln(i)
                    ! generate logarithm of 2D random field for Young's modulus (assumes that it's lognormally-distributed)
                    call sim2d(efld2D,nxe,nye,nxew,nyew,5*nye/4,dz,dz,kseed,MXM,MXK, &
                    bC0(2,:),bCT(2,:,:),bCC(2,:,:,:),bCS(2,:,:,:),bCI(2,:,:),bAT(2,:,:,:),bAC(2,:,:,:,:),bAS(2,:,:,:,:),bAI(2,:,:,:), bM, bk1, bk2, bkk,sdata_temp,distribution)
                    do j = 1,goodcases
                    ! calculate weighted geometric average of soil properties around pile
                        CKprops(iter,i,goodpiles(j)) = exp(sum(efld2D(extents(j,1,1):extents(j,1,2),extents(j,2,1):extents(j,2,2)) * sdist(j,extents(j,1,1):extents(j,1,2),extents(j,2,1):extents(j,2,2))) / sumweights(j))
                    end do
                end do
            end if

            !Store the soil in memory
            call soil_layers(xyi)                             !get layer boundaries
            fullbfld(iter,:,:,:) = bfld(2:nlayer,:,:)


        !             bfldave = bfldave + bfld(2,:,:)

            !Get layer heights at each pile for true, full, original soil (later used to determine actual differential settlement)
            if(multitype == 1) then !just use the layer boundaries within the pile
                do j = 1,goodcases
                    do i = 1,nlayer-1

                        CKheight(iter,i,goodpiles(j)) = dz * sum(bfld(i+1,plocation(goodpiles(j),1):plocation(goodpiles(j),1)+prad(1)-1,plocation(goodpiles(j),2):plocation(goodpiles(j),2)+prad(2)-1))/product(prad)

                    end do
                end do

            else !get the effective heights for each layer based on inverse distance weighted average of surrounding soil
                do j = 1,goodcases
                    do i = 1,nlayer-1
                        CKheight(iter,i,goodpiles(j)) = dz * sum(bfld(i+1,extents(j,1,1):extents(j,1,2),extents(j,2,1):extents(j,2,2)) * sdist(j,extents(j,1,1):extents(j,1,2),extents(j,2,1):extents(j,2,2))) / sumweights(j)
                    end do
                end do
            end if
        end do

        write(*,*)  !add newline

    end if

        !call CPU_TIME(ender)
        !write(*,*) start-ender
            !read(*,*)


            ! 		open(71805,file='bfldave.txt')
        ! 		do j = 1,nyew
        ! 			write(71805,'(1000000(F6.2,X))') bfldave(:,j)/nrep_MC
        ! 		end do
        !
        ! 		close(71805)
        ! 		stop
        !

        !open(505,file='truelayers_3L.dat',access='stream')
        !write(505) fullbfld(:,1:2,:,:)
        !close(505)
        !
        !stop



    end subroutine

    !This subroutine is deprecated, but leave it in for the sake of other deprecated features
    subroutine getheights(bfld,indices,plocation,xdim,ydim,npl,pdim,power,radius,nlayer,height)

    !get soil layer depths at individual piles. Specifically, the 'effective' layer depths,
    !calculated as an inverse-distance weighted average, with distances to an arbitrary power.
    !This reflects how soil further away from a pile has a smaller impact on its performance.
    !A contant height for the layers at each pile is required for both the 2D axisymmetric analysis
    !as well as the Mylonakis and Gazetas method for pile settlement.


    implicit none

    integer, intent(in) :: xdim,ydim,npl,power
    real, intent(in) :: bfld(:,:,:) !(nlayer-1,xdim,ydim)
    integer, intent(in) :: indices(2,xdim,ydim)
    integer, intent(in) :: plocation(npl,2),radius(2)

    integer, intent(in) :: pdim(2)
    integer, intent(in) :: nlayer

    integer :: xlow,xhigh,ylow,yhigh
    real :: dist(xdim,ydim)
    integer :: p,i,j
    real sumweights

    real, intent(out) :: height(nlayer-1,npl)

    real start,finish




    !Loop through the piles
    do p = 1,npl
        !Get the upper and lower bounds of soil indices in each dimension. I.e. a square of height information around the pile.

        call CPU_TIME(start)
        !do j = 1,1000

        xlow = plocation(p,1) - radius(1)
        ylow = plocation(p,2) - radius(2)
        xhigh= plocation(p,1) + radius(1) + pdim(1) - 1
        yhigh= plocation(p,2) + radius(2) + pdim(2) - 1

        dist = sqrt((indices(1,:,:)-plocation(p,1)-float(pdim(1))/4)**2 + (indices(2,:,:)-plocation(p,2)-float(pdim(2))/4)**2)  !Calculate the distances of each soil element from the pile using pythagoras
        dist = dist**(-power) !Multiply them by an arbitrary power (negative for inverse).
        dist(plocation(p,1):plocation(p,1)+pdim(1)-1,plocation(p,2):plocation(p,2)+pdim(2)-1) = maxval(dist)/10 !The elements within the pile radius must near-zero, as they are replaced by the pile itself.
        sumweights = sum(dist(xlow:xhigh,ylow:yhigh)) !get the sum of the weights as part of the weighted average calculation

        !end do
        call cpu_time(finish)
        !write(*,*) 'prep: ',(finish - start)/1000

        call CPU_TIME(start)
        !do j = 1,1000

        do i = 1,nlayer-1 !get the effective heights for each layer.
            height(i,p) = sum(bfld(i+1,xlow:xhigh,ylow:yhigh) * dist(xlow:xhigh,ylow:yhigh)) / sumweights
        end do

        !end do
        call cpu_time(finish)
        !write(*,*) 'height: ',(finish - start)/1000

        !read(*,*)

    end do

    end subroutine








    
    



end module




