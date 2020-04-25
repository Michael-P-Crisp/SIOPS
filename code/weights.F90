
!c  *********************************************************************
!c  *                                                                   *
!c  *                         subroutine getweights                     *
!c  *                                                                   *
!c  *********************************************************************
!c  Mixed Precision Version 1.01
!c  Written by D. Vaughan Griffiths and Gordon A. Fenton, Apr 9, 2011
!c  Latest Update: Jun 14, 2011
!c
!c  PURPOSE  Get soil weights used for PIE method using finite element analysis.
!c
!c  DESCRIPTION
!c
!   This routine performs a linear finite element analysis of the settlement
!   of a 3-D soil mass with a centrally-located, rigid footing. The soil weights
!	can then replace FEM completely when used with variable single-layer soil profiles.
!	See the fem_3d subroutine for a description of the finite element analysis code used.

!	The output is a file giving a list of normalised foundation settlements in mm:
!	1 kN applied load, Young's modulus 1 MPa. The settlement can then subsequently be
!	linearly scaled by both Young's modulus (inversely-proportional) and applied load (proportional).

!	Also, an additional file is produced for each specified foundation depth containing a list
!	of the soil weights for that particular pile design, most likely it's going in the order of (x, z, y), 
!	i.e y changing fastest and x changing slowest.


!	For more information on the PIE method, see:
!	"Effective Young's modulus of a spatially variable soil mass under a footing"
!	by Ching et al. (2018)

!	Program assembled by Michael Crisp, August 2018
    
module weights

contains


!c-------------------------------------------------------------------------

subroutine getweights(esize,srad,nzef,prad,npdepths,pdepths,cg_tol,cg_limit,startstress,datafolder)

 use fem_stuff
 USE fem_geom
 use fem_main
 use variables

 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15),npri=1,nstep=1,nod=8,nip=8
 INTEGER::fixed_freedoms,i,iel,k,loaded_nodes,   &
   nels,neq,nlen,nn,nr,nxef,nyef,nzef,ndepth,npl,ntol
 integer,parameter:: ndof=24,nip1=1,nprops=2,nodof=3,nst=6,ndim=3
 integer j !loop counter
 real(8) :: etime
 integer ::  nxp(1),nyp(1),pelm,pmod(1),tn,it,jr,ploop,pd
 REAL(8)::alpha,beta,det,dtim,one=1.0_iwp,penalty=1.0e20_iwp,up,      &
   zero=0.0_iwp
   CHARACTER(LEN=100)::argv,element='hexahedron'
   LOGICAL::cg_converged,solid=.TRUE. !whether the FEM solver has converged. Element solidity.

   integer, intent(in) :: prad(2) !width of pile in 2D (elements)
   integer,intent(in) :: npdepths !number of pile depths to test
   real :: nxpf,nypf	!float versions of the pile radius for checking
   
   integer, intent(in) :: srad !soil radius around pile

   real(8) :: coord(nod,ndim) !local element coordinates
   real(8),intent(in) :: cg_tol !convergence tolerance
   integer,intent(in) :: cg_limit !max number of iterations, 
   integer cg_iters !current number of iterations
   real(8), intent(in) :: esize !size of elements; length of cube side
   CHARACTER(LEN=15) ::type_2d !2d element type. Not used here.
   real start,finish
   integer xi,yi,zi,counter

!-----------------------dynamic arrays------------------------------------
 INTEGER, allocatable :: g(:),g_g(:,:),g_num(:,:),nf(:,:),num(:) !arrays used for the FEM
 REAL(4),ALLOCATABLE :: loads(:),prop(:,:),x(:),xnew(:)	!loads = applied load, xnew = displacement

 REAL(4),ALLOCATABLE :: soilweight(:)	!soil weights used for PIE method.
 real(4), allocatable :: sw(:,:,:,:) !3D array for the above data
 !real(4), allocatable :: sigma_z(:,:,:,:) !3D array of vertical stress

 real maxstress !maximum vertical stress in field
 integer, allocatable, intent(out) :: startstress(:,:)

 real(4),ALLOCATABLE :: evals(:),vvals(:) !Young's modulus and poisson's ratio for soil field.
 !real(8) g_pmul(ndof,nels),g_utemp(ndof,nels)
 real(4),ALLOCATABLE :: g_coord(:,:)	!global node coordinates
 integer :: istat
 logical :: regmesh	!whether the mesh has a single element size. Allows for lower RAM usage than otherwise.
integer, intent(in) :: pdepths(npdepths)	!pile depths to analyse (elements)
integer, allocatable :: tiednode(:) !equation number of the pile degree of freedom. Used for getting/setting displacement/load.
 REAL(8),ALLOCATABLE :: h_coords(:),z_coords(:)	!horizontal and vertical mesh coordinates
 real,ALLOCATABLE :: detdisp(:) !displacement for all pile settlement depths
 real(8),ALLOCATABLE :: sigma(:,:),strain(:,:) !stress and strain for the soilweight calculation

!-----------------------input and initialisation--------------------------
 character(1000) str2
 character(1000), intent(in) :: datafolder




regmesh = .true. !all elements must be same size for soilweight determination; regular mesh
 npl = 1			!one pile

!soil radius is pile position
nxp = srad
nyp = srad
!width of mesh in x and y directions. Note; code assumes pile is square for now. Non-square support is incomplete currently. Goes off x data.
nxef = 2 * srad + prad(1) 
nyef = 2 * srad + prad(2)

!do some checks for bogus input
if(nxef <= 0 .or. nzef <= 0 .or. prad(1) <= 0 .or. prad(2) <=0 .or. npdepths <= 0) then
	open(700,file='error.txt')
	write(700,*) 'Error, make sure that the following variables are integers greater than zero:'
	write(700,*) 'nxef, nzef, prad (x and y), npdepths'
	write(700,*) "You have written:",nxef,nzef,prad,npdepths
	close(700)
	stop
end if
if(cg_tol <= 0 .or. cg_tol >= 1) then
	open(700,file='error.txt')
	write(700,*) "Error, make sure that cg_tol is between 0 and 1 exclusive."
	write(700,*) "Recommend 0.0005"
	write(700,*) "You have written:",cg_tol
	close(700)
	stop
end if
if(cg_limit <= 0) then
	open(700,file='error.txt')
	write(700,*) 'Error, cg_limit must be an integer greater than 0'
	close(700)
	stop
end if
do i = 1,npdepths
	if(pdepths(i) < 0 .or. pdepths(i) > nzef) then
		open(700,file='error.txt')
		write(700,*) "Error, foundation depths must be positive integers less than the soil depth (elements)"
		write(700,*) "You have given pile depth and soil depth as:", pdepths(i),nzef
		close(700)
		stop
	end if
end do




nels = nxef*nyef*nzef
nn = ((nxef+1)*(nyef+1)*(nzef+1))







 
 allocate(evals(nels),vvals(nels),detdisp(npdepths),sw(npdepths,nxef,nyef,nzef)) !,sigma_z(6,nxef,nyef,nzef))
 allocate(g(ndof),g_g(ndof,nels),g_num(nod,nels),nf(nodof,nn),num(nod),g_coord(ndim,nn))
 allocate(sigma(nst,nels), strain(nst,nels),soilweight(nels)) !nst = no. stresses, e.g x,y,z,xy,xz,yz
 allocate(tiednode(npl))
 allocate(h_coords(nxef+1),z_coords(nzef+1))
 
 allocate(startstress(3,npdepths))
 startstress = 0

 vvals= 0.3 !set poisson's ratio constant for all elements. Hardcode for now.
 evals = 1.0 !unit Young's modulus. Can scale settlement afterwards according to pile force and soil stiffness.

h_coords(1)=0
z_coords(1)=0

!get the coordinate vectors
   do i=1,nxef
      h_coords(i+1)=esize + h_coords(i)
  end do
  
  do i=1,nzef
      z_coords(i+1)=esize + z_coords(i)
  end do

!do analysis once for each specified pile depth
write(*,*) 'Analysing Pile Designs'
write(*,*) 'Length (m)   settlement (mm)   time (s)'
do pd = 1,npdepths
     call cpu_time(start)

     !set boundary restraints
     call bc3d(nxef,nyef,nzef,nf,nodof,nn)


     !set equation numbers and sort out tied pile nodes
     call nfsort(nf,nxef,nyef,nzef,nn,tiednode,nxp,nyp,npl,prad,pdepths(pd))
     
     !number of equations
      neq=MAXVAL(nf)
       ALLOCATE(loads(0:neq),x(0:neq),xnew(0:neq))


 	    !get coordinate matrix for nodes, as well as the connectivity matrix
        DO iel=1,nels
           CALL hexahedron_xz(iel,h_coords,h_coords,z_coords,coord,num)
           CALL num_to_g(num,nf,g)
           g_num(:,iel)=num
           g_coord(:,num)=TRANSPOSE(coord)
           g_g(:,iel)=g
       END DO

    loads = 0
    loads(nf(3,tiednode(1))) =  -1

    !conduct 3D finite element analysis with iterative conjugate gradient solver

    call femit(xnew,cg_iters,cg_converged, & !output
	    cg_tol,cg_limit,g_coord,g_num,g_g,evals,vvals,nels,element,nn,nf,loads,neq,nod,ndof,regmesh)

    
	    
	!Exit with an error if the FEM solver hasn't converged.
	!if(.not. cg_converged) then
	!	open(700,file='error.txt')
	!	write(700,*) "Error, Iterative FEM subroutine did not converge."
	!	write(700,*) "Number of iterations was:",cg_iters
	!	write(700,*) "Suggest:"
	!	write(700,*) "Decreasing the tolerance,"
	!	write(700,*) "Decreasing the model complexity, or"
	!	write(700,*) "Increasing the allowable number of iterations."
	!	write(700,*) "(The last option is usually best)."
	!	close(700)
	!	stop
	!end if 
		

    type_2d='ignore'

    !get stress and strain information
    call get_stress(evals,vvals,xnew,g_coord,g_g,g_num,neq,nels,nod,ndof,ndim,nst,nn,type_2d,element,sigma, strain)

    !calculate soil weight
    do i =1,nels
	    soilweight(i) = sigma(1,i)*strain(1,i) + sigma(2,i)*strain(2,i) + sigma(3,i)*strain(3,i) + &
	             2*sigma(4,i)*strain(4,i) + 2*sigma(5,i)*strain(5,i) + 2*sigma(6,i)*strain(6,i)
    end do
    soilweight = soilweight/maxval(soilweight) !normalize so that largest weight is 1

    detdisp(pd) = xnew(nf(3,tiednode(1)))



    !output pile settlement curve, and soil weights for each pile depth
            !repack into 3D x,y,z array
			counter = 0
			do yi = 1,srad*2+prad(2)
				do zi = 1,nzef
					do xi = 1,srad*2+prad(1)
						counter = counter+1
						sw(pd,xi,yi,zi) = soilweight(counter)
                        !sigma_z(:,xi,yi,zi) = sigma(:,counter)
					end do
				end do
            end do
            
            
    
    !write(str2,'(A,I0,A)') 'soilweight',pd,'.txt'
    !open(668,file=str2)
    !write(str2,'(A,I0,A)') 'sigma_z',pd,'.txt'
    !open(6688,file=str2)
    !do xi = 1,nxef
    !    do yi = 1,nyef
    !        do zi = 1,nzef
    !            write(668,*) sw(xi,yi,zi)
	   !         write(6688,*) (sigma_z(j,xi,yi,zi),j=1,6)
    !        end do
    !    end do
    !end do
    !close(668)
    !close(6688)
    

    maxstress = maxval(sw(pd,:,:,:))
    
    do xi = srad,1,-1
        startstress(1,pd) = startstress(1,pd) + 1
        if(maxval(sw(pd,xi,:,:))/maxstress < 0.0001) exit
    end do
    
    do yi = srad,1,-1
        startstress(2,pd) = startstress(2,pd) + 1
        if(maxval(sw(pd,:,yi,:))/maxstress < 0.0001) exit
    end do
    
    startstress(3,pd) = pdepths(pd)
    do zi = pdepths(pd)+1,nzef
        startstress(3,pd) = startstress(3,pd) + 1
        if(maxval(sw(pd,:,:,zi))/maxstress < 0.0001) exit
    end do
    
    
    DEALLOCATE(loads,x,xnew)
    
   
    
    call cpu_time(finish)
     write(*,'(F10.3,3X,F15.10,3X,F8.1)') pdepths(pd)*esize,detdisp(pd),finish-start
    

end do



    write(str2,'(A,A,I0,A,F4.2,A)') trim(datafolder),'soilweights_prad-',prad(1),'_esize-',dz,'.dat'
    open(668,file=str2,access='stream')
	write(668) sw(:,srad-maxval(startstress(1,:))+1:srad+prad(1)+maxval(startstress(1,:)),srad-maxval(startstress(2,:))+1:srad+prad(2)+maxval(startstress(2,:)),:maxval(startstress(3,:)))
    close(668)
    


!save settlement curve to disk
write(str2,'(A,A,I0,A,F4.2,A)') trim(datafolder),'settlement_prad-',prad(1),'_esize-',dz,'.txt'
open(667,file=str2)
write(667,*) 'Pile depth (m), Normalised settlement (MPa*mm/kN)'
do i = 1,npdepths
	write(667,*) esize*pdepths(i),detdisp(i)
end do
close(667)

write(str2,'(A,A,I0,A,F4.2,A)') trim(datafolder),'model_bounds_prad-',prad(1),'_esize-',dz,'.txt'
open(667,file=str2)
write(667,*) 'Pile depth (m), Normalised settlement (MPa*mm/kN)'
do i = 1,npdepths
	write(667,*) startstress(:,i)
end do
close(667)


end subroutine

end module