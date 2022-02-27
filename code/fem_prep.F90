
!c  *********************************************************************
!c  *                                                                   *
!c  *                         subroutine fem3ds                         *
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
!   of a 3-D soil mass with a centrally-located, rigid footing.
!   Has the option of using variable element sizes (to allow for coarser resolution
!   away from the pile), to dramatically increase speed. LAS elements within a larger
!   FEM element will be geometrically averaged.

!	Program assembled by Michael Crisp, August 2018

!c-------------------------------------------------------------------------

module fem_prep

 use fem_stuff
 USE fem_geom
 use fem_main
 use edivide
 
 implicit none

contains

subroutine prep_fem3D(efld,esize,nxeo,nzeo,prad,pd,cg_tol,cg_limit,femvech,femvecv,disp,regmesh)

 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15),npri=1,nstep=1,nod=8,nip=8
 INTEGER::fixed_freedoms,i,iel,k,loaded_nodes,   &
   nels,neq,nlen,nn,nr,nxe,nze,ndepth,npl,ntol
 integer, intent(in) :: nxeo,nzeo
 integer,parameter:: ndof=24,nip1=1,nprops=2,nodof=3,nst=6,ndim=3
 integer j !loop counter
 real(8) :: etime
 integer ::  nxp(1),nyp(1),pelm,pmod(1),tn,it,jr,ploop
 REAL(8)::alpha,beta,det,dtim,one=1.0_iwp,penalty=1.0e20_iwp,up,      &
   zero=0.0_iwp
   CHARACTER(LEN=100)::argv,element='hexahedron'
   LOGICAL::cg_converged,solid=.TRUE. !whether the FEM solver has converged. Element solidity.

   integer :: prad(2) !width of pile in 2D (elements)
   real :: nxpf,nypf	!float versions of the pile radius for checking

   real(8) :: coord(nod,ndim) !local element coordinates
   real(8) :: cg_tol !convergence tolerance
   integer :: cg_limit,cg_iters !max number of iterations, current number of iterations
   real(8), intent(in) :: esize !size of elements; length of cube side
   CHARACTER(LEN=15) ::type_2d !2d element type. Not used here.
   
   integer, intent(in) :: pd
   real, intent(out) :: disp
   real, intent(in) :: efld(:,:,:)
   integer a,b,c
   integer count

   real, intent(in) :: femvech(:),femvecv(:) !size of FEM elements in horizontal and vertical directions
   integer :: snumh(size(femvech)),snumv(size(femvecv)) !number of soil elements per fem element in horizontal and vertical directions
!-----------------------dynamic arrays------------------------------------
 INTEGER, allocatable :: g(:),g_g(:,:),g_num(:,:),nf(:,:),num(:) !arrays used for the FEM
 REAL(4),ALLOCATABLE :: loads(:),prop(:,:),x(:),xnew(:)	!loads = applied load, xnew = displacement

 real(4),ALLOCATABLE :: evals(:),vvals(:) !Young's modulus and poisson's ratio for soil field.
 !real(8) g_pmul(ndof,nels),g_utemp(ndof,nels)
 real(4),ALLOCATABLE :: g_coord(:,:)	!global node coordinates
 integer :: istat
 logical,intent(in) :: regmesh	!whether the mesh has a single element size. Allows for lower RAM usage than otherwise.
integer, allocatable :: tiednode(:) !equation number of the pile degree of freedom. Used for getting/setting displacement/load.
 REAL(8),ALLOCATABLE :: h_coords(:),z_coords(:)	!horizontal and vertical mesh coordinates

!-----------------------input and initialisation--------------------------


!If regmesh is true, it will create a mesh where each soil element has its own fem element
!and use a lower-memory solver based on equal element size assumption
!If regmesh is false, it will create a varaible-element mesh and map the soil to it using
!a geometric average and the "femvec" vectors




 npl = 1			!one pile
 
if(regmesh) then
    nxe = nxeo
    nze = nzeo
else
    nxe = size(femvech)
    nze = size(femvecv)
end if
 

nxpf = (nxe-prad(1))/2
nypf = (nxe-prad(2))/2
nxp = nxpf
nyp = nypf



allocate(h_coords(nxe+1),z_coords(nze+1))
nels = nxe*nxe*nze
nn = ((nxe+1)*(nxe+1)*(nze+1))

 h_coords(1)=0
 z_coords(1)=0

 allocate(evals(nels),vvals(nels))
 allocate(g(ndof),g_g(ndof,nels),g_num(nod,nels),nf(nodof,nn),num(nod),g_coord(ndim,nn))
 allocate(tiednode(npl))
 
   





 vvals= 0.3 !set poisson's ratio constant for all elements. Hardcode for now.

 
    if (regmesh) then


          !get the coordinate vectors
		  do i=1,nxe
              h_coords(i+1)=esize + h_coords(i)
          end do
  
          do i=1,nze
              z_coords(i+1)=esize + z_coords(i)
          end do
          
          !repack soil
          count = 0
          do b = 1,nxe
          	do c = 1,nze
          		do a = 1,nxe
          			count = count + 1
          			evals(count) = efld(a,b,c)
          		end do
          	end do
          end do

          
    else
    
        !get number of soil elements per fem element from fem size vector
		snumh = (/ (nint(femvech(i)/esize), i=1,nxe) /)
		snumv = (/ (nint(femvecv(i)/esize), i=1,nze) /)
 
		!get young's modulus values for optimised mesh by geometric averaging
         call optimesh(efld,evals,2,snumh,snumv)
         
		 !get the coordinate vectors
		  do i=1,nxe
              h_coords(i+1)=femvech(i) + h_coords(i)
          end do
  
          do i=1,nze
              z_coords(i+1)=femvecv(i) + z_coords(i)
          end do
         
         
	end if


     !set boundary restraints
     call bc3d(nxe,nxe,nze,nf,nodof,nn,nod)


     !set equation numbers and sort out tied pile nodes
     call nfsort(nf,nxe,nxe,nze,nn,tiednode,nxp,nyp,npl,prad,pd, nod)
     
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
    
    write(*,*) 'here'

    !conduct 3D finite element analysis with iterative conjugate gradient solver
    call femit(xnew,cg_iters,cg_converged, & !output
	    cg_tol,cg_limit,g_coord,g_num,g_g,evals,vvals,nels,element,nn,nf,loads,neq,nod,ndof,regmesh)
	    
	!Exit with an error if the FEM solver hasn't converged.
	if(.not. cg_converged) then
		open(700,file='error.txt')
		write(700,*) "Error, Iterative FEM subroutine did not converge."
		write(700,*) "Number of iterations was:",cg_iters
		write(700,*) "Suggest:"
		write(700,*) "Decreasing the tolerance,"
		write(700,*) "Decreasing the model complexity, or"
		write(700,*) "Increasing the allowable number of iterations."
		write(700,*) "(The last option is usually best)."
		close(700)
	end if 


    disp = xnew(nf(3,tiednode(1)))

    
    DEALLOCATE(loads,x,xnew)


end subroutine




subroutine prep_fem2d(evals,vval,femvech,femvecv,pradin,disp,pd,esize)


INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15) !,nip=1,nodof=2,ndim=2,nod=8
 INTEGER::i,iel,k,loaded_nodes,neq,nlen,nn,np_types,nr
 REAL(iwp)::det
 real, parameter :: one=1.0_iwp,penalty=1.0e20_iwp,zero=0.0_iwp
 integer :: nels,nxe,nze
 integer,parameter:: nod=4,ndof=8,nip1=1,nprops=1,nodof=2,nst=3,ndim=2
 CHARACTER(LEN=15) :: type_2d
 CHARACTER(LEN=20) :: element
  CHARACTER(LEN=1) ::dir
  integer  prad!pile radius, original and axisymmetric version
  integer, intent(in) :: pradin
 integer, allocatable :: g_num(:,:),nf(:,:),g_g(:,:) !make this 
 real,allocatable :: g_coord(:,:)           !make this 
 real,allocatable :: loads(:),xnew(:)    !loads and output displacement
 !-------------uncomment if you want fixed freedoms-------------
 !integer, parameter :: fixed_freedoms = 0
 !integer:: node(fixed_freedoms),sense(fixed_freedoms),no(fixed_freedoms)                   
 !real(8) :: value(fixed_freedoms)
!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::g(:),num(:)   
 real(8) :: coord(nod,ndim)
 real(8), allocatable ::h_coords(:),z_coords(:)
  integer istat
 !logical, intent(in) :: debug
 real(4), intent(in) :: evals(:,:),vval
 real, allocatable :: tempfld(:)
 integer j !loop counter
 
 real, intent(in) :: femvech(:),femvecv(:)
 real, allocatable :: femvech_ax(:)
 integer snumv(size(femvecv)) !number of soil elements per FEM element in the vertical direction
 real, intent(out) :: disp !pile settlement
 integer pd !pile depth (fem elements)
 integer actwidth
 integer nodenum !node freedom associated with the force
 real(8), intent(in) :: esize


element = "quadrilateral"
type_2d='axisymmetric'
dir = 'z' 
nodenum = 1 !first node, vertical direction
 !vvals=v !set poisson's ratio constant for all elements



     !get various array sizes


     

    
    !get pile radius for axisymmetric mesh
if(mod(size(femvech),2)==0) then
    actwidth = size(femvech)/2
else
    actwidth = size(femvech)/2 + 1
end if
allocate(femvech_ax(actwidth))

if(mod(pradin,2) == 0) then
    prad = pradin/2 !If an even number, then use half
    femvech_ax = femvech(actwidth+1:)
else
    prad = pradin !otherwise, use the original amount and scale half of the first element
    femvech_ax = femvech(actwidth:)
    femvech_ax(1) = femvech_ax(1)/2

end if    
    nze = size(femvecv)
    nxe = size(femvech_ax)
    
    allocate(h_coords(nxe+1),z_coords(nze+1))
    nels = nxe*nze
    nn = ((nxe+1)*(nze+1))
    
    h_coords(1)=0
    z_coords(1)=0
    

 
    !get number of soil elements per fem element from fem size vector
    do i=1,nze
        snumv(i) = nint(femvecv(i)/esize)
    end do
		!snumv = (/ (nint(femvecv(i)/esize), i=1,nze) /)
 
		!get young's modulus values for optimised mesh by geometric averaging
        allocate(tempfld(nels))
         call optimesh2D(evals,tempfld,2,snumv,nxe)
		 !get the coordinate vectors
		  do i=1,nxe
              h_coords(i+1)=femvech_ax(i) + h_coords(i)
          end do
  
          do i=1,nze
              z_coords(i+1)=femvecv(i) + z_coords(i)
          end do
        
 
    
    
     allocate(g(ndof),g_g(ndof,nels),g_num(nod,nels),nf(nodof,nn),num(nod),g_coord(ndim,nn))

 !set equation numbers and sort out tied pile nodes
     !set boundary restraints
     call bc2d(nxe,nze,nf,nn,nodof)
     
     !get node numbering
     call nf2d(nf,nxe,nze,nn,nodenum,prad,pd)
     
    neq=MAXVAL(nf)
    ALLOCATE(loads(0:neq),xnew(0:neq))
    loads = 0
    loads(nodenum) =  0.1591549 !1/(2*pi) to account for circumference
    
    

     
    !get coordinate matrix for nodes, as well as the connectivity matrix
    DO iel=1,nels
        call geom_rect(element,iel,h_coords,z_coords,coord,num,dir)
        CALL num_to_g(num,nf,g)
        g_num(:,iel)=num
        g_coord(:,num)=TRANSPOSE(coord)
        g_g(:,iel)=g
    END DO
     
     
     call femaxi(xnew,nxe,nze,g_coord,g_num,g_g,tempfld,vval,type_2d,element,nels,nn,nf,loads,neq,nod,ndof)
     
     disp=xnew(nodenum)
     
     
     deallocate(femvech_ax,h_coords,z_coords,xnew,loads)
    deallocate(g,g_g,g_num,nf,num,g_coord,tempfld)


end subroutine

end module







