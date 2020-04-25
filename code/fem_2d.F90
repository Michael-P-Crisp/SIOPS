
!
!   This code has been adapted to serve as an axisymmetric pile for the 
!   site investigation analysis. Right now both direct stiffness and PCG
!   solution algorithms have been coded in, but the PCG code is giving
!   incorrect results, so the direct stiffness method must be used.
!   They have been found to take the same time to run, so the choice is
!   fairly inconsequential.

!   Note that the pile width is changed to account for it being a circle
!   instead of a square (the volume of concrete must be identical)

 !type_2d = 'plane' or 'axisymmetric'
 !element= 'quadrilateral' or 'triangle'
!    
!-------------------------------------------------------------------------
! Program 5.1 Plane or axisymmetric strain analysis of an elastic solid
!             using 3-, 6-, 10- or 15-node right-angled triangles or
!             4-, 8- or 9-node rectangular quadrilaterals. Mesh numbered
!             in x(r)- or y(z)- direction.
!-------------------------------------------------------------------------
subroutine femaxi(xnew,nxe,nye,g_coord,g_num,g_g,evals,vval,type_2d,element,nels,nn,nf,loads,neq,nod,ndof) !, &
!type_2d='axisymmetric',element= 'quadrilateral',fixed_freedoms=0,node,sense,no,value)


 USE fem_main 
 USE fem_geom 
 !use fem_stuff
 
 IMPLICIT NONE
 
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15),nip=9,nodof=2,ndim=2
 INTEGER::i,iel,k,loaded_nodes,neq,nlen,nn,np_types,nr,nst=3
 REAL(iwp)::det,one=1.0_iwp,penalty=1.0e20_iwp,zero=0.0_iwp
 integer, intent(in) :: nod,ndof,nels,nxe,nye
 CHARACTER(LEN=15), intent(in) :: type_2d
 CHARACTER(LEN=20), intent(in) :: element
 integer, intent(in) :: g_num(nod,nels),nf(nodof,nn),g_g(ndof,nels)
 real, intent(in) :: g_coord(ndim,nn)
 real :: loads(0:neq)
 real, intent(out) :: xnew(0:neq)	!output displacement
 real(8),allocatable :: bee(:,:),dee(:,:)
 !-------------uncomment if you want fixed freedoms-------------
 !integer, parameter :: fixed_freedoms = 0
 !integer:: node(fixed_freedoms),sense(fixed_freedoms),no(fixed_freedoms)                   
 !real(8) :: value(fixed_freedoms)
!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::g(:),kdiag(:),nfrepac(:,:,:),num(:)    
 REAL(8),ALLOCATABLE::coord(:,:),der(:,:),deriv(:,:), &
   eld(:),fun(:),gc(:),jac(:,:),      &
   points(:,:),prop(:,:),sigma(:),weights(:),        &
   x(:),diag_precon(:),d(:),store(:)  
 real, allocatable :: storkm(:,:,:),u(:),p(:),km(:,:),kv(:)
 real pi
 logical solid,cg_converged,direct
 REAL(iwp) :: alpha,beta,cg_tol,up
 integer cg_iters,cg_limit,count
  integer istat
 !logical, intent(in) :: debug
 real(4) :: evals(nels),vval
 integer j !loop counter
 character(1000) fname,dumploc
 integer :: ntol,npri=1,nstep=1
!-----------------------input and initialisation--------------------------

nlen=6;solid=.true.

 IF(type_2d=='axisymmetric') nst=4
 
 ALLOCATE(g(ndof),fun(nod),coord(nod,ndim),jac(ndim,ndim),der(ndim,nod),          &
   deriv(ndim,nod),km(ndof,ndof),eld(ndof),num(nod),gc(ndim),sigma(nst))
 ALLOCATE(kdiag(neq)) 
 allocate(bee(nst,ndof),dee(nst,nst))
 !allocate(storkm(ndof,ndof,nels),nfrepac(2,nxe+1,nze+1))

 
 kdiag=0
 !fixed_freedoms=0
  ALLOCATE(points(nip,ndim),weights(nip))

 
!############################################ use stiffness matrix approach ######################################
     
!-----------------------loop the elements to find global arrays sizes-----
 elements_1: DO iel=1,nels
   CALL fkdiag(kdiag,g_g(:,iel))
 END DO elements_1
 !CALL mesh(g_coord,g_num,argv,nlen,12)
 DO i=2,neq 
   kdiag(i)=kdiag(i)+kdiag(i-1) 
 END DO 
 ALLOCATE(kv(kdiag(neq)))
 !WRITE(11,'(2(A,I12))')                                                    &
 !  " There are",neq," equations and the skyline storage is ",kdiag(neq)
!-----------------------element stiffness integration and assembly--------
 !write(*,*) nip
CALL sample(element,points,weights) 
 kv=zero 
 gc=one
 elements_2: DO iel=1,nels
   CALL deemat(dee,evals(iel),vval)
   num=g_num(:,iel)
   coord=TRANSPOSE(g_coord(:,num)) 
   g=g_g(:,iel) 
   km=zero
   int_pts_1: DO i=1,nip
     CALL shape_fun(fun,points,i) 
     CALL shape_der(der,points,i)
     jac=MATMUL(der,coord) 
     det=determinant(jac) 
     CALL invert(jac)
     deriv=MATMUL(jac,der) 
     CALL beemat(bee,deriv)
     IF(type_2d=='axisymmetric')THEN
       gc=MATMUL(fun,coord) 
       bee(4,1:ndof-1:2)=fun(:)/gc(1)
     END IF
     km=km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)*gc(1)
   END DO int_pts_1
   CALL fsparv(kv,km,g,kdiag)
 END DO elements_2   
 

 
!uncomment if you want to have fixed freedoms
! IF(fixed_freedoms/=0)THEN
!   DO i=1,fixed_freedoms 
!    no(i)=nf(sense(i),node(i)) 
!   END DO
!   kv(kdiag(no))=kv(kdiag(no))+penalty 
!   loads(no)=kv(kdiag(no))*value
! END IF
!-----------------------equation solution---------------------------------
 CALL sparin(kv,kdiag) 
 CALL spabac(kv,loads,kdiag) 
 


 loads(0)=zero
 xnew=loads

	return !send loads (displacements) back
 

END subroutine
    
    
    
   
