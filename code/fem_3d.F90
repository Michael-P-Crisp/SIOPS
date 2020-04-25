!c  *********************************************************************
!c  *                                                                   *
!c  *                         subroutine fem3ds                         *
!c  *                                                                   *
!c  *********************************************************************
!c  Mixed Precision Version 1.01
!c  Written by D. Vaughan Griffiths and Gordon A. Fenton, Apr 9, 2011
!c  Latest Update: Jun 14, 2011
!c
!c  PURPOSE  finite element settlement analysis of one or two piles in a
!c           3-D soil mass.
!c
!c  DESCRIPTION!c  *********************************************************************
!c  *                                                                   *
!c  *                         subroutine fem3ds                         *
!c  *                                                                   *
!c  *********************************************************************
!c  Mixed Precision Version 1.01
!c  Written by D. Vaughan Griffiths and Gordon A. Fenton, Apr 9, 2011
!c  Latest Update: Jun 14, 2011
!c
!c  PURPOSE  finite element settlement analysis of one or two piles in a
!c           3-D soil mass.
!c
!c  DESCRIPTION
!c
!c  This routine performs a linear finite element analysis of the settlement
!c  of a 3-D soil mass in which one or two piles are founded. The analysis
!c  proceeds by using a Conjugate Gradient Iterative solver that avoids the
!c  need to assemble the entire stiffness matrix (a 50-element cube
!c  problem requires a stiffness matrix of about 50^5 for a total memory
!c  requirement in double precision of about 8*50^5 = 2.5 Gbytes). At the
!c  moment, 8-node brick elements are being used.
!c
!c  The soil mass has the following coordinate system;
!c
!c
!c                        --------------------------|
!c                  z    /                         /|
!c                  ^   /                         / |
!c                  |  y                         /  |
!c                  | /                         /   |
!c                  |/                         /    |
!c                  |-------x->---------------|     |
!c                  |                         |     |
!c                  |                         |     |
!c                  |                         |     |
!c                  |                         |     /
!c                  |                         |    /
!c                  |                         |   /
!c                  |                         |  /
!c                  |                         | /
!c                  |                         |/
!c                  ---------------------------
!c
!c  ARGUMENTS
!c
!c   istat	unit number connected to the file to which information,
!c		warning, and error messages are to be printed. (input)
!c
!c   loads	real*8 array used for the finite element analysis
!c g_coord	real*8 array used for the finite element analysis
!c
!c      ym	real*8 array of size at least (nxe x nye x nze) to contain
!c		the double precision version of the elastic modulus field.
!c		(ym = efld). (output)
!c
!c    efld	real*4 array of size at least (nxe x nye x nze) which contains
!c		the (optionally) random (locally averaged) elastic modulus
!c		field. (input)
!c
!c    nrfx	integer giving the leading dimension of the efld array.
!c		(input)
!c
!c    nrfy	integer giving the second dimension of the efld array.
!c		(input)
!c
!c       v	real*8 value containing Poisson's ratio. (input)
!c
!c       p	real*8 array used for the finite element analysis
!c       r	real*8 array used for the finite element analysis
!c       x	real*8 array used for the finite element analysis
!c    xnew	real*8 array used for the finite element analysis
!c       u	real*8 array used for the finite element analysis
!c  precon	real*8 array used for the finite element analysis
!c       d	real*8 array used for the finite element analysis
!c
!c     nxe	integer giving the number of elements describing the soil
!c		mass in the x-direction (width). (input)
!c
!c     nye	integer giving the number of elements describing the soil
!c		mass in the y-direction (depth). (input)
!c
!c     nze	integer giving the number of elements describing the soil
!c		mass in the z-direction (height). (input)
!c
!c      dx	real*8 value giving the x-direction dimension of each element.
!c		(input)
!c
!c      dy	real*8 value giving the y-direction dimension of each element.
!c		(input)
!c
!c      dz	real*8 value giving the z-direction dimension of each element.
!c		(input)
!c
!c     npl	integer containing the number of piles supported by
!c		the soil mass. (input)
!c
!c     nxp	integer vector of length at least npl containing the
!c		x-node index of the near edge of each pile (counting from
!c		node 1 at x = 0 in the x-direction). Each pile is assumed
!c		to be one element in plan size. (input)
!c
!c     nyp	integer vector of length at least npl containing the
!c		y-node index of the near edge of each pile (counting from
!c		node 1 at y = 0 in the y-direction). (input)
!c
!c       f	real*8 vector of length at least npl containing the load
!c		applied to each footing. (input)
!c
!c  cg_tol	real*8 value containing the convergence tolerance for the
!c		conjugate gradient iterations. (input)
!c
!c   limit	integer value containing the maximum number of iterations
!c		allowed during the conjugate gradient solution of the
!c		system of equations. If this limit is hit, the solution
!c		did not converge. (input)
!c
!c      nf	integer array used for the finite element analysis
!c   kdiag	integer array used for the finite element analysis
!c   ntied	integer array used for the finite element analysis
!c   g_num	integer array used for the finite element analysis
!c     g_g	integer array used for the finite element analysis
!c
!c  iloads	maximum number of degrees-of-freedom in the system. (input)
!c     inf	maximum number of nodes in the system. (input)
!c    inel	maximum number of elements. (input)
!c     inx	maximum number of elements in the x-direction. (input)
!c     iny	maximum number of elements in the y-direction. (input)
!c     inz	maximum number of elements in the z-direction. (input)
!c
!c  MAXTIE	maximum number of tied nodes per pile (should be 4). (input)
!c
!c    disp	real*8 vector of length at least equal to the number of
!c		piles which on output will contain the displacement
!c		(positive downwards) of each pile. (output)
!c
!c   iters	integer containing the actual number of iterations taken
!c		by the Conjugate Gradient solver to converge. (output)
!c
!c   debug	logical flag which is true if debug information is to be
!c		sent to the *.stt file. (input)
!c
!c  PARAMETERS
!c
!c	MAXFT	maximum number of piles. This is used to set up temporary
!c		arrays in this routine but must be at least as large as MAXFT
!c		defined in mrpile3s.f
!c
!c  REVISION HISTORY:
!c  1.01	streamlined nsfix computation below slightly (Jun 14/11)
!c-------------------------------------------------------------------------



      !   Pile FEM subroutine mostly replaced by Program 5.6 of Programming the Finite
!   Element Method, 5th edition. This leads to more succinct code, variable element size,
!   and partial support for additional mesh settings.
    
!-------------------------------------------------------------------------
! Program 5.6 Three-dimensional strain of an elastic solid using
!             8-, 14- or 20-node brick hexahedra. Mesh numbered in x-z
!             planes then in the y-direction. No global stiffness matrix
!             assembly. Diagonally preconditioned conjugate gradient solver.
!-------------------------------------------------------------------------

!Note: This code includes options for meshes where element size and stiffness is variable
!	(most flexible case, but with higher memory consumption), a case where the size is constant
!	and stiffness is variable (allows reduced memory consumption), and a case where the element
!	size and stiffness is constant (most rigid case, allows low memory + faster computation).
!	These settings can be specified by setting the regmesh variable to 3,2 or 1 respectively.


subroutine femit(xnew,cg_iters,cg_converged, & !output
	cg_tol,cg_limit,g_coord,g_num,g_g,evals,vvals,nels,element,nn,nf,loads,neq,nod,ndof,regmesh)
          !femaxi(xnew,nxe,nye,g_coord,g_num,g_g,evals,vvals,nels,nn,nf,loads,neq,nod,ndof)

USE fem_main 
 USE fem_geom 
 !use fem_stuff
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15),npri=1,nstep=1,nip=8
 INTEGER::fixed_freedoms,i,iel,k,loaded_nodes,   &
   nlen,nr,ndepth,npl,ntol
 integer,parameter:: nip1=1,nodof=3,nst=6,ndim=3
 integer j !loop counter
 REAL(8)::alpha,beta,det,dtim,one=1.0_iwp,penalty=1.0e20_iwp,up,      &
   zero=0.0_iwp 
   CHARACTER(LEN=20)::element !'hexahedron'
   LOGICAL:: solid=.TRUE.,paraview
   real(8) conval !convergence value
   logical output
   real cumdiff !cumulative average difference
   real amax !historial maximum absolute displacement throughout iteration
   
   real start,finish !timing variables


!-------------------dynamic arrays---------------------- 
 REAL(8),ALLOCATABLE::points(:,:),weights(:)                                                        
 INTEGER,ALLOCATABLE:: no(:),sense(:),node(:)
 real(4), allocatable :: storkm(:,:,:),value(:),store(:)
 !real(8) g_pmul(ndof,nels),g_utemp(ndof,nels)
 
 
!---------------input arrays and values------------------ 
 integer, intent(in) :: nels,neq,nn,cg_limit,nod,ndof
 integer, intent(in) :: g_g(ndof,nels),g_num(nod,nels),nf(nodof,nn)
 real(4), intent(in) :: evals(nels),vvals(nels),g_coord(ndim,nn)
 real(4) :: loads(0:neq)
 logical, intent(in) :: regmesh
 real(8), intent(in) :: cg_tol !, cg_cumtol !conv. tol., cumulative average conv. tol.
 !integer, intent(in) :: tol_limit !number of iterations that must be consequtively below tolerance
 !real(4) :: tlvec(tol_limit)	  !vector storing convergences to check for tol_limit
 
!---------------output arrays and values------------------ 
 logical, intent(out) :: cg_converged
 real(4), intent(out) :: xnew(0:neq)
 !real(4), intent(out) :: maxset(0:cg_limit)
 !real(4), intent(out) :: maxdiff(0:cg_limit)
 integer, intent(out) :: cg_iters
 
 !----------static arrays--------------
 INTEGER::g(ndof),num(nod)
 REAL(8)::bee(nst,ndof),coord(nod,ndim),dee(nst,nst),der(ndim,nod),       &       
   deriv(ndim,nod),eld(ndof),fun(nod),gc(ndim),jac(ndim,ndim),sigma(nst)
 real(4) :: km(ndof,ndof),utemp(24),base_km(ndof,ndof)
 REAL(4) ::d(0:neq),diag_precon(0:neq),x(0:neq),p(0:neq),u(0:neq) 
 real, allocatable :: g_pmul(:,:),g_utemp(:,:)
 
 
 !-----------external functions-------------
 !real, external :: sdot

!-----------------------input and initialisation--------------------------

	!maxset = 0 !maximum settlement each realization (theoretically pile settlement)
	!cumdiff = 0 !set cumulative average difference to zero
	

	 ALLOCATE(points(nip,ndim),weights(nip))
	 

	 CALL sample(element,points,weights) 
	 diag_precon=zero
	 x=zero 
	 cg_iters = 0
	 !tlvec = 10
	 
	 !output = fname /= ''
	 
	 !if(output) open(666,file=fname)
	 

	 
	 	!----------element stiffness integration, storage and preconditioner------ 
	!if(.false.) then !elements have same size and material properties (this doesn't seem to work accurately, and is slower too
	!
	!		allocate(g_pmul(ndof,nels),g_utemp(ndof,nels))
	!
	!	   CALL deemat(dee,1.0,sngl(vvals(1)))
	!	   base_km=zero
	!	   DO i=1,nip
	!		 CALL shape_der(der,points,i)
	!		 jac=MATMUL(der,TRANSPOSE(g_coord(:,g_num(:,1)))) 
	!		 det=determinant(jac) 
	!		 CALL invert(jac)
	!		 deriv=MATMUL(jac,der) 
	!		 CALL beemat(bee,deriv)
	!		 base_km=base_km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
 !          END DO 
 !
 !          km = base_km * evals(1)
	!	 
	!	   !preconditioner 
	!	   DO iel=1,nels
	!		   DO k=1,ndof 
	!			 diag_precon(g_g(k,iel))=diag_precon(g_g(k,iel))+km(k,k)
	!		   END DO
	!	   END DO
	!	   
	!	   	 diag_precon(1:)=one/diag_precon(1:) 
	!		 diag_precon(0)=zero 	 
	!		 d=diag_precon*loads 
	!		 p=d 
 !            
	!	   
	!	   write(*,*) 'done'
	!	   
	!	   !low memory, fast mesh solving
	!	  DO
	!	   cg_iters=cg_iters+1 
	!	   u=zero 
	!	   g_utemp=zero
	!	   DO iel=1,nels 
	!		 g_pmul(:,iel)=p(g_g(:,iel))
	!	   END DO
 !          g_utemp = matmul(km,g_pmul)
	!	   !CALL ssymm('L','l',ndof,nels,one,km,ndof,g_pmul,ndof,one,g_utemp,ndof)
	!	   DO iel=1,nels 
	!		 u(g_g(:,iel))=u(g_g(:,iel))+g_utemp(:,iel)
	!	   END DO
	!	   !IF(fixed_freedoms/=0)u(no)=p(no)*store 
	!	   !up=SDOT(neq,loads,1,d,1)
	!	   !alpha=up/SDOT(neq,p,1,u,1) 
 !          up = dot_product(loads,d)
 !          alpha = up/dot_product(p,u)
	!	   !CALL saxpy(neq,alpha,p,1,xnew,1)
 !          xnew = alpha*p + x
	!	   !alpha=-alpha 
	!	   !CALL saxpy(neq,alpha,u,1,loads,1) 
 !          loads = loads - u * alpha
	!	   d=diag_precon*loads
	!	   !beta=SDOT(neq,loads,1,d,1)/up 
 !          beta = dot_product(loads,d)/up
	!	   p=d+p*beta
 !          
	!	   CALL checon(xnew,x,cg_tol,cg_converged)
 !          x=xnew
	!	   IF(cg_converged.OR.cg_iters==cg_limit)EXIT
	!	END DO  
    !end if

		

	if(regmesh) then	!elements have same size, different properties
	
	
		   CALL deemat(dee,1.0,sngl(vvals(1)))	!element stiffness integration single element (e=1)
		   base_km=zero
		   DO i=1,nip
			 CALL shape_der(der,points,i)
			 jac=MATMUL(der,TRANSPOSE(g_coord(:,g_num(:,1)))) 
			 det=determinant(jac) 
			 CALL invert(jac)
			 deriv=MATMUL(jac,der) 
			 CALL beemat(bee,deriv)
			 base_km=base_km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
		   END DO 

		 
		   !preconditioner 
		   DO iel=1,nels
			   DO k=1,ndof 
				 diag_precon(g_g(k,iel))=diag_precon(g_g(k,iel))+base_km(k,k)*evals(iel) 
			   END DO
		   END DO  
		   
			 diag_precon(1:)=one/diag_precon(1:) 
			 diag_precon(0)=zero 	 
			 d=diag_precon*loads 
			 p=d 
			 
			 
		 !call cpu_time(start)
		   
		 !low-memory regular mesh solving
		 DO 
		   cg_iters=cg_iters+1 
		   !write(*,*) cg_iters
		   u=zero
		   DO iel=1,nels
			 g=g_g(:,iel)  
			 u(g)=u(g)+MATMUL(base_km*evals(iel),p(g)) 
		   END DO 
		   !uncomment if you want fixed freedoms
		!   IF(fixed_freedoms/=0)u(no)=p(no)*store 
		   up=DOT_PRODUCT(loads,d)
		   alpha=up/DOT_PRODUCT(p,u) 
		   xnew=x+p*alpha 
		   loads=loads-u*alpha
		   d=diag_precon*loads
		   beta=DOT_PRODUCT(loads,d)/up
		   p=d+p*beta
		   
		   CALL checon(xnew,x,cg_tol,cg_converged)
		   !maxset(cg_iters) = maxval(abs(xnew))
		   !maxdiff(cg_iters) = MAXVAL(ABS(xnew-x))

		   !write(*,*) 100*cg_tol/conval
		   !if(output) write(666,*) maxval(xnew)
		   x=xnew
		   IF(cg_converged.OR.cg_iters==cg_limit)EXIT
           !call cpu_time(finish)
           !write(*,*) cg_iters,finish-start
         END DO 
    
         
	else !elements have different size and properties
		
		allocate(storkm(ndof,ndof,nels))
		
		 DO iel=1,nels
		   CALL deemat(dee,evals(iel),vvals(iel))
		   num=g_num(:,iel) 
		   g=g_g(:,iel) 
		   coord=TRANSPOSE(g_coord(:,num))
		   km=zero
		   gauss_pts_1: DO i=1,nip
			 CALL shape_der(der,points,i)
			 jac=MATMUL(der,coord) 
			 det=determinant(jac) 
			 CALL invert(jac)
			 deriv=MATMUL(jac,der) 
			 CALL beemat(bee,deriv)
			 km=km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
		   END DO gauss_pts_1
		   storkm(:,:,iel)=km
		   DO k=1,ndof 
			 diag_precon(g(k))=diag_precon(g(k))+km(k,k) 
		   END DO
		 END DO  
		 
			 diag_precon(1:)=one/diag_precon(1:) 
			 diag_precon(0)=zero 	 
			 d=diag_precon*loads 
			 p=d  
		 
		 !higher memory variable element size mesh
		pcg: DO 
		   cg_iters=cg_iters+1 
		   u=zero
		   elements_2: DO iel=1,nels
			 g=g_g(:,iel) 
			 !km= storkm(:,:,iel)
			 u(g)=u(g)+MATMUL(storkm(:,:,iel),p(g)) 
		   END DO elements_2
		   !uncomment if you want fixed freedoms
		!   IF(fixed_freedoms/=0)u(no)=p(no)*store 
		   up=DOT_PRODUCT(loads,d)
		   alpha=up/DOT_PRODUCT(p,u) 
		   xnew=x+p*alpha 
		   loads=loads-u*alpha
		   d=diag_precon*loads 
		   beta=DOT_PRODUCT(loads,d)/up 
		   p=d+p*beta
		   CALL checon(xnew,x,cg_tol,cg_converged)
		   !maxset(cg_iters) = maxval(xnew)
		   !maxdiff(cg_iters) = MAXVAL(ABS(xnew-x))
		   !write(*,*) amax
		   !write(*,*) 100*cg_tol/conval
		   !if(output) write(666,*) maxval(xnew)
		   x=xnew
		   IF(cg_converged.OR.cg_iters==cg_limit)EXIT
		 END DO pcg
		 
		 
! 	 else 
! 		write(*,*) 'Error, incorrect mesh property setting selected'
! 		write(*,*) '1 = constant shape and properties'
! 		write(*,*) '2 = constant shape, variable properties'
! 		write(*,*) '3 = variable shape and properties'
! 		stop
	 end if
	 
	 
	 !if(output) close(666)


	!uncomment for fixed freedoms, place in preconditioning code
	 !fixed_freedoms = 0
	 !IF(fixed_freedoms/=0)THEN
	 !  ALLOCATE(node(fixed_freedoms),sense(fixed_freedoms),                   &
	 !    value(fixed_freedoms),no(fixed_freedoms),store(fixed_freedoms))
	 !  READ(10,*)(node(i),sense(i),value(i),i=1,fixed_freedoms)
	 !  DO  i=1,fixed_freedoms 
	 !    no(i)=nf(sense(i),node(i)) 
	 !  END DO
	 !  diag_precon(no)=diag_precon(no)+penalty 
	 !  loads(no)=diag_precon(no)*value
	 !  store=diag_precon(no)
	 !END IF
	 
	!if(paraview) then
	! CALL mesh_ensi(argv,nlen,g_coord,g_num,element,evals,nf,dble(loads(1:)),      &
	!                nstep,npri,dtim,solid,atype,jr,ploop,slash)
	!end if



      return
 end subroutine
      


!c
!c  This routine performs a linear finite element analysis of the settlement
!c  of a 3-D soil mass in which one or two piles are founded. The analysis
!c  proceeds by using a Conjugate Gradient Iterative solver that avoids the
!c  need to assemble the entire stiffness matrix (a 50-element cube
!c  problem requires a stiffness matrix of about 50^5 for a total memory
!c  requirement in double precision of about 8*50^5 = 2.5 Gbytes). At the
!c  moment, 8-node brick elements are being used.
!c
!c  The soil mass has the following coordinate system;
!c
!c
!c                        --------------------------|
!c                  z    /                         /|
!c                  ^   /                         / |
!c                  |  y                         /  |
!c                  | /                         /   |
!c                  |/                         /    |
!c                  |-------x->---------------|     |
!c                  |                         |     |
!c                  |                         |     |
!c                  |                         |     |
!c                  |                         |     /
!c                  |                         |    /
!c                  |                         |   /
!c                  |                         |  /
!c                  |                         | /
!c                  |                         |/
!c                  ---------------------------
!c
!c  ARGUMENTS
!c
!c   istat	unit number connected to the file to which information,
!c		warning, and error messages are to be printed. (input)
!c
!c   loads	real*8 array used for the finite element analysis
!c g_coord	real*8 array used for the finite element analysis
!c
!c      ym	real*8 array of size at least (nxe x nye x nze) to contain
!c		the double precision version of the elastic modulus field.
!c		(ym = efld). (output)
!c
!c    efld	real*4 array of size at least (nxe x nye x nze) which contains
!c		the (optionally) random (locally averaged) elastic modulus
!c		field. (input)
!c
!c    nrfx	integer giving the leading dimension of the efld array.
!c		(input)
!c
!c    nrfy	integer giving the second dimension of the efld array.
!c		(input)
!c
!c       v	real*8 value containing Poisson's ratio. (input)
!c
!c       p	real*8 array used for the finite element analysis
!c       r	real*8 array used for the finite element analysis
!c       x	real*8 array used for the finite element analysis
!c    xnew	real*8 array used for the finite element analysis
!c       u	real*8 array used for the finite element analysis
!c  precon	real*8 array used for the finite element analysis
!c       d	real*8 array used for the finite element analysis
!c
!c     nxe	integer giving the number of elements describing the soil
!c		mass in the x-direction (width). (input)
!c
!c     nye	integer giving the number of elements describing the soil
!c		mass in the y-direction (depth). (input)
!c
!c     nze	integer giving the number of elements describing the soil
!c		mass in the z-direction (height). (input)
!c
!c      dx	real*8 value giving the x-direction dimension of each element.
!c		(input)
!c
!c      dy	real*8 value giving the y-direction dimension of each element.
!c		(input)
!c
!c      dz	real*8 value giving the z-direction dimension of each element.
!c		(input)
!c
!c     npl	integer containing the number of piles supported by
!c		the soil mass. (input)
!c
!c     nxp	integer vector of length at least npl containing the
!c		x-node index of the near edge of each pile (counting from
!c		node 1 at x = 0 in the x-direction). Each pile is assumed
!c		to be one element in plan size. (input)
!c
!c     nyp	integer vector of length at least npl containing the
!c		y-node index of the near edge of each pile (counting from
!c		node 1 at y = 0 in the y-direction). (input)
!c
!c       f	real*8 vector of length at least npl containing the load
!c		applied to each footing. (input)
!c
!c  cg_tol	real*8 value containing the convergence tolerance for the
!c		conjugate gradient iterations. (input)
!c
!c   limit	integer value containing the maximum number of iterations
!c		allowed during the conjugate gradient solution of the
!c		system of equations. If this limit is hit, the solution
!c		did not converge. (input)
!c
!c      nf	integer array used for the finite element analysis
!c   kdiag	integer array used for the finite element analysis
!c   ntied	integer array used for the finite element analysis
!c   g_num	integer array used for the finite element analysis
!c     g_g	integer array used for the finite element analysis
!c
!c  iloads	maximum number of degrees-of-freedom in the system. (input)
!c     inf	maximum number of nodes in the system. (input)
!c    inel	maximum number of elements. (input)
!c     inx	maximum number of elements in the x-direction. (input)
!c     iny	maximum number of elements in the y-direction. (input)
!c     inz	maximum number of elements in the z-direction. (input)
!c
!c  MAXTIE	maximum number of tied nodes per pile (should be 4). (input)
!c
!c    disp	real*8 vector of length at least equal to the number of
!c		piles which on output will contain the displacement
!c		(positive downwards) of each pile. (output)
!c
!c   iters	integer containing the actual number of iterations taken
!c		by the Conjugate Gradient solver to converge. (output)
!c
!c   debug	logical flag which is true if debug information is to be
!c		sent to the *.stt file. (input)
!c
!c  PARAMETERS
!c
!c	MAXFT	maximum number of piles. This is used to set up temporary
!c		arrays in this routine but must be at least as large as MAXFT
!c		defined in mrpile3s.f
!c
!c  REVISION HISTORY:
!c  1.01	streamlined nsfix computation below slightly (Jun 14/11)
!c-------------------------------------------------------------------------



      !   Pile FEM subroutine mostly replaced by Program 5.6 of Programming the Finite
!   Element Method, 5th edition. This leads to more succinct code, variable element size,
!   and partial support for additional mesh settings.
    
!-------------------------------------------------------------------------
! Program 5.6 Three-dimensional strain of an elastic solid using
!             8-, 14- or 20-node brick hexahedra. Mesh numbered in x-z
!             planes then in the y-direction. No global stiffness matrix
!             assembly. Diagonally preconditioned conjugate gradient solver.
!-------------------------------------------------------------------------

!Note: This code includes options for meshes where element size and stiffness is variable
!	(most flexible case, but with higher memory consumption), a case where the size is constant
!	and stiffness is variable (allows reduced memory consumption), and a case where the element
!	size and stiffness is constant (most rigid case, allows low memory + faster computation).
!	These settings can be specified by setting the regmesh variable to 3,2 or 1 respectively.


subroutine femitold(xnew,cg_iters,cg_converged, & !output
	cg_tol,cg_limit,g_coord,g_num,g_g,evals,vvals,nels,element,nn,nf,loads,neq,nod,ndof,regmesh)
          !femaxi(xnew,nxe,nye,g_coord,g_num,g_g,evals,vvals,nels,nn,nf,loads,neq,nod,ndof)

USE fem_main 
 USE fem_geom 
 !use fem_stuff
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15),npri=1,nstep=1,nip=8
 INTEGER::fixed_freedoms,i,iel,k,loaded_nodes,   &
   nlen,nr,ndepth,npl,ntol
 integer,parameter:: nip1=1,nodof=3,nst=6,ndim=3
 integer j !loop counter
 REAL(8)::alpha,beta,det,dtim,one=1.0_iwp,penalty=1.0e20_iwp,up,      &
   zero=0.0_iwp 
   CHARACTER(LEN=20)::element !'hexahedron'
   LOGICAL:: solid=.TRUE.,paraview
   real(8) conval !convergence value
   logical output
   real cumdiff !cumulative average difference
   real amax !historial maximum absolute displacement throughout iteration
   
   real start,finish !timing variables


!-------------------dynamic arrays---------------------- 
 REAL(8),ALLOCATABLE::points(:,:),weights(:)                                                        
 INTEGER,ALLOCATABLE:: no(:),sense(:),node(:)
 real(4), allocatable :: storkm(:,:,:),value(:),store(:)
 !real(8) g_pmul(ndof,nels),g_utemp(ndof,nels)
 
 
!---------------input arrays and values------------------ 
 integer, intent(in) :: nels,neq,nn,cg_limit,nod,ndof
 integer, intent(in) :: g_g(ndof,nels),g_num(nod,nels),nf(nodof,nn)
 real(4), intent(in) :: evals(nels),vvals(nels),g_coord(ndim,nn)
 real(4) :: loads(0:neq)
 logical, intent(in) :: regmesh
 real(8), intent(in) :: cg_tol !, cg_cumtol !conv. tol., cumulative average conv. tol.
 !integer, intent(in) :: tol_limit !number of iterations that must be consequtively below tolerance
 !real(4) :: tlvec(tol_limit)	  !vector storing convergences to check for tol_limit
 
!---------------output arrays and values------------------ 
 logical, intent(out) :: cg_converged
 real(4), intent(out) :: xnew(0:neq)
 !real(4), intent(out) :: maxset(0:cg_limit)
 !real(4), intent(out) :: maxdiff(0:cg_limit)
 integer, intent(out) :: cg_iters
 
 !----------static arrays--------------
 INTEGER::g(ndof),num(nod)
 REAL(8)::bee(nst,ndof),coord(nod,ndim),dee(nst,nst),der(ndim,nod),       &       
   deriv(ndim,nod),eld(ndof),fun(nod),gc(ndim),jac(ndim,ndim),sigma(nst)
 real(4) :: km(ndof,ndof),utemp(24),base_km(ndof,ndof)
 REAL(4) ::d(0:neq),diag_precon(0:neq),x(0:neq),p(0:neq),u(0:neq) 
 real, allocatable :: g_pmul(:,:),g_utemp(:,:)
 
 !-----------external functions-------------
 real, external :: sdot

!-----------------------input and initialisation--------------------------

	!maxset = 0 !maximum settlement each realization (theoretically pile settlement)
	!cumdiff = 0 !set cumulative average difference to zero
	

	 ALLOCATE(points(nip,ndim),weights(nip))
	 

	 CALL sample(element,points,weights) 
	 diag_precon=zero
	 x=zero 
	 cg_iters = 0
	 !tlvec = 10
	 
	 !output = fname /= ''
	 
	 !if(output) open(666,file=fname)
	 

	 
	 	!----------element stiffness integration, storage and preconditioner------ 
! 	if(regmesh==1) then !elements have same size and material properties
! 	
! 			allocate(g_pmul(ndof,nels),g_utemp(ndof,nels))
! 	
! 		   CALL deemat(dee,1.0,sngl(vvals(1)))
! 		   base_km=zero
! 		   DO i=1,nip
! 			 CALL shape_der(der,points,i)
! 			 jac=MATMUL(der,TRANSPOSE(g_coord(:,g_num(:,1)))) 
! 			 det=determinant(jac) 
! 			 CALL invert(jac)
! 			 deriv=MATMUL(jac,der) 
! 			 CALL beemat(bee,deriv)
! 			 base_km=base_km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
! 		   END DO 
! 
! 		 
! 		   !preconditioner 
! 		   DO iel=1,nels
! 			   DO k=1,ndof 
! 				 diag_precon(g_g(k,iel))=diag_precon(g_g(k,iel))+base_km(k,k)*evals(1) 
! 			   END DO
! 		   END DO
! 		   
! 		   	 diag_precon(1:)=one/diag_precon(1:) 
! 			 diag_precon(0)=zero 	 
! 			 d=diag_precon*loads 
! 			 p=d 
! 		   
! 		   write(*,*) 'done'
! 		   
! 		   !low memory, fast mesh solving
! 		  DO
! 		   cg_iters=cg_iters+1 
! 		   u=zero 
! 		   g_utemp=zero
! 		   write(*,*) 'done1'
! 		   DO iel=1,nels 
! 			 g_pmul(:,iel)=p(g_g(:,iel))
! 		   END DO
! 		   write(*,*) 'done2'
! 		   CALL ssymm('L','l',ndof,nels,one,km,ndof,g_pmul,ndof,one,g_utemp,ndof)
! 		   write(*,*) 'done3'
! 		   DO iel=1,nels 
! 			 u(g_g(:,iel))=u(g_g(:,iel))+g_utemp(:,iel)
! 		   END DO
! 		   write(*,*) 'done4'
! 		   !IF(fixed_freedoms/=0)u(no)=p(no)*store 
! 		   up=SDOT(neq,loads,1,d,1)
! 		   alpha=up/SDOT(neq,p,1,u,1) 
! 		   CALL saxpy(neq,alpha,p,1,xnew,1)
! 		   write(*,*) 'done5'
! 		   alpha=-alpha 
! 		   CALL saxpy(neq,alpha,u,1,loads,1) 
! 		   write(*,*) 'done6'
! 		   d=diag_precon*loads
! 		   beta=SDOT(neq,loads,1,d,1)/up 
! 		   p=d+p*beta
! 		   CALL checon(xnew,x,cumdiff,cg_tol,cg_cumtol,cg_converged,conval,cg_iters,tlvec)
! 		   maxset(cg_iters) = maxval(xnew)
! 		   if(output) write(*,*) 100*cg_tol/conval
! 		   IF(cg_converged.OR.cg_iters==cg_limit)EXIT
! 		END DO  
		

		

	if(regmesh) then	!elements have same size, different properties
	
	
		   CALL deemat(dee,1.0,sngl(vvals(1)))	!element stiffness integration single element (e=1)
		   base_km=zero
		   DO i=1,nip
			 CALL shape_der(der,points,i)
			 jac=MATMUL(der,TRANSPOSE(g_coord(:,g_num(:,1)))) 
			 det=determinant(jac) 
			 CALL invert(jac)
			 deriv=MATMUL(jac,der) 
			 CALL beemat(bee,deriv)
			 base_km=base_km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
		   END DO 

		 
		   !preconditioner 
		   DO iel=1,nels
			   DO k=1,ndof 
				 diag_precon(g_g(k,iel))=diag_precon(g_g(k,iel))+base_km(k,k)*evals(iel) 
			   END DO
		   END DO  
		   
			 diag_precon(1:)=one/diag_precon(1:) 
			 diag_precon(0)=zero 	 
			 d=diag_precon*loads 
			 p=d 
			 
			 
		 !call cpu_time(start)
		   
		 !low-memory regular mesh solving
		 DO 
		   cg_iters=cg_iters+1 
		   !write(*,*) cg_iters
		   u=zero
		   DO iel=1,nels
			 g=g_g(:,iel)  
			 u(g)=u(g)+MATMUL(base_km*evals(iel),p(g)) 
		   END DO 
		   !uncomment if you want fixed freedoms
		!   IF(fixed_freedoms/=0)u(no)=p(no)*store 
		   up=DOT_PRODUCT(loads,d)
		   alpha=up/DOT_PRODUCT(p,u) 
		   xnew=x+p*alpha 
		   loads=loads-u*alpha
		   d=diag_precon*loads
		   beta=DOT_PRODUCT(loads,d)/up
		   p=d+p*beta
		   
		   CALL checon(xnew,x,cg_tol,cg_converged)
		   !maxset(cg_iters) = maxval(abs(xnew))
		   !maxdiff(cg_iters) = MAXVAL(ABS(xnew-x))

		   !write(*,*) 100*cg_tol/conval
		   !if(output) write(666,*) maxval(xnew)
		   x=xnew
		   IF(cg_converged.OR.cg_iters==cg_limit)EXIT
           !call cpu_time(finish)
           !write(*,*) cg_iters,finish-start
		END DO 
    
	else !elements have different size and properties
		
		allocate(storkm(ndof,ndof,nels))
		
		 DO iel=1,nels
		   CALL deemat(dee,evals(iel),vvals(iel))
		   num=g_num(:,iel) 
		   g=g_g(:,iel) 
		   coord=TRANSPOSE(g_coord(:,num))
		   km=zero
		   gauss_pts_1: DO i=1,nip
			 CALL shape_der(der,points,i)
			 jac=MATMUL(der,coord) 
			 det=determinant(jac) 
			 CALL invert(jac)
			 deriv=MATMUL(jac,der) 
			 CALL beemat(bee,deriv)
			 km=km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
		   END DO gauss_pts_1
		   storkm(:,:,iel)=km
		   DO k=1,ndof 
			 diag_precon(g(k))=diag_precon(g(k))+km(k,k) 
		   END DO
		 END DO  
		 
			 diag_precon(1:)=one/diag_precon(1:) 
			 diag_precon(0)=zero 	 
			 d=diag_precon*loads 
			 p=d  
		 
		 !higher memory variable element size mesh
		pcg: DO 
		   cg_iters=cg_iters+1 
		   u=zero
		   elements_2: DO iel=1,nels
			 g=g_g(:,iel) 
			 !km= storkm(:,:,iel)
			 u(g)=u(g)+MATMUL(storkm(:,:,iel),p(g)) 
		   END DO elements_2
		   !uncomment if you want fixed freedoms
		!   IF(fixed_freedoms/=0)u(no)=p(no)*store 
		   up=DOT_PRODUCT(loads,d)
		   alpha=up/DOT_PRODUCT(p,u) 
		   xnew=x+p*alpha 
		   loads=loads-u*alpha
		   d=diag_precon*loads 
		   beta=DOT_PRODUCT(loads,d)/up 
		   p=d+p*beta
		   CALL checon(xnew,x,cg_tol,cg_converged)
		   !maxset(cg_iters) = maxval(xnew)
		   !maxdiff(cg_iters) = MAXVAL(ABS(xnew-x))
		   !write(*,*) amax
		   !write(*,*) 100*cg_tol/conval
		   !if(output) write(666,*) maxval(xnew)
		   x=xnew
		   IF(cg_converged.OR.cg_iters==cg_limit)EXIT
		 END DO pcg
		 
		 
! 	 else 
! 		write(*,*) 'Error, incorrect mesh property setting selected'
! 		write(*,*) '1 = constant shape and properties'
! 		write(*,*) '2 = constant shape, variable properties'
! 		write(*,*) '3 = variable shape and properties'
! 		stop
	 end if
	 
	 
	 !if(output) close(666)


	!uncomment for fixed freedoms, place in preconditioning code
	 !fixed_freedoms = 0
	 !IF(fixed_freedoms/=0)THEN
	 !  ALLOCATE(node(fixed_freedoms),sense(fixed_freedoms),                   &
	 !    value(fixed_freedoms),no(fixed_freedoms),store(fixed_freedoms))
	 !  READ(10,*)(node(i),sense(i),value(i),i=1,fixed_freedoms)
	 !  DO  i=1,fixed_freedoms 
	 !    no(i)=nf(sense(i),node(i)) 
	 !  END DO
	 !  diag_precon(no)=diag_precon(no)+penalty 
	 !  loads(no)=diag_precon(no)*value
	 !  store=diag_precon(no)
	 !END IF
	 
	!if(paraview) then
	! CALL mesh_ensi(argv,nlen,g_coord,g_num,element,evals,nf,dble(loads(1:)),      &
	!                nstep,npri,dtim,solid,atype,jr,ploop,slash)
	!end if



      return
 end subroutine
      

