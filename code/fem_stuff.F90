module fem_stuff
    
! This module contains various subroutines related to
! Programming the Finite Element Method, 5th edition 
! These serve to adapt the FEM subroutine specifically to
! the pile settlement problem.
    
    implicit none
    
    contains
    
    
    !convert an integer to a string
function itoa(i) result(res)
  character(100) :: res
  integer,intent(in) :: i
  character(range(i)+2) :: tmp
  write(tmp,'(i0)') i
  res = trim(tmp)
end function
   
  subroutine  rcoords( x_coords, y_coords, z_coords,nxp,nyp,nxe,nze,nesh,nesv)
  
  real(8), intent(out) :: x_coords(:),y_coords(:),z_coords(:)
  integer nh,nv,i,j,nxe,nze
  real(4) :: nesh(:),nesv(:)
  integer nxp(:),nyp(:)
  logical readvals
  
  readvals=.false.
  
  if(readvals) then
    !nesh= (/ (1.0,j=1,10),(0.5,j=1,41),(1.0,j=1,10) /)
    nesh= (/ 4.0,3.0,2.5,2.0,1.5,1.5,1.0,1.0,1.0,(0.5,j=1,11),1.0,1.0,1.0,1.5,1.5,2.0,2.5,3.0,4.0 /)
    !nesh= (/ (1.0,j=1,19),(0.5,j=1,5),(1.0,j=1,19) /)
    nesv= (/ (1.0,j=1,25),1.5,1.5,1.5,2.0,2.0,2.5,3.0,4.0,5.0,6.0,6.0 /)
    !nesv= (/ (1.0,j=1,60) /)
    !nesv= (/ (0.5,j=1,nv) /)
    !nesh= (/ (0.5,j=1,nh) /)
    nxp=nint(real(nh)/2); nyp=nxp
  end if
      
  nesh=nesh
  nesv=nesv

  x_coords(1)=0.0
  y_coords(1)=0.0
  z_coords(1)=0.0
  
  
  do i=1,nxe
      x_coords(i+1)=nesh(i)+x_coords(i)
      y_coords(i+1)=nesh(i)+y_coords(i)
  end do
  
  do i=1,nze
      z_coords(i+1)=nesv(i)+z_coords(i)
  end do
  
  end subroutine
  
  subroutine nfsort(nf,nxe,nye,nze,nn,nodenum,nxp,nyp,npl,prad,pdepth)
  
  !This sorts out the equation numbers, and ties the node
  !freedoms at the top of the pile.
  
  integer ic,i,L,j,ns,nn,nodenum(:),prad(2),nprad,npdepths,k
  integer nf(:,:)
  integer nyp(:),nxp(:),nfix(npl)
  integer, intent(in) :: nxe,nye,nze,npl
  integer, allocatable :: ntied(:,:)
  integer, intent(in) :: pdepth
  


    nfix = (prad(1)+1)*(prad(2)+1)*(pdepth + 1) !number of tied nodes per pile
    
      
    allocate(ntied(nfix(npl),npl))

        
      do L=1,npl
           ic = 0
         do i = nyp(L), nyp(L)+prad(2)
            ns = (i)*(nxe+1)*(nze+1)
            do j = nxp(L), nxp(L)+prad(1)
            	do k = 0,pdepth
               		ic = ic + 1
               		ntied(ic,L) = ns + (j + 1) + (nxe+1)*k 
               	end do
            end do
         end do
         nf(1:2,ntied(1,L)) = 0 !fix node in horizontal displacement
         do i = 2, nfix(L)			! and set corresponding
            nf(1:3,ntied(i,L)) = 0		! equation numbers to 0 for now
         end do
      end do
      
      

      
!c						number equations via nf-array
      ic = 0  !<------------------------------------------------------<<
      do i = 1, nn
         do j = 1, 3
            if( nf(j,i) .ne. 0 ) then
               ic = ic + 1
               nf(j,i) = ic
            endif
         end do
      end do
!c						and modify nf-array to account
!c						for "tied" freedoms
      do L = 1, npl  !<------------------------------------------------------<<
         do i=2,nfix(L)
             nf(3,ntied(i,L))=nf(3,ntied(1,L)) !tie vertical freedoms together
         end do
      end do   
         
      !save tied equation numbers to work with loads
        nodenum(1:npl)=ntied(1,1:npl)
         
  end subroutine
  
  !equation numbering for 2D axisymmetric case
   subroutine nf2d(nf,nxe,nye,nn,nodenum,prad,pdepth)
  
  integer ic,i,L,j,k,ns,nn,nsfix,nodenum
  integer prad,pdepth
  integer nf(:,:)
  integer nxp(1),nfix
  integer, intent(in) :: nxe,nye
  integer ntied((prad+1)*(pdepth + 1),1),npl

  
  L=1
  npl=1
  nxp=1
 
         nfix = (prad+1)*(pdepth + 1) !number of tied nodes per pile
         
         
         ic = 0
         do i = 0, prad
            do k = 0,pdepth
                ic = ic + 1
                ntied(ic,1) = (i)*(nye+1) + k  + 1
            end do
         end do
         nf(1,ntied(1,L)) = 0 !fix node in horizontal displacement nf(1,ntied(1,1))
         do i = 2, nfix			! and set corresponding
            nf(1:2,ntied(i,L)) = 0		! equation numbers to 0 for now
         end do

!c						number equations via nf-array
      ic = 0  !<------------------------------------------------------<<
      do i = 1, nn
         do j = 1, 2
            if( nf(j,i) .ne. 0 ) then
               ic = ic + 1
               nf(j,i) = ic
            endif
         end do
      end do
!c						and modify nf-array to account
!c						for "tied" freedoms

         do i=2,nfix
               nf(2,ntied(i,L))=nf(2,ntied(1,L))
         end do

      
     ! nf(2,559)
         
         
               !nfl(L) = nf(ntied(1,L),3)
      !loads(nfl(L)) =  -1250.0 !-f(L)
      nodenum=ntied(1,1)
      
         
         
  end subroutine
  
  subroutine bc3d(nxe,nye,nze,nfo,nodof,nn)
  
  !This subroutine sets the restraints on the sides and
  !bottom of the soil, for the 8 node brick case.
  !It shouldn't be too hard to write new code for the 14
  !and 20 cases.
  
      integer nfo(:,:),nxe,nye,nze,nodof,nn
      integer nm,j,k,i
      integer nf(nn,nodof)
      
      nm=0
      do j=1,nze
        nm=nm+1
        nf(nm,1)=0
        nf(nm,2)=0
        nf(nm,3) = 1
        do k=2,nxe
          nm=nm+1
          nf(nm,1) = 1
          nf(nm,2)=0
          nf(nm,3) = 1
        end do
        nm=nm+1
        nf(nm,1)=0
        nf(nm,2)=0
        nf(nm,3) = 1
    end do
      do  k=1,nxe+1
        nm=nm+1
        nf(nm,1)=0
        nf(nm,2)=0
        nf(nm,3)=0
      end do
      do i=2,nye
        do j=1,nze
          nm=nm+1
          nf(nm,1)=0
          nf(nm,2) = 1
          nf(nm,3) = 1
          do k=2,nxe
            nm=nm+1
            nf(nm,1) = 1
            nf(nm,2) = 1
            nf(nm,3) = 1
         end do
          nm=nm+1
          nf(nm,1)=0
          nf(nm,2) = 1
          nf(nm,3) = 1
        end do
        do k=1,nxe+1
          nm=nm+1
          nf(nm,1)=0
          nf(nm,2)=0
          nf(nm,3)=0
        end do
    end do
      do j=1,nze
        nm=nm+1
        nf(nm,1)=0
        nf(nm,2)=0
        nf(nm,3) = 1
        do k=2,nxe
          nm=nm+1
          nf(nm,1) = 1
          nf(nm,2)=0
          nf(nm,3) = 1
        end do
        nm=nm+1
        nf(nm,1)=0
        nf(nm,2)=0
        nf(nm,3) = 1
    end do
      do k=1,nxe+1
        nm=nm+1
        nf(nm,1)=0
        nf(nm,2)=0
        nf(nm,3)=0
      end do
      
      !the FEM subroutine from the book has this indexing
      !reversed, so swap the indexes.
      nfo(1,:)=nf(:,1)
      nfo(2,:)=nf(:,2)
      nfo(3,:)=nf(:,3)

      
      return
  end subroutine
  
  
    subroutine bc2d(nxe,nze,nfo,nn,nodof)
!c
!c      set up boundary condition data for 2D axisymmetric case
!c
      integer, intent(in) :: nxe,nze,nn,nodof
      integer nf(nn,nodof)
      integer, intent(out) :: nfo(:,:)
      integer nm,i,j,k
      
      nm=0
      nf=1
       !set boundaries to the left
        do j=1,nze
          nm=nm+1
          nf(nm,1) = 0
          nf(nm,2) = 1
        end do
        !set boundaries to the bottom
        do k=1,nxe+1
          nm=k*(nze+1)
          nf(nm,1) = 0
          nf(nm,2) = 0
        end do
        !set boundaries to the right
        nm=(nze+1)*nxe
        do j=1,nze
          nm=nm+1
          nf(nm,1) = 0
          nf(nm,2) = 1
        end do        

      nfo(1,:)=nf(:,1)
      nfo(2,:)=nf(:,2)
      
      return
      end subroutine
  
  
  !fill in pile concrete for a multi-element-thick pile
    
      subroutine setpil(efld,nxe,nye,nze,nxp,nyp,ndepth,pmod,prad,npl)
      
      implicit none
    
      real*4 efld(:), pmod
      integer nxp, nyp,npl,ndepth,nze,nye,nxe,prad
      integer i,kk,k,j,startloc,ix,iy,L
      !do L=1,npl
      
      !write(*,*) 'inside setpil'
      !read(*,*)
      
          do ix=1,prad
              do iy=1,prad
                  startloc=nze*nxe*(nyp+iy-1)+nxp+ix
      
                  do i = 1,ndepth
                      efld(startloc)=pmod
                      startloc=startloc+nxe
                  end do
              end do
          end do
      !end do
      
      
      return
      end subroutine
      
      subroutine setpil2d(efld,nxe,nye,nze,nxp,nyp,ndepth,pmod,prad,npl)
      
      implicit none
    
      real*4 efld(:), pmod
      integer nxp, nyp,npl,ndepth,nze,nye,nxe,prad
      integer i,kk,k,j,startloc,ix,iy,L
      !do L=1,npl
      
      !write(*,*) 'inside setpil'
      !read(*,*)
      
          do ix=1,prad
                  startloc=nye*(ix-1)
                  do i = 1,ndepth
                      efld(startloc+i)=pmod
                  end do
          end do
      !end do
      
      
      return
    end subroutine
    
 
 
end module

