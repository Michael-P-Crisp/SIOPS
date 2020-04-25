!c  *********************************************************************
!c  *                                                                   *
!c  *                         subroutine edvide                         *
!c  *                                                                   *
!c  *********************************************************************
!c  Single Precision Version 1.0
!c  Written by Michael P Crisp
!c  Latest Update: September 12, 2016
!c
!c  PURPOSE Divide the LAS field by a specified factor in order to reduce the number
!   of elements
!   
!	For SI, this generates the subset itself by merging in arrays of constant values
!	based on the fieldmerge subroutine.

!   NOTE: ckdivide subroutine is deprecated owing to the merging of elements being
!   handled by the variable element size averging code (see optimesh subroutine)
!c
!c-------------------------------------------------------------------------

    module edivide
    contains

	
	
    !create a subset of the complete knowledge field for FEM
    
      subroutine ckdivide(efld,tempfld,nxew,nyew,nzew,femfactor,redtype)
      
      real(4), intent(in) :: efld(:,:,:)
      real(4), intent(out) :: tempfld(:,:,:)
      integer, intent(in) :: redtype, femfactor
      real produc,prod
      integer fem
      
      
      integer, intent(in) :: nxew,nyew,nzew
      integer ix,iy,iz,ic,x,y,z
      real ffs
      

      
      !if(femfactor==1) then !if there's no scaling, just return the original
      !    tempfld=efld
      !    return
      !end if
      
      fem=femfactor-1
      
      if(redtype==1) then !standard arithmetic average
          ffs=femfactor**3
          x=0
	      do ix=1,nxew-femfactor+1,femfactor
              x=x+1
              y=0
              do iy=1,nyew-femfactor+1,femfactor
                  y=y+1
                  z=0
                  do iz=1,nzew-femfactor+1,femfactor
                      z=z+1
                      tempfld(x,y,z) = sum(efld(ix:ix+fem,iy:iy+fem,iz:iz+fem))/ffs
                  end do
              end do
          end do
      end if
      
      if(redtype==2) then !geometric average
          ffs=femfactor**3
          x=0
	      do ix=1,nxew-femfactor+1,femfactor
              x=x+1
              y=0
              do iy=1,nyew-femfactor+1,femfactor
                  y=y+1
                  z=0
                  do iz=1,nzew-femfactor+1,femfactor
                      z=z+1
                      prod = exp(sum(log(efld(ix:ix+fem,iy:iy+fem,iz:iz+fem)))/ffs)
                      tempfld(x,y,z) =  prod
                  end do
              end do
          end do
      end if
      !
      
      end subroutine

! %%%%%%%%%%%%%%% subroutine optimesh %%%%%%%%%%%%%%%%%%%%%%%%%
!   
!   Purpose: Averages the LAS field to fit into a FEM mesh with variable element
!   sizes.
!
!   Simultaneously repacks the materials into a vector
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
subroutine optimesh(efld,tempfld,redtype,fesh,fesv)
      
      real(4), intent(in) :: efld(:,:,:)
      real(4), intent(out) :: tempfld(:)
      integer, intent(in) :: redtype
      real produc,prod
      integer fem,f5
      
      integer, intent(in) :: fesh(:),fesv(:)
      integer ix,iy,iz,ic,x,y,z,i,j,k
      real ffs
      logical ck !true if complete knowledge run
      real vsum,weight
      integer nxe,nye,nze,count
      real(4) t1(size(tempfld))
      real(4) t2(size(tempfld))
      integer xL,xU,yL,yU,zL,zU !bounds
      integer nefh,nefv !number of elements in vertical and horizontal directions
      
      nxe=size(efld,1)
      nye=size(efld,2)
      nze=size(efld,3)
      
      nefh = size(fesh)
      nefv = size(fesv)
      
      tempfld=0.0
      t1=0.0
      t2=0.0
      f5=1074
      
      count=0     
      if(redtype==1) then !standard arithmetic average
	      do iy=1,nefh !number of elements in FEM horizontal
              do iz=1,nefv
                    do ix=1,nefh
                      count=count+1
                      xL=sum(fesh(1:ix-1));    xU=sum(fesh(1:ix))
                      yL=sum(fesh(1:iy-1));    yU=sum(fesh(1:iy))
                      zL=sum(fesv(1:iz-1));    zU=sum(fesv(1:iz))
                      weight=(xU-xL)*(yU-yL)*(zU-zL) !number of elements being averaged
                      tempfld(count) = sum(efld(xL+1:xU,yL+1:yU,zL+1:zU))/weight
                  end do
              end do
          end do
      else if(redtype==2)   then      !geometric average
	      do iy=1,nefh !number of elements in FEM horizontal
              do iz=1,nefv
                    do ix=1,nefh
                      count=count+1
                      xL=sum(fesh(1:ix-1));    xU=sum(fesh(1:ix))
                      yL=sum(fesh(1:iy-1));    yU=sum(fesh(1:iy))
                      zL=sum(fesv(1:iz-1));    zU=sum(fesv(1:iz))
                      weight=(xU-xL)*(yU-yL)*(zU-zL) !number of elements being averaged
                      tempfld(count) = exp(sum(log(efld(xL+1:xU,yL+1:yU,zL+1:zU)))/weight)
                  end do
              end do
          end do
      else  !average of geometric and arithmetic
	      do iy=1,nefh !number of elements in FEM horizontal
              do iz=1,nefv
                    do ix=1,nefh
                      count=count+1
                      xL=sum(fesh(1:ix-1));    xU=sum(fesh(1:ix))
                      yL=sum(fesh(1:iy-1));    yU=sum(fesh(1:iy))
                      zL=sum(fesv(1:iz-1));    zU=sum(fesv(1:iz))
                      weight=(xU-xL)*(yU-yL)*(zU-zL) !number of elements being averaged
                      t1(count) = exp(sum(log(efld(xL+1:xU,yL+1:yU,zL+1:zU)))/weight)
                      t2(count) = sum(efld(xL+1:xU,yL+1:yU,zL+1:zU))/weight
                  end do
              end do
          end do
          tempfld=sqrt(t1*t2)
      end if
      
      
end subroutine

!averages blocks of the mesh like the subroutine above, but fits it BACK
!into the original high-res mesh rather than an a coarser, optimized one.
!Used to test the effect of E averaging on error.
subroutine simpmesh(efld,tempfld,nefh,nefv,redtype,fesh,fesv,nxp,nyp,nxew,nyew,nzew,onlye)
      
      real(4), intent(inout) :: efld(:,:,:)
      real(4) :: tempfld(:)
      integer, intent(in) :: redtype,nxp(:),nyp(:)
      real produc,prod
      integer fem,f5
      
      integer, intent(in) :: fesh(:),fesv(:),nefh,nefv,nxew,nyew,nzew
      integer ix,iy,iz,ic,x,y,z,i,j,k
      real ffs
      logical ck !true if complete knowledge run
      real vsum,weight
      integer nxe,nye,nze,count
      real(4) t1,t2
      integer xL,xU,yL,yU,zL,zU
      logical onlye
      
      nxe=size(efld,1)
      nye=size(efld,2)
      nze=size(efld,3)
      
      tempfld=0.0
      t1=0.0
      t2=0.0
      f5=1074
      
      if(onlye) then
      count=0     
      if(redtype==1) then !standard arithmetic average
	      do iy=1,nefh !number of elements in FEM horizontal
              do iz=1,nefv
                    do ix=1,nefh
                      count=count+1
                      xL=sum(fesh(1:ix-1));    xU=sum(fesh(1:ix))
                      yL=sum(fesh(1:iy-1));    yU=sum(fesh(1:iy))
                      zL=sum(fesv(1:iz-1));    zU=sum(fesv(1:iz))
                      weight=(xU-xL)*(yU-yL)*(zU-zL) !number of elements being averaged
                      efld(xL+1:xU,yL+1:yU,zL+1:zU) = sum(efld(xL+1:xU,yL+1:yU,zL+1:zU))/weight
                  end do
              end do
          end do
      else if(redtype==2)   then      !geometric average
	      do iy=1,nefh !number of elements in FEM horizontal
              do iz=1,nefv
                    do ix=1,nefh
                      count=count+1
                      xL=sum(fesh(1:ix-1));    xU=sum(fesh(1:ix))
                      yL=sum(fesh(1:iy-1));    yU=sum(fesh(1:iy))
                      zL=sum(fesv(1:iz-1));    zU=sum(fesv(1:iz))
                      weight=(xU-xL)*(yU-yL)*(zU-zL) !number of elements being averaged
                      efld(xL+1:xU,yL+1:yU,zL+1:zU) = exp(sum(log(efld(xL+1:xU,yL+1:yU,zL+1:zU)))/weight)
                  end do
              end do
          end do
      else  !average of geometric and arithmetic
	      do iy=1,nefh !number of elements in FEM horizontal
              do iz=1,nefv
                    do ix=1,nefh
                      count=count+1
                      xL=sum(fesh(1:ix-1));    xU=sum(fesh(1:ix))
                      yL=sum(fesh(1:iy-1));    yU=sum(fesh(1:iy))
                      zL=sum(fesv(1:iz-1));    zU=sum(fesv(1:iz))
                      weight=(xU-xL)*(yU-yL)*(zU-zL) !number of elements being averaged
                      t1 = exp(sum(log(efld(xL+1:xU,yL+1:yU,zL+1:zU)))/weight)
                      t2 = sum(efld(xL+1:xU,yL+1:yU,zL+1:zU))/weight
                      efld(xL+1:xU,yL+1:yU,zL+1:zU)=sqrt(t1*t2)
                  end do
              end do
          end do
      end if
      end if
      
      
      count=0
      do iy=1,nyew
          do iz=1,nzew
              do ix=1,nxew
                  count=count+1
                  tempfld(count)=efld(ix,iy,iz)
              end do
          end do
      end do

end subroutine

!Provide mesh averaging for 2D axisymmetric case: assumes properties are uniform in the horizontal direction
subroutine optimesh2D(efld,tempfld,redtype,fesv,nefh)
      
      real(4), intent(in) :: efld(:,:)
      real(4), intent(out) :: tempfld(:)
      integer, intent(in) :: redtype
      real produc,prod
      integer fem,f5
      
      integer, intent(in) :: fesv(:)
      integer ix,iy,iz,ic,x,y,z,i,j,k
      real ffs
      logical ck !true if complete knowledge run
      real vsum,weight
      integer nxe,nye,nze,cou
      integer zL,zU !bounds
      integer nefh,nefv !number of elements in vertical and horizontal directions
      integer counter 

      
      nxe=size(efld,1)
      nze=size(efld,2)
      
      nefv = size(fesv)
      
      tempfld=0.0

      
      counter=0     
      if(redtype==1) then !standard arithmetic average
              do ix=1,nefh
                  do iz=1,nefv
                      zL=sum(fesv(1:iz-1));    zU=sum(fesv(1:iz))
                      weight=(zU-zL)!number of elements being averaged
                      vsum = sum(efld(1,zL+1:zU))/weight
                  
                      counter=counter+1
                      tempfld(counter) = vsum
                  end do
              end do
      else if(redtype==2)   then      !geometric average
              
                do ix=1,nefh
                    do iz=1,nefv
                      zL=sum(fesv(1:iz-1));    zU=sum(fesv(1:iz))
                      weight=(zU-zL) !number of elements being averaged
                      vsum = exp(sum(log(efld(1,zL+1:zU)))/weight)
                
                      counter=counter+1
                      tempfld(counter) = vsum
                  end do
              end do
      end if
      
      
end subroutine


      
      end module
      