!c  *********************************************************************
!c  *                                                                   *
!c  *                         subroutine soilstat                        *
!c  *                                                                   *
!c  *********************************************************************
!c  Single Precision Version 1.0
!c  Written by Michael P Crisp
!c  Latest Update: March, 2017
!c
!c  PURPOSE calculate various statistical moments and autocorrelations for the soil
!c
!c  DESCRIPTION: Get means and standard devation of final soil.
!   calculate autocorrelations for various lags. This is defined in 
!   metres rather than elements. As this can be a time-consuming process
!   for high-res soils, it's much faster to do it at discrete steps rather
!   than continuously. 
!
!   WARNING:
!
!    Note that if the calculations aren't continuous there will
!   be a slight error in the autocorrelation, as the mean and standard deviation is
!   is based on the full soil rather than the subset that's being analysed. 
!   This may change in the future.
!
!   Also, ensure that the step size in metres in each direction isn't smaller than
!   the dimension of the element.
!c
!c
!c  REVISION HISTORY:
!c-------------------------------------------------------------------------
      module soilstat


    contains

    !calculate max differential settlement between any two piles (function of distance and ddmax)
    subroutine outputstats(nxew,nyew,nzew,dx,dy,dz,ckcfldl,cordist,corstep)
    
    implicit none
    
    integer count,count2,xlength,ylength,zlength,nel,i,x,y,z,step
    integer, intent(in) :: nxew,nyew,nzew
    real mean,sd,var
    real :: ckcfldl(:,:,:)
    real(8), intent(in) :: dx,dy,dz
    real corx(0:nxew),cory(0:nyew),corz(0:nzew)
    real cordist, corstep !maximum lag to calculate correlation, and lag increments (m)
    
		open(1713,file='soilstats.txt')
    
     	corx=0; cory=0; corz=0
        corx(0)=1; cory(0)=1; corz(0)=1
         xlength=floor(cordist/dx)-1
         ylength=floor(cordist/dy)-1
         zlength=floor(cordist/dz)-1
         nel=nxew*nyew*nzew
         mean=sum(ckcfldl)/nel
         ckcfldl = ckcfldl - mean
         sd=sqrt(sum((ckcfldl)**2)/nel)
         !z correlation
         var=sd**2*nxew*nyew
         count=0 ; step=floor(corstep/dz)
         do i=step,zlength,step
             count=count+1; count2=0
             do z=1,nzew-i,step
                 count2=count2+1
                 corz(count)=sum((ckcfldl(:,:,z))*(ckcfldl(:,:,z+i)))+corz(count)
             end do
             corz(count)=corz(count)/(var*count2)
         end do
         write(1713,'(8096(F6.4,1X))') corz(0:count)
         !write(*,*) 'done z'
         var=sd**2*nxew*nzew
         count=0 ; step=floor(corstep/dy)
         do i=step,ylength,step
             count=count+1; count2=0
             do y=1,nyew-i,step
                 count2=count2+1
                 cory(count)=sum((ckcfldl(:,y,:))*(ckcfldl(:,y+i,:)))+cory(count)
             end do
             cory(count)=cory(count)/(var*count2)
         end do
         write(1713,'(8096(F6.4,1X))') cory(0:count)
         !write(*,*) 'done y'
         var=sd**2*nzew*nyew
         count=0 ; step=floor(corstep/dx)
         do i=step,xlength,step
             count=count+1; count2=0
             do x=1,nxew-i
                 count2=count2+1
                 corx(count)=sum((ckcfldl(x,:,:))*(ckcfldl(x+i,:,:)))+corx(count)
             end do
             corx(count)=corx(count)/(var*count2)
         end do
         write(1713,'(8096(F6.4,1X))') corx(0:count)
         !write(*,*) 'done x'
          ckcfldl = ckcfldl + mean
        
             nel=nxew*nyew*nzew
             write(1709,'(10(F8.3,1X))') mean,exp(sum(log(dble(ckcfldl)))/nel),sd !mean, geometric mean, sd
             
             !cycle !move onto next realisation
             
         close(1713)

      end subroutine

    end module
