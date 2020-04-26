module piecewise_CMD
    
    !This module contains subroutines for a method for generating 3D random fields using
    !a piecewise covariance decomposition method. The advantage over LAS is that it allows
    !for accurate anisotropic soils. The disadvantage is that it has a streaked appearance,
    !although this streaking does not interfere with the soil statistics.

    
    !This code is adapted from the 1D random field covariance decomposition generator which
    !was in turn adapted from LAS2D. Note that any of the dimensions cannot exceed the size
    !of the size 0 LAS matrix, i.e. the MXK parameter, however this should be large enough.
    
    !The theory was developed by:
    !Li DQ, Xiao T, Zhang LM, & Cao ZJ. (2019). Stepwise covariance
    !    matrix decomposition for efficient simulation of multivariate
    !    large-scale three-dimensional random fields. Applied Mathematical
    !    Modelling, 68, 169-181. DOI: 10.1016/j.apm.2018.11.011
    
    contains
    
    
        !Initial processing 
      subroutine piecewise_init(nxe,nye,nze,dx,thh,thz,varfnc,R0x,R0y,R0z)

      !real cfld(*), c(7)
      real(8) thh, thz
      integer, intent(in) :: nxe
      integer kseed
      character(6) varfnc
      logical debug, liid, shofld, lxfld
      real*8 dvar, dpb, dx
      
      real(8), external :: dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      real*8 dmin1, dble
      real(4), intent(out) :: R0x(nxe,nxe), R0y(nye,nye), R0z(nze,nze)
      real(8) H, G ,da,db !for some of the more obscure correlation models
      !save XL, YL, ienter, liid
      !external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      !common/dparam/ dvar, dpb, dthx, dthy
      real(8) XL, YL, ZL
      
      data zero/0.0/, half/0.5/, one/1.0/, onept5/1.5/, twopt5/2.5/
      data ienter/0/, twopi/6.283185/
      data tol/1.e-3/

   1  format(a,a)

!c-------------------------------------- initialize ---------------------

         !liid   = ((thx .eq. zero) .and. (thy .eq. zero))

         XL     = float(nxe)*dx
         YL     = float(nye)*dx
         ZL     = float(nze)*dx

         !if( .not. liid ) then
         
         !c						initialize internal generator
      kseed = iseed( kseed )

            dvar = 1.d0
            if( varfnc .eq. 'dlafr2' .or. varfnc .eq. 'dlsfr2' ) then
               dpb = dx
               if( dx .lt. dx ) dpb = dx
            endif
            
            !if(min(thx,thy) < min(dx,dy) then
			!	write(*,*) 'warning, SOF is smaller than element. Will generate random noise.'
			!	return
				
			!else
            

				if( varfnc .eq. 'dlavx2' ) then
                   call piecewise_prep( dlavx2, nxe, nye, nze, XL, YL, ZL, R0x, R0y, R0z , tol, thh, thz )

				elseif( varfnc .eq. 'dlsep2' ) then
				   call piecewise_prep( dlsep2, nxe, nye, nze, XL, YL, ZL, R0x, R0y, R0z , tol, thh, thz )
                   
				elseif( varfnc .eq. 'dlspx2' ) then
				   call piecewise_prep( dlspx2, nxe, nye, nze, XL, YL, ZL, R0x, R0y, R0z , tol, thh, thz )
                   
				elseif( varfnc .eq. 'dlafr2' ) then
				   call piecewise_prep( dlafr2, nxe, nye, nze, XL, YL, ZL, R0x, R0y, R0z , tol, thh, thz )
                   
				elseif( varfnc .eq. 'dlsfr2' ) then
				   call piecewise_prep( dlsfr2, nxe, nye, nze, XL, YL, ZL, R0x, R0y, R0z , tol, thh, thz )
                   
				endif
            !end if
            

      
      return
      end
    

      
      !Build the covariance decomposition matrices in each direction
      
      subroutine piecewise_prep( dvarfn, N1x, N1y, N1z, XL, YL, ZL, Routx, Routy, Routz , tol, thx, thz )
      real*8 Rx(9,9,2), Bx(4,4), Sx(9,4)
      real*8 Ry(9,9,2), By(4,4), Sy(9,4)
      real*8 Rz(9,9,2), Bz(4,4), Sz(9,4)
      real*8 R0x(N1x*N1x), R0y(N1y*N1y), R0z(N1z*N1z)
      real*4 Routx(N1x,N1x), Routy(N1y,N1y), Routz(N1z,N1z)
      real(8), intent(in) :: XL,YL,ZL
      integer, intent(in) :: N1x, N1y, N1z
      real(8), intent(in) :: thx, thz !SOF in x and y directions
      real*8 T1, T2, dvarfn, dble
      logical lformR
      real(8) var, pb, H, G !field variance, averaging length, hurst components in horizontal and vertical directions
      real(8) da, db 
      external dvarfn
      integer x,y,z !loop counters
    
      !character dvarfn
   1  format(a,a,a)
   2  format(a,e13.4)
3     format(a,i4,a,i4,a,i4,a)

!c                    initialize internal generator
      !kseed = iseed( kseed ) !already initialized
!c                    form initial covariance matrix
      var = 1.d0
      R0x=0.d0
      R0y=0.d0
      R0z=0.d0
      
      !do x direction
      T1 = dble(XL)/dble(N1x)
      T2 = T1
      call dcvit2( dvarfn, R0x, N1x, Rx, 9, N1x, 1, T1, T2 , var, pb, H, G, da, db , thx, thx)
      

!c                    and compute its cholesky decomp
      call dchol2( R0x, N1x, N1x, rerr )
      if( rerr .gt. tol ) then
         write(*,1)'Warning: Cholesky decomposition of stage 0 covariance matrix'
         write(*,2)'         has maximum relative error of ',rerr
      endif

      
      
      
      
                  !do y direction
      T1 = dble(YL)/dble(N1y)
      T2 = T1
      call dcvit2( dvarfn, R0y, N1y, Ry, 9, N1y, 1, T1, T2 , var, pb, H, G, da, db , thx, thx)
      
      !c                    and compute its cholesky decomp
      call dchol2( R0y, N1y, N1y, rerr )
      if( rerr .gt. tol ) then
         write(*,1)'Warning: Cholesky decomposition of stage 0 covariance matrix'
         write(*,2)'         has maximum relative error of ',rerr
      endif

      
      
      
      
            !do z direction
      T1 = dble(ZL)/dble(N1z)
      T2 = T1
      call dcvit2( dvarfn, R0z, N1z, Rz, 9, N1z, 1, T1, T2 , var, pb, H, G, da, db , thz, thz)
      
      !c                    and compute its cholesky decomp
      call dchol2( R0z, N1z, N1z, rerr )
      if( rerr .gt. tol ) then
         write(*,1)'Warning: Cholesky decomposition of stage 0 covariance matrix'
         write(*,2)'         has maximum relative error of ',rerr
      endif


      
      !repack into single precision 2D arrays
      do x = 1,N1x
        Routx(x,:) = R0x((x-1)*N1x+1:x*N1x)
      end do
      do y = 1,N1y
        Routy(y,:) = R0y((y-1)*N1y+1:y*N1y)
      end do
      do z = 1,N1z
        Routz(z,:) = R0z((z-1)*N1z+1:z*N1z)
      end do

      return
      end subroutine
      
      
      
      subroutine matmul2(soil,s1,s2,R,r1,r2,U,u1,u2)
      
      integer, intent(in) :: s1,s2,r1,r2,u1,u2
      real, intent(in) :: R(r1,r2),U(u1,u2)
      real, intent(out) :: soil(s1,s2)
      
      soil = matmul(R,U)
      
      end subroutine
      
      
      
      
      !Generate the piecewise 3D random field through nested do loops.
      !This may be faster for very large soils due to hardware cache interaction.

    subroutine cmd_piecewise_loop( soil, nxe,nye,nze ,R0x, R0y, R0z)
      !use, intrinsic :: ISO_C_BINDING
      real(4), intent(out) :: soil(nxe,nye,nze)
      real(4) :: U(nxe,nye,nze)
      integer, intent(in) :: nxe,nye,nze
      real(4), intent(in) :: R0x(nxe,nxe), R0y(nye,nye), R0z(nze,nze)
      integer init
      integer L,i,j,x,y,z,x2 !loop counters
      real start, finish

!c-------------------------------------- generate realization -----------------
      iz = 0
!c                     generate field
      
      

      write(*,*) 'in loop'
      
      call cpu_time(start)
      
      ! get random numbers
      call vnorm( U, nxe*nye*nze )
    
    
      
    !process x direction

    soil = 0
    do x2 = 1,nxe
        do z = 1,nze
            do y = 1,nye
                do x = 1,nxe
                !call matmul2(soil2(:,y,z),nye,nze,Rx(x,:),nxe,U2(:,y,z),nye,nze)
                  soil(x,y,z) = soil(x,y,z) + R0x(x,x2) * U(x2,y,z)
              end do
          end do
        end do
    end do
      
      !process y direction
      U = 0
    do x2 = 1,nye
        do z = 1,nze
            do y = 1,nye
                do x = 1,nxe
                !call matmul2(soil2(:,y,z),nye,nze,Rx(x,:),nxe,U2(:,y,z),nye,nze)
                  U(x,y,z) = U(x,y,z) + R0y(y,x2) * soil(x,x2,z)
              end do
          end do
        end do
    end do
      
      !process z direction
      
    soil = 0
    
    do x2 = 1,nze
        do z = 1,nze
            do y = 1,nye
                do x = 1,nxe
                !call matmul2(soil2(:,y,z),nye,nze,Rx(x,:),nxe,U2(:,y,z),nye,nze)
                  soil(x,y,z) = soil(x,y,z) + R0z(z,x2) * U(x,y,x2)
              end do
          end do
        end do
    end do
      
     
      call cpu_time(finish)
      write(*,*) finish-start
      
    !open(505,file='piecewise.dat',access='stream')
    !write(505) soil
    !close(505)

    
      
    
    end subroutine
    
    
       !Generate the piecewise 3D random field through matmul.
      !This may be faster for smaller fields.
    
       subroutine cmd_piecewise_matmul( soil, nxe,nye,nze ,R0x, R0y, R0z)
      !use, intrinsic :: ISO_C_BINDING
      real(4), intent(out) :: soil(nxe,nye,nze)
      real(4) :: U(nxe,nye,nze)
      integer, intent(in) :: nxe,nye,nze
      real(4), intent(in) :: R0x(nxe,nxe), R0y(nye,nye), R0z(nze,nze)
      integer init
      integer L,i,j,x,y,z
      real start, finish

!c-------------------------------------- generate realization -----------------
      iz = 0
!c                     generate field
      
      
        write(*,*) 'in matmul'
      
      
      call cpu_time(start)
      
      ! get random numbers
      call vnorm( U, nxe*nye*nze )
      

      

      
    !process x direction
      
      do z = 1,nze
          do y = 1,nye
            !call matmul2(soil2(:,y,z),nye,nze,Rx(x,:),nxe,U2(:,y,z),nye,nze)
              soil(:,y,z) = matmul(R0x,U(:,y,z))
          end do
      end do
      
      !process y direction
      
      do z = 1,nze
          do x = 1,nxe
            !call matmul2(soil2(:,y,:),nxe,nze,Ry(y,:),nye,U2(:,y,:),nxe,nze)
            U(x,:,z) = matmul(R0y,soil(x,:,z))
          end do
      end do
      
      !process z direction
      
      do y = 1,nye
          do x = 1,nxe
            !call matmul2(soil2(:,:,z),nxe,nye,Rz(z,:),nze,U2(:,:,z),nxe,nye)
            soil(x,y,:) = matmul(R0z,U(x,y,:))
          end do
      end do
      
     
      call cpu_time(finish)
      write(*,*) finish-start
      
    read(*,*)
    !  
    ! 
    !
    !open(505,file='piecewise_fast.dat',access='stream')
    !write(505) soil
    !close(505)
    !stop
      
    
    end subroutine

    
    
end module
