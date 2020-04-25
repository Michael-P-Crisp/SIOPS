module LAS1D
    
    
    contains
    

    !1D random field setup and generation. See LAS2I and LAS2G for details.
    !This version is stripped down to only work with CMD, no actual LAS.
    !Although the stage 0 fields are typically big enough that this won't be an issue.
    
    
    !c  *********************************************************************
!!c  *                                                                   *
!c  *                         subroutine sim2sd                         *
!c  *                                                                   *
!c  *********************************************************************
      subroutine sim1sd_init(nxe,dx,thx,varfnc,C0)

      !real cfld(*), c(7)
      real(8) thx, thy
      integer, intent(in) :: nxe
      integer kseed
      character(6) varfnc
      logical debug, liid, shofld, lxfld
      real*8 dvar, dpb, dx
      
      real(8), external :: dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      real*8 dmin1, dble
      real, intent(out) :: C0((nxe*(nxe + 1))/2)
      real(8) H, G ,da,db !for some of the more obscure correlation models
      !save XL, YL, ienter, liid
      !external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      !common/dparam/ dvar, dpb, dthx, dthy
      
      data zero/0.0/, half/0.5/, one/1.0/, onept5/1.5/, twopt5/2.5/
      data ienter/0/, twopi/6.283185/
      data tol/1.e-3/

   1  format(a,a)

!c-------------------------------------- initialize ---------------------

         !liid   = ((thx .eq. zero) .and. (thy .eq. zero))

         XL     = float(nxe)*dx
         YL     = float(1)*dx

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
                   call las1i( dlavx2, nxe, XL,C0 , tol, thx )

				elseif( varfnc .eq. 'dlsep2' ) then
				   call las1i( dlsep2, nxe, XL,C0 , tol, thx )
                   
				elseif( varfnc .eq. 'dlspx2' ) then
				   call las1i( dlspx2, nxe, XL,C0 , tol, thx )
                   
				elseif( varfnc .eq. 'dlafr2' ) then
				   call las1i( dlafr2, nxe, XL,C0 , tol, thx )
                   
				elseif( varfnc .eq. 'dlsfr2' ) then
				   call las1i( dlsfr2, nxe, XL,C0 , tol, thx )
                   
				endif
            !end if
            

      
      return
      end
    

      
      !Generate 1D normally-distributed, zero-mean, unit-variance random field.
      !Must have proper processing done before calling
      
      subroutine las1i( dvarfn, N1, XL,C0 , tol, thx )
      real C0((N1*(N1 + 1))/2)
      real*8 R(9,9,2), B(4,4), S(9,4)
      real*8 R0(N1*N1)
      real*8 T1, T2, dvarfn, dble
      logical lformR
      real(8) var, pb, H, G !field variance, averaging length, hurst components in horizontal and vertical directions
      real(8) da, db 
      real(8) thx !SOF in x and y directions
      external dvarfn
      !character dvarfn
   1  format(a,a,a)
   2  format(a,e13.4)
   3  format(a,i4,a,i4,a,i4,a)

!c                    initialize internal generator
      !kseed = iseed( kseed ) !already initialized
!c                    form initial covariance matrix
      T1 = dble(XL)/dble(N1)
      T2 = T1
      var = 1.d0
      call dcvit2( dvarfn, R0, N1, R, 9, N1, 1, T1, T2 , var, pb, H, G, da, db , thx, thx)

!c                    and compute its cholesky decomp
      call dchol2( R0, N1, N1, rerr )
      if( rerr .gt. tol ) then
         write(*,1)'Warning: Cholesky decomposition of stage 0 covariance matrix'
         write(*,2)'         has maximum relative error of ',rerr
      endif
!c                    save in real*4 for LAS2G
      L = 0
      do j = 1, N1
		  do i = 1, j
			 L = L + 1
			 C0(L) = R0(i+(j-1)*N1)
		  end do
      end do

      return
      end subroutine
      
      
      
      
      

    subroutine las1g( Z, N1 ,C0)
      real Z(N1)
      integer init
      real C0((N1*(N1 + 1))/2), U(N1)
      integer L,i,j

!c-------------------------------------- generate realization -----------------
      iz = 0
!c                     generate stage 0 field
      call vnorm( U, N1 )
      L = 1
      do i = 1, N1
         Z(iz+i) = C0(L)*U(1)
         do j = 2, i
            Z(iz+i) = Z(iz+i) + C0(L+j-1)*U(j)
         end do
         L = L + i
      end do
    
      end subroutine
    
    
    
end module