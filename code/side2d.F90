!c  *********************************************************************
!c  *                                                                   *
!c  *                        subroutine side2d                          *
!c  *                                                                   *
!c  *********************************************************************
!c  Mixed Precision Version 2.31
!c  Written by Gordon A. Fenton, TUNS, Aug. 26, 1992
!c  Latest Update: Jun 9, 1999
!c
!c  PURPOSE  creates the parameter matrices required for the side cell
!c           subdivisions of LAS2G.
!c
!c  Requires:
!c   1) from libGAFsim:	DSIFA, DSISL, DCHOL2, DAXPY, DSWAP, DDOT, IDAMAX
!c
!c    n    is the rank of [CC][CC^T]. It can be either n = 3 for the 2-D case or
!c         n = 7 for the 3-D case.
!c    ms   is the mapping from the global covariance matrix `R' into the local
!c         `side' covariance matrix
!c
!c  REVISION HISTORY:
!c  2.31	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!c-----------------------------------------------------------------------------
      subroutine side2d(R,ir,B,ib,S,is,CS,ic,n,AS,ms,iout,tol)
      real CS(ic,*), AS(6,n,*)
      real*8 R(ir,*), B(ib,*), S(is,*)
      real*8 RS(6,6), DA(6), BB(7,7)
      integer ms(6,*), indx(6)

   1  format(a,a,a)
   2  format(a,e13.4)
!c							for each side...
      do ns = 1, 4
!c							extract R
         do j = 1, 6
            do i = 1, j
               RS(i,j) = R(ms(i,ns), ms(j,ns))
            end do
         end do
!c							factorize R
         call dsifa( RS, 6, 6, indx, ierr )
         if( ierr .ne. 0 ) then
            write(iout,1)'Error: unable to factorize side covariance matrix in SIDE2D.'
            stop
         endif
!c							make a copy of S
         do j = 1, n
            do i = 1, 6
               DA(i) = S(ms(i,ns),j)
            end do
!c							and solve for A
            call dsisl( RS, 6, 6, indx, DA )
!c							store in real*4
            do i = 1, 6
               AS(i,j,ns) = DA(i)
            end do
!c							update B
            do i = 1, j
               BB(i,j) = B(i,j) - S(ms(1,ns),i)*DA(1) - S(ms(2,ns),i)*DA(2) &
              - S(ms(3,ns),i)*DA(3) - S(ms(4,ns),i)*DA(4) - S(ms(5,ns),i)*DA(5) - S(ms(6,ns),i)*DA(6)
            end do
         end do
         
!c							Cholesky Decomposition
         call dchol2( BB, 7, n, rerr )
         if( rerr .gt. tol ) then
            write(iout,1)'Side2d: Cholesky decomposition of side covariance matrix BB'
            write(iout,2)'        has maximum relative error of ',rerr
         endif
!c							store in real*4
         ii = 0
         do j = 1, n
         do i = 1, j
            ii = ii + 1
            CS(ii,ns) = BB(i,j)
         end do
         end do
      end do

      return
      end
