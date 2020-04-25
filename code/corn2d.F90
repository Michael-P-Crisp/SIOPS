!c  *********************************************************************
!c  *                                                                   *
!c  *                        subroutine corn2d                          *
!c  *                                                                   *
!c  *********************************************************************
!c  Mixed Precision Version 2.21
!c  Written by Gordon A. Fenton, TUNS, Aug. 26, 1992
!c  Latest Update: Jun 9, 1999
!c
!c  PURPOSE  creates the parameter matrices required for the corner cell
!c           subdivisions of LAS2G (and LAS3G special case).
!c
!c  Requires:
!c   1) from libGAFsim:	DSIFA, DSISL, DCHOL2, DAXPY, DSWAP, DDOT, IDAMAX
!c
!c    n    is the rank of [CC][CC^T]. It can be either n = 3 for the 2-D case or
!c         n = 7 for the 3-D case.
!c    mc   is the mapping from the global covariance matrix `R' into the local
!c         `corner' covariance matrix
!c
!c  REVISION HISTORY:
!c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!c-----------------------------------------------------------------------------
      subroutine corn2d(R,ir,B,ib,S,is,CC,ic,n,AC,mc,iout,tol)
      real CC(ic,*), AC(4,n,*)
      real*8 R(ir,*), B(ib,*), S(is,*)
      real*8 RC(4,4), DA(4), BB(7,7)
      integer mc(4,*), indx(4)

   1  format(a,a,a)
   2  format(a,e13.4)
!c							extract R
      do j = 1, 4
         do i = 1, j
            RC(i,j) = R(mc(i,1), mc(j,1))
         end do
      end do
!c							factorize R
      call dsifa(RC,4,4,indx,ierr)
      if( ierr .ne. 0 ) then
         write(iout,1)'Error: unable to factorize corner covariance matrix in CORN2D.'
         stop
      endif
!c							for each corner...
      do nc = 1, 4
!c							make a copy of S
         do j = 1, n
            do i = 1, 4
               DA(i) = S(mc(i,nc),j)
            end do
!c							and solve for A
            call dsisl( RC, 4, 4, indx, DA )
!c							store in real*4
            do i = 1, 4
               AC(i,j,nc) = DA(i)
            end do
!c							update B
            do i = 1, j
               BB(i,j) = B(i,j) - S(mc(1,nc),i)*DA(1) - S(mc(2,nc),i)*DA(2) - S(mc(3,nc),i)*DA(3) - S(mc(4,nc),i)*DA(4)
            end do
         end do
!c							Cholesky Decomposition
         call dchol2( BB, 7, n, rerr )
         if( rerr .gt. tol ) then
            write(iout,1)'Corn2d: Cholesky decomposition of corner covariance matrix BB'
            write(iout,2)'        has maximum relative error of ',rerr
         endif
!c							store in real*4
         ii = 0
         do j = 1, n
             do i = 1, j
                ii = ii + 1
                CC(ii,nc) = BB(i,j)
             end do
         end do
      end do

      return
      end
