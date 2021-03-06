!c  *********************************************************************
!c  *                                                                   *
!c  *                        subroutine intr3d                          *
!c  *                                                                   *
!c  *********************************************************************
!c  Mixed Precision Version 1.1
!c  Written by Gordon A. Fenton, TUNS, Aug. 2, 1993
!c  Latest Update: Feb. 22, 1994
!c
!c  PURPOSE  creates the parameter matrices required for the interior cell
!c           subdivision of LAS3G.
!c
!c  Notes:
!c   1) arrays B and R are destroyed in this routine.
!c
!c  Requires:
!c   1) from libGAFsim:	DSIFA, DSISL, DCHOL2, DAXPY, DSWAP, DDOT, IDAMAX
!c-----------------------------------------------------------------------------
      subroutine intr3d( R, ir, B, ib, S, is, CI, AI, iout, tol )
      real*8 R(ir,*), B(ib,*), S(is,*), DA(27)
      real CI(*), AI(27,*)
      integer indx(27)

   1  format(a,a,a)
   2  format(a,e13.4)
!c							factorize R
      call dsifa(R,ir,27,indx,ierr)
      if( ierr .ne. 0 ) then
         write(iout,1)'Error: unable to factorize interior covariance matrix in INTR3D.'
         stop
      endif
!c							make a copy of S
      do j = 1, 7
         do i = 1, 27
            DA(i) = S(i,j)
          end do
!c							and solve for A
         call dsisl(R,ir,27,indx,DA)
!c							store in real*4
         do i = 1, 27
            AI(i,j) = DA(i)
         end do
!c							update B
         do i = 1, j
            do k = 1, 27
               B(i,j) = B(i,j) - S(k,i)*DA(k)
            end do
         end do
      end do
      
!c							Cholesky Decomposition
      call dchol2( B, ib, 7, rerr )
      if( rerr .gt. tol ) then
         write(iout,1)'Warning: Cholesky decomposition of interior covariance matrix B'
         write(iout,2)'         has maximum relative error of ',rerr
      endif
!c							store in real*4
      ii = 0
      do j = 1, 7
      do i = 1, j
         ii = ii + 1
         CI(ii) = B(i,j)
      end do
      end do

      return
      end
