!c  ***********************************************************************
!c  *                                                                     *
!c  *                          subroutine dchol2                          *
!c  *                                                                     *
!c  ***********************************************************************
!c  Double Precision Version 3.21
!c  Written by Gordon A. Fenton, TUNS, Feb. 20, 1994
!c  Latest Update: Jun 9, 1999
!c
!c  PURPOSE  to compute the LL' factorization of a symmetric covariance matrix
!c           [A]. [A] need only be non-negative definite.
!c
!c  This function reduces a matrix A into its LL' decomposition (where
!c  prime indicates transpose) in place. No scaling or pivoting is performed
!c  during the decomposition. This algorithm is specialized for the
!c  decomposition of a covariance matrix which may or may not be positive
!c  definite (but is at least non-negative definite). If a diagonal element
!c  is found to be (algorithmically) zero, then the entire associated column
!c  of L is set to zero. This allows the decomposition of covariance matrices
!c  such as
!c                 [A]                 [L]            [L']
!c           __          __      __          __ __          __
!c           | 1  1  1  1 |      | 1  0  0  0 | | 1  1  1  1 |
!c           | 1  1  1  1 |      | 1  0  0  0 | | 0  0  0  0 |
!c           | 1  1  1  1 |  =   | 1  0  0  0 | | 0  0  0  0 |
!c           | 1  1  1  1 |      | 1  0  0  0 | | 0  0  0  0 |
!c           `-          -'      `-          -' `-          -'
!c
!c
!c  An estimate of the maximum relative error is computed herein by comparing
!c  the lowermost diagonal element of A and L*L'.
!c
!c  Arguments to the routine are as follows;
!c
!c     A    real array of size at least n x n containing on input the
!c          matrix of coefficients in its upper triangle at least. On
!c          ouput, A will contain L' written in the upper triangle.
!c          (input/output)
!c
!c    ia    leading dimension of A exactly as specified in the calling routine.
!c          (input)
!c
!c     n    integer giving the size of the matrix A. (input)
!c
!c  rerr    estimate of the maximum relative error between L*L' and A. This
!c          is the relative error between the lower-rightmost diagonal element
!c          of A and L*L'. (output)
!c
!c  NOTE: the value of `tol' is chosen rather empirically using a fairly
!c        pathological 256 x 256 matrix. It was found that other values
!c        (such as tol = 0) could lead to large errors in [L]*[L^T] when
!c        the matrix is large and nearly singular.
!c
!c  REVISION HISTORY:
!c  3.21	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
!c-----------------------------------------------------------------------------
      subroutine dchol2( A, ia, n, rerr )
      implicit real*8 (a-h,o-z)
      dimension A(ia,*)
      real*4 rerr
      data zero/0.d0/
      sqrt(z) = dsqrt(z)
      abs(z)  = dabs(z)
!c					preliminaries
      A1  = A(n,n)
      if( n .gt. 3 ) then
         A2  = A(n-1,n-1)
         A3  = A(n-2,n-2)
      endif
      rerr = 0.
      if( n .gt. 128 ) then
         tol = 1.d-5
      else
         tol = 1.d-20
      endif

      do i = 1, n
!c					find diagonal elements
         s = zero
         do j = 1, i-1
            s = s + A(j,i)*A(j,i)
         end do
         t = A(i,i) - s

         if( t .le. tol ) then
!c					diagonal element is algorithmically 0
            A(i,i) = zero
            do j = i+1, n
               A(i,j) = zero
            end do
         else
!c					else find off-diagonal elements
            A(i,i) = sqrt(t)
            do j = i+1, n
               s = zero
               do k = 1, i-1
                  s = s + A(k,i)*A(k,j)
               end do
               A(i,j) = (A(i,j)/A(i,i)) - (s/A(i,i))
            end do
         endif
      end do
!c					estimate error
      if( n .gt. 3 ) then
         t1 = A(1,n)*A(1,n)+A(n-1,n)*A(n-1,n) + A(n,n)*A(n,n)
         t2 = A(1,n-1)*A(1,n-1) + A(n-1,n-1)*A(n-1,n-1)
         t3 = A(1,n-2)*A(1,n-2)
         do i = 2, n-2
            t1 = t1 + A(i,n)*A(i,n)
            t2 = t2 + A(i,n-1)*A(i,n-1)
            t3 = t3 + A(i,n-2)*A(i,n-2)
         end do
         r1 = abs((A1 - t1)/A1)
         r2 = abs((A2 - t2)/A2)
         if( r1 .gt. r2 ) r2 = r1
         r3 = abs((A3 - t3)/A3)
         if( r2 .gt. r3 ) r3 = r2
         rerr = r3
      else
         t = A(1,n)*A(1,n)
         do i = 2, n
            t = t + A(i,n)*A(i,n)
         end do
         rerr = abs((A1 - t)/A1)
      endif

      return
      end
