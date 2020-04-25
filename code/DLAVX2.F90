!c  *******************************************************************
!c  *                                                                 *
!c  *                     Function dlavx2                             *
!c  *                                                                 *
!c  *******************************************************************
!c  Double Precision Version 1.52
!c  Written by Gordon A. Fenton, TUNS, 1991
!c  Latest Update: Jun 7, 2001
!c
!c  PURPOSE  returns the covariance between two points in a 2-D Markovian
!c           random field.
!c
!c  This function returns the covariance between two points in a 2-D
!c  Markov random field, separated by lag vector {X,Y}.
!c  The covariance function for this process is given by
!c
!c         B(x,y) = var * exp{ -sqrt[(2*x/dthx)^2 + (2*y/dthy)^2] }
!c
!c  where var is the point variance, and dthx and dthy are the scales of
!c  fluctuation in the x- and y-directions, respectively. The parameters
!c  var, dthx, and dthy, are brought in via the common block /dparam/.
!c
!c  If var < 0, then this function returns the variance of a local average
!c  of the process, |var|*V(X,Y), averaged over the domain X x Y. This variance
!c  is obtained by 5-pt Gauss quadrature integration of the covariance function.
!c  The 4-fold integral is reduced to a 2-fold integration by taking advantage
!c  of the quadrant symmetry of the covariance function.
!c
!c  See "Random Fields" by E. Vanmarcke and G.A. Fenton's thesis
!c  for more details.
!c
!c  The arguments to this routine are just the X and Y components of the
!c  lag vector between the points (or the dimensions of the averaging area,
!c  if var < 0).
!c
!c  NOTE: when the element non-square, extremely small scales of fluctuation
!c	 tend to yield inaccurate local average variances (either directly
!c	 here using var < 0, or through dcvaa2). More Gauss points are
!c	 needed in the integration. See F77/sim/utils for a version of
!c	 dlavx2 with 20-pt Gaussian quadrature.
!c
!c  REVISION HISTORY:
!c  1.1	brought the pi/2 factor into the function directly, rather than
!c	making it a parameter (was `dum1').
!c  1.2	added lvarfn flag to end of /dparam/ common block. Set it to
!c	true if the variance function is to be returned. Otherwise this
!c	function returns the covariance. (Apr 29/00)
!c  1.3	eliminated lvarfn -- now return covariances only if var < 0 (Mar 27/01)
!c  1.4	properly handled sign on var for covariances (Apr 5/01)
!c  1.5	reversed default - now return covariances if var > 0.
!c	Variance function values now computed by Gauss quadrature (Apr 11/01)
!c  1.51	revised above documentation to reflect revision 1.5 (May 9/01)
!c  1.52	added comment on accuracy above (Jun 7/01)
!c---------------------------------------------------------------------------

    
      real*8 function dlavx2( X, Y , var, dpb, H, G, da, db ,dthx,dthy)
      parameter (NG = 5)
      implicit real*8 (a-h,o-z)
      dimension w(NG), z(NG)
      !common/dparam/ var, dpb, dthx, dthy, dthz
      data zero/0.d0/, one/1.d0/, two/2.d0/, three/3.d0/
      data quart/0.25d0/, half/0.5d0/
      data w/ .236926885056189d0,.478628670499366d0,.568888888888889d0,.478628670499366d0,.236926885056189d0/
      data z/-.906179845938664d0,-.538469310105683d0,.000000000000000d0,.538469310105683d0, .906179845938664d0/
      exp(zz)  = dexp(zz)
      abs(zz)  = dabs(zz)
      sqrt(zz) = dsqrt(zz)

      aY = abs(Y)
      aX = abs(X)
      if( var .lt. zero ) then
         if( (dthx .eq. zero) .and. (dthy .eq. zero) ) then
            if( (X .eq. zero) .and. (Y .eq. zero) ) then
               dlavx2 = -var
            else
               dlavx2 = zero
            endif
         elseif( dthx .eq. zero ) then
            if( X .eq. zero ) then
               r2 = half*aY
               ty = two/dthy
               d1 = zero
               do j = 1, NG
                  yj = r2*(one + z(j))
                  a2 = ty*yj
                  d1 = d1 + w(j)*(one-z(j))*exp(-a2)
               end do
               dlavx2 = -half*var*d1
            else
               dlavx2 = zero
            endif
         elseif( dthy .eq. zero ) then
            if( Y .eq. zero ) then
               r1 = half*aX
               tx = two/dthx
               d1 = zero
               do i = 1, NG
                  xi = r1*(one + z(i))
                  a1 = tx*xi
                  d1 = d1 + w(i)*(one-z(i))*exp(-a1)
               end do
               dlavx2 = -half*var*d1
            else
               dlavx2 = zero
            endif
         else
            r1 = half*aX
            r2 = half*aY
            tx = two/dthx
            ty = two/dthy
            d1 = zero
            do i = 1, NG
               xi = r1*(one + z(i))
               a1 = tx*xi
               d2 = zero
               do j = 1, NG
                  yj = r2*(one + z(j))
                  a2 = ty*yj
                  T  = sqrt(a1*a1 + a2*a2)
                  d2 = d2 + w(j)*(one-z(j))*exp(-T)
               end do
               d1 = d1 + w(i)*(one-z(i))*d2
            end do
            dlavx2 = -quart*var*d1
         endif
      else !here
         if( (dthx .eq. zero) .and. (dthy .eq. zero) ) then
            if( (X .eq. zero) .and. (Y .eq. zero) ) then
               dlavx2 = var
            else
               dlavx2 = zero
            endif
         elseif( dthx .eq. zero ) then
            if( X .eq. zero ) then
               dlavx2 = var*exp(-two*aY/dthy)
            else
               dlavx2 = zero
            endif
         elseif( dthy .eq. zero ) then
            if( Y .eq. zero ) then
               dlavx2 = var*exp(-two*aX/dthx)
            else
               dlavx2 = zero
            endif
         else
            a1 = two*aX/dthx
            a2 = two*aY/dthy
            T  = sqrt(a1*a1 + a2*a2)
            dlavx2 = var*exp( -T )
         endif
      endif

      return
      end function
      
