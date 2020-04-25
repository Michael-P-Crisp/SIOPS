!c  **********************************************************************
!c  *                                                                    *
!c  *                   Integer Function iseed                           *
!c  *                                                                    *
!c  **********************************************************************
!c  Single Precision Version 1.1
!c  Written by Gordon A. Fenton, Princeton, Dec. 8, 1988.
!c  Latest Update: Oct 14, 1996
!c
!c  PURPOSE  initializes the system pseudo-random number generated using
!c           process ID as default seed.
!c
!c  Initializes the local random number generator RANDF. If the argument integer
!c  seed (kseed) is zero, a random seed is calculated based on the process ID
!c  of the parent process. The function returns the actual seed used. This
!c  routine is system specific.
!c
!c  Notes:
!c   1) any function which is reasonably certain to return a different integer
!c      on each invocation of the calling process can be used to generate a
!c      pseudo-random seed here. Some possibilities might be wall clock time,
!c      Process ID number (UNIX based)...
!c
!c  Requires:
!c   1) from Fortran lib: GETPID
!c
!c  REVISION HISTORY:
!c  1.1	now using new randu function (RAN2 from Numerical Recipes, 2nd Ed)
!   1.2 Fixed bug where rfirst would always be zero due to integer nature Dec 2016
!c	(Oct 14/96)
!c--------------------------------------------------------------------------
      integer function iseed( kseed ) result(rfirst)
      integer kseed, getpid, jseed
!c
      jseed = kseed
      if( kseed .eq. 0 ) then
!c                                avoid seeds of 0
	  jseed = getpid() + 2
      endif
!c                                initialize the generator
      rfirst = randu( jseed )*1234567890 !I added in this multiplication, otherwise kseed would always be zero as it's an integer and 0<kseed<1

      end
