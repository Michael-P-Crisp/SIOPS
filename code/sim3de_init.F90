!c  ********************************************************************
!c  *                                                                  *
!c  *                        subroutine sim3de                         *
!c  *                                                                  *
!c  ********************************************************************
!c  Single Precision Version 1.1
!c  Written by Gordon A. Fenton, Dalhousie University, Feb 18, 2003
!c  Modified by Michael Crisp August 17, 2016
!c  Latest Update: August 17, 2016
!c
!c  PURPOSE  to generate a realization of a 3-D elastic modulus field
!c
!c  This routine creates a realization of an elastic modulus (E) field in 3-D.
!c  The generation proceeds as follows;
!c
!c   1) generate a zero-mean, unit-variance, Gaussian random field in 3-D
!c      with prescribed covariance structure.
!c
!c   2) transform the field into a lognormally distributed conductivity field
!c      on a point-wise basis by setting Y = exp(emn + esd*G) at each field
!c      coordinate (emn and esd being the point mean and standard deviation
!c      of the log-E field and G is the standard Gaussian field locally
!c      averaged over each element domain).
!c
!c   3) map discrete E values into the finite element mesh, so that
!c      each element takes on a particular E value
!c
!c  Arguments are described as follows;
!c
!c    efld	real array of size at least nxe x nye x nze which on output
!c		will contain the desired realization of E.
!c		The leading dimensions of efld are assumed to be nrfx and nrfy.
!c		(output)
!c
!c    nrfx	number of random field cells in the X direction. (input)
!c
!c    nrfy	number of random field cells in the Y direction. (input)
!c
!c    nrfz	number of random field cells in the Z direction (vertical).
!c		(input)
!c
!c      dx	real*8 value containing the physical size of an element in the
!c		X direction. (input)
!c
!c      dy	real*8 value containing the physical size of an element in the
!c		Y direction. (input)
!c
!c      dz	real*8 value containing the physical size of an element in the
!c		Z direction (vertical). (input)
!c
!c     thx	real*8 value containing the scale of fluctuation of log-E in
!c		the X direction.
!c		Note that if any of thx, thy, or thz are zero, the entire
!c		field is assumed to be made up of independent values (ie a
!c		white noise process). In this case, the mean and variance
!c		of each cell is determined directly from emn and esd, that
!c		is the affect of local averaging is not considered. (input)
!c
!c     thy	real*8 value containing the scale of fluctuation of log-E in
!c		the Y direction.
!c
!c     thz	real*8 value containing the scale of fluctuation of log-E in
!c		the Z direction.
!c
!c     emn	(real*8) target mean of E at a point. (input)
!c
!c     esd	(real*8) target standard deviation of E at a point. (input)
!c
!c   lunif	logical flag which is true if the E field is
!c		uniform, corresponding to an infinite scale of fluctuation
!c		in both directions. In this case, the E field
!c		is set equal to a single random value having mean emn and
!c		standard deviation esd. Local averaging is not performed.
!c		(input)
!c
!c   kseed	integer seed to be used for the pseudo-random number generator.
!c		If KSEED = 0, then a random seed will be used (based on the
!c		clock time when LAS2D is called for the first time).
!c		On output, KSEED is set to the value of the actual seed used.
!c		(input/output)
!c
!c  shofld	logical flag which is true if the generated random field
!c		is to be dumped to a DISPLAY format file. Only one field is
!c		dumped, so shofld is set to false after dumping. (input/output)
!c
!c   ieplt	integer vector of length at least 3 containing the element
!c		indices corresponding to the plane over which the random field
!c		plot is to be drawn. The elements of ieplt correspond
!c		to the (x,y,z) element indices and zero values in this vector
!c		correspond to all elements in that direction. Only one of the
!c		components is expected to be non-zero. So, for example, the
!c		vector (0,26,0) means that the x-z plane at y-element = 26 is
!c		to be drawn. If a 50 x 50 x 30 element domain has two
!c		footings, one centered at x-node = 26 and y-node = 15 and
!c		the other centered at x-node = 26 and y-node = 35, then
!c		normally the y-z plane through x-node = 26 would be drawn
!c		to see the random field under the two footings. Thus, in this
!c		case, idplt might be (25,0,0) or (26,0,0) (depending on which
!c		side of the nodal plane you want to be on). The z-indices are
!c		measured from the top surface of the soil mass, thus z = 1 is
!c		the top layer of the soil mass and z = nze is the bottom
!c		layer, etc. (input)
!c
!c
!c   iefld	unit number to which the first field of log-E is
!c		written. (input)
!c
!c     job	character string containing the title of the run. (input)
!c
!c    sub1	character string containing the subtitle of the run. (input)
!c
!c    sub2	character string containing the sub-subtitle of the run.
!c		(input)
!c
!c  varfnc	character string containing the name of the covariance
!c		function controlling the random fields. Possible covariance
!c		functions are (input)
!c		`dlavx3' - 3-D exponentially decaying (Markov) model
!c			   requires X-, Y-, and Z-direction scales of
!c			   fluctuation
!c		`dlafs3' - 3-D partially separable fractional Gaussian noise
!c			   model. Requires (H,G,delta) as parameters. In this
!c			   case, thx is H, thy is G, and delta is the minimum
!c			   element dimension. G is assumed to apply in the
!c			   Y-Z plane.
!c		`dlsep3' - 3-D separable (1D x 1D x 1D) Markov model
!c			   requires X-, Y-, and Z-direction scales of
!c			   fluctuation
!c		`dlsfr3' - 3-D separable fractional Gaussian noise model
!c			   requires (H_x,H_y,H_z,delta) as parameters. In this
!c			   case, thx is H_x, thy is H_y, thz is H_z, and
!c			   delta is the minimum element dimension.
!c		`dlspx3' - 3-D separable Gaussian decaying model
!c			   requires X-, Y-, and Z-direction scales of
!c			   fluctuation.
!c
!c   debug	logical flag which is true if debugging information is to be
!c		output to unit istat. (input)
!c
!c   istat	unit number to which error and warning messages are issued.
!c		(input)
!c
!c  dcheck	logical flag which is true if this is a data-check run only, in
!c		which case only the required field size and target statistics
!c		are actually computed herein. (input)
!c
!c     nxe	number of finite elements in the x-direction (this can be
!c		different than nrfx due to simulation requirements). (input)
!c
!c     nye	number of finite elements in the y-direction (this can be
!c		different than nrfy due to simulation requirements). (input)
!c
!c     nze	number of finite elements in the z-direction (this can be
!c		different than nrfz due to simulation requirements). (input)
!c
!c    eavg	arithmetic average of the E field for this realization. This
!c		is averaged over the finite element mesh, rather than the
!c		entire random field (in the event that nrfx does not equal
!c		nxe, etc). (output)
!c
!c    egeo	geometric average of the E field for this realization over
!c		the finite element mesh. (output)
!c
!c    ehrm	harmonic average of the E field for this realization over
!c		the finite element mesh. (output)
!c
!c  REVISION HISTORY:
!c  1.1	corrrected real*8 declarations of dlavx2, .. to dlavx3, .. (Jan 10/05)
!c  1.2    Add in option for a 2D normally distributed output for boundary between layers
!c---------------------------------------------------------------------------
      
      
      subroutine sim3d_init(nrfx,nrfy,nrfz,dx,dy,dz,dthx,dthy,dthz,varfnc,MXM,MXK, &
                      C0,CC,CE,CS,CI,AC,AE,AS,AI,ATC, ATS, ATI,CTC,CTS,CTI,M,k1,k2,k3,kk)
      integer MXK,MXM
      real*8, intent(in) :: dx,dy,dz,dthx,dthy,dthz
      integer, intent(in) :: nrfx,nrfy,nrfz
      integer kseed
      real*8 dvar, dpb !common variables
      real*8 dlavx3, dlsep3, dlspx3, dlafs3, dlsfr3
      real*8 dmin1, dble
      real*8 zerod, oned, dexp
      real*4 XL, YL, ZL
      real(8) H,G,F
      character(6), intent(in) :: varfnc
      external dlavx3, dlsep3, dlspx3, dlafs3, dlsfr3
      !save pmn, psd, liid, XL, YL, ZL, nn, xd, yd, zd
!c				export parameters to variance function
      !common/dparam/ dvar, dpb, dthx, dthy, dthz
      real, intent(out) :: CC(28,8,MXM), CE(28,12,MXM), CS(28,6,MXM), CI(28,MXM) 
      real, intent(out) :: C0(MXK*(MXK + 1)/2)
      real, intent(out) :: AC(8,7,8,MXM), AE(12,7,12,MXM), AS(18,7,6,MXM), AI(27,7,MXM)
      real, intent(out) :: ATC(4,7,4), ATS(6,7,4), ATI(9,7)
      real, intent(out) :: CTC(28,4), CTS(28,4), CTI(28)
      
      
      
      integer, intent(out) :: M,k1,k2,k3,kk
      data zerod/0.d0/, half/0.5/, one/1.0/, two/2.0/, oned/1.d0/
      data tol/1.e-3/
      
      
      XL     = float(nrfx)*sngl(dx)	! size of random field
      YL     = float(nrfy)*sngl(dy)
      ZL     = float(nrfz)*sngl(dz)
      dvar=oned
      
!c						initialize internal generator
      kseed = iseed( kseed )


!c-------------------------------------- initialize -------------------------

		!if(min(dthx,dthy) < min(dx,dy) then
		!		write(*,*) 'warning, SOF is smaller than element. Will generate random noise.'
		!		return
				
		!else

!c						get correlation matrices
            if( varfnc .eq. 'dlavx3' ) then
            
			   call las3i( dlavx3, nrfx,nrfy,nrfz, XL, YL, ZL, MXM, MXK, C0, CC, CE, CS, CI, AC, AE, AS, AI, &
                    ATC, ATS, ATI, CTC, CTS, CTI,M, k1, k2, k3, kk, iout, tol,dvar, dpb, dthx, dthy, dthz , H, G, F)  
                             
            elseif( varfnc .eq. 'dlsep3' ) then
			   call las3i( dlsep3, nrfx,nrfy,nrfz, XL, YL, ZL, MXM, MXK, C0, CC, CE, CS, CI, AC, AE, AS, AI, &
                    ATC, ATS, ATI, CTC, CTS, CTI,M, k1, k2, k3, kk, iout, tol,dvar, dpb, dthx, dthy, dthz , H, G, F)
           
            elseif( varfnc .eq. 'dlspx3' ) then
				call las3i( dlspx3, nrfx,nrfy,nrfz, XL, YL, ZL, MXM, MXK, C0, CC, CE, CS, CI, AC, AE, AS, AI, &
                    ATC, ATS, ATI, CTC, CTS, CTI,M, k1, k2, k3, kk, iout, tol,dvar, dpb, dthx, dthy, dthz , H, G, F)

            elseif( varfnc .eq. 'dlafs3' ) then
                call las3i( dlafs3, nrfx,nrfy,nrfz, XL, YL, ZL, MXM, MXK, C0, CC, CE, CS, CI, AC, AE, AS, AI, &
                    ATC, ATS, ATI, CTC, CTS, CTI,M, k1, k2, k3, kk, iout, tol,dvar, dpb, dthx, dthy, dthz , H, G, F)

            elseif( varfnc .eq. 'dlsfr3' ) then
                call las3i( dlsfr3, nrfx,nrfy,nrfz, XL, YL, ZL, MXM, MXK, C0, CC, CE, CS, CI, AC, AE, AS, AI, &
                    ATC, ATS, ATI, CTC, CTS, CTI,M, k1, k2, k3, kk, iout, tol,dvar, dpb, dthx, dthy, dthz , H, G, F)

                   
            else
				write(*,*) 'incorrect correlation function chosen, try "dlavx3"'
			end if
	    !end if


     


      return
      end subroutine
      
