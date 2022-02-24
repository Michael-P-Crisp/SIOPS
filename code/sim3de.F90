! This documentation is deprecated, but kept for historical reasons.

!c  ********************************************************************
!c  *                                                                  *
!c  *                        subroutine sim3de                         *
!c  *                                                                  *
!c  ********************************************************************
!c  Single Precision Version 1.1
!c  Written by Gordon A. Fenton, Dalhousie University, Feb 18, 2003
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
!c		The leading dimensions of efld are assumed to be nxe and nye.
!c		(output)
!c
!c    nxe	number of random field cells in the X direction. (input)
!c
!c    nye	number of random field cells in the Y direction. (input)
!c
!c    nze	number of random field cells in the Z direction (vertical).
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
!c     nxew	number of finite elements in the x-direction (this can be
!c		different than nxe due to simulation requirements). (input)
!c
!c     nyew	number of finite elements in the y-direction (this can be
!c		different than nye due to simulation requirements). (input)
!c
!c     nzew	number of finite elements in the z-direction (this can be
!c		different than nze due to simulation requirements). (input)
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
!c---------------------------------------------------------------------------



module sim3de 
      
      contains



      subroutine transform_soil(efld,efldsize,sdata,stype)
            !transform the original generated soil from zero-mean, unit-variance, normal distribution
            !to a field of the desired statistics and distribution

            implicit none

            integer, intent(in) :: efldsize
            real, intent(inout) :: efld(efldsize)
            real, intent(in) :: sdata(4)
            character, intent(in) :: stype

            if(stype == 'n') then
            
                  efld = sdata(1) + efld * sdata(2)


                  
            else if(stype == 'l') then
            
                        efld = exp(sdata(3) + efld * sdata(4))
                        
            else if(stype == 'b') then
            
                        efld = sdata(1) + 0.5*(sdata(2)-sdata(1))*(1.0 + tanh((sdata(3)+sdata(4)*efld)/6.2831853))
                        
            else 
                        write(*,*) 'Incorrect distribution selected. n = normal, l = lognormal, b = bounded.'
                        
            end if

      end subroutine

      
      subroutine sim3d(efldgd,nxe,nye,nze,zroom,nxew,nyew,nzew,dx,dy,dz,kseed,MXM,MXK, &
                      C0,CC,CE,CS,CI,AC,AE,AS,AI,ATC, ATS, ATI,CTC,CTS,CTI,M,k1,k2,k3,kk,sdata,stype,anisotropy)
      
      implicit none
      
      integer MXK,MXM,NGS
      integer, intent(in) :: nxe,nye,nze !original dimensions of field
      integer, intent(in) :: nxew,nyew,nzew !dimensions of desired field, subset of full field
      integer, intent(in) :: zroom !some extra depth to store previous stage
      real(8), intent(in) :: dx,dy,dz
      real*4 :: efld(nxe,nye,zroom) !ceiling(real(9*nrfz)/8.0)
      real(4), intent(out) :: efldgd(nxew,nyew,nzew)
      integer, intent(in) :: kseed
      real*8 dvar, dpb !common variables
      real*8 dmin1, dble
      real*8 zerod, oned, dexp
      integer xpos,ypos,zpos !values for random subset
      real*4 XL, YL, ZL
      real(4), intent(in) :: sdata(4) !mean and sd for (log)normal distributions, lower,upper,m,s parameters for bounded
      character, intent(in) :: stype !type of distribution n = normal, l = lognormal, b = bounded
      integer, intent(in) :: anisotropy
      !save pmn, psd, liid, XL, YL, ZL, nn, xd, yd, zd
!c				export parameters to variance function
      !common/dparam/ dvar, dpb, dthx, dthy, dthz
      real, intent(in) :: CC(28,8,MXM), CE(28,12,MXM), CS(28,6,MXM), CI(28,MXM) 
      real, intent(in) :: C0(MXK*(MXK + 1)/2)
      real, intent(in) :: AC(8,7,8,MXM), AE(12,7,12,MXM), AS(18,7,6,MXM), AI(27,7,MXM)
      real, intent(in) :: ATC(4,7,4), ATS(6,7,4), ATI(9,7)
      real, intent(in) :: CTC(28,4), CTS(28,4), CTI(28)
      
      real(8) mean,sd !mean, standard deviation
      real randu !random generator function
      integer n,j !loop counters
      integer istat
      
      
      
      integer, intent(in) :: M,k1,k2,k3,kk
      real half,one,two
      data zerod/0.d0/, half/0.5/, one/1.0/, two/2.0/, oned/1.d0/
      
      NGS=max(ceiling(real(7*max(nxe,nye,nze)/2)),MXK)
      
!initialize internal generator.
      !kseed = iseed( kseed )
      
      
      XL     = float(nxe)*sngl(dx)	! size of random field
      YL     = float(nye)*sngl(dy)
      ZL     = float(nze)*sngl(dz)
      dvar=oned
      
		 !if(min(dthx,dthy,dthz) < min(dx,dy,dz) then
		 
		 !	  vnorm(efld,nxe*nrfy*nrfz)

		 !else
		 

!----------- create zero-mean, unit-variance normally-distributed field --------

			   call las3g(efld,nxe,nye,nze,XL,YL,ZL,kseed,istat &
			 ,C0, CC, CE, CS, CI, AC, AE, AS, AI, &
						ATC, ATS, ATI, CTC, CTS, CTI,M, k1, k2, k3, kk, MXM,MXK,NGS)

!----------- create random subset of desired size ------------
!c  PURPOSE Take a subset of the full array, determined by
!   the difference between the original dimensions and working dimensions.
!   The location of this subset is randomly generated to avoid the banding 
!   problem that can be seen when averaging the fields from a large number 
!   of realisations, as encountered in Goldsworthy's thesis (2006).
!   This subset also works around the limitation of the LAS dimensions being
!   aN**b.
      
						
		!generate random offsets
		!kseed=randu(kseed)*1234567890
		xpos=NINT(randu(0)*(nxe-nxew))
		if(xpos==0) xpos=1
		!kseed=randu(kseed)*1234567890
		ypos=NINT(randu(0)*(nye-nyew))
		if(ypos==0) ypos=1
		zpos=NINT(randu(0)*(nze - nzew))
		if(zpos==0) zpos=1
		
		
		efldgd = efld(xpos:xpos+nxew-1,ypos:ypos+nyew-1,zpos:zpos+nzew-1)
		
		
        
						
!----------- scale subset to be exactly zero mean and unit variance -------
            mean = sum(efldgd)/size(efldgd)
            sd = sqrt(sum((efldgd-mean)**2)/(size(efldgd)-1))
            
            efldgd = (efldgd-mean)/sd

           call transform_soil(efldgd,size(efldgd),sdata,stype)

           
                !mean = sum(efldgd)/size(efld)
                !sd = sqrt(sum((efldgd)**2)/(n-1))


      return
    end subroutine
      
    
    
!Basically the same subroutine as above, except instead of outputting a random subset, it exports the full soil.
!This is useful when generating a very large soil from which to take random subsets later, as is the case with the
!"superset" option set to true in the soil input file.
subroutine sim3d_nosubset(efld,nxe,nye,nze,zroom,nxew,nyew,nzew,dx,dy,dz,kseed,MXM,MXK, &
                      C0,CC,CE,CS,CI,AC,AE,AS,AI,ATC, ATS, ATI,CTC,CTS,CTI,M,k1,k2,k3,kk,sdata,stype,anisotropy)
      
      implicit none
      
      integer MXK,MXM,NGS
      integer, intent(in) :: nxe,nye,nze !original dimensions of field
      integer, intent(in) :: nxew,nyew,nzew !dimensions of desired field, subset of full field
      integer, intent(in) :: zroom !some extra depth to store previous stage
      real(8), intent(in) :: dx,dy,dz
      real(4), intent(out) :: efld(nxe,nye,zroom) !ceiling(real(9*nrfz)/8.0)

      integer, intent(in) :: kseed
      real*8 dvar, dpb !common variables
      real*8 dmin1, dble
      real*8 zerod, oned, dexp
      integer xpos,ypos,zpos !values for random subset
      real*4 XL, YL, ZL
      real(4), intent(in) :: sdata(4) !mean and sd for (log)normal distributions, lower,upper,m,s parameters for bounded
      character, intent(in) :: stype !type of distribution n = normal, l = lognormal, b = bounded
      integer, intent(in) :: anisotropy
      !save pmn, psd, liid, XL, YL, ZL, nn, xd, yd, zd
!c				export parameters to variance function
      !common/dparam/ dvar, dpb, dthx, dthy, dthz
      real, intent(in) :: CC(28,8,MXM), CE(28,12,MXM), CS(28,6,MXM), CI(28,MXM) 
      real, intent(in) :: C0(MXK*(MXK + 1)/2)
      real, intent(in) :: AC(8,7,8,MXM), AE(12,7,12,MXM), AS(18,7,6,MXM), AI(27,7,MXM)
      real, intent(in) :: ATC(4,7,4), ATS(6,7,4), ATI(9,7)
      real, intent(in) :: CTC(28,4), CTS(28,4), CTI(28)
      
      real(8) mean,sd
      real mean2,sd2 !mean, standard deviation
      real randu !random generator function
      integer n,j,x,y,z !loop counters
      integer istat
      
      
      
      integer, intent(in) :: M,k1,k2,k3,kk
      real half,one,two
      data zerod/0.d0/, half/0.5/, one/1.0/, two/2.0/, oned/1.d0/
      
      NGS=max(ceiling(real(7*max(nxe,nye,nze)/2)),MXK)
      
!initialize internal generator.
      !kseed = iseed( kseed )
      
      
      XL     = float(nxe)*sngl(dx)	! size of random field
      YL     = float(nye)*sngl(dy)
      ZL     = float(nze)*sngl(dz)
      dvar=oned
      
		 !if(min(dthx,dthy,dthz) < min(dx,dy,dz) then
		 
		 !	  vnorm(efld,nxe*nrfy*nrfz)

		 !else
		 

!----------- create zero-mean, unit-variance normally-distributed field --------

			   call las3g(efld,nxe,nye,nze,XL,YL,ZL,kseed,istat &
			 ,C0, CC, CE, CS, CI, AC, AE, AS, AI, &
						ATC, ATS, ATI, CTC, CTS, CTI,M, k1, k2, k3, kk, MXM,MXK,NGS)

!----------- create random subset of desired size ------------
!c  PURPOSE Take a subset of the full array, determined by
!   the difference between the original dimensions and working dimensions.
!   The location of this subset is randomly generated to avoid the banding 
!   problem that can be seen when averaging the fields from a large number 
!   of realisations, as encountered in Goldsworthy's thesis (2006).
!   This subset also works around the limitation of the LAS dimensions being
!   aN**b.
               
               
               
      
            mean = sum(efld(:,:,:nze))/size(efld(:,:,:nze))
            sd = sqrt(sum((efld(:,:,:nze)-mean)**2)/(size(efld(:,:,:nze))-1))
 
            efld(:,:,:nze) = (efld(:,:,:nze)-mean)/sd
            

     

         
            
            
                !open(505,file='soil.dat',access='stream')
				!write(505) efld(:,:,:nze)
				!close(505)
            
      



      return
    end subroutine
      
    
end module
    
    
    
    

