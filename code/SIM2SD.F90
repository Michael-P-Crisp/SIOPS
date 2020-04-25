!c  *********************************************************************
!c  *                                                                   *
!c  *                         subroutine sim2sd                         *
!c  *                                                                   *
!c  *********************************************************************
    module sim2sd
    
    contains
    
      subroutine sim2d(cfld,nxe,nye,nxew,nyew,zroom,dx,dy,kseed,MXM,MXK,C0,CT,CC,CS,CI,AT,AC,AS,AI, M, k1, k2, kk,sdata,stype)

      integer nxe, nye, kseed, zroom
      real :: efld(nxe,zroom) !ceiling(real(5*nye)/4.0)
      real, intent(out) :: cfld(:,:) 
      character(6) varfnc
      integer, intent(in) :: MXM,MXK
      integer NGS
      real*8 dx, dy
      real*8 dmin1, dble
      integer xpos, ypos
      real, intent(in) :: C0((MXK*(MXK + 1))/2)
      real, intent(in) :: CT(6,2), CC(6,4,MXM), CS(6,4,MXM), CI(6,MXM)
      real, intent(in) :: AT(3,3,2), AC(4,3,4,MXM), AS(6,3,4,MXM), AI(9,3,MXM)
      integer :: M,K1,K2,KK
      real(4), intent(in) :: sdata(4) !mean and sd for (log)normal distributions, lower,upper,m,s parameters for bounded
      character, intent(in) :: stype !type of distribution n = normal, l = lognormal, b = bounded
      !save XL, YL, ienter, liid
      !external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      !common/dparam/ dvar, dpb, dthx, dthy, dthz
      data zero/0.0/, half/0.5/, one/1.0/, onept5/1.5/, twopt5/2.5/
      data ienter/0/, twopi/6.283185/
      
      real(8) mean,sd

   1  format(a,a)
   

	  NGS=max(3*(9*2**(MXM-1)),MXK) !This value could potentially be invalid if more than 2147483647 are required
	  
	 !initialize internal generator.
      kseed = iseed( kseed )

      XL     = float(nxe)*dx
      YL     = float(nye)*dy
      
      dvar = 1.d0
		if( varfnc .eq. 'dlafr2' .or. varfnc .eq. 'dlsfr2' ) then
		   dpb = dx
		   if( dy .lt. dx ) dpb = dy
		endif
      
      !c------------------------------------ generate the random field(s) -----------

	 !if(min(thx,thy) < min(dx,dy) then
		 
	 !	vnorm(efld,nxe*nye)

	 !else
	 
	 !if( sdata(2) < dy/100 .and. (stype == 'n' .or. stype == 'l')) then !
	 
	 !       cfld = 0.0 !if the standard deviation is 100x smaller than the element size, just set it to zero
	        
	 !else

	  !----------- create zero-mean, unit-variance normally-distributed field --------

			call las2g(efld,nxe,nye,XL,YL,kseed,istat,C0,CT,CC,CS,CI,AT,AC,AS,AI,M,K1,K2,KK,NN, MXM,MXK,NGS)
			
	  !----------- create random subset of desired size ------------		
	  	!generate random offsets
		jseed=randu(jseed)*1234567890
		xpos=NINT(randu(jseed)*(nxe-nxew))
		if(xpos==0) xpos=1
		jseed=randu(jseed)*1234567890
		ypos=NINT(randu(jseed)*(nye-nyew))
		if(ypos==0) ypos=1

		!subset layer roughness
		cfld = efld(xpos:xpos+nxew-1,ypos:ypos+nyew-1) 
	
	  !----------- scale subset to be exactly zero mean and unit variance -------
            
            !mean = sum(cfld)/size(cfld)
            !sd = sqrt(sum((cfld-mean)**2)/(size(cfld)-1))
            
            
            
            !cfld = (cfld-mean) !/sd

	  !---------------- apply neccessary transformations ----------------
	              
           if(stype == 'n') then
           
				cfld = sdata(1) + cfld * sdata(2)
				
		   else if(stype == 'l') then
		   
				cfld = exp(sdata(3) + cfld * sdata(4))
				
		   else if(stype == 'b') then
		   
				cfld = sdata(1) + 0.5*(sdata(2)-sdata(1))*(1.0 + tanh((sdata(3)+sdata(4)*cfld)/6.2831853))
				
		   else 
				write(*,*) 'Incorrect distribution selected. n = normal, l = lognormal, b = bounded.'
				
           end if
           
           

           
           
		   
	!end if

      
      return
      end

      end module