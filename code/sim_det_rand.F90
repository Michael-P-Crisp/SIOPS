!c  *********************************************************************
!c  *                                                                   *
!c  *                         subroutine sim_det_rand                   *
!c  *                                                                   *
!c  *********************************************************************
      subroutine sim_det_rand(cfld,nxe,nye,nze,kseed,sdata,stype,whitenoise)

      real, intent(out) :: cfld(nxe,nye,nze) !, c(7)
      real temp(1)
      integer nxe, nye, kseed
      real*8 dmin1, dble
      real(4), intent(in) :: sdata(4) !mean and sd for (log)normal distributions, lower,upper,m,s parameters for bounded
      character, intent(in) :: stype !type of distribution n = normal, l = lognormal, b = bounded
      logical, intent(in) :: whitenoise !true for white noise, else constant deterministic 
      data zero/0.0/, half/0.5/, one/1.0/, onept5/1.5/, twopt5/2.5/
      data ienter/0/, twopi/6.283185/


      !c------------------------------------ generate the random field(s) -----------

	  kseed = iseed( kseed )

	  if(whitenoise) then !get random noise
		 
			call vnorm(cfld,nxe*nye*nze)
			
			
			
		   if(stype == 'n') then
           
				cfld = sdata(1) + cfld * sdata(2)
				
		   else if(stype == 'l') then
		   
				cfld = exp(sdata(3) + cfld * sdata(4))
				
		   else if(stype == 'b') then
		   
				cfld = sdata(1) + 0.5*(sdata(2)-sdata(1))*(1.0 + tanh((sdata(3)+sdata(4)*cfld)/6.2831853))
				
		   else 
				write(*,*) 'Incorrect distribution selected. n = normal, l = lognormal, b = bounded.'
				
		   end if
			
			
			
	  else !set uniform deterministic value
	  
		   call vnorm(temp,1)
	  
	  
		   if(stype == 'n') then
           
				temp = sdata(1) + temp * sdata(2)
				
		   else if(stype == 'l') then
		   
				temp = exp(sdata(1) + temp * sdata(2))
				
		   else if(stype == 'b') then
		   
				temp = sdata(1) + 0.5*(sdata(2)-sdata(1))*(1.0 + tanh((sdata(3)+sdata(4)*temp)/6.2831853))
				
		   else 
				write(*,*) 'Incorrect distribution selected. n = normal, l = lognormal, b = bounded.'
				
		   end if
	  
			
			cfld = temp(1)
	  end if

	  
	              


      
      return
      end
