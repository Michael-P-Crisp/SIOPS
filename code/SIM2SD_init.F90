!c  *********************************************************************
!!c  *                                                                   *
!c  *                         subroutine sim2sd                         *
!c  *                                                                   *
!c  *********************************************************************
      subroutine sim2sd_init(nxe,nye,dx,dy,thx,thy,varfnc,MXM,MXK,C0,CT,CC,CS,CI,AT,AC,AS,AI, M, k1, k2, kk)

      !real cfld(*), c(7)
      real(8) thx, thy
      integer MXM,MXK
      integer, intent(in) :: nxe, nye
      integer kseed
      character(6) varfnc
      logical debug, liid, shofld, lxfld
      real*8 dvar, dpb, dx, dy
      real(8), external :: dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      real*8 dmin1, dble
      real, intent(out) :: C0((MXK*(MXK + 1))/2)
      real, intent(out) :: CT(6,2), CC(6,4,MXM), CS(6,4,MXM), CI(6,MXM)
      real, intent(out) :: AT(3,3,2), AC(4,3,4,MXM), AS(6,3,4,MXM), AI(9,3,MXM)
      integer, intent(out) :: M,K1,K2,KK
      real(8) H, G ,da,db !for some of the more obscure correlation models
      !save XL, YL, ienter, liid
      !external dlavx2, dlsep2, dlspx2, dlafr2, dlsfr2
      !common/dparam/ dvar, dpb, dthx, dthy
      
      data zero/0.0/, half/0.5/, one/1.0/, onept5/1.5/, twopt5/2.5/
      data ienter/0/, twopi/6.283185/
      data tol/1.e-3/

   1  format(a,a)

!c-------------------------------------- initialize ---------------------

         !liid   = ((thx .eq. zero) .and. (thy .eq. zero))

         XL     = float(nxe)*dx
         YL     = float(nye)*dy

         !if( .not. liid ) then
         
         !c						initialize internal generator
      kseed = iseed( kseed )

            dvar = 1.d0
            if( varfnc .eq. 'dlafr2' .or. varfnc .eq. 'dlsfr2' ) then
               dpb = dx
               if( dy .lt. dx ) dpb = dy
            endif
            
            !if(min(thx,thy) < min(dx,dy) then
			!	write(*,*) 'warning, SOF is smaller than element. Will generate random noise.'
			!	return
				
			!else
            

				if( varfnc .eq. 'dlavx2' ) then
				   call las2i( dlavx2, nxe, nye, XL, YL, MXM, MXK ,C0, CT, CC, CS, CI, AT, AC, AS, AI, M,k1, k2, kk, iout, tol , dvar, dpb, H, G, da, db,thx,thy )

				elseif( varfnc .eq. 'dlsep2' ) then
				   call las2i( dlsep2, nxe, nye, XL, YL, MXM, MXK ,C0, CT, CC, CS, CI, AT, AC, AS, AI, M,k1, k2, kk, iout, tol , dvar, dpb, H, G, da, db,thx,thy )
				   
				elseif( varfnc .eq. 'dlspx2' ) then
				   call las2i( dlspx2, nxe, nye, XL, YL, MXM, MXK ,C0, CT, CC, CS, CI, AT, AC, AS, AI, M,k1, k2, kk, iout, tol , dvar, dpb, H, G, da, db,thx,thy )
				   
				elseif( varfnc .eq. 'dlafr2' ) then
				   call las2i( dlafr2, nxe, nye, XL, YL, MXM, MXK ,C0, CT, CC, CS, CI, AT, AC, AS, AI, M,k1, k2, kk, iout, tol , dvar, dpb, H, G, da, db,thx,thy )
				   
				elseif( varfnc .eq. 'dlsfr2' ) then
				   call las2i( dlsfr2, nxe, nye, XL, YL, MXM, MXK ,C0, CT, CC, CS, CI, AT, AC, AS, AI, M,k1, k2, kk, iout, tol , dvar, dpb, H, G, da, db,thx,thy )
				   
				endif
            !end if
            

      
      return
      end
