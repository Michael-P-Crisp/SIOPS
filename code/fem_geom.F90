module fem_geom
    
! This module contains various subroutines from the book;
! Programming the Finite Element Method, 5th edition 
! Related to the FEM subroutine
    
contains


  character(len=100) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str
  
    character(len=100) function stri(k)
!   "Convert an integer to string."
    real(8), intent(in) :: k
    write (stri, *) k
    stri = adjustl(stri)
end function stri

SUBROUTINE hexahedron_xz(iel,x_coords,y_coords,z_coords,coord,num)
!
! This subroutine generates nodal coordinates and numbering for
! 8, 14 or 20-node "bricks" counting x-z planes in the y-direction.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::iel
 REAL(8),INTENT(IN)::x_coords(:),y_coords(:),z_coords(:)  !here
 REAL(8),INTENT(OUT)::coord(:,:)                          !here
 INTEGER,INTENT(OUT)::num(:)    
 REAL(8)::pt5=0.5_iwp                                     !here
 INTEGER::fac1,fac2,ip,iq,is,iplane,nod,nxe,nze
 nxe=UBOUND(x_coords,1)-1
 nze=UBOUND(z_coords,1)-1
 nod=UBOUND(num,1)
 iq=(iel-1)/(nxe*nze)+1
 iplane=iel-(iq-1)*nxe*nze
 is=(iplane-1)/nxe+1 
 ip=iplane-(is-1)*nxe
 SELECT CASE(nod)
 CASE(8)
   fac1=(nxe+1)*(nze+1)*(iq-1)
   num(1)=fac1+is*(nxe+1)+ip
   num(2)=fac1+(is-1)*(nxe+1)+ip
   num(3)=num(2)+1
   num(4)=num(1)+1
   num(5)=(nxe+1)*(nze+1)*iq+is*(nxe+1)+ip
   num(6)=(nxe+1)*(nze+1)*iq+(is-1)*(nxe+1)+ip
   num(7)=num(6)+1
   num(8)=num(5)+1
!
   coord(1:2,1)=x_coords(ip)
   coord(5:6,1)=x_coords(ip)
   coord(3:4,1)=x_coords(ip+1)
   coord(7:8,1)=x_coords(ip+1)
!
   coord(1:4,2)=y_coords(iq)
   coord(5:8,2)=y_coords(iq+1)
!
   coord(2:3,3)=z_coords(is)
   coord(6:7,3)=z_coords(is)
   coord(1:4:3,3)=z_coords(is+1)
   coord(5:8:3,3)=z_coords(is+1)
!
   
 CASE(14)
   fac1=(2*nxe+1)*(2*nze+1)*(iq-1)
   fac2=(2*nxe+1)*(2*nze+1)*iq
   num(1)=fac1+is*(2*nxe+1)+ip
   num(2)=num(1)-(2*nxe+1)
   num(3)=num(2)+1
   num(4)=num(1)+1
   num(5)=num(2)+nxe+1
   num(6)=fac1+(nxe+1)*(nze+1)+nxe*nze+(is-1)*(2*nxe+1)+nxe+ip
   num(7)=num(6)-nxe
   num(8)=num(6)+1
   num(9)=num(8)+nxe
   num(10)=fac2+is*(2*nxe+1)+ip
   num(11)=num(10)-(2*nxe+1)
   num(12)=num(11)+1
   num(13)=num(10)+1
   num(14)=num(11)+nxe+1
!
   coord(1:2,1)=x_coords(ip)
   coord(6,1)=x_coords(ip)
   coord(10:11,1)=x_coords(ip)
   coord(5:9:2,1)=pt5*(x_coords(ip)+x_coords(ip+1))
   coord(14,1)=pt5*(x_coords(ip)+x_coords(ip+1))
   coord(3:4,1)=x_coords(ip+1)
   coord(8,1)=x_coords(ip+1)
   coord(12:13,1)=x_coords(ip+1)
!
   coord(1:5,2)=y_coords(iq)
   coord(6:9,2)=pt5*(y_coords(iq)+y_coords(iq+1))
   coord(10:14,2)=y_coords(iq+1)
!
   coord(2:3,3)=z_coords(is)
   coord(7,3)=z_coords(is)
   coord(11:12,3)=z_coords(is)
   coord(5:6,3)=pt5*(z_coords(is)+z_coords(is+1))
   coord(8:14:6,3)=pt5*(z_coords(is)+z_coords(is+1))
   coord(1:4:3,3)=z_coords(is+1)
   coord(9,3)=z_coords(is+1)
   coord(10:13:3,3)=z_coords(is+1)
!
 CASE(20)
   fac1=((2*nxe+1)*(nze+1)+(2*nze+1)*(nxe+1))*(iq-1)
   fac2=((2*nxe+1)*(nze+1)+(2*nze+1)*(nxe+1))*iq
   num(1)=fac1+(3*nxe+2)*is+2*ip-1
   num(2)=fac1+(3*nxe+2)*is-nxe+ip-1
   num(3)=num(1)-3*nxe-2
   num(4)=num(3)+1
   num(5)=num(4)+1
   num(6)=num(2)+1
   num(7)=num(1)+2
   num(8)=num(1)+1
   num(9)=fac2-(nxe+1)*(nze+1)+(nxe+1)*is+ip
   num(10)=num(9)-nxe-1
   num(11)=num(10)+1
   num(12)=num(9)+1
   num(13)=fac2+(3*nxe+2)*is+2*ip-1
   num(14)=fac2+(3*nxe+2)*is-nxe+ip-1
   num(15)=num(13)-3*nxe-2
   num(16)=num(15)+1
   num(17)=num(16)+1
   num(18)=num(14)+1
   num(19)=num(13)+2
   num(20)=num(13)+1 
!
   coord(1:3,1)=x_coords(ip)
   coord(9:10,1)=x_coords(ip)
   coord(13:15,1)=x_coords(ip)
   coord(5:7,1)=x_coords(ip+1)
   coord(11:12,1)=x_coords(ip+1)
   coord(17:19,1)=x_coords(ip+1)
   coord(4:8:4,1)=pt5*(x_coords(ip)+x_coords(ip+1))
   coord(16:20:4,1)=pt5*(x_coords(ip)+x_coords(ip+1))
!
   coord(1:8,2)=y_coords(iq)
   coord(13:20,2)=y_coords(iq+1)
   coord(9:12,2)=pt5*(y_coords(iq)+y_coords(iq+1))
!
   coord(1,3)=z_coords(is+1)
   coord(7:9,3)=z_coords(is+1)
   coord(12:13,3)=z_coords(is+1)
   coord(19:20,3)=z_coords(is+1)
   coord(3:5,3)=z_coords(is)
   coord(10:11,3)=z_coords(is)
   coord(15:17,3)=z_coords(is)
   coord(2:6:4,3)=pt5*(z_coords(is)+z_coords(is+1))
   coord(14:18:4,3)=pt5*(z_coords(is)+z_coords(is+1))
!
 CASE DEFAULT
   WRITE(11,'(a)')"Wrong number of nodes for hexahedral element"
   STOP
 END SELECT
RETURN
END SUBROUTINE hexahedron_xz


SUBROUTINE mesh_size(element,nod,nels,nn,nxe,nye,nze)
!
!  This subroutine returns the number of elements (nels) and the number
!  of nodes (nn) in a 2-d geometry-created mesh.
!
 IMPLICIT NONE
 CHARACTER(LEN=15),INTENT(IN)::element
 INTEGER,INTENT(IN)::nod,nxe,nye
 INTEGER,INTENT(IN),OPTIONAL::nze
 INTEGER,INTENT(OUT)::nels,nn
 IF(element=="triangle")THEN
   nels=nxe*nye*2
   IF(nod==3)nn=(nxe+1)*(nye+1)
   IF(nod==6)nn=(2*nxe+1)*(2*nye+1)
   IF(nod==10)nn=(3*nxe+1)*(3*nye+1)
   IF(nod==15)nn=(4*nxe+1)*(4*nye+1)
 ELSE IF(element=="quadrilateral")THEN
   nels=nxe*nye
   IF(nod==4)nn=(nxe+1)*(nye+1)
   IF(nod==5)nn=(nxe+1)*(nye+1)+nxe*nye
   IF(nod==8)nn=(2*nxe+1)*(nye+1)+(nxe+1)*nye
   IF(nod==9)nn=(2*nxe+1)*(2*nye+1)
 ELSE IF(element=="hexahedron")THEN
   nels=nxe*nye*nze
   IF(nod==8)nn=(nxe+1)*(nye+1)*(nze+1)
   IF(nod==14)nn=4*nxe*nye*nze+2*(nxe*nye+nye*nze+nze*nxe)+nxe+nye+nze+1
   IF(nod==20)nn=((2*nxe+1)*(nze+1)+(nxe+1)*nze)*(nye+1)+                 &
     (nxe+1)*(nze+1)*nye
 END IF
RETURN
END SUBROUTINE mesh_size

SUBROUTINE formnf(nf)
!
! This subroutine forms the nf matrix.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN OUT)::nf(:,:)
 INTEGER::i,j,m
 m=0
 DO j=1,UBOUND(nf,2)
   DO i=1,UBOUND(nf,1)
     IF(nf(i,j)/=0)THEN
       m=m+1
       nf(i,j)=m
     END IF
   END DO
 END DO
RETURN
END SUBROUTINE formnf 

SUBROUTINE mesh_ensi(argv,nlen,g_coord,g_num,element,property,nf,loads,          &
                     nstep,npri,dtim,solid,it,jr,ploop,slash)
!
! This subroutine outputs a set of files in the Ensight gold format.
! Models in this format can be viewed in the visualisation tool ParaView.
!
! Element types supported:                Tested with:
!
! 2-node bar
! 3-node triangle                         p51  (4th edition p51_1.dat)
! 6-node triangle
! 4-node quadrilateral                    p115 (4th edition)
! 8-node quadrilateral                    p116 (4th edition)
! 4-node tetrahedron                      p54  (4th edition p54_2.dat)
! 8-node hexahedron                       p86  (4th edition)
! 20-node hexahedron                      p55  (4th edition)

  IMPLICIT none

  INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
  INTEGER,   INTENT(IN)         :: nlen,nstep,npri,it,jr,ploop
  INTEGER,   INTENT(IN)         :: g_num(:,:),nf(:,:)
  INTEGER                       :: i,j,k,l,m,n,nfe,nod,nels,ndim,nn
  INTEGER                       :: prnwidth,remainder
  REAL(8), INTENT(IN)         :: loads(:),dtim   !herex
  real(4), INTENT(IN)           :: g_coord(:,:)
  CHARACTER(LEN=100), INTENT(IN) :: argv,element
  LOGICAL, INTENT(IN)           :: solid
  real(4), INTENT(IN)           :: property(:)
  character                     slash

!------------------------------------------------------------------------------
! 1. Initialisation
!------------------------------------------------------------------------------

  nn   = UBOUND(g_coord,2) ; ndim = UBOUND(g_coord,1)
  nels = UBOUND(g_num,2)   ; nod  = UBOUND(g_num,1)

!------------------------------------------------------------------------------
! 2. Write case file
!------------------------------------------------------------------------------



  OPEN(12,FILE='paraview'//slash//trim(argv)//trim(str(it))//'_it_'//trim(str(jr))//'_pd_'//trim(str(ploop))//'.ensi.case')

  WRITE(12,'(A/A)')    "#", "# Post-processing file generated by subroutine &
                             &WRITE_ENSI in "
  WRITE(12,'(A,A,/A)') "#", " Smith, Griffiths and Margetts, 'Programming the &
                             &Finite Element Method',","# Wiley, 2013."
  WRITE(12,'(A/A/A)')  "#","# Ensight Gold Format","#"
  WRITE(12,'(2A/A/A/A/A/A/A)')   "# Problem name: ",trim(argv)//trim(str(it))//'_it_'//trim(str(jr))//'_pd_'//trim(str(ploop))//"#"
  WRITE(12,'(A/A/A)')  "FORMAT","type:  ensight gold","GEOMETRY"
  WRITE(12,'(2A/A)')   "model: 1  ",trim(argv)//'.ensi.geo',"VARIABLE"
  WRITE(12,'(2A)')     "scalar per element:  material      ",                 &
                        trim(argv)//trim(str(it))//'_it_'//trim(str(jr))//'_pd_'//trim(str(ploop))//'.ensi.matid'
  IF(solid) THEN
    WRITE(12,'(2A)')   "scalar per node:     restraint     ",                 &
                        trim(argv)//'.ensi.ndbnd'
    WRITE(12,'(2A)')   "vector per node:     displacement  ",                 &
                        trim(argv)//trim(str(it))//'_it_'//trim(str(jr))//'_pd_'//trim(str(ploop))//'.ensi.displ-*****'
  ELSE
    WRITE(12,'(2A)')   "scalar per node:     pressure      ",                 &
                        trim(argv)//'.ensi.pressure-*****'
  END IF
  WRITE(12,'(2A)')     "vector per node:     load          ",                 &
                        trim(argv)//'.ensi.ndlds'
  WRITE(12,'(A/A)')     "TIME","time set:     1"
  WRITE(12,'(A,I5)')    "number of steps:",nstep/npri
  WRITE(12,'(A,I5)')    "filename start number:",npri
  WRITE(12,'(A,I5)')    "filename increment:",npri
  WRITE(12,'(A)')       "time values:"
  prnwidth  = 5
  remainder = mod(nstep/npri,prnwidth)
  n         = ((nstep/npri) - remainder)/prnwidth
  IF(nstep/npri<=prnwidth) THEN
    DO i=1,nstep,npri
      IF(i==nstep) THEN
        WRITE(12,'(E12.5)') i*dtim
      ELSE
        WRITE(12,'(E12.5)',ADVANCE='no') i*dtim
      END IF
    END DO
  ELSE
    IF(remainder==0) THEN
      DO j=1,n
        m = ((j-1)*prnwidth)+1
        l = ((j-1)*prnwidth)+prnwidth
        WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
      END DO
    ELSE
!     DO j=1,n-1
      DO j=1,n
        m = ((j-1)*prnwidth)+1
        l = ((j-1)*prnwidth)+prnwidth
        WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
      END DO
      m = (n*prnwidth)+1
      l = (n*prnwidth)+remainder
      DO i=m,l
        IF(i==l) THEN
          WRITE(12,'(E12.5)') dtim*i
        ELSE
          WRITE(12,'(E12.5)',ADVANCE='no') dtim*i
        END IF
      END DO
    END IF
  END IF

  CLOSE(12)

!------------------------------------------------------------------------------
! 3. Write geometry file
!------------------------------------------------------------------------------
if(it==1 .or. it==2) then
  OPEN(13,FILE='paraview'//slash//trim(argv)//'.ensi.geo')
  WRITE(13,'(/2A)')   "Problem name: ", trim(argv)
  WRITE(13,'(A/A/A)') "Geometry files","node id given","element id given"
  WRITE(13,'(A/A)')   "part","      1"
  IF(ndim==2) WRITE(13,'(A)') "2d-mesh"
  IF(ndim==3) WRITE(13,'(A)') "Volume Mesh"
  WRITE(13,'(A)')     "coordinates"

  WRITE(13,'(I10)') nn
  DO j=1,ndim
    DO i=1,nn
      WRITE(13,'(E12.5)') g_coord(j,i)
    END DO
  END DO

  IF(ndim==2) THEN ! ensight requires zeros for the z-ordinate
    DO i=1,nn
      WRITE(13,'(A)') " 0.00000E+00"
    END DO
  END IF

  SELECT CASE(element)
    CASE('triangle')
      SELECT CASE(nod)
        CASE(3)
          WRITE(13,'(A/I10)') "tria3", nels
          DO i = 1,nels
            WRITE(13,'(3I10)')g_num(3,i),g_num(2,i),g_num(1,i)
          END DO
        CASE DEFAULT
          WRITE(13,'(A)')   "# Element type not recognised"
      END SELECT
    CASE('quadrilateral')
      SELECT CASE(nod)
        CASE(4)
          WRITE(13,'(A/I10)') "quad4", nels
          DO i = 1,nels
            WRITE(13,'(4I10)')g_num(1,i),g_num(4,i),g_num(3,i),g_num(2,i)
          END DO
        CASE(8)
          WRITE(13,'(A/I10)') "quad8", nels
          DO i = 1,nels
            WRITE(13,'(8I10)')g_num(1,i),g_num(7,i),g_num(5,i),g_num(3,i),    &
                              g_num(8,i),g_num(6,i),g_num(4,i),g_num(2,i)
          END DO
        CASE DEFAULT
          WRITE(13,'(A)')   "# Element type not recognised"
      END SELECT
    CASE('hexahedron')
      SELECT CASE(nod)
        CASE(8)
          WRITE(13,'(A/I10)') "hexa8", nels
          DO i = 1,nels
            WRITE(13,'(8I10)') g_num(1,i),g_num(4,i),g_num(8,i),g_num(5,i),   &
                               g_num(2,i),g_num(3,i),g_num(7,i),g_num(6,i)
          END DO
        CASE(20)
          WRITE(13,'(A/I10)') "hexa20", nels
          DO i = 1,nels
            WRITE(13,'(20I10)')                                               &
              g_num(1,i), g_num(7,i), g_num(19,i),g_num(13,i),g_num(3,i),     &
              g_num(5,i), g_num(17,i),g_num(15,i),g_num(8,i), g_num(12,i),    &
              g_num(20,i),g_num(9,i), g_num(4,i), g_num(11,i),g_num(16,i),    &
              g_num(10,i),g_num(2,i), g_num(6,i), g_num(18,i),g_num(14,i)
          END DO
        CASE DEFAULT
          WRITE(13,'(A)')   "# Element type not recognised"
      END SELECT
    CASE('tetrahedron')
      SELECT CASE(nod)
        CASE(4)
          WRITE(13,'(A/I10)') "tetra4", nels
          DO i = 1,nels
            WRITE(13,'(4I10)') g_num(1,i),g_num(3,i),g_num(2,i),g_num(4,i)
          END DO
        CASE DEFAULT
          WRITE(13,'(A)')   "# Element type not recognised"
      END SELECT
    CASE DEFAULT
      WRITE(13,'(A)')       "# Element type not recognised"
  END SELECT

  CLOSE(13)
end if
!------------------------------------------------------------------------------
! 4. Write file containing material IDs
!------------------------------------------------------------------------------

  OPEN(14,FILE='paraview'//slash//trim(argv)//trim(str(it))//'_it_'//trim(str(jr))//'_pd_'//trim(str(ploop))//'.ensi.matid')
  WRITE(14,'(A)') "Alya Ensight Gold --- Scalar per-element variable file"
  WRITE(14,'(A/A)') "part", "      1"

  SELECT CASE(element)
    CASE('triangle')
      SELECT CASE(nod)
        CASE(3)
          WRITE(14,'(A)') "tria3"
        CASE DEFAULT
          WRITE(14,'(A)') "# Element type not recognised"
      END SELECT
    CASE('quadrilateral')
      SELECT CASE(nod)
        CASE(4)
          WRITE(14,'(A)') "quad4"
        CASE(8)
          WRITE(14,'(A)') "quad8"
        CASE DEFAULT
          WRITE(14,'(A)') "# Element type not recognised"
      END SELECT
    CASE('hexahedron')
      SELECT CASE(nod)
        CASE(8)
          WRITE(14,'(A)') "hexa8"
        CASE(20)
          WRITE(14,'(A)') "hexa20"
        CASE DEFAULT
          WRITE(14,'(A)') "# Element type not recognised"
      END SELECT
    CASE('tetrahedron')
      SELECT CASE(nod)
        CASE(4)
          WRITE(14,'(A)') "tetra4"
        CASE DEFAULT
        WRITE(14,'(A)') "# Element type not recognised"
      END SELECT
    CASE DEFAULT
      WRITE(14,'(A)')   "# Element type not recognised"
  END SELECT

  DO i=1,nels; WRITE(14,'(F12.5)') property(i); END DO
      !'(I10)'

  WRITE(14,'(A)')

  CLOSE(14)

!------------------------------------------------------------------------------
! 5. Write boundary conditions. Encoded using formula: 4z + 2y + 1x
!
!    110 = 1   010 = 2   100 = 3   011 = 4   101 = 5   001 = 6   000 = 7
!------------------------------------------------------------------------------
if(it==1) then
  IF(solid) THEN
    OPEN(15,FILE='paraview'//slash//trim(argv)//'.ensi.ndbnd')
    WRITE(15,'(A)')     "Alya Ensight Gold --- Scalar per-node variable file"
    WRITE(15,'(A/A/A)') "part", "      1","coordinates"
    IF(ndim==3) THEN
      DO i=1,UBOUND(g_coord,2)
        nfe=0
        IF(nf(1,i)==0) nfe=nfe+1
        IF(nf(2,i)==0) nfe=nfe+2
        IF(nf(3,i)==0) nfe=nfe+4
        WRITE(15,'(I2)') nfe
      END DO
    ELSE IF(ndim==2) THEN
      DO i=1,nn
        nfe=0
        IF(nf(1,i)==0) nfe=nfe+1
        IF(nf(2,i)==0) nfe=nfe+2
        WRITE(15,'(I2)') nfe
      END DO
    ELSE
      PRINT *, "Wrong number of dimensions in mesh_ensi"
    END IF
  END IF

  CLOSE(15)

!------------------------------------------------------------------------------
! 6. Write loaded nodes
!------------------------------------------------------------------------------

  OPEN(16,FILE='paraview'//slash//trim(argv)//'.ensi.ndlds')
  WRITE(16,'(A)')     "Alya Ensight Gold --- Vector per-node variable file"
  WRITE(16,'(A/A/A)') "part", "      1","coordinates"
  DO j=1,UBOUND(nf,1)
    DO i=1, UBOUND(nf,2)
        if(nf(j,i)==0) then
            WRITE(16,'(E12.5)') 0.0
        else
            WRITE(16,'(E12.5)') loads(nf(j,i))
        end if
    END DO
  END DO
  CLOSE(16)
  
  end if
  

  RETURN

  END SUBROUTINE mesh_ensi
  

  
  
  

  SUBROUTINE checon(loads,oldlds,tol,converged)
!
! This subroutine sets converged to .FALSE. if relative change in loads
! and oldlds is greater than tol and updates oldlds.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(4),INTENT(IN)::loads(0:)            !herex
 real(8),INTENT(IN):: tol
 real(8) :: temp
 !integer, intent(in) :: iloop
 REAL(4),INTENT(IN OUT)::oldlds(0:)      !herex
 LOGICAL,INTENT(OUT)::converged
 !real(4), intent(in out) :: tlvec(:)
 CONVERGED=.TRUE.
 !tlvec(:size(tlvec)) = tlvec(2:) !increment everything down
 temp=MAXVAL(ABS(loads-oldlds))/MAXVAL(ABS(loads))
 !tlvec(size(tlvec)) = temp
 CONVERGED=(temp <= tol)
 !oldlds=loads
RETURN
  END SUBROUTINE checon
  
  
  SUBROUTINE dismsh_ensi(argv,nlen,step,nf,loads,neq,it,jr,ploop)
!
! This subroutine outputs displacements in the Ensight gold format for
! visualization in ParaView. ParaView also requires the output of subroutine
! mesh_ensi for the geometry.
!
  IMPLICIT none

  INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
  INTEGER,   INTENT(IN)         :: nlen,step,nf(:,:),neq,it,jr,ploop
  INTEGER                       :: i,j
  REAL(4), INTENT(IN)         :: loads(0:ubound(nf,2)*ubound(nf,1))
  CHARACTER(LEN=15), INTENT(IN) :: argv
  CHARACTER(LEN=5)              :: ch

  WRITE(ch,'(I5.5)') step ! convert integer to string using internal file

  OPEN(17,FILE='paraview\'//trim(argv)//trim(str(it))//'_it_'//trim(str(jr))//'_pd_'//trim(str(ploop))//'.ensi.displ-'//ch)

  WRITE(17,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
  WRITE(17,'(A/A/A)') "part", "      1","coordinates"

 !write(*,*) ubound(nf,2)
  
    DO j=1,UBOUND(nf,1)
    DO i=1, UBOUND(nf,2)
!        if(nf(j,i)==0) then
!            WRITE(17,'(E12.5)') 0.0
!        else
            WRITE(17,'(E12.5)') loads(nf(j,i))
!        end if
    END DO
    END DO
  

  IF(UBOUND(nf,1)==2) THEN ! ensight requires zeros for the z-ordinate
    DO i=1,UBOUND(nf,2)
      WRITE(17,'(A)') " 0.00000E+00"
    END DO
  END IF

  CLOSE(17)

  RETURN

END SUBROUTINE dismsh_ensi

!subroutine dismsh_ensi(argv,nlen,step,nf,loads,it,jr,ploop)
!!
!! this subroutine outputs displacements in the ensight gold format for
!! visualization in paraview. paraview also requires the output of subroutine
!! mesh_ensi for the geometry.
!!
!  implicit none
!
!  integer,parameter             :: iwp=selected_real_kind(15)
!  integer,   intent(in)         :: nlen,step,nf(:,:),it,jr,ploop
!  integer                       :: i,j
!  real(8), intent(in)         :: loads(:)     !here
!  character(len=100), intent(in) :: argv
!  character(len=5)              :: ch           !here
!
!  write(ch,'(i5.5)') step ! convert integer to string using internal file
!
!  open(17,file='paraview\'//trim(argv)//trim(str(it))//'_it_'//trim(str(jr))//'_pd_'//trim(str(ploop))//'.ensi.displ-'//ch)
!
!  write(17,'(a)') "alya ensight gold --- vector per-node variable file"
!  write(17,'(a/a/a)') "part", "      1","coordinates"
!
! 
!    do j=1,ubound(nf,1)
!    do i=1, ubound(nf,2)
!        if(nf(j,i)==0) then
!            write(17,'(e12.5)') 0.0
!        else
!            write(17,'(e12.5)') loads(nf(j,i))
!        end if
!    end do
!    end do
!  
!
!  if(ubound(nf,1)==2) then ! ensight requires zeros for the z-ordinate
!    do i=1,ubound(nf,2)
!      write(17,'(a)') " 0.00000e+00"
!    end do
!  end if
!
!  close(17)
!
!  return
!
!end subroutine dismsh_ensi

SUBROUTINE geom_rect(element,iel,x_coords,y_coords,coord,num,dir)
!
! This subroutine forms the coordinates and connectivity for a
! rectangular mesh of rt. angled triangular elements (3, 6, 10 or 15-node)
! or quadrilateral elements (4, 8 or 9-node) counting in the
! x- or y-dir. 
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x_coords(:),y_coords(:)
 REAL(iwp),INTENT(OUT)::coord(:,:)
 CHARACTER(LEN=15),INTENT(IN)::element
 CHARACTER(LEN=1),INTENT(IN)::dir
 INTEGER,INTENT(IN)::iel
 INTEGER,INTENT(OUT)::num(:)
 INTEGER::ip,iq,jel,fac1,nod,nxe,nye
 REAL(iwp)::pt5=0.5_iwp,two=2.0_iwp,d3=3.0_iwp 
 nxe=UBOUND(x_coords,1)-1
 nod=UBOUND(num,1)
 IF(element=='triangle')THEN
   nye=(UBOUND(y_coords,1)-1)*2
   IF(dir=='x'.OR.dir=='r')THEN
     jel=2*nxe*((iel-1)/(2*nxe))
     ip=(iel-jel+1)/2
     iq=2*((iel-1)/(2*nxe)+1)-1+((iel/2)*2)/iel
   ELSE  
     jel=(iel-1)/nye
     ip=jel+1
     iq=iel-nye*jel
   END IF
   SELECT CASE(nod)
   CASE(3)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(nxe+1)*(iq-1)/2+ip
         num(2)=num(1)+1              
         num(3)=(nxe+1)*(iq+1)/2+ip
       ELSE
         num(1)=(ip-1)*(nye+2)/2+(iq+1)/2
         num(2)=num(1)+(nye+2)/2
         num(3)=num(1)+1
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(1,2)=y_coords((iq+1)/2)
       coord(2,1)=x_coords(ip+1)   
       coord(2,2)=y_coords((iq+1)/2)
       coord(3,1)=x_coords(ip)   
       coord(3,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(nxe+1)*iq/2+ip+1     
         num(2)=num(1)-1               
         num(3)=(nxe+1)*(iq-2)/2+ip+1
       ELSE
         num(1)=ip*(nye+2)/2+(iq+2)/2
         num(2)=(ip-1)*(nye+2)/2+(iq+1)/2+1
         num(3)=num(1)-1
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(1,2)=y_coords((iq+2)/2)
       coord(2,1)=x_coords(ip)   
       coord(2,2)=y_coords((iq+2)/2)
       coord(3,1)=x_coords(ip+1) 
       coord(3,2)=y_coords(iq/2)
     END IF
   CASE(6)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(iq-1)*(2*nxe+1)+2*ip-1
         num(2)=num(1)+1 
         num(3)=num(1)+2 
         num(4)=(iq-1)*(2*nxe+1)+2*nxe+2*ip+1
         num(5)=(iq+1)*(2*nxe+1)+2*ip-1
         num(6)=num(4)-1 
       ELSE
         num(1)=2*(nye+1)*(ip-1)+iq
         num(2)=2*(nye+1)*(ip-1)+nye+1+iq
         num(3)=2*(nye+1)*ip+iq
         num(4)=num(2)+1
         num(5)=num(1)+2 
         num(6)=num(1)+1
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(1,2)=y_coords((iq+1)/2)
       coord(3,1)=x_coords(ip+1)   
       coord(3,2)=y_coords((iq+1)/2)
       coord(5,1)=x_coords(ip)   
       coord(5,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=iq*(2*nxe+1)+2*ip+1
         num(2)=num(1)-1 
         num(3)=num(1)-2 
         num(4)=(iq-2)*(2*nxe+1)+2*nxe+2*ip+1
         num(5)=(iq-2)*(2*nxe+1)+2*ip+1
         num(6)=num(4)+1 
       ELSE 
         num(1)=2*(nye+1)*ip+iq+1 
         num(2)=2*(nye+1)*(ip-1)+nye+iq+2
         num(3)=2*(nye+1)*(ip-1)+iq+1
         num(4)=num(2)-1 
         num(5)=num(1)-2
         num(6)=num(1)-1
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(1,2)=y_coords((iq+2)/2)
       coord(3,1)=x_coords(ip)   
       coord(3,2)=y_coords((iq+2)/2)
       coord(5,1)=x_coords(ip+1) 
       coord(5,2)=y_coords(iq/2)
     END IF
     coord(2,:)=pt5*(coord(1,:)+coord(3,:))
     coord(4,:)=pt5*(coord(3,:)+coord(5,:))
     coord(6,:)=pt5*(coord(5,:)+coord(1,:))
   CASE(10)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(iq-1)/2*(3*nxe+1)*3+3*ip-2
         num(2)=num(1)+1
         num(3)=num(1)+2
         num(4)=num(1)+3
         num(5)=(iq-1)/2*(3*nxe+1)*3+3*nxe+1+3*ip
         num(6)=(iq-1)/2*(3*nxe+1)*3+6*nxe+2+3*ip-1
         num(7)=(iq-1)/2*(3*nxe+1)*3+9*nxe+3+3*ip-2
         num(8)=num(6)-1
         num(9)=num(5)-2
         num(10)=num(9)+1
       ELSE
         num(1)=(9*(nye-2)/2+12)*(ip-1)+3*(iq-1)/2+1
         num(2)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)/2+4+3*(iq-1)/2+1
         num(3)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)+8+3*(iq-1)/2+1
         num(4)=(9*(nye-2)/2+12)*(ip-1)+9*(nye-2)/2+12+3*(iq-1)/2+1
         num(5)=num(3)+1 
         num(6)=num(2)+2
         num(7)=num(1)+3
         num(8)=num(1)+2
         num(9)=num(1)+1
         num(10)=num(2)+1
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(2,1)=x_coords(ip)+(x_coords(ip+1)-x_coords(ip))/d3
       coord(3,1)=x_coords(ip)+two*(x_coords(ip+1)-x_coords(ip))/d3
       coord(4,1)=x_coords(ip+1)
       coord(4,2)=y_coords((iq+1)/2)
       coord(5,2)=y_coords((iq+1)/2)+                                     &
         (y_coords((iq+3)/2)-y_coords((iq+1)/2))/d3
       coord(6,2)=y_coords((iq+1)/2)+                                     &
         two*(y_coords((iq+3)/2)-y_coords((iq+1)/2))/d3
       coord(7,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(iq-2)/2*(3*nxe+1)*3+9*nxe+3+3*ip+1
         num(2)=num(1)-1
         num(3)=num(1)-2
         num(4)=num(1)-3
         num(5)=(iq-2)/2*(3*nxe+1)*3+6*nxe+2+3*ip-1
         num(6)=(iq-2)/2*(3*nxe+1)*3+3*nxe+1+3*ip
         num(7)=(iq-2)/2*(3*nxe+1)*3+3*ip+1
         num(8)=num(6)+1
         num(9)=num(5)+2
         num(10)=num(9)-1
       ELSE
         num(1)=(9*(nye-2)/2+12)*(ip-1)+9*(nye-2)/2+12+3*iq/2+1
         num(2)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)+8+3*iq/2+1
         num(3)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)/2+4+3*iq/2+1
         num(4)=(9*(nye-2)/2+12)*(ip-1)+3*iq/2+1
         num(5)=num(3)-1
         num(6)=num(2)-2
         num(7)=num(1)-3
         num(8)=num(1)-2
         num(9)=num(1)-1
         num(10)=num(2)-1
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(2,1)=x_coords(ip+1)-(x_coords(ip+1)-x_coords(ip))/d3
       coord(3,1)=x_coords(ip+1)-two*(x_coords(ip+1)-x_coords(ip))/d3
       coord(4,1)=x_coords(ip)
       coord(4,2)=y_coords((iq+2)/2)
       coord(5,2)=y_coords((iq+2)/2)-(y_coords((iq+2)/2)-y_coords(iq/2))/d3
       coord(6,2)=y_coords((iq+2)/2)-                                     &
         two*(y_coords((iq+2)/2)-y_coords(iq/2))/d3
       coord(7,2) =y_coords(iq/2)
     END IF
     coord(5,1)=coord(3,1)
     coord(6,1)=coord(2,1)
     coord(7,1)=coord(1,1)
     coord(8,1)=coord(1,1)
     coord(9,1)=coord(1,1)
     coord(10,1)=coord(2,1)
     coord(1,2)=coord(4,2)
     coord(2,2)=coord(4,2)
     coord(3,2)=coord(4,2)
     coord(8,2)=coord(6,2)
     coord(9,2)=coord(5,2)
     coord(10,2)=coord(5,2)
   CASE(15)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
       fac1=4*(4*nxe+1)*(iq-1)/2
         num(1)=fac1+4*ip-3
         num(2)=num(1)+1
         num(3)=num(1)+2
         num(4)=num(1)+3
         num(5)=num(1)+4
         num(6)=fac1+ 4*nxe+1+4*ip
         num(7)=fac1+ 8*nxe+1+4*ip
         num(8)=fac1+12*nxe+1+4*ip
         num(9)=fac1+16*nxe+1+4*ip
         num(10)=num(8)-1
         num(11)=num(7)-2
         num(12)=num(6)-3
         num(13)=num(12)+1
         num(14)=num(12)+2
         num(15)=num(11)+1
       ELSE
         fac1=4*(2*nye+1)*(ip-1)+2*iq-1 
         num(1)=fac1
         num(2)=fac1+2*nye+1
         num(3)=fac1+4*nye+2 
         num(4)=fac1+6*nye+3 
         num(5)=fac1+8*nye+4
         num(6)=fac1+6*nye+4 
         num(7)=fac1+4*nye+4 
         num(8)=fac1+2*nye+4
         num(9)=fac1+4 
         num(10)=fac1+3 
         num(11)=fac1+2 
         num(12)=fac1+1
         num(13)=fac1+2*nye+2 
         num(14)=fac1+4*nye+3
         num(15)=fac1+2*nye+3  
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(1,2)=y_coords((iq+1)/2)
       coord(5,1)=x_coords(ip+1)   
       coord(5,2)=y_coords((iq+1)/2)
       coord(9,1)=x_coords(ip)   
       coord(9,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         fac1=4*(4*nxe+1)*(iq-2)/2
         num(1)=fac1+16*nxe+5+4*ip
         num(2)=num(1)-1
         num(3)=num(1)-2
         num(4)=num(1)-3
         num(5)=num(1)-4
         num(6)=fac1+12*nxe+1+4*ip
         num(7)=fac1+8*nxe+1+4*ip
         num(8)=fac1+4*nxe+1+4*ip
         num(9)=fac1+4*ip+1
         num(10)=num(8)+1
         num(11)=num(7)+2
         num(12)=num(6)+3
         num(13)=num(12)-1
         num(14)=num(12)-2
         num(15)=num(11)-1
       ELSE
         fac1=4*(2*nye+1)*(ip-1)+2*iq+8*nye+5 
         num(1)=fac1 
         num(2)=fac1-2*nye-1
         num(3)=fac1-4*nye-2 
         num(4)=fac1-6*nye-3 
         num(5)=fac1-8*nye-4
         num(6)=fac1-6*nye-4  
         num(7)=fac1-4*nye-4 
         num(8)=fac1-2*nye-4
         num(9)=fac1-4
         num(10)=fac1-3 
         num(11)=fac1-2 
         num(12)=fac1-1
         num(13)=fac1-2*nye-2  
         num(14)=fac1-4*nye-3
         num(15)=fac1-2*nye-3 
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(1,2)=y_coords((iq+2)/2)
       coord(5,1)=x_coords(ip)   
       coord(5,2)=y_coords((iq+2)/2)
       coord(9,1)=x_coords(ip+1) 
       coord(9,2)=y_coords(iq/2)
     END IF
     coord(3,:)=pt5*(coord(1,:)+coord(5,:))
     coord(7,:)=pt5*(coord(5,:)+coord(9,:))
     coord(11,:)=pt5*(coord(9,:)+coord(1,:))
     coord(2,:)=pt5*(coord(1,:)+coord(3,:))
     coord(4,:)=pt5*(coord(3,:)+coord(5,:))
     coord(6,:)=pt5*(coord(5,:)+coord(7,:))
     coord(8,:)=pt5*(coord(7,:)+coord(9,:))
     coord(10,:)=pt5*(coord(9,:)+coord(11,:))
     coord(12,:)=pt5*(coord(11,:)+coord(1,:))
     coord(15,:)=pt5*(coord(7,:)+coord(11,:))
     coord(14,:)=pt5*(coord(3,:)+coord(7,:))
     coord(13,:)=pt5*(coord(2,:)+coord(15,:))
   CASE DEFAULT
     WRITE(11,'(a)')"Wrong number of nodes for triangular element"
     STOP
   END SELECT
 ELSE !quadrilateral
   nye=UBOUND(y_coords,1)-1
   IF(dir=='x'.OR.dir=='r')THEN
     iq=(iel-1)/nxe+1
     ip=iel-(iq-1)*nxe
   ELSE
     ip=(iel-1)/nye+1
     iq=iel-(ip-1)*nye
   END IF
   SELECT CASE(nod)
   CASE(4)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(nxe+1)+ip		        		
       num(2)=(iq-1)*(nxe+1)+ip				
       num(3)=num(2)+1					
       num(4)=num(1)+1					
     ELSE
       num(1)=(ip-1)*(nye+1)+iq+1
       num(2)=num(1)-1
       num(3)=ip*(nye+1)+iq
       num(4)=num(3)+1
     END IF
!
     coord(1:2,1)=x_coords(ip)
     coord(3:4,1)=x_coords(ip+1)
     coord(1,2)=y_coords(iq+1)
     coord(2:3,2)=y_coords(iq)
     coord(4,2)=coord(1,2)
   CASE(5)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(2*nxe+1)+ip
       num(2)=(iq-1)*(2*nxe+1)+ip
       num(3)=num(2)+1
       num(4)=num(1)+1
       num(5)=iq*(2*nxe+1)+ip-nxe
     ELSE
       num(1)=(ip-1)*(2*nye+1)+iq+1
       num(2)=num(1)-1
       num(3)=ip*(2*nye+1)+iq
       num(4)=num(3)+1
       num(5)=ip*(2*nye+1)+iq-nye
     END IF
!
     coord(1:2,1)=x_coords(ip)
     coord(3:4,1)=x_coords(ip+1)
     coord(1,2)=y_coords(iq+1)
     coord(2:3,2)=y_coords(iq)
     coord(4,2)=coord(1,2)
     coord(5,:)=0.25_iwp*(coord(1,:)+coord(2,:)+coord(3,:)+coord(4,:))
   CASE(8)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(3*nxe+2)+2*ip-1                 
       num(2)=iq*(3*nxe+2)+ip-nxe-1		  
       num(3)=(iq-1)*(3*nxe+2)+2*ip-1		   
       num(4)=num(3)+1
       num(5)=num(4)+1
       num(6)=num(2)+1
       num(7)=num(1)+2
       num(8)=num(1)+1
     ELSE
       num(1)=(ip-1)*(3*nye+2)+2*iq+1
       num(2)=num(1)-1
       num(3)=num(1)-2
       num(4)=(ip-1)*(3*nye+2)+2*nye+iq+1
       num(5)=ip*(3*nye+2)+2*iq-1
       num(6)=num(5)+1
       num(7)=num(5)+2
       num(8)=num(4)+1
     END IF
!
     coord(1:3,1)=x_coords(ip)
     coord(5:7,1)=x_coords(ip+1)
     coord(4,1)=pt5*(coord(3,1)+coord(5,1))
     coord(8,1)=pt5*(coord(7,1)+coord(1,1))
     coord(1,2)=y_coords(iq+1)
     coord(7:8,2)=y_coords(iq+1)
     coord(3:5,2)=y_coords(iq)
     coord(2,2)=pt5*(coord(1,2)+coord(3,2))
     coord(6,2)=pt5*(coord(5,2)+coord(7,2))
   CASE(9)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(4*nxe+2)+2*ip-1
       num(2)=iq*(4*nxe+2)+2*ip-nxe-4
       num(3)= (iq-1)*(4*nxe+2)+2*ip-1
       num(4)=num(3)+1
       num(5)=num(4)+1
       num(6)=num(2)+2
       num(7)=num(1)+2
       num(8)=num(1)+1
       num(9)=num(2)+1
     ELSE
       num(1)=(ip-1)*2*(2*nye+1)+2*iq+1
       num(2)=num(1)-1
       num(3)=num(1)-2
       num(4)=(ip-1)*2*(2*nye+1)+2*nye+2*iq
       num(5)=ip*2*(2*nye+1)+2*iq-1
       num(6)=num(5)+1
       num(7)=num(5)+2
       num(8)=num(4)+2
       num(9)=num(4)+1
     END IF
!
     coord(1:3,1)=x_coords(ip)
     coord(5:7,1)=x_coords(ip+1)
     coord(4,1)=pt5*(coord(3,1)+coord(5,1))
     coord(8,1)=pt5*(coord(7,1)+coord(1,1))
     coord(1,2)=y_coords(iq+1)
     coord(7:8,2)=y_coords(iq+1)
     coord(3:5,2)=y_coords(iq)
     coord(2,2)=pt5*(coord(1,2)+coord(3,2))
     coord(6,2)=pt5*(coord(5,2)+coord(7,2))
     coord(9,:)=pt5*(coord(4,:)+coord(8,:))
   CASE DEFAULT
     WRITE(11,'(a)')"Wrong number of nodes for quadrilateral element"
     STOP
   END SELECT
 END IF
RETURN
END SUBROUTINE geom_rect
    
    end module
    
