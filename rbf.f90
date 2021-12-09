!
!*** END MLS
!  
!---------------------------
!*** RADIUS AND DERIVATIVES
!*** DETERMINED
!*** CHK0
!---------------------------
  SUBROUTINE RBF_RADIUS(NDI,X,XI,R,DRDX,D2RDX2)
    IMPLICIT REAL(8) (A-H,O-Z)
    REAL(8),DIMENSION(NDI)::X,XI
    REAL(8),DIMENSION(NDI)::DRDX
    REAL(8),DIMENSION(NDI,NDI)::D2RDX2
    R=RNORM2(NDI,X-XI)
    DRDX=0.0D00
    D2RDX2=0.0D00
    IF(R.GT.1.0D-30)THEN
       DO ID=1,NDI
          DRDX(ID)=(X(ID)-XI(ID))/R
          DO JD=1,NDI
             D2RDX2(ID,JD)=DELTAK(ID,JD)/R-(1.0D00/(R**3.0D00))*(X(ID)-XI(ID))*(X(JD)-XI(JD))
          END DO
       END DO
    END IF
  END SUBROUTINE RBF_RADIUS
  
!--------------------------
!*** RADIAL BASIS FUNCTION
!*** BASIS FUNCTION
!*** chk0
!--------------------------
  SUBROUTINE RBF_C2(RMAX,R,A,DADR,D2ADR2)
    IMPLICIT REAL(8) (A-H,O-Z)
    REAL(8)::RMAX,R,A,DADR,D2ADR2
    A=0.0D00
    DADR=0.0D00
    D2ADR2=0.0D00
    IF(R.LE.RMAX)THEN
       A=((1.0d00-(R/RMAX))**4)*(1.0D00+4.0d00*(R/RMAX))
       DADR=(20.0d00*R*(R-RMAX)**3.0D00)/(RMAX**5.0D00)
       D2ADR2=20.0d00*(R-RMAX)**2.0D00*(4.0D00*R-RMAX)/(RMAX**5.0D00)
    END IF
  END SUBROUTINE RBF_C2

!------------------------------------
!*** RADIAL AND DERIVATIVES
!*** XI IS THE "NODE" COORDINATE
!*** X IS THE GAUSS POINT COORDINATE
!*** chk0
!------------------------------------
  SUBROUTINE RBF_ONE(RMAX,NDI,X,XI,A,DADX,D2ADX2)
    IMPLICIT REAL(8) (A-H,O-Z)
    REAL(8)::RMAX,R,A,DADR,D2ADR2
    REAL(8),DIMENSION(NDI)::X,XI
    REAL(8),DIMENSION(NDI)::DADX,DRDX
    REAL(8),DIMENSION(NDI,NDI)::D2ADX2,D2RDX2
    CALL RBF_RADIUS(NDI,X,XI,R,DRDX,D2RDX2)
    CALL RBF_C2(RMAX,R,A,DADR,D2ADR2)
    DO ID=1,NDI
       DADX(ID)=DADR*DRDX(ID)
       DO JD=1,NDI
          D2ADX2(ID,JD)=D2ADR2*DRDX(ID)*DRDX(JD)+DADR*D2RDX2(ID,JD)
       END DO
    END DO
  END SUBROUTINE RBF_ONE
  
!------------------------------------------
!*** OBTAIN SHAPE FUNCTIONS AND DERIVATIVES
!*** FROM A GIVEN LIST OF NODES
!*** AT ONE POINT WITH COORDINATES X
!*** CHK0         
!------------------------------------------
  SUBROUTINE RBF_SHAPEFUNCTIONS(ISTYLE,RMAX,NDI,N,XN,X,FF,DFF,DFF2)
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER,PARAMETER::MPOL=10
    REAL(8)::RMAX
    REAL(8),DIMENSION(N,N)::AMATRIX,INVMATRIX,GMATRIX
    REAL(8),DIMENSION(N,MPOL)::BMATRIX
    REAL(8),DIMENSION(MPOL,N)::HMATRIX
    REAL(8),DIMENSION(NDI,NDI,N)::DFF2
    REAL(8),DIMENSION(NDI,N)::DFF
    REAL(8),DIMENSION(N)::FF,A
    REAL(8),DIMENSION(NDI,N)::XN
    REAL(8),DIMENSION(NDI)::X,XAV
    REAL(8),DIMENSION(MPOL)::B
    REAL(8),DIMENSION(NDI,MPOL)::DBDX
    REAL(8),DIMENSION(NDI,NDI,MPOL)::D2BDX2
    REAL(8),DIMENSION(MPOL,MPOL)::B1B,INVB1B
    REAL(8),DIMENSION(NDI,N)::DADX
    REAL(8),DIMENSION(NDI,NDI,N)::D2ADX2
    REAL(8),DIMENSION(MPOL)::POLYN
    REAL(8),DIMENSION(NDI,MPOL)::DPOLYN
!------------------------
!*** AVERAGE COORDINATES
!------------------------
    DO ID=1,NDI
       XAV(ID)=0.0D00
       DO IN=1,N
          XAV(ID)=XAV(ID)+XN(ID,IN)
       END DO
       XAV(ID)=XAV(ID)/(1.0D00*N)
    END DO
!------------------
!*** FORMS AMATRIX
!------------------
    DO JN=1,N
       DO IN=1,N
!-----------------------------------
!*** DADX AND D2ADX2 ARE TRASH here
!-----------------------------------
          CALL RBF_ONE(RMAX,NDI,XN(1:NDI,IN),XN(1:NDI,JN),AMATRIX(IN,JN),DADX(1:NDI,IN),D2ADX2(1:NDI,1:NDI,IN))
       END DO
    END DO
!-------------------------
!*** INVMATRIX=AMATRIX^-1
!-------------------------
    CALL INVMAT(N,DET,AMATRIX,INVMATRIX)
!------------------------    
!*** STANDARD POLYNOMIALS
!------------------------
!    CALL MLS_POLYNQBASE(NDI,1.0d00,X,XAV,POLYN,DPOLYN,MP)
!------------------
!*** FORMS BMATRIX
!------------------
    DO IN=1,N
       IF(ISTYLE.EQ.2)THEN
          CALL MLS_POLYNQBASE(NDI,1.0d00,XN(1:NDI,IN),XAV,POLYN,DPOLYN,MP)
       ELSE
          CALL MLS_POLYNLINBASE(NDI,1.0d00,XN(1:NDI,IN),XAV,POLYN,DPOLYN,MP)
       END IF
       DO IP=1,MP
          BMATRIX(IN,IP)=POLYN(IP)
       END DO
    END DO
!-----------------------
!*** FORMS BIG MATRICES
!-----------------------
    DO JD=1,MP
       DO ID=1,MP
          B1B(ID,JD)=0.0D00
          DO IN=1,N
             DO JN=1,N
                B1B(ID,JD)=B1B(ID,JD)+BMATRIX(IN,ID)*INVMATRIX(IN,JN)*BMATRIX(JN,JD)
             END DO
          END DO
       END DO
    END DO
!----------------
!*** INVERT B1B
!*** INTO INVB1B
!----------------
    CALL INVMAT(MP,DET,B1B(1:MP,1:MP),INVB1B(1:MP,1:MP))
!------------
!*** HMATRIX
!------------
    DO IN=1,N
       DO ID=1,MP
          HMATRIX(ID,IN)=0.0D00          
          DO JN=1,N
             DO JD=1,MP
                HMATRIX(ID,IN)=HMATRIX(ID,IN)+INVB1B(ID,JD)*BMATRIX(JN,JD)*INVMATRIX(JN,IN)
             END DO
          END DO
       END DO
    END DO
!------------
!*** GMATRIX
!------------
    DO JN=1,N
       DO IN=1,N
          GMATRIX(IN,JN)=INVMATRIX(IN,JN)
          DO ID=1,MP
             DO KN=1,N
                GMATRIX(IN,JN)=GMATRIX(IN,JN)-INVMATRIX(IN,KN)*BMATRIX(KN,ID)*HMATRIX(ID,JN)
             END DO
          END DO
       END DO
    END DO
!------------------
!*** FORMS VECTORS
!------------------
    DO IN=1,N
       CALL RBF_ONE(RMAX,NDI,X(1:NDI),XN(1:NDI,IN),A(IN),DADX(1:NDI,IN),D2ADX2(1:NDI,1:NDI,IN))
       IF(ISTYLE.EQ.2)THEN
          CALL MLS_POLYNQBASE(NDI,1.0d00,X(1:NDI),XAV,B,DBDX(1:NDI,1:MP),MP)
       ELSE
          CALL MLS_POLYNLINBASE(NDI,1.0d00,X(1:NDI),XAV,B,DBDX(1:NDI,1:MP),MP)
       END IF
    END DO
    D2BDX2=0.0D00
!------------------------------------------------------
!*** NOW DETERMINES RB SHAPE FUNCTIONS AND DERIVATIVES
!------------------------------------------------------
    IF(ISTYLE.EQ.0)THEN
       DO IN=1,N
          FF(IN)=DOTPROD(N,A(1:N),INVMATRIX(1:N,IN))
          DO ID=1,NDI
             DFF(ID,IN)=DOTPROD(N,DADX(ID,1:N),INVMATRIX(1:N,IN))
             DO JD=1,NDI
                DFF2(ID,JD,IN)=DOTPROD(N,D2ADX2(ID,JD,1:N),INVMATRIX(1:N,IN))
             ENDDO
          ENDDO
       ENDDO
    ELSE
       DO IN=1,N
          FF(IN)=DOTPROD(N,A(1:N),GMATRIX(1:N,IN))+DOTPROD(MP,B(1:MP),HMATRIX(1:MP,IN))
          DO ID=1,NDI
             DFF(ID,IN)=DOTPROD(N,DADX(ID,1:N),GMATRIX(1:N,IN))+DOTPROD(MP,DBDX(ID,1:MP),HMATRIX(1:MP,IN))
             DO JD=1,NDI
                DFF2(ID,JD,IN)=DOTPROD(N,D2ADX2(ID,JD,1:N),GMATRIX(1:N,IN))+DOTPROD(MP,D2BDX2(ID,JD,1:MP),HMATRIX(1:MP,IN))
             ENDDO
          ENDDO
       ENDDO
    END IF
  END SUBROUTINE RBF_SHAPEFUNCTIONS

!-------------------------------------------------
!*** FROM AN ELEMENT SUPPORT ESTIMATE, DETERMINES
!*** GAUSS POINT SUPPORT RADIUS 
!-------------------------------------------------
  SUBROUTINE RBF_GAUSSLEVELSHAPEFUNCTIONS(ISTYLE,RCENTROID,NDI,XG,XC,NTOT,XTOT,FF,DFF,DFF2)
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER::ISTYLE
    REAL(8)::RCENTROID,RGAUSS
    REAL(8),DIMENSION(NDI)::XG,XC
    REAL(8),DIMENSION(NTOT)::FF
    REAL(8),DIMENSION(NDI,NTOT)::DFF
    REAL(8),DIMENSION(NDI,NDI,NTOT)::DFF2
    REAL(8),DIMENSION(NDI,NTOT)::XTOT
    RGAUSS=RCENTROID-RNORM2(NDI,XG-XC)
    CALL RBF_TESTINGALLREQUIRED(ISTYLE,RGAUSS,NDI,NTOT,XTOT,XG,FF,DFF,DFF2)    
  END SUBROUTINE RBF_GAUSSLEVELSHAPEFUNCTIONS

!--------------------------------------------------------------------
!*** FOR REPRESENTATION PURPOSES ONLY:
!*** GET ALL NODES, MAKES A SELECTION
!*** AND RETURNS ALL QUANTITIES FOR ALL NODES, NOT ONLY THE SELECTION
!--------------------------------------------------------------------
  SUBROUTINE RBF_TESTINGALLREQUIRED(ISTYLE,RMAX,NDI,NTOT,XTOT,X,FF,DFF,DFF2)
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER::ISTYLE
    INTEGER,PARAMETER::MPOL=10
    INTEGER,DIMENSION(:),ALLOCATABLE::LISTN
    REAL(8)::RMAX
    REAL(8),DIMENSION(NTOT)::FF,FFLOC
    REAL(8),DIMENSION(NDI,NTOT)::DFF,DFFLOC
    REAL(8),DIMENSION(NDI,NDI,NTOT)::DFF2,DFF2LOC
    REAL(8),DIMENSION(NDI,NTOT)::XTOT,XN
    REAL(8),DIMENSION(NDI)::X,XMED
    REAL(8),DIMENSION(MPOL)::POLYN
    REAL(8),DIMENSION(NDI,MPOL)::DPOLYN
    REAL(8),DIMENSION(MPOL,NTOT)::U2
    !-----------------------
    !*** GETS CLOSEST NODES
    !-----------------------
    CALL MLS_GETSNODES(RMAX,X,NDI,NTOT,XTOT,N,LISTN)
    DO I=1,N
       DO ID=1,NDI
          XN(ID,I)=XTOT(ID,LISTN(I))
       END DO
    END DO
    !-----------------------------------
    !*** CHECK FOR REPEATED COORDINATES
    !-----------------------------------
    DO IN=1,N
       DO JN=IN+1,N
          IF(RNORM2(NDI,XN(1:NDI,IN)-XN(1:NDI,JN)).LE.1.0D-20)THEN
             WRITE(*,*)"NODE IN AND JN",IN,JN," ARE THE SAME"
             WRITE(*,*)"IN,X=",IN,XN(1:NDI,IN)
             WRITE(*,*)"JN,X=",JN,XN(1:NDI,JN)
          END IF
       END DO
    END DO
    !-------------------------------
    !*** SHAPE FUNCTION DERIVATIVES
    !-------------------------------
    IF(ISTYLE.GE.0)THEN
       CALL RBF_SHAPEFUNCTIONS(ISTYLE,RMAX,NDI,N,XN,X,FFLOC,DFFLOC,DFF2LOC)
    ELSE
       !--------   
       !**** MLS
       !--------
       TOL=1.0D-2    
       XMED=0.0D00
       DO I=1,N
          XMED=XMED+XN(1:NDI,I)
       END DO
       XMED=XMED/(1.0D00*N)
       CALL MLS_U2ATACOORDINATE(RMAX,TOL,NDI,N,XMED,XN,X,U2(1:MPOL,1:N))
       CALL MLS_POLYNQBASE(NDI,RMAX,X(1:NDI),XMED(1:NDI),POLYN(1:MPOL),DPOLYN(1:NDI,1:MPOL),M)
       CALL MLS_SFDER(NDI,N,M,POLYN(1:M),DPOLYN(1:NDI,1:M),U2(1:M,1:N),FFLOC,DFFLOC)
    END IF
    FF=0.0D00
    DFF=0.0D00
    DFF2=0.0D00
    DO I=1,N
       FF(LISTN(I))=FFLOC(I)
       DFF(1:NDI,LISTN(I))=DFFLOC(1:NDI,I)
       DFF2(1:NDI,1:NDI,LISTN(I))=DFF2LOC(1:NDI,1:NDI,I)
    END DO
    DEALLOCATE(LISTN)
  END SUBROUTINE RBF_TESTINGALLREQUIRED
!-------------------------
!*** GETS INFLUENCE NODES
!*** FOR A GIVEN POINT
!*** chk0
!-------------------------
  SUBROUTINE MLS_GETSNODES(RMAX,X,NDI,NTOT,XTOT,N,LISTN)
    IMPLICIT REAL(8)(A-H,O-Z)
    REAL(8)::RMAX
    REAL(8),DIMENSION(NDI,NTOT)::XTOT
    REAL(8),DIMENSION(NDI)::X
    INTEGER::NDI
    INTEGER,DIMENSION(:),ALLOCATABLE::LISTN
    N=0
    DO ITOT=1,NTOT
       IF(RNORM2(NDI,XTOT(1:NDI,ITOT)-X(1:NDI)).LE.RMAX)THEN
          N=N+1
       END IF
    END DO
    ALLOCATE(LISTN(N))
!*** NOW INSERTS THE NODES IN THE LIST
    N=0
    DO ITOT=1,NTOT
       IF(RNORM2(NDI,XTOT(1:NDI,ITOT)-X(1:NDI)).LE.RMAX)THEN
          N=N+1
          LISTN(N)=ITOT
       END IF
    END DO
  END SUBROUTINE MLS_GETSNODES

!-------------------------------------------------
!*** DETERMINES THE DEFORMATION GRADIENT
!*** AND GREEN-LAGRANGE STRAIN IN ENGINEERING FORM
!*** MLS STUFF
!*** CHK0,1
!-------------------------------------------------
  SUBROUTINE MLS_STRAIN(NDI,N,XTOT,XDEF,DFF,F,STRAIN)
    IMPLICIT REAL(8) (a-h,o-z)
    REAL(8),DIMENSION(NDI,N)::XTOT,XDEF
    REAL(8),DIMENSION(NDI,NDI)::C,F,FO,CO
    REAL(8),DIMENSION(NDI,N)::DFF
    REAL(8),DIMENSION(NDI*(NDI+1)/2)::STRAIN,STRAINO
!-------------------------
!*** DEFORMATION GRADIENT
!*** BEFORE CORRECTION
!-------------------------
    F=0.0D00
    DO IN=1,N
       DO JD=1,NDI
          DO ID=1,NDI
             F(ID,JD)=F(ID,JD)+DFF(JD,IN)*XDEF(ID,IN)
          END DO
       END DO
    END DO    
!------------------------------
!*** RIGHT CAUCHY-GREEN TENSOR
!------------------------------
    CALL MATMAT(NDI,NDI,NDI,F,F,C,2)
    DO ID=1,NDI*(NDI+1)/2
       CALL APOMAT(I1,I2,ID,NDI)
       STRAIN(ID)=0.5D00*(C(I1,I2)-DELTAK(I1,I2))
    END DO
    DO ID=NDI+1,NDI*(NDI+1)/2
       STRAIN(ID)=2.0D00*STRAIN(ID)
    END DO
!------------------------------------
!*** UNDEFORMED DEFORMATION GRADIENT
!------------------------------------
!!$    FO=0.0D00
!!$    DO IN=1,N
!!$       DO JD=1,NDI
!!$          DO ID=1,NDI
!!$             FO(ID,JD)=FO(ID,JD)+DFF(JD,IN)*XTOT(ID,IN)
!!$          END DO
!!$       END DO
!!$    END DO
!!$!------------------------------
!!$!*** RIGHT CAUCHY-GREEN TENSOR
!!$!------------------------------
!!$    CALL MATMAT(NDI,NDI,NDI,FO,FO,CO,2)
!!$    DO ID=1,NDI*(NDI+1)/2
!!$       CALL APOMAT(I1,I2,ID,NDI)
!!$       STRAINO(ID)=0.5D00*(CO(I1,I2)-DELTAK(I1,I2))
!!$    END DO
!!$    DO ID=NDI+1,NDI*(NDI+1)/2
!!$       STRAINO(ID)=2.0D00*STRAINO(ID)
!!$    END DO
!!$!--------------
!!$!*** NOW DO IT
!!$!--------------
!!$    STRAIN=STRAIN-STRAINO
!!$    CALL MATINV(NDI,DET,FO,FO,.FALSE.)
!!$    CALL MATMAT(NDI,NDI,NDI,F,FO,F,0)
!!$    WRITE(*,*)"f=",F
  END SUBROUTINE MLS_STRAIN

  SUBROUTINE MLS_2DFORCE(NK,F,S,NUCLEUS)
    IMPLICIT NONE
    DOUBLE PRECISION V(27),NK(2),F(2,2),S(3),NUCLEUS(2)
    V(22)=NK(2)*S(2)
    V(21)=NK(1)*S(1)
    NUCLEUS(1)=(F(1,2)*NK(1)+F(1,1)*NK(2))*S(3)+F(1,1)*V(21)+F(1,2)*V(22)
    NUCLEUS(2)=(F(2,2)*NK(1)+F(2,1)*NK(2))*S(3)+F(2,1)*V(21)+F(2,2)*V(22)
  END SUBROUTINE MLS_2DFORCE

  SUBROUTINE MLS_2D(NK,NL,F,S,DS,KERNEL)
    IMPLICIT NONE
    DOUBLE PRECISION V(67),NK(2),NL(2),F(2,2),S(3),DS(3,3),KERNEL(2,2)
    V(53)=NK(2)*(NL(2)*S(2)+NL(1)*S(3))+NK(1)*(NL(1)*S(1)+NL(2)*S(3))
    V(52)=F(2,2)*NK(1)+F(2,1)*NK(2)
    V(51)=F(1,2)*NK(1)+F(1,1)*NK(2)
    V(50)=F(2,2)*NK(2)
    V(49)=F(1,2)*NK(2)
    V(48)=F(2,1)*NK(1)
    V(47)=F(1,1)*NK(1)
    V(46)=F(2,2)*NL(1)+F(2,1)*NL(2)
    V(45)=F(1,2)*NL(1)+F(1,1)*NL(2)
    V(44)=F(2,2)*NL(2)
    V(60)=DS(2,2)*V(44)
    V(59)=DS(1,2)*V(44)+DS(1,3)*V(46)
    V(43)=F(1,2)*NL(2)
    V(54)=DS(1,2)*V(43)+DS(1,3)*V(45)
    V(42)=F(2,1)*NL(1)
    V(62)=DS(3,1)*V(42)+DS(3,2)*V(44)+DS(3,3)*V(46)
    V(61)=DS(2,1)*V(42)+DS(2,3)*V(46)
    V(58)=DS(1,1)*V(42)
    V(41)=F(1,1)*NL(1)
    V(56)=DS(3,1)*V(41)+DS(3,2)*V(43)+DS(3,3)*V(45)
    V(55)=DS(2,1)*V(41)+DS(2,3)*V(45)
    V(38)=V(47)*V(58)
    V(39)=V(49)*V(60)
    V(57)=V(38)+V(39)
    KERNEL(1,1)=V(53)+V(47)*(DS(1,1)*V(41)+V(54))+V(49)*(DS(2,2)*V(43)+V(55))+V(51)*V(56)
    KERNEL(1,2)=V(57)+V(47)*V(59)+V(49)*V(61)+V(51)*V(62)
    KERNEL(2,1)=V(48)*V(54)+V(50)*V(55)+V(52)*V(56)+V(57)
    KERNEL(2,2)=V(53)+V(48)*(V(58)+V(59))+V(50)*(V(60)+V(61))+V(52)*V(62)
  END SUBROUTINE MLS_2D

  SUBROUTINE MLS_3DFORCE(NK,F,S,NUCLEUS)
    IMPLICIT NONE
    DOUBLE PRECISION V(51),NK(3),F(3,3),S(6),NUCLEUS(3)
    V(46)=NK(3)*S(3)
    V(45)=NK(2)*S(2)
    V(44)=NK(1)*S(1)
    NUCLEUS(1)=(F(1,2)*NK(1)+F(1,1)*NK(2))*S(4)+(F(1,3)*NK(1)+F(1,1)*NK(3))*S(5)+(F(1,3)*NK(2)+F(1,2)*NK(3))*S(6)+F(1,1)*V&
         &(44)+F(1,2)*V(45)+F(1,3)*V(46)
    NUCLEUS(2)=(F(2,2)*NK(1)+F(2,1)*NK(2))*S(4)+(F(2,3)*NK(1)+F(2,1)*NK(3))*S(5)+(F(2,3)*NK(2)+F(2,2)*NK(3))*S(6)+F(2,1)*V&
         &(44)+F(2,2)*V(45)+F(2,3)*V(46)
    NUCLEUS(3)=(F(3,2)*NK(1)+F(3,1)*NK(2))*S(4)+(F(3,3)*NK(1)+F(3,1)*NK(3))*S(5)+(F(3,3)*NK(2)+F(3,2)*NK(3))*S(6)+F(3,1)*V&
         &(44)+F(3,2)*V(45)+F(3,3)*V(46)
  END SUBROUTINE MLS_3DFORCE
  
  SUBROUTINE MLS_3D(NK,NL,F,S,DS,KERNEL)
    IMPLICIT NONE
    DOUBLE PRECISION V(173),NK(3),NL(3),F(3,3),S(6),DS(6,6),KERNEL(3,3)
    V(144)=NK(1)*(NL(1)*S(1)+NL(2)*S(4)+NL(3)*S(5))+NK(3)*(NL(3)*S(3)+NL(1)*S(5)+NL(2)*S(6))+NK(2)*(NL(2)*S(2)+NL(1)*S(4)&
         &+NL(3)*S(6))
    V(143)=F(3,3)*NK(2)+F(3,2)*NK(3)
    V(142)=F(2,3)*NK(2)+F(2,2)*NK(3)
    V(141)=F(1,3)*NK(2)+F(1,2)*NK(3)
    V(140)=F(3,3)*NK(1)+F(3,1)*NK(3)
    V(139)=F(2,3)*NK(1)+F(2,1)*NK(3)
    V(138)=F(1,3)*NK(1)+F(1,1)*NK(3)
    V(137)=F(3,2)*NK(1)+F(3,1)*NK(2)
    V(136)=F(2,2)*NK(1)+F(2,1)*NK(2)
    V(135)=F(1,2)*NK(1)+F(1,1)*NK(2)
    V(134)=F(3,3)*NK(3)
    V(133)=F(2,3)*NK(3)
    V(132)=F(1,3)*NK(3)
    V(131)=F(3,2)*NK(2)
    V(130)=F(2,2)*NK(2)
    V(129)=F(1,2)*NK(2)
    V(128)=F(3,1)*NK(1)
    V(127)=F(2,1)*NK(1)
    V(126)=F(1,1)*NK(1)
    V(125)=F(3,3)*NL(2)+F(3,2)*NL(3)
    V(124)=F(2,3)*NL(2)+F(2,2)*NL(3)
    V(123)=F(1,3)*NL(2)+F(1,2)*NL(3)
    V(122)=F(3,3)*NL(1)+F(3,1)*NL(3)
    V(121)=F(2,3)*NL(1)+F(2,1)*NL(3)
    V(120)=F(1,3)*NL(1)+F(1,1)*NL(3)
    V(119)=F(3,2)*NL(1)+F(3,1)*NL(2)
    V(118)=F(2,2)*NL(1)+F(2,1)*NL(2)
    V(117)=F(1,2)*NL(1)+F(1,1)*NL(2)
    V(116)=F(3,3)*NL(3)
    V(115)=F(2,3)*NL(3)
    V(168)=DS(3,3)*V(115)
    V(114)=F(1,3)*NL(3)
    V(147)=DS(3,3)*V(114)
    V(113)=F(3,2)*NL(2)
    V(160)=DS(1,2)*V(113)+DS(1,3)*V(116)+DS(1,4)*V(119)+DS(1,5)*V(122)+DS(1,6)*V(125)
    V(112)=F(2,2)*NL(2)
    V(167)=DS(2,2)*V(112)
    V(154)=DS(1,2)*V(112)+DS(1,3)*V(115)+DS(1,4)*V(118)+DS(1,5)*V(121)+DS(1,6)*V(124)
    V(111)=F(1,2)*NL(2)
    V(148)=DS(1,2)*V(111)+DS(1,3)*V(114)+DS(1,4)*V(117)+DS(1,5)*V(120)+DS(1,6)*V(123)
    V(146)=DS(2,2)*V(111)
    V(110)=F(3,1)*NL(1)
    V(165)=DS(6,1)*V(110)+DS(6,2)*V(113)+DS(6,3)*V(116)+DS(6,4)*V(119)+DS(6,5)*V(122)+DS(6,6)*V(125)
    V(164)=DS(5,1)*V(110)+DS(5,2)*V(113)+DS(5,3)*V(116)+DS(5,4)*V(119)+DS(5,5)*V(122)+DS(5,6)*V(125)
    V(163)=DS(4,1)*V(110)+DS(4,2)*V(113)+DS(4,3)*V(116)+DS(4,4)*V(119)+DS(4,5)*V(122)+DS(4,6)*V(125)
    V(162)=DS(3,1)*V(110)+DS(3,2)*V(113)+DS(3,4)*V(119)+DS(3,5)*V(122)+DS(3,6)*V(125)
    V(161)=DS(2,1)*V(110)+DS(2,3)*V(116)+DS(2,4)*V(119)+DS(2,5)*V(122)+DS(2,6)*V(125)
    V(109)=F(2,1)*NL(1)
    V(166)=DS(1,1)*V(109)
    V(159)=DS(6,1)*V(109)+DS(6,2)*V(112)+DS(6,3)*V(115)+DS(6,4)*V(118)+DS(6,5)*V(121)+DS(6,6)*V(124)
    V(158)=DS(5,1)*V(109)+DS(5,2)*V(112)+DS(5,3)*V(115)+DS(5,4)*V(118)+DS(5,5)*V(121)+DS(5,6)*V(124)
    V(157)=DS(4,1)*V(109)+DS(4,2)*V(112)+DS(4,3)*V(115)+DS(4,4)*V(118)+DS(4,5)*V(121)+DS(4,6)*V(124)
    V(156)=DS(3,1)*V(109)+DS(3,2)*V(112)+DS(3,4)*V(118)+DS(3,5)*V(121)+DS(3,6)*V(124)
    V(155)=DS(2,1)*V(109)+DS(2,3)*V(115)+DS(2,4)*V(118)+DS(2,5)*V(121)+DS(2,6)*V(124)
    V(108)=F(1,1)*NL(1)
    V(153)=DS(6,1)*V(108)+DS(6,2)*V(111)+DS(6,3)*V(114)+DS(6,4)*V(117)+DS(6,5)*V(120)+DS(6,6)*V(123)
    V(152)=DS(5,1)*V(108)+DS(5,2)*V(111)+DS(5,3)*V(114)+DS(5,4)*V(117)+DS(5,5)*V(120)+DS(5,6)*V(123)
    V(151)=DS(4,1)*V(108)+DS(4,2)*V(111)+DS(4,3)*V(114)+DS(4,4)*V(117)+DS(4,5)*V(120)+DS(4,6)*V(123)
    V(150)=DS(3,1)*V(108)+DS(3,2)*V(111)+DS(3,4)*V(117)+DS(3,5)*V(120)+DS(3,6)*V(123)
    V(149)=DS(2,1)*V(108)+DS(2,3)*V(114)+DS(2,4)*V(117)+DS(2,5)*V(120)+DS(2,6)*V(123)
    V(145)=DS(1,1)*V(108)
    V(104)=V(127)*V(145)+V(130)*V(146)+V(133)*V(147)
    V(105)=V(128)*V(145)+V(131)*V(146)+V(134)*V(147)
    V(106)=V(128)*V(166)+V(131)*V(167)+V(134)*V(168)
    KERNEL(1,1)=V(144)+V(126)*(V(145)+V(148))+V(129)*(V(146)+V(149))+V(132)*(V(147)+V(150))+V(135)*V(151)+V(138)*V(152)+V&
         &(141)*V(153)
    KERNEL(1,2)=V(104)+V(126)*V(154)+V(129)*V(155)+V(132)*V(156)+V(135)*V(157)+V(138)*V(158)+V(141)*V(159)
    KERNEL(1,3)=V(105)+V(126)*V(160)+V(129)*V(161)+V(132)*V(162)+V(135)*V(163)+V(138)*V(164)+V(141)*V(165)
    KERNEL(2,1)=V(104)+V(127)*V(148)+V(130)*V(149)+V(133)*V(150)+V(136)*V(151)+V(139)*V(152)+V(142)*V(153)
    KERNEL(2,2)=V(144)+V(136)*V(157)+V(139)*V(158)+V(142)*V(159)+V(127)*(V(154)+V(166))+V(130)*(V(155)+V(167))+V(133)*(V&
         &(156)+V(168))
    KERNEL(2,3)=V(106)+V(127)*V(160)+V(130)*V(161)+V(133)*V(162)+V(136)*V(163)+V(139)*V(164)+V(142)*V(165)
    KERNEL(3,1)=V(105)+V(128)*V(148)+V(131)*V(149)+V(134)*V(150)+V(137)*V(151)+V(140)*V(152)+V(143)*V(153)
    KERNEL(3,2)=V(106)+V(128)*V(154)+V(131)*V(155)+V(134)*V(156)+V(137)*V(157)+V(140)*V(158)+V(143)*V(159)
    KERNEL(3,3)=V(144)+V(128)*(DS(1,1)*V(110)+V(160))+V(131)*(DS(2,2)*V(113)+V(161))+V(134)*(DS(3,3)*V(116)+V(162))+V(137&
         &)*V(163)+V(140)*V(164)+V(143)*V(165)
  END SUBROUTINE MLS_3D
