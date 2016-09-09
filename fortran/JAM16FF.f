
      !PROGRAM JAM16FF
      !IMPLICIT NONE
      !INTEGER ipos
      !CHARACTER(len=20)::flav,hadron
      !REAL*8 z,Q2,dp,get_zF

      !ipos=0
      !flav='dp'
      !hadron='pion'
      !z=0.5
      !Q2=10.0

      !print *, get_zF(ipos,hadron,flav,z,Q2)

      !END PROGRAM


***************************************************************************
*
*     JAM16 Fragmentation Functions
*     
*     Returns z*D_q+(z,Q2) where q+ = q + qbar
*
*     ipos (integer): posterior number (0 - 200)
*     flav (string): 'up','dp','sp','cp','bp','gl'
*     hadron (string): 'pion','kaon'
*  
*
***************************************************************************

      FUNCTION get_zF(ipos,hadron,flav,z,Q2)
      IMPLICIT NONE
      INTEGER nq2,nz,ipos
      PARAMETER (nq2=30,nz=100)
      INTEGER i,iz,iq2
      REAL*8 ZA(nz),Q2A(nq2),up(nq2,nz),dp(nq2,nz),sp(nq2,nz),cp(nq2,nz)
      REAL*8 bp(nq2,nz),gl(nq2,nz)
      REAL*8 z,Q2,zF,error,RINTERP2D,get_zF
      CHARACTER(len=20)::filename,folder,hadron,flav

      !ipos = 0
      !z=0.5
      !Q2=10.0
      !hadron='pion'
      !flav='up'

      if (hadron.eq.'pion') then
         folder='FFpion'
      elseif (hadron.eq.'kaon') then
         folder='FFkaon'
      endif

      if (ipos < 10) then
         WRITE(filename,'(A3,I1)') 'zF-',ipos
      elseif (ipos>10.and.ipos<100) then
         WRITE(filename,'(A3,I2)') 'zF-',ipos
      elseif (ipos>100) then
         WRITE(filename,'(A3,I3)') 'zF-',ipos
      endif

      OPEN(10,FILE=trim(folder)//'/'//trim(filename)//'.dat',
     &                                            STATUS='old')
      DO i=1,3,1
         READ(10,*)
      ENDDO
      READ(10,*) ZA
      READ(10,*) Q2A
      READ(10,*) 
      DO iq2=1,nq2
         DO iz=1,nz
            READ(10,*) dp(iq2,iz),up(iq2,iz),sp(iq2,iz),cp(iq2,iz),
     &                 bp(iq2,iz),gl(iq2,iz)
         ENDDO
      ENDDO
      CLOSE(10)


      if (flav.eq.'dp') then
         zF= RINTERP2D(Q2A,ZA,dp,Q2,z,nq2,nz)
!         Call interp2D(Q2A,ZA,dp,nQ2,nz,Q2,z,zF,error)
      elseif (flav.eq.'up') then
         zF = RINTERP2D(Q2A,ZA,up,Q2,z,nq2,nz)
!         Call interp2D(Q2A,ZA,up,nQ2,nz,Q2,z,zF,error)
      elseif (flav.eq.'sp') then
         zF = RINTERP2D(Q2A,ZA,sp,Q2,z,nq2,nz)
!         Call interp2D(Q2A,ZA,sp,nQ2,nz,Q2,z,zF,error)
      elseif (flav.eq.'cp') then
         zF = RINTERP2D(Q2A,ZA,cp,Q2,z,nq2,nz)
!         Call interp2D(Q2A,ZA,cp,nQ2,nz,Q2,z,zF,error)
      elseif (flav.eq.'bp') then
         zF = RINTERP2D(Q2A,ZA,bp,Q2,z,nq2,nz)
!         Call interp2D(Q2A,ZA,bp,nQ2,nz,Q2,z,zF,error)
      elseif (flav.eq.'gl') then
         zF = RINTERP2D(Q2A,ZA,gl,Q2,z,nq2,nz)
!         Call interp2D(Q2A,ZA,gl,nQ2,nz,Q2,z,zF,error)
      endif

      !print *,zF
      get_zF = zF

      RETURN
      END

**************************************************************************

      SUBROUTINE interp2D(x1a,x2a,ya,m,n,x1,x2,y,dy)
      INTEGER m,n,NMAX,MMAX
      REAL*8 dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
      PARAMETER (NMAX=100,MMAX=30)
      INTEGER j,k
      REAL*8 ymtmp(MMAX),yntmp(NMAX)
      DO j=1,m
         DO k=1,n
            yntmp(k)=ya(j,k)
         ENDDO
         CALL POLINT4(x2a,yntmp,x2,ymtmp(j))
!          CALL interp1D(x2a,yntmp,n,x2,ymtmp(j),dy)
      ENDDO
      CALL interp1D(x1a,ymtmp,m,x1,y,dy)
      RETURN
      END


**************************************************************************


      SUBROUTINE interp1D(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=100)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do i=1,n 
         dift=abs(x-xa(i))
         if (dift.lt.dif) then
            ns=i
            dif=dift
         endif
         c(i)=ya(i)
         d(i)=ya(i)
      enddo
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
         do i=1,n-m
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
            if(den.eq.0.) then 
               print *, 'failure in polint'
               exit
            endif
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
         enddo
         if (2*ns.lt.n-m)then
            dy=c(ns+1)
         else
            dy=d(ns)
            ns=ns-1
         endif
         y=y+dy
      enddo
      return
      END

**************************************************************************

        FUNCTION RINTERP2D (X,Y,F,X0,Y0,NX,NY)

        IMPLICIT NONE
        INTEGER NX,NY
        INTEGER I_EXTRAP,I,J,I1,I2,J1,J2
        REAL*8  X(NX),Y(NY),F(NX,NY)
        REAL*8  X0,Y0
        REAL*8  RINTX1,RINTX2,RINTERP2D

        I_EXTRAP = 0

        IF ((X0.GT.X(NX)).OR.(X0.LT.X(1))
     &       .OR.(Y0.GT.Y(NY)).OR.(Y0.LT.Y(1))) THEN
          I_EXTRAP = I_EXTRAP + 1
          RINTERP2D = 0.D0
          RETURN
        ENDIF

        DO I=2,NX
          DO J=2,NY
            IF ((X0.LE.X(I)).AND.(Y0.LE.Y(J))) THEN
              GOTO 1
            ENDIF
          ENDDO
        ENDDO

 1      I1 = I - 1
        J1 = J - 1
        RINTX1 = ((X0-X(I1))/(X(I1+1)-X(I1)))*(F(I1+1,J1)-F(I1,J1))
     &         + F(I1,J1)
        RINTX2 = ((X0-X(I1))/(X(I1+1)-X(I1)))*(F(I1+1,J1+1)-F(I1,J1+1))
     &         + F(I1,J1+1)

        RINTERP2D = ((Y0-Y(J1))/(Y(J1+1)-Y(J1)))*(RINTX2-RINTX1)+RINTX1

        RETURN
        END
***************************************************************************

      SUBROUTINE POLINT4 (XA,YA,X,Y)
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  The POLINT4 routine is based on the POLINT routine from "Numerical Recipes",
C  but assuming N=4, and ignoring the error estimation.
C  suggested by Z. Sullivan. 
      DIMENSION XA(*),YA(*)

      H1=XA(1)-X
      H2=XA(2)-X
      H3=XA(3)-X
      H4=XA(4)-X

      W=YA(2)-YA(1)
      DEN=W/(H1-H2)
      D1=H2*DEN
      C1=H1*DEN
      
      W=YA(3)-YA(2)
      DEN=W/(H2-H3)
      D2=H3*DEN
      C2=H2*DEN

      W=YA(4)-YA(3)
      DEN=W/(H3-H4)
      D3=H4*DEN
      C3=H3*DEN

      W=C2-D1
      DEN=W/(H1-H3)
      CD1=H3*DEN
      CC1=H1*DEN

      W=C3-D2
      DEN=W/(H2-H4)
      CD2=H4*DEN
      CC2=H2*DEN

      W=CC2-CD1
      DEN=W/(H1-H4)
      DD1=H4*DEN
      DC1=H1*DEN

      If((H3+H4).lt.0D0) Then
         Y=YA(4)+D3+CD2+DD1
      Elseif((H2+H3).lt.0D0) Then
         Y=YA(3)+D2+CD1+DC1
      Elseif((H1+H2).lt.0D0) Then
         Y=YA(2)+C2+CD1+DC1
      ELSE
         Y=YA(1)+C1+CC1+DC1
      ENDIF
      RETURN
      END
