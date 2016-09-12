
      PROGRAM JAM16FF
      IMPLICIT NONE
      
      INTEGER ipos,iz
      CHARACTER(len=20)::flav,hadron
      REAL*8 z,Q2,dp,get_zF,TEST

      ipos=0
      flav='dp'
      hadron='pion'
      !z=0.5D0
      Q2=1.D0

      CALL GRID_INIT(hadron,ipos)

      DO iz=5,100
         z = FLOAT(iz)/100.D0
         dp = get_zF(flav,z,Q2)
         print *,z,Q2,dp
      ENDDO

      END PROGRAM


***************************************************************************
*
*     JAM16 Fragmentation Functions
*     
*     Returns z*D_q+(z,Q2) where q+ = q + qbar
*
*     flav (string): 'up','dp','sp','cp','bp','gl'
*  
*     GRID_INIT must be called only once before using FUNCTION get_zF!
*
*     SUBROUTINE GRID_INIT requires two inputs:
*
*       - ipos (integer): posterior number (0 - 200)
*         ipos=0 is central value
*         ipos=1-200 are replicas from MC analysis
*
*       - hadron (string): 'pion','kaon'
*
***************************************************************************

      FUNCTION get_zF(flav,z,Q2)
      IMPLICIT NONE
      INTEGER nq2,nz,ipos,INIT
      PARAMETER (nq2=30,nz=100)
      INTEGER i,iz,iq2,J,J1,J2
      REAL*8 ZA(nz),QA(nq2),up(nz,nq2),dp(nz,nq2),sp(nz,nq2),cp(nz,nq2)
      REAL*8 bp(nz,nq2),gl(nz,nq2)
      REAL*8 z,Q2,zF,error,RINTERP2D,A1,A2,S1,S2,logQ,ANS,get_zF,get_fz
      REAL*8 S,SA,lam,Qini,Q,DQP(nz,nq2),Q2A(nq2)
      CHARACTER(len=20)::flav
      COMMON/FF_INIT/INIT
      COMMON/FF_GRIDS/up,dp,sp,cp,bp,gl
      COMMON/Z_Q2_GRIDS/ZA,QA

      IF (INIT.eq.0) then
         print *,'Grid was not initialized! Must call GRID_INIT first.'
         stop
      ENDIF

      Qini=1.D0
      q = DSQRT(Q2)
      lam=0.2268D0
      S = DLOG(DLOG(q/lam)/DLOG(Qini/lam))

      DO 2 J=1,nq2
      SA=DLOG(DLOG(QA(J)/lam)/DLOG(Qini/lam))
      IF(S.LT.SA)THEN
         J2=J
         IF(J2.EQ.1)J2=2
         J1=J2-1
         S2=DLOG(DLOG(QA(J2)/lam)/DLOG(Qini/lam))
         S1=DLOG(DLOG(QA(J1)/lam)/DLOG(Qini/lam))
         GOTO 1
      ENDIF
 2    CONTINUE
 1    CONTINUE

      A1=get_fz(flav,z,J1)
      A2=get_fz(flav,z,J2)
      ANS=A1*(S-S2)/(S1-S2)+A2*(S-S1)/(S2-S1)

      get_zF=ANS
      RETURN
      END

**************************************************************************

      FUNCTION get_fz(flav,z,J)
      IMPLICIT NONE
      INTEGER nz,nq2,J,In1,I
      PARAMETER (nz = 100, nq2 = 30)
      REAL*8 zz(4),fz(4),up(nz,nq2),dp(nz,nq2),cp(nz,nq2),sp(nz,nq2)
      REAL*8 bp(nz,nq2),gl(nz,nq2),ZA(nz),QA(nq2),DQP(nz,nq2)
      REAL*8 get_fz,z,ans,error
      CHARACTER(len=20)::flav
      COMMON/FF_GRIDS/up,dp,sp,cp,bp,gl
      COMMON/Z_Q2_GRIDS/ZA,QA
      

      IF (flav.eq.'up')then
         DQP=up
      ELSEIF(flav.eq.'dp')then
         DQP=dp
      ELSEIF(flav.eq.'sp')then
         DQP=sp
      ELSEIF(flav.eq.'cp')then
         DQP=cp
      ELSEIF(flav.eq.'bp')then
         DQP=bp
      ELSEIF(flav.eq.'gl')then
         DQP=gl
      ENDIF

      DO 1 I=1,nz 
      IF(z.LT.ZA(I))GOTO 2
 1    CONTINUE
 2    I=I-2
      If(I.LE.0.d0)I=1
      If(I.GT.(nz-3))I=nz-3
      zz(1)=ZA(I)
      zz(2)=ZA(I+1)
      zz(3)=ZA(I+2)
      zz(4)=ZA(I+3)
      fz(1)=DQP(I,J)*z
      fz(2)=DQP(I+1,J)*z
      fz(3)=DQP(I+2,J)*z
      fz(4)=DQP(I+3,J)*z
      call POLINT4(zz,fz,z,ans)
      get_fz=ans/z

      RETURN
      END


**************************************************************************

      SUBROUTINE GRID_INIT(hadron,ipos)
      IMPLICIT NONE
      INTEGER ipos,INIT
      INTEGER nq2,nz
      PARAMETER (nq2=30,nz=100)
      INTEGER i,iz,iq2
      CHARACTER(len=20)::hadron,filename,folder
      REAL*8 ZA(nz),QA(nq2),up(nz,nq2),dp(nz,nq2),sp(nz,nq2),cp(nz,nq2)
      REAL*8 bp(nz,nq2),gl(nz,nq2)
      COMMON/FF_INIT/INIT
      COMMON/FF_GRIDS/up,dp,sp,cp,bp,gl
      COMMON/Z_Q2_GRIDS/ZA,QA


      if (hadron.eq.'pion') then
!         folder='JAM16FF_pi'
         folder='FFpion'
      elseif (hadron.eq.'kaon') then
!         folder='JAM16FF_K'
         folder='FFkaon'
      endif

      if (ipos < 10) then
         WRITE(filename,'(A3,I1)') 'xF-',ipos
      elseif (ipos>10.and.ipos<100) then
         WRITE(filename,'(A3,I2)') 'xF-',ipos
      elseif (ipos>100) then
         WRITE(filename,'(A3,I3)') 'xF-',ipos
      endif

      print *,filename

      OPEN(10,FILE=trim(folder)//'/'//trim(filename)//'.tab',
     &                                            STATUS='old')
      READ(10,*) ZA
      READ(10,*) QA
      DO iz=1,nz
         DO iq2=1,nq2
            READ(10,*) dp(iz,iq2),up(iz,iq2),sp(iz,iq2),cp(iz,iq2),
     &                 bp(iz,iq2),gl(iz,iq2)
         ENDDO
      ENDDO
      CLOSE(10)
      INIT = 1
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

        FUNCTION RINTERP2D(X,Y,F,X0,Y0,NX,NY)

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


      SUBROUTINE ratint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n),TINY
      PARAMETER (NMAX=4,TINY=1.e-25)
      INTEGER i,m,ns
      REAL*8 dd,h,hh,t,w,c(NMAX),d(NMAX) 
      ns=1
      hh=abs(x-xa(1))
      do i=1,n 
         h=abs(x-xa(i)) 
         if (h.eq.0.)then
            y=ya(i)
            dy=0.0
            return
         else if (h.lt.hh) then 
            ns=i
            hh=h 
         endif
         c(i)=ya(i)
         d(i)=ya(i)+TINY 
      enddo
      y=ya(ns) 
      ns=ns-1
      do m=1,n-1
         do i=1,n-m
            w=c(i+1)-d(i)
            h=xa(i+m)-x 
            t=(xa(i)-x)*d(i)/h
            dd=t-c(i+1)
            if(dd.eq.0.)then
               print *,'failure in ratint'
               stop
            endif
            dd=w/dd
            d(i)=c(i+1)*dd
            c(i)=t*dd 
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
