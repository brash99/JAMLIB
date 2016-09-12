
c      PROGRAM JAM16FF
c      IMPLICIT NONE
      
c      INTEGER ipos,iz
c      CHARACTER(len=20)::flav,hadron
c      REAL*8 z,Q2,dp,get_zF,TEST,up

c      ipos=0
c      flav='up'
c      hadron='pion'
c      z=0.05D0
c      Q2=1.D0

c      CALL GRID_INIT(hadron,ipos)

c      up = get_zF(flav,z,Q2)

c      DO iz=5,100
c         z = FLOAT(iz)/100.D0
c         dp = get_zF(flav,z,Q2)
c         print *,z,Q2,dp
c      ENDDO

c      END PROGRAM


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
      INTEGER nq2,nz,ipos,INIT,ier,iwrk(2),lwrk,kwrk
      PARAMETER (nq2=30,nz=100,lwrk=8,kwrk=2)
      INTEGER i,iz,iq2,J,J1,J2
      REAL*8 ZA(nz),QA(nq2),up(nz,nq2),dp(nz,nq2),sp(nz,nq2),cp(nz,nq2)
      REAL*8 bp(nz,nq2),gl(nz,nq2),zarr(1),qarr(1),zF(1,1)
      REAL*8 z,Q2,error,RINTERP2D,A1,A2,S1,S2,logQ,ANS,get_zF,get_fz
      REAL*8 S,SA,lam,Qini,Q,DQP(nz,nq2),Q2A(nq2),CUP(3000),TZ(104)
      REAL*8 TQ2(34),wrk(lwrk)
      CHARACTER(len=20)::flav
      COMMON/FF_INIT/INIT
      COMMON/FF_GRIDS/up,dp,sp,cp,bp,gl,CUP
      COMMON/Z_Q2_GRIDS/ZA,QA,TZ,TQ2

      IF (INIT.eq.0) then
         print *,'Grid was not initialized! Must call GRID_INIT first.'
         stop
      ENDIF

c      Qini=1.D0
c      q = DSQRT(Q2)
c      lam=0.2268D0
c      S = DLOG(DLOG(q/lam)/DLOG(Qini/lam))

c      DO 2 J=1,nq2
c      SA=DLOG(DLOG(QA(J)/lam)/DLOG(Qini/lam))
c      IF(S.LT.SA)THEN
c         J2=J
c         IF(J2.EQ.1)J2=2
c         J1=J2-1
c         S2=DLOG(DLOG(QA(J2)/lam)/DLOG(Qini/lam))
c         S1=DLOG(DLOG(QA(J1)/lam)/DLOG(Qini/lam))
c         GOTO 1
c      ENDIF
c 2    CONTINUE
c 1    CONTINUE

c      A1=get_fz(flav,z,J1)
c      A2=get_fz(flav,z,J2)
c      ANS=A1*(S-S2)/(S1-S2)+A2*(S-S1)/(S2-S1)

      zarr(1)=z
      qarr(1)=Q2

      CALL bispev(TZ,104,TQ2,34,CUP,3,3,zarr,1,qarr,1,zF,wrk,lwrk,
     &            iwrk,kwrk,ier)

c      print *, TZ,TQ2!zarr(1),qarr(1),zF(1,1),ier

      get_zF=zF(1,1)
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
c      COMMON/FF_GRIDS/up,dp,sp,cp,bp,gl
c      COMMON/Z_Q2_GRIDS/ZA,QA
      

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
      REAL*8 bp(nz,nq2),gl(nz,nq2),CUP(3000),TZ(104),TQ2(34)
      COMMON/FF_INIT/INIT
      COMMON/FF_GRIDS/up,dp,sp,cp,bp,gl,CUP
      COMMON/Z_Q2_GRIDS/ZA,QA,TZ,TQ2


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

      OPEN(10,FILE=trim(folder)//'/'//trim(filename)//'.tab',
     &                                            STATUS='old')
      READ(10,*) TZ
      READ(10,*) TQ2
      READ(10,*) CUP
!      DO iz=1,nz
!         DO iq2=1,nq2
!            READ(10,*) dp(iz,iq2),up(iz,iq2),sp(iz,iq2),cp(iz,iq2),
!     &                 bp(iz,iq2),gl(iz,iq2)
!         ENDDO
!      ENDDO
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

***********************************************************************
      subroutine bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
     * iwrk,kwrk,ier)
      integer nx,ny,kx,ky,mx,my,lwrk,kwrk,ier
      integer iwrk(kwrk)
      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),
     * wrk(lwrk)
      integer i,iw,lwest
      ier = 10
      lwest = (kx+1)*mx+(ky+1)*my
      if(lwrk.lt.lwest) go to 100
      if(kwrk.lt.(mx+my)) go to 100
      if(mx-1) 100,30,10
 10     do 20 i=2,mx
        if(x(i).lt.x(i-1)) go to 100
 20       continue
 30         if(my-1) 100,60,40
 40           do 50 i=2,my
        if(y(i).lt.y(i-1)) go to 100
 50       continue
 60         ier = 0
      iw = mx*(kx+1)+1
      call fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk(1),wrk(iw),
     * iwrk(1),iwrk(mx+1))
 100    return
      end
*************************************************************************

      subroutine fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wx,wy,lx,ly)
      integer nx,ny,kx,ky,mx,my
      integer lx(mx),ly(my)
      real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),
     * wx(mx,kx+1),wy(my,ky+1)
      integer kx1,ky1,l,l1,l2,m,nkx1,nky1
      real*8 arg,sp,tb,te
      real*8 h(6)
      kx1 = kx+1
      nkx1 = nx-kx1
      tb = tx(kx1)
      te = tx(nkx1+1)
      l = kx1
      l1 = l+1
      do 40 i=1,mx
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
 10         if(arg.lt.tx(l1) .or. l.eq.nkx1) go to 20
        l = l1
        l1 = l+1
        go to 10
 20         call fpbspl(tx,nx,kx,arg,l,h)
        lx(i) = l-kx1
        do 30 j=1,kx1
          wx(i,j) = h(j)
 30           continue
 40             continue
      ky1 = ky+1
      nky1 = ny-ky1
      tb = ty(ky1)
      te = ty(nky1+1)
      l = ky1
      l1 = l+1
      do 80 i=1,my
        arg = y(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
 50         if(arg.lt.ty(l1) .or. l.eq.nky1) go to 60
        l = l1
        l1 = l+1
        go to 50
 60         call fpbspl(ty,ny,ky,arg,l,h)
        ly(i) = l-ky1
        do 70 j=1,ky1
          wy(i,j) = h(j)
 70           continue
 80             continue
      m = 0
      do 130 i=1,mx
        l = lx(i)*nky1
        do 90 i1=1,kx1
          h(i1) = wx(i,i1)
 90           continue
        do 120 j=1,my
          l1 = l+ly(j)
          sp = 0.
          do 110 i1=1,kx1
            l2 = l1
            do 100 j1=1,ky1
              l2 = l2+1
              sp = sp+c(l2)*h(i1)*wy(j,j1)
 100                  continue
            l1 = l1+nky1
 110              continue
          m = m+1
          z(m) = sp
 120          continue
 130            continue
      return
      end

***********************************************************

      subroutine fpbspl(t,n,k,x,l,h)
      real*8 x
      integer n,k,l
      real*8 t(n),h(6)
      real*8 f,one
      integer i,j,li,lj
      real*8 hh(5)
      one = 0.1e+01
      h(1) = one
      do 20 j=1,k
        do 10 i=1,j
          hh(i) = h(i)
 10           continue
        h(1) = 0.
        do 20 i=1,j
          li = l+i
          lj = li-j
          f = hh(i)/(t(li)-t(lj))
          h(i) = h(i)+f*(t(li)-x)
          h(i+1) = f*(x-t(lj))
 20         continue
      return
      end
