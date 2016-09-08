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
      REAL*8 z,Q2,zF,error,get_zF
      CHARACTER(len=20)::filename,folder,hadron,flav

      if (hadron.eq.'pion') then
         folder='FFpion'
      elseif (hadron.eq.'kaon') then
         folder='FFkaon'
      endif

      ipos=0
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
         CALL interp2D(Q2A,ZA,dp,nq2,nz,Q2,z,zF,error)
      elseif (flav.eq.'up') then
         CALL interp2D(Q2A,ZA,up,nq2,nz,Q2,z,zF,error)
      elseif (flav.eq.'sp') then
         CALL interp2D(Q2A,ZA,sp,nq2,nz,Q2,z,zF,error)
      elseif (flav.eq.'cp') then
         CALL interp2D(Q2A,ZA,cp,nq2,nz,Q2,z,zF,error)
      elseif (flav.eq.'bp') then
         CALL interp2D(Q2A,ZA,bp,nq2,nz,Q2,z,zF,error)
      elseif (flav.eq.'gl') then
         CALL interp2D(Q2A,ZA,gl,nq2,nz,Q2,z,zF,error)      
      endif

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
         CALL interp1D(x2a,yntmp,n,x2,ymtmp(j),dy)
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
