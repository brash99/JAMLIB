************************************************************************
*
*     JAM16 Pion and Kaon Fragmentation Functions
*
*     Authors: N. Sato, J.J. Ethier, W. Melnitchouk, M. Hirai,
*              S. Kumano, A. Accardi
*     arXiv:1609.00899
*     
*     Returns z*D_q+ and z*D_g where q+ = q + qbar
*
*     flav (character*2): 'up','dp','sp','cp','bp','gl'
*  
*     GRID_INIT must be called only once before using FUNCTION get_zF!
*
*     SUBROUTINE GRID_INIT requires two inputs:
*
*       - ipos (integer): posterior number (0 to 199) from MC analysis
*       - hadron (character*10): 'pion','kaon'
*
************************************************************************

      FUNCTION get_zF(flav,z,Q2)
      IMPLICIT NONE
      INTEGER nz,nq2,nc,ipos,INIT,ier,lwrk,kwrk
      PARAMETER (nz=104,nq2=34,nc=3000,lwrk=8,kwrk=2)
      INTEGER iwrk(kwrk)
      REAL*8 TZ(nz),TQ2(nq2),CUP(nc),CDP(nc),CSP(nc),CCP(nc)
      REAL*8 CBP(nc),CG(nc),CQP(nc)
      REAL*8 z,Q2,zarr(1),qarr(1),zF(1,1),get_zF
      REAL*8 wrk(lwrk)
      CHARACTER flav*2
      COMMON/FF_INIT/INIT
      COMMON/FF_COEFFS/CUP,CDP,CSP,CCP,CBP,CG
      COMMON/Z_Q2_KNOTS/TZ,TQ2

      IF (INIT.eq.0) then
         print *,'Grid was not initialized! Must call GRID_INIT first.'
         stop
      ENDIF

      zarr(1)=z
      qarr(1)=Q2

      IF (flav.eq.'up') THEN
         CQP = CUP
      ELSEIF (flav.eq.'dp') THEN
         CQP = CDP
      ELSEIF (flav.eq.'sp') THEN
         CQP = CSP
      ELSEIF (flav.eq.'cp') THEN
         CQP = CCP
      ELSEIF (flav.eq.'bp') THEN
         CQP = CBP
      ELSEIF (flav.eq.'gl') THEN
         CQP = CG
      ENDIF

      CALL bispev(TZ,nz,TQ2,nq2,CQP,3,3,zarr,1,qarr,1,zF,wrk,lwrk,
     &            iwrk,kwrk,ier)

      get_zF=zF(1,1)

      RETURN
      END

************************************************************************

      SUBROUTINE GRID_INIT(hadron,ipos)
      IMPLICIT NONE
      INTEGER ipos,INIT
      INTEGER nq2,nz,nc
      PARAMETER (nz=104,nq2=34,nc=3000)
      CHARACTER*10 hadron,filename,folder
      REAL*8 TZ(nz),TQ2(nq2),CUP(nc),CDP(nc),CSP(nc),CCP(nc)
      REAL*8 CBP(nc),CG(nc)
      COMMON/FF_INIT/INIT
      COMMON/FF_COEFFS/CUP,CDP,CSP,CCP,CBP,CG
      COMMON/Z_Q2_KNOTS/TZ,TQ2


      if (hadron.eq.'pion') then
         folder='FFpion'
      elseif (hadron.eq.'kaon') then
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
      READ(10,*) CDP
      READ(10,*) CSP
      READ(10,*) CCP
      READ(10,*) CBP
      READ(10,*) CG
      CLOSE(10)
      INIT = 1
      RETURN
      END

************************************************************************
*
*     BIVARIATE SPLINE ROUTINE FROM FITPACK
*
************************************************************************
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

************************************************************************

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

************************************************************************

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

************************************************************************
