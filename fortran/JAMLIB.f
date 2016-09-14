************************************************************************
*
*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
*          _   _    __  __ _     ___ ____                        
*         | | / \  |  \/  | |   |_ _| __ )                       
*         | |/ _ \ | |\/| | |    | ||  _ \                       
*     | |_| / ___ \| |  | | |___ | || |_) |                      
*      \___/_/   \_\_|  |_|_____|___|____/                       
*                                                           
*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
*                                                           
*     Authors:                                                   
*     Nobuo Sato         (Jefferson Lab)                         
*     Jake Ethier        (College of William and Mary)           
*     Wally Melnitchouk  (Jefferson Lab)                         
*     Alberto Accardi    (Hampton University and Jefferson Lab)
*
*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
*
*     FUNCTION get_xF(x,Q2,flav)
*
*     INPUT: x,Q^2, and flavor
*
*     flav (character*2): 'up','dp','sp','cp','bp','gl','u','d','s'
*                         'c','b','ub','db','sb','cb','bb','p','n'
*
*     Flavors 'qp' = q + qb (i.e. 'up' = u + ubar)
*             'p'  = proton (for T4PPDF and structure functions)
*             'n'  = neutron (for T4PPDF and structure functions)
*
*     Returns x times the function (PPDF,FF,etc.)
*
*     GRID_INIT must be called only once before using FUNCTION get_xF!
*
*     SUBROUTINE GRID_INIT requires three inputs:
*
*     - lib (character*10): library (JAM15,JAM16,etc.)
*     - dist (character*10): distribution type (PPDF,FFpion,FFkaon,etc.)
*     - ipos (integer): posterior number (0 to 199) from MC analysis
*
************************************************************************

      FUNCTION get_xF(x,Q2,flav)
      IMPLICIT NONE
      INTEGER nx,nq2,nc,ipos,INIT,ier,lwrk,kwrk
      PARAMETER (nx=104,nq2=34,nc=3000,lwrk=8,kwrk=2)
      INTEGER iwrk(kwrk)
      REAL*8 TX(nx),TQ2(nq2),CUP(nc),CDP(nc),CSP(nc),CCP(nc),CBP(nc)
      REAL*8 CG(nc),CU(nc),CD(nc),CS(nc),CC(nc),CB(nc),CUB(nc),CDB(nc)
      REAL*8 CSB(nc),CCB(nc),CBB(nc),CP(nc),CN(nc),CQ(nc)
      REAL*8 x,Q2,xarr(1),qarr(1),xF(1,1),get_xF
      REAL*8 wrk(lwrk)
      CHARACTER flav*2
      COMMON/GRID_INI/INIT
      COMMON/BSPLINE_COEFFS/CUP,CDP,CSP,CCP,CBP,CG,CU,CD,CS,CC,CB,CUB,
     &                      CDB,CSB,CCB,CBB,CP,CN
      COMMON/X_Q2_KNOTS/TX,TQ2

      IF (INIT.eq.0) then
         print *,'Grid was not initialized! Must call GRID_INIT first.'
         stop
      ENDIF

      IF (flav.eq.'up') THEN
         CQ = CUP
      ELSEIF (flav.eq.'dp') THEN
         CQ = CDP
      ELSEIF (flav.eq.'sp') THEN
         CQ = CSP
      ELSEIF (flav.eq.'cp') THEN
         CQ = CCP
      ELSEIF (flav.eq.'bp') THEN
         CQ = CBP
      ELSEIF (flav.eq.'gl') THEN
         CQ = CG
      ELSEIF (flav.eq.'u') THEN
         CQ = CU
      ELSEIF (flav.eq.'d') THEN
         CQ = CD
      ELSEIF (flav.eq.'s') THEN
         CQ = CU
      ELSEIF (flav.eq.'c') THEN
         CQ = CU
      ELSEIF (flav.eq.'b') THEN
         CQ = CU
      ELSEIF (flav.eq.'ub') THEN
         CQ = CU
      ELSEIF (flav.eq.'db') THEN
         CQ = CU
      ELSEIF (flav.eq.'sb') THEN
         CQ = CU
      ELSEIF (flav.eq.'cb') THEN
         CQ = CU
      ELSEIF (flav.eq.'bb') THEN
         CQ = CU
      ELSEIF (flav.eq.'p') THEN
         CQ = CP
      ELSEIF (flav.eq.'n') THEN
         CQ = CN
      ENDIF

      xarr(1)=x
      qarr(1)=Q2

      CALL bispev(TX,nx,TQ2,nq2,CQ,3,3,xarr,1,qarr,1,xF,wrk,lwrk,
     &            iwrk,kwrk,ier)

      get_xF=xF(1,1)

      RETURN
      END

************************************************************************

      SUBROUTINE GRID_INIT(lib,dist,ipos)
      IMPLICIT NONE
      INTEGER ipos,INIT
      INTEGER nq2,nx,nc
      PARAMETER (nx=104,nq2=34,nc=3000)
      CHARACTER*10 lib,dist,filename
      REAL*8 TX(nx),TQ2(nq2),CUP(nc),CDP(nc),CSP(nc),CCP(nc),CBP(nc)
      REAL*8 CG(nc),CU(nc),CD(nc),CS(nc),CC(nc),CB(nc),CUB(nc),CDB(nc)
      REAL*8 CSB(nc),CCB(nc),CBB(nc),CP(nc),CN(nc)
      LOGICAL PLUS
      COMMON/GRID_INI/INIT
      COMMON/BSPLINE_COEFFS/CUP,CDP,CSP,CCP,CBP,CG,CU,CD,CS,CC,CB,CUB,
     &                      CDB,CSB,CCB,CBB,CP,CN
      COMMON/X_Q2_KNOTS/TX,TQ2

      !! Need to update for future JAM libraries
      IF (lib.eq.'JAM15'.and.dist.eq.'PPDF') THEN
         PLUS = .TRUE.
      ELSEIF (lib.eq.'JAM16'.and.dist(1:2).eq.'FF') THEN
         PLUS = .TRUE.
      ELSE
         PLUS = .FALSE.
      ENDIF

      if (ipos < 10) then
         WRITE(filename,'(A3,I1)') 'xF-',ipos
      elseif (ipos>10.and.ipos<100) then
         WRITE(filename,'(A3,I2)') 'xF-',ipos
      elseif (ipos>100) then
         WRITE(filename,'(A3,I3)') 'xF-',ipos
      endif

      OPEN(10,FILE=trim(lib)//'/'//trim(dist)//'/'//
     &               trim(filename)//'.tab',STATUS='old')

      READ(10,*) TX
      READ(10,*) TQ2

      !! Need to update for future JAM libraries
      IF (dist.eq.'PPDF'.or.dist(1:2).eq.'FF') THEN
         IF (PLUS.eqv..TRUE.) THEN
            READ(10,*) CUP
            READ(10,*) CDP
            READ(10,*) CSP
            READ(10,*) CCP
            READ(10,*) CBP
            READ(10,*) CG
         ELSE
            READ(10,*) CSB
            READ(10,*) CDB
            READ(10,*) CUB
            READ(10,*) CU
            READ(10,*) CD
            READ(10,*) CS
            READ(10,*) CC
            READ(10,*) CB
            READ(10,*) CG
         ENDIF
      ELSEIF (dist.eq.'T3PPDF') THEN
         READ(10,*) CU
         READ(10,*) CD
      ELSEIF (dist.eq.'T4PPDF') THEN
         READ(10,*) CP
         READ(10,*) CN
      ENDIF
      CLOSE(10)
      INIT = 1
      RETURN
      END

************************************************************************
*
*     BIVARIATE SPLINE ROUTINE FROM FITPACK
*     http://www.netlib.org/dierckx/
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
