*****************************************************************
*
*      2D Interpolation 
*      Taken from Numerical Recipes in Fortran77
*
*****************************************************************

      SUBROUTINE polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
      INTEGER m,n,NMAX,MMAX
      REAL dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
      PARAMETER (NMAX=20,MMAX=20)
      INTEGER j,k
      REAL ymtmp(MMAX),yntmp(NMAX) 
      do 12 j=1,m
         do 11 k=1,n 
            yntmp(k)=ya(j,k)
         enddo 11
         call polint(x2a,yntmp,n,x2,ymtmp(j),dy) 
      enddo 12
      call polint(x1a,ymtmp,m,x1,y,dy)
      return
      END
