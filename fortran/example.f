***************************************************************************                                                                                   
      PROGRAM EXAMPLE
      IMPLICIT NONE
      INTEGER ipos
      CHARACTER flav*2,lib*10,dist*10
      REAL*8 x,Q2,up,get_xF

      lib='JAM15'
      dist='PPDF'
      ipos=0

      x=0.5
      Q2=10.0
      flav='up'

      CALL GRID_INIT(lib,dist,ipos)
      up = get_xF(x,Q2,flav)
      print *, up

      END PROGRAM

**************************************************************************                                                                                    
