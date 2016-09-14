***************************************************************************                                                                                   
      PROGRAM EXAMPLE
      IMPLICIT NONE
      INTEGER ipos
      CHARACTER flav*2,hadron*10
      REAL*8 z,Q2,up,get_zF

      ipos=0
      flav='up'
      hadron='pion'
      z=0.5
      Q2=10.0

      CALL GRID_INIT(hadron,ipos)
      up = get_zF(flav,z,Q2)
      print *, up

      END PROGRAM

**************************************************************************                                                                                    
