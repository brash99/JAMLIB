***************************************************************************                                                                                   
      PROGRAM EXAMPLE
      IMPLICIT NONE
      INTEGER ipos
      CHARACTER(len=20)::flav,hadron
      REAL*8 z,Q2,dp,get_zF

      ipos=0
      flav='dp'
      hadron='pion'
      z=0.5
      Q2=10.0

      print *, get_zF(ipos,hadron,flav,z,Q2)

      END PROGRAM

**************************************************************************                                                                                    
