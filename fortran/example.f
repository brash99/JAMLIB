***************************************************************************                                                                                   
      PROGRAM EXAMPLE
      IMPLICIT NONE
      INTEGER ipos,i
      CHARACTER flav*2,lib*10,dist*10
      REAL*8 x,Q2,up,get_xF

      !!Current flag choices (as of 09/16/16):
      !! For lib='JAM15': dist='PPDF','T3PPDF','T4PPDF'
      !! For lib='JAM16': dist='FFpion','FFkaon'

      lib='JAM16'
      dist='FFpion'
      ipos=0

      x=0.5
      Q2=10.0
      flav='up'

      !! To obtain value from single posterior:
      CALL GRID_INIT(lib,dist,ipos)
      up = get_xF(x,Q2,flav)
      print *, 'xF(z=0.5,Q2=10) from single posterior (ipos=0):',up

      !!To obtain average value:
      up = 0.D0
      DO i=0,199
         CALL GRID_INIT(lib,dist,i)
         up=up+get_xF(x,Q2,flav)  !!Sum values from each posterior
      ENDDO
      up = up/200.D0            !!Divide by total number of posteriors

      print *,'Average xF(z=0.5,Q2=10):',up

      END PROGRAM

**************************************************************************                                                                                    
