      REAL FUNCTION  GETP (AK, BK, P, PTOP)
C
      IMPLICIT NONE
C
C     Dummy Arguments
C
      REAL   :: AK, BK, P, PTOP

      GETP = AK + BK * (P - PTOP)

      RETURN
      END FUNCTION GETP
