C
C     ----------------------------------------------------------------
C
      SUBROUTINE COPYRE(PA,PB,KLEN)
C
C**   *SUBROUTINE* *COPYRE*  COPY *KLEN* REAL VALUES FROM
C                            *PA* TO *PB*.
C
      IMPLICIT NONE
C
C     Dummy Arguments
C
      INTEGER, INTENT(IN)    :: KLEN
      REAL,    INTENT(IN)    :: PA(KLEN)
      REAL,    INTENT(INOUT) :: PB(KLEN)
C
C     Local Variables
C
      INTEGER :: JK
C     ----------------------------------------------------------------
C
C*        1.       COPY PA INTO PB.
C                  ---- -- ---- ---
      DO JK=1,KLEN
         PB(JK)=PA(JK)
      ENDDO
C
C     ----------------------------------------------------------------
C
      RETURN
      END SUBROUTINE COPYRE
