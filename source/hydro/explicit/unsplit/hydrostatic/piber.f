      SUBROUTINE PIBER(PS, NX, AK, BK, AKH, BKH, DAK, DBK,
     &     HYDROP, HYDRODP, HYDROMP)
C
      USE mo_grid, only: ie,je,ke,ke1, ptop
C
      IMPLICIT NONE

C
C**** GEOPOT   -   UP:BERECHNUNG DES HYDROSTATISCHEN DRUCKS AN DEN NEBENFL.
C**   AUFRUF   :   CALL PIBER(PS, NX)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BERECHNUNG DES HYDROSTATISCHEN DRUCKS AN DEN NEBENFL.
C**   VERSIONS-
C**   DATUM    :   15.06.
C**                2016
C**
C**   EINGABE-
C**   PARAMETER:   PS: HYDROSTATISCHER BODENDRUCK
C**                NX: ZEITINDEX
C**   AUSGABE-
C**   PARAMETER:   HYDROP:  HYDROSTATISCHER DRUCK AN DEN NEBENFLÄCHEN
C**                HYDRODP: DIFFERENZ HYDROSTATISCHER DRUCK
C**                HYDROMP: MITTELWERT HYDROSTATISCHER DRUCK
C**
C**   VERFASSER:   K.SIECK
      !
      ! Dummy Arguments
      !
      INTEGER, INTENT(IN) :: NX
      !
      REAL, INTENT(IN) :: AK(KE1), BK(KE1), AKH(KE), BKH(KE),
     &                    DAK(KE), DBK(KE)
      !
      REAL, INTENT(IN) :: PS(IE,JE,3)
      !
      REAL, INTENT(OUT) :: HYDROP(IE,JE,KE1,3), HYDRODP(IE,JE,KE,3),
     &                     HYDROMP(IE,JE,KE,3)
!
!     Local variables
!
      REAL, EXTERNAL :: GETP
      INTEGER :: I, J, K
!
      DO K  = 1, KE
        DO J = 1, JE
          DO I = 1, IE
            HYDROP(I, J, K, NX)  =
     &           GETP(AK(K),  BK(K),  PS(I, J, NX), PTOP)
            HYDRODP(I, J, K, NX) =
     &           GETP(DAK(K), DBK(K), PS(I, J, NX), PTOP)
            HYDROMP(I, J, K, NX) =
     &           GETP(AKH(K), BKH(K), PS(I, J, NX), PTOP)
          ENDDO
        ENDDO
      ENDDO

      DO J = 1, JE
        DO I = 1, IE
          HYDROP(I, J, KE1, NX) = PS(I, J, NX)
        ENDDO
      ENDDO

      END SUBROUTINE PIBER
