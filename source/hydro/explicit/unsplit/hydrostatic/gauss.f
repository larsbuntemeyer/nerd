C
C     SUBROUTINE GAUSS
C
C@(#) BIBLIOTHEK EM: EUROPA-MODELL / DEUTSCHLAND-MODELL
C@(#) MODUL GAUSS.F, V2.12 VOM 4/27/94, EXTRAHIERT AM 4/28/94
C
C**** GAUSS    -   UP:LOESUNG TRIDIAGONALER GLEICHUNGSSYSTEME
C**   AUFRUF   :   CALL GAUSS (IANF,IEND,M,IAGA)
C**   ENTRIES  :      ---
C**   ZWECK    :   PARALLELE LOESUNG VON IANF-IEND+1 SYSTEMEN VON KE
C**                SIMULTANEN DREIPUNKT-DIFFERENZENGLEICHUNGEN (EIN-
C**                DIMENSIONALE RANDWERTAUFGABE)
C**   VERSIONS-
C**   DATUM    :   20. 2. 89
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   IANF: INDEX DES ERSTEN TRIDIAGONALEN SYSTEMS
C**                IEND: INDEX DES LETZTEN     ,,         ,,
C**                M   : ( = 1,3) BEZEICHNET EINEN BESTIMMTEN SATZ VON
C**                      GLEICHUNGSSYSTEMEN
C**
C**   PARAMETER:      ---
C**
C**   COMMON-
C**   BLOECKE  :   PARAM
C**
C**   METHODE  :   ZEILENWEISE INVERSION ( ELIMINATIONSVERFAHREN VON
C**                GAUSS )
C**   FEHLERBE-
C**   HANDLUNG :      ---
C**   VERFASSER:   G. DOMS
C
      SUBROUTINE GAUSS (IANF,IEND,M,AGA,AGB,AGC,AGD,AGE)
C
      USE MO_PARORG
      IMPLICIT NONE
C
C     ERWARTET WERDEN DIE KOEFFIZIENTENMATRIZEN AGA(I,K,M), AGB(I,K,M)
C     UND AGC(I,K,M) SOWIE DIE INHOMOGENEN VEKTOREN AGD(I,K,M).
C     DIE LOESUNGSVEKTOREN ERSCHEINEN AUF AGE(I,K,M).
C
      INTEGER, INTENT(IN)    :: IANF,IEND,M
C
C     DIMENSIONIERUNG DER ERFORDERLICHEN FELDER
C
      REAL,    INTENT(IN)    :: AGA(IE,KE,4), AGB(IE,KE,4)
      REAL,    INTENT(INOUT) :: AGC(IE,KE,4), AGD(IE,KE,4), AGE(IE,KE,4)

C
C     Local Variables
C
      INTEGER :: I,K
      REAL    :: ZZZ
C
C     TRANSFORMATION DER KOEFFIZINENTENMATRIZEN IN EINE OBERE DREIECKS-
C     FORM MIT 1 IN DER DIAGONALEN

      DO I  = IANF , IEND
         AGC(I,1,M) = AGC(I,1,M)/AGB(I,1,M)
         AGD(I,1,M) = AGD(I,1,M)/AGB(I,1,M)
      ENDDO
      DO K  =   2 , KE-1
         DO I  = IANF, IEND
            ZZZ        = 1. / ( AGB(I,K,M) - AGA(I,K,M)*AGC(I,K-1,M) )
            AGC(I,K,M) = AGC(I,K,M)*ZZZ
            AGD(I,K,M) = ( AGD(I,K,M) - AGA(I,K,M)*AGD(I,K-1,M) )*ZZZ
         ENDDO
      ENDDO
      DO I  = IANF, IEND
         AGE(I,KE,M) = ( AGD(I,KE,M) - AGA(I,KE,M)*AGD(I,KE-1,M) ) /
     &                 ( AGB(I,KE,M) - AGA(I,KE,M)*AGC(I,KE-1,M) )
      ENDDO

C     LOESUNG DURCH RUECKWAERTSSUBSTITUTION

      DO K  = KE-1 , 1 , -1
         DO I  = IANF , IEND
            AGE(I,K,M) = AGD(I,K,M) - AGC(I,K,M)*AGE(I,K+1,M)
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE GAUSS