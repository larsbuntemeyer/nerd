C
      MODULE MO_PARKON
C
      IMPLICIT NONE
C
C@(#) BIBLIOTHEK EM: EUROPA-MODELL / DEUTSCHLAND-MODELL
C@(#) MODUL PARKON.H, V2.12 VOM 4/27/94, EXTRAHIERT AM 4/28/94
C
      REAL    ::          B1    , B2W   , B2E  , B3    , B4W  , B4E   ,
     &                    UC1   , UC2   , UCL  , RHDE  ,
     &                    AKS2  , AKS4  , AKT
C
C**** PARKON   -   CB:AUSGEWAEHLTE PARAMETRISIERUNGSKONSTANTEN
C**
C**   BESCHREIBUNG DER VARIABLEN:
C**   B1       :   PARAMETER ZUR BERECHNUNG DES SAETTIGUNGS-  (PA)
C**   B2W      :   DAMPFDRUCKES UEBER WASSER (W) UND EIS (E)  (--)
C**   B2E      :                 "                            (--)
C**   B3       :                 "                            (K)
C**   B4W      :                 "                            (K)
C**   B4E      :                 "                            (K)
C**   UC1      :   PARAMETER ZUR BERECHNUNG DES BEDECKUNGS-   (--)
C**   UC2      :   GRADES MIT WOLKEN IM UNGESAETTIGTEN FALL   (--)
C**   UCL      :                 "
C**   RHDE     :   PARAMETER ZUR BERECHNUNG DES SCHNEEBEDECK- (--)
C**                TEN BODENS
C**   AKS2     :   PARAMETER BEI DER HORIZONTALDIFFUSION 2.O. (--)
C**   AKS4     :   PARAMETER BEI DER HORIZONTALDIFFUSION 4.O. (--)
C**   AKT      :   VON KARMAN-KONSTANTE                       (--)
C**
C**   VERFASSER:   D.MAJEWSKI
C     17.11.15: CONVERTED TO MODULE (LARS BUNTEMEYER)
C
      END MODULE MO_PARKON
