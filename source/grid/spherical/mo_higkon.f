C
      MODULE MO_HIGKON
C
      IMPLICIT NONE
C
      REAL    ::          EDDLAM , EDDPHI , EDADPHI, DLADDPH, DPHDDLA,
     &                    PI     , ED2DT  , DT2    , DTDEH  ,
     &                    GQ     , GH     , EDG    , DEGRAD , RADDEG ,
     &                    RDRD   , EMRDRD , RDDRM1 , TCRITH , TCRITL ,
     &                    WCPR   , RDWCP  , B234W  , B234E  , NEHDDT
C
C**** HIGKON   -   CB:HILFSGROESSEN
C**
C**   BESCHREIBUNG DER VARIABLEN:
C**   EDDLAM   :   INVERSE MASCHENWEITE IN LAM-RICHTUNG      (1/RAD)
C**   EDDPHI   :   INVERSE MASCHENWEITE IN PHI-RICHTUNG      (1/RAD)
C**   EDADPHI  :   INVERS.(MASCHENWEITE*RERD) IN PHI-RICHTUNG(1/RAD*M)
C**   DLADDPH  :   DLAM/DPHI                                 (--)
C**   DPHDDLA  :   DPHI/DLAM                                 (--)
C**   PI       :   ZAHL PI                                   (--)
C**   ED2DT    :   1/(2*ZEITSCHRITT)                         (1/S)
C**   DT2      :   2*ZEITSCHRITT                             (S)
C**   DTDEH    :   ZEITSCHRITT/3600 S                        (--)
C**   NEHDDT   :   3600 S/ZEITSCHRITT                        (--)
C**   GQ       :   QUADRAT DER ERDBESCHLEUNIGUNG             (M/S**2)**2
C**   GH       :   HAELFTE DER ERDBESCHLEUNIGUNG : G/2       (M/S**2)
C**   EDG      :   INVERSES DER ERDBESCHLEUNIGUNG            (S**2)/M)
C**   DEGRAD   :   UMRECHNUNGSFAKTOR DEGREE IN RAD           (--)
C**   RADDEG   :   UMRECHNUNGSFAKTOR RAD IN DEGREE           (--)
C**   RDRD     :   R/RD                                      (--)
C**   EMRDRD   :   1 - R/RD                                  (--)
C**   RDDRM1   :   RD/R - 1                                  (--)
C**   TCRITH   :   273.16                                    (K)
C**   TCRITL   :   268.16                                    (K)
C**   WCPR     :   1/WCP                                     (KG*K/J)
C**   RDWCP    :   R/WCP                                     (--)
C**   B234W    :   B2W*(B3 - B4W)                            (K**2)
C**   B234E    :   B2E*(B3 - B4E)                            (K**2)
C**
C**   VERFASSER:   D.MAJEWSKI
C     17.11.15: CONVERTED TO MODULE (LARS BUNTEMEYER)
C
      END MODULE MO_HIGKON
