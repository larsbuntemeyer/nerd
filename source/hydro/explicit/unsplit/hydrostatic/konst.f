      SUBROUTINE KONST
      !
      USE MO_GRID
      USE MO_ORG
      USE MO_PARORG
      USE MO_COMDYN
      USE MO_PHYKON
      USE MO_HIGKON
      USE MO_PARKON
      USE MO_FAKTINF
      !
      IMPLICIT NONE
      !
C
C**** KONST    -   UP:BESETZEN VON KONSTANTEN IN DEN COMMON-BLOECKEN
C****                 *PARAM* , *ORG*, *COMDYN*,*PHYKON*, *HIGKON*,
C****                 *PARKON*, *COMDIA*, *COMTER*
C**   AUFRUF   :   CALL KONST IN HP *EMORG*
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BESETZEN VON KONSTANTEN IN DEN COMMON-BLOECKEN
C**               *PARAM*, *ORG*, *COMDYN*,*PHYKON*, *HIGKON*, *PARKON*,
C**               *COMDIA*, *COMTER*
C**   VERSIONS-
C**   DATUM    :   01.02.89
C**                2007
C**
C**   EXTERNALS:   KEINE
C**
C**   EINGABE-
C**   PARAMETER:   KEINE
C**   AUSGABE-
C**   PARAMETER:   KEINE
C**
C**   COMMON-
C**   BLOECKE  :   GRID, ORG, COMDYN, PHYKON, HIGKON, PARKON,
C**                DATTER, COMTER, UNITCH, UNITNR
C**
C**   METHODE  :   VORBESETZEN MIT ARITHMETISCHEN STATEMENTS
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
C
C     KONSTANTEN AUS CB *ORG*
C     WEITERE VARIABLEN AUS CB *ORG* WERDEN IM HP *EMORG* UND IM
C     UP *MEDEA* BESETZT
      IAA    =   1
      IEA    =  IE
      JAA    =   1
      JEA    =  JE

      IAH    = IAA + 2
      JAH    = JAA + 2
      IEH    = IEA - 2
      JEH    = JEA - 2

      IAHGG  = MYPOSGRD(5) + 1
      JAHGG  = MYPOSGRD(6) + 1
      IEHGG  = MYPOSGRD(7) - 1
      JEHGG  = MYPOSGRD(8) - 1

      IAU    = IAA + 2
      JAU    = JAA + 2
      IF (NEIGHBOR(3) .EQ. -1) THEN
         IEU = IEA - 3
      ELSE
         IEU = IEA - 2
      ENDIF
      JEU    = JEA - 2

      IAV    = IAA + 2
      JAV    = JAA + 2
      IEV    = IEA - 2
      IF (NEIGHBOR(2) .EQ. -1) THEN
         JEV = JEA - 3
      ELSE
         JEV = JEA - 2
      ENDIF

C     BESTIMMUNG DER BERECHNUNGSINDICES
      IF (NEIGHBOR(1) .EQ. -1) THEN
         IAHCOMP = 1
      ELSE
         IAHCOMP = IAH
      ENDIF
      IF (NEIGHBOR(2) .EQ. -1) THEN
         JEHCOMP = JE
      ELSE
         JEHCOMP = JEH
      ENDIF
      IF (NEIGHBOR(3) .EQ. -1) THEN
         IEHCOMP = IE
      ELSE
         IEHCOMP = IEH
      ENDIF
      IF (NEIGHBOR(4) .EQ. -1) THEN
         JAHCOMP = 1
      ELSE
         JAHCOMP = JAH
      ENDIF

      FAKRMY = 1.0

C     KONSTANTEN AUS CB *COMDYN*
      SITR   = 300.0
      SIPSR  = 900.0 E2

C     KONSTANTEN AUS CB *PHYKON*
      T0    =   273.16
      R     =   287.05
      RD    =   461.51
      WCP   =  1005.0
      WLK   =     2.501 E6
      WLF   =     0.334 E6
      WLS   =     2.835 E6
      G     =     9.80665
      RERD  =  6371.229 E3
      STAG  = 86164.09054
      RHF   =  1000.0
      SIGMA =     5.6697 E-8
      SOKO  =  1368.0

C     KONSTANTEN AUS CB *PARKON*
      B1    =   610.78
      B2W   =    17.2693882
      B2E   =    21.8745584
      B3    =   273.16
      B4W   =    35.86
      B4E   =     7.66

      UC1   =     0.8
      UC2   =  SQRT(3.0)
      UCL   =     1.00

      RHDE  =     1.0/5.0*1000.0
      AKT   =     0.4

C     KONSTANTEN AUS CB *HIGKON*
      PI      =   4.0*ATAN(1.0)
      EDDLAM  =   1.0/(DLAM*PI/180.0)
      EDDPHI  =   1.0/(DPHI*PI/180.0)
      EDADPHI =   EDDPHI/RERD
      DLADDPH =   DLAM/DPHI
      DPHDDLA =   DPHI/DLAM

C     DIE KONSTANTEN ED2DT, DT2, DTDEH, NEHDDT WERDEN IM UPS *PROGEXP*
C     BZW. *PROGSLA* BERECHNET, WEIL SIE VOM ZEITSCHRITT DT ABHAENGEN.
      GQ      =   G**2
      GH      =   G*0.5
      EDG     =   1.0/G
      DEGRAD  =   PI/180.0
      RADDEG  =   180.0/PI
      RDRD    =   R/RD
      EMRDRD  =   1.0 - RDRD
      RDDRM1  =   RD/R - 1.0
      WCPR    =   1.0/WCP
      RDWCP   =   R  /WCP
      B234W   =   B2W*(B3 - B4W)
      B234E   =   B2E*(B3 - B4E)
C
C     KONSTANTEN AUS CB *FAKTINF*
C
      FAKINF   = 10000.
      EDFAKINF = 1./FAKINF
      NFRAC    = 3
C
      RETURN
      END SUBROUTINE KONST
