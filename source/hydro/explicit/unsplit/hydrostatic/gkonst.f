      SUBROUTINE GKONST
     &     (PHI   , RLA  ,
     &      FC    , CPHI , ACPHIR, RMY  , AK    , BK    ,
     &      SISTM, SIVMT , SIVMTI, SICQ  ,
     &      SIGAM , SITAU, SINUE , A1T  , A2T   , VVFH  , TRIGSI,
     &      TRIGSJ, RZ1I , RZ2I  , RZ1J , RZ2J  , IFAXI , IFAXJ,
     &      GCPHI , GACPHIR,DDMPUJ,DDMPVJ)
C
C                  SUBROUTINE GKONST
C
C**** GKONST   -   UP:BERECHNUNG DER KONSTANTEN (FC, PHI, RLA, RMY)
C****                 FELDER, DIE NUR VON LAENGE/BREITE DER GITTERPUNKTE
C****                 ABHAENGEN. AUSSERDEM BERECHNUNG DER MATRIZEN
C****                 (SISTM, SIVMT, SIVMTI, SIGAM, SITAU) UND VEKTOREN
C****                 (SINUE, SICQ) FUER DAS SI-ZEITSCHRITT-VERFAHREN.
C**
C**   AUFRUF   :   CALL GKONST
C**
C**   ENTRIES  :   KEINE
C**
C**   ZWECK    :   BERECHNUNG DER KONSTANTEN (FC, PHI, RLA, RMY)
C**                FELDER, DIE NUR VON LAENGE/BREITE DER GITTERPUNKTE
C**                ABHAENGEN. AUSSERDEM BERECHNUNG DER MATRIZEN
C**                (SISTM, SIVMT, SIVMTI, SIGAM, SITAU) UND VEKTOREN
C**                (SINUE, SICQ) FUER DAS SI-ZEITSCHRITT-VERFAHREN.
C**   VERSIONS-
C**   DATUM    :   09.02.89
C**                06.09.93  STOP DURCH CALL EMABORT AUSGETAUSCHT   GERT
C**                23.09.93  NAMENSAENDERUNG DER NAGLIBROUTINEN     GERT
C**                2007
C**
C**   EXTERNALS:   PHSTOPH, RLSTORL, F02AGE, F04AEE (NAGLIB),
C**                PASSKO, AEROCOM
C**   EINGABE-
C**   PARAMETER:   KEINE
C**   AUSGABE-
C**   PARAMETER:   KEINE
C**
C**   COMMON-
C**   BLOECKE  :   PARAM, ORG, GRID, COMDYN, PHYKON, HIGKON, PROGCHK
C**                PARKON, COMDIA, COMPCST
C**
C**   METHODE  :   BERECHNUNG VON PHI, RLA DURCH RUECKTRANSFORMATION
C**                DER EM-TRANSFORMIERTEN KOORDINATEN;
C**                NORMAL MODES SIND DIE EIGENVEKTOREN DER MATRIX
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
C
      USE MO_COMDIA
      USE MO_PARORG
      USE MO_ORG
      USE MO_GRID
      USE MO_COMDYN
      USE MO_PHYKON
      USE MO_PARKON
      USE MO_HIGKON
      USE MO_PROGCHK
C
      IMPLICIT NONE
C
C     Dummy Arguments
C
      REAL, INTENT(INOUT) ::  PHI   (IE,JE), RLA   (IE,JE),
     &                        FC    (IE,JE), CPHI  (JE, 2),
     &                        ACPHIR(JE, 2), RMY   (IE, JE,3)

      REAL, INTENT(IN)    ::  AK    (KE1  ), BK    (KE1  ),
     &                        SIVMTI(KE,KE)
      REAL, INTENT(INOUT) ::  SICQ  (KE   ), SIVMT (KE,KE),
     &                        SISTM (KE,KE),
     &                        SIGAM (KE,KE), SITAU (KE,KE),
     &                        SINUE (KE   )
      REAL                ::  SIVMTITMP(KE,KE)

      REAL, INTENT(INOUT) ::  A1T   (KE1  ), A2T   (KE1  )

      REAL, INTENT(INOUT) ::  VVFH  (KE   )

      REAL, INTENT(INOUT) ::  RZ1I  ( 6), RZ2I  ( 6),
     &                        RZ1J  ( 6), RZ2J  ( 6)
      INTEGER, INTENT(INOUT) ::  IFAXI (11), IFAXJ (11)
C
      REAL, INTENT(INOUT) ::  TRIGSI(5*(MOIE-1)/2), TRIGSJ(5*(MOJE-1)/2)
      REAL, INTENT(INOUT) ::  GCPHI(MOJE,2), GACPHIR(MOJE,2)

      REAL, INTENT(OUT) :: DDMPUJ(JE), DDMPVJ(JE)
C
C     Local Variables
C
      REAL    :: ZEINH  (KE,KE), ZADUM  (KE,KE ),
     &           ZBDUM  (KE,KE), ZRDUM  (KE    ),
     &           ZIDUM  (KE   ), ZBETAR (KE    ),
     &           ZDPR   (KE   ), ZALAR  (KE    ),
     &           ZA1I   (2* JE), ZA1J   (2* IE ),
     &           ZWIDUM (KE   ), ZVLDUM (KE,KE ),
     &           ZRCONDE(KE   ), ZRCONDV(KE    ),
     &           ZSCALE (KE   ),
     &           ZWORK  (KE*KE+6*KE+100),
     &           ZAF    (KE,KE), 
     &           ZR     (KE   ), ZC     (KE    ),
     &           ZFERR  (KE   ), ZBERR  (KE    ),
     &           ZWORK2 (4*KE )
C
      INTEGER :: IWORK  (2*KE-2), IWORK2 (KE    )
C
      INTEGER :: IPIV   (KE    )
C
      CHARACTER*1 :: YEQUED
C
      INTEGER :: K4,K3,K2,K1,JEVG,JEUG,JEHG,JG,K,JAVG,JAUG,JAHG,J,ILO,
     &           IG,IHI,IEVG,IEUG,IEHG,IAVG,IAUG,IAHG,I,IERR
      REAL    :: ZABNRM,ZDTDDTR,ZDRV,ZDRU,ZDRH,ZBETAMX,ZALOPUR,ZPNU,
     &           ZALOPOR,ZAKS4MX,ZAKS2MX,ZWEIOM,ZRLAS,ZRCOND,ZPUR,ZPOR,
     &           ZPNO,ZPHIS,ZDXMIN
CKS
      REAL :: DXJ(JE), ACDT, CDDAMP
C
C     LWORK >= KE*(KE+6) (= 520 FUER KE=20)
      INTEGER :: LWORK
C      PARAMETER (LWORK=1000)
C
C     External Functions
C     
      REAL, EXTERNAL :: RLSTORL,PHSTOPH,GETP
C
      LWORK = KE*(KE+6) + 100
C
C      FELDER INITIALISIEREN
C
      CALL SETRA(TRIGSI,5*(MOIE-1)/2,0.0)
      CALL SETRA(TRIGSJ,5*(MOJE-1)/2,0.0)
      CALL SETRA(RZ1I,6,0.0)
      CALL SETRA(RZ1J,6,0.0)
      CALL SETRA(RZ2I,6,0.0)
      CALL SETRA(RZ2J,6,0.0)
      CALL SETIA(IFAXI,11,0)
      CALL SETIA(IFAXJ,11,0)
      CALL SETRA(ZEINH,KE*KE,0.0)
      CALL SETRA(ZADUM,KE*KE,0.0)
      CALL SETRA(ZBDUM,KE*KE,0.0)
      CALL SETRA(ZRDUM,KE,0.0)
      CALL SETRA(ZIDUM,KE,0.0)
      CALL SETRA(ZWIDUM,KE,0.0)
      CALL SETRA(ZBETAR,KE,0.0)
      CALL SETRA(ZDPR,KE,0.0)
      CALL SETRA(ZALAR,KE,0.0)
      CALL SETRA(ZA1I,2*JE,0.0)
      CALL SETRA(ZA1J,2*IE,0.0)
      CALL SETRA(ZVLDUM, KE*KE, 0.0)
      CALL SETRA(ZSCALE, KE, 0.0)
      CALL SETRA(ZRCONDE, KE, 0.0)
      CALL SETRA(ZRCONDV, KE, 0.0)
      CALL SETRA(ZWORK, LWORK, 0.0)
      CALL SETIA(IWORK, 2*KE-2, 0)
      CALL SETRA(ZAF, KE*KE, 0.0)
      CALL SETIA(IPIV, KE, 0)
      CALL SETRA(ZR, KE, 0.0)
      CALL SETRA(ZC, KE, 0.0)
      CALL SETRA(ZFERR, KE, 0.0)
      CALL SETRA(ZBERR, KE, 0.0)
      CALL SETRA(ZWORK2, 4*KE, 0.0)
      CALL SETIA(IWORK2, KE, 0)
      CALL SETRA(SIVMTITMP,KE*KE, 0.0)
      CALL SETRA(DDMPUJ, JE, 0.0)
      CALL SETRA(DDMPVJ, JE, 0.0)
C
      ZWEIOM = 4.0*PI/STAG

C     BERECHNUNG VON PHI, RLA, FC, CPHI, ACPHIR UND RMY
C     JG IST DER GLOBALE J-INDEX
      DO J     = 1,JE
         JG          = MYPOSGRD(2) - 3 + J
         IF (JG .LT. MYPOSGRD(6)) JG = MYPOSGRD(6)
         IF (JG .GT. MYPOSGRD(8)) JG = MYPOSGRD(8)
         ZPHIS       = PHILU + (JG - 1)*DPHI
         CPHI  (J,1) = COS(ZPHIS*DEGRAD)
         CPHI  (J,2) = COS((ZPHIS + 0.5*DPHI)*DEGRAD)
         ACPHIR(J,1) = 1.0/(RERD*CPHI(J,1))
         ACPHIR(J,2) = 1.0/(RERD*CPHI(J,2))
CKS
         DXJ(J) = RERD * DLAM * DEGRAD * CPHI(J,1)
CKS
         DO I    = 1,IE
            IG          = MYPOSGRD(1) - 3 + I
            IF (IG .LT. MYPOSGRD(5)) IG = MYPOSGRD(5)
            IF (IG .GT. MYPOSGRD(7)) IG = MYPOSGRD(7)
            ZRLAS       = RLALU + (IG - 1)*DLAM
            IF(ZRLAS.GT.180.0) ZRLAS = ZRLAS - 360.0

            PHI(I,J)   = PHSTOPH(ZPHIS, ZRLAS, POLPHI)*DEGRAD
            RLA(I,J)   = RLSTORL(ZPHIS, ZRLAS, POLPHI, POLLAM)*DEGRAD
         ENDDO
      ENDDO

C     BERECHNUNG VON GCPHI UND GACPHIR IM GESAMTGEBIET
      DO J=1,MOJE
         ZPHIS        = PHILU + (J - 1)*DPHI
         GCPHI  (J,1) = COS(ZPHIS*DEGRAD)
         GCPHI  (J,2) = COS((ZPHIS + 0.5*DPHI)*DEGRAD)
         GACPHIR(J,1) = 1.0/(RERD*GCPHI(J,1))
         GACPHIR(J,2) = 1.0/(RERD*GCPHI(J,2))
      ENDDO
CKS   DIVERGENCE DAMPING COEFFICIENTS
      ACDT = DT * SQRT(MIN( RERD*GCPHI( 1,1)*DLAM*DEGRAD,
     &     RERD*GCPHI(MOJE,1)*DLAM*DEGRAD)**2 + (RERD*DPHI*DEGRAD)**2)
      CDDAMP = CODAMP * ACDT
      DDMPUJ(:) = CDDAMP/DXJ(:)
      DDMPVJ(:) = CDDAMP/(RERD*DPHI*DEGRAD)
CKS

      DO J  = 1,JE
         DO I  = 1,IE
            FC(I,J)  = ZWEIOM*SIN(PHI(I,J))
         ENDDO
      ENDDO

C     BERECHNUNG VON RMY:
C     RMY(I,J,1) IST AM PS, H, QDW-GITTERPUNKT DEFINIERT
C     RMY(I,J,2) IST AM U         -GITTERPUNKT DEFINIERT
C     RMY(I,J,3) IST AM V         -GITTERPUNKT DEFINIERT
      ZDTDDTR = DT/300.0*0.5/MAX( DLAM, DPHI )
C     GLOBALE ANFANGS- UND ENDPUNKTE DER BERECHNUNGEN
      IAHG = 2
      JAHG = 2
      IEHG = MYPOSGRD(7) - 1
      JEHG = MYPOSGRD(8) - 1
      IAUG = 2
      JAUG = 2
      IEUG = MYPOSGRD(7) - 2
      JEUG = MYPOSGRD(8) - 1
      IAVG = 2
      JAVG = 2
      IEVG = MYPOSGRD(7) - 1
      JEVG = MYPOSGRD(8) - 2
      DO J = 1,JE
         JG = MYPOSGRD(2) - 3 + J
         IF (JG .LT. MYPOSGRD(6)) JG = MYPOSGRD(6)
         IF (JG .GT. MYPOSGRD(8)) JG = MYPOSGRD(8)
         DO I = 1,IE
            IG      = MYPOSGRD(1) - 3 + I
            IF (IG .LT. MYPOSGRD(5)) IG = MYPOSGRD(5)
            IF (IG .GT. MYPOSGRD(7)) IG = MYPOSGRD(7)

            ZDRH    = MIN(IG - IAHG + 0.25 , IEHG - IG + 0.25 ,
     &           JG - JAHG + 0.25 , JEHG - JG + 0.25 )
            IF(ZDRH.LE.0.0) THEN
               RMY(I,J,1) = 1.0*ZDTDDTR
            ELSE
               RMY(I,J,1) = (1.0 - TANH(0.5*ZDRH))/TANH(0.5*ZDRH)*
     &              ZDTDDTR
            ENDIF

            ZDRU    = MIN(IG - IAUG + 0.75 , IEUG - IG + 0.75 ,
     &           JG - JAUG + 0.25 , JEUG - JG + 0.25 )
            IF(ZDRU.LE.0.0) THEN
               RMY(I,J,2) = 1.0*ZDTDDTR
            ELSE
               RMY(I,J,2) = (1.0 - TANH(0.5*ZDRU))/TANH(0.5*ZDRU)*
     &              ZDTDDTR
            ENDIF

            ZDRV    = MIN(IG - IAVG + 0.25 , IEVG - IG + 0.25 ,
     &           JG - JAVG + 0.75 , JEVG - JG + 0.75 )
            IF(ZDRV.LE.0.0) THEN
               RMY(I,J,3) = 1.0*ZDTDDTR
            ELSE
               RMY(I,J,3) = (1.0 - TANH(0.5*ZDRV))/TANH(0.5*ZDRV)*
     &                  ZDTDDTR
            ENDIF
         ENDDO
      ENDDO

C     BERECHNUNG DER STRUKTURMATRIX UND DER MATRIX DER VERTIKALEN
C     NORMAL MODES FUER SI-ZEITSCHRITT-VERFAHREN
C     EINHEITSMATRIX
      DO K1     = 1,KE
         DO K2     = 1,KE
            ZEINH(K1,K2) = 0.0
            IF(K1.EQ.K2) THEN
               ZEINH(K1,K2) = 1.0
            ENDIF
         ENDDO
      ENDDO

C     REFERENZATMOSPHAERE

      ZBETAMX   = ALOG(2.0)
      ZPOR      = PTOP
      ZPUR      = 0.0
      IF (LPTOP0) THEN
        ZALOPOR   = 0.0
      ELSE
        ZALOPOR   = ALOG(PTOP)
      END IF
      ZALOPUR   = 0.0

      DO K   = 1,KE
         ZPUR      = GETP(AK(K+1), BK(K+1), SIPSR, PTOP)
         ZALOPUR   = ALOG(ZPUR)
         ZDPR(K)   = ZPUR    - ZPOR
         ZALAR(K)  = ZALOPUR - ZALOPOR
         ZBETAR(K) = MIN(ZBETAMX,1.0 - ZPOR*ZALAR(K)/ZDPR(K))
         ZPOR      = ZPUR
         ZALOPOR   = ZALOPUR
      ENDDO

C     DIE MATRIZEN SIGAM, SITAU UND DER VEKTOR SINUE
      DO K1     = 1,KE
         SINUE(K1)    = ZDPR(K1)
         DO K2     = 1,KE

            IF(K1.LT.K2) THEN
               SIGAM(K1,K2) = 0.0
               SITAU(K1,K2) = RDWCP*SITR/ZDPR(K2)*ZALAR(K2)*ZDPR(K1)
            ELSE IF(K1.EQ.K2) THEN
               SIGAM(K1,K2) = R*ZBETAR(K1)
               SITAU(K1,K2) = RDWCP*SITR*ZBETAR(K1)
            ELSE IF(K1.GT.K2) THEN
               SIGAM(K1,K2) = R*ZALAR(K1)
               SITAU(K1,K2) = 0.0
            ENDIF

         ENDDO
      ENDDO

      K3=0
C     STRUKTURMATRIX UND EIGENVEKTOREN
      DO K1     = 1,KE
         DO K2     = 1,KE
            K3           = MAX(K1,K2) + 1
            SISTM(K1,K2) = ZDPR(K1)/SIPSR

            DO K4     = K3,KE
               SISTM(K1,K2) = SISTM(K1,K2) + RDWCP*ZALAR(K4)**2 *
     &              ZDPR(K1)/ZDPR(K4)
            ENDDO

            IF(K1.LT.K2) THEN
              SISTM(K1,K2) = R*SITR*(SISTM(K1,K2) + RDWCP*ZBETAR(K2)*
     &              ZALAR(K2)*ZDPR(K1)/ZDPR(K2))
            ELSE IF(K1.EQ.K2) THEN
              SISTM(K1,K2) = R*SITR*(SISTM(K1,K2) + RDWCP*ZBETAR(K2)**2)
            ELSE IF(K1.GT.K2) THEN
              SISTM(K1,K2) = R*SITR*(SISTM(K1,K2) + RDWCP*ZBETAR(K1)*
     &              ZALAR(K1))
            ENDIF

            ZBDUM(K1,K2)    = SISTM(K1,K2)
         ENDDO
      ENDDO

      CALL DGEEVX('B','V','V','B',KE,ZBDUM,KE,SICQ,ZWIDUM,ZVLDUM,KE,
     &     SIVMT,KE,ILO,IHI,ZSCALE,ZABNRM,ZRCONDE,ZRCONDV,
     &     ZWORK,LWORK,IWORK,IERR)

      IF (IERR .NE. 0) THEN
         CALL REMARK('GKONST: FEHLER IN SGEEVX')
         CALL EMABORT
      END IF

      CALL BUBBSORT(SICQ,SIVMT,KE,KE)

      DO K2=1,KE
         DO K1=1,KE
            SIVMTITMP(K1,K2)=SIVMT(K1,K2)
         ENDDO
      ENDDO

      CALL DGESVX('E','N',KE,KE,SIVMTITMP,KE,ZAF,KE,IPIV,YEQUED,ZR,
     &     ZC,ZEINH,KE,SIVMTI,KE,ZRCOND,ZFERR,ZBERR,
     &     ZWORK2,IWORK2,IERR)

      IF (IERR .NE. 0) THEN
         CALL REMARK('GKONST: FEHLER IN SGESVX')
         CALL EMABORT
      END IF

C     VERTIKAL VARIIERENDER IMPLIZITHEITSGRAD DER V-DIFFUSION
      DO K  = 1 , KE1
         A1T   (K) = 0.75
      ENDDO
      A1T(KE1)  = 1.2
      A1T(KE)   = 1.2
      A1T(KE-1) = 1.1
      A1T(KE-2) = 1.0
      A1T(KE-3) = 0.90
      A1T(KE-4) = 0.80

      DO K  = 1 , KE1
         A2T   (K) = 1.0 - A1T(K)
      ENDDO

C     VERTIKALKOORDINATEN-PARAMETER: MITTELWERT UND DIFFERENZ
C     EVT. VERTIKAL VARIIERENDER VERSTAERKUNGSFAKTOR DER H-DIFFUSION
CKS   MOVED TO ECAPREP!
      DO K  = 1 , KE
CKS         AKH   (K) = 0.5*(AK(K) + AK(K+1))
C         BKH   (K) = 0.5*(BK(K) + BK(K+1))
C         DAK   (K) =    - AK(K) + AK(K+1)
CKS         DBK   (K) =    - BK(K) + BK(K+1)
         VVFH  (K) = 1.0
      ENDDO

C     DIE HORIZONTALDIFFUSION WIRD SO DIMENSIONIERT, DASS WELLEN DER
C     WELLENLAENGE 2*DX IN 2*DT AUF 1/E IN DER AMPLITUDE REDUZIERT
C     WERDEN; DER ZWEITE WERT IN DER MIN-FUNKTION IST DER STABILITAETS-
C     KRITISCHE WERT.
      ZAKS2MX   =  1.0/( 16.0*DT)
      ZAKS4MX   =  1.0/(128.0*DT)
CKS   CHANGED THE FORMULATION BACK TO OLD REMO STYLE. SHOULD BE REVISITED
C     ONCE REMO-NH IS WORKING!
C      AKS2      =  MIN ( 1.0/(2.0*DT*PI**2) * ZDTDDTR, ZAKS2MX )
C      AKS4      =  MIN ( 1.0/(2.0*DT*PI**4) * ZDTDDTR, ZAKS4MX )
      AKS2      =  MIN ( 1.0/(2.0*DT*PI**2), ZAKS2MX )
      AKS4      =  MIN ( 1.0/(2.0*DT*PI**4), ZAKS4MX )
CKS
CTR  VVFH is always limited by AKS4, even for 2nd order-only diffusion
      DO K  = 1 , KE
         VVFH  (K) =  MIN ( VVFH(K)*AKS4, ZAKS4MX )/AKS4
      ENDDO

C     KFL850, KFL800, KFL500, KFL400, KFL300 BESTIMMEN
      DO K = 1,KE
         ZPNO   = GETP(AK(K  ), BK(K  ), 1.0 E5, PTOP)
         ZPNU   = GETP(AK(K+1), BK(K+1), 1.0 E5, PTOP)
         IF(850.0 E2 .GE. ZPNO .AND. 850.0 E2 .LT. ZPNU) KFL850 = K
         IF(800.0 E2 .GE. ZPNO .AND. 800.0 E2 .LT. ZPNU) KFL800 = K
         IF(500.0 E2 .GE. ZPNO .AND. 500.0 E2 .LT. ZPNU) KFL500 = K
         IF(400.0 E2 .GE. ZPNO .AND. 400.0 E2 .LT. ZPNU) KFL400 = K
         IF(300.0 E2 .GE. ZPNO .AND. 300.0 E2 .LT. ZPNU) KFL300 = K
      ENDDO

C     FFT'S IN I- UND J-RICHTUNG INITIALISIEREN
      CALL PASSKO(TRIGSI, 5*(MOIE-1)/2, MOIE-1, PI, RZ1I, RZ2I, IFAXI)
      CALL PASSKO(TRIGSJ, 5*(MOJE-1)/2, MOJE-1, PI, RZ1J, RZ2J, IFAXJ)

C     MAXIMALE WINDGESCHWINDIGKEIT VBCFL BERECHNEN, DIE BEI DEM
C     GEWAEHLTEN ZEITSCHRITT DT WEGEN DES CFL-KRITERIUMS NICHT
C     UEBERSCHRITTEN WERDEN DARF.
C     ZUR BERECHNUNG MUSS GCPHI BENUTZT WERDEN.
      ZDXMIN = MIN ( RERD*GCPHI( 1,1)*DLAM*DEGRAD,
     &     RERD*GCPHI(MOJE,1)*DLAM*DEGRAD, RERD*DPHI*DEGRAD )
      IF ( LSITS ) THEN
         VBCFL = ZDXMIN/(SQRT(2.0)*DT)
      ELSE
         VBCFL = ZDXMIN/(2.0*SQRT(2.0)*DT) - SQRT(SICQ(1))
      ENDIF


      RETURN
      END SUBROUTINE GKONST
