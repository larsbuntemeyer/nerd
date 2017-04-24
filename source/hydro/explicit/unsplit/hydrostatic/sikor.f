      SUBROUTINE SIKOR
     &   (U     , V    , T    , PS   , ACPHIR , CPHI  ,
     &    SISTM , SIGAM, SITAU, SINUE, SIVMT  , SIVMTI, SICQ,
     &    TRIGSI, RZ1I , RZ2I , IFAXI, GACPHIR, GCPHI ,
     &    PINT  , DWDT , AK   , BK   , AKH    , BKH   , DAK,
     &    DBK   , HYDROP,HYDRODP,HYDROMP)
C
C
C@(#) BIBLIOTHEK EM: EUROPA-MODELL / DEUTSCHLAND-MODELL
C@(#) MODUL SIKOR.F, V2.12 VOM 4/27/94, EXTRAHIERT AM 4/28/94
C
C**** SIKOR    -   BERECHNUNG DER SEMI-IMPLIZITEN KORREKTUREN
C**   AUFRUF   :   CALL SIKOR
C**   ENTRIES  :      ---
C**   ZWECK    :   BERECHNUNG DER SEMI-IMPLIZITEN KORREKTURTERME FUER
C**                U, V, T UND PS.
C**   VERSIONS-
C**   DATUM    :   29.3.89
C**
C**   EXTERNALS:   GAELKO , SOLVE , SETRA ,
C**                MXM, MXV, SGEMM, SGEMV
C**
C**   EINGABE-
C**   PARAMETER:      ---
C**   AUSGABE-
C**   PARAMETER:      ---
C**
C**   COMMON-
C**   BLOECKE  :   ORG, HIGKON, COMDYN, PHYKON
C**
C**   METHODE  :   ENTKOPPLUNG DER HELMHOLTZGLEICHUNG DES SEMI-IMPLIZI-
C**                TEN ZEITSCHRITTVERFAHRENS DURCH DIAGONALISIERUNG DER
C**                STRUKTURMATRIX; LOESUNG DER ENTKOPPELTEN GLEICHUNGEN
C**                MIT FFA-ALGORITHMUS; BERECHNUNG DER KORREKTURTERME
C**                UND ENTSPRECHENDE MODIFIKATION FUER U, V, T UND PS.
C**   FEHLERBE-
C**   HANDLUNG :      ---
C**   VERFASSER:   G.DOMS

C     VERSION MIT CRAY-ROUTINEN FUER MATRIXOPERATIONEN
C     (AUSKOMMENTIERT) UND FTN77-VERSION DES FFT-SOLVERS
C
      USE MO_PARORG
      USE MO_ORG
      USE MO_HIGKON
      USE MO_COMDYN
      USE MO_PHYKON
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C     Dummy Arguments
C
      REAL,    INTENT(INOUT) ::  U(IE,JE,KE,3),  V(IE,JE,KE,3),
     &                           T(IE,JE,KE,3), PS(IE,JE,   3),
     &                           PINT(IE,JE,KE1,3), DWDT(IE,JE,KE,3)

      REAL,    INTENT(IN)    ::  ACPHIR (  JE,2), CPHI (  JE,2)
      REAL,    INTENT(IN)    ::  GACPHIR(MOJE,2), GCPHI(MOJE,2)

      REAL,    INTENT(IN)    ::  SISTM(KE,KE),  SIGAM(KE,KE),
     &                           SITAU(KE,KE),  SINUE(KE   ),
     &                           SIVMT(KE,KE), SIVMTI(KE,KE),
     &                           SICQ (   KE)
      REAL,    INTENT(IN)    ::  AK(KE1), BK(KE1), AKH(KE)  ,
     &                           BKH(KE), DAK(KE), DBK(KE)

      REAL,    INTENT(IN)    ::  TRIGSI(5*(MOIE-1)/2), RZ1I(6), RZ2I(6)
      INTEGER, INTENT(IN)    ::  IFAXI(11)
CKS
      REAL, INTENT(INOUT) :: HYDROP(IE,JE,KE1,3), HYDRODP(IE,JE,KE,3),
     &                       HYDROMP(IE,JE,KE,3)
C
C-----------------------------------------------------------------------
C     Local Declarations

C     ZREBUF MUSS NUR MIT DEN LOKALEN DIMENSIONEN DIMENSIONIERT WERDEN (UND
C     NICHT ETWA MIT IALLOGI, IALLOGJ), DA EINE REGULAERE GEBIETSZERLEGUNG
C     VORAUSGESETZT WIRD. (IE*JE*KE,7) WAERE AUCH AUSREICHEND, WENN ZREBUF(1,3)
C     NACH DEM KOPIEREN DES INHALTS AUF U BZW. V GLEICH WIEDER VERWENDET WIRD.
C
      REAL :: ZREBUF(IE*JE*KE,8)

      REAL :: ZDEXP (IE,JE,KE),
     &        ZDEXD (IE,JE,KE),
     &        ZDTTD (IE,JE,KE), ZDEXDS(IE,JE),
     &        ZF    (IE,JE,KE),
     &        ZSIGAT(IE,JE)

      REAL ::   ZFAKTU(MOJE)
      REAL ::   ZAP   (MOJE), ZBP(MOJE,KE), ZBS(MOJE),
     &          ZCP   (MOJE), ZDI(MOIE)   , ZRLAM(KE)

      REAL ::   ZDPSI(IE,JE), ZDTSI(IE,JE)
      REAL ::   ZPS  (IE,JE)
      REAL ::   PSDT (IE,JE)
      REAL :: CYK, ZAPL, ZBSL, ZCPL, ZD2X, ZD2Y, ZDLAMQ, ZDLDDPQ,  
     &        ZDTTPS, ZDTTT, ZDTTUX, ZDTTVO, ZDTTVU, ZDTTVY, ZDUSI,  
     &        ZFAKT1, ZFAKT3, ZFAKTV, ZDLDDTQ, ZDVSI
      REAL :: ZFAKTZ, ZFKPS
C
      INTEGER :: I, IL, IMESLEN, J, K, KS, MGAUSS, MT, NFFT
      INTEGER :: KP1
C
C-----------------------------------------------------------------------
C     SENDEN DER INNEREN RANDZEILE VON U UND V

      IF (NEIGHBOR(3) .NE. -1) THEN
         IMESLEN = 0
         DO K = 1 ,KE
            DO J = JAH , JEH
               IMESLEN = IMESLEN + 1
               ZREBUF(IMESLEN,1) = U(IEH,J,K,NE)
            ENDDO
         ENDDO
         TYPE  = 80
         COUNT = IMESLEN*1
         DEST  = NEIGHBOR(3)
         CALL PSENDR(ZREBUF(1,1))
      ENDIF

      IF (NEIGHBOR(2) .NE. -1) THEN
         IMESLEN = 0
         DO K = 1 ,KE
            DO I = IAH , IEH
               IMESLEN = IMESLEN + 1
               ZREBUF(IMESLEN,2) = V(I,JEH,K,NE)
            ENDDO
         ENDDO
         TYPE  = 81
         COUNT = IMESLEN*1
         DEST  = NEIGHBOR(2)
         CALL PSENDR(ZREBUF(1,2))
      ENDIF

C     1. BERECHNUNG DES 'EXPLIZITEN' POTENTIALS DTT(P),EX (ZDEXP)
C        UND DER 'EXPLIZITEN' DIVERGENZ DTT(D),EX (ZDEXD)

C     POTENTIAL- UND DIVERGENZFELDER MIT NULL VORBESETZEN
C     CALL SETRA(ZDEXD,IEJEKE,0.0)
C     CALL SETRA(ZDTTD,IEJEKE,0.0)
C     CALL SETRA(ZDEXP,IEJEKE,0.0)
      DO K = 1 , KE
         ZRLAM(K) = SIGAM(K,1)
         DO J = 1, JE
            DO I = 1, IE
               ZDEXD(I,J,K) = 0.0
               ZDTTD(I,J,K) = 0.0
               ZDEXP(I,J,K) = 0.0
            ENDDO
         ENDDO
      ENDDO

C     EMPFANG DER INNEREN RANDZEILE VON U UND V
      CYK=2
      TAGCOUNT=2
      TAGTABLE(1)=80
      TAGTABLE(2)=81
      IF (NEIGHBOR(1) .EQ. -1) CYK = CYK - 1
      IF (NEIGHBOR(4) .EQ. -1) CYK = CYK - 1
      IL=0
      DO  WHILE (IL.LT.CYK)
         CALL PTEST
         CALL PRECVR(ZREBUF(1,3))

         IF (SOURCE.EQ.NEIGHBOR(1)) THEN
            IMESLEN = 0
            DO K = 1 , KE
               DO J = JAH , JEH
                  IMESLEN = IMESLEN + 1
                  U(IAH-1,J,K,NE) = ZREBUF(IMESLEN,3)
               ENDDO
            ENDDO
            IL = IL+1
         ENDIF
         IF (SOURCE.EQ.NEIGHBOR(4)) THEN
            IMESLEN = 0
            DO K = 1 , KE
               DO I = IAH , IEH
                  IMESLEN = IMESLEN + 1
                  V(I,JAH-1,K,NE) = ZREBUF(IMESLEN,3)
               ENDDO
            ENDDO
            IL = IL+1
         ENDIF
      ENDDO

      CALL PSTOP

C     BERECHNUNG DES EXPLIZITEN POTENTIALS
      ZFKPS = R*SITR/SIPSR
      CALL SETRA(ZSIGAT,IEJE,0.0)

      DO J = JAH , JEH
        DO I = IAH , IEH
          ZPS (I,J) = PS (I,J,NE)
        ENDDO
      ENDDO

      DO K = KE, 1, -1
         DO J = JAH , JEH
            DO I = IAH , IEH
               ZDTTPS       = PS(I,J,  NE) - 2.*PS(I,J,  NJ) +
     &              PS(I,J,  NA)
               ZDTTT        = T (I,J,K,NE) - 2.*T (I,J,K,NJ) +
     &              T (I,J,K,NA)
               ZDEXP(I,J,K) = ZSIGAT(I,J)   + ZDTTT*SIGAM(K,K)
               ZDEXP(I,J,K) = ZDEXP (I,J,K) + ZFKPS*ZDTTPS
               ZSIGAT(I,J)  = ZSIGAT(I,J)   + ZDTTT*ZRLAM(K)
C              BERECHNUNG DER EXPLIZITEN DIVERGENZ
               ZDTTUX = U(I  ,J,K,NE) -2.*U(I  ,J,K,NJ) + U(I  ,J,K,NA)
     &                - U(I-1,J,K,NE) +2.*U(I-1,J,K,NJ) - U(I-1,J,K,NA)
               ZDTTVO = V(I,J  ,K,NE) -2.*V(I,J  ,K,NJ) + V(I,J  ,K,NA)
               ZDTTVU = V(I,J-1,K,NE) -2.*V(I,J-1,K,NJ) + V(I,J-1,K,NA)
               ZDTTVY = ZDTTVO*CPHI(J,2) - ZDTTVU*CPHI(J-1,2)
               ZDEXD(I,J,K) = ACPHIR(J,1)*
     &              ( EDDLAM*ZDTTUX + EDDPHI*ZDTTVY )
            ENDDO
         ENDDO
      ENDDO

C     AUSTAUSCH DER INNEREN RANDZEILEN VON ZDEXP

      IF (NEIGHBOR(1) .NE. -1) THEN
         IMESLEN = 0
         DO K = 1 , KE
            DO J = JAH , JEH
               IMESLEN = IMESLEN + 1
               ZREBUF(IMESLEN,4) = ZDEXP(IAH,J,K)
            ENDDO
         ENDDO
         TYPE  = 84
         COUNT = IMESLEN*1
         DEST  = NEIGHBOR(1)
         CALL PSENDR(ZREBUF(1,4))
      ENDIF

      IF (NEIGHBOR(2) .NE. -1) THEN
         IMESLEN = 0
         DO K = 1 , KE
            DO I = IAH , IEH
               IMESLEN = IMESLEN + 1
               ZREBUF(IMESLEN,5) = ZDEXP(I,JEH,K)
            ENDDO
         ENDDO
         TYPE  = 85
         COUNT = IMESLEN*1
         DEST  = NEIGHBOR(2)
         CALL PSENDR(ZREBUF(1,5))
      ENDIF

      IF (NEIGHBOR(3) .NE. -1) THEN
         IMESLEN = 0
         DO K = 1 , KE
            DO J = JAH , JEH
               IMESLEN = IMESLEN + 1
               ZREBUF(IMESLEN,6) = ZDEXP(IEH,J,K)
            ENDDO
         ENDDO
         TYPE  = 86
         COUNT = IMESLEN*1
         DEST  = NEIGHBOR(3)
         CALL PSENDR(ZREBUF(1,6))
      END IF

      IF (NEIGHBOR(4) .NE. -1) THEN
         IMESLEN = 0
         DO K = 1 , KE
            DO I = IAH , IEH
               IMESLEN = IMESLEN + 1
               ZREBUF(IMESLEN,7) = ZDEXP(I,JAH,K)
            ENDDO
         ENDDO
         TYPE  = 87
         COUNT = IMESLEN*1
         DEST  = NEIGHBOR(4)
         CALL PSENDR(ZREBUF(1,7))
      END IF

C     2. KOEFFIZIENTEN UND QUELLFUNKTION DER DISKRETISIERTEN UND MIT
C        DEM FAKTOR (DLAM*A*CPHI)**2 MULTIPLIZIERTEN HELMHOLTZGLEICHUNG
C        BEREITSTELLEN

      ZDLDDPQ = DLADDPH**2
      ZDLDDTQ = (1./(EDDLAM*DT))**2
!DIR$ UNROLL
      DO J = JAHGG , JEHGG
         ZAP(J) = ZDLDDPQ*GCPHI(J,1)*GCPHI(J-1,2)
         ZCP(J) = ZDLDDPQ*GCPHI(J,1)*GCPHI(J  ,2)
         ZBS(J) = - ( 2. + ZAP(J) + ZCP(J) )
      ENDDO

C     EMPFANG DER INNEREN RANDZEILEN VON ZDEXP
      CYK         = 4
      TAGCOUNT    = 4
      TAGTABLE(1) = 84
      TAGTABLE(2) = 85
      TAGTABLE(3) = 86
      TAGTABLE(4) = 87
      IF (NEIGHBOR(4) .EQ. -1) CYK = CYK - 1
      IF (NEIGHBOR(3) .EQ. -1) CYK = CYK - 1
      IF (NEIGHBOR(2) .EQ. -1) CYK = CYK - 1
      IF (NEIGHBOR(1) .EQ. -1) CYK = CYK - 1
      IL = 0
      DO  WHILE (IL.LT.CYK)
         CALL PTEST
         CALL PRECVR(ZREBUF(1,8))

         IF (SOURCE.EQ.NEIGHBOR(1)) THEN
            IMESLEN = 0
            DO K = 1 , KE
               DO J = JAH , JEH
                  IMESLEN = IMESLEN + 1
                  ZDEXP(IAH-1,J,K) = ZREBUF(IMESLEN,8)
               ENDDO
            ENDDO
            IL = IL+1
         ENDIF
         IF (SOURCE.EQ.NEIGHBOR(2)) THEN
            IMESLEN = 0
            DO K = 1 ,KE
               DO I = IAH , IEH
                  IMESLEN = IMESLEN + 1
                  ZDEXP(I,JEH+1,K) = ZREBUF(IMESLEN,8)
               ENDDO
            ENDDO
            IL = IL+1
         ENDIF

         IF (SOURCE.EQ.NEIGHBOR(3)) THEN
            IMESLEN = 0
            DO K = 1 ,KE
               DO J = JAH , JEH
                  IMESLEN = IMESLEN + 1
                  ZDEXP(IEH+1,J,K) = ZREBUF(IMESLEN,8)
               ENDDO
            ENDDO
            IL = IL+1
         ENDIF

         IF (SOURCE.EQ.NEIGHBOR(4)) THEN
            IMESLEN = 0
            DO K = 1 ,KE
               DO I = IAH , IEH
                  IMESLEN = IMESLEN + 1
                  ZDEXP(I,JAH-1,K) = ZREBUF(IMESLEN,8)
               ENDDO
            ENDDO
            IL = IL+1
         ENDIF
      END DO

      CALL PSTOP

C     QUELLFUNKTION BERECHNEN UND AUF ZDTTD ZWISCHENSPEICHERN

      ZDLAMQ = (1./EDDLAM)**2
      DO K = 1 , KE
         DO J = JAH , JEH
            ZFAKT3 = ZDLAMQ/(ACPHIR(J,1)**2)
            ZAPL   = ZDLDDPQ * CPHI(J,1) * CPHI(J-1,2)
            ZCPL   = ZDLDDPQ * CPHI(J,1) * CPHI(J  ,2)
            ZBSL   = - ( 2. + ZAPL + ZCPL )
            DO I = IAH , IEH
               ZD2X = ZDEXP(I+1,J,K) + ZBSL*ZDEXP(I,J,K) +
     &              ZDEXP(I-1,J,K)
               ZD2Y = ZAPL*ZDEXP(I,J-1,K) + ZCPL*ZDEXP(I,J+1,K)
               ZDTTD(I,J,K) = DT*(ZD2X+ZD2Y) - ZFAKT3*ZDEXD(I,J,K)
            ENDDO
         ENDDO
      ENDDO


C     3. TRANSFORMATION DER QUELLFUNKTION (-> ZDEXDS), RECHTE SEITEN DER
C        ENTKOPPELTEN HELMHOLTZGLEICHUNGEN BEREITSTELLEN (-> ZF), LOESEN
C        UND ERGEBNISSE (->ZX) FUER TRANSFORMIERTE DIVERGENZ DTT(D*)
C        ABSPEICHERN (->ZDEXD)

      MGAUSS = MYGAUSS_IUP - MYGAUSS_ILO + 3
      NFFT   = MYFFT_JUP   - MYFFT_JLO   + 3
      MT     = 5*(MOIE - 1)/2
      CALL SETRA(ZF,IEJEKE,0.0)

      DO K = 1 , KE
         CALL SETRA(ZDEXDS,IEJE,0.0)
         DO KS = 1 , KE
            DO J = JAH, JEH
               DO I = IAH, IEH
                  ZDEXDS(I,J) = ZDEXDS(I,J) + ZDTTD(I,J,KS)*SIVMT(KS,K)
               ENDDO
            ENDDO
         ENDDO
C        RECHTE SEITE DER HELMHOLTZGLEICHUNG (LOKAL)
         ZFAKTZ = 1./(SICQ(K)*DT**2)

         DO J = JAH , JEH
            DO I = IAH , IEH
               ZF(I,J,K) = ZFAKTZ*ZDEXDS(I,J)
            ENDDO
         ENDDO

C        LETZTER KOEFFIZIENT DER HELMHOLTZGLEICHUNG (GLOBAL)
C        DIESER KOEFFIZIENT IST ABHAENGIG VON DER SCHICHT K
         ZFAKT1 = ZDLDDTQ/SICQ(K)
         DO J = JAHGG , JEHGG
            ZFAKT3 = ZFAKT1/(GACPHIR(J,1)**2)
            ZBP(J,K) = ZBS(J) - ZFAKT3
         ENDDO

      ENDDO

C     HELMHOLTZGLEICHUNG LOESEN UND ERGEBNISSE FUER DTT(D*) AUF ZDEXD
C     ZWISCHENSPEICHERN
C     (KOEFFIZIENTEN DER GAUSS-ELIMINATION (COEF, BPT: UP GAELKO)
C     WERDEN JETZT IN SOLVE3D BERECHNET)

C      CALL GAELKO( ZCOEF,ZBPT,IEG,JEG,ZAP,ZBP,ZCP,ZDI,IEG1,JEG1,PI )
C      MT = 5*IEG1/2

      CALL SOLVE3D (ZDEXD,ZF,ZAP,ZBP,ZCP,ZDI,PI,TRIGSI,RZ1I,RZ2I,IFAXI,
     &              NFFT,MGAUSS,MT)

!DIR$ TASK

C     4. RUECKTRANSFORMATION DER LOESUNG (->DTT(D)) UND BERECHNUNG DES
C        POTENTIALTERMS DTT(P)

      CALL SETRA(ZDTTD,IEJEKE,0.0)
      DO K = 1 , KE
         DO KS = 1 , KE
            DO J  = JAH , JEH
               DO I  = IAH , IEH
                  ZDTTD(I,J,K) = ZDTTD(I,J,K) +
     &                 ZDEXD(I,J,KS)*SIVMTI(KS,K)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO K = 1 , KE
         DO KS = 1 , KE
            DO J = JAH , JEH
               DO I = IAH , IEH
                  ZDEXP(I,J,K) = ZDEXP(I,J,K) -
     &                 DT*ZDTTD(I,J,KS)*SISTM(KS,K)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

C     SENDEN DER RANDZEILEN VON ZDEXP AN LINKEN UND UNTEREN NACHBARN
      IF (NEIGHBOR(1) .NE. -1) THEN
         IMESLEN = 0
         DO K = 1 , KE
            DO J = JAH , JEH
               IMESLEN = IMESLEN + 1
               ZREBUF(IMESLEN,1) = ZDEXP(IAH,J,K)
            ENDDO
         ENDDO
         TYPE  = 80
         COUNT = IMESLEN*1
         DEST  = NEIGHBOR(1)
         CALL PSENDR(ZREBUF(1,1))
      END IF

      IF (NEIGHBOR(4) .NE. -1) THEN
         IMESLEN = 0
         DO K = 1 , KE
            DO I = IAH , IEH
               IMESLEN = IMESLEN + 1
               ZREBUF(IMESLEN,2) = ZDEXP(I,JAH,K)
            ENDDO
         ENDDO
         TYPE  = 81
         COUNT = IMESLEN*1
         DEST  = NEIGHBOR(4)
         CALL PSENDR(ZREBUF(1,2))
      END IF

C     5. SEMI-IMPLIZITE KORREKTUREN BERECHNEN UND PROGNOSTISCHE
C        VARIABLE ENTSPRECHEND AUFDATIEREN

C     A.) TEMPERATUR-KORREKTUR

      DO K = 1 , KE

         CALL SETRA(ZDTSI,IEJE,0.0)
         DO KS = 1 , K
            DO J = JAH , JEH
               DO I = IAH , IEH
                  ZDTSI(I,J) = ZDTSI(I,J) + ZDTTD(I,J,KS)*SITAU(KS,K)
               ENDDO
            ENDDO
         ENDDO

         DO J = JAH , JEH
            DO I = IAH , IEH
               T(I,J,K,NE) = T(I,J,K,NE) -DT*ZDTSI(I,J)
            ENDDO
         ENDDO

      ENDDO

C     B.) BODENDRUCK-KORREKTUR

      CALL SETRA(ZDPSI,IEJE,0.0)
      DO KS = 1 , KE
         DO J = JAH,JEH
            DO I = IAH,IEH
               ZDPSI(I,J) = ZDPSI(I,J) + ZDTTD(I,J,KS)*SINUE(KS)
            ENDDO
         ENDDO
      ENDDO
C
      DO J = JAH,JEH
         DO I = IAH,IEH
            PS(I,J,NE) = PS(I,J,NE) - DT*ZDPSI(I,J)
         ENDDO
      ENDDO

C     EMPFANG DER RANDZEILEN VON ZDEXP VOM RECHTEN UND OBEREN NACHBARN

      CYK         = 2
      TAGCOUNT    = 2
      TAGTABLE(1) = 80
      TAGTABLE(2) = 81
      IF (NEIGHBOR(3) .EQ. -1) CYK = CYK - 1
      IF (NEIGHBOR(2) .EQ. -1) CYK = CYK - 1
      IL = 0
      DO  WHILE (IL.LT.CYK)
         CALL PTEST
         CALL PRECVR(ZREBUF(1,3))

         IF (SOURCE.EQ.NEIGHBOR(2)) THEN
            IMESLEN = 0
            DO K = 1 , KE
               DO I = IAH , IEH
                  IMESLEN = IMESLEN + 1
                  ZDEXP(I,JEH+1,K) = ZREBUF(IMESLEN,3)
               ENDDO
            ENDDO
            IL = IL+1
         ENDIF

         IF (SOURCE.EQ.NEIGHBOR(3)) THEN
         IMESLEN = 0
         DO K = 1 ,KE
            DO J = JAH , JEH
               IMESLEN = IMESLEN + 1
               ZDEXP(IEH+1,J,K) = ZREBUF(IMESLEN,3)
            ENDDO
         ENDDO
         IL = IL+1
      ENDIF
      ENDDO

C     C.) U UND V - KORREKTUR

      ZFAKTV = DT*EDADPHI

      DO K = 1 , KE
         DO J = JAU , JEU
            ZFAKTU(J) = DT*EDDLAM*ACPHIR(J,1)
            DO I = IAU , IEU
               ZDUSI       = ZFAKTU(J)*( ZDEXP(I+1,J,K) - ZDEXP(I,J,K) )
               U(I,J,K,NE) = U(I,J,K,NE) - ZDUSI
            ENDDO
         ENDDO

         DO J = JAV , JEV
            DO I = IAV , IEV
               ZDVSI       = ZFAKTV*( ZDEXP(I,J+1,K) - ZDEXP(I,J,K) )
               V(I,J,K,NE) = V(I,J,K,NE) - ZDVSI
            ENDDO
         ENDDO
      ENDDO

C============================================================
C  Code from remo-nh
C============================================================
C
C     COMPUTE HYDROSTATIC PRESSURES
C
      CALL PIBER(PS, NE, AK, BK, AKH, BKH, DAK, DBK,
     &     HYDROP, HYDRODP, HYDROMP)

      DO J  = JAH , JEH
        DO I  = IAH , IEH
          PINT(I,J,1,NE) = HYDROP(I,J,1,NE)
          PSDT(I,J)      = PS(I,J,NE) - ZPS(I,J)
        ENDDO
      ENDDO


      DO K  = 2 , KE1
        DO J  = JAH , JEH
          DO I  = IAH , IEH
            PINT(I,J,K,NE) = PINT(I,J,K,NE) + BK(K)*PSDT(I,J)
          ENDDO
        ENDDO
      ENDDO

      DO K  = 1 , KE
        KP1 = K + 1
        DO J  = JAH , JEH
          DO I  = IAH , IEH
            DWDT(I,J,K,NE) = (PINT  (I,J,KP1,NE) - PINT  (I,J,K,NE))
     &                      /(HYDROP(I,J,KP1,NE) - HYDROP(I,J,K,NE))
          ENDDO
        ENDDO
      ENDDO

      CALL PWAIT

      END SUBROUTINE SIKOR
