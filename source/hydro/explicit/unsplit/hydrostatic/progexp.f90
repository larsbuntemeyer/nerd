!
!
!
SUBROUTINE PROGEXP                                                 &
!
! WARNING! Here are 3 dummy argument names that
! differ from actual parameter names.
! The rest if fine.
!
!  OM850M(1), OM500M(1), OM300M(1),
!     |        |        |
!     v        v        v
(ZOM850M, ZOM500M, ZOM300M,                                     &
 AK    , BK    , A1T    ,                                       &
 A2T   , VVFH  , FIB    , FC     , ACPHIR , CPHI , PS    ,      &
 TGL   , TGW   , TGI    , QDBL   , QDBW   , QDBI , BFLHSL,      &
 BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, TG   , QDB   ,      &
 BFLHS , BFLQDS, BFLUS  , BFLVS  , U      , V    , T     ,      &
 QD    , QW    , FI     , TMKVMH , TMCHL  , TMCHW,              &
 TMCHI , TMCM  , TMCH   , SOTHDT , TTK    , QDTK , UVTK  ,      &
 TTS   , QDTS  , QWTS   , INFRL  , INFRW  , INFRI, QI    ,      &
 QITS  , PINT  , DWDT   , ETAS   , DDMPUJ , DDMPVJ,HYDRODP)    
!
!
!**** PROGEXP  -   UP: PROGNOSE FUER EINEN ZEITSCHRITT
!**   AUFRUF   :   CALL PROGEXP ( JATPROG, JETPROG,
!**               1               ZOM850M, ZOM500M, ZOM300M )
!**   ENTRIES  :      ---
!**   ZWECK    :   BERECHNUNG ALLER PROGNOSEVARIABLEN FUER
!**                DEN ZEITPUNKT NZT+1 (EXPLIZITE PROGNOSE)
!**   VERSIONS-
!**   DATUM    :   22.7.1991
!**
!**   EXTERNALS:   GAUSS  , HQTOTQ , COPYRE
!**
!**   EINGABE-
!**   PARAMETER:   JATPROG: ANFANGS-J-INDEX DES PROGNOSEBEREICHS;
!**                JETPROG: END    -J-INDEX DES PROGNOSEBEREICHS;
!**                         DIE FESTLEGUNG ERFOLGT IN UP *PROGORG*.
!**   AUSGABE-
!**   PARAMETER    ZOM850M: OMEGA-MITTELWERT IN ETWA 850 HPA (KONTROLLE)
!**                ZOM500M: OMEGA-MITTELWERT IN ETWA 500 HPA (KONTROLLE)
!**                ZOM300M: OMEGA-MITTELWERT IN ETWA 300 HPA (KONTROLLE)
!**
!**   COMMON-
!**   BLOECKE  :   PARAM  , ORG, COMDYN, COMPYH, PHYKON, HIGKON, PARKON
!**                PROGCHK, COMDIA, COMPCST, COMPMLF, COMPSLF, COMPGP3
!**
!**   METHODE  :   ZEITLICH:  SEMI-IMPLIZIT/LEAP-FROG
!**                RAEUMLICH: FINITE DIFFERENZEN 2.ORDNUNG
!**                           ARAKAWA-C-GITTER (HORIZONTAL)
!**                           ETA-KOORDINATE (VERTIKAL)
!**   FEHLERBE-
!**   HANDLUNG :      ---
!**   VERFASSER:   G. DOMS UND D. MAJEWSKI
!
USE MO_PARORG,               ONLY:   JE, IEKE, NEIGHBOR
USE MO_ORG,                  ONLY:   LAISTEP, IAA, IEA, NJ, JAH,       &
                                     LPTOP0, IEH, NA, NA2, NE,          &
                                     IAV, IEV, NJ2, IAH, JEH, IEU,      &
                                     PTOP, IAU, JEV, IEHGG, JEHGG,      &
                                     IAHGG, JAHGG                        
USE MO_COMDYN,               ONLY:   ALCNVA, VBMXV, VBCFL,              &
                                     LDIVDAMP, LHDIFF2                   
USE MO_PHYKON,               ONLY:   R, RERD, WCP, WLK, WLS              
USE MO_HIGKON,               ONLY:   ED2DT, EDDPHI, RDDRM1, WCPR,       &
                                     EDDLAM, EDADPHI, DT2, EDG,         &
                                     DTDEH, RDRD, EMRDRD                 
USE MO_PROGCHK,              ONLY:   KFL850, KFL500, KFL300              
USE MO_FAKTINF,              ONLY:   EDFAKINF                            
USE MO_DYNAMICS,             ONLY:   INDEX1, INDEX2, INDEX3,            &
                                     INDEX4, JM2, JM3, JM4,             &
                                     JN5, JO2, JO3, JO4, JO5, JS4,      &
                                     JS5, JSP, JU3, JU4, JM5,           &
                                     JU5, JZM, JZN, JZO, JZS, JZU,      &
                                     ZX2, ZTKVZ, ZALOG2, Z4DRERD,       &
                                     ZTKVL, ZA1A, ZA2A, LMASSF           
USE MO_DYNAMICS,             ONLY:   DYN_ADVECT_DIFFUSE,                &
                                     DYN_DIFFUSE_HORIZONTAL,            &
                                     DYN_INIT_DIFFUSION,                &
                                     DYN_INIT_INDICES,                  &
                                     DYN_SWAP_INDICES,                  &
                                     DYN_PRECOMPUTE
USE MO_MAGNUS,               ONLY:   FGEW, FGQD 
USE MO_MEMORY_DYNAMICS
!
IMPLICIT NONE
!
! WCPR = 1/c_p : inverse of specific heat of dry air (4.7)
! WLK  = L_v   : latent heat of evaporation (4.7)
! T    = h - 1/c_p *
!
 
!-----------------------------------------------------------------------
 
! ERFORDERLICHE EM-FELDER AUS DEM LANGZEITSPEICHER DIMENSIONIEREN
! ---------------------------------------------------------------
! VERTIKAL-KOORDINATEN-PARAMETER
REAL, INTENT(IN) :: AK(KE1), BK(KE1)
 
! VERTIKAL VARIIERENDER IMPLIZITHEITSGRAD DER V-DIFFUSION
REAL, INTENT(IN) :: A1T(KE1), A2T(KE1)
 
! EVT. VERTIKAL VARIIERENDER VERSTAERKUNGSFAKTOR DER H-DIFFUSION
REAL, INTENT(IN) :: VVFH(KE)
 
! EXTERNE PARAMETER
REAL, INTENT(IN) :: FIB(IE,JE), FC(IE,JE), ACPHIR(JE,2),CPHI(JE,2)
 
! PROGNOSTISCHE BODENFELDER
 
REAL, INTENT(OUT) ::  PS(IE,JE,3)
REAL, INTENT(IN)  ::  TGL (IE,JE,3), TGW (IE,JE,3), TGI (IE,JE,3),      &
       QDB(IE,JE,3), QDBL(IE,JE,3), QDBW(IE,JE,3), QDBI(IE,JE,3),       &
       TG (IE,JE,3)                                                      
!DIAGNOSTISCHE BODENFELDER                                                 
REAL, INTENT(OUT) ::   BFLHS  (IE,JE),                                  &
         BFLHSL (IE,JE), BFLHSW (IE,JE), BFLHSI (IE,JE),                &
         BFLQDS (IE,JE),                                                &
         BFLQDSL(IE,JE), BFLQDSW(IE,JE), BFLQDSI(IE,JE),                &
         BFLUS  (IE,JE), BFLVS  (IE,JE)                                  
!                                                                          
INTEGER, INTENT(IN) :: INFRL  (IE,JE), INFRW  (IE,JE),                   &
                       INFRI  (IE,JE)                                     
                                                                           
!ATMOSPHAEREN-FELDER                                                       
REAL, INTENT(OUT) ::   U    (IE,JE,KE,3), V (IE,JE,KE,3),                &
                       T    (IE,JE,KE,3), QD(IE,JE,KE,3),                             &
                       QW   (IE,JE,KE,3)                                               
REAL, INTENT(IN) ::   FI(IE,JE,KE,2)                                      
                                                                           
!FELDER DER KOEFFIZIENTEN UND FLUESSE (PHYSIKALISCHE UPS)                  
REAL, INTENT(IN) :: TMKVMH(IE*(KE-1),JE,2)                                
                                                                          
REAL, INTENT(IN) ::   TMCM(IE,JE), TMCH(IE,JE), TMCHL(IE,JE), TMCHW(IE,JE), TMCHI(IE,JE)           
                                                                           
REAL, INTENT(IN) ::   SOTHDT(IEKE,JE,2)                                   
                                                                           
REAL, INTENT(IN) ::   TTK(IEKE,JE), QDTK(IEKE,JE), UVTK(IEKE,JE,2)        
                                                                           
REAL, INTENT(IN) ::   TTS(IEKE,JE), QDTS(IEKE,JE), QWTS(IEKE,JE)          
                                                                           
REAL, INTENT(IN)  :: DWDT(IE,JE,KE ,3)                                     
REAL, INTENT(OUT) :: PINT(IE,JE,KE1,3)                                    
REAL, INTENT(OUT) :: ETAS(IE,JE,KE1)                                      
REAL, INTENT(IN)  :: DDMPUJ(JE), DDMPVJ(JE)                                
!                                                                          
REAL, INTENT(IN)  :: HYDRODP(IE,JE,KE,3)                                   
REAL, INTENT(OUT) ::   QI(IE,JE,KE,3), QITS  (IEKE,JE)                   
                                                                           
                                                                           
                                                                           
! LOKALE FELDER DIMENSIONIEREN                                             
! MOVED TO MO_MEMORY_DYNAMICS                                              
                                                                           
CHARACTER  YTXT*26                                                        
!                                                                          
!                                                                          
REAL :: Z1, Z2, Z3, Z4, ZA1, ZA2,                                        &
        ZA3, ZAGAM, ZAGAT, ZAGCM, ZAGCT, ZDELTF, ZEPSRAY,                &
        ZFADVX                                                            
REAL :: ZFADVY, ZFDIVO, ZFDIVU, ZFDIVX, ZFKIO, ZFKIU, ZFVZO,             &
        ZOM300M, ZOM500M, ZOM850M, ZQD1, ZQD2, ZQD3,                     &
        ZQD4, ZQW1, ZQW2, ZT1, ZFVZU                                      
REAL :: ZT2, ZT3, ZT4, ZTMKVHM, ZTRC, ZX1, ZXO,                          &
        ZXU, ZZGEW
INTEGER :: I, J, K, JP1, KP1, KM1
REAL :: DDMPU, DDMPV
REAL, EXTERNAL :: GETP
!
INTEGER :: I_TRACER
!
!-----------------------------------------------------------------------
!
! PARAMETER ZUR BERECHNUNG DER DIFFUSION
! UND LN(2) SETZEN.
!
za1a = alcnva
IF ( za1a<0.5 ) za1a = 0.5
za2a = (1.-za1a)
zalog2 = ALOG(2.)
z4drerd = 4.0/rerd
!
! WENN DIE MAXIMALE WINDGESCHWINDIGKEIT 95 % DES STABILITAETS-
! KRITISCHEN WERTES UEBERSCHREITET, WIRD IN DER U- UND V-GLEICHUNG
! RAYLEIGH-REIBUNG EINGEFUEHRT, UM DEN WIND ABZUBREMSEN.
!
IF ( vbmxv>0.95*vbcfl ) THEN
  zepsray = 0.0005*ed2dt*(vbmxv-0.95*vbcfl)/(0.05*vbcfl)
  IF ( vbmxv>vbcfl*1.05 ) THEN
     WRITE (ytxt,'(A,F5.1,A)') 'WARNING: VBMAX=' , vbmxv , ' M/S'
     CALL REMARK(ytxt)
     IF ( vbmxv>250.0 ) THEN
        PRINT * , 'VBMAX EXCEEDS 250 M/S'
        STOP 1
     ENDIF
  ENDIF
ELSE
  zepsray = 0.0
ENDIF
! 
! OMEGA FLAECHENMITTELWERTE VORBEREITEN
!
Zom850m = 0.
Zom500m = 0.
Zom300m = 0.
! 
! SETZEN DES STEUERPARAMETERS FUER TROCKENE (QD=QW=0) BZW. FEUCHTE
! (QD>0, QW>0) PROGNOSE; WENN EIN ADIABATISCHER INITIALISIERUNGS-
! ZEITSCHRITT GERECHNET WIRD, IST ZTRC=0.0, SONST ZTRC=1.0
IF ( laistep ) THEN
  ztrc = 0.0
ELSE
  ztrc = 1.0
ENDIF
! 
!-----------------------------------------------------------------------
! 
!     EXPLIZITE PROGNOSE OHNE RANDRELAXATION
!     DIE PROGNOSE ERFOLGT 'SCHEIBENWEISE' VON J = JATPROG BIS JETPROG.
!     IN EINEM TASK WERDEN JEWEILS ZUSAMMENHAENGENDE BEREICHE VON
!     J = JATPROG BIS JETPROG BEHANDELT; DIE TASK-EINTEILUNG UND
!     STEUERUNG ERFOLGT IM UP *PROGORG*; MAXIMAL KOENNEN VIER TASKS
!     PARALLEL LAUFEN; D.H. DER BEREICH JAH BIS JEH WIRD DANN IN VIER
!     UNTERBEREICHE GETEILT.
 
!***********************************************************************
!*                                                                     *
!*    ZUORDNUNG VON SCHEIBENINDIZES  'SUED'--------->-----------'NORD' *
!*    -----------------------------                                    *
!*                                          J-2 J-1      J      J+1 J+2*
!*    U, T, QD, QW, PHY.UP'S                     +       +       +     *
!*    V                                              +       +         *
!*                                                                     *
!*    ZALOPN(K+1) ALOG( P(K+1/2) )              JU3     JM3     JO3    *
!*    ZPHF(K)                               JS5 JU5     JM5     JO5 JN5*
!*    ZDP (K)                                   JU3     JM3     JO3    *
!*    ZPNF(K)                               JS5 JU5     JM5     JO5 JN5*
!*    ZTV (K)                               JS5 JU5     JM5     JO5 JN5*
!*    ZGQD(K)                                   JU3     JM3     JO3    *
!*    ZALOPH(K)   ALOG( P(K) )                  JU3     JM3     JO3    *
!*    ZBETAK(K)   HILFSGROESSE F. FIH, ALOPN    JU3     JM3     JO3    *
!*    ZFIHF(K)    FI(K) AN HAUPTFLAECHEN    JS5 JU3     JM3     JO3 JN5*
!*                FUER PHYS. PARAMETRISIER.                            *
!*    ZFIH(K)     FI(K) AN HAUPTFLAECHEN        JU3     JM3     JO3    *
!*                FUER DRUCKGRADIENTTERM                               *
!*    ZGU(K)      .5*(DP(I+1)+DP(I))*U                  JM2     JO2    *
!*    ZGV(K)      .5*(DP(J+1)+DP(J))*V              JU3     JM3     JO3*
!*    ZZETA(K)    POTENTIELLE VORTICITY             JM2     JO2        *
!*    ZEKIN(K)    KINETISCHE ENERGIE                    JM2     JO2    *
!*    ZPPHI(K)    R*TV*GRADY(LN(P))                 JM2     JO2        *
!*    ZSDIV(K+1)  VERT. SUMME DER DIVERGENZEN           JM2     JO2    *
!*    ZETAS(K+1)  MODIF. VERTIKALGESCHWINDIGK.          JM2     JO2    *
!*    ZSKHH(K)    DP*SKH/DLAM**2 AM H-GP.               JM2     JO2    *
!*    ZSKHU(K)    DP*SKH/DLAM**2 AM U-GP.               JM2     JO2    *
!*    ZSKHV(K)    DP*SKH/DLAM**2 AM V-GP.           JM2     JO2        *
!*    ZSKHZ(K)    DP*SKH/DLAM**2 AM Z-GP.           JM2     JO2        *
!*    ZHJ(K)      H(NJ)                         JU3     JM3     JO3    *
!*    ZQDWJ(K)    QDW(NJ)                       JU3     JM3     JO3    *
!*    ZTMKVM(K)   TMKVM(K)                      JU3     JM3     JO3    *
!*    ZTMKVH(K)   TMKVH(K)                      JU3     JM3     JO3    *
!*    ZUTK(K)     U-TENDENZ (KONVEKTIV)                 JM2     JO2    *
!*    ZVTK(K)     V-TENDENZ (KONVEKTIV)                 JM2     JO2    *
!*                                                                     *
!*    ZPLAM(K)    R*TV*GRADX(LN(P))                      +             *
!*    ZALPOM(K)   ALPHA*OMEGA                            +             *
!*    ZTPA (K)    HP (NA)                  JS5  JU5     JM5     JO5 JN5*
!*    ZQDWA(K)    QDW(NA)                  JS5  JU5     JM5     JO5 JN5*
!*    ZPLAM(K)    DRUCKGRADIENT                          +             *
!*    SOWIE WEITERE LOKALE FELDER OHNE SCHEI-            +             *
!*    BENINDEX SIND IN DER J-FLAECHE DEFINIERT           +             *
!*                                                                     *
!*                                         J-2  J-1      J      J+1 J+2*
!*                                                                     *
!***********************************************************************
! 
!     VORBESETZUNG LOKALER FELDER AM OBER- UND UNTERRAND
!     --------------------------------------------------
DO i = iaa , iea
  ZPNF(i,1,1) = Pint(i,jah-2,1,nj)
  ZPNF(i,1,2) = Pint(i,jah-1,1,nj)
  ZPNF(i,1,3) = Pint(i,jah,1,nj)
  ZPNF(i,1,4) = Pint(i,jah+1,1,nj)
  ZPNF(i,1,5) = Pint(i,jah+2,1,nj)
  ZSDIV(i,1,1) = 0.
  ZSDIV(i,1,2) = 0.
  ZETAS(i,1,1) = 0.
  ZETAS(i,1,2) = 0.
  ZETAS(i,ke1,1) = 0.
  ZETAS(i,ke1,2) = 0.
  AGA(i,1,1) = 0.
  AGA(i,1,2) = 0.
  AGA(i,1,3) = 0.
  AGC(i,ke,1) = 0.
  AGC(i,ke,2) = 0.
  AGC(i,ke,3) = 0.
  AGA(i,1,4) = 0.
  AGC(i,ke,4) = 0.
ENDDO
! 
!KS introduction of pressure top not equal 0.0
!
IF ( lptop0 ) THEN
  zalopn(iaa:iea,1,1) = 0.0
  zalopn(iaa:iea,1,2) = 0.0
  zalopn(iaa:iea,1,3) = 0.0
ELSE
  zalopn(iaa:iea,1,1) = ALOG(ZPNF(iaa:iea,1,1))
  zalopn(iaa:iea,1,2) = ALOG(ZPNF(iaa:iea,1,2))
  zalopn(iaa:iea,1,3) = ALOG(ZPNF(iaa:iea,1,3))
ENDIF
! 
Etas(iaa:iea,:,1:ke1) = 0.
zdampu(1:ie,1:ke) = 0.
zdampv(1:ie,1:ke) = 0.
zdextu(1:ie) = 0.
zdextv(1:ie) = 0.
! 
! ANFANGSINDIZES FUER ZYKLISCHES UMSPEICHERN SETZEN
CALL DYN_INIT_INDICES
!
! LOKALE FELDER AM 'SUEDRAND' (J=JATPROG, JATPROG-1) VORBESETZEN
! --------------------------------------------------------------
 
! ACHTUNG:
! UM BEIM AUTOTASKING (PARALLER DURCHLAUF VON UP *PROGEXP*) IDENTI-
! SCHE ERGEBNISSE ZU ERZIELEN, MUSS DARAUF GEACHTET WERDEN, DASS
! DIE FORMULIERUNG DER TERME FUER DEN 'SUEDRAND' UND FUER DAS INNE-
! RE DES PROGNOSEGEBIETES EXAKT IDENTISCH SIND. SONST KANN DURCH DIE
! UNVERMEIDLICHEN RUNDUNGSFEHLER DAS ERGEBNIS JE NACH ANZAHL DER
! TASKS VERSCHIEDEN SEIN.
! 
jzs = jah - 2
jzu = jah - 1
jzm = jah
jzo = jah + 1
! 
! GEPACKTE FELDER DER PHYSIKALISCHEN TENDENZEN UND KOEFFIZIENTEN
! AUSPACKEN
CALL COPYRE(Tmkvmh(1,jzu,1),ZTMKVM(1,2,ju3),ie*(ke-1))
CALL COPYRE(Tmkvmh(1,jzu,2),ZTMKVH(1,2,ju3),ie*(ke-1))
CALL COPYRE(Tmkvmh(1,jzm,1),ZTMKVM(1,2,jm3),ie*(ke-1))
CALL COPYRE(Tmkvmh(1,jzm,2),ZTMKVH(1,2,jm3),ie*(ke-1))
CALL COPYRE(Uvtk(1,jzm,1),ZUTK(1,1,jm2),ieke)
CALL COPYRE(Uvtk(1,jzm,2),ZVTK(1,1,jm2),ieke)
! 
CALL DYN_PRECOMPUTE
!
! LOKALE HILFSFELDER FUER HORIZONTALDIFFUSION
CALL DYN_INIT_DIFFUSION(j)
!
! RAENDER DER SCHEIBENFELDER FUER DIE HQ-ZERLEGUNG VORBESETZEN;
! (BELIEBIGE KONSISTENTE WERTE)
!VDIR NOVECTOR
DO i = 1 , 2
  DO k = 1 , ke
     HE(i,k) = wcp*T(1,jzu,k,nj) + wlk*Qd(1,jzu,k,nj)
     QDWE(i,k) = Qd(1,jzu,k,nj) + Qw(1,jzu,k,nj) + Qi(1,jzu,k,nj)
     TSTART(i,k) = T(1,jzu,k,nj)
     GQDSTA(i,k) = ZGQD(1,k,ju3)
     PHFSTA(i,k) = ZPHF(1,k,ju5)
     PHFE(i,k) = ZPHF(1,k,ju5)
  ENDDO
ENDDO
 
!VDIR NOVECTOR
DO i = ie - 1 , ie
  DO k = 1 , ke
     HE(i,k) = HE(1,k)
     QDWE(i,k) = QDWE(1,k)
     TSTART(i,k) = TSTART(1,k)
     GQDSTA(i,k) = GQDSTA(1,k)
     PHFSTA(i,k) = PHFSTA(1,k)
     PHFE(i,k) = PHFE(1,k)
  ENDDO
ENDDO
!---------------------------------------------------------------------------- 
! BEGINN DER SCHLEIFE IN J-RICHTUNG
! ---------------------------------
! ACHTUNG:
! UM BEIM AUTOTASKING (PARALLER DURCHLAUF VON UP *PROGEXP*) IDENTI-
! SCHE ERGEBNISSE ZU ERZIELEN, MUSS DARAUF GEACHTET WERDEN, DASS
! DIE FORMULIERUNG DER TERME FUER DEN 'SUEDRAND' UND FUER DAS INNE-
! RE DES PROGNOSEGEBIETES EXAKT IDENTISCH SIND. SONST KANN DURCH DIE
! UNVERMEIDLICHEN RUNDUNGSFEHLER DAS ERGEBNIS JE NACH ANZAHL DER
! TASKS VERSCHIEDEN SEIN.
DO j = jah , jeh
 
  jzn = j + 2
 
  ! LOKALE FELDER FUER J   BESETZEN
  ! GEPACKTE FELDER DER PHYSIKALISCHEN TENDENZEN UND KOEFFIZIENTEN
  ! AUSPACKEN
  CALL COPYRE(Ttk(1,j),zttk,ieke)
  CALL COPYRE(Qdtk(1,j),zqdtk,ieke)
  CALL COPYRE(Tts(1,j),ztts,ieke)
  CALL COPYRE(Qdts(1,j),zqdts,ieke)
  CALL COPYRE(Qwts(1,j),zqwts,ieke)
  CALL COPYRE(Sothdt(1,j,1),zsodta,ieke)
  CALL COPYRE(Sothdt(1,j,2),zthdta,ieke)
  CALL COPYRE(Qits(1,j),zqits,ieke)
  ! LOKALE FELDER FUER J+1 BESETZEN
  ! GEPACKTE FELDER DER PHYSIKALISCHEN TENDENZEN UND KOEFFIZIENTEN
  ! AUSPACKEN
  CALL COPYRE(Tmkvmh(1,j+1,1),ZTMKVM(1,2,jo3),ie*(ke-1))
  CALL COPYRE(Tmkvmh(1,j+1,2),ZTMKVH(1,2,jo3),ie*(ke-1))
  CALL COPYRE(Uvtk(1,j+1,1),ZUTK(1,1,jo2),ieke)
  CALL COPYRE(Uvtk(1,j+1,2),ZVTK(1,1,jo2),ieke)
  ! LOKALE FELDER FUER J+2 BESETZEN
  zx1 = .5*r*Acphir(j,1)*eddlam
  zx2 = .5*rerd*(Cphi(j,1)+Cphi(j+1,1))
  zxo = Cphi(j+1,1)*eddphi
  zxu = Cphi(j,1)*eddphi
  ! 
  DO k = 1 , ke
     kp1 = k + 1
     DO i = iaa , iea
        jp1 = j + 1
        ZDP(i,k,jo3) = Hydrodp(i,jp1,k,nj)
        !LM ADDED ZDPINT FOR OMEGA-ALPHA TERM
        ZDPINT(i,k,jo3) = Pint(i,jp1,kp1,nj) - Pint(i,jp1,k,nj)
        ! FUER GEWICHTETE HDIFF
        !==============================================================
        ! !!!! Should ZDP5 be based on PINT?
        !==============================================================
        ZDP5(i,k,jn5) = Hydrodp(i,jzn,k,na)
        ZDPV(i,k,jo4) = 0.5*(ZDP5(i,k,jo5)+ZDP5(i,k,jn5))
        ! pressure at half levels
        ZPHF(i,k,jn5) = 0.5*(Pint(i,jzn,k,nj)+Pint(i,jzn,kp1,nj))
        ZPNF(i,kp1,jn5) = Pint(i,jzn,kp1,nj)
        zalopn(i,kp1,jo3) = ALOG(ZPNF(i,kp1,jo5))
        ZTV(i,k,jn5) = T(i,jzn,k,nj)*(1.0+rddrm1*Qd(i,jzn,k,nj)-(Qw(i,jzn,k,nj)+Qi(i,jzn,k,nj)))
        zzgew = FGEW(T(i,jp1,k,nj))
        ZGQD(i,k,jo3) = FGQD(zzgew,ZPHF(i,k,jo5))
     ENDDO
  ENDDO
  !
  ! FUER GEWICHTETE HDIFF
  !
  DO k = 1 , ke
     ! FEHLERKORREKTUR FUER ZDPU
     DO i = iaa , ieh + 1
        ZDPU(i,k,jo3) = 0.5*(ZDP5(i,k,jo5)+ZDP5(i+1,k,jo5))
     ENDDO
  ENDDO
  !
  !KS SPLITTED LOOP TO REMOVE IF INSIDE LOOP
  !
  !KS FIRST PART K<KE
  DO k = 1 , ke - 1
     DO i = iaa , iea
        ZFIHF(i,k,jn5) = Fi(i,jzn,k+1,na2) + r*ZTV(i,k,jn5)*ALOG(ZPNF(i,k+1,jn5)/ZPHF(i,k,jn5))    &
                       & /Dwdt(i,jzn,k,nj)
        ! ZTPA (I,JN5,K)    = T(I,JZN,K,NA) + ZFIHF(I,K,JN5)*WCPR
        ZTPA(i,jzn,k) = T(i,jzn,k,na) + ZFIHF(i,k,jn5)*wcpr
     ENDDO
  ENDDO
  !KS      SECOND PART K=KE
  DO i = iaa , iea
     ZFIHF(i,ke,jn5) = Fib(i,jzn) + r*ZTV(i,ke,jn5)*ALOG(ZPNF(i,ke+1,jn5)/ZPHF(i,ke,jn5))          &
                     & /Dwdt(i,jzn,ke,nj)
     ! ZTPA (I,JN5,K)    = T(I,JZN,KE,NA) + ZFIHF(I,KE,JN5)*WCPR
     ZTPA(i,jzn,k) = T(i,jzn,ke,na) + ZFIHF(i,ke,jn5)*wcpr
  ENDDO
  !
  !KS introduction of pressure top not equal 0.0
  IF ( lptop0 ) THEN
     zbetak(iaa:iea,1,jo3) = zalog2
  ELSE
     DO i = iaa , iea
        zbetak(i,1,jo3) = 1. - ZPNF(i,1,jo5)/ZDP(i,1,jo3)*(zalopn(i,2,jo3)-zalopn(i,1,jo3))
     ENDDO
  ENDIF
  !
  !KS SPLITTED LOOP TO REMOVE IF INSIDE LOOP
  !
  !KS FIRST PART K=1
  DO i = iaa , iea
     ZALOPH(i,1,jo3) = zalopn(i,2,jo3) - zbetak(i,1,jo3)
     ZFIH(i,1,jo3) = Fi(i,j+1,2,nj2) + r*ZTV(i,1,jo5)*zbetak(i,1,jo3)/Dwdt(i,j+1,1,nj)
  ENDDO
  !KS SECOND PART K>1
  DO k = 2 , ke - 1
     DO i = iaa , iea
        zbetak(i,k,jo3) = 1. - ZPNF(i,k,jo5)/ZDP(i,k,jo3)*(zalopn(i,k+1,jo3)-zalopn(i,k,jo3))
        ZFIH(i,k,jo3) = Fi(i,j+1,k+1,nj2) + r*ZTV(i,k,jo5)*zbetak(i,k,jo3)/Dwdt(i,j+1,k,nj)
        ZALOPH(i,k,jo3) = zalopn(i,k+1,jo3) - zbetak(i,k,jo3)
     ENDDO
  ENDDO
  !
  ! CORIOLISPARAMETER AM ZETA-PUNKT BEREITSTELLEN
  !
  DO i = iaa , ieh + 1
     ZFCZET(i) = 0.25*(Fc(i,j)+Fc(i+1,j)+Fc(i,j+1)+Fc(i+1,j+1))
  ENDDO
 
  DO i = iaa , iea
     zbetak(i,ke,jo3) = 1. - ZPNF(i,ke,jo5)/ZDP(i,ke,jo3)*(zalopn(i,ke1,jo3)-zalopn(i,ke,jo3))
     ZALOPH(i,ke,jo3) = zalopn(i,ke1,jo3) - zbetak(i,ke,jo3)
     ZFIH(i,ke,jo3) = Fib(i,j+1) + r*ZTV(i,ke,jo5)*zbetak(i,ke,jo3)/Dwdt(i,j+1,ke,nj)
  ENDDO
  ! 
  ! WEITERE LOKALE FELDER AN U/ZETA-GITTERPUNKTEN (HF) BESETZEN , JO
  !
  DO k = 1 , ke
     DO i = iaa , ieh + 1
        ZGU(i,k,jo2) = .5*(ZDP(i,k,jo3)+ZDP(i+1,k,jo3))*U(i,j+1,k,nj)
        ZPLAM(i,k) = zx1*(ZTV(i+1,k,jm5)+ZTV(i,k,jm5))*(ZALOPH(i+1,k,jm3)-ZALOPH(i,k,jm3))
        z1 = zx2*ZFCZET(i)
        z2 = eddlam*(V(i+1,j,k,nj)-V(i,j,k,nj))
        z3 = zxo*U(i,j+1,k,nj) - zxu*U(i,j,k,nj)
        z4 = (ZDP(i,k,jm3)+ZDP(i+1,k,jm3))*Cphi(j,1) + (ZDP(i,k,jo3)+ZDP(i+1,k,jo3))*Cphi(j+1,1)
        ZZETA(i,k,jo2) = z4drerd*(z1+z2-z3)/z4
     ENDDO
  ENDDO
  !
  ! WEITERE LOKALE FELDER AN H/V-GITTERPUNKTEN (HF) BESETZEN
  ! 
  zx1 = 0.5*r*edadphi
  IF ( (j<jeh) .OR. (NEIGHBOR(2)/=-1) ) THEN
     zfkio = Cphi(j+1,2)/Cphi(j+1,1)
     zfkiu = Cphi(j,2)/Cphi(j+1,1)
     zfdivx = Acphir(j+1,1)*eddlam
     zfdivo = Acphir(j+1,1)*eddphi*Cphi(j+1,2)
     zfdivu = Acphir(j+1,1)*eddphi*Cphi(j,2)
 
     DO k = 1 , ke
        DO i = iah - 1 , ieh + 1
           ZGV(i,k,jo3) = .5*(ZDP(i,k,jo3)+Hydrodp(i,j+2,k,nj))*V(i,j+1,k,nj)
           ZPPHI(i,k,jo2) = zx1*(ZTV(i,k,jo5)+ZTV(i,k,jm5))*(ZALOPH(i,k,jo3)-ZALOPH(i,k,jm3))
           ZEKIN(i,k,jo2) = .25*(U(i-1,j+1,k,nj)**2+U(i,j+1,k,nj)**2+zfkiu*V(i,j,k,nj)             &
                          & **2+zfkio*V(i,j+1,k,nj)**2)
        ENDDO
     ENDDO
 
     DO k = 1 , ke
        DO i = iah - 1 , ieh + 1
           ZSDIV(i,k+1,jo2) = ZSDIV(i,k,jo2) + zfdivx*(ZGU(i,k,jo2)-ZGU(i-1,k,jo2))                &
                            & + zfdivo*ZGV(i,k,jo3) - zfdivu*ZGV(i,k,jm3)
        ENDDO
     ENDDO
     !
     ! DIVERGENCE DAMPING
     ! ------------------
     !
     IF ( ldivdamp ) THEN
        ddmpu = Ddmpuj(j)
        ddmpv = Ddmpvj(j)
        !
        ! HORIZONTAL MODE
        !
        DO k = 1 , ke
           DO i = iah - 1 , ieu
              zdampu(i,k) = ((ZSDIV(i+1,k+1,jm2)-ZSDIV(i,k+1,jm2))                                 &
                          & -(ZSDIV(i+1,k,jm2)-ZSDIV(i,k,jm2)))*2./(ZDP(i,k,jm3)+ZDP(i+1,k,jm3))   &
                          & *ddmpu
              zdampv(i,k) = ((ZSDIV(i,k+1,jo2)-ZSDIV(i,k+1,jm2))-(ZSDIV(i,k,jo2)-ZSDIV(i,k,jm2)))  &
                          & *2./(ZDP(i,k,jo3)+ZDP(i,k,jm3))*ddmpv
           ENDDO
        ENDDO
        !
        ! EXTERNAL MODE
        !
        DO i = iah - 1 , ieu
           zdextu(i) = (ZSDIV(i+1,ke1,jm2)-ZSDIV(i,ke1,jm2))*2./(Ps(i,j,nj)+Ps(i+1,j,nj))*ddmpu
           zdextv(i) = (ZSDIV(i,ke1,jo2)-ZSDIV(i,ke1,jm2))*2./(Ps(i,j+1,nj)+Ps(i,j,nj))*ddmpv
        ENDDO
     ENDIF
     ! 
     !
     !
     DO k = 2 , ke
        DO i = iah - 1 , ieh + 1
           ZETAS(i,k,jo2) = Bk(k)*ZSDIV(i,ke1,jo2) - ZSDIV(i,k,jo2)
        ENDDO
        Etas(iah-1:ieh+1,j+1,k) = ZETAS(iah-1:ieh+1,k,jo2)
     ENDDO
     !
  ELSE
     DO k = 1 , ke
        DO i = iah , ieh
           ZPPHI(i,k,jo2) = zx1*(ZTV(i,k,jo5)+ZTV(i,k,jm5))*(ZALOPH(i,k,jo3)-ZALOPH(i,k,jm3))
        ENDDO
     ENDDO
  ENDIF
  ! 
  !CCCLBC  Compute explizit horizontal diffusion
  CALL DYN_DIFFUSE_HORIZONTAL(j)
  !CCCLB
  !
  ! EXPLIZITE PROGNOSE OHNE RANDRELAXATIONSTERME
  ! --------------------------------------------
  !
  zx2 = rerd*Acphir(j,1)
  zfadvx = .5*Acphir(j,1)*eddlam
  zfadvy = .5*Acphir(j,1)*eddphi
  ! 
  ! 1. BODENDRUCK + DRUCK
  ! ---------------------
  DO i = iah , ieh
     PSDT(i) = ZSDIV(i,ke1,jm2)
     Ps(i,j,ne) = Ps(i,j,na) - PSDT(i)*dt2
     Pint(i,j,1,ne) = GETP(Ak(1),Bk(1),Ps(i,j,ne),ptop)
  ENDDO
  ! 
  DO k = 2 , ke + 1
     km1 = k - 1
     DO i = iah , ieh
        Pint(i,j,k,ne) = -(Bk(km1)+Bk(k))*dt2*PSDT(i)*Dwdt(i,j,km1,nj) + Pint(i,j,km1,na)          &
                       & + Pint(i,j,k,na) - Pint(i,j,km1,ne)
     ENDDO
  ENDDO
  ! 
  ! 2. H - QDW - PROGNOSE
  ! ---------------------
  !KS
  ! Coefficients of the tridiagonal system (Eq. 5.3.17)
  !
  ! AGA_K (PSI^t+dt)_k-1 +
  ! AGB_K (PSI^t+dt)_k   +
  ! AGC_K (PSI^t+dt)_k+1 =
  ! AGD_K
  !
  ! PSI stands for T, QW, QD or QI depending on the last index
  ! of AGA, AGB, AGC and AGD
  !KS
  !
  ! HORIZONTALADVEKTION UND QUELLTERME
  !
  CALL DYN_ADVECT_DIFFUSE(j)
  ! 
  DO k = 1 , ke
     DO i = iah , ieh
        HE(i,k) = AGE(i,k,1)*wcp + wlk*AGE(i,k,2)*ztrc
        QDWE(i,k) = (AGE(i,k,2)+AGE(i,k,3)+AGE(i,k,4))*ztrc
        QDLE(i,k) = (AGE(i,k,2)+AGE(i,k,3))*ztrc
        QDIE(i,k) = AGE(i,k,4)*ztrc
        !KS  TSTART(I,K) = T(I,J,K,NJ)
        !  GQDSTA(I,K) = ZGQD(I,K,JM3)
        !  PHFSTA(I,K) = ZPHF(I,K,JM5)
        !KS  PHFE  (I,K) = AKH(K) + BKH(K)*PS(I,J,NE)
     ENDDO
  ENDDO
  !
  ! BERECHNUNG DES FLUESSES LATENTER WAERME AM BODEN;
  ! DIESER FLUESS WIRD AUFSUMMIERT
  !
  DO i = iah , ieh
     IF ( Infrl(i,j)>0 ) Bflqdsl(i,j) = Bflqdsl(i,j) - Tmchl(i,j)                                  &
                                      & *edg*(wlk*(A2t(ke1)*(Qdbl(i,j,na)-(Qd(i,j,ke,na)           &
                                      & +(Qw(i,j,ke,na))))+A1t(ke1)*(Qdbl(i,j,ne)-QDLE(i,ke)))     &
                                      & -wls*(A2t(ke1)*Qi(i,j,ke,na)+A1t(ke1)*QDIE(i,ke)))*dtdeh
     IF ( Infrw(i,j)>0 ) Bflqdsw(i,j) = Bflqdsw(i,j) - Tmchw(i,j)                                  &
                                      & *edg*(wlk*(A2t(ke1)*(Qdbw(i,j,na)-(Qd(i,j,ke,na)           &
                                      & +(Qw(i,j,ke,na))))+A1t(ke1)*(Qdbw(i,j,ne)-QDLE(i,ke)))     &
                                      & -wls*(A2t(ke1)*Qi(i,j,ke,na)+A1t(ke1)*QDIE(i,ke)))*dtdeh
     IF ( Infri(i,j)>0 ) Bflqdsi(i,j) = Bflqdsi(i,j) - Tmchi(i,j)                                  &
                                      & *edg*(wlk*(A2t(ke1)*(Qdbi(i,j,na)-(Qd(i,j,ke,na)           &
                                      & +(Qw(i,j,ke,na))))+A1t(ke1)*(Qdbi(i,j,ne)-QDLE(i,ke)))     &
                                      & -wls*(A2t(ke1)*Qi(i,j,ke,na)+A1t(ke1)*QDIE(i,ke)))*dtdeh
     Bflqds(i,j) = (FLOAT(Infrl(i,j))*Bflqdsl(i,j)+FLOAT(Infrw(i,j))*Bflqdsw(i,j)+FLOAT(Infri(i,j))&
                 & *Bflqdsi(i,j))*edfakinf
  ENDDO
  !
  ! MASSENFLUSS-KORREKTURSCHEMA  (FUER QDW<0)
  !
  IF ( lmassf ) THEN
     DO i = iah , ieh
        ZQKOR(i) = 0.
     ENDDO
     DO k = 1 , ke
        DO i = iah , ieh
           ZQKOR(i) = QDWE(i,k) + ZQKOR(i)*ZEDDPQ(i,k)
           IF ( ZQKOR(i)<0. ) THEN
              QDWE(i,k) = 0.
              ZQKOR(i) = ZQKOR(i)*ZDP(i,k,jm3)
           ELSE
              QDWE(i,k) = ZQKOR(i)
              ZQKOR(i) = 0.
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  ! 
  ! BEI ADIABATISCHER INITIALISIERUNG: T MIT UNVERAENDERTEM QD BERECH-
  ! NEN; SONST IST DIE SKALIGE KONDENSATIONSRATE UNGLEICH NULL.
  !
  IF ( laistep ) THEN
     DO k = 1 , ke
        DO i = iah , ieh
           Qd(i,j,k,ne) = Qd(i,j,k,na)
           Qw(i,j,k,ne) = Qw(i,j,k,na)
           T(i,j,k,ne) = wcpr*HE(i,k)
           Qi(i,j,k,ne) = Qi(i,j,k,na)
        ENDDO
     ENDDO
 
  ELSE
     DO k = 1 , ke
        DO i = iah , ieh
           Qd(i,j,k,ne) = AGE(i,k,2)
           Qw(i,j,k,ne) = AGE(i,k,3)
           T(i,j,k,ne) = AGE(i,k,1)
           Qi(i,j,k,ne) = AGE(i,k,4)
        ENDDO
     ENDDO
     DO i = iah , ieh
        ZQKOR(i) = 0.
     ENDDO
     DO k = 1 , ke
        DO i = iah , ieh
           zqw1 = Qw(i,j,k,ne)        ! water content
           zqw2 = Qi(i,j,k,ne)        ! ice content
           ! This is described in Part I, Section 4.8
           IF ( zqw1<0. ) THEN         ! 4.8.5
                     ! WCPR = 1/c_p : inverse of specific heat of dry air (4.7)
                     ! WLK  = L_v   : latent heat of evaporation (4.7)
                     ! T = h - 1/c_p *
              T(i,j,k,ne) = T(i,j,k,ne) - wlk*Qw(i,j,k,ne)*wcpr
              Qd(i,j,k,ne) = Qd(i,j,k,ne) + Qw(i,j,k,ne)
              Qw(i,j,k,ne) = 0.          ! outside clouds
           ENDIF
           IF ( zqw2<0. ) THEN
              T(i,j,k,ne) = T(i,j,k,ne) - wls*Qi(i,j,k,ne)*wcpr
              Qd(i,j,k,ne) = Qd(i,j,k,ne) + Qi(i,j,k,ne)
              Qi(i,j,k,ne) = 0.
           ENDIF
           ZQKOR(i) = Qd(i,j,k,ne) + ZQKOR(i)*ZEDDPQ(i,k)
           IF ( ZQKOR(i)<0. ) THEN
              Qd(i,j,k,ne) = 0.
              Qw(i,j,k,ne) = 0.
              Qi(i,j,k,ne) = 0.
              ZQKOR(i) = ZQKOR(i)*ZDP(i,k,jm3)
           ELSE
              Qd(i,j,k,ne) = ZQKOR(i)
              ZQKOR(i) = 0.
           ENDIF
        ENDDO
     ENDDO
     !
     ! BERECHNUNG DES FLUESSES SENSIBLER WAERME AM BODEN;
     !
     DO i = iah , ieh
        IF ( Infrl(i,j)>0 ) Bflhsl(i,j) = Bflhsl(i,j) - Tmchl(i,j)                                 &
           & *edg*wcp*(A2t(ke1)*(Tgl(i,j,na)-T(i,j,ke,na))+A1t(ke1)*(Tgl(i,j,ne)-T(i,j,ke,ne))     &
           & +wcpr*(Fib(i,j)-ZFIHF(i,ke,jm5)))*dtdeh
        IF ( Infrw(i,j)>0 ) Bflhsw(i,j) = Bflhsw(i,j) - Tmchw(i,j)                                 &
           & *edg*wcp*(A2t(ke1)*(Tgw(i,j,na)-T(i,j,ke,na))+A1t(ke1)*(Tgw(i,j,ne)-T(i,j,ke,ne))     &
           & +wcpr*(Fib(i,j)-ZFIHF(i,ke,jm5)))*dtdeh
        IF ( Infri(i,j)>0 ) Bflhsi(i,j) = Bflhsi(i,j) - Tmchi(i,j)                                 &
           & *edg*wcp*(A2t(ke1)*(Tgi(i,j,na)-T(i,j,ke,na))+A1t(ke1)*(Tgi(i,j,ne)-T(i,j,ke,ne))     &
           & +wcpr*(Fib(i,j)-ZFIHF(i,ke,jm5)))*dtdeh
        Bflhs(i,j) = (FLOAT(Infrl(i,j))*Bflhsl(i,j)+FLOAT(Infrw(i,j))*Bflhsw(i,j)+FLOAT(Infri(i,j))&
                   & *Bflhsi(i,j))*edfakinf
     ENDDO
     !
     ! OMEGA-WERTE BESTIMMTER MODELLFLAECHEN SUMMIEREN (PROGCHK)
     !
     DO i = iah , ieh
        Zom850m = Zom850m + ABS(ZALPOM(i,kfl850)*ZPHF(i,kfl850,jm5)/(r*ZTV(i,kfl850,jm5)))
        Zom500m = Zom500m + ABS(ZALPOM(i,kfl500)*ZPHF(i,kfl500,jm5)/(r*ZTV(i,kfl500,jm5)))
        Zom300m = Zom300m + ABS(ZALPOM(i,kfl300)*ZPHF(i,kfl300,jm5)/(r*ZTV(i,kfl300,jm5)))
     ENDDO
  ENDIF        ! LAISTEP
  !
  ! 3. U - PROGNOSE
  ! ---------------
  ! 
  zx1 = Acphir(j,1)*eddlam
  zfvzo = 0.125*Cphi(j,2)/Cphi(j,1)
  zfvzu = 0.125*Cphi(j-1,2)/Cphi(j,1)
 
  DO k = 1 , ke
     DO i = iau , ieu
        ZGRAD(i,k) = -ZPLAM(i,k)                                                                   &
                   & - zx1*((ZFIH(i+1,k,jm3)-ZFIH(i,k,jm3))*0.5*(Dwdt(i,j,k,nj)+Dwdt(i+1,j,k,nj))  &
                   & +ZEKIN(i+1,k,jm2)-ZEKIN(i,k,jm2))
        ZVZEQ(i,k) = (ZZETA(i,k,jo2)+ZZETA(i,k,jm2))                                               &
                   & *(zfvzo*(ZGV(i+1,k,jm3)+ZGV(i,k,jm3))+zfvzu*(ZGV(i+1,k,ju3)+ZGV(i,k,ju3)))
        ZEDDPQ(i,k) = 2./(ZDP(i+1,k,jm3)+ZDP(i,k,jm3))
        AGB(i,k,1) = ed2dt
        AGD(i,k,1) = (ed2dt-zepsray)*U(i,j,k,na) + ZGRAD(i,k) + ZVZEQ(i,k) + ZUDIFH(i,k)           &
                   & + 0.5*(ZUTK(i+1,k,jm2)+ZUTK(i,k,jm2)) + zdampu(i,k) + zdextu(i)
     ENDDO
  ENDDO
  ! 
  ! VERTIKALDIFFUSION UND VERTIKALADVEKTION
  !
  DO i = iau , ieu
     zagcm = .25*(ZETAS(i,2,jm2)+ZETAS(i+1,2,jm2))*ZEDDPQ(i,1)
     zagct = -0.5*(ZTMKVM(i,2,jm3)+ZTMKVM(i+1,2,jm3))*ZEDDPQ(i,1)
     AGC(i,1,1) = zagcm*za1a + zagct*A1t(2)
     AGB(i,1,1) = AGB(i,1,1) - AGC(i,1,1)
     AGD(i,1,1) = AGD(i,1,1) - (za2a*zagcm+A2t(2)*zagct)*(U(i,j,2,na)-U(i,j,1,na))
  ENDDO
  ! 
  DO k = 2 , ke - 1
     DO i = iau , ieu
        zagam = -0.25*(ZETAS(i,k,jm2)+ZETAS(i+1,k,jm2))*ZEDDPQ(i,k)
        zagcm = 0.25*(ZETAS(i,k+1,jm2)+ZETAS(i+1,k+1,jm2))*ZEDDPQ(i,k)
        zagat = -0.5*(ZTMKVM(i,k,jm3)+ZTMKVM(i+1,k,jm3))*ZEDDPQ(i,k)
        zagct = -0.5*(ZTMKVM(i,k+1,jm3)+ZTMKVM(i+1,k+1,jm3))*ZEDDPQ(i,k)
        AGA(i,k,1) = zagam*za1a + zagat*A1t(k)
        AGC(i,k,1) = zagcm*za1a + zagct*A1t(k+1)
        AGB(i,k,1) = AGB(i,k,1) - AGA(i,k,1) - AGC(i,k,1)
        AGD(i,k,1) = AGD(i,k,1) - (za2a*zagam+A2t(k)*zagat)*(U(i,j,k-1,na)-U(i,j,k,na))            &
                   & - (za2a*zagcm+A2t(k+1)*zagct)*(U(i,j,k+1,na)-U(i,j,k,na))
     ENDDO
  ENDDO
  !
  DO i = iau , ieu
     zagam = -0.25*(ZETAS(i,ke,jm2)+ZETAS(i+1,ke,jm2))*ZEDDPQ(i,ke)
     zagat = -0.5*(ZTMKVM(i,ke,jm3)+ZTMKVM(i+1,ke,jm3))*ZEDDPQ(i,ke)
     zagct = -0.5*(Tmcm(i,j)+Tmcm(i+1,j))*ZEDDPQ(i,ke)
     AGA(i,ke,1) = zagam*za1a + zagat*A1t(ke)
     AGB(i,ke,1) = AGB(i,ke,1) - AGA(i,ke,1) - zagct*A1t(ke1)
     AGD(i,ke,1) = AGD(i,ke,1) - (za2a*zagam+A2t(ke)*zagat)*(U(i,j,ke-1,na)-U(i,j,ke,na))          &
                 & + A2t(ke1)*zagct*U(i,j,ke,na)
  ENDDO
  ! 
  ! GAUSS-ELIMINATION UND U-ENDE-WERTE AUF U(I,J,K,NE) ABLEGEN
  !
  CALL GAUSS(iau,ieu,index1,AGA,AGB,AGC,AGD,AGE)
  !
  DO k = 1 , ke
     DO i = iau , ieu
        U(i,j,k,ne) = AGE(i,k,1)
     ENDDO
  ENDDO
  !
  ! BERECHNUNG DES U-IMPULSFLUSSES AM BODEN;
  ! DIESER FLUESS WIRD AUFSUMMIERT
  DO i = iau , ieu
     Bflus(i,j) = Bflus(i,j) + 0.5*(Tmcm(i,j)+Tmcm(i+1,j))                                         &
                & *edg*(A2t(ke1)*U(i,j,ke,na)+A1t(ke1)*U(i,j,ke,ne))*dtdeh
  ENDDO
  !
  ! 4. V - PROGNOSE
  ! ---------------
  !
  IF ( j<=jev ) THEN
 
     DO k = 1 , ke
        DO i = iav , iev
           ZGRAD(i,k) = -ZPPHI(i,k,jo2) - edadphi*((ZFIH(i,k,jo3)-ZFIH(i,k,jm3))                   &
                      & *0.5*(Dwdt(i,j,k,nj)+Dwdt(i,j+1,k,nj))+ZEKIN(i,k,jo2)-ZEKIN(i,k,jm2))
           ZUZEQ(i,k) = -0.125*(ZZETA(i-1,k,jo2)+ZZETA(i,k,jo2))                                   &
                      & *(ZGU(i-1,k,jo2)+ZGU(i,k,jo2)+ZGU(i-1,k,jm2)+ZGU(i,k,jm2))
           ZEDDPQ(i,k) = 2./(ZDP(i,k,jo3)+ZDP(i,k,jm3))
           AGB(i,k,2) = ed2dt
           AGD(i,k,2) = (ed2dt-zepsray)*V(i,j,k,na) + ZGRAD(i,k) + ZUZEQ(i,k) + ZVDIFH(i,k)        &
                      & + 0.5*(ZVTK(i,k,jo2)+ZVTK(i,k,jm2)) + zdampv(i,k) + zdextv(i)
        ENDDO
     ENDDO
     !
     ! VERTIKALDIFFUSION UND VERTIKALADVEKTION
     !
     DO i = iav , iev
        zagcm = 0.25*(ZETAS(i,2,jm2)+ZETAS(i,2,jo2))*ZEDDPQ(i,1)
        zagct = -0.5*(ZTMKVM(i,2,jm3)+ZTMKVM(i,2,jo3))*ZEDDPQ(i,1)
        AGC(i,1,2) = zagcm*za1a + zagct*A1t(2)
        AGB(i,1,2) = AGB(i,1,2) - AGC(i,1,2)
        AGD(i,1,2) = AGD(i,1,2) - (za2a*zagcm+A2t(2)*zagct)*(V(i,j,2,na)-V(i,j,1,na))
     ENDDO
 
     DO k = 2 , ke - 1
        DO i = iav , iev
           zagam = -0.25*(ZETAS(i,k,jm2)+ZETAS(i,k,jo2))*ZEDDPQ(i,k)
           zagcm = 0.25*(ZETAS(i,k+1,jm2)+ZETAS(i,k+1,jo2))*ZEDDPQ(i,k)
           zagat = -0.5*(ZTMKVM(i,k,jm3)+ZTMKVM(i,k,jo3))*ZEDDPQ(i,k)
           zagct = -0.5*(ZTMKVM(i,k+1,jm3)+ZTMKVM(i,k+1,jo3))*ZEDDPQ(i,k)
           AGA(i,k,2) = zagam*za1a + zagat*A1t(k)
           AGC(i,k,2) = zagcm*za1a + zagct*A1t(k+1)
           AGB(i,k,2) = AGB(i,k,2) - AGA(i,k,2) - AGC(i,k,2)
           AGD(i,k,2) = AGD(i,k,2) - (za2a*zagam+A2t(k)*zagat)*(V(i,j,k-1,na)-V(i,j,k,na))         &
                      & - (za2a*zagcm+A2t(k+1)*zagct)*(V(i,j,k+1,na)-V(i,j,k,na))
        ENDDO
     ENDDO
 
     DO i = iav , iev
        zagam = -0.25*(ZETAS(i,ke,jm2)+ZETAS(i,ke,jo2))*ZEDDPQ(i,ke)
        zagat = -0.5*(ZTMKVM(i,ke,jm3)+ZTMKVM(i,ke,jo3))*ZEDDPQ(i,ke)
        zagct = -0.5*(Tmcm(i,j)+Tmcm(i,j+1))*ZEDDPQ(i,ke)
        AGA(i,ke,2) = zagam*za1a + zagat*A1t(ke)
        AGB(i,ke,2) = AGB(i,ke,2) - AGA(i,ke,2) - zagct*A1t(ke1)
        AGD(i,ke,2) = AGD(i,ke,2) - (za2a*zagam+A2t(ke)*zagat)*(V(i,j,ke-1,na)-V(i,j,ke,na))       &
                    & + A2t(ke1)*zagct*V(i,j,ke,na)
     ENDDO
     ! 
     ! GAUSS ELIMINATION UND V-ENDE-WERTE AUF V(I,J,K,NE) ABLEGEN
     !
     CALL GAUSS(iav,iev,index2,AGA,AGB,AGC,AGD,AGE)
 
     DO k = 1 , ke
        DO i = iav , iev
           V(i,j,k,ne) = AGE(i,k,2)
        ENDDO
     ENDDO
     ! 
     ! BERECHNUNG DES V-IMPULSFLUSSES AM BODEN;
     ! DIESER FLUESS WIRD AUFSUMMIERT
     !
     DO i = iav , iev
        Bflvs(i,j) = Bflvs(i,j) + 0.5*(Tmcm(i,j)+Tmcm(i,j+1))                                      &
                   & *edg*(A2t(ke1)*V(i,j,ke,na)+A1t(ke1)*V(i,j,ke,ne))*dtdeh
     ENDDO
 
  ENDIF ! IF(J.LE.JEV)
  ! 
  ! ZYKLISCHE VERTAUSCHUNG DER SCHEIBENINDIZES
  !
  CALL DYN_SWAP_INDICES
  !
ENDDO  ! J = JAH, JEH
! 
! DIESER FLUESS WIRD AUFSUMMIERT
! OMEGA-MITTELWERTE BESTIMMTER MODELLFLAECHEN BILDEN
!
zdeltf = FLOAT((iehgg-iahgg+1)*(jehgg-jahgg+1))
Zom850m = Zom850m/zdeltf
Zom500m = Zom500m/zdeltf
Zom300m = Zom300m/zdeltf
! 
!-----------------------------------------------------------------------
END SUBROUTINE PROGEXP
