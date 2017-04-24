
MODULE mo_memory_dynamics

USE MO_PARORG, ONLY: IE,JE,KE,KE1

IMPLICIT NONE


INTEGER, PARAMETER :: N_TRACER = 4

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
 REAL, ALLOCATABLE, DIMENSION(:,:)    ::   ZQITS,  ZQIDIH
 REAL, ALLOCATABLE, DIMENSION(:)      ::   PSDT
 REAL, ALLOCATABLE, DIMENSION(:,:,:)  ::   ZTPA, ZGQD, ZTV, ZFIHF, ZLAPQI          
                                                                 
!NEUE FELDER FUER GEWICHTETE SIGMA-HORIZONTALDIFFUSION:          
 REAL, ALLOCATABLE, DIMENSION(:,:,:)  ::   ZPHF, ZPNF, ZDP, ZDP5, ZDPU, ZDPV,         &
!LM  ZDPINT: NONHYDROSTATIC VERTICAL PRESSURE DIFFERENCE
                                           ZDPINT
                                                                     
 REAL, ALLOCATABLE, DIMENSION(:,:,:)  ::   ZTMKVM, ZTMKVH, ZUTK, ZVTK 
 REAL, ALLOCATABLE, DIMENSION(:,:)    ::   ZTTK, ZQDTK, ZTTS, ZQDTS, ZQWTS,           &
                                           ZSODTA, ZTHDTA, ZPLAM
                                                                     
 REAL, ALLOCATABLE, DIMENSION(:,:,:)  ::   ZALOPN, ZALOPH, ZBETAK, ZFIH, ZPPHI,       &
                                           ZGU, ZGV, ZEKIN, ZZETA, ZETAS, ZSDIV                     
                                                                 
 REAL, ALLOCATABLE, DIMENSION(:,:,:)  ::   ZLAPT, ZLAPQD, ZLAPQW, ZLAPU, ZLAPV
                                                                     
 REAL, ALLOCATABLE, DIMENSION(:,:)    ::   ZTADV,  ZTDIFH, ZQDDIH,                    &
                                           ZQWDIH, ZALPOM, ZQDADV           
                                                                     
 REAL, ALLOCATABLE, DIMENSION(:,:)    ::   ZGRAD,  ZVZEQ, ZEDDPQ, ZUDIFH,             &
                                           ZVDIFH, ZUZEQ, ZDAMPU, ZDAMPV,             & 
!LM  ZEDDPQINT: INVERSE OF NONHYDROSTATIC VERTICAL                  &
                                           ZEDDPQINT                                      
                                                                     
!    PRESSURE DIFFERENCE FOR OMEGA-ALPHA TERM                        
 REAL, ALLOCATABLE, DIMENSION(:)      ::   ZQKOR, ZFCZET, ZDEXTU, ZDEXTV
                                                                     
!    FELDER FUER SUBROUTINE GAUSS                                    
 REAL, ALLOCATABLE, DIMENSION(:,:,:)  ::   AGA, AGB, AGC, AGD, AGE
                                                                     
!    FELDER FUER SUBROUTINE HQTOTQ                                   
 REAL, ALLOCATABLE, DIMENSION(:,:)    ::   HE, QDWE, PHFE, TSTART, GQDSTA, PHFSTA,    &
                                           QDLE, QDIE  

CONTAINS


SUBROUTINE allocate_dynamics
IMPLICIT NONE
INTEGER :: NJ
NJ = 1
      ALLOCATE (ZQITS(IE,KE  ), ZLAPQI(IE,KE,NJ+2), ZQIDIH(IE,KE) )
      ALLOCATE (PSDT(IE) )
      ALLOCATE (ZTPA(IE,JE,KE), ZGQD (IE,KE,NJ+2),                      &
               ZTV (IE,KE,NJ+4), ZFIHF(IE,KE,NJ+4)   )
                                                                     
!    NEUE FELDER FUER GEWICHTETE SIGMA-HORIZONTALDIFFUSION:          
     ALLOCATE (      ZPHF(IE,KE,NJ+4), ZPNF(IE,KE1,NJ+4), ZDP(IE,KE,NJ+2),   &
               ZDP5(IE,KE,NJ+4), ZDPU(IE,KE ,NJ+2), ZDPV(IE,KE,NJ+3),        & 
!LM  ZDPINT: NONHYDROSTATIC VERTICAL PRESSURE DIFFERENCE            
               ZDPINT(IE,KE,3)     )                                  
                                                                     
     ALLOCATE (ZTMKVM(IE,2:KE,NJ+2),                                   &
               ZTMKVH(IE,2:KE,NJ+2),                                   &
               ZTTK  (IE,  KE  ),                                   &
               ZQDTK (IE,  KE  ),                                   &
               ZUTK  (IE,  KE,NJ+1),                                   &
               ZVTK  (IE,  KE,NJ+1),                                   &
               ZTTS  (IE,  KE  ),                                   &
               ZQDTS (IE,  KE  ),                                   &
               ZQWTS (IE,  KE  ),                                   &
               ZSODTA(IE,  KE  ),                                   &
               ZTHDTA(IE,  KE  )  )
                                                                     
     ALLOCATE (ZALOPN(IE,KE1,NJ+2), ZALOPH(IE,KE ,NJ+3),                  &
               ZBETAK(IE,KE ,NJ+2), ZFIH  (IE,KE ,NJ+2),                  &
               ZPPHI (IE,KE ,NJ+1), ZPLAM (IE,KE   ),                  &
               ZGU   (IE,KE ,NJ+1), ZGV   (IE,KE ,NJ+2),                  &
               ZEKIN (IE,KE ,NJ+1), ZZETA (IE,KE ,NJ+1),                  &
               ZETAS (IE,KE1,NJ+1), ZSDIV (IE,KE1,NJ+1)  )                  
                                                                     
     ALLOCATE (ZLAPT(IE,KE,NJ+2), ZLAPQD(IE,KE,NJ+2), ZLAPQW(IE,KE,NJ+2),    &
               ZLAPU(IE,KE,NJ+2), ZLAPV (IE,KE,NJ+2)    )       
                                                                     
     ALLOCATE (ZTADV (IE,KE), ZTDIFH(IE,KE), ZQDDIH(IE,KE),         &
               ZQWDIH(IE,KE), ZALPOM(IE,KE), ZQDADV(IE,KE)   )
                                                                     
     ALLOCATE (ZGRAD (IE,KE), ZVZEQ (IE,KE),                        &
               ZEDDPQ(IE,KE), ZUDIFH(IE,KE),                        &
               ZVDIFH(IE,KE), ZUZEQ (IE,KE),                        &
               ZDAMPU (IE,KE), ZDAMPV (IE,KE),                      & 
!LM  ZEDDPQINT: INVERSE OF NONHYDROSTATIC VERTICAL                  &
               ZEDDPQINT(IE,KE)        )                              
                                                                     
!    PRESSURE DIFFERENCE FOR OMEGA-ALPHA TERM                        
                                                                     
     ALLOCATE (ZQKOR(IE), ZFCZET(IE),                               &
               ZDEXTU(IE), ZDEXTV(IE) )                               
                                                                     
!    FELDER FUER SUBROUTINE GAUSS                                    
     ALLOCATE (AGA(IE,KE,N_TRACER), AGB(IE,KE,N_TRACER),            &
               AGC(IE,KE,N_TRACER), AGD(IE,KE,N_TRACER),            &
               AGE(IE,KE,N_TRACER)  )
                                                                     
!    FELDER FUER SUBROUTINE HQTOTQ                                   
     ALLOCATE (HE    (IE,KE), QDWE  (IE,KE),                        &
               PHFE  (IE,KE),                                       &
               TSTART(IE,KE), GQDSTA(IE,KE),                        &
               PHFSTA(IE,KE),QDLE  (IE,KE),                         &
               QDIE  (IE,KE)   )

END SUBROUTINE allocate_dynamics  
 
SUBROUTINE deallocate_dynamics
DEALLOCATE(  ZQITS,  ZQIDIH, PSDT, ZTPA, ZGQD, ZTV, ZFIHF, ZLAPQI,  & 
             ZPHF, ZPNF, ZDP, ZDP5, ZDPU, ZDPV, ZDPINT,             &
             ZTMKVM, ZTMKVH, ZUTK, ZVTK, ZTTK, ZQDTK, ZTTS, ZQDTS,  &
             ZQWTS, ZSODTA, ZTHDTA, ZPLAM,                          &
             ZALOPN, ZALOPH, ZBETAK, ZFIH, ZPPHI,                   &
             ZGU, ZGV, ZEKIN, ZZETA, ZETAS, ZSDIV,                  &  
             ZLAPT, ZLAPQD, ZLAPQW, ZLAPU, ZLAPV,                   &
             ZTADV,  ZTDIFH, ZQDDIH, ZQWDIH, ZALPOM, ZQDADV,        &           
             ZGRAD,  ZVZEQ, ZEDDPQ, ZUDIFH, ZVDIFH, ZUZEQ, ZDAMPU,  &
             ZDAMPV, ZEDDPQINT, ZQKOR, ZFCZET, ZDEXTU, ZDEXTV,      &
             AGA, AGB, AGC, AGD, AGE, HE, QDWE, PHFE, TSTART,       &
             GQDSTA, PHFSTA, QDLE, QDIE  )
END SUBROUTINE deallocate_dynamics


END MODULE mo_memory_dynamics
