!
!----------------------------------------------------------------------------------------------
!
subroutine hydro_solve(dt)
 !
 use mo_grid, only: ib,je,ie,nx,ibg,ieg,jbg,jeg,kbg,keg, ke, ke1,ieke=>nxnz
 use mo_database
 use mo_parameters, only: nvar,ndim
 use mo_namelist, only: bc
 use mo_hydro, only: fill_guardcells_1D
 !
 implicit none
 !
 real, intent(inout) :: dt
 integer :: i,j,k
 !
 real :: zom850m, zom500m, zom300m
 REAL  ::  TGL (IE,JE,3), TGW (IE,JE,3), TGI (IE,JE,3),      &
       QDBL(IE,JE,3), QDBW(IE,JE,3), QDBI(IE,JE,3)
 REAL  ::   BFLHS  (IE,JE),                                  &
         BFLHSL (IE,JE), BFLHSW (IE,JE), BFLHSI (IE,JE),                &
         BFLQDS (IE,JE),                                                &
         BFLQDSL(IE,JE), BFLQDSW(IE,JE), BFLQDSI(IE,JE),                &
         BFLUS  (IE,JE), BFLVS  (IE,JE)                                  
 REAL :: TMKVMH(IE*(KE-1),JE,2)                                
 REAL :: TMCM(IE,JE), TMCHL(IE,JE), TMCHW(IE,JE), TMCHI(IE,JE)           
 REAL :: SOTHDT(IEKE,JE,2)                                   
 REAL ::   TTK(IEKE,JE), QDTK(IEKE,JE), UVTK(IEKE,JE,2)        
                                                               
 REAL ::   TTS(IEKE,JE), QDTS(IEKE,JE), QWTS(IEKE,JE)          
 INTEGER :: INFRL  (IE,JE), INFRW  (IE,JE),                   &
                        INFRI  (IE,JE)                                     
 REAL  ::   QITS  (IEKE,JE)                   
 REAL  :: DDMPUJ(JE), DDMPVJ(JE)                                
 !
 ! dummy variables
 !
 !
 !call fill_guardcells
 !
 !
 call PROGEXP                                                 &
   (ZOM850M, ZOM500M, ZOM300M,                                     &
     FIB    , FC     , PS    ,      &
    TGL   , TGW   , TGI    , QDBL   , QDBW   , QDBI , BFLHSL,      &
    BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, TG   , QDB   ,      &
    BFLHS , BFLQDS, BFLUS  , BFLVS  , U      , V    , T     ,      &
    QD    , QW    , FI     , TMKVMH , TMCHL  , TMCHW,              &
    TMCHI , TMCM  , TMCH   , SOTHDT , TTK    , QDTK , UVTK  ,      &
    TTS   , QDTS  , QWTS   , INFRL  , INFRW  , INFRI, QI    ,      &
    QITS  , PINT  , DWDT   , ETAS   , DDMPUJ , DDMPVJ,HYDRODP)    
 !
end subroutine hydro_solve
