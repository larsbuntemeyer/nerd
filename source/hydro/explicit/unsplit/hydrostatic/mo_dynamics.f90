!*==mo_dynoldmics.spg  processed by SPAG 6.72Rc at 17:28 on 30 Mar 2017
!
!
MODULE MO_DYNAMICS
USE mo_grid,      only: iah=>ib, ieh=>ie, jah=>jb, jeh=>je, iaa=>ibg, iea=>ieg, &
                        ke=>nz, ke1=>nz1, eddlam, edadphi, eddphi
USE mo_driver,    only: nnow, nold, nold2, nnow2, ed2dt 
USE mo_constants,         ONLY: R, RERD=>r_earth, wcpr, RDDRM1
USE MO_MEMORY_EC4,        ONLY: CPHI, A1T, A2T, TG, ACPHIR, VVFH, HYDRODP, FC, &
                                ETAS, BK 
USE MO_MEMORY_MAIN,       ONLY: T, QD, QW, QI, DWDT, TMCH, FIB, QDB, U, V,   &
                                PINT, FI
!USE MO_MEMORY_DYNAMICS,   ONLY: ZEDDPQ, ZDP, ZEDDPQINT, ZDPINT, ZGV,         &
!                                ZTADV, ZQDADV, ZGU, ZPLAM, ZPPHI, ZTV,       &
!                                ZALOPN, ZSDIV, ZBETAK, ZALPOM, AGB, AGD,     &
!                                ZTDIFH, ZTTK, ZTTS, ZTHDTA, ZSODTA, ZETAS,   &
!                                AGC, ZTMKVH, AGA, ZFIHF, AGE, N_TRACER,      &
!                                ZQDDIH, ZQDTK, ZQDTS, ZQWDIH, ZQWTS,         &
!                                ZQIDIH, ZQITS, ZLAPT, ZLAPQD, ZLAPQW, ZLAPU, &
!                                ZLAPV, ZUDIFH, ZVDIFH, ZLAPQI, ZDPU, ZDPV,   &
!                                ZDP5, ZTPA, ZPHF
USE MO_MEMORY_DYNAMICS 
USE MO_MAGNUS,            ONLY:   FGEW, FGQD 
!
! cphi    = cos(phi)
! ed2dt   = 1/(2*dt)
! zlap(f) = div(grad(f*dp))/dp : why is it weighted with pressure?
! c_p     = 1005.0   : specific heat of dry air (J/(KG*K))
! wcpr    = 1/c_p    : inverse of specific heat of dry air (4.7)
! wlk     = L_v      : latent heat of evaporation (4.7)
! R       = 287.05   : gas constant dry air
! Rd      = 461.51   : gas constant water vapour
! rddrm1  = Rd/R - 1
! T       = h - 1/c_p *  
!
! grid descprition
!
!   boundary  |                     |  boundary
! |  1  |  2  |  3  |  ... |  IE-2  | IE-1 |  IE  | 
!   IAA       | IAH           IEH   |        IEA
!
!  do i=iah,ieh : loop inside the subdomain
!  do i=1  ,iea : loop with boundary cells
!
IMPLICIT NONE
!
! PARAMETER definitions
INTEGER, PARAMETER :: INDEX1=1, INDEX2=2, INDEX3=3, INDEX4=4
! ZTKVZ , ZTKVL : GEWICHTSFAKTOREN ZUR HORIZONTALEN MITTELUNG DER
!                 VERTIKALEN TRANSPORTKOEFFIZIENTEN
REAL,    PARAMETER :: ZTKVZ = 0.9
REAL,    PARAMETER :: ZTKVL = (1.-ZTKVZ)*0.25
!
!BERECHNUNG MIT ODER OHNE MASSENFLUSS-KORREKTURSCHEMA
LOGICAL, PARAMETER :: LMASSF = .TRUE.
!
! local variables
INTEGER :: JM2, JM3, JM4, JN5, JO2, JO3, JO4, JO5, JS4, JS5, JSP, JU3, JU4, JM5 
INTEGER :: JU5, JZM, JZN, JZO, JZS, JZU
REAL    :: ZX2, ZA1A, ZA2A, ZALOG2, Z4DRERD
!
!
REAL :: LHDIFF2, LPTOP0, AKS2, AKS4
!
!
INTERFACE diffuse_horizontal_sigma
!   MODULE PROCEDURE diffuse_horizontal_sigma_2d
   MODULE PROCEDURE diffuse_horizontal_sigma_3d
   MODULE PROCEDURE diffuse_horizontal_sigma_3d1
END INTERFACE
!
CONTAINS
!
SUBROUTINE dyn_precompute
IMPLICIT NONE
INTEGER :: i,j,k,kp1
REAL    :: zzgew, z1, z2, z3, z4, zxo, zxu, zx1, &
           zfkio, zfkiu, zfdivx, zfdivo, zfdivu
ZX2 = .5*RERD*(CPHI(JZU,1)+CPHI(JZM,1))
ZXO = CPHI(JZM,1)*EDDPHI
ZXU = CPHI(JZU,1)*EDDPHI

!JZS = JAH - 2
!JZU = JAH - 1
!JZM = JAH
!JZO = JAH + 1
DO k = 1 , ke
  kp1 = k + 1
  DO i = iaa , iea
     zdp(i,k,ju3) = Hydrodp(i,jzu,k,nnow)
     zdp(i,k,jm3) = Hydrodp(i,jzm,k,nnow)
!LM   ADDED ZDPINT FOR OMEGA-ALPHA TERM
     ! --------------------------------------------
     ! pressure difference between half levels 
     ! "pressure thickness" of full level
     zdpint(i,k,ju3) = Pint(i,jzu,kp1,nnow) - Pint(i,jzu,k,nnow)
     zdpint(i,k,jm3) = Pint(i,jzm,kp1,nnow) - Pint(i,jzm,k,nnow)
!     DEFINITION VON ZDP5 FUER GEWICHTETE SIGMA HDIFFUSION:
!==============================================================
!    !!!! Should ZDP5 be based on PINT?
!==============================================================
     ! hydrostatic pressure
     zdp5(i,k,js5) = Hydrodp(i,jzs,k,nold)
     zdp5(i,k,ju5) = Hydrodp(i,jzu,k,nold)
     zdp5(i,k,jm5) = Hydrodp(i,jzm,k,nold)
     zdp5(i,k,jo5) = Hydrodp(i,jzo,k,nold)
     ! --------------------------------------------
     ! non-hydrostatic pressure at full levels
     zphf(i,k,js5) = 0.5*(Pint(i,jzs,k,nnow)+Pint(i,jzs,kp1,nnow))
     zphf(i,k,ju5) = 0.5*(Pint(i,jzu,k,nnow)+Pint(i,jzu,kp1,nnow))
     zphf(i,k,jm5) = 0.5*(Pint(i,jzm,k,nnow)+Pint(i,jzm,kp1,nnow))
     zphf(i,k,jo5) = 0.5*(Pint(i,jzo,k,nnow)+Pint(i,jzo,kp1,nnow))
     ! --------------------------------------------
     ! non-hydrostatic pressure at half levels
     zpnf(i,kp1,js5) = Pint(i,jzs,kp1,nnow)
     zpnf(i,kp1,ju5) = Pint(i,jzu,kp1,nnow)
     zpnf(i,kp1,jm5) = Pint(i,jzm,kp1,nnow)
     zpnf(i,kp1,jo5) = Pint(i,jzo,kp1,nnow)
 
     zalopn(i,kp1,ju3) = ALOG(zpnf(i,kp1,ju5))
     zalopn(i,kp1,jm3) = ALOG(zpnf(i,kp1,jm5))
     ! --------------------------------------------
     ! virtual temperature
     ! Part I, Eq. 4.4.7
     ! Tv = T * (1 + (Rd/R - 1) * Qd - Qw + Qi)
     ztv(i,k,js5) = T(i,jzs,k,nnow)*(1.0+rddrm1*Qd(i,jzs,k,nnow)-(Qw(i,jzs,k,nnow)+Qi(i,jzs,k,nnow)))
     ztv(i,k,ju5) = T(i,jzu,k,nnow)*(1.0+rddrm1*Qd(i,jzu,k,nnow)-(Qw(i,jzu,k,nnow)+Qi(i,jzu,k,nnow)))
     ztv(i,k,jm5) = T(i,jzm,k,nnow)*(1.0+rddrm1*Qd(i,jzm,k,nnow)-(Qw(i,jzm,k,nnow)+Qi(i,jzm,k,nnow)))
     ztv(i,k,jo5) = T(i,jzo,k,nnow)*(1.0+rddrm1*Qd(i,jzo,k,nnow)-(Qw(i,jzo,k,nnow)+Qi(i,jzo,k,nnow)))
     ! Magnus Equation
     zzgew = FGEW(T(i,jzu,k,nnow))
     zgqd(i,k,ju3) = FGQD(zzgew,zphf(i,k,ju5))
     zzgew = FGEW(T(i,jzm,k,nnow))
     zgqd(i,k,jm3) = FGQD(zzgew,zphf(i,k,jm5))
  ENDDO
ENDDO
 
! DEFINITION VON ZDPU UND ZDPV FUER GEWICHTETE SIGMA HDIFFUSION
DO k = 1 , ke
  ! FEHLERKORREKTUR FUER ZDPU-DEFINITION
  DO i = iaa , ieh + 1
     zdpu(i,k,ju3) = 0.5*(zdp5(i,k,ju5)+zdp5(i+1,k,ju5))
     zdpu(i,k,jm3) = 0.5*(zdp5(i,k,jm5)+zdp5(i+1,k,jm5))
  ENDDO
ENDDO
!
!KS   SPLITTED THE FOLLOWING LOOP INTO TWO TO REMOVE THE IF INSIDE A LOOP
!
!KS   FIRST PART WITH FOR ALL K<KE
DO k = 1 , ke - 1
  DO i = iaa , iea
     ! pressure at half levels between cells (i+1/2,j+1/2) (cell corners)
     zdpv(i,k,js4) = 0.5*(zdp5(i,k,js5)+zdp5(i,k,ju5))
     zdpv(i,k,ju4) = 0.5*(zdp5(i,k,ju5)+zdp5(i,k,jm5))
     zdpv(i,k,jm4) = 0.5*(zdp5(i,k,jm5)+zdp5(i,k,jo5))
!           ***************
     zalpom(i,k) = 0.0
     ! potential energy 
     ! FI: geoptential height [m]
     ! log(p_hl/p_fl) = log p_hl - log p_fl
     ! [zfihf] = J/kg = m^2 / s^2  : specific energy
     zfihf(i,k,js5) = Fi(i,jzs,k+1,nold2) + r*ztv(i,k,js5)*ALOG(zpnf(i,k+1,js5)/zphf(i,k,js5))       &
                    & /Dwdt(i,jzs,k,nnow)
     zfihf(i,k,ju5) = Fi(i,jzu,k+1,nold2) + r*ztv(i,k,ju5)*ALOG(zpnf(i,k+1,ju5)/zphf(i,k,ju5))       &
                    & /Dwdt(i,jzu,k,nnow)
     zfihf(i,k,jm5) = Fi(i,jzm,k+1,nold2) + r*ztv(i,k,jm5)*ALOG(zpnf(i,k+1,jm5)/zphf(i,k,jm5))       &
                    & /Dwdt(i,jzm,k,nnow)
     zfihf(i,k,jo5) = Fi(i,jzo,k+1,nold2) + r*ztv(i,k,jo5)*ALOG(zpnf(i,k+1,jo5)/zphf(i,k,jo5))       &
                    & /Dwdt(i,jzo,k,nnow)
     ! ztpa = T + zfihf / c_p
     ! heating due to potential height? adiabatic heating? 
     ztpa(i,jzs,k) = T(i,jzs,k,nold) + zfihf(i,k,js5)*wcpr
     ztpa(i,jzu,k) = T(i,jzu,k,nold) + zfihf(i,k,ju5)*wcpr
     ztpa(i,jzm,k) = T(i,jzm,k,nold) + zfihf(i,k,jm5)*wcpr
     ztpa(i,jzo,k) = T(i,jzo,k,nold) + zfihf(i,k,jo5)*wcpr
  ENDDO
ENDDO
!
!KS   SECOND PART FOR K=KE
!
DO i = iaa , iea
  ! pressure at half levels between cells (i+1/2,j+1/2) (cell corners)
  zdpv(i,ke,js4) = 0.5*(zdp5(i,ke,js5)+zdp5(i,ke,ju5))
  zdpv(i,ke,ju4) = 0.5*(zdp5(i,ke,ju5)+zdp5(i,ke,jm5))
  zdpv(i,ke,jm4) = 0.5*(zdp5(i,ke,jm5)+zdp5(i,ke,jo5))
  !        ***************
  zalpom(i,ke) = 0.0
  ! potential energy 
  ! FI: geoptential height [m]
  ! log(p_hl/p_fl) = log p_hl - log p_fl
  ! [zfihf] = J/kg = m^2 / s^2  : specific energy
  zfihf(i,ke,js5) = Fib(i,jzs) + r*ztv(i,ke,js5)*ALOG(zpnf(i,ke+1,js5)/zphf(i,ke,js5))             &
                  & /Dwdt(i,jzs,ke,nnow)
  zfihf(i,ke,ju5) = Fib(i,jzu) + r*ztv(i,ke,ju5)*ALOG(zpnf(i,ke+1,ju5)/zphf(i,ke,ju5))             &
                  & /Dwdt(i,jzu,ke,nnow)
  zfihf(i,ke,jm5) = Fib(i,jzm) + r*ztv(i,ke,jm5)*ALOG(zpnf(i,ke+1,jm5)/zphf(i,ke,jm5))             &
                  & /Dwdt(i,jzm,ke,nnow)
  zfihf(i,ke,jo5) = Fib(i,jzo) + r*ztv(i,ke,jo5)*ALOG(zpnf(i,ke+1,jo5)/zphf(i,ke,jo5))             &
                  & /Dwdt(i,jzo,ke,nnow)
  ! ztpa = T + zfihf / c_p
  ! heating due to potential height? adiabatic heating? 
  ztpa(i,jzs,ke) = T(i,jzs,ke,nold) + zfihf(i,ke,js5)*wcpr
  ztpa(i,jzu,ke) = T(i,jzu,ke,nold) + zfihf(i,ke,ju5)*wcpr
  ztpa(i,jzm,ke) = T(i,jzm,ke,nold) + zfihf(i,ke,jm5)*wcpr
  ztpa(i,jzo,ke) = T(i,jzo,ke,nold) + zfihf(i,ke,jo5)*wcpr
ENDDO
!
!KS   introduction of pressure top not equal 0.0
IF ( lptop0 ) THEN
  zbetak(iaa:iea,1,ju3) = zalog2
  zbetak(iaa:iea,1,jm3) = zalog2
  zbetak(iaa:iea,1,jo3) = zalog2
ELSE
  DO i = iaa , iea
     zbetak(i,1,ju3) = 1. - zpnf(i,1,ju5)/zdp(i,1,ju3)*(zalopn(i,2,ju3)-zalopn(i,1,ju3))
     zbetak(i,1,jm3) = 1. - zpnf(i,1,jm5)/zdp(i,1,jm3)*(zalopn(i,2,jm3)-zalopn(i,1,jm3))
     zbetak(i,1,jo3) = 1. - zpnf(i,1,jo5)/zdp(i,1,jo3)*(zalopn(i,2,jo3)-zalopn(i,1,jo3))
  ENDDO
ENDIF
!
!KS   SPLITTED THE FOLLOWING LOOP INTO TWO TO REMOVE THE IF INSIDE A LOOP
!
!KS   FIRST PART FOR K=1
DO i = iaa , iea
  zaloph(i,1,ju3) = zalopn(i,2,ju3) - zbetak(i,1,ju3)
  zaloph(i,1,jm3) = zalopn(i,2,jm3) - zbetak(i,1,jm3)
  zfih(i,1,ju3) = Fi(i,jzu,2,nnow2) + r*ztv(i,1,ju5)*zbetak(i,1,ju3)/Dwdt(i,jzu,1,nnow)
  zfih(i,1,jm3) = Fi(i,jzm,2,nnow2) + r*ztv(i,1,jm5)*zbetak(i,1,jm3)/Dwdt(i,jzm,1,nnow)
ENDDO
!
!KS   SECOND PART FOR K>1
!
DO k = 2 , ke - 1
  DO i = iaa , iea
     zbetak(i,k,ju3) = 1. - zpnf(i,k,ju5)/zdp(i,k,ju3)*(zalopn(i,k+1,ju3)-zalopn(i,k,ju3))
     zbetak(i,k,jm3) = 1. - zpnf(i,k,jm5)/zdp(i,k,jm3)*(zalopn(i,k+1,jm3)-zalopn(i,k,jm3))
     zfih(i,k,ju3) = Fi(i,jzu,k+1,nnow2) + r*ztv(i,k,ju5)*zbetak(i,k,ju3)/Dwdt(i,jzu,k,nnow)
     zfih(i,k,jm3) = Fi(i,jzm,k+1,nnow2) + r*ztv(i,k,jm5)*zbetak(i,k,jm3)/Dwdt(i,jzm,k,nnow)
     zaloph(i,k,ju3) = zalopn(i,k+1,ju3) - zbetak(i,k,ju3)
     zaloph(i,k,jm3) = zalopn(i,k+1,jm3) - zbetak(i,k,jm3)
  ENDDO
ENDDO
!KS
!     CORIOLISPARAMETER AM ZETA-PUNKT BEREITSTELLEN
!
DO i = iaa , ieh + 1
  zfczet(i) = 0.25*(Fc(i,jzu)+Fc(i+1,jzu)+Fc(i,jzm)+Fc(i+1,jzm))
ENDDO
 
DO i = iaa , iea
  zbetak(i,ke,ju3) = 1. - zpnf(i,ke,ju5)/zdp(i,ke,ju3)*(zalopn(i,ke1,ju3)-zalopn(i,ke,ju3))
  zbetak(i,ke,jm3) = 1. - zpnf(i,ke,jm5)/zdp(i,ke,jm3)*(zalopn(i,ke1,jm3)-zalopn(i,ke,jm3))
  zaloph(i,ke,ju3) = zalopn(i,ke1,ju3) - zbetak(i,ke,ju3)
  zaloph(i,ke,jm3) = zalopn(i,ke1,jm3) - zbetak(i,ke,jm3)
  zfih(i,ke,ju3) = Fib(i,jzu) + r*ztv(i,ke,ju5)*zbetak(i,ke,ju3)/Dwdt(i,jzu,ke,nnow)
  zfih(i,ke,jm3) = Fib(i,jzm) + r*ztv(i,ke,jm5)*zbetak(i,ke,jm3)/Dwdt(i,jzm,ke,nnow)
ENDDO
!
! WEITERE LOKALE FELDER AN U/ZETA-GITTERPUNKTEN (HF) VORBESETZEN
!
DO k = 1 , ke
  DO i = iaa , ieh + 1
     zgu(i,k,jm2) = .5*(zdp(i,k,jm3)+zdp(i+1,k,jm3))*U(i,jzm,k,nnow)
     z1 = zx2*zfczet(i)
     z2 = eddlam*(V(i+1,jzu,k,nnow)-V(i,jzu,k,nnow))
     z3 = zxo*U(i,jzm,k,nnow) - zxu*U(i,jzu,k,nnow)
     z4 = (zdp(i,k,ju3)+zdp(i+1,k,ju3))*Cphi(jzu,1) + (zdp(i,k,jm3)+zdp(i+1,k,jm3))*Cphi(jzm,1)
     zzeta(i,k,jm2) = z4drerd*(z1+z2-z3)/z4
  ENDDO
 
!     WEITERE LOKALE FELDER AN H/V-GITTERPUNKTEN VORBESETZEN
  zx1 = 0.5*r*edadphi
  zfkio = Cphi(jzm,2)/Cphi(jzm,1)
  zfkiu = Cphi(jzu,2)/Cphi(jzm,1)
  zfdivx = Acphir(jzm,1)*eddlam
  zfdivo = Acphir(jzm,1)*eddphi*Cphi(jzm,2)
  zfdivu = Acphir(jzm,1)*eddphi*Cphi(jzu,2)
 
  DO i = iah - 1 , ieh + 1
     zgv(i,k,ju3) = .5*(zdp(i,k,ju3)+Hydrodp(i,jzm,k,nnow))*V(i,jzu,k,nnow)
     zgv(i,k,jm3) = .5*(zdp(i,k,jm3)+Hydrodp(i,jzm+1,k,nnow))*V(i,jzm,k,nnow)
     zpphi(i,k,jm2) = zx1*(ztv(i,k,jm5)+ztv(i,k,ju5))*(zaloph(i,k,jm3)-zaloph(i,k,ju3))
     zekin(i,k,jm2) = .25*(U(i-1,jzm,k,nnow)**2+U(i,jzm,k,nnow)**2+zfkiu*V(i,jzu,k,nnow)                 &
                    & **2+zfkio*V(i,jzm,k,nnow)**2)
  ENDDO
ENDDO
 
DO k = 1 , ke
  DO i = iah - 1 , ieh + 1
     zsdiv(i,k+1,jm2) = zsdiv(i,k,jm2) + zfdivx*(zgu(i,k,jm2)-zgu(i-1,k,jm2)) + zfdivo*zgv(i,k,jm3)&
                      & - zfdivu*zgv(i,k,ju3)
  ENDDO
ENDDO
 
DO k = 2 , ke
  DO i = iah - 1 , ieh + 1
     zetas(i,k,jm2) = Bk(k)*zsdiv(i,ke1,jm2) - zsdiv(i,k,jm2)
  ENDDO
  Etas(iah-1:ieh+1,jzm,k) = zetas(iah-1:ieh+1,k,jm2)
ENDDO

END SUBROUTINE dyn_precompute
!
!
!
SUBROUTINE dyn_init_diffusion(j)
IMPLICIT NONE
INTEGER, INTENT(IN) :: j
INTEGER :: i,k
!LOKALE HILFSFELDER FUER HORIZONTALDIFFUSION
 DO k = 1 , ke
    DO i = iah - 1 , ieh + 1
       !GEWICHTUNG VON ZLAPT, ZLAPQD, ZLAPQW MIT GEWICHT DER SCHICHT
       ! diffusion with sigma weighting at cell centers
       ! last four digits give the coordinates for weighting
       zlapt (i,k,ju3) = diffuse_horizontal_sigma(ztpa          , i,jzu,k, ju3,ju4,js4,jzu)
       zlapqd(i,k,ju3) = diffuse_horizontal_sigma(QD(:,:,:,nold), i,jzu,k, ju3,ju4,js4,jzu)
       zlapqw(i,k,ju3) = diffuse_horizontal_sigma(QW(:,:,:,nold), i,jzu,k, ju3,ju4,js4,jzu)
       zlapqi(i,k,ju3) = diffuse_horizontal_sigma(QI(:,:,:,nold), i,jzu,k, ju3,ju4,js4,jzu)
       zlapt (i,k,jm3) = diffuse_horizontal_sigma(ztpa          , i,jm5,k, jm3,jm4,ju4,jm5)
       zlapqd(i,k,jm3) = diffuse_horizontal_sigma(QD(:,:,:,nold), i,jm5,k, jm3,jm4,ju4,jm5)
       zlapqw(i,k,jm3) = diffuse_horizontal_sigma(QW(:,:,:,nold), i,jm5,k, jm3,jm4,ju4,jm5)
       zlapqi(i,k,jm3) = diffuse_horizontal_sigma(QI(:,:,:,nold), i,jm5,k, jm3,jm4,ju4,jm5)
       ! diffusion without sigma weighting at cell faces
       ! last digit indicates left or right cell face
       zlapu (i,k,ju3) = diffuse_horizontal_zeta(U(:,:,:,nold), i,jzu,k, 0)   
       zlapv (i,k,ju3) = diffuse_horizontal_zeta(V(:,:,:,nold), i,jzu,k, 1)   
       zlapu (i,k,jm3) = diffuse_horizontal_zeta(U(:,:,:,nold), i,jzm,k, 0)   
       zlapv (i,k,jm3) = diffuse_horizontal_zeta(V(:,:,:,nold), i,jzm,k, 1)   
    ENDDO
 ENDDO

END SUBROUTINE dyn_init_diffusion
!
SUBROUTINE dyn_advect_diffuse(j)
!*--DYN_ADVECT_AND_DIFFUSE33
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: j
!
! Local variables
!
REAL    :: za1, za2, za3, zagcm, ztmkvhm, zagam, zagat, zagct, zgvco, zgvcu    
INTEGER :: i, k, i_tracer 
!
!
!        2. H - QDW - PROGNOSE
!        ---------------------
!KS
!        Coefficients of the tridiagonoldl system (Eq. 5.3.17)
!
!        AGA_K (PSI^t+dt)_k-1 +
!        AGB_K (PSI^t+dt)_k   + 
!        AGC_K (PSI^t+dt)_k+1 =
!        AGD_K
!
!        PSI stands for T, QW, QD or QI depending on the last index
!        of AGA, AGB, AGC and AGD
!KS
!
DO k = 1 , ke
  DO i = iah , ieh
     ZEDDPQ(i,k) = 1./ZDP(i,k,jm3)
     !LM ADDED INVERSE OF NONHYDROSTATIC PRESSURE DIFFERENCE
     ZEDDPQINT(i,k) = 1./ZDPINT(i,k,jm3)
     ! horizontal advection
     ZTADV(i,k)     = advect_horizontal(T,i,j,k,nnow) 
     ZQDADV(i,k)    = advect_horizontal(QD,i,j,k,nnow)
     ! 
     ! alpha*omega 
     ZALPOM(i,k)    = alpha_omega(i,j,k) !(.5*(za1+za2)*ZEDDPQ(i,k)+za3)
     !
     ! Source terms and explizit solutions go into the RHS of the equations (AGD)
     !
     ! TODO: At this point, the sum depends on the order of the terms, so we do
     ! not clean up here right now. We have to check out this for different
     ! compuiler options and version.
     AGB(i,k,1) = ed2dt
     !!                                  ! horiz. diffusion   phy. tend.    horiz. advection
     !AGD(i,k,1) = ed2dt*T (i,j,k,nold) + ZTDIFH(i,k)      + ZTTS(i,k)   + ZTADV(i,k)      +               &
     !!                                  ! alpha*omega        phy. tend.    phy. tend.        phy. tend.
     !                                  + ZALPOM(i,k)*wcpr + ZTTK(i,k)   + ZTHDTA(i,k)     + ZSODTA(i,k)
     AGD(i,k,1) = ed2dt*T (i,j,k,nold) + ZTADV(i,k)      + ZTDIFH(i,k)   + ZALPOM(i,k)*wcpr      +               &
     !                                  ! alpha*omega        phy. tend.    phy. tend.        phy. tend.
                                       + ZTTK(i,k) + ZTTS(i,k)   + ZTHDTA(i,k)     + ZSODTA(i,k)
     !                                  !
     !                                  ! horiz. diffusion   phy. tend.    horiz. advection
     !AGD(i,k,2) = ed2dt*QD(i,j,k,nold) + ZQDDIH(i,k)      + ZQDTS(i,k)  + ZQDADV(i,k)      + ZQDTK(i,k)
     AGD(i,k,2) = ed2dt*QD(i,j,k,nold) + ZQDADV(i,k)      + ZQDDIH(i,k)  + ZQDTK(i,k)      + ZQDTS(i,k)
                                       ! 
                                       ! horiz. diffusion   phy. tend.
     AGD(i,k,3) = ed2dt*QW(i,j,k,nold) + ZQWDIH(i,k)      + ZQWTS(i,k)
                                       !
                                       ! horiz. diffusion   phy. tend.
     AGD(i,k,4) = ed2dt*QI(i,j,k,nold) + ZQIDIH(i,k)      + ZQITS(i,k)
     ! 
     !  OKT. 2015: VERTIKAL INTEGRIERTER FEUCHTETRANSPORT
     !         QJX(I,J) = QJX(I,J) +
     !      ZFADVX*ZGU(I,K,JM2)*(QD(I+1,J,K,NJ)+QD(I,J,K,NJ))/G*DT
     !         QJY(I,J) = QJY(I,J) +
     !      ZFADVY*ZGVCO*(QD(I,J+1,K,NJ)+QD(I,J,K,NJ))/G*DT
     !
  ENDDO
ENDDO
!
!  VERTIKALADVEKTION UND -DIFFUSION
!
DO i = iah , ieh
  zagcm   = .5*ZETAS(i,2,jm2)*ZEDDPQ(i,1)
  ztmkvhm = ztkvz*ZTMKVH(i,2,jm3) + ztkvl*(ZTMKVH(i+1,2,jm3)+ZTMKVH(i,2,jo3)+ZTMKVH(i-1,2,jm3)     &
            +ZTMKVH(i,2,ju3))
  zagct   = -ztmkvhm*ZEDDPQ(i,1)
  AGC(i,1,1) = zagcm*za1a + zagct*A1T(2)
  AGB(i,1,1) = AGB(i,1,1) - AGC(i,1,1)
  AGD(i,1,1) = AGD(i,1,1) - (za2a*zagcm+A2T(2)*zagct)*(T(i,j,2,nold)-T(i,j,1,nold))                    &
             & - wcpr*zagct*(ZFIHF(i,2,jm5)-ZFIHF(i,1,jm5))
  AGD(i,1,2) = AGD(i,1,2) - (za2a*zagcm+A2T(2)*zagct)*(QD(i,j,2,nold)-QD(i,j,1,nold))
  AGD(i,1,3) = AGD(i,1,3) - (za2a*zagcm+A2T(2)*zagct)*(QW(i,j,2,nold)-QW(i,j,1,nold))
  AGD(i,1,4) = AGD(i,1,4) - (za2a*zagcm+A2T(2)*zagct)*(QI(i,j,2,nold)-QI(i,j,1,nold))
  DO i_tracer = 2 , N_TRACER
     AGB(i,1,i_tracer) = AGB(i,1,1)
     AGC(i,1,i_tracer) = AGC(i,1,1)
  ENDDO
ENDDO
!
!  Equations from Part II of EM/DM
!
DO k = 2 , ke - 1
  DO i = iah , ieh
!        1.RHS Eq. 5.4.5
     zagam = -.5*ZETAS(i,k,jm2)*ZEDDPQ(i,k)
!        1.RHS Eq. 5.4.6
     zagcm = .5*ZETAS(i,k+1,jm2)*ZEDDPQ(i,k)
!        2.RHS Eq. 5.4.5
     ztmkvhm = ztkvz*ZTMKVH(i,k,jm3) + ztkvl*(ZTMKVH(i+1,k,jm3)+ZTMKVH(i,k,jo3)+ZTMKVH(i-1,k,jm3)  &
             & +ZTMKVH(i,k,ju3))
     zagat = -ztmkvhm*ZEDDPQ(i,k)
!        2.RHS Eq. 5.4.6
     ztmkvhm = ztkvz*ZTMKVH(i,k+1,jm3) + ztkvl*(ZTMKVH(i+1,k+1,jm3)+ZTMKVH(i,k+1,jo3)              &
             & +ZTMKVH(i-1,k+1,jm3)+ZTMKVH(i,k+1,ju3))
     zagct = -ztmkvhm*ZEDDPQ(i,k)
!        Eq. 5.4.5
     AGA(i,k,1) = zagam*za1a + zagat*A1T(k)
!        Eq. 5.4.6
     AGC(i,k,1) = zagcm*za1a + zagct*A1T(k+1)
!        Eq. 5.4.7
     AGB(i,k,1) = AGB(i,k,1) - AGA(i,k,1) - AGC(i,k,1)
!        Eq. 5.4.8
     AGD(i,k,1) = AGD(i,k,1) - (za2a*zagam+A2T(k)*zagat)*(T(i,j,k-1,nold)-T(i,j,k,nold))               &
                & - (za2a*zagcm+A2T(k+1)*zagct)*(T(i,j,k+1,nold)-T(i,j,k,nold))                        &
                & - wcpr*zagat*(ZFIHF(i,k-1,jm5)-ZFIHF(i,k,jm5))                                   &
                & - wcpr*zagct*(ZFIHF(i,k+1,jm5)-ZFIHF(i,k,jm5))
     AGD(i,k,2) = AGD(i,k,2) - (za2a*zagam+A2T(k)*zagat)*(QD(i,j,k-1,nold)-QD(i,j,k,nold))             &
                & - (za2a*zagcm+A2T(k+1)*zagct)*(QD(i,j,k+1,nold)-QD(i,j,k,nold))
     AGD(i,k,3) = AGD(i,k,3) - (za2a*zagam+A2T(k)*zagat)*(QW(i,j,k-1,nold)-QW(i,j,k,nold))             &
                & - (za2a*zagcm+A2T(k+1)*zagct)*(QW(i,j,k+1,nold)-QW(i,j,k,nold))
     AGD(i,k,4) = AGD(i,k,4) - (za2a*zagam+A2T(k)*zagat)*(QI(i,j,k-1,nold)-QI(i,j,k,nold))             &
                & - (za2a*zagcm+A2T(k+1)*zagct)*(QI(i,j,k+1,nold)-QI(i,j,k,nold))
     DO i_tracer = 2 , N_TRACER
        AGA(i,k,i_tracer) = AGA(i,k,1)
        AGB(i,k,i_tracer) = AGB(i,k,1)
        AGC(i,k,i_tracer) = AGC(i,k,1)
     ENDDO
  ENDDO
ENDDO
!
DO i = iah , ieh
  zagam = -.5*ZETAS(i,ke,jm2)*ZEDDPQ(i,ke)
  ztmkvhm = ztkvz*ZTMKVH(i,ke,jm3) + ztkvl*(ZTMKVH(i+1,ke,jm3)+ZTMKVH(i,ke,jo3)+ZTMKVH(i-1,ke,jm3) &
          & +ZTMKVH(i,ke,ju3))
  zagat = -ztmkvhm*ZEDDPQ(i,ke)
  zagct = -TMCH(i,j)*ZEDDPQ(i,ke)
  AGA(i,ke,1) = zagam*za1a + zagat*A1T(ke)
  AGB(i,ke,1) = AGB(i,ke,1) - AGA(i,ke,1) - zagct*A1T(ke1)
  AGD(i,ke,1) = AGD(i,ke,1) - (za2a*zagam+A2T(ke)*zagat)*(T(i,j,ke-1,nold)-T(i,j,ke,nold)) - A2T(ke1)  &
              & *zagct*(TG(i,j,nold)-T(i,j,ke,nold)) - A1T(ke1)*zagct*TG(i,j,nold)                       &
              & - wcpr*zagat*(ZFIHF(i,ke-1,jm5)-ZFIHF(i,ke,jm5))                                   &
              & - wcpr*zagct*(FIB(i,j)-ZFIHF(i,ke,jm5))
  AGD(i,ke,2) = AGD(i,ke,2) - (za2a*zagam+A2T(ke)*zagat)*(QD(i,j,ke-1,nold)-QD(i,j,ke,nold)) - A2T(ke1)&
              & *zagct*(QDB(i,j,nold)-QD(i,j,ke,nold)) - A1T(ke1)*zagct*QDB(i,j,nold)
  AGD(i,ke,3) = AGD(i,ke,3) - (za2a*zagam+A2T(ke)*zagat)*(QW(i,j,ke-1,nold)-QW(i,j,ke,nold)) + A2T(ke1)&
              & *zagct*QW(i,j,ke,nold)
  AGD(i,ke,4) = AGD(i,ke,4) - (za2a*zagam+A2T(ke)*zagat)*(QI(i,j,ke-1,nold)-QI(i,j,ke,nold)) + A2T(ke1)&
              & *zagct*QI(i,j,ke,nold)
  !
  ! the matrix coefficients for h-qdw are all the same
  ! since they are all transported by the same scheme:
  ! implicit VERTICAL advection and turbulent fluxes.
  ! HORIZONTAL advection is treated explicitily, that
  ! is why the advection part is directly added onto
  ! the RHS vector AGD.
  !
  DO i_tracer = 2 , N_TRACER
     AGA(i,ke,i_tracer) = AGA(i,ke,1)
     AGB(i,ke,i_tracer) = AGB(i,ke,1)
  ENDDO
ENDDO
 
!  GAUSS - ELIMINATION UND H-QDW-ZERLEGUNG
DO i_tracer = 1 , N_TRACER
  CALL GAUSS(iah,ieh,i_tracer,AGA,AGB,AGC,AGD,AGE)
ENDDO
 
 
END SUBROUTINE dyn_advect_diffuse

SUBROUTINE prepare_RHS 
END SUBROUTINE prepare_RHS

SUBROUTINE dyn_diffuse_horizontal(j)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: j
  INTEGER :: i,k
  REAL    :: div4
  ! LOKALE HILFSFELDER FUER HORIZONTALDIFFUSION
  DO k = 1 , ke
     DO i = iah - 1 , ieh + 1
        ! GEWICHTETE HDIFF
        !GEWICHTUNG VON ZLAPT, ZLAPQD, ZLAPQW MIT GEWICHT DER SCHICHT
        ! diffusion with sigma weighting at cell centers
        ! last four digits give the coordinates for weighting
        zlapt (i,k,jo3) = diffuse_horizontal_sigma(ztpa,           i,j+1,k, jo3,jo4,jm4,jo5) 
        zlapqd(i,k,jo3) = diffuse_horizontal_sigma(QD(:,:,:,nold), i,j+1,k, jo3,jo4,jm4,jo5) 
        zlapqw(i,k,jo3) = diffuse_horizontal_sigma(QW(:,:,:,nold), i,j+1,k, jo3,jo4,jm4,jo5) 
        zlapqi(i,k,jo3) = diffuse_horizontal_sigma(QI(:,:,:,nold), i,j+1,k, jo3,jo4,jm4,jo5) 
        ! diffusion without sigma weighting at cell faces
        ! last digit indicates left or right cell face
        zlapu(i,k,jo3)  = diffuse_horizontal_zeta(U(:,:,:,nold), i,j+1,k, 0) 
        zlapv(i,k,jo3)  = diffuse_horizontal_zeta(V(:,:,:,nold), i,j+1,k, 1)  
     ENDDO
  ENDDO
 
  !AM RECHTEN UND LINKEN RAND H-DIFFUSION 2. ORDNUNG
  DO k = 1 , ke
     DO i = iah , ieh
        IF(lhdiff2 .OR. at_boundary(i,j,k)) THEN
           ! aks2 = 1.0 / (2*dt*pi**2)
           ! aks4 = 1.0 / (2*dt*pi**4)
           ztdifh(i,k) = Vvfh(k)*aks2*zlapt (i,k,jm3)
           zqddih(i,k) = Vvfh(k)*aks2*zlapqd(i,k,jm3)
           zqwdih(i,k) = Vvfh(k)*aks2*zlapqw(i,k,jm3)
           zudifh(i,k) = Vvfh(k)*aks2*zlapu (i,k,jm3)
           zvdifh(i,k) = Vvfh(k)*aks2*zlapv (i,k,jm3)
           zqidih(i,k) = Vvfh(k)*aks2*zlapqi(i,k,jm3)
        ELSE
           ! IM INNEREN DES MODELLGEBIETES H-DIFFUSION 4. ORDNUNG
           ! GEWICHTETE HDIFF
           ! TODO: The diffusion factor of the 4th order diffusion has
           ! to be included in the function. Multiplying it here afterwards
           ! will give different results! I don't understand that yet.
           ztdifh(i,k)  = diffuse_4th_order(zlapt ,i,j,k)
           zqddih(i,k)  = diffuse_4th_order(zlapqd,i,j,k)
           zqwdih(i,k)  = diffuse_4th_order(zlapqw,i,j,k)
           zqidih(i,k)  = diffuse_4th_order(zlapqi,i,j,k)
           ! diffuse velocities at cell face               
           zudifh(i,k) = -Vvfh(k)*aks4 * diffuse_horizontal_zeta_slice(zlapu,i,j,ju3,jm3,jo3,k,0)        
           zvdifh(i,k) = -Vvfh(k)*aks4 * diffuse_horizontal_zeta_slice(zlapv,i,j,ju3,jm3,jo3,k,1)        
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE dyn_diffuse_horizontal
!
!
!> explizit horizontal diffusion at zeta point
!
! CPHI  (J,1) = COS(ZPHIS*DEGRAD)              : cos(phi) at cell centers
! CPHI  (J,2) = COS((ZPHIS + 0.5*DPHI)*DEGRAD) : cos(phi) at cell faces
! Diffusion operator is dimensionless.
REAL FUNCTION diffuse_horizontal_zeta1(var,i,jm1,j,jp1,k,nt,np)
  IMPLICIT NONE
  REAL, INTENT(IN) :: var(:,:,:,:)
  INTEGER, INTENT(IN) :: i,jm1,j,jp1,k,nt,np
  ! nt: time level
  ! np: coordinate of cos(phi)
  diffuse_horizontal_zeta1 =                                                     & 
                   var(i+1,j,k,nt) + var(i-1,j,k,nt) - 2.0*var(i,j,k,nt)  &
 + ( Cphi(j+np,2-np) *(var(i,  jp1,k,nt) - var(i,  j,k,nt))                 &
 -   Cphi(jm1+np,  2-np) *(var(i,  j,k,nt) - var(i,    jm1,k,nt))  )              &
 /   Cphi(j,1+np)
END FUNCTION diffuse_horizontal_zeta1
!
!
!
REAL FUNCTION diffuse_horizontal_zeta(var,i,j,k,np)
  IMPLICIT NONE
  REAL,    INTENT(IN) :: var(:,:,:)
  INTEGER, INTENT(IN) :: i,j,k,np
  ! nt: time level
  ! np: coordinate of cos(phi)
  diffuse_horizontal_zeta =                                            & 
                        var(i+1,j,  k) + var(i-1,j,k) - 2.0*var(i,j,k) &
 + ( Cphi(j+np  ,2-np)*(var(i,  j+1,k) - var(i,  j,k))                 &
 -   Cphi(j-1+np,2-np)*(var(i,  j,  k) - var(i,j-1,k))  )              &
 /   Cphi(j     ,1+np)
END FUNCTION diffuse_horizontal_zeta
!
!
!
REAL FUNCTION diffuse_horizontal_zeta_slice(var,i,j,jm1,jn,jp1,k,np)
  IMPLICIT NONE
  REAL,    INTENT(IN) :: var(:,:,:)
  INTEGER, INTENT(IN) :: i,j,jm1,jn,jp1,k,np
  ! nt: time level
  ! np: coordinate of cos(phi)
  diffuse_horizontal_zeta_slice =                     & 
           (             var(i+1,k,jn ) + var(i-1,k,jn ) - 2.0*var(i,k,jn) &
 + ( Cphi(j+np  ,2-np)*(var(i,  k,jp1) - var(i,  k,jn ))                  &
 -   Cphi(j-1+np,2-np)*(var(i,  k,jn ) - var(i,  k,jm1))  )              &
 /   Cphi(j     ,1+np))
END FUNCTION diffuse_horizontal_zeta_slice
!
!
!> explizit horizontal diffusion with sigma weighting on a latitude slice
!
REAL FUNCTION diffuse_horizontal_sigma_2d(var,i,j,k,nt)
  IMPLICIT NONE
  REAL,    INTENT(IN) :: var(:,:,:)
  INTEGER, INTENT(IN) :: i,j,k,nt
  ! GEWICHTETE HDIFF
  diffuse_horizontal_sigma_2d =                                       &
    (              (var(i+1,k,jo5)-var(i,  k,jo5))*zdpu(i,  k,jo3)    &
     -             (var(i  ,k,jo5)-var(i-1,k,jo5))*zdpu(i-1,k,jo3)    &
     +(Cphi(j+1,2)*(var(i,  k,jn5)-var(i,  k,jo5))*zdpv(i,  k,jo4)    &
     - Cphi(j,  2)*(var(i,  k,jo5)-var(i,  k,jm5))*zdpv(i,  k,jm4))   & 
     / Cphi(j+1,1)) / zdp5(i,k,jo5)
END FUNCTION diffuse_horizontal_sigma_2d
!
!
!> explizit horizontal diffusion with sigma weighting at grid point
!
REAL FUNCTION diffuse_horizontal_sigma_3d(var,i,j,k,jw1,jw2,jw3,jwn)
  IMPLICIT NONE
  REAL,    INTENT(IN) :: var(:,:,:)
  INTEGER, INTENT(IN) :: i,j,k,jw1,jw2,jw3,jwn
  !
  ! GEWICHTETE HDIFF
  diffuse_horizontal_sigma_3d =                                       &
     (            (var(i+1,j  ,k)-var(i,  j  ,k))*zdpu(i,  k,jw1)     &
   -              (var(i,  j  ,k)-var(i-1,j  ,k))*zdpu(i-1,k,jw1)     &
   +( Cphi(j,2)*  (var(i,  j+1,k)-var(i,  j  ,k))*zdpv(i,  k,jw2)     &
     -Cphi(j-1,2)*(var(i,  j  ,k)-var(i,  j-1,k))*zdpv(i,  k,jw3))    & 
   /  Cphi(j,1)) / zdp5(i,k,jwn)
  !
END FUNCTION diffuse_horizontal_sigma_3d
!
!
!
REAL FUNCTION diffuse_4th_order(var,i,j,k)
  IMPLICIT NONE
  REAL,    INTENT(IN) :: var(:,:,:)
  INTEGER, INTENT(IN) :: i,j,k
  diffuse_4th_order = -Vvfh(k)*aks4*                                & 
        (            (var(i+1,k,jm3)-var(i,  k,jm3))*zdpu(i,  k,jm3)   &
        -            (var(i,  k,jm3)-var(i-1,k,jm3))*zdpu(i-1,k,jm3)   &
       +(Cphi(j,  2)*(var(i,  k,jo3)-var(i,  k,jm3))*zdpv(i,  k,jm4)   &
       - Cphi(j-1,2)*(var(i,  k,jm3)-var(i,  k,ju3))*zdpv(i,  k,ju4))  & 
       / Cphi(j,  1)) / zdp5(i,k,jm5)
END FUNCTION diffuse_4th_order
!
!
!
REAL FUNCTION diffuse_horizontal_sigma_3d1(var,i,jm1,j,jp1,k)
  IMPLICIT NONE
  REAL, INTENT(IN) :: var(:,:,:)
  INTEGER, INTENT(IN) :: i,jm1,j,jp1,k
  diffuse_horizontal_sigma_3d1 =                                         &
       (              (var(i+1,j  ,k)-var(i,  j  ,k))*zdpu(i,  k,ju3)    &
      -               (var(i,  j  ,k)-var(i-1,j  ,k))*zdpu(i-1,k,ju3)    & 
      + ( Cphi(j  ,2)*(var(i,  jp1,k)-var(i,  j  ,k))*zdpv(i,  k,ju4)    &
        - Cphi(jm1,2)*(var(i,  j  ,k)-var(i,  jm1,k))*zdpv(i,  k,js4))   &
      /   Cphi(j  ,1)) / zdp5(i,k,j)
END FUNCTION diffuse_horizontal_sigma_3d1
!
!
!
REAL FUNCTION diffuse_qd1(var,i,j,k,nt,j1,j2,j3)
  IMPLICIT NONE
  REAL,    INTENT(IN) :: var(:,:,:)
  INTEGER, INTENT(IN) :: i,j,k,nt, j1,j2,j3
  ! GEWICHTETE HDIFF
  diffuse_qd1 = -Vvfh(k)*aks4*    &
 (                     (var(i+1,k,jm3)-var(i,  k,jm3))*zdpu(i,  k,jm3)   &
    -                  (var(i,  k,jm3)-var(i-1,k,jm3))*zdpu(i-1,k,jm3)   &
    + (  Cphi(j,  2) * (var(i,  k,jo3)-var(i,  k,jm3))*zdpv(i,  k,jm4)   &
        -Cphi(j-1,2) * (var(i,  k,jm3)-var(i,  k,ju3))*zdpv(i,  k,ju4))  &
       / Cphi(j  ,1)) / zdp5(i,k,jm5)
END FUNCTION diffuse_qd1


!
!
!> explizit horizontal advection 
!
REAL FUNCTION advect_horizontal(var,i,j,k,nt)
  IMPLICIT NONE
  REAL,    INTENT(IN) :: var(:,:,:,:)
  INTEGER, INTENT(IN) :: i,j,k,nt
  REAL :: zgvco, zgvcu, zt1, zt2, zt3, zt4
  REAL :: zfadvx, zfadvy
  !
  ZFADVX= .5*ACPHIR(J,1)*EDDLAM
  ZFADVY= .5*ACPHIR(J,1)*EDDPHI
  !
  zgvco = ZGV(i,k,jm3)*CPHI(j  ,2)
  zgvcu = ZGV(i,k,ju3)*CPHI(j-1,2)
  !
  zt1 = ZGU(i,  k,jm2) * (var(i+1,  j,k,nnow)-var(i,    j,k,nnow))
  zt2 = ZGU(i-1,k,jm2) * (var(i,    j,k,nnow)-var(i-1,  j,k,nnow))
  zt3 =         zgvco  * (var(i,  j+1,k,nnow)-var(i,    j,k,nnow))
  zt4 =         zgvcu  * (var(i,    j,k,nnow)-var(i,  j-1,k,nnow))
  !
  advect_horizontal = -(zfadvx*(zt1+zt2) + zfadvy*(zt3+zt4)) * ZEDDPQ(i,k)
  !
END FUNCTION advect_horizontal


REAL FUNCTION alpha_omega(i,j,k)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j,k
  REAL :: za1, za2, za3, zgvco, zgvcu
  zgvco = ZGV(i,k,jm3)*CPHI(j,2)
  zgvcu = ZGV(i,k,ju3)*CPHI(j-1,2)
  za1   = ZGU(i,k,jm2)*ZPLAM(i,k) + ZGU(i-1,k,jm2)*ZPLAM(i-1,k)
  za2   = (zgvco*ZPPHI(i,k,jo2)+zgvcu*ZPPHI(i,k,jm2))*zx2
  !LM HERE ZEDDPQINT IS NOW USED INSTEAD OF ZEDDPQ
  !LM AND ZA3 IS MULTIPLIED WITH DWDT (SEE GOETTEL 4.1.22)
  za3 = -r*ZTV(i,k,jm5)*ZEDDPQINT(i,k)*DWDT(i,j,k,nnow)                                           &
      & *((ZALOPN(i,k+1,jm3)-ZALOPN(i,k,jm3))*ZSDIV(i,k,jm2)+(ZSDIV(i,k+1,jm2)-ZSDIV(i,k,jm2))  &
      & *ZBETAK(i,k,jm3))
  alpha_omega = (.5*(za1+za2)*ZEDDPQ(i,k)+za3)
END FUNCTION alpha_omega
 
SUBROUTINE dyn_swap_indices
  ! ZYKLISCHE VERTAUSCHUNG DER SCHEIBENINDIZES
  IMPLICIT NONE
  JM2 = 3 - JM2; JO2 = 3 - JO2
  JSP = JU3; JU3 = JM3; JM3 = JO3; JO3 = JSP
  JSP = JS4; JS4 = JU4; JU4 = JM4; JM4 = JO4; JO4 = JSP
  JSP = JS5; JS5 = JU5; JU5 = JM5; JM5 = JO5; JO5 = JN5; JN5 = JSP
END SUBROUTINE dyn_swap_indices

SUBROUTINE dyn_init_indices
  ! ANFANGSINDIZES FUER ZYKLISCHES UMSPEICHERN SETZEN
  IMPLICIT NONE
  JM2 = 1; JO2 = 2
  JU3 = 1; JM3 = 2; JO3 = 3
  JS4 = 1; JU4 = 2; JM4 = 3; JO4 = 4
  JS5 = 1; JU5 = 2; JM5 = 3; JO5 = 4; JN5 = 5
END SUBROUTINE dyn_init_indices

LOGICAL FUNCTION at_boundary(i,j,k)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j,k
  LOGICAL :: boundary_check(4)
  at_boundary = .FALSE.
  boundary_check(1) = (i==iah .AND. neighbor(1)==-1)
  boundary_check(2) = (j==jeh .AND. neighbor(2)==-1)
  boundary_check(3) = (i==ieh .AND. neighbor(3)==-1)
  boundary_check(4) = (j==jah .AND. neighbor(4)==-1)
  IF(ANY(boundary_check)) at_boundary = .TRUE.
END FUNCTION at_boundary


 
END MODULE MO_DYNAMICS
