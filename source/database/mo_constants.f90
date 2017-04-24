!
!C**   BESCHREIBUNG DER VARIABLEN:
!C**   EDDLAM   :   INVERSE MASCHENWEITE IN LAM-RICHTUNG      (1/RAD)
!C**   EDDPHI   :   INVERSE MASCHENWEITE IN PHI-RICHTUNG      (1/RAD)
!C**   EDADPHI  :   INVERS.(MASCHENWEITE*RERD) IN PHI-RICHTUNG(1/RAD*M)
!C**   DLADDPH  :   DLAM/DPHI                                 (--)
!C**   DPHDDLA  :   DPHI/DLAM                                 (--)
!C**   PI       :   ZAHL PI                                   (--)
!C**   ED2DT    :   1/(2*ZEITSCHRITT)                         (1/S)
!C**   DT2      :   2*ZEITSCHRITT                             (S)
!C**   DTDEH    :   ZEITSCHRITT/3600 S                        (--)
!C**   NEHDDT   :   3600 S/ZEITSCHRITT                        (--)
!C**   GQ       :   QUADRAT DER ERDBESCHLEUNIGUNG             (M/S**2)**2
!C**   GH       :   HAELFTE DER ERDBESCHLEUNIGUNG : G/2       (M/S**2)
!C**   EDG      :   INVERSES DER ERDBESCHLEUNIGUNG            (S**2)/M)
!C**   DEGRAD   :   UMRECHNUNGSFAKTOR DEGREE IN RAD           (--)
!C**   RADDEG   :   UMRECHNUNGSFAKTOR RAD IN DEGREE           (--)
!C**   RDRD     :   R/RD                                      (--)
!C**   EMRDRD   :   1 - R/RD                                  (--)
!C**   RDDRM1   :   RD/R - 1                                  (--)
!C**   TCRITH   :   273.16                                    (K)
!C**   TCRITL   :   268.16                                    (K)
!C**   WCPR     :   1/WCP                                     (KG*K/J)
!C**   RDWCP    :   R/WCP                                     (--)
!C**   B234W    :   B2W*(B3 - B4W)                            (K**2)
!C**   B234E    :   B2E*(B3 - B4E)                            (K**2)
!C**
!
module mo_constants 
   !
   implicit none
   !
   real, parameter      ::   speed_of_light     = 2.99792458e0
   real, parameter      ::   Boltzmann_constant = 1.3806488e-16
   !real, parameter      :: gas_constant       = 5.189e19
   real, parameter      ::   gas_constant       = 8.3145119843E7
   !
   real, parameter      ::   r_earth  = 6371.229E3
   real, parameter      ::   PI       =   4.0*ATAN(1.0)
   real, parameter      ::   T0       =   273.16
   real, parameter      ::   R        =   287.05
   real, parameter      ::   RD       =   461.51
   real, parameter      ::   WCP      =  1005.0
   real, parameter      ::   WLK      =     2.501E6
   real, parameter      ::   WLF      =     0.334E6
   real, parameter      ::   WLS      =     2.835E6
   real, parameter      ::   G        =     9.80665
   real, parameter      ::   RERD     =  6371.229E3
   real, parameter      ::   STAG     = 86164.09054
   real, parameter      ::   RHF      =  1000.0
   real, parameter      ::   SIGMA    =     5.6697E-8
   real, parameter      ::   SOKO     =  1368.0
   !
   real, parameter      ::   B1    =   610.78
   real, parameter      ::   B2W   =    17.2693882
   real, parameter      ::   B2E   =    21.8745584
   real, parameter      ::   B3    =   273.16
   real, parameter      ::   B4W   =    35.86
   real, parameter      ::   B4E   =     7.66

   real, parameter      ::   UC1   =     0.8
   real, parameter      ::   UC2   =  SQRT(3.0)
   real, parameter      ::   UCL   =     1.00

   real, parameter      ::   RHDE  =     1.0/5.0*1000.0
   real, parameter      ::   AKT   =     0.4

   real, parameter      ::   GQ      =   G**2
   real, parameter      ::   GH      =   G*0.5
   real, parameter      ::   EDG     =   1.0/G
   real, parameter      ::   DEGRAD  =   PI/180.0
   real, parameter      ::   RADDEG  =   180.0/PI
   real, parameter      ::   RDRD    =   R/RD
   real, parameter      ::   EMRDRD  =   1.0 - RDRD
   real, parameter      ::   RDDRM1  =   RD/R - 1.0
   real, parameter      ::   WCPR    =   1.0/WCP
   real, parameter      ::   RDWCP   =   R  /WCP
   real, parameter      ::   B234W   =   B2W*(B3 - B4W)
   real, parameter      ::   B234E   =   B2E*(B3 - B4E)
   !
   real, parameter      ::   TCRITH  =   273.16
   real, parameter      ::   TCRITL  =   268.16
   !
end module mo_constants 
