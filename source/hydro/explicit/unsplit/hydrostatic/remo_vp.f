      ! remo for the "vertical plane"
      ! author: lennart.marien@hzg.de
      PROGRAM REMOVP
      USE netcdf
      USE VARDEF
      USE IO
      IMPLICIT NONE

      ! read namelist input
      CALL READ_NML
      ! set some derived parameters
      KE1 = KE + 1
      DPHI = DLAM
      NANF = 0
      ! first output timestep
      NOUT = CEILING((STOUT - DT/2.)/DT)
      ! last timestep
      NENDE = CEILING((TENDE - DT/2.)/DT + 2.)

      ! allocate temporary and output fields
      CALL ALLOCATIONS
      ! retrieve aks and bks from file
      CALL READ_AKBK
      ! now initialize various constants and control parameters
      CALL GKONST
      
      ! call test case initializer according to given namelist
      ! add more initializers at the end of this file and
      ! extend this conditional to add more cases
      IF (LSTRAKA) THEN
        CALL INIT_STRAKA_BUBBLE
      ELSE IF (LMOUNTWAV) THEN
        CALL INIT_MOUNTAIN_WAVE
      ELSE IF (LGRAVWAVE) THEN
        CALL INIT_GRAVWAVE
      ELSE
        STOP "ERROR! PLEASE DEFINE TESTCASE."
      END IF

      TG(:,1) = T(:,KE,1)

      ! make sure all initial timesteps are equal 
      U(:,:,2) = U(:,:,1)
      U(:,:,3) = U(:,:,1)
      W(:,:,2) = W(:,:,1)
      W(:,:,3) = W(:,:,1)
      PINT(:,:,2) = PINT(:,:,1)
      PINT(:,:,3) = PINT(:,:,1)
      T(:,:,2) = T(:,:,1)
      T(:,:,3) = T(:,:,1)
      PS(:,2) = PS(:,1)
      PS(:,3) = PS(:,1)
      TG(:,2) = TG(:,1)
      TG(:,3) = TG(:,1)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !
      NZT = NANF
      NZTANF = NANF
      IF (NANF>0) NZTANF = NANF + 1
      !
 
      ! create output file
      CALL NC_CREATE
      ! create dimensions
      CALL NC_DEFDIM
      ! create coordinate and output variables
      CALL NC_DEFVAR
         
C     ZEITSCHRITT HALBIEREN, WENN NZT = 0 (VORWAERTSZEITSCHRITT)
C
      IF (NZT.EQ.0) THEN
          DT = DT*0.5
      ENDIF
      
      ! main loop starts here
      DO WHILE (NZT+1.LE.NENDE)

      DT2 = 2.0*DT
      ED2DT = 1.0/DT2

C     ZYKLISCHE VERTAUSCHUNG DER ZEIT-INDIZES
C
      NSP = NA
      NA  = NJ
      NJ  = NE
      NE  = NSP
      
      ZFDIVX = ACPHIR(1)*EDDLAM
      ZX1 = .5*R*ACPHIR(1)*EDDLAM
        
      ! compute half-level pressures
      CALL CALC_P_HL(NJ)
      IF (LHYDRO) THEN
         ! in contrast to standard REMO ZPNF is based on the hydrostatic
         ! pressure but this is just a cosmetic change because where it
         ! should be the non-hydrostatic pressure PINT is used directly
         ! in this version instead
         PINT(:,:,NJ) = ZPNF(:,:)
         DWDT(:,:,NJ) = 1.
      END IF
      ! compute logarithm of half-level pressures
      CALL CALC_ALOGP_HL
      ! compute pressure thicknesses
      CALL CALC_DP(NA)
      CALL CALC_DP(NJ)
      ! compute zbetak
      CALL CALC_BETAK
      ! compute logarithm of full level pressures
      CALL CALC_ALOGP_FL(NJ)
      ! compute full-level pressures by averaging (for diffusion)
      CALL CALC_P_FL_DIF
      ! compute full-level geopotential (for diffusion)
      ! by calling this before computing the geopotential
      ! for time step NJ the last timestep will be used for
      ! diffusion just as it is right now in REMO
      CALL CALC_FI_FL_DIF
      ! compute geopotential on half levels
      ! in REMO "FI" is updated at the end of the timestep
      ! instead but this should make no difference
      CALL CALC_FI_HL(NJ)
      ! compute geopotential on full levels
      CALL CALC_FI_FL

C     ZTKVZ , ZTKVL : GEWICHTSFAKTOREN ZUR HORIZONTALEN MITTELUNG DER
C                     VERTIKALEN TRANSPORTKOEFFIZIENTEN
      ZTKVZ = 0.9
      ZTKVL = (1.-ZTKVZ)*0.5
C
C     EDDY DIFFUSION
C
      DO K = KE,2,-1
        DO I = IAA, IEA
          KM1 = K - 1
          ZDPN2     = ZDP(I,K,NA) + ZDP(I,KM1,NA)
          ZTN       = (ZDP(I,K,NA)*T(I,KM1,NA) +
     &         ZDP(I,KM1,NA)*T(I,K,NA))/ZDPN2
          HYDROP    = AK(K) + BK(K)*(PS(I,NA)-PTOP)
          ZRHO      = HYDROP/(R*ZTN)
          ZGROP     = 2.0*G*G*ZRHO*ZRHO/ZDPN2
          ZTMKVM(I,K) = Z2ZM(I,KM1)*ZGROP
          ZTMKVH(I,K) = Z2ZH(I,KM1)*ZGROP
        ENDDO
      ENDDO
      
      DO K = 1, KE
         DO I = IAA, IEA
            ZTPA(I,K) = T(I,K,NA) + ZFIHF(I,K)*WCPR
         END DO

         ! "mass" variable
         DO I = IAA , IEH + 1
            ZGU(I,K) = .5*(ZDP(I,K,NJ)+
     &           ZDP(I+1,K,NJ))*U(I,K,NJ)
            ZPLAM(I,K)   = ZX1*( T(I+1,K,NJ) +    T(I,K,NJ) )
     &           *(ZALOPH(I+1,K) - ZALOPH(I,K))
            ZDPU(I,K) = 0.5*(ZDP(I,K,NA) + ZDP(I+1,K,NA))
         END DO
 
         ! vertical sum of divergence
         DO I = IAH - 1, IEH + 1
            ZSDIV(I,K+1) = ZSDIV(I,K) +
     &           ZFDIVX*( ZGU(I,K)-ZGU(I-1,K) )
            ZEKIN(I,K) = 
     &           .25*(  U(I-1,K,NJ)**2 + U(I,K,NJ)**2 )
         ENDDO
      ENDDO

      DO K = 1, KE
        DO I = IAA, IEA
          TPOT(I,K) =  T(I,K,NJ)*(P0/EXP(ZALOPH(I,K)))**(R/CP)
          TPOTPERT(I,K) = TPOT(I,K) - TPOTINI(I,K)
        END DO
      END DO

      DO K = 1,KE
         DO I = IAH - 1, IEH + 1
            ZLAPU(I,K) = U(I+1,K,NA) + U(I-1,K,NA) - 2.0*U(I,K,NA)
            ZLAPW(I,K) = W(I+1,K,NA) + W(I-1,K,NA) - 2.0*W(I,K,NA)
            ZLAPDWDT(I,K) = DWDT(I+1,K,NA) + DWDT(I-1,K,NA) - 
     &                      2.0*DWDT(I,K,NA)
!            ZLAPT(I,K) = 
!     &             ((ZTPA(I+1,K)-ZTPA(I  ,K))*ZDPU(I  ,K)-
!     &              (ZTPA(I  ,K)-ZTPA(I-1,K))*ZDPU(I-1,K))/
!     &             ZDP(I,K,NA)
            ! right now no pressure weighting of temperature diffusion
            ! because it's unclear whether it should be based on PINT or
            ! not
            ZLAPT(I,K) = T(I+1,K,NA) + T(I-1,K,NA) - 2.0*T(I,K,NA)

         END DO
      END DO

      DO K = 1, KE
        DO I = IAH, IEH
!         Make sure to use correct diffusion for periodic boundaries
          IF ( LHDIFF2
     &         .OR.(I.EQ.IAH.AND..NOT.LPERIODIC)
     &         .OR.(I.EQ.IEH.AND..NOT.LPERIODIC)) THEN
              ZUDIFH(I,K) = AKS2*ZLAPU(I,K)
              ZTDIFH(I,K) = AKS2*ZLAPT(I,K)
            ELSE
              ZUDIFH(I,K) = -AKS4*( ZLAPU(I+1,K) + ZLAPU(I-1,K) -
     &                    2.0*ZLAPU(I,K) )
              ZTDIFH(I,K) = -AKS4*( ZLAPT(I+1,K) + ZLAPT(I-1,K) -
     &                    2.0*ZLAPT(I,K) )
            END IF
         END DO
      END DO

      DO K = 2 , KE
        DO I = IAH - 1, IEH + 1
          ZETAS(I,K) = BK(K)*ZSDIV(I,KE1) - ZSDIV(I,K)
        ENDDO
      ENDDO
  
      ZFADVX= .5*ACPHIR(1)*EDDLAM

C     EXPLIZITE PROGNOSE OHNE RANDRELAXATIONSTERME
C     --------------------------------------------
C     1. BODENDRUCK + DRUCK
C     ---------------------

      DO I = IAH, IEH
        PSDT(I) = ZSDIV(I,KE1)
        PS(I,NE) = PS(I,NA) - PSDT(I)*DT2
        PINT(I,1,NE) = PTOP
      ENDDO

      DO K = 2, KE1
         KM1 = K - 1
         DO I = IAH, IEH
          PINT(I,K,NE) = -(BK(KM1) + BK(K)) * DT2 * PSDT(I) *
     &            DWDT(I,KM1,NJ) +
     &            PINT(I,KM1,NA) + PINT(I,K,NA) - 
     &            PINT(I,KM1,NE)
         END DO
      END DO

C     2. TEMPERATUR
C     ---------------------

C     HORIZONTALADVEKTION UND QUELLTERME
      DO K = 1 , KE
         DO I = IAH , IEH
            ZEDDPQ(I,K) = 1./ZDP(I,K,NJ)
            ZT1         = ZGU(I  ,K)*
     &              ( T(I+1,K,NJ) - T(I  ,K,NJ) )
            ZT2         = ZGU(I-1,K)*
     &              ( T(I  ,K,NJ) - T(I-1,K,NJ) )
            ZTADV(I,K)  = -( ZFADVX*(ZT1+ZT2) )
     &              *ZEDDPQ(I,K)

C           ALPHA*OMEGA
            ZA1 =  ZGU(I  ,K)*ZPLAM(I  ,K) +
     &           ZGU(I-1,K)*ZPLAM(I-1,K)
            ZA3 = -2*R*T(I,K,NJ)*(1/(PINT(I,K,NJ) + PINT(I,K+1,NJ))
     &            * DWDT(I,K,NJ) * 0.5 * (ZSDIV (I,K)
     &         + ZSDIV (I,K+1) ) )

            ZALPOM(I,K) = ( .5*( ZA1 )*ZEDDPQ(I,K) + ZA3 )

            AGB(I,K) = ED2DT
            AGD(I,K) = ED2DT*T(I,K,NA) + ZTADV(I,K) +
     &              ZTDIFH(I,K) + ZALPOM(I,K) *WCPR
         ENDDO
      ENDDO

C     VERTIKALADVEKTION UND -DIFFUSION
      DO I = IAH , IEH
         ZAGCM    = .5*ZETAS(I,2)*ZEDDPQ(I,1)
         ZTMKVHM  = ZTKVZ*ZTMKVH(I  ,2) + ZTKVL*
     &        ( ZTMKVH(I+1,2) + ZTMKVH(I-1,2) )
         ZAGCT      = - ZTMKVHM*ZEDDPQ(I,1)
         AGC(I,1) = ZAGCM*ZA1A + ZAGCT*A1T(2)
         AGB(I,1) = AGB(I,1) - AGC(I,1)
         AGD(I,1) = AGD(I,1) - ( ZA2A*ZAGCM + A2T(2)*ZAGCT )*
     &        ( T(I,2,NA) - T(I,1,NA) )
     &        -WCPR* ZAGCT*( ZFIHF(I,2) - ZFIHF(I,1) )
      ENDDO
!     Equations from Part II of EM/DM
      DO K = 2 , KE-1
         DO I  = IAH , IEH
!           1.RHS Eq. 5.4.5
            ZAGAM      = -.5*ZETAS(I,K)*ZEDDPQ(I,K)
!           1.RHS Eq. 5.4.6
            ZAGCM      =  .5*ZETAS(I,K+1)*ZEDDPQ(I,K)
!           2.RHS Eq. 5.4.5
            ZTMKVHM    = ZTKVZ*ZTMKVH(I  ,K) + ZTKVL*
     &           ( ZTMKVH(I+1,K) + ZTMKVH(I-1,K) )
            ZAGAT      = - ZTMKVHM*ZEDDPQ(I,K)
!           2.RHS Eq. 5.4.6
            ZTMKVHM    = ZTKVZ*ZTMKVH(I  ,K+1) + ZTKVL*
     &           ( ZTMKVH(I+1,K+1) + ZTMKVH(I-1,K+1) )
            ZAGCT      = - ZTMKVHM*ZEDDPQ(I,K)
!           Eq. 5.4.5
            AGA(I,K) = ZAGAM*ZA1A + ZAGAT*A1T(K)
!           Eq. 5.4.6
            AGC(I,K) = ZAGCM*ZA1A + ZAGCT*A1T(K+1)
!           Eq. 5.4.7
            AGB(I,K) = AGB(I,K) - AGA(I,K) - AGC(I,K)
!           Eq. 5.4.8
            AGD(I,K) = AGD(I,K) - ( ZA2A*ZAGAM + A2T(K)*ZAGAT )*
     &              ( T(I,K-1,NA) - T(I,K,NA) )
     &              - ( ZA2A*ZAGCM + A2T(K+1)*ZAGCT )*
     &              ( T(I,K+1,NA) - T(I,K,NA) )
     &              -WCPR* ZAGAT*( ZFIHF(I,K-1) - ZFIHF(I,K))
     &              -WCPR* ZAGCT*( ZFIHF(I,K+1) - ZFIHF(I,K))
         ENDDO
      ENDDO
C
      DO I = IAH , IEH
         ZAGAM     = -.5*ZETAS(I,KE)*ZEDDPQ(I,KE)
         ZTMKVHM   = ZTKVZ*ZTMKVH(I  ,KE) + ZTKVL*
     &           ( ZTMKVH(I+1,KE) + ZTMKVH(I-1,KE) )
         ZAGAT     = -  ZTMKVHM *ZEDDPQ(I,KE)
         ZAGCT     = - TMCH(I)*ZEDDPQ(I,KE)
         AGA(I,KE) = ZAGAM*ZA1A + ZAGAT*A1T(KE)
         AGB(I,KE) = AGB(I,KE) - AGA(I,KE) - ZAGCT*A1T(KE1)
         AGD(I,KE) = AGD(I,KE) - ( ZA2A*ZAGAM + A2T(KE)*ZAGAT)*
     &           ( T(I,KE-1,NA) - T(I,KE,NA) )
     &           - A2T(KE1)*ZAGCT*( TG(I,NA) - T(I,KE,NA) )
     &           - A1T(KE1)*ZAGCT*  TG(I,NA)
     &           -WCPR* ZAGAT*( ZFIHF(I,KE-1) - ZFIHF(I,KE) )
     &           -WCPR* ZAGCT*( FIB  (I)      - ZFIHF(I,KE) )
      ENDDO
      
C     GAUSS - ELIMINATION UND H-QDW-ZERLEGUNG

      CALL GAUSS ( IAH , IEH ,AGA,AGB,AGC,AGD,AGE,IE,KE)

      DO K = 1 , KE
         DO I = IAH , IEH
            T (I,K,NE) = AGE(I,K)
         ENDDO
      ENDDO

C     3. U - PROGNOSE
C        ---------------

      ZX1   = ACPHIR(1)*EDDLAM

      DO K = 1 , KE
         DO I = IAU , IEU
            ZGRAD(I,K) = -ZPLAM(I,K) - ZX1 * 0.5 * (DWDT(I,K,NJ)
     &             + DWDT(I+1,K,NJ)) * ( (ZFIH(I+1,K) - ZFIH(I,K)) 
     &             + ZEKIN(I+1,K) - ZEKIN(I,K) )
            ZEDDPQ(I,K) =  2./( ZDP(I+1,K,NJ) + ZDP(I,K,NJ) )
            AGB(I,K) = ED2DT
            AGD(I,K) = ED2DT*U(I,K,NA) + ZGRAD(I,K) + ZUDIFH(I,K)
         ENDDO
      ENDDO

C     VERTIKALDIFFUSION UND VERTIKALADVEKTION
      DO I = IAU , IEU
         ZAGCM = .25*( ZETAS (I,2) + ZETAS(I+1,2) )*
     &       ZEDDPQ(I,1)
         ZAGCT = -0.5*(ZTMKVM(I,2)+ ZTMKVM(I+1,2) )*
     &       ZEDDPQ(I,1)
         AGC(I,1) = ZAGCM*ZA1A + ZAGCT*A1T(2)
         AGB(I,1) = AGB(I,1) - AGC(I,1)
         AGD(I,1) = AGD(I,1) - ( ZA2A*ZAGCM + A2T(2)*ZAGCT )*
     &           ( U(I,2,NA) - U(I,1,NA ) )
      ENDDO

      DO K = 2 , KE-1
         DO I = IAU , IEU
            ZAGAM = - 0.25*(ZETAS(I,K)+ZETAS(I+1,K))*
     &              ZEDDPQ(I,K)
            ZAGCM =   0.25*(ZETAS(I,K+1)+ZETAS(I+1,K+1))*
     &              ZEDDPQ(I,K)
            ZAGAT = - 0.5 *(ZTMKVM(I,K  )+ZTMKVM(I+1,K  ))*
     &              ZEDDPQ(I,K)
            ZAGCT = - 0.5 *(ZTMKVM(I,K+1)+ZTMKVM(I+1,K+1))*
     &              ZEDDPQ(I,K)
            AGA(I,K) = ZAGAM*ZA1A + ZAGAT*A1T(K)
            AGC(I,K) = ZAGCM*ZA1A + ZAGCT*A1T(K+1)
            AGB(I,K) = AGB(I,K) - AGA(I,K) - AGC(I,K)
            AGD(I,K) = AGD(I,K) - ( ZA2A*ZAGAM + A2T(K)*ZAGAT )*
     &              ( U(I,K-1,NA) - U(I,K,NA) )
     &              - ( ZA2A*ZAGCM + A2T(K+1)*ZAGCT )*
     &              ( U(I,K+1,NA) - U(I,K,NA) )
         ENDDO
      ENDDO

      DO I = IAU , IEU
         ZAGAM = - 0.25*(ZETAS  (I,KE)+ZETAS (I+1,KE))*
     &           ZEDDPQ(I,KE)
         ZAGAT = - 0.5 *( ZTMKVM(I,KE)+ZTMKVM(I+1,KE))*
     &           ZEDDPQ(I,KE)
         ZAGCT = - 0.5 *( TMCM (I) + TMCM (I+1) )    *
     &           ZEDDPQ(I,KE)
         AGA(I,KE) = ZAGAM*ZA1A + ZAGAT*A1T(KE)
         AGB(I,KE) = AGB(I,KE) - AGA(I,KE) - ZAGCT*A1T(KE1)
         AGD(I,KE) = AGD(I,KE) - ( ZA2A*ZAGAM + A2T(KE)*ZAGAT )*
     &        ( U(I,KE-1,NA) - U(I,KE,NA) )
     &        + A2T(KE1)*ZAGCT*U(I,KE,NA)
      ENDDO

C     GAUSS-ELIMINATION UND U-ENDE-WERTE AUF U(I,J,K,NE) ABLEGEN
      CALL GAUSS ( IAU , IEU , AGA,AGB,AGC,AGD,AGE,IE,KE)

      DO K = 1 , KE
         DO I = IAU , IEU
            U(I,K,NE) = AGE(I,K)
         ENDDO
      ENDDO
!     RAYLEIGH DAMPING
      IF ( LRAYDAMP ) THEN
        CALL RAYLEIGH(1)
      ENDIF
C     PERIODIC BOUNDARY CONDITIONS
      IF ( LPERIODIC ) THEN
        U    (1    :IAU-1,:,NE) = U    (IEU-2:IEU-1,:,NE)
        U    (IEU+1:IE   ,:,NE) = U    (IAU+1:IAU+2,:,NE)
        T    (1    :IAH-1,:,NE) = T    (IEH-2:IEH-1,:,NE)
        T    (IEH+1:IE   ,:,NE) = T    (IAH+1:IAH+2,:,NE)
        PS   (1    :IAH-1,  NE) = PS   (IEH-2:IEH-1,  NE)
        PS   (IEH+1:IE   ,  NE) = PS   (IAH+1:IAH+2,  NE)
        PINT (1    :IAH-1,:,NE) = PINT (IEH-2:IEH-1,:,NE)
        PINT (IEH+1:IE   ,:,NE) = PINT (IAH+1:IAH+2,:,NE)
        ZETAS(1    :IAH-1,:   ) = ZETAS(IEH-2:IEH-1,:   )
        ZETAS(IEH+1:IE   ,:   ) = ZETAS(IAH+1:IAH+2,:   )
      ENDIF

      IF (.NOT.LHYDRO) THEN
        CALL CALC_P_HL(NE)
        CALL CALC_DP(NE)
         CALL CDZDT
!       RAYLEIGH DAMPING
        IF ( LRAYDAMP ) THEN
          CALL RAYLEIGH(2)
        ENDIF
!       Periodic boundary conditions
        IF ( LPERIODIC ) THEN
          W    (1    :IAH-1,:,NE) = W    (IEH-2:IEH-1,:,NE)
          W    (IEH+1:IE   ,:,NE) = W    (IAH+1:IAH+2,:,NE)
        END IF
        IF (NZT > 3) THEN
            CALL CDWDT
          CALL FINAL
        END IF
        DO K = 1, KE
          DO I = IAH, IEH
            DWDT(I,K,NE) = (PINT(I,K+1,NE) - PINT(I,K,NE))/ZDP(I,K,NE)
          END DO
        END DO
        IF ( LPERIODIC ) THEN
          U    (1    :IAU-1,:,NE) = U    (IEU-2:IEU-1,:,NE)
          U    (IEU+1:IE   ,:,NE) = U    (IAU+1:IAU+2,:,NE)
          T    (1    :IAH-1,:,NE) = T    (IEH-2:IEH-1,:,NE)
          T    (IEH+1:IE   ,:,NE) = T    (IAH+1:IAH+2,:,NE)
          PS   (1    :IAH-1,  NE) = PS   (IEH-2:IEH-1,  NE)
          PS   (IEH+1:IE   ,  NE) = PS   (IAH+1:IAH+2,  NE)
          PINT (1    :IAH-1,:,NE) = PINT (IEH-2:IEH-1,:,NE)
          PINT (IEH+1:IE   ,:,NE) = PINT (IAH+1:IAH+2,:,NE)
          W    (1    :IAH-1,:,NE) = W    (IEH-2:IEH-1,:,NE)
          W    (IEH+1:IE   ,:,NE) = W    (IAH+1:IAH+2,:,NE)
          DWDT (1    :IAH-1,:,NE) = DWDT (IEH-2:IEH-1,:,NE)
          DWDT (IEH+1:IE   ,:,NE) = DWDT (IAH+1:IAH+2,:,NE)
        END IF
      END IF

      ! write output for this timestep
      IF ((NZT==0).OR.(NZT.EQ.NOUT+1)) THEN
          CALL NC_WRITE_TIMESTEP(NJ)
          ! compute next output time step
          NOUT = MIN(NENDE-2, CEILING((TAC + STOUT - DT)/DT))
      END IF

      ! asselin filtering
      DO I = 1, IE
         PS(I,NJ) = PS(I,NJ)
     &            + EPSASS*(PS(I,NE)- 2.*PS(I,NJ) + PS(I,NA))
         DO K = 1, KE
            U(I,K,NJ) = U(I,K,NJ)
     &                + EPSASS*(U(I,K,NE)- 2.*U(I,K,NJ) + U(I,K,NA))
            T(I,K,NJ) = T(I,K,NJ)
     &                + EPSASS*(T(I,K,NE)- 2.*T(I,K,NJ) + T(I,K,NA))
            DWDT(I,K,NJ) = DWDT(I,K,NJ)
     &                + EPSASS*(DWDT(I,K,NE)- 2.*DWDT(I,K,NJ) + 
     &                DWDT(I,K,NA))
         END DO
         DO K = 1, KE1
            W(I,K,NJ) = W(I,K,NJ)
     &                + EPSASS*(W(I,K,NE)- 2.*W(I,K,NJ) + W(I,K,NA))
            PINT(I,K,NJ) = PINT(I,K,NJ)
     &                + EPSASS*(PINT(I,K,NE)- 2.*PINT(I,K,NJ) + 
     &                PINT(I,K,NA))
         END DO
      END DO

C     NAECHSTEN ZEITSCHRITT BERECHNEN
C
      TAC = TAC + DT
      IF (NZT.EQ.0) THEN
        DT = DT*2.0
      ENDIF
      NZT    = NZT + 1

      ! main loop ends here
      ENDDO 
      ! write latest timestep NE, write time dimension, close file
      CALL NC_CLOSE
      ! deallocate fields
      CALL DEALLOCATIONS
      
      CONTAINS

        SUBROUTINE CALC_P_HL(NT)
          INTEGER, INTENT(IN) :: NT

          ZPNF(:,1) = PTOP
          DO K = 1, KE
             DO I = IAA, IEA
                ! based on PINT in real REMO but respective places have
                ! been adjusted so that it makes no difference 
                ! basically oversight and later laziness on my part
                ZPNF(I,K+1) = AK(K+1) + BK(K+1)*(PS(I,NT)-PTOP)
             END DO
          END DO
        END SUBROUTINE CALC_P_HL

        SUBROUTINE CALC_ALOGP_HL

          IF (LPTOP0) THEN
             ZALOPN(:,1) = ALOG(2.0)
          ELSE
             ZALOPN(:,1) = ALOG(PTOP)
          END IF

          DO K = 1, KE
            DO I = IAA, IEA
                ZALOPN(I,K+1) = ALOG(PINT(I,K+1,NJ))
             END DO
          END DO
        END SUBROUTINE CALC_ALOGP_HL

        SUBROUTINE CALC_DP(NT)
          INTEGER, INTENT(IN) :: NT

          DO K = 1, KE
             DO I = IAA, IEA
                ! note that this has a time dimension unlike ZDP in real
                ! REMO - there instead ZDP, ZDP5 and HYDRODP are used
                ZDP(I,K,NT) = DAK(K) + DBK(K)*(PS(I,NT)-PTOP) 
                ZDPINT(I,K) = PINT(I,K+1,NJ) - PINT(I,K,NJ)
             END DO
          END DO

        END SUBROUTINE CALC_DP

        SUBROUTINE CALC_BETAK

          DO I = IAA, IEA
             IF (LPTOP0) THEN
                ZBETAK(I,1) = ALOG(2.0)
             ELSE 
                ZBETAK(I,1) = 1. - PINT(I,1,NJ)/ZDPINT(I,1)
     &                    * (ZALOPN(I,2)-ZALOPN(I,1))
             END IF
          END DO

         DO K = 2, KE-1
            DO I = IAA , IEA
               ZBETAK(I,K) = 1. - PINT(I,K,NJ)/ZDPINT(I,K)
     &                       * (ZALOPN(I,K+1)-ZALOPN(I,K))
            END DO
         END DO
      
         DO I = IAA, IEA
            ZBETAK(I,KE) = 1. - PINT(I,KE,NJ) / ZDPINT(I,KE)
     &                   * (ZALOPN(I,KE1) - ZALOPN(I,KE))
         END DO

        END SUBROUTINE CALC_BETAK

        SUBROUTINE CALC_ALOGP_FL(NT)
          INTEGER, INTENT(IN) :: NT

          DO K = 1, KE
             KP1 = K + 1
             DO I = IAA , IEA
                ZALOPH(I,K) = ALOG(0.5*(PINT(I,K,NT) + PINT(I,KP1,NT)))
             END DO
          END DO

        END SUBROUTINE CALC_ALOGP_FL

        SUBROUTINE CALC_P_FL_DIF
      
          DO K = 1, KE
             DO I = IAA, IEA
                ZPHF(I,K) = 0.5*(ZPNF(I,K) + ZPNF(I,K+1))
             END DO
          END DO

        END SUBROUTINE CALC_P_FL_DIF

        SUBROUTINE CALC_T_FROM_TPOT(NT)
          INTEGER, INTENT(IN) :: NT

          DO I = IAA, IEA
             DO K = 1,KE
                T(I,K,NT) =  TPOT(I,K)*(EXP(ZALOPH(I,K))/P0)**(R/CP)
             END DO
          END DO

        END SUBROUTINE CALC_T_FROM_TPOT

        SUBROUTINE CALC_FI_HL(NT)          
          REAL :: DZ
          INTEGER, INTENT(IN)  :: NT

          DO I = IAA, IEA
             FI(I,KE1)    = FIB(I)
          END DO

          DO K = KE , 1  , -1
          KP1 = K + 1
            DO I = IAA , IEA
              DZ = T(I,K,NT)*R*2.
     &             /((PINT(I,KP1,NT) + PINT(I,K,NT)))
     &             * ZDP(I,K,NT)
              FI(I,K) = DZ + FI(I,KP1)
            END DO
          END DO
         
        END SUBROUTINE CALC_FI_HL

        SUBROUTINE CALC_FI_FL

          DO K = 1, KE
             KP1 = K + 1
             DO I = IAA , IEA
                ZFIH(I,K) = 0.5*(FI(I,K) + FI(I,KP1))
             END DO
          END DO
      
        END SUBROUTINE CALC_FI_FL

        SUBROUTINE CALC_FI_FL_DIF

          DO K = 1, KE
             KP1 = K + 1
             DO I = IAA , IEA
                ZFIHF(I,K) = 0.5*(FI(I,K) + FI(I,KP1))
             END DO
          END DO

        END SUBROUTINE CALC_FI_FL_DIF

        SUBROUTINE CDZDT
          REAL :: Z(IE,KE1), DZ, DZ2, ZPW(IE,KE), TTB(IE), ADVZ(IE,KE),
     &            VADV(IE,KE)
          INTEGER :: KP1

          Z(:,:) = 0.
          ZPW(:,:) = 0.
          TTB(:) = 0.
          ADVZ(:,:) = 0.
          VADV(:,:) = 0.
         !INITIALIZE W AND Z AT THE SURFACE. HERE W IS DEFINED ON HALF LEVELS!
          DO I = IAA, IEA
             Z(I,KE1)    = FIB(I)/G
             W(I,KE1,NE) = 0. 
          END DO

          DO K = KE , 1  , -1
          KP1 = K + 1
            DO I = IAA , IEA
              DZ = T(I,K,NE)*R*2.
     &             /((PINT(I,KP1,NE) + PINT(I,K,NE))*G)
     &             * ZDP(I,K,NE)
              DZ2 = T(I,K,NA)*R*2.
     &              /((PINT(I,KP1,NA) + PINT(I,K,NA))*G)
     &              * ZDP(I,K,NA)
              Z(I,K) = DZ + Z(I,KP1)
              W(I,K,NE) = W(I,KP1,NE) + DZ - DZ2
            END DO
        END DO


         !COMPUTE THE HEIGHT CHANGE (ZPW) AND THE LOCAL CHANGE (W) ON FULL LEVELS.
         !NOW W IS DEFINED ON FULL LEVELS!
         DO K = 1 , KE
            KP1 = K + 1
            DO I = IAA , IEA
               ZPW(I,K)  = (Z(I,K) + Z(I,KP1))*0.5
               W(I,K,NE) = (W(I,K,NE) + W(I,KP1,NE))*0.5/DT2
            END DO
         END DO

C--------------VERTICAL ADVECTION---------------------------------------
         DO I = IAH , IEH
            TTB(I) = 0.
         END DO

         DO K = 1 , KE-1
            DO I = IAH , IEH
               TTAL = (ZPW(I,K+1) - ZPW(I,K))*ZETAS(I,K+1)*0.5
               W(I,K,NE) = W(I,K,NE)
     &             + (TTAL + TTB(I))/ZDP(I,K,NE)
               VADV(I,K) = (TTAL + TTB(I))/ZDP(I,K,NE)
               TTB(I)    = TTAL
            END DO
         END DO

         DO I = IAH , IEH
            W(I,KE,NE) = 
     &           TTB(I)/ZDP(I,KE,NJ) + W(I,KE,NE)
            VADV(I,KE) = TTB(I)/ZDP(I,KE,NE)
         END DO

!     CALCULATES DIAGNOSTICALLY THE HORIZONTAL ADVECTION OF HEIGHT
         DO K = 1 , KE
            DO I = IAH , IEH
            ZT1         = ZGU(I,K)*
     &           ( ZPW(I+1,K) - ZPW(I,K) )
            ZT2         = ZGU(I-1,K)*
     &           ( ZPW(I,K) - ZPW(I-1,K) )
            ADVZ(I,K) = ( ZFADVX*(ZT1+ZT2) )
     &           /ZDP(I,K,NE)
            END DO
         END DO

         DO K = 1, KE
            DO I = IAH, IEH
               W(I,K,NE) = W(I,K,NE) + ADVZ(I,K)
            END DO
         END DO

        END SUBROUTINE CDZDT

        SUBROUTINE CDWDT
          REAL :: TDWDT(IE,KE), ADVW(IE,KE), VADV(IE,KE), TWB(IE), TWAL

          DO K = 1, KE
            DO I = IAA , IEA
               DWDT(I,K,NE) = 0.
               TDWDT(I,K) = (W(I,K,NE) - W(I,K,NA))/DT2
            END DO
          END DO

C--------------VERTICAL ADVECTION---------------------------------------
         TWB(:) = 0.

         DO K = 1 , KE-1
            KP1 = K + 1
            DO I = IAH , IEH
               TWAL = (W(I,KP1,NE) - W(I,K,NE))*ZETAS(I,KP1)*0.5
               DWDT(I,K,NE) = DWDT(I,K,NE) +
     &                  (TWAL + TWB(I))/ZDP(I,K,NE)
               VADV(I,K) = (TWAL + TWB(I))/ZDP(I,K,NE)
               TWB(I)    = TWAL
            END DO
         END DO

         DO I = IAH , IEH
            DWDT(I,KE,NE) = TWB(I)/ZDP(I,KE,NE) + DWDT(I,KE,NE)
            VADV(I,KE) = TWB(I)/ZDP(I,KE,NE)
         END DO

         DO K = 1 , KE
            DO I = IAH , IEH
            ZT1         = ZGU(I,K)*
     &           ( W(I+1,K,NE) - W(I,K,NE) )
            ZT2         = ZGU(I-1,K)*
     &           ( W(I,K,NE) - W(I-1,K,NE) )
            ADVW(I,K) = ( ZFADVX*(ZT1+ZT2) )
     &           /ZDP(I,K,NE)
            END DO
         END DO

         DO K = 1, KE
            DO I = IAH, IEH
               DWDT(I,K,NE) = DWDT(I,K,NE) + ADVW(I,K) + TDWDT(I,K)
            END DO
         END DO
!        Periodic boundary conditions
         IF ( LPERIODIC ) THEN
           DWDT (1    :IAH-1,:,NE) = DWDT (IEH-2:IEH-1,:,NE)
           DWDT (IEH+1:IE   ,:,NE) = DWDT (IAH+1:IAH+2,:,NE)
         ELSE
           DWDT (1,:,NE) = 1.0
           DWDT (IE,:,NE) = 1.0
         END IF
C        Spatial filter for DWDT
C        (horizontal)
         DO K = 1 , KE
           DO I = IAH , IEH
             DIFD(I,K) = DWDT(I-1,K,NE) - 2. * DWDT(I,K,NE)
     &            + DWDT(I+1,K,NE)
           END DO
         END DO
C        (vertical)
         DO K = 2 , KE-1
           DO I = IAH , IEH
             DIFD(I,K) = DWDT(I,K-1,NE) - 2. * DWDT(I,K,NE)
     &            + DWDT(I,K+1,NE) + DIFD(I,K)
           END DO
         END DO
C        (apply)
         DO K = 1 , KE
           DO I = IAH , IEH
             DWDT(I,K,NE) = DWDT(I,K,NE) + DDIFW * DIFD(I,K)
             END DO
         END DO
C        Time filter for DWDT
         DO K = 1, KE
            DO I = IAH, IEH
              DWDT(I,K,NE) = ((DWDT(I,K,NE)/G) + 1.) * DWP1 +
     &             DWDT(I,K,NA)*DWP
            END DO
         END DO
         
        END SUBROUTINE CDWDT
        
        SUBROUTINE FINAL
!         FINAL UPDATE OF NON-HYDROSTATIC CORRECTION

!
!         LOCAL VARIABLES
!
          REAL    :: CAPPA, GDT, GDT2, FCC, WGHT
          REAL    :: DPTL, DPTU, DELP, DPSTR
          REAL    :: PP1, RDPLDN, RDPLUP
          REAL    :: TTFC(KE)
          REAL    :: PDR(KE,IE), TEMP(KE,IE), PRES(KE1,IE), RTOP(KE,IE)
          REAL    :: WDOT(KE,IE), WVEL(KE1,IE), WRK(KE,IE)
          REAL ::
     &         B1(KE)   , B2(KE)   , B3(KE)   , C0(KE)    ,
     &         PONE(KE1), PSTR(KE1), PNP1(KE1), COFF(KE1) ,
     &         CHI(KE1)

          REAL ::  P1(KE1,IE)
          REAL :: DP, TMP

          CAPPA = R/CP
          GDT   = G*2.*DT
          GDT2  = GDT*GDT
          FCC   = -R/GDT2
CKS       FCC IS MULTIPLIED BY 4 IN ETA MODEL!

CKS       WGHT  = CAPPA
CKS       ETA MODEL AND WRF USE 0.35 INSTEAD OF CAPPA=287.05/1005.=0.2856!
          WGHT = 0.35

          DO I = IAH, IEH
             PDR(:,I) = ZDP(I,:,NE)
             TTFC(:) = 1. - CAPPA 
             RTOP(:,I) = T(I,:,NA)*CAPPA
     &                /(PINT(I,1:KE,NA) + PINT(I,2:KE1,NA))
             TEMP(:,I) = T(I,:,NE)
             PRES(:,I) = PINT(I,:,NE)
             WDOT(:,I) = DWDT(I,:,NE)
             WVEL(:,I) = W(I,:,NE)
             WRK(:,I)  = TTFC(:)*PDR(:,I)*FCC
          END DO

          DO I = IAH, IEH

!            UPPER BOUNDARY CONDITION
             CHI(1)  = 0.

             PONE(1) = AK(1) + BK(1) * (PS(I,NE)-PTOP)
             PSTR(1) = PONE(1)
             PNP1(1) = PONE(1)
             P1(1,I) = PONE(1)

             P1(2:KE1,I) = PRES(2:KE1,I)
             PONE(2:KE1) = PRES(2:KE1,I)
             DO K = 2 , KE1
                KM1 = K - 1
                DPSTR   = WDOT(KM1,I)*PDR(KM1,I)
                PSTR(K) = PSTR(KM1) + DPSTR
                PP1     = PNP1(KM1) + DPSTR
                DP      = (PP1 - PONE(K))*WGHT
                PNP1(K) = PONE(K) + DP
             END DO

CKS          WHY IS COFF ALWAYS MULTIPLIED WITH 0.5? IN WRF AND ETA
CKS          THIS IS NOT THE CASE!
CTR          MOVED FACTOR 0.5 TO THE DEFINITION OF COFF
CKS          Commented 0.5 out because it seems not make sense

             COFF(1:KE) = TEMP(:,I)*WRK(:,I)
     &                 /((PNP1(1:KE) + PNP1(2:KE1)))**2 !* 0.5

             DO K = 2 , KE
                KM1 = K - 1
                KP1 = K + 1
                RDPLDN = 1./PDR(K,I)
                RDPLUP = 1./PDR(KM1,I)

                B1(K) = COFF(KM1) + RDPLUP
                B2(K) = (COFF(KM1) + COFF(K)) - (RDPLUP + RDPLDN)
                B3(K) = COFF(K) + RDPLDN

                C0(K) = -((PSTR(KM1)+PSTR(K)-(PONE(KM1)+PONE(K)))*
     &                   COFF(KM1) + (PSTR(K)+PSTR(KP1) - 
     &                   (PONE(K)+PONE(KP1)))*COFF(K))
             END DO

             B2(KE) = B2(KE) + B3(KE)

             DO K = 3 , KE
                KM1 = K - 1
                TMP   = -B1(K)/B2(KM1)
                B2(K) = B3(KM1)*TMP + B2(K)
                C0(K) = C0(KM1)*TMP + C0(K)
             END DO

!            LOWER BOUNDARY CONDITION
             CHI(KE)  = C0(KE)/B2(KE)
             CHI(KE1) = CHI(KE)

             DO K = KE-1 , 2 , -1
                CHI(K) = (-B3(K)*CHI(K+1) + C0(K))/B2(K)
             END DO

             PNP1(1:KE1)   = CHI(1:KE1) + PSTR(1:KE1)
             PRES(1:KE1,I) = PNP1(1:KE1)
          END DO 

C----------BACKSUBSTITUTION---------------------------------------------
          DO I = IAH, IEH
             DPTU = 0.
             DO K = 1 , KE
                KP1  = K + 1
                DPTL = PRES(KP1,I) - P1(KP1,I)
                TEMP(K,I) = TEMP(K,I) + (DPTU + DPTL)*RTOP(K,I)
                DELP      = (PRES(KP1,I) - PRES(K,I))/PDR(K,I)
                WVEL(K,I) = WVEL(K,I) + (DELP - WDOT(K,I))*GDT
                WDOT(K,I)   = DELP 
                DPTU = DPTL
             END DO
          END DO

C-----------------------------------------------------------------------
C     offload working variables
C-----------------------------------------------------------------------
          DO I = IAH, IEH
             T(I,:,NE)    = TEMP(:,I)
             DWDT(I,:,NE) = WDOT(:,I)
             PINT(I,:,NE) = PRES(:,I)
             W(I,:,NE)    = WVEL(:,I)
          END DO

        END SUBROUTINE FINAL

      SUBROUTINE RAYLEIGH(key_ray)

      INTEGER :: key_ray
      
      REAL :: ZTOP, ZLEV, ZVFK
      REAL :: ZWGTVD(IE,KE)
      REAL :: ZURE, ZWR, DINIT, TIME
      
      ZWGTVD(:,:) = 0.
      DINIT = 0.
!
!     INCREASED RAYLEIGH DAMPING TO DAMP STARTUP TRANSIENTS
!
      IF (DINIT0.GT.0.) THEN
        TIME = NZT * DT
        DINIT=DINIT0*EXP(-TIME/TINIT)
        ZWGTVD(:,:) = DINIT/DT
      ENDIF
!
!     Damping at the top
!
      ZVFK = 1. / ( ZNRNEST*DT )  
      DO I = 1,IE
        ZTOP = ZFIH(I,1)/G

        DO K = 1,KE
          ZLEV = ZFIH(I,K)/G

          IF ( ZDAMP.GT.ZLEV ) EXIT
          ZWGTVD(I,K) = ZWGTVD(I,K) + ZVFK *
     &        (COS(0.5*PI*( ZTOP-ZLEV )/( ZTOP-ZDAMP )))**2
        ENDDO
      ENDDO
      
      SELECT CASE (key_ray)
      CASE(1)
        DO K = 1,KE
          DO I = IAA, IEA
            ZURE = UM
            U(I,K,NE)  = U(I,K,NE) -
     &           ZWGTVD(I,K)*( U(I,K,NE) - ZURE )
            U(I,K,NE)  = U(I,K,NE) -
     &           ALPHABOUND(I,2)*( U(I,K,NE) - ZURE )
          ENDDO
        ENDDO
      CASE(2)
        DO K = 1,KE
          DO I = IAA, IEA
            ZWR = 0.0
            W(I,K,NE) = W(I,K,NE) -
     &           ZWGTVD(I,K)*(W(I,K,NE) - ZWR)
            W(I,K,NE) = W(I,K,NE) -
     &           ALPHABOUND(I,1)*(W(I,K,NE) - ZWR)
          ENDDO
        ENDDO
      END SELECT

      END SUBROUTINE RAYLEIGH

      SUBROUTINE GETALPHA

      DO I = 1, IE
        ALPHABOUND(I,1) = RMY(I,1)/(1. + RMY(I,1))
        ALPHABOUND(I,2) = RMY(I,2)/(1. + RMY(I,2))
      END DO

      ALPHABOUND(   1:2 ,:) = 1.0
      ALPHABOUND(IE-1:IE,1) = 1.0
      ALPHABOUND(IE-2:IE,2) = 1.0
      
      END SUBROUTINE GETALPHA

      
      SUBROUTINE GKONST
!     initialize control parameters of time stepping
!     time indices
      NA = 3
      NJ = 1
      NE = 2
!     inner domain (without boundary zone) indices
      IAA = 1
      IAH = IAA + 2
      IEH = IE - 2
      IEA = IE
      IAU = IAA + 2
      IEU = IEA - 2

      DO I=IAA, IEA
        RLON(I) = LALU + (I-1) * DLAM
      END DO

!     initialize various constants
      P0      =  101325.
      PI      =  4.0*ATAN(1.0)
      EDDLAM  =  1.0/(DLAM*PI/180.0)
      G       =  9.80665
      R       =  287.05   
      CP      =  1005.0
      WCPR    =  1.0/CP
      RERD    =  6371.229 E3      
      ZPHIS   =  PHILU
      DEGRAD  =  PI/180.0
      CPHI    =  (/ COS(ZPHIS*DEGRAD), COS((ZPHIS + 0.5*DPHI)*
     &                DEGRAD) /)
      ACPHIR  =  (/ 1.0/(RERD*CPHI(1)), 1.0/(RERD*CPHI(2)) /)
      ZA1A    =  0.5
      ZA2A    =  (1.-ZA1A)
      TAC     =  0.
      AKS2MX  =  1.0/(16.0*DT)
      AKS4MX  =  1.0/(128.0*DT)
      DIFFAC  =  6.*DT/300.0*0.5/MAX(DLAM,DPHI)
!     right now this is set to high diffusion, longterm this should
!     revert to the old formulation 
      DIFFAC  =  1.
      AKS2    =  MIN( 1.0/(2.0*DT*PI**2)*DIFFAC, AKS2MX )
      AKS4    =  MIN( 1.0/(2.0*DT*PI**4)*DIFFAC, AKS4MX )
      DWP1    =  1.0 - DWP

!     VERTIKAL VARIIERENDER IMPLIZITHEITSGRAD DER V-DIFFUSION

      A1T(KE1)  = 1.2
      A1T(KE)   = 1.2
      A1T(KE-1) = 1.1
      A1T(KE-2) = 1.0
      A1T(KE-3) = 0.90
      A1T(KE-4) = 0.80

      DO K  = 1 , KE1
        A2T   (K) = 1.0 - A1T(K)
      ENDDO
C     BERECHNUNG VON RMY:
C     RMY(I,1) IST AM PS, H, QDW-GITTERPUNKT DEFINIERT
C     RMY(I,2) IST AM U         -GITTERPUNKT DEFINIERT
      ZDTDDTR = 1.0 !DT/300.0*0.5/MAX( DLAM, DPHI )
      RMY(:,:) = 0.0
      IAHG = 2
      IAUG = 2
      IEHG = IE - 1
      IEUG = IE - 2
      DO I = 1,IE
        ZDRH = MIN(I - IAHG + 0.25 , IEHG - I + 0.25)
        IF(ZDRH.LE.0.0) THEN
          RMY(I,1) = 1.0*ZDTDDTR
        ELSE
          RMY(I,1) = (1.0 - TANH(0.5*ZDRH))/TANH(0.5*ZDRH)*
     &         ZDTDDTR
        ENDIF
        ZDRU = MIN(I - IAUG + 0.75 , IEUG - I + 0.75)
        IF(ZDRU.LE.0.0) THEN
          RMY(I,2) = 1.0*ZDTDDTR
        ELSE
          RMY(I,2) = (1.0 - TANH(0.5*ZDRU))/TANH(0.5*ZDRU)*
     &         ZDTDDTR
        END IF
      END DO
      CALL GETALPHA

      END SUBROUTINE GKONST

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! please add new initializers only hereafter and other routines
        ! only before
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE INIT_STRAKA_BUBBLE
!     first define hydrostatically balanced background state with 
!     constant potential temperature of 300 K
      U(:,:,1) = 0.
      W(:,:,1) = 0.
      T(:,:,1) = 288.15
      PS(:,1) = 101325.
      FIB(:) = 0.
      CALL CALC_P_HL(1)
      PINT(:,:,1) = ZPNF(:,:)
      CALL CALC_ALOGP_HL
      CALL CALC_DP(1)
      CALL CALC_BETAK
      CALL CALC_ALOGP_FL(1)
      TPOT(:,:) = TPOT0
      TPOTINI(:,:) = TPOT(:,:)
      CALL CALC_T_FROM_TPOT(1)
      CALL CALC_FI_HL(1)
      CALL CALC_FI_FL

!     now add straka-like bubble
!     distance in m per degree longitude
      DISX = RERD*PI/180.*COS(PHILU) 

      DO K = 1, KE
        DO I = IAA, IEA
          BFAC = SQRT(((I-BCI)*DLAM*DISX/BRX)**2 + ((ZFIH(I,K)/G -
     &         BCZ)/BRZ)**2)
          IF (BFAC <= 1.) THEN
            T(I,K,1) = T(I,K,1) + BVAL*(COS(PI*BFAC) + 1.)/2.
          END IF
        END DO
      END DO
      END SUBROUTINE INIT_STRAKA_BUBBLE

      
      SUBROUTINE INIT_MOUNTAIN_WAVE
      
!     Initialize the mountain test case
!     from MPAS 2D

      REAL :: XC, AMPL, vgradT, DX
      REAL :: X(IE), HX(IE), HXPL(IE)
      INTEGER :: I

      ! first define hydrostatically balanced background state with 
      ! constant potential temperature of 300 K
      U(:,:,1) = 0.
      W(:,:,1) = 0.
      T(:,:,1) = 288.15
      PS(:,1) = 101325.
      FIB(:) = 0.
      TPOT(:,:) = TPOT0
      
!     adjust FIB
!      HM = 500.
!      XA = 2000.
      AMPL = 1.
      DX = RERD * CPHI(1) * DLAM * DEGRAD
      XC = 0.5 * IE * DX
      
      DO I=1,IE
        X   (I) = (I-1)*DX
        HX  (I) = HM/(1.+((X(I)-XC)/XA)**2)
        HXPL(I) = AMPL*HX(I)
      END DO
      FIB(:) = HXPL(:)*G
!     Periodic boundary conditions
      IF ( LPERIODIC ) THEN
        FIB(1    :IAH-1) = FIB(IEH-2:IEH-1)
        FIB(IEH+1:IE   ) = FIB(IAH+1:IAH+2)
      END IF
!     adjust U
      U(:,:,:) = UM
!     adjust PS
      vgradT = 0.0065           ! [K/m]
!      PS(:,1) = PS(:,1)*(1 - FIB(:) * vgradT / TPOT(:,KE))**(G/R/vgradT)
      PS(:,1) = PS(:,1)*EXP(-G*FIB(:)/(R*TPOT0))

      CALL CALC_P_HL(1)
      PINT(:,:,1) = ZPNF(:,:)
      CALL CALC_ALOGP_HL
      CALL CALC_DP(1)
      CALL CALC_BETAK
      CALL CALC_ALOGP_FL(1)
      CALL CALC_T_FROM_TPOT(1)
      CALL CALC_FI_HL(1)
      CALL CALC_FI_FL

!     Set stability for mountain wave test case. Recalculate temperature and
!     geopotential
      TPOT(:,:) = TPOT0 * EXP(BRUNT**2./G*ZFIH(:,:)/G)
      TPOTINI(:,:) = TPOT(:,:)
      CALL CALC_T_FROM_TPOT(1)
      CALL CALC_FI_HL(1)
      CALL CALC_FI_FL
      
      END SUBROUTINE INIT_MOUNTAIN_WAVE


      SUBROUTINE INIT_GRAVWAVE
!
!     INTERTIA-GRAVITY WAVES AFTER SKAMEROCK AND KLEMP
!
      REAL :: DTPOT(IE,KE)
      REAL :: ZTOP, XC, DX, X
      INTEGER :: I, K
!     DEFINE BACKGROUND STATE
      U(:,:,1) = UM
      W(:,:,1) = 0.
      T(:,:,1) = 288.15
      PS(:,1) = 101325.
      FIB(:) = 0.
      TPOT(:,:) = TPOT0
      TPOTINI(:,:) = TPOT0

      CALL CALC_P_HL(1)
      PINT(:,:,1) = ZPNF(:,:)
      CALL CALC_ALOGP_HL
      CALL CALC_DP(1)
      CALL CALC_BETAK
      CALL CALC_ALOGP_FL(1)
      CALL CALC_T_FROM_TPOT(1)
      CALL CALC_FI_HL(1)
      CALL CALC_FI_FL

      DO I = 1,3
        TPOT(:,:) = TPOT0 * EXP(BRUNTGW**2./G*ZFIH(:,:)/G)
        TPOTINI(:,:) = TPOT(:,:)
        CALL CALC_T_FROM_TPOT(1)
        CALL CALC_FI_HL(1)
        CALL CALC_FI_FL
      ENDDO
!     ADD POTENTIAL TEMPERATURE PERTURBATION
      DX = RERD * CPHI(1) * DLAM * DEGRAD
      XC = 0.5*IE * DX
      DO K = 1, KE
        DO I = IAA, IEA
          ZTOP = ZFIH(I,1)/G
          X = (I-1)*DX
          DTPOT(I,K) = DTPOT0 * SIN(PI*ZFIH(I,K)/G/ZTOP) /
     &         (1 + (X-XC)**2/PHW**2)
        END DO
      END DO

      TPOT(:,:) = TPOT(:,:) + DTPOT(:,:)
      CALL CALC_T_FROM_TPOT(1)
      CALL CALC_FI_HL(1)
      CALL CALC_FI_FL

      END SUBROUTINE INIT_GRAVWAVE

      END PROGRAM REMOVP
