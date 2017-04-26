
subroutine init_grid
   ! 
   use mo_grid
   use mo_namelist
   use mo_parameters 
   use mo_constants
   !
   !
   implicit none
   !
   integer        :: i,j,k,kp1
   real :: zphis
   !
   !ndim    = 1
   nx      = nxb 
   ny      = nyb
   nz      = nzb
   nz1     = nzb+k3d
   !
   ! Precalculate loop indices as following:
   ! -----------------------------------------------------------------------------------------------
   !  ibg|              |    ib    |                 |     ie    |                   |     ieg     |
   ! -----------------------------------------------------------------------------------------------
   ! 
   ! | 1 | ... | nguard | nguard+1 | ... | ... | ... | nguard+nx | nguard+nx+1 | ... | 2*nguard+nx |
   ! -----------------------------------------------------------------------------------------------
   ! |   guard cells    |             valid cells                |          guard cells            |
   !                   xmin                                     xmax
   ! 
   !
   ! loop ranges including guardcells
   !
   ibg     = 1
   ieg     = nx+2*nguard 
   jbg     = 1
   jeg     = ny+2*k2d*nguard 
   kbg     = 1
   ! no guard cells in the layer dimension 
   keg     = nz
   keg1    = keg+k3d
   !
   ! loop ranges excluding guardcells
   !
   ib      = nguard+1
   ie      = nx+nguard
   jb      = k2d*nguard+1
   je      = ny+k2d*nguard
   ! no guard cells in the layer dimension 
   kb      = 1
   ke      = nz
   ke1     = ke+k3d
   !
   nxny    = ieg*jeg
   nxnz    = ieg*keg
   !
   !
   write(*,*) '----- Grid_init -------------------'
   write(*,*) 'grid parameters:'
   write(*,*) 'ndim:', ndim
   write(*,*) 'k2d, k3d:', k2d, k3d
   write(*,*) 'nx,ny,nz:', nx,ny,nz
   write(*,*) 'nguard:', nguard
   write(*,*) 'ib,jb,kb:', ib,jb,kb
   write(*,*) 'ie,je,ke:', ie,je,ke
   write(*,*) 'ibg,jbg,kbg:', ibg,jbg,kbg
   write(*,*) 'ieg,jeg,keg:', ieg,jeg,keg
   !
   !  allocate coordinate fields
   !
   allocate(xcCoord(ibg:ieg))
   allocate(ycCoord(jbg:jeg))
   allocate(zcCoord(kbg:keg))
   allocate(xlCoord(ibg:ieg))
   allocate(ylCoord(jbg:jeg))
   allocate(zlCoord(kbg:keg))
   allocate(xrCoord(ibg:ieg))
   allocate(yrCoord(jbg:jeg))
   allocate(zrCoord(kbg:keg))
   allocate(GCPHI(ny+2*k2d*nguard,2))
   ALLOCATE(GACPHIR(ny+2*k2d*nguard,2))
   ALLOCATE(ACPHIR(ny+2*k2d*nguard,2))
   ALLOCATE(CPHI(ny+2*k2d*nguard,2)) 
   ALLOCATE(A1T(nz1), A2T(nz1))
   ALLOCATE(AK(nz1), BK(nz1))
   ALLOCATE(AKH(nz), BKH(nz))
   ALLOCATE(DAK(nz), DBK(nz))
   ALLOCATE(VVFH(nz))
   ALLOCATE(LVL(nz))
   ALLOCATE(LVL2(nz1))
   !
   !  inititialize the grid
   !
   if(xmax-xmin<0.d0) then
      write(*,*) 'ERROR in Grind_init:' 
      write(*,*) 'xmax < xmin' 
   endif
   if(ndim>1.and.ymax-ymin<0.d0) then
      write(*,*) 'ERROR in Grind_init:' 
      write(*,*) 'ymax < ymin' 
   endif
   if(ndim==3.and.zmax-zmin<0.d0) then
      write(*,*) 'ERROR in Grind_init:' 
      write(*,*) 'zmax < zmin' 
   endif
   !
   dx = (xmax-xmin)/nx
   dy = (ymax-ymin)/ny
   dz = (zmax-zmin)/nz
   !
   dlam = dx
   dphi = dy 
   EDDLAM  =   1.0/(DLAM*PI/180.0)
   EDDPHI  =   1.0/(DPHI*PI/180.0)
   EDADPHI =   EDDPHI/R_EARTH
   DLADDPH =   DLAM/DPHI
   DPHDDLA =   DPHI/DLAM
   !
   write(*,*) 'dx:', dx
   write(*,*) 'dy:', dy
   write(*,*) 'dz:', dz
   !
   write(*,*) 'initializing grid-coordinates'
   !
   do i=1,nx+2*nguard
      do j=1,ny+2*k2d*nguard
            !
            ! The cell center-coordinates
            !
            xcCoord(i) = xmin - 1.0*nguard*dx + 0.5*dx + (i-1)*dx
            ycCoord(j) = ymin - 1.0*nguard*dy + 0.5*dy + (j-1)*dy
            !
            ! The left cell-face coordinates
            !
            xlCoord(i) = xmin - 1.0*nguard*dx + (i-1)*dx
            ylCoord(j) = ymin - 1.0*nguard*dy + (j-1)*dy
            !
            ! The right cell-face coordinates
            !
            xrCoord(i) = xmin - 1.0*nguard*dx  + i*dx
            yrCoord(j) = ymin - 1.0*nguard*dy + j*dy
            !
      enddo
   enddo
   do k=1,ke
     zcCoord(k) = k
   enddo
!     BERECHNUNG VON PHI, RLA, FC, CPHI, ACPHIR UND RMY
!     JG IST DER GLOBALE J-INDEX
   DO J=1,JEG
      ZPHIS       = ycCoord(j)!PHILU + (JG - 1)*DPHI
      CPHI  (J,1) = COS(ZPHIS*DEGRAD)
      CPHI  (J,2) = COS((ZPHIS + 0.5*DPHI)*DEGRAD)
      ACPHIR(J,1) = 1.0/(RERD*CPHI(J,1))
      ACPHIR(J,2) = 1.0/(RERD*CPHI(J,2))
   ENDDO

!  BERECHNUNG VON GCPHI UND GACPHIR IM GESAMTGEBIET
   DO J=1,JEG
      ZPHIS        = ycCoord(j)!PHILU + (J - 1)*DPHI
      GCPHI  (J,1) = COS(ZPHIS*DEGRAD)
      GCPHI  (J,2) = COS((ZPHIS + 0.5*DPHI)*DEGRAD)
      GACPHIR(J,1) = 1.0/(RERD*GCPHI(J,1))
      GACPHIR(J,2) = 1.0/(RERD*GCPHI(J,2))
   ENDDO
   !
   call read_akbk
   !
   DO K  = 1 , KE
     KP1 = K + 1
     AKH   (K) = 0.5*(AK(K) + AK(KP1))
     BKH   (K) = 0.5*(BK(K) + BK(KP1))
     DAK   (K) =    - AK(K) + AK(KP1)
     DBK   (K) =    - BK(K) + BK(KP1)
   ENDDO
   !
   PTOP = AK(1)
   !KS   Check if the model top is different from 0 Pa.
   IF (PTOP < 1.0E-10) THEN
     LPTOP0 = .TRUE.
   ELSE
     LPTOP0 = .FALSE.
   END IF
   !
   write(*,*) '----- Grid_init done --------------'
   !
end subroutine init_grid
