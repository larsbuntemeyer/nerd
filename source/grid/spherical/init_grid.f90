
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
   integer        :: i,j,k
   !
   !ndim    = 1
   nx      = nxb 
   ny      = nyb
   nz      = nzb
   nz      = nzb+1
   !k2d     = 0
   !k3d     = 0
   !nguard  = 2
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
   keg     = nz+2*k3d*nguard 
   !
   ! loop ranges excluding guardcells
   !
   ib      = nguard+1
   ie      = nx+nguard
   jb      = k2d*nguard+1
   je      = ny+k2d*nguard
   kb      = k3d*nguard+1
   ke      = nz+k3d*nguard
   !
   !
   !
   write(*,*) '----- Grid_init -------------------'
   write(*,*) 'grid parameters:'
   write(*,*) 'ndim:', ndim
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
   ALLOCATE(GCPHI(ny,2))
   ALLOCATE(GACPHIR(ny,2))
   ALLOCATE(ACPHIR(ny,2))
   ALLOCATE(CPHI(ny,2)) 
   ALLOCATE(A1T(nz1), A2T(nz1))
   ALLOCATE(AK(nz1), BK(nz1))
   ALLOCATE(AKH(nz), BKH(nz))
   ALLOCATE(DAK(nz), DBK(nz))
   ALLOCATE(VVFH(nz))
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
         do k=1,nz+2*k3d*nguard
            !
            ! The cell center-coordinates
            !
            xcCoord(i) = xmin - 1.0*nguard*dx + 0.5*dx + (i-1)*dx
            ycCoord(j) = ymin - 1.0*nguard*dy + 0.5*dy + (j-1)*dy
            zcCoord(k) = zmin - 1.0*nguard*dz + 0.5*dz + (k-1)*dz
            !
            ! The left cell-face coordinates
            !
            xlCoord(i) = xmin - 1.0*nguard*dx + (i-1)*dx
            ylCoord(j) = ymin - 1.0*nguard*dy + (j-1)*dy
            zlCoord(k) = zmin - 1.0*nguard*dz + (k-1)*dz
            !
            ! The right cell-face coordinates
            !
            xrCoord(i) = xmin - 1.0*nguard*dx  + i*dx
            yrCoord(j) = ymin - 1.0*nguard*dy + j*dy
            zrCoord(k) = zmin - 1.0*nguard*dz + k*dz
            !
         enddo
      enddo
   enddo
   !
   write(*,*) '----- Grid_init done --------------'
   !
end subroutine init_grid
