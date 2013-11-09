!
subroutine Grid_init
   !
   use Grid_data
   use RuntimeParameters
   ! 
   implicit none
   !
   integer        :: i,j,k
   !
   ndim    = 1
   nx      = 100
   ny      = 1
   nz      = 1
   k2d     = 0
   k3d     = 0
   nguard  = 1
   !
   write(*,*) '----- Grid_init -------------------'
   write(*,*) 'grid parameters:'
   write(*,*) 'ndim:', ndim
   write(*,*) 'nx,ny,nz:', nx,ny,nz
   write(*,*) 'nguard:', nguard
   !
   !  allocate coordinate fields
   !
   allocate(xcCoord(nx+2*nguard))
   allocate(ycCoord(ny+2*nguard))
   allocate(zcCoord(nz+2*nguard))
   allocate(xlCoord(nx+2*nguard))
   allocate(ylCoord(ny+2*nguard))
   allocate(zlCoord(nz+2*nguard))
   allocate(xrCoord(nx+2*nguard))
   allocate(yrCoord(ny+2*nguard))
   allocate(zrCoord(nz+2*nguard))
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
            xcCoord(i) = xmin - 1.d0*nguard*dx + 0.5d0*dx + (i-1)*dx
            ycCoord(j) = ymin - 1.d0*nguard*dy + 0.5d0*dy + (j-1)*dy
            zcCoord(k) = zmin - 1.d0*nguard*dz + 0.5d0*dz + (k-1)*dz
            !
            ! The left cell-face coordinates
            !
            xlCoord(i) = xmin - 1.d0*nguard*dx + (i-1)*dx
            ylCoord(j) = ymin - 1.d0*nguard*dy + (j-1)*dy
            zlCoord(k) = zmin - 1.d0*nguard*dz + (k-1)*dz
            !
            ! The right cell-face coordinates
            !
            xrCoord(i) = xmin - 1.d0*nguard*dx  + i*dx
            yrCoord(j) = ymin - 1.d0*nguard*dy + j*dy
            zrCoord(k) = zmin - 1.d0*nguard*dz + k*dz
            !
         enddo
      enddo
   enddo
   !
   write(*,*) '----- Grid_init done --------------'
   !
end subroutine Grid_init
!
