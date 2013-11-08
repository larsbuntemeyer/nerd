!
subroutine Database_init
   !
   use Database
   use Grid
   !
   implicit none
   !
   write(*,*) '----- Database_init ---------------'
   !
   nvar = 5 
   !
   write(*,*) 'allocating'
   allocate(dens(nx+2*nguard,ny+2*k2d*nguard,nz+2*k3d*nguard))
   allocate(eint(nx+2*nguard,ny+2*k2d*nguard,nz+2*k3d*nguard))
   !
   allocate(vx(nx+2*nguard,ny+k2d*2*nguard,nz+k3d*2*nguard))
   allocate(vy(nx+2*nguard,ny+k2d*2*nguard,nz+k3d*2*nguard))
   allocate(vz(nx+2*nguard,ny+k2d*2*nguard,nz+k3d*2*nguard))
   !
   allocate(vxf(nx+1+2*nguard,ny+k2d*(1+2*nguard),nz+k3d*(1+2*nguard)))
   allocate(vyf(nx+1+2*nguard,ny+k2d*(1+2*nguard),nz+k3d*(1+2*nguard)))
   allocate(vzf(nx+1+2*nguard,ny+k2d*(1+2*nguard),nz+k3d*(1+2*nguard)))
   !
   allocate(ener(nx+2*nguard,ny+2*k2d*nguard,nz+2*k3d*nguard))
   allocate(pres(nx+2*nguard,ny+2*k2d*nguard,nz+2*k3d*nguard))
   allocate(temp(nx+2*nguard,ny+2*k2d*nguard,nz+2*k3d*nguard))
   write(*,*) '----- Database_init done ----------'
   !
end subroutine Database_init
!
