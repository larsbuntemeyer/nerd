!
subroutine Hydro_solve(dt)
   !
   use Grid
   use Database
   !
   implicit none
   !
   double precision, intent(in) :: dt
   !
   ! the current fluxes:
   ! the flux is defined at the cell interfaces, so we
   ! for an x-grid with nx cells, we have nx+2*nguard+1
   ! cell interfaces 
   !
   double precision, dimension(nvar,nx+2*nguard+1,ny+k2d*(1+2*nguard),nz+k3d*(1+2*nguard)) :: flux
   double precision, dimension(nvar,nx+2*nguard,ny+k2d*2*nguard,nz+k3d*2*nguard) :: state
   !
   double precision :: rho,ei,ekin,vtot2
   double precision :: fx,fy,fz 
   !
   integer :: i,j,k,ivar
   !
   if(ndim>1) then
      write(*,*) '------------------------'
      write(*,*) 'ERROR in Hydro_solve:'
      write(*,*) 'only 1D possilbe!'
      write(*,*) '------------------------'
      stop
   endif
   !
   ! Fill guardcells 
   !
   do i=1,nguard
      do j=1,ny+2*k2d*nguard
         do k=1,nz+2*k3d*nguard
            !
            dens(i,j,k) = dens(nguard+1,j,k)
            eint(i,j,k) = eint(nguard+1,j,k)
            vx(i,j,k) = vx(nguard+1,j,k)
            vy(i,j,k) = vx(nguard+1,j,k)
            vz(i,j,k) = vx(nguard+1,j,k)
            !
         enddo
      enddo
   enddo 
   do i=1+nx+nguard,nx+2*nguard
      do j=1,ny+2*k2d*nguard
         do k=1,nz+2*k3d*nguard
            !
            dens(i,j,k) = dens(nx+nguard,j,k)
            eint(i,j,k) = eint(nx+nguard,j,k)
            vx(i,j,k) = vx(nx+nguard,j,k)
            vy(i,j,k) = vx(nx+nguard,j,k)
            vz(i,j,k) = vx(nx+nguard,j,k)
            !
         enddo
      enddo
   enddo 
   do j=1,1+k2d*nguard
      do i=1,nx+2*nguard
         do k=1,nz+2*k3d*nguard
            !
            dens(i,j,k) = dens(i,k2d*nguard+1,k)
            eint(i,j,k) = eint(i,k2d*nguard+1,k)
            vx(i,j,k) = vx(i,k2d*nguard+1,k)
            vy(i,j,k) = vy(i,k2d*nguard+1,k)
            vz(i,j,k) = vz(i,k2d*nguard+1,k)
            !
         enddo
      enddo
   enddo 
   do j=ny+k2d*(nguard+1),ny+2*k2d*nguard
      do i=1,nx+2*nguard
         do k=1,nz+2*k3d*nguard
            !
            dens(i,j,k) = dens(i,ny+k2d*nguard,k)
            eint(i,j,k) = eint(i,ny+k2d*nguard,k)
            vx(i,j,k) = vx(i,ny+k2d*nguard,k)
            vy(i,j,k) = vy(i,ny+k2d*nguard,k)
            vz(i,j,k) = vz(i,ny+k2d*nguard,k)
            !
         enddo
      enddo
   enddo 
   do k=1,1+k3d*nguard
      do i=1,nx+2*nguard
         do j=1,ny+2*k2d*nguard
            !
            dens(i,j,k) = dens(i,j,k3d*nguard+1)
            eint(i,j,k) = eint(i,j,k3d*nguard+1)
            vx(i,j,k) = vx(i,j,k3d*nguard+1)
            vy(i,j,k) = vy(i,j,k3d*nguard+1)
            vz(i,j,k) = vz(i,j,k3d*nguard+1)
            !
         enddo
      enddo
   enddo 
   do k=nz+k3d*(nguard+1),nz+2*k3d*nguard
      do i=1,nx+2*nguard
         do j=1,ny+2*k3d*nguard
            !
            dens(i,j,k) = dens(i,j,nz+k3d*nguard)
            eint(i,j,k) = eint(i,j,nz+k3d*nguard)
            vx(i,j,k) = vx(i,j,nz+k3d*nguard)
            vy(i,j,k) = vy(i,j,nz+k3d*nguard)
            vz(i,j,k) = vz(i,j,nz+k3d*nguard)
            !
         enddo
      enddo
   enddo 
   !
   ! Compute conserved variables
   !
   do i=1,nx+2*nguard
      do j=1,ny+2*k2d*nguard
         do k=1,nz+2*k3d*nguard
            !
            vtot2 = vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
            rho   = dens(i,j,k)
            ekin  = 0.5d0 * vtot2
            !
            state(1,i,j,k) = rho 
            state(2,i,j,k) = rho * vx(i,j,k)
            state(3,i,j,k) = rho * vy(i,j,k)
            state(4,i,j,k) = rho * vz(i,j,k)
            state(5,i,j,k) = rho * (eint(i,j,k) + ekin)
            !
         enddo
      enddo
   enddo 
   !
   ! Interpolate velocities at cell interfaces
   !
   do i=nguard+1,nx+nguard+1
      do j=nguard+k2d,ny+k2d*(nguard+1)
         do k=nguard+k3d,nz+k3d*(nguard+1)
            !
            vxf(i,j,k) = 0.5d0*(state(2,i,j,k)   / state(1,i,j,k) +  &
                                state(2,i-1,j,k) / state(1,i-1,j,k)) 
            vyf(i,j,k) = 0.5d0*(state(3,i,j,k)   / state(1,i,j,k) +  &
                                state(3,i,j-1,k) / state(1,i,j-1,k)) 
            vyf(i,j,k) = 0.5d0*(state(4,i,j,k)   / state(1,i,j,k) +  &
                                state(4,i,j,k-1) / state(1,i,j,k-1)) 
            !
         enddo
      enddo
   enddo 
   !
   ! construct the fluxes at cell interfaces       
   !
   do i=nguard+1,nx+nguard+1
      do j=nguard+k2d,ny+k2d*(nguard+1)
         do k=nguard+k3d,nz+k3d*(nguard+1)
            !
            if(vxf(i,j,k) < 0.d0) then
               flux(1,i,j,k) = vxf(i,j,k) *  state(1,i,j,k) 
               flux(2,i,j,k) = vxf(i,j,k) *  state(2,i,j,k) + pres(i,j,k)
               flux(3,i,j,k) = vxf(i,j,k) *  state(3,i,j,k)
               flux(4,i,j,k) = vxf(i,j,k) *  state(4,i,j,k)
               flux(5,i,j,k) = vxf(i,j,k) * (state(5,i,j,k) + pres(i,j,k))
            else 
               flux(1,i,j,k) = vxf(i,j,k) *  state(1,i-1,j,k) 
               flux(2,i,j,k) = vxf(i,j,k) *  state(2,i-1,j,k) + pres(i,j,k)
               flux(3,i,j,k) = vxf(i,j,k) *  state(3,i-1,j,k)
               flux(4,i,j,k) = vxf(i,j,k) *  state(4,i-1,j,k)
               flux(5,i,j,k) = vxf(i,j,k) * (state(5,i-1,j,k) + pres(i,j,k))
            endif
            !
         enddo
      enddo
   enddo 
   !
   ! Now advect
   !
   do i=nguard+1,nx+nguard
      do j=nguard+k2d,ny+k2d*nguard
         do k=nguard+k3d,nz+k3d*nguard
            do ivar=1,nvar
               state(ivar,i,j,k) = state(ivar,i,j,k) -         &
                                   dt/dx * (flux(ivar,i+1,j,k) - flux(ivar,i,j,k))
            enddo
         enddo
      enddo
   enddo 
   !
   ! Update simple variables
   !
   do i=nguard+1,nx+nguard
      do j=nguard+k2d,ny+k2d*nguard
         do k=nguard+k3d,nz+k3d*nguard
            do ivar=1,nvar
               rho = state(1,i,j,k)
               dens(i,j,k) = rho 
               vx(i,j,k)   = state(2,i,j,k) / rho
            !   vy(i,j,k)   = state(3,i,j,k) / rho
            !   vz(i,j,k)   = state(4,i,j,k) / rho
               vtot2       = vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2
               ekin        = 0.5d0 * vtot2
               ener(i,j,k) = state(5,i,j,k) / rho
               eint(i,j,k) = ener(i,j,k) - ekin
            enddo
         enddo
      enddo
   enddo 
   !
end subroutine Hydro_solve
