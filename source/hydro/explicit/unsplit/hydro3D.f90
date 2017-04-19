subroutine hydro3D(dt)
      !
      use Grid, only: ib,ie,kbg,keg,jbg,jeg,ibg,ieg,jb,je,kb,ke, &
                      xcCoord,ycCoord,zcCoord,       &
                      xlCoord,xrCoord,zlCoord,zrCoord, &
                      ylCoord,yrCoord,k2d,k3d
      use Database, only: dens,pres,u,v,w,eint,nvar,uf,vf,wf,ener
      use Eos, only: Eos_gamma
      !
      implicit none
      !
      real, intent(in) :: dt
      !
      ! the current fluxes:
      ! the flux is defined at the cell interfaces, so we
      ! for an x-grid with nx cells, we have nx+2*nguard+1
      ! cell interfaces 
      !
      real, dimension(nvar,ibg:ieg+1,jbg:jeg+1,kbg:keg+1) :: xflux,yflux,zflux
      real, dimension(nvar,ibg:ieg  ,jbg:jeg  ,kbg:keg  ) :: q
      real, dimension(nvar,ibg:ieg  ,jbg:jeg  ,kbg:keg  ) :: dq
      !
      real :: rho,ei,ekin,vtot2,dp,dx,dy,dz,rhoinv
      real :: fx,fy,fz
      !
      integer :: i,j,k,ivar
      !
      ! Compute conserved variables
      !
      do k=kbg,keg
         do j=jbg,jeg
            do i=ibg,ieg
               !
               vtot2 = u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2
               rho   = dens(i,j,k)
               ekin  = 0.5d0 * vtot2
               !
               q(1,i,j,k) = rho 
               q(2,i,j,k) = rho * u(i,j,k)
               q(3,i,j,k) = rho * v(i,j,k)
               q(4,i,j,k) = rho * w(i,j,k)
               q(5,i,j,k) = rho * (eint(i,j,k) + ekin)
               !
            enddo
         enddo
      enddo 
      !
      ! Interpolate velocities at cell interfaces
      !
      do k=kb,ke+k3d
         do j=jb,je+k2d
            do i=ib,ie+1
               !
               uf(i,j,k) = 0.5d0 * (u(i,j,k) + u(i-1,j,k))
               if(k2d==1) vf(i,j,k) = 0.5d0 * (v(i,j,k) + v(i,j-1,k))
               if(k3d==1) wf(i,j,k) = 0.5d0 * (w(i,j,k) + w(i,j,k-1))
               !
            enddo
         enddo
      enddo 
      !
      ! construct the fluxes at cell interfaces       
      !
      do k=kb,ke+k3d
         do j=jb,je+k2d
            do i=ib,ie+1
               !
               do ivar=1,nvar
                 xflux(ivar,i,j,k) = interface_flux(q(ivar,i-2,j,k),q(ivar,i-1,j,k), &
                                                q(ivar,i,j,k),  q(ivar,i+1,j,k), &
                                                uf(i,j,k),dt)
                 if(k2d==1) then
                 yflux(ivar,i,j,k) = interface_flux(q(ivar,i,j-2,k),q(ivar,i,j-1,k), &
                                                q(ivar,i,j,k),  q(ivar,i,j+1,k), &
                                                vf(i,j,k),dt)
                 endif
                 if(k3d==1) then
                 zflux(ivar,i,j,k) = interface_flux(q(ivar,i,j,k-2),q(ivar,i,j,k-1), &
                                                q(ivar,i,j,k),  q(ivar,i,j,k+1), &
                                                wf(i,j,k),dt)
                 endif
               enddo
               !
            enddo
         enddo
      enddo 
      !
      ! Now advect
      !
      do i=ib,ie
         do j=jb,je
            do k=kb,ke
               do ivar=1,nvar
                  q(ivar,i,j,k) = q(ivar,i,j,k) +         &
                       dt/(xrCoord(i)-xlCoord(i)) * (xflux(ivar,i,j,k) - xflux(ivar,i+1,j,k))
                  if(k2d==1) then
                  q(ivar,i,j,k) = q(ivar,i,j,k) +         &
                       dt/(yrCoord(i)-ylCoord(i)) * (yflux(ivar,i,j,k) - yflux(ivar,i,j+1,k))
                  endif
                  if(k3d==1) then
                  q(ivar,i,j,k) = q(ivar,i,j,k) +         &
                       dt/(zrCoord(i)-zlCoord(i)) * (zflux(ivar,i,j,k) - zflux(ivar,i,j,k+1))
                  endif
               enddo
            enddo
         enddo
      enddo 
      !
      ! Update simple variables
      !
      do i=ib,ie
         do j=jb,je
            do k=kb,ke
               rho = q(1,i,j,k)
               dens(i,j,k) = rho 
               u(i,j,k)   = q(2,i,j,k) / rho
               v(i,j,k)   = q(3,i,j,k) / rho
               w(i,j,k)   = q(4,i,j,k) / rho
               vtot2       = u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2
               ekin        = 0.5d0 * vtot2
               ener(i,j,k) = q(5,i,j,k) / rho
               eint(i,j,k) = ener(i,j,k) - ekin
            enddo
         enddo
      enddo 
      !
      call Eos_gamma 
      !
      ! Now add pressure force
      !
      do i=ib,ie
         do j=jb,je
            do k=kb,ke
               rhoinv = 1.0/dens(i,j,k)
               dp = pres(i+1,j,k) - pres(i-1,j,k)
               dx = xcCoord(i+1)  - xcCoord(i-1)
               u(i,j,k)   = u(i,j,k) - dt * dp/dx * rhoinv
               !if(k2d==1) then 
               !  dp = pres(i,j+1,k) - pres(i,j-1,k)
               !  dy = ycCoord(i-1)  - ycCoord(i+1)
               !  v(i,j,k)   = v(i,j,k) - dt/dy * dp * rhoinv
               !endif 
               !if(k3d==1) then 
               !  dp = pres(i,j+1,k) - pres(i,j-1,k)
               !  dz = zcCoord(i-1)  - zcCoord(i+1)
               !  w(i,j,k)   = w(i,j,k) - dt/dz * dp * rhoinv
               !endif 
               ! 
               !vtot2       = u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2
               !ekin        = 0.5d0 * vtot2
               !ener(i,j,k) = q(5,i,j,k) / rho
               !eint(i,j,k) = ener(i,j,k) - ekin
            enddo
         enddo
      enddo
      call Eos_gamma 
      !
end subroutine hydro3D
