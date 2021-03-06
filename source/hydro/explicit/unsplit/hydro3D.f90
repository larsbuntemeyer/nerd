subroutine hydro3D(dt)
      !
      use mo_grid, only: ib,ie,kbg,keg,jbg,jeg,ibg,ieg,jb,je,kb,ke, &
                      xcCoord,ycCoord,zcCoord,       &
                      xlCoord,xrCoord,zlCoord,zrCoord, &
                      ylCoord,yrCoord
      use mo_parameters, only: nvar,ndim,k2d,k3d
      use mo_database, only: dens,pres,u,v,w,eint,uf,vf,wf,ener
      use mo_hydro, only: interface_flux
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
      real :: fx,fy,fz,accl,dts,idt
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
               uf(i,j,k) = 0.5 * (u(i,j,k) + u(i-1,j,k))
               if(ndim>1)  vf(i,j,k) = 0.5 * (v(i,j,k) + v(i,j-1,k))
               if(ndim==3) wf(i,j,k) = 0.5 * (w(i,j,k) + w(i,j,k-1))
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
                 dx = xcCoord(i) - xcCoord(i-1)
                 xflux(ivar,i,j,k) = interface_flux(q(ivar,i-2,j,k),q(ivar,i-1,j,k), &
                                                q(ivar,i,j,k),  q(ivar,i+1,j,k), &
                                                uf(i,j,k),dt,dx)
                 if(ndim>1) then
                 dy = ycCoord(j) - ycCoord(j-1)
                 yflux(ivar,i,j,k) = interface_flux(q(ivar,i,j-2,k),q(ivar,i,j-1,k), &
                                                q(ivar,i,j,k),  q(ivar,i,j+1,k), &
                                                vf(i,j,k),dt,dy)
                 endif
                 if(ndim==3) then
                 dz = zcCoord(k) - zcCoord(k-1)
                 zflux(ivar,i,j,k) = interface_flux(q(ivar,i,j,k-2),q(ivar,i,j,k-1), &
                                                q(ivar,i,j,k),  q(ivar,i,j,k+1), &
                                                wf(i,j,k),dt,dz)
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
                  if(ndim>1) then 
                    q(ivar,i,j,k) = q(ivar,i,j,k) +         &
                         dt/(yrCoord(j)-ylCoord(j)) * (yflux(ivar,i,j,k) - yflux(ivar,i,j+1,k))
                  endif
                  if(ndim==3) then
                  q(ivar,i,j,k) = q(ivar,i,j,k) +         &
                       dt/(zrCoord(k)-zlCoord(k)) * (zflux(ivar,i,j,k) - zflux(ivar,i,j,k+1))
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
      call eos
      !
      ! Now add pressure force
      !
      do k=kb,ke
        do j=jb,je
          do i=ib,ie
               rhoinv = 1.0/dens(i,j,k)
               dp = pres(i+1,j,k) - pres(i-1,j,k)
               dx = xcCoord(i+1)  - xcCoord(i-1)
               !dp = 0.5*(pres(i+1,j,k)+pres(i,j,k)) - 0.5*(pres(i-1,j,k)+pres(i,j,k))
               !dx = xrCoord(i)  - xlCoord(i)
               u(i,j,k)   = u(i,j,k) - dt * dp/dx * rhoinv
               if(ndim>1) then 
                 dp = pres(i,j+1,k) - pres(i,j-1,k) 
                 dy = ycCoord(j+1)  - ycCoord(j-1)
                 !dp = 0.5*(pres(i,j,k)+pres(i,j+1,k)) - 0.5*(pres(i,j,k)+pres(i,j-1,k))
                 !dy = yrCoord(j)  - ylCoord(j)
                 accl = -dt*dp/dy*rhoinv -dt*9.81
                 !v(i,j,k)   = v(i,j,k) - dt * dp/dy * rhoinv
                 !v(i,j,k)   = v(i,j,k) - dt * 9.81
                 v(i,j,k)   = v(i,j,k) + accl
                 !if(abs(accl).gt.1.e-4) then
                 !  print*, 'accl > 0', accl, xcCoord(i),ycCoord(j),i,j,pres(i,j+1,k),pres(i,j,k),pres(i,j-1,k)
                 !  stop
                 !endif 
               endif 
               if(ndim==3) then 
                 dp = pres(i,j,k+1) - pres(i,j,k-1)
                 dz = zcCoord(k+1)  - zcCoord(k-1)
                 w(i,j,k)   = w(i,j,k) - dt * dp/dz * rhoinv
               endif 
               ! 
               !vtot2       = u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2
               !ekin        = 0.5d0 * vtot2
               !ener(i,j,k) = q(5,i,j,k) / rho
               !eint(i,j,k) = ener(i,j,k) - ekin
            enddo
         enddo
      enddo
      !
end subroutine hydro3D
