!
!
!
subroutine sweep1D(dt,dx,dir,rho,u,v,w,eint,pres,n,nvar)
   !
   use mo_namelist, only: gamma,fl
   use mo_grid, only: ibg,ieg,nx,ib,ie
   use mo_hydro
   !
   implicit none
   !
   real,    intent(in) :: dx
   real,    intent(in) :: dt
   integer, intent(in) :: dir,n,nvar
   real, dimension(n), intent(in)    :: pres 
   real, dimension(n), intent(inout) :: rho,u,v,w,eint
   !
   ! the current fluxes:
   ! the flux is defined at the cell interfaces, so we
   ! for an x-grid with nx cells, we have nx+2*nguard+1
   ! cell interfaces...
   !
   real, dimension(nvar,n)   :: q
   real, dimension(n)        :: enth
   ! Roe averages at cell interfaces
   real, dimension(n+1)      :: velx_a,vely_a,velz_a
   real, dimension(n+1)      :: enth_a,vel_a,cs_a
   ! average eigenvalues at cell interfaces
   real, dimension(nvar,n+1) :: dq,t,p,r,l,a,K1,K2,K3,K4,K5
   real, dimension(nvar,n+1) :: fl_left,fl_right,fl_diss,fl_roe
   !
   real :: rho2,rho2m1,ei,ekin,vtot2,dq5a,ener,dtdx,pres_a,rho_a
   real :: fx,fy,fz
   !
   integer :: i,j,k,ivar
   !
   ! State Vector
   ! contains conserved variables and enthalpy (including guard cells)
   !
   do i=ibg,ieg
      !
      vtot2 = u(i)**2+v(i)**2+w(i)**2
      ekin  = 0.5 * vtot2
      !
      q(1,i)  = rho(i)
      q(2,i)  = rho(i) * u(i)
      q(3,i)  = rho(i) * v(i)
      q(4,i)  = rho(i) * w(i)
      q(5,i)  = rho(i) * (eint(i) + ekin)
      enth(i) = eint(i) + ekin + pres(i)/rho(i) 
      !
   enddo
   !
   ! Compute jumps in conserved quantities at cell interfaces
   !
   do ivar=1,nvar
     do i=ib,ie+1
       !
       dq(ivar,i) = q(ivar,i) - q(ivar,i-1)
       !
     enddo
     !dq(ivar,1)   = 0.0!dq(ivar,2)
     !dq(ivar,n+1) = 0.0!dq(ivar,n)
   enddo
   !
   ! Compute Roe averages at cell interfaces
   !
   do i=ib,ie+1
     !
     rho2m1 = sqrt(rho(i-1))
     rho2   = sqrt(rho(i))
     pres_a = (pres(i-1)*rho2m1 + pres(i)*rho2) / (rho2m1+rho2)
     rho_a = (rho(i-1)*rho2m1 + rho(i)*rho2) / (rho2m1+rho2)
     velx_a(i) = (u(i-1)*rho2m1 + u(i)*rho2) / (rho2m1+rho2)
     vely_a(i) = (v(i-1)*rho2m1 + v(i)*rho2) / (rho2m1+rho2)
     velz_a(i) = (w(i-1)*rho2m1 + w(i)*rho2) / (rho2m1+rho2)
     enth_a(i) = (enth(i-1)*rho2m1 + enth(i)*rho2) / (rho2m1+rho2)
     vel_a(i) = sqrt(velx_a(i)**2+vely_a(i)**2+velz_a(i)**2)
     !cs_a(i) = sqrt((gamma-1.0)*(enth_a(i)-0.5*vel_a(i)**2))
     !cs_a(i) = (gamma-1.0)*sqrt(enth_a(i)-0.5*vel_a(i)**2)
     cs_a(i) = sqrt(gamma*pres_a/rho_a) 
     !
   enddo
   !
   do i=ib,ie+1
     !
     ! Eigenvalues 
     !
     l(1,i) = velx_a(i) - cs_a(i)
     l(2,i) = velx_a(i)
     l(3,i) = velx_a(i)
     l(4,i) = velx_a(i)
     l(5,i) = velx_a(i) + cs_a(i)
     !
     ! K-Vectors (Eigenvectors)
     !
     K1(1,i) = 1.0
     K1(2,i) = velx_a(i) - cs_a(i)
     K1(3,i) = vely_a(i)
     K1(4,i) = velz_a(i)
     K1(5,i) = enth_a(i) - velx_a(i)*cs_a(i)
     !
     K2(1,i) = 1.0
     K2(2,i) = velx_a(i)
     K2(3,i) = vely_a(i)
     K2(4,i) = velz_a(i)
     K2(5,i) = 0.5*vel_a(i)**2
     !
     K3(1,i) = 0.0
     K3(2,i) = 0.0 
     K3(3,i) = 1.0 
     K3(4,i) = 0.0 
     K3(5,i) = vely_a(i)
     !
     K4(1,i) = 0.0
     K4(2,i) = 0.0 
     K4(3,i) = 0.0 
     K4(4,i) = 1.0 
     K4(5,i) = velz_a(i)
     !
     K5(1,i) = 1.0
     K5(2,i) = velx_a(i) + cs_a(i)
     K5(3,i) = vely_a(i)
     K5(4,i) = velz_a(i)
     K5(5,i) = enth_a(i) + velx_a(i)*cs_a(i)
     !
     ! wavestrengths
     !
     a(3,i) = dq(3,i) - vely_a(i)*dq(1,i)
     a(4,i) = dq(4,i) - velz_a(i)*dq(1,i)
     dq5a   = dq(5,i) - (dq(3,i)-vely_a(i)*dq(1,i))*vely_a(i) - (dq(4,i)-velz_a(i)*dq(1,i))*velz_a(i)
     a(2,i) = (gamma-1.0)/(cs_a(i)**2) * (dq(1,i)*(enth_a(i)-velx_a(i)**2)+velx_a(i)*dq(2,i)-dq5a)
     a(1,i) = 1.0/(2.0*cs_a(i)) * (dq(1,i)*(velx_a(i)+cs_a(i))-dq(2,i)-cs_a(i)*a(2,i))
     a(5,i) = dq(1,i)-(a(1,i)+a(2,i))
     !
   enddo
   !
   ! timestep
   !
   !dt = cfl_timestep(nx+1,nvar,l(:,ib:ie+1),dx)
   !
   ! compute left and right flux vectors at cell interfaces
   !
   do i=ib,ie+1
     fl_left(1,i)  =  q(1,i-1)*u(i-1)
     fl_left(2,i)  = (q(2,i-1)*u(i-1)+pres(i-1))
     fl_left(3,i)  =  q(3,i-1)*u(i-1)
     fl_left(4,i)  =  q(4,i-1)*u(i-1)
     fl_left(5,i)  = (q(5,i-1)+pres(i-1))*u(i-1)
     fl_right(1,i) =  q(1,i)*u(i)
     fl_right(2,i) = (q(2,i)*u(i)+pres(i))
     fl_right(3,i) =  q(3,i)*u(i)
     fl_right(4,i) = q(4,i)*u(i)
     fl_right(5,i) = (q(5,i)+pres(i))*u(i)
   enddo
   !
   ! compute slopes for the flux limiter
   !
   call compute_slope(a,l,n,nvar,r)
   !
   ! compute flip flop function and flux limiter
   !
   do ivar=1,nvar
     do i=ib,ie+1
       ! flip flop
       t(ivar,i) = sign(1.0,l(ivar,i)) 
       ! flux limiter
       p(ivar,i) = phi(fl,r(ivar,i))
     enddo
   enddo
   !
   ! compute the dissipative flux
   !
   dtdx = dt/dx
   do ivar=1,nvar
     do i=ib,ie+1
       fl_diss(ivar,i) = a(1,i)*l(1,i)*K1(ivar,i) * (t(1,i)+p(1,i)*(l(1,i)*dtdx-t(1,i))) &
                       + a(2,i)*l(2,i)*K2(ivar,i) * (t(2,i)+p(2,i)*(l(2,i)*dtdx-t(2,i))) &
                       + a(3,i)*l(3,i)*K3(ivar,i) * (t(3,i)+p(3,i)*(l(3,i)*dtdx-t(3,i))) &
                       + a(4,i)*l(4,i)*K4(ivar,i) * (t(4,i)+p(4,i)*(l(4,i)*dtdx-t(4,i))) &
                       + a(5,i)*l(5,i)*K5(ivar,i) * (t(5,i)+p(5,i)*(l(5,i)*dtdx-t(5,i)))
     enddo
   enddo
   !
   ! compute Roe flux at cell interfaces 
   !
   do ivar=1,nvar
     do i=ib,ie+1
       fl_roe(ivar,i) = 0.5*(fl_left(ivar,i)+fl_right(ivar,i)) - 0.5*fl_diss(ivar,i)
     enddo
   enddo
   !
   ! advect in flux conserving form
   !
   do ivar=1,nvar
     do i=ib,ie+1
       q(ivar,i) = q(ivar,i) + dtdx*(fl_roe(ivar,i)-fl_roe(ivar,i+1))
     enddo
   enddo 
   !
   ! Update simple variables
   !
   do i=ib,ie
      rho(i) = q(1,i)
      u(i)   = q(2,i) / q(1,i)
      v(i)   = q(3,i) / q(1,i)
      w(i)   = q(4,i) / q(1,i)
      ener   = q(5,i) / q(1,i)
      vtot2  = u(i)**2+v(i)**2+w(i)**2
      ekin   = 0.5 * vtot2
      eint(i) = ener - ekin
   enddo
   !
   contains
     !
     ! flip flop function
     !
     real function ff(x)
     implicit none
     real, intent(in) :: x
     ff = sign(1.0,x)
     end function ff
     !
end subroutine sweep1D
