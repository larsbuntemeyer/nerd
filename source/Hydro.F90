!
! Hydro module
!
module Hydro
   !
   use Grid
   use Database
   !
   implicit none
   !
   public  :: Hydro_init, Hydro_solve
   !
   private :: hydro1D, sweep1D, fill_guardcells, minmod, &
              interface_flux
   !
contains
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine Hydro_init
      !
      implicit none
      !
      write(*,*) '----- Hydro_init ------------------'
      write(*,*) 'Hydro is not doing anything yet,   '
      write(*,*) 'because i have to admit, i was lazy'
      write(*,*) '----- Hydro_init done -------------'
      !
   end subroutine Hydro_init
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine sweep1D(dt,dir)
      !
      implicit none
      !
      real,    intent(in) :: dt
      integer, intent(in) :: dir
      !
      ! the current fluxes:
      ! the flux is defined at the cell interfaces, so we
      ! for an x-grid with nx cells, we have nx+2*nguard+1
      ! cell interfaces
      !
      real, dimension(nvar,ibg:ieg+1,jbg:jeg+1,kbg:keg+1) :: flux,dq
      real, dimension(nvar,ibg:ieg  ,jbg:jeg  ,kbg:keg  ) :: q
      !
      real :: rho,ei,ekin,vtot2
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
      ! Compute jumps in conserved quantities
      !
      do ivar=1,nvar
        do k=kbg,keg
          do j=jbg,jeg
            do i=ibg+1,ieg-1
              !
              dq(ivar,i,j,k) = q(ivar,i,j,k) - q(ivar,i-1,j,k)
              !
            enddo
          enddo
        enddo
      enddo
      !
   end subroutine sweep1D
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine Hydro_solve(dt)
      !
      implicit none
      !
      real, intent(in) :: dt
      !
      if(ndim>1) then
         write(*,*) '------------------------'
         write(*,*) 'ERROR in Hydro_solve:'
         write(*,*) 'only 1D possilbe!'
         write(*,*) '------------------------'
         stop
      endif
      !
      call fill_guardcells
      call hydro1D(dt)
      !
   end subroutine Hydro_solve
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine hydro1D(dt)
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
      real, dimension(nvar,ibg:ieg+1,jbg:jeg+1,kbg:keg+1) :: flux
      real, dimension(nvar,ibg:ieg  ,jbg:jeg  ,kbg:keg  ) :: q,dq
      !
      real :: rho,ei,ekin,vtot2
      real :: fx,fy,fz
      !
      integer :: i,j,k,ivar
      !
      ! Compute conserved variables
      !
      do i=ibg,ieg
         do j=jbg,jeg
            do k=kbg,keg
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
      ! Compute jumps conserved variables
      !
      do i=ibg,ieg
         do j=jbg,jeg
            do k=kbg,keg
               do ivar=1,nvar
                  !
                  dq(ivar,i,j,k) = q(ivar,i,j,k)-q(ivar,i-1,j,k)
                  !
               enddo
            enddo
         enddo
      enddo 
      !
      ! Interpolate velocities at cell interfaces
      !
      do i=ib,ie+1
         do j=jb,je+k2d
            do k=kb,ke+k3d
               !
               !uf(i,j,k) = 0.5d0*(q(2,i,j,k)   / q(1,i,j,k) +  &
               !                    q(2,i-1,j,k) / q(1,i-1,j,k)) 
               !vf(i,j,k) = 0.5d0*(q(3,i,j,k)   / q(1,i,j,k) +  &
               !                    q(3,i,j-1,k) / q(1,i,j-1,k)) 
               !wf(i,j,k) = 0.5d0*(q(4,i,j,k)   / q(1,i,j,k) +  &
               !                    q(4,i,j,k-1) / q(1,i,j,k-1)) 
               uf(i,j,k) = 0.5d0 * (u(i,j,k) + u(i-1,j,k))
               !
            enddo
         enddo
      enddo 
      !
      ! construct the fluxes at cell interfaces       
      !
      do i=ib,ie+1
         do j=jb,je+k2d
            do k=kb,ke+k3d
               !
               flux(1,i,j,k) = interface_flux(q(1,i-2,j,k),q(1,i-1,j,k), &
                                              q(1,i,j,k),  q(1,i+1,j,k), &
                                              uf(i,j,k),dt)
             !  if(uf(i,j,k) < 0.d0) then
             !     flux(1,i,j,k) = uf(i,j,k) *  q(1,i,j,k) 
             !     flux(2,i,j,k) = uf(i,j,k) *  q(2,i,j,k) + pres(i,j,k)
             !     flux(3,i,j,k) = uf(i,j,k) *  q(3,i,j,k)
             !     flux(4,i,j,k) = uf(i,j,k) *  q(4,i,j,k)
             !     flux(5,i,j,k) = uf(i,j,k) * (q(5,i,j,k) + pres(i,j,k))
             !  else 
             !     flux(1,i,j,k) = uf(i,j,k) *  q(1,i-1,j,k) 
             !     flux(2,i,j,k) = uf(i,j,k) *  q(2,i-1,j,k) + pres(i,j,k)
             !     flux(3,i,j,k) = uf(i,j,k) *  q(3,i-1,j,k)
             !     flux(4,i,j,k) = uf(i,j,k) *  q(4,i-1,j,k)
             !     flux(5,i,j,k) = uf(i,j,k) * (q(5,i-1,j,k) + pres(i,j,k))
             !  endif
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
                       dt/(xrCoord(i)-xlCoord(i)) * (flux(ivar,i,j,k) - flux(ivar,i+1,j,k))
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
               !u(i,j,k)   = q(2,i,j,k) / rho
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
   end subroutine hydro1D
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine fill_guardcells
      !
      implicit none
      !
      integer :: i,j,k,ivar
      !
      ! Fill guardcells 
      !
      do i=ibg,nguard
         do j=jbg,jeg
            do k=jbg,keg
               !
               dens(i,j,k) = dens(ib,j,k)
               eint(i,j,k) = eint(ib,j,k)
               u(i,j,k) = u(ib,j,k)
               v(i,j,k) = u(ib,j,k)
               w(i,j,k) = u(ib,j,k)
               !
            enddo
         enddo
      enddo 
      do i=ie+1,ieg
         do j=jbg,jeg
            do k=kbg,keg
               !
               dens(i,j,k) = dens(ie,j,k)
               eint(i,j,k) = eint(ie,j,k)
               u(i,j,k) = u(ie,j,k)
               v(i,j,k) = u(ie,j,k)
               w(i,j,k) = u(ie,j,k)
               !
            enddo
         enddo
      enddo 
      !
   end subroutine fill_guardcells
   !
   !----------------------------------------------------------------------------------------------
   !
   real function interface_flux(q1,q2,q3,q4,v_face,dt)
      !
      implicit none
      !
      real :: q1,q2,q3,q4
      real :: v_face,dt
      real :: r,phi,flux,theta
      !
      theta = sign(1.e0,v_face)
      !
      if(abs(q3-q2).gt.0.e0) then
         if(v_face.ge.0e0) then
            r = (q2-q1)/(q3-q2)
         else 
            r = (q4-q3)/(q3-q2)
         endif
      else
         r = 0.d0
      endif
      !
      select case(fl)
         !
         case('donor-cell')
            phi = 0.d0
         case('Lax-Wendroff')
            phi = 1.d0
         case('Beam-Warming')
            phi = r
         case('Fromm')
            phi = 0.5e0*(1.e0+r)
         case('minmod')
            phi = minmod(1.e0,r)
         case('superbee')
            phi = max(0.e0,min(1.e0,2.e0*r),min(2.e0,r))
         case default
            phi = 0.e0
      end select
      !
      flux = 0.5d0*v_face*((1.e0+theta)*q2+(1.e0-theta)*q3) +  &
             0.5d0*abs(v_face)*(1.e0-abs(v_face*dt/dx))*phi*(q3-q2)
             !
      interface_flux = flux
      !
   end function interface_flux
   !
   !----------------------------------------------------------------------------------------------
   !
   real function phi(fl,r)
   !
   implicit none
   !
   real,         intent(in) :: r
   character(*), intent(in) :: fl
   !
   reaL :: limiter
   !
   select case(fl)
     !
     case('donor-cell')
       limiter = 0.e0
     case('Lax-Wendroff')
       limiter = 1.e0
     case('Beam-Warming')
       limiter = r
     case('Fromm')
       limiter = 0.5e0*(1.e0+r)
     case('minmod')
       limiter = minmod(1.e0,r)
     case('superbee')
       limiter = max(0.e0,min(1.e0,2.e0*r),min(2.e0,r))
     case('hyperbee')
       limiter = 1.0!hyperbee(r)
     case('MC')
       limiter = max(0.e0,min(0.5e0*(1.e0+r),2.e0,2.e0*r))
     case('van Leer')
       limiter = (r+abs(r))/(1.e0+abs(r))
     case('van Albada 1')
       limiter = (r*r+r)/(r*r+1.e0)
     case('van Albada 2')
       limiter = (2.e0*r)/(r*r+1.e0)
     case default
       limiter = 0.e0
   end select
   !
   phi = limiter
   !
   end function phi
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine roe_average(q,qa,w,n)
      !
      ! This subroutine computes Roe's average qa
      ! of a quantity a using the weights w.
      !
      implicit none
      !
      integer,               intent(in)  :: n
      real,    dimension(n), intent(in)  :: q,w
      real,    dimension(n), intent(out) :: qa
      !
      integer :: i
      !
      do i=2,n-1
        qa(i) = (q(i-1)*w(i-1)+q(i)*w(i))/(w(i-1)+w(i))
      enddo
      !
      qa(1) = qa(2)
      qa(n) = qa(n-1)
      !
   end subroutine roe_average
   !
   !----------------------------------------------------------------------------------------------
   !
   real function slope_limiter(slope,q1,q2,q3,q4,v)
      !
      implicit none
      !
      character(80)    :: slope
      real :: q1,q2,q3,q4
      real :: v
      !
      select case(slope)
         !
      end select
      !
      slope_limiter = 0.d0
      !
   end function slope_limiter
   !
   !----------------------------------------------------------------------------------------------
   !
   real function minmod(a,b)
      !
      implicit none
      !
      real :: a,b,c
      !
      if(a*b.gt.0.d0) then
         if (abs(a).lt.abs(b)) then 
             c = a
         else
             c = b
         endif
      else
         c = 0
      endif
      !
      minmod = c 
      !
   end function minmod
   !
   !----------------------------------------------------------------------------------------------
   !
end module Hydro
