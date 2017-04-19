!
! Hydro module
!
module mo_hydro
   !
   implicit none
   !
contains
   !
   !
   !
   subroutine fill_guardcells_1D(dens,pres,eint,u,v,w,ib,ie,n,bc)
     use mo_grid, only: nguard
     use mo_namelist, only: outflow, periodic
     implicit none
     real ,intent(inout), dimension(n) :: dens,pres,eint,u,v,w
     integer, intent(in) :: ib,ie,n 
     integer, intent(in) :: bc
     integer :: i,j,k 
     select case(bc)
       case(outflow) 
         do i=ib,ib+nguard-1
           dens(i) = dens(ib+nguard)
           pres(i) = pres(ib+nguard)
           eint(i) = eint(ib+nguard)
           u(i)    = u(ib+nguard)
           v(i)    = v(ib+nguard)
           w(i)    = w(ib+nguard)
         enddo 
         do i=ie-nguard+1,ie
           dens(i) = dens(ie-nguard)
           pres(i) = pres(ie-nguard)
           eint(i) = eint(ie-nguard)
           u(i)    = u(ie-nguard)
           v(i)    = v(ie-nguard)
           w(i)    = w(ie-nguard)
         enddo
       case(periodic) 
         do i=ib,ib+nguard-1
           dens(i) = dens(ie-nguard-i+1)
           pres(i) = pres(ie-nguard-i+1)
           eint(i) = eint(ie-nguard-i+1)
           u(i)    = u(ie-nguard-i+1)
           v(i)    = v(ie-nguard-i+1)
           w(i)    = w(ie-nguard-i+1)
         enddo 
         do i=ie-nguard+1,ie
           dens(i) = dens(i-n+nguard+1)
           pres(i) = pres(i-n+nguard+1)
           eint(i) = eint(i-n+nguard+1)
           u(i)    = u(i-n+nguard+1)
           v(i)    = v(i-n+nguard+1)
           w(i)    = w(i-n+nguard+1)
         enddo
       case default
     end select 
   end subroutine fill_guardcells_1D
   !
   subroutine fill_guardcells
      !
      use mo_grid
      use mo_database
      !
      implicit none
      !
      integer :: i,j,k,ivar
      !
      ! Fill guardcells x-direction
      !
      do i=ibg,nguard
        do j=jb,je
          do k=kb,ke
            !
            dens(i,j,k) = dens(ib,j,k)
            pres(i,j,k) = pres(ib,j,k)
            eint(i,j,k) = eint(ib,j,k)
            u(i,j,k)    = u(ib,j,k)
            !
          enddo
        enddo
      enddo 
      do i=ie+1,ieg
        do j=jb,je
          do k=kb,ke
            !
            dens(i,j,k) = dens(ie,j,k)
            pres(i,j,k) = pres(ie,j,k)
            eint(i,j,k) = eint(ie,j,k)
            u(i,j,k)    = u(ie,j,k)
            !
          enddo
        enddo
      enddo 
      !
      ! Fill guardcells z-direction
      !
      do k=kbg,nguard
        do j=jb,je
          do i=ib,ie
            !
            dens(i,j,k) = dens(i,j,ke)
            pres(i,j,k) = pres(i,j,ke)
            eint(i,j,k) = eint(i,j,ke)
            w(i,j,k)    = w(i,j,kb)
            !
          enddo
        enddo
      enddo 
      do k=ke+1,keg
        do j=jb,je
          do i=ib,ie
            !
            dens(i,j,k) = dens(i,j,ke)
            pres(i,j,k) = pres(i,j,ke)
            eint(i,j,k) = eint(i,j,ke)
            w(i,j,k)    = w(ie,j,k)
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
      use mo_grid
      use mo_namelist
      !
      implicit none
      !
      real :: q1,q2,q3,q4
      real :: v_face,dt
      real :: r,limiter,flux,theta
      !
      theta = sign(1.0,v_face)
      !
      if(abs(q3-q2).gt.0.0) then
         if(v_face.ge.0.0) then
            r = (q2-q1)/(q3-q2)
         else 
            r = (q4-q3)/(q3-q2)
         endif
      else
         r = 0.0
      endif
      !
      limiter = phi(fl,r)
      !
      flux = 0.5d0*v_face*((1.e0+theta)*q2+(1.e0-theta)*q3) +  &
             0.5d0*abs(v_face)*(1.e0-abs(v_face*dt/dx))*limiter*(q3-q2)
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
       limiter = 0.0
     case('Lax-Wendroff')
       limiter = 1.0
     case('Beam-Warming')
       limiter = r
     case('Fromm')
       limiter = 0.5*(1.0+r)
     case('minmod')
       limiter = minmod(1.0,r)
     case('superbee')
       limiter = max(0.0,min(1.0,2.0*r),min(2.0,r))
     case('hyperbee')
       limiter = 1.0!hyperbee(r)
     case('MC')
       limiter = max(0.0,min(0.5*(1.0+r),2.0,2.0*r))
     case('van Leer')
       limiter = (r+abs(r))/(1.0+abs(r))
     case('van Albada 1')
       limiter = (r*r+r)/(r*r+1.0)
     case('van Albada 2')
       limiter = (2.0*r)/(r*r+1.0)
     case default
       limiter = 0.0
   end select
   !
   phi = limiter
   !
   end function phi
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine compute_slope(a,l,n,nvar,r)
     !
     implicit none
     !
     integer, intent(in) :: n,nvar
     real, dimension(nvar,n), intent(in)  :: a,l
     real, dimension(nvar,n), intent(out) :: r
     !
     integer :: i,ivar
     !
     do ivar=1,nvar
       do i=2,n-1
         if(abs(a(ivar,i)) > 0.0) then
            if(l(ivar,i) >= 0.0) then
              r(ivar,i) = a(ivar,i-1)/a(ivar,i)
            else
              r(ivar,i) = a(ivar,i+1)/a(ivar,i)
            endif
         else
            r(ivar,i) = 0.0
         endif
       enddo
       r(ivar,1) = r(ivar,3)
       r(ivar,2) = r(ivar,3)
       r(ivar,n) = r(ivar,n-2)
       r(ivar,n-1) = r(ivar,n-2)
     enddo
     !
   end subroutine compute_slope
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
   real function cfl_timestep(n,nvar,l,dx)
   !
   use mo_namelist 
   
   implicit none
   !
   integer, intent(in) :: n,nvar
   real, intent(in), dimension(nvar,n) :: l
   real, intent(in)    :: dx
   !
   integer :: i,ivar
   real :: dt,lmax,lmin
   real, dimension(n-1) :: dti
   !
   do i=1,n-1
     lmax = maxval(l(:,i))
     if(lmax.lt.0.0) lmax=0.0
     lmin = minval(l(:,i+1))
     if(lmin.gt.0.0) lmin=0.0
     if(lmax-lmin > 0.0) then
       dti(i) = dx/(lmax-lmin)
     else
       dti(i) = dtmax
     endif
   enddo 
   !
   dt = cfl * minval(dti)
   !
   if(dt < dtmin) then
      write(*,*) 'WARNING: cfl timestep is less than minimum timestep'
      write(*,*) 'using dtmin'
      dt = dtmin
   endif
   !
   cfl_timestep = dt
   !
   end function cfl_timestep
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
end module mo_hydro 
