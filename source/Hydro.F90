!
!
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
   private :: hydro1D, fill_guardcells
   !
contains
   !
   !
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
   !
   !
   subroutine Hydro_solve(dt)
      !
      implicit none
      !
      double precision, intent(in) :: dt
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
   !
   !
   subroutine hydro1D(dt)
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
      double precision, dimension(nvar,ibg:ieg+1,jbg:jeg+1,kbg:keg+1) :: flux
      double precision, dimension(nvar,ibg:ieg  ,jbg:jeg  ,kbg:keg  ) :: state
      !
      double precision :: rho,ei,ekin,vtot2
      double precision :: fx,fy,fz 
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
               state(1,i,j,k) = rho 
               state(2,i,j,k) = rho * u(i,j,k)
               state(3,i,j,k) = rho * v(i,j,k)
               state(4,i,j,k) = rho * w(i,j,k)
               state(5,i,j,k) = rho * (eint(i,j,k) + ekin)
               !
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
               !uf(i,j,k) = 0.5d0*(state(2,i,j,k)   / state(1,i,j,k) +  &
               !                    state(2,i-1,j,k) / state(1,i-1,j,k)) 
               !vf(i,j,k) = 0.5d0*(state(3,i,j,k)   / state(1,i,j,k) +  &
               !                    state(3,i,j-1,k) / state(1,i,j-1,k)) 
               !wf(i,j,k) = 0.5d0*(state(4,i,j,k)   / state(1,i,j,k) +  &
               !                    state(4,i,j,k-1) / state(1,i,j,k-1)) 
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
               if(uf(i,j,k) < 0.d0) then
                  flux(1,i,j,k) = uf(i,j,k) *  state(1,i,j,k) 
                  flux(2,i,j,k) = uf(i,j,k) *  state(2,i,j,k) + pres(i,j,k)
                  flux(3,i,j,k) = uf(i,j,k) *  state(3,i,j,k)
                  flux(4,i,j,k) = uf(i,j,k) *  state(4,i,j,k)
                  flux(5,i,j,k) = uf(i,j,k) * (state(5,i,j,k) + pres(i,j,k))
               else 
                  flux(1,i,j,k) = uf(i,j,k) *  state(1,i-1,j,k) 
                  flux(2,i,j,k) = uf(i,j,k) *  state(2,i-1,j,k) + pres(i,j,k)
                  flux(3,i,j,k) = uf(i,j,k) *  state(3,i-1,j,k)
                  flux(4,i,j,k) = uf(i,j,k) *  state(4,i-1,j,k)
                  flux(5,i,j,k) = uf(i,j,k) * (state(5,i-1,j,k) + pres(i,j,k))
               endif
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
                  state(ivar,i,j,k) = state(ivar,i,j,k) +         &
                                      dt/dx * (flux(ivar,i,j,k) - flux(ivar,i+1,j,k))
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
               rho = state(1,i,j,k)
               dens(i,j,k) = rho 
               !u(i,j,k)   = state(2,i,j,k) / rho
               v(i,j,k)   = state(3,i,j,k) / rho
               w(i,j,k)   = state(4,i,j,k) / rho
               vtot2       = u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2
               ekin        = 0.5d0 * vtot2
               ener(i,j,k) = state(5,i,j,k) / rho
               eint(i,j,k) = ener(i,j,k) - ekin
            enddo
         enddo
      enddo 
      !
   end subroutine hydro1D
   !
   !
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
  !    do j=ibg,ieg
  !       do i=jbg,nguard
  !          do k=kbg,keg
  !             !
  !             dens(i,j,k) = dens(i,jb,k)
  !             eint(i,j,k) = eint(i,jb,k)
  !             u(i,j,k) = u(i,jb,k)
  !             v(i,j,k) = v(i,jb,k)
  !             w(i,j,k) = w(i,jb,k)
  !             !
  !          enddo
  !       enddo
  !    enddo 
  !    do j=ny+k2d*(nguard+1),ny+2*k2d*nguard
  !       do i=1,nx+2*nguard
  !          do k=1,nz+2*k3d*nguard
  !             !
  !             dens(i,j,k) = dens(i,ny+k2d*nguard,k)
  !             eint(i,j,k) = eint(i,ny+k2d*nguard,k)
  !             u(i,j,k) = u(i,ny+k2d*nguard,k)
  !             v(i,j,k) = v(i,ny+k2d*nguard,k)
  !             w(i,j,k) = w(i,ny+k2d*nguard,k)
  !             !
  !          enddo
  !       enddo
  !    enddo 
  !    do k=1,1+k3d*nguard
  !       do i=1,nx+2*nguard
  !          do j=1,ny+2*k2d*nguard
  !             !
  !             dens(i,j,k) = dens(i,j,k3d*nguard+1)
  !             eint(i,j,k) = eint(i,j,k3d*nguard+1)
  !             u(i,j,k) = u(i,j,k3d*nguard+1)
  !             v(i,j,k) = v(i,j,k3d*nguard+1)
  !             w(i,j,k) = w(i,j,k3d*nguard+1)
  !             !
  !          enddo
  !       enddo
  !    enddo 
  !    do k=nz+k3d*(nguard+1),nz+2*k3d*nguard
  !       do i=1,nx+2*nguard
  !          do j=1,ny+2*k3d*nguard
  !             !
  !             dens(i,j,k) = dens(i,j,nz+k3d*nguard)
  !             eint(i,j,k) = eint(i,j,nz+k3d*nguard)
  !             u(i,j,k) = u(i,j,nz+k3d*nguard)
  !             v(i,j,k) = v(i,j,nz+k3d*nguard)
  !             w(i,j,k) = w(i,j,nz+k3d*nguard)
  !             !
  !          enddo
  !       enddo
  !    enddo 
   end subroutine fill_guardcells
   !
end module Hydro
