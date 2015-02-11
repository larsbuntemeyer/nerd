!
module Database
   !
   use Grid
   !
   implicit none
   !
   public :: Database_init
   !
   !  number of state-variables
   !
   integer, save :: nvar
   !
   !  the state-variables
   !
   real, dimension(:,:,:), allocatable :: dens
   real, dimension(:,:,:), allocatable :: eint
   real, dimension(:,:,:), allocatable :: u
   real, dimension(:,:,:), allocatable :: v
   real, dimension(:,:,:), allocatable :: w
   !
   !  derived variables
   !
   real, dimension(:,:,:), allocatable :: ener
   real, dimension(:,:,:), allocatable :: pres
   real, dimension(:,:,:), allocatable :: temp
   !
   !  the velocity field at the cell interfaces
   !
   real, dimension(:,:,:), allocatable :: uf
   real, dimension(:,:,:), allocatable :: vf
   real, dimension(:,:,:), allocatable :: wf
   !
   !
   !--------------------------------------------------------
contains
   !
   !
   !
   subroutine Database_init
      !
      implicit none
      !
      write(*,*) '----- Database_init ---------------'
      !
      nvar = 5 
      !
      write(*,*) 'allocating'
      allocate(dens(ibg:ieg,jbg:jeg,kbg:keg))
      allocate(eint(ibg:ieg,jbg:jeg,kbg:keg))
      !
      allocate(u(ibg:ieg,jbg:jeg,kbg:keg))
      allocate(v(ibg:ieg,jbg:jeg,kbg:keg))
      allocate(w(ibg:ieg,jbg:jeg,kbg:keg))
      !
      allocate(uf(ibg:ieg+1,jbg:jeg+1,kbg:keg+1))
      allocate(vf(jbg:ieg+1,jbg:jeg+1,kbg:keg+1))
      allocate(wf(ibg:ieg+1,jbg:jeg+1,kbg:keg+1))
      !
      allocate(ener(nx+2*nguard,ny+2*k2d*nguard,nz+2*k3d*nguard))
      allocate(pres(nx+2*nguard,ny+2*k2d*nguard,nz+2*k3d*nguard))
      allocate(temp(nx+2*nguard,ny+2*k2d*nguard,nz+2*k3d*nguard))
      !
      write(*,*) '----- Database_init done ----------'
      !
   end subroutine Database_init
   !
end module Database
!
