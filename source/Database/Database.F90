module Database
!
   use Grid
!
   implicit none
!
   double precision, dimension(:,:,:), allocatable :: rho
   double precision, dimension(:,:,:), allocatable :: temp
   double precision, dimension(:,:,:), allocatable :: eint
   double precision, dimension(:,:,:), allocatable :: ener
   double precision, dimension(:,:,:), allocatable :: pres
!
!
!
   contains
!
   subroutine Database_init
!
      implicit none
!
      write(*,*) '-----------------------'
      write(*,*) 'database_init'
      write(*,*) 'allocating'
      allocate(rho(nx,ny,nz))
      allocate(temp(nx,ny,nz))
      allocate(eint(nx,ny,nz))
      allocate(ener(nx,ny,nz))
      allocate(pres(nx,ny,nz))
      write(*,*) 'success'
      write(*,*) '-----------------------'
!
      
!
   end subroutine Database_init
!
end module Database
