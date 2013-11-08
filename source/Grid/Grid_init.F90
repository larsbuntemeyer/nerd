!
subroutine grid_init
!
   use grid_data
! 
   implicit none
!
   ndim = 1
   nx   = 100
   ny   = 1
   nz   = 1
!
   write(*,*) '-----------------------'
   write(*,*) 'grid_init'
   write(*,*) 'grid parameters:'
   write(*,*) 'ndim:', ndim
   write(*,*) 'nx,ny,nz:', nx,ny,nz
   write(*,*) '-----------------------'
!
end subroutine grid_init
!
