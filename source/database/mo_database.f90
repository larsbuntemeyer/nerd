!
module mo_database
   !
   implicit none
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
   real, dimension(:,:,:,:), allocatable :: t 
   real, dimension(:,:,:,:), allocatable :: ux
   real, dimension(:,:,:,:), allocatable :: uy
   real, dimension(:,:,:,:), allocatable :: uz
   real, dimension(:,:,:,:), allocatable :: qd
   real, dimension(:,:,:,:), allocatable :: qw
   real, dimension(:,:,:,:), allocatable :: qi
   real, dimension(:,:,:,:), allocatable :: dwdt
   real, dimension(:,:,:,:), allocatable :: tmch
   real, dimension(:,:,:,:), allocatable :: fi
   real, dimension(:,:,:,:), allocatable :: fib
   real, dimension(:,:,:,:), allocatable :: qdb
   real, dimension(:,:,:,:), allocatable :: pint

   !
   !--------------------------------------------------------
end module mo_database
!
