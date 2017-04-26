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
   !  fields for stratified model
   !
   real, dimension(:,:,:,:), allocatable :: t 
   real, dimension(:,:,:,:), allocatable :: vx
   real, dimension(:,:,:,:), allocatable :: vy
   real, dimension(:,:,:,:), allocatable :: vz
   real, dimension(:,:,:,:), allocatable :: qd
   real, dimension(:,:,:,:), allocatable :: qw
   real, dimension(:,:,:,:), allocatable :: qi
   real, dimension(:,:,:,:), allocatable :: dwdt
   real, dimension(:,:,:,:), allocatable :: fi
   real, dimension(:,:,:,:), allocatable :: pint
   real, dimension(:,:,:,:), allocatable :: hydrop
   real, dimension(:,:,:,:), allocatable :: hydrodp
   real, dimension(:,:,:,:), allocatable :: hydromp
   !
   real, dimension(:,:,:), allocatable :: etas 
   real, dimension(:,:,:), allocatable :: tg 
   real, dimension(:,:,:), allocatable :: qdb
   real, dimension(:,:,:), allocatable :: ps
   !
   real, dimension(:,:)  , allocatable :: fc 
   real, dimension(:,:)  , allocatable :: fib
   real, dimension(:,:)  , allocatable :: tmch
   !
   !--------------------------------------------------------
end module mo_database
!
