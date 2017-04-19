!
module mo_database
   !
   !
   implicit none
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
end module mo_database
!
