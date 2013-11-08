!
module Database
!
   implicit none
!
!  number of state-variables
!
   integer, save :: nvar
!
!  the state-variables
!
   double precision, dimension(:,:,:), allocatable :: dens
   double precision, dimension(:,:,:), allocatable :: eint
   double precision, dimension(:,:,:), allocatable :: vx
   double precision, dimension(:,:,:), allocatable :: vy
   double precision, dimension(:,:,:), allocatable :: vz
!
!  derived variables
!
   double precision, dimension(:,:,:), allocatable :: ener
   double precision, dimension(:,:,:), allocatable :: pres
   double precision, dimension(:,:,:), allocatable :: temp
!
!  the velocity field at the cell interfaces
!
   double precision, dimension(:,:,:), allocatable :: vxf
   double precision, dimension(:,:,:), allocatable :: vyf
   double precision, dimension(:,:,:), allocatable :: vzf
!
   interface
      subroutine Database_init
         implicit none
      end subroutine
   end interface
!
end module Database
!
