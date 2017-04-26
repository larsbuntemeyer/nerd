!
!
!
module mo_namelist 
   !
   implicit none
   !
   !  time
   !
   integer          , save :: n_max
   real , save :: t_max
   !
   !  more
   !
   real        :: xmin,xmax
   real        :: ymin,ymax
   real        :: zmin,zmax
   !
   real        :: cfl
   real        :: dtmin
   real        :: dtmax
   real        :: dtini
   !
   real        :: gamma
   real        :: mu_mol
   !
   character(80)           :: fl
   !
   integer, parameter :: outflow=1, periodic=2, reflective=3
   integer, parameter :: bc = reflective
   !
   integer     :: output_interval
   logical     :: lhdiff2
   logical     :: laistep
   logical     :: ldivdamp
   !
   character(256) :: akbk_file
   !
end module mo_namelist
