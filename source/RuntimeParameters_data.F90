
module RuntimeParameters_data
!
   implicit none
!
!  time
!
   integer          , save :: n_max
   double precision , save :: t_max
!
!  more
!
   double precision        :: xmin,xmax
   double precision        :: ymin,ymax
   double precision        :: zmin,zmax
!
   double precision        :: cfl
!
   double precision        :: gamma
   double precision        :: mu_mol
!
end module RuntimeParameters_data
