!
!
!
module RuntimeParameters
!
   implicit none
!
   public :: RuntimeParameters_init
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
contains
!
!
!
subroutine RuntimeParameters_init
   !
   write(*,*) '----- Runtime_parameters_init -----'
   write(*,*) 'initializing runtime parameters'
   !
   xmin=0.d0
   ymin=0.d0
   zmin=0.d0
   xmax=1.d0
   ymax=0.d0
   zmax=0.d0
   !
   cfl = 0.8
   n_max = 10
   gamma = 5.d0/3.d0
   mu_mol = 1.0
   !
   write(*,*) '----- Runtime_parameters_init done-'
   !
end subroutine RuntimeParameters_init
!
end module RuntimeParameters
