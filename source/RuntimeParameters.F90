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
   integer, parameter :: outflow=1, periodic=2
   integer, parameter :: bc = periodic
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
   xmin=0.0
   ymin=0.0
   zmin=0.0
   xmax=1.0
   ymax=0.0
   zmax=1.0
   !
   cfl = 0.8
   dtmin = 1.d-10
   dtmax = 0.1
   t_max = 2.0
   dtini = 1.e-5
   n_max = 1000
   gamma = 7.d0/5.d0
   mu_mol = 1.0
   !
   fl  = 'donor-cell'
   !
   write(*,*) '----- Runtime_parameters_init done-'
   !
end subroutine RuntimeParameters_init
!
end module RuntimeParameters
