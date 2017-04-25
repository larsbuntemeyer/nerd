
subroutine init_namelist 
   !
   use mo_namelist
   !
   implicit none
   !
   write(*,*) '----- Runtime_parameters_init -----'
   write(*,*) 'initializing runtime parameters'
   !
   xmin=-10.
   ymin=-10.
   zmin=0.
   xmax=10.0
   ymax=10.0
   zmax=10.0
   !
   cfl = 0.8
   dtmin = 1.d-10
   dtmax = 5.e-2
   t_max = 0. !0.1
   dtini = 1.e-5
   n_max = 100000
   gamma = 7.d0/5.d0
   mu_mol = 1.0
   lptop0 = .false.
   !
   output_interval = 10
   !
   fl  = 'donor-cell'
   !
   akbk_file = 'akbk_bubble2.txt'
   !
   lptop0 = .TRUE.
   lhdiff2 = .TRUE.
   laistep = .FALSE.
   ldivdamp = .TRUE.
   !
   write(*,*) '----- Runtime_parameters_init done-'
   !
end subroutine init_namelist 
