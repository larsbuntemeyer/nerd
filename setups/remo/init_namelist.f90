
subroutine init_namelist 
   !
   use mo_namelist
   !
   implicit none
   !
   write(*,*) '----- Runtime_parameters_init -----'
   write(*,*) 'initializing runtime parameters'
   !
   xmin=-1.0
   ymin=-0.01
   zmin=0.
   xmax=3.0
   ymax=0.01
   zmax=0.
   !
   cfl = 0.8
   dtmin = 1.0 
   dtmax = 1.0 
   t_max = 1.0 !0.1
   dtini = 1.0
   n_max = 100000
   gamma = 7.d0/5.d0
   mu_mol = 1.0
   !
   output_interval = 10
   !
   fl  = 'donor-cell'
   !
   akbk_file = 'akbk_27.txt'
   !
   lhdiff2 = .TRUE.
   laistep = .FALSE.
   ldivdamp = .TRUE.
   !
   write(*,*) '----- Runtime_parameters_init done-'
   !
end subroutine init_namelist 
