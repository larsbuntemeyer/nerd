
subroutine init_namelist 
   !
   use mo_namelist
   !
   implicit none
   !
   write(*,*) '----- Runtime_parameters_init -----'
   write(*,*) 'initializing runtime parameters'
   !
   xmin=0.0
   ymin=0.0
   zmin=0.0
   xmax=1.0
   ymax=0.0
   zmax=0.0
   !
   cfl = 0.8
   dtmin = 1.d-10
   dtmax = 0.5
   t_max = 0.2 !0.1
   dtini = 1.e-3
   n_max = 100000
   gamma = 7.d0/5.d0
   mu_mol = 1.0
   !
   fl  = 'superbee'
   !
   write(*,*) '----- Runtime_parameters_init done-'
   !
end subroutine init_namelist 
