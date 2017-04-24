subroutine init_driver
   !
   use mo_driver
   !
   implicit none
   !
   !write(*,*) '----- Driver_init_nerd-------------'
   !!
   call dump_parameters
   call init_namelist 
   call init_grid
   call init_database
   call init_hydro
   call init_domain
   call init_io
   !!
   dt2 = 2.0*dt
   ed2dt = 1.0/dt2
   !write(*,*) '----- Driver_init_nerd done--------'
   !
end subroutine init_driver
