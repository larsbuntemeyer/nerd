subroutine init_driver
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
   call init_io
   call init_domain
   !!
   !write(*,*) '----- Driver_init_nerd done--------'
   !
end subroutine init_driver
