
subroutine init_driver
 !
 implicit none
 !
 write(*,*) '----- Driver_init_nerd-------------'
 !
 !call RuntimeParameters_init
 !call Grid_init
 !call Database_init
 !call Hydro_init
 !call Io_init
 !!call Simulation_init_domain
 !call Simulation_init_advect
 !call Eos_gamma 
 !
 !write(*,*) 'driver parameters:'
 !write(*,*) 'n_max:', n_max
 !write(*,*) 't_max:', t_max
 !
 write(*,*) '----- Driver_init_nerd done--------'
   !
end subroutine init_driver
