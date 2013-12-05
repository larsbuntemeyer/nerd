!
!
!
module Driver
!
   implicit none
!
   public :: Driver_init, Driver_evolve,   &
             Driver_finish
!
contains
!
!
!
subroutine Driver_init
   !
   use RuntimeParameters
   use Grid
   use Database
   use Hydro
   use Simulation
   use Io
   use Eos
   !
   implicit none
   !
   write(*,*) '----- Driver_init_nerd-------------'
   write(*,*) 'driver parameters:'
   write(*,*) 'n_max:', n_max
   write(*,*) 't_max:', t_max
   !
   call RuntimeParameters_init
   call Grid_init
   call Database_init
   call Hydro_init
   call Io_init
   call Simulation_init_domain
   call Eos_gamma 
   !
   write(*,*) '----- Driver_init_nerd done--------'
   !
end subroutine Driver_init
!
!
!
subroutine Driver_evolve
   !
   use RuntimeParameters, only: n_max,cfl
   use Hydro
   use Io 
   use Eos
   use Database 
   !
   implicit none
   !
   integer          :: step
   double precision :: current_time
   double precision :: dt, vmax 
   !
   write(*,*) '----- Driver_evolve_nerd ----------'
   !
   current_time = 0.d0
   !
   do step=1,n_max
      !
      vmax = maxval(u)
      if(vmax>0.d0) then
         dt = cfl*dx/vmax
      else
         dt = 0.d0
      endif
      !
      current_time = current_time + dt  
      !
      call Hydro_solve(dt) 
   !   call Eos_gamma 
      !
      write(*,*) 'step:', step,                    &
                 'current time:', current_time,    &
                 'timestep:',dt
   enddo
   !
   write(*,*) '----- Driver_evolve_nerd done -----'
   !
end subroutine Driver_evolve
!
!
!
subroutine Driver_finish
   !
   use Io 
   !
   implicit none
   !
   write(*,*) '----- Driver_finish_nerd ----------'
   write(*,*) 'writing to file...'
   !
   call Io_write_to_file
   !
   write(*,*) 'writing to file... done'
   write(*,*) '----- Driver_finish_nerd done------'
   !
end subroutine Driver_finish
!
end module Driver
