!
!
!
module Driver
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
   public  :: Driver_init, Driver_evolve,   &
              Driver_finish
   private :: get_cfl_timestep
   !
contains
!
!
!
subroutine Driver_init
   !
   implicit none
   !
   write(*,*) '----- Driver_init_nerd-------------'
   !
   call RuntimeParameters_init
   call Grid_init
   call Database_init
   call Hydro_init
   call Io_init
   !call Simulation_init_domain
   call Simulation_init_advect
   call Eos_gamma 
   !
   write(*,*) 'driver parameters:'
   write(*,*) 'n_max:', n_max
   write(*,*) 't_max:', t_max
   !
   write(*,*) '----- Driver_init_nerd done--------'
   !
end subroutine Driver_init
!
!
!
subroutine Driver_evolve
   !
   implicit none
   !
   integer   :: step
   real      :: current_time
   real      :: dt, vmax
   !
   write(*,*) '----- Driver_evolve_nerd ----------'
   !
   !
   current_time = 0.d0
   step = 0
   dt = dtini
   !
   do while(current_time<t_max)
      !
      write(*,'(I5,4D18.8)') step,current_time,dt, sum(ener(ib:ie,jb:je,kb:ke)), &
                                                   sum(dx*dens(ib:ie,jb:je,kb:ke))
      !
      call Eos_gamma
      ! 
      call Hydro_solve(dt)
      !
      !write(*,'(I5,4D18.8)') step,current_time,dt, sum(ener), sum(dx*dens)
      ! 
      current_time = current_time + dt  
      step         = step+1
      ! 
      !
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
!
!
real function get_cfl_timestep()
   !
   implicit none
   !
   real :: dt
   real :: vmax
   !
   vmax = maxval(abs(u))
   if(ndim.ge.2) vmax = max(vmax,maxval(abs(v)))
   if(ndim.ge.3) vmax = max(vmax,maxval(abs(w)))
   !
   if(vmax>0.d0) then
      dt = cfl*dx/vmax
   else
      dt = dtmax
   endif
   !
   if(dt.lt.dtmin) then
      write(*,*) 'WARNING: cfl timestep is less than minimum timestep'
      write(*,*) 'using dtmin'
      dt = dtmin
   endif
   !
   get_cfl_timestep = dt
   !
end function get_cfl_timestep
!
end module Driver
