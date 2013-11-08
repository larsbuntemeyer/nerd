!
subroutine Driver_evolve_nerd
   !
   use RuntimeParameters, only: n_max,cfl
   use Hydro
   use Io 
   use Eos 
   !
   implicit none
   !
   integer          :: step
   double precision :: current_time
   double precision :: dt 
   !
   write(*,*) '----- Driver_evolve_nerd ----------'
   !
   do step=1,n_max
      !
      dt = 1.e-6!hYDRO_CFL_TIMEstEP(CFL)
      current_time = current_time + dt  
      !
      call Hydro_solve(dt) 
      call Eos_gamma 
      !
      write(*,*) 'step:', step,                    &
                 'current time:', current_time,    &
                 'timestep:',dt
   enddo
   !
   write(*,*) '----- Driver_evolve_nerd done -----'
   !
end subroutine Driver_evolve_nerd
!
