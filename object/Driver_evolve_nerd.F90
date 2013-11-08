!
subroutine Driver_evolve_nerd
!
   use RuntimeParameters, only: n_max
!
   implicit none
!
   integer          :: timestep
   double precision :: current_time
!
   do timestep=0,n_max
      write(*,*) 'timestep:', timestep,                          &
                 'current time:', current_time
   enddo
!
end subroutine Driver_evolve_nerd
!
