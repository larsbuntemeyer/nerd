!
subroutine Driver_init_nerd
   !
   use Driver_data
   use RuntimeParameters
   use Grid
   use Database
   use Hydro
   !
   implicit none
   !
   write(*,*) '----- Driver_init_nerd ------------'
   write(*,*) 'driver parameters:'
   write(*,*) 'n_max:', n_max
   write(*,*) 't_max:', t_max
   call RuntimeParameters_init
   call Grid_init
   call Database_init
   call Hydro_init
   !
   write(*,*) '----- Driver_init_nerd done--------'
   !
end subroutine Driver_init_nerd
!
