module Hydro
!
   use Hydro_data
!
   implicit none
!
   interface
      subroutine Hydro_init
      end subroutine Hydro_init
   end interface
   interface
      subroutine Hydro_solve(dt)
         implicit none
         double precision, intent(in) :: dt
      end subroutine Hydro_solve
   end interface
   interface
      double precision function Hydro_cfl_timestep(cfl)
         implicit none
         double precision :: cfl
      end function Hydro_cfl_timestep
   end interface
!
end module Hydro
