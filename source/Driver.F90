module Driver
!
   use Driver_data
!
   implicit none
!
   interface
      subroutine Driver_init_nerd
         implicit none
      end subroutine Driver_init_nerd
   end interface
!
   interface
      subroutine Driver_evolve_nerd
         implicit none
      end subroutine Driver_evolve_nerd
   end interface
!
   interface
      subroutine Driver_finish_nerd
         implicit none
      end subroutine Driver_finish_nerd
   end interface
!
end module Driver
