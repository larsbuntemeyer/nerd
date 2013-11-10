!
!
!
program nerd
   !
   use Driver
   !
   implicit none
   !
   call write_banner
   !
   call Driver_init
   call Driver_evolve
   call Driver_finish
   !
   write(*,*) ''   
   write(*,*) 'finished NERD'   
   write(*,*) ''   
   !
contains 
   !
   subroutine write_banner
      !
      implicit none 
      !
      write(*,*) ''
      write(*,*) ''
      write(*,*) '====================================================='
      write(*,*) '               _   ____________  ____                ' 
      write(*,*) '              / | / / ____/ __ \/ __ \               '
      write(*,*) '             /  |/ / __/ / /_/ / / / /               '
      write(*,*) '            / /|  / /___/ _, _/ /_/ /                '
      write(*,*) '           /_/ |_/_____/_/ |_/_____/                 '
      write(*,*) '                                                     '
      write(*,*) '====================================================='
      write(*,*) ''
      write(*,*) ''
      !
   end subroutine write_banner
   !
end program nerd
!
