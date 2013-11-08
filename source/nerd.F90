!
program nerd
   !
   use Driver
   !
   implicit none
   !
   real :: x,y,z,a,b
   integer :: i
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
   !
   call Driver_init_nerd 
   call Driver_evolve_nerd 
   call Driver_finish_nerd 
   !
   write(*,*) ''   
   write(*,*) 'finished NERD'   
   write(*,*) ''   
   !
end program nerd
!
