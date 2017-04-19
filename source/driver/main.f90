program main
 !
 implicit none
 !
 call write_banner
 !
 call init_driver
 call evolve
 call finish_driver
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
end program main
