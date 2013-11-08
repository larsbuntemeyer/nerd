!
subroutine Io_init
   !
   use Io_data
   !
   implicit none
   !
   write(*,*) '----- Io_init ---------------------'
   !
   outdata_unit = 2 
   outdata_base = 'output.dat'
   open(unit=outdata_unit,file=outdata_base)
   !
   write(*,*) 'output data file:',outdata_base
   write(*,*) '----- Io_init done ----------------'
   ! 
end subroutine Io_init
