subroutine evolve
!
use mo_namelist
use mo_database
use mo_grid, only: ib,ie,jb,je,kb,ke,dx
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
   call eos
   ! 
   call hydro_solve(dt)
   !
   current_time = current_time + dt  
   step         = step+1
   ! 
   !
enddo
!
write(*,*) '----- Driver_evolve_nerd done -----'
!
end subroutine evolve
