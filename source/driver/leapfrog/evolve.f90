subroutine evolve
!
use mo_namelist
use mo_database
use mo_parameters, only: ndim
use mo_grid, only: ib,ie,jb,je,kb,ke,dx,dy
use mo_driver
!
implicit none
!
integer :: nsp
!
write(*,*) '----- Driver_evolve_nerd ----------'
!
!
current_time = 0.d0
step = 0
dt = dtini
!
nold = 3; nnow=1; nnew = 2
!
do while(current_time<t_max)
   !
   nsp=nold; nold=nnow; nnow=nnew; nnew=nsp;
   !
   write(*,'(I5,6D18.8)') step,current_time,dt, sum(ener(ib:ie,jb:je,kb:ke)), &
                                                sum(dx*dens(ib:ie,jb:je,kb:ke)), &
                                                maxval(abs(u(ib:ie,jb:je,kb:ke))), &
                                                maxval(abs(v(ib:ie,jb:je,kb:ke)))
   !
   call eos
   ! 
   call hydro_solve(dt)
   !
   current_time = current_time + dt  
   step         = step+1
   ! 
   dt = cfl*minval(dx/((abs(u))))
   !if(ndim>1) dt = min(dt,cfl*dy/(maxval(abs(u))))
   if(dt.gt.dtmax) dt = dtmax
   !
   if(mod(step,output_interval)==0) then
      write(*,*) '-------- writing to file ----------'
      call io_write_to_file
   endif
   !
enddo
!
write(*,*) '----- Driver_evolve_nerd done -----'
!
end subroutine evolve
