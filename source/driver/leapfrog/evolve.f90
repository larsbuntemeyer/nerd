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
current_time = 0.d0
step = 0
dt = dtini
dt2 = 2.0*dt
ed2dt = 1.0/dt2
dtdeh = dt/3600.0 
!
nold = 3; nnow=1; nnew = 2
nold2 = 2; nnow2 = 1;
!
!
do while(current_time<t_max)
   !
   nsp=nold; nold=nnow; nnow=nnew; nnew=nsp;
   nold2 = 3-nold2; nnow2=3-nnow2;
   !
   write(*,'(I5,6D18.8)') step,current_time,dt, sum(ener(ib:ie,jb:je,kb:ke)), &
                                                sum(dx*dens(ib:ie,jb:je,kb:ke)), &
                                                maxval(abs(u(ib:ie,jb:je,kb:ke))), &
                                                maxval(abs(v(ib:ie,jb:je,kb:ke)))
   !
   call eos
   ! 
   !call hydro_solve(dt)
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
