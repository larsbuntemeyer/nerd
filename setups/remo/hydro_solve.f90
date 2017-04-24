!
!----------------------------------------------------------------------------------------------
!
subroutine hydro_solve(dt)
 !
 use mo_grid, only: ib,ie,nx,ibg,ieg,jbg,jeg,kbg,keg
 use mo_database
 use mo_parameters, only: nvar,ndim
 use mo_namelist, only: bc
 use mo_hydro, only: fill_guardcells_1D
 !
 implicit none
 !
 real, intent(inout) :: dt
 integer :: i,j,k
 real,dimension(nvar,nx+1) :: F_l,F_r
 !
 if(ndim>2) then
    write(*,*) '------------------------'
    write(*,*) 'ERROR in Hydro_solve:'
    write(*,*) 'only 1D and 2D possilbe!'
    write(*,*) '------------------------'
    stop
 endif
 !
 !call fill_guardcells
 !
 do k=kbg,keg
   do j=jbg,jeg
      call fill_guardcells_1D(dens(:,j,k),pres(:,j,k),eint(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),ibg,ieg,ieg,2)
   enddo  
 enddo 
 if(ndim>1) then 
   do k=kbg,keg
     do i=ibg,ieg
        call fill_guardcells_1D(dens(i,:,k),pres(i,:,k),eint(i,:,k),u(i,:,k),v(i,:,k),w(i,:,k),jbg,jeg,jeg,3)
     enddo  
   enddo
 endif
 !
 !
 call hydro3D(dt)
 !
 !
end subroutine hydro_solve
