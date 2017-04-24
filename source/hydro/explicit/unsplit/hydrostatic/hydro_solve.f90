!
!----------------------------------------------------------------------------------------------
!
subroutine hydro_solve(dt)
 !
 use mo_grid, only: ndim,ib,ie,nx,ibg,ieg,jbg,jeg,kbg,keg
 use mo_database
 use mo_parameters, only: nvar
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
      call fill_guardcells_1D(dens(:,j,k),pres(:,j,k),eint(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),ibg,ieg,ieg,bc)
   enddo  
 enddo 
 if(ndim>1) then 
   do k=kbg,keg
     do i=ibg,ieg
        call fill_guardcells_1D(dens(i,:,k),pres(i,:,k),eint(i,:,k),u(i,:,k),v(i,:,k),w(i,:,k),jbg,jeg,jeg,bc)
     enddo  
   enddo
 endif
 call hydro3D(dt)
 !do k=kb,ke
 !  do j=jb,je
 !    !call compute_xflux(dens(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),eint(:,j,k),pres(:,j,k), &
 !    !                   nx+2*nguard,nvar,ib,ie,F_l,F_r)
 !    call fill_guardcells_1D(dens(:,j,k),pres(:,j,k),eint(:,j,k),u(:,j,k),ibg,ieg,ieg,2)
 !    call sw
 !  enddo  
 !enddo  
 !
 !
end subroutine hydro_solve
!
