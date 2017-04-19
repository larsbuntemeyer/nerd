!
!----------------------------------------------------------------------------------------------
!
subroutine hydro_solve(dt)
 !
 use Grid, only: ndim,ib,ie
 use Database
 use RuntimeParameters, only: bc
 !
 implicit none
 !
 real, intent(inout) :: dt
 integer :: i,j,k
 real,dimension(nvar,nx+1) :: F_l,F_r
 !
 if(ndim>1) then
    write(*,*) '------------------------'
    write(*,*) 'ERROR in Hydro_solve:'
    write(*,*) 'only 1D possilbe!'
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
 !do j=jbg,jeg
 !  do i=ibg,ieg
 !     call fill_guardcells_1D(dens(i,j,:),pres(i,j,:),eint(i,j,:),u(i,j,:),v(i,j,:),w(i,j,:),kbg,keg,keg,2)
 !  enddo  
 !enddo  
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
