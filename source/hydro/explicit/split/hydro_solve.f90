!
!----------------------------------------------------------------------------------------------
!
subroutine hydro_solve(dt)
 !
 use mo_grid, only: ndim,ib,ie,jb,je,kb,ke,ibg,ieg,jbg,jeg,kbg,keg,dx,dy,dz,nx,ny,nz, nguard
 use mo_database, only: dens,pres,eint,u,v,w,eint,nvar
 use mo_hydro, only: fill_guardcells_1D
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
 do k=kb,ke
   do j=jb,je
     !call compute_xflux(dens(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),eint(:,j,k),pres(:,j,k), &
     !                   nx+2*nguard,nvar,ib,ie,F_l,F_r)
     call fill_guardcells_1D(dens(:,j,k),pres(:,j,k),eint(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),ibg,ieg,ieg,1)
     call sweep1D(0.5*dt,dx,1,dens(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),eint(:,j,k), &
                  pres(:,j,k),nx+2*nguard,nvar) 
   enddo  
 enddo  
 !
 !do j=jb,je
 !  do i=ib,ie
 !    call fill_guardcells_1D(dens(i,j,:),pres(i,j,:),eint(i,j,:),w(i,j,:),v(i,j,:),w(i,j,:),kbg,keg,keg,2)
 !    call sweep1D(dt,dz,1,dens(i,j,:),w(i,j,:),u(i,j,:),v(i,j,:),eint(i,j,:), &
 !                 pres(i,j,:),nz+2*nguard,nvar) 
 !  enddo  
 !enddo
 ! 
 do k=kb,ke
   do j=jb,je
     !call compute_xflux(dens(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),eint(:,j,k),pres(:,j,k), &
     !                   nx+2*nguard,nvar,ib,ie,F_l,F_r)
     call fill_guardcells_1D(dens(:,j,k),pres(:,j,k),eint(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),ibg,ieg,ieg,1)
     call sweep1D(0.5*dt,dx,1,dens(:,j,k),u(:,j,k),v(:,j,k),w(:,j,k),eint(:,j,k), &
                  pres(:,j,k),nx+2*nguard,nvar) 
   enddo  
 enddo  
 !
end subroutine hydro_solve
!
