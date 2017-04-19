!
!
!
subroutine init_database
   !
   use mo_database
   use mo_grid
   !
   implicit none
   !
   write(*,*) '----- Database_init ---------------'
   !
   nvar = 5 
   !
   write(*,*) 'allocating'
   allocate(dens(ibg:ieg,jbg:jeg,kbg:keg))
   allocate(eint(ibg:ieg,jbg:jeg,kbg:keg))
   !
   allocate(u(ibg:ieg,jbg:jeg,kbg:keg))
   allocate(v(ibg:ieg,jbg:jeg,kbg:keg))
   allocate(w(ibg:ieg,jbg:jeg,kbg:keg))
   !
   allocate(uf(ibg:ieg+1,jbg:jeg+1,kbg:keg+1))
   allocate(vf(jbg:ieg+1,jbg:jeg+1,kbg:keg+1))
   allocate(wf(ibg:ieg+1,jbg:jeg+1,kbg:keg+1))
   !
   allocate(ener(nx+2*nguard,ny+2*k2d*nguard,nz+2*k3d*nguard))
   allocate(pres(nx+2*nguard,ny+2*k2d*nguard,nz+2*k3d*nguard))
   allocate(temp(nx+2*nguard,ny+2*k2d*nguard,nz+2*k3d*nguard))
   !
   write(*,*) '----- Database_init done ----------'
   !
end subroutine init_database
!
