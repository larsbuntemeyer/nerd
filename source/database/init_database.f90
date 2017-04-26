!
!
!
subroutine init_database
   !
   use mo_database
   use mo_grid, only: ibg,ieg,jbg,jeg,kbg,keg,ie=>ieg,je=>jeg,ke=>keg,ke1=>keg1, &
                      nx,ny,nz
   use mo_parameters, only: nguard,k2d,k3d
   !
   implicit none
   !
   write(*,*) '----- Database_init ---------------'
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
   allocate(ps(ie,je,3))
   ALLOCATE(T(IE,JE,KE,3))
   ALLOCATE(VX(IE,JE,KE,3))
   ALLOCATE(VY(IE,JE,KE,3))
   ALLOCATE(VZ(IE,JE,KE1,3))
   ALLOCATE(QD(IE,JE,KE,3))
   ALLOCATE(QW(IE,JE,KE,3))
   ALLOCATE(QI(IE,JE,KE,3))
   ALLOCATE(FI(IE,JE,KE,2))
   ALLOCATE(DWDT(IE,JE,KE,3))
   ALLOCATE(PINT(IE,JE,KE1,3))
   ALLOCATE(ETAS(IE,JE,KE1))
   ALLOCATE(FIB(IE,JE))
   ALLOCATE(TG(IE,JE,3))
   ALLOCATE(QDB(IE,JE,3))
   ALLOCATE(TMCH(IE,JE))
   ALLOCATE(FC(IE,JE))
   ALLOCATE(HYDROP(IE,JE,KE1,3), HYDRODP(IE,JE,KE,3),HYDROMP(IE,JE,KE,3))
   write(*,*) '----- Database_init done ----------'
   !
end subroutine init_database
!
