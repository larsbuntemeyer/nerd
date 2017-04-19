!
module Simulation
   !
   implicit none
   !
   ! Here one can add his own runtime-parameters
   !
   real    :: radius
   real    :: rho_left,rho_right
   real    :: vx_left,vx_right
   !
contains
!
subroutine Simulation_init_domain
   !
   use Grid
   use RuntimeParameters
   use Database
   !
   implicit none
   !
   integer             :: i,j,k 
   real    :: distance
   real    :: xctr,yctr,zctr
   real    :: xsize,ysize,zsize
   real    :: rho_r,rho_l,vx_r,vx_l,eint_l,eint_r,pres_r,pres_l
   real    :: ek,ei,e, border
   integer,parameter :: seed = 86456
   !
   ! Toro Test 1
   !
   rho_l  = 1.0
   vx_l   = 0.75
   pres_l = 1.0
   rho_r  = 0.125
   vx_r   = 0.125
   pres_r = 0.1 
   border = 0.3
   !!
   !! Toro Test 2
   !!
   !rho_l  = 1.0
   !vx_l   = -2.0 
   !pres_l = 0.4
   !rho_r  = 1.0
   !vx_r   = 2.0
   !pres_r = 0.4 
   !border = 0.3
   !
   ! Toro Test 3
   !
   !rho_l  = 1.0
   !vx_l   = 0.0 
   !pres_l = 1000.0
   !rho_r  = 1.0
   !vx_r   = 0.0
   !pres_r = 0.01 
   !border = 0.5
   !
   ! Make internal energy consistent with pressure and density
   !
   eint_l = pres_l/(rho_l*(gamma-1.0))
   eint_r = pres_r/(rho_r*(gamma-1.0))
   !
   xctr  = 0.5*(xmax-xmin)
   yctr  = 0.5*(ymax-ymin)
   zctr  = 0.5*(zmax-zmin)
   xsize = (xmax-xmin)
   ysize = (ymax-ymin)
   zsize = (zmax-zmin)
   !
   write(*,*) '----- init_domain -----------------'
   !
   !do i=1+nguard,nx+nguard
   do i=ib,ie
      do j=jb,je
         do k=kb,ke
            !
            ! Here one should init the physics 
            !
            !distance = (xcCoord(i) - xctr)**2 
            !if(ndim>2) then 
            !   distance = distance + (ycCoord(j) - yctr)**2 
            !endif
            !if(ndim==3) then 
            !   distance = distance + (zcCoord(k) - zctr)**2 
            !endif
            !
            !distance = sqrt(distance)
            !
            !--------------------------------------
            !distance = (xcCoord(i) - xctr)**2 
            !distance = distance + (zcCoord(k) - zctr)**2
            !distance = sqrt(distance)
            !!
            !if(distance < 0.5*border) then
            !  dens(i,j,k) = rho_l
            !  eint(i,j,k) = eint_l
            !else
            !  dens(i,j,k) = rho_r
            !  eint(i,j,k) = eint_r
            !endif
            !v(i,j,k) = 0.0 
            !u(i,j,k) = 0.0 
            !w(i,j,k) = 0.0 
            !ek = 0.0 
            !ener(i,j,k) = eint(i,j,k) + ek
            !--------------------------------------
            !if(zcCoord(k) < border) then
            !  dens(i,j,k) = rho_l
            !  w(i,j,k)    = vx_l
            !  eint(i,j,k) = eint_l
            !else
            !  dens(i,j,k) = rho_r
            !  w(i,j,k)    = vx_r
            !  eint(i,j,k) = eint_r
            !endif
            !!
            !v(i,j,k) = 0.0 
            !u(i,j,k) = 0.0 
            !ek = 0.5*w(i,j,k)**2
            !ener(i,j,k) = eint(i,j,k) + ek
            !--------------------------------------
            if(xcCoord(i) < border) then
              dens(i,j,k) = rho_l
              u(i,j,k)    = vx_l
              eint(i,j,k) = eint_l
            else
              dens(i,j,k) = rho_r
              u(i,j,k)    = vx_r
              eint(i,j,k) = eint_r
            endif
            !
            v(i,j,k) = 0.0 
            w(i,j,k) = 0.0 
            ek = 0.5*u(i,j,k)**2
            ener(i,j,k) = eint(i,j,k) + ek
            !--------------------------------------
            !if(zcCoord(k) > 0.7 .or. zcCoord(k) < 0.3) then
            !  u(i,j,k)    = -0.5 + 1.d-3*(2.0*rand(seed)-1.0)  
            !  dens(i,j,k) = 1.0
            !else
            !  u(i,j,k)    =  0.5 + 1.d-3*(2.0*rand(seed)-1.0)  
            !  dens(i,j,k) = 2.0
            !endif
            !v(i,j,k) = 0.0 
            !w(i,j,k) = 0.0 + 0.1*(2.0*rand(seed)-1.0)  
            !pres(i,j,k) = 2.5
            !eint(i,j,k) = pres(i,j,k)/(dens(i,j,k)*(gamma-1.0))
            !ek = 0.5*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
            !ener(i,j,k) = eint(i,j,k) + ek
            !
         enddo
      enddo
   enddo
   !
   write(*,*) '----- init_domain done ------------'
   !
end subroutine Simulation_init_domain
   !
   !
   !
subroutine Simulation_init_advect
   !
   use Grid
   use RuntimeParameters
   use Database
   !
   implicit none
   !
   integer             :: i,j,k 
   real    :: distance
   real    :: xctr,yctr,zctr
   real    :: xsize,ysize,zsize
   real    :: rho_r,rho_l,vx_r,vx_l,eint_l,eint_r,pres_r,pres_l
   real    :: ek,ei,e, border
   integer,parameter :: seed = 86456
   !
   ! Toro Test 1
   !
   rho_l  = 1.0
   vx_l   = 0.75
   pres_l = 1.0
   rho_r  = 0.125
   vx_r   = 0.0
   pres_r = 0.1 
   border = 0.3
   !!
   !! Toro Test 2
   !!
   !rho_l  = 1.0
   !vx_l   = -2.0 
   !pres_l = 0.4
   !rho_r  = 1.0
   !vx_r   = 2.0
   !pres_r = 0.4 
   !border = 0.3
   !
   ! Toro Test 3
   !
   !rho_l  = 1.0
   !vx_l   = 0.0 
   !pres_l = 1000.0
   !rho_r  = 1.0
   !vx_r   = 0.0
   !pres_r = 0.01 
   !border = 0.5
   !
   ! Make internal energy consistent with pressure and density
   !
   eint_l = pres_l/(rho_l*(gamma-1.0))
   eint_r = pres_r/(rho_r*(gamma-1.0))
   !
   xctr  = 0.5*(xmax-xmin)
   yctr  = 0.5*(ymax-ymin)
   zctr  = 0.5*(zmax-zmin)
   xsize = (xmax-xmin)
   ysize = (ymax-ymin)
   zsize = (zmax-zmin)
   !
   write(*,*) '----- init_domain -----------------'
   !
   !do i=1+nguard,nx+nguard
   do i=ib,ie
      do j=jb,je
         do k=kb,ke
            !
            ! Here one should init the physics 
            !
            !distance = (xcCoord(i) - xctr)**2 
            !if(ndim>2) then 
            !   distance = distance + (ycCoord(j) - yctr)**2 
            !endif
            !if(ndim==3) then 
            !   distance = distance + (zcCoord(k) - zctr)**2 
            !endif
            !! 
            !distance = sqrt(distance)
            !
            !--------------------------------------
            !distance = (xcCoord(i) - xctr)**2 
            !!distance = distance + (zcCoord(k) - zctr)**2
            !distance = sqrt(distance)
            !!!
            !if(distance < 0.5*border) then
            !  dens(i,j,k) = rho_l
            !  eint(i,j,k) = eint_l
            !else
            !  dens(i,j,k) = rho_r
            !  eint(i,j,k) = eint_r
            !endif
            !u(i,j,k) = 0.0
            !v(i,j,k) = 0.0
            !w(i,j,k) = 0.0
            !--------------------------------------
            if(xcCoord(i) < border) then
              dens(i,j,k) = rho_l
              u(i,j,k)    = vx_l
              eint(i,j,k) = eint_l
            else
              dens(i,j,k) = rho_r
              u(i,j,k)    = vx_r
              eint(i,j,k) = eint_r
            endif
            !
            v(i,j,k) = 0.0 
            w(i,j,k) = 0.0 
            ek = 0.5*u(i,j,k)**2
            ener(i,j,k) = eint(i,j,k) + ek
            !--------------------------------------
            !if(zcCoord(k) < border) then
            !  dens(i,j,k) = rho_l
            !  w(i,j,k)    = 0.1!vx_l
            !  eint(i,j,k) = eint_l
            !else
            !  dens(i,j,k) = rho_r
            !  w(i,j,k)    = 0.1!vx_r
            !  eint(i,j,k) = eint_r
            !endif
            !!
            !v(i,j,k) = 0.0 
            !u(i,j,k) = 0.0 
            !ek = 0.5*u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2
            !ener(i,j,k) = eint(i,j,k) + ek
            !
         enddo
      enddo
   enddo
   !
   write(*,*) '----- init_domain done ------------'
   !
end subroutine Simulation_init_advect
!
end module Simulation
