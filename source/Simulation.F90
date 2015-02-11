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
   real    :: rho_r,rho_l,vx_r,vx_l
   real    :: ek,ei,e
   !
   rho_r = rho_right
   rho_l = rho_left
   vx_l  = vx_left
   vx_r  = vx_right
   rho_r = 1.d-1
   rho_l = 1.d0
   vx_l  = 1.d-1
   vx_r  = 1.d-1
   !
   xctr  = 0.5d0*(xmax-xmin)
   yctr  = 0.5d0*(ymax-ymin)
   zctr  = 0.5d0*(zmax-zmin)
   xsize = (xmax-xmin)
   ysize = (ymax-ymin)
   zsize = (zmax-zmin)
   !
   write(*,*) '----- init_domain -----------------'
   !
   do i=1+nguard,nx+nguard
      do j=1+k2d*nguard,ny+k2d*nguard
         do k=1+k3d*nguard,nz+k3d*nguard
            !
            ! Here one should init the physics 
            !
            distance = (xcCoord(i) - xctr)**2 
            if(ndim>2) then 
               distance = distance + (ycCoord(j) - yctr)**2 
            endif
            if(ndim==3) then 
               distance = distance + (zcCoord(k) - zctr)**2 
            endif
            !
            !distance = sqrt(distance)
            !
            if(xcCoord(i) < xctr) then
              dens(i,j,k) = rho_l
              u(i,j,k) = vx_l
            else
              dens(i,j,k) = rho_r
              u(i,j,k) = vx_r
            endif
            !
            ek = 0.5d0*u(i,j,k)**2
            e = 1.0 + ek
            ener(i,j,k) = e
            eint(i,j,k) = e - ek
            !
         enddo
      enddo
   enddo
   !
   write(*,*) '----- init_domain done ------------'
   !
end subroutine Simulation_init_domain
!
end module Simulation
