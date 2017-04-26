subroutine init_domain
   !
   use mo_grid
   use mo_namelist
   use mo_database, only: ps,t,u=>vx,w=>vz,fib, fi, qw, qd, qi
   use mo_constants
   use mo_parameters, only: ndim
   !
   implicit none
   !
   integer             :: i,j,k 
   real    :: distance
   real    :: xctr,yctr,zctr
   real    :: xsize,ysize,zsize
   real    :: rho_in,rho_out,eint_in,eint_out,pres_in,pres_out
   real    :: ek,ei,e, radius, disx, bfac
   integer,parameter :: seed = 86456
   !
   real, parameter :: BCZ=3000.0,BCI=150.0,BRZ=2000.0,BRX=4000.0,BVAL=6.6
   !
   real    :: rho_left,rho_right, eint_env, zfih
   real    :: vx_left,vx_right, p,pt, dp, rho_t
   !
   !
   !eint_env = 300.0/(mu_mol*(gamma-1.0)) * gas_constant
   eint_env = R*300.0/(gamma-1.0) 
   rho_t = 0.1
   pt = 52000.0
   rho_in  = 1.0
   pres_in = 0.1
   rho_out  = 0.1
   pres_out = 0.1 
   radius = 1500.
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
   QW(:,:,:,1) = 0.
   QD(:,:,:,1) = 0.
   QI(:,:,:,1) = 0.
   U(:,:,:,1) = 0.
   W(:,:,:,1) = 0.
   T(:,:,:,1) = 288.15
   PS(:,:,1) = 101325.
   FIB = 0.
   FI = 0.
!  now add straka-like bubble
!  distance in m per degree longitude
!   DISX = RERD*PI/180.*COS(ycCoord(1)) 
!   do k=kb,ke-1
!      do j=je,jb
!         do i=ib,ie
!            !
!            zfih = 0.5*(fi(i,j,k,1)+fi(i,j,k+1,1))
!            BFAC = SQRT(((I-BCI)*DLAM*DISX/BRX)**2 + ((ZFIH/G - BCZ)/BRZ)**2)
!            IF (BFAC <= 1.) THEN
!              T(I,J,K,1) = T(I,J,K,1) + BVAL*(COS(PI*BFAC) + 1.)/2.
!            END IF
!         enddo
!      enddo
!   enddo
   !do k=kb,ke
   !   do j=je-1,jb,-1
   !      do i=ib,ie
   !         !
   !         distance = (xcCoord(i) - xctr)**2 
   !         if(ndim>1) then 
   !            distance = distance + (ycCoord(j) - yctr)**2 
   !         endif
   !         if(ndim==3) then 
   !            distance = distance + (zcCoord(k) - zctr)**2 
   !         endif
   !         !
   !         distance = sqrt(distance)
   !         !
   !         if(distance < radius) then
   !           dens(i,j,k) = dens(i,j,k)*1.1
   !           !pres(i,j,k) = dens(i,j,k)*(R*300.0)
   !           eint(i,j,k) = pres(i,j,k)/(dens(i,j,k)*(gamma-1.0))
   !         endif
   !         u(i,j,k) = 0.
   !         v(i,j,k) = 0.
   !         w(i,j,k) = 0.
   !         ek = 0.5*(u(i,j,k)**2+v(i,j,k)**2+w(i,j,k)**2)
   !         ener(i,j,k) = eint(i,j,k) + ek
   !         !
   !      enddo
   !   enddo
   !enddo
   !
   write(*,*) '----- init_domain done ------------'
   !
end subroutine init_domain
