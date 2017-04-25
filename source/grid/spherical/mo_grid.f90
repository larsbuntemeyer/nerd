!
!
!
module mo_grid
   !
   !use mo_parameters, only: nguard, ndim, k2d, k3d
   !
   implicit none
   !
   !integer,          save :: ndim
   integer,          save :: nx,ny,nz
   integer,          save :: nz1,nxny
   !integer,          save :: k2d,k3d 
   !integer,          save :: nguard
   integer,          save :: ib,jb,kb
   integer,          save :: ie,je,ke
   integer,          save :: ibg,jbg,kbg
   integer,          save :: ieg,jeg,keg
   real, save :: dx,dy,dz
   REAL    ::          EDDLAM , EDDPHI , EDADPHI, DLADDPH, DPHDDLA
   real    :: dlam,dphi
   !
   !  The cell coordinates
   !  c: center coordinate of the cell
   !  l: left cell-face coordinate
   !  r: right cell-face coordinate
   !
   real, dimension(:),      allocatable :: xcCoord
   real, dimension(:),      allocatable :: ycCoord
   real, dimension(:),      allocatable :: zcCoord
   real, dimension(:),      allocatable :: xlCoord
   real, dimension(:),      allocatable :: ylCoord
   real, dimension(:),      allocatable :: zlCoord
   real, dimension(:),      allocatable :: xrCoord
   real, dimension(:),      allocatable :: yrCoord
   real, dimension(:),      allocatable :: zrCoord
   !
   REAL, ALLOCATABLE, DIMENSION(:,:) :: GCPHI, GACPHIR, ACPHIR, CPHI  
   !
   REAL, ALLOCATABLE, DIMENSION(:)   :: AK, BK, AKH, BKH, DAK, DBK, A1T, A2T, VVFH
   !
   real :: ptop
   !
   !
end module mo_grid
