!
!
!
module mo_grid
   !
   implicit none
   !
   integer,          save :: ndim
   integer,          save :: nx,ny,nz
   integer,          save :: k2d,k3d 
   integer,          save :: nguard
   integer,          save :: ib,jb,kb
   integer,          save :: ie,je,ke
   integer,          save :: ibg,jbg,kbg
   integer,          save :: ieg,jeg,keg
   real, save :: dx,dy,dz
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
!
!
end module mo_grid
