!
module Grid_data
   !
   implicit none
   !
   integer,          save :: ndim
   integer,          save :: nx,ny,nz
   integer,          save :: k2d,k3d 
   integer,          save :: nguard
   double precision, save :: dx,dy,dz 
   !
   !  The cell coordinates
   !  c: center coordinate of the cell
   !  l: left cell-face coordinate
   !  r: right cell-face coordinate
   !
   double precision, dimension(:),      allocatable :: xcCoord
   double precision, dimension(:),      allocatable :: ycCoord
   double precision, dimension(:),      allocatable :: zcCoord
   double precision, dimension(:),      allocatable :: xlCoord
   double precision, dimension(:),      allocatable :: ylCoord
   double precision, dimension(:),      allocatable :: zlCoord
   double precision, dimension(:),      allocatable :: xrCoord
   double precision, dimension(:),      allocatable :: yrCoord
   double precision, dimension(:),      allocatable :: zrCoord
   !
end module Grid_data
!
