!
module mo_io
!
use netcdf
use mo_parameters, only: ndim
!
implicit none
!
character(len=80),save   :: outdata_base
character(len=80),save   :: indata_base
integer          ,save   :: outdata_unit
!
integer :: ncid
integer :: x_dimid,y_dimid,z_dimid
integer :: rec_dimid
integer :: x_varid,y_varid,z_varid
!
integer :: dens_id, pres_id, eint_id 
integer :: hydrop_id, t_id, vx_id, vy_id, vz_id
integer :: rec_varid
integer :: u_id, v_id, w_id 
!
integer :: dimids(4)
!
integer :: nt
integer :: count(4)
!
integer, parameter :: nt_close = 10
!
!
contains
!
subroutine check(status)
integer, intent ( in) :: status

if(status /= nf90_noerr) then 
  print *, trim(nf90_strerror(status))
  stop "Stopped"
end if
end subroutine check  
   !
end module mo_io
!
