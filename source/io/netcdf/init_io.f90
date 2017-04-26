subroutine init_io
!
use netcdf
use mo_io
use mo_database
use mo_parameters, only: ndim,k2d,k3d
use mo_grid!, only: nx,ny,nz,ib,ie
!
implicit none
!

!
write(*,*) '----- Io_init ---------------------'
!
outdata_base = 'output.nc'

! Create the file. 
call check( nf90_create(outdata_base, nf90_clobber, ncid) )

! Define the dimensions. The record dimension is defined to have
! unlimited length - it can grow as needed. In this example it is
! the time dimension.
call check( nf90_def_dim(ncid, 'x', nx, x_dimid) )
if(k2d==1)  call check( nf90_def_dim(ncid, 'y', nx, y_dimid) )
if(k3d==1)  call check( nf90_def_dim(ncid, 'z', nz, z_dimid) )
call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, rec_dimid) )

! Define the coordinate variables. We will only define coordinate
! variables for lat and lon.  Ordinarily we would need to provide
! an array of dimension IDs for each variable's dimensions, but
! since coordinate variables only have one dimension, we can
! simply provide the address of that dimension ID (lat_dimid) and
! similarly for (lon_dimid).
call check( nf90_def_var(ncid, 'x', NF90_DOUBLE, x_dimid, x_varid) )
if(k2d==1) call check( nf90_def_var(ncid, 'y', NF90_DOUBLE, y_dimid, y_varid) )
if(k3d==1) call check( nf90_def_var(ncid, 'z', NF90_DOUBLE, z_dimid, z_varid) )
call check( nf90_def_var(ncid, 'time', NF90_DOUBLE, rec_dimid, rec_varid) )

! The dimids array is used to pass the dimids of the dimensions of
! the netCDF variables. Both of the netCDF variables we are creating
! share the same four dimensions. In Fortran, the unlimited
! dimension must come last on the list of dimids.
if(ndim==1) then
  dimids = (/ x_dimid, rec_dimid, -1, -1 /)
elseif(ndim==2 .and. k3d==0) then
  dimids = (/ x_dimid, y_dimid, rec_dimid, -1 /)
elseif(ndim==2 .and. k2d==0) then
  dimids = (/ x_dimid, z_dimid, rec_dimid, -1 /)
elseif(ndim==3) then
  dimids = (/ x_dimid, y_dimid, z_dimid, rec_dimid /)
endif
!
print*, 'defining variables'
!call check( nf90_def_var(ncid, 'dens', NF90_REAL, reshape(dimids,(/ndim+1/)), dens_id) )
call check( nf90_def_var(ncid, 'dens', NF90_REAL, dimids(1:ndim+1), dens_id) )
call check( nf90_def_var(ncid, 'eint', NF90_REAL, dimids(1:ndim+1), eint_id) )
call check( nf90_def_var(ncid, 'pres', NF90_REAL, dimids(1:ndim+1), pres_id) )
call check( nf90_def_var(ncid, 'u'   , NF90_REAL, dimids(1:ndim+1), u_id) )
if(k2d==1)  call check( nf90_def_var(ncid, 'v', NF90_REAL, dimids(1:ndim+1), v_id) )
if(k3d==1)  call check( nf90_def_var(ncid, 'w', NF90_REAL, dimids(1:ndim+1), w_id) )
call check( nf90_def_var(ncid, 'hydrop'   , NF90_REAL, dimids(1:ndim+1), hydrop_id) )
call check( nf90_def_var(ncid, 't'   , NF90_REAL, dimids(1:ndim+1), t_id) )
!
! End define mode.
call check( nf90_enddef(ncid) )
!
nt = 0
if(ndim==2 .and. k2d==0) then
  count = (/ nx, nz, 1 , 1/)
else
  count = (/ nx, ny, nz, 1 /)
endif

! Write the coordinate variable data. This will put the latitudes
! and longitudes of our data grid into the netCDF file.
call check( nf90_put_var(ncid, x_varid, xcCoord(ib:ie)) )
if(k2d==1)  call check( nf90_put_var(ncid, y_varid, ycCoord(jb:je)) )
if(k3d==1)  call check( nf90_put_var(ncid, z_varid, zcCoord(kb:ke)) )

call check( nf90_close(ncid) )
!
call io_write_to_file
!
write(*,*) 'output data file:',outdata_base
write(*,*) '----- Io_init done ----------------'
! 
end subroutine init_io
