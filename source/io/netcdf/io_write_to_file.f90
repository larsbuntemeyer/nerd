
subroutine io_write_to_file
!
use mo_grid
use mo_database
use mo_io
use mo_parameters, only: k2d,k3d
use mo_driver
!
implicit none
integer :: i,j,k, odim
integer :: start(3)
!
if (mod(nt,nt_close)==0) call check( nf90_open(outdata_base, nf90_write, ncid) )

nt = nt+1
!
! Write the coordinate variable data. This will put the latitudes
! and longitudes of our data grid into the netCDF file.
call check( nf90_put_var(ncid, rec_varid, current_time, start=(/nt/)) )

! These settings tell netcdf to write one timestep of data. (The
! setting of start(4) inside the loop below tells netCDF which
! timestep to write.)
start = 1
odim  = ndim
!if(ndim==2 .and. k2d==0) odim=3

! Write the pretend data. This will write our surface pressure and
! surface temperature data. The arrays only hold one timestep worth
! of data. We will just rewrite the same data for each timestep. In
! a real :: application, the data would change between timesteps.
!
print*, 'count: ', count
call check( nf90_put_var(ncid, dens_id, dens(ib:ie,jb:je,kb:ke), start = (/start(1:odim),nt/), &
                            count = count) )
call check( nf90_put_var(ncid, pres_id, pres(ib:ie,jb:je,kb:ke), start = (/start(1:odim),nt/), &
                            count = count) )
call check( nf90_put_var(ncid, eint_id, eint(ib:ie,jb:je,kb:ke), start = (/start(1:odim),nt/), &
                            count = count) )
call check( nf90_put_var(ncid, u_id, u(ib:ie,jb:je,kb:ke), start = (/start(1:odim),nt/), &
                            count = count) )
if(k2d==1) call check( nf90_put_var(ncid, v_id, v(ib:ie,jb:je,kb:ke), start = (/start(1:odim),nt/), &
                            count = count) )
if(k3d==1) call check( nf90_put_var(ncid, w_id, w(ib:ie,jb:je,kb:ke), start = (/start(1:odim),nt/), &
                            count = count) )
call check( nf90_put_var(ncid, hydrop_id, hydrop(ib:ie,jb:je,kb:ke,1), start = (/1,1,1,nt/), &
                            count = count) )
call check( nf90_put_var(ncid, t_id, t(ib:ie,jb:je,kb:ke,1), start = (/1,1,1,nt/), &
                            count = count) )

if (mod(nt,nt_close)==0) call check( nf90_close(ncid) )
!       
end subroutine io_write_to_file
!
