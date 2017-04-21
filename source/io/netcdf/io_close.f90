subroutine io_close

use mo_io, only: ncid, check
use netcdf

implicit none
!
! Close the file. This causes netCDF to flush all buffers and make
! sure your data are really written to disk.
call check( nf90_close(ncid) )

end subroutine io_close
