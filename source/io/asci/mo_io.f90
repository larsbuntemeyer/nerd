!
module mo_io
   !
   !use grid 
   !use database
   !!
   !implicit none
   !!
   !public :: Io_init, Io_write_to_file
   !!
   character(len=80),save   :: outdata_base
   character(len=80),save   :: indata_base
   integer          ,save   :: outdata_unit
   !!
   !!
   !contains 
   !!
   !!
   !subroutine Io_write_to_file
   !   !
   !   implicit none
   !   integer :: i,j,k
   !   !
   !   do i=ibg,ieg
   !    do j=jb,je
   !     do k=kbg,keg
   !      !
   !      write(outdata_unit,'(5F18.8)') xcCoord(i),dens(i,j,k), &
   !                                     u(i,j,k),pres(i,j,k),eint(i,j,k)
   !      !
   !     enddo  
   !    enddo  
   !   enddo  
   !   !       
   !end subroutine Io_write_to_file
   !
end module mo_io
!
