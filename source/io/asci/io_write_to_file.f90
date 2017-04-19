
subroutine io_write_to_file
   !
   use mo_grid, only: ib,ie,jb,je,kb,ke,xcCoord
   use mo_database
   use mo_io, only: outdata_unit
   !
   implicit none
   integer :: i,j,k
   !
   do i=ib,ie
    do j=jb,je
     do k=kb,ke
      !
      write(outdata_unit,'(5F18.8)') xcCoord(i),dens(i,j,k), &
                                     u(i,j,k),pres(i,j,k),eint(i,j,k)
      !
     enddo  
    enddo  
   enddo  
   !       
end subroutine io_write_to_file
!
