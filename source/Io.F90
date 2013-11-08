!
module Io
   !
   use Io_data
   use Grid 
   use Database
   !
   implicit none
   !
   integer :: i,j,k
   !
   interface 
      subroutine Io_init
         implicit none
      end subroutine Io_init
   end interface
   !
   contains 
   !
   subroutine Io_write_to_file
      !
      implicit none
      !
     ! do i=1+nguard,nx+nguard
     !    do j=1+k2d*nguard,ny+k2d*nguard
     !       do k=1+k3d*nguard,nz+k3d*nguard
     !          !
     !          write(outdata_unit,*) i,j,k,dens(i,j,k)
     !          !
     !       enddo  
     !    enddo  
     ! enddo  
      do i=1,nx+2*nguard
         do j=1,ny+k2d*2*nguard
            do k=1,nz+k3d*2*nguard
               !
               write(outdata_unit,*) i,j,k,dens(i,j,k),vx(i,j,k),pres(i,j,k)
               !
            enddo  
         enddo  
      enddo  
      !       
   end subroutine Io_write_to_file
   !
end module Io
!
