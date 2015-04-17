!
!
!
module Eos  
   !
   use Grid 
   use Database
   use RuntimeParameters
   use PhysicalConstants 
   !
   implicit none
   !
   public :: Eos_gamma
   !
   integer :: i,j,k
   !
   contains 
   !
   !
   !
   subroutine Eos_gamma
      !
      implicit none
      !
      do k=kb,ke
         do j=jb,je
            do i=ib,ie
               !
               pres(i,j,k) = eint(i,j,k)*dens(i,j,k)*(gamma-1.d0)      
               temp(i,j,k) = eint(i,j,k)*mu_mol*(gamma-1.d0)/gas_constant      
               !
            enddo  
         enddo  
      enddo  
      !       
   end subroutine Eos_gamma 
   !
end module Eos
!
