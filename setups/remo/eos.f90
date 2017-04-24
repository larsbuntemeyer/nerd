
subroutine eos 
   !
   use mo_constants
   use mo_database
   use mo_namelist
   use mo_grid, only: ib,ie,kb,ke,jb,je
   !
   implicit none
   integer :: i,j,k
   real, parameter :: R=287.05
   !
   do k=kb,ke
      do j=jb,je
         do i=ib,ie
            !
            pres(i,j,k) = eint(i,j,k)*dens(i,j,k)*(gamma-1.0)      
            !temp(i,j,k) = eint(i,j,k)*mu_mol*(gamma-1.0)/gas_constant
            temp(i,j,k) = pres(i,j,k)/(R*dens(i,j,k))
            !
         enddo  
      enddo  
   enddo
   !       
end subroutine eos 
