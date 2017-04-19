!
!
!
subroutine sweep1D(dt,dx,dir,rho,u,v,w,eint,pres,n,nvar)
   !
   !
   implicit none
   !
   real,    intent(in) :: dx
   real,    intent(in) :: dt
   integer, intent(in) :: dir,n,nvar
   real, dimension(n), intent(in)    :: pres 
   real, dimension(n), intent(inout) :: rho,u,v,w,eint
   !
   !
end subroutine sweep1D
