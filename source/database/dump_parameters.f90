subroutine dump_parameters

use mo_parameters

implicit none


write(*,*) ''
write(*,*) 'parameters from setup:'
write(*,*) ''
write(*,100) 'ndim', ndim
write(*,*) ''
write(*,100) 'nxb', nxb
write(*,100) 'nyb', nyb
write(*,100) 'nzb', nzb
write(*,*) ''
write(*,100) 'nguard', nguard
write(*,*) ''
write(*,100) 'k2d', k2d
write(*,100) 'k3d', k3d
write(*,*) ''
write(*,101) 'il_bnd, iu_bnd', il_bnd,iu_bnd
write(*,101) 'jl_bnd, ju_bnd', jl_bnd,ju_bnd
write(*,101) 'kl_bnd, ku_bnd', kl_bnd,ku_bnd
write(*,*) ''
write(*,*) ''


100 format (A10,I5)
101 format (A21,2I5)

end subroutine dump_parameters
