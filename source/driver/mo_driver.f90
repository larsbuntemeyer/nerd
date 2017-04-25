module mo_driver

implicit none

integer   :: step
real      :: current_time
real      :: dt,dt2,ed2dt
real      :: dtdeh
integer   :: nold, nnow, nnew
integer   :: nold2, nnow2
real      :: vbmxv,vbcfl
real, parameter   :: alcnva=0.5

end module mo_driver
