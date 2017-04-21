!!****h* /source/database/amr/paramesh2.0/physicaldata
!!
!! NAME
!!
!!    physicaldata
!!
!!
!! SYNOPSIS
!!
!!   include file with physical data 
!!   
!!
!! DESCRIPTION
!!
!!   This file has physical data needed for the code.
!!   Description of data is included in the file close to 
!!   the declarations
!!  
!!***


module mo_parameters 

!
!
!
!Pre-Processor Control
!
!#ifndef CRAY
!#define shmem_real_get shmem_get8
!#define SHMEM_REAL_GET SHMEM_GET8
!#define shmem_real_put shmem_put8
!#define SHMEM_REAL_PUT SHMEM_PUT8
!#endif
!
! set pre-processor variable to control use of different timesteps on
! different grid blocks
!#define VAR_DT                                                 !<<< USER EDIT
!
! Does the algorithm use predictor-corrector type timestepping?
!#define PRED_CORR                                              !<<< USER EDIT
!
! Will any grid blocks represent obstacles to flow?
!#define EMPTY_CELLS                                             !<<< USER EDIT
!
! set the model dimension here. If running 2.5D also edit l2p5d below.
! Skip this.  We will set N_DIM directly from compilation command line.
!#define N_DIM 2
!
  integer k1d
  parameter(k1d=1)
!
!-----------------------------------------------------------------
! physicaldata.fh

!!inserted for compatability with paramesh3.0..  i want this always
!!defined for the mesh
      logical, save :: NPGFlag = .false.

!
! set physical dimension of model and number of edges on each grid block
  integer ndim,nbedges
#if N_DIM == 3
  parameter(ndim=3)
#endif
#if N_DIM == 2
  parameter(ndim=2)
#endif
#if N_DIM == 1
  parameter(ndim=1)
#endif
!
! l2p5d needs to be declared not matter what the dimension of the
! problem is set to !
  integer l2p5d
  parameter(l2p5d=0)                                        !<<< USER EDIT
                                                                ! 1 if 2.5D
                                                                ! 0 if 2 or 3D
!
  parameter(nbedges=ndim*2**(ndim-1))
!
!
!
! an increment variable for the z dimension to enable the same code to
! work for 2D or 3D models.
  integer k3d
  parameter(k3d=(ndim-1)/2)
!
  integer k2d
  parameter(k2d=ndim/2)
!
!
! set size of grid blocks
! The spatial relationship between grid points in parents
! and their children is subtly different depending on whether
! the block dimensions (ie nxb,nyb,nzb) are even or odd. This has 
! significant consequences when defining interpolation within
! restriction or prolongation operations.
!
  integer nxb,nyb,nzb,maxdim

#ifdef NXB
  parameter(nxb = NXB)
#else
  parameter(nxb=8)                                        !<<< USER EDIT  
#endif

#ifdef NYB
  parameter(nyb = NYB)
#else
  #if N_DIM > 1
     parameter(nyb=8)                                    !<<< USER EDIT
  #else
    parameter(nyb=1)
  #endif
#endif

#ifdef NZB
  parameter(nzb = NZB)
#else
  #if N_DIM > 2   
    parameter(nzb=8)                                     !<<< USER EDIT
  #else
    parameter(nzb=1)
  #endif
#endif

!

!
!
!
  parameter(maxdim=max(nxb,nyb,nzb))
!
!
! these guard cell offsets are required to accomodate differences
! in cases when block dimensions are odd or even
  integer gc_off_x,gc_off_y,gc_off_z
  parameter(gc_off_x=mod(nxb,2))
  parameter(gc_off_y=mod(nyb,2))
  parameter(gc_off_z=mod(nzb,2))
!
! set the maximum number of blocks per processor
  integer maxblocks
#ifdef MAXBLOCKS
  parameter (maxblocks = MAXBLOCKS)
#else
#if N_DIM == 3
  parameter (maxblocks = 200)              !<<< USER EDIT
#else /* N_DIM < 3 */
  parameter (maxblocks = 1000)             !<<< USER EDIT 
#endif /*N_DIM*/
#endif /*MAXBLOCKS*/
!
!
!
!
!..this include file is very important at this location, as it sets a 
!..parameter (ionmax) that determines the array sizes and do-loop limits 
!..of the mesh, hydro, eos and burn modules. it touches just about everything.
!

#include "nerd_defines.fh"

!
! set number of unknowns associated with each grid cell
  integer nvar
  !parameter(nvar=NERD_NUMBER_OF_VARIABLES)
  parameter(nvar=5)
  integer nvar2
  parameter(nvar2=5)
  integer nvarsm
  parameter (nvarsm=2)
!
! set the number of guard cell layers at each boundary
  integer nguard
  parameter(nguard=NERD_NUMBER_OF_GUARD_CELLS)
!
!
! common block storing the solution for cell-centered quantities.
! unksm stores copies of global variables which DO NOT need guard cells
! AND do not need to be saved from one timestep to the next !!!
!
  integer il_bnd,iu_bnd
  integer jl_bnd,ju_bnd
  integer kl_bnd,ku_bnd
  parameter(il_bnd=1, iu_bnd=nxb+2*nguard)
  parameter(jl_bnd=1, ju_bnd=nyb+2*nguard*k2d)
  parameter(kl_bnd=1, ku_bnd=nzb+2*nguard*k3d)
  integer nxlo,nylo,nzlo,nxhi,nyhi,nzhi
  parameter (nxlo=nguard+1,nylo=nguard*k2d+1,nzlo=nguard*k3d+1)
  parameter (nxhi=nguard+nxb)
  parameter (nyhi=nguard*k2d+nyb)
  parameter (nzhi=nguard*k3d+nzb)
  
!  real, DIMENSION(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,                   &
!       &     kl_bnd:ku_bnd,maxblocks), TARGET :: unk
!  real, DIMENSION(nvar2,il_bnd:iu_bnd,maxblocks),                     &
!       &     TARGET :: unk2
!  real, DIMENSION(nvar2,jl_bnd:ju_bnd,maxblocks),                     &
!       &     TARGET :: unk3
!  real, DIMENSION(nvar2,kl_bnd:ku_bnd,maxblocks),                     &
!       &     TARGET :: unk4
!  real, DIMENSION(nvarsm,nxlo:nxhi,                                   &
!       &     nylo:nyhi,nzlo:nzhi,maxblocks), TARGET :: unksm
!  
!  
! 
!!
!!
!!
!! The convention for relating variables associated with cell faces to the
!! variables defined at cell centers is as follows:
!!
!! If iface_off=0 :
!!         the array facevarx(:,i,j,k,:) for example defines data
!!         on the x(i-1/2) face of the (i,j,k)-th mesh cell.
!! If iface_off=-1 :
!!         the array facevarx(:,i,j,k,:) for example defines data
!!         on the x(i+1/2) face of the (i,j,k)-th mesh cell.
!!
!  integer iface_off
!  parameter(iface_off=0)                                  !<<< USER EDIT
!!
!!
!! The number of data words needed on a cell face is set by nfacevar.
!!
!!
!  integer nfacevar
!! 2 added to store strong fields at faces for all components of B
!  parameter(nfacevar=0)   !<<< USER EDIT
!!
!  integer nbndvar
!  parameter(nbndvar=max(1,nfacevar))
!!
!  integer maxblocksf
!  parameter(maxblocksf= 1+(maxblocks-1)*min(1,nfacevar) )
!!
!! common block storing the solution for cell-face-centered quantities.
!  real ::                                                           &
!       &     facevarx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,              &
!       &     kl_bnd:ku_bnd,maxblocksf)                                    &
!       &     ,facevary(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,           &
!       &     kl_bnd:ku_bnd,maxblocksf)                                    &
!       &     ,facevarz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,               &
!       &     kl_bnd:ku_bnd+k3d,maxblocksf)
!  
!!
!! set data length of grid blocks
!  integer len_block
!  parameter(len_block=iu_bnd*ju_bnd*ku_bnd*nvar)
!!
!  integer len_blockfx,len_blockfy,len_blockfz
!  parameter(len_blockfx=(nxb+2*nguard+1)*(nyb+2*nguard*k2d)*        &
!       &                          (nzb+2*nguard*k3d))
!  parameter(len_blockfy=(nxb+2*nguard)*(nyb+2*nguard+1)*            &
!       &                          (nzb+2*nguard*k3d))
!  parameter(len_blockfz=(nxb+2*nguard)*(nyb+2*nguard)*              &
!       &                          ((nzb+2*nguard)*k3d+1))
!!
!!
!! common block for timestep control
!  integer maxlevels
!  parameter(maxlevels=20)
!  common/timecntl/                                                  &
!       & time_loc(maxblocks),dtlevel(maxlevels),ldtcomplete(maxblocks)
!  real time_loc,dtlevel
!  logical ldtcomplete
!!
!!
!#if defined(VAR_DT) || defined(PRED_CORR)
!!      parameter(maxblocks_t=(maxblocks-1)*ivar_dt+1)
!!      parameter(nvar_t=(nvar-1)*ivar_dt+1)
!  common/tsolution/                                                 &
!       &      t_unk(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,       &
!       &       maxblocks)
!  real t_unk
!#endif
!!
!#ifdef PRED_CORR
!  real ::                                                         &
!       &           tfacevarx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,       &
!       &                          kl_bnd:ku_bnd,maxblocksf)               &
!       &          ,tfacevary(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,     &
!       &                          kl_bnd:ku_bnd,maxblocksf)               &
!       &          ,tfacevarz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,         &
!       &                          kl_bnd:ku_bnd+k3d,maxblocksf)
!  
!#endif
!!
!! To average fluxes set red_f = .25.
!! To sum fluxes set red_f = 1.0.
!!
!! changed -- 2-24-00 
!! we are now converting the flux densities to fluxes before the call to
!! amr_flux_conserve.  This is to get the proper geometry factors included
!! in non-cartesian cases.  The fluxes are converted back after the call.
!! 
!! red_f is set according to dimension, to sum instead of average
!  real red_f
!!      parameter(red_f = 0.25)   
!!
!#if N_DIM == 1
!  parameter (red_f = 0.25)
!#elif N_DIM == 2
!  parameter (red_f = 0.5)
!#elif N_DIM == 3
!  parameter (red_f = 1.0)
!#endif
!!
!!-----------------------------------------------------------------
!! include header file defining data structure on cell faces
!!
!!-----------------------------------------------------------------
!!
!! This file defines a data structure to be used for quantities
!! which may need to be defined at grid block interfaces, eg fluxes,
!! pressures.
!!
!!
!!
!! storage used for fluxes at block boundaries. This is used when conservation
!! constraints need to be imposed.
!!
!! updated 2-15-00 -- allocate storage for the internal energy flux

#include "nerd_defines.fh"
  
  integer nfluxvar
  parameter(nfluxvar=NERD_NUMBER_OF_FLUXES)
!
  integer nfluxes
  parameter(nfluxes=max(1,nfluxvar))
!
  integer maxblocksfl
  parameter(maxblocksfl= 1+(maxblocks-1)*min(1,nfluxvar) )
!!
!!
!!
!!..in 1d the flux_y, flux_z, tflux_y, and tflux_z arrays are not used,
!!..but do need to be declared. thus, in 1d the parameter maxblocksfl
!!..has been replaced with a 1. this significantly reduces the
!!..memory footprint for 1d simulations.
!!
!!..in 2d the flux_z and tflux_z arrays are not used,
!!..but do need to be declared. thus, in 2d the parameter maxblocksfl
!!..has been replaced with a 1. this significantly reduces the
!!..memory footprint for 2d simulations.
!!
!#if N_DIM == 1
!  real, target ::                                                  &
!       flux_x(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocksfl),&
!       flux_y(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd,1),          &
!       flux_z(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2,1)
!!
!  real, target ::                                                   &
!       tflux_x(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocksfl),&
!       tflux_y(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd,1),          &
!       tflux_z(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2,1)
!
!#endif
!!
!!
!#if N_DIM == 2
!  real, target ::                                                   &
!       flux_x(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocksfl), &
!       flux_y(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd,maxblocksfl), &
!       flux_z(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2,1)
!  
!!
!  real, target ::                                                    &
!       tflux_x(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocksfl), &
!       tflux_y(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd,maxblocksfl), &
!       tflux_z(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2,1)
!#endif
!!
!!
!#if N_DIM == 3
!  real, target ::                                                  &
!       flux_x(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocksfl),&
!       flux_y(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd,maxblocksfl),&
!       flux_z(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2,maxblocksfl)
!!
!  real, target ::                                                    &
!       tflux_x(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocksfl), &
!       tflux_y(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd,maxblocksfl), &
!       tflux_z(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2,maxblocksfl)
!#endif
!!
!!
!!
!! storage used for cell edges at block boundaries. 
!! This is used when quantities located at cell edge centers need to
!! be used consistently at the boundaries between blocks at different
!! refinement levels.
!!
!  integer nedgevar
!  parameter(nedgevar=1)                                     !<<< USER EDIT
!!
!  integer nedges
!  parameter(nedges=max(1,nedgevar))
!!
!  integer maxblockse
!  parameter(maxblockse= 1+(maxblocks-1)*min(1,nedgevar) )
!!
!!
!!..flash does not presently use these edge storage variables
!!..but it might when magnetic fields are included. until then,
!!..the maxblockse size declaration has been replaced with 1 in order
!!..to reduce the memory footprint.
!!
!!      common/edges/
!!     .  bedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1,
!!     .    kl_bnd:ku_bnd+1,maxblockse),
!!     .  bedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1,
!!     .    kl_bnd:ku_bnd+1,maxblockse),
!!     .  bedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2,
!!     .    kl_bnd:ku_bnd+1,maxblockse),
!!     .  bedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2,
!!     .    kl_bnd:ku_bnd+1,maxblockse),
!!     .  bedge_facez_x(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,
!!     .    1:2,maxblockse),
!!     .  bedge_facez_y(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,
!!     .    1:2,maxblockse),
!!
!!
!  real ::                                                         &
!       &  bedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1,                       &
!       &    kl_bnd:ku_bnd+1,1),                                           &
!       &  bedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1,                       &
!       &    kl_bnd:ku_bnd+1,1),                                           &
!       &  bedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2,                       &
!       &    kl_bnd:ku_bnd+1,1),                                           &
!       &  bedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2,                       &
!       &    kl_bnd:ku_bnd+1,1),                                           &
!       &  bedge_facez_x(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,           &
!       &    1:2,1),                                                       &
!       &  bedge_facez_y(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,           &
!       &    1:2,1),                                                       &
!       &  recvarx1e(nedges,1:2,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd+1),          &
!       &  recvary1e(nedges,il_bnd:iu_bnd+1,1:2,kl_bnd:ku_bnd+1),          &
!       &  recvarz1e(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,1:2),          &
!       &  recvarx2e(nedges,1:2,jl_bnd:ju_bnd+1,kl_bnd:ku_bnd+1),          &
!       &  recvary2e(nedges,il_bnd:iu_bnd+1,1:2,kl_bnd:ku_bnd+1),          &
!       &  recvarz2e(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,1:2)
!!
!!
!!..flash does not presently use these edge storage variables
!!..but it might when mhd is added in.
!!..the maxblockse size declaration has been reduced to 1 in order
!!..to help reduce the memory footprint.
!!..28nov99 fxt
!!
!!
!!      common/tedges/
!!     .  tbedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1,
!!     .    kl_bnd:ku_bnd+1,maxblockse),
!!     .  tbedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1,
!!     .    kl_bnd:ku_bnd+1,maxblockse),
!!     .  tbedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2,
!!     .    kl_bnd:ku_bnd+1,maxblockse),
!!     .  tbedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2,
!!     .    kl_bnd:ku_bnd+1,maxblockse),
!!     .  tbedge_facez_x(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,
!!     .    1:2,maxblockse),
!!     .  tbedge_facez_y(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,
!!     .    1:2,maxblockse)
!!
!  real ::                                                         &
!       &  tbedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1,                      &
!       &    kl_bnd:ku_bnd+1,1),                                           &
!       &  tbedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1,                      &
!       &    kl_bnd:ku_bnd+1,1),                                           &
!       &  tbedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2,                      &
!       &    kl_bnd:ku_bnd+1,1),                                           &
!       &  tbedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2,                      &
!       &    kl_bnd:ku_bnd+1,1),                                           &
!       &  tbedge_facez_x(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,          &
!       &    1:2,1),                                                       &
!       &  tbedge_facez_y(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,          &
!       &    1:2,1)
!!
!!
!!
!! workspace arrays used for inter-block communications
!  integer nbndmax
!  parameter(nbndmax=max(nbndvar,nfluxes))
!  real ::                                                         &
!       &     recvarx1(nbndmax,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd)            &
!       &    ,recvary1(nbndmax,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd)            &
!       &    ,recvarz1(nbndmax,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2)            &
!       &    ,bndtempx1(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd)           &
!       &    ,bndtempy1(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd)           &
!       &    ,bndtempz1(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2)
!  
!!
!!
!!
!!
!! parameters used in communication calls
!  integer len_block_bndx,len_block_bndy,len_block_bndz
!  parameter(len_block_bndx=2*(nyb+2*nguard*k2d)*                    &
!       &                                 (nzb+2*nguard*k3d))
!  parameter(len_block_bndy=2*(nxb+2*nguard*k2d)*                    &
!       &                                 (nzb+2*nguard*k3d))
!  parameter(len_block_bndz=2*(nxb+2*nguard)*(nyb+2*nguard))
!!
!  integer len_block_ex,len_block_ey,len_block_ez
!  parameter(len_block_ex=2*(nyb+k2d+2*nguard*k2d)*                  &
!       &                                 (nzb+k3d+2*nguard*k3d))
!  parameter(len_block_ey=2*(nxb+1+2*nguard)*                        &
!       &                                 (nzb+k3d+2*nguard*k3d))
!  parameter(len_block_ez=2*(nxb+1+2*nguard)*                        &
!       &                                 (nyb+k2d+2*nguard))
!  

end module mo_parameters 






