!!****h* source/database/amr/dBaseDeclarations
!!
!! NAME
!!  dBaseDeclarations
!!
!! SYNOPSIS
!!  Database data that does not belong to paramesh.  
!!
!! SEE ALSO
!!  dBaseInclude.F90 for data that is shared with paramesh
!!
!! DESCRIPTION
!!  Amongst all the data here, there is also:
!!
!!  dbase_defines.fh :  created by the setup script
!!    These are problem specific parameters, specified at compile time.
!!    All preprocessor definitions:
!!      DBASE_VAR_NAMES_*, DBASE_FLUX_NAMES_*, 
!!      DBASE_VAR_ADVECT,  DBASE_VAR_RENORM, DBASE_VAR_CONSERVE,
!!      DBASE_NUMBER_OF_ADVECT_VARIABLES,
!!      DBASE_NUMBER_OF_RENORM_VARIABLES,
!!      DBASE_NUMBER_OF_CONSERVE_VARIABLES,
!!      DBASE_NUMBER_OF_SPECIES
!!      DBASE_NUMBER_OF_MASS_SCALARS
!!
!! SOURCE
!!  dbase_defines.fh
!!  
!! ATTRIBUTES
!!   ionmax, varnam,
!!   dbFluxNames, dbNamedVarAdvect, dbNamedVarRenorm,
!!   dbNamedVarConserve, dbNumVarAdvect, dbNumVarRenorm,  
!!   dbNumVarConserve, dbVarAdvect, dbVarRenorm, dbVarConserve,  
!!   irhoflx, iuflx, ipflx, iutflx, iuttflx,  
!!   ieflx, ieintflx, inuc_begin, inucflx_begin,  
!!   iref1, iref2, iref3, iref4,
!!   MAXCELLS, 
!!   dbX1, dbY1, dbZ1, dbX2, dbY2, dbZ2,  
!!   dbXn, dbYn, dbZn,  
!!   dbNumVariables, dbNumSpecies, dbNumFluxes,  
!!   dbNumVarNamed, dbNumFluxNamed, 
!!   dbSpeciesBegin, dbSpeciesEnd,  
!!   dbMaxBlocks, dbNumZones,  
!!   dbiPoint,  
!!   dbixVector, dbiyVector, dbizVector,  
!!   dbixyPlane, dbixzPlane, dbiyzPlane,  
!!   dbiyxPlane, dbizxPlane, dbizyPlane,  
!!   dbixyzCube, dbiAllVariables, dbiAllSpecies,  
!!   dbixCoord, dbiyCoord, dbizCoord,  
!!   dbiznl, dbizn, dbiznr,
!!   dbiugrid,
!!   ixznl, ixzn, ixznr, 
!!   iugridx, iyznl, iyzn, iyznr, 
!!   iugridy, izznl, izzn, izznr, 
!!   iugridz,  
!!   bnd_reflect, bnd_outflow, bnd_periodic,  
!!   bnd_user, bnd_hydrostat,
!!   CARTESIAN, CYLINDRICAL, SPHERICAL, POLAR, 
!!   geom_cartesian, geom_planar, geom_cylrad, geom_sphrad,  
!!   geom_cylang, geom_sphtheta, geom_sphphi,  
!!   sweep_all, sweep_x, sweep_y, sweep_z,  
!!   sweep_order_xyz, sweep_order_zyx,  
!!   itemp_old, ishockt, 
!!   time, dt, dtold,
!!   redshift, oldredshift,
!!   cpuseconds,
!!   nstep, nbegin, MyPE, MasterPE, NumPEs
!!***   

module dBaseDeclarations

  use mo_parameters 

  implicit none

!---------------------------------------------------------------------
!
! Global variables
!
!---------------------------------------------------------------------

  real time, dt, dtold, redshift, oldredshift, cpuseconds
  integer nstep, nbegin, MyPE, MasterPE, NumPEs

!--------------------------------------------------------------------
!
! Parameters set by setup script
!
!--------------------------------------------------------------------

#include "dbase_defines.fh"

  integer, parameter     :: dbNumVariables = nvar
  integer, parameter     :: dbNumFluxes    = nfluxes
  integer, parameter     :: dbNumSpecies   = DBASE_NUMBER_OF_SPECIES
  integer, parameter     :: dbNumMassScalars= DBASE_NUMBER_OF_MASS_SCALARS
  integer, parameter     :: ionmax         = DBASE_NUMBER_OF_SPECIES
  integer, parameter     :: nMassScalars   = DBASE_NUMBER_OF_MASS_SCALARS

  integer, parameter     :: dbNumVarNamed  =                               &
                             dbNumVariables - dbNumSpecies - dbNumMassScalars

  integer, parameter     :: dbNumFluxNamed =                               &
                             dbNumFluxes - dbNumSpecies - dbNumMassScalars

  character (len=4)      :: varnam(dbNumVarNamed)

  character (len=16)     :: dbFluxNames(max(dbNumFluxNamed,1))

  integer                :: dbNamedVarAdvect(dbNumVarNamed)
  integer                :: dbNamedVarRenorm(dbNumVarNamed)
  integer                :: dbNamedVarConserve(dbNumVarNamed)

!  data varnam & 
!    DBASE_VAR_NAMES_1
!    DBASE_VAR_NAMES_2
!    DBASE_VAR_NAMES_3
!    DBASE_VAR_NAMES_4
!    DBASE_VAR_NAMES_5
!    DBASE_VAR_NAMES_6
!    DBASE_VAR_NAMES_7
!    DBASE_VAR_NAMES_8
!    DBASE_VAR_NAMES_9
!    DBASE_VAR_NAMES_10
!    DBASE_VAR_NAMES_11
!    DBASE_VAR_NAMES_12
!
!#ifdef DBASE_FLUX_NAMES_1
!  data dbFluxNames & 
!    DBASE_FLUX_NAMES_1
!    DBASE_FLUX_NAMES_2
!    DBASE_FLUX_NAMES_3
!    DBASE_FLUX_NAMES_4
!    DBASE_FLUX_NAMES_5
!    DBASE_FLUX_NAMES_6
!    DBASE_FLUX_NAMES_7
!    DBASE_FLUX_NAMES_8
!    DBASE_FLUX_NAMES_9
!    DBASE_FLUX_NAMES_10
!    DBASE_FLUX_NAMES_11
!    DBASE_FLUX_NAMES_12
!#endif
!
!  data dbNamedVarAdvect &
!      DBASE_VAR_ADVECT_1
!      DBASE_VAR_ADVECT_2
!      DBASE_VAR_ADVECT_3
!      DBASE_VAR_ADVECT_4
!      DBASE_VAR_ADVECT_5
!      DBASE_VAR_ADVECT_6
!      DBASE_VAR_ADVECT_7
!      DBASE_VAR_ADVECT_8
!      DBASE_VAR_ADVECT_9
!      DBASE_VAR_ADVECT_10
!      DBASE_VAR_ADVECT_11
!      DBASE_VAR_ADVECT_12
!
!  data dbNamedVarRenorm &
!      DBASE_VAR_RENORM_1
!      DBASE_VAR_RENORM_2
!      DBASE_VAR_RENORM_3
!      DBASE_VAR_RENORM_4
!      DBASE_VAR_RENORM_5
!      DBASE_VAR_RENORM_6
!      DBASE_VAR_RENORM_7
!      DBASE_VAR_RENORM_8
!      DBASE_VAR_RENORM_9
!      DBASE_VAR_RENORM_10
!      DBASE_VAR_RENORM_11
!      DBASE_VAR_RENORM_12
!
!  data dbNamedVarConserve &
!      DBASE_VAR_CONSERVE_1
!      DBASE_VAR_CONSERVE_2
!      DBASE_VAR_CONSERVE_3
!      DBASE_VAR_CONSERVE_4
!      DBASE_VAR_CONSERVE_5
!      DBASE_VAR_CONSERVE_6
!      DBASE_VAR_CONSERVE_7
!      DBASE_VAR_CONSERVE_8
!      DBASE_VAR_CONSERVE_9
!      DBASE_VAR_CONSERVE_10
!      DBASE_VAR_CONSERVE_11
!      DBASE_VAR_CONSERVE_12

  integer, parameter :: dbNumVarAdvect =   DBASE_NUMBER_OF_ADVECT_VARIABLES + &
                                           DBASE_NUMBER_OF_SPECIES
  integer, parameter :: dbNumVarRenorm =   DBASE_NUMBER_OF_RENORM_VARIABLES
  integer, parameter :: dbNumVarConserve = DBASE_NUMBER_OF_CONSERVE_VARIABLES+&
                                           DBASE_NUMBER_OF_SPECIES

!---------------------------------------------------------------------------
!
! Variables to be set in dBaseInit 
! avoid using them, get local copy through dBaseKeyNumber instead
!
!---------------------------------------------------------------------------

  integer, target :: dbVarAdvect(dbNumVarAdvect)
  integer, target :: dbVarRenorm(dbNumVarRenorm)
  integer, target :: dbVarConserve(dbNumVarConserve)

  integer :: irhoflx
  integer :: iuflx
  integer :: ipflx
  integer :: iutflx
  integer :: iuttflx
  integer :: ieflx
  integer :: ieintflx

  integer, parameter :: inuc_begin = dbNumVarNamed+1
  integer, parameter :: inucflx_begin = dbNumFluxNamed+1

  integer :: iref1
  integer :: iref2
  integer :: iref3
  integer :: iref4


!---------------------------------------------------------------------------
!
! Local copied of some global parameters and their derivatives
!
!---------------------------------------------------------------------------

!!parameter used to make it easier for the code to deal with
!!different ways of defining the iu_bnd, etc parameters.
!!Especially if these are not set as parameters.

  integer, parameter :: MAXCELLS = max(iu_bnd, ju_bnd, ku_bnd)

  integer, parameter :: GC  = 1
  integer, parameter :: NGC = 0

! lower bounds with guardcells

  integer, parameter :: iLo_gc=1
  integer, parameter :: jLo_gc=1
  integer, parameter :: kLo_gc=1

! upper bounds with guardcells

  integer, parameter :: iHi_gc=nxb+2*nguard
  integer, parameter :: jHi_gc=nyb+2*nguard*k2d
  integer, parameter :: kHi_gc=nzb+2*nguard*k3d

! lower bounds with out guardcells

  integer, parameter :: iLo=nxlo
  integer, parameter :: jLo=nylo
  integer, parameter :: kLo=nzlo

! upper bounds with out guardcells

  integer, parameter :: iHi=nxhi
  integer, parameter :: jHi=nyhi
  integer, parameter :: kHi=nzhi

! block size

  integer, parameter :: dbXn = iHi_gc - iLo_gc + 1
  integer, parameter :: dbYn = jHi_gc - jLo_gc + 1
  integer, parameter :: dbZn = kHi_gc - kLo_gc + 1

! everything else

  integer, parameter :: dbSpeciesBegin = dbNumVarNamed + 1 
  integer, parameter :: dbSpeciesEnd   = dbNumVarNamed + dbNumSpecies

  integer, parameter :: dbMaxBlocks    = maxblocks

  integer, parameter :: dbNumZones     = nvar2

!---------------------------------------------------------------------------
!
! Numerical constants local to database (all private)
!
!---------------------------------------------------------------------------


!  Numerical constants to specify direction in dBase Get/Put Data

   integer, parameter :: dbiPoint   = 9000

   integer, parameter :: dbixVector = 9001
   integer, parameter :: dbiyVector = 9002
   integer, parameter :: dbizVector = 9003

   integer, parameter :: dbixyPlane = 9012
   integer, parameter :: dbixzPlane = 9013
   integer, parameter :: dbiyzPlane = 9023

   integer, parameter :: dbiyxPlane = 9021
   integer, parameter :: dbizxPlane = 9031
   integer, parameter :: dbizyPlane = 9032

   integer, parameter :: dbixyzCube = 9123

   integer, parameter :: dbiAllVariables = 9901
   integer, parameter :: dbiAllSpecies   = 9902


!  Numerical constants to specify direction in dBase Get/Put Coords

   integer, parameter :: dbixCoord  = 8001
   integer, parameter :: dbiyCoord  = 8002
   integer, parameter :: dbizCoord  = 8003

   integer, parameter :: dbiLowerFace = -1
   integer, parameter :: dbiUpperFace =  1

!  Numerical constants to specify variable in dBase Get/Put Coords

   integer, parameter :: dbiznl     = 1
   integer, parameter :: dbizn      = 2
   integer, parameter :: dbiznr     = 3
   integer, parameter :: dbiugrid   = 4
   integer, parameter :: dbidz      = 5

!---------------------------------------------------------------------------
!
! Numerical constants can be accessed through use dBase.  
!
!---------------------------------------------------------------------------

!  Pointers for the unk2 array (x coordinate information)

   integer, parameter :: ixznl      = dbiznl
   integer, parameter :: ixzn       = dbizn
   integer, parameter :: ixznr      = dbiznr
   integer, parameter :: iugridx    = dbiugrid
   integer, parameter :: idx        = dbidz

!   Pointers for the unk3 array (y coordinate information)

   integer, parameter :: iyznl      = dbiznl
   integer, parameter :: iyzn       = dbizn
   integer, parameter :: iyznr      = dbiznr
   integer, parameter :: iugridy    = dbiugrid
   integer, parameter :: idy        = dbidz

!  Pointers for the unk4 array (z coordinate information)

   integer, parameter :: izznl      = dbiznl
   integer, parameter :: izzn       = dbizn
   integer, parameter :: izznr      = dbiznr
   integer, parameter :: iugridz    = dbiugrid
   integer, parameter :: idz        = dbidz

!-------------------------------------------------------------------------------
!
!  Constants for use by the boundary condition code

   integer, parameter :: bnd_reflect    = -20
   integer, parameter :: bnd_outflow    = -21
   integer, parameter :: bnd_periodic   = -22
   integer, parameter :: bnd_user       = -23
   integer, parameter :: bnd_hydrostat  = -41

!-------------------------------------------------------------------------------
!
!  Constants for use by geometry-dependent coordinate type

   integer, parameter :: geom_cartesian = 0
   integer, parameter :: geom_planar   = 0
   integer, parameter :: geom_cylrad   = 1
   integer, parameter :: geom_sphrad   = 2
   integer, parameter :: geom_cylang   = 3
   integer, parameter :: geom_sphtheta = 4
   integer, parameter :: geom_sphphi   = 5

! the above constants are used for the different coordinate
! axes.  Define some parameters that describe the entire mesh
! geometry.
   integer, parameter :: CARTESIAN = 0
   integer, parameter :: CYLINDRICAL = 1
   integer, parameter :: SPHERICAL = 2
   integer, parameter :: POLAR = 3

   integer, save      :: meshGeometry

!-------------------------------------------------------------------------------
!
!  Constants to indicate which sweep we're doing, for
!  dimensionally split hydro solvers

   integer, parameter :: sweep_all       = 0
   integer, parameter :: sweep_x         = 1
   integer, parameter :: sweep_y         = 2
   integer, parameter :: sweep_z         = 3
   integer, parameter :: sweep_order_xyz = 1
   integer, parameter :: sweep_order_zyx = 2


!---------------------------------------------------------------------------
!
! Pointers for the unksm array (?)

  integer, parameter :: itemp_old = 1
  integer, parameter :: ishockt   = 2


  public
!!The data list is soo long.. I am making it all public for now
!----------------------------------------------------------------------
!
!  private
!  public ::ionmax, varnam, &
!       dbFluxNames, dbNamedVarAdvect, dbNamedVarRenorm, &
!       dbNamedVarConserve, dbNumVarAdvect, dbNumVarRenorm,   &
!       dbNumVarConserve, dbVarAdvect, dbVarRenorm, dbVarConserve,   &
!       irhoflx, iuflx, ipflx, iutflx, iuttflx,   &
!       ieflx, ieintflx, inuc_begin, inucflx_begin,   &
!       iref1, iref2, iref3, iref4, &
!       dbX1, dbY1, dbZ1, dbX2, dbY2, dbZ2,   &
!       dbXn, dbYn, dbZn,   &
!       dbNumVariables, dbNumSpecies, dbNumFluxes,   &
!       dbNumVarNamed. dbNumFluxNamed,   &
!       dbMaxBlocks, dbNumZones,   &
!       dbiPoint,   &
!       dbixVector, dbiyVector, dbizVector,   &
!       dbixyPlane, dbixzPlane, dbiyzPlane,   &
!       dbiyxPlane, dbizxPlane, dbizyPlane,   &
!       dbixyzCube, dbiAllVariables, dbiAllSpecies,   &
!       dbixCoord, dbiyCoord, dbizCoord,   &
!       dbiznl, dbizn, dbiznr, &
!       dbiugrid, &
!       ixznl, ixzn, ixznr, &
!       iugridx, iyznl, iyzn, iyznr, &
!       iugridy, izznl, izzn, izznr, &
!       iugridz,   &
!       bnd_reflect, bnd_outflow, bnd_periodic,   &
!       bnd_user, bnd_hydrostat,   &
!       CARTESIAN, CYLINDRICAL, SPHERICAL, POLAR, &
!       geom_cartesian, geom_planar, geom_cylrad, geom_sphrad,   &
!       geom_cylang, geom_sphtheta, geom_sphphi,   &
!       sweep_all, sweep_x, sweep_y, sweep_z,   &
!       sweep_order_xyz, sweep_order_zyx,   &
!       itemp_old, ishockt, 
!Stuff we actually create.. setting it from the code
!       time, dt, dtold, &
!       cpuseconds,
!       redshift, oldredshift, &
!       nstep, nbegin, MyPE, MasterPE, NumPEs


end module dBaseDeclarations
