module flashModules

  implicit none

  integer, PARAMETER :: NUM_MODULES = 8

  public :: NUM_MODULES, getFlashModules

contains

!!****f* object/flashModules
!!
!! NAME
!!
!!  getFlashModules
!!
!!
!! SYNOPSIS
!!
!!
!!  getFlashModules(module_names)
!!
!!  getFlashModules(character())
!!
!!
!! DESCRIPTION
!!
!!  Return a character array of size NUM_MODULES containing
!!  the names of all of the FLASH modules used to assemble
!!  the current executable
!!
!!  The module_names variable should be declared as
!!
!!    use flashModules
!!
!!  #include "flash_defines.fh"
!!    character (len=MAX_STRING_LENGTH) :: flash_modules(NUM_MODULES) 
!!
!!
!!  The length of each character string is set to MAX_STRING_LENGTH,
!!  which is defined in the automatically generated flash_defines.fh
!!
!!***
  subroutine getFlashModules(module_names)

    implicit none

    integer, parameter                :: MAX_STRING_LENGTH = 80
    character (len=MAX_STRING_LENGTH) :: module_names(NUM_MODULES)

!! temporary holds the result of the cat/cut from the setup_mod -- it is
!! dimensioned to be the same size as the result from the cut so we do
!! not overflow.  
    character (len=80) :: temporary

    temporary =  "Database"
    module_names(1) =  temporary(1:min(MAX_STRING_LENGTH,80))
  
    temporary =  "Driver"
    module_names(2) =  temporary(1:min(MAX_STRING_LENGTH,80))
  
    temporary =  "Grid"
    module_names(3) =  temporary(1:min(MAX_STRING_LENGTH,80))
  
    temporary =  "Main"
    module_names(4) =  temporary(1:min(MAX_STRING_LENGTH,80))
  
    temporary =  "RuntimeParameters"
    module_names(5) =  temporary(1:min(MAX_STRING_LENGTH,80))
  
    temporary =  "physics"
    module_names(6) =  temporary(1:min(MAX_STRING_LENGTH,80))
  
    temporary =  "physics/Hydro"
    module_names(7) =  temporary(1:min(MAX_STRING_LENGTH,80))
  
    temporary =  "../setups/example"
    module_names(8) =  temporary(1:min(MAX_STRING_LENGTH,80))
  

    return

  end subroutine getFlashModules
end module flashModules

