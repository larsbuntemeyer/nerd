
# name of the exectuable
EXE       = nerd 

# compiler and flags
FCOMP      = gfortran 
F90FLAGS  = -ffree-form
F77FLAGS  = -ffixed-form
FFLAGS_OPT = -cpp -c -fconvert=big-endian
LFLAGS_OPT = -o
LINK       = $(FCOMP)

# netcdf flags 
NCINCLUDE  = 
NCLIBS     = 


ECHO = echo
