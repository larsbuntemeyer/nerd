
# name of the exectuable
EXE       = nerd 

# compiler and flags
FCOMP      = gfortran 
F90FLAGS  = -ffree-form
F77FLAGS  = -ffixed-form
FFLAGS_OPT = -cpp -c -fconvert=big-endian
LFLAGS_OPT = -o
LINK       = $(FCOMP)
FFLAGS_OPT  = -c -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -finit-real=nan -fdefault-real-8 

# netcdf flags 
NCINCLUDE  = 
NCLIBS     = 


ECHO = echo
