
# name of the exectuable
EXE       = nerd 

# compiler and flags
FCOMP      = gfortran 
F90FLAGS  = -ffree-form
F77FLAGS  = -ffixed-form
FFLAGS_OPT = -cpp -O3 -c -fconvert=big-endian -fdefault-real-8 -I/usr/include
LFLAGS_OPT = -L/usr/lib -lnetcdff -lnetcdf  -o 
LINK       = $(FCOMP)
#FFLAGS_OPT  = -cpp -c -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -finit-real=nan -fdefault-real-8 

# netcdf flags 
INCLUDE  =  -I/usr/include
LIB     = -L/usr/lib -lnetcdff -lnetcdf


ECHO = echo
