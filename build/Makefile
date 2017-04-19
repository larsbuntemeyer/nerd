#
# makefile for NERD 
#

F90 = gfortran
#F90 = ifort
EXE = nerd

SOURCE  = ../source
INCLUDE = ../include

# compiler flags
#CFLAGS  = -g -check all -fpe0 -warn -traceback -debug extended
#CFLAGS = -fbacktrace -fbounds-check -fcheck-array-temporaries -ffpe-trap=invalid,zero,overflow -Wall -Wextra
#CFLAGS = -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace \
         -ffree-line-length-0 -fcheck=all -finit-real=nan -ffpe-trap=zero,overflow,underflow -fdefault-real-8  
#CFLAGS = -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace \
#         -ffree-line-length-0 -fcheck=all -finit-real=nan -fdefault-real-8 
CFLAGS = -O3



# linker flags
LFLAGS  = -lpgplot -lX11

vpath %.F90 $(SOURCE) 
vpath %.h   $(INCLUDE)

# Link
$(EXE): $(EXE).o $(OBJECTS) 
	$(F90) $(EXE).o $(OBJECTS) $(LFLAGS) -o $(EXE)

include $(SOURCE)/Makefile.Database
include $(SOURCE)/Makefile.Eos
include $(SOURCE)/Makefile.Driver
include $(SOURCE)/Makefile.Grid
include $(SOURCE)/Makefile.Hydro
include $(SOURCE)/Makefile.RuntimeParameters
include $(SOURCE)/Makefile.Simulation
include $(SOURCE)/Makefile.Io
include $(SOURCE)/Makefile.PhysicalConstants

OBJECTS = $(Database) $(Driver) $(Grid) $(Hydro) \
          $(RuntimeParameters) $(Simulation) $(Io) \
          $(Eos) $(PhysicalConstants) 

# Compile
$(EXE).o : $(OBJECTS)

%.o: %.F90 
	$(F90) -c $(CFLAGS) $<

clean:
	rm $(OBJECTS) $(EXE) $(EXE).o *.mod 
