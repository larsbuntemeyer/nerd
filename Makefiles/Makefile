#
# makefile for NERD 
#

#F90 = gfotran
F90 = ifort

SOURCE  = ./source
INCLUDE = ./include

# compiler flags
CFLAGS  = -O3

# linker flags
LFLAGS  =

vpath %.F90 $(SOURCE) 
vpath %.h   $(INCLUDE)

# Link
nerd:   nerd.o $(OBJECTS) 
	$(F90) nerd.o $(OBJECTS) $(LFLAGS) -o nerd

include $(SOURCE)/Makefile.runtime_parameters
include $(SOURCE)/Makefile.driver
include $(SOURCE)/Makefile.database
include $(SOURCE)/Makefile.grid
include $(SOURCE)/Makefile.hydro

OBJECTS = $(driver) $(grid) $(database) $(hydro) \
          $(runtime_parameters)  

# Compile
nerd.o : $(OBJECTS)

%.o: %.F90 
	$(F90) -c $(CFLAGS) $<

clean:
	rm $(OBJECTS) nerd nerd.o *.mod 
