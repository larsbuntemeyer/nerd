# FLASH makefile definitions for the Intel ifc compilers on Linux

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

HDF5_PATH = $(HDF5_DIR)#$(HOME)/hdf5-1.6.9/hdf5
#HDF5_PATH = $(HDF5_BASE)
SZ_PATH = /usr/local/szip

# MPICH_PATH = /usr/mpi/gcc/mvapich2-1.0.3
# MPICH_PATH = /usr/mpi/gcc/openmpi-1.2.6
#MPICH_PATH = /usr/mpi/intel-9.1/openmpi-1.2.6

BOOST_PATH = /usr/include/boost
#SAPPORO_PATH = /usr/src/sapporo_v1.5
SAPPORO_PATH = /usr/src/sapporo_light
CUDA_PATH = /usr/local/cuda


#----------------------------------------------------------------------------
# Compiler and linker commands
#
#  We use the f90 compiler as the linker, so some C libraries may explicitly
#  need to be added into the link line.
#----------------------------------------------------------------------------

FCOMP   = ftn 
CCOMP   = cc
CPPCOMP = CC
LINK    = ftn 

PP     = -D


#----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized
#  code, one for testing, and one for debugging.  The default is to use the 
#  _OPT version.  Specifying -debug to setup will pick the _DEBUG version,
#  these should enable bounds checking.  Specifying _TEST is used for 
#  flash_test, and is set for quick code generation, and (sometimes) 
#  profiling.  The Makefile generated by setup will assign the generic token 
#  (ex. FFLAGS) to the proper set of flags (ex. FFLAGS_OPT).
#----------------------------------------------------------------------------

#FFLAGS_OPT   =  -c -i4 -r8 -mcmodel=large -i-dynamic -O3 -I/usr/local/include
FFLAGS_OPT   =  -c -O3 -i4 -r8 -mcmodel=large -dynamic -I/usr/local/include 
#FFLAGS_OPT   = -c -O3 -mcmodel=large -I/usr/local/include
FFLAGS_DEBUG =  -c -i4 -r8
FFLAGS_TEST  =  -c -i4 -r8


CFLAGS_HDF5 = -I$(HDF5_PATH)/include

#CFLAGS_OPT   = -c -O3 -mcmodel=large -i-dynamic 
CFLAGS_OPT   = -c -O3 -i-dynamic
CFLAGS_DEBUG = -c -g
CFLAGS_TEST  = -c -02


#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT, 
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -mcmodel=large -dynamic -o 
#LFLAGS_OPT   = -o -mcmodel=large
LFLAGS_DEBUG = -o
LFLAGS_TEST  = -o


#----------------------------------------------------------------------------
# Library specific linking
#
#  If a FLASH module has a 'LIBRARY xxx' line in its Config file, we need to
#  create a macro in this Makefile.h for LIB_xxx, which will be added to the
#  link line when FLASH is built.  This allows us to switch between different
#  (incompatible) libraries.  We also create a _OPT, _DEBUG, and _TEST
#  library macro to add any performance-minded libraries (like fast math),
#  depending on how FLASH was setup.
#----------------------------------------------------------------------------

LIB_HDF5 = -L$(HDF5_PATH)/lib -lhdf5 -lz -L$(SZ_PATH)/lib
LIB_GPU  = -I$(SAPPORO_PATH)  -L$(SAPPORO_PATH) -lsapporo -L/usr/lib64 /usr/lib64/libboost_thread.a -L$(CUDA_PATH)/lib -lcuda -lcudart

LIB_OPT = 
#LIB_OPT = -L$(MPICH_PATH)/lib
LIB_DEBUG = 
LIB_TEST =


#----------------------------------------------------------------------------
# Additional machine-dependent object files
#
#  Add any machine specific files here -- they will be compiled and linked
#  when FLASH is built.
#----------------------------------------------------------------------------

MACHOBJ = 


#----------------------------------------------------------------------------
# Additional commands
#---------------------------------------------------------------------------- 

MV = mv -f
AR = ar -r
RM = rm -f
CD = cd
RL = ranlib
ECHO = echo