
# name of the exectuable
EXE       = nerd

# compiler
FCOMP      = ifort
LINK       = $(FCOMP)

################### flag setting ###########################

# netcdf flags (if necessary, please uncomment)
NCINCLUDE  = -I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-parallel-bullxmpi-intel14/include -DNETCDF 
NCLIBS     = -L/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-parallel-bullxmpi-intel14/lib -lnetcdff -L/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.0-bullxmpi-intel14/lib -L/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-parallel-bullxmpi-intel14/lib -Wl,-rpath,/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-parallel-bullxmpi-intel14/lib -lnetcdf -L/sw/rhel6-x64/hdf5/hdf5-1.8.14-parallel-bullxmpi-intel14/lib -Wl,-rpath,/sw/rhel6-x64/hdf5/hdf5-1.8.14-parallel-bullxmpi-intel14/lib -lhdf5 -lhdf5_hl -L/sw/rhel6-x64/sys/libaec-0.3.2-intel14/lib -Wl,-rpath,/sw/rhel6-x64/sys/libaec-0.3.2-intel14/lib -lsz -lz -lcurl -Wl,-rpath,/sw/rhel6-x64/netcdf/parallel_netcdf-1.6.0-bullxmpi-intel14/lib -lnetcdf


################### flag setting ###########################

# some special flags for free and fixed format
F90FLAGS   = -free
F77FLAGS   = -fixed

# general fortran compiler flags
#FFLAGS_OPT = -c -O3 -convert big_endian -r8 -fpp $(NCINCLUDE) #-check all -debug all -traceback 
FFLAGS_OPT = -c -fpp -check all -debug all -traceback $(NCINCLUDE) 

# linker flags
LFLAGS_OPT = $(NCLIBS) -o 

ECHO = echo
