#
# makefile for driver module
#

Io += Io.o Io_data.o Io_init.o 

Io.o : Io_data.o

Io_init.o : Io_data.o

