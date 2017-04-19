# Makefile for the FLASH code (driver module)

driver = main.o init_driver.o evolve.o finish_driver.o init_namelist.o mo_namelist.o init_domain.o 

evolve.o: mo_namelist.o mo_database.o mo_grid.o

init_namelist.o: mo_namelist.o
