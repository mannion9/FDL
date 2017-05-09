CF 			= gfortran
objects = Main.o io.o Constants.o Init.o Boundary.o  Timestep.o
files   = Main.f90 io.f90 Init.f90 Boundary.f90  Timestep.f90

exe : $(objects)
	@$(CF) -o exe $(objects)
Constants.mod : Constants.o Constants.f90
	@$(CF) -c Constants.f90
Constants.o : Constants.f90
	@$(CF) -c Constants.f90
Main.o : $(files) Constants.mod
	@$(CF) -c Main.f90
io.o : io.f90 Constants.mod
	@$(CF) -c io.f90
Init.o : Init.f90 Constants.mod
	@$(CF) -c Init.f90
Boundary.o : Boundary.f90 Constants.mod
	@$(CF) -c Boundary.f90
Timestep.o : Timestep.f90 Constants.mod
	@$(CF) -c Timestep.f90
