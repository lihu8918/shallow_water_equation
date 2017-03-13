SRC = var.f90 proc.f90 main.f90
OBJ = *.mod
FC  = ifort

ocean    : $(SRC)
	$(FC) -O2 -o ocean $(SRC)

clean:
	rm $(OBJ) ocean

