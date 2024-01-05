
PACKAGE= -ljemalloc  -larpack -llapack -lblas
Optionl=  -g -static -ffpe-trap=invalid
openmp=
PEPSOption= -fno-var-tracking-assignments


FC=mpif90
TensorDIRName=Tensor-4.1.3
TENDIR=../Tensor/lib$(TensorDIRName)
TensorDIRLink= -I$(TENDIR)
TensorlibLink= -L$(TENDIR)   -l$(TensorDIRName)    $(PACKAGE)



GM := density.o TargetElement.o TargetFunction.o  main.o

all: run 

run:	$(GM)
	$(FC)  $(openmp) -o $@ $^ $(TensorlibLink)

./%.o :	./%.f90
	$(FC) $(openmp)  -o $@ -c $(<)  $(TensorDIRLink)  $(Optionl)
clean:
	rm *.o *.mod $(all) run
