include ../../make.inc

LBFGSTool := LinearSearchTools.o LineSearchEnergy.o  OptimizationTools.o GeneralOptimizationElement.o GeneralLineSearchEnergy.o GeneralOptimizationEngine.o OptimizationTensorInterface.o




all: $(LBFGSTool) 

$(LBFGSTool):%.o : ./%.f90

./%.o :	./%.f90
	$(FC) $(openmp)  -o $@ -c $(<)  $(Optionl) -I../Tensor-$(TensorVersion) 
	
clean:
	rm -rf *.o *.mod *.a

