include $(HOME)/makefileConf

Tensor := Tools.o  pointer_Tools.o  Basic_Tools.o  data_base.o    QuantumNumber.o dimension.o \
			 U1.o Parity.o Tensor.o parameter.o Fun_for_Tensor.o



all: $(Tensor) 

$(Tensor):%.o : ./%.f90

./%.o :	./%.f90
	$(FC) $(openmp)  -o $@ -c $(<) -cpp -c -Ofast $(Optionl) $(PACKAGE)
	
clean:
	rm -rf *.o *.mod *.a

