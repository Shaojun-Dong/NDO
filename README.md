

This collection of FORTRAN files is aimed at trainning the Neural Density Operators (NDO) for quantum tomography straightforward, simple and fast. The most important files or  documents are:

------------------------------

Tensor/               : The basic code of tensor and optimization tools.

NDO/density.f90       : define the reduced density matrix from NDOs.

NDO/TargetElement.f90 : define the basic operations that are used in Tensor/src/optimization-2.1.9/GeneralOptimizationElement.f90.

NDO/TargetFunction.f90: define the cost function of our problems.

Datasets.nb           : Generation of synthetic datasets for different types of quantum walk systems.

------------------------------

To compile the code, do as following:

1. Make the tensor lib with the code below

cd Tensor 

make

2. Go to NDO and modify the file of makefile, set the value "PACKAGE" to link with lapack(or MKL)

PACKAGE:  links with  lapack and blas

and then input

make


3. Then you can run the NDO with "./run"

