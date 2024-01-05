####################################################################
#  TNSP make include file.                                         #
####################################################################
#
#
#    Below is the fortran compiler
FC=mpif90
#
#
#    Below is the version of Tensor
#       If the version of the Tensor is older than 2.2,one should use 
#    the code "make Old" to compile the TNSP.Other version use "make"
#    to compile the TNSP
version=4.1.3
#
#
#    Below is the version of eigen_value
#       This package should link to arpack when used
#       If the version of the Tensor is older than 2.2, ignore this
#    line 
eigvalue_version=4.1
#
#
#    compiler options for your machine. It can not nothing. 
#        for complex number Optionl=-xavx
#        for real number Optionl=-O3
#     or
#        Optionl=
#
#    Below is the version of optimization Tools
optimization_version=2.1.9
#   if you want to backtrace the error, add the code below to Optionl
#   Optionl= -O3  -g -static -ffpe-trap=invalid
#
#   if you compile with ifort, you should use -traceback -g compiler option,
#   that is
#   Optionl= -g -traceback
#
Optionl= -O3   -g -static -ffpe-trap=invalid
#
# openmp options
openmp= 
#
#
#
#
# The name for output package
outputNameold= libTensor-$(version)-$(Symmetry_version)
outputName= libTensor-$(version)
#
#
IFIFORT=1
#
####################################################################
#                                                                  #
# Do not modify the code below                                     #
#                                                                  #
####################################################################
TensorVersion=$(version)
TNSDIR=./src/Tensor-$(version)
SYMM=./src/Symmetry-$(Symmetry_version)
EIG=./src/eigvalue-$(eigvalue_version)
OPTI=./src/optimization-$(optimization_version)
MPIDIR=./src/mpi
IFORTDIR=./src/ifort
ECHO    = echo


