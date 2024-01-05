#define NDivideNFUNCNAME NumDivideNumi
#define DATATYPE integer
#include "templet/NumDivideNum0.f90"
#undef NDivideNFUNCNAME
#undef DATATYPE

#define NDivideNFUNCNAME NumDivideNums
#define DATATYPE real*4
#include "templet/NumDivideNum0.f90"
#undef NDivideNFUNCNAME
#undef DATATYPE

#define NDivideNFUNCNAME NumDivideNumd
#define DATATYPE real*8
#include "templet/NumDivideNum0.f90"
#undef NDivideNFUNCNAME
#undef DATATYPE

#define NDivideNFUNCNAME NumDivideNumc
#define DATATYPE complex*8
#include "templet/NumDivideNum0.f90"
#undef NDivideNFUNCNAME
#undef DATATYPE

#define NDivideNFUNCNAME NumDivideNumz
#define DATATYPE complex*16
#include "templet/NumDivideNum0.f90"
#undef NDivideNFUNCNAME
#undef DATATYPE



	subroutine NumDivideNum(Res,A,B)
		class(*),target,intent(inout)::Res
		class(*),target,intent(in)::A
		class(*),target,intent(in)::B
		select type(Res)
			type is (integer)
				call NumDivideNumi(Res,A,B)
			type is (real(kind=4))
				call NumDivideNums(Res,A,B)
			type is (real(kind=8))
				call NumDivideNumd(Res,A,B)
			type is (complex(kind=4))
				call NumDivideNumc(Res,A,B)
			type is (complex(kind=8))
				call NumDivideNumz(Res,A,B)
			class default
				call num_divide_num(Res,A,B)
		end select
	end subroutine



	subroutine default_num_divide_num(R,A,B)
		class(*),target,intent(inout)::R
		class(*),target,intent(in)::A
		class(*),target,intent(in)::B
		call writemess('ERROR (/) for class data',-1)
		call writemess('DO NOT setting the default subroutine of num_divide_num',-1)
		call error_stop
	end subroutine

	

