#define NumMinusNumFUNCNAME NumMinusNumi
#define DATATYPE integer
#include "templet/NumMinusNum0.f90"
#undef NumMinusNumFUNCNAME
#undef DATATYPE

#define NumMinusNumFUNCNAME NumMinusNums
#define DATATYPE real*4
#include "templet/NumMinusNum0.f90"
#undef NumMinusNumFUNCNAME
#undef DATATYPE

#define NumMinusNumFUNCNAME NumMinusNumd
#define DATATYPE real*8
#include "templet/NumMinusNum0.f90"
#undef NumMinusNumFUNCNAME
#undef DATATYPE

#define NumMinusNumFUNCNAME NumMinusNumc
#define DATATYPE complex*8
#include "templet/NumMinusNum0.f90"
#undef NumMinusNumFUNCNAME
#undef DATATYPE

#define NumMinusNumFUNCNAME NumMinusNumz
#define DATATYPE complex*16
#include "templet/NumMinusNum0.f90"
#undef NumMinusNumFUNCNAME
#undef DATATYPE








	subroutine NumMinusNum(R,A,B)
		class(*),target,intent(inout)::R
		class(*),target,intent(in)::A
		class(*),target,intent(in)::B
		select type(R)
			type is (integer)
				call NumMinusNumi(R,A,B)
			type is (real(kind=4))
				call NumMinusNums(R,A,B)
			type is (real(kind=8))
				call NumMinusNumd(R,A,B)
			type is (complex(kind=4))
				call NumMinusNumc(R,A,B)
			type is (complex(kind=8))
				call NumMinusNumz(R,A,B)
			class default
				call Num_Minus_num(R,A,B)
		end select
	end subroutine

	subroutine default_num_Minus_num(R,A,B)
		class(*),target,intent(inout)::R
		class(*),target,intent(in)::A
		class(*),target,intent(in)::B
		call writemess('ERROR (-) for class data',-1)
		call writemess('DO NOT setting the default subroutine of Num_Minus_num',-1)
		call error_stop
	end subroutine
