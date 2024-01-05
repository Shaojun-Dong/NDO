#define NumPlusNumFUNCNAME NumPlusNumi
#define DATATYPE integer
#include "templet/NumPlusNum0.f90"
#undef NumPlusNumFUNCNAME
#undef DATATYPE

#define NumPlusNumFUNCNAME NumPlusNums
#define DATATYPE real*4
#include "templet/NumPlusNum0.f90"
#undef NumPlusNumFUNCNAME
#undef DATATYPE

#define NumPlusNumFUNCNAME NumPlusNumd
#define DATATYPE real*8
#include "templet/NumPlusNum0.f90"
#undef NumPlusNumFUNCNAME
#undef DATATYPE

#define NumPlusNumFUNCNAME NumPlusNumc
#define DATATYPE complex*8
#include "templet/NumPlusNum0.f90"
#undef NumPlusNumFUNCNAME
#undef DATATYPE

#define NumPlusNumFUNCNAME NumPlusNumz
#define DATATYPE complex*16
#include "templet/NumPlusNum0.f90"
#undef NumPlusNumFUNCNAME
#undef DATATYPE


#include "Plus/NumPlusNuma.f90"






	subroutine NumPlusNum(R,A,B)
			class(*),target,intent(inout)::R
			class(*),target,intent(in)::A
			class(*),target,intent(in)::B
		select type(R)
			type is (integer)
				call NumPlusNumi(R,A,B)
			type is (real(kind=4))
				call NumPlusNums(R,A,B)
			type is (real(kind=8))
				call NumPlusNumd(R,A,B)
			type is (complex(kind=4))
				call NumPlusNumc(R,A,B)
			type is (complex(kind=8))
				call NumPlusNumz(R,A,B)
			type is (character(len=*))
				call NumPlusNuma(R,A,B)
			class default
				call Num_Plus_num(R,A,B)
		end select
	end subroutine

	subroutine default_num_plus_num(R,A,B)
			class(*),target,intent(inout)::R
			class(*),target,intent(in)::A
			class(*),target,intent(in)::B
		call writemess('ERROR (+) for class data',-1)
		call writemess('DO NOT setting the default subroutine of Num_Plus_num',-1)
		call error_stop
	end subroutine
