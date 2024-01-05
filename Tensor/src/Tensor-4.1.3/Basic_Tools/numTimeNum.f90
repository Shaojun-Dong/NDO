#define NumTimenumFUNCNAME NumTimenumi
#define DATATYPE integer
#include "templet/NumTimeNum0.f90"
#undef NumTimenumFUNCNAME
#undef DATATYPE

#define NumTimenumFUNCNAME NumTimenums
#define DATATYPE real*4
#include "templet/NumTimeNum0.f90"
#undef NumTimenumFUNCNAME
#undef DATATYPE

#define NumTimenumFUNCNAME NumTimenumd
#define DATATYPE real*8
#include "templet/NumTimeNum0.f90"
#undef NumTimenumFUNCNAME
#undef DATATYPE

#define NumTimenumFUNCNAME NumTimenumc
#define DATATYPE complex*8
#include "templet/NumTimeNum0.f90"
#undef NumTimenumFUNCNAME
#undef DATATYPE

#define NumTimenumFUNCNAME NumTimenumz
#define DATATYPE complex*16
#include "templet/NumTimeNum0.f90"
#undef NumTimenumFUNCNAME
#undef DATATYPE



	subroutine NumTimenum(Res,A,B)
		class(*),target,intent(inout)::Res
		class(*),target,intent(in)::A
		class(*),target,intent(in)::B
		select type(Res)
			type is (integer)
				call NumTimenumi(Res,A,B)
			type is (real(kind=4))
				call NumTimenums(Res,A,B)
			type is (real(kind=8))
				call NumTimenumd(Res,A,B)
			type is (complex(kind=4))
				call NumTimenumc(Res,A,B)
			type is (complex(kind=8))
				call NumTimenumz(Res,A,B)
			class default
				call Num_Time_num(Res,A,B)
		end select
	end subroutine

	subroutine default_Num_Time_num(R,A,B)
		class(*),target,intent(inout)::R
		class(*),target,intent(in)::A
		class(*),target,intent(in)::B
		call writemess('ERROR (*) for class data',-1)
		call writemess('DO NOT setting the default subroutine of Num_Time_num',-1)
		call error_stop
	end subroutine

