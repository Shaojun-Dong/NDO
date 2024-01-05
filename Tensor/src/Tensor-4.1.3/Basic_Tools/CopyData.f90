#define copyDataFUNCNAME copyDatai
#define DATATYPE integer
#include "templet/CopyData0.f90"
#undef copyDataFUNCNAME
#undef DATATYPE

#define copyDataFUNCNAME copyDatas
#define DATATYPE real*4
#include "templet/CopyData0.f90"
#undef copyDataFUNCNAME
#undef DATATYPE

#define copyDataFUNCNAME copyDatad
#define DATATYPE real*8
#include "templet/CopyData0.f90"
#undef copyDataFUNCNAME
#undef DATATYPE

#define copyDataFUNCNAME copyDatac
#define DATATYPE complex*8
#include "templet/CopyData0.f90"
#undef copyDataFUNCNAME
#undef DATATYPE

#define copyDataFUNCNAME copyDataz
#define DATATYPE complex*16
#include "templet/CopyData0.f90"
#undef copyDataFUNCNAME
#undef DATATYPE

#define copyDataFUNCNAME copyDatal
#define DATATYPE logical
#include "templet/CopyData0.f90"
#undef copyDataFUNCNAME
#undef DATATYPE

#define copyDataFUNCNAME copyDataa
#define DATATYPE character(len=*)
#include "templet/CopyData0.f90"
#undef copyDataFUNCNAME
#undef DATATYPE






	subroutine copyData(outA,inA)
		class(*),intent(in)::inA
		class(*),intent(inout)::outA
		integer::i
		select type(outA)
			type is (integer)
				call copyDatai(outA,inA)
			type is (real(kind=4))
				call copyDatas(outA,inA)
			type is (real(kind=8))
				call copyDatad(outA,inA)
			type is (complex(kind=4))
				call copyDatac(outA,inA)
			type is (complex(kind=8))
				call copyDataz(outA,inA)
			type is (logical)
				call copyDatal(outA,inA)
			type is (character(len=*))
				call copyDataa(outA,inA)
			class default
				call copy_Class_Data(outA,inA)
		end select
	end subroutine

	subroutine default_copy_Class_Data(outA,inA)
		class(*),target,intent(in)::inA
		class(*),target,intent(inout)::outA
		call writemess('ERROR copy class data',-1)
		call writemess('DO NOT setting the default subroutine of copy_Class_Data',-1)
		call error_stop
	end subroutine
