#define SetDataFUNCNAME SetDatai
#define DATATYPE integer
#include "templet/SetData0.f90"
#undef SetDataFUNCNAME
#undef DATATYPE

#define SetDataFUNCNAME SetDatas
#define DATATYPE real*4
#include "templet/SetData0.f90"
#undef SetDataFUNCNAME
#undef DATATYPE

#define SetDataFUNCNAME SetDatad
#define DATATYPE real*8
#include "templet/SetData0.f90"
#undef SetDataFUNCNAME
#undef DATATYPE

#define SetDataFUNCNAME SetDatac
#define DATATYPE complex*8
#include "templet/SetData0.f90"
#undef SetDataFUNCNAME
#undef DATATYPE

#define SetDataFUNCNAME SetDataz
#define DATATYPE complex*16
#include "templet/SetData0.f90"
#undef SetDataFUNCNAME
#undef DATATYPE

#define SetDataFUNCNAME SetDatal
#define DATATYPE logical
#include "templet/SetData0.f90"
#undef SetDataFUNCNAME
#undef DATATYPE

#define SetDataFUNCNAME SetDataa
#define DATATYPE character(len=*)
#include "templet/SetData0.f90"
#undef SetDataFUNCNAME
#undef DATATYPE






	subroutine SetClassData(outA,inA,length)
		class(*),target,intent(inout)::outA(:)
		class(*),target,intent(in)::inA
		integer,intent(in)::length
		integer::i
		select type(outA)
			type is (integer)
				call SetDatai(outA,inA,length)
			type is (real(kind=4))
				call SetDatas(outA,inA,length)
			type is (real(kind=8))
				call SetDatad(outA,inA,length)
			type is (complex(kind=4))
				call SetDatac(outA,inA,length)
			type is (complex(kind=8))
				call SetDataz(outA,inA,length)
			type is (logical)
				call SetDatal(outA,inA,length)
			type is (character(len=*))
				call SetDataa(outA,inA,length)
			class default
				call set_Class_Data(outA,inA,length)
		end select
	end subroutine

	subroutine default_Set_Class_Data(outA,inA,length)
		class(*),target,intent(inout)::outA(:)
		class(*),target,intent(in)::inA
		integer,intent(in)::length
		integer::i
		class(*),pointer::outp
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default  Set_Class_Data is called    #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do i=1,length
			call ClassPointer1DFunc(outA,i,outp)
			call copy_Class_Data(outp,inA)
		end do
		return
	end subroutine
