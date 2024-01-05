#define ArrayDivideArrayFUNCNAME ArrayDivideArrayi
#define DATATYPE integer
#include "templet/ArrayDivideArray0.f90"
#undef ArrayDivideArrayFUNCNAME
#undef DATATYPE

#define ArrayDivideArrayFUNCNAME ArrayDivideArrays
#define DATATYPE real*4
#include "templet/ArrayDivideArray0.f90"
#undef ArrayDivideArrayFUNCNAME
#undef DATATYPE

#define ArrayDivideArrayFUNCNAME ArrayDivideArrayd
#define DATATYPE real*8
#include "templet/ArrayDivideArray0.f90"
#undef ArrayDivideArrayFUNCNAME
#undef DATATYPE

#define ArrayDivideArrayFUNCNAME ArrayDivideArrayc
#define DATATYPE complex*8
#include "templet/ArrayDivideArray0.f90"
#undef ArrayDivideArrayFUNCNAME
#undef DATATYPE

#define ArrayDivideArrayFUNCNAME ArrayDivideArrayz
#define DATATYPE complex*16
#include "templet/ArrayDivideArray0.f90"
#undef ArrayDivideArrayFUNCNAME
#undef DATATYPE



	subroutine ArrayDivideArray(Res,A,B,length)
		class(*),intent(inout)::Res(:)
		class(*),intent(in)::A(:),B(:)
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in ArrayDivideArray, length<0',-1)
			call error_stop
		end if
		select type(Res)
			type is (integer)
				call ArrayDivideArrayi(Res,A,B,length)
			type is (real(kind=4))
				call ArrayDivideArrays(Res,A,B,length)
			type is (real(kind=8))
				call ArrayDivideArrayd(Res,A,B,length)
			type is (complex(kind=4))
				call ArrayDivideArrayc(Res,A,B,length)
			type is (complex(kind=8))
				call ArrayDivideArrayz(Res,A,B,length)
			class default
				call Array_Divide_Array(Res,A,B,length)
		end select
	end subroutine

	subroutine default_Array_Divide_Array(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:),B(:)
		integer,intent(in)::length
		integer::i
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default Array_Time_Array is called   #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do i=1,length
			call numDividenum(Res(i),A(i),B(i))
		end do
		return
	end subroutine
	

