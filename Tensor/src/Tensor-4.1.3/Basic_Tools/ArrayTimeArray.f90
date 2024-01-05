#define ArrayTimeArrayFUNCNAME ArrayTimeArrayi
#define DATATYPE integer
#include "templet/ArrayTimeArray0.f90"
#undef ArrayTimeArrayFUNCNAME
#undef DATATYPE

#define ArrayTimeArrayFUNCNAME ArrayTimeArrays
#define DATATYPE real*4
#include "templet/ArrayTimeArray0.f90"
#undef ArrayTimeArrayFUNCNAME
#undef DATATYPE

#define ArrayTimeArrayFUNCNAME ArrayTimeArrayd
#define DATATYPE real*8
#include "templet/ArrayTimeArray0.f90"
#undef ArrayTimeArrayFUNCNAME
#undef DATATYPE

#define ArrayTimeArrayFUNCNAME ArrayTimeArrayc
#define DATATYPE complex*8
#include "templet/ArrayTimeArray0.f90"
#undef ArrayTimeArrayFUNCNAME
#undef DATATYPE

#define ArrayTimeArrayFUNCNAME ArrayTimeArrayz
#define DATATYPE complex*16
#include "templet/ArrayTimeArray0.f90"
#undef ArrayTimeArrayFUNCNAME
#undef DATATYPE



	subroutine ArrayTimeArray(Res,A,B,length)
		class(*),intent(inout)::Res(:)
		class(*),intent(in)::A(:),B(:)
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in ArrayTimeArray, length<0',-1)
			call error_stop
		end if
		select type(Res)
			type is (integer)
				call ArrayTimeArrayi(Res,A,B,length)
			type is (real(kind=4))
				call ArrayTimeArrays(Res,A,B,length)
			type is (real(kind=8))
				call ArrayTimeArrayd(Res,A,B,length)
			type is (complex(kind=4))
				call ArrayTimeArrayc(Res,A,B,length)
			type is (complex(kind=8))
				call ArrayTimeArrayz(Res,A,B,length)
			class default
				call Array_Time_Array(Res,A,B,length)
		end select
	end subroutine

	subroutine default_Array_Time_Array(Res,A,B,length)
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
			call numTimenum(Res(i),A(i),B(i))
		end do
		return
	end subroutine
	

