#include "Minus/ClassMinusClass11.f90"
#include "Minus/ClassMinusClass12.f90"
#include "Minus/ClassMinusClass13.f90"
#include "Minus/ClassMinusClass14.f90"
#include "Minus/ClassMinusClass15.f90"
#include "Minus/ClassMinusClass21.f90"
#include "Minus/ClassMinusClass22.f90"
#include "Minus/ClassMinusClass23.f90"
#include "Minus/ClassMinusClass24.f90"
#include "Minus/ClassMinusClass25.f90"
#include "Minus/ClassMinusClass31.f90"
#include "Minus/ClassMinusClass32.f90"
#include "Minus/ClassMinusClass33.f90"
#include "Minus/ClassMinusClass34.f90"
#include "Minus/ClassMinusClass35.f90"
#include "Minus/ClassMinusClass41.f90"
#include "Minus/ClassMinusClass42.f90"
#include "Minus/ClassMinusClass43.f90"
#include "Minus/ClassMinusClass44.f90"
#include "Minus/ClassMinusClass45.f90"
#include "Minus/ClassMinusClass51.f90"
#include "Minus/ClassMinusClass52.f90"
#include "Minus/ClassMinusClass53.f90"
#include "Minus/ClassMinusClass54.f90"
#include "Minus/ClassMinusClass55.f90"

#define ClassMinusClassFUNCName ClassMinusClassi
#define ClassMinusClassName ClassMinusClassix
#define DATATYPE integer
#include "templet/Minus0.f90"
#undef ClassMinusClassFUNCName
#undef ClassMinusClassName
#undef DATATYPE

#define ClassMinusClassFUNCName ClassMinusClasss
#define ClassMinusClassName ClassMinusClasssx
#define DATATYPE real*4
#include "templet/Minus0.f90"
#undef ClassMinusClassFUNCName
#undef ClassMinusClassName
#undef DATATYPE

#define ClassMinusClassFUNCName ClassMinusClassd
#define ClassMinusClassName ClassMinusClassdx
#define DATATYPE real*8
#include "templet/Minus0.f90"
#undef ClassMinusClassFUNCName
#undef ClassMinusClassName
#undef DATATYPE

#define ClassMinusClassFUNCName ClassMinusClassc
#define ClassMinusClassName ClassMinusClasscx
#define DATATYPE complex*8
#include "templet/Minus0.f90"
#undef ClassMinusClassFUNCName
#undef ClassMinusClassName
#undef DATATYPE

#define ClassMinusClassFUNCName ClassMinusClassz
#define ClassMinusClassName ClassMinusClasszx
#define DATATYPE complex*16
#include "templet/Minus0.f90"
#undef ClassMinusClassFUNCName
#undef ClassMinusClassName
#undef DATATYPE

	subroutine ClassMinusClass(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:),B(:)
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in ClassMinusClass, length<0',-1)
			call error_stop
		end if
		select type(Res)
			type is (integer)
				call ClassMinusClassi(Res,A,B,length)
			type is (real(kind=4))
				call ClassMinusClasss(Res,A,B,length)
			type is (real(kind=8))
				call ClassMinusClassd(Res,A,B,length)
			type is (complex(kind=4))
				call ClassMinusClassc(Res,A,B,length)
			type is (complex(kind=8))
				call ClassMinusClassz(Res,A,B,length)
			class default
				call Class_Minus_Class(Res,A,B,length)
		end select
	end subroutine


	subroutine default_ClassMinusClass(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:),B(:)
		integer,intent(in)::length
		integer::i
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default Class_Minus_Class is called  #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do i=1,length
			call num_Minus_num(Res(i),A(i),B(i))
		end do
		return
	end subroutine