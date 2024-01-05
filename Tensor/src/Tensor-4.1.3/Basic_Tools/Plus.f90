#include "Plus/ClassPlusClass11.f90"
#include "Plus/ClassPlusClass12.f90"
#include "Plus/ClassPlusClass13.f90"
#include "Plus/ClassPlusClass14.f90"
#include "Plus/ClassPlusClass15.f90"
#include "Plus/ClassPlusClass21.f90"
#include "Plus/ClassPlusClass22.f90"
#include "Plus/ClassPlusClass23.f90"
#include "Plus/ClassPlusClass24.f90"
#include "Plus/ClassPlusClass25.f90"
#include "Plus/ClassPlusClass31.f90"
#include "Plus/ClassPlusClass32.f90"
#include "Plus/ClassPlusClass33.f90"
#include "Plus/ClassPlusClass34.f90"
#include "Plus/ClassPlusClass35.f90"
#include "Plus/ClassPlusClass41.f90"
#include "Plus/ClassPlusClass42.f90"
#include "Plus/ClassPlusClass43.f90"
#include "Plus/ClassPlusClass44.f90"
#include "Plus/ClassPlusClass45.f90"
#include "Plus/ClassPlusClass51.f90"
#include "Plus/ClassPlusClass52.f90"
#include "Plus/ClassPlusClass53.f90"
#include "Plus/ClassPlusClass54.f90"
#include "Plus/ClassPlusClass55.f90"
#include "Plus/ClassPlusClass71.f90"
#include "Plus/ClassPlusClass72.f90"
#include "Plus/ClassPlusClass73.f90"
#include "Plus/ClassPlusClass74.f90"
#include "Plus/ClassPlusClass75.f90"
#include "Plus/ClassPlusClass76.f90"
#include "Plus/ClassPlusClass77.f90"

#define ClassPlusClassFUNCName ClassPlusClassi
#define ClassPlusClassName ClassPlusClassix
#define DATATYPE integer
#include "templet/Plus0.f90"
#undef ClassPlusClassFUNCName
#undef ClassPlusClassName
#undef DATATYPE

#define ClassPlusClassFUNCName ClassPlusClasss
#define ClassPlusClassName ClassPlusClasssx
#define DATATYPE real*4
#include "templet/Plus0.f90"
#undef ClassPlusClassFUNCName
#undef ClassPlusClassName
#undef DATATYPE

#define ClassPlusClassFUNCName ClassPlusClassd
#define ClassPlusClassName ClassPlusClassdx
#define DATATYPE real*8
#include "templet/Plus0.f90"
#undef ClassPlusClassFUNCName
#undef ClassPlusClassName
#undef DATATYPE

#define ClassPlusClassFUNCName ClassPlusClassc
#define ClassPlusClassName ClassPlusClasscx
#define DATATYPE complex*8
#include "templet/Plus0.f90"
#undef ClassPlusClassFUNCName
#undef ClassPlusClassName
#undef DATATYPE

#define ClassPlusClassFUNCName ClassPlusClassz
#define ClassPlusClassName ClassPlusClasszx
#define DATATYPE complex*16
#include "templet/Plus0.f90"
#undef ClassPlusClassFUNCName
#undef ClassPlusClassName
#undef DATATYPE
	
	subroutine ClassPlusClassa(Res,A,B,length)
		character(len=*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:)
		class(*),target,intent(in)::B(:)
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in ClassPlusClass, length<0',-1)
			call error_stop
		end if
		select type(A)
			type is (integer)
				call ClassPlusClassax(Res,A,B,length)
			type is (real(kind=4))
				call ClassPlusClassax(Res,A,B,length)
			type is (real(kind=8))
				call ClassPlusClassax(Res,A,B,length)
			type is (complex(kind=4))
				call ClassPlusClassax(Res,A,B,length)
			type is (complex(kind=8))
				call ClassPlusClassax(Res,A,B,length)
			type is (logical)
				call ClassPlusClassax(Res,A,B,length)
			type is (character(len=*))
				call ClassPlusClassax(Res,A,B,length)
			class default
				call Class_Plus_Class(Res,A,B,length)
		end select
	end subroutine

	subroutine ClassPlusClass(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:),B(:)
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in ClassPlusClass, length<0',-1)
			call error_stop
		end if
		select type(Res)
			type is (integer)
				call ClassPlusClassi(Res,A,B,length)
			type is (real(kind=4))
				call ClassPlusClasss(Res,A,B,length)
			type is (real(kind=8))
				call ClassPlusClassd(Res,A,B,length)
			type is (complex(kind=4))
				call ClassPlusClassc(Res,A,B,length)
			type is (complex(kind=8))
				call ClassPlusClassz(Res,A,B,length)
			type is (character(len=*))
				call ClassPlusClassa(Res,A,B,length)
			class default
				call Class_Plus_Class(Res,A,B,length)
		end select
	end subroutine


	subroutine default_ClassPlusClass(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:),B(:)
		integer,intent(in)::length
		integer::i
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default Class_Plus_Class is called   #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do i=1,length
			call num_plus_num(Res(i),A(i),B(i))
		end do
		return
	end subroutine