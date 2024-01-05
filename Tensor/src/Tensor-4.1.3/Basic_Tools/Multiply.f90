#include "Multiply/multiply11.f90"
#include "Multiply/multiply12.f90"
#include "Multiply/multiply13.f90"
#include "Multiply/multiply14.f90"
#include "Multiply/multiply15.f90"
#include "Multiply/multiply21.f90"
#include "Multiply/multiply22.f90"
#include "Multiply/multiply23.f90"
#include "Multiply/multiply24.f90"
#include "Multiply/multiply25.f90"
#include "Multiply/multiply31.f90"
#include "Multiply/multiply32.f90"
#include "Multiply/multiply33.f90"
#include "Multiply/multiply34.f90"
#include "Multiply/multiply35.f90"
#include "Multiply/multiply41.f90"
#include "Multiply/multiply42.f90"
#include "Multiply/multiply43.f90"
#include "Multiply/multiply44.f90"
#include "Multiply/multiply45.f90"
#include "Multiply/multiply51.f90"
#include "Multiply/multiply52.f90"
#include "Multiply/multiply53.f90"
#include "Multiply/multiply54.f90"
#include "Multiply/multiply55.f90"

#define ClassTimenumFUNCName ClassTimenumi
#define ClassTimenum2FUNCName ClassTimenum2i
#define ClassTimenumName ClassTimenumix
#define DATATYPE integer
#include "templet/multiply0.f90"
#undef ClassTimenumFUNCName
#undef ClassTimenum2FUNCName
#undef ClassTimenumName
#undef DATATYPE

#define ClassTimenumFUNCName ClassTimenums
#define ClassTimenum2FUNCName ClassTimenum2s
#define ClassTimenumName ClassTimenumsx
#define DATATYPE real*4
#include "templet/multiply0.f90"
#undef ClassTimenumFUNCName
#undef ClassTimenum2FUNCName
#undef ClassTimenumName
#undef DATATYPE

#define ClassTimenumFUNCName ClassTimenumd
#define ClassTimenum2FUNCName ClassTimenum2d
#define ClassTimenumName ClassTimenumdx
#define DATATYPE real*8
#include "templet/multiply0.f90"
#undef ClassTimenumFUNCName
#undef ClassTimenum2FUNCName
#undef ClassTimenumName
#undef DATATYPE

#define ClassTimenumFUNCName ClassTimenumc
#define ClassTimenum2FUNCName ClassTimenum2c
#define ClassTimenumName ClassTimenumcx
#define DATATYPE complex*8
#include "templet/multiply0.f90"
#undef ClassTimenumFUNCName
#undef ClassTimenum2FUNCName
#undef ClassTimenumName
#undef DATATYPE

#define ClassTimenumFUNCName ClassTimenumz
#define ClassTimenum2FUNCName ClassTimenum2z
#define ClassTimenumName ClassTimenumzx
#define DATATYPE complex*16
#include "templet/multiply0.f90"
#undef ClassTimenumFUNCName
#undef ClassTimenum2FUNCName
#undef ClassTimenumName
#undef DATATYPE

	subroutine classTimeNum(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:)
		class(*),target,intent(in)::B
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in multiply, length<0',-1)
			call error_stop
		end if
		select type(Res)
			type is (integer)
				call ClassTimenumi(Res,A,B,length)
			type is (real(kind=4))
				call ClassTimenums(Res,A,B,length)
			type is (real(kind=8))
				call ClassTimenumd(Res,A,B,length)
			type is (complex(kind=4))
				call ClassTimenumc(Res,A,B,length)
			type is (complex(kind=8))
				call ClassTimenumz(Res,A,B,length)
			class default
				call Class_Time_num(Res,A,B,length)
		end select
	end subroutine

	subroutine NumTimeClass(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A
		class(*),target,intent(in)::B(:)
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in multiply, length<0',-1)
			call error_stop
		end if
		select type(Res)
			type is (integer)
				call ClassTimenum2i(Res,A,B,length)
			type is (real(kind=4))
				call ClassTimenum2s(Res,A,B,length)
			type is (real(kind=8))
				call ClassTimenum2d(Res,A,B,length)
			type is (complex(kind=4))
				call ClassTimenum2c(Res,A,B,length)
			type is (complex(kind=8))
				call ClassTimenum2z(Res,A,B,length)
			class default
				call num_time_Class(Res,A,B,length)
		end select
	end subroutine

	subroutine default_Class_Time_num(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:)
		class(*),target,intent(in)::B
		integer,intent(in)::length
		integer::i
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default Class_Time_num is called   #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do i=1,length
			call Num_Time_num(Res(i),A(i),B)
		end do
	end subroutine

	subroutine default_num_time_Class(R,A,B,length)
		class(*),target,intent(inout)::R(:)
		class(*),target,intent(in)::A
		class(*),target,intent(in)::B(:)
		integer,intent(in)::length
		integer::i
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default num_time_Class is called     #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do i=1,length
			call Num_time_num(R(i),A,B(i))
		end do
	end subroutine