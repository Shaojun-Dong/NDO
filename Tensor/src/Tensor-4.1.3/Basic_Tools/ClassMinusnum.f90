#define ClassMinusnumFUNCNAME ClassMinusnumi
#define DATATYPE integer
#include "templet/ClassMinusnum0.f90"
#undef ClassMinusnumFUNCNAME
#undef DATATYPE

#define ClassMinusnumFUNCNAME ClassMinusnums
#define DATATYPE real*4
#include "templet/ClassMinusnum0.f90"
#undef ClassMinusnumFUNCNAME
#undef DATATYPE

#define ClassMinusnumFUNCNAME ClassMinusnumd
#define DATATYPE real*8
#include "templet/ClassMinusnum0.f90"
#undef ClassMinusnumFUNCNAME
#undef DATATYPE

#define ClassMinusnumFUNCNAME ClassMinusnumc
#define DATATYPE complex*8
#include "templet/ClassMinusnum0.f90"
#undef ClassMinusnumFUNCNAME
#undef DATATYPE

#define ClassMinusnumFUNCNAME ClassMinusnumz
#define DATATYPE complex*16
#include "templet/ClassMinusnum0.f90"
#undef ClassMinusnumFUNCNAME
#undef DATATYPE



	subroutine ClassMinusnum(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:)
		class(*),target,intent(in)::B
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in ClassMinusnum, length<0',-1)
			call error_stop
		end if
		select type(Res)
			type is (integer)
				call ClassMinusnumi(Res,A,B,length)
			type is (real(kind=4))
				call ClassMinusnums(Res,A,B,length)
			type is (real(kind=8))
				call ClassMinusnumd(Res,A,B,length)
			type is (complex(kind=4))
				call ClassMinusnumc(Res,A,B,length)
			type is (complex(kind=8))
				call ClassMinusnumz(Res,A,B,length)
			class default
				call Class_Minus_num(Res,A,B,length)
		end select
	end subroutine


	subroutine default_Class_Minus_num(R,A,B,length)
		class(*),target,intent(inout)::R(:)
		class(*),target,intent(in)::A(:)
		class(*),target,intent(in)::B
		integer,intent(in)::length
		integer::i
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default Class_Minus_num is called   #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do i=1,length
			call num_Minus_num(R(i),A(i),B)
		end do
	end subroutine


	

