#define ClassPlusnumFUNCNAME ClassPlusnumi
#define DATATYPE integer
#include "templet/ClassPlusnum0.f90"
#undef ClassPlusnumFUNCNAME
#undef DATATYPE

#define ClassPlusnumFUNCNAME ClassPlusnums
#define DATATYPE real*4
#include "templet/ClassPlusnum0.f90"
#undef ClassPlusnumFUNCNAME
#undef DATATYPE

#define ClassPlusnumFUNCNAME ClassPlusnumd
#define DATATYPE real*8
#include "templet/ClassPlusnum0.f90"
#undef ClassPlusnumFUNCNAME
#undef DATATYPE

#define ClassPlusnumFUNCNAME ClassPlusnumc
#define DATATYPE complex*8
#include "templet/ClassPlusnum0.f90"
#undef ClassPlusnumFUNCNAME
#undef DATATYPE

#define ClassPlusnumFUNCNAME ClassPlusnumz
#define DATATYPE complex*16
#include "templet/ClassPlusnum0.f90"
#undef ClassPlusnumFUNCNAME
#undef DATATYPE


#define ClassPlusnumFUNCNAME ClassPlusnuma
#define DATATYPE character(len=*)
#include "templet/ClassPlusnum0.f90"
#undef ClassPlusnumFUNCNAME
#undef DATATYPE


	subroutine ClassPlusnum(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:)
		class(*),target,intent(in)::B
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in ClassPlusnum, length<0',-1)
			call error_stop
		end if
		select type(Res)
			type is (integer)
				call ClassPlusnumi(Res,A,B,length)
			type is (real(kind=4))
				call ClassPlusnums(Res,A,B,length)
			type is (real(kind=8))
				call ClassPlusnumd(Res,A,B,length)
			type is (complex(kind=4))
				call ClassPlusnumc(Res,A,B,length)
			type is (complex(kind=8))
				call ClassPlusnumz(Res,A,B,length)
			type is (character(len=*))
				call ClassPlusnuma(Res,A,B,length)
			class default
				call Class_Plus_num(Res,A,B,length)
		end select
	end subroutine

	subroutine default_Class_Plus_num(R,A,B,length)
		class(*),target,intent(inout)::R(:)
		class(*),target,intent(in)::A(:)
		class(*),target,intent(in)::B
		integer,intent(in)::length
		integer::i
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default Num_Plus_num is called   #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do i=1,length
			call Num_Plus_num(R(i),A(i),B)
		end do
	end subroutine

	

