#define numMinusClassFUNCNAME numMinusClassi
#define DATATYPE integer
#include "templet/numMinusClass0.f90"
#undef numMinusClassFUNCNAME
#undef DATATYPE

#define numMinusClassFUNCNAME numMinusClasss
#define DATATYPE real*4
#include "templet/numMinusClass0.f90"
#undef numMinusClassFUNCNAME
#undef DATATYPE

#define numMinusClassFUNCNAME numMinusClassd
#define DATATYPE real*8
#include "templet/numMinusClass0.f90"
#undef numMinusClassFUNCNAME
#undef DATATYPE

#define numMinusClassFUNCNAME numMinusClassc
#define DATATYPE complex*8
#include "templet/numMinusClass0.f90"
#undef numMinusClassFUNCNAME
#undef DATATYPE

#define numMinusClassFUNCNAME numMinusClassz
#define DATATYPE complex*16
#include "templet/numMinusClass0.f90"
#undef numMinusClassFUNCNAME
#undef DATATYPE


	subroutine numMinusClass(Res,A,B,length)
		class(*),intent(inout)::Res(:)
		class(*),intent(in)::A,B(:)
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in numMinusClass, length<0',-1)
			call error_stop
		end if
		select type(Res)
			type is (integer)
				call numMinusClassi(Res,A,B,length)
			type is (real(kind=4))
				call numMinusClasss(Res,A,B,length)
			type is (real(kind=8))
				call numMinusClassd(Res,A,B,length)
			type is (complex(kind=4))
				call numMinusClassc(Res,A,B,length)
			type is (complex(kind=8))
				call numMinusClassz(Res,A,B,length)
			class default
				call num_Minus_Class(Res,A,B,length)
		end select
	end subroutine


	subroutine default_num_Minus_Class(R,A,B,length)
		class(*),target,intent(inout)::R(:)
		class(*),target,intent(in)::A
		class(*),target,intent(in)::B(:)
		integer,intent(in)::length
		integer::i
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default num_Minus_Class is called   #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do i=1,length
			call Num_Minus_num(R(i),A,B(i))
		end do
	end subroutine


	

