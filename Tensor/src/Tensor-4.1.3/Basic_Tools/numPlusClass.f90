#define numPlusClassFUNCNAME numPlusClassi
#define DATATYPE integer
#include "templet/numPlusClass0.f90"
#undef numPlusClassFUNCNAME
#undef DATATYPE

#define numPlusClassFUNCNAME numPlusClasss
#define DATATYPE real*4
#include "templet/numPlusClass0.f90"
#undef numPlusClassFUNCNAME
#undef DATATYPE

#define numPlusClassFUNCNAME numPlusClassd
#define DATATYPE real*8
#include "templet/numPlusClass0.f90"
#undef numPlusClassFUNCNAME
#undef DATATYPE

#define numPlusClassFUNCNAME numPlusClassc
#define DATATYPE complex*8
#include "templet/numPlusClass0.f90"
#undef numPlusClassFUNCNAME
#undef DATATYPE

#define numPlusClassFUNCNAME numPlusClassz
#define DATATYPE complex*16
#include "templet/numPlusClass0.f90"
#undef numPlusClassFUNCNAME
#undef DATATYPE



#define numPlusClassFUNCNAME numPlusClassa
#define DATATYPE character(len=*)
#include "templet/numPlusClass0.f90"
#undef numPlusClassFUNCNAME
#undef DATATYPE

	

	subroutine numPlusClass(Res,A,B,length)
		class(*),target,intent(inout)::Res(:)
		class(*),target,intent(in)::A
		class(*),target,intent(in)::B(:)
		integer,intent(in)::length
		if(length.eq.0)return
		if(length.lt.0)then
			call writemess('ERROR in numPlusClass, length<0',-1)
			call error_stop
		end if
		select type(Res)
			type is (integer)
				call numPlusClassi(Res,A,B,length)
			type is (real(kind=4))
				call numPlusClasss(Res,A,B,length)
			type is (real(kind=8))
				call numPlusClassd(Res,A,B,length)
			type is (complex(kind=4))
				call numPlusClassc(Res,A,B,length)
			type is (complex(kind=8))
				call numPlusClassz(Res,A,B,length)
			type is (character(len=*))
				call numPlusClassa(Res,A,B,length)
			class default
				call num_Plus_Class(Res,A,B,length)
		end select
	end subroutine

	subroutine default_num_Plus_Class(R,A,B,length)
		class(*),target,intent(inout)::R(:)
		class(*),target,intent(in)::A
		class(*),target,intent(in)::B(:)
		integer,intent(in)::length
		integer::i
		logical,save::FLAG=.true.
		if(FLAG)then
			call writemess('##############  WARNING  ####################')
			call writemess('#  the default num_Plus_Class is called   #')
			call writemess('#############################################')
			FLAG=.false.
		end if
		do i=1,length
			call Num_Plus_num(R(i),A,B(i))
		end do
	end subroutine


	

