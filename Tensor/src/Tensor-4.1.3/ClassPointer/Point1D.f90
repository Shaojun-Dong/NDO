#define Point1D_All_ValueFUNCNAME Point1D_All_Valuei
#define Point1D_Some_ValueFUNCNAME Point1D_Some_Valuei
#define Point1D_a_ValueFUNCNAME Point1D_a_Valuei
#define DATATYPE integer
#define DATATYPE2 integer
#include "templet/Point1D0.f90"
#undef Point1D_All_ValueFUNCNAME
#undef Point1D_Some_ValueFUNCNAME
#undef Point1D_a_ValueFUNCNAME
#undef DATATYPE
#undef DATATYPE2

#define Point1D_All_ValueFUNCNAME Point1D_All_Values
#define Point1D_Some_ValueFUNCNAME Point1D_Some_Values
#define Point1D_a_ValueFUNCNAME Point1D_a_Values
#define DATATYPE real(kind=4)
#define DATATYPE2 real(kind=4)
#include "templet/Point1D0.f90"
#undef Point1D_All_ValueFUNCNAME
#undef Point1D_Some_ValueFUNCNAME
#undef Point1D_a_ValueFUNCNAME
#undef DATATYPE
#undef DATATYPE2

#define Point1D_All_ValueFUNCNAME Point1D_All_Valued
#define Point1D_Some_ValueFUNCNAME Point1D_Some_Valued
#define Point1D_a_ValueFUNCNAME Point1D_a_Valued
#define DATATYPE real(kind=8)
#define DATATYPE2 real(kind=8)
#include "templet/Point1D0.f90"
#undef Point1D_All_ValueFUNCNAME
#undef Point1D_Some_ValueFUNCNAME
#undef Point1D_a_ValueFUNCNAME
#undef DATATYPE
#undef DATATYPE2

#define Point1D_All_ValueFUNCNAME Point1D_All_Valuec
#define Point1D_Some_ValueFUNCNAME Point1D_Some_Valuec
#define Point1D_a_ValueFUNCNAME Point1D_a_Valuec
#define DATATYPE complex(kind=4)
#define DATATYPE2 complex(kind=4)
#include "templet/Point1D0.f90"
#undef Point1D_All_ValueFUNCNAME
#undef Point1D_Some_ValueFUNCNAME
#undef Point1D_a_ValueFUNCNAME
#undef DATATYPE
#undef DATATYPE2

#define Point1D_All_ValueFUNCNAME Point1D_All_Valuez
#define Point1D_Some_ValueFUNCNAME Point1D_Some_Valuez
#define Point1D_a_ValueFUNCNAME Point1D_a_Valuez
#define DATATYPE complex(kind=8)
#define DATATYPE2 complex(kind=8)
#include "templet/Point1D0.f90"
#undef Point1D_All_ValueFUNCNAME
#undef Point1D_Some_ValueFUNCNAME
#undef Point1D_a_ValueFUNCNAME
#undef DATATYPE
#undef DATATYPE2


#define Point1D_All_ValueFUNCNAME Point1D_All_Valuel
#define Point1D_Some_ValueFUNCNAME Point1D_Some_Valuel
#define Point1D_a_ValueFUNCNAME Point1D_a_Valuel
#define DATATYPE logical
#define DATATYPE2 logical
#include "templet/Point1D0.f90"
#undef Point1D_All_ValueFUNCNAME
#undef Point1D_Some_ValueFUNCNAME
#undef Point1D_a_ValueFUNCNAME
#undef DATATYPE
#undef DATATYPE2

#define Point1D_All_ValueFUNCNAME Point1D_All_Valuea
#define Point1D_Some_ValueFUNCNAME Point1D_Some_Valuea
#define Point1D_a_ValueFUNCNAME Point1D_a_Valuea
#define DATATYPE character(len=characterlen)
#define DATATYPE2 character(len=*)
#include "templet/Point1D0.f90"
#undef Point1D_All_ValueFUNCNAME
#undef Point1D_Some_ValueFUNCNAME
#undef Point1D_a_ValueFUNCNAME
#undef DATATYPE
#undef DATATYPE2


	subroutine Point1D_All_Valuea2(length,inA,outp,chalen)
		integer,intent(in)::length,chalen
		class(*),target::inA(:)
		character(len=chalen),pointer::outp(:)
		select type(inA)
			type is (character(len=*))
				if(size(inA).lt.length)then
					call writemess('ERROR in pointing data,length',-1)
					call error_stop
				end if
				if(chalen.ne.len(inA))then
					call writemess('ERROR in pointing data,chalen',-1)
					call error_stop
				end if
				outp(1:length)=>inA(1:length)
			class default
				call writemess('ERROR type  in pointting',-1)
				call error_stop
		end select
		return
	end subroutine
	subroutine Point1D_Some_Valuea2(inA,ith,jth,outp,chalen)
		integer::ith,jth,chalen
		class(*),target::inA(:)
		character(len=chalen),pointer::outp(:)
		select type(inA)
			type is (character(len=*))
				if(size(inA).lt.jth)then
					call writemess('ERROR in pointing data,length',-1)
					call error_stop
				end if
				if(chalen.ne.len(inA))then
					call writemess('ERROR in pointing data,chalen',-1)
					call error_stop
				end if
				outp(1:(jth-ith+1))=>inA(ith:jth)
			class default
				call writemess('ERROR type in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine Point1D_a_Valuea2(inA,ith,outp,chalen)
		integer::ith,chalen
		class(*),target::inA(:)
		character(len=chalen),pointer::outp
		select type(inA)
			type is (character(len=*))
				if(size(inA).lt.ith)then
					call writemess('ERROR in pointing data,length',-1)
					call error_stop
				end if
				if(chalen.ne.len(inA))then
					call writemess('ERROR in pointing data,chalen',-1)
					call error_stop
				end if
				outp=>inA(ith)
			class default
				call writemess('ERROR type in pointting',-1)
				call error_stop
		end select
		return
	end subroutine

	!***************************************

	subroutine Point1D_All_Value(length,inA,outp)
		class(*),target::inA(:)
		class(*),pointer::outp(:)
		integer,intent(in)::length
		if(ClassSize(inA).lt.length)then
			call writemess('ERROR in pointing data,length',-1)
			call error_stop
		end if
		select type(inA)
			type is (character(len=*))
				call point1D_All_Value_char(length,inA,outp)
			class default
				outp(1:length)=>inA(1:length)
		end select
		return
	end subroutine
	subroutine point1D_All_Value_char(length,inA,outp)
		character(len=*),target::inA(:)
		class(*),pointer::outp(:)
		character(len=len(inA)),pointer::inACha(:)
		integer,intent(in)::length
		inACha=>inA
		outp=>inACha
		return
	end subroutine

	subroutine Point1D_Some_Value(inA,ith,jth,outp)
		integer::ith,jth
		class(*),target::inA(:)
		class(*),pointer::outp(:)
		if(ClassSize(inA).lt.jth)then
			call writemess('ERROR in pointing data,length',-1)
			call error_stop
		end if
		select type(inA)
			type is (character(len=*))
				call Point1D_Some_Value_char(inA,ith,jth,outp)
			class default
				outp(1:(jth-ith+1))=>inA(ith:jth)
		end select
		return
	end subroutine
	subroutine Point1D_Some_Value_char(inA,ith,jth,outp)
		integer::ith,jth
		character(len=*),target::inA(:)
		class(*),pointer::outp(:)
		character(len=len(inA)),pointer::inACha(:)
		inACha=>inA
		outp(1:(jth-ith+1))=>inACha(ith:jth)
		return
	end subroutine


	subroutine Point1D_a_Value(inA,ith,outp)
		integer::ith
		class(*),target::inA(:)
		class(*),pointer::outp
		if(ClassSize(inA).lt.ith)then
			call writemess('ERROR in pointing data,length',-1)
			call error_stop
		end if
		select type(inA)
			type is (character(len=*))
				call Point1D_a_Value_char(inA,ith,outp)
			class default
				outp=>inA(ith)
		end select
		return
	end subroutine

	subroutine Point1D_a_Value_char(inA,ith,outp)
		integer::ith
		character(len=*),target::inA(:)
		class(*),pointer::outp
		character(len=len(inA)),pointer::inACha(:)
		inACha=>inA
		outp=>inACha(ith)
		return
	end subroutine


	subroutine pointerERROR_test()
		character(len=20),target::chaarray(4)
		character(len=20),pointer::ap(:)
		class(*),pointer::clp(:),clp2(:)
		chaarray=['aaaaaa','bbbbbb','cccccc','ssssss']
		clp=>chaarray
		clp2=>clp(2:3)
		select type (clp2)
			type is (character(len=*))
				write(*,*)len(clp2),clp2
		end select
	end subroutine
