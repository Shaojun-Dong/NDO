#define getaDataValueFUNCNAME getValuei
#define getAllDataValueFUNCNAME getAllValuei
#define DATATYPE integer
#include "templet/getaValue0.f90"
#undef getaDataValueFUNCNAME
#undef getAllDataValueFUNCNAME
#undef DATATYPE

#define getaDataValueFUNCNAME getValues
#define getAllDataValueFUNCNAME getAllValues
#define DATATYPE real*4
#include "templet/getaValue0.f90"
#undef getaDataValueFUNCNAME
#undef getAllDataValueFUNCNAME
#undef DATATYPE

#define getaDataValueFUNCNAME getValued
#define getAllDataValueFUNCNAME getAllValued
#define DATATYPE real*8
#include "templet/getaValue0.f90"
#undef getaDataValueFUNCNAME
#undef getAllDataValueFUNCNAME
#undef DATATYPE

#define getaDataValueFUNCNAME getValuec
#define getAllDataValueFUNCNAME getAllValuec
#define DATATYPE complex*8
#include "templet/getaValue0.f90"
#undef getaDataValueFUNCNAME
#undef getAllDataValueFUNCNAME
#undef DATATYPE

#define getaDataValueFUNCNAME getValuez
#define getAllDataValueFUNCNAME getAllValuez
#define DATATYPE complex*16
#include "templet/getaValue0.f90"
#undef getaDataValueFUNCNAME
#undef getAllDataValueFUNCNAME
#undef DATATYPE

#define getaDataValueFUNCNAME getValuel
#define getAllDataValueFUNCNAME getAllValuel
#define DATATYPE logical
#include "templet/getaValue0.f90"
#undef getaDataValueFUNCNAME
#undef getAllDataValueFUNCNAME
#undef DATATYPE


	function getValuea(Da,ith)result(Res)
		character(len=:),allocatable::Res
		class(DataArray),intent(in)::Da
		integer,intent(in)::ith
		if(Da%getType().eq.8)then
			allocate(character(len=Da%DataCharacterLen)::Res)
		else
			allocate(character(len=characterlen)::Res)
		end if
		call getavalue(Res,Da%ClassData,ith)
		return
	end function
	function getAllValuea(Da)result(Res)
		character(len=:),allocatable::Res(:)
		class(DataArray),intent(in)::Da
		if(Da%getType().eq.8)then
			allocate(character(len=Da%DataCharacterLen)::Res(Da%TotalData))
		else
			allocate(character(len=characterlen)::Res(Da%TotalData))
		end if
		call FastCopyArray(Res,Da%ClassData,Da%TotalData)
		return
	end function

	subroutine GetClassValue(Da,outVal,ith)
		class(DataArray),intent(in)::Da
		integer,intent(in)::ith
		class(*),intent(inout)::outVal
		call getavalue(outVal,Da%ClassData,ith)
	end subroutine
	subroutine GetSomeClassValue(Da,outVal,ith,jth)
		class(DataArray),intent(in)::Da
		integer,intent(in)::ith,jth
		class(*),intent(inout)::outVal(:)
		class(*),pointer::p(:)
		call ClassPointer1DFunc(Da%ClassData,ith,jth,p)
		call FastCopyArray(outVal,p,jth-ith+1)
	end subroutine
	subroutine GetAllClassValue(Da,outVal)
		class(DataArray),intent(in)::Da
		class(*),intent(inout)::outVal(:)
		call FastCopyArray(outVal,Da%ClassData,Da%getTotalData())
	end subroutine
