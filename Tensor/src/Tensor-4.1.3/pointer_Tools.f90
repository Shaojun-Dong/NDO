module Pointer_Tools
	use Tools
	implicit none
	public


	interface Pointer1DFunc
		module procedure Point1D_All_Valuei
		module procedure Point1D_All_Values
		module procedure Point1D_All_Valued
		module procedure Point1D_All_Valuec
		module procedure Point1D_All_Valuez
		module procedure Point1D_All_Valuel
		module procedure Point1D_All_Valuea
		module procedure Point1D_All_Valuea2

		module procedure Point1D_Some_Valuei
		module procedure Point1D_Some_Values
		module procedure Point1D_Some_Valued
		module procedure Point1D_Some_Valuec
		module procedure Point1D_Some_Valuez
		module procedure Point1D_Some_Valuel
		module procedure Point1D_Some_Valuea
		module procedure Point1D_Some_Valuea2

		module procedure Point1D_a_Valuei
		module procedure Point1D_a_Values
		module procedure Point1D_a_Valued
		module procedure Point1D_a_Valuec
		module procedure Point1D_a_Valuez
		module procedure Point1D_a_Valuel
		module procedure Point1D_a_Valuea
		module procedure Point1D_a_Valuea2
	end interface

	interface Pointer2DFunc
		module procedure Point2D_All_Valuei
		module procedure Point2D_All_Values
		module procedure Point2D_All_Valued
		module procedure Point2D_All_Valuec
		module procedure Point2D_All_Valuez
		module procedure Point2D_All_Valuel
		module procedure Point2D_All_Valuea
		module procedure Point2D_All_Valuea2

		module procedure Point2D_Some_Valuei
		module procedure Point2D_Some_Values
		module procedure Point2D_Some_Valued
		module procedure Point2D_Some_Valuec
		module procedure Point2D_Some_Valuez
		module procedure Point2D_Some_Valuel
		module procedure Point2D_Some_Valuea
		module procedure Point2D_Some_Valuea2

		module procedure Point2D_a_Valuei
		module procedure Point2D_a_Values
		module procedure Point2D_a_Valued
		module procedure Point2D_a_Valuec
		module procedure Point2D_a_Valuez
		module procedure Point2D_a_Valuel
		module procedure Point2D_a_Valuea
		module procedure Point2D_a_Valuea2


		module procedure Point2D_All_Value2i
		module procedure Point2D_All_Value2s
		module procedure Point2D_All_Value2d
		module procedure Point2D_All_Value2c
		module procedure Point2D_All_Value2z
		module procedure Point2D_All_Value2l
		module procedure Point2D_All_Value2a
		module procedure Point2D_All_Value2a2

		module procedure Point2D_Some_Value2i
		module procedure Point2D_Some_Value2s
		module procedure Point2D_Some_Value2d
		module procedure Point2D_Some_Value2c
		module procedure Point2D_Some_Value2z
		module procedure Point2D_Some_Value2l
		module procedure Point2D_Some_Value2a
		module procedure Point2D_Some_Value2a2

		module procedure Point2D_a_Value2i
		module procedure Point2D_a_Value2s
		module procedure Point2D_a_Value2d
		module procedure Point2D_a_Value2c
		module procedure Point2D_a_Value2z
		module procedure Point2D_a_Value2l
		module procedure Point2D_a_Value2a
		module procedure Point2D_a_Value2a2
	end interface

	interface Pointer3DFunc
		module procedure Point3D_All_Valuei
		module procedure Point3D_All_Values
		module procedure Point3D_All_Valued
		module procedure Point3D_All_Valuec
		module procedure Point3D_All_Valuez
		module procedure Point3D_All_Valuel
		module procedure Point3D_All_Valuea
		module procedure Point3D_All_Valuea2

		module procedure Point3D_Some_Valuei
		module procedure Point3D_Some_Values
		module procedure Point3D_Some_Valued
		module procedure Point3D_Some_Valuec
		module procedure Point3D_Some_Valuez
		module procedure Point3D_Some_Valuel
		module procedure Point3D_Some_Valuea
		module procedure Point3D_Some_Valuea2

		module procedure Point3D_a_Valuei
		module procedure Point3D_a_Values
		module procedure Point3D_a_Valued
		module procedure Point3D_a_Valuec
		module procedure Point3D_a_Valuez
		module procedure Point3D_a_Valuel
		module procedure Point3D_a_Valuea
		module procedure Point3D_a_Valuea2

		module procedure Point3D_All_Value2i
		module procedure Point3D_All_Value2s
		module procedure Point3D_All_Value2d
		module procedure Point3D_All_Value2c
		module procedure Point3D_All_Value2z
		module procedure Point3D_All_Value2l
		module procedure Point3D_All_Value2a
		module procedure Point3D_All_Value2a2

		module procedure Point3D_Some_Value2i
		module procedure Point3D_Some_Value2s
		module procedure Point3D_Some_Value2d
		module procedure Point3D_Some_Value2c
		module procedure Point3D_Some_Value2z
		module procedure Point3D_Some_Value2l
		module procedure Point3D_Some_Value2a
		module procedure Point3D_Some_Value2a2

		module procedure Point3D_a_Value2i
		module procedure Point3D_a_Value2s
		module procedure Point3D_a_Value2d
		module procedure Point3D_a_Value2c
		module procedure Point3D_a_Value2z
		module procedure Point3D_a_Value2l
		module procedure Point3D_a_Value2a
		module procedure Point3D_a_Value2a2
	end interface

	interface Pointer4DFunc
		module procedure Point4D_All_Valuei
		module procedure Point4D_All_Values
		module procedure Point4D_All_Valued
		module procedure Point4D_All_Valuec
		module procedure Point4D_All_Valuez
		module procedure Point4D_All_Valuel
		module procedure Point4D_All_Valuea
		module procedure Point4D_All_Valuea2

		module procedure Point4D_Some_Valuei
		module procedure Point4D_Some_Values
		module procedure Point4D_Some_Valued
		module procedure Point4D_Some_Valuec
		module procedure Point4D_Some_Valuez
		module procedure Point4D_Some_Valuel
		module procedure Point4D_Some_Valuea
		module procedure Point4D_Some_Valuea2

		module procedure Point4D_a_Valuei
		module procedure Point4D_a_Values
		module procedure Point4D_a_Valued
		module procedure Point4D_a_Valuec
		module procedure Point4D_a_Valuez
		module procedure Point4D_a_Valuel
		module procedure Point4D_a_Valuea
		module procedure Point4D_a_Valuea2

		module procedure Point4D_All_Value2i
		module procedure Point4D_All_Value2s
		module procedure Point4D_All_Value2d
		module procedure Point4D_All_Value2c
		module procedure Point4D_All_Value2z
		module procedure Point4D_All_Value2l
		module procedure Point4D_All_Value2a
		module procedure Point4D_All_Value2a2	

		module procedure Point4D_Some_Value2i
		module procedure Point4D_Some_Value2s
		module procedure Point4D_Some_Value2d
		module procedure Point4D_Some_Value2c
		module procedure Point4D_Some_Value2z
		module procedure Point4D_Some_Value2l
		module procedure Point4D_Some_Value2a
		module procedure Point4D_Some_Value2a2

		module procedure Point4D_a_Value2i
		module procedure Point4D_a_Value2s
		module procedure Point4D_a_Value2d
		module procedure Point4D_a_Value2c
		module procedure Point4D_a_Value2z
		module procedure Point4D_a_Value2l
		module procedure Point4D_a_Value2a
		module procedure Point4D_a_Value2a2
	end interface

	interface ClassPointer1DFunc
		module procedure Point1D_All_Value
		module procedure Point1D_Some_Value
		module procedure Point1D_a_Value
	end interface
	
	interface ClassPointer2DFunc
		module procedure Point2D_All_Value
		module procedure Point2D_Some_Value
		module procedure Point2D_a_Value
		module procedure Point2D_All_Value2
		module procedure Point2D_Some_Value2
		module procedure Point2D_a_Value2
	end interface

	interface ClassPointer3DFunc
		module procedure Point3D_All_Value
		module procedure Point3D_Some_Value
		module procedure Point3D_a_Value
		module procedure Point3D_All_Value2
		module procedure Point3D_Some_Value2
		module procedure Point3D_a_Value2
	end interface

	interface ClassPointer4DFunc
		module procedure Point4D_All_Value
		module procedure Point4D_Some_Value
		module procedure Point4D_a_Value
		module procedure Point4D_All_Value2
		module procedure Point4D_Some_Value2
		module procedure Point4D_a_Value2
	end interface


	interface classSize
		module procedure classSize1
		module procedure classSize2
		module procedure classSize3
		module procedure classSize4
	end interface

contains

	integer function classSize1(A)result(classSize)
		class(*),intent(in)::A(:)
		select type(A)
			type is (character(len=*))
				classSize=ClassSize_char1(A)
			class default
				classSize=size(A)
		end select
		return
	end function
	integer function ClassSize_char1(A)result(Res)
		character(len=*),target,intent(in)::A(:)
		character(len=len(A)),pointer::inACha(:)
		inACha=>A
		Res=size(inACha)
		return
	end function
	integer function classSize2(A,ith)result(classSize)
		class(*),intent(in)::A(:,:)
		integer,optional,intent(in)::ith
		select type(A)
			type is (character(len=*))
				classSize=ClassSize_char2(A)
			class default
				classSize=size(A,ith)
		end select
		return
	end function
	integer function ClassSize_char2(A)result(Res)
		character(len=*),target,intent(in)::A(:,:)
		character(len=len(A)),pointer::inACha(:,:)
		inACha=>A
		Res=size(inACha)
		return
	end function
	integer function classSize3(A,ith)result(classSize)
		class(*),intent(in)::A(:,:,:)
		integer,optional,intent(in)::ith
		select type(A)
			type is (character(len=*))
				classSize=ClassSize_char3(A)
			class default
				classSize=size(A,ith)
		end select
		return
	end function
	integer function ClassSize_char3(A)result(Res)
		character(len=*),target,intent(in)::A(:,:,:)
		character(len=len(A)),pointer::inACha(:,:,:)
		inACha=>A
		Res=size(inACha)
		return
	end function
	integer function classSize4(A,ith)result(classSize)
		class(*),intent(in)::A(:,:,:,:)
		integer,optional,intent(in)::ith
		select type(A)
			type is (character(len=*))
				classSize=ClassSize_char4(A)
			class default
				classSize=size(A,ith)
		end select
		return
	end function
	integer function ClassSize_char4(A)result(Res)
		character(len=*),target,intent(in)::A(:,:,:,:)
		character(len=len(A)),pointer::inACha(:,:,:,:)
		inACha=>A
		Res=size(inACha)
		return
	end function


#include "ClassPointer/Point1D.f90"

#include "ClassPointer/Point2D.f90"

#include "ClassPointer/Point2DChar.f90"

#include "ClassPointer/Point3D.f90"

#include "ClassPointer/Point3DChar.f90"

#include "ClassPointer/Point4D.f90"

#include "ClassPointer/Point4DChar.f90"

end module
