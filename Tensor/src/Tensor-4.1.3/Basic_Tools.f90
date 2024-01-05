module Basic_Tools
	use Tools
	use mpi
	use Pointer_Tools
	implicit none
	public
	
#include "Basic_Tools/Minus_interface.f90"

#include "Basic_Tools/Plus_interface.f90"

#include "Basic_Tools/Multiply_interface.f90"
	
#include "Basic_Tools/default_interface.f90"

	
	interface ModifyAllValue
		module procedure SetClassData
	end interface

	interface ClassWritemess
		module procedure ClassWritemess0
		module procedure ClassWritemess_form
		module procedure ClassWritemess_array_form
		module procedure ClassWritemess_array
	end interface


	interface CopyArray
		module procedure FastcopyARRAY
		module procedure copyARRAYSameType2D
		module procedure copyARRAYSameType3D
	end interface

contains

#include "Basic_Tools/CopyArray.f90" 
!A(:)=B(:)

#include "Basic_Tools/CopyData.f90"  
!A=B

#include "Basic_Tools/SetData.f90"  
!A(:)=B

!*****************

#include "Basic_Tools/Plus.f90"      
!R(:)=A(:)+B(:)

#include "Basic_Tools/PlusNum.f90"   
!R=A+B

#include "Basic_Tools/ClassPlusnum.f90" 
! R(:)=A(:)+B

#include "Basic_Tools/numPlusClass.f90" 
! R(:)=A+B(:)

!*****************

#include "Basic_Tools/Minus.f90" 
!R(:)=A(:)-B(:)

#include "Basic_Tools/MinusNum.f90" 
!R=A-B

#include "Basic_Tools/ClassMinusnum.f90" 
! R(:)=A(:)-B

#include "Basic_Tools/numMinusClass.f90" 
! R(:)=A-B(:)

!*****************

#include "Basic_Tools/numTimeNum.f90" 
!R=A*B

#include "Basic_Tools/Multiply.f90" 
!R(:)=A(:)*B
!R(:)=A*B(:)

#include "Basic_Tools/ArrayTimeArray.f90" 
!R(i)=A(i)*B(i)

!*****************

#include "Basic_Tools/numDivideNum.f90" 
!R=A/B

#include "Basic_Tools/divide.f90" 
!R(:)=A(:)/B

#include "Basic_Tools/ArrayDivideArray.f90" 
!R(i)=A(i)/B(i)

!*****************

#include "Basic_Tools/Modify.f90" 

#include "Basic_Tools/inout.f90" 

#include "Basic_Tools/MPI.f90" 

#include "Basic_Tools/ClassWritemess.f90" 

#include "Basic_Tools/getaValue.f90" 




	integer function select_Combine_type(classtype1,classtype2)result(Res)
		integer,intent(in)::classtype1,classtype2
		integer::flag
		if(classtype1.eq.classtype2)then
			Res=classtype1
			return
		end if
		if(classtype1.eq.7)then
			Res=7
			return
		end if
		if(classtype2.eq.7)then
			Res=7
			return
		end if
		flag=10*classtype1+classtype2
		select case(flag)
		!int+classtype ---> classtype
			case (12)!int,real4
				Res=2
			case (13)!int,real8
				Res=3
			case (14)!int,compelx(kind=4)
				Res=4
			case (15)!int,compelx(kind=8)
				Res=5
			case (17)!int,character
				Res=7
		!real4+classtype ---> max{2,classtype}
			case (21)!real(kind=4),int
				Res=2
			case (23)!real(kind=4),real8
				Res=3
			case (24)!real(kind=4),compelx(kind=4)
				Res=4
			case (25)!real(kind=4),compelx(kind=8)
				Res=5
			case (27)!real(kind=4),character
				Res=7
		!depend on 		classtype
			case (31)!real(kind=8),int
				Res=3
			case (32)!real(kind=8),real4
				Res=3
			case (34)!real(kind=8),compelx(kind=4)
				Res=5
			case (35)!real(kind=8),compelx(kind=8)
				Res=5
			case (37)!real(kind=8),character
				Res=7
		!depend on 		classtype		
			case (41)!compelx(kind=4),int
				Res=4
			case (42)!compelx(kind=4),real4
				Res=4
			case (43)!compelx(kind=4),real8
				Res=5
			case (45)!compelx(kind=4),compelx(kind=8)
				Res=5
			case (47)!compelx(kind=4),character 
				Res=7
		!complex*16+classtype ---> max{5,classtype}	
			case (51)!compelx(kind=8),int
				Res=5
			case (52)!compelx(kind=8),real4
				Res=5
			case (53)!compelx(kind=8),real8
				Res=5
			case (54)!compelx(kind=8),compelx(kind=4)
				Res=5
			case (57)!compelx(kind=8),character
				Res=7
		!character+classtype ---> max{7,classtype}		
			case (71)!character,int
				Res=7
			case (72)!character,real4
				Res=7
			case (73)!character),real8
				Res=7
			case (74)!character,real4
				Res=7
			case (75)!character),real8
				Res=7
			case (76)!character,logical 
				Res=7
			case default
				write(*,*)"ERROR, no such type in select"
				write(*,*)"flag=",flag
				call error_stop
		end select
		return
	end function
end module
