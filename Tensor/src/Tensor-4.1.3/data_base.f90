!
!                   _ooOoo_
!                  o8888888o
!                  88" . "88
!                  (| -_- |)
!                  O\  =  /O
!               ____/`---'\____
!             .'  \\|     |//  `.
!            /  \\|||  :  |||//  \
!           /  _||||| -:- |||||-  \
!           |   | \\\  -  /// |   |
!           | \_|  ''\---/''  |   |
!           \  .-\__  `-`  ___/-. /
!         ___`. .'  /--.--\  `. . __
!      ."" '<  `.___\_<|>_/___.'  >'"".
!     | | :  `- \`.;`\ _ /`;.`/ - ` : | |
!     \  \ `-.   \_ __\ /__ _/   .-` /  /
!======`-.____`-.___\_____/___.-`____.-'======
!                   `=---='
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!       Buddha blessed , no BUG 
! Report bugs of the package to sj.dong@outlook.com
!

! allocate(A,source=B):allocate A as the same shape and datatype of B, and copy B to A
! allocate(A,mold=B)  :allocate A as the same shape and datatype of B, But DO NOT copy B to A


module data_base_Tools
	use Pointer_Tools
	use Basic_Tools 
	use Tools
	use memory_type
	use mpi
	implicit none
	private
	type(memory),private::WorkingMemory


	
	public::DataArray
	type DataArray
		private
		class(*),public,pointer::ClassData(:)
		logical::DynamicClass=.true.
		integer::classType=0
						           !if classType=1, integer
						           !if classType=2, real(4)
						           !if classType=3, real(8)
						           !if classType=4, complex(4)
						           !if classType=5, complex(8)
						           !if classType=6, logical
						           !if classType=7, character(len=characterlen)
						           !if classType=8, character(len=:)
		integer::totalData=0
		integer::MomeryLength=0
		integer,public::DataCharacterLen=0
		!MomeryLength=size(classData)
		!            if MomeryLength = totalBlock, then the dataarray is full
		integer::totalBlock=0
		integer,allocatable::starti(:)
		integer,allocatable::endi(:)

		integer,allocatable:: iData(:)
		real(kind=4),allocatable::sData(:)
		real(kind=8),allocatable::dData(:)
		complex(kind=4),allocatable::cData(:)
		complex(kind=8),public,allocatable::zData(:)
		logical,allocatable:: ldata(:)
		character(len=characterlen),allocatable:: adata(:)
		character(len=:),allocatable::adata2(:)
		class(*),allocatable::clsData(:)

	contains
		procedure,public::empty=>emptyDataArray
		procedure,public::deallocate=>deallocateDataArray
		procedure::allocateDataArray1,allocateDataArray2,allocateDataArrayChar1,allocateDataArrayChar2,&
					allocateDataArray3,allocateDataArray4,allocateDataArray5,allocateDataArray6,allocateDataArrayChar4
		generic,public::allocate=>allocateDataArray1,allocateDataArray2,allocateDataArrayChar1,&
								allocateDataArrayChar2,allocateDataArray3,allocateDataArray4,&
								allocateDataArray5,allocateDataArray6,allocateDataArrayChar4
		procedure::allocateDataArray1Check,allocateDataArray2Check
		procedure::allocateDataArrayChar1Check,allocateDataArrayChar2Check
		procedure::allocateDataArray3Check,allocateDataArray4Check
		generic,public::allocateData=>allocateDataArray1Check,allocateDataArray2Check,&
								allocateDataArrayChar1Check,allocateDataArrayChar2Check,&
								allocateDataArray3Check,allocateDataArray4Check,&
								allocateDataArray5,allocateDataArray6
		procedure::allocateDataArrayMomerySource1
		generic,public::allocateMomeryClassType=>allocateDataArrayMomerySource1
		procedure::allocateDataArrayMomerySource1Check
		generic,public::allocateDataMomeryClassType=>allocateDataArrayMomerySource1Check
		

		procedure::allocate_From_source1,allocate_From_source2,allocate_From_source3
		generic,public::allocateClassType=>allocate_From_source1,allocate_From_source2,allocate_From_source3
		procedure::allocateDataArrayMomery2,allocateDataArrayMomery1,allocateDataArrayMomery3
		generic,public::allocateDataArrayMomery=>allocateDataArrayMomery2,allocateDataArrayMomery1,allocateDataArrayMomery3
		generic,public::allocateMomery=>allocateDataArrayMomery2,allocateDataArrayMomery1,allocateDataArrayMomery3
		procedure::allocateDataArrayMomery1Check,allocateDataArrayMomery2Check
		generic,public::allocateDataMomery=>allocateDataArrayMomery1Check,allocateDataArrayMomery2Check,&
												allocateDataArrayMomery3
		procedure::Set_block_momery1,Set_block_momery2
		generic,public::Set_block_momery=>Set_block_momery1,Set_block_momery2

		procedure,public::getType
		procedure,public::getClassType
		procedure,public::getTotalBlock
		procedure,public::getBlockIndex
		procedure::getFlagAll,getFlagith
		generic,public::getFlag=>getFlagAll,getFlagith
		procedure::getTotalDataith,getTotalDataAll
		generic,public::getTotalData=>getTotalDataith,getTotalDataAll
		procedure::getBlockLengthAll
		generic,public::getBlockTotalData=>getBlockLengthAll,getTotalDataith
		generic,public::getBlocklength=>getBlockLengthAll,getTotalDataith

		procedure::DataArrayithPointeri,DataArrayithPointers,DataArrayithPointerd,DataArrayithPointerc
		procedure::DataArrayithPointerz,DataArrayithPointerl,DataArrayithPointera
		procedure::DataArrayAllPointeri,DataArrayAllPointers,DataArrayAllPointerd,DataArrayAllPointerc
		procedure::DataArrayAllPointerz,DataArrayAllPointerl,DataArrayAllPointera
		generic,public::pointer=>DataArrayithPointeri,DataArrayithPointers,DataArrayithPointerd,DataArrayithPointerc,&
								DataArrayithPointerz,DataArrayithPointerl,DataArrayithPointera,&
								DataArrayAllPointeri,DataArrayAllPointers,DataArrayAllPointerd,DataArrayAllPointerc,&
								DataArrayAllPointerz,DataArrayAllPointerl,DataArrayAllPointera
		generic,public::pointAllData=>DataArrayAllPointeri,DataArrayAllPointers,DataArrayAllPointerd,DataArrayAllPointerc,&
								DataArrayAllPointerz,DataArrayAllPointerl,DataArrayAllPointera
		procedure::pointStarti1D,pointStarti2D,pointStarti3D,pointStarti4D
		procedure::pointEndi1D,pointEndi2D,pointEndi3D,pointEndi4D
		generic,public::pointStarti=>pointStarti1D,pointStarti2D,pointStarti3D,pointStarti4D
		generic,public::pointEndi=>pointEndi1D,pointEndi2D,pointEndi3D,pointEndi4D

		procedure::getValuei,getValues,getValued,getValuec,getValuez,getValuel,getValuea
		procedure::getAllValuei,getAllValues,getAllValued,getAllValuec,getAllValuez
		procedure::getAllValuel,getAllValuea
		generic,public::ii=>getValuei,getAllValuei
		generic,public::si=>getValues,getAllValues
		generic,public::di=>getValued,getAllValued
		generic,public::ci=>getValuec,getAllValuec
		generic,public::zi=>getValuez,getAllValuez
		generic,public::li=>getValuel,getAllValuel
		generic,public::ai=>getValuea,getAllValuea
		procedure::GetClassValue,GetSomeClassValue,GetAllClassValue
		generic,public::getValue=>GetClassValue,GetSomeClassValue,GetAllClassValue

		procedure::setValueith,setValueAll,setValueAll2,setSomeValue
		generic,public::setValue=>setValueith,setValueAll,setValueAll2,setSomeValue
		procedure::classTypepointAllData
		generic,public::classpointAllData=>classTypepointAllData
		procedure::classpointerith
		generic,public::classPointer=>classpointerith,classTypepointAllData
		procedure,public::setDynamic
		procedure,public::unsetDynamic
		procedure::setType1,setType2
		generic,public::setType=>setType1,setType2
		generic,public::setClassType=>setType1,setType2
		procedure,public::ifDynamicClass
		procedure,public::write=>writeExternalData
		procedure,public::read=>readExternalData
		procedure,public::print=>printData
		procedure,public::pointTotalData
		procedure,public::CutOffDataArray
		procedure,public::resetTotalData
		procedure,public::isnanData
		procedure,public::getCharacterLen
	end type


	public::assignment(=)
	interface assignment(=)
		module procedure assignmentDataArray
		module procedure assignmentDataArrayi
		module procedure assignmentDataArrays
		module procedure assignmentDataArrayd
		module procedure assignmentDataArrayc
		module procedure assignmentDataArrayz
		module procedure assignmentDataArrayl
		module procedure assignmentDataArraya
		module procedure assignmentDataArray2i
		module procedure assignmentDataArray2s
		module procedure assignmentDataArray2d
		module procedure assignmentDataArray2c
		module procedure assignmentDataArray2z
		module procedure assignmentDataArray2l
		module procedure assignmentDataArray2a
		module procedure assignmentDataArrayVali
		module procedure assignmentDataArrayVals
		module procedure assignmentDataArrayVald
		module procedure assignmentDataArrayValc
		module procedure assignmentDataArrayValz
		module procedure assignmentDataArrayVall
		module procedure assignmentDataArrayVala
		module procedure assignmentDataArrayVal2i
		module procedure assignmentDataArrayVal2s
		module procedure assignmentDataArrayVal2d
		module procedure assignmentDataArrayVal2c
		module procedure assignmentDataArrayVal2z
		module procedure assignmentDataArrayVal2l
		module procedure assignmentDataArrayVal2a
	end interface


	public::paste
	interface paste
		module procedure pasteTwoDataArray
	end interface

	public::pasteTwoDataArrayRoutine
	

	public::reorderBlock,check_input_order,send_DataArray,BCAST_DataArra
	public::MPI_SUM_DataArra,MPI_MAX_DataArray,MPI_MIN_DataArray

	public::set_allocate_Class_data_subroutine,select_data_type_char
	interface
		subroutine allocate_Class_data_interface(A,length,iType)
		class(*),intent(inout),allocatable::A(:)
		integer,intent(in)::length,iType
		end subroutine allocate_Class_data_interface
	end interface
	procedure(allocate_Class_data_interface),public,pointer::allocate_Class_data=>default_allocate_Class_data

contains

	!*********************************************************
	!
	! procedure pointer
	!
	!*********************************************************

	subroutine set_allocate_Class_data_subroutine(Func)
		procedure(allocate_Class_data_interface)::Func
		allocate_Class_data=>Func
	end subroutine

	subroutine default_allocate_Class_data(A,length,iType)
		class(*),intent(inout),allocatable::A(:)
		integer,intent(in)::length,iType
		call writemess('ERROR class type in allocate,iType='+iType,-1)
		call error_stop
	end subroutine

	!*********************************************************
	!
	! set value
	!
	!*********************************************************

#include "data_base_include/getaValue.f90"

	subroutine setValueith(Da,ith,inData)
		class(DataArray),intent(inout)::Da
		class(*),intent(in)::inData
		integer,intent(in)::ith
		if(.not.Da%getFlag())then
			call writemess('DO NOT allocate the DataArray yet',-1)
			call error_stop
		end if
		if(ith.gt.Da%getTotalData())then
			call writemess('ERROR in setting the ith element of the dataarray',-1)
			call writemess('ith='+ith,-1)
			call writemess('Da%getTotalData()='+Da%getTotalData(),-1)
			call error_stop
		end if 
		call ModifyaValue(Da%ClassData,inData,ith)
		return
	end subroutine
	subroutine setValueAll(Da,inData)
		class(DataArray),intent(inout)::Da
		class(*),intent(in)::inData
		if(.not.Da%getFlag())then
			call writemess('DO NOT allocate the DataArray yet',-1)
			call error_stop
		end if
		call ModifyAllValue(Da%ClassData,inData,Da%getTotalData())
		return
	end subroutine
	subroutine setSomeValue(Da,ith,inData)
		class(DataArray),intent(inout)::Da
		class(*),intent(in)::inData(:)
		integer,intent(in)::ith(2)
		if(.not.Da%getFlag())then
			call writemess('DO NOT allocate the DataArray yet',-1)
			call error_stop
		end if
		if(ith(2).gt.Da%getTotalData())then
			call writemess('ERROR in setting the ith element of the dataarray',-1)
			call writemess('ith(2)='+ith(2),-1)
			call writemess('Da%getTotalData()='+Da%getTotalData(),-1)
			call error_stop
		end if 
		if(ith(1).le.0)then
			call writemess('ERROR in setting the ith element of the dataarray',-1)
			call writemess('ith(1)='+ith(1),-1)
			call error_stop
		end if 
		if(ClassSize(inData).ne.(ith(2)-ith(1)+1))then
			call writemess('ERROR in setting the ith element of the dataarray',-1)
			call writemess('ClassSize(inData)='+ClassSize(inData),-1)
			call writemess('(ith(2)-ith(1)+1)='+(ith(2)-ith(1)+1),-1)
			call error_stop
		end if
		call ModifySomeValue(Da%ClassData,inData,ith(1),ith(2))
		return
	end subroutine
	subroutine setValueAll2(Da,inData)
		class(DataArray),intent(inout)::Da
		class(*),intent(in)::inData(:)
		if(.not.Da%getFlag())then
			call writemess('DO NOT allocate the DataArray yet',-1)
			call error_stop
		end if
		if(ClassSize(inData).ne.Da%getTotalData())then
			call writemess('ERROR in setting the ith element of the dataarray',-1)
			call writemess('ClassSize(inData)='+ClassSize(inData),-1)
			call writemess('Da%getTotalData()='+Da%getTotalData(),-1)
			call error_stop
		end if
		call fastCopyArray(Da%ClassData,inData,Da%getTotalData())
		return
	end subroutine


	!*********************************************************
	!
	! allocate
	!
	!*********************************************************

#include "data_base_include/allocate.f90"

#include "data_base_include/allocate_source.f90"

#include "data_base_include/momery.f90"	

	subroutine emptyDataArray(Da)
		class(DataArray),intent(inout)::Da
		integer::i
		if(deallocate_memory_flag)then
			call deallocateDataArray(Da)
			return
		end if
		if(allocated(Da%Endi))Da%Endi=0
		if(Da%DynamicClass)then
			Da%classType=0
			Da%DataCharacterLen=0
		end if
		Da%totalData=0
		Da%totalBlock=0
		Da%ClassData=>null()
		return
	end subroutine

	subroutine deallocateDataArray(Da)
		class(DataArray),intent(inout)::Da
		integer::i
		if(allocated(Da%ClsData))deallocate(Da%ClsData)
		if(allocated(Da%iData))deallocate(Da%iData)
		if(allocated(Da%sData))deallocate(Da%sData)
		if(allocated(Da%dData))deallocate(Da%dData)
		if(allocated(Da%cData))deallocate(Da%cData)
		if(allocated(Da%zData))deallocate(Da%zData)
		if(allocated(Da%lData))deallocate(Da%lData)
		if(allocated(Da%aData))deallocate(Da%aData)
		if(allocated(Da%aData2))deallocate(Da%aData2)
		Da%ClassData=>null()
		if(Da%DynamicClass)then
			Da%classType=0
		end if
		Da%DataCharacterLen=0
		Da%totalData=0
		Da%totalBlock=0
		Da%MomeryLength=0
		return
	end subroutine




	!***************************************************
	!
	!		assignmentDataArray
	!
	!****************************************************


#include "data_base_include/assignment.f90"

	subroutine assignmentDataArray(inoutDa,inDa)
		type(DataArray),intent(inout)::inoutDa
		type(DataArray),intent(in)::inDa
		call inoutDa%empty()
		if(.not.inDa%getFlag())return
			
		if((inoutDa%DynamicClass).or.(inoutDa%classType.eq.0))then
			inoutDa%classType=inDa%getType()
			inoutDa%DataCharacterLen=inDa%DataCharacterLen
		end if

		inoutDa%TotalData=inDa%getTotalData()
		call allocateDataArrayClassData(inoutDa,inoutDa%TotalData,inoutDa%classType,inDa%DataCharacterLen)
		call FastcopyArray(inoutDa%ClassData,inDa%ClassData,inoutDa%totalData)

		inoutDa%totalBlock=inDa%getTotalBlock()
		call allocateCheck(inoutDa%starti,inoutDa%totalBlock)
		call allocateCheck(inoutDa%endi,inoutDa%totalBlock)
		inoutDa%starti(1:inoutDa%totalBlock)=inDa%starti(1:inoutDa%totalBlock)
		inoutDa%endi(1:inoutDa%totalBlock)=inDa%endi(1:inoutDa%totalBlock)
		return
	end subroutine


	!*********************************************************
	!
	! pointer function for DataArray
	!
	!*********************************************************

#include "data_base_include/pointer.f90"


	subroutine pointTotalData(Da,p)
		class(DataArray),target,intent(in)::Da
		integer,pointer::p
		p=>Da%TotalData
		return
	end subroutine

	subroutine classTypepointAllData(Da,p)
		class(DataArray),target,intent(in)::Da
		class(*),pointer,intent(inout)::p(:)
		call ClassPointer1DFunc(Da%TotalData,Da%ClassData,p)
		return
	end subroutine
	subroutine classpointerith(Da,p,ith)
		class(DataArray),target,intent(in)::Da
		integer,intent(in)::ith
		class(*),pointer,intent(inout)::p(:)
		if(.not.Da%getFlag())then
			call writemess('The DataArray is empty',-1)
			call error_Stop
		end if
		if(Da%getFlag(ith))then
			call ClassPointer1DFunc(Da%ClassData,Da%starti(ith),Da%endi(ith),p)
		else
			call writemess('ERROR in pointting to the data of DataArray',-1)
			call writemess('the Block ith is empty',-1)
			call writemess('ith='+ith,-1)
			call error_stop
		end if
		return
	end subroutine

	!*********************************************************
	!
	! setting
	!
	!*********************************************************

	subroutine setDynamic(Da)
		class(DataArray),intent(inout)::Da
		Da%DynamicClass=.true.
		return
	end subroutine
	subroutine unsetDynamic(Da)
		class(DataArray),intent(inout)::Da
		Da%DynamicClass=.false.
		return
	end subroutine
	subroutine setType1(Da,classType,chalen)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::classType
		integer,optional,intent(in)::chalen
		if(Da%getFlag())then
			call writemess('CAN not set type to a non-empty data',-1)
			call error_stop
		end if
		Da%classType=classType
		if(Da%classType.eq.8)then
			if(present(chalen))then
				Da%DataCharacterLen=chalen
			else
				Da%DataCharacterLen=CharacterLen
			end if
		else if(Da%classType.eq.7)then
			Da%DataCharacterLen=CharacterLen
		else
			Da%DataCharacterLen=0
		end if
		Da%DynamicClass=.false.
		return
	end subroutine
	subroutine setType2(Da,classType,chalen)
		class(DataArray),intent(inout)::Da
		character(len=*),intent(in)::classType
		integer,optional,intent(in)::chalen
		if(Da%getFlag())then
			call writemess('CAN not set type to a non-empty data',-1)
			call error_stop
		end if
		Da%classType=select_data_type_char(ClassType)
		if(Da%classType.eq.8)then
			if(present(chalen))then
				Da%DataCharacterLen=chalen
			else
				Da%DataCharacterLen=CharacterLen
			end if
		else if(Da%classType.eq.7)then
			Da%DataCharacterLen=CharacterLen
		else
			Da%DataCharacterLen=0
		end if
		Da%DynamicClass=.false.
		return
	end subroutine
	logical function ifDynamicClass(Da)
		class(DataArray),intent(in)::Da
		ifDynamicClass=Da%DynamicClass
		return
	end function
	function getType(Da)
		integer::getType
		class(DataArray),intent(in)::Da
		getType=Da%classType
		return
	end function
	function getClassType(Da)
		character(len=50)::getClassType
		class(DataArray),intent(in)::Da
		getClassType=out_data_class_type(Da%classType)
		return
	end function
	function getTotalBlock(Da)
		integer::getTotalBlock
		class(DataArray),intent(in)::Da
		getTotalBlock=Da%totalBlock
		return
	end function
	function getFlagAll(Da)
		logical::getFlagAll
		class(DataArray),intent(in)::Da
		getFlagAll=Da%totalData.gt.0
		return
	end function
	function getFlagith(Da,ith)
		logical::getFlagith
		class(DataArray),intent(in)::Da
		integer,intent(in)::ith
		if(Da%totalData.le.0)then
			getFlagith=.false.
			return
		end if
		if(ith.gt.Da%getTotalBlock())then
			call writemess('ERROR in getFlag(ith),ith>totalBlock',-1)
			call writemess('ith='+ith)
			call writemess('totalBlock='+Da%getTotalBlock())
			call error_stop
		end if
		getFlagith=Da%endi(ith).gt.0
		return
	end function
	function getTotalDataAll(Da)
		integer::getTotalDataAll
		class(DataArray),intent(in)::Da
		getTotalDataAll=Da%totalData
		return
	end function
	function getTotalDataith(Da,ith)
		integer::getTotalDataith
		class(DataArray),intent(in)::Da
		integer,intent(in)::ith
		if(ith.gt.Da%getTotalBlock())then
			call writemess('ERROR in getTotalDataith(ith),ith>totalBlock',-1)
			call writemess('ith='+ith)
			call writemess('totalBlock='+Da%getTotalBlock())
			call error_stop
		end if
		if(Da%endi(ith).eq.0)then
			getTotalDataith=0
		else
			getTotalDataith=Da%endi(ith)-Da%starti(ith)+1
		end if
		return
	end function

	function getBlockLengthAll(Da)
		integer,allocatable::getBlockLengthAll(:)
		class(DataArray),intent(in)::Da
		integer::i,TotalBlock
		TotalBlock=Da%getTotalBlock()
		allocate(getBlockLengthAll(TotalBlock))
		do i=1,TotalBlock
			getBlockLengthAll(i)=getTotalDataith(Da,i)
		end do
		return
	end function

	function getBlockIndex(Da,ith)
		integer,allocatable::getBlockIndex(:)
		class(DataArray),intent(in)::Da
		integer,intent(in)::ith
		if(ith.gt.Da%getTotalBlock())then
			call writemess('ERROR in getBlockIndex(ith),ith>totalBlock',-1)
			call writemess('ith='+ith)
			call writemess('totalBlock='+Da%getTotalBlock())
			call error_stop
		end if
		allocate(getBlockIndex(2))
		getBlockIndex=[Da%starti(ith),Da%endi(ith)]
		return
	end function

	function out_data_class_type(classType)result(classtypechar)
		character(len=30)::classtypechar
		integer,intent(in) ::classType
		select case(classType)
			case (1)
				classtypechar='integer'
			case (2)
				classtypechar='real(kind=4)'
			case (3)
				classtypechar='real(kind=8)'
			case (4)
				classtypechar='complex(kind=4)'
			case (5)
				classtypechar='complex(kind=8)'
			case (6)
				classtypechar='logical'
			case (7)
				classtypechar='character'
			case (8)
				classtypechar='character(len=:)'
			case (0)
				classtypechar='null'
			case default 
				classtypechar='other'
		end 	select
		return
	end function
	function select_data_type_char(indata)result(select_data_type)
		integer::select_data_type
		character(len=*),intent(in) ::indata
		if(indata.equ.'integer') then
			select_data_type=1
			return
		end if
		if((indata.equ.'real*4').or.(indata.equ.'real(kind=4)').or.(indata.equ.'real')) then
			select_data_type=2
			return
		end if
		if((indata.equ.'real*8').or.(indata.equ.'real(kind=8)').or.(indata.equ.'double')) then
			select_data_type=3
			return
		end if
		if((indata.equ.'complex*8').or.(indata.equ.'complex(kind=4)').or.(indata.equ.'complex')) then
			select_data_type=4
			return
		end if
		if((indata.equ.'complex*16').or.(indata.equ.'complex(kind=8)')) then
			select_data_type=5
			return
		end if
		if(indata.equ.'logical') then
			select_data_type=6
			return
		end if
		if(indata.equ.'character') then
			select_data_type=7
			return
		end if
		if(indata.equ.'character(len=:)') then
			select_data_type=8
			return
		end if
		call writemess('ERROR type, type='+indata,-1)
		call error_stop()
		return
	end function

	logical function isnanData(A)
		class(DataArray),intent(in)::A
		integer::i,totaldata
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		totalData=A%totalData
		if(totalData.eq.0)then
			call writemess(" There is no data in Tensor, when checking if there is element is NAN",-1)
			call error_stop
		end if
		isnanData=.false.
		select case(A%classType)
			case(1)
				call A%pointer(ip)
				do i=1,totalData
					if(isnan(real(ip(i))) )then
						isnanData=.true.
						return
					end if
				end do
			case(2)
				call A%pointer(sp)
				do i=1,totalData
					if(isnan(sp(i)) )then
						isnanData=.true.
						return
					end if
				end do
			case(3)
				call A%pointer(dp)
				do i=1,totalData
					if(isnan(dp(i)) )then
						isnanData=.true.
						return
					end if
				end do
			case(4)
				call A%pointer(cp)
				do i=1,totalData
					if(isnan(real(cp(i),kind=4)).or.isnan(aimag(cp(i))) )then
						isnanData=.true.
						return
					end if
				end do
			case(5)
				call A%pointer(zp)
				do i=1,totalData
					if(isnan(dble(zp(i))).or.isnan(dimag(zp(i))) )then
						isnanData=.true.
						return
					end if
				end do
			case default
				call writemess(" ERROR in isnan",-1)
				call error_stop
		end select
		return
	end function

	function getCharacterLen(A)
		integer::getCharacterLen
		class(DataArray),intent(in)::A
		getCharacterLen=A%DataCharacterLen
		return
	end function

	!*********************************************************
	!
	! Paste two DataArrqy
	!
	!*********************************************************

	function pasteTwoDataArray(A,B)result(C)
		type(DataArray)::C
		type(DataArray),intent(in)::A,B
		integer::i,TotalDataA,TotalDataB,TotalDataC
		integer::classType,TotalBlockA,TotalBlockB,TotalBlockC
		class(*),pointer::TMP(:)
		TotalDataA=A%getTotalData()
		TotalBlockA=A%getTotalBlock()

		TotalDataB=B%getTotalData()
		TotalBlockB=B%getTotalBlock()

		TotalDataC=TotalDataA+TotalDataB
		TotalBlockC=TotalBlockA+TotalBlockB

		classType=select_Combine_type(A%getType(),B%getType())
		call C%setType(classType)
		C%DataCharacterLen=max(A%DataCharacterLen,B%DataCharacterLen)

		C%totalData=TotalDataC
		C%totalBlock=TotalBlockC
		call allocateDataArrayClassData(C,TotalDataC,C%classType,C%DataCharacterLen)
		call ClassPointer1DFunc(C%ClassData,1,TotalDataA,TMP)
		call FastCopyArray(TMP,A%ClassData,TotalDataA)
		call ClassPointer1DFunc(C%ClassData,TotalDataA+1,TotalDataC,TMP)
		call FastCopyArray(TMP,B%ClassData,TotalDataB)

		call allocateCheck(C%starti,TotalBlockC)
		C%starti(1:TotalBlockA)=A%starti(1:TotalBlockA)
		C%starti(TotalBlockA+1:TotalBlockC)=TotalDataA+B%starti(1:TotalBlockB)

		call allocateCheck(C%endi,TotalBlockC)
		C%endi(1:TotalBlockA)=A%endi(1:TotalBlockA)
		C%endi(TotalBlockA+1:TotalBlockC)=TotalDataA+B%endi(1:TotalBlockB)

		return
	end function

	subroutine pasteTwoDataArrayRoutine(C,A,B)
		type(DataArray),intent(inout)::C
		type(DataArray),intent(in)::A
		type(DataArray),intent(in)::B
		class(*),pointer::TMP(:)
		integer::i,TotalDataA,TotalDataB,TotalDataC
		integer::classType,TotalBlockA,TotalBlockB,TotalBlockC
		TotalDataA=A%getTotalData()
		TotalBlockA=A%getTotalBlock()

		TotalDataB=B%getTotalData()
		TotalBlockB=B%getTotalBlock()

		TotalDataC=TotalDataA+TotalDataB
		TotalBlockC=TotalBlockA+TotalBlockB

		
		if((C%DynamicClass).or.(C%classType.eq.0))then
			classType=select_Combine_type(A%getType(),B%getType())
			C%classType=classType
			C%DataCharacterLen=max(A%DataCharacterLen,B%DataCharacterLen)
		end if

		
		C%totalData=TotalDataC
		C%totalBlock=TotalBlockC
		call allocateDataArrayClassData(C,TotalDataC,C%classType,C%DataCharacterLen)
		call ClassPointer1DFunc(C%ClassData,1,TotalDataA,TMP)
		call FastCopyArray(TMP,A%ClassData,TotalDataA)
		call ClassPointer1DFunc(C%ClassData,TotalDataA+1,TotalDataC,TMP)
		call FastCopyArray(TMP,B%ClassData,TotalDataB)

		call allocateCheck(C%starti,TotalBlockC)
		C%starti(1:TotalBlockA)=A%starti(1:TotalBlockA)
		C%starti(TotalBlockA+1:TotalBlockC)=TotalDataA+B%starti(1:TotalBlockB)

		call allocateCheck(C%endi,TotalBlockC)
		C%endi(1:TotalBlockA)=A%endi(1:TotalBlockA)
		C%endi(TotalBlockA+1:TotalBlockC)=TotalDataA+B%endi(1:TotalBlockB)
		return
	end subroutine


	!*********************************************************
	!
	! write and read
	!
	!*********************************************************

	subroutine writeExternalData(A,uni)
		class(DataArray),intent(in)::A
		integer,intent(in)::uni
		integer::i
		write(uni,*)'readable_data',A%getFlag()
		if(.not.A%getFlag())then
			write(uni,*)'Empty DataArray'
			return
		end if
		write(uni,*)A%classType,A%totalData,A%totalBlock,A%DataCharacterLen
		write(uni,*)A%starti(1:A%totalBlock)
		write(uni,*)A%endi(1:A%totalBlock)
		call write_out_data(A%ClassData,A%totalData,uni)
		write(uni,*)'End_data'
		return
	end subroutine
	subroutine readExternalData(A,uni)
		class(DataArray),intent(inout)::A
		integer,intent(in)::uni
		character(len=50)::Noused
		integer::i
		logical::TMPLOgi
		read(uni,*)Noused,TMPLOgi
		if(.not.TMPLOgi)then
			read(uni,*)Noused
			return
		end if
		read(uni,*)A%classType,A%totalData,A%totalBlock,A%DataCharacterLen
		A%MomeryLength=A%TotalData
		call allocateCheck(A%starti,A%totalBlock)
		call allocateCheck(A%endi,A%totalBlock)
		read(uni,*)(A%starti(i),i=1,A%totalBlock)
		read(uni,*)(A%endi(i),i=1,A%totalBlock)
		call allocateDataArrayClassData(A,A%totalData,A%classType,A%DataCharacterLen)
		call read_externl_data(A%ClassData,A%totalData,uni)
		read(uni,*)Noused
	end subroutine
	subroutine printData(A_,PrintType)
		class(DataArray),intent(in)::A_
		integer,optional::PrintType
		class(*),pointer::clp(:)
		integer::i
		do i=1,A_%getTotalBlock()
			if(A_%getFlag(i))then
				call A_%ClassPointer(Clp,i)
				call ClassWritemess(Clp,PrintType)
			else
				call writemess('Empty Block')
			end if
		end do
		return
	end subroutine

	!*********************************************************
	!
	! permutation
	!
	!*********************************************************



	subroutine reorderBlock(inoutDa,newOrder)
		type(DataArray),intent(inout)::inoutDa
		integer,intent(in)::newOrder(:)
		integer,pointer::oldindex(:)
		integer::lenOrder,i
		lenOrder=inoutDa%getTotalBlock()
		if(size(NewOrder).ne.lenOrder)then
			call writemess('ERROR in reorderBlock,size(newOrder)',-1)
			call writemess('size(newOrder)='+size(newOrder),-1)
			call writemess('inoutDa%getTotalBlock()='+inoutDa%getTotalBlock(),-1)
			call error_stop
		end if
		call WorkingMemory%Check()
		call WorkingMemory%allocate(1,lenOrder)
		call WorkingMemory%get_memory(oldindex,lenOrder)
		oldindex=inoutDa%starti(1:lenOrder)
		do i=1,lenOrder
			inoutDa%starti(i)=oldindex(newOrder(i))
		end do

		oldindex=inoutDa%endi(1:lenOrder)
		do i=1,lenOrder
			inoutDa%endi(i)=oldindex(newOrder(i))
		end do

		call WorkingMemory%free()
		if(check_same_name_flag)call check_input_order(newOrder)
		return
	end subroutine

	subroutine check_input_order(order)
		integer,intent(in)::order(:)
		integer::check1,check2
		check1=sum(order)
		check2=(1+size(order))/2
		if(check1.ne.check2)then
			call writemess('ERROR in in put order',-1)
			call writemess(order,-1)
			call error_Stop
		end if
	end subroutine


	!*********************************************************
	!
	!  trucation
	!
	!*********************************************************

	subroutine CutOffDataArray(Da,BlockNum,outDa)
		class(DataArray),intent(in)::Da
		type(DataArray),intent(inout)::outDa
		integer,intent(in)::BlockNum(:)
		integer::NewTotalData,classType,i,TotalBlock,ii
		class(*),pointer::ip(:),Newip(:)
		TotalBlock=Da%getTotalBlock()
		if(size(BlockNum).ne.TotalBlock)then
			call writemess('ERROR in CutOffDataArray, TotalBlock',-1)
			call writemess('size(BlockNum)='+size(BlockNum),-1)
			call writemess('Da%getTotalData()='+Da%getTotalBlock(),-1)
			call error_stop
		end if
		NewTotalData=sum(BlockNum)
		
		if((outDa%DynamicClass).or.(outDa%classType.eq.0))then
			classType=Da%getType()
			outDa%classType=classType
			outDa%DataCharacterLen=Da%DataCharacterLen
		end if

		call allocateDataArray_no_empty_Block(outDa,BlockNum,outDa%classType)
		ii=0
		do i=1,TotalBlock
			if(BlockNum(i).gt.0)then
				ii=ii+1
				call Da%ClassPointer(ip,i)
				call outDa%ClassPointer(Newip,ii)
				call FastCopyArray(Newip,ip,BlockNum(i))
			end if
		end do

		return
	end subroutine


	!*********************************************************
	!
	!  code for MPI
	!
	!*********************************************************

	subroutine send_DataArray(Da1,Da2,ID1,ID2,MPIcommon)
		type(DataArray),intent(in)::Da1
		type(DataArray),intent(inout)::Da2
		integer,intent(in)::ID1,ID2
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,len1,len2,istatus(MPI_STATUS_SIZE),mpi_comm
		integer::sendDatai(5)
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		logical,pointer::lp(:)
		character(len=characterlen),pointer::ap(:)

		tag=1
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
		if(present(MPIcommon))then
			if((ID1.ge.proNum).or.(ID2.ge.proNum))return
		end if
		
		if(ID1.eq.ID2) return !The some cpu, do nothing
		
		if((proID.ne.ID1).and.(proID.ne.ID2)) return!The proID do not sent or recv, return

		if(proID.eq.ID1) then
			call mpi_send(Da1%DynamicClass,1,MPI_logical,ID2,tag,MPI_Comm,ierr)
		end if
		if(proID.eq.ID2) then
			call mpi_recv(Da2%DynamicClass,1,MPI_logical,ID1,tag,MPI_Comm,istatus,ierr)
		end if

		!**************************** integer **********************************************		
		if(proID.eq.ID1) then
			sendDatai=[Da1%classType,Da1%totalData,Da1%MomeryLength,Da1%totalBlock,Da1%DataCharacterLen]
			call mpi_send(sendDatai,4,MPI_integer,ID2,tag,MPI_Comm,ierr)
			if(Da1%totalData.le.0)return
		end if
		if(proID.eq.ID2) then
			call mpi_recv(sendDatai,4,MPI_integer,ID1,tag,MPI_Comm,istatus,ierr)
			Da2%classType=sendDatai(1)
			Da2%totalData=sendDatai(2)
			Da2%MomeryLength=sendDatai(3)
			Da2%totalBlock=sendDatai(4)
			Da2%DataCharacterLen=sendDatai(5)
			if(Da2%totalData.le.0)return
		end if

		!**************************** starti **********************************************		
		if(proID.eq.ID1) then
			call mpi_send(Da1%starti,Da1%totalBlock,MPI_integer,ID2,tag,MPI_Comm,ierr)
		end if
		if(proID.eq.ID2) then
			call allocateCheck(Da2%starti,Da2%totalBlock)
			call mpi_recv(Da2%starti,Da2%totalBlock,MPI_integer,ID1,tag,MPI_Comm,istatus,ierr)
		end if

		!**************************** endi **********************************************		
		if(proID.eq.ID1) then
			call mpi_send(Da1%endi,Da1%totalBlock,MPI_integer,ID2,tag,MPI_Comm,ierr)
		end if
		if(proID.eq.ID2) then
			call allocateCheck(Da2%endi,Da2%totalBlock)
			call mpi_recv(Da2%endi,Da2%totalBlock,MPI_integer,ID1,tag,MPI_Comm,istatus,ierr)
		end if

		if(proID.eq.ID1) then
			call ClassMPISend(Da1%ClassData,Da1%totalData,ID2,tag,MPI_Comm,ierr)
		end if
		if(proID.eq.ID2) then
			call allocateDataArrayClassData(Da2,Da2%totalData,Da2%classType,Da2%DataCharacterLen)
			call ClassMpiRecv(Da2%ClassData,Da2%totalData,ID1,tag,MPI_Comm,ierr)
		end if
		return
	end subroutine
	subroutine BCAST_DataArra(Ten1,ID,MPIcommon)
		type(DataArray),intent(inout)::Ten1
		integer,intent(in)::ID
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,len1,len2,istatus(MPI_STATUS_SIZE),mpi_comm
		integer::sendDatai(5)
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		logical,pointer::lp(:)
		character(len=characterlen),pointer::ap(:)
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		
		tag=1
		call mpi_comm_rank(mpi_comm,proID,ierr)
		call mpi_comm_size(mpi_comm,proNum,ierr )
		
		if(present(MPIcommon))then
			if(ID.ge.proNum)return
		end if
		
		call MPI_BCAST(Ten1%DynamicClass,1,MPI_logical,ID,mpi_comm,ierr)

		if(proID.eq.ID) then
			sendDatai=[Ten1%classType,Ten1%totalData,Ten1%MomeryLength,Ten1%totalBlock,Ten1%DataCharacterLen]
		end if
		!************************** integer *********************************************		
		call MPI_BCAST(sendDatai,4,MPI_integer,ID,mpi_comm,ierr)
		Ten1%classType=sendDatai(1)
		Ten1%totalData=sendDatai(2)
		Ten1%MomeryLength=sendDatai(3)
		Ten1%totalBlock=sendDatai(4)
		Ten1%DataCharacterLen=sendDatai(5)
		if(Ten1%totalData.le.0)return

		if(proID.ne.ID) then
			call allocateDataArrayClassData(Ten1,Ten1%totalData,Ten1%classType,Ten1%DataCharacterLen)
			call allocateCheck(Ten1%endi,Ten1%totalBlock)
			call allocateCheck(Ten1%starti,Ten1%totalBlock)
		end if
		call MPI_BCAST(Ten1%starti(1:Ten1%totalBlock),Ten1%totalBlock,MPI_integer,ID,mpi_comm,ierr)
		call MPI_BCAST(Ten1%endi(1:Ten1%totalBlock),Ten1%totalBlock,MPI_integer,ID,mpi_comm,ierr)
		call ClassMpiBcast(Ten1%ClassData,Ten1%totalData,ID,mpi_comm,ierr)
		return
	end subroutine
	subroutine MPI_SUM_DataArra(inTData,outTData,MPIcommon)!Do not check input data and do no allocate memery
		type(DataArray),intent(inout)::outTData
		type(DataArray)::inTData
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,len1,len2,istatus(MPI_STATUS_SIZE),mpi_comm
		integer::totalData
		integer,pointer::inip(:),outip(:)
		real*4,pointer::insp(:),outsp(:)
		real*8,pointer::indp(:),outdp(:)
		complex*8,pointer::incp(:),outcp(:)
		complex*16,pointer::inzp(:),outzp(:)
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		tag=1
		totalData=outTData%TotalData
		call ClassMpiAllreduce(inTData%ClassData,outTData%ClassData,inTData%totalData,MPI_SUM,mpi_comm,ierr)
		return
	end subroutine
	
	subroutine MPI_MAX_DataArray(inTData,outTData,MPIcommon)!Do not check input data and do no allocate memery
		type(DataArray),intent(inout)::outTData
		type(DataArray)::inTData
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,len1,len2,istatus(MPI_STATUS_SIZE),mpi_comm,totalData
		integer,pointer::inip(:),outip(:)
		real*4,pointer::insp(:),outsp(:)
		real*8,pointer::indp(:),outdp(:)
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		tag=1
		totalData=outTData%TotalData
		call ClassMpiAllreduce(inTData%ClassData,outTData%ClassData,inTData%totalData,MPI_MAX,mpi_comm,ierr)
		return
	end subroutine
	
	subroutine MPI_MIN_DataArray(inTData,outTData,MPIcommon)!Do not check input data and do no allocate memery
		type(DataArray),intent(inout)::outTData
		type(DataArray)::inTData
		integer::ierr
		integer,optional,intent(in)::MPIcommon
		integer::proID,proNum,tag,len1,len2,istatus(MPI_STATUS_SIZE),mpi_comm,totalData
		integer,pointer::inip(:),outip(:)
		real*4,pointer::insp(:),outsp(:)
		real*8,pointer::indp(:),outdp(:)
		if(present(MPIcommon))then
			mpi_comm=MPIcommon
		else
			mpi_comm=mpi_comm_world
		end if
		tag=1
		totalData=outTData%TotalData
		call ClassMpiAllreduce(inTData%ClassData,outTData%ClassData,inTData%totalData,MPI_MIN,mpi_comm,ierr)
		return
	end subroutine

end module
