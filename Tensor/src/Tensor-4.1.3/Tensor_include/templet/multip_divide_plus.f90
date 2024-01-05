	subroutine multiplyFuncName(A,B)
		class(Tensor),intent(inout)::A
		DATATYPE,intent(in)::B
		integer::classType
		integer,pointer::Aip(:),Rip(:)
		real*4,pointer::Asp(:),Rsp(:)
		real*8,pointer::Adp(:),Rdp(:)
		complex*8,pointer::Acp(:),Rcp(:)
		complex*16,pointer::Azp(:),Rzp(:)
		character(len=characterlen),pointer::Aap(:),Rap(:)
		classType=A%getType()
		select case(classType)
			case(1)
				call A%Data%pointAllData(Aip)
				Aip=Aip*B
			case(2)
				call A%Data%pointAllData(Asp)
				if(DATATYPENumber.eq.2)then
					call sscal (A%getTotalData(), B, Asp, 1)
				else
					call sscal (A%getTotalData(), real(B,kind=4), Asp, 1)
				end if
			case(3)
				call A%Data%pointAllData(Adp)
				if(DATATYPENumber.eq.3)then
					call dscal (A%getTotalData(), B, Adp, 1)
				else
					call dscal (A%getTotalData(), dble(B), Adp, 1)
				end if
			case(4)
				call A%Data%pointAllData(Acp)
				if(DATATYPENumber.eq.4)then
					call cscal (A%getTotalData(), B, Acp, 1)
				else
					call cscal (A%getTotalData(), cmplx(B,kind=4), Acp, 1)
				end if
			case(5)
				call A%Data%pointAllData(Azp)
				if(DATATYPENumber.eq.5)then
					call zscal (A%getTotalData(), B, Azp, 1)
				else
					call zscal (A%getTotalData(), dcmplx(B), Azp, 1)
				end if
			case default
				call writemess('ERROR type, datatpye='+classType,-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine divideFuncName(A,B)
		class(Tensor),intent(inout)::A
		DATATYPE,intent(in)::B
		integer::classType
		integer,pointer::Aip(:)
		real*4,pointer::Asp(:)
		real*8,pointer::Adp(:)
		complex*8,pointer::Acp(:)
		complex*16,pointer::Azp(:)
		character(len=characterlen),pointer::Aap(:)
		classType=A%getType()
		select case(classType)
			case(1)
				call A%Data%pointAllData(Aip)
				Aip=Aip*(1/B)
			case(2)
				call A%Data%pointAllData(Asp)
				call sscal (A%getTotalData(), real(1./B,kind=4), Asp, 1)
			case(3)
				call A%Data%pointAllData(Adp)
				call dscal (A%getTotalData(), dble(1d0/B), Adp, 1)
			case(4)
				call A%Data%pointAllData(Acp)
				call cscal (A%getTotalData(), cmplx(1./B,kind=4), Acp, 1)
			case(5)
				call A%Data%pointAllData(Azp)
				call zscal (A%getTotalData(), dcmplx(1d0/B), Azp, 1)
			case default
				call writemess('ERROR type, datatpye='+classType,-1)
				call error_stop
		end select
		return
	end subroutine

	subroutine TplusNumSubroutineFuncName(A,B,Res)
		class(Tensor),intent(in)::A
		DATATYPE,intent(in)::B
		type(Tensor),intent(inout)::Res
		integer::classType
		character(len=characterlen),pointer::Aap(:),Rap(:)
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		classType=select_type_in_add_minu_class_type(A%getType(),DATATYPENumber)
		call Res%Data%AllocateData(A%Data,classType)
		Res%Dimension=A%Dimension
		call ClassPlusnum(Res%Data%ClassData,A%Data%ClassDAta,B,A%getTotalData())
		return
	end subroutine
	subroutine TminusNumSubroutineFuncName(A,B,Res)
		class(Tensor),intent(in)::A
		DATATYPE,intent(in)::B
		type(Tensor),intent(inout)::Res
		integer::classType
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		classType=select_type_in_add_minu_class_type(A%getType(),DATATYPENumber)
		call Res%Data%AllocateData(A%Data,classType)
		Res%Dimension=A%Dimension
		call ClassMinusnum(Res%Data%ClassData,A%Data%ClassData,B,A%getTotalData())
		return
	end subroutine
	subroutine NumminusTSubroutineFuncName(A,B,Res)
		class(Tensor),intent(in)::A
		DATATYPE,intent(in)::B
		type(Tensor),intent(inout)::Res
		integer::classType
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		classType=select_type_in_add_minu_class_type(A%getType(),DATATYPENumber)
		call Res%Data%AllocateData(A%Data,classType)
		Res%Dimension=A%Dimension
		call numMinusClass(Res%Data%ClassData,B,A%Data%ClassData,A%getTotalData())
		return
	end subroutine
	function TplusNumFuncName(A,B)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		DATATYPE,intent(in)::B
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		call Res%Allocate(A,select_type_in_add_minu_class_type(A%getType(),DATATYPENumber))
		call ClassPlusnum(Res%Data%ClassData,A%Data%ClassData,B,A%getTotalData())
		return
	end function
	function NumplusTFuncName(B,A)result(Res)
		type(Tensor)::Res
		DATATYPE,intent(in)::B
		type(Tensor),intent(in)::A
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		call Res%Allocate(A,select_type_in_add_minu_class_type(A%getType(),DATATYPENumber))
		call numPlusClass(Res%Data%ClassData,B,A%Data%ClassData,A%getTotalData())
		return
	end function
	function TMinusNumFuncName(A,B)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		DATATYPE,intent(in)::B
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		call Res%Allocate(A,select_type_in_add_minu_class_type(A%getType(),DATATYPENumber))
		call ClassMinusnum(Res%Data%ClassData,A%Data%ClassData,B,A%getTotalData())
		return
	end function
	function NumMinusTFuncName(B,A)result(Res)
		type(Tensor)::Res
		DATATYPE,intent(in)::B
		type(Tensor),intent(in)::A
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		call Res%Allocate(A,select_type_in_add_minu_class_type(A%getType(),DATATYPENumber))
		call numMinusClass(Res%Data%ClassData,B,A%Data%ClassData,A%getTotalData())
		return
	end function

	function TmultiplyNumFuncName(A,B)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		DATATYPE,intent(in)::B
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		call Res%Allocate(A,select_type_in_add_minu_class_type(A%getType(),DATATYPENumber))
		call ClassTimenum(Res%Data%ClassData,A%Data%ClassData,B,A%getTotalData())
	end function

	function NumMultiplyTFuncName(B,A)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		DATATYPE,intent(in)::B
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		call Res%Allocate(A,select_type_in_add_minu_class_type(A%getType(),DATATYPENumber))
		call ClassTimenum(Res%Data%ClassData,A%Data%ClassData,B,A%getTotalData())
	end function

	function TdivideNumFuncName(A,B)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		DATATYPE,intent(in)::B
		if(.not.A%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		call Res%Allocate(A,select_type_in_add_minu_class_type(A%getType(),DATATYPENumber))
		call Classdividenum(Res%Data%ClassData,A%Data%ClassData,B,A%getTotalData())
	end function

	function NumdivideTFuncName(A,B)result(Res)
		type(Tensor)::Res
		DATATYPE,intent(in)::A
		type(Tensor),intent(in)::B
		real*4,pointer::Rsp(:)
		real*8,pointer::Rdp(:)
		complex*8,pointer::Rcp(:)
		complex*16,pointer::Rzp(:)
		class(*),pointer::Bp,Rp
		if(.not.B%Data%getFlag())then
			call writemess('ERROR: There is no data in the input tensor',-1)
			call error_stop
		end if
		if(B%getTotalData().ne.1)then
			call writemess('ERROR in num/Tensor, the tensor shold be regards as a number',-1)
			call error_stop
		end if
		call Res%Allocate(B,select_type_in_add_minu_class_type(B%getType(),DATATYPENumber))
		call ClassPointer1DFunc(Res%Data%ClassData,1,Rp)
		call ClassPointer1DFunc(B%Data%ClassData,1,Bp)
		call NumdivideNum(Rp,A,Bp)
		return
	end function