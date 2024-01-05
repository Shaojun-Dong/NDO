module element_Tools
	use General_Optimization_Element
	use Tensor_Tools
	use Tools
	implicit none

	private
	real*8::adam_smallest_delta=1d-10
	public::OptData
	type, extends(OptimElement_structure) :: OptData
		type(Tensor),allocatable::Alldata(:)
	contains
		procedure,public::pointer=>PointerFunc
		procedure,public::storeData
		procedure,public::write=>writeData
		procedure,public::read=>readData
		procedure,public::readOldData
		procedure,public::initialData
		procedure,public::allocate=>allocatememory
	end type

	public::assignment(=)
	interface assignment(=)
		module procedure OptDataassignment0
	end interface

	public::initialElementInterFace,PointToOptData,permute,OptDatanorm2
contains
	subroutine allocatememory(OT)
		class(OptData),intent(inout)::OT
		if(allocated(OT%Alldata))deallocate(OT%Alldata)
		allocate(OT%Alldata(9))
	end subroutine
	subroutine initialData(OT,ParameterLen,phy_length,elementScal)
		class(OptData),intent(inout)::OT
		integer,intent(in)::ParameterLen,phy_length
		real*8,intent(in)::elementScal
		allocate(OT%Alldata(9))
		call OT%Alldata(1)%allocate([ParameterLen,phy_length],'real*8')
		call OT%Alldata(1)%setName(1,'A.i')
		call OT%Alldata(1)%setName(2,'A.j')
		call OT%Alldata(2)%allocate([ParameterLen,phy_length],'real*8')
		call OT%Alldata(2)%setName(1,'A.i')
		call OT%Alldata(2)%setName(2,'A.j')
		call OT%Alldata(3)%allocate([ParameterLen,phy_length],'real*8')
		call OT%Alldata(3)%setName(1,'A.i')
		call OT%Alldata(3)%setName(2,'A.j')
		call OT%Alldata(4)%allocate([ParameterLen,phy_length],'real*8')
		call OT%Alldata(4)%setName(1,'A.i')
		call OT%Alldata(4)%setName(2,'A.j')
		call OT%Alldata(5)%allocate([ParameterLen],'real*8')
		call OT%Alldata(5)%setName(1,'A.i')
		call OT%Alldata(6)%allocate([ParameterLen],'real*8')
		call OT%Alldata(6)%setName(1,'A.i')
		call OT%Alldata(7)%allocate([phy_length],'real*8')
		call OT%Alldata(7)%setName(1,'A.i')
		call OT%Alldata(8)%allocate([phy_length],'real*8')
		call OT%Alldata(8)%setName(1,'A.i')
		call OT%Alldata(9)%allocate([ParameterLen],'real*8')
		call OT%Alldata(9)%setName(1,'A.i')
		call random(OT,elementScal)
	end subroutine
	subroutine permute(A)
		type(OptData), intent(inout) :: A
			INTEGER::I
		if(.not.allocated(A%Alldata))then
			call writemess('ERROR in permute')
			call error_stop
		end if
		do i=1,9
			if(A%Alldata(i)%getRank().eq.2)then
				call A%Alldata(i)%permute(['A.i','A.j'])
			else if(A%Alldata(i)%getRank().gt.2)then
				call writemess('ERROR in permute')
				call error_stop
			end if
		end do
		return
	end subroutine
	subroutine random(OT,elementScal)
		type(OptData), intent(inout) :: OT
		real*8,intent(in)::elementScal
		integer::i
		if(.not.allocated(OT%Alldata))then
			call writemess('ERROR in random')
			call error_stop
		end if
		do i=1,9
			call OT%Alldata(i)%random([-elementScal,elementScal])
		end do
		return
	end subroutine

	subroutine writeData(OT,uni)
		class(OptData),intent(in)::OT
		integer,intent(in)::uni
		integer::i
		if(.not.allocated(OT%Alldata))then
			call writemess('ERROR ')
			call error_stop
		end if
		do i=1,9
			call OT%Alldata(i)%write(uni)
		end do
		return
	end subroutine
	subroutine readData(OT,uni)
		class(OptData),intent(inout)::OT
		integer,intent(in)::uni
		integer::i
		allocate(OT%Alldata(9))
			
		do i=1,9
			call OT%Alldata(i)%read(uni)
		end do
		return
	end subroutine
	subroutine readOldData(OT,uni)
		class(OptData),intent(inout)::OT
		integer,intent(in)::uni
		integer::i
		allocate(OT%Alldata(9))
		do i=1,9
			call OT%Alldata(i)%read(uni)
		end do
		call OT%Alldata(1)%setName(1,'A.i')
		call OT%Alldata(1)%setName(2,'A.j')
		call OT%Alldata(2)%setName(1,'A.i')
		call OT%Alldata(2)%setName(2,'A.j')
		call OT%Alldata(3)%setName(1,'A.i')
		call OT%Alldata(3)%setName(2,'A.j')
		call OT%Alldata(4)%setName(1,'A.i')
		call OT%Alldata(4)%setName(2,'A.j')
		call OT%Alldata(5)%setName(1,'A.i')
		call OT%Alldata(6)%setName(1,'A.i')
		call OT%Alldata(7)%setName(1,'A.i')
		call OT%Alldata(8)%setName(1,'A.i')
		call OT%Alldata(9)%setName(1,'A.i')
		return
	end subroutine


	subroutine OptDataassignment0(Rp,Ap)
		type(OptData),intent(inout)::Rp
		type(OptData),intent(in)::Ap
			integer::i
		if(.not.allocated(Rp%Alldata))then
			allocate(Rp%Alldata(9))
		end if
		if(.not.allocated(Ap%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		do i=1,9
			Rp%Alldata(i)=Ap%Alldata(i)
		end do
		return
	end subroutine

	subroutine PointerFunc(OT,W1,W2,U1,U2,C1,C2,B1,B2,D)
		class(OptData),target::OT
		type(Tensor),pointer::W1
		type(Tensor),pointer::W2
		type(Tensor),pointer::U1
		type(Tensor),pointer::U2
		type(Tensor),pointer::C1
		type(Tensor),pointer::C2
		type(Tensor),pointer::B1
		type(Tensor),pointer::B2
		type(Tensor),pointer::D
		if(.not.allocated(OT%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		W1=>OT%Alldata(1)
		W2=>OT%Alldata(2)
		U1=>OT%Alldata(3)
		U2=>OT%Alldata(4)
		C1=>OT%Alldata(5)
		C2=>OT%Alldata(6)
		B1=>OT%Alldata(7)
		B2=>OT%Alldata(8)
		D=>OT%Alldata(9)
		return
	end subroutine
	subroutine storeData(OT,inputData)
		class(OptData),target::OT
		type(Tensor)::inputData(9)
			integer::i
		if(.not.allocated(OT%Alldata))then
			allocate(OT%Alldata(9))
		end if
		do i=1,9
			OT%Alldata(i)=inputData(i)
		end do
		return
	end subroutine

	!******************************************
	!
	! Define the function used in OptimElement_structure
	!
	!******************************************

	subroutine PointToOptData(OT,T)
		class(OptimElement_structure),target::OT
		type(OptData),pointer::T
		select type(OT)
		type is (OptData)
			T=>OT
		class default
			call writemess('ERROR in PointToOptData',-1)
			call writemess('Input class is not OptData',-1)
			call error_stop
		end select
		return
	end subroutine

	subroutine initialElementInterFace()
		call writemess('Initial the OptData InterFace for OptimElement_structure')
		allocateOptimElementArray=>OptDataallocateArray
		MinusFunc=>OptDataMinus
		assignmentFunc=>OptDataassignment
		multiplyFunc=>OptDatamultiply
		A_minus_t_time_B=>OptDataA_minus_t_time_B
		!A_minus_t_time_B:Res=A-(t*B) or A=A-(t*B) or 
		norm2Func=>OptDatanorm2
		DotFunc=>OptDataDot
		ZeroFunc=>OptDataZero
		ErrorFunc=>OptDataErrorFunc
		RandomElement=>OptDataRandom
		printFunc=>OptDataPrintString
		A_over_sqrt_B=>OptData_A_over_sqrt_B
		element_product=>OptData_element_product
	end subroutine
	subroutine OptDataallocateArray(Res,len)
		class(OptimElement_structure),allocatable,target::Res(:)
		integer,intent(in)::len
		type(OptData)::TMP
		allocate(Res(len),mold=TMP)
		return
	end subroutine



	subroutine OptDataMinus(Res,A,B)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		type(OptData),pointer::Ap,Bp,Rp
		integer::i
		call PointToOptData(A,Ap)
		call PointToOptData(B,Bp)
		call PointToOptData(Res,Rp)
		if(.not.allocated(Rp%Alldata))then
			allocate(Rp%Alldata(9))
		end if
		if(.not.allocated(Ap%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		if(.not.allocated(Bp%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		call permute(Ap)
		call permute(Bp)
		do i=1,9
			Rp%Alldata(i)=Ap%Alldata(i)-Bp%Alldata(i)
		end do
		return
	end subroutine

	subroutine OptDataassignment(Res,A)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		type(OptData),pointer::Ap,Rp
		integer::i
		call PointToOptData(A,Ap)
		call PointToOptData(Res,Rp)
		if(.not.allocated(Rp%Alldata))then
			allocate(Rp%Alldata(9))
		end if
		if(.not.allocated(Ap%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		do i=1,9
			Rp%Alldata(i)=Ap%Alldata(i)
		end do
		return
	end subroutine

	subroutine OptDatamultiply(Res,num)
		class(OptimElement_structure),target::Res
		real*8::num
		type(OptData),pointer::Rp
		integer::i
		call PointToOptData(Res,Rp)
		if(.not.allocated(Rp%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		do i=1,9
			Rp%Alldata(i)=Rp%Alldata(i)*num
		end do
		return
	end subroutine

	subroutine OptDataA_minus_t_time_B(Res,A,t,B)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		real*8::t
		type(OptData),pointer::Ap,Bp,Rp
		integer::i
		call PointToOptData(A,Ap)
		call PointToOptData(B,Bp)
		call PointToOptData(Res,Rp)
		if(.not.allocated(Rp%Alldata))then
			allocate(Rp%Alldata(9))
		end if
		if(.not.allocated(Ap%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		if(.not.allocated(Bp%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		call permute(Ap)
		call permute(Bp)
		do i=1,9
			Rp%Alldata(i)=Ap%Alldata(i)-(Bp%Alldata(i)*t)
		end do
		return
	end subroutine

	subroutine OptDatanorm2(A,Res)
		real*8,target::Res
		class(OptimElement_structure),target::A
		type(OptData),pointer::Ap
		integer::i
		call PointToOptData(A,Ap)
		if(.not.allocated(Ap%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		Res=0d0
		do i=1,9
			Res=Res+Ap%Alldata(i)%dnorm2()
		end do
		return
	end subroutine

	subroutine OptDataDot(Res,A,B)
		real*8::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		type(OptData),pointer::Ap,Bp
		integer::i
		call PointToOptData(A,Ap)
		call PointToOptData(B,Bp)
		if(.not.allocated(Ap%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		if(.not.allocated(Bp%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		call permute(Ap)
		call permute(Bp)
		Res=0d0
		do i=1,9
			Res=Res+(Ap%Alldata(i).dot.Bp%Alldata(i))
		end do
	end subroutine


	subroutine OptDataZero(Res)
		class(OptimElement_structure),target::Res
		type(OptData),pointer::Rp
		integer::i
		call PointToOptData(Res,Rp)
		if(.not.allocated(Rp%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		do i=1,9
			call Rp%Alldata(i)%zero()
		end do
		return
	end subroutine

	subroutine OptDataErrorFunc()
		call error_stop
	end subroutine

	subroutine OptDataRandom(Res)
		class(OptimElement_structure),target::Res
		type(OptData),pointer::Rp
		integer::i
		call PointToOptData(Res,Rp)
		if(.not.allocated(Rp%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		do i=1,9
			call Rp%Alldata(i)%random([-1d0,1d0])
		end do
		return
	end subroutine
	subroutine OptDataPrintString(w)
		class(*),intent(in)::w
		select type(w)
			type is (character(len=*))
				call writemess(w)
			type is (real(kind=8))
				call writemess(w)
			type is (real(kind=4))
				call writemess(w)
			type is (complex(kind=8))
				call writemess(w)
			type is (complex(kind=4))
				call writemess(w)
			type is (integer)
				call writemess(w)
			type is (logical)
				call writemess(w)
			type is (Tensor)
				call writemess(w)
			class default
				call writemess('ERROR class in auxDataPrintString')
				call error_stop
		end select
	end subroutine

	subroutine OptData_element_product(Res,A,B)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		type(OptData),pointer::Rp,Ap,Bp
		integer::datai,i,totalData
		call PointToOptData(Res,Rp)
		call PointToOptData(A,Ap)
		call PointToOptData(B,Bp)
		if(.not.allocated(Rp%Alldata))then
			allocate(Rp%Alldata(9))
		end if
		if(.not.allocated(Ap%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		if(.not.allocated(Bp%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		call permute(Ap)
		call permute(Bp)
		do datai=1,9
			totalData=Ap%Alldata(datai)%getTotalData()
			if(totalDAta.ne.Bp%Alldata(datai)%getTotalData())then
				call writemess('ERROR in element_product',-1)
				call error_stop
			end if
			Rp%Alldata(datai)=Ap%Alldata(datai)
			do i=1,TotalData
				call Rp%Alldata(datai)%setValue(i,Ap%Alldata(datai)%i(i)*Bp%Alldata(datai)%i(i))
			end do
		end do

		
		return
	end subroutine

	subroutine OptData_A_over_sqrt_B(Res,A,B)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		type(OptData),pointer::Rp,Ap,Bp
		character(len=characterlen)::w
		integer::datai,i,totalData
		call PointToOptData(Res,Rp)
		call PointToOptData(A,Ap)
		call PointToOptData(B,Bp)
		if(.not.allocated(Rp%Alldata))then
			allocate(Rp%Alldata(9))
		end if
		if(.not.allocated(Ap%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		if(.not.allocated(Bp%Alldata))then
			call writemess('ERROR')
			call error_stop
		end if
		call permute(Ap)
		call permute(Bp)
		do datai=1,9
			totalData=Ap%Alldata(datai)%getTotalData()
			if(totalDAta.ne.Bp%Alldata(datai)%getTotalData())then
				call writemess('ERROR in element_product',-1)
				call error_stop
			end if
			if(Ap%Alldata(datai)%getType().ne.3)then
				call writemess('ERROR in adam, data type should be real*8')
				call error_Stop
			end if
			Rp%Alldata(datai)=Ap%Alldata(datai)
			do i=1,TotalData
				if(dsqrt(Bp%Alldata(datai)%di(i)).lt.adam_smallest_delta)then
					if((Ap%Alldata(datai)%di(i).gt.adam_smallest_delta))then
						w='**** WARNING in adam: -A/srqt(B)='+(Ap%Alldata(datai)%di(i)/dsqrt(Bp%Alldata(datai)%di(i)))
						w=w+', A='+Ap%Alldata(datai)%di(i)
						w=w+' ,sqrt(B)='+sqrt(Bp%Alldata(datai)%di(i))+'  ****'
						call writemess(w,-1)
					end if
					call Rp%Alldata(datai)%setValue(i,0d0)
				else
					call Rp%Alldata(datai)%setValue(i,Ap%Alldata(datai)%di(i)/dsqrt(Bp%Alldata(datai)%di(i)))
				end if
			end do
		end do
		return
	end subroutine

end module
