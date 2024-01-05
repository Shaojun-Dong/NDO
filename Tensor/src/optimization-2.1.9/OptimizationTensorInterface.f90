module Optimization_Tensor_interface
	use General_Optimization_Element
	use General_Optimization_Tools
	use Tensor_Tools
	use Tools
	implicit none
	private
	real*8::adam_smallest_delta=1d-8
	public::OptimElement_structure
	public::OptimTensor
	type, extends(OptimElement_structure) :: OptimTensor
		type(Tensor)::data
	contains
		procedure,public::pointer=>PointToData
	end type OptimTensor


	public::Tensor2OptimTensor,initialTensorInterFace

	public::assignment(=)
	interface assignment(=)
		module procedure OptimTensorToTensor
		module procedure TensorToOptimTensor
	end interface
contains

	subroutine TensorToOptimTensor(OT,T)
		type(OptimTensor),intent(inout)::OT
		type(Tensor),intent(in)::T
		OT%Data=T
		return
	end subroutine
	subroutine OptimTensorToTensor(T,OT)
		type(Tensor),intent(inout)::T
		type(OptimTensor),intent(in)::OT
		T=OT%Data
		return
	end subroutine

	subroutine PointToData(OT,T)
		class(OptimTensor),target::OT
		type(Tensor),pointer::T
		T=>OT%Data
		return
	end subroutine
	subroutine Tensor2OptimTensor(OT,T)
		class(OptimElement_structure),target::OT
		type(Tensor),pointer::T
		select type(OT)
		type is (OptimTensor)
			T=>OT%Data
		class default
			call writemess('ERROR in Tensor2OptimTensor',-1)
			call writemess('Input class is not OptimTensor',-1)
			call error_stop
		end select
		return
	end subroutine
	subroutine initialTensorInterFace()
		call writemess('Initial the Tensor InterFace for OptimElement_structure')
		allocateOptimElementArray=>defineallocateArray
		printFunc=>definePrintString
		MinusFunc=>defineMinus
		assignmentFunc=>defineassignment
		multiplyFunc=>definemultiply
		A_minus_t_time_B=>defineA_minus_t_time_B
		!A_minus_t_time_B:Res=A-(t*B) or A=A-(t*B) or 
		norm2Func=>definenorm2
		DotFunc=>defineDot
		YDotFunc=>defineYDot
		writeoutElement=>defineprint
		ZeroFunc=>defineZero
		ErrorFunc=>defineErrorFunc
		RandomElement=>defineRandom
		A_over_sqrt_B=>define_A_over_sqrt_B
		element_product=>define_element_product
	end subroutine
	subroutine defineallocateArray(Res,len)
		class(OptimElement_structure),allocatable,target::Res(:)
		integer,intent(in)::len
		type(OptimTensor)::TMP
		if(allocated(Res))deallocate(Res)
		allocate(Res(len),mold=TMP)
		return
	end subroutine

	subroutine definePrintString(w)
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

	subroutine defineprint(A)
		class(OptimElement_structure),allocatable,target::A
		type(Tensor),pointer::p
		call Tensor2OptimTensor(A,p)
		call p%print()
	end subroutine

	subroutine defineMinus(Res,A,B)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		type(Tensor),pointer::Ap,Bp,Rp
		call Tensor2OptimTensor(A,Ap)
		call Tensor2OptimTensor(B,Bp)
		call Tensor2OptimTensor(Res,Rp)
			Rp=Ap-Bp
		return
	end subroutine

	subroutine defineassignment(Res,A)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		type(Tensor),pointer::Ap,Rp
		call Tensor2OptimTensor(A,Ap)
		call Tensor2OptimTensor(Res,Rp)
			Rp=Ap
		return
	end subroutine

	subroutine definemultiply(Res,num)
		class(OptimElement_structure),target::Res
		real*8::num
		type(Tensor),pointer::Rp
		call Tensor2OptimTensor(Res,Rp)
		Rp=Rp*num
		return
	end subroutine

	subroutine defineA_minus_t_time_B(Res,A,t,B)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		real*8::t
		type(Tensor),pointer::Ap,Bp,Rp
		call Tensor2OptimTensor(A,Ap)
		call Tensor2OptimTensor(B,Bp)
		call Tensor2OptimTensor(Res,Rp)
		Rp=Ap-(Bp*t)
		return
	end subroutine

	subroutine definenorm2(A,Res)
		real*8,target::Res
		class(OptimElement_structure),target::A
		type(Tensor),pointer::Ap
		call Tensor2OptimTensor(A,Ap)
		Res=Ap%dnorm2()
		return
	end subroutine

	subroutine defineDot(Res,A,B)
		real*8::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		type(Tensor),pointer::Ap,Bp
		call Tensor2OptimTensor(A,Ap)
		call Tensor2OptimTensor(B,Bp)
		Res=(Ap.dot.Bp)
	end subroutine

	subroutine defineYDot(Res,A,B)
		real*8::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		type(Tensor),pointer::Ap,Bp
		call Tensor2OptimTensor(A,Ap)
		call Tensor2OptimTensor(B,Bp)
		Res=(Ap.dot.Bp)
	end subroutine

	subroutine defineZero(Res)
		class(OptimElement_structure),target::Res
		type(Tensor),pointer::Rp
		call Tensor2OptimTensor(Res,Rp)
		call Rp%zero()
	end subroutine

	subroutine defineErrorFunc()
		call error_stop
	end subroutine


	subroutine defineRandom(Res)
		class(OptimElement_structure),target::Res
		type(Tensor),pointer::Rp
		call Tensor2OptimTensor(Res,Rp)
		call Rp%random([-1d0,1d0])
	end subroutine
	
	subroutine define_element_product(Res,A,B)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		type(Tensor),pointer::Rp,Ap,Bp
			integer::i,totalData
		call Tensor2OptimTensor(Res,Rp)
		call Tensor2OptimTensor(A,Ap)
		call Tensor2OptimTensor(B,Bp)
		
		totalData=Ap%getTotalData()
		if(totalDAta.ne.Bp%getTotalData())then
			call writemess('ERROR in element_product',-1)
			call error_stop
		end if
		call Rp%allocate(Ap%dim(),Ap%getType())
		do i=1,TotalData
			call Rp%setValue(i,Ap%i(i)*Bp%i(i))
		end do
		return
	end subroutine

	subroutine define_A_over_sqrt_B(Res,A,B)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		type(Tensor),pointer::Rp,Ap,Bp
		character(len=characterlen)::w
		integer::i,totalData
		call Tensor2OptimTensor(Res,Rp)
		call Tensor2OptimTensor(A,Ap)
		call Tensor2OptimTensor(B,Bp)
		
		totalData=Ap%getTotalData()
		if(totalDAta.ne.Bp%getTotalData())then
			call writemess('ERROR in element_product',-1)
			call error_stop
		end if
		if(Ap%getType().ne.3)then
			call writemess('ERROR in adam, data type should be real*8')
			call error_Stop
		end if
		call Rp%allocate(Ap%dim(),Ap%getType())
		do i=1,TotalData
			if(dsqrt(Bp%di(i)).lt.adam_smallest_delta)then
				if((Ap%di(i).gt.adam_smallest_delta))then
					w='**** WARNING in adam: -A/srqt(B)='+(Ap%di(i)/dsqrt(Bp%di(i)))+', A='+Ap%di(i)
					w=w+' ,sqrt(B)='+sqrt(Bp%di(i))+'  ****'
					call writemess(w,-1)
				end if
				call Rp%setValue(i,0d0)
			else
				call Rp%setValue(i,Ap%di(i)/dsqrt(Bp%di(i)-adam_smallest_delta))
			end if
		end do
		return
	end subroutine


end module
