module General_Optimization_Element
	implicit none

	public::OptimElement_structure
	type :: OptimElement_structure
	contains
	end type OptimElement_structure

	interface
		subroutine Minus_interface(Res,A,B)
			import :: OptimElement_structure
			class(OptimElement_structure),target::Res
			class(OptimElement_structure),target::A
			class(OptimElement_structure),target::B
		end subroutine Minus_interface
	end interface


	interface
		subroutine multiply_interface(Res,num)
			import :: OptimElement_structure
			class(OptimElement_structure),target::Res
			real*8::num
		end subroutine multiply_interface
	end interface


	interface
		subroutine assignment_interface(Res,A)
			import :: OptimElement_structure
			class(OptimElement_structure),target::Res
			class(OptimElement_structure),target::A
		end subroutine assignment_interface
	end interface

	interface
		subroutine A_minus_t_time_B_interface(Res,A,t,B)
			import :: OptimElement_structure
			class(OptimElement_structure),target::Res
			class(OptimElement_structure),target::A
			class(OptimElement_structure),target::B
			real*8::t
		end subroutine A_minus_t_time_B_interface
	end interface

	interface
		subroutine norm2_interface(A,Res)
			import :: OptimElement_structure
			class(OptimElement_structure),target::A
			real*8,target::Res
		end subroutine norm2_interface
	end interface


	interface
		subroutine PrintString_interface(w)
			class(*),intent(in)::w
		end subroutine PrintString_interface
	end interface

	interface
		subroutine allocateOptimElementArray_interface(Res,len)
			import :: OptimElement_structure
			class(OptimElement_structure),allocatable,target::Res(:)
			integer,intent(in)::len
		end subroutine allocateOptimElementArray_interface
	end interface


	interface
		subroutine Dot_interface(Res,A,B)
			import :: OptimElement_structure
			real*8::Res
			class(OptimElement_structure),target::A
			class(OptimElement_structure),target::B
		end subroutine Dot_interface
	end interface

	interface
		subroutine Error_interface()
		end subroutine Error_interface
	end interface

	interface
		subroutine print_interface(A)
			import :: OptimElement_structure
			class(OptimElement_structure),target::A
		end subroutine print_interface
	end interface
	interface
		subroutine zeroFunc_interface(A)
			import :: OptimElement_structure
			class(OptimElement_structure),target::A
		end subroutine zeroFunc_interface
	end interface

	interface
		subroutine YDot_interface(Res,A,B)
			import :: OptimElement_structure
			real*8::Res
			class(OptimElement_structure),target::A
			class(OptimElement_structure),target::B
		end subroutine YDot_interface
	end interface

	interface
		subroutine RandomElement_interface(Res)
			import :: OptimElement_structure
			class(OptimElement_structure),target::Res
		end subroutine RandomElement_interface
	end interface

	interface
		subroutine element_product_interface(Res,A,B)
			import :: OptimElement_structure
			class(OptimElement_structure),target::Res
			class(OptimElement_structure),target::A
			class(OptimElement_structure),target::B
		end subroutine element_product_interface
	end interface

	interface
		subroutine A_over_sqrt_B_interface(Res,A,B)
			import :: OptimElement_structure
			class(OptimElement_structure),target::Res
			class(OptimElement_structure),target::A
			class(OptimElement_structure),target::B
		end subroutine A_over_sqrt_B_interface
	end interface

	procedure(allocateOptimElementArray_interface),pointer::allocateOptimElementArray=>defalutallocateArray
	procedure(Minus_interface),pointer::MinusFunc=>defalutMinus
	procedure(assignment_interface),pointer::assignmentFunc=>defalutassignment
	procedure(multiply_interface),pointer::multiplyFunc=>defalutmultiply
	procedure(A_minus_t_time_B_interface),pointer::A_minus_t_time_B=>defalutA_minus_t_time_B
	!A_minus_t_time_B: A=A-(t*B) or Res=A-(t*B)
	procedure(norm2_interface),pointer::norm2Func=>defalutnorm2
	procedure(PrintString_interface),pointer::printFunc=>defalutPrintString
	procedure(Dot_interface),pointer::DotFunc=>defalutDot
	procedure(Error_interface),pointer::ErrorFunc=>defalutError_step
	procedure(print_interface),pointer::writeoutElement=>defalutwriteoutElement
	procedure(zeroFunc_interface),pointer::zeroFunc=>defalutzeroFunc
	procedure(YDot_interface),pointer::YDotFunc=>defalutYDot
	procedure(RandomElement_interface),pointer::RandomElement=>defalutRandomElement
	procedure(element_product_interface),pointer::element_product=>default_element_product
	procedure(A_over_sqrt_B_interface),pointer::A_over_sqrt_B=>default_A_over_sqrt_B



contains

	subroutine defalutzeroFunc(Res)
		class(OptimElement_structure),target::Res
		call printFunc('Do NOT set the zeroFunc yet in the module General_Optimization_Element')
		call ErrorFunc
	end subroutine

	subroutine defalutallocateArray(Res,len)
		class(OptimElement_structure),allocatable,target::Res(:)
		integer,intent(in)::len
		type(OptimElement_structure)::TMP
		call printFunc('Do NOT set the allocateArrayFunc yet in the module General_Optimization_Element')
		call ErrorFunc

		allocate(Res(len),mold=TMP)
		return
	end subroutine

	subroutine defalutPrintString(w)
		class(*),intent(in)::w
		select type(w)
			type is (character(len=*))
				write(*,*)trim(w)
			type is (real(kind=8))
				write(*,*)w
			type is (integer)
				write(*,*)w
		end select
	end subroutine

	subroutine defalutwriteoutElement(A)
		class(OptimElement_structure),allocatable,target::A
		call printFunc('Do NOT set the printFunc yet in the module General_Optimization_Element')
		call ErrorFunc
	end subroutine

	subroutine defalutError_step( )
		stop
	end subroutine

	subroutine defalutMinus(Res,A,B)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		call printFunc('Do NOT set the MinusFunc yet in the module General_Optimization_Element')
		call ErrorFunc
		return
	end subroutine

	subroutine defalutassignment(Res,A)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		call printFunc('Do NOT set the assignmentFunc yet in the module General_Optimization_Element')
		call ErrorFunc
		return
	end subroutine


	subroutine defalutmultiply(Res,num)
		class(OptimElement_structure),target::Res
		real*8::num
		call printFunc('Do NOT set the divisionFunc yet in the module General_Optimization_Element')
		call ErrorFunc
		return
	end subroutine

	subroutine defalutA_minus_t_time_B(Res,A,t,B)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		real*8::t
		call printFunc('Do NOT set the A_minus_t_time_BFunc yet in the module General_Optimization_Element')
		call ErrorFunc
		return
	end subroutine

	subroutine defalutnorm2(A,Res)
		real*8,target::Res
		class(OptimElement_structure),target::A
		call printFunc('Do NOT set the norm2Func yet in the module General_Optimization_Element')
		call ErrorFunc
		return
	end subroutine

	subroutine defalutDot(Res,A,B)
		real*8::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		call printFunc('Do NOT set the DotFunc yet in the module General_Optimization_Element')
		call ErrorFunc
	end subroutine

	subroutine defalutYDot(Res,A,B)
		real*8::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		call printFunc('Do NOT set the YDotFunc yet in the module General_Optimization_Element')
		call ErrorFunc
	end subroutine

	subroutine defalutRandomElement(Res)
		class(OptimElement_structure),target::Res
		call printFunc('Do NOT set the RandomElement yet in the module General_Optimization_Element')
		call ErrorFunc
	end subroutine


	subroutine default_element_product(Res,A,B)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		call printFunc('Do NOT set the element_product yet in the module General_Optimization_Element')
		call ErrorFunc
	end subroutine

	subroutine default_A_over_sqrt_B(Res,A,B)
		class(OptimElement_structure),target::Res
		class(OptimElement_structure),target::A
		class(OptimElement_structure),target::B
		call printFunc('Do NOT set the A_over_sqrt_B yet in the module General_Optimization_Element')
		call ErrorFunc
	end subroutine

	

end module
