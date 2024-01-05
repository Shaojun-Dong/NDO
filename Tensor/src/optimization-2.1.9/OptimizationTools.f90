module  Optimization_Tools
	use tensor_Tools
	use LinearSearchTools
	use LineSearchEnergy
	use Tools
	implicit none
	private
	integer,private,parameter::write_time_num=20
	integer,private,parameter::max_stop_counter=5
	!***************************************************
	!       abstract  definitation of OptimEngine_structure
	!***************************************************

	public::OptimEngine_structure
	type, abstract :: OptimEngine_structure
		logical,private::printFlag=.false.
		type(Tensor),allocatable::st(:)
		type(Tensor),allocatable::yt(:)
		type(Tensor),private,allocatable::workingMemory(:)
		type(Tensor),pointer::metric!In Newton method:
							 !   metric is the  Hessian matrix
							 !In Natural Gradient method:
							 !   metric is  the metric matrix
		real*8::SolveLinearEquationRCOND=-1d0
		integer::CGLLSmaxRunning=100
		real*8::CGLLSerror=1d-16

		integer::length=0
		logical::FullFlag=.false.
		integer::endindex=0
		real*8::Adam_rho1=0.9
		real*8::Adam_rho2=0.999
		real*8::Adam_rho1t=0.9
		real*8::Adam_rho2t=0.999
		real*8::Adam_delta=1d-8
		logical::Adam_warning=.false.
		integer::Adam_resetStep=-1
		integer::Adam_resetcounter=0

		integer::CG_direction_Flag=2
		real*8::max_step_in_Linear_search=-1d0
		real*8::first_step_in_Linear_search=0.1d0
		real*8::zero_gradient_in_RGM=-1!in RGM, gradient will regard as 0 when gradient< zero_gradient_in_RGM
		real*8::stop_gradient=-1 !it will stop if norm(gradient)<stop_gradient
		real*8::stop_step=1d-10 !it will stop if x<stop_step for max_stop_counter times, x is the search step in every step
		character(len=50)::method='null'
		real*8::penaltyfunctionFactor=10
		integer::solve_linear_equation_method=1
		procedure(LinearSearch_subroutine),pointer::LinearSearch=>LinearSearch1
		procedure(Metric_func_interface),pointer::Metric_func=>null()
		procedure(SJ_interface),pointer::SJ=>null()
		procedure(SJPointer_interface),pointer::SJPointer=>null()
		procedure(matrix_time_vector_interface),pointer::LLSMV=>null()
		procedure(matrix_matrix_time_vector_interface),pointer::LLSMMV=>null()
		procedure(Step1Subroutine_interface),pointer::inStep1=>null()
		procedure(Step2Subroutine_interface),pointer::inStep2=>null()
		procedure(Step3Subroutine_interface),pointer::inStep3=>null()
		procedure(Step4Subroutine_interface),pointer::inStep4=>null()
		procedure(BeforeStepSubroutine_interface),pointer::BeforeStep=>null()
		procedure(EndStepSubroutine_interface),pointer::EndStep=>null()
		procedure(RescalDirection_interface),pointer::RescalDirection=>DefaultRescalDirection
		procedure(AllocateMemory_interface),pointer::AllocateWorkingMemory=>DefaultAllocateWorkingMemory
		procedure(set_DefaultDirection_interface),pointer::set_Direction=>set_DefaultDirection
		procedure(Energy_FunctionInterface),pointer::Energy_Function=>defaultEnergy_Function
		procedure(GradientFunctionInterface),pointer::GradientFunction=>defaultGradientFunction
		procedure(InititalCGLLSStartPointInterface),pointer::InititalCGLLSStartPoint=>defaultInititalCGLLSStartPoint
	contains
		procedure(targetFunc), deferred :: target_Function
		procedure::Optimization1,Optimization2,Optimization3
		generic,public::Optim=>Optimization1,Optimization2,Optimization3

		procedure,public::setprintFlag=>LBFGSsetprintFlag
		procedure,public::NewElement=>pointNewElement
		procedure,public::i=>element_i
		procedure,public::allocate=>allocatememory
		procedure,public::deallocate=>deallocatememory
		procedure,public::getLength=>datalength
		procedure,public::resetEndpoint
		procedure,public::dataSize
		procedure::check_stop
		procedure,public::pointMemory

		procedure,public::set_CG_direction_flag
		procedure,public::set_stop_error
		procedure,public::set_stop_gradient=>set_stop_error
		procedure,public::set_stop_Step
		procedure::set_method1,set_method0
		generic,public::set_method=>set_method1,set_method0
		procedure,public::get_method
		procedure,public::set_penaltyFactor
		procedure,public::set_Hessian_func
		procedure,public::set_Metric_func
		procedure,public::set_SJ_func
		procedure,public::set_SJPointer_func
		procedure,public::set_CGLLS_parameter
		procedure,public::set_SJ_MV_func
		procedure,public::set_SJ_MMV_func
		

		procedure,public::set_max_step_in_Linear_search
		procedure,public::set_first_step_in_Linear_search
		procedure,public::set_linear_search_type
		procedure,public::set_linear_search_function
		procedure,public::set_zero_gradient_in_RGM
		generic,public::set_linear_search=>set_linear_search_type,set_linear_search_function
		procedure,public::Set_inStep1Func
		procedure,public::Set_inStep2Func
		procedure,public::Set_inStep3Func
		procedure,public::Set_inStep4Func
		procedure,public::Set_beforeStepFunc
		procedure,public::Set_EndStepFunc
		procedure,public::unSet_inStep1Func
		procedure,public::unSet_inStep2Func
		procedure,public::unSet_inStep3Func
		procedure,public::unSet_inStep4Func
		procedure,public::unSet_beforeStepFunc
		procedure,public::unSet_EndStepFunc
		procedure,public::Set_adam_reset_Step
		procedure,public::Set_adam_rhos
		procedure,public::Set_adam_delta
		procedure,public::set_RCOND_in_SolveLinearEquation
		procedure,public::set_solve_linear_equation_method
		procedure,public::allocateInternalMemory

		procedure,public::setEnergyFunction
		procedure,public::setGradientFunction
		procedure,public::setInititalCGLLSStartPoint
		procedure::runCGLLS2
		procedure::runCGLLS1
		procedure::runCGLLS_MMx
	end type OptimEngine_structure

	abstract interface
		subroutine targetFunc(A,outVal,outGradient,point)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: A
			real*8,target,intent(inout)::outVal
			type(Tensor),target,intent(inout)::outGradient
			type(Tensor),target,intent(in)::point
		end subroutine targetFunc
	end interface

	

	interface
		subroutine Energy_FunctionInterface(OE,outVal,point)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::outVal
			type(Tensor),target,intent(in)::point
		end subroutine Energy_FunctionInterface
	end interface
	interface
		subroutine GradientFunctionInterface(OE,Gradient,point)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			type(Tensor),target,intent(in)::point
			type(Tensor),target,intent(inout)::Gradient
		end subroutine GradientFunctionInterface
	end interface

	abstract interface
		subroutine LinearSearch_subroutine(OE,max_running,point,dir,x,outValue,inoutgra)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure), intent(inout) :: OE
			real*8,intent(inout)::x,outValue
			type(Tensor),intent(inout)::point
			type(Tensor),intent(inout)::inoutgra
			type(Tensor),intent(in)::dir
			integer,intent(in)::max_running
		end subroutine LinearSearch_subroutine
	end interface

	abstract interface
		subroutine Step1Subroutine_interface(OE,Value,direction,Gradient,point,ith,t)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Value
			type(Tensor),target,intent(inout)::Gradient
			type(Tensor),target,intent(inout)::point
			type(Tensor),target,intent(inout)::direction
			integer,target,intent(in)::ith
			real*8,target::t
		end subroutine Step1Subroutine_interface
	end interface
	abstract interface
		subroutine Step2Subroutine_interface(OE,Value,direction,Gradient,point,ith,t)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Value
			type(Tensor),target,intent(inout)::Gradient
			type(Tensor),target,intent(inout)::point
			type(Tensor),target,intent(inout)::direction
			integer,target,intent(in)::ith
			real*8,target::t
		end subroutine Step2Subroutine_interface
	end interface
	abstract interface
		subroutine Step3Subroutine_interface(OE,Value,direction,Gradient,point,ith,t)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Value
			type(Tensor),target,intent(inout)::Gradient
			type(Tensor),target,intent(inout)::point
			type(Tensor),target,intent(inout)::direction
			integer,target,intent(in)::ith
			real*8,target::t
		end subroutine Step3Subroutine_interface
	end interface

	abstract interface
		subroutine Step4Subroutine_interface(OE,Value,direction,Gradient,point,ith,t)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Value
			type(Tensor),target,intent(inout)::Gradient
			type(Tensor),target,intent(inout)::point
			type(Tensor),target,intent(inout)::direction
			integer,target,intent(in)::ith
			real*8,target::t
		end subroutine Step4Subroutine_interface
	end interface

	abstract interface
		subroutine BeforeStepSubroutine_interface(OE,Value,Gradient,point)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Value
			type(Tensor),target,intent(inout)::Gradient
			type(Tensor),target,intent(inout)::point
		end subroutine BeforeStepSubroutine_interface
	end interface

	abstract interface
		subroutine EndStepSubroutine_interface(OE,Value,Gradient,point)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Value
			type(Tensor),target,intent(inout)::Gradient
			type(Tensor),target,intent(inout)::point
		end subroutine EndStepSubroutine_interface
	end interface


	abstract interface
		subroutine RescalDirection_interface(OE,Dir)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure), intent(inout) :: OE
			type(Tensor),intent(inout)::Dir
		end subroutine RescalDirection_interface
	end interface

	abstract interface
		subroutine AllocateMemory_interface(OE,Gradient,direction)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure), intent(inout) :: OE
			type(Tensor),pointer,intent(inout)::Gradient,direction
		end subroutine AllocateMemory_interface
	end interface

	abstract interface
		subroutine set_DefaultDirection_interface(OE,direction,Gradient,point)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure), intent(inout) :: OE
			type(Tensor),intent(in)::point,Gradient
			type(Tensor),intent(inout)::direction
		end subroutine set_DefaultDirection_interface
	end interface


	abstract interface
		subroutine Metric_func_interface(OE,Hessian,direction,Gradient,point)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			type(Tensor),target,intent(in)::point,Gradient
			type(Tensor),target,intent(in)::direction
			type(Tensor),target,intent(inout)::Hessian
		end subroutine Metric_func_interface
	end interface

	abstract interface
		subroutine SJ_interface(OE,Y,G,direction,Gradient,point)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			type(Tensor),target,intent(in)::point,Gradient
			type(Tensor),target,intent(in)::direction
			type(Tensor),target,intent(inout)::Y,G
		end subroutine SJ_interface
	end interface

	abstract interface
		subroutine SJPointer_interface(OE,Yp,Gp,direction,Gradient,point)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			type(Tensor),target,intent(in)::point,Gradient
			type(Tensor),target,intent(in)::direction
			type(Tensor),pointer,intent(inout)::Yp,Gp
		end subroutine SJPointer_interface
	end interface



	interface
		subroutine matrix_time_vector_interface(OE,outAp,inp)!outAp=A*inp
			import :: Tensor
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			type(Tensor),target,intent(inout)::outAp
			type(Tensor),target,intent(in)::inp
		end subroutine matrix_time_vector_interface
	end interface
	interface
		subroutine matrix_matrix_time_vector_interface(OE,outAp,inp,MVFlag)!if MVFlag=true, outAp=A^T*A*inp, else outAp=A^T*inp
			import :: Tensor
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			type(Tensor),target,intent(inout)::outAp
			type(Tensor),target,intent(in)::inp
			logical,intent(in)::MVFlag
		end subroutine matrix_matrix_time_vector_interface
	end interface

	!***************************************************
	!        definitation of OptimRunner
	!***************************************************

	public::OptimRunner
	type, extends(OptimEngine_structure) :: OptimRunner
		procedure(external_target_Function_interface),pointer,NOPASS,private::externalFunc=>null()
	contains
		procedure::target_Function=>LBFGSFunc
		procedure,public::set_target_function
	end type

	interface
		subroutine external_target_Function_interface(outVal,outGradient,point)
			use tensor_Tools
			real*8,target,intent(inout)::outVal
			type(Tensor),target,intent(inout)::outGradient
			type(Tensor),target,intent(in)::point
		end subroutine external_target_Function_interface
	end interface


	interface
		subroutine InititalCGLLSStartPointInterface(OE,inoutx,b)
			use tensor_Tools
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			type(Tensor),target, intent(inout) :: inoutx
			type(Tensor),target, intent(in)::b
		end subroutine InititalCGLLSStartPointInterface
	end interface

contains


	subroutine defaultEnergy_Function(OE,outVal,point)
		class(OptimEngine_structure),target, intent(inout) :: OE
		real*8,target,intent(inout)::outVal
		type(Tensor),target,intent(in)::point
		type(Tensor),target::outGradient
		call OE%target_Function(outVal,outGradient,point)
		return
	end subroutine

	subroutine defaultGradientFunction(OE,outGradient,point)
		class(OptimEngine_structure),target, intent(inout) :: OE
		type(Tensor),target,intent(in)::point
		type(Tensor),target,intent(inout)::outGradient
		real*8,target::outVal
		call OE%target_Function(outVal,outGradient,point)
		return
	end subroutine

	subroutine setEnergyFunction(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(Energy_FunctionInterface)::Func
		call writemess('Set the Energy_Function as external function')
		call writemess('Set LinearSearch type as type4  in OptimEngine_structure') 
		call writemess('Use only the Energy in the Searching') 
		OE%LinearSearch=>LinearSearch4
		OE%Energy_Function=>Func
		return
	end subroutine

	subroutine setGradientFunction(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(GradientFunctionInterface)::Func
		call writemess('Set the Energy_Function as external function')
		OE%GradientFunction=>Func
		return
	end subroutine

	!***************************************************
	!           Basic function for OptimEngine_structure
	!***************************************************

	subroutine LBFGSpointerFunc(p,PointTarget,ith)
		type(Tensor),pointer,intent(inout)::p
		type(Tensor),target,intent(in)::PointTarget(:)
		integer::ith
		p=>PointTarget(ith)
		return
	end subroutine

	subroutine set_RCOND_in_SolveLinearEquation(OE,RCOND)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::RCOND
		OE%SolveLinearEquationRCOND=RCOND
		return
	end subroutine

	subroutine set_solve_linear_equation_method(OE,Flag)
		class(OptimEngine_structure), intent(inout) :: OE
		integer,intent(in)::Flag
		call writemess('Set the solve_linear_equation_method as solve_linear_equation_method='+Flag)
		OE%solve_linear_equation_method=Flag
	end subroutine

	

	subroutine point(BTool,st,yt,ith)
		class(OptimEngine_structure),intent(in)::BTool
		type(Tensor),pointer,intent(inout)::st,yt
		integer,intent(in)::ith
		if(ith.gt.BTool%length)then
			call writemess('ERROR in point to LBFGSTool',-1)
			call writemess('ith='+ith,-1)
			call writemess('BTool%length='+BTool%length,-1)
			call error_stop
		end if
		call LBFGSpointerFunc(st,BTool%st,ith)
		call LBFGSpointerFunc(yt,BTool%yt,ith)
		return
	end subroutine

	subroutine pointNewElement(BTool,st,yt)
		class(OptimEngine_structure),intent(inout)::BTool
		type(Tensor),pointer,intent(inout)::st,yt
		BTool%endindex=BTool%endindex+1
		if(BTool%endindex.gt.BTool%length)then
			BTool%endindex=1
			BTool%FullFlag=.true.
		end if
		call point(BTool,st,yt,BTool%endindex)
		return
	end subroutine

	subroutine element_i(BTool,st,yt,ith)
		class(OptimEngine_structure),intent(in)::BTool
		type(Tensor),pointer,intent(inout)::st,yt
		integer,intent(in)::ith
		integer::i
		if(ith.gt.BTool%length)then
			call writemess('ERROR in element_i,1',-1)
			call error_stop
		end if
		if(BTool%FullFlag)then
			i=BTool%endindex+ith
			if(i.gt.BTool%length)then
				i=i-BTool%length
			end if
		else
			i=ith
		end if
		call point(BTool,st,yt,i)
		return
	end subroutine

	subroutine allocatememory(BTool,length)
		class(OptimEngine_structure),intent(inout)::BTool
		integer,intent(in)::length
		if(length.le.0)then
			call writemess('ERROR in allocatememory for LBFGSTool',-1)
			call error_stop
		end if
		allocate(BTool%st(length))
		allocate(BTool%yt(length))
		BTool%length=length
		BTool%FullFlag=.false.
		BTool%endindex=0
		return
	end subroutine

	subroutine allocateInternalMemory(A,length)
		class(OptimEngine_structure),intent(inout)::A
		integer,intent(in)::length
		if(allocated(A%workingMemory))then
			if(size(A%workingMemory).ne.length)then
				deallocate(A%workingMemory)
				allocate(A%workingMemory(length))
			end if
		else
			allocate(A%workingMemory(length))
		end if
		return
	end subroutine

	subroutine deallocatememory(BTool)
		class(OptimEngine_structure),intent(inout)::BTool
		integer::i
		if(BTool%length.eq.0)return
		do i=1,BTool%length
			call BTool%st(i)%deallocate()
			call BTool%yt(i)%deallocate()
		end do
		deallocate(BTool%yt)
		deallocate(BTool%st)
		BTool%length=0
		BTool%FullFlag=.false.
		BTool%endindex=0
		return
	end subroutine

	subroutine resetEndpoint(BTool)
		class(OptimEngine_structure),intent(inout)::BTool
		BTool%endindex=0
		BTool%FullFlag=.false.
		return
	end subroutine
	function datalength(BTool)
		integer::datalength
		class(OptimEngine_structure),intent(in)::BTool
		if(BTool%FullFlag)then
			datalength=BTool%length
		else
			datalength=BTool%endindex
		end if
		return
	end function

	function dataSize(BTool)
		integer::dataSize
		class(OptimEngine_structure),intent(in)::BTool
		dataSize=BTool%length
		return
	end function

	subroutine LBFGSsetprintFlag(OE,printFlag)
		class(OptimEngine_structure), intent(inout) :: OE
		logical,intent(in)::printFlag
		OE%printFlag=printFlag
		return
	end subroutine

	subroutine pointMemory(OE,p,ith)
		class(OptimEngine_structure),target, intent(inout) :: OE
		type(Tensor),pointer::p
		integer,intent(in)::ith
		if(.not.allocated(OE%workingMemory))then
			call writemess('ERROR in pointMemory, DO NOT allocate memory yet',-1)
			call error_stop
		end if
		if(ith.gt.size(OE%workingMemory))then
			call writemess('ERROR in pointMemory, ith> size(memory)',-1)
			call error_stop
		end if
		p=>OE%workingMemory(ith)
		return
	end subroutine

	function get_method(OE)
		character(len=200)::get_method
		class(OptimEngine_structure), intent(inout) :: OE
		get_method=OE%method
		return
	end function

	!***************************************************
	!           initial function for OptimEngine_structure
	!***************************************************


	subroutine set_Metric_func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(Metric_func_interface)::Func
		call writemess('Set the Metric_func in OptimEngine_structure')
		OE%Metric_func=>Func
		return
	end subroutine

	subroutine set_SJ_func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(SJ_interface)::Func
		call writemess('Set the SJ_func in OptimEngine_structure')
		call writemess(' One should calculat the differential matrix of intermediated space Y')
		call writemess(' And the metric G in the space Y')
		call writemess(' The metric in the problem is Y^T * G * Y')
		OE%SJ=>Func
		return
	end subroutine
	subroutine set_SJPointer_func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(SJPointer_interface)::Func
		call writemess('Set the SJPointer_func in OptimEngine_structure')
		call writemess(' One should calculat the differential matrix of intermediated space Y')
		call writemess(' And the metric G in the space Y')
		call writemess(' The metric in the problem is Y^T * G * Y')
		call writemess(' SJPointer_func to make the pointer Yp =>Y and Gp => G')
		OE%SJPointer=>Func
		return
	end subroutine
	subroutine set_SJ_MV_func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(matrix_time_vector_interface)::Func
		call writemess('Set the matrix times vector function used in SJ in OptimEngine_structure')
		call writemess(' The metric in the problem is M=Y^T * G * Y')
		call writemess(' Function is to output the p, p=M*x')
		OE%LLSMV=>Func
		return
	end subroutine
	subroutine set_SJ_MMV_func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(matrix_matrix_time_vector_interface)::Func
		call writemess('Set the transpose_matrix times matrix times vector function used in SJ in OptimEngine_structure')
		call writemess(' The metric in the problem is M=Y^T * G * Y')
		call writemess(' Function is to output the p from a giving input tensor x and logical MVFlag,')
		call writemess('  if MVFlag=.false.')
		call writemess('     p=M^T*M*x ')
		call writemess('  if MVFlag=.true.')
		call writemess('     p=M^T*x ')
		OE%LLSMMV=>Func
		return
	end subroutine
	subroutine set_CGLLS_parameter(OE,CGLLSmaxRunning,CGLLSerror)
		class(OptimEngine_structure), intent(inout) :: OE
		integer,intent(in)::CGLLSmaxRunning
		real*8,intent(in)::CGLLSerror
		call writemess('Set the parameter in CGLLS, which will use only in SJ')
		call writemess(' maxRunning for the root of Ax=b is CGLLSmaxRunning='+CGLLSmaxRunning)
		call writemess(' step error for the root of Ax=b is CGLLSerror='+CGLLSerror)
		OE%CGLLSmaxRunning=CGLLSmaxRunning
		OE%CGLLSerror=CGLLSerror
		return
	end subroutine

	subroutine setInititalCGLLSStartPoint(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(InititalCGLLSStartPointInterface)::Func
		call writemess('Set the InititalCGLLSStartPoint function, it will use for CGLLS')
		OE%InititalCGLLSStartPoint=>Func
	end subroutine

	subroutine set_Hessian_func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(Metric_func_interface)::Func
		call writemess('Set the Metric_func (Hessian) in OptimEngine_structure')
		OE%Metric_func=>Func
		return
	end subroutine


	subroutine set_CG_direction_flag(OE,Flag)
		class(OptimEngine_structure), intent(inout) :: OE
		integer,intent(in)::Flag
		call writemess('Set the CG_direction_flag((only use in method=CG) in OptimEngine_structure')
		call writemess(' CG_direction_flag='+Flag)
		call writemess('     1: Crowder_Wolfe')
		call writemess('     2: Fletcher_Reeves')
		call writemess('     3: Dixon')
		call writemess('     4: Polak_Ribiere_Polyak')
		call writemess('     5: Dai_Yuan')
		OE%CG_direction_Flag=Flag
		return
	end subroutine

	subroutine set_stop_error(OE,error)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::error
		call writemess('Set the step error(stop_gradient) in OptimEngine_structure')
		OE%stop_gradient=error
		return
	end subroutine
	subroutine set_stop_Step(OE,step)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::step
		call writemess('Set the stop length of the step(stop_step) in OptimEngine_structure')
		OE%stop_step=step
		return
	end subroutine
	subroutine set_max_step_in_Linear_search(OE,step)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::step
		call writemess('Set the max step in Linear search in OptimEngine_structure')
		OE%max_step_in_Linear_search=step
		return
	end subroutine
	subroutine set_first_step_in_Linear_search(OE,step)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::step
		call writemess('Set the first step in Linear search in OptimEngine_structure')
		OE%first_step_in_Linear_search=step
		return
	end subroutine
	subroutine set_zero_gradient_in_RGM(OE,zero_gradient)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::zero_gradient
		call writemess('Set zero_gradient_in_RGM(only use in method=RGM) in OptimEngine_structure')
		OE%zero_gradient_in_RGM=zero_gradient
		return
	end subroutine

	function check_stop(OE,value,gradient,step)
		logical::check_stop
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(in)::gradient
		real*8,intent(in)::step,value
		real*8::norm
		character(len=characterlen)::w
		integer,save::stop_counter=0

		check_stop=.false.
		if(OE%stop_gradient.gt.0) then
			norm=gradient%dnorm()
			if(abs(value).gt.1d0)then
				check_stop=abs(norm/value).le.OE%stop_gradient
				if(check_stop)then
					w='Search is going to stop, norm of the |gradient/value|='+abs(norm/value)
					w=w+', The giving stop error='+OE%stop_gradient
					call writemess(w)
					return
				end if
			else
				check_stop=norm.le.OE%stop_gradient
				if(check_stop)then
					w='Search is going to stop, norm of the gradient='+norm
					w=w+', The giving stop error='+OE%stop_gradient
					call writemess(w)
					return
				end if
			end if
		end if

		if(step.le.OE%stop_step)then
			stop_counter=stop_counter+1
			if(stop_counter.ge.max_stop_counter)then
				check_stop=.true.
				w='Search is going to stop, the search step='+step
				call writemess(w)
				return
			end if
		else
			stop_counter=0
		end if
		return
	end function

	subroutine set_linear_search_type(OE,flag)
		class(OptimEngine_structure), intent(inout) :: OE
		integer,intent(in)::flag
		if(flag.eq.1)then
			call writemess('Set LinearSearch type as type1  in OptimEngine_structure') 
			OE%LinearSearch=>LinearSearch1
			return
		end if
		if(flag.eq.2)then
			call writemess('Set LinearSearch type as type2  in OptimEngine_structure') 
			call writemess('Alway keep the point with min value')
			OE%LinearSearch=>LinearSearch2
			return
		end if
		if(flag.eq.3)then
			call writemess('Set LinearSearch type as type3  in OptimEngine_structure') 
			call writemess('Use the penalty function to avoid too large output of x')
			OE%LinearSearch=>LinearSearch3
			return
		end if
		if(flag.eq.4)then
			call writemess('Set LinearSearch type as type4  in OptimEngine_structure') 
			call writemess('Use only the Energy in the Searching') 
			OE%LinearSearch=>LinearSearch4
			return
		end if
		call writemess('No such case',-1)
		call error_stop
	end subroutine

	subroutine set_linear_search_function(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(LinearSearch_subroutine)::Func
		call writemess('Set the LinearSearch subroutine as the external function in OptimEngine_structure')
		OE%LinearSearch=>Func
		return
	end subroutine

	subroutine Set_inStep1Func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(Step1Subroutine_interface)::Func
		call writemess('Set the inStep1Func in OptimEngine_structure')
		call writemess(' the inStep1Func is called at the begining of the loop and before calling the targetfunction')
		OE%inStep1=>Func
		return
	end subroutine
	subroutine Set_inStep2Func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(Step2Subroutine_interface)::Func
		call writemess('Set the inStep2Func in OptimEngine_structure')
		call writemess(' the inStep2Func is after calling the targetfunction and before setting the searching direction')
		OE%inStep2=>Func
		return
	end subroutine
	subroutine Set_inStep3Func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(Step3Subroutine_interface)::Func
		call writemess('Set the inStep3Func in OptimEngine_structure')
		call writemess(' the inStep3Func is after calling the searching_direction and before going on a step')
		OE%inStep3=>Func
		return
	end subroutine
	subroutine Set_inStep4Func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(Step4Subroutine_interface)::Func
		call writemess('Set the inStep4Func in OptimEngine_structure')
		call writemess(' the inStep4Func is after going on a step and at the end of the loop ')
		OE%inStep4=>Func
		return
	end subroutine
	subroutine unSet_inStep1Func(OE)
		class(OptimEngine_structure), intent(inout) :: OE
		call writemess('unSet the inStep1Func in OptimEngine_structure')
		OE%inStep1=>null()
		return
	end subroutine
	subroutine unSet_inStep2Func(OE)
		class(OptimEngine_structure), intent(inout) :: OE
		call writemess('inSet the inStep2Func in OptimEngine_structure')
		OE%inStep2=>null()
		return
	end subroutine
	subroutine unSet_inStep3Func(OE)
		class(OptimEngine_structure), intent(inout) :: OE
		call writemess('Set the inStep3Func in OptimEngine_structure')
		OE%inStep3=>null()
		return
	end subroutine
	subroutine unSet_inStep4Func(OE)
		class(OptimEngine_structure), intent(inout) :: OE
		call writemess('Set the inStep3Func in OptimEngine_structure')
		OE%inStep3=>null()
		return
	end subroutine
	subroutine Set_BeforeStepFunc(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(BeforeStepSubroutine_interface)::Func
		call writemess('Set the BeforeStepFunc in OptimEngine_structure')
		OE%BeforeStep=>Func
		return
	end subroutine
	subroutine unSet_BeforeStepFunc(OE)
		class(OptimEngine_structure), intent(inout) :: OE
		call writemess('unSet the BeforeStepFunc in OptimEngine_structure')
		OE%BeforeStep=>null()
		return
	end subroutine
	subroutine Set_EndStepFunc(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(EndStepSubroutine_interface)::Func
		call writemess('Set the EndStepFun in OptimEngine_structure')
		OE%EndStep=>Func
		return
	end subroutine
	subroutine unSet_EndStepFunc(OE)
		class(OptimEngine_structure), intent(inout) :: OE
		call writemess('unSet the EndStepFun in OptimEngine_structure')
		OE%EndStep=>null()
		return
	end subroutine

	subroutine set_method0(OE,method)
		class(OptimEngine_structure), intent(inout) :: OE
		character(len=*),intent(in)::method
		OE%method=method
		call writemess('Set the optimal method in OptimEngine_structure as:'+method)
		select case(method)
			case ('GM')
				OE%RescalDirection=>DefaultRescalDirection
				OE%AllocateWorkingMemory=>GMMemory
				OE%set_Direction=>set_GMDirection
			case ('RGM')
				OE%RescalDirection=>RandomGradientRescalDirection
				OE%AllocateWorkingMemory=>RGMMemory
				OE%set_Direction=>set_RGMDirection
			case ('CG')
				OE%RescalDirection=>DefaultRescalDirection
				OE%AllocateWorkingMemory=>CGMemory
				OE%set_Direction=>set_CGDirection
			case ('BFGS')
				OE%RescalDirection=>DefaultRescalDirection
				OE%AllocateWorkingMemory=>BFGSMemory
				OE%set_Direction=>set_BFGSDirection
			case ('LBFGS')
				OE%RescalDirection=>DefaultRescalDirection
				OE%AllocateWorkingMemory=>LBFGSMemory
				OE%set_Direction=>set_LBFGSDirection
			case ('adam')
				OE%RescalDirection=>adamRescalDirection
				OE%AllocateWorkingMemory=>AdamMemory
				OE%set_Direction=>set_AdamDirection
			case ('Newton')
				OE%RescalDirection=>DefaultRescalDirection
				OE%AllocateWorkingMemory=>Natural_Memory
				OE%set_Direction=>set_natural_direction
			case ('NG')
				OE%RescalDirection=>DefaultRescalDirection
				OE%AllocateWorkingMemory=>Natural_Memory
				OE%set_Direction=>set_natural_direction
			case ('SJ')
				OE%RescalDirection=>DefaultRescalDirection
				OE%AllocateWorkingMemory=>SJ_Memory
				OE%set_Direction=>set_SJ_direction
			case ('GM2')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>GMMemory
				OE%set_Direction=>set_GMDirection
				call writemess(' DO NOT rescal the searching direction')
			case ('CG2')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>CGMemory
				OE%set_Direction=>set_CGDirection
				call writemess(' DO NOT rescal the searching direction')
			case ('BFGS2')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>BFGSMemory
				OE%set_Direction=>set_BFGSDirection
				call writemess(' DO NOT rescal the searching direction')
			case ('LBFGS2')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>LBFGSMemory
				OE%set_Direction=>set_LBFGSDirection
				call writemess(' DO NOT rescal the searching direction')
			case ('Newton2')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>Natural_Memory
				OE%set_Direction=>set_natural_direction
				call writemess(' DO NOT rescal the searching direction')
			case ('NG2')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>Natural_Memory
				OE%set_Direction=>set_natural_direction
				call writemess(' DO NOT rescal the searching direction')
			case ('SJ2')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>SJ_Memory
				OE%set_Direction=>set_SJ_direction
				call writemess(' DO NOT rescal the searching direction')
			case default
				call writemess('NO such case of method, the default method are:')
				call writemess(' GM     :  Gradient Method')
				call writemess(' RGM    : Random Gradient Method')
				call writemess(' CG     : Conjugate Gradient method ')
				call writemess(' BFGS   : BFGS method ')
				call writemess(' LBFGS  : LBFGS method ')
				call writemess(' adam   : Adam method ')
				call writemess(' Newton : Newton method ')
				call writemess(' NG     : natural gradient method ')
				call writemess(' ****2  : the same method, but Do not rescall the direction')
				call writemess(' If Newton method is in used, the subroutine getting the  Hessian should be add')
				call writemess(' If natural gradient method is in used, the subroutine getting the metric should be add')
				call error_stop
		end select
		return
	end subroutine

	subroutine set_method1(OE,Rescal_Direction,Allocate_Memory,set_Direction)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(RescalDirection_interface)::Rescal_Direction
		procedure(AllocateMemory_interface)::Allocate_Memory
		procedure(set_DefaultDirection_interface)::set_Direction
		OE%method='other'
		call writemess('Set the optimal method in OptimEngine_structure as user define method')
		OE%RescalDirection=>Rescal_Direction
		OE%AllocateWorkingMemory=>Allocate_Memory
		OE%set_Direction=>set_Direction
		return
	end subroutine

	subroutine set_penaltyFactor(OE,factor)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::factor
		call writemess('Set the penaltyfunctionFactor='+factor)
		OE%penaltyfunctionFactor=factor
		return
	end subroutine

	!***************************************************
	!           parameter for adam
	!***************************************************

	subroutine Set_adam_reset_Step(OE,Nstep)
		class(OptimEngine_structure), intent(inout) :: OE
		integer,intent(in)::Nstep
		OE%Adam_resetStep=Nstep
		call writemess('Set Adam_resetStep='+Nstep)
		return
	end subroutine
	subroutine Set_adam_rhos(OE,rho1,rho2)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::rho1,rho2
		OE%Adam_rho1=rho1
		OE%Adam_rho2=rho2
		call writemess('Set Adam_rho1='+rho1)
		call writemess('Set Adam_rho2='+rho2)
		return
	end subroutine

	subroutine Set_adam_delta(OE,delta,message)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::delta
		logical,optional,intent(in)::message
		OE%Adam_delta=delta
		call writemess('Set Adam_delta='+delta)
		if(present(message))then
			OE%Adam_warning=message
		end if
		if(message)then
			call writemess('If gradient(i)^2 < delta, than regard it as 0, and print the warning')
		else
			call writemess('If gradient(i)^2 < delta, than regard it as 0, But do not print the warning')
		end if
		return
	end subroutine
	!***************************************************
	!           Basic function for OptimRunner
	!***************************************************

	subroutine set_target_function(LBFGSer,Func)
		class(OptimRunner),intent(inout)::LBFGSer
		procedure(external_target_Function_interface)::Func
		call writemess('Set the target_function in OptimRunner')
		LBFGSer%externalFunc=>Func
		return
	end subroutine

	subroutine LBFGSFunc(A,outVal,outGradient,point)
		class(OptimRunner),target, intent(inout) :: A
		real*8,target,intent(inout)::outVal
		type(Tensor),target,intent(inout)::outGradient
		type(Tensor),target,intent(in)::point
		if(associated(A%externalFunc))then
			call A%externalFunc(outVal,outGradient,point)
		else
			call writemess('DO not set the target function yet',-1)
			call error_stop
		end if
		return
	end subroutine

	!***************************************************
	!           default function
	!***************************************************

	subroutine DefaultRescalDirection(OE,Dir)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(inout)::Dir
		Dir=Dir/Dir%dnorm()
		return
	end subroutine

	subroutine DefaultNotRescalDirection(OE,Dir)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(inout)::Dir
		return
	end subroutine

	subroutine DefaultAllocateWorkingMemory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),pointer,intent(inout)::Gradient,direction
		call printERRORMessage()
	end subroutine
	subroutine set_DefaultDirection(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(in)::point,Gradient
		type(Tensor),intent(inout)::direction
		call printERRORMessage()
	end subroutine
	subroutine printERRORMessage()
		call writemess('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		call writemess('% DO NOT set the optimal method yet     %')
		call writemess('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		call error_stop
	end subroutine

	
	!***************************************************
	!            optimization
	!***************************************************


	subroutine Optimization1(OE,T0,tau0,numStep,Point,outValue)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::T0,tau0
		integer,intent(in)::numStep
		type(Tensor),intent(inout)::Point
		real*8,intent(inout)::outValue
		integer::i
		real*8::t
		type(Tensor),pointer::Gradient,direction

		if((OE%method.equ.'LBFGS').and.(OE%dataSize().le.0))then
			call writemess('ERROR in LBFGS',-1)
			call writemess(' DO NOT allocate memory for running',-1)
			call error_stop
		end if
		if(OE%printFlag)call reset_time_calculator(numStep,write_time_num) 
		call OE%resetEndpoint()
		call OE%AllocateWorkingMemory(Gradient,direction)
		if(associated(OE%BeforeStep))call OE%BeforeStep(outValue,Gradient,point)
		do i=1,numStep
			t=T0*tau0/(T0+dble(i))
			if(associated(OE%inStep1))call OE%inStep1(outValue,direction,Gradient,point,i,t)
			call OE%target_Function(outValue,Gradient,point)
			if(associated(OE%inStep2))call OE%inStep2(outValue,direction,Gradient,point,i,t)
			call OE%set_direction(direction,Gradient,point)
			call OE%RescalDirection(direction)
			if(associated(OE%inStep3))call OE%inStep3(outValue,direction,Gradient,point,i,t)
			point=point+(direction*t)
			if(associated(OE%inStep4))call OE%inStep4(outValue,direction,Gradient,point,i,t)
			if(OE%printFlag)call time_calculator()
			if(OE%check_stop(outValue,Gradient,t))exit
		end do
		if(associated(OE%EndStep))call OE%EndStep(outValue,Gradient,point)
		return
	end subroutine

	subroutine Optimization2(OE,t,numStep,Point,outValue)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::t
		integer,intent(in)::numStep
		type(Tensor),intent(inout)::Point
		real*8,intent(inout)::outValue
		integer::i
		type(Tensor),pointer::Gradient,direction
		if((OE%method.equ.'LBFGS').and.(OE%dataSize().le.0))then
			call writemess('ERROR in LBFGS',-1)
			call writemess(' DO NOT allocate memory for running',-1)
			call error_stop
		end if
		if(OE%printFlag)call reset_time_calculator(numStep,write_time_num) 
			
		call OE%AllocateWorkingMemory(Gradient,direction)
		if(associated(OE%BeforeStep))call OE%BeforeStep(outValue,Gradient,point)
		do i=1,numStep
			if(associated(OE%inStep1))call OE%inStep1(outValue,direction,Gradient,point,i,t)
			call OE%target_Function(outValue,Gradient,point)
			if(associated(OE%inStep2))call OE%inStep2(outValue,direction,Gradient,point,i,t)
			call OE%set_direction(direction,Gradient,point)
			call OE%RescalDirection(direction)
			if(associated(OE%inStep3))call OE%inStep3(outValue,direction,Gradient,point,i,t)
			point=point+(direction*t)
			if(associated(OE%inStep4))call OE%inStep4(outValue,direction,Gradient,point,i,t)
			if(OE%printFlag)call time_calculator()
			if(OE%check_stop(outValue,Gradient,t))exit
		end do
		if(associated(OE%EndStep))call OE%EndStep(outValue,Gradient,point)
		return
	end subroutine

	subroutine Optimization3(OE,NlinearStep,numStep,Point,outValue)
		class(OptimEngine_structure), intent(inout) :: OE
		integer,intent(in)::NlinearStep,numStep
		type(Tensor),intent(inout)::Point
		real*8,intent(inout)::outValue
		real*8::x
		integer::i
		type(Tensor),pointer::Gradient,direction
		if((OE%method.equ.'LBFGS').and.(OE%dataSize().le.0))then
			call writemess('ERROR in LBFGS',-1)
			call writemess(' DO NOT allocate memory for running',-1)
			call error_stop
		end if
		if(OE%printFlag)call reset_time_calculator(numStep,write_time_num) 
		call OE%AllocateWorkingMemory(Gradient,direction)
		if(associated(OE%BeforeStep))call OE%BeforeStep(outValue,Gradient,point)
		x=OE%first_step_in_Linear_search
		do i=1,numStep
			if(associated(OE%inStep1))call OE%inStep1(outValue,direction,Gradient,point,i,x)
			call OE%target_Function(outValue,Gradient,point)
			if(associated(OE%inStep2))call OE%inStep2(outValue,direction,Gradient,point,i,x)
			call OE%set_direction(direction,Gradient,point)
			call OE%RescalDirection(direction)
			if(associated(OE%inStep3))call OE%inStep3(outValue,direction,Gradient,point,i,x)
			call OE%LinearSearch(NlinearStep,point,direction,x,outValue,Gradient)
			if(associated(OE%inStep4))call OE%inStep4(outValue,direction,Gradient,point,i,x)
			if(OE%printFlag)call time_calculator()
			if(OE%check_stop(outValue,Gradient,x))exit
		end do
		if(associated(OE%EndStep))call OE%EndStep(outValue,Gradient,point)
		return
	end subroutine

	!***************************************************
	!           Adam optimization
	!***************************************************

	subroutine adamRescalDirection(OE,Dir)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(inout)::Dir
	end subroutine

	subroutine AdamMemory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),pointer,intent(inout)::Gradient,direction
		integer::memorylen
		type(Tensor),pointer::s,g
		memorylen=4
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.memorylen)then
				deallocate(OE%workingMemory)
				allocate(OE%workingMemory(memorylen))
			end if
		else
			allocate(OE%workingMemory(memorylen))
		end if
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		call OE%pointMemory(s,3)
		call OE%pointMemory(g,4)
		call s%empty()
		call g%empty()
		return
	end subroutine

	subroutine set_AdamDirection(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(in)::point,Gradient
		type(Tensor),intent(inout)::direction
		type(Tensor),pointer::s,r
		real*8::a,b,snorm,rnorm
		real*8::rho,rhot
		call OE%pointMemory(s,3)
		call OE%pointMemory(r,4)

		OE%Adam_resetcounter=OE%Adam_resetcounter+1
		if(OE%Adam_resetcounter.eq.OE%Adam_resetStep)then
			call s%empty()
			call r%empty()
			OE%Adam_resetcounter=0
		end if

		if(s%getFlag())then
			rho=OE%Adam_rho1
			OE%Adam_rho1t=OE%Adam_rho1t*rho
			rhot=OE%Adam_rho1t
			a=rho!/(1d0-rhot)
			b=(1d0-rho)!/(1d0-rhot)
			s=(a*s)+(b*Gradient)
		else
			s=Gradient
			OE%Adam_rho1t=OE%Adam_rho1
		end if


		if(r%getFlag())then
			rho=OE%Adam_rho2
			OE%Adam_rho2t=OE%Adam_rho2t*rho
			rhot=OE%Adam_rho2t
			a=rho!/(1d0-rhot)
			b=(1d0-rho)!/(1d0-rhot)
			r=(a*r)+(b*element_product(Gradient,Gradient))
		else
			r=element_product(Gradient,Gradient)
			OE%Adam_rho2t=OE%Adam_rho2
		end if

		snorm=1d0/(1d0-OE%Adam_rho1t)
		rnorm=1d0/(1d0-OE%Adam_rho2t)
		direction=element_minu_A_over_sqrt_B(s*snorm,r*rnorm,OE%adam_delta,OE%Adam_warning)
		!call direction%print
		!call s%print
		!call r%print
		!call Gradient%print
		!read(*,*)
		return
	end subroutine


	function element_product(A,B)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		type(Tensor),intent(in)::B
		integer::i,itype,totalData
		real*8,pointer::Ap(:),Bp(:),Resp(:)
		iType=A%getType()
		if(iType.eq.3)then
			totalData=A%getTotalData()
			if(totalDAta.ne.B%getTotalData())then
				call writemess('ERROR in element_product',-1)
				call error_stop
			end if
			call Res%allocate(A%dim(),itype)
			call Res%pointer(Resp)
			call A%pointer(Ap)
			call B%pointer(Bp)
			do i=1,TotalData
				Resp(i)=Ap(i)*Bp(i)
			end do
			return
		end if
		call writemess('ERROR in element_product',-1)
		call error_stop
	end function

	function element_minu_A_over_sqrt_B(A,B,delta,printmessage)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		type(Tensor),intent(in)::B
		real*8,intent(in)::delta
		logical,intent(in)::printmessage
		integer::i,itype,totalData
		real*8,pointer::Ap(:),Bp(:),Resp(:)
		
		iType=A%getType()
		if(iType.eq.3)then
			totalData=A%getTotalData()
			if(totalDAta.ne.B%getTotalData())then
				call writemess('ERROR in element_product',-1)
				call error_stop
			end if
			call Res%allocate(A%dim(),itype)
			call Res%pointer(Resp)
			call A%pointer(Ap)
			call B%pointer(Bp)
			do i=1,TotalData
				if(sqrt(Bp(i)).lt.delta)then
					if((Ap(i).gt.delta).and.printmessage)then
						call writemess('**** WARNING in adam: -A/srqt(B)='+(Ap(i)/sqrt(Bp(i)))+', A='+Ap(i)+' ,sqrt(B)='+sqrt(Bp(i))+'  ****',-1)
					end if
					Resp(i)=-0d0
				else
					Resp(i)=-1d0*Ap(i)/sqrt(Bp(i))
				end if
			end do
			return
		end if
		call writemess('ERROR in element_product',-1)
		call error_stop
	end function

	!***************************************************
	!            Gradient optimization
	!***************************************************

	subroutine GMMemory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),pointer,intent(inout)::Gradient,direction
		integer::memorylen
		memorylen=2
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.memorylen)then
				deallocate(OE%workingMemory)
				allocate(OE%workingMemory(memorylen))
			end if
		else
			allocate(OE%workingMemory(memorylen))
		end if
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		return
	end subroutine
	subroutine set_GMDirection(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(in)::point,Gradient
		type(Tensor),intent(inout)::direction
		direction=(-1d0)*Gradient
		return
	end subroutine

	!***************************************************
	!           random Gradient optimization
	!***************************************************

	subroutine RandomGradient(OE,Direction,gradient)!dir=-1*gra
		class(OptimEngine_structure), intent(in) :: OE
		type(Tensor),intent(in)::gradient
		type(Tensor),intent(inout)::Direction	
		integer::i
		real*8::numr,numi
		real*8,pointer::dirp(:),grap(:)
		complex*16,pointer::zdirp(:),zgrap(:)
		if(gradient%getType().eq.3)then
			call Direction%empty()
			call Direction%allocate([gradient%getTotalData()],'real*8')
			call gradient%pointer(grap)
			call Direction%pointer(dirp)
			do i=1,gradient%getTotalData()
				if(abs(grap(i)).gt.OE%zero_gradient_in_RGM)then
					if(grap(i).gt.0)then
						dirp(i)=-1d0*randomnumber()
					else
						dirp(i)=randomnumber()
					end if
				else
					dirp(i)=0
				end if
			end do
		else if(gradient%getType().eq.5)then
			call Direction%empty()
			call Direction%allocate([gradient%getTotalData()],'complex*16')
			call gradient%pointer(zgrap)
			call Direction%pointer(zdirp)
			do i=1,gradient%getTotalData()
				numr=dreal(zgrap(i))
				numi=aimag(zgrap(i))
				if(abs(numr).gt.OE%zero_gradient_in_RGM)then
					if(numr.gt.0)then
						numr=-1d0*randomnumber()
					else
						numr=randomnumber()
					end if
				else
					numr=0
				end if

				if(abs(numi).gt.OE%zero_gradient_in_RGM)then
					if(numi.gt.0)then
						numi=-1d0*randomnumber()
					else
						numi=randomnumber()
					end if
				else
					numi=0
				end if

				dirp(i)=dcmplx(numr,numi)
			end do
		else
			call writemess('ERROR in data type of the gradient',-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine RandomGradientRescalDirection(OE,Dir)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(inout)::Dir
		return
	end subroutine

	subroutine RGMMemory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),pointer,intent(inout)::Gradient,direction
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.2)then
				deallocate(OE%workingMemory)
				allocate(OE%workingMemory(2))
			end if
		else
			allocate(OE%workingMemory(2))
		end if
		
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		return
	end subroutine
	subroutine set_RGMDirection(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(in)::point,Gradient
		type(Tensor),intent(inout)::direction
		call RandomGradient(OE,direction,Gradient)
		return
	end subroutine

	!***************************************************
	!           Conjugate Gradient optimization
	!***************************************************

	subroutine CG_direction(OE,dir,Gra,priorgra)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(in)::Gra
		type(Tensor),intent(inout)::dir,priorgra
		real*8::direction_corr
		if(dir%getFlag()) then
			direction_corr=correctPara(OE,dir,gra,priorGra)
			dir=(direction_corr*dir)-gra
		else
			dir=(-1d0)*gra
		end if
		priorgra=Gra
		return
	end subroutine
	
	function correctPara(OE,dir,newgra,gra)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::correctPara
		type(Tensor),intent(in)::newgra,gra,dir
		real*8::nor
		select case(OE%CG_direction_Flag)
			case (1)
				correctPara=Crowder_Wolfe(OE,dir,newgra,gra)
			case (2)
				correctPara=Fletcher_Reeves(OE,dir,newgra,gra)
			case (3)
				correctPara=Dixon(OE,dir,newgra,gra)
			case (4)
				correctPara=Polak_Ribiere_Polyak(OE,dir,newgra,gra)
			case (5)
				correctPara=Dai_Yuan(OE,dir,newgra,gra)
			case default
				call writemess('ERROR CG_direction_Flag='+OE%CG_direction_Flag,-1)
				call error_stop
		end select
		return
	end function

	function Crowder_Wolfe(OE,dir,newgra,gra)result(Res)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::Res
		type(Tensor),intent(in)::dir,newgra,gra
		real*8::TMPr
		Res=newgra.dot.(newgra-gra)
		TMPr=dir.dot.(newgra-gra)
		Res=Res/TMPr
		return
	end function
	function Fletcher_Reeves(OE,dir,newgra,gra)result(Res)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::Res
		type(Tensor),intent(in)::dir,newgra,gra
		real*8::TMPr
		Res=newgra%dnorm2()
		TMPr=gra%dnorm2()
		Res=Res/TMPr
		return
	end function
	function Dixon(OE,dir,newgra,gra)result(Res)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::Res
		type(Tensor),intent(in)::dir,newgra,gra
		real*8::TMPr
		Res=newgra%dnorm2()
		TMPr=dir.dot.gra
		Res=-1d0*Res/TMPr
		return
	end function
	function Polak_Ribiere_Polyak(OE,dir,newgra,gra)result(Res)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::Res
		type(Tensor),intent(in)::dir,newgra,gra
		real*8::TMPr
		Res=newgra.dot.(newgra-gra)
		TMPr=gra%dnorm2()
		Res=Res/TMPr
		return
	end function
	function Dai_Yuan(OE,dir,newgra,gra)result(Res)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::Res
		type(Tensor),intent(in)::dir,newgra,gra
		real*8::TMPr
		Res=newgra%dnorm2()
		TMPr=dir.dot.(newgra-gra)
		Res=Res/TMPr
		return
	end function

	subroutine CGMemory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),pointer,intent(inout)::Gradient,direction
		integer::memorylen
		memorylen=4
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.memorylen)then
				deallocate(OE%workingMemory)
				allocate(OE%workingMemory(memorylen))
			end if
		else
			allocate(OE%workingMemory(memorylen))
		end if
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		return
	end subroutine
	subroutine set_CGDirection(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(in)::point,Gradient
		type(Tensor),intent(inout)::direction
		type(Tensor),pointer::gra0,dir0
		call OE%pointMemory(gra0,3)
		call OE%pointMemory(dir0,4)
		call CG_direction(OE,dir0,Gradient,gra0)
		direction=dir0
		return
	end subroutine


	!***************************************************
	!            BFGS optimization
	!***************************************************

	subroutine BFGS_direction(OE,C,y,s)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(inout)::C,y,s
		type(Tensor)::one
		integer::dimn
		type(Tensor)::sy,ss,tempC
		real*8::rho
		if(.not.y%getFlag())return
		dimn=y%getTotalData()
		call one%eye(dimn,dimn,y%getType())
		if(.not.C%getFlag())then
			call writemess('ERROR in BFGS_direction',-1)
			call error_stop
		end if
		rho=y.dot.s
		rho=1d0/rho
		sy=s.kron.y
		ss=s.kron.s
		tempC=one-(rho*sy)
		C=tempC * C
		C=C*( tempC.p.[2,1] )
		C=C+(rho*ss)
		!C=( (one-(rho*sy)) * C * ( one-(rho*(.T.sy)) ) )+(rho*ss)
		return
	end subroutine

	subroutine BFGSMemory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),pointer,intent(inout)::Gradient,direction
		type(Tensor),pointer::Ct,yt,st,SavePoint,SaveGradient
		integer::memorylen
		memorylen=7
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.memorylen)then
				deallocate(OE%workingMemory)
				allocate(OE%workingMemory(memorylen))
			end if
		else
			allocate(OE%workingMemory(memorylen))
		end if
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		call OE%pointMemory(Ct,3)
		call OE%pointMemory(yt,4)
		call OE%pointMemory(st,5)
		call OE%pointMemory(SavePoint,6)
		call OE%pointMemory(SaveGradient,7)
		call yt%empty()
		call St%empty()
		call Ct%empty()
		call SavePoint%empty()
		call SaveGradient%empty()
		return
	end subroutine

	subroutine set_BFGSDirection(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(in)::point,Gradient
		type(Tensor),intent(inout)::direction
		type(Tensor),pointer::Ct,yt,st,SavePoint,SaveGradient
		call OE%pointMemory(Ct,3)
		call OE%pointMemory(yt,4)
		call OE%pointMemory(st,5)
		call OE%pointMemory(SavePoint,6)
		call OE%pointMemory(SaveGradient,7)
		st=point-SavePoint
		yt=Gradient-SaveGradient
		SavePoint=point
		SaveGradient=Gradient

		if(.not.Ct%getFlag())then
			call Ct%eye(Point%getTotalData(),Point%getTotalData(),Point%getType())
			Ct=((st.dot.yt)/(yt.dot.yt))*Ct
		else
			call BFGS_direction(OE,Ct,yt,st)
		end if
		direction=Ct*Gradient*(-1d0)
		
		return
	end subroutine


	!***************************************************
	!            LBFGS optimization
	!***************************************************

	! p is \Delta f(x_i)
	! s_i= x_{i+1} - x_i
	! y_i= p_{i+1} - p_i

	subroutine LBFGS_direction(LBFGS,p,inputp)
		class(OptimEngine_structure),intent(inout)::LBFGS
		type(Tensor),intent(inout)::p
		type(Tensor),intent(in)::inputp
		real*8::temp,temp2,beta
		real*8,allocatable::alpha(:)
		integer::i,length
		type(Tensor),pointer::s,y
		length=LBFGS%getLength()
		p=inputp*(-1d0)
		if(length.eq.0)return

		allocate(alpha(length))
		do i=length,1,-1
			call LBFGS%i(s,y,i)
			temp=y.dot.s
			alpha(i)=s.dot.p
			alpha(i)=alpha(i)/temp
			p=p-(y*alpha(i))
		end do

		call LBFGS%i(s,y,length)
		temp=y.dot.s
		temp=temp/y%dnorm2()
		p=p*temp
	

		do i=1,length
			call LBFGS%i(s,y,i)
			temp=y.dot.s
			beta=y.dot.p
			beta=beta/temp
			p=p+(s*(alpha(i)-beta))
		end do
		return
	end subroutine

	subroutine LBFGSMemory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),pointer,intent(inout)::Gradient,direction
		integer::memorylen
		type(Tensor),pointer::SavePoint,SaveGradient
		memorylen=4
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.memorylen)then
				deallocate(OE%workingMemory)
				allocate(OE%workingMemory(memorylen))
			end if
		else
			allocate(OE%workingMemory(memorylen))
		end if
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		call OE%pointMemory(SavePoint,3)
		call OE%pointMemory(SaveGradient,4)
		call SavePoint%empty()
		call SaveGradient%empty()
		return
	end subroutine
	subroutine set_LBFGSDirection(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(in)::point,Gradient
		type(Tensor),intent(inout)::direction
		type(Tensor),pointer::yt,st,SavePoint,SaveGradient
		call OE%pointMemory(SavePoint,3)
		call OE%pointMemory(SaveGradient,4)
		call OE%NewElement(st,yt)
		yt=Gradient-SaveGradient
		st=point-SavePoint
		SaveGradient=Gradient
		SavePoint=point
		call LBFGS_direction(OE,Direction,Gradient)
		return
	end subroutine
	!***************************************************
	!            Newton of natural gradient optimization
	!***************************************************


	subroutine Natural_Memory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),pointer,intent(inout)::Gradient,direction
		integer::memorylen
		memorylen=3
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.memorylen)then
				deallocate(OE%workingMemory)
				allocate(OE%workingMemory(memorylen))
			end if
		else
			allocate(OE%workingMemory(memorylen))
		end if
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		call OE%pointMemory(OE%metric,3)
		call OE%metric%empty()
		return
	end subroutine


	subroutine set_natural_direction(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(in)::point,Gradient
		type(Tensor),intent(inout)::direction
		type(Tensor),pointer::gra0,dir0
		if(associated(OE%Metric_func))then
			call OE%Metric_func(OE%metric,direction,Gradient,point)
		else
			if(.not.OE%metric%getFlag())then
				call writemess('DO not set the Metric_func yet',-1)
				call error_stop
			end if
		end if

		if(OE%solve_linear_equation_method.eq.2)then
			if(OE%SolveLinearEquationRCOND.gt.0)then
				CALL direction%Solvelinequ(OE%metric,Gradient,OE%SolveLinearEquationRCOND)
			else
				CALL direction%Solvelinequ(OE%metric,Gradient)
			end if
		else if(OE%solve_linear_equation_method.eq.1)then
			CALL direction%SolveLLS(OE%metric,Gradient,OE%SolveLinearEquationRCOND)
		else
			call writemess('ERROR solve_linear_equation_method,solve_linear_equation_method='+OE%solve_linear_equation_method)
			call error_stop
		end if
		direction=(-1d0)*direction
		return
	end subroutine

	!***************************************************
	!           NoName optimization
	!***************************************************


	subroutine SJ_Memory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),pointer,intent(inout)::Gradient,direction
		integer::memorylen
		type(Tensor),pointer::Y,G
		memorylen=8
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.memorylen)then
				deallocate(OE%workingMemory)
				allocate(OE%workingMemory(memorylen))
			end if
		else
			allocate(OE%workingMemory(memorylen))
		end if
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		return
	end subroutine


	subroutine set_SJ_direction(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor),intent(in)::point,Gradient
		type(Tensor),intent(inout)::direction
		type(Tensor),pointer::gra0,dir0,Y,G,Ap,r,p,Ab
		call OE%pointMemory(Ap,3)
		call OE%pointMemory(r,4)
		call OE%pointMemory(p,5)
		call OE%pointMemory(Ab,6)
		if(associated(OE%SJ))then
			call OE%pointMemory(Y,7)
			call OE%pointMemory(G,8)
			call OE%SJ(Y,G,direction,Gradient,point)
			call OE%runCGLLS1(direction,Gradient,Y,G,Ap,r,p)
		else if(associated(OE%SJPointer)) then
			call OE%SJPointer(Y,G,direction,Gradient,point)
			call OE%runCGLLS1(direction,Gradient,Y,G,Ap,r,p)
		else if(associated(OE%LLSMV)) then
			call OE%runCGLLS2(direction,Gradient,Ap,r,p)
		else if(associated(OE%LLSMMV)) then
			call OE%runCGLLS_MMx(direction,Gradient,Ap,r,p,Ab)
		else
			call writemess('DO not set the SJ yet',-1)
			call error_stop
		end if
		
		direction=(-1d0)*direction
		return
	end subroutine


	!*************************************************************
	!*************************************************************
	!    linear search1 use the information of Energy and Gradient
	!*************************************************************

	subroutine LinearSearch1(OE,max_running,point,dir,x,outValue,inoutgra)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(inout)::x,outValue
		type(Tensor),intent(inout)::point
		type(Tensor),intent(inout)::inoutgra
		type(Tensor),intent(in)::dir
		integer,intent(in)::max_running
		logical::stopflag
		real*8::x1,x2,x3,outx
		real*8::f1,f2,f3,g1,g2,g3
		type(Tensor)::gradient,NewPoint
		integer::i
		type(Tensor)::Savepoint,Savegradient
		real*8::savex,saveValue
		logical::max_step_flag
		if(max_running.lt.2)then
			call writemess('ERROR in LinearSearch,input max_running should >=2',-1)
			call error_stop
		end if
		max_step_flag=OE%max_step_in_Linear_search.gt.0
		x1=0d0
		g1=inoutgra.dot.dir
		f1=outValue

		x2=x
		if(x2.le.OE%stop_step) then
			x2=OE%first_step_in_Linear_search
			x=x2
		end if

		NewPoint=point+(x2*dir)
		call OE%target_Function(f2,gradient,NewPoint)
		g2=gradient.dot.dir
		if(max_step_flag)then
			Savepoint=NewPoint
			Savegradient=gradient
			savex=x2
			saveValue=f2
		end if


		if(dabs(g2).le.OE%stop_gradient)then
			Point=NewPoint
			x=x2
			outValue=f2
			inoutgra=gradient
			return
		end if
		call LinearSearch_third_point(outx,f1,g1,f2,g2,x2)
		x3=outx
		NewPoint=point+(x3*dir)
		call OE%target_Function(f3,gradient,NewPoint)
		g3=gradient.dot.dir

		if(dabs(g3).le.OE%stop_gradient)then
			Point=NewPoint
			x=x3
			outValue=f3
			inoutgra=gradient
			return
		end if

		stopflag=.false.

		do i=3,max_running
			call LinearSearch_fouth_point(outx,stopflag,x1,f1,g1,x2,f2,g2,x3,f3,g3,x)
			if(isnan(outx).or.(outx.lt.0d0))then
				write(*,*)"ERROR in LinearSearch"
				write(*,*)"x1,x2"
				write(*,*)x1
				write(*,*)x2
				stop
			end if
			NewPoint=Point+(outx*dir)
			call OE%target_Function(f3,gradient,NewPoint)
			g3=gradient.x.dir
			x3=outx
			if((dabs(g3).le.OE%stop_gradient).or. stopflag) exit
		end do
		if(max_step_flag.and.(outx.gt.OE%max_step_in_Linear_search))then
			Point=Savepoint
			x=Savex
			outValue=SaveValue
			inoutgra=Savegradient
		else
			Point=NewPoint
			x=outx
			outValue=f3
			inoutgra=gradient
		end if
		return
	end subroutine


	!*************************************************************
	!*************************************************************
	!  linear search2 use the information of Energy and Gradient
	!*************************************************************

	!LinearSearch2
		!Alway keep the point with min value

	subroutine LinearSearch2(OE,max_running,point,dir,x,outValue,inoutgra)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(inout)::x,outValue
		type(Tensor),intent(inout)::point
		type(Tensor),intent(inout)::inoutgra
		type(Tensor),intent(in)::dir
		integer,intent(in)::max_running
		logical::stopflag
		real*8::x1,x2,x3,outx
		real*8::f1,f2,f3,g1,g2,g3
		type(Tensor)::gradient,NewPoint
		integer::i
		type(Tensor)::minpoint,mingradient
		real*8::minx,minValue
		logical::max_step_flag
		if(max_running.lt.2)then
			call writemess('ERROR in LinearSearch,input max_running should >=2',-1)
			call error_stop
		end if
		max_step_flag=OE%max_step_in_Linear_search.gt.0
		x1=0d0
		g1=inoutgra.dot.dir
		f1=outValue

		x2=x
		if(x2.le.OE%stop_step) then
			x2=OE%first_step_in_Linear_search
			x=x2
		end if

		NewPoint=point+(x2*dir)
		call OE%target_Function(f2,gradient,NewPoint)
		g2=gradient.dot.dir

		minx=x2
		minValue=f2
		minpoint=NewPoint
		mingradient=gradient


		if(dabs(g2).le.OE%stop_gradient)then
			Point=NewPoint
			x=x2
			outValue=f2
			inoutgra=gradient
			return
		end if

		
		call LinearSearch_third_point(outx,f1,g1,f2,g2,x2)
		x3=outx
		NewPoint=point+(x3*dir)
		call OE%target_Function(f3,gradient,NewPoint)
		g3=gradient.dot.dir

		if(dabs(g3).le.OE%stop_gradient)then
			Point=NewPoint
			x=x3
			outValue=f3
			inoutgra=gradient
			return
		end if
		if(f3.lt.minValue)then
			if(max_step_flag)then
				if(x3.le.OE%max_step_in_Linear_search)then
					minx=x3
					minValue=f3
					minpoint=NewPoint
					mingradient=gradient
				end if
			else
				minx=x3
				minValue=f3
				minpoint=NewPoint
				mingradient=gradient
			end if
		end if

		stopflag=.false.

		do i=3,max_running
			call LinearSearch_fouth_point(outx,stopflag,x1,f1,g1,x2,f2,g2,x3,f3,g3,x)
			if(isnan(outx).or.(outx.lt.0d0))then
				write(*,*)"ERROR in LinearSearch"
				write(*,*)"x1,x2"
				write(*,*)x1
				write(*,*)x2
				stop
			end if
			NewPoint=Point+(outx*dir)
			call OE%target_Function(f3,gradient,NewPoint)
			g3=gradient.x.dir
			x3=outx
			if(f3.lt.minValue)then
				if(max_step_flag)then
					if(x3.le.OE%max_step_in_Linear_search)then
						minx=x3
						minValue=f3
						minpoint=NewPoint
						mingradient=gradient
					end if
				else
					minx=x3
					minValue=f3
					minpoint=NewPoint
					mingradient=gradient
				end if
			end if
			if((dabs(g3).le.OE%stop_gradient).or. stopflag) exit
		end do
		Point=minpoint
		x=minx
		outValue=minValue
		inoutgra=mingradient
		return
	end subroutine

	!*************************************************************
	!*************************************************************
	!  linear search3 use the information of Energy and Gradient
	!*************************************************************

	!LinearSearch3
		!Use the penalty function to avoid too large output of x

	subroutine penalty_function(OE,outf,outg,gradient,NewPoint,Point,dir,x)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(inout)::outf,outg
		type(Tensor),intent(inout)::gradient,NewPoint
		type(Tensor),intent(in)::Point,dir
		real*8,intent(in)::x
		real*8::x_max,tempx,factor
		x_max=OE%max_step_in_Linear_search
		NewPoint=point+(x*dir)
		call OE%target_Function(outf,gradient,NewPoint)
		outg=gradient.dot.dir
		if(x_max.gt.0)then
			if(x.gt.x_max)then
				tempx=(x-x_max)
				factor=OE%penaltyfunctionFactor
				outf=outf+(tempx*tempx*tempx*tempx*factor)
				outg=outg+(4*factor*tempx*tempx*tempx)
			end if
		end if
		return
	end subroutine

	subroutine LinearSearch3(OE,max_running,point,dir,x,outValue,inoutgra)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(inout)::x,outValue
		type(Tensor),intent(inout)::point
		type(Tensor),intent(inout)::inoutgra
		type(Tensor),intent(in)::dir
		integer,intent(in)::max_running
		logical::stopflag
		real*8::x1,x2,x3,outx
		real*8::f1,f2,f3,g1,g2,g3
		type(Tensor)::gradient,NewPoint
		integer::i
		type(Tensor)::Savepoint,Savegradient
		real*8::savex,saveValue
		if(max_running.lt.2)then
			call writemess('ERROR in LinearSearch,input max_running should >=2',-1)
			call error_stop
		end if
		x1=0d0
		g1=inoutgra.dot.dir
		f1=outValue

		x2=x
		if(x2.le.OE%stop_step) then
			x2=OE%first_step_in_Linear_search
			x=x2
		end if

		call penalty_function(OE,f2,g2,gradient,NewPoint,Point,dir,x2)


		if(dabs(g2).le.OE%stop_gradient)then
			Point=NewPoint
			x=x2
			outValue=f2
			inoutgra=gradient
			return
		end if
		call LinearSearch_third_point(outx,f1,g1,f2,g2,x2)
		x3=outx

		call penalty_function(OE,f3,g3,gradient,NewPoint,Point,dir,x3)
		if(dabs(g3).le.OE%stop_gradient)then
			Point=NewPoint
			x=x3
			outValue=f3
			inoutgra=gradient
			return
		end if

		stopflag=.false.

		do i=3,max_running
			call LinearSearch_fouth_point(outx,stopflag,x1,f1,g1,x2,f2,g2,x3,f3,g3,x)
			if(isnan(outx).or.(outx.lt.0d0))then
				write(*,*)"ERROR in LinearSearch"
				write(*,*)"x1,x2"
				write(*,*)x1
				write(*,*)x2
				stop
			end if
			call penalty_function(OE,f3,g3,gradient,NewPoint,Point,dir,outx)
			x3=outx
			if((dabs(g3).le.OE%stop_gradient).or. stopflag) exit
		end do
		Point=NewPoint
		x=outx
		outValue=f3
		inoutgra=gradient
		return
	end subroutine


	!*************************************************************
	!*************************************************************
	!  linear search3 use the information of Energy only
	!      call OE%Energy_Function do not call OE%target_Function
	!      DO not use and modify inoutgra
	!*************************************************************

	subroutine LinearSearch4(OE,max_running,point,dir,x,outValue,inoutgra)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(inout)::x,outValue
		type(Tensor),intent(inout)::point
		type(Tensor),intent(inout)::inoutgra
		type(Tensor),intent(in)::dir
		integer,intent(in)::max_running
		logical::max_step_flag
		real*8::allx(3),allf(3),newx,newf,g0
		type(Tensor)::NewPoint,workingA,workingabc,workingy
		integer::ith,i
		if(max_running.lt.2)then
			call writemess('ERROR in LinearSearch,input max_running should >=2',-1)
			call error_stop
		end if
		max_step_flag=OE%max_step_in_Linear_search.gt.0
		allx(1)=0d0
		allf(1)=outValue

		allx(2)=x
		NewPoint=point+(allx(2)*dir)
		call OE%Energy_Function(allf(2),NewPoint)

		g0=inoutgra.dot.dir
		allx(3)=find_the_second_point(allx(2),allf(2),outValue,g0)
		NewPoint=point+(allx(3)*dir)
		call OE%Energy_Function(allf(3),NewPoint)

		if(max_running.eq.2)then
			if(allf(2).lt.allf(3))then
				x=allx(2)
				outValue=allf(2)
			else
				x=allx(3)
				outValue=allf(3)
			end if
			point=point+(x*dir)
			return
		end if

		do i=3,max_running
			call reorderdata(allx,allf)
			newx=find_the_next_point(allx,allf)
			NewPoint=point+(newx*dir)
			call OE%Energy_Function(newf,NewPoint)
			call select_new_point(allx,allf,newx,newf)
		end do

		ith=Find_min_out_index(allf)
		if(allx(ith).le.0)then
			allf(ith)=allf(ith+1)
			allx(ith)=allx(ith+1)
			ith=Find_min_out_index(allf)
		end if
		x=allx(ith)
		outValue=allf(ith)
		point=point+(x*dir)
		return
	end subroutine




	!use for the solve LLS with CG method
	!
	! AX=b, A is a n times n positive definite matrix
	!
	! Use for the two situations below
	!
	!
	!   A=Y' * G * Y,  Y' is the transposition of Y
	!
	!    Y is N time n matrix, N may not be equation to n
	!



	subroutine YGYp(Y,G,p,res)
		type(Tensor),target,intent(inout)::Y,G,p,res
		res=(G*(Y*p))*Y
		return

		res=Y*p
		res=G*res
		res=res*Y
		return
	end subroutine

	subroutine runCGLLS1(OE,inoutx,b,Y,G,memoryAp,memory_r,memory_p)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor)::inoutx,b,Y,G,memoryAp,memory_r,memory_p
		integer::i
		real*8::rr,alpha,beta
		call OE%InititalCGLLSStartPoint(inoutx,b)
		call YGYp(Y,G,inoutx,memory_r)
		memory_r=b-memory_r
		memory_p=memory_r
		do i=1,OE%CGLLSmaxRunning
			rr=memory_r%dnorm2()
			if(rr.le.OE%CGLLSerror)exit
			call YGYp(Y,G,memory_p,memoryAp)
			alpha=rr/(memory_p.dot.memoryAp)
			inoutx=inoutx+(alpha*memory_p)
			memory_r=memory_r-(alpha*memoryAp)
			beta=memory_r%dnorm2()/rr
			memory_p=memory_r+(beta*memory_p)
			
		end do
		if(OE%CGLLSerror.gt.0d0)then
			if(rr.gt.OE%CGLLSerror)then
				call writemess('*******  CGLLS WRONNING in OptimizationTools.f90 ******* ')
				call writemess('    The iteration may fail in convergence, the error is: ')
				call writemess(rr)
				call writemess('    The norm2 of gra is: ')
				call writemess(b%dnorm2())
				call writemess('********************************************************* ')
			end if
		end if
		return
	end subroutine

	subroutine runCGLLS2(OE,inoutx,b,memoryAp,memory_r,memory_p)
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor)::inoutx,b,memoryAp,memory_r,memory_p
		integer::i
		real*8::rr,alpha,beta
		call OE%InititalCGLLSStartPoint(inoutx,b)
		call OE%LLSMV(memory_r,inoutx)
		memory_r=b-memory_r
		memory_p=memory_r
		do i=1,OE%CGLLSmaxRunning
			rr=memory_r%dnorm2()
			if(rr.le.OE%CGLLSerror)exit
			call OE%LLSMV(memoryAp,memory_p)
			alpha=rr/(memory_p.dot.memoryAp)
			inoutx=inoutx+(alpha*memory_p)
			memory_r=memory_r-(alpha*memoryAp)
			beta=memory_r%dnorm2()/rr
			memory_p=memory_r+(beta*memory_p)
		end do
		if(OE%CGLLSerror.gt.0d0)then
			if(rr.gt.OE%CGLLSerror)then
				call writemess('*******  CGLLS WRONNING in OptimizationTools.f90 ******* ')
				call writemess('    The iteration may fail in convergence, the error is: ')
				call writemess(rr)
				call writemess('    The norm2 of gra is: ')
				call writemess(b%dnorm2())
				call writemess('********************************************************* ')
			end if
		end if
		return
	end subroutine

	subroutine runCGLLS_MMx(OE,inoutx,inb,memoryAp,memory_r,memory_p,b)! solve the x in A^T * A * x = A^T * b
		class(OptimEngine_structure), intent(inout) :: OE
		type(Tensor)::inoutx,inb,memoryAp,memory_r,memory_p,b
		integer::i
		real*8::norm2,rr,alpha,beta
		call OE%InititalCGLLSStartPoint(inoutx,inb)
		call OE%LLSMMV(b,inb,.true.)

		call OE%LLSMMV(memory_r,inoutx,.false.)
		memory_r=b-memory_r
		memory_p=memory_r
		do i=1,OE%CGLLSmaxRunning
			rr=memory_r%dnorm2()
			call OE%LLSMMV(memoryAp,memory_p,.false.)
			alpha=rr/(memory_p.dot.memoryAp)
			inoutx=inoutx+(alpha*memory_p)
			memory_r=memory_r-(alpha*memoryAp)
			beta=memory_r%dnorm2()/rr
			memory_p=memory_r+(beta*memory_p)
			norm2=memory_p%dnorm2()
			if(norm2.le.OE%CGLLSerror)exit
		end do
		if(OE%CGLLSerror.gt.0d0)then
			if(rr.gt.OE%CGLLSerror)then
				call writemess('*******  CGLLS WRONNING in OptimizationTools.f90 ******* ')
				call writemess('    The iteration may fail in convergence, the error is: ')
				call writemess(rr)
				call writemess('    The norm2 of gra is: ')
				call writemess(b%dnorm2())
				call writemess('********************************************************* ')
			end if
		end if
		return
	end subroutine

	subroutine defaultInititalCGLLSStartPoint(OE,inoutx,b)
		class(OptimEngine_structure),target, intent(inout) :: OE
		type(Tensor),target, intent(inout) :: inoutx
		type(Tensor),target, intent(in)::b
		call inoutx%allocate([b%getTotalData()],b%getType())
		call inoutx%zero()
		return
	end subroutine
end module
