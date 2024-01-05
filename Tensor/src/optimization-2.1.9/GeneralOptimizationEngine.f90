module General_Optimization_Tools
	use General_Optimization_Element
	use GeneralLineSearchEnergy
	implicit none
	private

	integer,private,parameter::max_stop_counter=5
	integer,private,parameter::LBFGS_default_memory_length=10
	real*8,private::zero_gradient=1d-16

	public::OptimElement_structure
	public::OptimEngine_structure
	type  :: OptimEngine_structure
		real*8::stop_gradient=-1 !it will stop if norm(gradient)<stop_gradient
		real*8::stop_step=1d-10 !it will stop if x<stop_step for max_stop_counter times, x is the search step in every step
		integer::CG_direction_Flag=1
		integer::CGLLSmaxRunning=300
		real*8::CGLLSerror=1d-16
		real*8::line_search_min_step=-1d0
		real*8::first_step_in_Line_search=0.1d0
		real*8::search_restart_step=-1d-5
		real*8::line_search_restart_max_num=-1
		character(len=50)::method='null'

		!use for adam
		real*8::Adam_rho1=0.9
		real*8::Adam_rho2=0.999
		real*8::Adam_rho1t=0.9
		real*8::Adam_rho2t=0.999
		integer::Adam_resetStep=-1
		integer::Adam_resetcounter=0

		!use for LBFGS method
		! st,yt,length,FullFlag,endindex
		class(OptimElement_structure),allocatable::st(:)
		class(OptimElement_structure),allocatable::yt(:)
		integer::length=0
		logical::FullFlag=.false.
		integer::endindex=0
		!use for NL method
		! st,yt,length,FullFlag,endindex...
		class(OptimElement_structure),allocatable::nl_YTy(:)    	!YparX^T * y 
		class(OptimElement_structure),allocatable::nl_YTgama(:)    	!YparX^T * gama 
		class(OptimElement_structure),allocatable::nl_yg(:)    		!ygradient in y space
		class(OptimElement_structure),allocatable::nl_yx(:)    		!ypoint in y space
		class(OptimElement_structure),allocatable::nl_yo(:)    		!yother in y space
		real*8,allocatable::nl_sy(:,:)								!sy(i,j)=si.dot.yj=(nl_yx(i+1)-nl_yx(i))*(nl_yg(j+1)-nl_yg(j))
		real*8,allocatable::nl_sgama(:,:)							!sgama(i,j)=si.dot.gamaj
		real*8,allocatable::nl_yxg(:,:)								!yxg(i,j)=yxi.dot.ygj
		real*8,allocatable::nl_yxx(:,:)								!yxx(i,j)=yxi.dot.yxj
		real*8::nl_episilo=0
		real*8::nl_diag_value=0
		integer::nl_length=0		
		integer::flag1=0
		integer::flag2=0
		integer::flag3=0
		integer::saveflag3=0
		real*8::line_search_final_step_ratio=1d0
		procedure(nl_YTYV_interface),pointer::NLYTYV=>null()
		procedure(nl_YTV_interface),pointer::NLYTV=>null()
		procedure(nl_Y_point_gradient_interface),pointer::NLYPG=>null()
		procedure(YDotY_interface),pointer::YDotY=>default_YDotY
		procedure(LparY_Y_interface),pointer::LparY_Y=>default_LparY_Y
		procedure(Y_LparY_interface),pointer::Y_LparY=>default_Y_LparY
		procedure(YparX_Y_interface),pointer::YparX_Y=>default_YparX_Y
		procedure(YparX_LparY_interface),pointer::YparX_LparY=>default_YparX_LparY
		procedure(store_other_interface),pointer::store_other=>default_store_other

		class(OptimElement_structure),private,allocatable::workingMemory(:)
		procedure(LinearSearch_subroutine),pointer::LinearSearch=>LinearSearch4
		procedure(Step1Subroutine_interface),pointer::inStep1=>null()
		procedure(Step2Subroutine_interface),pointer::inStep2=>null()
		procedure(Step3Subroutine_interface),pointer::inStep3=>null()
		procedure(Step4Subroutine_interface),pointer::inStep4=>null()
		procedure(BeforeStepSubroutine_interface),pointer::BeforeStep=>null()
		procedure(EndStepSubroutine_interface),pointer::EndStep=>null()
		procedure(matrix_time_vector_interface),pointer::LLSMV=>null()
		procedure(externalCGLLS_interface),pointer::externalCGLLS=>null()
		procedure(RescalDirection_interface),pointer::RescalDirection=>DefaultRescalDirection
		procedure(AllocateMemory_interface),pointer::AllocateWorkingMemory=>DefaultAllocateWorkingMemory
		procedure(set_DefaultDirection_interface),pointer::set_Direction=>set_DefaultDirection
		procedure(targetFunc_interface), pointer :: target_Function=>Defaulttarget_Function
		procedure(InititalCGLLSStartPointInterface),pointer::InititalCGLLSStartPoint=>defaultInititalCGLLSStartPointFromZero
		procedure(Energy_FunctionInterface),pointer::Energy_Function=>defaultEnergy_Function
		procedure(PointToOptimData_interface),pointer::PointToOptimData=>DefaultPointToOptimData
		procedure(restart_line_search_subroutone_interface),pointer::restart_line_search_subroutone=>defult_restart_line_search_subroutone
	contains
		procedure,public::deallocate=>deallocate_engine
		procedure,public::pointMemory
		procedure::check_stop
		procedure::set_method1,set_method2
		generic,public::set_method=>set_method1,set_method2
		procedure::runCGLLS
		procedure,public::set_MV_func
		procedure,public::set_CGLLS_parameter
		procedure::setInititalCGLLSStartPoint1,setInititalCGLLSStartPoint2
		generic,public::setInititalCGLLSStartPoint=>setInititalCGLLSStartPoint1,setInititalCGLLSStartPoint2
		procedure,public::set_CG_direction_flag
		procedure,public::set_target_function
		procedure,public::set_Energy_Function
		procedure,public::set_stop_error
		procedure,public::set_stop_Step
		procedure,public::set_stop_gradient=>set_stop_error
		procedure,public::set_line_search_restart_step
		procedure,public::set_line_search_restart_max_num
		procedure,public::set_line_search_final_step_ratio
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

		!use for LBFGS method
		procedure::pointNewElement
		generic,public::NewElement=>pointNewElement
		procedure::element_i
		generic,public::i=>element_i
		procedure,public::allocate_LBFGS_memory
		procedure,public::deallocate_LBFGS_memory
		procedure,public::LBFGS_data_length
		procedure,public::reset_LBFGS_Endpoint
		procedure,public::LBFGS_data_Size

		!use for NL method
		procedure,public::allocate_NL_memory
		procedure,public::deallocate_NL_memory
		procedure,public::NL_data_Size
		procedure,public::set_NL_length
		procedure,public::set_NL_episilo
		procedure,public::set_NL_function
		procedure,public::NLMV
		procedure,public::set_NL_Idetity_function
		procedure,public::set_NL_parameter
		procedure,public::set_LparY_Y
		procedure,public::set_Y_LparY
		procedure,public::set_YparX_Y
		procedure,public::set_YparX_LparY
		procedure,public::set_YDotY
		procedure,public::set_store_other
		procedure,public::default_LparY_Y
		procedure,public::default_Y_LparY
		procedure,public::default_YparX_Y
		procedure,public::default_YparX_LparY
		procedure,public::default_YDotY
		procedure,public::default_store_other


		procedure::Optimization1,Optimization2,Optimization3,Optimization4,Optimization5
		generic,public::Optim=>Optimization1,Optimization2,Optimization3,Optimization4,Optimization5
		procedure,public::set_PointToOptimData_function
		procedure,public::set_line_search_min_step
		procedure,public::set_first_step_in_Line_search
		procedure,public::set_line_search_type
		procedure,public::set_LLS_func
		procedure,public::set_line_search_function
		generic,public::set_line_search=>set_line_search_type,set_line_search_function
		procedure,public::set_restart_line_search_subroutone
		procedure,public::set_LBFGS_length
	end type OptimEngine_structure

	interface
		subroutine targetFunc_interface(A,outVal,outGradient,point)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: A
			real*8,target,intent(inout)::outVal
			class(OptimElement_structure),target,intent(inout)::outGradient
			class(OptimElement_structure),target,intent(in)::point
		end subroutine targetFunc_interface
	end interface

	interface
		subroutine Energy_FunctionInterface(A,outVal,point)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: A
			real*8,target,intent(inout)::outVal
			class(OptimElement_structure),target,intent(in)::point
		end subroutine Energy_FunctionInterface
	end interface

	interface
		subroutine AllocateMemory_interface(OE,Gradient,direction)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure), intent(inout) :: OE
			class(OptimElement_structure),pointer,intent(inout)::Gradient,direction
		end subroutine AllocateMemory_interface
	end interface

	interface
		subroutine Step1Subroutine_interface(OE,Value,direction,Gradient,point,ith,t)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Value
			class(OptimElement_structure),target,intent(inout)::Gradient
			class(OptimElement_structure),target,intent(inout)::point
			class(OptimElement_structure),target,intent(inout)::direction
			integer,target,intent(in)::ith
			real*8,target::t
		end subroutine Step1Subroutine_interface
	end interface

	interface
		subroutine Step2Subroutine_interface(OE,Value,direction,Gradient,point,ith,t)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Value
			class(OptimElement_structure),target,intent(inout)::Gradient
			class(OptimElement_structure),target,intent(inout)::point
			class(OptimElement_structure),target,intent(inout)::direction
			integer,target,intent(in)::ith
			real*8,target::t
		end subroutine Step2Subroutine_interface
	end interface

	interface
		subroutine Step3Subroutine_interface(OE,Value,direction,Gradient,point,ith,t)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Value
			class(OptimElement_structure),target,intent(inout)::Gradient
			class(OptimElement_structure),target,intent(inout)::point
			class(OptimElement_structure),target,intent(inout)::direction
			integer,target,intent(in)::ith
			real*8,target::t
		end subroutine Step3Subroutine_interface
	end interface

	interface
		subroutine Step4Subroutine_interface(OE,Value,direction,Gradient,point,ith,t)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Value
			class(OptimElement_structure),target,intent(inout)::Gradient
			class(OptimElement_structure),target,intent(inout)::point
			class(OptimElement_structure),target,intent(inout)::direction
			integer,target,intent(in)::ith
			real*8,target::t
		end subroutine Step4Subroutine_interface
	end interface

	interface
		subroutine BeforeStepSubroutine_interface(OE,Value,Gradient,point)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Value
			class(OptimElement_structure),target,intent(inout)::Gradient
			class(OptimElement_structure),target,intent(inout)::point
		end subroutine BeforeStepSubroutine_interface
	end interface

	interface
		subroutine EndStepSubroutine_interface(OE,Value,Gradient,point)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Value
			class(OptimElement_structure),target,intent(inout)::Gradient
			class(OptimElement_structure),target,intent(inout)::point
		end subroutine EndStepSubroutine_interface
	end interface

	interface
		subroutine RescalDirection_interface(OE,Dir)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure), intent(inout) :: OE
			class(OptimElement_structure),intent(inout)::Dir
		end subroutine RescalDirection_interface
	end interface

	interface
		subroutine set_DefaultDirection_interface(OE,direction,Gradient,point)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure), intent(inout) :: OE
			class(OptimElement_structure),intent(in)::point,Gradient
			class(OptimElement_structure),intent(inout)::direction
		end subroutine set_DefaultDirection_interface
	end interface

	interface
		subroutine InititalCGLLSStartPointInterface(OE,inoutx,b)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			class(OptimElement_structure),target, intent(inout) :: inoutx
			class(OptimElement_structure),target, intent(in)::b
		end subroutine InititalCGLLSStartPointInterface
	end interface

	interface
		subroutine matrix_time_vector_interface(OE,outAp,inp)!outAp=A*inp
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			class(OptimElement_structure),target,intent(inout)::outAp
			class(OptimElement_structure),target,intent(in)::inp
		end subroutine matrix_time_vector_interface
	end interface


	interface
		subroutine PointToOptimData_interface(A,point)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: A
			class(OptimElement_structure),pointer,intent(inout)::point
		end subroutine PointToOptimData_interface
	end interface

	interface
		subroutine LinearSearch_subroutine(OE,max_running,point,dir,x,outValue,inoutgra)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,intent(inout),target::x,outValue
			class(OptimElement_structure),target,intent(inout)::point
			class(OptimElement_structure),target,intent(inout)::inoutgra
			class(OptimElement_structure),target,intent(in)::dir
			integer,target,intent(in)::max_running
		end subroutine LinearSearch_subroutine
	end interface

	interface
		subroutine restart_line_search_subroutone_interface(OE,max_running,point,dir,x,outValue,inoutgra)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,intent(inout),target::x,outValue
			class(OptimElement_structure),target,intent(inout)::point
			class(OptimElement_structure),target,intent(inout)::inoutgra
			class(OptimElement_structure),target,intent(in)::dir
			integer,target,intent(in)::max_running
		end subroutine restart_line_search_subroutone_interface
	end interface

	interface
		subroutine externalCGLLS_interface(OE,direction,Gradient)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure), intent(inout) :: OE
			class(OptimElement_structure),intent(inout)::direction
			class(OptimElement_structure),intent(in)::Gradient
		end subroutine externalCGLLS_interface
	end interface

	interface
		subroutine nl_YTYV_interface(OE,inoutp,inp) !outAp=A^T*A*inp
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			class(OptimElement_structure),target,intent(inout)::inoutp
			class(OptimElement_structure),target,intent(in)::inp
		end subroutine nl_YTYV_interface
	end interface

	interface
		subroutine nl_YTV_interface(OE,inoutp,inp) !outAp=A^T*inp
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			class(OptimElement_structure),target,intent(inout)::inoutp
			class(OptimElement_structure),target,intent(in)::inp
		end subroutine nl_YTV_interface
	end interface
	
	interface
		subroutine nl_Y_point_gradient_interface(OE,yPoint,yGradient,point) 
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			class(OptimElement_structure),target,intent(inout)::yPoint
			class(OptimElement_structure),target,intent(inout)::yGradient
			class(OptimElement_structure),target,intent(in)::point
		end subroutine nl_Y_point_gradient_interface
	end interface

	interface
		subroutine YDotY_interface(OE,Res,point,gradient,B,Gb,Ob)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Res
			class(OptimElement_structure),target,intent(in)::point
			class(OptimElement_structure),target,intent(in)::gradient
			class(OptimElement_structure),target,intent(in)::B
			class(OptimElement_structure),target,intent(in)::Gb
			class(OptimElement_structure),target,intent(in)::Ob
		end subroutine YDotY_interface
	end interface

	interface
		subroutine LparY_Y_interface(OE,Res,point,gradient,B,Gb,Ob)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Res
			class(OptimElement_structure),target,intent(in)::point
			class(OptimElement_structure),target,intent(in)::gradient
			class(OptimElement_structure),target,intent(in)::B
			class(OptimElement_structure),target,intent(in)::Gb
			class(OptimElement_structure),target,intent(in)::Ob
		end subroutine LparY_Y_interface
	end interface

	interface
		subroutine Y_LparY_interface(OE,Res,point,gradient,B,Gb,Ob)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			real*8,target,intent(inout)::Res
			class(OptimElement_structure),target,intent(in)::point
			class(OptimElement_structure),target,intent(in)::gradient
			class(OptimElement_structure),target,intent(in)::B
			class(OptimElement_structure),target,intent(in)::Gb
			class(OptimElement_structure),target,intent(in)::Ob
		end subroutine Y_LparY_interface
	end interface

	interface
		subroutine YparX_Y_interface(OE,Res,point,gradient,B,Gb,Ob)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			class(OptimElement_structure),target,intent(inout)::Res
			class(OptimElement_structure),target,intent(in)::point
			class(OptimElement_structure),target,intent(in)::gradient
			class(OptimElement_structure),target,intent(in)::B
			class(OptimElement_structure),target,intent(in)::Gb
			class(OptimElement_structure),target,intent(in)::Ob
		end subroutine YparX_Y_interface
	end interface

	interface
		subroutine YparX_LparY_interface(OE,Res,point,gradient,B,Gb,Ob)
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			class(OptimElement_structure),target,intent(inout)::Res
			class(OptimElement_structure),target,intent(in)::point
			class(OptimElement_structure),target,intent(in)::gradient
			class(OptimElement_structure),target,intent(in)::B
			class(OptimElement_structure),target,intent(in)::Gb
			class(OptimElement_structure),target,intent(in)::Ob
		end subroutine YparX_LparY_interface
	end interface

	interface
		subroutine store_other_interface(OE,Other) 
			use General_Optimization_Element
			import :: OptimEngine_structure
			class(OptimEngine_structure),target, intent(inout) :: OE
			class(OptimElement_structure),target,intent(inout)::Other
		end subroutine store_other_interface
	end interface

contains
	subroutine deallocate_engine(OE)
		class(OptimEngine_structure),target, intent(inout) :: OE
		OE%stop_gradient=-1 !it will stop if norm(gradient)<stop_gradient
		OE%stop_step=1d-10 !it will stop if x<stop_step for max_stop_counter times, x is the search step in every step
		OE%CG_direction_Flag=1
		OE%CGLLSmaxRunning=300
		OE%CGLLSerror=1d-16
		OE%line_search_min_step=-1d0
		OE%first_step_in_Line_search=0.1d0
		OE%search_restart_step=-1d-5
		OE%line_search_restart_max_num=-1
		OE%method='null'

		!use for adam
		OE%Adam_rho1=0.9
		OE%Adam_rho2=0.999
		OE%Adam_rho1t=0.9
		OE%Adam_rho2t=0.999
		OE%Adam_resetStep=-1
		OE%Adam_resetcounter=0

		!use for LBFGS method
		! st,yt,length,FullFlag,endindex
		if(allocated(OE%st))then
			deallocate(OE%st)
			deallocate(OE%yt)
		end if
		OE%length=0
		OE%FullFlag=.false.
		OE%endindex=0
		!use for NL method
		! st,yt,length,FullFlag,endindex...
		if(allocated(OE%nl_YTy))then
			deallocate(OE%nl_YTy)
			deallocate(OE%nl_YTgama)
			deallocate(OE%nl_yg)
			deallocate(OE%nl_yx)
			deallocate(OE%nl_yo)
			deallocate(OE%nl_sy)
			deallocate(OE%nl_sgama)
			deallocate(OE%nl_yxg)
			deallocate(OE%nl_yxx)
		end if
		OE%nl_episilo=0
		OE%nl_diag_value=0
		OE%nl_length=0		
		OE%flag1=0
		OE%flag2=0
		OE%flag3=0
		OE%saveflag3=0
		OE%line_search_final_step_ratio=1d0
		

		OE%NLYTYV=>null()
		OE%NLYTV=>null()
		OE%NLYPG=>null()
		OE%YDotY=>default_YDotY
		OE%LparY_Y=>default_LparY_Y
		OE%Y_LparY=>default_Y_LparY
		OE%YparX_Y=>default_YparX_Y
		OE%YparX_LparY=>default_YparX_LparY
		OE%store_other=>default_store_other
		if(allocated(OE%workingMemory))deallocate(OE%workingMemory)
		OE%LinearSearch=>LinearSearch4
		OE%inStep1=>null()
		OE%inStep2=>null()
		OE%inStep3=>null()
		OE%inStep4=>null()
		OE%BeforeStep=>null()
		OE%EndStep=>null()
		OE%LLSMV=>null()
		OE%externalCGLLS=>null()
		OE%RescalDirection=>DefaultRescalDirection
		OE%AllocateWorkingMemory=>DefaultAllocateWorkingMemory
		OE%set_Direction=>set_DefaultDirection
		OE%target_Function=>Defaulttarget_Function
		OE%InititalCGLLSStartPoint=>defaultInititalCGLLSStartPointFromZero
		OE%Energy_Function=>defaultEnergy_Function
		OE%PointToOptimData=>DefaultPointToOptimData
		OE%restart_line_search_subroutone=>defult_restart_line_search_subroutone
		return
	end subroutine
	subroutine default_LparY_Y(OE,Res,point,gradient,B,Gb,Ob)
		class(OptimEngine_structure),target, intent(inout) :: OE
		real*8,target,intent(inout)::Res
		class(OptimElement_structure),target,intent(in)::point
		class(OptimElement_structure),target,intent(in)::gradient
		class(OptimElement_structure),target,intent(in)::B
		class(OptimElement_structure),target,intent(in)::Gb
		class(OptimElement_structure),target,intent(in)::Ob
		call printFunc(' DO NOT set the LparY_Y subroutine yet')
		call ErrorFunc
		return
	end subroutine
	subroutine default_Y_LparY(OE,Res,point,gradient,B,Gb,Ob)
		class(OptimEngine_structure),target, intent(inout) :: OE
		real*8,target,intent(inout)::Res
		class(OptimElement_structure),target,intent(in)::point
		class(OptimElement_structure),target,intent(in)::gradient
		class(OptimElement_structure),target,intent(in)::B
		class(OptimElement_structure),target,intent(in)::Gb
		class(OptimElement_structure),target,intent(in)::Ob
		call printFunc(' DO NOT set the Y_LparY subroutine yet')
		call ErrorFunc
		return
	end subroutine
	subroutine default_YparX_Y(OE,Res,point,gradient,B,Gb,Ob)
		class(OptimEngine_structure),target, intent(inout) :: OE
		class(OptimElement_structure),target,intent(inout)::Res
		class(OptimElement_structure),target,intent(in)::point
		class(OptimElement_structure),target,intent(in)::gradient
		class(OptimElement_structure),target,intent(in)::B
		class(OptimElement_structure),target,intent(in)::Gb
		class(OptimElement_structure),target,intent(in)::Ob
		call printFunc(' DO NOT set the YparX_Y subroutine yet')
		call ErrorFunc
		return
	end subroutine
	subroutine default_YparX_LparY(OE,Res,point,gradient,B,Gb,Ob)
		class(OptimEngine_structure),target, intent(inout) :: OE
		class(OptimElement_structure),target,intent(inout)::Res
		class(OptimElement_structure),target,intent(in)::point
		class(OptimElement_structure),target,intent(in)::gradient
		class(OptimElement_structure),target,intent(in)::B
		class(OptimElement_structure),target,intent(in)::Gb
		class(OptimElement_structure),target,intent(in)::Ob
		call printFunc(' DO NOT set the YparX_LparY subroutine yet')
		call ErrorFunc
		return
	end subroutine
	subroutine default_YDotY(OE,Res,point,gradient,B,Gb,Ob)
		class(OptimEngine_structure),target, intent(inout) :: OE
		real*8,target,intent(inout)::Res
		class(OptimElement_structure),target,intent(in)::point
		class(OptimElement_structure),target,intent(in)::gradient
		class(OptimElement_structure),target,intent(in)::B
		class(OptimElement_structure),target,intent(in)::Gb
		class(OptimElement_structure),target,intent(in)::Ob
		call printFunc(' DO NOT set the YDotY subroutine yet')
		call ErrorFunc
		return
	end subroutine
	subroutine default_store_other(OE,Other) 
		class(OptimEngine_structure),target, intent(inout) :: OE
		class(OptimElement_structure),target,intent(inout)::Other
		call printFunc(' DO NOT set the store_other subroutine yet')
		call ErrorFunc
	end subroutine 
	subroutine set_LparY_Y(OE,Func)
		class(OptimEngine_structure),target, intent(inout) :: OE
		procedure(LparY_Y_interface)::Func
		call printFunc(' set the LparY_Y subroutine')
		OE%LparY_Y=>Func
		return
	end subroutine
	subroutine set_Y_LparY(OE,Func)
		class(OptimEngine_structure),target, intent(inout) :: OE
		procedure(Y_LparY_interface)::Func
		call printFunc(' set the Y_LparY subroutine')
		OE%Y_LparY=>Func
		return
	end subroutine
	subroutine set_YparX_Y(OE,Func)
		class(OptimEngine_structure),target, intent(inout) :: OE
		procedure(YparX_Y_interface)::Func
		call printFunc(' set the YparX_Y subroutine')
		OE%YparX_Y=>Func
		return
	end subroutine
	subroutine set_YparX_LparY(OE,Func)
		class(OptimEngine_structure),target, intent(inout) :: OE
		procedure(YparX_LparY_interface)::Func
		call printFunc(' set the YparX_LparY subroutine')
		OE%YparX_LparY=>Func
		return
	end subroutine
	subroutine set_YDotY(OE,Func)
		class(OptimEngine_structure),target, intent(inout) :: OE
		procedure(YDotY_interface)::Func
		call printFunc(' set the YDotY subroutine')
		OE%YDotY=>Func
		return
	end subroutine
	subroutine set_store_other(OE,Func)
		class(OptimEngine_structure),target, intent(inout) :: OE
		procedure(store_other_interface)::Func
		call printFunc(' set the store_other subroutine')
		OE%store_other=>Func
		return
	end subroutine


	subroutine defaultEnergy_Function(OE,outVal,point)
		class(OptimEngine_structure),target, intent(inout) :: OE
		real*8,target,intent(inout)::outVal
		class(OptimElement_structure),target,intent(in)::point
		class(OptimElement_structure),allocatable,target::outGradient
		allocate(outGradient,source=point)
		call OE%target_Function(outVal,outGradient,point)
		return
	end subroutine

	subroutine set_Energy_Function(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(Energy_FunctionInterface)::Func
		call printFunc('Set the Energy_Function as external function')
		call printFunc('Set LinearSearch type as type4  in OptimEngine_structure') 
		call printFunc('Use only the Energy in the Searching') 
		OE%LinearSearch=>LinearSearch4
		OE%Energy_Function=>Func
		return
	end subroutine

	subroutine set_target_function(OE,Func)
		class(OptimEngine_structure),intent(inout)::OE
		procedure(targetFunc_interface)::Func
		call printFunc('Set the target_function in OptimEngine_structure')
		OE%target_Function=>Func
		return
	end subroutine

	subroutine set_PointToOptimData_function(OE,Func)
		class(OptimEngine_structure),intent(inout)::OE
		procedure(PointToOptimData_interface)::Func
		call printFunc('Set the PointToOptimData in OptimEngine_structure')
		OE%PointToOptimData=>Func
		return
	end subroutine

	!***************************************************
	!           initial function for OptimEngine_structure
	!***************************************************

	subroutine set_MV_func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(matrix_time_vector_interface)::Func
		call printFunc('Set the matrix times vector function used in Natural Gradient Method in OptimEngine_structure')
		call printFunc(' The metric in the problem is M=Y^T * G * Y')
		call printFunc(' Function is to output the p, p=M*x')
		OE%LLSMV=>Func
		return
	end subroutine

	subroutine set_LLS_func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(externalCGLLS_interface)::Func
		call printFunc('Set the externalCGLLS function used in Natural Gradient Method in OptimEngine_structure')
		call printFunc(' The metric in the problem is M=Y^T * G * Y')
		call printFunc(' Function is to output the d, the root the equation M*d=g, where g is the gradient')
		OE%externalCGLLS=>Func
		return
	end subroutine

	subroutine set_CGLLS_parameter(OE,CGLLSmaxRunning,CGLLSerror)
		class(OptimEngine_structure), intent(inout) :: OE
		integer,intent(in)::CGLLSmaxRunning
		real*8,intent(in)::CGLLSerror
		call printFunc('Set the parameter in CGLLS, which will use only in Natural Gradient Method')
		call printFunc(' maxRunning for the root of Ax=b is CGLLSmaxRunning=')
		call printFunc(CGLLSmaxRunning)
		call printFunc(' step error for the root of Ax=b is CGLLSerror=')
		call printFunc(CGLLSerror)
		OE%CGLLSmaxRunning=CGLLSmaxRunning
		OE%CGLLSerror=CGLLSerror
		return
	end subroutine

	subroutine setInititalCGLLSStartPoint1(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(InititalCGLLSStartPointInterface)::Func
		call printFunc('Set the InititalCGLLSStartPoint function, it will use for CGLLS')
		OE%InititalCGLLSStartPoint=>Func
	end subroutine
	
	subroutine setInititalCGLLSStartPoint2(OE,Flag)
		class(OptimEngine_structure), intent(inout) :: OE
		character(len=*)::Flag
		call printFunc('Set the InititalCGLLSStartPoint function, it will use for CGLLS')
		select case(Flag)
			case ('zero')
				call printFunc('InititalCGLLSStartPoint is 0')
				OE%InititalCGLLSStartPoint=>defaultInititalCGLLSStartPointFromZero
			case ('gradient')
				call printFunc('InititalCGLLSStartPoint is gradient')
				OE%InititalCGLLSStartPoint=>defaultInititalCGLLSStartPointFromGra
			case ('zero_first')
				call printFunc('the first InititalCGLLSStartPoint is 0, and then it use last output point')
				OE%InititalCGLLSStartPoint=>defaultInititalCGLLSStartPointFromZero2
			case ('gradient_first')
				call printFunc('the first InititalCGLLSStartPoint is gradient, and then it use last output point')
				OE%InititalCGLLSStartPoint=>defaultInititalCGLLSStartPointFromGra2
			case default
				call printFunc('ERROR InititalCGLLSStartPoint Flag')
				call ErrorFunc
		end select
		return
	end subroutine

	subroutine set_CG_direction_flag(OE,Flag)
		class(OptimEngine_structure), intent(inout) :: OE
		integer,intent(in)::Flag
		call printFunc('Set the CG_direction_flag((only use in method=CG) in OptimEngine_structure')
		call printFunc(' CG_direction_flag=')
		call printFunc(Flag)
		call printFunc('     1: Crowder_Wolfe')
		call printFunc('     2: Fletcher_Reeves')
		call printFunc('     3: Dixon')
		call printFunc('     4: Polak_Ribiere_Polyak')
		call printFunc('     5: Dai_Yuan')
		OE%CG_direction_Flag=Flag
		return
	end subroutine

	subroutine set_stop_error(OE,error)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::error
		call printFunc('Set the step error(stop_gradient) in OptimEngine_structure')
		OE%stop_gradient=error
		return
	end subroutine
	subroutine set_stop_Step(OE,step)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::step
		call printFunc('Set the stop length of the step(stop_step) in OptimEngine_structure')
		OE%stop_step=step
		return
	end subroutine
	subroutine set_line_search_restart_step(OE,step)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::step
		call printFunc('Set the search_restart_step  in OptimEngine_structure,search_restart_step=')
		call printFunc(step)
		call printFunc(' If line search step < line_search_min_step, the program will reset the step as')
		call printFunc('  step = search_restart_step')
		OE%search_restart_step=step
		return
	end subroutine
	
	subroutine set_line_search_min_step(OE,step)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::step
		call printFunc('Set the line_search_min_step in OptimEngine_structure, line_search_min_step=')
		call printFunc(step)
		call printFunc(' If line search step < line_search_min_step, and search_restart_step>0')
		call printFunc(' the program will reset the step as')
		call printFunc('  step = search_restart_step')
		OE%line_search_min_step=step
		return
	end subroutine
	subroutine set_first_step_in_Line_search(OE,step)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::step
		call printFunc('Set the first step in Linear search in OptimEngine_structure')
		OE%first_step_in_Line_search=step
		return
	end subroutine
	subroutine set_line_search_restart_max_num(OE,line_search_restart_max_num)
		class(OptimEngine_structure), intent(inout) :: OE
		integer,intent(in)::line_search_restart_max_num
		call printFunc('Set the line_search_restart_max_num in OptimEngine_structure, line_search_restart_max_num=')
		call printFunc(line_search_restart_max_num)
		OE%line_search_restart_max_num=line_search_restart_max_num
		return
	end subroutine
	

	function check_stop(OE,value,gradient,step)
		logical::check_stop
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),intent(inout)::gradient
		real*8,intent(in)::step,value
		real*8::norm,tmp
		integer,save::stop_counter=0

		check_stop=.false.
		if(OE%stop_gradient.gt.0) then
			call norm2Func(gradient,norm)
			if(abs(value).gt.1d0)then
				check_stop=abs(norm/value).le.OE%stop_gradient
				if(check_stop)then
					call printFunc('Search is going to stop, norm of the |gradient/value|=')
					tmp=abs(norm/value)
					call printFunc(tmp)
					call printFunc('The giving stop error=')
					call printFunc(OE%stop_gradient)
					return
				end if
			else
				check_stop=norm.le.OE%stop_gradient
				if(check_stop)then
					call printFunc('Search is going to stop, norm of the gradient=')
					call printFunc(norm)
					call printFunc('The giving stop error=')
					call printFunc(OE%stop_gradient)
					return
				end if
			end if
		end if

		if(step.le.OE%stop_step)then
			stop_counter=stop_counter+1
			if(stop_counter.ge.max_stop_counter)then
				check_stop=.true.
				call printFunc('Search is going to stop, the search step=')
				call printFunc(step)
				return
			end if
		else
			stop_counter=0
		end if
		return
	end function

	subroutine set_line_search_type(OE,flag)
		class(OptimEngine_structure), intent(inout) :: OE
		integer,intent(in)::flag
		if(flag.eq.4)then
			call printFunc('Set LinearSearch type as type4  in OptimEngine_structure') 
			call printFunc('Use only the Energy in the Searching') 
			OE%LinearSearch=>LinearSearch4
			return
		end if
		call printFunc('No such case')
		call ErrorFunc
	end subroutine

	subroutine set_line_search_function(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(LinearSearch_subroutine)::Func
		call printFunc('Set the LinearSearch subroutine as the external function in OptimEngine_structure')
		OE%LinearSearch=>Func
		return
	end subroutine

	subroutine set_restart_line_search_subroutone(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(LinearSearch_subroutine)::Func
		call printFunc('Set the restart_line_search_subroutone as the external function in OptimEngine_structure')
		OE%restart_line_search_subroutone=>Func
		return
	end subroutine

	subroutine set_line_search_final_step_ratio(OE,ratio)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::ratio
		call printFunc('Set line_search_final_step_ratio as ratio=')
		call printFunc(  ratio )
		call printFunc('The ratio is used in line search, reset the output step as x=x*ratio')
		OE%line_search_final_step_ratio=ratio
		return
	end subroutine


	subroutine pointMemory(OE,p,ith)
		class(OptimEngine_structure),target, intent(inout) :: OE
		class(OptimElement_structure),pointer::p
		integer,intent(in)::ith
		if(.not.allocated(OE%workingMemory))then
			call printFunc('ERROR in pointMemory, DO NOT allocate memory yet')
			call ErrorFunc
		end if
		if(ith.gt.size(OE%workingMemory))then
			call printFunc('ERROR in pointMemory, ith> size(memory)')
			call ErrorFunc
		end if
		p=>OE%workingMemory(ith)
		return
	end subroutine

	subroutine set_method1(OE,method)
		class(OptimEngine_structure), intent(inout) :: OE
		character(len=*),intent(in)::method
		call printFunc('Set the optimal method in OptimEngine_structure as:')
		call printFunc(method)
		select case(method)
			case ('GM')
				OE%RescalDirection=>DefaultRescalDirection
				OE%AllocateWorkingMemory=>GMMemory
				OE%set_Direction=>set_GMDirection
			case ('RGM')
				OE%RescalDirection=>RandomDirection
				OE%AllocateWorkingMemory=>GMMemory
				OE%set_Direction=>set_GMDirection
			case ('adam')
				call printFunc('adam will not rescal the searching direction')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>AdamMemory
				OE%set_Direction=>set_AdamDirection
			case ('CG')
				OE%RescalDirection=>DefaultRescalDirection
				OE%AllocateWorkingMemory=>CGMemory
				OE%set_Direction=>set_CGDirection
			case ('RCG')
				OE%RescalDirection=>RandomDirection
				OE%AllocateWorkingMemory=>CGMemory
				OE%set_Direction=>set_CGDirection
			case ('LBFGS')
				OE%RescalDirection=>DefaultRescalDirection
				OE%AllocateWorkingMemory=>LBFGSMemory
				OE%set_Direction=>set_LBFGSDirection
			case ('RLBFGS')
				OE%RescalDirection=>RandomDirection
				OE%AllocateWorkingMemory=>LBFGSMemory
				OE%set_Direction=>set_LBFGSDirection
			case ('SJ')
				OE%RescalDirection=>DefaultRescalDirection
				OE%AllocateWorkingMemory=>SJ_Memory
				OE%set_Direction=>set_SJ_direction
			case ('NG')
				OE%RescalDirection=>DefaultRescalDirection
				OE%AllocateWorkingMemory=>SJ_Memory
				OE%set_Direction=>set_SJ_direction
			case ('RSJ')
				OE%RescalDirection=>RandomDirection
				OE%AllocateWorkingMemory=>SJ_Memory
				OE%set_Direction=>set_SJ_direction
			case ('RNG')
				OE%RescalDirection=>RandomDirection
				OE%AllocateWorkingMemory=>SJ_Memory
				OE%set_Direction=>set_SJ_direction
			case ('NL')
				OE%RescalDirection=>DefaultRescalDirection
				OE%LLSMV=>NLMV
				OE%AllocateWorkingMemory=>NL_Memory
				OE%set_Direction=>set_NL_direction
			case ('GM2')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>GMMemory
				OE%set_Direction=>set_GMDirection
				call printFunc(' DO NOT rescal the searching direction')
			case ('adam2')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>AdamMemory
				OE%set_Direction=>set_AdamDirection
				call printFunc(' DO NOT rescal the searching direction')
			case ('CG2')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>CGMemory
				OE%set_Direction=>set_CGDirection
				call printFunc(' DO NOT rescal the searching direction')
			case ('LBFGS2')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>LBFGSMemory
				OE%set_Direction=>set_LBFGSDirection
				call printFunc(' DO NOT rescal the searching direction')
			case ('SJ2')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>SJ_Memory
				OE%set_Direction=>set_SJ_direction
				call printFunc(' DO NOT rescal the searching direction')
			case ('NG2')
				OE%RescalDirection=>DefaultNotRescalDirection
				OE%AllocateWorkingMemory=>SJ_Memory
				OE%set_Direction=>set_SJ_direction
				call printFunc(' DO NOT rescal the searching direction')
			case default
				call printFunc('NO such case of method, the default method are:')
				call printFunc(' GM     :  Gradient Method')
				call printFunc(' RGM    : Random Gradient Method')
				call printFunc(' CG     : Conjugate Gradient method ')
				call printFunc(' LBFGS  : LBFGS method ')
				call printFunc(' NG     : natural gradient method ')
				call printFunc(' ****2  : the same method, but Do not rescall the direction')
				call printFunc(' R****  : the same method, but randomly the direction')
				call printFunc(' If natural gradient method is in used, the subroutine getting the metric should be add')
				call ErrorFunc
		end select
		OE%method=method
		return
	end subroutine

	subroutine set_method2(OE,set_Direction)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(set_DefaultDirection_interface)::set_Direction
		OE%method='BFGS'
		call printFunc('Set the optimal method in OptimEngine_structure as user define method')
		call printFunc('BFGS')
		call printFunc('BFGS')
		call printFunc('BFGS')
		call printFunc('BFGS')
		OE%RescalDirection=>DefaultRescalDirection
		OE%AllocateWorkingMemory=>GMMemory
		OE%set_Direction=>set_Direction
		return
	end subroutine


	subroutine Set_inStep1Func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(Step1Subroutine_interface)::Func
		call printFunc('Set the inStep1Func in OptimEngine_structure')
		call printFunc(' the inStep1Func is called at the begining of the loop and before calling the targetfunction')
		OE%inStep1=>Func
		return
	end subroutine
	subroutine Set_inStep2Func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(Step2Subroutine_interface)::Func
		call printFunc('Set the inStep2Func in OptimEngine_structure')
		call printFunc(' the inStep2Func is after calling the targetfunction and before setting the searching direction')
		OE%inStep2=>Func
		return
	end subroutine
	subroutine Set_inStep3Func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(Step3Subroutine_interface)::Func
		call printFunc('Set the inStep3Func in OptimEngine_structure')
		call printFunc(' the inStep3Func is after calling the searching_direction and before going on a step')
		OE%inStep3=>Func
		return
	end subroutine
	subroutine Set_inStep4Func(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(Step4Subroutine_interface)::Func
		call printFunc('Set the inStep4Func in OptimEngine_structure')
		call printFunc(' the inStep4Func is after going on a step and at the end of the loop ')
		OE%inStep4=>Func
		return
	end subroutine
	subroutine unSet_inStep1Func(OE)
		class(OptimEngine_structure), intent(inout) :: OE
		call printFunc('unSet the inStep1Func in OptimEngine_structure')
		OE%inStep1=>null()
		return
	end subroutine
	subroutine unSet_inStep2Func(OE)
		class(OptimEngine_structure), intent(inout) :: OE
		call printFunc('inSet the inStep2Func in OptimEngine_structure')
		OE%inStep2=>null()
		return
	end subroutine
	subroutine unSet_inStep3Func(OE)
		class(OptimEngine_structure), intent(inout) :: OE
		call printFunc('Set the inStep3Func in OptimEngine_structure')
		OE%inStep3=>null()
		return
	end subroutine
	subroutine unSet_inStep4Func(OE)
		class(OptimEngine_structure), intent(inout) :: OE
		call printFunc('Set the inStep3Func in OptimEngine_structure')
		OE%inStep3=>null()
		return
	end subroutine
	subroutine Set_BeforeStepFunc(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(BeforeStepSubroutine_interface)::Func
		call printFunc('Set the BeforeStepFunc in OptimEngine_structure')
		OE%BeforeStep=>Func
		return
	end subroutine
	subroutine unSet_BeforeStepFunc(OE)
		class(OptimEngine_structure), intent(inout) :: OE
		call printFunc('unSet the BeforeStepFunc in OptimEngine_structure')
		OE%BeforeStep=>null()
		return
	end subroutine
	subroutine Set_EndStepFunc(OE,Func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(EndStepSubroutine_interface)::Func
		call printFunc('Set the EndStepFun in OptimEngine_structure')
		OE%EndStep=>Func
		return
	end subroutine
	subroutine unSet_EndStepFunc(OE)
		class(OptimEngine_structure), intent(inout) :: OE
		call printFunc('unSet the EndStepFun in OptimEngine_structure')
		OE%EndStep=>null()
		return
	end subroutine

	!nl_YV_func:
	!   d Y
 	!  -----  *  Pi
 	!   d xi
	!
	!
	!nl_YTV_func:
	!              T
	!       /      \
	!      |  d Yj  |
 	!  Qj* | ----- |
 	!       \ d x  /
	!
	!nl_YPG_func:
	!                d E
	!    Yj  and    -----
	!                d Yj

	subroutine set_NL_function(OE,nl_YTV_func,nl_YPG_func,episilo,diag_value)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(nl_YTV_interface)::nl_YTV_func
		procedure(nl_Y_point_gradient_interface)::nl_YPG_func
		real*8,intent(in)::episilo
		real*8,intent(in)::diag_value
		call printFunc('Set the episilo in NL method as episilo=')
		call printFunc(episilo)
		call printFunc('Set the diag_value in NL method as episilo=')
		call printFunc(diag_value)
		OE%nl_episilo = episilo
		OE%nl_diag_value = diag_value
		call printFunc('Set the matrix times vector function used in NL in OptimEngine_structure')
		call printFunc(' nl_YTV_func is to output the p, p = p * dY/dx')
		call printFunc(' nl_YPG_func is to output the ypoint,ygradient')
		OE%NLYTV=>nl_YTV_func
		OE%NLYPG=>nl_YPG_func
		return
	end subroutine

	subroutine set_NL_parameter(OE,episilo,diag_value)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::episilo
		real*8,intent(in)::diag_value
		call printFunc('Set the episilo in NL method as episilo=')
		call printFunc(episilo)
		call printFunc('Set the diag_value in NL method as episilo=')
		call printFunc(diag_value)
		OE%nl_episilo = episilo
		OE%nl_diag_value = diag_value
		return
	end subroutine

	subroutine set_NL_Idetity_function(OE,nl_YTYV_func)
		class(OptimEngine_structure), intent(inout) :: OE
		procedure(nl_YTYV_interface)::nl_YTYV_func
		call printFunc('Set the matrix times vector function used in NL in OptimEngine_structure')
		call printFunc(' nl_YTV_func is to output the p, p=p * dY/dx')
		call printFunc(' nl_YPG_func is to output the ypoint,ygradient')
		OE%NLYTYV=>nl_YTYV_func
		return
	end subroutine

	!***************************************************
	!           use in LBFGS
	!***************************************************

	subroutine OEpoint_St_Yt(OE,st,yt,ith)
		class(OptimEngine_structure),target,intent(in)::OE
		class(OptimElement_structure),pointer,intent(inout)::st,yt
		integer,intent(in)::ith
		if(ith.gt.OE%length)then
			call printFunc('ERROR in point to LBFGSTool')
			!call writemess('ith='+ith)
			!call writemess('OE%length='+OE%length)
			call ErrorFunc
		end if
		st=>OE%st(ith)
		yt=>OE%yt(ith)
		return
	end subroutine

	subroutine pointNewElement(OE,st,yt)
		class(OptimEngine_structure),intent(inout)::OE
		class(OptimElement_structure),pointer,intent(inout)::st,yt
		if(.not.allocated(OE%st))then
			call printFunc('DO NOT allocate st yet')
			call ErrorFunc
		end if
		OE%endindex=OE%endindex+1
		if(OE%endindex.gt.OE%length)then
			OE%endindex=1
			OE%FullFlag=.true.
		end if
		call OEpoint_St_Yt(OE,st,yt,OE%endindex)
		return
	end subroutine

	subroutine element_i(OE,st,yt,ith)
		class(OptimEngine_structure),intent(in)::OE
		class(OptimElement_structure),pointer,intent(inout)::st,yt
		integer,intent(in)::ith
		integer::i
		if(ith.gt.OE%length)then
			call printFunc('ERROR in element_i,1')
			call ErrorFunc
		end if
		if(OE%FullFlag)then
			i=OE%endindex+ith
			if(i.gt.OE%length)then
				i=i-OE%length
			end if
		else
			i=ith
		end if
		call OEpoint_St_Yt(OE,st,yt,i)
		return
	end subroutine


	subroutine allocate_LBFGS_memory(OE,length)
		class(OptimEngine_structure),intent(inout)::OE
		integer,intent(in)::length
		if(length.le.0)then
			call printFunc('ERROR in allocatememory for allocate_LBFGS_memory')
			call ErrorFunc
		end if
		call allocateOptimElementArray(OE%st,length)
		call allocateOptimElementArray(OE%yt,length)
		OE%length=length
		OE%FullFlag=.false.
		OE%endindex=0
		return
	end subroutine

	subroutine deallocate_LBFGS_memory(OE)
		class(OptimEngine_structure),intent(inout)::OE
		integer::i
		if(OE%length.eq.0)return
		deallocate(OE%yt)
		deallocate(OE%st)
		OE%length=0
		OE%FullFlag=.false.
		OE%endindex=0
		return
	end subroutine

	subroutine reset_LBFGS_Endpoint(OE)
		class(OptimEngine_structure),intent(inout)::OE
		OE%endindex=0
		OE%FullFlag=.false.
		return
	end subroutine

	function LBFGS_data_length(OE)
		integer::LBFGS_data_length
		class(OptimEngine_structure),intent(in)::OE
		if(OE%FullFlag)then
			LBFGS_data_length=OE%length
		else
			LBFGS_data_length=OE%endindex
		end if
		return
	end function

	function LBFGS_data_Size(OE)
		integer::LBFGS_data_Size
		class(OptimEngine_structure),intent(in)::OE
		LBFGS_data_Size=OE%length
		return
	end function

	subroutine set_LBFGS_length(OE,length)
		class(OptimEngine_structure),intent(inout)::OE
		integer,intent(in)::length
		call printFunc('Set the length in LBFGS method as length=')
		call printFunc(length)
		OE%length = length
		return
	end subroutine


	!***************************************************
	!           use in NLBFGS
	!***************************************************
	
	subroutine allocate_NL_memory(OE,length)
		class(OptimEngine_structure),intent(inout)::OE
		integer,intent(in)::length
		if(length.le.0)then
			call printFunc('ERROR in allocatememory for allocate_LBFGS_memory')
			call ErrorFunc
		end if
		call allocateOptimElementArray(OE%nl_YTy,length)
		call allocateOptimElementArray(OE%nl_YTgama,length)
		call allocateOptimElementArray(OE%nl_yg,length+1)
		call allocateOptimElementArray(OE%nl_yx,length+1)
		call allocateOptimElementArray(OE%nl_yo,length+1)
		allocate(OE%nl_sy(length,length))
		allocate(OE%nl_sgama(length,length))
		allocate(OE%nl_yxg(length+1,length+1))
		allocate(OE%nl_yxx(length+1,length+1))
		OE%nl_length=length
		return
	end subroutine

	subroutine deallocate_NL_memory(OE)
		class(OptimEngine_structure),intent(inout)::OE
		integer::i
		if(OE%nl_length.eq.0)return
		deallocate(OE%nl_YTy)
		deallocate(OE%nl_YTgama)
		deallocate(OE%nl_yg)
		deallocate(OE%nl_yx)
		deallocate(OE%nl_yo)
		deallocate(OE%nl_yxg)
		deallocate(OE%nl_yxx)
		deallocate(OE%nl_sy)
		deallocate(OE%nl_sgama)

		OE%nl_length=0
		OE%flag1=0
		OE%flag2=0
		OE%flag3=0
		OE%saveflag3=0
		return
	end subroutine

	function NL_data_Size(OE)
		integer::NL_data_Size
		class(OptimEngine_structure),intent(in)::OE
		NL_data_Size=OE%nl_length
		return
	end function

	subroutine set_NL_length(OE,length)
		class(OptimEngine_structure),intent(inout)::OE
		integer,intent(in)::length
		call printFunc('Set the length in NL method as length=')
		call printFunc(length)
		OE%nl_length = length
		return
	end subroutine

	subroutine set_NL_episilo(OE,episilo)
		class(OptimEngine_structure),intent(inout)::OE
		real*8,intent(in)::episilo
		call printFunc('Set the episilo in NL method as episilo=')
		call printFunc(episilo)
		OE%nl_episilo = episilo
		return
	end subroutine

	subroutine set_NL_diag_value(OE,diag_value)
		class(OptimEngine_structure),intent(inout)::OE
		real*8,intent(in)::diag_value
		call printFunc('Set the length in NL method as length=')
		call printFunc(diag_value)
		OE%nl_diag_value = diag_value
		return
	end subroutine

	subroutine NLMV(OE,inoutp,inp)
		! Set the matrix times vector function used in SJ in OptimEngine_structure
		! The metric in the problem is M=Y^T * G * Y
		! Function is to output the p, p=Y^T * G * Y * x
		class(OptimEngine_structure),target, intent(inout) :: OE
		class(OptimElement_structure),target,intent(inout)::inoutp
		class(OptimElement_structure),target,intent(in)::inp
		integer::index
		integer::i,j
		real*8,allocatable::yp(:),gamap(:)
		real*8::norm2
		class(OptimElement_structure),pointer::v1,v2,v3,v4,v5
		class(OptimElement_structure),pointer::nl_point,nl_gradient

		call OE%pointMemory(v1,6)
		call OE%pointMemory(v2,7)
		call OE%pointMemory(v3,8)
		call OE%pointMemory(v4,9)
		call OE%pointMemory(v5,10)
		call OE%pointMemory(nl_point,11)
		call OE%pointMemory(nl_gradient,12)

		call assignmentFunc(v1,inp)
		call assignmentFunc(v2,inp)
		call assignmentFunc(v3,inp)
		call assignmentFunc(v4,inp)
		call assignmentFunc(v5,inp)

		call OE%NLYTYV(inoutp,inp)

		allocate(yp(OE%nl_length))
		allocate(gamap(OE%nl_length))

		if((OE%flag1.gt.0).or.(OE%flag2.eq.0))then
			if(OE%flag1.gt.0)then
				call OE%Y_LparY(OE%nl_yxg(1,1),nl_point,nl_gradient,OE%nl_yx(1),OE%nl_yg(1),OE%nl_yo(1))
				call OE%YDotY(OE%nl_yxx(1,1),nl_point,nl_gradient,OE%nl_yx(1),OE%nl_yg(1),OE%nl_yo(1))
				OE%flag1 = -1	
			end if
		else if((OE%flag2.gt.0).or.(OE%flag3.eq.0))then
			if(OE%flag2.gt.0)then
				call OE%Y_LparY(OE%nl_yxg(2,2),nl_point,nl_gradient,OE%nl_yx(2),OE%nl_yg(2),OE%nl_yo(2))
				call OE%Y_LparY(OE%nl_yxg(2,1),nl_point,nl_gradient,OE%nl_yx(1),OE%nl_yg(1),OE%nl_yo(1))
				call OE%LparY_Y(OE%nl_yxg(1,2),nl_point,nl_gradient,OE%nl_yx(1),OE%nl_yg(1),OE%nl_yo(1))
				OE%nl_sy(1,1) = OE%nl_yxg(2,2) - OE%nl_yxg(2,1) - OE%nl_yxg(1,2) + OE%nl_yxg(1,1)
				call OE%YDotY(OE%nl_yxx(2,2),nl_point,nl_gradient,OE%nl_yx(2),OE%nl_yg(2),OE%nl_yo(2))
				call OE%YDotY(OE%nl_yxx(2,1),nl_point,nl_gradient,OE%nl_yx(1),OE%nl_yg(1),OE%nl_yo(1))
				OE%nl_yxx(1,2) = OE%nl_yxx(2,1)
				OE%nl_sgama(1,1) = OE%nl_yxx(2,2) - OE%nl_yxx(2,1) - OE%nl_yxx(1,2) + OE%nl_yxx(1,1)
				OE%nl_sy(1,1) = OE%nl_sy(1,1) - (OE%nl_episilo*OE%nl_sgama(1,1))
				call OE%YparX_LparY(v2,nl_point,nl_gradient,OE%nl_yx(2),OE%nl_yg(2),OE%nl_yo(2))
				call MinusFunc(v4,v2,nl_gradient)
				call norm2Func(v4,norm2)
				call OE%YparX_LparY(v3,nl_point,nl_gradient,OE%nl_yx(1),OE%nl_yg(1),OE%nl_yo(1))
				call MinusFunc(OE%nl_YTy(1),v2,v3)
				call OE%YparX_Y(v2,nl_point,nl_gradient,OE%nl_yx(2),OE%nl_yg(2),OE%nl_yo(2))
				call OE%YparX_Y(v3,nl_point,nl_gradient,OE%nl_yx(1),OE%nl_yg(1),OE%nl_yo(1))
				call MinusFunc(OE%nl_YTgama(1),v2,v3)
				call A_minus_t_time_B(OE%nl_YTy(1),OE%nl_YTy(1),OE%nl_episilo,OE%nl_YTgama(1))
			end if
			call DotFunc(yp(1),OE%nl_YTy(1),inp)
			call DotFunc(gamap(1),OE%nl_YTgama(1),inp)
			call A_minus_t_time_B(inoutp,inoutp,-yp(1)/OE%nl_sy(1,1),OE%nl_YTy(1))
			call A_minus_t_time_B(inoutp,inoutp,gamap(1)/OE%nl_sgama(1,1),OE%nl_YTgama(1))
			call A_minus_t_time_B(inoutp,inoutp,-OE%nl_diag_value,inp)
			OE%flag2 = -1
		else 
			index = OE%flag3+1
			if(OE%flag3.ne.OE%saveflag3)then
				OE%saveflag3=OE%flag3
				call OE%LparY_Y(OE%nl_yxg(index,index+1),nl_point,nl_gradient,OE%nl_yx(index),OE%nl_yg(index),OE%nl_yo(index))
				call OE%YDotY(OE%nl_yxx(index,index+1),nl_point,nl_gradient,OE%nl_yx(index),OE%nl_yg(index),OE%nl_yo(index))
				do i=1,index+1
					call OE%Y_LparY(OE%nl_yxg(index+1,i),nl_point,nl_gradient,OE%nl_yx(i),OE%nl_yg(i),OE%nl_yo(i))
					call OE%YDotY(OE%nl_yxx(index+1,i),nl_point,nl_gradient,OE%nl_yx(i),OE%nl_yg(i),OE%nl_yo(i))
				end do
				do i=1,index	
					OE%nl_sy(index,i) = OE%nl_yxg(index+1,i+1) - OE%nl_yxg(index,i+1) - OE%nl_yxg(index+1,i) + OE%nl_yxg(index,i)
					OE%nl_sgama(index,i) = OE%nl_yxx(index+1,i+1) - OE%nl_yxx(index,i+1) - OE%nl_yxx(index+1,i) + OE%nl_yxx(index,i)
					OE%nl_sy(index,i) = OE%nl_sy(index,i) - (OE%nl_episilo*OE%nl_sgama(index,i))
					do j=1,i-1
						OE%nl_sgama(index,i) = OE%nl_sgama(index,i) + (OE%nl_sy(index,j)*OE%nl_sy(i,j)/OE%nl_sy(j,j))
						OE%nl_sgama(index,i) = OE%nl_sgama(index,i) - (OE%nl_sgama(index,j)*OE%nl_sgama(i,j)/OE%nl_sgama(j,j))
					end do
				end do
				do i=1,index
					if(i.eq.1)then
						call OE%YparX_LparY(v2,nl_point,nl_gradient,OE%nl_yx(1),OE%nl_yg(1),OE%nl_yo(1))
						call OE%YparX_Y(v4,nl_point,nl_gradient,OE%nl_yx(1),OE%nl_yg(1),OE%nl_yo(1))
					end if
					call OE%YparX_LparY(v3,nl_point,nl_gradient,OE%nl_yx(i+1),OE%nl_yg(i+1),OE%nl_yo(i+1))
					call OE%YparX_Y(v5,nl_point,nl_gradient,OE%nl_yx(i+1),OE%nl_yg(i+1),OE%nl_yo(i+1))
					call MinusFunc(OE%nl_YTy(i),v3,v2)
					call MinusFunc(OE%nl_YTgama(i),v5,v4)
					call A_minus_t_time_B(OE%nl_YTy(i),OE%nl_YTy(i),OE%nl_episilo,OE%nl_YTgama(i))
					call assignmentFunc(v2,v3)
					call assignmentFunc(v4,v5)
				end do
				do i=2,index
					do j=1,i-1
						call A_minus_t_time_B(OE%nl_YTgama(i),OE%nl_YTgama(i),-OE%nl_sy(i,j)/OE%nl_sy(j,j),OE%nl_YTy(j))
						call A_minus_t_time_B(OE%nl_YTgama(i),OE%nl_YTgama(i),OE%nl_sgama(i,j)/OE%nl_sgama(j,j),OE%nl_YTgama(j))
					end do
				end do
			end if
			do i=1,index
				call DotFunc(yp(i),OE%nl_YTy(i),inp)
				call DotFunc(gamap(i),OE%nl_YTgama(i),inp)
			end do
			do i=1,index
				call A_minus_t_time_B(inoutp,inoutp,-yp(i)/OE%nl_sy(i,i),OE%nl_YTy(i))
				call A_minus_t_time_B(inoutp,inoutp,gamap(i)/OE%nl_sgama(i,i),OE%nl_YTgama(i))
			end do
			call A_minus_t_time_B(inoutp,inoutp,-OE%nl_diag_value,inp)
		end if	
		if(OE%flag3.eq.(OE%nl_length-1))then
			call OE%deallocate_NL_memory()
			call OE%allocate_NL_memory(LBFGS_default_memory_length)
			return
		end if
		return
	end subroutine

	!***************************************************
	!           default function
	!***************************************************

	subroutine DefaultRescalDirection(OE,Dir)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),intent(inout)::Dir
		real*8::norm2
		call norm2Func(Dir,norm2)
		norm2=dsqrt(norm2)
		if(norm2.le.zero_gradient)then
			call printFunc('*******  DefaultRescalDirection WRONNING in OptimizationTools.f90 ******* ')
			call printFunc('   |Dir| is too small(<1d-16), |Dir|=')
			call printFunc(norm2)
			call printFunc('   DO NOT rescal the direction')
			call printFunc('************************************************************************* ')
		else
			call multiplyFunc(Dir,1d0/norm2)
		end if
		return
	end subroutine

	subroutine RandomDirection(OE,Dir)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),intent(inout)::Dir
		call RandomElement(Dir)
		return
	end subroutine

	subroutine DefaultNotRescalDirection(OE,Dir)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),intent(inout)::Dir
		return
	end subroutine

	subroutine DefaultAllocateWorkingMemory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),pointer,intent(inout)::Gradient,direction
		call printERRORMessage()
	end subroutine
	subroutine set_DefaultDirection(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),intent(in)::point,Gradient
		class(OptimElement_structure),intent(inout)::direction
		call printERRORMessage()
	end subroutine
	subroutine printERRORMessage()
		call printFunc('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		call printFunc('% DO NOT set the optimal method yet     %')
		call printFunc('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		call ErrorFunc
	end subroutine
	subroutine Defaulttarget_Function(A,outVal,outGradient,point)
		class(OptimEngine_structure),target, intent(inout) :: A
		real*8,target,intent(inout)::outVal
		class(OptimElement_structure),target,intent(inout)::outGradient
		class(OptimElement_structure),target,intent(in)::point
		call printFunc('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		call printFunc('% DO NOT set the target_Function yet    %')
		call printFunc('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		call ErrorFunc
	end subroutine

	subroutine DefaultPointToOptimData(A,point)
		class(OptimEngine_structure),target, intent(inout) :: A
		class(OptimElement_structure),pointer,intent(inout)::point
		call printFunc('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		call printFunc('% DO NOT set the PointToOptimData yet    %')
		call printFunc('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		call ErrorFunc
	end subroutine

	


	!***************************************************
	!            optimization
	!***************************************************

	subroutine Optimization1(OE,T0,tau0,numStep,Point,outValue)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::T0,tau0
		integer,intent(in)::numStep
		class(OptimElement_structure),intent(inout)::Point
		real*8,intent(inout)::outValue
		integer::i
		real*8::t
		class(OptimElement_structure),pointer::Gradient,direction

		call OE%AllocateWorkingMemory(Gradient,direction)
		call assignmentFunc(Gradient,point)
		if(associated(OE%BeforeStep))call OE%BeforeStep(outValue,Gradient,point)
		do i=1,numStep
			t=T0*tau0/(T0+dble(i))
			if(associated(OE%inStep1))call OE%inStep1(outValue,direction,Gradient,point,i,t)
			call OE%target_Function(outValue,Gradient,point)
			if(associated(OE%inStep2))call OE%inStep2(outValue,direction,Gradient,point,i,t)
			call OE%set_direction(direction,Gradient,point)
			call OE%RescalDirection(direction)
			if(associated(OE%inStep3))call OE%inStep3(outValue,direction,Gradient,point,i,t)
			call A_minus_t_time_B(point,point,t,direction)
			if(associated(OE%inStep4))call OE%inStep4(outValue,direction,Gradient,point,i,t)
			if(OE%check_stop(outValue,Gradient,t))exit
		end do
		if(associated(OE%EndStep))call OE%EndStep(outValue,Gradient,point)
		return
	end subroutine

	subroutine Optimization2(OE,t,numStep,Point,outValue)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::t
		integer,intent(in)::numStep
		class(OptimElement_structure),intent(inout)::Point
		real*8,intent(inout)::outValue
		integer::i
		class(OptimElement_structure),pointer::Gradient,direction
			
		call OE%AllocateWorkingMemory(Gradient,direction)
		call assignmentFunc(Gradient,point)
		if(associated(OE%BeforeStep))call OE%BeforeStep(outValue,Gradient,point)
		do i=1,numStep
			if(associated(OE%inStep1))call OE%inStep1(outValue,direction,Gradient,point,i,t)
			call OE%target_Function(outValue,Gradient,point)
			if(associated(OE%inStep2))call OE%inStep2(outValue,direction,Gradient,point,i,t)
			call OE%set_direction(direction,Gradient,point)
			call OE%RescalDirection(direction)
			if(associated(OE%inStep3))call OE%inStep3(outValue,direction,Gradient,point,i,t)
			call A_minus_t_time_B(point,point,t,direction)
			if(associated(OE%inStep4))call OE%inStep4(outValue,direction,Gradient,point,i,t)
			if(OE%check_stop(outValue,Gradient,t))exit
		end do
		if(associated(OE%EndStep))call OE%EndStep(outValue,Gradient,point)
		return
	end subroutine

	subroutine Optimization3(OE,NlinearStep,numStep,Point,outValue)
		class(OptimEngine_structure), intent(inout) :: OE
		integer,intent(in)::NlinearStep,numStep
		class(OptimElement_structure),intent(inout)::Point
		real*8,intent(inout)::outValue
		real*8::x
		integer::i
		class(OptimElement_structure),pointer::Gradient,direction

		call OE%AllocateWorkingMemory(Gradient,direction)
		call assignmentFunc(Gradient,point)
		if(associated(OE%BeforeStep))call OE%BeforeStep(outValue,Gradient,point)
		x=OE%first_step_in_Line_search
		do i=1,numStep
			if(associated(OE%inStep1))call OE%inStep1(outValue,direction,Gradient,point,i,x)
			call OE%target_Function(outValue,Gradient,point)
			if(associated(OE%inStep2))call OE%inStep2(outValue,direction,Gradient,point,i,x)
			call OE%set_direction(direction,Gradient,point)
			call OE%RescalDirection(direction)
			if(associated(OE%inStep3))call OE%inStep3(outValue,direction,Gradient,point,i,x)
			call OE%LinearSearch(NlinearStep,point,direction,x,outValue,Gradient)
			if(associated(OE%inStep4))call OE%inStep4(outValue,direction,Gradient,point,i,x)
			if(OE%check_stop(outValue,Gradient,x))exit
		end do
		if(associated(OE%EndStep))call OE%EndStep(outValue,Gradient,point)
		return
	end subroutine

	subroutine Optimization4(OE,T0,tau0,numStep,outValue)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8,intent(in)::T0,tau0
		integer,intent(in)::numStep
		class(OptimElement_structure),pointer::Point
		real*8,intent(inout)::outValue
		integer::i
		real*8::t
		class(OptimElement_structure),pointer::Gradient,direction

		call OE%PointToOptimData(Point)

		call OE%AllocateWorkingMemory(Gradient,direction)
		call assignmentFunc(Gradient,point)
		if(associated(OE%BeforeStep))call OE%BeforeStep(outValue,Gradient,point)
		do i=1,numStep
			t=T0*tau0/(T0+dble(i))
			if(associated(OE%inStep1))call OE%inStep1(outValue,direction,Gradient,point,i,t)
			call OE%target_Function(outValue,Gradient,point)
			if(associated(OE%inStep2))call OE%inStep2(outValue,direction,Gradient,point,i,t)
			call OE%set_direction(direction,Gradient,point)
			call OE%RescalDirection(direction)
			if(associated(OE%inStep3))call OE%inStep3(outValue,direction,Gradient,point,i,t)
			call A_minus_t_time_B(point,point,t,direction)
			if(associated(OE%inStep4))call OE%inStep4(outValue,direction,Gradient,point,i,t)
			if(OE%check_stop(outValue,Gradient,t))exit
		end do
		if(associated(OE%EndStep))call OE%EndStep(outValue,Gradient,point)
		return
	end subroutine

	subroutine Optimization5(OE,t,numStep,outValue)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::t
		integer,intent(in)::numStep
		class(OptimElement_structure),pointer::Point
		real*8,intent(inout)::outValue
		integer::i
		class(OptimElement_structure),pointer::Gradient,direction

		call OE%PointToOptimData(Point)
			
		call OE%AllocateWorkingMemory(Gradient,direction)
		call assignmentFunc(Gradient,point)
		if(associated(OE%BeforeStep))call OE%BeforeStep(outValue,Gradient,point)
		do i=1,numStep
			if(associated(OE%inStep1))call OE%inStep1(outValue,direction,Gradient,point,i,t)
			call OE%target_Function(outValue,Gradient,point)
			if(associated(OE%inStep2))call OE%inStep2(outValue,direction,Gradient,point,i,t)
			call OE%set_direction(direction,Gradient,point)
			call OE%RescalDirection(direction)
			if(associated(OE%inStep3))call OE%inStep3(outValue,direction,Gradient,point,i,t)
			call A_minus_t_time_B(point,point,t,direction)
			if(associated(OE%inStep4))call OE%inStep4(outValue,direction,Gradient,point,i,t)
			if(OE%check_stop(outValue,Gradient,t))exit
		end do
		if(associated(OE%EndStep))call OE%EndStep(outValue,Gradient,point)
		return
	end subroutine

	!Optimization6 no this type subroutine, line search should modify the point
	! and call the energy function. the data of point should be a input parameter
	
	!***************************************************
	!            LBFGS optimization
	!***************************************************

	!***************************************************
	!            LBFGS optimization
	!***************************************************

	! p is \Delta f(x_i)
	! s_i= x_{i+1} - x_i
	! y_i= p_{i+1} - p_i

	subroutine LBFGS_direction(OE,p,inputp)
		class(OptimEngine_structure),intent(inout)::OE
		class(OptimElement_structure),intent(inout)::p
		class(OptimElement_structure),intent(in)::inputp

		real*8::temp,temp2,beta
		real*8,allocatable::alpha(:)
		integer::i,length

		class(OptimElement_structure),pointer::s,y

		call assignmentFunc(p,inputp)

		length = OE%LBFGS_data_length()
		allocate(alpha(length))
		do i=length,1,-1
			call OE%i(s,y,i)
			call DotFunc(temp,y,s)
			call DotFunc(alpha(i),s,p)
			alpha(i)=alpha(i)/temp
			call A_minus_t_time_B(p,p,alpha(i),y)
		end do

		!p = Hk * p

		do i=1,length
			call OE%i(s,y,i)
			call DotFunc(temp,y,s)
			call DotFunc(beta,y,p)
			beta=beta/temp
			call A_minus_t_time_B(p,p,beta-alpha(i),s)
		end do



		return
	end subroutine

	subroutine LBFGSMemory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),pointer,intent(inout)::Gradient,direction
		integer::memorylen
		memorylen=4
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.memorylen)then
				deallocate(OE%workingMemory)
				call allocateOptimElementArray(OE%workingMemory,memorylen)
			end if
		else
			call allocateOptimElementArray(OE%workingMemory,memorylen)
		end if
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		if((OE%LBFGS_data_Size().le.0))then
			call printFunc(' DO not allocate momery for LBFGS method')
			call printFunc(' set the momery to defaut length=')
			call printFunc(LBFGS_default_memory_length)
			call OE%allocate_LBFGS_memory(LBFGS_default_memory_length)
			OE%length=LBFGS_default_memory_length
		else
			call OE%allocate_LBFGS_memory(OE%length)
		end if
		return
	end subroutine

	subroutine set_LBFGSDirection(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),intent(in)::point,Gradient
		class(OptimElement_structure),intent(inout)::direction
		class(OptimElement_structure),pointer::yt,st,SavePoint,SaveGradient
		logical,save::first=.true.
		call OE%pointMemory(SavePoint,3)
		call OE%pointMemory(SaveGradient,4)
	
		if(first)then
			call assignmentFunc(direction,Gradient)
			call assignmentFunc(SaveGradient,Gradient)
			call assignmentFunc(SavePoint,point)
			first=.false.
			return		
		end if
		
		call OE%NewElement(st,yt)
		call MinusFunc(yt,Gradient,SaveGradient)	
		call MinusFunc(st,point,SavePoint)
		!call A_minus_t_time_B(yt,Gradient,1d0,SaveGradient)	
		!call A_minus_t_time_B(st,point,1d0,SavePoint)
		
		call LBFGS_direction(OE,Direction,Gradient)

		call assignmentFunc(SaveGradient,Gradient)
		call assignmentFunc(SavePoint,point)

		return
	end subroutine

	!***************************************************
	!           Adam optimization
	!***************************************************


	subroutine AdamMemory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),pointer,intent(inout)::Gradient,direction
		!type(OptimElement_structure),pointer::s,g
		integer::memorylen
		memorylen=6
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.memorylen)then
				deallocate(OE%workingMemory)
				call allocateOptimElementArray(OE%workingMemory,memorylen)
			end if
		else
			call allocateOptimElementArray(OE%workingMemory,memorylen)
		end if
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		!call OE%pointMemory(s,3)
		!call OE%pointMemory(g,4)
		!call s%empty()
		!call g%empty()

		return
	end subroutine

	subroutine set_AdamDirection(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),intent(in)::point,Gradient
		class(OptimElement_structure),intent(inout)::direction
		class(OptimElement_structure),pointer::s,r,TMP,TMP2
		real*8::a,b,snorm,rnorm
		real*8::rho,rhot
		call OE%pointMemory(s,3)
		call OE%pointMemory(r,4)
		call OE%pointMemory(TMP,5)
		call OE%pointMemory(TMP2,6)
		OE%Adam_resetcounter=OE%Adam_resetcounter+1
		if(OE%Adam_resetcounter.eq.OE%Adam_resetStep)then
			OE%Adam_resetcounter=1
		end if

		if(OE%Adam_resetcounter.gt.1)then
			rho=OE%Adam_rho1
			OE%Adam_rho1t=OE%Adam_rho1t*rho
			rhot=OE%Adam_rho1t
			a=rho 
			b=(1d0-rho) 
			!s=(rho*s)+((1-rho)*Gradient)
			call assignmentFunc(tmp,s)
			call multiplyFunc(tmp,a)
			call A_minus_t_time_B(s,tmp,-b,Gradient) !s=(a*s)+(b*Gradient)
			
		else
			call assignmentFunc(s,Gradient)!s=Gradient
			OE%Adam_rho1t=OE%Adam_rho1
		end if


		if(OE%Adam_resetcounter.gt.1)then
			rho=OE%Adam_rho2
			OE%Adam_rho2t=OE%Adam_rho2t*rho
			rhot=OE%Adam_rho2t
			a=rho 
			b=(1d0-rho) 

			!r=(rho*r)+((1-rho)*Gradient*Gradient)
			call assignmentFunc(tmp2,r)
			call multiplyFunc(tmp2,a)
			call element_product(tmp,Gradient,Gradient)
			call A_minus_t_time_B(r,tmp2,-b,tmp)
			!r=(a*r)+(b*element_product(Gradient,Gradient))
		else
			call element_product(r,Gradient,Gradient)
			!r=element_product(Gradient,Gradient)
			OE%Adam_rho2t=OE%Adam_rho2
		end if

		snorm=1d0/(1d0-OE%Adam_rho1t)
		rnorm=1d0/(1d0-OE%Adam_rho2t)

		call assignmentFunc(tmp,s)
		call multiplyFunc(tmp,snorm)
		call assignmentFunc(tmp2,r)
		call multiplyFunc(tmp2,rnorm)
		call A_over_sqrt_B(direction,tmp,tmp2)
		!call direction%print
		!call s%print
		!call r%print
		!call Gradient%print
		!read(*,*)
		return
	end subroutine

	!***************************************************
	!            Gradient optimization
	!***************************************************

	subroutine GMMemory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),pointer,intent(inout)::Gradient,direction
		integer::memorylen
		memorylen=2
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.memorylen)then
				deallocate(OE%workingMemory)
				call allocateOptimElementArray(OE%workingMemory,memorylen)
			end if
		else
			call allocateOptimElementArray(OE%workingMemory,memorylen)
		end if
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		return
	end subroutine
	subroutine set_GMDirection(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),intent(in)::point,Gradient
		class(OptimElement_structure),intent(inout)::direction
		call assignmentFunc(direction,Gradient)
		return
	end subroutine

	!***************************************************
	!           Conjugate Gradient optimization
	!***************************************************

	subroutine CG_direction(OE,dir,Gra,priorgra)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),intent(in)::Gra
		class(OptimElement_structure),intent(inout)::dir,priorgra
		real*8::direction_corr
		logical,save::First=.true.
		if(First) then
			call assignmentFunc(dir,gra)
			First=.false.
		else
			direction_corr=correctPara(OE,dir,gra,priorGra)
			call A_minus_t_time_B(dir,gra,direction_corr,dir)
		end if
		call assignmentFunc(priorgra,Gra)
		return
	end subroutine
	
	function correctPara(OE,dir,newgra,gra)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::correctPara
		class(OptimElement_structure),intent(in)::newgra,gra,dir
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
				call defalutPrintString('ERROR CG_direction_Flag=')
				call ErrorFunc
		end select
		return
	end function

	function Crowder_Wolfe(OE,dir,newgra,gra)result(Res)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::Res
		class(OptimElement_structure),intent(in)::dir,newgra,gra
		real*8::TMPr
		class(OptimElement_structure),allocatable::TMP
		allocate(TMP,source=dir)
		call MinusFunc(TMP,newgra,gra)
		call DotFunc(Res,newgra,TMP)
		call DotFunc(TMPr,dir,TMP)
		Res=Res/TMPr
		return
	end function

	function Fletcher_Reeves(OE,dir,newgra,gra)result(Res)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::Res
		class(OptimElement_structure),intent(in)::dir,newgra,gra
		real*8::TMPr
		call norm2Func(newgra,Res)
		call norm2Func(gra,TMPr)
		Res=Res/TMPr
		return
	end function
	function Dixon(OE,dir,newgra,gra)result(Res)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::Res
		class(OptimElement_structure),intent(in)::dir,newgra,gra
		real*8::TMPr
		call norm2Func(newgra,Res)
		call DotFunc(TMPr,dir,gra)
		Res=-1d0*Res/TMPr
		return
	end function
	function Polak_Ribiere_Polyak(OE,dir,newgra,gra)result(Res)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::Res
		class(OptimElement_structure),intent(in)::dir,newgra,gra
		real*8::TMPr
		class(OptimElement_structure),allocatable::TMP
		allocate(TMP,source=dir)
		call MinusFunc(TMP,newgra,gra)
		call DotFunc(Res,newgra,TMP)
		call norm2Func(gra,TMPr)
		!Res=newgra.dot.(newgra-gra)
		!TMPr=gra%dnorm2()
		Res=Res/TMPr
		return
	end function
	function Dai_Yuan(OE,dir,newgra,gra)result(Res)
		class(OptimEngine_structure), intent(inout) :: OE
		real*8::Res
		class(OptimElement_structure),intent(in)::dir,newgra,gra
		real*8::TMPr
		class(OptimElement_structure),allocatable::TMP
		allocate(TMP,source=dir)
		call norm2Func(newgra,Res)
		call MinusFunc(TMP,newgra,gra)
		call DotFunc(TMPr,dir,TMP)
		!Res=newgra%dnorm2()
		!TMPr=dir.dot.(newgra-gra)
		Res=Res/TMPr
		return
	end function

	subroutine CGMemory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),pointer,intent(inout)::Gradient,direction
		integer::memorylen
		memorylen=4
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.memorylen)then
				deallocate(OE%workingMemory)
				call allocateOptimElementArray(OE%workingMemory,memorylen)
			end if
		else
			call allocateOptimElementArray(OE%workingMemory,memorylen)
		end if
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		return
	end subroutine
	subroutine set_CGDirection(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),intent(in)::point,Gradient
		class(OptimElement_structure),intent(inout)::direction
		class(OptimElement_structure),pointer::gra0,dir0
		call OE%pointMemory(gra0,3)
		call OE%pointMemory(dir0,4)
		call CG_direction(OE,dir0,Gradient,gra0)
		call assignmentFunc(direction,dir0)
		return
	end subroutine

	!***************************************************
	!           NoName optimization
	!***************************************************


	subroutine SJ_Memory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),pointer,intent(inout)::Gradient,direction
		integer::memorylen
		class(OptimElement_structure),pointer::Y,G
		memorylen=8
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.memorylen)then
				deallocate(OE%workingMemory)
				call allocateOptimElementArray(OE%workingMemory,memorylen)
			end if
		else
			call allocateOptimElementArray(OE%workingMemory,memorylen)
		end if
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		return
	end subroutine


	subroutine set_SJ_direction(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),intent(in)::point,Gradient
		class(OptimElement_structure),intent(inout)::direction
		class(OptimElement_structure),pointer::gra0,dir0,Y,G,Ap,r,p
		call OE%pointMemory(Ap,3)
		call OE%pointMemory(r,4)
		call OE%pointMemory(p,5)
		call assignmentFunc(Ap,Gradient)
		call assignmentFunc(r,Gradient)
		call assignmentFunc(p,Gradient)
		if(associated(OE%LLSMV))then
			call OE%runCGLLS(direction,Gradient,Ap,r,p)
		else if(associated(OE%externalCGLLS))then
			call OE%externalCGLLS(direction,Gradient)
		else
			call printFunc('DO not set the SJ yet')
			call errorFunc
		end if
		
		return
	end subroutine

	!***************************************************
	!           NL optimization
	!***************************************************


	subroutine NL_Memory(OE,Gradient,direction)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),pointer,intent(inout)::Gradient,direction
		integer::memorylen
		class(OptimElement_structure),pointer::Y,G
		memorylen=15
		if(allocated(OE%workingMemory))then
			if(size(OE%workingMemory).ne.memorylen)then
				deallocate(OE%workingMemory)
				call allocateOptimElementArray(OE%workingMemory,memorylen)
			end if
		else
			call allocateOptimElementArray(OE%workingMemory,memorylen)
		end if
		call OE%pointMemory(Gradient,1)
		call OE%pointMemory(direction,2)
		if((OE%NL_data_Size().le.0))then
			call printFunc(' DO not allocate momery for NL method')
			call printFunc(' set the momery to defaut length=')
			call printFunc(LBFGS_default_memory_length)
			call OE%allocate_NL_memory(LBFGS_default_memory_length)
			OE%nl_length=LBFGS_default_memory_length
		else
			call OE%allocate_NL_memory(OE%nl_length)
		end if
		return
	end subroutine

	subroutine set_NL_direction(OE,direction,Gradient,point)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure),intent(in)::point,Gradient
		class(OptimElement_structure),intent(inout)::direction
		class(OptimElement_structure),pointer::Ap,r,p
		class(OptimElement_structure),pointer::nl_point,nl_gradient

		call OE%pointMemory(Ap,3)
		call OE%pointMemory(r,4)
		call OE%pointMemory(p,5)
		call OE%pointMemory(nl_point,11)
		call OE%pointMemory(nl_gradient,12)

		call assignmentFunc(Ap,Gradient)		
		call assignmentFunc(r,Gradient)	
		call assignmentFunc(p,Gradient)
		call assignmentFunc(nl_point,point)	
		call assignmentFunc(nl_gradient,Gradient)
	

		if(associated(OE%YDotY))then
			if(OE%flag1.eq.0)then
				call assignmentFunc(OE%nl_yx(1),point)	
				call assignmentFunc(OE%nl_yg(1),Gradient)
				call OE%store_other(OE%nl_yo(1))
				!call OE%NLYPG(OE%nl_yx(1),OE%nl_yg(1),point)
				OE%flag1=OE%flag1+1
			else if(OE%flag2.eq.0)then
				call assignmentFunc(OE%nl_yx(2),point)	
				call assignmentFunc(OE%nl_yg(2),Gradient)
				call OE%store_other(OE%nl_yo(2))
				!call OE%NLYPG(OE%nl_yx(2),OE%nl_yg(2),point)
				OE%flag2=OE%flag2+1
			else
				OE%flag3=OE%flag3+1
				!call OE%NLYPG(OE%nl_yx(OE%flag3+2),OE%nl_yg(OE%flag3+2),point)
				call assignmentFunc(OE%nl_yx(OE%flag3+2),point)	
				call assignmentFunc(OE%nl_yg(OE%flag3+2),Gradient)
				call OE%store_other(OE%nl_yo(OE%flag3+2))
			end if
		else
			call printFunc('DO not set the NL yet')
			call errorFunc
		end if

		if(associated(OE%YparX_Y))then
			call OE%runCGLLS(direction,Gradient,Ap,r,p)
		else if(associated(OE%externalCGLLS))then
			call OE%externalCGLLS(direction,Gradient)
		else
			call printFunc('DO not set the NL yet')
			call errorFunc
		end if
		return
	end subroutine

	!use for the solve LLS with CG method
	!
	! AX=b, A is a n times n positive definite matrix
	!
	!


	subroutine runCGLLS(OE,inoutx,b,memoryAp,memory_r,memory_p)
		class(OptimEngine_structure), intent(inout) :: OE
		class(OptimElement_structure)::inoutx,b,memoryAp,memory_r,memory_p
		integer::i
		real*8::rr,alpha,beta
		call OE%InititalCGLLSStartPoint(inoutx,b)
		call OE%LLSMV(memory_p,inoutx)
		call MinusFunc(memory_r,b,memory_p)
		call assignmentFunc(memory_p,memory_r)
		do i=1,OE%CGLLSmaxRunning
			call norm2Func(memory_r,rr)
			if(rr.le.OE%CGLLSerror)exit
			call OE%LLSMV(memoryAp,memory_p)
			call DotFunc(alpha,memory_p,memoryAp)
			alpha=rr/alpha
			call A_minus_t_time_B(inoutx,inoutx,-alpha,memory_p)
			call A_minus_t_time_B(memory_r,memory_r,alpha,memoryAp)
			call norm2Func(memory_r,beta)
			beta=beta/rr
			call A_minus_t_time_B(memory_p,memory_r,-beta,memory_p)
		end do

		if(OE%CGLLSerror.gt.0d0)then
			if(rr.gt.OE%CGLLSerror)then
				call printFunc('*******  CGLLS WRONNING in OptimizationTools.f90 ******* ')
				call printFunc('    The iteration may fail in convergence, the error is: ')
				call printFunc(rr)
				call printFunc('    The norm2 of gra is: ')
				call norm2Func(b,beta)
				call printFunc(beta)
				call printFunc('********************************************************* ')
			end if
		end if
		return
	end subroutine


	subroutine defaultInititalCGLLSStartPointFromZero(OE,inoutx,b)
		class(OptimEngine_structure),target, intent(inout) :: OE
		class(OptimElement_structure),target, intent(inout) :: inoutx
		class(OptimElement_structure),target, intent(in)::b
		call assignmentFunc(inoutx,b)
		call zeroFunc(inoutx)
		return
	end subroutine

	subroutine defaultInititalCGLLSStartPointFromGra(OE,inoutx,b)
		class(OptimEngine_structure),target, intent(inout) :: OE
		class(OptimElement_structure),target, intent(inout) :: inoutx
		class(OptimElement_structure),target, intent(in)::b
		call assignmentFunc(inoutx,b)
		return
	end subroutine
	subroutine defaultInititalCGLLSStartPointFromZero2(OE,inoutx,b)
		class(OptimEngine_structure),target, intent(inout) :: OE
		class(OptimElement_structure),target, intent(inout) :: inoutx
		class(OptimElement_structure),target, intent(in)::b
		logical,save::First=.true.
		if(First)then
			call assignmentFunc(inoutx,b)
			call zeroFunc(inoutx)
			First=.false.
		end if
		return
	end subroutine

	subroutine defaultInititalCGLLSStartPointFromGra2(OE,inoutx,b)
		class(OptimEngine_structure),target, intent(inout) :: OE
		class(OptimElement_structure),target, intent(inout) :: inoutx
		class(OptimElement_structure),target, intent(in)::b
		logical,save::First=.true.
		if(First)then
			call assignmentFunc(inoutx,b)
			First=.false.
		end if
		return
	end subroutine

	!*************************************************************
	!*************************************************************
	!  linear search4 use the information of Energy only
	!      call OE%Energy_Function do not call OE%target_Function
	!      DO not use and modify inoutgra
	!*************************************************************

	subroutine LinearSearch4(OE,max_running,point,diretion,x,outValue,inoutgra)
		class(OptimEngine_structure),target, intent(inout) :: OE
		real*8,intent(inout),target::x,outValue
		class(OptimElement_structure),target,intent(inout)::point
		class(OptimElement_structure),target,intent(inout)::inoutgra
		class(OptimElement_structure),target,intent(in)::diretion
		integer,target,intent(in)::max_running
		real*8::allx(3),allf(3),newx,newf,g0
		class(OptimElement_structure),allocatable::NewPoint
		class(OptimElement_structure),pointer::dir
		integer,save::restart_num=0
		integer::ith,i
		if(max_running.lt.2)then
			call printFunc('ERROR in LineSearch,input max_running should >=2')
			call ErrorFunc
		end if
		dir=>diretion
		allocate(NewPoint,mold=point)
		call assignmentFunc(NewPoint,point)
		allx(1)=0d0
		allf(1)=outValue
		call DotFunc(g0,inoutgra,dir)
		if(g0.lt.0)then
			call printFunc('   ')
			call printFunc('               WORNING in LineSearch')
			call printFunc('diretion*gradient<0 in LineSearch, reset the direction as the gradient')
			call printFunc('   ')
			call OE%RescalDirection(inoutgra)
			dir=>inoutgra
			call DotFunc(g0,inoutgra,dir)
		end if
		g0=-g0


		allx(2)=x
		call A_minus_t_time_B(NewPoint,point,allx(2),dir)
		call OE%Energy_Function(allf(2),NewPoint)

		
		allx(3)=find_the_second_point(allx(2),allf(2),outValue,g0)

		call A_minus_t_time_B(NewPoint,point,allx(3),dir)
		call OE%Energy_Function(allf(3),NewPoint)

		if(max_running.eq.2)then
			if(allf(2).lt.allf(3))then
				x=allx(2)
				outValue=allf(2)
			else
				x=allx(3)
				outValue=allf(3)
			end if
			call A_minus_t_time_B(point,point,x,dir)
			return
		end if

		do i=3,max_running
			call reorderdata(allx,allf)
			newx=find_the_next_point(allx,allf)
			call A_minus_t_time_B(NewPoint,point,newx,dir)
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
		x=x*OE%line_search_final_step_ratio
		
		if((OE%search_restart_step.gt.0).and.(x.le.OE%line_search_min_step)&
			   .and.(restart_num.lt.OE%line_search_restart_max_num))then
			call OE%restart_line_search_subroutone(max_running,point,diretion,x,outValue,inoutgra)
			restart_num=restart_num+1
		else
			call A_minus_t_time_B(point,point,x,dir)
		end if
		
		return
	end subroutine

	subroutine defult_restart_line_search_subroutone(OE,max_running,point,diretion,x,outValue,inoutgra)
		class(OptimEngine_structure),target, intent(inout) :: OE
		real*8,intent(inout),target::x,outValue
		class(OptimElement_structure),target,intent(inout)::point
		class(OptimElement_structure),target,intent(inout)::inoutgra
		class(OptimElement_structure),target,intent(in)::diretion
		integer,target,intent(in)::max_running
		call OE%RescalDirection(inoutgra)
		x=OE%search_restart_step
		call A_minus_t_time_B(point,point,x,inoutgra)
		call OE%Energy_Function(outValue,point)
		return
	end subroutine

end module
