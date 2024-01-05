#define subSymTensorFuncName subSymTensori
#define subTensorFuncName subTensori
#define DATATYPE integer
#include "templet/subTensor0.f90"
#undef subSymTensorFuncName
#undef subTensorFuncName
#undef DATATYPE

#define subSymTensorFuncName subSymTensors
#define subTensorFuncName subTensors
#define DATATYPE real*4
#include "templet/subTensor0.f90"
#undef subSymTensorFuncName
#undef subTensorFuncName
#undef DATATYPE

#define subSymTensorFuncName subSymTensord
#define subTensorFuncName subTensord
#define DATATYPE real*8
#include "templet/subTensor0.f90"
#undef subSymTensorFuncName
#undef subTensorFuncName
#undef DATATYPE

#define subSymTensorFuncName subSymTensorc
#define subTensorFuncName subTensorc
#define DATATYPE complex*8
#include "templet/subTensor0.f90"
#undef subSymTensorFuncName
#undef subTensorFuncName
#undef DATATYPE

#define subSymTensorFuncName subSymTensorz
#define subTensorFuncName subTensorz
#define DATATYPE complex*16
#include "templet/subTensor0.f90"
#undef subSymTensorFuncName
#undef subTensorFuncName
#undef DATATYPE

#define subSymTensorFuncName subSymTensorl
#define subTensorFuncName subTensorl
#define DATATYPE logical
#include "templet/subTensor0.f90"
#undef subSymTensorFuncName
#undef subTensorFuncName
#undef DATATYPE

#define subSymTensorFuncName subSymTensora
#define subTensorFuncName subTensora
#define DATATYPE character(len=characterlen)
#include "templet/subTensor0.f90"
#undef subSymTensorFuncName
#undef subTensorFuncName
#undef DATATYPE

	function SubTensor1Func(T,Legi,QNi,degi,keepQN_)result(Res)
		type(Tensor)::Res
		class(Tensor),intent(in) ::T
		integer,intent(in)::Legi,QNi,degi
		logical,optional,intent(in)::keepQN_
		if(T%getSymmetryFlag())then
			select case(T%getType())
				case(1)
					call subSymTensori(Res,T,Legi,QNi,degi,keepQN_)
				case(2)
					call subSymTensors(Res,T,Legi,QNi,degi,keepQN_)
				case(3)
					call subSymTensord(Res,T,Legi,QNi,degi,keepQN_)
				case(4)
					call subSymTensorc(Res,T,Legi,QNi,degi,keepQN_)
				case(5)
					call subSymTensorz(Res,T,Legi,QNi,degi,keepQN_)
				case(6)
					call subSymTensorl(Res,T,Legi,QNi,degi,keepQN_)
				case(7)
					call subSymTensora(Res,T,Legi,QNi,degi,keepQN_)
			end select
			return
		end if
		if(degi.ne.1)then
			call writemess('ERROR in subTensor, the tensor is non-symmetry, degi should be 1',-1)
			call error_stop
		end if
		select case(T%getType())
			case(1)
				call subTensori(Res,T,Legi,QNi,keepQN_)
			case(2)
				call subTensors(Res,T,Legi,QNi,keepQN_)
			case(3)
				call subTensord(Res,T,Legi,QNi,keepQN_)
			case(4)
				call subTensorc(Res,T,Legi,QNi,keepQN_)
			case(5)
				call subTensorz(Res,T,Legi,QNi,keepQN_)
			case(6)
				call subTensorl(Res,T,Legi,QNi,keepQN_)
			case(7)
				call subTensora(Res,T,Legi,QNi,keepQN_)
		end select
		return
	end function

	function SubTensor2Func(T,Legi,dimi,keepQN_)result(Res)
		type(Tensor)::Res
		class(Tensor),intent(in) ::T
		integer,intent(in)::Legi,dimi
		logical,optional,intent(in)::keepQN_
		if(T%getSymmetryFlag())then
			call writemess('ERROR in subTensor, the tensor is of symmetry',-1)
			call writemess(' NO NOT input the degi yet',-1)
			call error_stop
		end if
		select case(T%getType())
			case(1)
				call subTensori(Res,T,Legi,dimi,keepQN_)
			case(2)
				call subTensors(Res,T,Legi,dimi,keepQN_)
			case(3)
				call subTensord(Res,T,Legi,dimi,keepQN_)
			case(4)
				call subTensorc(Res,T,Legi,dimi,keepQN_)
			case(5)
				call subTensorz(Res,T,Legi,dimi,keepQN_)
			case(6)
				call subTensorl(Res,T,Legi,dimi,keepQN_)
			case(7)
				call subTensora(Res,T,Legi,dimi,keepQN_)
		end select
		return
	end function