module targetFunction_Tools
	use Density_type
	use Tensor_Tools
	use dimension_tools
	use Tools
	use General_Optimization_Tools
	use General_Optimization_Element
	use Optimization_Tensor_interface
	use element_Tools
	use Parameter_Tools
	implicit none
	logical,private::MLE_matrix_flag=.true.
	type, extends(OptimEngine_structure) :: OptimEngine
		type(Tensor),allocatable::DataValue(:)
		type(Tensor),allocatable::Umatrix(:)
		type(Tensor),allocatable::density(:)
		type(OptData)::DeltaA
		integer::ParameterLen=0
		integer::phy_length=0
		integer::num_of_density=0
		integer::delta_step_for_save=0
		real*8::diag_small_number=1d-6
		!data define below for MLE
		type(Tensor)::L
		!data define below for MLE_numerical
		real*8::deltax=0.001

	end type
	
contains

	function divergencesFunc(OE)
		real*8::divergencesFunc
		class(OptimEngine_structure),intent(in)::OE
		real*8::PlogP
		integer::i,j
		select type(OE)
		type is (OptimEngine)
			divergencesFunc=0d0
			do i=1,OE%num_of_density
				do j=1,OE%phy_length
					if(OE%DataValue(i)%di(j).gt.1d-16)then
						PlogP=OE%DataValue(i)%di(j)*log(OE%DataValue(i)%di(j)/OE%density(i)%di([j,j]))
						divergencesFunc=divergencesFunc+PlogP
					end if
				end do
			end do
		end select
		return
	end function

	subroutine MVFunc(OE,outAp_,inp_)!outAp=A*inp
		class(OptimEngine_structure),target, intent(inout) :: OE
		class(OptimElement_structure),target,intent(inout)::outAp_
		class(OptimElement_structure),target,intent(in)::inp_
		type(OptData),pointer::inp,inoutp
		type(Tensor)::TMPDelta
		type(Tensor)::tmp,allname,temp
		integer::i
		complex*16,pointer::zp(:)
		call PointToOptData(inp_,inp)
		call PointToOptData(outAp_,inoutp)
		select type(OE)
		type is (OptimEngine)

		if(.not.allocated(inoutp%Alldata))allocate(inoutp%Alldata(9))
		do i=1,9
			tmp=contract(OE%DeltaA%Alldata(i),'A.S1','A.S2')
			tmp=tmp.kron.OE%density(1)
			call tmp%permute(OE%DeltaA%Alldata(i)%getName())
			TMPDelta=OE%DeltaA%Alldata(i) - tmp
			
			temp=TMPDelta
			call temp%pointer(zp)
			zp=conjg(zp)
			allname=inp%Alldata(i)%getName()
			inoutp%Alldata(i)=contract(TMPDelta,allname%ai(),inp%Alldata(i),allname%ai())
			call inoutp%Alldata(i)%contract(temp,['A.S1','A.S2'],['A.S1','A.S2'])

			call inoutp%Alldata(i)%permute(allname%ai())
			if(OE%diag_small_number.lt.0)then
				inoutp%Alldata(i)=dble(inoutp%Alldata(i)-(OE%diag_small_number*inp%Alldata(i)%dmax('abs')*inp%Alldata(i)))
			else
				inoutp%Alldata(i)=dble(inoutp%Alldata(i)+(OE%diag_small_number*inp%Alldata(i)))
			end if
		end do
		end select
		return
	end subroutine 

	subroutine targetFunction(OE,outVal,outGradient_,point_)
		class(OptimEngine_structure),target, intent(inout) :: OE
		real*8,target,intent(inout)::outVal
		class(OptimElement_structure),target,intent(inout)::outGradient_
		class(OptimElement_structure),target,intent(in)::point_
		type(OptData),pointer::outGradient,point
		type(Tensor),pointer::W1,W2,U1,U2,C1,C2,B1,B2,D
		type(Tensor)::ULRToolTensor,temp
		type(OptData)::traceGradientData
		integer::i,j
		real*8::PlogP
		type(Tensor)::rhoAndGra(10),ULTensor,URTensor
		call PointToOptData(outGradient_,outGradient)
		call PointToOptData(point_,point)
		select type(OE)
		type is (OptimEngine)
			if(OE%num_of_density.eq.0)then
				call writemess('ERROR in gradientFunc')
				call error_stop
			end if
			call point%pointer(W1,W2,U1,U2,C1,C2,B1,B2,D)
			rhoAndGra=DensityMatrixAndGradient(W1,W2,U1,U2,C1,C2,B1,B2,D)
			call OE%DeltaA%storeData(rhoAndGra(2:10))

			call traceGradientData%allocate()
			do i=1,9
				traceGradientData%Alldata(i)=contract(OE%DeltaA%Alldata(i),'A.S1','A.S2')*OE%num_of_density
			end do
			

			do i=1,OE%num_of_density
				OE%density(i)=OE%Umatrix(i)*rhoAndGra(1)*(.H.OE%Umatrix(i))
				call OE%density(i)%setName(1,'A.S1')
				call OE%density(i)%setName(2,'A.S2')
			end do
			ULTensor=transform2Tensor(OE%DataValue,OE%density,OE%Umatrix)
			URTensor=AddColH(OE%Umatrix)
			call ULRToolTensor%contract(ULTensor,['U.L','U.1'],URTensor,['U+.L','U+.1'])

			call outGradient%allocate()
			do i=1,9
				call temp%contract(ULRToolTensor,['U.2 ','U+.2'],OE%DeltaA%Alldata(i),['A.S1','A.S2'])
				call temp%permute(traceGradientData%Alldata(i)%getName())
				outGradient%Alldata(i)=dble(traceGradientData%Alldata(i)-temp)
			end do
			

			outVal=0d0
			do i=1,OE%num_of_density
				do j=1,OE%phy_length
					PlogP=OE%DataValue(i)%di(j)*log(OE%density(i)%di([j,j]))
					outVal=outVal-PlogP
					!PlogP=OE%DataValue(i)%di(j)*log(OE%DataValue(i)%di(j)/OE%density(i)%di([j,j]))
					!outVal=outVal+PlogP
				end do
			end do
		class default
			call writemess('ERROR in targetFunction',-1)
			call error_stop
		end select
	end subroutine 

	subroutine EnergyFunction(OE,outVal,point_)
		class(OptimEngine_structure),target, intent(inout) :: OE
		real*8,target,intent(inout)::outVal
		class(OptimElement_structure),target,intent(in)::point_
		type(OptData),pointer::point
		type(Tensor),pointer::W1,W2,U1,U2,C1,C2,B1,B2,D
		integer::i,j
		real*8::PlogP
		type(Tensor)::rhoAndGra(10)
		call PointToOptData(point_,point)
		select type(OE)
		type is (OptimEngine)
			if(OE%num_of_density.eq.0)then
				call writemess('ERROR in gradientFunc')
				call error_stop
			end if
			call point%pointer(W1,W2,U1,U2,C1,C2,B1,B2,D)
			rhoAndGra=DensityMatrixAndGradient(W1,W2,U1,U2,C1,C2,B1,B2,D)
			do i=1,OE%num_of_density
				OE%density(i)=OE%Umatrix(i)*rhoAndGra(1)*(.H.OE%Umatrix(i))
				call OE%density(i)%setName(1,'A.S1')
				call OE%density(i)%setName(2,'A.S2')
			end do
			outVal=0d0
			do i=1,OE%num_of_density
				do j=1,OE%phy_length
					PlogP=OE%DataValue(i)%di(j)*log(OE%density(i)%di([j,j]))
					outVal=outVal-PlogP
					!PlogP=OE%DataValue(i)%di(j)*log(OE%DataValue(i)%di(j)/OE%density(i)%di([j,j]))
					!outVal=outVal+PlogP
				end do
			end do
		class default
			call writemess('ERROR in targetFunction',-1)
			call error_stop
		end select
	end subroutine 




	function PsigmaOverRhoSigmeU(DataValue,Rho,Umatrix)Result(Res)
		type(Tensor)::Res
		type(Tensor)::DataValue,Rho,Umatrix
		complex*16,pointer::zp(:,:)
		integer::i
		Res=Umatrix
		call Res%pointer(zp)
		do i=1,Res%Dim(1)
			zp(i,:)=zp(i,:)*(DataValue%di(i)/Rho%di([i,i]))
		end do
		return
	end function

	function transform2Tensor(DataValue,Rho,Umatrix)Result(Res)
		type(Tensor)::Res
		type(Tensor)::DataValue(:),Rho(:),Umatrix(:)
		integer::i
		type(Tensor)::temp
		temp=PsigmaOverRhoSigmeU(DataValue(1),Rho(1),Umatrix(1))
		Res=temp
		do i=2,size(Rho)
			temp=PsigmaOverRhoSigmeU(DataValue(i),Rho(i),Umatrix(i))
			Res=addTensorCol(Res,temp)
		end do
		call Res%setName(Res%getRank(),'U.L')
		return
	end function
	function addCol(A)Result(Res)
		type(Tensor)::Res
		type(Tensor)::A(:)
		integer::i
		type(Tensor)::temp
		Res=A(1)
		do i=2,size(A)
			Res=addTensorCol(Res,A(i))
		end do
		call Res%setName(Res%getRank(),'U.L')
		return
	end function
	function addColH(A)Result(Res)
		type(Tensor)::Res
		type(Tensor)::A(:)
		integer::i
		type(Tensor)::temp
		Res=.con.A(1)
		do i=2,size(A)
			Res=addTensorCol(Res,.con.A(i))
		end do
		call Res%setName(Res%getRank(),'U+.L')
		call Res%setName('U+')
		return
	end function

	subroutine inputcheck(point)
		type(OptData),intent(inout)::point
		type(Tensor),pointer::W1,W2,U1,U2,C1,C2,B1,B2,D
		real*8::norm2
		call point%pointer(W1,W2,U1,U2,C1,C2,B1,B2,D)
		call find_exp_factor(W1,W2,U1,U2,C1,C2,B1,B2,D)
		return
	end subroutine


	subroutine initialOptim(A,DataValue,Umatrix,parameter)
		class(OptimEngine), intent(inout) :: A
		type(List)::Parameter(3)
		type(Tensor),intent(in)::DataValue(:),Umatrix(:)
		integer::ParameterLen
		logical::alive
		integer::i
		if(Parameter(3)%ai('model').equ.'RBM')then
			call initialElementInterFace()
			call A%set_target_function(targetFunction)
			call A%set_Energy_Function(EnergyFunction)
			call A%set_MV_func(MVFunc)
		else if(Parameter(3)%ai('model').equ.'MLE') then
			call initialTensorInterFace()
			call A%set_target_function(numerical_differentiation_MLE)
			call A%set_Energy_Function(numerical_Energy_MLE)
		else
			call writemess('ERROR model')
			call error_stop
		end if
		call A%set_stop_error(Parameter(2)%di('stop_error'))
		call A%Set_inStep3Func(inStepTimeCounter)
		call A%set_CG_direction_flag(Parameter(1)%ii('CG_direction_flag'))
		call A%set_method(Parameter(3)%ai('search_method'))
		call A%set_CGLLS_parameter(Parameter(1)%ii('CGLLSnum'),Parameter(2)%di('CGLLSerror'))
		call A%set_LBFGS_length(Parameter(1)%ii('LBFGS_length'))
		A%diag_small_number=Parameter(2)%di('diag_small_number')
		call writemess('run the method:'+Parameter(3)%ai('search_method'))
		if(Parameter(1)%ii('lineSearchNum').gt.0)then
			call writemess('Going to run the line search, lineSearchNum='+Parameter(1)%ii('lineSearchNum'))
		end if
		parameterLen=Parameter(1)%ii('parameterLen')
		if(parameterLen.eq.0)then
			parameterLen=size(DataValue)*2
		else if(parameterLen.lt.0)then
			parameterLen=size(DataValue)*(parameterLen)*(-1)
		end if

		A%ParameterLen=ParameterLen
		A%num_of_density=size(DataValue)
		A%phy_length=DataValue(1)%getTotalData()

		A%delta_step_for_save=Parameter(1)%ii('delta_step_for_save')
		allocate(A%DataValue(A%num_of_density))
		allocate(A%Umatrix(A%num_of_density))
		allocate(A%density(A%num_of_density))
		A%DataValue=DataValue
		A%Umatrix=Umatrix
		do i=1,A%num_of_density
			call A%Umatrix(i)%setName('U')
		end do

		if(Parameter(3)%ai('read_address').equ.'null')then
			call open_file(1111,'output/divergences.dat','replace')
			close(1111)
		else
			inquire(file='output/divergences.dat',exist=alive)
			if(.not.alive)then
				call open_file(1111,'output/divergences.dat','replace')
				close(1111)
			end if
		end if
		return
	end subroutine
	subroutine inputDataCheck(DataValue,Umatrix)
		type(Tensor),intent(in)::DataValue(:),Umatrix(:)
		integer::dataLength,i
		real*8::testValue
		type(Tensor)::UU
		call writemess('Check the input Data...')
		dataLength=size(DataValue)
		if(dataLength.ne.size(Umatrix))then
			call writemess('ERROR in inputDataCheck')
			call error_stop
		end if
		do i=1,dataLength
			call checkDataValue(DataValue(i))
			testValue=DataValue(i)%dsum()
			if(testValue.nequ.1d0)then
				call writemess('ERROR in DataValue('+i+')')
				call writemess('sum(DataValue('+i+'))='+testValue)
				call error_stop
			end if

			UU=Umatrix(i)*(.H.Umatrix(i))
			if(.not.equal_identity(UU))then
				call writemess('ERROR in Umatrix('+i+')')
				call writemess('Umatrix('+i+')*Umatrix('+i+')^+ !=identity')
				call UU%print
				call error_stop
			end if
			UU=(.H.Umatrix(i))*Umatrix(i)
			if(.not.equal_identity(UU))then
				call writemess('ERROR in Umatrix('+i+')')
				call writemess('Umatrix('+i+')^+ * Umatrix('+i+') !=identity')
				call UU%print
				call error_stop
			end if
		end do
		call writemess('Check the input Data...done')
		return
	end subroutine

	function equal_identity(A)
		type(Tensor)::A
		logical::equal_identity
		integer::i,j,M,N
		complex*16::va
		equal_identity=.true.
		M=A%dim(1)
		N=A%dim(2)
		do i=1,M
			do j=1,N
				va=A%zi([i,j])
				if(i.eq.j)then
					if(abs(va-dcmplx(1d0,0d0)).gt.1d-8)then
						equal_identity=.false.
						return
					end if
				else
					if(abs(va).gt.1d-8)then
						equal_identity=.false.
						return
					end if
				end if
			end do
		end do
		return
	end function

	subroutine checkDataValue(A)
		type(Tensor)::A
		integer::i
		real*8::va
		do i=1,A%getTotalData()
			va=A%di(i)
			if(va.lt.0)then
				if(abs(va).gt.1d-16)then
					call writemess('ERROR in checkDataValue,va='+va)
					call error_stop
				end if
			end if
		end do
		return
	end subroutine
	

	subroutine inStepTimeCounter(OE,Value,direction,Gradient,point,ith,t)
		class(OptimEngine_structure),target, intent(inout) :: OE
		real*8,target,intent(inout)::Value
		class(OptimElement_structure),target,intent(inout)::Gradient
		class(OptimElement_structure),target,intent(inout)::point
		class(OptimElement_structure),target,intent(inout)::direction
		integer,target,intent(in)::ith
		real*8,target::t
		integer,save::counter=0
		real*8::normGra,normP,diver
		counter=counter+1
		call time_calculator()
		call open_file(1111,'output/divergences.dat','old','end')
		call norm2Func(Gradient,normGra)
		call norm2Func(point,normP)
		diver=divergencesFunc(OE)
		write(1111,*)Value,t,normGra,normP,diver
		close(1111)
		select type(OE)
			type is (OptimEngine)
				select type(point)
					type is (OptData)
						if(counter.eq.OE%delta_step_for_save)then
							call open_file(114,'save/temp_point.dat','replace')
							call point%write(114)
							close(114)
						end if
						if(counter.gt.OE%delta_step_for_save)counter=0
					type is (OptimTensor)
						if(counter.eq.OE%delta_step_for_save)then
							call open_file(114,'save/temp_point.dat','replace')
							call point%data%write(114)
							close(114)
						end if
						if(counter.gt.OE%delta_step_for_save)counter=0
				end select
		end select
	end subroutine 



	!************************************************************************
	!
	!          code below for MLE
	!
	!************************************************************************


	subroutine PointToOptimTensor(OT,T)
		class(OptimElement_structure),target::OT
		type(Tensor),pointer::T
		select type(OT)
		type is (OptimTensor)
			T=>OT%Data
		class default
			call writemess('ERROR in PointToOptimTensor',-1)
			call writemess('Input class is not OptData',-1)
			call error_stop
		end select
		return
	end subroutine



	function PhysicalRho(L)result(Rho)
		type(Tensor)::Rho
		type(Tensor),intent(in)::L
		real*8::norm2
		Rho=(.H.L)*L
		norm2=Rho%dtrace()
		if(norm2.ne.0d0)then
			Rho=Rho/Rho%dtrace()
		end if
		call Rho%setName('rho')
		return
	end function



	subroutine initialOptimTensor(OT,phy_length,elementScal)
		class(OptimTensor),intent(inout)::OT
		integer,intent(in)::phy_length
		real*8,intent(in)::elementScal
		integer::i,j,total
		total=0
		do i=1,phy_length
			total=total+1
			do j=1,i-1
				total=total+2
			end do
		end do
		call OT%data%allocate([total],'real*8')
		call OT%data%random([-elementScal,elementScal])
		return
	end subroutine

	function transform2L_MEL(OData,Dim)result(L)
		type(Tensor)::L
		type(Tensor),intent(in)::OData
		integer,intent(in)::dim
		integer::i,j,k
		complex*16,pointer::Ldata(:,:)
		real*8,pointer::Op(:)
		real*8::tmpr,tmpi

		call L%allocate([dim,dim],'complex*16')
		call L%zero()
		call L%pointer(LData)
		call OData%pointer(Op)
		k=0
		do i=1,dim
			do j=1,i
				if(i.eq.j)then
					k=k+1
					if(k.gt.OData%getTotalData())then
						call writemess('ERROR in transform2L_MEL',-1)
						call error_stop
					end if
					LData(i,j)=Op(k)
				else
					k=k+1
					if(k.gt.OData%getTotalData())then
						call writemess('ERROR in transform2L_MEL',-1)
						call error_stop
					end if
					tmpr=Op(k)

					k=k+1
					if(k.gt.OData%getTotalData())then
						call writemess('ERROR in transform2L_MEL',-1)
						call error_stop
					end if
					tmpi=Op(k)
					LData(i,j)=dcmplx(tmpr,tmpi)
				end if
			end do
		end do
		call L%setName('L')
		return
	end function

	subroutine numerical_differentiation_MLE(OE,outVal,outGradient_,point_)
		class(OptimEngine_structure),target, intent(inout) :: OE
		real*8,target,intent(inout)::outVal
		class(OptimElement_structure),target,intent(inout)::outGradient_
		class(OptimElement_structure),target,intent(in)::point_
		type(Tensor),pointer::point,outGradient
		type(Tensor)::deltaPoint
		real*8::TMPVal,gra
		integer::i,j
		select type(OE)
		type is (OptimEngine)

			call PointToOptimTensor(outGradient_,outGradient)
			call PointToOptimTensor(point_,point)
			call Energy_for_numerical_differentiation_MLE(OE,outVal,point)
			call outGradient%allocate(point)
			call outGradient%zero()
			do i=1,point%getTotalData()
				deltaPoint=point
				call deltaPoint%setValue(i,point%di(i)+OE%deltax)
				call Energy_for_numerical_differentiation_MLE(OE,TMPVal,deltaPoint)
				gra=(TMPVal-outVal)/OE%deltax
				call outGradient%setValue(i,gra)
			end do
		class default
			call writemess('ERROR in numerical_differentiation_MLE',-1)
			call error_stop
		end select
	end subroutine

	subroutine numerical_Energy_MLE(OE,outVal,point_)
		class(OptimEngine_structure),target, intent(inout) :: OE
		real*8,target,intent(inout)::outVal
		class(OptimElement_structure),target,intent(in)::point_
		type(Tensor),pointer::point
		type(Tensor)::deltaPoint
		select type(OE)
		type is (OptimEngine)

			call PointToOptimTensor(point_,point)
			call Energy_for_numerical_differentiation_MLE(OE,outVal,point)

		class default
			call writemess('ERROR in numerical_differentiation_MLE',-1)
			call error_stop
		end select
	end subroutine

	subroutine Energy_for_numerical_differentiation_MLE(OE,outVal,point)
		class(OptimEngine_structure),target, intent(inout) :: OE
		real*8,target,intent(inout)::outVal
		type(Tensor),intent(in)::point
		integer::i,j
		real*8::PlogP
		type(Tensor)::rho
		select type(OE)
		type is (OptimEngine)
			if(OE%num_of_density.eq.0)then
				call writemess('ERROR in gradientFunc')
				call error_stop
			end if

			OE%L=transform2L_MEL(point,OE%phy_length)
			rho=PhysicalRho(OE%L)
			do i=1,OE%num_of_density
				OE%density(i)=OE%Umatrix(i)*rho*(.H.OE%Umatrix(i))
			end do
			outVal=0d0
			do i=1,OE%num_of_density
				do j=1,OE%phy_length
					PlogP=OE%DataValue(i)%di(j)*log(OE%density(i)%di([j,j]))
					outVal=outVal-PlogP
				end do
			end do
		class default
			call writemess('ERROR in targetFunction',-1)
			call error_stop
		end select
	end subroutine


end module
