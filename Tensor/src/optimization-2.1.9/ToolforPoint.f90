module Tool_For_Point
	use Tools
	use Basic_Tools
	use Tensor_Tools
	use dimension_Tools
	implicit none
contains

	function TensorArrayInfo(A)result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A(:)
		integer::i,lenA,si,ei
		integer,pointer::ip(:,:)
		lenA=size(A)
		call Res%allocate([lenA,2],'integer')
		call Res%pointer(ip)
		call Res%zero()
		ei=0
		si=0
		do i=1,lenA
			if(A(i)%getFlag())then
				si=ei+1
				ei=ei+A(i)%getTotalData()
				ip(i,:)=[si,ei]
			else
				ip(i,:)=[0,0]
			end if
		end do
		return
	end function
	function TensorArrayLegInfo(A)result(Res)
		type(dimension),allocatable::Res(:)
		type(Tensor),intent(in)::A(:)
		integer::i
		allocate(Res(size(A)))
		do i=1,size(A)
			Res(i)=A(i)%Dimension
		end do
		return
	end function

	subroutine ArrayToTensor(outT,A,Info,classType)
		type(Tensor),intent(inout)::outT
		type(Tensor),intent(in)::A(:),Info
		character(len=*),optional::classType
		integer::i,maxei,si,ei
		integer,pointer::ip(:,:)
		class(*),pointer::Ap(:),AllData(:)
		maxei=info%imax()
		if(present(classType))then
			call outT%allocate([maxei],classType)
		else
			if(maxei.ne.outT%getTotalData())then
				call writemess('ERROR in ArrayToTensor, total length')
				call error_stop
			end if
		end if
		call outT%zero()
		call outT%ClassPointer(AllData)
		call Info%pointer(ip)
		do i=1,Info%dim(1)
			si=ip(i,1)
			ei=ip(i,2)
			if(ei.gt.0)then
				call A(i)%ClassPointer(Ap)
				if(ClassSize(Ap).ne.(ei-si+1))then
					call writemess('ERROR in ArrayToTensor',-1)
					call error_stop
				end if
				call ModifySomeValue(AllData,Ap,si,ei)
			end if
		end do
		return
	end subroutine

	subroutine TensorToArray(outA,T,Info,DimInf)
		type(Tensor),intent(inout)::outA(:)
		type(Tensor),intent(in)::T,Info
		type(Tensor),optional,intent(in)::DimInf(:)
		integer::i,totaldata,si,ei,checki
		integer,pointer::ip(:,:)
		class(*),pointer::Ap(:),AllData(:),TMPp(:)
		call T%ClassPointer(AllData)
		call Info%pointer(ip)
		do i=1,Info%dim(1)
			si=ip(i,1)
			ei=ip(i,2)
			if(ei.gt.0)then
				if(.not.outA(i)%getFlag())then
					call writemess('ERROR in TensorToArray',-1)
					call error_stop
				end if
				if(present(DimInf))then
					do checki=1,outA(i)%getRank()
						if(DimInf(i)%getName(checki).nequ.outA(i)%getName(checki))then
							call writemess('ERROR in TensorToArray,legs info',-1)
							call error_stop
						end if
					end do
				end if
				call outA(i)%ClassPointer(Ap)
				totaldata=outA(i)%getTotalData()
				if(ClassSize(Ap).ne.(ei-si+1))then
					call writemess('ERROR in TensorToArray',-1)
					call error_stop
				end if
				call ClassPointer1DFunc(AllData,si,ei,TMPp)
				call CopyArray(Ap,TMPp,totaldata)
			else
				if(outA(i)%getFlag())then
					call writemess('ERROR in TensorToArray',-1)
					call error_stop
				end if
			end if
		end do
		return
	end subroutine
end module
