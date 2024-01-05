	subroutine get_a_Data(A,ith,outData)
		class(Tensor),intent(in)::A
		integer,intent(in)::ith
		class(*),intent(inout)::outData
		class(*),pointer::p
		if(.not.A%Data%getFlag())then
			call writemess('ERROR, the tensor is empty',-1)
			call error_stop
		end if
		call ClassPointer1DFunc(A%Data%ClassData,ith,p)
		call CopyData(outData,p)
	end subroutine
	subroutine get_a_Data_Vec(A,vec,Res)
		class(Tensor),intent(in)::A
		integer,intent(in)::vec(:)
		class(*),intent(inout)::Res
		integer,pointer::dim(:)
		class(*),pointer::p2(:,:),p3(:,:,:),p4(:,:,:,:),p0
		integer::rank,index
		if(.not.A%Data%getFlag())then
			call writemess('ERROR, the tensor is empty',-1)
			call error_stop
		end if
		if(A%getSymmetryFlag())then
			call writemess('ERROR in getting the element in the tensor, the tensor is of symmetry',-1)
			call error_stop
		end if
		rank=A%getRank()
		if(size(vec).ne.rank)then
			call writemess('ERROR in getting the element in the tensor',-1)
			call writemess('size(vec)='+size(vec),-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call error_stop
		end if
		call A%pointDim(dim)
		call check_indices_and_dim(dim,vec)
		select case(rank)
			case(1)
				call ClassPointer1DFunc(A%Data%ClassData,vec(1),p0)
			case(2)
				call A%ClassPointer(p2)
				call ClassPointer2DFunc(p2,vec(1),vec(2),p0)
			case(3)
				call A%ClassPointer(p3)
				call ClassPointer3DFunc(p3,vec(1),vec(2),vec(3),p0)
			case(4)
				call A%ClassPointer(p4)
				call ClassPointer4DFunc(p4,vec(1),vec(2),vec(3),vec(4),p0)
			case default
				index=addressToIndes(vec,dim)
				call ClassPointer1DFunc(A%Data%ClassData,index,p0)
		end select
		call CopyData(Res,p0)
		return
	end subroutine
	subroutine get_a_Block_Data(A,blocki,i,Res)
		class(Tensor),intent(in)::A
		integer,intent(in)::blocki,i
		class(*),intent(inout)::Res
		class(*),pointer::p(:)
		class(*),pointer::p0
		if(.not.A%Data%getFlag())then
			call writemess('ERROR, the tensor is empty',-1)
			call error_stop
		end if
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in getting the element in the tensor, the tensor is not of symmetry',-1)
			call error_stop
		end if
		if(blocki.gt.A%getTotalBlock())then
			call writemess('ERROR in getting the element in the tensor',-1)
			call writemess('blocki='+blocki,-1)
			call error_Stop
		end if
		if(i.gt.A%getTotalData())then
			call writemess('ERROR in getting the element in the tensor',-1)
			call writemess('blocki='+blocki,-1)
			call writemess('i='+i,-1)
			call writemess('A%getTotalData(blocki)='+A%getTotalData(blocki),-1)
			call error_Stop
		end if
		if(.not.A%getFlag(blocki))then
			call writemess('ERROR, the block in the tensor is empty',-1)
			call error_stop
		end if
		call A%ClassPointer(p,blocki)
		if(.not.associated(p))then
			call writemess('ERROR in getValueSymFUNCNAME',-1)
			call error_stop
		end if

		call ClassPointer1DFunc(p,i,p0)
		call CopyData(Res,p0)
		return
	end subroutine
	subroutine get_all_Data(A,Res)
		class(Tensor),intent(in)::A
		class(*),intent(inout)::Res(:)
		if(.not.A%Data%getFlag())then
			call writemess('ERROR, the tensor is empty',-1)
			call error_stop
		end if
		call FastCopyArray(Res,A%Data%ClassData,A%getTotalData())
		return
	end subroutine 
	subroutine get_Blocki_data(A,ith,Res)
		class(Tensor),intent(in)::A
		integer,intent(in)::ith
		class(*),intent(inout)::Res(:)
		class(*),pointer::p(:)
		if(.not.A%Data%getFlag())then
			call writemess('ERROR, the tensor is empty',-1)
			call error_stop
		end if
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in getting the element in the tensor, the tensor is not of symmetry',-1)
			call error_stop
		end if
		if(.not.A%getFlag(ith))then
			call writemess('ERROR, the block in the tensor is empty',-1)
			call error_stop
		end if
		call A%ClassPointer(p,ith)
		if(.not.associated(p))then
			call writemess('ERROR in getValueSymFUNCNAME',-1)
			call error_stop
		end if
		call FastCopyArray(Res,p,A%getTotalData(ith))
		return
	end subroutine
	subroutine get_Blockvec_data(A,vec,Res)
		class(Tensor),intent(in)::A
		integer,intent(in)::vec(:)
		class(*),intent(inout)::Res(:)
		class(*),pointer::p(:)
		if(.not.A%Data%getFlag())then
			call writemess('ERROR, the tensor is empty',-1)
			call error_stop
		end if
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in getting the element in the tensor, the tensor is not of symmetry',-1)
			call error_stop
		end if
		if(.not.A%getFlag(vec))then
			call writemess('ERROR, the block in the tensor is empty',-1)
			call error_stop
		end if
		call A%ClassPointer(p,vec)
		if(.not.associated(p))then
			call writemess('ERROR in getValueSymFUNCNAME',-1)
			call error_stop
		end if
		call FastCopyArray(Res,p,classSize(p))
		return
	end subroutine
