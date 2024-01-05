	function getValue0FUNCNAME(A,i)Result(Res)
		class(Tensor),intent(in)::A
		DATATYPE::Res
		integer,intent(in)::i
		if(.not.A%Data%getFlag())then
			call writemess('ERROR, the tensor is empty',-1)
			call error_stop
		end if
		if(i.gt.A%getTotalData())then
			call writemess('ERROR in getting the element in the tensor',-1)
			call writemess('i='+i,-1)
			call writemess('A%getTotalData()='+A%getTotalData(),-1)
			call error_Stop
		end if
		call getaValue(Res,A%Data%ClassData,i)
		return
	end function
	function getValueSymFUNCNAME(A,blocki,i)Result(Res)
		class(Tensor),intent(in)::A
		DATATYPE::Res
		integer,intent(in)::blocki,i
		class(*),pointer::p(:)
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
		call getaValue(Res,p,i)
		return
	end function

	function getValueVecFUNCNAME(A,vec)Result(Res)
		class(Tensor),intent(in)::A
		DATATYPE::Res
		integer,intent(in)::vec(:)
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
				call getaValue(Res,A%Data%ClassData,vec(1)) 
			case(2)
				call A%ClassPointer(p2)
				call ClassPointer2DFunc(p2,vec(1),vec(2),p0)
				call CopyData(Res,p0)
			case(3)
				call A%ClassPointer(p3)
				call ClassPointer3DFunc(p3,vec(1),vec(2),vec(3),p0)
				call CopyData(Res,p0)
			case(4)
				call A%ClassPointer(p4)
				call ClassPointer4DFunc(p4,vec(1),vec(2),vec(3),vec(4),p0)
				call CopyData(Res,p0)
			case default
				index=addressToIndes(vec,dim)
				call getaValue(Res,A%Data%ClassData,index) 
		end select
		return
	end function

	function getValueSymVecFUNCNAME(A,Blockvec,vec)Result(Res)
		class(Tensor),intent(in)::A
		DATATYPE::Res
		integer,intent(in)::Blockvec(:),vec(:)
		integer,pointer::dim(:)
		class(*),pointer::p(:),p2(:,:),p3(:,:,:),p4(:,:,:,:),p0
		integer::rank,index
		integer,pointer::BlockDim(:)
		if(.not.A%Data%getFlag())then
			call writemess('ERROR, the tensor is empty',-1)
			call error_stop
		end if
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in getting the element in the tensor, the tensor is not of symmetry',-1)
			call error_stop
		end if
		rank=A%getRank()
		if(size(vec).ne.rank)then
			call writemess('ERROR in getting the element in the tensor',-1)
			call writemess('size(vec)='+size(vec),-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call error_stop
		end if
		if(size(Blockvec).ne.rank)then
			call writemess('ERROR in getting the element in the tensor',-1)
			call writemess('size(Blockvec)='+size(Blockvec),-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call error_stop
		end if
		if(.not.A%getFlag(Blockvec))then
			call writemess('ERROR, the block in the tensor is empty',-1)
			call error_stop
		end if
		call A%pointDim(dim)
		call check_indices_and_dim(dim,Blockvec)
		select case(rank)
			case(1)
				call A%ClassPointer(p,Blockvec)
				if(.not.associated(p))then
					call writemess('ERROR in getValueSymVecFUNCNAME',-1)
					call error_stop
				end if
				call getaValue(Res,p,vec(1)) 
			case(2)
				call A%ClassPointer(p2,Blockvec)
				if(.not.associated(p2))then
					call writemess('ERROR in getValueSymVecFUNCNAME',-1)
					call error_stop
				end if
				call ClassPointer2DFunc(p2,vec(1),vec(2),p0)
				call CopyData(Res,p0)
			case(3)
				call A%ClassPointer(p3,Blockvec)
				if(.not.associated(p3))then
					call writemess('ERROR in getValueSymVecFUNCNAME',-1)
					call error_stop
				end if
				call ClassPointer3DFunc(p3,vec(1),vec(2),vec(3),p0)
				call CopyData(Res,p0)
			case(4)
				call A%ClassPointer(p4,Blockvec)
				if(.not.associated(p4))then
					call writemess('ERROR in getValueSymVecFUNCNAME',-1)
					call error_stop
				end if
				call ClassPointer4DFunc(p4,vec(1),vec(2),vec(3),vec(4),p0)
				call CopyData(Res,p0)
			case default
				call A%ClassPointer(p,Blockvec)
				if(.not.associated(p))then
					call writemess('ERROR in getValueSymVecFUNCNAME',-1)
					call error_stop
				end if
				call WorkingMemory%Check()
				call WorkingMemory%allocate(1,rank)
				call WorkingMemory%get_memory(BlockDim,rank)
				BlockDim=A%getBlockDim(Blockvec)
				call check_indices_and_dim(BlockDim,vec,size(p))
				index=addressToIndes(vec,BlockDim)
				call WorkingMemory%free()
				if(index.gt.size(p))then
					call writemess('ERROR in getValueSymVecFUNCNAME',-1)
					call error_stop
				end if
				call getaValue(Res,p,index) 
		end select
		return
	end function

	function getValueAllFUNCNAME(A)Result(Res)
		class(Tensor),intent(in)::A
		DATATYPE2,allocatable::Res(:)
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		logical,pointer::lp(:)
		character(len=characterlen),pointer::ap(:)
		integer::TotalData,i
		if(.not.A%Data%getFlag())then
			call writemess('ERROR, the tensor is empty',-1)
			call error_stop
		end if
		TotalData=A%getTotalData()
		allocate(Res(TotalData))
		call FastCopyArray(Res,A%Data%ClassData,TotalData)
		return
	end function