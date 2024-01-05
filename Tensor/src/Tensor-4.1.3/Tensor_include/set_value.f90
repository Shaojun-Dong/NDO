
	!-------------use eigher in tensor or symtensor------------------------------------

	subroutine set_tensor_Some_value0(A,ith,val,jth)!use for tensor only
		class(Tensor),intent(inout)::A
		type(Tensor),intent(in)::val
		integer,intent(in)::ith(2)
		integer,intent(in)::jth(2)
		class(*),pointer::TMP(:)
		if(A%getSymmetryFlag())then
			call writemess('Can NOT setting the data to the tensor, the tensor is of symmetry',-1)
			call error_stop
		end if
		call ClassPointer1DFunc(Val%Data%ClassData,jth(1),jth(2),TMP)
		call ModifySomeValue(A%Data%ClassData,TMP,ith(1),ith(2))
		return
	end subroutine
	subroutine set_tensor_Some_value1(A,ith,val,jth)!use for tensor only
		class(Tensor),intent(inout)::A
		class(*),target,intent(in)::val(:)
		integer,intent(in)::ith(2)
		integer,optional,intent(in)::jth(2)
		class(*),pointer::TMP(:)
		if(A%getSymmetryFlag())then
			call writemess('Can NOT setting the data to the tensor, the tensor is of symmetry',-1)
			call error_stop
		end if
		select type(val)
			type is (Tensor)
				call writemess('ERROR in setting value to the ith element of the tensor',-1)
				call error_stop
			class default
				if(present(jth))then
					call ClassPointer1DFunc(Val,jth(1),jth(2),TMP)
				else
					TMP=>Val
				end if
				call ModifySomeValue(A%Data%ClassData,TMP,ith(1),ith(2))
		end select
		return
	end subroutine

	subroutine set_tensor_All_value0(A,val)!use for tensor only
		class(Tensor),intent(inout)::A
		class(*),intent(in)::val
		integer::totalData
		class(*),pointer::TMP0
		if(A%getSymmetryFlag())then
			call writemess('Can NOT setting the data to the tensor, the tensor is of symmetry',-1)
			call error_stop
		end if
		select type(val)
			type is (Tensor)
				totalData=val%getTotalData()
				if(totalData.eq.A%getTotalDAta())then
					call FastCopyArray(A%Data%ClassData,val%Data%ClassData,totalData)
				else if(totalData.eq.1)then
					call ClassPointer1DFunc(val%Data%ClassData,1,TMP0)
					call ModifyAllValue(A%Data%ClassData,TMP0,A%getTotalDAta())
				else
					call writemess('ERROR in setting value to the ith element of the tensor',-1)
					call writemess('size(block)='+A%getTotalDAta(),-1)
					call writemess('size(input)='+totalData,-1)
					call error_stop
				end if
			class default
				call ModifyAllValue(A%Data%ClassData,val,A%getTotalData())
		end select
		return
	end subroutine

	subroutine set_tensor_value0(A,ith,val)!use for tensor only
		class(Tensor),intent(inout)::A
		integer,intent(in)::ith
		class(*),intent(in)::val
		integer::totalData
		class(*),pointer::TMP0
		if(A%getSymmetryFlag())then
			call writemess('Can NOT setting the data to the tensor, the tensor is of symmetry',-1)
			call error_stop
		end if
		select type(val)
			type is (Tensor)
				totalData=val%getTotalData()
				if(totalData.ne.1)then
					call writemess('ERROR in setting value to the ith element of the tensor',-1)
					call error_stop
				end if
				call ClassPointer1DFunc(val%Data%ClassData,1,TMP0)
				call ModifyaValue(A%Data%ClassData,TMP0,ith)
			class default
				call ModifyaValue(A%Data%ClassData,val,ith)
		end select
		return
	end subroutine
	subroutine set_tensor_vec_value0(A,vec,val)!use for tensor only
		class(Tensor),intent(inout)::A
		integer,intent(in)::vec(:)
		class(*),intent(in)::val
		class(*),pointer::p,TMP0
		integer,pointer::dim(:)
		integer::totalData,index
		if(A%getSymmetryFlag())then
			call writemess('Can NOT setting the data to the tensor, the tensor is of symmetry',-1)
			call error_stop
		end if
		call A%pointDim(dim)
		call check_indices_and_dim(dim,vec)
		select case(size(vec))
			case(1)
				call ClassPointer1DFunc(A%Data%ClassData,vec(1),p)
			case(2)
				call ClassPointer2DFunc(A%getTotalData(),A%Data%ClassData,dim(1),dim(2),vec(1),vec(2),p)
			case(3)
				call ClassPointer3DFunc(A%getTotalData(),A%Data%ClassData,dim(1),dim(2),dim(3),vec(1),vec(2),vec(3),p)
			case(4)
				call ClassPointer4DFunc(A%getTotalData(),A%Data%ClassData,dim(1),dim(2),dim(3),dim(4),vec(1),vec(2),vec(3),vec(4),p)
			case default
				index=addressToIndes(vec,Dim)
				call ClassPointer1DFunc(A%Data%ClassData,index,p)
		end select 
		select type(val)
			type is (Tensor)
				totalData=val%getTotalData()
				if(totalData.ne.1)then
					call writemess('ERROR in setting value to the ith element of the tensor',-1)
					call error_stop
				end if
				call ClassPointer1DFunc(val%Data%ClassData,1,TMP0)
				call CopyData(p,TMP0)
			class default
				call CopyData(p,val)
		end select
		return
	end subroutine

	subroutine set_tensor_block_value1(A,ith,val)!use for symtensor only
		class(Tensor),intent(inout)::A
		integer,intent(in)::ith
		class(*),intent(in)::val(:)
		class(*),pointer::p(:)
		integer::totalData
		if(.not.A%getSymmetryFlag())then
			call writemess('Can NOT setting the data the symmetry tensor, the tensor is not of symmetry',-1)
			call error_stop
		end if
		select type(val)
			type is (Tensor)
				call writemess('ERROR in setting value to the ith block of the tensor',-1)
				call error_stop
			class default
				call A%Classpointer(p,ith)
				if(.not.associated(p))then
					call writemess('The block in the tensor is empty',-1)
					call writemess('Does this block obey the symmetry rule of if this block has been allocated?')
					call error_stop
				end if
				totalData=ClassSize(val)
				if(ClassSize(p).ne.totalData)then
					call writemess('ERROR in setting value to the ith block of the tensor',-1)
					call writemess('size(block)='+size(p),-1)
					call writemess('size(input)='+totalData,-1)
					call error_stop
				end if
				call FastCopyArray(p,val,totalData)
		end select
		return
	end subroutine

	subroutine set_tensor_block_value0(A,ith,val)!use for symtensor only
		class(Tensor),intent(inout)::A
		integer,intent(in)::ith
		class(*),intent(in)::val
		class(*),pointer::p(:),TMP0
		integer::totalData
		if(.not.A%getSymmetryFlag())then
			call writemess('Can NOT setting the data the symmetry tensor, the tensor is not of symmetry',-1)
			call error_stop
		end if
		select type(val)
			type is (Tensor)
				call A%Classpointer(p,ith)
				if(.not.associated(p))then
					call writemess('The block in the tensor is empty',-1)
					call writemess('Does this block obey the symmetry rule of if this block has been allocated?')
					call error_stop
				end if
				totalData=val%getTotalDAta()
				if(size(p).ne.totalData)then
					call writemess('ERROR in setting value to the ith block of the tensor',-1)
					call writemess('size(block)='+size(p),-1)
					call writemess('size(input)='+totalData,-1)
					call error_stop
				end if
				call FastCopyArray(p,val%Data%classData,totalData)
			class default
				call A%Classpointer(p,ith)
				if(.not.associated(p))then
					call writemess('The block in the tensor is empty',-1)
					call writemess('Does this block obey the symmetry rule of if this block has been allocated?')
					call error_stop
				end if
				totalData=1
				if(ClassSize(p).ne.totalData)then
					call writemess('ERROR in setting value to the ith block of the tensor',-1)
					call writemess('size(block)='+size(p),-1)
					call writemess('size(input)='+totalData,-1)
					call error_stop
				end if
				call ClassPointer1DFunc(p,1,TMP0)
				call CopyData(TMP0,val)
		end select
		return
	end subroutine

	subroutine set_tensor_blockvec_value0(A,vec,val)!use for tensor only
		class(Tensor),intent(inout)::A
		integer,intent(in)::vec(:)
		class(*),intent(in)::val
		class(*),pointer::p(:),TMP0
		integer,pointer::dim(:)
		integer::totalData
		if(.not.A%getSymmetryFlag())then
			call writemess('Can NOT setting the data the symmetry tensor, the tensor is not of symmetry',-1)
			call error_stop
		end if
		call A%pointDim(dim)
		call check_indices_and_dim(dim,vec)
		call A%ClassPointer(p,vec) 
		if(.not.associated(p))then
			call writemess('The block in the tensor is empty',-1)
			call writemess('Does this block obey the symmetry rule of if this block has been allocated?')
			call error_stop
		end if
		select type(val)
			type is (Tensor)
				totalData=val%getTotalDAta()
				if(size(p).ne.totalData)then
					call writemess('ERROR in setting value to the ith block of the tensor',-1)
					call writemess('size(block)='+size(p),-1)
					call writemess('size(input)='+totalData,-1)
					call error_stop
				end if
				call FastCopyArray(p,val%Data%classData,totalData)
			class default
				totalData=1
				if(ClassSize(p).ne.totalData)then
					call writemess('ERROR in setting value to the ith block of the tensor',-1)
					call writemess('size(block)='+size(p),-1)
					call writemess('size(input)='+totalData,-1)
					call error_stop
				end if
				call ClassPointer1DFunc(p,1,TMP0)
				call CopyData(TMP0,val)
		end select
		return
	end subroutine
	subroutine set_tensor_blockvec_value1(A,vec,val)!use for tensor only
		class(Tensor),intent(inout)::A
		integer,intent(in)::vec(:)
		class(*),intent(in)::val(:)
		class(*),pointer::p(:)
		integer,pointer::dim(:)
		integer::totalData
		if(.not.A%getSymmetryFlag())then
			call writemess('Can NOT setting the data the symmetry tensor, the tensor is not of symmetry',-1)
			call error_stop
		end if
		call A%pointDim(dim)
		call check_indices_and_dim(dim,vec)
		call A%ClassPointer(p,vec) 
		if(.not.associated(p))then
			call writemess('The block in the tensor is empty',-1)
			call writemess('Does this block obey the symmetry rule of if this block has been allocated?')
			call error_stop
		end if
		select type(val)
			type is (Tensor)
				call writemess('ERROR in setting value to the ith block of the tensor',-1)
				call error_stop
			class default
				totalData=Classsize(val)
				if(ClassSize(p).ne.totalData)then
					call writemess('ERROR in setting value to the ith block of the tensor',-1)
					call writemess('size(block)='+size(p),-1)
					call writemess('size(input)='+totalData,-1)
					call error_stop
				end if
				call FastCopyArray(p,val,totalData)
		end select
		return
	end subroutine

	subroutine check_indices_and_dim(dim,indices,TotalData)
		integer,intent(in)::dim(:),indices(:)
		integer,optional,intent(in)::TotalData
		integer::i
		if(size(dim).ne.size(indices))then
			call writemess('ERROR in input indices and dimenison',-1)
			call writemess('size(dim)='+size(dim),-1)
			call writemess('size(indices)='+size(indices),-1)
			call error_stop
		end if
		do i=1,size(dim)
			if(indices(i).le.0)then
				call writemess('ERROR in input indices',-1)
				call writemess('indices(i)='+indices(i),-1)
				call writemess('indices:',-1)
				call writemess(indices,-1)
				call error_stop
			end if
			if(indices(i).gt.dim(i))then
				call writemess('ERROR in input indices and dimenison',-1)
				call writemess('indices(i)='+indices(i),-1)
				call writemess('dim(i)='+dim(i),-1)
				call writemess('indices:',-1)
				call writemess(indices,-1)
				call writemess('dim:',-1)
				call writemess(dim,-1)
				call error_stop
			end if
		end do
		if(present(TotalData))then
			if(product(dim).ne.TotalData)then
				call writemess('ERROR in input indices',-1)
				call writemess('TotalData='+TotalData,-1)
				call writemess('product(dim)='+product(dim),-1)
				call writemess('dim:',-1)
				call writemess(dim,-1)
				call error_stop
			end if
		end if
		return
	end subroutine
	!-------------use both in tensor and symtensor------------------------------------

	subroutine set_All_value(A,val)
		class(Tensor),intent(inout)::A
		class(*),intent(in)::val
		class(*),pointer::TMP0
		integer::totalData
		select type(val)
			type is (Tensor)
				totalData=val%getTotalData()
				if(totalData.eq.A%getTotalDAta())then
					call FastCopyArray(A%Data%ClassData,val%Data%ClassData,totalData)
				else if(totalData.eq.1)then
					call ClassPointer1DFunc(val%Data%ClassData,1,TMP0)
					call ModifyAllValue(A%Data%ClassData,TMP0,A%getTotalData())
				else
					call writemess('ERROR in setting value to the ith element of the tensor',-1)
					call writemess('size(block)='+A%getTotalDAta(),-1)
					call writemess('size(input)='+totalData,-1)
					call error_stop
				end if
				
			class default
				call ModifyAllValue(A%Data%ClassData,val,A%getTotalDAta())
		end select
		return
	end subroutine
	subroutine set_All_value2(A,val)
		class(Tensor),intent(inout)::A
		class(*),intent(in)::val(:)
		integer::totalData
		totalData=A%getTotalData()
		if(ClassSize(val).ne.totalData)then
			call writemess('ERROR in setData',-1)
			call writemess('ClassSize(val)='+ClassSize(val),-1)
			call writemess('A%getTotalData()='+A%getTotalData(),-1)
			call error_stop
		end if
		call FastCopyArray(A%Data%ClassData,val,totalData)
		return
	end subroutine

	!set_ith_value
	!  if A is symmetry type, set the blocki
	!  if A is not symmetry type, set the elementi

	subroutine set_ith_value0(A,ith,val)
		class(Tensor),intent(inout)::A
		integer,intent(in)::ith
		class(*),intent(in)::val
		integer::totalData
		if(A%getSymmetryFlag())then
			call set_tensor_block_value0(A,ith,val)
		else
			call set_tensor_value0(A,ith,val)
		end if
		return
	end subroutine

	subroutine set_ith_value1(A,ith,val)
		class(Tensor),intent(inout)::A
		integer,intent(in)::ith
		class(*),intent(in)::val(:)
		integer::totalData
		if(A%getSymmetryFlag())then
			call set_tensor_block_value1(A,ith,val)
		else
			call writemess('ERROR in setting value for non-symmetry tensor',-1)
			call error_stop
		end if
		return
	end subroutine

	subroutine set_vec_value0(A,vec,val)
		class(Tensor),intent(inout)::A
		integer,intent(in)::vec(:)
		class(*),intent(in)::val
		if(A%getSymmetryFlag())then
			call set_tensor_blockvec_value0(A,vec,val)
		else
			call set_tensor_vec_value0(A,vec,val)
		end if
		return
	end subroutine

	subroutine set_vec_value1(A,vec,val)
		class(Tensor),intent(inout)::A
		integer,intent(in)::vec(:)
		class(*),intent(in)::val(:)
		if(A%getSymmetryFlag())then
			call set_tensor_blockvec_value1(A,vec,val)
		else
			call writemess('ERROR in setting value for non-symmetry tensor',-1)
			call error_stop
		end if
		return
	end subroutine

	subroutine set_some_value0(A,ith,val,jth)
		class(Tensor),intent(inout)::A
		type(Tensor),intent(in)::val
		integer,intent(in)::ith(2)
		integer,intent(in)::jth(2)
		class(*),pointer::TMP(:)
		call ClassPointer1DFunc(Val%Data%ClassData,jth(1),jth(2),TMP)
		call ModifySomeValue(A%Data%ClassData,TMP,ith(1),ith(2))
	end subroutine

	subroutine set_some_value1(A,ith,val,jth)!use for tensor only
		class(Tensor),intent(inout)::A
		class(*),intent(in)::val(:)
		integer,intent(in)::ith(2)
		integer,intent(in)::jth(2)
		class(*),pointer::TMP(:)
		select type(val)
			type is (Tensor)
				call writemess('ERROR in setting value to the ith element of the tensor',-1)
				call error_stop
			class default
				call ClassPointer1DFunc(Val,jth(1),jth(2),TMP)
				call ModifySomeValue(A%Data%ClassData,TMP,ith(1),ith(2))
		end select
		return
	end subroutine