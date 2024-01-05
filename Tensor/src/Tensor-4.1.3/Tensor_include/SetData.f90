	subroutine set_tensor_Some_Data(A,ith1,ith2,val,jth1,jth2)!use for tensor only
		class(Tensor),intent(inout)::A
		class(*),target,intent(in)::val(:)
		integer,intent(in)::ith1,ith2
		integer,optional,intent(in)::jth1,jth2
		class(*),pointer::TMP(:)
		if(present(jth1).and.present(jth2))then
			call ClassPointer1DFunc(Val,jth1,jth2,TMP)
		else
			TMP=>val
		end if
		call ModifySomeValue(A%Data%ClassData,TMP,ith1,ith2)
		return
	end subroutine

	subroutine set_tensor_All_data(A,val)!use for tensor only
		class(Tensor),intent(inout)::A
		class(*),intent(in)::val
		call ModifyAllValue(A%Data%ClassData,val,A%getTotalData())
		return
	end subroutine
	subroutine set_tensor_All_data2(A,val)!use for tensor only
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

	subroutine set_tensor_data(A,ith,val)!use for tensor only
		class(Tensor),intent(inout)::A
		integer,intent(in)::ith
		class(*),intent(in)::val
		call ModifyaValue(A%Data%ClassData,val,ith)
		return
	end subroutine
	subroutine set_tensor_vec_data(A,vec,val)!use for tensor only
		class(Tensor),intent(inout)::A
		integer,intent(in)::vec(:)
		class(*),intent(in)::val
		class(*),pointer::p
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
		call CopyData(p,val)
		return
	end subroutine

	subroutine set_tensor_block_data(A,ith,val)!use for symtensor only
		class(Tensor),intent(inout)::A
		integer,intent(in)::ith
		class(*),intent(in)::val(:)
		class(*),pointer::p(:)
		integer::totalData
		if(.not.A%getSymmetryFlag())then
			call writemess('Can NOT setting the data the symmetry tensor, the tensor is not of symmetry',-1)
			call error_stop
		end if
		call A%Classpointer(p,ith)
		if(.not.associated(p))then
			call writemess('The block in the tensor is empty',-1)
			call writemess('Does this block obey the symmetry rule of if this block has been allocated?')
			call error_stop
		end if
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
		return
	end subroutine
	subroutine set_tensor_blockvec_data(A,vec,val)!use for tensor only
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
		totalData=ClassSize(val)
		if(ClassSize(p).ne.totalData)then
			call writemess('ERROR in setting value to the ith block of the tensor',-1)
			call writemess('size(block)='+size(p),-1)
			call writemess('size(input)='+totalData,-1)
			call error_stop
		end if
		call FastCopyArray(p,val,totalData)
		return
	end subroutine

	subroutine set_tensor_blockvec_data2(A,vec,val)!use for tensor only
		class(Tensor),intent(inout)::A
		integer,intent(in)::vec(:)
		class(*),intent(in)::val(:,:)
		class(*),pointer::p(:,:)
		integer,pointer::dim(:)
		integer::dim1,dim2
		if(.not.A%getSymmetryFlag())then
			call writemess('Can NOT setting the data the symmetry tensor, the tensor is not of symmetry',-1)
			call error_stop
		end if
		if(size(vec).ne.2)then
			call writemess('ERROR in set_tensor_blockvec_data2',-1)
			call writemess('size(vec)='+size(vec),-1)
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
		dim1=ClassSize(val,1)
		dim2=ClassSize(val,2)
		if(ClassSize(p,1).ne.dim1)then
			call writemess('ERROR in setting value to the ith block of the tensor',-1)
			call writemess('ClassSize(p,1)='+ClassSize(p,1),-1)
			call writemess('ClassSize(val,1)='+ClassSize(val,1),-1)
			call error_stop
		end if
		if(ClassSize(p,2).ne.dim2)then
			call writemess('ERROR in setting value to the ith block of the tensor',-1)
			call writemess('ClassSize(p,2)='+ClassSize(p,2),-1)
			call writemess('ClassSize(val,2)='+ClassSize(val,2),-1)
			call error_stop
		end if
		call CopyArray(p,val,dim1,dim2)
		return
	end subroutine
	subroutine set_tensor_blockvec_data3(A,vec,val)!use for tensor only
		class(Tensor),intent(inout)::A
		integer,intent(in)::vec(:)
		class(*),intent(in)::val(:,:,:)
		class(*),pointer::p(:,:,:)
		integer,pointer::dim(:)
		integer::dim1,dim2,dim3
		if(.not.A%getSymmetryFlag())then
			call writemess('Can NOT setting the data the symmetry tensor, the tensor is not of symmetry',-1)
			call error_stop
		end if
		if(size(vec).ne.3)then
			call writemess('ERROR in set_tensor_blockvec_data2',-1)
			call writemess('size(vec)='+size(vec),-1)
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
		dim1=ClassSize(val,1)
		dim2=ClassSize(val,2)
		dim3=ClassSize(val,3)
		if(ClassSize(p,1).ne.dim1)then
			call writemess('ERROR in setting value to the ith block of the tensor',-1)
			call writemess('ClassSize(p,1)='+ClassSize(p,1),-1)
			call writemess('ClassSize(val,1)='+ClassSize(val,1),-1)
			call error_stop
		end if
		if(ClassSize(p,2).ne.dim2)then
			call writemess('ERROR in setting value to the ith block of the tensor',-1)
			call writemess('ClassSize(p,2)='+ClassSize(p,2),-1)
			call writemess('ClassSize(val,2)='+ClassSize(val,2),-1)
			call error_stop
		end if
		if(ClassSize(p,3).ne.dim3)then
			call writemess('ERROR in setting value to the ith block of the tensor',-1)
			call writemess('ClassSize(p,3)='+ClassSize(p,3),-1)
			call writemess('ClassSize(val,3)='+ClassSize(val,3),-1)
			call error_stop
		end if
		call CopyArray(p,val,dim1,dim2,dim3)
		return
	end subroutine

