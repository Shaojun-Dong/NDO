



	subroutine eyeQN2(A,QN,leftQN)
		class(Tensor),intent(inout)::A
		type(QuanNum),intent(in)::QN
		logical,intent(in)::leftQN
		type(DataArray)::AData
		integer::iType,nonSymDimi,i,ii,jj,ith
		class(*),pointer::mold,Rp(:,:),Ap(:),pin,pout
		integer,pointer::ARROW(:),dim(:),deg(:)
		real*4,pointer::sp(:)
		if(.not.A%Data%getFlag())then
			call writemess('ERROR in eye, there is no data in the input tensor',-1)
			call error_stop
		end if
		if(A%getDimFlag())then
			if(A%getRank().ne.1)then
				call writemess('ERROR in eye, the input  tensor should be a vector',-1)
				call error_stop
			end if
		end if
		AData=A%Data
		call A%Empty()
		A%Dimension=[QN,QN]
		if((.not.LeftQN))then
			ith=1
		else
			ith=2
		end if
		if(A%getFermiFlag())then
			call A%pointArrow(arrow)
			arrow(ith)=-arrow(ith)
		end if
		call hermitian_conjugate_dimension(A%Dimension,ith)
		call A%pointDim(dim)
		call A%pointdeg(deg,1)
		nonSymDimi=sum(deg)
		iType=AData%getType()
		if((iType.ge.0).and.(iType.le.8))then
			call A%Data%allocateDataMomery(product(dim),nonSymDimi*nonSymDimi,AData%getType())
		else
			call ClassPointer1DFunc(AData%ClassData,1,mold)
			call A%Data%allocateDataMomeryClassType(product(dim),nonSymDimi*nonSymDimi,mold,AData%getType())
		end if
		do i=1,dim(1)
			if(AData%getFlag(i))then
				call A%setBlockMomery([i,i])
				call A%ClassPointer(Rp,[i,i])
				call AData%ClassPointer(Ap,i)
				do jj=1,deg(i)
					do ii=1,deg(i)
						call ClassPointer2DFunc(Rp,ii,jj,pout)
						if(ii.eq.jj)then
							call ClassPointer1DFunc(Ap,ii,pin)
							call copyData(pout,pin)
						else
							call copyData(pout,0d0)
						end if
					end do
				end do
			end if
		end do
		return
	end subroutine

	subroutine eyeQN3(A,QN,classType,leftQN)
		class(Tensor),intent(inout)::A
		type(QuanNum),intent(in)::QN
		integer,intent(in)::classType
		logical,intent(in)::leftQN
		integer::i,ii,jj,ith
		class(*),pointer::mold,Rp(:,:),Ap(:),pin,pout
		integer,pointer::ARROW(:),dim(:),deg(:)
		type(Dimension)::dimen
		if(A%Data%getFlag())then
			call writemess('ERROR in eye, there IS data in the input tensor',-1)
			call error_stop
		end if
		dimen=[QN,QN]
		if((.not.LeftQN))then
			ith=1
		else
			ith=2
		end if
		if(dimen%getFermiFlag())then
			call dimen%pointArrow(arrow)
			arrow(ith)=-arrow(ith)
		end if
		call hermitian_conjugate_dimension(dimen,ith)
		call A%allocate(dimen,ClassType)
		call A%pointDim(dim)
		call A%pointdeg(deg,1)
		do i=1,dim(1)
			if(A%getFlag([i,i]))then
				call A%ClassPointer(Rp,[i,i])
				do jj=1,deg(i)
					do ii=1,deg(i)
						call ClassPointer2DFunc(Rp,ii,jj,pout)
						if(ii.eq.jj)then
							call copyData(pout,1d0)
						else
							call copyData(pout,0d0)
						end if
					end do
				end do
			end if
		end do
		return
	end subroutine
	subroutine eyeQN4(A,QN,classType,leftQN)
		class(Tensor),intent(inout)::A
		type(QuanNum),intent(in)::QN
		character(len=*),intent(in)::classType
		logical,intent(in)::leftQN
		integer::iclassType
		iclassType=select_data_type_char(classType)
		call eyeQN3(A,QN,iclassType,leftQN)
		return
	end subroutine
	subroutine eyeQN5(A,QN,source,classType,leftQN)
		class(Tensor),intent(inout)::A
		type(QuanNum),intent(in)::QN
		class(*),intent(in)::source
		integer,intent(in)::classType
		logical,intent(in)::leftQN
		integer::i,ii,jj,ith
		class(*),pointer::mold,Rp(:,:),Ap(:),pin,pout
		integer,pointer::ARROW(:),dim(:),deg(:)
		real*4,pointer::sp(:)
		type(Dimension)::dimen
		if(A%Data%getFlag())then
			call writemess('ERROR in eye, there IS data in the input tensor',-1)
			call error_stop
		end if
		dimen=[QN,QN]
		if((.not.LeftQN))then
			ith=1
		else
			ith=2
		end if
		if(A%getFermiFlag())then
			call A%pointArrow(arrow)
			arrow(ith)=-arrow(ith)
		end if
		call hermitian_conjugate_dimension(dimen,ith)
		call A%allocateClassType(dimen,source,ClassType)
		call A%pointDim(dim)
		call A%pointdeg(deg,1)
		do i=1,dim(1)
			if(A%getFlag([i,i]))then
				call A%ClassPointer(Rp,[i,i])
				do jj=1,deg(i)
					do ii=1,deg(i)
						call ClassPointer2DFunc(Rp,ii,jj,pout)
						if(ii.eq.jj)then
							call copyData(pout,1d0)
						else
							call copyData(pout,0d0)
						end if
					end do
				end do
			end if
		end do
		return
	end subroutine
	function eyeQNTensor1(QN,classType,leftQN)result(Res)
		type(Tensor)::Res
		type(QuanNum),intent(in)::QN
		integer,intent(in)::classType
		logical,intent(in)::leftQN
		call eyeQN3(Res,QN,classType,leftQN)
		return
	end function
	function eyeQNTensor2(QN,classType,leftQN)result(Res)
		type(Tensor)::Res
		type(QuanNum),intent(in)::QN
		character(len=*),intent(in)::classType
		logical,intent(in)::leftQN
		call eyeQN4(Res,QN,classType,leftQN)
		return
	end function
	function eyeQNTensor3(QN,source,classType,leftQN)result(Res)
		type(Tensor)::Res
		type(QuanNum),intent(in)::QN
		class(*),intent(in)::source
		integer,intent(in)::classType
		logical,intent(in)::leftQN
		call eyeQN5(Res,QN,source,classType,leftQN)
		return
	end function

	!***************************************

	subroutine eyeQN6(A,QN,leftQN)
		class(Tensor),intent(inout)::A
		type(QuanNum),intent(in)::QN
		character(len=*),intent(in)::leftQN
		if(leftQN.equ.'left')then
			call eyeQN2(A,QN,.true.)
		else if(leftQN.equ.'leftQN')then
			call eyeQN2(A,QN,.true.)
		else if(leftQN.equ.'right')then
			call eyeQN2(A,QN,.false.)
		else  if(leftQN.equ.'rightQN')then
			call eyeQN2(A,QN,.false.)
		else
			call writemess('ERROR value leftQN='+leftQN,-1)
			call writemess('leftQN should be one of:',-1)
			call writemess(' left',-1)
			call writemess(' leftQN',-1)
			call writemess(' right',-1)
			call writemess(' rightQN',-1)
			call error_stop
		end if
		return
	end subroutine

	subroutine eyeQN7(A,QN,classType,leftQN)
		class(Tensor),intent(inout)::A
		type(QuanNum),intent(in)::QN
		integer,intent(in)::classType
		character(len=*),intent(in)::leftQN
		if(leftQN.equ.'left')then
			call eyeQN3(A,QN,classType,.true.)
		else if(leftQN.equ.'leftQN')then
			call eyeQN3(A,QN,classType,.true.)
		else if(leftQN.equ.'right')then
			call eyeQN3(A,QN,classType,.false.)
		else  if(leftQN.equ.'rightQN')then
			call eyeQN3(A,QN,classType,.false.)
		else
			call writemess('ERROR value leftQN='+leftQN,-1)
			call writemess('leftQN should be one of:',-1)
			call writemess(' left',-1)
			call writemess(' leftQN',-1)
			call writemess(' right',-1)
			call writemess(' rightQN',-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine eyeQN8(A,QN,classType,leftQN)
		class(Tensor),intent(inout)::A
		type(QuanNum),intent(in)::QN
		character(len=*),intent(in)::classType
		character(len=*),intent(in)::leftQN
		if(leftQN.equ.'left')then
			call eyeQN4(A,QN,classType,.true.)
		else if(leftQN.equ.'leftQN')then
			call eyeQN4(A,QN,classType,.true.)
		else if(leftQN.equ.'right')then
			call eyeQN4(A,QN,classType,.false.)
		else  if(leftQN.equ.'rightQN')then
			call eyeQN4(A,QN,classType,.false.)
		else
			call writemess('ERROR value leftQN='+leftQN,-1)
			call writemess('leftQN should be one of:',-1)
			call writemess(' left',-1)
			call writemess(' leftQN',-1)
			call writemess(' right',-1)
			call writemess(' rightQN',-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine eyeQN9(A,QN,source,classType,leftQN)
		class(Tensor),intent(inout)::A
		type(QuanNum),intent(in)::QN
		class(*),intent(in)::source
		integer,intent(in)::classType
		character(len=*),intent(in)::leftQN
		if(leftQN.equ.'left')then
			call eyeQN5(A,QN,source,classType,.true.)
		else if(leftQN.equ.'leftQN')then
			call eyeQN5(A,QN,source,classType,.true.)
		else if(leftQN.equ.'right')then
			call eyeQN5(A,QN,source,classType,.false.)
		else  if(leftQN.equ.'rightQN')then
			call eyeQN5(A,QN,source,classType,.false.)
		else
			call writemess('ERROR value leftQN='+leftQN,-1)
			call writemess('leftQN should be one of:',-1)
			call writemess(' left',-1)
			call writemess(' leftQN',-1)
			call writemess(' right',-1)
			call writemess(' rightQN',-1)
			call error_stop
		end if
		return
	end subroutine
	function eyeQNTensor4(QN,classType,leftQN)result(Res)
		type(Tensor)::Res
		type(QuanNum),intent(in)::QN
		integer,intent(in)::classType
		character(len=*),intent(in)::leftQN
		call eyeQN7(Res,QN,classType,leftQN)
		return
	end function
	function eyeQNTensor5(QN,classType,leftQN)result(Res)
		type(Tensor)::Res
		type(QuanNum),intent(in)::QN
		character(len=*),intent(in)::classType
		character(len=*),intent(in)::leftQN
		call eyeQN8(Res,QN,classType,leftQN)
		return
	end function
	function eyeQNTensor6(QN,source,classType,leftQN)result(Res)
		type(Tensor)::Res
		type(QuanNum),intent(in)::QN
		class(*),intent(in)::source
		integer,intent(in)::classType
		character(len=*),intent(in)::leftQN
		call eyeQN9(Res,QN,source,classType,leftQN)
		return
	end function




	!


	subroutine eyeMN1(A,M,N)
		class(Tensor),intent(inout)::A
		integer,intent(in)::M,N
		class(*),pointer::Rp(:,:),Ap(:),pout,pin
		integer::i,j
		type(DataArray)::AData
		if(.not.A%Data%getFlag())then
			call writemess('ERROR in eye, there is no data in the input tensor',-1)
			call error_stop
		end if
		if(A%getDimFlag())then
			if(A%getRank().ne.1)then
				call writemess('ERROR in eye, the input  tensor should be a vector',-1)
				call error_stop
			end if
		end if
		AData=A%Data
		call A%allocate([M,N],A%getType())
		call A%ClassPointer(Rp)
		call AData%ClassPointer(Ap)
		do j=1,N
			do i=1,M
				call ClassPointer2DFunc(Rp,i,j,pout)
				if(i.eq.j)then
					call ClassPointer1DFunc(Ap,i,pin)
					call CopyData(pout,pin)
				else
					call CopyData(pout,0d0)
				end if
			end do
		end do
		return
	end subroutine
	subroutine eyeMN2(A,M,N,classType)
		class(Tensor),intent(inout)::A
		integer,intent(in)::M,N
		integer,intent(in)::classType
		class(*),pointer::Rp(:,:),pout
		integer::i,j
		if(A%Data%getFlag())then
			call writemess('ERROR in eye, there IS data in the input tensor',-1)
			call error_stop
		end if
		call A%allocate([M,N],classType)
		call A%ClassPointer(Rp)
		do j=1,N
			do i=1,M
				call ClassPointer2DFunc(Rp,i,j,pout)
				if(i.eq.j)then
					call CopyData(pout,1d0)
				else
					call CopyData(pout,0d0)
				end if
			end do
		end do
		return
	end subroutine
	subroutine eyeMN3(A,M,N,classType)
		class(Tensor),intent(inout)::A
		integer,intent(in)::M,N
		character(len=*),intent(in)::classType
		integer::iclassType
		iclassType=select_data_type_char(classType)
		call eyeMN2(A,M,N,iclassType)
		return
	end subroutine
	subroutine eyeMN4(A,M,N,source,classType)
		class(Tensor),intent(inout)::A
		integer,intent(in)::M,N
		class(*),intent(in)::source
		integer,intent(in)::classType
		class(*),pointer::Rp(:,:),pout
		integer::i,j
		if(A%Data%getFlag())then
			call writemess('ERROR in eye, there IS data in the input tensor',-1)
			call error_stop
		end if
		call A%allocateClassType([M,N],source,classType)
		call A%ClassPointer(Rp)
		do j=1,N
			do i=1,M
				call ClassPointer2DFunc(Rp,i,j,pout)
				if(i.eq.j)then
					call CopyData(pout,1d0)
				else
					call CopyData(pout,0d0)
				end if
			end do
		end do
		return
	end subroutine
	function eyeMNTensor1(M,N,classType)result(Res)
		type(Tensor)::Res
		integer,intent(in)::M,N
		integer,intent(in)::classType
		call eyeMN2(Res,M,N,classType)
	end function
	function eyeMNTensor2(M,N,classType)result(Res)
		type(Tensor)::Res
		integer,intent(in)::M,N
		character(len=*),intent(in)::classType
		call eyeMN3(Res,M,N,classType)
	end function
	function eyeMNTensor3(M,N,source,classType)result(Res)
		type(Tensor)::Res
		integer,intent(in)::M,N
		class(*),intent(in)::source
		integer,intent(in)::classType
		call eyeMN4(Res,M,N,source,classType)
	end function


	subroutine eye0(A)
		class(Tensor),intent(inout)::A
		class(*),pointer::Ap(:,:),pout
		integer::i,ii
		if(.not.A%getFlag())then
			call writemess('Allocate Tensor first',-1)
			call error_stop
		end if
		if(.not.A%DAta%getFlag())then
			call writemess('Allocate Tensor first',-1)
			call error_stop
		end if
		if(A%getRank().ne.2)then
			call writemess('error in call A%eye(), A should be a matrix',-1)
			call error_stop
		end if
		if(A%dim(1).ne.A%dim(2))then
			call writemess('error in call A%eye(), A should be a matrix with the same dimension',-1)
			call error_stop
		end if
		call A%zero()
		if(A%getSymmetryFlag())then
			do i=1,A%dim(1)
				call A%Classpointer(Ap,[i,i])
				if(.not.associated(Ap))then
					call writemess('The block in the tensor is empty',-1)
					call writemess('Does this block obey the symmetry rule of if this block has been allocated?')
					call error_stop
				end if
				do ii=1,ClassSize(Ap,1)
					call ClassPointer2DFunc(Ap,ii,ii,pout)
					call copyData(pout,1)
				end do
			end do
		else
			call A%Classpointer(Ap)
			if(.not.associated(Ap))then
				call writemess('The block in the tensor is empty',-1)
				call writemess('Does this block obey the symmetry rule of if this block has been allocated?')
				call error_stop
			end if
			do i=1,A%dim(1)
				call ClassPointer2DFunc(Ap,i,i,pout)
				call copyData(pout,1)
			end do
		end if
		return
	end subroutine

