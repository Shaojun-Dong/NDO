	subroutine DATAName2T0(A,B)
		type(Tensor),intent(inout)::A
		DATATYPE,intent(in)::B
		A%Data=B
		A%Dimension=[1]
		return
	end subroutine
	subroutine DATAName2T1(A,B)
		type(Tensor),intent(inout)::A
		DATATYPE,intent(in)::B(:)
		A%Data=B
		A%Dimension=[size(B)]
		return
	end subroutine

	subroutine DATAName2T2(A,B)
		type(Tensor),intent(inout)::A
		DATATYPE,intent(in)::B(:,:)
		integer::m,n
		m=size(B,1)
		n=size(B,2)
		call A%Data%allocateData(size(B),DATATYPENumber)
		call ReshapCopyDATAName(A%Data%ClassData,B,size(B))
		A%Dimension=[m,n]
		return
	end subroutine

	subroutine DATAName2T3(A,B)
		type(Tensor),intent(inout)::A
		DATATYPE,target,intent(in)::B(:,:,:)
		integer::m,n,o
		m=size(B,1)
		n=size(B,2)
		o=size(B,3)
		call A%Data%allocateData(size(B),DATATYPENumber)
		call ReshapCopyDATAName(A%Data%ClassData,B,size(B))
		A%Dimension=[m,n,o]
		return
	end subroutine

	subroutine ReshapCopyDATAName(A,B,length)
		integer,intent(in)::length
		class(*)::A(:)
		DATATYPE::B(length)
		call FastCopyArray(A,B,length)
		return
	end subroutine

	subroutine DATAName2T4(A,B)
		type(Tensor),intent(inout)::A
		DATATYPE,intent(in)::B(:,:,:,:)
		integer::m,n,o,p
		m=size(B,1)
		n=size(B,2)
		o=size(B,3)
		p=size(B,4)
		call A%Data%allocateData(size(B),DATATYPENumber)
		call ReshapCopyDATAName(A%Data%ClassData,B,size(B))
		A%Dimension=[m,n,o,p]
		return
	end subroutine


	subroutine T2DATAName0(A,B)
		DATATYPE,intent(inout)::A
		type(Tensor),intent(in)::B
		class(*),pointer::p0
		if(.not.B%getFlag())then
			call writemess('The tensor is empty in assignment',-1)
			call error_stop
		end if
		if(.not.B%Data%getFlag())then
			if(B%getSymmetryFlag())then
				A=0
				return
			end if
			call writemess('ERROR in assignment for val=Tensor, totalData',-1)
			call writemess('B%getTotalData()='+B%getTotalData())
			call error_stop
		end if
		if(B%getTotalData().ne.1)then
			call writemess('ERROR in assignment for val=Tensor, totalData',-1)
			call writemess('B%getTotalData()='+B%getTotalData())
			call error_stop
		end if
		call ClassPointer1DFunc(B%Data%ClassData,1,p0)
		call copyData(A,p0)
		return
	end subroutine

	subroutine T2DATAName1(A,B)
		DATATYPE,intent(inout)::A(:)
		type(Tensor),intent(in)::B
		if(.not.B%getFlag())then
			call writemess('The tensor is empty in assignment',-1)
			call error_stop
		end if
		if(.not.B%Data%getFlag())then
			call writemess('ERROR in assignment for val=Tensor, totalData',-1)
			call writemess('B%getTotalData()='+B%getTotalData())
			call error_stop
		end if
		A=B%Data
		return
	end subroutine

	subroutine T2DATAName2(A,B)
		DATATYPE,intent(inout)::A(:,:)
		type(Tensor),intent(in)::B
		integer::m,n
		integer,pointer::dim(:)
		if(.not.B%getFlag())then
			call writemess('The tensor is empty in assignment',-1)
			call error_stop
		end if
		if(.not.B%Data%getFlag())then
			call writemess('ERROR in assignment for val=Tensor, totalData',-1)
			call writemess('B%getTotalData()='+B%getTotalData())
			call error_stop
		end if
		m=size(A,1)
		n=size(A,2)
		if(B%getRank().ne.2)then
			call writemess('ERROR in assignment for DATATYPE',-1)
			call writemess('B%getRank()='+B%getRank())
			call error_Stop
		end if
		call B%pointDim(dim)
		if(dim(1).ne.m)then
			call writemess('ERROR in assignment for DATATYPE',-1)
			call writemess('size(A,1)='+size(A,1),-1)
			call writemess('dim(1)='+dim(1),-1)
			call error_Stop
		end if
		if(dim(2).ne.n)then
			call writemess('ERROR in assignment for DATATYPE',-1)
			call writemess('size(A,2)='+size(A,2),-1)
			call writemess('dim(2)='+dim(2),-1)
			call error_Stop
		end if
		call ReshapCopy2DATAName(A,B%Data%ClassData,B%getTotalData())
		return
	end subroutine

	subroutine ReshapCopy2DATAName(A,B,length)
		integer,intent(in)::length
		DATATYPE::A(length)
		class(*)::B(:)
		call FastCopyArray(A,B,length)
		return
	end subroutine

	subroutine T2DATAName3(A,B)
		DATATYPE,intent(inout)::A(:,:,:)
		type(Tensor),intent(in)::B
		integer::m,n,o
		integer,pointer::dim(:)
		if(.not.B%getFlag())then
			call writemess('The tensor is empty in assignment',-1)
			call error_stop
		end if
		if(.not.B%Data%getFlag())then
			call writemess('ERROR in assignment for val=Tensor, totalData',-1)
			call writemess('B%getTotalData()='+B%getTotalData())
			call error_stop
		end if
		m=size(A,1)
		n=size(A,2)
		o=size(A,3)
		if(B%getRank().ne.3)then
			call writemess('ERROR in assignment for DATATYPE',-1)
			call error_Stop
		end if
		call B%pointDim(dim)
		if(dim(1).ne.m)then
			call writemess('ERROR in assignment for DATATYPE',-1)
			call writemess('size(A,1)='+size(A,1),-1)
			call writemess('dim(1)='+dim(1),-1)
			call error_Stop
		end if
		if(dim(2).ne.n)then
			call writemess('ERROR in assignment for DATATYPE',-1)
			call writemess('size(A,2)='+size(A,2),-1)
			call writemess('dim(2)='+dim(2),-1)
			call error_Stop
		end if
		if(dim(3).ne.o)then
			call writemess('ERROR in assignment for DATATYPE',-1)
			call writemess('size(A,3)='+size(A,3),-1)
			call writemess('dim(3)='+dim(3),-1)
			call error_Stop
		end if
		call ReshapCopy2DATAName(A,B%Data%ClassData,B%getTotalData())
		return
	end subroutine

	subroutine T2DATAName4(A,B)
		DATATYPE,intent(inout)::A(:,:,:,:)
		type(Tensor),intent(in)::B
		integer::m,n,o,p
		integer,pointer::dim(:)
		if(.not.B%getFlag())then
			call writemess('The tensor is empty in assignment',-1)
			call error_stop
		end if
		if(.not.B%Data%getFlag())then
			call writemess('ERROR in assignment for val=Tensor, totalData',-1)
			call writemess('B%getTotalData()='+B%getTotalData())
			call error_stop
		end if
		m=size(A,1)
		n=size(A,2)
		o=size(A,3)
		p=size(A,4)
		if(B%getRank().ne.4)then
			call writemess('ERROR in assignment for DATATYPE',-1)
			call error_Stop
		end if
		call B%pointDim(dim)
		if(dim(1).ne.m)then
			call writemess('ERROR in assignment for DATATYPE',-1)
			call writemess('size(A,1)='+size(A,1),-1)
			call writemess('dim(1)='+dim(1),-1)
			call error_Stop
		end if
		if(dim(2).ne.n)then
			call writemess('ERROR in assignment for DATATYPE',-1)
			call writemess('size(A,2)='+size(A,2),-1)
			call writemess('dim(2)='+dim(2),-1)
			call error_Stop
		end if
		if(dim(3).ne.o)then
			call writemess('ERROR in assignment for DATATYPE',-1)
			call writemess('size(A,3)='+size(A,3),-1)
			call writemess('dim(3)='+dim(3),-1)
			call error_Stop
		end if
		if(dim(4).ne.p)then
			call writemess('ERROR in assignment for DATATYPE',-1)
			call writemess('size(A,4)='+size(A,4),-1)
			call writemess('dim(4)='+dim(4),-1)
			call error_Stop
		end if
		call ReshapCopy2DATAName(A,B%Data%ClassData,B%getTotalData())
		return
	end subroutine
