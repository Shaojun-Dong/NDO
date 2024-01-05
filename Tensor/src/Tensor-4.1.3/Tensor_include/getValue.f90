#define getValue0FUNCNAME getValue0i
#define getValueSymFUNCNAME getValueSymi
#define getValueVecFUNCNAME getValueVeci
#define getValueSymVecFUNCNAME getValueSymVeci
#define getValueAllFUNCNAME getValueAlli
#define DATATYPE integer
#define DATATYPE2 integer
#include "templet/getValue0.f90"
#undef getValue0FUNCNAME
#undef getValueSymFUNCNAME
#undef getValueVecFUNCNAME
#undef getValueSymVecFUNCNAME
#undef getValueAllFUNCNAME
#undef DATATYPE
#undef DATATYPE2

#define getValue0FUNCNAME getValue0s
#define getValueSymFUNCNAME getValueSyms
#define getValueVecFUNCNAME getValueVecs
#define getValueSymVecFUNCNAME getValueSymVecs
#define getValueAllFUNCNAME getValueAlls
#define DATATYPE real*4
#define DATATYPE2 real*4
#include "templet/getValue0.f90"
#undef getValue0FUNCNAME
#undef getValueSymFUNCNAME
#undef getValueVecFUNCNAME
#undef getValueSymVecFUNCNAME
#undef getValueAllFUNCNAME
#undef DATATYPE
#undef DATATYPE2

#define getValue0FUNCNAME getValue0d
#define getValueSymFUNCNAME getValueSymd
#define getValueVecFUNCNAME getValueVecd
#define getValueSymVecFUNCNAME getValueSymVecd
#define getValueAllFUNCNAME getValueAlld
#define DATATYPE real*8
#define DATATYPE2 real*8
#include "templet/getValue0.f90"
#undef getValue0FUNCNAME
#undef getValueSymFUNCNAME
#undef getValueVecFUNCNAME
#undef getValueSymVecFUNCNAME
#undef getValueAllFUNCNAME
#undef DATATYPE
#undef DATATYPE2

#define getValue0FUNCNAME getValue0c
#define getValueSymFUNCNAME getValueSymc
#define getValueVecFUNCNAME getValueVecc
#define getValueSymVecFUNCNAME getValueSymVecc
#define getValueAllFUNCNAME getValueAllc
#define DATATYPE complex*8
#define DATATYPE2 complex*8
#include "templet/getValue0.f90"
#undef getValue0FUNCNAME
#undef getValueSymFUNCNAME
#undef getValueVecFUNCNAME
#undef getValueSymVecFUNCNAME
#undef getValueAllFUNCNAME
#undef DATATYPE
#undef DATATYPE2

#define getValue0FUNCNAME getValue0z
#define getValueSymFUNCNAME getValueSymz
#define getValueVecFUNCNAME getValueVecz
#define getValueSymVecFUNCNAME getValueSymVecz
#define getValueAllFUNCNAME getValueAllz
#define DATATYPE complex*16
#define DATATYPE2 complex*16
#include "templet/getValue0.f90"
#undef getValue0FUNCNAME
#undef getValueSymFUNCNAME
#undef getValueVecFUNCNAME
#undef getValueSymVecFUNCNAME
#undef getValueAllFUNCNAME
#undef DATATYPE
#undef DATATYPE2

#define getValue0FUNCNAME getValue0l
#define getValueSymFUNCNAME getValueSyml
#define getValueVecFUNCNAME getValueVecl
#define getValueSymVecFUNCNAME getValueSymVecl
#define getValueAllFUNCNAME getValueAlll
#define DATATYPE logical
#define DATATYPE2 logical
#include "templet/getValue0.f90"
#undef getValue0FUNCNAME
#undef getValueSymFUNCNAME
#undef getValueVecFUNCNAME
#undef getValueSymVecFUNCNAME
#undef getValueAllFUNCNAME
#undef DATATYPE
#undef DATATYPE2

#define getValue0FUNCNAME getValue0a
#define getValueSymFUNCNAME getValueSyma
#define getValueVecFUNCNAME getValueVeca
#define getValueSymVecFUNCNAME getValueSymVeca
#define getValueAllFUNCNAME getValueAlla
#define DATATYPE character(len=A%data%DataCharacterLen)
#define DATATYPE2 character(len=A%data%DataCharacterLen)
#include "templet/getValue0.f90"
#undef getValue0FUNCNAME
#undef getValueSymFUNCNAME
#undef getValueVecFUNCNAME
#undef getValueSymVecFUNCNAME
#undef getValueAllFUNCNAME
#undef DATATYPE
#undef DATATYPE2


	function getValue0T(A,i)Result(Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
		integer,intent(in)::i
		class(*),pointer::p0
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
		call Res%allocate([1],A%getType())
		call ClassPointer1DFunc(A%Data%ClassData,i,p0)
		call ModifyaValue(Res%Data%ClassData,p0,1)
		return
	end function
	function getValueSymT(A,blocki,i)Result(Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
		integer,intent(in)::blocki,i
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

		call Res%allocate([1],A%getType())
		call ClassPointer1DFunc(p,i,p0)
		call ModifyaValue(Res%Data%ClassData,p0,1)
		return
	end function

	function getValueVecT(A,vec)Result(Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
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
		call Res%allocate([1],A%getType())
		call ModifyaValue(Res%Data%ClassData,p0,1)
		return
	end function

	function getValueSymVecT(A,Blockvec,vec)Result(Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
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
				call ClassPointer1DFunc(p,vec(1),p0)
			case(2)
				call A%ClassPointer(p2,Blockvec)
				if(.not.associated(p2))then
					call writemess('ERROR in getValueSymVecFUNCNAME',-1)
					call error_stop
				end if
				call ClassPointer2DFunc(p2,vec(1),vec(2),p0)
			case(3)
				call A%ClassPointer(p3,Blockvec)
				if(.not.associated(p3))then
					call writemess('ERROR in getValueSymVecFUNCNAME',-1)
					call error_stop
				end if
				call ClassPointer3DFunc(p3,vec(1),vec(2),vec(3),p0)
			case(4)
				call A%ClassPointer(p4,Blockvec)
				if(.not.associated(p4))then
					call writemess('ERROR in getValueSymVecFUNCNAME',-1)
					call error_stop
				end if
				call ClassPointer4DFunc(p4,vec(1),vec(2),vec(3),vec(4),p0)
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
				call ClassPointer1DFunc(p,index,p0)
		end select

		call Res%allocate([1],A%getType())
		call ModifyaValue(Res%Data%ClassData,p0,1)
		return
	end function

	function getValueAllT(A)Result(Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
		integer::TotalData,i
		if(.not.A%Data%getFlag())then
			call writemess('ERROR, the tensor is empty',-1)
			call error_stop
		end if
		TotalData=A%getTotalData()
		call Res%allocate([TotalData],A%getType())
		call FastCopyArray(Res%Data%ClassData,A%Data%ClassData,TotalData)
		return
	end function

	function getithTensor(A,ith)Result(Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
		integer,intent(in)::ith
		class(*),pointer::p(:),p0
		integer,pointer::Blockdim(:),index(:)
		integer::rank
		if(.not.A%Data%getFlag())then
			call Res%empty()
			return
		end if
		rank=A%getRank()
		if(A%getSymmetryFlag())then
			call A%ClassPointer(p,ith)
			if(associated(p))then
				call WorkingMemory%check()
				call WorkingMemory%allocate(1,rank+rank)
				call WorkingMemory%get_memory(Blockdim,rank)
				call WorkingMemory%get_memory(index,rank)
				call IndesToaddress(A%dim(),index,ith)
				Blockdim=A%getBlockDim(index)
				call Res%allocate(Blockdim,A%getType())
				call WorkingMemory%free()
				call FastCopyArray(Res%Data%ClassData,p,Res%getTotalData())
			else
				call Res%empty()
			end if 
		else
			call Res%allocate([1],A%getType())
			call ClassPointer1DFunc(A%Data%ClassData,ith,p0)
			call ModifyaValue(Res%Data%ClassData,p0,1)
		end if
		return
	end function
	function getvecTensor(A,vec)Result(Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
		integer,intent(in)::vec(:)
		class(*),pointer::p(:),p0
		integer,pointer::Blockdim(:),dim(:)
		integer::index,rank
		if(.not.A%Data%getFlag())then
			call Res%empty()
			return
		end if
		rank=A%getRank()
		call A%pointDim(Dim)
		call check_indices_and_dim(dim,vec)
		if(A%getSymmetryFlag())then
			call A%ClassPointer(p,vec)
			if(associated(p))then
				call WorkingMemory%check()
				call WorkingMemory%allocate(1,rank)
				call WorkingMemory%get_memory(Blockdim,rank)
				Blockdim=A%getBlockDim(vec)
				call Res%allocate(Blockdim,A%getType())
				call WorkingMemory%free()
				call FastCopyArray(Res%Data%ClassData,p,Res%getTotalData())
			else
				call Res%empty()
			end if 
		else
			call Res%allocate([1],A%getType())
			index=addressToIndes(vec,Dim)
			call ClassPointer1DFunc(A%Data%ClassData,index,p0)
			call ModifyaValue(Res%Data%ClassData,p0,1)
		end if
	end function
	function SymgetithTensor(A,ith)Result(Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
		integer,intent(in)::ith
		if(.not.A%Data%getFlag())then
			call Res%empty()
			return
		end if
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in getting the ith block , the input tensor is not of symmetry',-1)
			call A%diminfo(.true.)
			call error_stop
		end if
		Res=getithTensor(A,ith)
		return
	end function
	function SymgetvecTensor(A,vec)Result(Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
		integer,intent(in)::vec(:)
		if(.not.A%Data%getFlag())then
			call Res%empty()
			return
		end if
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in getting the ith block , the input tensor is not of symmetry',-1)
			call A%diminfo(.true.)
			call error_stop
		end if
		Res=getvecTensor(A,vec)
	end function
	function getValueAllTensor(A)Result(Res)
		class(Tensor),intent(in)::A
		type(Tensor)::Res
		call Res%allocate([A%getTotalData()],A%getType())
		call FastCopyArray(Res%Data%ClassData,A%Data%ClassData,A%getTotalData())
		return
	end function


	type(Tensor) function fermi_element_vec(fH,QNindex,Degindex)result(subfH)
		class(Tensor),intent(in)::fH
		integer,intent(in)::QNindex(:),Degindex(:)
		integer::rank,i
		type(QuanNum),allocatable::QN(:)
		type(Tensor)::TMP
		if(.not.fH%getSymmetryFlag())then
			call writemess('The input tensor should be a symmetry tensor',-1)
			call error_stop
		end if
		rank=fH%getRank()
		allocate(QN(rank))
		do i=1,rank
			call QN(i)%setQN([fH%getQN(i,QNindex(i))])
			call QN(i)%setDeg([1])
			if(fH%getFermiFlag())then
				call QN(i)%setFermiArrow(fH%getFermiArrow(i))
			end if
		end do
		call subfH%allocateMomery(QN,1,fH%getType())
		call subfH%setBlockMomery(1,1)
		TMP=fH%Blocki(QNindex)
		call subfH%setValue(1,TMP%i(Degindex))
		if(fH%getNameFlag())then
			do i=1,rank
				call subfH%setName(i,fH%getName(i))
			end do
		end if
		return
	end function

	type(Tensor) function fermi_element_vec2(fH,index)result(subfH)
		class(Tensor),intent(in)::fH
		integer,intent(in)::index(:)
		integer,allocatable::QNindex(:),Degindex(:),dimen(:)
		integer::rank,i
		type(QuanNum),allocatable::QN(:)
		rank=fH%getRank()
		if(fH%getSymmetryFlag())then
			allocate(QNindex(rank))
			allocate(Degindex(rank))
			call fH%index2QNinfo(QNindex,Degindex,index)
			subfH=fermi_element_vec(fH,QNindex,Degindex)
			return
		end if
		subfH=fH%Ti(index)
		allocate(dimen(rank))
		dimen=1
		call subfH%resetdim(dimen)
		if(fH%getNameFlag())then
			do i=1,rank
				call subfH%setName(i,fH%getName(i))
			end do
		end if
		return
	end function


	subroutine store_All_Data(A,AllData)
		class(Tensor),intent(inout)::A
		type(Tensor),intent(in)::AllData(:)
		integer::i,TotalData,BlockLength,si,ei
		integer::NumAllData,Datai
		NumAllData=size(AllData)
		TotalData=0
		do i=1,NumAllData
			TotalData=TotalData+AllData(i)%getTotalData()
		end do
		call A%empty()
		call A%allocate([TotalData],AllData(1)%getType(),AllData(1)%Data%DataCharacterLen)
		si=0
		ei=0
		do Datai=1,NumAllData
			if(AllData(Datai)%getFlag())then
				BlockLength=AllData(Datai)%getTotalData()
				ei=si+BlockLength
				si=si+1
				call ModifySomeValue(A%Data%ClassData,AllData(Datai)%Data%ClassData,si,ei)
				si=ei
			end if
		end do
		if(ei.ne.TotalData)then
			call writemess(' Internal error in store_All_Data',-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine store_a_Tensor(Res,T)
		class(Tensor),intent(inout)::Res
		type(Tensor),intent(in)::T
		integer::TotalData
		if(.not.T%getFlag())then
			call Res%empty()
			return
		end if
		TotalData=T%getTotalData()
		call Res%empty()
		call Res%allocate([TotalData],T%getType(),T%Data%DataCharacterLen)
		call FastcopyARRAY(Res%Data%ClassData,T%Data%ClassData,TotalData)
		return
	end subroutine

	subroutine distribute_a_tensor(T,Res)
		class(Tensor),intent(in)::T
		type(Tensor),intent(inout)::Res
		integer::TotalData
		if(.not.T%getFlag())then
			call Res%empty()
			return
		end if
		TotalData=T%getTotalData()
		if(Res%getTotalData().ne.TotalData)then
			call writemess('ERROR in distribute,TotalData',-1)
			call error_stop
		end if
		call FastcopyARRAY(Res%Data%ClassData,T%Data%ClassData,TotalData)
		return
	end subroutine

	subroutine distribute_All_Data(A,AllData)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::AllData(:)
		integer::i,TotalData,BlockLength,si,ei
		integer::NumAllData,Datai
		class(*),pointer::Ap(:)
		NumAllData=size(AllData)
		TotalData=A%getTotalData()
		si=0
		ei=0
		do Datai=1,NumAllData
			if(AllData(Datai)%getFlag())then
				BlockLength=AllData(Datai)%getTotalData()
				ei=si+BlockLength
				si=si+1
				if(ei.gt.TotalData)then
					call writemess(' error in distribute_All_Data',-1)
					call error_stop
				end if
				call ClassPointer1DFunc(A%Data%ClassData,si,ei,Ap)
				call FastcopyARRAY(AllData(Datai)%Data%ClassData,Ap,BlockLength)
				si=ei
			end if
		end do
		if(ei.ne.TotalData)then
			call writemess(' error in distribute_All_Data',-1)
			call error_stop
		end if
		return
	end subroutine
