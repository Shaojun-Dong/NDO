	!ArrowSignCorrect
		!fermi-arrow=1     mean the leg is <a| :  --->----[A]
		!fermi-arrow=-1 	mean the leg is |b> :  ---<----[B]
		!
		!   the contraction should be    <a|b>
		!   if not do the Correction
		!      __         __
		!     |  |       |  |
		!     |A |--->---|B |
		!     |  |--->---|  | 
		!     |  |---<---|  |
		!     |__|       |__|
		!
		!the fermi-arrow of A are   -1,  -1,  1
		!the fermi-arrow of B are    1,   1, -1
		!The first two leg will not cause the SignCorrect
		! the last leg, fermi-arrow of A is 1 and fermi-arrow of B is -1, it cause SignCorrect
		!
		!
		!If fermi-arrow=0 , will not monitor the order
	
	subroutine BackwardForContract(Res,A,vecA)
		type(Tensor),intent(inout)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::vecA(:)
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),Wlegs(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:),invervec(:)
		integer::rank,i,counter,vec_len,Arrow
		logical::goon,Flag
		rank=A%getRank()
		vec_len=size(vecA)
		Flag=Check_not_permute_backward_case1(vecA,vec_len,rank)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		call WorkingMemory%get_memory(invervec,vec_len)
		counter=0
		newOrder(rank-vec_len+1:rank)=vecA
		do i=1,rank-vec_len
			counter=counter+1
			goon=if_In_the_array(counter,vecA,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vecA,vec_len)
			end do
			if(counter.gt.rank)then
				call writemess('ERROR in PermuteBackWard',-1)
				call writemess(newOrder)
				call writemess(vecA)
				call writemess("i="+i+',counter='+counter+',rank='+rank)
				call error_stop
			end if
			newOrder(i)=counter
		end do
		counter=0
		do i=vec_len,1,-1
			counter=counter+1
			invervec(counter)=vecA(i)
		end do

		Res%Dimension=A%Dimension
		if(Flag)then
			Res%Data=A%Data
		else
			call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		end if
		do i=1,vec_len
			Arrow=Res%getFermiArrow(vecA(i))
			if(Arrow.gt.0)then
				call Contract_fix_Sgin(Res,vecA(i))
			end if
		end do
		call ArrayData_backward_sign(Res,invervec,legCheckith,workingindex,Wlegs,vec_len)
		if(Flag)then
			call WorkingMemory%free()
			return
		end if
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine BackwardOnelegForContract(Res,A,ith)
		type(Tensor),intent(inout)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),legCheckith(:)
		integer,pointer::ei(:),dim_i(:)
		integer::rank,i,counter,vec_len,Arrow
		logical::goon,Flag
		rank=A%getRank()
		Flag=ith.eq.rank
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank-ith+1)
		counter=0
		newOrder(rank)=ith
		do i=1,rank-1
			counter=counter+1
			if(counter.eq.ith)then
				counter=counter+1
			end if
			if(counter.gt.rank)then
				call writemess('ERROR in PermuteBackWard_ith',-1)
				call error_stop
			end if
			newOrder(i)=counter
		end do
		Res%Dimension=A%Dimension
		if(Flag)then
			Res%Data=A%Data
		else
			call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		end if
		Arrow=Res%getFermiArrow(ith)
		if(Arrow.gt.0)then
			call Contract_fix_Sgin(Res,ith)
		end if
		call ArrayData_backward_one_sign(Res,ith,legCheckith,workingindex)
		if(Flag)then
			call WorkingMemory%free()
			return
		end if
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine ForwardForContract(Res,A,vecA)
		type(Tensor),intent(inout)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::vecA(:)
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),Wlegs(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:),invervec(:)
		integer::rank,i,counter,vec_len,Arrow
		logical::goon,Flag
		rank=A%getRank()
		vec_len=size(vecA)
		Flag=Check_not_permute_forward_case1(vecA,vec_len)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		call WorkingMemory%get_memory(invervec,vec_len)
		counter=0
		newOrder(1:vec_len)=vecA
		do i=vec_len+1,rank
			counter=counter+1
			goon=if_In_the_array(counter,vecA,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vecA,vec_len)
			end do
			newOrder(i)=counter
		end do

		counter=0
		do i=vec_len,1,-1
			counter=counter+1
			invervec(counter)=vecA(i)
		end do

		Res%Dimension=A%Dimension
		if(Flag)then
			Res%Data=A%Data
		else
			call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		end if
		do i=1,vec_len
			Arrow=Res%getFermiArrow(vecA(i))
			if(Arrow.lt.0)then
				call Contract_fix_Sgin(Res,vecA(i))
			end if
		end do
		call ArrayData_forard_sign(Res,invervec,legCheckith,workingindex,Wlegs,vec_len)
		if(Flag)then
			call WorkingMemory%free()
			return
		end if
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine ForwardOneLegForContract(Res,A,ith)
		type(Tensor),intent(inout)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		integer::rank,i,counter,vec_len,Arrow
		logical::goon,Flag
		rank=A%getRank()
		Flag=ith.eq.1
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+ith)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,ith)
		counter=0
		newOrder(1)=ith
		do i=2,rank
			counter=counter+1
			if(counter.eq.ith)then
				counter=counter+1
			end if
			newOrder(i)=counter
		end do

		Res%Dimension=A%Dimension
		if(Flag)then
			Res%Data=A%Data
		else
			call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		end if
		Arrow=Res%getFermiArrow(ith)
		if(Arrow.lt.0)then
			call Contract_fix_Sgin(Res,ith)
		end if
		call ArrayData_forard_one_sign(Res,ith,legCheckith,workingindex)
		if(Flag)then
			call WorkingMemory%free()
			return
		end if
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine ContractCheck(A,vecA,B,vecB)
		type(Tensor),intent(in)::A,B
		integer,intent(in)::vecA(:),vecB(:)
		integer::lenVec,i
		integer,pointer::ip1(:),ip2(:)
		lenVec=size(vecA)
		if(size(vecB).ne.lenVec)then
			call writemess('ERROR in contract,length of legs',-1)
			call writemess(vecA,-1)
			call writemess(vecB,-1)
			call error_stop
		end if
		if(A%getFermiFlag().neqv.B%getFermiFlag())then
			call writemess('ERROR in getFermiFlag',-1)
			call A%diminfo(.true.)
			call B%diminfo(.true.)
			call error_stop
		end if
		if(A%getSymmetryFlag().neqv.B%getSymmetryFlag())then
			call writemess('ERROR in getSymmetryFlag',-1)
			call A%diminfo(.true.)
			call B%diminfo(.true.)
			call error_stop
		end if
		if(A%getFermiFlag())then
			call A%pointArrow(ip1)
			call B%pointArrow(ip2)
			do i=1,lenVec
				if(mointer_order_flag)then
					!if(ip1(vecA(i))*ip2(vecB(i)).ge.0)then
					!	call writemess('ERROR in fermi-arrow',-1)
					!	call writemess(vecA)
					!	call writemess(vecB)
					!	call A%diminfo(.true.)
					!	call B%diminfo(.true.)
					!	call error_stop
					!end if
					if((ip1(vecA(i)).ge.0).and.(ip2(vecB(i)).ge.0))then
						call writemess('ERROR in fermi-arrow',-1)
						call writemess(vecA)
						call writemess(vecB)
						call A%diminfo(.true.)
						call B%diminfo(.true.)
						call error_stop
					end if
					if((ip1(vecA(i)).le.0).and.(ip2(vecB(i)).le.0))then
						call writemess('ERROR in fermi-arrow',-1)
						call writemess(vecA)
						call writemess(vecB)
						call A%diminfo(.true.)
						call B%diminfo(.true.)
						call error_stop
					end if

				else
					!if(ip1(vecA(i))*ip2(vecB(i)).gt.0)then
					!	call writemess('ERROR in fermi-arrow',-1)
					!	call writemess(vecA)
					!	call writemess(vecB)
					!	call A%diminfo(.true.)
					!	call B%diminfo(.true.)
					!	call error_stop
					!end if
					if((ip1(vecA(i)).gt.0).and.(ip2(vecB(i)).gt.0))then
						call writemess('ERROR in fermi-arrow',-1)
						call writemess(vecA)
						call writemess(vecB)
						call A%diminfo(.true.)
						call B%diminfo(.true.)
						call error_stop
					end if
					if((ip1(vecA(i)).gt.0).and.(ip2(vecB(i)).gt.0))then
						call writemess('ERROR in fermi-arrow',-1)
						call writemess(vecA,-1)
						call writemess(vecB,-1)
						call A%diminfo(.true.)
						call B%diminfo(.true.)
						call error_stop
					end if
				end if
			end do
		end if

		if(A%getSymmetryFlag())then
			call checkSymmetryRule(A%Dimension,vecA,B%Dimension,vecB)
		end if
		return
	end subroutine

	subroutine ArrowSignCorrect(A,inA,vecA,B,inB,vecB)
		type(Tensor),target,intent(in)::inA,inB
		type(Tensor),target,intent(inout)::A,B
		integer,intent(in)::vecA(:),vecB(:)
		class(Tensor),pointer::Tp,inTp
		integer::TotalDataA,TotalDataB
		Tp=>A
		inTp=>inA
		if(associated(Tp,inTp))then
			call writemess('input Tensors can not be the same variable',-1)
			call error_stop
		end if
		Tp=>B
		inTp=>inB
		if(associated(Tp,inTp))then
			call writemess('input Tensors can not be the same variable',-1)
			call error_stop
		end if
		
		if(inA%getFermiFlag())then
			TotalDataA=inA%getTotalData()
			TotalDataB=inB%getTotalData()
			if(TotalDataA.le.TotalDataB)then
				call BackwardForContract(A,inA,vecA)
				call inB%forward(vecB,B)
			else
				call inA%backward(vecA,A)
				call ForwardForContract(B,inB,vecB)
			end if
		else
			call inA%backward(vecA,A)
			call inB%forward(vecB,B)
		end if
		return
	end subroutine

	subroutine ArrowSignCorrectOneLeg(A,inA,ithA,B,inB,ithB)
		type(Tensor),target,intent(in)::inA,inB
		type(Tensor),target,intent(inout)::A,B
		integer,intent(in)::ithA,ithB
		class(Tensor),pointer::Tp,inTp
		integer::TotalDataA,TotalDataB
		Tp=>A
		inTp=>inA
		if(associated(Tp,inTp))then
			call writemess('input Tensors can not be the same variable',-1)
			call error_stop
		end if
		Tp=>B
		inTp=>inB
		if(associated(Tp,inTp))then
			call writemess('input Tensors can not be the same variable',-1)
			call error_stop
		end if
		if(inA%getFermiFlag())then
			TotalDataA=inA%getTotalData()
			TotalDataB=inB%getTotalData()
			if(TotalDataA.le.TotalDataB)then
				call BackwardOneLegForContract(A,inA,ithA)
				call inB%forward(ithB,B)
			else
				call inA%backward(ithA,A)
				call ForwardOneLegForContract(B,inB,ithB)
			end if
		else
			call inA%backward(ithA,A)
			call inB%forward(ithB,B)
		end if
		return
	end subroutine
	
	subroutine un_set_mointer_order_flag()
		mointer_order_flag=.false.
	end subroutine
	subroutine set_mointer_order_flag()
		mointer_order_flag=.true.
	end subroutine

	!******************  contract  *********************
		!	T1:[i1,i2,i3,i4,i5,i6,i7,i8]
		!	T2:[j1,j2,j3,j4,j5,j6,j7,j8,j9,j10]
		!	i1=(/5,1,2/)
		!	i2=(/10,3,5/)
		!	then the result will be T1'*T2'
		!	where ,
		!	T1'=[i3,i4,i6,i7,i8,(i5*i1*i2)]
		!	T2'=[(j10*j3*j5),j1,j2,j4,j6,j7,j8,j9]
		!  The sign of T1' will be determine by T1, T1 permute to (/2,1,5/) with fermi rule **********NOTE here not (/5,1,2/) *************
		!         And the permute to [i3,i4,i6,i7,i8,(i5*i1*i2)] non-fermi rule
		!  T2' will be determine by T1, T1 permute to (/10,3,5/)  with fermi rule 
		!  example   
		!       C1*C2 |11> = C1*C2 * C1^+ * C2^+ |0> = -C2*C1 * C1^+ * C2^+ |0> = -|00>
		!     The order of C1*C2 shoud reoder as C2*C1
		!     but the data store in memery is C1*C2
		!     so the function permute to C2*C1 with fermi rule, and then permute to C1*C2 with non-fermi rule
		!
		! 	input Tensor should be in its original dimenison,there is no contract on it
		!	if present len_of_contract, len_of_contract specify the length of  i1, and i2

	function contract_vec(A,vecA,B,vecB) result(T)
		type(Tensor)::T
		type(Tensor),intent(in) :: A,B
		integer,intent(in) :: vecA(:),vecB(:)
		if(.not.A%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		call ContractCheck(A,vecA,B,vecB)
		call ArrowSignCorrect(TMPContract1,A,vecA,TMPContract2,B,vecB)
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0,size(vecA))
		return
	end function

	function contract_ith(A,ithA,B,ithB) result(T)
		type(Tensor)::T
		type(Tensor),intent(in) :: A,B
		integer,intent(in) :: ithA,ithB
		if(.not.A%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		call ContractCheck(A,[ithA],B,[ithB])
		call ArrowSignCorrectOneLeg(TMPContract1,A,ithA,TMPContract2,B,ithB)
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0)
		return
	end function


	function contract_ChaVec(A,wA,B,wB) result(T)
		type(Tensor)::T
		type(Tensor),intent(in) :: A,B
		character(len=*),intent(in) :: wA(:),wB(:)
		integer::len_vec
		if(.not.A%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		len_vec=size(wA)
		call allocateCheck(TMPi1,len_vec)
		call allocateCheck(TMPi2,len_vec)
		TMPi1=A%FindOrder(wA)
		TMPi2=B%FindOrder(wB)
		call ContractCheck(A,TMPi1(1:len_vec),B,TMPi2(1:len_vec))
		call ArrowSignCorrect(TMPContract1,A,TMPi1(1:len_vec),TMPContract2,B,TMPi2(1:len_vec))
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0,len_vec)
		return
	end function

	function contract_Cha(A,wA,B,wB) result(T)
		type(Tensor)::T
		type(Tensor),intent(in) :: A,B
		character(len=*),intent(in) :: wA,wB
		integer::iA,iB
		if(.not.A%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		iA=A%FindOrder(wA)
		iB=B%FindOrder(wB)
		call ContractCheck(A,[iA],B,[iB])
		call ArrowSignCorrectOneLeg(TMPContract1,A,iA,TMPContract2,B,iB)
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0)
		return
	end function

	function contract_Same_name(A,B) result(T)
		type(Tensor)::T
		type(Tensor),intent(in) :: A,B
		character(len=len_of_Name),pointer::name1(:),name2(:)
		integer::lenofname,rank1,rank2
		if(.not.A%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		rank1=A%getRank()
		rank2=B%getRank()
		call allocateCheck(TMPa1,max(rank1,rank2))
		call A%pointName(name1)
		call B%pointName(name2)
		call find_same_name(name1,name2,TMPa1,lenofname) 

		call allocateCheck(TMPi1,lenofname)
		call allocateCheck(TMPi2,lenofname)
		TMPi1=A%FindOrder(TMPa1(1:lenofname))
		TMPi2=B%FindOrder(TMPa1(1:lenofname))

		call ContractCheck(A,TMPi1(1:lenofname),B,TMPi2(1:lenofname))
		call ArrowSignCorrect(TMPContract1,A,TMPi1(1:lenofname),TMPContract2,B,TMPi2(1:lenofname))
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0,lenofname)
		return
	end function

	subroutine find_same_name(inname1,inname2,SameName,lenofname)
		character(len=*),intent(in)::inname1(:),inname2(:)
		character(len=*),intent(inout)::SameName(:)
		integer,intent(inout)::lenofname
		integer::i,j,k,sizeSameName
		k=0
		lenofname=0
		sizeSameName=size(SameName)
		do i=1,size(inname1)
			do j=1,size(inname2)
				if(inname1(i).equ.inname2(j))then
					k=k+1
					if(k.gt.sizeSameName)then
						call writemess('ERROR in finding the same name, maybe something wrong in the tensorname',-1)
						call writemess('You can diminfo to see if the tensor names are right')
						call writemess('It is not allow to have two or more same name in one tensor')
						call error_stop()
					end if
					SameName(k)=inname1(i)
					lenofname=lenofname+1
				end if
			end do
		end do
		return
	end subroutine



	!******************  contract subroutine *********************

	subroutine contract_vecT(T,A,vecA,B,vecB)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(in) :: A,B
		integer,intent(in) :: vecA(:),vecB(:)
		if(.not.A%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		call ContractCheck(A,vecA,B,vecB)
		call ArrowSignCorrect(TMPContract1,A,vecA,TMPContract2,B,vecB)
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0,size(vecA))
		return
	end subroutine

	subroutine contract_vecA(T,vecA,B,vecB)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(in) :: B
		integer,intent(in) :: vecA(:),vecB(:)
		if(.not.T%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		call ContractCheck(T,vecA,B,vecB)
		call ArrowSignCorrect(TMPContract1,T,vecA,TMPContract2,B,vecB)
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0,size(vecA))
		return
	end subroutine

	subroutine contract_vecB(T,A,vecA,vecB)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(in) :: A
		integer,intent(in) :: vecA(:),vecB(:)
		if(.not.A%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		call ContractCheck(A,vecA,T,vecB)
		call ArrowSignCorrect(TMPContract1,A,vecA,TMPContract2,T,vecB)
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0,size(vecA))
		return
	end subroutine

	subroutine contract_ithT(T,A,ithA,B,ithB)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(in) :: A,B
		integer,intent(in) :: ithA,ithB
		if(.not.A%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		call ContractCheck(A,[ithA],B,[ithB])
		call ArrowSignCorrectOneLeg(TMPContract1,A,ithA,TMPContract2,B,ithB)
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0)
		return
	end subroutine

	subroutine contract_ithA(T,ithA,B,ithB)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(in) :: B
		integer,intent(in) :: ithA,ithB
		if(.not.T%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		call ContractCheck(T,[ithA],B,[ithB])
		call ArrowSignCorrectOneLeg(TMPContract1,T,ithA,TMPContract2,B,ithB)
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0)
		return
	end subroutine

	subroutine contract_ithB(T,A,ithA,ithB)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(in) :: A
		integer,intent(in) :: ithA,ithB
		if(.not.A%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		call ContractCheck(A,[ithA],T,[ithB])
		call ArrowSignCorrectOneLeg(TMPContract1,A,ithA,TMPContract2,T,ithB)
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0)
		return
	end subroutine

	subroutine contract_ChaVecT(T,A,wA,B,wB)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(in) :: A,B
		character(len=*),intent(in) :: wA(:),wB(:)
		integer::len_vec
		if(.not.A%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		len_vec=size(wA)
		call allocateCheck(TMPi1,len_vec)
		call allocateCheck(TMPi2,len_vec)
		TMPi1=A%FindOrder(wA)
		TMPi2=B%FindOrder(wB)
		call ContractCheck(A,TMPi1(1:len_vec),B,TMPi2(1:len_vec))
		call ArrowSignCorrect(TMPContract1,A,TMPi1(1:len_vec),TMPContract2,B,TMPi2(1:len_vec))
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0,len_vec)
		return
	end subroutine

	subroutine contract_ChaVecA(T,wA,B,wB)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(in) :: B
		character(len=*),intent(in) :: wA(:),wB(:)
		integer::len_vec
		if(.not.T%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		len_vec=size(wA)
		call allocateCheck(TMPi1,len_vec)
		call allocateCheck(TMPi2,len_vec)
		TMPi1=T%FindOrder(wA)
		TMPi2=B%FindOrder(wB)
		call ContractCheck(T,TMPi1(1:len_vec),B,TMPi2(1:len_vec))
		call ArrowSignCorrect(TMPContract1,T,TMPi1(1:len_vec),TMPContract2,B,TMPi2(1:len_vec))
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0,len_vec)
		return
	end subroutine

	subroutine contract_ChaVecB(T,A,wA,wB)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(in) :: A
		character(len=*),intent(in) :: wA(:),wB(:)
		integer::len_vec
		if(.not.A%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		len_vec=size(wA)
		call allocateCheck(TMPi1,len_vec)
		call allocateCheck(TMPi2,len_vec)
		TMPi1=A%FindOrder(wA)
		TMPi2=T%FindOrder(wB)
		call ContractCheck(A,TMPi1(1:len_vec),T,TMPi2(1:len_vec))
		call ArrowSignCorrect(TMPContract1,A,TMPi1(1:len_vec),TMPContract2,T,TMPi2(1:len_vec))
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0,len_vec)
		return
	end subroutine

	subroutine contract_ChaT(T,A,wA,B,wB)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(in) :: A,B
		character(len=*),intent(in) :: wA,wB
		integer::iA,iB
		if(.not.A%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		iA=A%FindOrder(wA)
		iB=B%FindOrder(wB)
		call ContractCheck(A,[iA],B,[iB])
		call ArrowSignCorrectOneLeg(TMPContract1,A,iA,TMPContract2,B,iB)
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0)
		return
	end subroutine

	subroutine contract_ChaA(T,wA,B,wB)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(in) :: B
		character(len=*),intent(in) :: wA,wB
		integer::iA,iB
		if(.not.T%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		iA=T%FindOrder(wA)
		iB=B%FindOrder(wB)
		call ContractCheck(T,[iA],B,[iB])
		call ArrowSignCorrectOneLeg(TMPContract1,T,iA,TMPContract2,B,iB)
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0)
		return
	end subroutine

	subroutine contract_ChaB(T,A,wA,wB)
		class(Tensor),intent(inout)::T
		class(Tensor),intent(in) :: A
		character(len=*),intent(in) :: wA,wB
		integer::iA,iB
		if(.not.A%getFlag())then
			call writemess('There is no data in the first fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		if(.not.T%getFlag())then
			call writemess('There is no data in the second fTensor, when contracting fTensor',-1)
			call error_stop()
		end if
		iA=A%FindOrder(wA)
		iB=T%FindOrder(wB)
		call ContractCheck(A,[iA],T,[iB])
		call ArrowSignCorrectOneLeg(TMPContract1,A,iA,TMPContract2,T,iB)
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0)
		return
	end subroutine

	subroutine contract_Same_nameA(T,B) 
		class(Tensor),intent(inout)::T
		class(Tensor),intent(in) :: B
		character(len=len_of_Name),pointer::name1(:),name2(:)
		integer::lenofname,rank1,rank2
		if(.not.T%getFlag())then
			call writemess('There is no data in the first Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		if(.not.B%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		rank1=T%getRank()
		rank2=B%getRank()
		call allocateCheck(TMPa1,max(rank1,rank2))
		call T%pointName(name1)
		call B%pointName(name2)
		call find_same_name(name1,name2,TMPa1,lenofname) 

		call allocateCheck(TMPi1,lenofname)
		call allocateCheck(TMPi2,lenofname)
		TMPi1=T%FindOrder(TMPa1(1:lenofname))
		TMPi2=B%FindOrder(TMPa1(1:lenofname))

		call ContractCheck(T,TMPi1(1:lenofname),B,TMPi2(1:lenofname))
		call ArrowSignCorrect(TMPContract1,T,TMPi1(1:lenofname),TMPContract2,B,TMPi2(1:lenofname))
		call T%empty()
		call ProductTensorRoutine(T,TMPContract1,TMPContract2,1,0,lenofname)
		return
	end subroutine






	subroutine contract_ownlegs_routine(T,ith1,ith2) 
		class(Tensor)::T
		type(Tensor),pointer::pT
		integer,intent(in)::ith1,ith2
		type(Dimension)::NewDimen
		integer::rank,classtype,dim1,i,k
		integer,pointer::idata(:,:,:),newidata(:)
		real*4,pointer::sdata(:,:,:),newsdata(:)
		real*8,pointer::ddata(:,:,:),newddata(:)
		complex*8,pointer::cdata(:,:,:),newcdata(:)
		complex*16,pointer::zdata(:,:,:),newzdata(:)
		if(.not.T%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		rank=T%getRank()
		if(rank.eq.2)then
			T=T%trace()
			return
		end if
		pT=>TMPContract1
		pT=T.pf.ith1
		call pT%forward(ith2)
		dim1=pT%dim(1)
		if(dim1.ne.pT%dim(2))then
			call writemess(' ERROR in contract(ith1,ith2), dimension')
			call error_stop
		end if
		NewDimen=pT%Dimension.subdim.[3,rank]
		classtype=T%getType()
		call T%empty()
		call T%allocate(NewDimen,classtype)
		call pT%fuse(3,rank)

		select case (classtype)
			case (1)
				call pT%pointer(idata)
				call T%pointer(newidata)
				do k=1,T%getTotalData()
					newidata(k)=0
					do i=1,dim1
						newidata(k)=newidata(k)+idata(i,i,k)
					end do
				end do
			case (2)
				call pT%pointer(sdata)
				call T%pointer(newsdata)
				do k=1,T%getTotalData()
					newsdata(k)=0
					do i=1,dim1
						newsdata(k)=newsdata(k)+sdata(i,i,k)
					end do
				end do
			case(3)
				call pT%pointer(ddata)
				call T%pointer(newddata)
				do k=1,T%getTotalData()
					newddata(k)=0
					do i=1,dim1
						newddata(k)=newddata(k)+ddata(i,i,k)
					end do
				end do
			case(4)
				call pT%pointer(cdata)
				call T%pointer(newcdata)
				do k=1,T%getTotalData()
					newcdata(k)=0
					do i=1,dim1
						newcdata(k)=newcdata(k)+cdata(i,i,k)
					end do
				end do
			case(5)
				call pT%pointer(zdata)
				call T%pointer(newzdata)
				do k=1,T%getTotalData()
					newzdata(k)=0
					do i=1,dim1
						newzdata(k)=newzdata(k)+zdata(i,i,k)
					end do
				end do
			case default
				call writemess(' ERROR in contract(name1,name2), clasee type')
				call error_stop
		end select
		return
	end subroutine

	type(Tensor) function contract_ownlegs(Tin,ith1,ith2) Result(Res)
		type(Tensor),intent(in)::Tin
		type(Tensor),pointer::pT
		integer,intent(in)::ith1,ith2
		type(Dimension)::NewDimen
		integer::rank,classtype,dim1,i,k
		integer,pointer::idata(:,:,:),newidata(:)
		real*4,pointer::sdata(:,:,:),newsdata(:)
		real*8,pointer::ddata(:,:,:),newddata(:)
		complex*8,pointer::cdata(:,:,:),newcdata(:)
		complex*16,pointer::zdata(:,:,:),newzdata(:)
		if(.not.Tin%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		rank=Tin%getRank()
		if(rank.eq.2)then
			Res=Tin%trace()
			return
		end if
		pT=>TMPContract1
		pT=Tin.pf.ith1
		call pT%forward(ith2)
		dim1=pT%dim(1)
		if(dim1.ne.pT%dim(2))then
			call writemess(' ERROR in contract(ith1,ith2), dimension')
			call error_stop
		end if
		NewDimen=pT%Dimension.subdim.[3,rank]
		classtype=Tin%getType()
		call Res%empty()
		call Res%allocate(NewDimen,classtype)
		call pT%fuse(3,rank)

		select case (classtype)
			case (1)
				call pT%pointer(idata)
				call Res%pointer(newidata)
				do k=1,Res%getTotalData()
					newidata(k)=0
					do i=1,dim1
						newidata(k)=newidata(k)+idata(i,i,k)
					end do
				end do
			case (2)
				call pT%pointer(sdata)
				call Res%pointer(newsdata)
				do k=1,Res%getTotalData()
					newsdata(k)=0
					do i=1,dim1
						newsdata(k)=newsdata(k)+sdata(i,i,k)
					end do
				end do
			case(3)
				call pT%pointer(ddata)
				call Res%pointer(newddata)
				do k=1,Res%getTotalData()
					newddata(k)=0
					do i=1,dim1
						newddata(k)=newddata(k)+ddata(i,i,k)
					end do
				end do
			case(4)
				call pT%pointer(cdata)
				call Res%pointer(newcdata)
				do k=1,Res%getTotalData()
					newcdata(k)=0
					do i=1,dim1
						newcdata(k)=newcdata(k)+cdata(i,i,k)
					end do
				end do
			case(5)
				call pT%pointer(zdata)
				call Res%pointer(newzdata)
				do k=1,Res%getTotalData()
					newzdata(k)=0
					do i=1,dim1
						newzdata(k)=newzdata(k)+zdata(i,i,k)
					end do
				end do
			case default
				call writemess(' ERROR in contract(name1,name2), clasee type')
				call error_stop
		end select
		return
	end function


	type(Tensor) function contract_name_ownlegs(Tin,name1,name2) Result(Res)
		type(Tensor),intent(in)::Tin
		type(Tensor),pointer::pT
		character(len=*),intent(in)::name1,name2
		type(Dimension)::NewDimen
		integer::rank,classtype,dim1,i,k
		integer,pointer::idata(:,:,:),newidata(:)
		real*4,pointer::sdata(:,:,:),newsdata(:)
		real*8,pointer::ddata(:,:,:),newddata(:)
		complex*8,pointer::cdata(:,:,:),newcdata(:)
		complex*16,pointer::zdata(:,:,:),newzdata(:)
		if(.not.Tin%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		rank=Tin%getRank()
		if(rank.eq.2)then
			Res=Tin%trace()
			return
		end if
		pT=>TMPContract1
		pT=Tin.pf.name1
		call pT%forward(name2)
		dim1=pT%dim(1)
		if(dim1.ne.pT%dim(2))then
			call writemess(' ERROR in contract(ith1,ith2), dimension')
			call error_stop
		end if
		NewDimen=pT%Dimension.subdim.[3,rank]
		classtype=Tin%getType()
		call Res%empty()
		call Res%allocate(NewDimen,classtype)
		call pT%fuse(3,rank)

		select case (classtype)
			case (1)
				call pT%pointer(idata)
				call Res%pointer(newidata)
				do k=1,Res%getTotalData()
					newidata(k)=0
					do i=1,dim1
						newidata(k)=newidata(k)+idata(i,i,k)
					end do
				end do
			case (2)
				call pT%pointer(sdata)
				call Res%pointer(newsdata)
				do k=1,Res%getTotalData()
					newsdata(k)=0
					do i=1,dim1
						newsdata(k)=newsdata(k)+sdata(i,i,k)
					end do
				end do
			case(3)
				call pT%pointer(ddata)
				call Res%pointer(newddata)
				do k=1,Res%getTotalData()
					newddata(k)=0
					do i=1,dim1
						newddata(k)=newddata(k)+ddata(i,i,k)
					end do
				end do
			case(4)
				call pT%pointer(cdata)
				call Res%pointer(newcdata)
				do k=1,Res%getTotalData()
					newcdata(k)=0
					do i=1,dim1
						newcdata(k)=newcdata(k)+cdata(i,i,k)
					end do
				end do
			case(5)
				call pT%pointer(zdata)
				call Res%pointer(newzdata)
				do k=1,Res%getTotalData()
					newzdata(k)=0
					do i=1,dim1
						newzdata(k)=newzdata(k)+zdata(i,i,k)
					end do
				end do
			case default
				call writemess(' ERROR in contract(name1,name2), clasee type')
				call error_stop
		end select
		return
	end function

	subroutine contract_name_ownlegs_routine(T,name1,name2) 
		class(Tensor)::T
		type(Tensor),pointer::pT
		character(len=*),intent(in)::name1,name2
		type(Dimension)::NewDimen
		integer::rank,classtype,dim1,i,k
		integer,pointer::idata(:,:,:),newidata(:)
		real*4,pointer::sdata(:,:,:),newsdata(:)
		real*8,pointer::ddata(:,:,:),newddata(:)
		complex*8,pointer::cdata(:,:,:),newcdata(:)
		complex*16,pointer::zdata(:,:,:),newzdata(:)
		if(.not.T%getFlag())then
			call writemess('There is no data in the second Tensor, when contracting Tensor',-1)
			call error_stop()
		end if
		rank=T%getRank()
		if(rank.eq.2)then
			T=T%trace()
			return
		end if
		pT=>TMPContract1
		pT=T.pf.name1
		call pT%forward(name2)
		dim1=pT%dim(1)
		if(dim1.ne.pT%dim(2))then
			call writemess(' ERROR in contract(name1,name2), dimension')
			call error_stop
		end if
		NewDimen=pT%Dimension.subdim.[3,rank]
		classtype=T%getType()
		call T%empty()
		call T%allocate(NewDimen,classtype)
		call pT%fuse(3,rank)

		select case (classtype)
			case (1)
				call pT%pointer(idata)
				call T%pointer(newidata)
				do k=1,T%getTotalData()
					newidata(k)=0
					do i=1,dim1
						newidata(k)=newidata(k)+idata(i,i,k)
					end do
				end do
			case (2)
				call pT%pointer(sdata)
				call T%pointer(newsdata)
				do k=1,T%getTotalData()
					newsdata(k)=0
					do i=1,dim1
						newsdata(k)=newsdata(k)+sdata(i,i,k)
					end do
				end do
			case(3)
				call pT%pointer(ddata)
				call T%pointer(newddata)
				do k=1,T%getTotalData()
					newddata(k)=0
					do i=1,dim1
						newddata(k)=newddata(k)+ddata(i,i,k)
					end do
				end do
			case(4)
				call pT%pointer(cdata)
				call T%pointer(newcdata)
				do k=1,T%getTotalData()
					newcdata(k)=0
					do i=1,dim1
						newcdata(k)=newcdata(k)+cdata(i,i,k)
					end do
				end do
			case(5)
				call pT%pointer(zdata)
				call T%pointer(newzdata)
				do k=1,T%getTotalData()
					newzdata(k)=0
					do i=1,dim1
						newzdata(k)=newzdata(k)+zdata(i,i,k)
					end do
				end do
			case default
				call writemess(' ERROR in contract(name1,name2), clasee type')
				call error_stop
		end select
		return
	end subroutine