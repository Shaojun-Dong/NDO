	
	!*********************************************************
	!
	! permutation tools
	!
	!*********************************************************

	function if_In_the_array(ith,vec,len_vec)result(Res)!***********************************************
		logical::Res
		integer,intent(in)::ith,len_vec
		integer,intent(in)::vec(len_vec)
		integer::i
		Res=.false.
		do i=1,len_vec
			if(ith.eq.vec(i))then
				Res=.true.
				return
			end if
		end do
		return
	end function

	function Check_not_permute_forward_case1(vec,lenvec)
		logical::Check_not_permute_forward_case1
		integer,intent(in)::vec(:),lenvec
		integer,pointer::check(:)
		integer::i
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,lenvec)
		call WorkingMemory%get_memory(check,lenvec)
		do i=1,lenvec
			check(i)=i
		end do
		Check_not_permute_forward_case1=vec.equ.check
		call WorkingMemory%free()
		return
	end function

	function Check_not_permute_forward_case2(vec,lenvec,check)
		logical::Check_not_permute_forward_case2
		integer,intent(in)::vec(:),lenvec
		integer,intent(inout)::check(:)
		integer::i
		do i=1,lenvec
			check(i)=i
		end do
		Check_not_permute_forward_case2=vec.equ.check
		return
	end function

	function Check_not_permute_backward_case1(vec,lenvec,rank)
		logical::Check_not_permute_backward_case1
		integer,intent(in)::vec(:),lenvec,rank
		integer,pointer::check(:)
		integer::i,ii
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,lenvec)
		call WorkingMemory%get_memory(check,lenvec)
		ii=0
		do i=rank-lenvec+1,rank
			ii=ii+1
			check(ii)=i
		end do
		Check_not_permute_backward_case1=vec.equ.check
		call WorkingMemory%free()
		return
	end function
	function Check_not_permute_backward_case2(vec,lenvec,rank,check)
		logical::Check_not_permute_backward_case2
		integer,intent(in)::vec(:),lenvec,rank
		integer,intent(inout)::check(:)
		integer::i,ii
		ii=0
		do i=rank-lenvec+1,rank
			ii=ii+1
			check(ii)=i
		end do
		Check_not_permute_backward_case2=vec.equ.check
		return
	end function

	subroutine forward_one_sign(Block,dimen,ithleg,legCheckith,indices,maxinde,FirstIndex)
		integer,intent(inout)::Block(:),indices(:)
		integer,intent(in)::maxinde(:),legCheckith(:),ithleg
		type(Dimension),intent(in)::dimen
		logical,intent(in)::FirstIndex
		integer::iQN,i
		logical::goon,ifp,testrule
		if(FirstIndex)return
		indices=1
		goon=.true.
		i=0
		do while (goon)
			i=i+1
			if(block(i).ne.0)then
				if(check_same_name_flag)then
					testrule = if_symmetry_Rule(dimen,indices)
					if(.not.testrule)then
						call writemess('ERRRO symmetry rule',-1)
						call writemess('indices=')
						call writemess(indices,'I4')
						call dimen%diminfo(.true.)
						call error_stop
					end if
				end if
				call QaunNumParity(iQN,dimen,ithleg,indices(ithleg) )
				if(iQN .eq.(-1))then
					call ifParity(ifp,dimen,indices,legCheckith)
					if(ifp) then
						block(i)=-block(i)
					end if
				end if
			end if
			goon=index_counter(indices,maxinde)
		end do
		return
	end subroutine

	subroutine backward_one_sign(Block,dimen,ithleg,legCheckith,indices,maxinde,LastIndex)
		integer,intent(inout)::Block(:),indices(:)
		integer,intent(in)::maxinde(:),legCheckith(:),ithleg
		type(Dimension),intent(in)::dimen
		logical,intent(in)::LastIndex
		integer::iQN,i
		logical::goon,ifp,testrule
		if(LastIndex)return
		indices=1
		goon=.true.
		i=0
		do while (goon)
			i=i+1
			if(block(i).ne.0)then
				if(check_same_name_flag)then
					testrule = if_symmetry_Rule(dimen,indices)
					if(.not.testrule)then
						call writemess('ERRRO symmetry rule',-1)
						call writemess('indices=')
						call writemess(indices,'I4')
						call dimen%diminfo(.true.)
						call error_stop
					end if
				end if
				call QaunNumParity(iQN,dimen,ithleg,indices(ithleg) )
				if(iQN .eq.(-1))then
					call ifParity(ifp,dimen,indices,legCheckith)
					if(ifp) then
						block(i)=-block(i)
					end if
				end if
			end if
			goon=index_counter(indices,maxinde)
		end do
		return
	end subroutine
 

	subroutine forward_fermi_sign(Block,dimen,legs,legCheckith,indices,maxinde,Wlegs,lenleg)
		integer,intent(inout)::Block(:),indices(:),Wlegs(:),legCheckith(:)
		integer,intent(in)::maxinde(:),legs(:),lenleg
		type(Dimension),intent(in)::dimen
		integer::i,j,ii,ii_in_use
		logical::FirstIndex,First
		Wlegs=legs
		ii=0
		First=.true.
		do i=lenleg,1,-1
			ii=ii+1
			do j=1,i-1
				if(Wlegs(j).lt.Wlegs(i))Wlegs(j)=Wlegs(j)+1
			end do
			legCheckith(ii)=legs(i)
			ii_in_use=ii
			do j=legs(i)-1,1,-1
				if(.not.if_In_the_array(j,legCheckith(1:ii),ii))then
					ii_in_use=ii_in_use+1
					legCheckith(ii_in_use)=j
				end if
			end do
			if(First)then
				FirstIndex=legs(i).eq.1
				First=.false.
			else
				FirstIndex=.false.
			end if
			call forward_one_sign(Block,dimen,legs(i),legCheckith(1:ii_in_use),indices,maxinde,FirstIndex)
		end do
		return
	end subroutine
	subroutine backward_fermi_sign(Block,dimen,legs,legCheckith,indices,maxinde,Wlegs,lenleg)
		integer,intent(inout)::Block(:),indices(:),Wlegs(:),legCheckith(:)
		integer,intent(in)::maxinde(:),legs(:),lenleg
		type(Dimension),intent(in)::dimen
		integer::i,j,ii,ii_in_use,rank
		logical::LastIndex,First
		Wlegs=legs
		ii=0
		rank=dimen%getRank()
		First=.true.
		do i=1,lenleg
			ii=ii+1
			do j=i+1,lenleg
				if(Wlegs(j).gt.Wlegs(i))Wlegs(j)=Wlegs(j)-1
			end do
			legCheckith(ii)=legs(i)
			ii_in_use=ii
			do j=legs(i)+1,rank
				if(.not.if_In_the_array(j,legCheckith(1:ii),ii))then
					ii_in_use=ii_in_use+1
					legCheckith(ii_in_use)=j
				end if
			end do
			if(First)then
				LastIndex=legs(i).eq.dimen%getRank()
				First=.false.
			else
				LastIndex=.false.
			end if
			call backward_one_sign(Block,dimen,legs(i),legCheckith(1:ii_in_use),indices,maxinde,LastIndex)
		end do
		return
	end subroutine


	subroutine ArrayData_forard_sign(A,vec,legCheckith,workingindex,Wlegs,vec_len)
		type(Tensor),intent(inout)::A
		integer,intent(inout)::legCheckith(:),workingindex(:),Wlegs(:)
		integer,intent(in)::vec_len,vec(:)
		integer,pointer::ei(:)
		integer::i
		if(A%getFermiFlag().and.A%Data%getFlag())then
			call A%Data%pointEndi(ei)
			call forward_fermi_sign(ei,A%Dimension,vec,legCheckith,&
						workingindex,A%dim(),Wlegs,vec_len)
			call fix_fermi_sign(A)
		end if
		return
	end subroutine
	subroutine ArrayData_forard_one_sign(A,ith,legCheckith,workingindex)
		type(Tensor),intent(inout)::A
		integer,intent(inout)::legCheckith(:),workingindex(:)
		integer,intent(in)::ith
		integer,pointer::ei(:)
		integer::i
		if(A%getFermiFlag().and.A%Data%getFlag())then
			call A%Data%pointEndi(ei)
			do i=1,ith
				if((i).gt.size(legCheckith))then
					call writemess('ERROR in ArrayData_forard_one_sign',-1)
					call error_stop
				end if
				legCheckith(i)=ith-i+1
			end do
			if((ith).ne.size(legCheckith))then
				call writemess('ERROR in ArrayData_forard_one_sign',-1)
				call error_stop
			end if
			call forward_one_sign(ei,A%Dimension,ith,legCheckith(1:ith),&
							workingindex,A%dim(),ith.eq.1)
			call fix_fermi_sign(A)
		end if
	end subroutine

	subroutine ArrayData_backward_sign(A,vec,legCheckith,workingindex,Wlegs,vec_len)
		type(Tensor),intent(inout)::A
		integer,intent(inout)::legCheckith(:),workingindex(:),Wlegs(:)
		integer,intent(in)::vec_len,vec(:)
		integer,pointer::ei(:)
		integer::i
		if(A%getFermiFlag().and.A%Data%getFlag())then
			call A%Data%pointEndi(ei)
			call backward_fermi_sign(ei,A%Dimension,vec,legCheckith,&
						workingindex,A%dim(),Wlegs,vec_len)
			call fix_fermi_sign(A)
		end if
	end subroutine

	subroutine ArrayData_backward_one_sign(A,ith,legCheckith,workingindex)
		type(Tensor),intent(inout)::A
		integer,intent(inout)::legCheckith(:),workingindex(:)
		integer,intent(in)::ith
		integer,pointer::ei(:)
		integer::i,rank
		if(A%getFermiFlag().and.A%Data%getFlag())then
			rank=A%getRank()
			call A%Data%pointEndi(ei)
			do i=ith,rank
				if((i-ith+1).gt.size(legCheckith))then
					call writemess('ERROR in ArrayData_forard_one_sign',-1)
					call error_stop
				end if
				legCheckith(i-ith+1)=i
			end do
			if((rank-ith+1).ne.size(legCheckith))then
				call writemess('ERROR in ArrayData_forard_one_sign',-1)
				call writemess('size(legCheckith)='+size(legCheckith),-1)
				call writemess('(rank-ith+1)='+((rank-ith+1)),-1)
				call error_stop
			end if
			call backward_one_sign(ei,A%Dimension,ith,legCheckith,workingindex,A%dim(),ith.eq.rank)
			call fix_fermi_sign(A)
		end if
	end subroutine

	

	subroutine permutationArrayData(outDa,inDa,plan,oldDim,WorkingDim,workingindex)
		type(DataArray),target,intent(in)::inDa
		type(DataArray),target,intent(inout)::outDa
		integer, intent(in) ::  plan(:)
		integer,intent(inout)::WorkingDim(:),workingindex(:)
		type(Dimension),intent(in)::oldDim
		integer::classType,TotalBlock,i,rank
		integer,pointer::inip(:),outip(:)
		real*4,pointer::insp(:),outsp(:)
		real*8,pointer::indp(:),outdp(:)
		complex*8,pointer::incp(:),outcp(:)
		complex*16,pointer::inzp(:),outzp(:)
		logical,pointer::inlp(:),outlp(:)
		character(len=characterlen),pointer::inap(:),outap(:)
		type(DataArray),pointer::inDap,outDap
		class(*),pointer::inclp(:),outclp(:)
		inDap=>inDa
		outDap=>outDa
		if(associated(inDap,outDap))then
			call writemess('input Tensors can not be the same variable in the permutation',-1)
			call error_stop
		end if
		call outDa%empty()
		call outDa%allocate(inDa)
		TotalBlock=inDa%getTotalBlock()
		classType=inDa%GetType()
		rank=oldDim%getRank()
		select case(classType)
			case (1)
				do i=1,TotalBlock
					if(inDa%getFlag(i))then
						call inDa%pointer(inip,i)
						call outDa%pointer(outip,i)
						call GetBlockDimensionRoutine(oldDim,WorkingDim,i,workingindex)
						call tensor_transpose(inip,outip,WorkingDim,plan,rank)
					end if
				end do
			case (2)
				do i=1,TotalBlock
					if(inDa%getFlag(i))then
						call inDa%pointer(insp,i)
						call outDa%pointer(outsp,i)
						call GetBlockDimensionRoutine(oldDim,WorkingDim,i,workingindex)
						call tensor_transpose(insp,outsp,WorkingDim,plan,rank)
					end if
				end do
			case (3)
				do i=1,TotalBlock
					if(inDa%getFlag(i))then
						call inDa%pointer(indp,i)
						call outDa%pointer(outdp,i)
						call GetBlockDimensionRoutine(oldDim,WorkingDim,i,workingindex)
						call tensor_transpose(indp,outdp,WorkingDim,plan,rank)
					end if
				end do
			case (4)
				do i=1,TotalBlock
					if(inDa%getFlag(i))then
						call inDa%pointer(incp,i)
						call outDa%pointer(outcp,i)
						call GetBlockDimensionRoutine(oldDim,WorkingDim,i,workingindex)
						call tensor_transpose(incp,outcp,WorkingDim,plan,rank)
					end if
				end do
			case (5)
				do i=1,TotalBlock
					if(inDa%getFlag(i))then
						call inDa%pointer(inzp,i)
						call outDa%pointer(outzp,i)
						call GetBlockDimensionRoutine(oldDim,WorkingDim,i,workingindex)
						call tensor_transpose(inzp,outzp,WorkingDim,plan,rank)
					end if
				end do
			case (6)
				do i=1,TotalBlock
					if(inDa%getFlag(i))then
						call inDa%pointer(inlp,i)
						call outDa%pointer(outlp,i)
						call GetBlockDimensionRoutine(oldDim,WorkingDim,i,workingindex)
						call tensor_transpose(inlp,outlp,WorkingDim,plan,rank)
					end if
				end do
			case (7)
				do i=1,TotalBlock
					if(inDa%getFlag(i))then
						call inDa%pointer(inap,i)
						call outDa%pointer(outap,i)
						call GetBlockDimensionRoutine(oldDim,WorkingDim,i,workingindex)
						call tensor_transpose(inap,outap,WorkingDim,plan,rank)
					end if
				end do
			case default
				do i=1,TotalBlock
					if(inDa%getFlag(i))then
						call inDa%ClassPointer(inclp,i)
						call outDa%ClassPointer(outclp,i)
						call GetBlockDimensionRoutine(oldDim,WorkingDim,i,workingindex)
						call tensor_transpose_class(inclp,outclp,WorkingDim,plan,rank)
					end if
				end do
		end select
		return
	end subroutine
	subroutine permutationArrayBlock(outDa,inDa,dim_i,plan,rank)
		type(DataArray),target,intent(in)::inDa
		type(DataArray),target,intent(inout)::outDa
		integer,intent(in)::plan(:),rank,dim_i(:)
		integer,pointer::insi(:),outsi(:),inei(:),outei(:)
		if(outDa%getTotalBlock().gt.1)then
			call inDa%pointStarti(insi)
			call outDa%pointStarti(outsi)
			call tensor_transpose(insi,outsi,dim_i,plan,rank)

			call inDa%pointEndi(inei)
			call outDa%pointEndi(outei)
			call tensor_transpose(inei,outei,dim_i,plan,rank)
		end if
	end subroutine
	subroutine fix_fermi_sign(A)
		type(Tensor),intent(in)::A
		integer::i,itype
		integer,pointer::ip(:)
		real*4,pointer::sp(:)
		real*8,pointer::dp(:)
		complex*8,pointer::cp(:)
		complex*16,pointer::zp(:)
		integer,pointer::Blocksign(:)
		itype=A%getType()
		call A%Data%pointEndi(Blocksign)
		select case (itype)
			case (1)
				do i=1,size(Blocksign)
					if(Blocksign(i).lt.0)then
						Blocksign(i)=-Blocksign(i)
						call A%pointer(ip,i)
						ip=-ip
					end if
				end do
			case (2)
				do i=1,size(Blocksign)
					if(Blocksign(i).lt.0)then
						Blocksign(i)=-Blocksign(i)
						call A%pointer(sp,i)
						sp=-sp
					end if
				end do
			case (3)
				do i=1,size(Blocksign)
					if(Blocksign(i).lt.0)then
						Blocksign(i)=-Blocksign(i)
						call A%pointer(dp,i)
						dp=-dp
					end if
				end do
			case (4)
				do i=1,size(Blocksign)
					if(Blocksign(i).lt.0)then
						Blocksign(i)=-Blocksign(i)
						call A%pointer(cp,i)
						cp=-cp
					end if
				end do
			case (5)
				do i=1,size(Blocksign)
					if(Blocksign(i).lt.0)then
						Blocksign(i)=-Blocksign(i)
						call A%pointer(zp,i)
						zp=-zp
					end if
				end do
			case default
				call writemess(' ERROR class type in fermi permutation',-1)
				call writemess('Call A%getType()='+A%getType(),-1)
				call error_stop
		end select
		return
	end subroutine

	!*********************************************************
	!
	! permute forward
	!
	!*********************************************************

	function PermuteForWard_ith(A,ith)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		integer::rank,i,counter
		if(ith.eq.1)then
			Res=A
			return
		end if
		rank=A%getRank()
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
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_one_sign(Res,ith,legCheckith,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	subroutine PermuteForWard_ith_routine2(A,ith,Res)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		integer,intent(in)::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		integer::rank,i,counter
		if(ith.eq.1)then
			Res=A
			return
		end if
		rank=A%getRank()
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
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_one_sign(Res,ith,legCheckith,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine PermuteForWard_ith_routine(A,ith)
		class(Tensor),intent(inout)::A
		integer,intent(in)::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),legCheckith(:),dim_i(:)
		integer::rank,i,counter
		if(ith.eq.1)then
			return
		end if
		rank=A%getRank()
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
		TMPTensor1%Data=A%Data
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_one_sign(A,ith,legCheckith,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	function PermuteForWard_cha(A,cha)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		character(len=*),intent(in)::cha
		integer::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		character(len=characterlen),pointer::indices(:)
		integer::rank,i,counter
		if(index(cha,indexsymbol).eq.0)then
			TMPpermutation_index=A%getAllName(cha)
			if(.not.TMPpermutation_index%getFlag())then
				call writemess('ERROR in PermuteForWard_cha',-1)
				call error_stop
			end if
			call TMPpermutation_index%pointer(indices)
			call PermuteForWard_veccha_routine2(A,indices,Res)
			return
		end if
		rank=A%getRank()
		ith=A%FindOrder(cha)
		if(ith.eq.1)then
			Res=A
			return
		end if
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
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_one_sign(Res,ith,legCheckith,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	subroutine PermuteForWard_cha_routine(A,cha)
		class(Tensor),intent(inout)::A
		character(len=*),intent(in)::cha
		integer::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),legCheckith(:),dim_i(:)
		character(len=characterlen),pointer::indices(:)
		integer::rank,i,counter
		if(index(cha,indexsymbol).eq.0)then
			TMPpermutation_index=A%getAllName(cha)
			if(.not.TMPpermutation_index%getFlag())then
				call writemess('ERROR in PermuteForWard_cha',-1)
				call error_stop
			end if
			call TMPpermutation_index%pointer(indices)
			call PermuteForWard_veccha_routine(A,indices)
			return
		end if

		rank=A%getRank()
		ith=A%FindOrder(cha)
		if(ith.eq.1)then
			return
		end if
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
		TMPTensor1%Data=A%Data
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_one_sign(A,ith,legCheckith,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine


	function PermuteForWard_vec(A,vec)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::vec(:)
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),Wlegs(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(vec)
		if(Check_not_permute_forward_case1(vec,vec_len))then
			Res=A
			return
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		counter=0
		newOrder(1:vec_len)=vec
		do i=vec_len+1,rank
			counter=counter+1
			goon=if_In_the_array(counter,vec,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vec,vec_len)
			end do
			newOrder(i)=counter
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_sign(Res,vec,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	subroutine PermuteForWard_vec_routine2(A,vec,Res)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		integer,intent(in)::vec(:)
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),Wlegs(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(vec)
		if(Check_not_permute_forward_case1(vec,vec_len))then
			Res=A
			return
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		counter=0
		newOrder(1:vec_len)=vec
		do i=vec_len+1,rank
			counter=counter+1
			goon=if_In_the_array(counter,vec,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vec,vec_len)
			end do
			newOrder(i)=counter
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_sign(Res,vec,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine PermuteForWard_vec_routine(A,vec)
		class(Tensor),intent(inout)::A
		integer,intent(in)::vec(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::newOrder(:),Wlegs(:),legCheckith(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(vec)
		if(Check_not_permute_forward_case1(vec,vec_len))then
			return
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		counter=0
		newOrder(1:vec_len)=vec
		do i=vec_len+1,rank
			counter=counter+1
			goon=if_In_the_array(counter,vec,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vec,vec_len)
			end do
			newOrder(i)=counter
		end do
		TMPTensor1%Data=A%Data
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_sign(A,vec,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	function PermuteForWard_veccha(A,cha)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		character(len=*),intent(in)::cha(:)
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),Wlegs(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:),vec(:),checkvec(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(cha)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len+vec_len+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		call WorkingMemory%get_memory(vec,vec_len)
		call WorkingMemory%get_memory(checkvec,vec_len)

		do i=1,vec_len
			newOrder(i)=A%FindOrder(cha(i))
			vec(i)=newOrder(i)
		end do

		if(Check_not_permute_forward_case2(vec,vec_len,checkvec))then
			Res=A
			call WorkingMemory%free()
			return
		end if

		counter=0
		do i=vec_len+1,rank
			counter=counter+1
			goon=if_In_the_array(counter,newOrder(1:vec_len),vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,newOrder(1:vec_len),vec_len)
			end do
			newOrder(i)=counter
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_sign(Res,vec,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	subroutine PermuteForWard_veccha_routine2(A,cha,Res)
		class(Tensor),intent(in)::A
		character(len=*),intent(in)::cha(:)
		type(Tensor),intent(inout)::Res
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),Wlegs(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:),vec(:),checkvec(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(cha)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len+vec_len+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		call WorkingMemory%get_memory(vec,vec_len)
		call WorkingMemory%get_memory(checkvec,vec_len)

		do i=1,vec_len
			newOrder(i)=A%FindOrder(cha(i))
			vec(i)=newOrder(i)
		end do

		if(Check_not_permute_forward_case2(vec,vec_len,checkvec))then
			Res=A
			call WorkingMemory%free()
			return
		end if

		counter=0
		do i=vec_len+1,rank
			counter=counter+1
			goon=if_In_the_array(counter,newOrder(1:vec_len),vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,newOrder(1:vec_len),vec_len)
			end do
			newOrder(i)=counter
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_sign(Res,vec,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine PermuteForWard_veccha_routine(A,cha)
		class(Tensor),intent(inout)::A
		character(len=*),intent(in)::cha(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:),vec(:)
		integer,pointer::newOrder(:),Wlegs(:),legCheckith(:),checkvec(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(cha)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len+vec_len+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		call WorkingMemory%get_memory(vec,vec_len)
		call WorkingMemory%get_memory(checkvec,vec_len)
		do i=1,vec_len
			newOrder(i)=A%FindOrder(cha(i))
			vec(i)=newOrder(i)
		end do
		if(Check_not_permute_forward_case2(vec,vec_len,checkvec))then
			call WorkingMemory%free()
			return
		end if
		counter=0
		do i=vec_len+1,rank
			counter=counter+1
			goon=if_In_the_array(counter,newOrder(1:vec_len),vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,newOrder(1:vec_len),vec_len)
			end do
			newOrder(i)=counter
		end do
		TMPTensor1%Data=A%Data
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_sign(A,vec,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine PermuteForWard_cha_routine2(A,cha,Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
		character(len=*),intent(in)::cha
		integer::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		character(len=characterlen),pointer::indices(:)
		integer::rank,i,counter
		if(index(cha,indexsymbol).eq.0)then
			TMPpermutation_index=A%getAllName(cha)
			if(.not.TMPpermutation_index%getFlag())then
				call writemess('ERROR in PermuteForWard_cha',-1)
				call error_stop
			end if
			call TMPpermutation_index%pointer(indices)
			call PermuteForWard_veccha_routine2(A,indices,Res)
			return
		end if
		rank=A%getRank()
		ith=A%FindOrder(cha)
		if(ith.eq.1)then
			Res=A
			return
		end if
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
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_one_sign(Res,ith,legCheckith,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	!*********************************************************
	!
	! permute backward
	!
	!*********************************************************

	function PermuteBackWard_ith(A,ith)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		integer::rank,i,counter
		rank=A%getRank()
		if(ith.eq.rank)then
			Res=A
			return
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank-ith+1)
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
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_backward_one_sign(Res,ith,legCheckith,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	subroutine PermuteBackWard_ith_routine2(A,ith,Res)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		integer,intent(in)::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		integer::rank,i,counter
		rank=A%getRank()
		if(ith.eq.rank)then
			Res=A
			return
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank-ith+1)
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
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_backward_one_sign(Res,ith,legCheckith,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine
	
	subroutine PermuteBackWard_ith_routine(A,ith)
		class(Tensor),intent(inout)::A
		integer,intent(in)::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),legCheckith(:),dim_i(:)
		integer::rank,i,counter
		rank=A%getRank()
		if(ith.eq.rank)then
			return
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank-ith+1)
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
		TMPTensor1%Data=A%Data
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_backward_one_sign(A,ith,legCheckith,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	function PermuteBackWard_cha(A,cha)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		character(len=*),intent(in)::cha
		integer::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		character(len=characterlen),pointer::indices(:)
		integer::rank,i,counter
		if(index(cha,indexsymbol).eq.0)then
			TMPpermutation_index=A%getAllName(cha)
			if(.not.TMPpermutation_index%getFlag())then
				call writemess('ERROR in PermuteForWard_cha',-1)
				call error_stop
			end if
			call TMPpermutation_index%pointer(indices)
			call PermuteBackWard_veccha_routine2(A,indices,Res)
			return
		end if
		rank=A%getRank()
		ith=A%FindOrder(cha)
		if(ith.eq.rank)then
			Res=A
			return
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank-ith+1)
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
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_backward_one_sign(Res,ith,legCheckith,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	subroutine PermuteBackWard_cha_routine(A,cha)
		class(Tensor),intent(inout)::A
		character(len=*),intent(in)::cha
		integer::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),legCheckith(:),dim_i(:)
		character(len=characterlen),pointer::indices(:)
		integer::rank,i,counter
		if(index(cha,indexsymbol).eq.0)then
			TMPpermutation_index=A%getAllName(cha)
			if(.not.TMPpermutation_index%getFlag())then
				call writemess('ERROR in PermuteForWard_cha',-1)
				call error_stop
			end if
			call TMPpermutation_index%pointer(indices)
			call PermuteBackWard_veccha_routine(A,indices)
			return
		end if
		rank=A%getRank()
		ith=A%FindOrder(cha)
		if(ith.eq.rank)then
			return
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank-ith+1)
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
		TMPTensor1%Data=A%Data
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_backward_one_sign(A,ith,legCheckith,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	function PermuteBackWard_vec(A,vec)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::vec(:)
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),Wlegs(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(vec)
		if(Check_not_permute_backward_case1(vec,vec_len,rank))then
			Res=A
			return
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		counter=0
		newOrder(rank-vec_len+1:rank)=vec
		do i=1,rank-vec_len
			counter=counter+1
			goon=if_In_the_array(counter,vec,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vec,vec_len)
			end do
			if(counter.gt.rank)then
				call writemess('ERROR in PermuteBackWard',-1)
				call writemess(newOrder)
				call writemess(vec)
				call writemess("i="+i+',counter='+counter+',rank='+rank)
				call error_stop
			end if
			newOrder(i)=counter
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_backward_sign(Res,vec,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	subroutine PermuteBackWard_vec_routine(A,vec)
		class(Tensor),intent(inout)::A
		integer,intent(in)::vec(:)
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),Wlegs(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(vec)
		if(Check_not_permute_backward_case1(vec,vec_len,rank))then
			return
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		counter=0
		newOrder(rank-vec_len+1:rank)=vec
		do i=1,rank-vec_len
			counter=counter+1
			goon=if_In_the_array(counter,vec,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vec,vec_len)
			end do
			if(counter.gt.rank)then
				call writemess('ERROR in PermuteBackWard',-1)
				call writemess(newOrder)
				call writemess(vec)
				call writemess("i="+i+',counter='+counter+',rank='+rank)
				call error_stop
			end if
			newOrder(i)=counter
		end do
		TMPTensor1%Data=A%Data
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_backward_sign(A,vec,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine PermuteBackWard_vec_routine2(A,vec,Res)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		integer,intent(in)::vec(:)
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),Wlegs(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(vec)
		if(Check_not_permute_backward_case1(vec,vec_len,rank))then
			Res=A
			return
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		counter=0
		newOrder(rank-vec_len+1:rank)=vec
		do i=1,rank-vec_len
			counter=counter+1
			goon=if_In_the_array(counter,vec,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vec,vec_len)
			end do
			if(counter.gt.rank)then
				call writemess('ERROR in PermuteBackWard',-1)
				call writemess(newOrder)
				call writemess(vec)
				call writemess("i="+i+',counter='+counter+',rank='+rank)
				call error_stop
			end if
			newOrder(i)=counter
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_backward_sign(Res,vec,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	function PermuteBackWard_veccha(A,cha)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		character(len=*),intent(in)::cha(:)
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),Wlegs(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:),vec(:),veccheck(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(cha)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len+vec_len+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		call WorkingMemory%get_memory(vec,vec_len)
		call WorkingMemory%get_memory(veccheck,vec_len)

		counter=0
		do i=1,vec_len
			vec(i)=A%FindOrder(cha(i))
			newOrder(rank-vec_len+i)=vec(i)
		end do
		if(Check_not_permute_backward_case2(vec,vec_len,rank,veccheck))then
			Res=A
			call WorkingMemory%free()
			return
		end if
		do i=1,rank-vec_len
			counter=counter+1
			goon=if_In_the_array(counter,vec,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vec,vec_len)
			end do
			if(counter.gt.rank)then
				call writemess('ERROR in PermuteBackWard',-1)
				call error_stop
			end if
			newOrder(i)=counter
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_backward_sign(Res,vec,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	subroutine PermuteBackWard_veccha_routine2(A,cha,Res)
		class(Tensor),intent(in)::A
		character(len=*),intent(in)::cha(:)
		type(Tensor),intent(inout)::Res
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),Wlegs(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:),vec(:),veccheck(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(cha)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len+vec_len+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		call WorkingMemory%get_memory(vec,vec_len)
		call WorkingMemory%get_memory(veccheck,vec_len)

		counter=0
		do i=1,vec_len
			vec(i)=A%FindOrder(cha(i))
			newOrder(rank-vec_len+i)=vec(i)
		end do
		if(Check_not_permute_backward_case2(vec,vec_len,rank,veccheck))then
			Res=A
			call WorkingMemory%free()
			return
		end if
		do i=1,rank-vec_len
			counter=counter+1
			goon=if_In_the_array(counter,vec,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vec,vec_len)
			end do
			if(counter.gt.rank)then
				call writemess('ERROR in PermuteBackWard',-1)
				call error_stop
			end if
			newOrder(i)=counter
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_backward_sign(Res,vec,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine


	subroutine PermuteBackWard_veccha_routine(A,cha)
		class(Tensor),intent(inout)::A
		character(len=*),intent(in)::cha(:)
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),Wlegs(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:),vec(:),veccheck(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(cha)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len+vec_len+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		call WorkingMemory%get_memory(vec,vec_len)
		call WorkingMemory%get_memory(veccheck,vec_len)


		counter=0
		do i=1,vec_len
			vec(i)=A%FindOrder(cha(i))
			newOrder(rank-vec_len+i)=vec(i)
		end do
		if(Check_not_permute_backward_case2(vec,vec_len,rank,veccheck))then
			call WorkingMemory%free()
			return
		end if
		do i=1,rank-vec_len
			counter=counter+1
			goon=if_In_the_array(counter,vec,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vec,vec_len)
			end do
			if(counter.gt.rank)then
				call writemess('ERROR in PermuteBackWard',-1)
				call error_stop
			end if
			newOrder(i)=counter
		end do
		TMPTensor1%Data=A%Data
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_backward_sign(A,vec,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine PermuteBackWard_cha_routine2(A,cha,Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
		character(len=*),intent(in)::cha
		integer::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:)
		integer,pointer::ei(:),legCheckith(:),dim_i(:)
		character(len=characterlen),pointer::indices(:)
		integer::rank,i,counter
		if(index(cha,indexsymbol).eq.0)then
			TMPpermutation_index=A%getAllName(cha)
			if(.not.TMPpermutation_index%getFlag())then
				call writemess('ERROR in PermuteForWard_cha',-1)
				call error_stop
			end if
			call TMPpermutation_index%pointer(indices)
			call PermuteBackWard_veccha_routine2(A,indices,Res)
			return
		end if
		rank=A%getRank()
		ith=A%FindOrder(cha)
		if(ith.eq.rank)then
			Res=A
			return
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank-ith+1)
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
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_backward_one_sign(Res,ith,legCheckith,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine


	!*********************************************************
	!
	! permutation
	!
	!*********************************************************

	function permutation_vec(A,newOrder)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::newOrder(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::Wlegs(:),legCheckith(:)
		integer::rank,i,counter,vec_len
		rank=A%getRank()
		vec_len=size(newOrder)
		if(rank.ne.vec_len)then
			call writemess('ERROR in permutation_cha_routine',-1)
			call writemess('rank='+rank,-1)
			call writemess('size(newOrder)='+size(newOrder),-1)
			call error_stop
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_sign(Res,newOrder,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	function permutation_cha(A,NameOrder)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		character(len=*),intent(in)::NameOrder(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::newOrder(:),Wlegs(:),legCheckith(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(NameOrder)
		if(rank.ne.vec_len)then
			call writemess('ERROR in permutation_cha_routine',-1)
			call writemess('rank='+rank,-1)
			call writemess('size(order)='+vec_len,-1)
			call error_stop
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		do i=1,rank
			newOrder(i)=A%FindOrder(NameOrder(i))
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_sign(Res,newOrder,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	subroutine permutation_vec_routine(A,newOrder)
		class(Tensor),intent(inout)::A
		integer,intent(in)::newOrder(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::Wlegs(:),legCheckith(:)
		integer::rank,i,counter,vec_len
		TMPTensor1%Data=A%Data
		rank=A%getRank()
		vec_len=size(newOrder)
		if(rank.ne.vec_len)then
			call writemess('ERROR in permutation_cha_routine',-1)
			call writemess('rank='+rank,-1)
			call writemess('size(order)='+vec_len,-1)
			call error_stop
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_sign(A,newOrder,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine permutation_cha_routine(A,NameOrder)
		class(Tensor),intent(inout)::A
		character(len=*),intent(in)::NameOrder(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::newOrder(:),Wlegs(:),legCheckith(:)
		integer::rank,i,counter,vec_len
		logical::goon
		TMPTensor1%Data=A%Data
		rank=A%getRank()
		vec_len=size(NameOrder)
		if(rank.ne.vec_len)then
			call writemess('ERROR in permutation_cha_routine',-1)
			call writemess('rank='+rank,-1)
			call writemess('size(order)='+vec_len,-1)
			call error_stop
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		do i=1,rank
			newOrder(i)=A%FindOrder(NameOrder(i))
		end do
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_sign(A,newOrder,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine permutation_vec_routine2(A,newOrder,Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
		integer,intent(in)::newOrder(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::Wlegs(:),legCheckith(:)
		integer::rank,i,counter,vec_len
		rank=A%getRank()
		vec_len=size(newOrder)
		if(rank.ne.vec_len)then
			call writemess('ERROR in permutation_cha_routine',-1)
			call writemess('rank='+rank,-1)
			call writemess('size(order)='+vec_len,-1)
			call error_stop
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_sign(Res,newOrder,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine permutation_cha_routine2(A,NameOrder,Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
		character(len=*),intent(in)::NameOrder(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::newOrder(:),Wlegs(:),legCheckith(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(NameOrder)
		if(rank.ne.vec_len)then
			call writemess('ERROR in permutation_cha_routine',-1)
			call writemess('rank='+rank,-1)
			call writemess('size(order)='+vec_len,-1)
			call error_stop
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		do i=1,rank
			newOrder(i)=A%FindOrder(NameOrder(i))
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call ArrayData_forard_sign(Res,newOrder,legCheckith,workingindex,Wlegs,vec_len)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine



	!*********************************************************
	!
	! permute forward DO NOT Fix the fermi sign, use in contract
	!
	!*********************************************************
	
	function NotfermiPermuteForWard_ith(A,ith)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:)
		integer,pointer::ei(:),dim_i(:)
		integer::rank,i,counter
		rank=A%getRank()
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
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
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	subroutine NotfermiPermuteForWard_ith_routine(A,ith)
		class(Tensor),intent(inout)::A
		integer,intent(in)::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),dim_i(:)
		integer::rank,i,counter
		rank=A%getRank()
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		counter=0
		newOrder(1)=ith
		do i=2,rank
			counter=counter+1
			if(counter.eq.ith)then
				counter=counter+1
			end if
			newOrder(i)=counter
		end do
		TMPTensor1%Data=A%Data
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	function NotfermiPermuteForWard_vec(A,vec)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::vec(:)
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),dim_i(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(vec)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		counter=0
		newOrder(1:vec_len)=vec
		do i=vec_len+1,rank
			counter=counter+1
			goon=if_In_the_array(counter,vec,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vec,vec_len)
			end do
			newOrder(i)=counter
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	subroutine NotfermiPermuteForWard_vec_routine(A,vec)
		class(Tensor),intent(inout)::A
		integer,intent(in)::vec(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::newOrder(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(vec)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		counter=0
		newOrder(1:vec_len)=vec
		do i=vec_len+1,rank
			counter=counter+1
			goon=if_In_the_array(counter,vec,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vec,vec_len)
			end do
			newOrder(i)=counter
		end do
		TMPTensor1%Data=A%Data
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	!*********************************************************
	!
	! permute backward DO NOT Fix the fermi sign, use in contract
	!
	!*********************************************************

	function NotfermiPermuteBackWard_ith(A,ith)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:)
		integer,pointer::dim_i(:)
		integer::rank,i,counter
		rank=A%getRank()
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
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
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function
	
	subroutine NotfermiPermuteBackWard_ith_routine(A,ith)
		class(Tensor),intent(inout)::A
		integer,intent(in)::ith
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:),dim_i(:)
		integer::rank,i,counter
		rank=A%getRank()
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
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
		TMPTensor1%Data=A%Data
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	function NotfermiPermuteBackWard_vec(A,vec)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::vec(:)
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:)
		integer,pointer::dim_i(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(vec)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		counter=0
		newOrder(rank-vec_len+1:rank)=vec
		do i=1,rank-vec_len
			counter=counter+1
			goon=if_In_the_array(counter,vec,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vec,vec_len)
			end do
			if(counter.gt.rank)then
				call writemess('ERROR in PermuteBackWard',-1)
				call writemess(newOrder)
				call writemess(vec)
				call writemess("i="+i+',counter='+counter+',rank='+rank)
				call error_stop
			end if
			newOrder(i)=counter
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	subroutine NotfermiPermuteBackWard_vec_routine(A,vec)
		class(Tensor),intent(inout)::A
		integer,intent(in)::vec(:)
		integer,pointer::workingindex(:),WorkingDim(:),newOrder(:)
		integer,pointer::dim_i(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(vec)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		counter=0
		newOrder(rank-vec_len+1:rank)=vec
		do i=1,rank-vec_len
			counter=counter+1
			goon=if_In_the_array(counter,vec,vec_len)
			do while(goon)
				counter=counter+1
				goon=if_In_the_array(counter,vec,vec_len)
			end do
			if(counter.gt.rank)then
				call writemess('ERROR in PermuteBackWard',-1)
				call writemess(newOrder)
				call writemess(vec)
				call writemess("i="+i+',counter='+counter+',rank='+rank)
				call error_stop
			end if
			newOrder(i)=counter
		end do
		TMPTensor1%Data=A%Data
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	!*********************************************************
	!
	! permutation DO NOT Fix the fermi sign
	!
	!*********************************************************

	function Notfermipermutation_vec(A,newOrder)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		integer,intent(in)::newOrder(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::Wlegs(:),legCheckith(:)
		integer::rank,i,counter,vec_len
		rank=A%getRank()
		vec_len=size(newOrder)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	function Notfermipermutation_cha(A,NameOrder)Result(Res)
		type(Tensor)::Res
		type(Tensor),intent(in)::A
		character(len=*),intent(in)::NameOrder(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::newOrder(:),Wlegs(:),legCheckith(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(NameOrder)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		do i=1,rank
			newOrder(i)=A%FindOrder(NameOrder(i))
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end function

	subroutine Notfermipermutation_vec_routine(A,newOrder)
		class(Tensor),intent(inout)::A
		integer,intent(in)::newOrder(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::Wlegs(:),legCheckith(:)
		integer::rank,i,counter,vec_len
		TMPTensor1%Data=A%Data
		rank=A%getRank()
		vec_len=size(newOrder)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine Notfermipermutation_cha_routine(A,NameOrder)
		class(Tensor),intent(inout)::A
		character(len=*),intent(in)::NameOrder(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::newOrder(:),Wlegs(:),legCheckith(:)
		integer::rank,i,counter,vec_len
		logical::goon
		TMPTensor1%Data=A%Data
		rank=A%getRank()
		vec_len=size(NameOrder)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		do i=1,rank
			newOrder(i)=A%FindOrder(NameOrder(i))
		end do
		call permutationArrayData(A%Data,TMPTensor1%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(A%Data,TMPTensor1%Data,dim_i,newOrder,rank)
		call permutationDimension(A%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine Notfermipermutation_vec_routine2(A,newOrder,Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
		integer,intent(in)::newOrder(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::Wlegs(:),legCheckith(:)
		integer::rank,i,counter,vec_len
		rank=A%getRank()
		vec_len=size(newOrder)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine

	subroutine Notfermipermutation_cha_routine2(A,NameOrder,Res)
		type(Tensor)::Res
		class(Tensor),intent(in)::A
		character(len=*),intent(in)::NameOrder(:)
		integer,pointer::workingindex(:),WorkingDim(:),dim_i(:)
		integer,pointer::newOrder(:),Wlegs(:),legCheckith(:)
		integer::rank,i,counter,vec_len
		logical::goon
		rank=A%getRank()
		vec_len=size(NameOrder)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank+rank+rank+rank+vec_len)
		call WorkingMemory%get_memory(WorkingDim,rank)
		call WorkingMemory%get_memory(workingindex,rank)
		call WorkingMemory%get_memory(newOrder,rank)
		call WorkingMemory%get_memory(legCheckith,rank)
		call WorkingMemory%get_memory(Wlegs,vec_len)
		do i=1,rank
			newOrder(i)=A%FindOrder(NameOrder(i))
		end do
		Res%Dimension=A%Dimension
		call permutationArrayData(Res%Data,A%Data,newOrder,A%Dimension,WorkingDim,workingindex)
		call A%pointDim(dim_i)
		call permutationArrayBlock(Res%Data,A%Data,dim_i,newOrder,rank)
		call permutationDimension(Res%Dimension,newOrder)
		call WorkingMemory%free()
		return
	end subroutine