	



	! T(:,ith,:)=P(ith)*T(:,ith,:) ,input Tensor
		!where P(ith)=+-1 , the parity of ith


	subroutine Contract_fix_Sgin(T,ithleg)!(block,lenofblock,dimen,Info)
		type(Tensor),intent(inout)::T
		integer,intent(in)::ithleg
		integer::i,rank,M,N,L,totalBlock
		integer,pointer::Blocksign(:),dim(:)
		
		if(.not.T%Data%getFlag())return

		rank=T%getRank()
		call T%Data%pointEndi(Blocksign)
		if(ithleg.eq.1)then
			M=T%dim(1)
			N=1
			do i=2,rank
				N=N*T%dim(i)
			end do
			call Contract_fix_Sgin_first_index(Blocksign,M,N,T%dimension)
			call fix_fermi_sign(T)
			return
		end if
		if(ithleg.eq.rank)then
			N=T%dim(rank)
			M=1
			do i=rank-1,1,-1
				M=M*T%dim(i)
			end do
			call Contract_fix_Sgin_Last_index(Blocksign,M,N,rank,T%dimension)
			call fix_fermi_sign(T)
			return
		end if
		call T%pointDim(dim)
		M=product(dim(1:(ithleg-1)))
		N=T%dim(ithleg)
		L=product(dim(ithleg+1:rank))
		!M=T%dim(1)
		!do i=2,ithleg-1
		!	M=M*T%dim(i)
		!end do
		!N=T%dim(ithleg)
		!L=T%dim(ithleg+1)
		!do i=ithleg+2,rank
		!	L=L*T%dim(i)
		!end do
		call Contract_fix_Sgin_ith_index(Blocksign,M,N,L,ithleg,T%dimension)
		call fix_fermi_sign(T)
		return
	end subroutine	


	subroutine Contract_fix_Sgin_first_index(ei,M,N,dimen)
		type(Dimension),intent(in)::dimen
		integer,intent(in)::M,N
		integer,intent(inout)::ei(M,N)
		integer::i,j,iQN
		do i=1,M
			call QaunNumParity(iQN,dimen,1,i)
			if(iQN.eq.(-1))then
				do j=1,N
					if(ei(i,j).ne.0)then
						ei(i,j)=-ei(i,j)
					end if
				end do
			end if
		end do
	end subroutine

	subroutine Contract_fix_Sgin_Last_index(ei,M,N,rank,dimen)
		type(Dimension),intent(in)::dimen
		integer,intent(in)::M,N,rank
		integer,intent(inout)::ei(M,N)
		integer::i,j,iQN
		do j=1,N
			call QaunNumParity(iQN,dimen,rank,j)
			if(iQN.eq.(-1))then
				do i=1,M
					if(ei(i,j).ne.0)then
						ei(i,j)=-ei(i,j)
					end if
				end do
			end if
		end do
	end subroutine

	subroutine Contract_fix_Sgin_ith_index(ei,M,N,L,ith,dimen)
		type(Dimension),intent(in)::dimen
		integer,intent(in)::M,N,L,ith
		integer,intent(inout)::ei(M,N,L)
		integer::i,j,k,iQN
		do j=1,N
			call QaunNumParity(iQN,dimen,ith,j)
			if(iQN.eq.(-1))then
				do k=1,L
					do i=1,M
						if(ei(i,j,k).ne.0)then
							ei(i,j,k)=-ei(i,j,k)
						end if
					end do
				end do
			end if
		end do
	end subroutine



	!Reverse_Fermi_Rule
		! [T1]---->---[T2]
		! change to
		! [T1]----<---[T2]
		!
		!modify the rule of last leg of T1 and the first leg of T2
		!
		!	

	subroutine Reverse_Fermi_Rule1(T1,T2)
		type(Tensor),intent(inout)::T1,T2
		integer::rank1
		rank1=T1%getRank()
		if(T1%getTotalBlock().gt.T2%getTotalBlock())then
			call Contract_fix_Sgin(T2,1)
		else
			call Contract_fix_Sgin(T1,rank1)
		end if
		call T1%setFermiArrow(rank1,-T1%getFermiArrow(rank1))
		call T2%setFermiArrow(1,-T2%getFermiArrow(1))
		return
	end subroutine
	
	
	subroutine Reverse_Fermi_Rule_specify2(T1,ith1,T2,ith2)
		integer,intent(in)::ith1,ith2
		type(Tensor),intent(inout)::T1,T2
		if(T1%getTotalBlock().gt.T2%getTotalBlock())then
			call Contract_fix_Sgin(T2,ith2)
		else
			call Contract_fix_Sgin(T1,ith1)
		end if
		call T1%setFermiArrow(ith1,-T1%getFermiArrow(ith1))
		call T2%setFermiArrow(ith2,-T2%getFermiArrow(ith2))
		return
	end subroutine
	subroutine Reverse_Fermi_Rule_specify3(T1,ith1)
		integer,intent(in)::ith1
		class(Tensor),intent(inout)::T1
		call Contract_fix_Sgin(T1,ith1)
		call T1%setFermiArrow(ith1,-T1%getFermiArrow(ith1))
		return
	end subroutine

	subroutine Reverse_Fermi_Rule_specify4(T1,name1,T2,name2)
		character(len=*),intent(in)::name1,name2
		type(Tensor),intent(inout)::T1,T2
		integer::ith1,ith2
		ith1=T1%FindOrder(name1)
		ith2=T2%FindOrder(name2)
		if(T1%getTotalBlock().gt.T2%getTotalBlock())then
			call Contract_fix_Sgin(T2,ith2)
		else
			call Contract_fix_Sgin(T1,ith1)
		end if
		call T1%setFermiArrow(ith1,-T1%getFermiArrow(ith1))
		call T2%setFermiArrow(ith2,-T2%getFermiArrow(ith2))
		return
	end subroutine

	subroutine Reverse_Fermi_Rule_specify5(T1,name1)
		character(len=*),intent(in)::name1
		class(Tensor),intent(inout)::T1
		integer::ith1
		ith1=T1%FindOrder(name1)
		call Contract_fix_Sgin(T1,ith1)
		call T1%setFermiArrow(ith1,-T1%getFermiArrow(ith1))
		return
	end subroutine

	subroutine Reverse_Fermi_Rule_specify6(T1,name1,T2,name2,T1orT2)
		character(len=*),intent(in)::name1,name2
		type(Tensor),intent(inout)::T1,T2
		logical,intent(in)::T1orT2
		integer::ith1,ith2
		ith1=T1%FindOrder(name1)
		ith2=T2%FindOrder(name2)
		if(T1orT2)then
			call Contract_fix_Sgin(T1,ith1)
		else
			call Contract_fix_Sgin(T2,ith2)
		end if
		call T1%setFermiArrow(ith1,-T1%getFermiArrow(ith1))
		call T2%setFermiArrow(ith2,-T2%getFermiArrow(ith2))
		return
	end subroutine

	subroutine Reverse_Fermi_Rule_specify7(T1,ith1,T2,ith2,T1orT2)
		integer,intent(in)::ith1,ith2
		type(Tensor),intent(inout)::T1,T2
		logical,intent(in)::T1orT2
		if(T1orT2)then
			call Contract_fix_Sgin(T1,ith1)
		else
			call Contract_fix_Sgin(T2,ith2)
		end if
		call T1%setFermiArrow(ith1,-T1%getFermiArrow(ith1))
		call T2%setFermiArrow(ith2,-T2%getFermiArrow(ith2))
		return
	end subroutine

