	subroutine eyeQNDATAName(A,QN,Res,legQN)
		type(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		type(QuanNum),intent(in)::QN
		logical,intent(in),optional::legQN
		EYEDATAType,pointer::Rp(:,:),Ap(:)
		integer,pointer::ARROW(:),dim(:),deg(:)
		real*4,pointer::sp(:)
		integer::nonSymDimi,i,ii,jj,ith
		Res%Dimension=[QN,QN]
		call Res%pointArrow(arrow)
		if(present(legQN).and.(.not.legQN))then
			ith=1
		else
			ith=2
		end if
		arrow(ith)=-arrow(ith)
		call reverseDimensionRule(Res%Dimension,ith)
		call Res%pointDim(dim)
		call Res%pointdeg(deg,1)
		nonSymDimi=sum(deg)
		call Res%Data%allocateDataArrayMomery(product(dim),nonSymDimi*nonSymDimi,A%getType())
		do i=1,dim(1)
			if(A%Data%getFlag(i))then
				call Res%setBlockMomery([i,i])
				call Res%pointer(Rp,[i,i])
				call A%Data%pointer(Ap,i)
				do jj=1,deg(i)
					do ii=1,deg(i)
						if(ii.eq.jj)then
							Rp(ii,jj)=Ap(ii)
						else
							Rp(ii,jj)=0
						end if
					end do
				end do
			end if
		end do
		return
	end subroutine

	subroutine eyeNDATAName(A,N,Res)
		type(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		integer,intent(in)::N
		EYEDATAType,pointer::Rp(:,:),Ap(:)
		integer::i,j
		call Res%allocate([N,N],A%getType())
		call Res%pointer(Rp)
		call A%Data%pointer(Ap)
		do j=1,N
			do i=1,N
				if(i.eq.j)then
					Rp(i,i)=Ap(i)
				else
					Rp(i,j)=0
				end if
			end do
		end do
		return
	end subroutine

	subroutine eyeDimDATAName(A,Dimen,Res)
		type(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		type(dimension),intent(in)::Dimen
		EYEDATAType,pointer::Rp(:,:),Ap(:)
		integer,pointer::ARROW(:),dim(:),deg(:)
		real*4,pointer::sp(:)
		integer::nonSymDimi,i,ii,jj,j
		character(len=len_of_name),pointer::names(:)
		if(Dimen%getRank().ne.1)then
			call writemess('ERROR in eye tensor',-1)
			call error_stop
		end if
		call pasteDimension(Res%Dimension,Dimen,1,1,Dimen,1,1)
		if(Res%getNameFlag())then
			call Res%pointName(Names)
			Names(2)='NoName.2'
		end if
		if(Res%getSymmetryFlag())then
			call reverseDimensionRule(Res%Dimension,2)
			call Res%pointDim(dim)
			call Res%pointdeg(deg,1)
			nonSymDimi=sum(deg)
			if(Res%getFermiFlag())then
				call Res%pointArrow(arrow)
				arrow(2)=-arrow(2)
			end if
			call Res%Data%allocateDataArrayMomery(product(dim),nonSymDimi*nonSymDimi,A%getType())
			do i=1,dim(1)
				if(A%Data%getFlag(i))then
					call Res%setBlockMomery([i,i])
					call Res%pointer(Rp,[i,i])
					call A%Data%pointer(Ap,i)
					do jj=1,deg(i)
						do ii=1,deg(i)
							if(ii.eq.jj)then
								Rp(ii,jj)=Ap(ii)
							else
								Rp(ii,jj)=0
							end if
						end do
					end do
				end if
			end do
			return
		end if
		call Res%pointDim(dim)
		call Res%Data%allocate(product(dim),A%getType())
		call Res%pointer(Rp)
		call A%Data%pointer(Ap)
		do j=1,dim(2)
			do i=1,dim(1)
				if(i.eq.j)then
					Rp(i,i)=Ap(i)
				else
					Rp(i,j)=0
				end if
			end do
		end do
		return
	end subroutine

