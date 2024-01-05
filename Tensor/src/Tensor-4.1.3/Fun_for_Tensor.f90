module Function_for_Tensor
	use Tensor_tools
	use Tools
	use U1_Tools
	use Quantumnumber_tools
	use Dimension_Tools
	use Basic_Tools
	use Pointer_Tools
	use Parity_Tools
	implicit none

	interface antiSubTensor
		module procedure antiSubTensorSubroutine1
		module procedure antiSubTensorSubroutine2
	end interface

	!antiSubTensor
		! use the SubTensor to split a symTensor T into same Tensor t1,t2,t3
		! after modify the data of t1,t2,t3, use antiSubTensor to store the
		! data back to T
		!


contains

	subroutine AntisubSymTensor8(T,inTi,Legi,QNi,degi)
		type(Tensor),intent(inout) ::T
		type(Tensor),intent(inout)::inTi
		real*8,pointer::Rp3(:,:,:),TP3(:,:,:)
		real*8,pointer::Rp2(:,:),TP2(:,:)
		real*8,pointer::Tp(:),Rp(:)
		integer,intent(in)::Legi,QNi,degi
		type(QuanNum)::NewQN
		integer::i,j,ii,rank,Blocki
		integer::Dim1,Dim2,Dim3,Deg1,Deg2,Deg3
		integer::RDim1,RDim2,RDim3,RDeg1,RDeg2,RDeg3
		integer,pointer::dim(:),Rdim(:)
		integer,allocatable::indices(:)
		real*4::QN
		rank=T%GetRank()
		if(legi.gt.rank)then
			call writemess('ERROR in subSymTensor',-1)
			call writemess('legi='+legi,-1)
			call writemess('rank='+rank,-1)
			call error_Stop
		end if
		allocate(indices(rank))


		call inTi%pointDim(Rdim)
		call T%pointDim(Dim)
		if(QNi.gt.dim(legi))then
			call writemess('ERROR in subSymTensor',-1)
			call writemess('QNi='+QNi,-1)
			call writemess('dim(legi)='+dim(legi),-1)
			call error_Stop
		end if
		if(degi.gt.T%getDeg(legi,QNi))then
			call writemess('ERROR in subSymTensor',-1)
			call writemess('degi='+QNi,-1)
			call writemess('T%getDeg(legi,QNi)='+T%getDeg(legi,QNi),-1)
			call error_Stop
		end if

			if((legi.gt.1).and.(legi.lt.rank))then
				dim1=product(dim(1:legi-1))
				dim2=dim(legi)
				dim3=product(dim(legi+1:rank))

				Rdim1=dim1
				Rdim2=1
				Rdim3=dim3
				do j=1,dim3
					do i=1,dim1
						Blocki=addressToIndes([i,QNi,j],[dim1,dim2,dim3])
						if(T%getFlag(Blocki))then
							call T%pointer(Tp,Blocki)
							call IndesToaddress(dim,indices,Blocki)
							RDeg1=inTi%getBlockDim(1,legi-1,indices(1:legi-1))
							RDeg2=1
							RDeg3=inTi%getBlockDim(legi+1,rank,indices(legi+1:rank))
							Deg1=RDeg1
							Deg2=T%getDeg(legi,QNi)
							Deg3=RDeg3
							call inTi%pointer([Rdim1,Rdim2,Rdim3],Rp3,[Rdeg1,Rdeg2,Rdeg3],[i,1,j])
							if(size(Tp).ne.(Deg1*Deg2*Deg3))then
								call writemess('ERROR in subSymTensor',-1)
								call error_stop
							end if
							Tp3(1:Deg1,1:Deg2,1:Deg3)=>Tp
							Tp3(:,degi,:)=Rp3(:,1,:)
						end if
					end do
				end do
			else if(legi.eq.1)then
				dim1=dim(1)
				dim2=product(dim(2:rank))

				Rdim1=1
				Rdim2=dim2
				do j=1,dim2
					Blocki=addressToIndes([QNi,j],[dim1,dim2])
					if(T%getFlag(Blocki))then
						call T%pointer(Tp,Blocki)
						call IndesToaddress(dim,indices,Blocki)
						RDeg1=1
						RDeg2=inTi%getBlockDim(2,rank,indices(2:rank))
						Deg1=T%getDeg(1,QNi)
						Deg2=RDeg2
						call inTi%pointer([Rdim1,Rdim2],Rp2,[Rdeg1,Rdeg2],[1,j])
						if(size(Tp).ne.(Deg1*Deg2))then
							call writemess('ERROR in subSymTensor',-1)
							call error_stop
						end if
						Tp2(1:Deg1,1:Deg2)=>Tp
						Tp2(degi,:)=Rp2(1,:)
					end if
				end do
			else if(legi.eq.rank)then
				dim1=product(dim(1:rank-1))
				dim2=dim(rank)

				Rdim1=dim1
				Rdim2=1
				do i=1,dim1
					Blocki=addressToIndes([i,QNi],[dim1,dim2])
					if(T%getFlag(Blocki))then
						call T%pointer(Tp,Blocki)
						call IndesToaddress(dim,indices,Blocki)
						RDeg1=inTi%getBlockDim(1,rank-1,indices(1:rank-1))
						RDeg2=1
						Deg1=RDeg1
						Deg2=T%getDeg(rank,QNi)
						call inTi%pointer([Rdim1,Rdim2],Rp2,[Rdeg1,Rdeg2],[i,1])
						if(size(Tp).ne.(Deg1*Deg2))then
							call writemess('ERROR in subSymTensor',-1)
							call error_stop
						end if
						Tp2(1:Deg1,1:Deg2)=>Tp
						Tp2(:,degi)=Rp2(:,1)
					end if
				end do
			end if
		
		return
	end subroutine

	subroutine antiSubTensorSubroutine1(T,Tin,legName)
		type(Tensor),intent(inout)::T
		type(Tensor),intent(inout)::Tin(:)
		character(len=*),intent(in)::legName
		integer::i,ith,rank,k,l
		type(TEnsor)::AllName
		type(QuanNum)::QN
		call T%backward(T%getName(legName))
		QN=T%quantumnumber(T%getRank())
		ith=0
		AllName=T%getName()
		rank=T%getRank()
		do k=1,QN%getQNlength()
			do l=1,T%getDeg(rank,k)
				if(ith.gt.size(Tin))then
					call writemess('ERROR in antiSubTensorSubroutine')
					call error_stop
				end if
				ith=ith+1
				call Tin(ith)%permute(AllName%ai())
				call AntisubSymTensor8(T,Tin(ith),rank,k,l)
			end do
		end do
		return
	end subroutine

	subroutine antiSubTensorSubroutine2(T,Tin)
		type(Tensor),intent(inout)::T
		type(Tensor),intent(inout)::Tin(:)
		integer::i,j,ith,rank
		type(QuanNum)::QN
		QN=T%quantumnumber(T%getRank())
		rank=T%getRank()
		ith=0
		do i=1,QN%getQNlength()
			do j=1,T%getDeg(rank,i)
				ith=ith+1
				if(ith.gt.size(Tin))then
					call error_stop
				end if
				call AntisubSymTensor8(T,Tin(ith),rank,i,j)
			end do
		end do
		return
	end subroutine


end module
