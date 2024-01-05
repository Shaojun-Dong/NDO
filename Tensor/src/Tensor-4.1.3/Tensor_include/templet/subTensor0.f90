	subroutine subSymTensorFuncName(Res,T,Legi,QNi,degi,keepQN_)
		type(Tensor),intent(inout)::Res
		type(Tensor),intent(in) ::T
		DATATYPE,pointer::Rp3(:,:,:),TP3(:,:,:)
		DATATYPE,pointer::Rp2(:,:),TP2(:,:)
		DATATYPE,pointer::Tp(:),Rp(:)
		integer,intent(in)::Legi,QNi,degi
		logical,optional,intent(in)::keepQN_
		logical::keepQN
		type(QuanNum)::NewQN
		integer::i,j,ii,rank,Blocki
		integer::Dim1,Dim2,Dim3,Deg1,Deg2,Deg3
		integer::RDim1,RDim2,RDim3,RDeg1,RDeg2,RDeg3
		integer,pointer::dim(:),indices(:),Rdim(:)
		real*4::QN
		if(present(keepQN_))then
			keepQN=keepQN_
		else
			keepQN=.false.
		end if
		rank=T%GetRank()
		if(legi.gt.rank)then
			call writemess('ERROR in subSymTensor',-1)
			call writemess('legi='+legi,-1)
			call writemess('rank='+rank,-1)
			call error_Stop
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank)
		call WorkingMemory%get_memory(indices,rank)

		if(keepQN)then
			QN=T%GetQN(legi,QNi)
			call NewQN%setRule(1)
			if(T%getFermiFlag())then
				call NewQN%setFermiArrow(T%getFermiArrow(legi))
			end if
			call NewQN%setQN([QN])
			call NewQN%setDeg(1,1)
			if((legi.gt.1).and.(legi.lt.rank))then
				call pasteDimension(Res%Dimension,T%dimension,1,legi-1,NewQN,T%dimension,legi+1,rank)
			else if(legi.eq.1)then
				call pasteDimension(Res%Dimension,NewQN,T%dimension,2,rank)
			else if(legi.eq.rank)then
				call pasteDimension(Res%Dimension,T%dimension,1,rank-1,NewQN)
			else
				call writemess('ERROR in subSymTensor',-1)
				call error_Stop
			end if
			if(T%getNameFlag())then
				call Res%setName(legi,T%getName(legi))
			end if
		else
			if((legi.gt.1).and.(legi.lt.rank))then
				call pasteDimension(Res%Dimension,T%dimension,1,legi-1,T%dimension,legi+1,rank)
			else if(legi.eq.1)then
				call getsubDimension(T%Dimension,[2,rank],Res%Dimension)
			else if(legi.eq.rank)then
				call getsubDimension(T%Dimension,[1,rank-1],Res%Dimension)
			else
				call writemess('ERROR in subSymTensor',-1)
				call error_Stop
			end if
			
		end if
		call Res%pointDim(Rdim)
		call Res%data%allocateDataArrayMomery(product(Rdim),T%getTotalData(),T%getType())
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
		if(keepQN)then
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
							RDeg1=Res%getBlockDim(1,legi-1,indices(1:legi-1))
							RDeg2=1
							RDeg3=Res%getBlockDim(legi+1,rank,indices(legi+1:rank))
							Deg1=RDeg1
							Deg2=T%getDeg(legi,QNi)
							Deg3=RDeg3
							call Res%setBlockMomery([Rdim1,Rdim2,Rdim3],[i,1,j],Deg1*Deg3)
							call Res%pointer([Rdim1,Rdim2,Rdim3],Rp3,[Rdeg1,Rdeg2,Rdeg3],[i,1,j])
							if(size(Tp).ne.(Deg1*Deg2*Deg3))then
								call writemess('ERROR in subSymTensor',-1)
								call error_stop
							end if
							Tp3(1:Deg1,1:Deg2,1:Deg3)=>Tp
							Rp3(:,1,:)=Tp3(:,degi,:)
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
						RDeg2=Res%getBlockDim(2,rank,indices(2:rank))
						Deg1=T%getDeg(1,QNi)
						Deg2=RDeg2
						call Res%setBlockMomery([Rdim1,Rdim2],[1,j],RDeg2)
						call Res%pointer([Rdim1,Rdim2],Rp2,[Rdeg1,Rdeg2],[1,j])
						if(size(Tp).ne.(Deg1*Deg2))then
							call writemess('ERROR in subSymTensor',-1)
							call error_stop
						end if
						Tp2(1:Deg1,1:Deg2)=>Tp
						Rp2(1,:)=Tp2(degi,:)
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
						RDeg1=Res%getBlockDim(1,rank-1,indices(1:rank-1))
						RDeg2=1
						Deg1=RDeg1
						Deg2=T%getDeg(rank,QNi)
						call Res%setBlockMomery([Rdim1,Rdim2],[i,1],RDeg1)
						call Res%pointer([Rdim1,Rdim2],Rp2,[Rdeg1,Rdeg2],[i,1])
						if(size(Tp).ne.(Deg1*Deg2))then
							call writemess('ERROR in subSymTensor',-1)
							call error_stop
						end if
						Tp2(1:Deg1,1:Deg2)=>Tp
						Rp2(:,1)=Tp2(:,degi)
					end if
				end do
			end if
		
		else
			if((legi.gt.1).and.(legi.lt.rank))then
				dim1=product(dim(1:legi-1))
				dim2=dim(legi)
				dim3=product(dim(legi+1:rank))

				Rdim1=dim1
				Rdim3=dim3
				do j=1,dim3
					do i=1,dim1
						Blocki=addressToIndes([i,QNi,j],[dim1,dim2,dim3])
						if(T%getFlag(Blocki))then
							call T%pointer(Tp,Blocki)
							call IndesToaddress(dim,indices,Blocki)
							RDeg1=Res%getBlockDim(1,legi-1,indices(1:legi-1))
							RDeg3=Res%getBlockDim(legi+1,rank,indices(legi+1:rank))
							Deg1=RDeg1
							Deg2=T%getDeg(legi,QNi)
							Deg3=RDeg3
							call Res%setBlockMomery([Rdim1,Rdim3],[i,1,j],Deg1*Deg3)
							call Res%pointer([Rdim1,Rdim3],Rp2,[Rdeg1,Rdeg3],[i,j])
							if(size(Tp).ne.(Deg1*Deg2*Deg3))then
								call writemess('ERROR in subSymTensor',-1)
								call error_stop
							end if
							Tp3(1:Deg1,1:Deg2,1:Deg3)=>Tp
							Rp2(:,:)=Tp3(:,degi,:)
						end if
					end do
				end do
			else if(legi.eq.1)then
				dim1=dim(1)
				dim2=product(dim(2:rank))

				Rdim2=dim2
				do j=1,dim2
					Blocki=addressToIndes([QNi,j],[dim1,dim2])
					if(T%getFlag(Blocki))then
						call T%pointer(Tp,Blocki)
						call IndesToaddress(dim,indices,Blocki)
						RDeg2=Res%getBlockDim(2,rank,indices(2:rank))
						Deg1=T%getDeg(1,QNi)
						Deg2=RDeg2
						call Res%setBlockMomery(j,RDeg2)
						call Res%pointer(Rp,j)
						if(size(Tp).ne.(Deg1*Deg2))then
							call writemess('ERROR in subSymTensor',-1)
							call error_stop
						end if
						Tp2(1:Deg1,1:Deg2)=>Tp
						Rp(:)=Tp2(degi,:)
					end if
				end do
			else if(legi.eq.rank)then
				dim1=product(dim(1:rank-1))
				dim2=dim(rank)

				Rdim1=dim1
				do i=1,dim1
					Blocki=addressToIndes([i,QNi],[dim1,dim2])
					if(T%getFlag(Blocki))then
						call T%pointer(Tp,Blocki)
						call IndesToaddress(dim,indices,Blocki)
						RDeg1=Res%getBlockDim(1,rank-1,indices(1:rank-1))
						Deg1=RDeg1
						Deg2=T%getDeg(rank,QNi)
						call Res%setBlockMomery(i,RDeg1)
						call Res%pointer(Rp,i)
						if(size(Tp).ne.(Deg1*Deg2))then
							call writemess('ERROR in subSymTensor',-1)
							call error_stop
						end if
						Tp2(1:Deg1,1:Deg2)=>Tp
						Rp(:)=Tp2(:,degi)
					end if
				end do
			end if
		
		end if
		call WorkingMemory%free()
		return
	end subroutine

	subroutine subTensorFuncName(Res,T,Legi,dimi,keepQN_)
		type(Tensor),intent(inout)::Res
		type(Tensor),intent(in) ::T
		DATATYPE,pointer::Rp3(:,:,:),TP3(:,:,:)
		DATATYPE,pointer::Rp2(:,:),TP2(:,:)
		DATATYPE,pointer::Tp(:),Rp(:)
		integer,intent(in)::Legi,dimi
		logical,optional,intent(in)::keepQN_
		logical::keepQN
		type(Dimension)::Newdimi
		integer::i,j,ii,rank,Blocki
		integer::Dim1,Dim2,Dim3,Deg1,Deg2,Deg3
		integer::RDim1,RDim2,RDim3,RDeg1,RDeg2,RDeg3
		integer,pointer::dim(:),indices(:),Rdim(:)
		if(present(keepQN_))then
			keepQN=keepQN_
		else
			keepQN=.false.
		end if
		rank=T%GetRank()
		if(legi.gt.rank)then
			call writemess('ERROR in subTensor',-1)
			call writemess('legi='+legi,-1)
			call writemess('rank='+rank,-1)
			call error_Stop
		end if
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank)
		call WorkingMemory%get_memory(indices,rank)

		if(keepQN)then
			
			if((legi.gt.1).and.(legi.lt.rank))then
				Res%Dimension=[1]
				call pasteDimension(Newdimi,T%dimension,1,legi-1,Res%Dimension,1,1)
				call pasteDimension(Res%Dimension,Newdimi,1,legi,T%dimension,legi+1,rank)
			else if(legi.eq.1)then
				Newdimi=[1]
				call pasteDimension(Res%Dimension,Newdimi,1,1,T%dimension,2,rank)
			else if(legi.eq.rank)then
				Newdimi=[1]
				call pasteDimension(Res%Dimension,T%dimension,1,rank-1,Newdimi,1,1)
			else
				call writemess('ERROR in subSymTensor',-1)
				call error_Stop
			end if
			if(T%getNameFlag())then
				call Res%setName(legi,T%getName(legi))
			end if
		else
			if((legi.gt.1).and.(legi.lt.rank))then
				call pasteDimension(Res%Dimension,T%dimension,1,legi-1,T%dimension,legi+1,rank)
			else if(legi.eq.1)then
				call getsubDimension(T%Dimension,[2,rank],Res%Dimension)
			else if(legi.eq.rank)then
				call getsubDimension(T%Dimension,[1,rank-1],Res%Dimension)
			else
				call writemess('ERROR in subSymTensor',-1)
				call error_Stop
			end if
			
		end if
		call Res%pointDim(Rdim)
		call Res%data%allocate(product(Rdim),T%getType())
		call T%pointDim(Dim)
		if(dimi.gt.dim(legi))then
			call writemess('ERROR in subTensor',-1)
			call writemess('dimi='+dimi,-1)
			call writemess('dim(legi)='+dim(legi),-1)
			call error_Stop
		end if
		if(keepQN)then
			if((legi.gt.1).and.(legi.lt.rank))then
				dim1=product(dim(1:legi-1))
				dim2=dim(legi)
				dim3=product(dim(legi+1:rank))
				call T%pointer(Tp)
				Tp3(1:dim1,1:dim2,1:dim3)=>Tp

				Rdim1=dim1
				Rdim2=1
				Rdim3=dim3
				call Res%pointer(Rp)
				Rp3(1:dim1,1:1,1:dim3)=>Rp
				Rp3(:,1,:)=Tp3(:,dimi,:)
			else if(legi.eq.1)then
				dim1=dim(1)
				dim2=product(dim(2:rank))
				call T%pointer(Tp)
				Tp2(1:dim1,1:dim2)=>Tp

				Rdim1=1
				Rdim2=dim2
				call Res%pointer(Rp)
				Rp2(1:1,1:dim2)=>Rp
				Rp2(1,:)=Tp2(dimi,:)
			else if(legi.eq.rank)then
				dim1=product(dim(1:rank-1))
				dim2=dim(rank)
				call T%pointer(Tp)
				Tp2(1:dim1,1:dim2)=>Tp

				Rdim1=dim1
				Rdim2=1
				call Res%pointer(Rp)
				Rp2(1:dim1,1:1)=>Rp
				Rp2(:,1)=Tp2(:,dimi)
				
			end if
		
		else
			if((legi.gt.1).and.(legi.lt.rank))then
				dim1=product(dim(1:legi-1))
				dim2=dim(legi)
				dim3=product(dim(legi+1:rank))
				call T%pointer(Tp)
				Tp3(1:dim1,1:dim2,1:dim3)=>Tp

				Rdim1=dim1
				Rdim3=dim3
				call Res%pointer(Rp)
				Rp2(1:dim1,1:dim3)=>Rp
				Rp2=Tp3(:,dimi,:)
			else if(legi.eq.1)then
				dim1=dim(1)
				dim2=product(dim(2:rank))
				call T%pointer(Tp)
				Tp2(1:dim1,1:dim2)=>Tp

				Rdim2=dim2
				call Res%pointer(Rp)
				Rp=Tp2(dimi,:)
			else if(legi.eq.rank)then
				dim1=product(dim(1:rank-1))
				dim2=dim(rank)
				call T%pointer(Tp)
				Tp2(1:dim1,1:dim2)=>Tp

				Rdim1=dim1
				call Res%pointer(Rp)
				Rp=Tp2(:,dimi)
				
			end if
		
		end if
		call WorkingMemory%free()
		return
	end subroutine