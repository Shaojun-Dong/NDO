
	subroutine FuncName1(A,ip)
		class(Tensor),target,intent(in)::A
		FuncDataType,pointer,intent(inout)::ip(:)
		if(.not.A%getFlag())then
			call writemess('The Tensor is empty',-1)
			call error_stop
		end if
		if(A%getType().ne.DATATYPENum)then
			call writemess('ERROR in pointer, datatype='+DATATYPENum,-1)
			call writemess('A%getType()='+A%getType(),-1)
			call error_stop
		end if
		if(.not.A%Data%getFlag())then
			ip=>null()
			return
		end if
		call A%Data%pointer(ip)
		return
	end subroutine

	subroutine FuncName2(A,ip)
		class(Tensor),target,intent(in)::A
		FuncDataType,pointer,intent(inout)::ip(:,:)
		integer,pointer::dim(:)
		if(.not.A%getFlag())then
			call writemess('The Tensor is empty',-1)
			call error_stop
		end if
		if(A%Data%getTotalBlock().ne.1)then
			call writemess(' The Tensor is Symmetry type, error pointting',-1)
			call error_stop
		end if
		if(A%getRank().ne.2)then
			call writemess('ERROR in Tensor%pointer, rank',-1)
			call writemess('A%getRank()='+A%getRank())
			call error_stop
		end if
		if(A%getType().ne.DATATYPENum)then
			call writemess('ERROR in pointer, datatype='+DATATYPENum,-1)
			call writemess('A%getType()='+A%getType(),-1)
			call error_stop
		end if
		if(.not.A%Data%getFlag())then
			ip=>null()
			return
		end if
		call A%pointDim(dim)
		call Pointer2DFunc(A%getTotalData(),A%Data%ClassData,dim(1),dim(2),ip)
		return
	end subroutine

	subroutine FuncName3(A,ip)
		class(Tensor),target,intent(in)::A
		FuncDataType,pointer,intent(inout)::ip(:,:,:)
		integer,pointer::dim(:)
		if(.not.A%getFlag())then
			call writemess('The Tensor is empty',-1)
			call error_stop
		end if
		if(A%Data%getTotalBlock().ne.1)then
			call writemess(' The Tensor is Symmetry type, error pointting',-1)
			call error_stop
		end if
		if(A%getRank().ne.3)then
			call writemess('ERROR in Tensor%pointer, rank',-1)
			call writemess('A%getRank()='+A%getRank())
			call error_stop
		end if
		if(A%getType().ne.DATATYPENum)then
			call writemess('ERROR in pointer, datatype='+DATATYPENum,-1)
			call writemess('A%getType()='+A%getType(),-1)
			call error_stop
		end if
		if(.not.A%Data%getFlag())then
			ip=>null()
			return
		end if
		call A%pointDim(dim)
		call Pointer3DFunc(A%getTotalData(),A%Data%ClassData,dim(1),dim(2),dim(3),ip)
		return
	end subroutine

	subroutine FuncName4(A,ip)
		class(Tensor),target,intent(in)::A
		FuncDataType,pointer,intent(inout)::ip(:,:,:,:)
		integer,pointer::dim(:)
		if(.not.A%getFlag())then
			call writemess('The Tensor is empty',-1)
			call error_stop
		end if
		if(A%Data%getTotalBlock().ne.1)then
			call writemess(' The Tensor is Symmetry type, error pointting',-1)
			call error_stop
		end if
		if(A%getRank().ne.4)then
			call writemess('ERROR in Tensor%pointer, rank',-1)
			call writemess('A%getRank()='+A%getRank())
			call error_stop
		end if
		if(A%getType().ne.DATATYPENum)then
			call writemess('ERROR in pointer, datatype='+DATATYPENum,-1)
			call writemess('A%getType()='+A%getType(),-1)
			call error_stop
		end if
		if(.not.A%Data%getFlag())then
			ip=>null()
			return
		end if
		call A%pointDim(dim)
		call Pointer4DFunc(A%getTotalData(),A%Data%ClassData,dim(1),dim(2),dim(3),dim(4),ip)
		return
	end subroutine


	subroutine FuncBlockName_vec(A,ip,ith)
		class(Tensor),target,intent(in)::A
		FuncDataType,pointer,intent(inout)::ip(:)
		integer,intent(in)::ith(:)
		integer::index
		integer,pointer::si4(:,:,:,:),ei4(:,:,:,:)
		integer,pointer::si3(:,:,:),ei3(:,:,:)
		integer,pointer::si2(:,:),ei2(:,:)
		integer,pointer::Dim(:)
		if(.not.A%getFlag())then
			call writemess('The Tensor is empty',-1)
			call error_stop
		end if
		if(.not.A%getSymmetryFlag())then
			call writemess(' The Tensor is not Symmetry type, error pointting',-1)
			call error_stop
		end if
		if(A%getType().ne.DATATYPENum)then
			call writemess('ERROR in pointer, datatype='+DATATYPENum,-1)
			call writemess('A%getType()='+A%getType(),-1)
			call error_stop
		end if
		if(.not.A%Data%getFlag())then
			ip=>null()
			return
		end if
		call A%pointDim(Dim)
		call check_indices_and_dim(dim,ith)
		select case(size(ith))
			case(1)
				if(A%getFlag(ith(1)))then
					call A%Data%Pointer(ip,ith(1))
				else
					ip=>null()
				end if
			case(2)
				call A%Data%pointStarti(dim(1),dim(2),si2)
				call A%Data%pointEndi(dim(1),dim(2),ei2)
				if(ei2(ith(1),ith(2)).gt.0)then
					call Pointer1DFunc(A%Data%ClassData, si2(ith(1),ith(2)) , ei2(ith(1),ith(2)),ip )
				else
					ip=>null()
				end if
			case(3)
				call A%Data%pointStarti(dim(1),dim(2),dim(3),si3)
				call A%Data%pointEndi(dim(1),dim(2),dim(3),ei3)
				if(ei3(ith(1),ith(2),ith(3)).gt.0)then
					call Pointer1DFunc(A%Data%ClassData, si3(ith(1),ith(2),ith(3)) , ei3(ith(1),ith(2),ith(3)), ip )
				else
					ip=>null()
				end if
			case(4)
				call A%Data%pointStarti(dim(1),dim(2),dim(3),dim(4),si4)
				call A%Data%pointEndi(dim(1),dim(2),dim(3),dim(4),ei4)
				if(ei4(ith(1),ith(2),ith(3),ith(4)).gt.0)then
					call Pointer1DFunc(A%Data%ClassData,&
						si4(ith(1),ith(2),ith(3),ith(4)),ei4(ith(1),ith(2),ith(3),ith(4)),&
														 ip )
				else
					ip=>null()
				end if
			case default
				index=addressToIndes(ith,Dim)
				if(A%getFlag(index))then
					call A%Data%Pointer(ip,index)
				else
					ip=>null()
				end if
		end select
		return 
	end subroutine

	subroutine FuncBlockName1(A,ip,ith)
		class(Tensor),target,intent(in)::A
		FuncDataType,pointer,intent(inout)::ip(:)
		integer,intent(in)::ith
		if(.not.A%getFlag())then
			call writemess('The Tensor is empty',-1)
			call error_stop
		end if
		if(A%getType().ne.DATATYPENum)then
			call writemess('ERROR in pointer, datatype='+DATATYPENum,-1)
			call writemess('A%getType()='+A%getType(),-1)
			call error_stop
		end if
		if(.not.A%Data%getFlag())then
			ip=>null()
			return
		end if
		if(A%getFlag(ith))then
			call A%Data%pointer(ip,ith)
		else
			ip=>null()
		end if
		return 
	end subroutine

	subroutine FuncBlockName2(A,ip,ith)
		class(Tensor),target,intent(in)::A
		FuncDataType,pointer,intent(inout)::ip(:,:)
		class(*),pointer::ip0(:)
		integer,intent(in)::ith(:)
		integer,pointer::si(:,:),ei(:,:),dim(:)
		integer::blockdim(2)
		if(.not.A%getFlag())then
			call writemess('The Tensor is empty',-1)
			call error_stop
		end if
		if(.not.A%getSymmetryFlag())then
			call writemess(' The Tensor is not Symmetry type, error pointting',-1)
			call error_stop
		end if
		if(size(ith).ne.2)then
			call writemess('ERROR in Tensor%pointer, size(ith)',-1)
			call writemess('size(ith)='+size(ith))
			call error_stop
		end if
		if(A%getType().ne.DATATYPENum)then
			call writemess('ERROR in pointer, datatype='+DATATYPENum,-1)
			call writemess('A%getType()='+A%getType(),-1)
			call error_stop
		end if
		if(.not.A%Data%getFlag())then
			ip=>null()
			return
		end if
		call A%pointDim(dim)
		call check_indices_and_dim(dim,ith)
		call A%Data%pointStarti(dim(1),dim(2),Si)
		call A%Data%pointEndi(dim(1),dim(2),Ei)

		blockdim=A%getBlockDim(ith)
		if(ei(ith(1),ith(2)).gt.0)then
			call classPointer1DFunc(A%Data%ClassData, si(ith(1),ith(2)) , ei(ith(1),ith(2)),ip0 )
			call Pointer2DFunc(ei(ith(1),ith(2))-si(ith(1),ith(2))+1,ip0,blockdim(1),blockdim(2),ip)
		else
			ip=>null()
		end if
		return 
	end subroutine

	subroutine FuncDimBlockName2(A,Adim,ip,Blockdim,ith)
		class(Tensor),target,intent(in)::A
		FuncDataType,pointer,intent(inout)::ip(:,:)
		class(*),pointer::ip0(:)
		integer,intent(in)::ith(:),Adim(:),Blockdim(:)
		integer,pointer::si(:,:),ei(:,:)
		integer::totalBlock
		if(.not.A%getFlag())then
			call writemess('The Tensor is empty',-1)
			call error_stop
		end if
		if(Blockdim(1).eq.0)then
			ip=>null()
			return
		end if
		if(Blockdim(2).eq.0)then
			ip=>null()
			return
		end if
		totalBlock=A%Data%getTotalBlock()
		if(.not.A%getSymmetryFlag())then
			call writemess(' The Tensor is not Symmetry type, error pointting',-1)
			call error_stop
		end if
		if(size(Adim).ne.2)then
			call writemess('ERROR  dim',-1)
			call writemess('size(Adim)='+size(Adim))
			call error_stop
		end if
		if(size(Blockdim).ne.2)then
			call writemess('ERROR  Blockdim',-1)
			call writemess('size(Blockdim)='+size(Blockdim))
			call error_stop
		end if
		if(A%getType().ne.DATATYPENum)then
			call writemess('ERROR in pointer, datatype='+DATATYPENum,-1)
			call writemess('A%getType()='+A%getType(),-1)
			call error_stop
		end if
		call check_indices_and_dim(Adim,ith,totalBlock)
		if(.not.A%Data%getFlag())then
			ip=>null()
			return
		end if
		call A%Data%pointStarti(Adim(1),Adim(2),Si)
		call A%Data%pointEndi(Adim(1),Adim(2),Ei)
		if(product(Blockdim).ne.(ei(ith(1),ith(2))-si(ith(1),ith(2))+1))then
			call writemess('ERROR in point to the block',-1)
			call writemess('product(Blockdim)='+product(Blockdim),-1)
			call writemess('ei(ith(1),ith(2))-si(ith(1),ith(2))+1='+(ei(ith(1),ith(2))-si(ith(1),ith(2))+1),-1)
			call writemess('si(ith(1),ith(2))='+si(ith(1),ith(2)),-1)
			call writemess('ei(ith(1),ith(2))='+ei(ith(1),ith(2)),-1)
			call writemess('ith(1)='+ith(1),-1)
			call writemess('ith(2)='+ith(2),-1)
			call error_stop
		end if
		if(ei(ith(1),ith(2)).gt.0)then
			call ClassPointer1DFunc(A%Data%ClassData, si(ith(1),ith(2)) , ei(ith(1),ith(2)),ip0 )
			call Pointer2DFunc(ei(ith(1),ith(2))-si(ith(1),ith(2))+1,ip0,blockdim(1),blockdim(2),ip)
		else
			ip=>null()
		end if
		return 
	end subroutine

	subroutine FuncBlockName3(A,ip,ith)
		class(Tensor),target,intent(in)::A
		FuncDataType,pointer,intent(inout)::ip(:,:,:)
		integer,intent(in)::ith(:)
		class(*),pointer::ip0(:)
		integer,pointer::si(:,:,:),ei(:,:,:),dim(:)
		integer::bd(3),totalBlock
		if(.not.A%getFlag())then
			call writemess('The Tensor is empty',-1)
			call error_stop
		end if
		totalBlock=A%Data%getTotalBlock()
		if(.not.A%getSymmetryFlag())then
			call writemess(' The Tensor is not Symmetry type, error pointting',-1)
			call error_stop
		end if
		if(A%getRank().ne.3)then
			call writemess('ERROR in Tensor%pointer, rank',-1)
			call writemess('A%getRank()='+A%getRank())
			call error_stop
		end if
		if(A%getType().ne.DATATYPENum)then
			call writemess('ERROR in pointer, datatype='+DATATYPENum,-1)
			call writemess('A%getType()='+A%getType(),-1)
			call error_stop
		end if
		call A%pointDim(dim)
		call check_indices_and_dim(dim,ith)
		if(.not.A%Data%getFlag())then
			ip=>null()
			return
		end if
		call A%Data%pointStarti(dim(1),dim(2),dim(3),Si)
		call A%Data%pointEndi(dim(1),dim(2),dim(3),Ei)

		bd=A%getBlockDim(ith)

		if(ei(ith(1),ith(2),ith(3)).gt.0)then
			call ClassPointer1DFunc(A%Data%ClassData, si(ith(1),ith(2),ith(3)) , ei(ith(1),ith(2),ith(3)),ip0 )
			call Pointer3DFunc(ei(ith(1),ith(2),ith(3))-si(ith(1),ith(2),ith(3))+1,&
							ip0,bd(1),bd(2),bd(3),ip)
		else
			ip=>null()
		end if
		return 
	end subroutine

	subroutine FuncDimBlockName3(A,Adim,ip,Blockdim,ith)
		class(Tensor),target,intent(in)::A
		FuncDataType,pointer,intent(inout)::ip(:,:,:)
		integer,intent(in)::ith(:),Adim(:),Blockdim(:)
		integer,pointer::si(:,:,:),ei(:,:,:)
		class(*),pointer::ip0(:)
		integer::endindex,startindex,totalBlock
		if(.not.A%getFlag())then
			call writemess('The Tensor is empty',-1)
			call error_stop
		end if
		if(Blockdim(1).eq.0)then
			ip=>null()
			return
		end if
		if(Blockdim(2).eq.0)then
			ip=>null()
			return
		end if
		if(Blockdim(3).eq.0)then
			ip=>null()
			return
		end if
		totalBlock=A%Data%getTotalBlock()
		if(.not.A%getSymmetryFlag())then
			call writemess(' The Tensor is not Symmetry type, error pointting',-1)
			call error_stop
		end if
		if(size(Adim).ne.3)then
			call writemess('ERROR  dim',-1)
			call writemess('size(Adim)='+size(Adim))
			call error_stop
		end if
		if(size(Blockdim).ne.3)then
			call writemess('ERROR  Blockdim',-1)
			call writemess('size(Blockdim)='+size(Blockdim))
			call error_stop
		end if
		if(A%getType().ne.DATATYPENum)then
			call writemess('ERROR in pointer, datatype='+DATATYPENum,-1)
			call writemess('A%getType()='+A%getType(),-1)
			call error_stop
		end if
		call check_indices_and_dim(Adim,ith,totalBlock)
		if(.not.A%Data%getFlag())then
			ip=>null()
			return
		end if
		call A%Data%pointStarti(Adim(1),Adim(2),Adim(3),Si)
		call A%Data%pointEndi(Adim(1),Adim(2),Adim(3),Ei)
		startindex=si(ith(1),ith(2),ith(3))
		endindex=ei(ith(1),ith(2),ith(3))
		if(product(Blockdim).ne.(endindex-startindex+1))then
			call writemess('ERROR in point to the block',-1)
			call writemess('product(Blockdim)='+product(Blockdim),-1)
			call writemess('endindex-startindex+1='+(endindex-startindex+1),-1)
			call writemess('startindex='+startindex,-1)
			call writemess('endindex='+endindex,-1)
			call writemess('ith(1)='+ith(1),-1)
			call writemess('ith(2)='+ith(2),-1)
			call writemess('ith(3)='+ith(3),-1)
			call error_stop
		end if
		if(endindex.gt.0)then
			call ClassPointer1DFunc(A%Data%ClassData, startindex , endindex,ip0 )
			call Pointer3DFunc(endindex-startindex+1,ip0,Blockdim(1),Blockdim(2),Blockdim(3),ip)
		else
			ip=>null()
		end if
		return 
	end subroutine


	subroutine FuncBlockName4(A,ip,ith)
		class(Tensor),target,intent(in)::A
		FuncDataType,pointer,intent(inout)::ip(:,:,:,:)
		integer,intent(in)::ith(:)
		integer,pointer::si(:,:,:,:),ei(:,:,:,:),dim(:)
		integer::totalBlock,bd(4),endindex,startindex
		class(*),pointer::ip0(:)
		if(.not.A%getFlag())then
			call writemess('The Tensor is empty',-1)
			call error_stop
		end if
		totalBlock=A%Data%getTotalBlock()
		if(.not.A%getSymmetryFlag())then
			call writemess(' The Tensor is not Symmetry type, error pointting',-1)
			call error_stop
		end if
		if(A%getRank().ne.4)then
			call writemess('ERROR in Tensor%pointer, rank',-1)
			call writemess('A%getRank()='+A%getRank())
			call error_stop
		end if
		if(A%getType().ne.DATATYPENum)then
			call writemess('ERROR in pointer, datatype='+DATATYPENum,-1)
			call writemess('A%getType()='+A%getType(),-1)
			call error_stop
		end if
		call A%pointDim(dim)
		call check_indices_and_dim(dim,ith)
		if(.not.A%Data%getFlag())then
			ip=>null()
			return
		end if
		call A%Data%pointStarti(dim(1),dim(2),dim(3),dim(4),Si)
		call A%Data%pointEndi(dim(1),dim(2),dim(3),dim(4),Ei)

		bd=A%getBlockDim(ith)
		startindex=si(ith(1),ith(2),ith(3),ith(4))
		endindex=ei(ith(1),ith(2),ith(3),ith(4))
		if(endindex.gt.0)then
			call ClassPointer1DFunc(A%Data%ClassData,startindex , endindex,ip0 )
			call Pointer4DFunc(endindex-startindex+1,ip0,bd(1),bd(2),bd(3),bd(4), ip)
		else
			ip=>null()
		end if
		return 
	end subroutine


	subroutine PrintName(ip,dim,w)
		FuncDataType,target,intent(in)::ip(:)
		integer,intent(in)::dim(:)
		character(len=*),optional,intent(in)::w
		FuncDataType2,pointer::ip2(:,:),ip3(:,:,:),ip4(:,:,:,:)
		integer::rank,i,j,k
		rank=size(dim)
		select case(rank)
			case (1)
				if(present(w))then
					call writemess(ip,w)
				else
					call writemess(ip)
				end if
			case (2)
				ip2(1:dim(1),1:dim(2))=>ip
				do i=1,dim(1)
					if(present(w))then
						call writemess(ip2(i,:),w)
					else
						call writemess(ip2(i,:))
					end if
				end do
			case (3)
				ip3(1:dim(1),1:dim(2),1:dim(3))=>ip
				do j=1,dim(3)
					call writemess('(*,*,'+j+')')
					do i=1,dim(1)
						if(present(w))then
							call writemess(ip3(i,:,j),w)
						else
							call writemess(ip3(i,:,j))
						end if
					end do
				end do
			case (4)
				ip4(1:dim(1),1:dim(2),1:dim(3),1:dim(4))=>ip
				do k=1,dim(4)
					do j=1,dim(3)
						call writemess('(*,*,'+j+','+k+')')
						do i=1,dim(1)
							if(present(w))then
								call writemess(ip4(i,:,j,k),w)
							else
								call writemess(ip4(i,:,j,k))
							end if
						end do
					end do
				end do
			case default
				if(present(w))then
					call writemess(ip,w)
				else
					call writemess(ip)
				end if
				call writemess('dimension:')
				call writemess(dim)
		end select
		return
	end subroutine



	subroutine StoreDataType(A,Block,indexstart,indexend,Dim)
		FuncDataType,intent(inout)::A(:)
		FuncDataType,intent(in)::Block(:)
		integer,intent(in)::indexstart(:),indexend(:),Dim(:)
		integer,allocatable::indices(:)
		logical::goon
		integer::ith,Blockith
		allocate(indices(size(indexstart)))

		indices=indexstart
		goon=.true.
		Blockith=0
		do while(goon)
			ith=addressToIndes(indices,dim)
			Blockith=Blockith+1
			if(Blockith.gt.size(Block))then
				call writemess("ERROR in StoreData,SymTensor",-1)
				call writemess('size(Block)='+size(Block),-1)
				call writemess('Blockith='+Blockith,-1)
				call writemess('indexstart, indexend and dim are',-1)
				call writemess(indexstart,-1)
				call writemess(indexend,-1)
				call writemess(Dim,-1)
				call error_Stop()
			end if
			A(ith)=Block(Blockith)
			goon=index_counter(indices,indexstart,indexend,1)
		end do
		return
	end subroutine
	subroutine StoreSomeDataType(A,dimA,Block,dimBlock,indexstart,indexend,Dims,dime)
		FuncDataType,intent(inout)::A(:)
		FuncDataType,intent(in)::Block(:)
		integer,intent(in)::indexstart(:),indexend(:),Dims(:),dime(:),dimA(:),dimBlock(:)
		integer,allocatable::indices(:),indices2(:)
		logical::goon,goon2
		integer::ith,Blockith
		allocate(indices(size(indexstart)))
		allocate(indices2(size(Dims)))

		indices=indexstart
		indices2=Dims
		goon=.true.
		Blockith=0
		do while(goon)
			ith=addressToIndes(indices,dimA)
			Blockith=addressToIndes(indices2,dimBlock)
			if(Blockith.gt.size(Block))then
				call writemess("ERROR in StoreData,SymTensor",-1)
				call writemess('size(Block)='+size(Block),-1)
				call writemess('Blockith='+Blockith,-1)
				call error_Stop()
			end if
			A(ith)=Block(Blockith)
			goon=index_counter(indices,indexstart,indexend,1)
			goon2=index_counter(indices2,Dims,dime,1)
		end do
		return
	end subroutine

