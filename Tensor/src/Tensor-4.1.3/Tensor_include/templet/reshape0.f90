	subroutine subTensorDimDEDATANAME(A,DimMin,DimMax,Res)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		integer,intent(in)::DimMin(:),DimMax(:)
		DEDATANAME,pointer::Resip(:),Aip(:)
		DEDATANAME,pointer::Resip2(:,:),Aip2(:,:)
		DEDATANAME,pointer::Resip3(:,:,:),Aip3(:,:,:)
		DEDATANAME,pointer::Resip4(:,:,:,:),Aip4(:,:,:,:)
		DEDATANAME,pointer::Resip5(:,:,:,:,:),Aip5(:,:,:,:,:)
		integer,pointer::index(:),dim(:)
		character(len=Len_of_Name),pointer::AName(:)
		integer::i,ii,rank,TotalData
		logical::goon
		if(A%getSymmetryFlag())then
			call writemess('ERROR in subTensorDim. The input Tensor is of Symmery',-1)
			call error_stop
		end if
		rank=A%getRank()
		if(size(DimMin).ne.rank)then
			call writemess('ERROR in subTensorDim,rank',-1)
			call writemess('size(DimMin)='+size(DimMin),-1)
			call writemess(DimMin,-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call error_stop
		end if
		if(size(DimMax).ne.rank)then
			call writemess('ERROR in subTensorDim,rank',-1)
			call writemess('size(DimMax)='+size(DimMax),-1)
			call writemess(DimMax,-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call error_stop
		end if
		if(rank.eq.1)then
			TotalData=DimMax(1)-DimMin(1)+1
			call Res%allocate([TotalData],A%getType())
			call Res%pointer(Resip)
			call A%pointer(Aip)
			Resip=Aip(DimMin(1):DimMax(1))
			if(A%getNameFlag())then
				call A%pointName(AName)
				call Res%setName(1,AName(1))
			end if
			return
		end if

		if(rank.eq.2)then
			call Res%allocate((DimMax-DimMin+1),A%getType())
			call Res%pointer(Resip2)
			call A%pointer(Aip2)
			TotalData=Res%getTotalData()
			Resip2=Aip2(DimMin(1):DimMax(1),DimMin(2):DimMax(2))
			if(A%getNameFlag())then
				call A%pointName(AName)
				do i=1,Res%getRank()
					call Res%setName(i,AName(i))
				end do
			end if
			return
		end if

		if(rank.eq.3)then
			call Res%allocate((DimMax-DimMin+1),A%getType())
			call Res%pointer(Resip3)
			call A%pointer(Aip3)
			TotalData=Res%getTotalData()
			Resip3=Aip3(DimMin(1):DimMax(1),DimMin(2):DimMax(2),DimMin(3):DimMax(3))
			if(A%getNameFlag())then
				call A%pointName(AName)
				do i=1,Res%getRank()
					call Res%setName(i,AName(i))
				end do
			end if
			return
		end if

		if(rank.eq.4)then
			call Res%allocate((DimMax-DimMin+1),A%getType())
			call Res%pointer(Resip4)
			call A%pointer(Aip4)
			TotalData=Res%getTotalData()
			Resip4=Aip4(DimMin(1):DimMax(1),DimMin(2):DimMax(2),DimMin(3):DimMax(3)&
									,DimMin(4):DimMax(4))
			if(A%getNameFlag())then
				call A%pointName(AName)
				do i=1,Res%getRank()
					call Res%setName(i,AName(i))
				end do
			end if
			return
		end if

		call Res%allocate((DimMax-DimMin+1),A%getType())
		call Res%pointer(Resip)
		call A%pointer(Aip)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank)
		call WorkingMemory%get_memory(index,rank) 
		index=DimMin
		call A%pointDim(dim)
		do i=1,Res%getTotalData()
			if(.not.goon)then
				call writemess('ERROR in subTensorDim',-1)
				call error_stop
			end if
			ii=addressToIndes(index,dim)
			Resip(i)=Aip(ii)
			goon=index_counter(index,DimMax)
		end do
		if(A%getNameFlag())then
			call A%pointName(AName)
			do i=1,Res%getRank()
				call Res%setName(i,AName(i))
			end do
		end if
		call WorkingMemory%free()
		return
	end subroutine

	subroutine subTensorMaxDimDEDATANAME(A,DimMax,Res)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		integer,intent(in)::DimMax(:)
		DEDATANAME,pointer::Resip(:),Aip(:)
		DEDATANAME,pointer::Resip2(:,:),Aip2(:,:)
		DEDATANAME,pointer::Resip3(:,:,:),Aip3(:,:,:)
		DEDATANAME,pointer::Resip4(:,:,:,:),Aip4(:,:,:,:)
		DEDATANAME,pointer::Resip5(:,:,:,:,:),Aip5(:,:,:,:,:)
		integer,pointer::index(:),dim(:)
		character(len=Len_of_Name),pointer::AName(:)
		integer::i,ii,rank,TotalData
		logical::goon
		if(A%getSymmetryFlag())then
			call writemess('ERROR in subTensorDim. The input Tensor is of Symmery',-1)
			call error_stop
		end if
		rank=A%getRank()
		if(size(DimMax).ne.rank)then
			call writemess('ERROR in subTensorDim,rank',-1)
			call writemess('size(DimMax)='+size(DimMax),-1)
			call writemess(DimMax,-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call error_stop
		end if
		if(rank.eq.1)then
			TotalData=DimMax(1)
			call Res%allocate([TotalData],A%getType())
			call Res%pointer(Resip)
			call A%pointer(Aip)
			Resip=Aip(1:DimMax(1))
			if(A%getNameFlag())then
				call A%pointName(AName)
				call Res%setName(1,AName(1))
			end if
			return
		end if

		if(rank.eq.2)then
			call Res%allocate(DimMax,A%getType())
			call Res%pointer(Resip2)
			call A%pointer(Aip2)
			TotalData=Res%getTotalData()
			Resip2=Aip2(1:DimMax(1),1:DimMax(2))
			if(A%getNameFlag())then
				call A%pointName(AName)
				do i=1,Res%getRank()
					call Res%setName(i,AName(i))
				end do
			end if
			return
		end if

		if(rank.eq.3)then
			call Res%allocate(DimMax,A%getType())
			call Res%pointer(Resip3)
			call A%pointer(Aip3)
			TotalData=Res%getTotalData()
			Resip3=Aip3(1:DimMax(1),1:DimMax(2),1:DimMax(3))
			if(A%getNameFlag())then
				call A%pointName(AName)
				do i=1,Res%getRank()
					call Res%setName(i,AName(i))
				end do
			end if
			return
		end if

		if(rank.eq.4)then
			call Res%allocate(DimMax,A%getType())
			call Res%pointer(Resip4)
			call A%pointer(Aip4)
			TotalData=Res%getTotalData()
			Resip4=Aip4(1:DimMax(1),1:DimMax(2),1:DimMax(3),1:DimMax(4))
			if(A%getNameFlag())then
				call A%pointName(AName)
				do i=1,Res%getRank()
					call Res%setName(i,AName(i))
				end do
			end if
			return
		end if

		call Res%allocate(DimMax,A%getType())
		call Res%pointer(Resip)
		call A%pointer(Aip)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,rank)
		call WorkingMemory%get_memory(index,rank) 
		index=1
		call A%pointDim(dim)
		do i=1,Res%getTotalData()
			if(.not.goon)then
				call writemess('ERROR in subTensorDim',-1)
				call error_stop
			end if
			ii=addressToIndes(index,dim)
			Resip(i)=Aip(ii)
			goon=index_counter(index,DimMax)
		end do
		if(A%getNameFlag())then
			call A%pointName(AName)
			do i=1,Res%getRank()
				call Res%setName(i,AName(i))
			end do
		end if
		call WorkingMemory%free()
		return
	end subroutine



	subroutine subTensorDeg1DDEDATANAME(A,DegMinMax,Res)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		integer,intent(in)::DegMinMax(:,:)
		DEDATANAME,pointer::Rip(:),Aip(:)
		integer,pointer::deg(:),dim(:),ii(:)
		character(len=Len_of_Name),pointer::AName(:)
		integer::i,rank,blocklength,iilen
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in subTensorDim. The input Tensor is not of Symmery',-1)
			call error_stop
		end if
		rank=A%getRank()
		if(rank.ne.1)then
			call writemess('ERROR in subTensorDeg1DSubroutine',-1)
			call error_stop
		end if
		if(size(DegMinMax,1).ne.2)then
			call writemess('ERROR in subTensorDim,rank',-1)
			call writemess('size(DegMinMax,1)='+size(DegMinMax,1),-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call error_stop
		end if
		if(size(DegMinMax,2).ne.A%dim(1))then
			call writemess('ERROR in subTensorDim,rank',-1)
			call writemess('size(DegMinMax,2)='+size(DegMinMax,2),-1)
			call writemess(DegMinMax(1,:),-1)
			call writemess(DegMinMax(2,:),-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call error_stop
		end if
		call Res%allocateMomery(A%Dimension,A%getTotalData(),A%getType())
		call Res%Dimension%setDeg(1,DegMinMax(2,:)-DegMinMax(1,:)+1)

		call Res%pointDim(Dim)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,Dim(1))
		call WorkingMemory%get_memory(ii,dim(1))

		call Res%pointDeg(Deg,1)
		iilen=0
		do i=1,Dim(1)
			if(Deg(i).le.0)then
				Deg(i)=0
			else
				iilen=iilen+1
				ii(iilen)=i
			end if
		end do
		if(iilen.eq.0)then
			call writemess(' The subdeg will delete one leg in the tensor',-1)
			call writemess('DegMinMax(1,:)=',-1)
			call writemess(DegMinMax(1,:),-1)
			call writemess('DegMinMax(2,:)=',-1)
			call writemess(DegMinMax(2,:),-1)
			call error_stop
		end if
		if(iiLen.ne.dim(1))then
			call Res%killZeroDeg()
			call Res%pointDim(Dim)
			call Res%Data%allocateDataArrayMomery(product(Dim),A%getTotalData(),A%getType())
		else
			call Res%pointDim(Dim)
		end if
		
		do i=1,Dim(1)
			if(A%getFlag(ii(i)))then
				blocklength=(DegMinMax(2,ii(i))-DegMinMax(1,ii(i))+1)
				if(blocklength.gt.0)then
					call A%pointer(Aip,ii(i))
					call Res%setBlockMomery([i])
					call REs%pointer(Rip,i)
					if(size(Rip).ne.blocklength)then
						call writemess('ERROR in subTensorDeg2DSubroutine',-1)
						call error_stop
					end if
					if(size(Aip).lt.blocklength)then
						call writemess('ERROR in subTensorDeg2DSubroutine',-1)
						call error_stop
					end if
					Rip=Aip(DegMinMax(1,ii(i)):DegMinMax(2,ii(i)))
				end if
			end if
		end do
		call WorkingMemory%free()
		return
	end subroutine

	subroutine subTensorMaxDeg1DDEDATANAME(A,DegMax,Res)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		integer,intent(in)::DegMax(:)
		DEDATANAME,pointer::Rip(:),Aip(:)
		integer,pointer::deg(:),dim(:),ii(:)
		character(len=Len_of_Name),pointer::AName(:)
		integer::i,rank,blocklength,iilen
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in subTensorDim. The input Tensor is not of Symmery',-1)
			call error_stop
		end if
		rank=A%getRank()
		if(rank.ne.1)then
			call writemess('ERROR in subTensorDeg1DSubroutine',-1)
			call error_stop
		end if
		if(size(DegMax).ne.A%dim(1))then
			call writemess('ERROR in subTensorDim,rank',-1)
			call writemess('size(DegMax)='+size(DegMax),-1)
			call writemess(DegMax(:),-1)
			call writemess('A%getRank()='+A%getRank(),-1)
			call error_stop
		end if
		call Res%allocateMomery(A%Dimension,A%getTotalData(),A%getType())
		call Res%Dimension%setDeg(1,DegMax)

		call Res%pointDim(Dim)
		call WorkingMemory%check()
		call WorkingMemory%allocate(1,Dim(1))
		call WorkingMemory%get_memory(ii,dim(1))


		call Res%pointDeg(Deg,1)
		iilen=0
		do i=1,Dim(1)
			if(Deg(i).le.0)then
				Deg(i)=0
			else
				iilen=iilen+1
				ii(iilen)=i
			end if
		end do
		if(iilen.eq.0)then
			call writemess(' The subdeg will delete one leg in the tensor',-1)
			call writemess('DegMax=',-1)
			call writemess(DegMax,-1)
			call error_stop
		end if
		if(iiLen.ne.dim(1))then
			call Res%killZeroDeg()
			call Res%pointDim(Dim)
			call Res%Data%allocateDataArrayMomery(product(Dim),A%getTotalData(),A%getType())
		else
			call Res%pointDim(Dim)
		end if
		
		do i=1,Dim(1)
			if(A%getFlag(ii(i)))then
				blocklength=DegMax(ii(i))
				if(blocklength.gt.0)then
					call A%pointer(Aip,ii(i))
					call Res%setBlockMomery([i])
					call REs%pointer(Rip,i)
					if(size(Rip).ne.blocklength)then
						call writemess('ERROR in subTensorDeg2DSubroutine',-1)
						call error_stop
					end if
					if(size(Aip).lt.blocklength)then
						call writemess('ERROR in subTensorDeg2DSubroutine',-1)
						call error_stop
					end if
					Rip=Aip(1:DegMax(ii(i)))
				end if
			end if
		end do
		call WorkingMemory%free()
		return
	end subroutine


	subroutine subTensorDeg2DDEDATANAME(A,DegMinMax1,DegMinMax2,Res)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		integer,intent(in)::DegMinMax1(:,:),DegMinMax2(:,:)
		DEDATANAME,pointer::Rip(:,:),Aip(:,:)
		integer,pointer::dim(:),deg(:),ii(:),jj(:)
		character(len=Len_of_Name),pointer::AName(:)
		integer::i,j,rank,blocklength,iilen,jjlen
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in subTensorDim. The input Tensor is not of Symmery',-1)
			call error_stop
		end if
		rank=A%getRank()
		if(rank.ne.2)then
			call writemess('ERROR in subTensorDeg1DSubroutine',-1)
			call error_stop
		end if
		call Res%allocateMomery(A%Dimension,A%getTotalData(),A%getType())
		call Res%Dimension%setDeg(1,DegMinMax1(2,:)-DegMinMax1(1,:)+1)
		call Res%Dimension%setDeg(2,DegMinMax2(2,:)-DegMinMax2(1,:)+1)

		call Res%pointDim(Dim)

		call WorkingMemory%check()
		call WorkingMemory%allocate(1,sum(Dim))
		call WorkingMemory%get_memory(ii,dim(1))
		call WorkingMemory%get_memory(jj,dim(2))


		call Res%pointDeg(Deg,1)
		iilen=0
		do i=1,Dim(1)
			if(Deg(i).le.0)then
				Deg(i)=0
			else
				iilen=iilen+1
				ii(iilen)=i
			end if
		end do
		if(iilen.eq.0)then
			call writemess(' The subdeg will delete one leg in the tensor',-1)
			call writemess('DegMinMax1(1,:)=',-1)
			call writemess(DegMinMax1(1,:),-1)
			call writemess('DegMinMax1(2,:)=',-1)
			call writemess(DegMinMax1(2,:),-1)
			call error_stop
		end if

		call Res%pointDeg(Deg,2)
		jjlen=0
		do i=1,Dim(2)
			if(Deg(i).le.0)then
				Deg(i)=0
			else
				jjlen=jjlen+1
				jj(jjlen)=i
			end if
		end do
		if(jjlen.eq.0)then
			call writemess(' The subdeg will delete one leg in the tensor',-1)
			call writemess('DegMinMax2(1,:)=',-1)
			call writemess(DegMinMax2(1,:),-1)
			call writemess('DegMinMax2(2,:)=',-1)
			call writemess(DegMinMax2(2,:),-1)
			call error_stop
		end if
		if((iiLen.ne.dim(1)).or.(jjLen.ne.dim(2)))then
			call Res%killZeroDeg()
			call Res%pointDim(Dim)
			call Res%Data%allocateDataArrayMomery(product(Dim),A%getTotalData(),A%getType())
		else
			call Res%pointDim(Dim)
		end if
		do j=1,dim(2)
			do i=1,dim(1)
				if(A%getFlag([ii(i),jj(j)]))then
					blocklength=(DegMinMax1(2,ii(i))-DegMinMax1(1,ii(i))+1)*(DegMinMax2(2,jj(j))-DegMinMax2(1,jj(j))+1)
					if(blocklength.gt.0)then
						call A%pointer(Aip,[ii(i),jj(j)])
						call Res%setBlockMomery([i,j])
						call REs%pointer(Rip,[i,j])
						if(size(Rip).ne.blocklength)then
							call writemess('ERROR in subTensorDeg2DSubroutine',-1)
							call error_stop
						end if
						if(size(Aip).lt.blocklength)then
							call writemess('ERROR in subTensorDeg2DSubroutine',-1)
							call error_stop
						end if
						Rip=Aip(DegMinMax1(1,ii(i)):DegMinMax1(2,ii(i)),&
									DegMinMax2(1,jj(j)):DegMinMax2(2,jj(j)))
					end if
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end subroutine

	subroutine subTensorMaxDeg2DDEDATANAME(A,DegMax1,DegMax2,Res)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		integer,intent(in)::DegMax1(:),DegMax2(:)
		DEDATANAME,pointer::Rip(:,:),Aip(:,:)
		integer,pointer::dim(:),deg(:),ii(:),jj(:)
		character(len=Len_of_Name),pointer::AName(:)
		integer::i,j,rank,blocklength,iilen,jjlen
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in subTensorDim. The input Tensor is not of Symmery',-1)
			call error_stop
		end if
		rank=A%getRank()
		if(rank.ne.2)then
			call writemess('ERROR in subTensorDeg1DSubroutine',-1)
			call error_stop
		end if
		call Res%allocateMomery(A%Dimension,A%getTotalData(),A%getType())
		call Res%Dimension%setDeg(1,DegMax1)
		call Res%Dimension%setDeg(2,DegMax2)


		call Res%pointDim(Dim)

		call WorkingMemory%check()
		call WorkingMemory%allocate(1,sum(Dim))
		call WorkingMemory%get_memory(ii,dim(1))
		call WorkingMemory%get_memory(jj,dim(2))

		call Res%pointDeg(Deg,1)
		iilen=0
		do i=1,Dim(1)
			if(Deg(i).le.0)then
				Deg(i)=0
			else
				iilen=iilen+1
				ii(iilen)=i
			end if
		end do
		if(iilen.eq.0)then
			call writemess(' The subdeg will delete one leg in the tensor',-1)
			call writemess('DegMax1=',-1)
			call writemess(DegMax1,-1)
			call error_stop
		end if

		call Res%pointDeg(Deg,2)
		jjlen=0
		do i=1,Dim(2)
			if(Deg(i).le.0)then
				Deg(i)=0
			else
				jjlen=jjlen+1
				jj(jjlen)=i
			end if
		end do
		if(jjlen.eq.0)then
			call writemess(' The subdeg will delete one leg in the tensor',-1)
			call writemess('DegMax2=',-1)
			call writemess(DegMax2,-1)
			call error_stop
		end if
		if((iiLen.ne.dim(1)).or.(jjLen.ne.dim(2)))then
			call Res%killZeroDeg()
			call Res%pointDim(Dim)
			call Res%Data%allocateDataArrayMomery(product(Dim),A%getTotalData(),A%getType())
		else
			call Res%pointDim(Dim)
		end if

		
		do j=1,dim(2)
			do i=1,dim(1)
				if(A%getFlag([ii(i),jj(j)]))then
					blocklength=DegMax1(ii(i))*DegMax2(jj(j))
					if(blocklength.gt.0)then
						call A%pointer(Aip,[ii(i),jj(j)])
						call Res%setBlockMomery([i,j])
						call REs%pointer(Rip,[i,j])
						if(size(Rip).ne.blocklength)then
							call writemess('ERROR in subTensorDeg2DSubroutine',-1)
							call writemess('i='+i+', j='+j+' ,blocklength='+blocklength)
							call writemess('ii(i)='+ii(i)+', jj(j)='+jj(j))
							call writemess('DegMax1(ii(i))='+DegMax1(ii(i)))
							call writemess('DegMax2(jj(j))='+DegMax2(jj(j)))
							call writemess('size(Rip)='+size(Rip))
							call Res%diminfo(.true.)
							call error_stop
						end if
						if(size(Aip).lt.blocklength)then
							call writemess('ERROR in subTensorDeg2DSubroutine',-1)
							call error_stop
						end if
						Rip=Aip(1:DegMax1(ii(i)),1:DegMax2(jj(j)))
					end if
				end if
			end do
		end do
		call WorkingMemory%free()
		return
	end subroutine

	subroutine subTensorDeg3DDEDATANAME(A,DegMinMax1,DegMinMax2,DegMinMax3,Res)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		integer,intent(in)::DegMinMax1(:,:),DegMinMax2(:,:),DegMinMax3(:,:)
		DEDATANAME,pointer::Rip(:,:,:),Aip(:,:,:)
		integer,pointer::dim(:),deg(:),ii(:),jj(:),kk(:)
		character(len=Len_of_Name),pointer::AName(:)
		integer::i,j,k,rank,blocklength,iilen,jjlen,kklen
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in subTensorDim. The input Tensor is not of Symmery',-1)
			call error_stop
		end if
		rank=A%getRank()
		if(rank.ne.3)then
			call writemess('ERROR in subTensorDeg1DSubroutine',-1)
			call error_stop
		end if
		call Res%allocateMomery(A%Dimension,A%getTotalData(),A%getType())
		call Res%Dimension%setDeg(1,DegMinMax1(2,:)-DegMinMax1(1,:)+1)
		call Res%Dimension%setDeg(2,DegMinMax2(2,:)-DegMinMax2(1,:)+1)
		call Res%Dimension%setDeg(3,DegMinMax3(2,:)-DegMinMax3(1,:)+1)

		call Res%pointDim(Dim)

		call WorkingMemory%check()
		call WorkingMemory%allocate(1,sum(Dim))
		call WorkingMemory%get_memory(ii,dim(1))
		call WorkingMemory%get_memory(jj,dim(2))

		call Res%pointDeg(Deg,1)
		iilen=0
		do i=1,Dim(1)
			if(Deg(i).le.0)then
				Deg(i)=0
			else
				iilen=iilen+1
				ii(iilen)=i
			end if
		end do
		if(iilen.eq.0)then
			call writemess(' The subdeg will delete one leg in the tensor',-1)
			call writemess('DegMinMax1(1,:)=',-1)
			call writemess(DegMinMax1(1,:),-1)
			call writemess('DegMinMax1(2,:)=',-1)
			call writemess(DegMinMax1(2,:),-1)
			call error_stop
		end if

		call Res%pointDeg(Deg,2)
		jjlen=0
		do i=1,Dim(2)
			if(Deg(i).le.0)then
				Deg(i)=0
			else
				jjlen=jjlen+1
				jj(jjlen)=i
			end if
		end do
		if(jjlen.eq.0)then
			call writemess(' The subdeg will delete one leg in the tensor',-1)
			call writemess('DegMinMax2(1,:)=',-1)
			call writemess(DegMinMax2(1,:),-1)
			call writemess('DegMinMax2(2,:)=',-1)
			call writemess(DegMinMax2(2,:),-1)
			call error_stop
		end if

		call Res%pointDeg(Deg,3)
		kklen=0
		do i=1,Dim(3)
			if(Deg(i).le.0)then
				Deg(i)=0
			else
				kklen=kklen+1
				kk(kklen)=i
			end if
		end do
		if(kklen.eq.0)then
			call writemess(' The subdeg will delete one leg in the tensor',-1)
			call writemess('DegMinMax3(1,:)=',-1)
			call writemess(DegMinMax3(1,:),-1)
			call writemess('DegMinMax3(2,:)=',-1)
			call writemess(DegMinMax3(2,:),-1)
			call error_stop
		end if

		if((iiLen.ne.dim(1)).or.(jjLen.ne.dim(2)).or.(kkLen.ne.dim(3)))then
			call Res%killZeroDeg()
			call Res%pointDim(Dim)
			call Res%Data%allocateDataArrayMomery(product(Dim),A%getTotalData(),A%getType())
		else
			call Res%pointDim(Dim)
		end if
		do k=1,dim(3)
			do j=1,dim(2)
				do i=1,dim(1)
					if(A%getFlag([ii(i),jj(j),kk(k)]))then
						blocklength= (DegMinMax1(2,ii(i))-DegMinMax1(1,ii(i))+1) &
									*(DegMinMax2(2,jj(j))-DegMinMax2(1,jj(j))+1)&
									*(DegMinMax3(2,kk(k))-DegMinMax3(1,kk(k))+1)
						if(blocklength.gt.0)then
							call A%pointer(Aip,[ii(i),jj(j),kk(k)])
							call Res%setBlockMomery([i,j,k])
							call REs%pointer(Rip,[i,j,k])
							if(size(Rip).ne.blocklength)then
								call writemess('ERROR in subTensorDeg2DSubroutine',-1)
								call error_stop
							end if
							if(size(Aip).lt.blocklength)then
								call writemess('ERROR in subTensorDeg2DSubroutine',-1)
								call error_stop
							end if
							Rip=Aip(DegMinMax1(1,ii(i)):DegMinMax1(2,ii(i)),&
												  DegMinMax2(1,jj(j)):DegMinMax2(2,jj(j)),&
												  DegMinMax3(1,kk(k)):DegMinMax3(2,kk(k)))
						end if
					end if
				end do
			end do
		end do
		call WorkingMemory%free()
		return
	end subroutine


	subroutine subTensorMaxDeg3DDEDATANAME(A,DegMax1,DegMax2,DegMax3,Res)
		class(Tensor),intent(in)::A
		type(Tensor),intent(inout)::Res
		integer,intent(in)::DegMax1(:),DegMax2(:),DegMax3(:)
		DEDATANAME,pointer::Rip(:,:,:),Aip(:,:,:)
		integer,pointer::dim(:),deg(:),ii(:),jj(:),kk(:)
		character(len=Len_of_Name),pointer::AName(:)
		integer::i,j,k,rank,blocklength,iilen,jjlen,kklen
		if(.not.A%getSymmetryFlag())then
			call writemess('ERROR in subTensorDim. The input Tensor is not of Symmery',-1)
			call error_stop
		end if
		rank=A%getRank()
		if(rank.ne.3)then
			call writemess('ERROR in subTensorDeg1DSubroutine',-1)
			call error_stop
		end if
		call Res%allocateMomery(A%Dimension,A%getTotalData(),A%getType())
		call Res%Dimension%setDeg(1,DegMax1)
		call Res%Dimension%setDeg(2,DegMax2)
		call Res%Dimension%setDeg(3,DegMax3)


		call Res%pointDim(Dim)

		call WorkingMemory%check()
		call WorkingMemory%allocate(1,sum(Dim))
		call WorkingMemory%get_memory(ii,dim(1))
		call WorkingMemory%get_memory(jj,dim(2))
		call WorkingMemory%get_memory(kk,dim(3))

		call Res%pointDeg(Deg,1)
		iilen=0
		do i=1,Dim(1)
			if(Deg(i).le.0)then
				Deg(i)=0
			else
				iilen=iilen+1
				ii(iilen)=i
			end if
		end do
		if(iilen.eq.0)then
			call writemess(' The subdeg will delete one leg in the tensor',-1)
			call writemess('DegMax1=',-1)
			call writemess(DegMax1,-1)
			call error_stop
		end if

		call Res%pointDeg(Deg,2)
		jjlen=0
		do i=1,Dim(2)
			if(Deg(i).le.0)then
				Deg(i)=0
			else
				jjlen=jjlen+1
				jj(jjlen)=i
			end if
		end do
		if(jjlen.eq.0)then
			call writemess(' The subdeg will delete one leg in the tensor',-1)
			call writemess('DegMax2=',-1)
			call writemess(DegMax2,-1)
			call error_stop
		end if

		call Res%pointDeg(Deg,3)
		kklen=0
		do i=1,Dim(3)
			if(Deg(i).le.0)then
				Deg(i)=0
			else
				kklen=kklen+1
				kk(kklen)=i
			end if
		end do
		if(kklen.eq.0)then
			call writemess(' The subdeg will delete one leg in the tensor',-1)
			call writemess('DegMax3=',-1)
			call writemess(DegMax3,-1)
			call error_stop
		end if

		if((iiLen.ne.dim(1)).or.(jjLen.ne.dim(2)).or.(kkLen.ne.dim(3)))then
			call Res%killZeroDeg()
			call Res%pointDim(Dim)
			call Res%Data%allocateDataArrayMomery(product(Dim),A%getTotalData(),A%getType())
		else
			call Res%pointDim(Dim)
		end if
		do k=1,dim(3)
			do j=1,dim(2)
				do i=1,dim(1)
					if(A%getFlag([ii(i),jj(j),kk(k)]))then
						blocklength= DegMax1(ii(i))*DegMax2(jj(j))*DegMax3(kk(k))
						if(blocklength.gt.0)then
							call A%pointer(Aip,[ii(i),jj(j),kk(k)])
							call Res%setBlockMomery([i,j,k])
							call REs%pointer(Rip,[i,j,k])
							if(size(Rip).ne.blocklength)then
								call writemess('ERROR in subTensorDeg2DSubroutine',-1)
								call error_stop
							end if
							if(size(Aip).lt.blocklength)then
								call writemess('ERROR in subTensorDeg2DSubroutine',-1)
								call error_stop
							end if
							Rip=Aip(1:DegMax1(ii(i)),&
												  1:DegMax2(jj(j)),&
												  1:DegMax3(kk(k)))
						end if
					end if
				end do
			end do
		end do
		call WorkingMemory%free()
		return
	end subroutine


