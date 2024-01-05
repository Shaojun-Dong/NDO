	!Tensor-3.5.8  SVDcutoff_kill_inData TO SPEED Up

	subroutine SVDDataSubroutineDATATYPE(A,U,S,V,m,n,min_MN,info,ClassType,cut)
		integer,intent(in)::min_MN
		DATATYPE,intent(inout)::A(:,:)
		DATATYPE,intent(inout)::U(:,:),V(:,:)
		DATATYPE2,intent(inout)::s(:)
		integer,intent(in)::ClassType,info
		integer,intent(in),optional::cut
		complex*16,pointer :: zwork(:),uzdata(:,:),vzdata(:,:)
		complex*8,pointer :: cwork(:),ucdata(:,:),vcdata(:,:)
		real*8,pointer :: dwork(:),uddata(:,:),vddata(:,:),sdata8(:)
		real*4,pointer :: swork(:),usdata(:,:),vsdata(:,:),sdata4(:)
		integer m,n,lw,rwdim,i,j
		integer,pointer :: iw(:)
		real*8,pointer :: drw(:)
		real*4,pointer :: srw(:)
		integer :: ms_max,max_MN,length,cut2
		integer::totallengthofmemory
		call WorkingMemory%check()
		length=size(A)
		max_MN=max(M,N)
		lw=min_MN*min_MN+2*min_MN+max_MN+1
		rwdim=min_MN*max(5*min_MN+7,2*max_MN+2*min_MN+1)+1
		select case(ClassType)
			case (5)
				totallengthofmemory=length+lw+m*min_MN+min_MN*n
				call WorkingMemory%allocate(5,totallengthofmemory)
				totallengthofmemory=rwdim+min_MN
				call WorkingMemory%allocate(3,totallengthofmemory)
				call WorkingMemory%allocate(1,8*min_MN)
				call WorkingMemory%get_memory(drw,rwdim)
				call WorkingMemory%get_memory(iw,8*min_MN)
				call WorkingMemory%get_memory(zwork,lw)
				if(present(cut))then
					cut2=min(cut,min_MN)
					call WorkingMemory%get_memory(uzdata,m,min_MN)
					call WorkingMemory%get_memory(vzdata,min_MN,n)
					call WorkingMemory%get_memory(sdata8,min_MN)
					call ZGESVD('S','S',m,n,A,m,sdata8,uzdata,m,vzdata,min_MN,zWORK,lw,drw,INFO)
					S=sdata8(1:cut2)
					U=uzdata(:,1:cut2)
					V=Vzdata(1:cut2,:)
				else
					call ZGESVD('S','S',m,n,A,m,S,U,m,V,min_MN,zWORK,lw,drw,INFO)
				end if
			case(4)
				if(WorkingMemory%ifDynamic())then
					totallengthofmemory=length+lw+m*min_MN+min_MN*n
					call WorkingMemory%allocate(4,totallengthofmemory)
					totallengthofmemory=rwdim+min_MN
					call WorkingMemory%allocate(2,totallengthofmemory)
					call WorkingMemory%allocate(1,8*min_MN)
				end if
				call WorkingMemory%get_memory(srw,rwdim)
				call WorkingMemory%get_memory(iw,8*min_MN)
				call WorkingMemory%get_memory(cwork,lw)
				if(present(cut))then
					cut2=min(cut,min_MN)
					call WorkingMemory%get_memory(ucdata,m,min_MN)
					call WorkingMemory%get_memory(vcdata,min_MN,n)
					call WorkingMemory%get_memory(sdata4,min_MN)
					call CGESVD('S','S',m,n,A,m,sdata4,ucdata,m,vcdata,min_MN,cWORK,lw,srw,INFO)
					S=sdata4(1:cut2)
					U=ucdata(:,1:cut2)
					V=vcdata(1:cut2,:)
				else
					call CGESVD('S','S',m,n,A,m,S,U,m,V,min_MN,cWORK,lw,srw,INFO)
				end if
			case(3)
				if(WorkingMemory%ifDynamic())then
					totallengthofmemory=length+lw+m*min_MN+min_MN*n+rwdim+min_MN
					call WorkingMemory%allocate(3,totallengthofmemory)
					call WorkingMemory%allocate(1,8*min_MN)
				end if
				call WorkingMemory%get_memory(drw,rwdim)
				call WorkingMemory%get_memory(iw,8*min_MN)
				call WorkingMemory%get_memory(dwork,lw)
				if(present(cut))then
					cut2=min(cut,min_MN)
					call WorkingMemory%get_memory(uddata,m,min_MN)
					call WorkingMemory%get_memory(vddata,min_MN,n)
					call WorkingMemory%get_memory(sdata8,min_MN)
					call DGESVD('S','S',m,n,A,m,sdata8,uddata,m,vddata,min_MN,dWORK,lw,INFO)
					S=sdata8(1:cut2)
					U=uddata(:,1:cut2)
					V=vddata(1:cut2,:)
				else
					call DGESVD('S','S',m,n,A,m,S,U,m,V,min_MN,dWORK,lw,INFO)
				end if
			case(2)
				if(WorkingMemory%ifDynamic())then
					totallengthofmemory=length+lw+m*min_MN+min_MN*n+rwdim+min_MN
					call WorkingMemory%allocate(2,totallengthofmemory)
					call WorkingMemory%allocate(1,8*min_MN)
				end if
				call WorkingMemory%get_memory(srw,rwdim)
				call WorkingMemory%get_memory(iw,8*min_MN)
				call WorkingMemory%get_memory(swork,lw)
				if(present(cut))then
					cut2=min(cut,min_MN)
					call WorkingMemory%get_memory(usdata,m,min_MN)
					call WorkingMemory%get_memory(vsdata,min_MN,n)
					call WorkingMemory%get_memory(sdata4,min_MN)
					call SGESVD('S','S',m,n,A,m,sdata4,usdata,m,vsdata,min_MN,sWORK,lw,INFO)
					S=sdata4(1:cut2)
					U=usdata(:,1:cut2)
					V=vsdata(1:cut2,:)
				else
					call SGESVD('S','S',m,n,A,m,S,U,m,V,min_MN,sWORK,lw,INFO)
				end if
		end select
		call WorkingMemory%free()
		return
	end subroutine



	subroutine SVDSavingSingleValueDegDATATYPE(S,Newdeg,inNumSave,AllSindex,AllSdata)
		type(Tensor),intent(in)::S
		integer,intent(inout)::Newdeg(:)
		integer,intent(inout)::AllSindex(:)
		integer,intent(in)::inNumSave
		DATATYPE2,pointer,intent(inout)::AllSdata(:)
		integer::i,j,k,NumSave,lenNewdeg,i1,i2
		NumSave=min(size(AllSdata),inNumSave)
		i1=0
		i2=0
		do i=1,S%getTotalBlock()
			i1=i2+1
			i2=i2+S%getTotalData(i)
			AllSindex(i1:i2)=i
		end do
		call reorderDataDATATYPE(AllSdata,AllSindex)
		Newdeg=0
		lenNewdeg=size(Newdeg)
		do i=1,NumSave
			if(AllSindex(i).gt.lenNewdeg)then
				call writemess('ERROR in SVDSavingSingleValueDeg',-1)
				call error_stop
			end if
			Newdeg(AllSindex(i))=Newdeg(AllSindex(i))+1
		end do
		return
	end subroutine

	subroutine SVDSavingSingleValueDegValueDATATYPE(S,Newdeg,minNumSave,maxNumSave,NumSave,maxValue,&
		VType,AllSindex,AllSdata,outmaxValue)
		integer,intent(in)::minNumSave,maxNumSave
		DATATYPE2,intent(inout)::AllSdata(:)
		integer,intent(inout)::AllSindex(:)
		type(Tensor),intent(in)::S
		integer,intent(inout)::NumSave
		character(len=*),intent(in)::VType
		real*8,intent(in)::maxValue
		real*8,intent(inout)::outmaxValue
		integer,intent(inout)::Newdeg(:)
		integer::i,j,k,Total,lenNewdeg
		real*8::sumValue,tempr
		integer::AllData,i1,i2
		Total=size(AllSdata)
		NumSave=min(maxNumSave,total)
		i1=0
		i2=0
		do i=1,S%getTotalBlock()
			i1=i2+1
			i2=i2+S%getTotalData(i)
			AllSindex(i1:i2)=i
		end do
		AllData=Total
		call reorderDataDATATYPE(AllSdata,AllSindex)
		Newdeg=0
		sumValue=0
		Total=0
		lenNewdeg=size(Newdeg)
		tempr=1d0/AllSdata(1)
		do i=1,NumSave
			Total=Total+1
			if(AllSindex(i).gt.lenNewdeg)then
				call writemess('ERROR in SVDSavingSingleValueDegValue',-1)
				call error_stop
			end if
			Newdeg(AllSindex(i))=Newdeg(AllSindex(i))+1
			if(VType.equ.'sum2')then
				sumValue=sumValue+(AllSdata(i)*AllSdata(i))
				if((sumValue.gt.maxValue).and.(i.ge.minNumSave))then
					sumValue=AllSdata(max(i-1,1))
					Total=Total-1
					Newdeg(AllSindex(i))=Newdeg(AllSindex(i))-1
					exit
				end if
			else if(VType.equ.'1_n') then
				sumValue=abs(AllSdata(i)*tempr)
				if( ( sumValue.le.maxValue).and.(i.ge.minNumSave)) then
					sumValue=AllSdata(max(i-1,1))
					Total=Total-1
					Newdeg(AllSindex(i))=Newdeg(AllSindex(i))-1
					exit
				end if
			else if(VType.equ.'1n') then
				sumValue=abs(AllSdata(i)*tempr)
				if( ( sumValue.le.maxValue).and.(i.ge.minNumSave)) then
					sumValue=AllSdata(max(i-1,1))
					Total=Total-1
					Newdeg(AllSindex(i))=Newdeg(AllSindex(i))-1
					exit
				end if
			else if(VType.equ.'sum')then
				sumValue=sumValue+AllSdata(i)
				if( (sumValue.gt.maxValue).and.(i.ge.minNumSave))then
					sumValue=AllSdata(max(i-1,1))
					Total=Total-1
					Newdeg(AllSindex(i))=Newdeg(AllSindex(i))-1
					exit
				end if
			else if(VType.equ.'min')then
				sumValue=AllSdata(i)
				if( (sumValue.le.maxValue).and.(i.ge.minNumSave))then
					sumValue=AllSdata(max(i-1,1))
					Total=Total-1
					Newdeg(AllSindex(i))=Newdeg(AllSindex(i))-1
					exit
				end if
			else
				call writemess('No such case in SVD',-1)
				call writemess('VType='+VType,-1)
				call error_stop
			end if
		end do
		if(Total.eq.AllData)then
			outmaxValue=0d0
		else
			outmaxValue=sumValue
		end if
		NumSave=Total
		return
	end subroutine

	subroutine SVDSavingSingleValueNonSymmetryDATATYPE(NewDimLen,minNumSave,maxNumSave,&
		NumSave,maxValue,VType,AllSdata,outmaxValue)
		integer,intent(in)::minNumSave,maxNumSave
		DATATYPE2,intent(inout)::AllSdata(:)
		integer,intent(inout)::NumSave
		character(len=*),intent(in)::VType
		real*8,intent(in)::maxValue
		real*8,intent(inout)::outmaxValue
		integer,intent(inout)::NewDimLen
		integer::i,j,k,Total
		real*8::sumValue,tempr
		integer::AllData,i1,i2
		Total=size(AllSdata)
		NumSave=min(maxNumSave,total)
		AllData=Total
		NewDimLen=0
		sumValue=0
		Total=0
		tempr=1d0/AllSdata(1)
		do i=1,NumSave
			Total=Total+1
			NewDimLen=NewDimLen+1
			if(VType.equ.'sum2')then
				sumValue=sumValue+(AllSdata(i)*AllSdata(i))
				if((sumValue.gt.maxValue).and.(i.ge.minNumSave))then
					sumValue=AllSdata(max(i-1,1))
					Total=Total-1
					NewDimLen=NewDimLen-1
					exit
				end if
			else if(VType.equ.'1_n') then
				sumValue=abs(AllSdata(i)*tempr)
				if( ( sumValue.le.maxValue).and.(i.ge.minNumSave)) then
					sumValue=AllSdata(max(i-1,1))
					Total=Total-1
					NewDimLen=NewDimLen-1
					exit
				end if
			else if(VType.equ.'1n') then
				sumValue=abs(AllSdata(i)*tempr)
				if( ( sumValue.le.maxValue).and.(i.ge.minNumSave)) then
					sumValue=AllSdata(max(i-1,1))
					Total=Total-1
					NewDimLen=NewDimLen-1
					exit
				end if
			else if(VType.equ.'sum')then
				sumValue=sumValue+AllSdata(i)
				if( (sumValue.gt.maxValue).and.(i.ge.minNumSave))then
					sumValue=AllSdata(max(i-1,1))
					Total=Total-1
					NewDimLen=NewDimLen-1
					exit
				end if
			else if(VType.equ.'min')then
				sumValue=AllSdata(i)
				if( (sumValue.le.maxValue).and.(i.ge.minNumSave))then
					sumValue=AllSdata(max(i-1,1))
					Total=Total-1
					NewDimLen=NewDimLen-1
					exit
				end if
			else
				call writemess('No such case in SVD',-1)
				call writemess('VType='+VType,-1)
				call error_stop
			end if
		end do
		if(Total.eq.AllData)then
			outmaxValue=0d0
		else
			outmaxValue=sumValue
		end if
		NumSave=Total
		return
	end subroutine

	subroutine reorderDataDATATYPE(a,inde)
		DATATYPE2,intent(inout) :: a(:)
		integer,intent(inout):: inde(:)
		DATATYPE2 :: temp
		integer :: i,j,n,tempi
		n=size(a)
		do i=1,n-1
			 do j=i+1,n
				  if (a(i) .lt. a(j)) then
						temp = a(i)
						tempi=inde(i)
						a(i) = a(j)
						inde(i)=inde(j)
						a(j) = temp
						inde(j)=tempi
				  endif
			 enddo
		enddo
		return
	end subroutine

	function ifNonZeroBlockDATATYPE(A,ith,jth)result(Res)
		logical::res
		type(Tensor),intent(in)::A
		integer,intent(in)::ith,jth
		DATATYPE,pointer::Ap(:)
		DATATYPE::zero
		integer::i
		Res=.false.
		zero=0d0
		call A%pointer(Ap,[ith,jth])
		do i=1,size(Ap)
			if(Ap(i).ne.zero)then
				Res=.true.
				return
			end if
		end do
		return
	end function

	subroutine SVDMatrixKillNumDATATYPE(A,inoutU,inoutS,inoutV,CUTOFF)
		class(Tensor),intent(inout)::A
		type(Tensor),target,intent(inout)::inoutU,inoutS,inoutV
		integer,optional,intent(in)::CUTOFF
		DATATYPE,pointer::Ap(:,:),Up(:,:),Vp(:,:)
		DATATYPE2,pointer::sp(:),AllSdata(:)
		integer::M,N,min_MN,classtype,info,Newdimi,testi,testj
		type(Dimension)::UDim,VDim
		type(QuanNum)::NewQN(2)
		integer::TotalData
		integer,pointer::dim(:)
		real*4,pointer::QN1(:),QN2(:)
		integer::i,j,indexs,STotalData,DegLen
		integer,pointer::AllSindex(:)
		type(Tensor),pointer::U,S,V
		classtype=A%getType()
		if(A%getRank().ne.2)then
			call writemess('ERROR in SVD, the input tensor should be a matrix',-1)
			call error_Stop
		end if
		call inoutS%empty()
		if(.not.A%getSymmetryFlag())then
			U=>inoutU
			S=>inoutS
			V=>inoutV
			m=A%dim(1)
			n=A%dim(2)
			min_MN=min(M,N)
			if(present(CUTOFF))then
				Newdimi=min(min_MN,CUTOFF)
			else
				Newdimi=min_MN
			end if
			call pasteDimension(UDim,A%Dimension,1,1,[Newdimi])
			call pasteDimension(VDim,[Newdimi],A%Dimension,2,2)
			call U%allocate(UDim,classtype)
			call V%Allocate(VDim,classtype)
			select case(classtype)
				case(5)
					call s%allocate([Newdimi],3)
				case(4)
					call s%allocate([Newdimi],2)
				case default
					call s%allocate([Newdimi],classtype)
			end select
			call A%pointer(Ap)
			call U%pointer(Up)
			call V%pointer(Vp)
			call S%pointer(Sp)
			call SVDDataSubroutineDATATYPE(Ap,Up,Sp,Vp,m,n,min_MN,info,classtype,CUTOFF)
			if(info.ne.0) then
				call writemess('Error in svd ,info='+info,-1)
				call writemess('output The data in ./_SVD_ERROR_LOG.err',-1)
				open(unit=9991,file='./_SVD_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
				call A%write(9991)
				call U%write(9991)
				call S%write(9991)
				call V%write(9991)
				close(9991)
				call error_stop()
			end if
			return
		end if

		if(present(CUTOFF))then
			U=>SVD_U
			S=>SVD_S
			V=>SVD_V
		else
			U=>inoutU
			S=>inoutS
			V=>inoutV
		end if
		call choosSameQN(NewQN,A%Dimension)
		call pasteDimension(U%Dimension,A%Dimension,1,1,NewQN(1))
		call pasteDimension(V%Dimension,NewQN(2),A%Dimension,2,2)
		TotalData=A%getTotalData()
		call U%pointDim(dim)
		call U%Data%allocateDataArrayMomery(product(dim),TotalData,A%getType())
		call V%pointDim(dim)
		call V%Data%allocateDataArrayMomery(product(dim),TotalData,A%getType())
		call V%pointDeg(dim,1)
		select case(classtype)
			case(5)
				call S%Data%allocateDataArrayMomery(NewQN(1)%getQNlength(),sum(NewQN(1)%getDeg()),3)
			case(4)
				call S%Data%allocateDataArrayMomery(NewQN(1)%getQNlength(),sum(NewQN(1)%getDeg()),2)
			case default
				call S%Data%allocateDataArrayMomery(NewQN(1)%getQNlength(),sum(NewQN(1)%getDeg()),A%getType())
		end select
		call A%pointDim(dim)
		call A%pointQN(QN1,1)
		call A%pointQN(QN2,2)
		testi=0
		testj=0
		do j=1,dim(2)
			do i=1,dim(1)
				if(A%getFlag([i,j]).and.ifNonZeroBlockDATATYPE(A,i,j))then
					if(i.lt.testi)then ! U(1), S1=1.5 S2=0.5, the diag mement is not (i,i) but (i+1,i)
						call writemess("ERROR in SVD on SymTensor",-1)
						call writemess("The quantum number are:"+QN1(i)+' , '+QN2(j),-1)
						call A%diminfo(.true.)
						call error_stop()
					end if
					if(j.lt.testj)then ! U(1), S1=1.5 S2=0.5, the diag mement is not (i,i) but (i+1,i)
						call writemess("ERROR in SVD on SymTensor",-1)
						call writemess("The quantum number are:"+QN1(i)+' , '+QN2(j),-1)
						call A%diminfo(.true.)
						call error_stop()
					end if
					testi=i
					testj=j
					indexs=NewQN(1)%getIndex(QN2(j))
					if(indexs.eq.0)then
						call writemess('ERRRO in SVD',-1)
						call NewQN(1)%print
						call writemess('j='+j)
						call writemess('QN2(j)='+QN2(j))
						call error_stop
					end if
					call A%pointer(Ap,[i,j])
					call U%setBlockMomery([i,indexs])
					call U%pointer(Up,[i,indexs])
					call V%setBlockMomery([indexs,j])
					call V%pointer(Vp,[indexs,j])
					m=size(Ap,1)
					n=size(Ap,2)
					min_MN=min(m,n)
					call S%Data%Set_block_momery(indexs,min_MN)
					call S%Data%pointer(Sp,indexs)
					call SVDDataSubroutineDATATYPE(Ap,Up,Sp,Vp,m,n,min_MN,info,classtype)
				end if
			end do
		end do
		if(present(CUTOFF))then
			DegLen=U%dim(2)
			STotalData=S%getTotalData()
			call WorkingMemory%allocate(1,STotalData)
			call WorkingMemory%allocate(classtype,STotalData)
			call WorkingMemory%get_memory(AllSindex,STotalData)
			call WorkingMemory%get_memory(AllSdata,STotalData)
			call S%Data%pointAllData(sp)
			AllSdata=Sp
			call allocateCheck(SVDCUTDeg,DegLen)
			call SVDSavingSingleValueDegDATATYPE(S,SVDCUTDeg(1:DegLen),CUTOFF,AllSindex,AllSdata)
			call WorkingMemory%free()
			call U%CutTensorDeg(2,SVDCUTDeg(1:DegLen),inoutU)
			call V%CutTensorDeg(1,SVDCUTDeg(1:DegLen),inoutV)
			call S%Data%CutOffDataArray(SVDCUTDeg(1:DegLen),inoutS%Data)
		end if
		return
	end subroutine

	subroutine SVDMatrixKillValueDATATYPE(A,inoutU,inoutS,inoutV,minNumSave,maxNumSave,&
		inmaxValue,VType,outNumSave,outmaxValue)
		class(Tensor),intent(inout)::A
		type(Tensor),target,intent(inout)::inoutU,inoutS,inoutV
		integer,intent(in)::minNumSave,maxNumSave
		real*8,intent(in)::inmaxValue
		character(len=*),intent(in)::VType
		integer,optional,intent(inout)::outNumSave
		real*8,optional,intent(inout)::outmaxValue
		real*8::maxValue
		DATATYPE,pointer::Ap(:,:),Up(:,:),Vp(:,:)
		DATATYPE2,pointer::sp(:),AllSdata(:)
		integer::M,N,min_MN,classtype,info,Newdimi,NumSave
		type(Dimension)::UDim,VDim
		type(QuanNum)::NewQN(2)
		integer::TotalData
		integer,pointer::dim(:)
		real*4,pointer::QN1(:),QN2(:)
		integer::i,j,indexs,STotalData,DegLen,NewDimLen,testi,testj
		integer,pointer::AllSindex(:)
		type(Tensor),pointer::U,S,V
		classtype=A%getType()
		if(A%getRank().ne.2)then
			call writemess('ERROR in SVD, the input tensor should be a matrix',-1)
			call error_Stop
		end if
		call inoutS%empty()
		U=>SVD_U
		S=>SVD_S
		V=>SVD_V
		if(.not.A%getSymmetryFlag())then
			m=A%dim(1)
			n=A%dim(2)
			min_MN=min(M,N)
			Newdimi=min_MN
			call pasteDimension(UDim,A%Dimension,1,1,[Newdimi])
			call pasteDimension(VDim,[Newdimi],A%Dimension,2,2)
			call U%allocate(UDim,classtype)
			call V%Allocate(VDim,classtype)
			select case(classtype)
				case(5)
					call s%allocate([Newdimi],3)
				case(4)
					call s%allocate([Newdimi],2)
				case default
					call s%allocate([Newdimi],classtype)
			end select
			call A%pointer(Ap)
			call U%pointer(Up)
			call V%pointer(Vp)
			call S%pointer(Sp)
			call SVDDataSubroutineDATATYPE(Ap,Up,Sp,Vp,m,n,min_MN,info,classtype)
			if(info.ne.0) then
				call writemess('Error in svd ,info='+info,-1)
				call writemess('output The data in ./_SVD_ERROR_LOG.err',-1)
				open(unit=9991,file='./_SVD_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
				call A%write(9991)
				call U%write(9991)
				call S%write(9991)
				call V%write(9991)
				close(9991)
				call error_stop()
			end if
			STotalData=S%getTotalData()
			call WorkingMemory%allocate(classtype,STotalData)
			call WorkingMemory%get_memory(AllSdata,STotalData)
			AllSdata=Sp
			call SVDSavingSingleValueNonSymmetryDATATYPE(NewDimLen,minNumSave,maxNumSave,&
											NumSave,inmaxValue,VType,AllSdata,maxValue)
			call WorkingMemory%free()
			call U%pointDim(dim)
			call U%CutTensorDim([dim(1),NewDimLen],inoutU)
			call V%pointDim(dim)
			call V%CutTensorDim([NewDimLen,dim(2)],inoutV)
			call S%CutTensorDim([NewDimLen],inoutS)
			if(present(outNumSave))outNumSave=NumSave
			if(present(outNumSave))outmaxValue=maxValue
			return
		end if

		call choosSameQN(NewQN,A%Dimension)
		call pasteDimension(U%Dimension,A%Dimension,1,1,NewQN(1))
		call pasteDimension(V%Dimension,NewQN(2),A%Dimension,2,2)
		TotalData=A%getTotalData()
		call U%pointDim(dim)
		call U%Data%allocateDataArrayMomery(product(dim),TotalData,A%getType())
		call V%pointDim(dim)
		call V%Data%allocateDataArrayMomery(product(dim),TotalData,A%getType())
		call V%pointDeg(dim,1)
		select case(classtype)
			case(5)
				call S%Data%allocateDataArrayMomery(NewQN(1)%getQNlength(),sum(NewQN(1)%getDeg()),3)
			case(4)
				call S%Data%allocateDataArrayMomery(NewQN(1)%getQNlength(),sum(NewQN(1)%getDeg()),2)
			case default
				call S%Data%allocateDataArrayMomery(NewQN(1)%getQNlength(),sum(NewQN(1)%getDeg()),A%getType())
		end select
		call A%pointDim(dim)
		call A%pointQN(QN1,1)
		call A%pointQN(QN2,2)
		testi=0
		testj=0
		do j=1,dim(2)
			do i=1,dim(1)
				if(A%getFlag([i,j]).and.ifNonZeroBlockDATATYPE(A,i,j))then
					if(i.lt.testi)then ! U(1), S1=1.5 S2=0.5, the diag mement is not (i,i) but (i+1,i)
						call writemess("ERROR in SVD on SymTensor",-1)
						call writemess("The quantum number are:"+QN1(i)+' , '+QN2(j),-1)
						call A%diminfo(.true.)
						call error_stop()
					end if
					if(j.lt.testj)then ! U(1), S1=1.5 S2=0.5, the diag mement is not (i,i) but (i+1,i)
						call writemess("ERROR in SVD on SymTensor",-1)
						call writemess("The quantum number are:"+QN1(i)+' , '+QN2(j),-1)
						call A%diminfo(.true.)
						call error_stop()
					end if
					testi=i
					testj=j
					indexs=NewQN(1)%getIndex(QN2(j))
					if(indexs.eq.0)then
						call writemess('ERRRO in SVD',-1)
						call error_stop
					end if
					call A%pointer(Ap,[i,j])
					call U%setBlockMomery([i,indexs])
					call U%pointer(Up,[i,indexs])
					call V%setBlockMomery([indexs,j])
					call V%pointer(Vp,[indexs,j])
					m=size(Ap,1)
					n=size(Ap,2)
					min_MN=min(m,n)
					call S%Data%Set_block_momery(indexs,min_MN)
					call S%Data%pointer(Sp,indexs)
					call SVDDataSubroutineDATATYPE(Ap,Up,Sp,Vp,m,n,min_MN,info,classtype)
				end if
			end do
		end do

		DegLen=U%dim(2)
		STotalData=S%getTotalData()
		call WorkingMemory%allocate(1,STotalData)
		call WorkingMemory%allocate(classtype,STotalData)
		call WorkingMemory%get_memory(AllSindex,STotalData)
		call WorkingMemory%get_memory(AllSdata,STotalData)
		call S%Data%pointAllData(sp)
		AllSdata=Sp
		call allocateCheck(SVDCUTDeg,DegLen)
		call SVDSavingSingleValueDegValueDATATYPE(S,SVDCUTDeg(1:DegLen),minNumSave,maxNumSave&
									,NumSave,inmaxValue,VType,AllSindex,AllSdata,maxValue)
		call WorkingMemory%free()
		call U%CutTensorDeg(2,SVDCUTDeg(1:DegLen),inoutU)
		call V%CutTensorDeg(1,SVDCUTDeg(1:DegLen),inoutV)
		call S%Data%CutOffDataArray(SVDCUTDeg(1:DegLen),inoutS%Data)
		if(present(outNumSave))outNumSave=NumSave
		if(present(outmaxValue))outmaxValue=maxValue
		return
	end subroutine

	subroutine SVDMatrixKillNum2DATATYPE(A,inoutU,inoutS,inoutV,WT1,WT2,WT3,CUTOFF)
		class(Tensor),intent(inout)::A
		type(Tensor),target,intent(inout)::inoutU,inoutS,inoutV,WT1,WT2,WT3
		integer,optional,intent(in)::CUTOFF
		DATATYPE,pointer::Ap(:,:),Up(:,:),Vp(:,:)
		DATATYPE2,pointer::sp(:),AllSdata(:)
		integer::M,N,min_MN,classtype,info,Newdimi,testi,testj
		type(Dimension)::UDim,VDim
		type(QuanNum)::NewQN(2)
		integer::TotalData
		integer,pointer::dim(:)
		real*4,pointer::QN1(:),QN2(:)
		integer::i,j,indexs,STotalData,DegLen
		integer,pointer::AllSindex(:)
		type(Tensor),pointer::U,S,V
		classtype=A%getType()
		if(A%getRank().ne.2)then
			call writemess('ERROR in SVD, the input tensor should be a matrix',-1)
			call error_Stop
		end if
		call inoutS%empty()
		if(.not.A%getSymmetryFlag())then
			U=>inoutU
			S=>inoutS
			V=>inoutV
			m=A%dim(1)
			n=A%dim(2)
			min_MN=min(M,N)
			if(present(CUTOFF))then
				Newdimi=min(min_MN,CUTOFF)
			else
				Newdimi=min_MN
			end if
			call pasteDimension(UDim,A%Dimension,1,1,[Newdimi])
			call pasteDimension(VDim,[Newdimi],A%Dimension,2,2)
			call U%allocate(UDim,classtype)
			call V%Allocate(VDim,classtype)
			select case(classtype)
				case(5)
					call s%allocate([Newdimi],3)
				case(4)
					call s%allocate([Newdimi],2)
				case default
					call s%allocate([Newdimi],classtype)
			end select
			call A%pointer(Ap)
			call U%pointer(Up)
			call V%pointer(Vp)
			call S%pointer(Sp)
			call SVDDataSubroutineDATATYPE(Ap,Up,Sp,Vp,m,n,min_MN,info,classtype,CUTOFF)
			if(info.ne.0) then
				call writemess('Error in svd ,info='+info,-1)
				call writemess('output The data in ./_SVD_ERROR_LOG.err',-1)
				open(unit=9991,file='./_SVD_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
				call A%write(9991)
				call U%write(9991)
				call S%write(9991)
				call V%write(9991)
				close(9991)
				call error_stop()
			end if
			return
		end if

		if(present(CUTOFF))then
			U=>WT1
			S=>WT2
			V=>WT3
		else
			U=>inoutU
			S=>inoutS
			V=>inoutV
		end if
		call choosSameQN(NewQN,A%Dimension)
		call pasteDimension(U%Dimension,A%Dimension,1,1,NewQN(1))
		call pasteDimension(V%Dimension,NewQN(2),A%Dimension,2,2)
		TotalData=A%getTotalData()
		call U%pointDim(dim)
		call U%Data%allocateDataArrayMomery(product(dim),TotalData,A%getType())
		call V%pointDim(dim)
		call V%Data%allocateDataArrayMomery(product(dim),TotalData,A%getType())
		call V%pointDeg(dim,1)
		select case(classtype)
			case(5)
				call S%Data%allocateDataArrayMomery(NewQN(1)%getQNlength(),sum(NewQN(1)%getDeg()),3)
			case(4)
				call S%Data%allocateDataArrayMomery(NewQN(1)%getQNlength(),sum(NewQN(1)%getDeg()),2)
			case default
				call S%Data%allocateDataArrayMomery(NewQN(1)%getQNlength(),sum(NewQN(1)%getDeg()),A%getType())
		end select
		call A%pointDim(dim)
		call A%pointQN(QN1,1)
		call A%pointQN(QN2,2)
		testi=0
		testj=0
		do j=1,dim(2)
			do i=1,dim(1)
				if(A%getFlag([i,j]).and.ifNonZeroBlockDATATYPE(A,i,j))then
					if(i.lt.testi)then ! U(1), S1=1.5 S2=0.5, the diag mement is not (i,i) but (i+1,i)
						call writemess("ERROR in SVD on SymTensor",-1)
						call writemess('The input matrix is not diag',-1)
						call writemess("The quantum number are:"+QN1(i)+' , '+QN2(j),-1)
						call A%diminfo(.true.)
						call error_stop()
					end if
					if(j.lt.testj)then ! U(1), S1=1.5 S2=0.5, the diag mement is not (i,i) but (i+1,i)
						call writemess("ERROR in SVD on SymTensor",-1)
						call writemess('The input matrix is not diag',-1)
						call writemess("The quantum number are:"+QN1(i)+' , '+QN2(j),-1)
						call A%diminfo(.true.)
						call error_stop()
					end if
					testi=i
					testj=j
					indexs=NewQN(1)%getIndex(QN2(j))
					if(indexs.eq.0)then
						call writemess('ERRRO in SVD',-1)
						call NewQN(1)%print
						call writemess('j='+j)
						call writemess('QN2(j)='+QN2(j))
						call error_stop
					end if
					call A%pointer(Ap,[i,j])
					call U%setBlockMomery([i,indexs])
					call U%pointer(Up,[i,indexs])
					call V%setBlockMomery([indexs,j])
					call V%pointer(Vp,[indexs,j])
					m=size(Ap,1)
					n=size(Ap,2)
					min_MN=min(m,n)
					call S%Data%Set_block_momery(indexs,min_MN)
					call S%Data%pointer(Sp,indexs)
					call SVDDataSubroutineDATATYPE(Ap,Up,Sp,Vp,m,n,min_MN,info,classtype)
				end if
			end do
		end do
		if(present(CUTOFF))then
			DegLen=U%dim(2)
			STotalData=S%getTotalData()
			call WorkingMemory%allocate(1,STotalData)
			call WorkingMemory%allocate(classtype,STotalData)
			call WorkingMemory%get_memory(AllSindex,STotalData)
			call WorkingMemory%get_memory(AllSdata,STotalData)
			call S%Data%pointAllData(sp)
			AllSdata=Sp
			call allocateCheck(SVDCUTDeg,DegLen)
			call SVDSavingSingleValueDegDATATYPE(S,SVDCUTDeg(1:DegLen),CUTOFF,AllSindex,AllSdata)
			call WorkingMemory%free()
			call U%CutTensorDeg(2,SVDCUTDeg(1:DegLen),inoutU)
			call V%CutTensorDeg(1,SVDCUTDeg(1:DegLen),inoutV)
			call S%Data%CutOffDataArray(SVDCUTDeg(1:DegLen),inoutS%Data)
		end if
		return
	end subroutine

	subroutine SVDMatrixKillValue2DATATYPE(A,inoutU,inoutS,inoutV,WT1,WT2,WT3,minNumSave,maxNumSave&
		,inmaxValue,VType,outNumSave,outmaxValue)
		class(Tensor),intent(inout)::A
		type(Tensor),target,intent(inout)::inoutU,inoutS,inoutV,WT1,WT2,WT3
		integer,intent(in)::minNumSave,maxNumSave
		real*8,intent(in)::inmaxValue
		character(len=*),intent(in)::VType
		integer,optional,intent(inout)::outNumSave
		real*8,optional,intent(inout)::outmaxValue
		real*8::maxValue
		DATATYPE,pointer::Ap(:,:),Up(:,:),Vp(:,:)
		DATATYPE2,pointer::sp(:),AllSdata(:)
		integer::M,N,min_MN,classtype,info,Newdimi,NumSave
		type(Dimension)::UDim,VDim
		type(QuanNum)::NewQN(2)
		integer::TotalData
		integer,pointer::dim(:)
		real*4,pointer::QN1(:),QN2(:)
		integer::i,j,indexs,STotalData,DegLen,NewDimLen,testi,testj
		integer,pointer::AllSindex(:)
		type(Tensor),pointer::U,S,V
		classtype=A%getType()
		if(A%getRank().ne.2)then
			call writemess('ERROR in SVD, the input tensor should be a matrix',-1)
			call error_Stop
		end if
		call inoutS%empty()
		U=>WT1
		S=>WT2
		V=>WT3
		if(.not.A%getSymmetryFlag())then
			m=A%dim(1)
			n=A%dim(2)
			min_MN=min(M,N)
			Newdimi=min_MN
			call pasteDimension(UDim,A%Dimension,1,1,[Newdimi])
			call pasteDimension(VDim,[Newdimi],A%Dimension,2,2)
			call U%allocate(UDim,classtype)
			call V%Allocate(VDim,classtype)
			select case(classtype)
				case(5)
					call s%allocate([Newdimi],3)
				case(4)
					call s%allocate([Newdimi],2)
				case default
					call s%allocate([Newdimi],classtype)
			end select
			call A%pointer(Ap)
			call U%pointer(Up)
			call V%pointer(Vp)
			call S%pointer(Sp)
			call SVDDataSubroutineDATATYPE(Ap,Up,Sp,Vp,m,n,min_MN,info,classtype)
			if(info.ne.0) then
				call writemess('Error in svd ,info='+info,-1)
				call writemess('output The data in ./_SVD_ERROR_LOG.err',-1)
				open(unit=9991,file='./_SVD_ERROR_LOG.err',STATUS='replace',POSITION='APPEND')
				call A%write(9991)
				call U%write(9991)
				call S%write(9991)
				call V%write(9991)
				close(9991)
				call error_stop()
			end if
			STotalData=S%getTotalData()
			call WorkingMemory%allocate(classtype,STotalData)
			call WorkingMemory%get_memory(AllSdata,STotalData)
			AllSdata=Sp
			call SVDSavingSingleValueNonSymmetryDATATYPE(NewDimLen,minNumSave,maxNumSave,&
									NumSave,inmaxValue,VType,AllSdata,maxValue)
			call WorkingMemory%free()
			call U%pointDim(dim)
			call U%CutTensorDim([dim(1),NewDimLen],inoutU)
			call V%pointDim(dim)
			call V%CutTensorDim([NewDimLen,dim(2)],inoutV)
			call S%CutTensorDim([NewDimLen],inoutS)
			if(present(outNumSave))outNumSave=NumSave
			if(present(outNumSave))outmaxValue=maxValue
			return
		end if

		call choosSameQN(NewQN,A%Dimension)
		call pasteDimension(U%Dimension,A%Dimension,1,1,NewQN(1))
		call pasteDimension(V%Dimension,NewQN(2),A%Dimension,2,2)
		TotalData=A%getTotalData()
		call U%pointDim(dim)
		call U%Data%allocateDataArrayMomery(product(dim),TotalData,A%getType())
		call V%pointDim(dim)
		call V%Data%allocateDataArrayMomery(product(dim),TotalData,A%getType())
		call V%pointDeg(dim,1)
		select case(classtype)
			case(5)
				call S%Data%allocateDataArrayMomery(NewQN(1)%getQNlength(),sum(NewQN(1)%getDeg()),3)
			case(4)
				call S%Data%allocateDataArrayMomery(NewQN(1)%getQNlength(),sum(NewQN(1)%getDeg()),2)
			case default
				call S%Data%allocateDataArrayMomery(NewQN(1)%getQNlength(),sum(NewQN(1)%getDeg()),A%getType())
		end select
		call A%pointDim(dim)
		call A%pointQN(QN1,1)
		call A%pointQN(QN2,2)
		testi=0
		testj=0
		do j=1,dim(2)
			do i=1,dim(1)
				if(A%getFlag([i,j]).and.ifNonZeroBlockDATATYPE(A,i,j))then
					if(i.lt.testi)then ! U(1), S1=1.5 S2=0.5, the diag mement is not (i,i) but (i+1,i)
						call writemess("ERROR in SVD on SymTensor",-1)
						call writemess('The input matrix is not diag',-1)
						call writemess("The quantum number are:"+QN1(i)+' , '+QN2(j),-1)
						call A%diminfo(.true.)
						call error_stop()
					end if
					if(j.lt.testj)then ! U(1), S1=1.5 S2=0.5, the diag mement is not (i,i) but (i+1,i)
						call writemess("ERROR in SVD on SymTensor",-1)
						call writemess('The input matrix is not diag',-1)
						call writemess("The quantum number are:"+QN1(i)+' , '+QN2(j),-1)
						call A%diminfo(.true.)
						call error_stop()
					end if
					testi=i
					testj=j
					indexs=NewQN(1)%getIndex(QN2(j))
					if(indexs.eq.0)then
						call writemess('ERRRO in SVD',-1)
						call error_stop
					end if
					call A%pointer(Ap,[i,j])
					call U%setBlockMomery([i,indexs])
					call U%pointer(Up,[i,indexs])
					call V%setBlockMomery([indexs,j])
					call V%pointer(Vp,[indexs,j])
					m=size(Ap,1)
					n=size(Ap,2)
					min_MN=min(m,n)
					call S%Data%Set_block_momery(indexs,min_MN)
					call S%Data%pointer(Sp,indexs)
					call SVDDataSubroutineDATATYPE(Ap,Up,Sp,Vp,m,n,min_MN,info,classtype)
				end if
			end do
		end do

		DegLen=U%dim(2)
		STotalData=S%getTotalData()
		call WorkingMemory%allocate(1,STotalData)
		call WorkingMemory%allocate(classtype,STotalData)
		call WorkingMemory%get_memory(AllSindex,STotalData)
		call WorkingMemory%get_memory(AllSdata,STotalData)
		call S%Data%pointAllData(sp)
		AllSdata=Sp
		call allocateCheck(SVDCUTDeg,DegLen)
		call SVDSavingSingleValueDegValueDATATYPE(S,SVDCUTDeg(1:DegLen),minNumSave,&
					maxNumSave,NumSave,inmaxValue,VType,AllSindex,AllSdata,maxValue)
		call WorkingMemory%free()
		call U%CutTensorDeg(2,SVDCUTDeg(1:DegLen),inoutU)
		call V%CutTensorDeg(1,SVDCUTDeg(1:DegLen),inoutV)
		call S%Data%CutOffDataArray(SVDCUTDeg(1:DegLen),inoutS%Data)
		if(present(outNumSave))outNumSave=NumSave
		if(present(outmaxValue))outmaxValue=maxValue
		return
	end subroutine

	