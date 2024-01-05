	subroutine allocateCheck_cha_len(A,lenA,chalen)
		character(len=:),allocatable,intent(inout)::A(:)
		integer,intent(in)::lenA,chalen
		if(allocated(A)) then
			if(chalen.eq.len(A))then
				if(size(A).ltne.lenA) then
					deallocate(A)
					allocate(character(len=chalen)::A(lenA))
				end if
			else
				deallocate(A)
				allocate(character(len=chalen)::A(lenA))
			end if
		else
			allocate(character(len=chalen)::A(lenA))
		end if
		return
	end subroutine
	subroutine allocateDataArrayClassData(A,length,iType,chalen)
		type(DataArray),intent(inout)::A
		integer,intent(in)::length,iType
		integer,optional,intent(in)::chalen
		if(iType.eq.8)then
			if(present(chalen))then
				A%DataCharacterLen=chalen
			else
				if(A%DataCharacterLen.le.0)then
					A%DataCharacterLen=characterlen
				end if
			end if
		else
			A%DataCharacterLen=0
		end if
		if(iType.eq.7)then
			A%DataCharacterLen=characterlen
		end if
		select case(iType)
			case(1)
				call allocateCheck(A%idata,length)
				call ClassPointer1DFunc(length,A%idata,A%ClassData)
			case(2)
				call allocateCheck(A%sdata,length)
				call ClassPointer1DFunc(length,A%sdata,A%ClassData)
			case(3)
				call allocateCheck(A%ddata,length)
				call ClassPointer1DFunc(length,A%ddata,A%ClassData)
			case(4)
				call allocateCheck(A%cdata,length)
				call ClassPointer1DFunc(length,A%cdata,A%ClassData)
			case(5)
				call allocateCheck(A%zdata,length)
				call ClassPointer1DFunc(length,A%zdata,A%ClassData)
			case(6)
				call allocateCheck(A%ldata,length)
				call ClassPointer1DFunc(length,A%ldata,A%ClassData)
			case(7)
				call allocateCheck(A%adata,length)
				call ClassPointer1DFunc(length,A%adata,A%ClassData)
			case(8)
				call allocateCheck_cha_len(A%adata2,length,A%DataCharacterLen)
				call ClassPointer1DFunc(length,A%adata2,A%ClassData)
			case(0)
				if(length.ne.0)then
					call writemess('ERROR in allocateDataArrayClassData, error type')
					call writemess('length='+length)
					call writemess('iType='+iType)
					call error_stop
				endif
			case default
				call allocate_Class_data(A%ClsData,length,iType)
				call ClassPointer1DFunc(length,A%ClsData,A%ClassData)
		end select
		A%MomeryLength=length
		return
	end subroutine

	subroutine allocateDataArray1(DA,TotalData,classType,chalen)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalData
		integer,intent(in)::classType
		integer,optional,intent(in)::chalen
		if(TotalData.lt.0)then
			call writemess('ERROR in allocate dataarray',-1)
			call writemess('length of the data='+TotalData,-1)
			call error_stop
		end if	
		if(.not.Da%DynamicClass)then
			if((DA%classType.ne.0).and.(classType.ne.Da%classType))then
				call writemess('Can NOT change the  class type of the tensor',-1)
				call writemess('inout classType='+classType,-1)
				call writemess('      classType='+Da%classType,-1)
				call error_stop
			end if
		end if
		DA%TotalData=TotalData
		DA%totalBlock=1
		DA%classType=classType
		call allocateCheck(Da%starti,1)
		call allocateCheck(Da%endi,1)
		Da%starti(1)=1
		Da%endi(1)=TotalData
		call allocateDataArrayClassData(Da,TotalData,classType,chalen)
		return
	end subroutine

	subroutine allocateDataArray2(DA,BlockTotalData,classType,chalen)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::BlockTotalData(:)
		integer,intent(in)::classType
		integer::TotalBlock,i,ith,jth,TotalData
		integer,optional,intent(in)::chalen
		TotalBlock=size(BlockTotalData)
		TotalData=sum(BlockTotalData)
		if(TotalData.lt.0)then
			call writemess('ERROR in allocate dataarray',-1)
			call writemess('length of the data='+TotalData,-1)
			call error_stop
		end if		
		if(.not.Da%DynamicClass)then
			if((DA%classType.ne.0).and.(classType.ne.Da%classType))then
				call writemess('Can NOT change the  class type of the tensor',-1)
				call writemess('inout classType='+classType,-1)
				call writemess('      classType='+Da%classType,-1)
				call error_stop
			end if
		end if
		DA%TotalData=TotalData
		DA%totalBlock=TotalBlock
		DA%classType=classType
		call allocateCheck(Da%starti,TotalBlock)
		call allocateCheck(Da%endi,TotalBlock)
		call allocateDataArrayClassData(Da,TotalData,classType,chalen)
		ith=0
		jth=0
		do i=1,TotalBlock
			if(BlockTotalData(i).gt.0)then
				ith=jth+1
				jth=jth+BlockTotalData(i)
				if(jth.gt.TotalData)then
					call writemess('ERROR in allocateDataArray2',-1)
					call error_stop
				end if
				Da%starti(i)=ith
				Da%endi(i)=jth
			else
				Da%starti(i)=1
				Da%endi(i)=0
			end if
		end do
		return
	end subroutine

	subroutine allocateDataArray3(DA,DaType)
		class(DataArray),intent(inout)::Da
		type(DataArray),intent(in)::DaType
		integer::chalen
		call Da%empty()
		if(.not.DaType%getFlag())return
		if(.not.Da%DynamicClass)then
			if((DA%classType.ne.0).and.(DaType%classType.ne.Da%classType))then
				call writemess('Can NOT change the  class type of the tensor',-1)
				call writemess('inout classType='+DaType%classType,-1)
				call writemess('      classType='+Da%classType,-1)
				call error_stop
			end if
		end if
		DA%TotalData=DaType%TotalData
		DA%totalBlock=DaType%totalBlock
		DA%classType=DaType%classType
		if(DaType%classType.eq.7)then
			chalen=DaType%getCharacterlen()
		else
			chalen=0
		end if
		call allocateCheck(Da%starti,DaType%totalBlock)
		call allocateCheck(Da%endi,DaType%totalBlock)
		Da%starti(1:DaType%totalBlock)=DaType%starti(1:DaType%totalBlock)
		Da%endi(1:DaType%totalBlock)=DaType%endi(1:DaType%totalBlock)
		call allocateDataArrayClassData(Da,DA%TotalData,DA%classType,chalen)
		return
	end subroutine
	subroutine allocateDataArray4(DA,DaType,ClassType,chalen)
		class(DataArray),intent(inout)::Da
		type(DataArray),intent(in)::DaType
		integer,intent(in)::ClassType
		integer,optional,intent(in)::chalen
		call Da%empty()
		if(.not.DaType%getFlag())return
		if(.not.Da%DynamicClass)then
			if((DA%classType.ne.0).and.(classType.ne.Da%classType))then
				call writemess('Can NOT change the  class type of the tensor',-1)
				call writemess('inout classType='+classType,-1)
				call writemess('      classType='+Da%classType,-1)
				call error_stop
			end if
		end if
		DA%TotalData=DaType%TotalData
		DA%totalBlock=DaType%totalBlock
		DA%classType=ClassType
		call allocateCheck(Da%starti,DaType%totalBlock)
		call allocateCheck(Da%endi,DaType%totalBlock)
		Da%starti(1:DaType%totalBlock)=DaType%starti(1:DaType%totalBlock)
		Da%endi(1:DaType%totalBlock)=DaType%endi(1:DaType%totalBlock)
		call allocateDataArrayClassData(Da,DA%TotalData,DA%classType,chalen)
		return
	end subroutine
	subroutine allocateDataArray5(DA,TotalData)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalData
		integer::chalen
		if(TotalData.lt.0)then
			call writemess('ERROR in allocate dataarray',-1)
			call writemess('length of the data='+TotalData,-1)
			call error_stop
		end if	
		if(DA%classType.eq.0)then
			call writemess('DO NOT set class type to the data array yet',-1)
			call error_stop
		end if
		if(DA%classType.eq.7)then
			chalen=DA%getCharacterlen()
		else
			chalen=0
		end if
		DA%TotalData=TotalData
		DA%totalBlock=1
		call allocateCheck(Da%starti,1)
		call allocateCheck(Da%endi,1)
		Da%starti(1)=1
		Da%endi(1)=TotalData
		call allocateDataArrayClassData(Da,TotalData,DA%classType,chalen)
		return
	end subroutine
	subroutine allocateDataArray6(DA,BlockTotalData)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::BlockTotalData(:)
		integer::TotalBlock,i,ith,jth,TotalData
		integer::chalen
		TotalBlock=size(BlockTotalData)
		TotalData=sum(BlockTotalData)
		if(TotalData.lt.0)then
			call writemess('ERROR in allocate dataarray',-1)
			call writemess('length of the data='+TotalData,-1)
			call error_stop
		end if		
		if(DA%classType.eq.0)then
			call writemess('DO NOT set class type to the data array yet',-1)
			call error_stop
		end if
		if(DA%classType.eq.7)then
			chalen=DA%getCharacterlen()
		else
			chalen=0
		end if
		DA%TotalData=TotalData
		DA%totalBlock=TotalBlock
		call allocateCheck(Da%starti,TotalBlock)
		call allocateCheck(Da%endi,TotalBlock)
		call allocateDataArrayClassData(Da,TotalData,DA%classType,chalen)
		ith=0
		jth=0
		do i=1,TotalBlock
			if(BlockTotalData(i).gt.0)then
				ith=jth+1
				jth=jth+BlockTotalData(i)
				if(jth.gt.TotalData)then
					call writemess('ERROR in allocateDataArray2',-1)
					call error_stop
				end if
				Da%starti(i)=ith
				Da%endi(i)=jth
			else
				Da%starti(i)=1
				Da%endi(i)=0
			end if
		end do
		return
	end subroutine

	subroutine allocateDataArrayChar1(DA,TotalData,classType_,chalen)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalData
		character(len=*),intent(in)::classType_
		integer,optional,intent(in)::chalen
		integer::classType
		ClassType=select_data_type_char(ClassType_)
		call allocateDataArray1(DA,TotalData,classType,chalen)
	end subroutine

	subroutine allocateDataArrayChar2(DA,TotalData,classType_,chalen)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalData(:)
		character(len=*),intent(in)::classType_
		integer,optional,intent(in)::chalen
		integer::classType
		ClassType=select_data_type_char(ClassType_)
		call allocateDataArray2(DA,TotalData,classType,chalen)
		return
	end subroutine

	subroutine allocateDataArrayChar4(DA,DaType,ClassType_,chalen)
		class(DataArray),intent(inout)::Da
		type(DataArray),intent(in)::DaType
		character(len=*),intent(in)::classType_
		integer,optional,intent(in)::chalen
		integer::classType
		ClassType=select_data_type_char(ClassType_)
		call allocateDataArray4(DA,DaType,classType,chalen)
		return
	end subroutine


	subroutine allocateDataArray1Check(DA,TotalData,classType,chalen)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalData
		integer,intent(in)::classType
		integer,optional,intent(in)::chalen
		if((Da%DynamicClass).or.(Da%classType.eq.0))then
			Da%classType=classType
		end if
		if(present(chalen))DA%DataCharacterLen=chalen
		call allocateDataArray5(DA,TotalData)
		return
	end subroutine
	subroutine allocateDataArray2Check(DA,BlockTotalData,classType,chalen)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::BlockTotalData(:)
		integer,intent(in)::classType
		integer,optional,intent(in)::chalen
		if((Da%DynamicClass).or.(Da%classType.eq.0))then
			Da%classType=classType
		end if
		if(present(chalen))DA%DataCharacterLen=chalen
		call allocateDataArray6(DA,BlockTotalData)
		return
	end subroutine
	subroutine allocateDataArray3Check(DA,DaType)
		class(DataArray),intent(inout)::Da
		type(DataArray),intent(in)::DaType
		integer::chalen
		call Da%empty()
		if(.not.DaType%getFlag())return
		if((Da%DynamicClass).or.(Da%classType.eq.0))then
			Da%classType=DaType%classType
		end if
		if(DaType%classType.eq.7)then
			chalen=DaType%getCharacterlen()
		else
			chalen=0
		end if
		DA%TotalData=DaType%TotalData
		DA%totalBlock=DaType%totalBlock
		call allocateCheck(Da%starti,DaType%totalBlock)
		call allocateCheck(Da%endi,DaType%totalBlock)
		Da%starti(1:DaType%totalBlock)=DaType%starti(1:DaType%totalBlock)
		Da%endi(1:DaType%totalBlock)=DaType%endi(1:DaType%totalBlock)
		call allocateDataArrayClassData(Da,DA%TotalData,DA%classType,chalen)
		return
	end subroutine
	subroutine allocateDataArray4Check(DA,DaType,ClassType,chalen)
		class(DataArray),intent(inout)::Da
		type(DataArray),intent(in)::DaType
		integer,intent(in)::ClassType
		integer,optional,intent(in)::chalen
		call Da%empty()
		if(.not.DaType%getFlag())return
		if((Da%DynamicClass).or.(Da%classType.eq.0))then
			Da%classType=ClassType
		end if
		DA%TotalData=DaType%TotalData
		DA%totalBlock=DaType%totalBlock
		call allocateCheck(Da%starti,DaType%totalBlock)
		call allocateCheck(Da%endi,DaType%totalBlock)
		Da%starti(1:DaType%totalBlock)=DaType%starti(1:DaType%totalBlock)
		Da%endi(1:DaType%totalBlock)=DaType%endi(1:DaType%totalBlock)
		call allocateDataArrayClassData(Da,DA%TotalData,DA%classType,chalen)
		return
	end subroutine
	subroutine allocateDataArrayChar1Check(DA,TotalData,classType_,chalen)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalData
		character(len=*),intent(in)::classType_
		integer,optional,intent(in)::chalen
		integer::classType
		ClassType=select_data_type_char(ClassType_)
		if((Da%DynamicClass).or.(Da%classType.eq.0))then
			Da%classType=classType
		end if
		if(present(chalen))Da%DataCharacterLen=chalen
		call allocateDataArray5(DA,TotalData)
	end subroutine

	subroutine allocateDataArrayChar2Check(DA,TotalData,classType_,chalen)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalData(:)
		character(len=*),intent(in)::classType_
		integer,optional,intent(in)::chalen
		integer::classType
		ClassType=select_data_type_char(ClassType_)
		if((Da%DynamicClass).or.(Da%classType.eq.0))then
			Da%classType=classType
		end if
		if(present(chalen))Da%DataCharacterLen=chalen
		call allocateDataArray6(DA,TotalData)
		return
	end subroutine

	subroutine allocateDataArray_no_empty_Block(DA,BlockTotalData,classType,chalen)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::BlockTotalData(:)
		integer,intent(in)::classType
		integer,optional,intent(in)::chalen
		integer::TotalBlock,i,ith,jth,TotalData,ii
		TotalBlock=0
		do i=1,size(BlockTotalData)
			if(BlockTotalData(i).ne.0)TotalBlock=TotalBlock+1
		end do
		TotalData=sum(BlockTotalData)
		if(TotalData.eq.0)then
			call Da%empty()
			return
		end if
		if(TotalData.lt.0)then
			call writemess('ERROR in allocate dataarray',-1)
			call writemess('length of the data='+TotalData,-1)
			call error_stop
		end if	
		if(.not.Da%DynamicClass)then
			if((DA%classType.ne.0).and.(classType.ne.Da%classType))then
				call writemess('Can NOT change the  class type of the tensor',-1)
				call writemess('inout classType='+classType,-1)
				call writemess('      classType='+Da%classType,-1)
				call error_stop
			end if
		end if	
		DA%TotalData=TotalData
		DA%totalBlock=TotalBlock
		DA%classType=classType
		call allocateCheck(Da%starti,TotalBlock)
		call allocateCheck(Da%endi,TotalBlock)
		call allocateDataArrayClassData(Da,TotalData,classType,chalen)
		ith=0
		jth=0
		ii=0
		do i=1,size(BlockTotalData)
			if(BlockTotalData(i).gt.0)then
				ii=ii+1
				ith=jth+1
				jth=jth+BlockTotalData(i)
				if(jth.gt.TotalData)then
					call writemess('ERROR in allocateDataArray2',-1)
					call error_stop
				end if
				Da%starti(ii)=ith
				Da%endi(ii)=jth
			end if
		end do
		return
	end subroutine
