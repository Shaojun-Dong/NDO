	subroutine setDefineClassType(Da,defineClassType)
		class(DataArray),intent(inout)::Da
		integer,optional,intent(in)::defineClassType
		if(present(defineClassType))then
			if((defineClassType.ge.(0)).and.(defineClassType.le.8))then
				call writemess('Can not label the datatype='+defineClassType,-1)
				call writemess('datatype=0~8 have been defined',-1)
				call error_stop
			end if
			DA%classType=defineClassType
		else
			DA%classType=-1
		end if
		return
	end subroutine
	subroutine allocate_From_source1(Da,TotalData,source,defineClassType)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalData
		class(*),intent(in)::source
		integer,optional,intent(in)::defineClassType
		if(TotalData.lt.0)then
			call writemess('ERROR in allocate dataarray',-1)
			call writemess('length of the data='+TotalData,-1)
			call error_stop
		end if	
		if(.not.Da%DynamicClass)then
			if((DA%classType.ne.0))then
				call writemess('Can NOT change the  class type of the tensor',-1)
				call writemess('      classType='+Da%classType,-1)
				call error_stop
			end if
		end if
		DA%TotalData=TotalData
		DA%totalBlock=1
		call setDefineClassType(Da,defineClassType)
		call allocateCheck(Da%starti,1)
		call allocateCheck(Da%endi,1)
		Da%starti(1)=1
		Da%endi(1)=TotalData
		if(allocated(Da%clsData))deallocate(Da%clsData)
		allocate(Da%clsData(TotalData),mold=source)
		call ClassPointer1DFunc(TotalData,Da%clsData,Da%ClassData)
		Da%MomeryLength=TotalData
		return
	end subroutine

	subroutine allocate_From_source2(DA,BlockTotalData,source,defineClassType)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::BlockTotalData(:)
		class(*),intent(in)::source
		integer,optional,intent(in)::defineClassType
		integer::TotalBlock,i,ith,jth,TotalData
		TotalBlock=size(BlockTotalData)
		TotalData=sum(BlockTotalData)
		if(TotalData.lt.0)then
			call writemess('ERROR in allocate dataarray',-1)
			call writemess('length of the data='+TotalData,-1)
			call error_stop
		end if		
		if(.not.Da%DynamicClass)then
			if((DA%classType.ne.0))then
				call writemess('Can NOT change the  class type of the tensor',-1)
				call writemess('      classType='+Da%classType,-1)
				call error_stop
			end if
		end if
		DA%TotalData=TotalData
		DA%totalBlock=TotalBlock
		call setDefineClassType(Da,defineClassType)
		call allocateCheck(Da%starti,TotalBlock)
		call allocateCheck(Da%endi,TotalBlock)
		if(allocated(Da%clsData))deallocate(Da%clsData)
		allocate(Da%clsData(TotalData),mold=source)
		call ClassPointer1DFunc(TotalData,Da%clsData,Da%ClassData)
		Da%MomeryLength=TotalData
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

	subroutine allocate_From_source3(DA,DaType,source,defineClassType)
		class(DataArray),intent(inout)::Da
		type(DataArray),intent(in)::DaType
		class(*),intent(in)::source
		integer,optional,intent(in)::defineClassType
		call Da%empty()
		if(.not.DaType%getFlag())return
		if(.not.Da%DynamicClass)then
			if((DA%classType.ne.0))then
				call writemess('Can NOT change the  class type of the tensor',-1)
				call writemess('      classType='+Da%classType,-1)
				call error_stop
			end if
		end if
		DA%TotalData=DaType%TotalData
		DA%totalBlock=DaType%totalBlock
		call setDefineClassType(Da,defineClassType)
		call allocateCheck(Da%starti,DaType%totalBlock)
		call allocateCheck(Da%endi,DaType%totalBlock)
		Da%starti(1:DaType%totalBlock)=DaType%starti(1:DaType%totalBlock)
		Da%endi(1:DaType%totalBlock)=DaType%endi(1:DaType%totalBlock)
		if(allocated(Da%clsData))deallocate(Da%clsData)
		allocate(Da%clsData(DA%TotalData),mold=source)
		call ClassPointer1DFunc(DA%TotalData,Da%clsData,Da%ClassData)

		Da%MomeryLength=DA%TotalData
		return
	end subroutine


	subroutine allocateDataArrayMomerySource1(DA,TotalBlock,TotalData,source,defineClassType)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalBlock,TotalData
		class(*),intent(in)::source
		integer,optional,intent(in)::defineClassType
		if(.not.Da%DynamicClass)then
			if((DA%classType.ne.0))then
				call writemess('Can NOT change the  class type of the tensor',-1)
				call writemess('      classType='+Da%classType,-1)
				call error_stop
			end if
		end if
		DA%TotalData=0
		DA%totalBlock=TotalBlock
		call setDefineClassType(Da,defineClassType)
		call allocateCheck(Da%starti,TotalBlock)
		call allocateCheck(Da%endi,TotalBlock)
		Da%endi(1:TotalBlock)=0
		if(allocated(Da%clsData))deallocate(Da%clsData)
		allocate(Da%clsData(TotalData),mold=source)
		call ClassPointer1DFunc(TotalData,Da%clsData,Da%ClassData)
		Da%MomeryLength=TotalData
		return
	end subroutine

	subroutine allocateDataArrayMomerySource1Check(DA,TotalBlock,TotalData,source,defineClassType)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalBlock,TotalData
		class(*),intent(in)::source
		integer,optional,intent(in)::defineClassType
		if((Da%DynamicClass).or.(Da%classType.eq.0))then
			if(present(defineClassType))then
				Da%classType=defineClassType
			else
				call writemess('ERROR in allocateDataMemory',-1)
				call error_stop
			end if
		end if
		if((Da%classType.gt.0).and.(Da%classType.le.8))then
			call allocateDataArrayMomery1(DA,TotalBlock,TotalData,Da%classType)
		else
			call allocateDataArrayMomerySource1(DA,TotalBlock,TotalData,source,defineClassType)
		end if
		return
	end subroutine