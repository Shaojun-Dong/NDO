	subroutine allocateDataArrayMomery1(DA,TotalBlock,TotalData,classType,chalen)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalBlock,TotalData
		integer,intent(in)::classType
		integer,optional,intent(in)::chalen
		if(.not.Da%DynamicClass)then
			if((DA%classType.ne.0).and.(classType.ne.Da%classType))then
				call writemess('Can NOT change the  class type of the tensor',-1)
				call writemess('inout classType='+classType,-1)
				call writemess('      classType='+Da%classType,-1)
				call error_stop
			end if
		end if
		DA%TotalData=0
		DA%totalBlock=TotalBlock
		DA%classType=classType
		call allocateCheck(Da%starti,TotalBlock)
		call allocateCheck(Da%endi,TotalBlock)
		Da%endi(1:TotalBlock)=0
		call allocateDataArrayClassData(Da,TotalData,DA%classType,chalen)
		return
	end subroutine
	subroutine allocateDataArrayMomery2(DA,TotalBlock,TotalData,classType,chalen)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalBlock,TotalData
		character(len=*),intent(in)::classType
		integer,optional,intent(in)::chalen
		if(.not.Da%DynamicClass)then
			if((DA%classType.ne.0).and.(select_data_type_char(classType).ne.Da%classType))then
				call writemess('Can NOT change the  class type of the tensor',-1)
				call writemess('inout classType='+classType,-1)
				call writemess('      classType='+Da%classType,-1)
				call error_stop
			end if
		end if
		DA%TotalData=0
		DA%totalBlock=TotalBlock
		DA%classType=select_data_type_char(classType)
		call allocateCheck(Da%starti,TotalBlock)
		call allocateCheck(Da%endi,TotalBlock)
		Da%endi(1:TotalBlock)=0
		call allocateDataArrayClassData(Da,TotalData,DA%classType,chalen)
		return
	end subroutine
	subroutine allocateDataArrayMomery3(DA,TotalBlock,TotalData)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalBlock,TotalData
		if(Da%classType.eq.0)then
			call writemess('DO not set class type to the tensor yet',-1)
			call error_stop
		end if
		DA%TotalData=0
		DA%totalBlock=TotalBlock
		call allocateCheck(Da%starti,TotalBlock)
		call allocateCheck(Da%endi,TotalBlock)
		Da%endi(1:TotalBlock)=0
		call allocateDataArrayClassData(Da,TotalData,DA%classType)
		return
	end subroutine

	subroutine allocateDataArrayMomery1Check(DA,TotalBlock,TotalData,classType)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalBlock,TotalData
		integer,intent(in)::classType
		if((Da%DynamicClass).or.(Da%classType.eq.0))then
			Da%classType=classType
		end if
		call allocateDataArrayMomery1(DA,TotalBlock,TotalData,Da%classType)
		return
	end subroutine
	subroutine allocateDataArrayMomery2Check(DA,TotalBlock,TotalData,classType)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::TotalBlock,TotalData
		character(len=*),intent(in)::classType
		if((Da%DynamicClass).or.(Da%classType.eq.0))then
			Da%classType=select_data_type_char(classType)
		end if
		call allocateDataArrayMomery1(DA,TotalBlock,TotalData,Da%classType)
		return
	end subroutine

	subroutine Set_block_momery1(Da,ith,length)
		class(DataArray),intent(inout)::Da
		integer,intent(in)::ith,length
		if((DA%TotalData+length).gt.Da%MomeryLength)then
			call writemess('ERROR in Set_block_momery, totalData',-1)
			call writemess('DA%TotalData='+DA%TotalData)
			call writemess('DA%TotalData+length='+(DA%TotalData+length))
			call writemess('Da%MomeryLength='+Da%MomeryLength)
			call error_stop
		end if
		if(ith.gt.Da%TotalBlock)then
			call writemess('ERROR in Set_block_momery',-1)
			call writemess('Da%TotalBlock='+Da%TotalBlock,-1)
			call writemess('ith='+ith,-1) 
			call error_stop
		end if
		if(Da%endi(ith).ne.0)then
			call writemess('The block has been allocated',-1)
			call error_stop
		end if
		Da%starti(ith)=DA%TotalData+1
		Da%endi(ith)=Da%TotalData+length
		Da%TotalData=Da%endi(ith)
		return
	end subroutine
	subroutine Set_block_momery2(Da,dim,vec,length)
		class(DataArray),target,intent(inout)::Da
		integer,intent(in)::dim(:),vec(:),length
		integer,pointer::ei2(:,:),si2(:,:)
		integer,pointer::ei3(:,:,:),si3(:,:,:)
		integer,pointer::ei4(:,:,:,:),si4(:,:,:,:)
		integer::index
		if((DA%TotalData+length).gt.Da%MomeryLength)then
			call writemess('ERROR in Set_block_momery, totalData',-1)
			call writemess('DA%TotalData='+(DA%TotalData),-1)
			call writemess('length='+(length),-1)
			call writemess('DA%TotalData+length='+(DA%TotalData+length),-1)
			call writemess('Da%MomeryLength='+Da%MomeryLength,-1)
			call error_stop
		end if
		select case(size(vec))
			case (1)
				if(Da%endi(vec(1)).ne.0)then
					call writemess('The block has been allocated',-1)
					call error_stop
				end if
				Da%starti(vec(1))=DA%TotalData+1
				Da%endi(vec(1))=Da%TotalData+length
				Da%TotalData=Da%endi(vec(1))
			case (2)
				si2(1:dim(1),1:dim(2))=>Da%starti(1:Da%totalBlock)
				ei2(1:dim(1),1:dim(2))=>Da%endi(1:Da%totalBlock)
				if(ei2(vec(1),vec(2)).ne.0)then
					call writemess('The block has been allocated',-1)
					call writemess('ei2(vec(1),vec(2))='+ei2(vec(1),vec(2)))
					call error_stop
				end if
				si2(vec(1),vec(2))=Da%TotalData+1
				ei2(vec(1),vec(2))=Da%TotalData+length
				Da%TotalData=ei2(vec(1),vec(2))
			case (3)
				si3(1:dim(1),1:dim(2),1:dim(3))=>Da%starti(1:Da%totalBlock)
				ei3(1:dim(1),1:dim(2),1:dim(3))=>Da%endi(1:Da%totalBlock)
				if(ei3(vec(1),vec(2),vec(3)).ne.0)then
					call writemess('The block has been allocated',-1)
					call error_stop
				end if
				si3(vec(1),vec(2),vec(3))=Da%TotalData+1
				ei3(vec(1),vec(2),vec(3))=Da%TotalData+length
				Da%TotalData=ei3(vec(1),vec(2),vec(3))
			case (4)
				si4(1:dim(1),1:dim(2),1:dim(3),1:dim(4))=>Da%starti(1:Da%totalBlock)
				ei4(1:dim(1),1:dim(2),1:dim(3),1:dim(4))=>Da%endi(1:Da%totalBlock)
				if(ei4(vec(1),vec(2),vec(3),vec(4)).ne.0)then
					call writemess('The block has been allocated',-1)
					call error_stop
				end if
				si4(vec(1),vec(2),vec(3),vec(4))=Da%TotalData+1
				ei4(vec(1),vec(2),vec(3),vec(4))=Da%TotalData+length
				Da%TotalData=ei4(vec(1),vec(2),vec(3),vec(4))
			case default
				index=addressToIndes(vec,dim)
				if(Da%endi(index).ne.0)then
					call writemess('The block has been allocated',-1)
					call error_stop
				end if
				Da%starti(index)=DA%TotalData+1
				Da%endi(index)=Da%TotalData+length
				Da%TotalData=Da%endi(index)
		end select
		return
	end subroutine

	subroutine resetTotalData(Da,TotalData)
		class(DataArray),target,intent(inout)::Da
		integer,intent(in)::TotalData
		if(Da%MomeryLength.lt.TotalData)then
			call writemess('ERROR in resetTotalData',-1)
			call writemess('Da%MomeryLength='+Da%MomeryLength,-1)
			call writemess('TotalData='+TotalData,-1)
			call error_stop
		end if
		DA%TotalData=TotalData
		DA%totalBlock=1
		Da%starti(1)=1
		Da%endi(1)=TotalData
		return
	end subroutine
