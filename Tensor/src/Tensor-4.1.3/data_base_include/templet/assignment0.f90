	subroutine assignmentDataArrayDATANAME(idata,inDa)
		DATATYPE,intent(inout)::idata(:)
		type(DataArray),intent(in)::inDa
		integer::totalData
		if(.not.inDa%getFlag())then
			call writemess('ERROR in assignmentDataArray, empty data',-1)
			call error_stop
		end if
		totalData=inDa%getTotalData()
		if(size(iData).lt.totalData)then
			call writemess('ERROR in assignmentDataArray, totalData',-1)
			call writemess('size(iData)='+size(iData))
			call writemess('totalData='+totalData)
			call error_stop
		end if
		call FastcopyARRAY(idata,inDa%ClassData,totalData)
		return
	end subroutine
	subroutine assignmentDataArray2DATANAME(inoutDa,idata)
		type(DataArray),intent(inout)::inoutDa
		DATATYPE,intent(IN)::idata(:)
		integer::totalData
		totalData=size(idata)
		if(inoutDa%DynamicClass)then
			call inoutDa%allocate(TotalData,DATACLASSTYPE)
		else
			call inoutDa%allocate(TotalData)
		end if
		call FastcopyARRAY(inoutDa%ClassData,idata,totalData)
		return
	end subroutine

	subroutine assignmentDataArrayValDATANAME(idata,inDa)
		DATATYPE,intent(inout)::idata
		type(DataArray),intent(in)::inDa
		if(.not.inDa%getFlag())then
			call writemess('ERROR in assignmentDataArray, empty data',-1)
			call error_stop
		end if
		if(inDa%getTotalData().gt.1)then
			call writemess('ERROR in assignmentDataArray to a value, totalData',-1)
			call writemess('totalData='+inDa%getTotalData())
			call error_stop
		end if
		call GetaValue(idata,inDa%ClassData,1)
		return
	end subroutine
	subroutine assignmentDataArrayVal2DATANAME(inoutDa,idata)
		type(DataArray),intent(inout)::inoutDa
		DATATYPE,intent(IN)::idata
		if(inoutDa%DynamicClass)then
			call inoutDa%allocate(1,DATACLASSTYPE)
		else
			call inoutDa%allocate(1)
		end if
		call ModifyaValue(inoutDa%ClassData,idata,1)
		return
	end subroutine
