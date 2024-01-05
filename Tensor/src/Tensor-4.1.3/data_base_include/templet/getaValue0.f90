	function getaDataValueFUNCNAME(Da,ith)result(Res)
		DATATYPE::Res
		class(DataArray),intent(in)::Da
		integer,intent(in)::ith
		if(ith.gt.Da%getTotalData())then
			call writemess('ERROR in getting the ith element of the dataarray',-1)
			call writemess('ith='+ith,-1)
			call writemess('Da%getTotalData()='+Da%getTotalData(),-1)
			call error_stop
		end if
		call getavalue(Res,Da%ClassData,ith)
		return
	end function

	function getAllDataValueFUNCNAME(Da)result(Res)
		DATATYPE,allocatable::Res(:)
		class(DataArray),intent(in)::Da
		allocate(REs(Da%TotalData))
		call FastCopyArray(Res,Da%ClassData,Da%TotalData)
		return
	end function