	subroutine DataArrayithPointerFuncName(Da,p,ith)
		class(DataArray),target,intent(in)::Da
		integer,intent(in)::ith
		DATATYPE,pointer,intent(inout)::p(:)
		if(Da%getTotalBlock().le.0)then
			call writemess('There is no block in the dataarray',-1)
			call error_stop
		end if
		if(ith.gt.Da%getTotalBlock())then
			call writemess('ERROR in pointting to the data of DataArray',-1)
			call writemess('ith='+ith,-1)
			call writemess('Da%getTotalBlock()='+Da%getTotalBlock(),-1)
			call error_stop
		end if
		if(Da%getType().eq.DATATYPENumber)then
			if(Da%getFlag(ith))then
				call Pointer1DFunc(Da%ClassData,Da%starti(ith),Da%endi(ith),p)
			else
				call writemess('ERROR in pointting to the data of DataArray',-1)
				call writemess('the Block ith is empty',-1)
				call writemess('ith='+ith,-1)
				call error_stop
			end if
		else
			call writemess('ERROR in pointting to the data of DataArray',-1)
			call writemess('Da%getType()='+Da%getType(),-1)
			call error_stop
		end if
		return
	end subroutine
	subroutine DataArrayAllPointerFuncName(Da,p)
		class(DataArray),target,intent(in)::Da
		DATATYPE,pointer,intent(inout)::p(:)
		if(Da%getTotalBlock().le.0)then
			call writemess('There is no block in the dataarray',-1)
			call error_stop
		end if
		if(Da%getType().eq.DATATYPENumber)then
			call Pointer1DFunc(Da%totalData,Da%ClassData,p)
		else
			call writemess('ERROR in pointting to the data of DataArray',-1)
			call writemess('Da%getType()='+Da%getType(),-1)
			call error_stop
		end if
		return
	end subroutine
