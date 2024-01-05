
	subroutine GetaValue(outA,inA,ith)
		class(*),target,intent(in)::inA(:)
		class(*),target,intent(inout)::outA
		integer,intent(in)::ith
		class(*),pointer::p
		if(classsize(inA).lt.ith)then
			call writemess('ERROR in GetaValue, ith>size',-1)
			call writemess('ith='+ith,-1)
			call writemess('size='+classsize(inA),-1)
			call error_stop
		end if
		call ClassPointer1DFunc(inA,ith,p)
		call copyData(outA,p)
		return
	end subroutine

	subroutine GetSameValue(outA,inA,ith,jth)
		class(*),target,intent(in)::inA(:)
		class(*),target,intent(inout)::outA(:)
		integer,intent(in)::ith,jth
		class(*),pointer::p(:)
		if(classsize(inA).lt.ith)then
			call writemess('ERROR in GetaValue, ith>size',-1)
			call writemess('ith='+ith,-1)
			call writemess('size='+classsize(inA),-1)
			call error_stop
		end if
		call ClassPointer1DFunc(inA,ith,jth,p)
		call FastCopyArray(outA,p,jth-ith+1)
		return
	end subroutine