	subroutine allocateFUNCNAME(A,length)
		class(*),intent(inout),allocatable::A(:)
		integer,intent(in)::length
		if(allocated(A)) then
			if(size(A).lt.length) then
				deallocate(A)
				allocate(DATATYPE::A(length))
			end if
		else
			if(length.gt.0)then
				allocate(DATATYPE::A(length))
			end if
		end if
		return
	end subroutine
