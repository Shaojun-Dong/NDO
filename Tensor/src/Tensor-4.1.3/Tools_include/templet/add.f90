	function FUNCNAME(A,B)result(R)
		character(len=characterlen),allocatable::R(:)
		DATATYPE_1,intent(in)::A(:)
		DATATYPE_2,intent(in)::B(:)
		integer::lenA,i
		lenA=size(A)
		if(lenA.ne.size(B))then
			call writemess('ERROR in (+) for array, total length',-1)
			call error_stop
		end if
		allocate(R(lenA))
		do i=1,lenA
			R(i)=A(i)+B(i)
		end do
		return
	end function
