	subroutine char2arrayNAME(A,B)
		DATATYPE,intent(inout)::A(:)
		character(len=*),intent(in)::B
		integer::i
		do i=1,size(A)
			A(i)=B
		end do
		return
	end subroutine
	subroutine array2charNAME(A,B)
		character(len=*),intent(inout)::A(:)
		DATATYPE,intent(in)::B
		integer::i
		do i=1,size(A)
			A(i)=B
		end do
		return
	end subroutine

	subroutine logi2arrayNAME(A,B)
		DATATYPE,intent(inout)::A(:)
		logical,intent(in)::B
		integer::i
		do i=1,size(A)
			A(i)=B
		end do
		return
	end subroutine
	subroutine array2logiNAME(A,B)
		logical,intent(inout)::A(:)
		DATATYPE,intent(in)::B
		integer::i
		do i=1,size(A)
			A(i)=B
		end do
		return
	end subroutine