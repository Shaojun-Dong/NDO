	subroutine ArrayDivideArrayFUNCNAME(Res,A,B,length)
		DATATYPE,intent(inout)::Res(:)
		class(*),intent(in)::A(:)
		class(*),intent(in)::B(:)
		integer,intent(in)::length
		integer::i
		if(size(Res).lt.length)then
			call writemess('ERROR in (/), length',-1)
			call error_stop
		end if
		if(size(B).lt.length)then
			call writemess('ERROR in (/), length',-1)
			call error_stop
		end if
		if(size(A).lt.length)then
			call writemess('ERROR in (/), length',-1)
			call error_stop
		end if
		select type(A)
			type is (integer)
				select type(B)
					type is (integer)
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (real(kind=4))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (real(kind=8))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (complex(kind=4))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (complex(kind=8))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					class default
						call Array_Divide_Array(Res,A,B,length)
				end select

			type is (real(kind=4))
				select type(B)
					type is (integer)
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (real(kind=4))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (real(kind=8))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (complex(kind=4))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (complex(kind=8))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					class default
						call Array_Divide_Array(Res,A,B,length)
				end select

			type is (real(kind=8))
				select type(B)
					type is (integer)
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (real(kind=4))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (real(kind=8))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (complex(kind=4))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (complex(kind=8))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					class default
						call Array_Divide_Array(Res,A,B,length)
				end select

			type is (complex(kind=4))
				select type(B)
					type is (integer)
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (real(kind=4))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (real(kind=8))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (complex(kind=4))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (complex(kind=8))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					class default
						call Array_Divide_Array(Res,A,B,length)
				end select

			type is (complex(kind=8))
				select type(B)
					type is (integer)
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (real(kind=4))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (real(kind=8))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (complex(kind=4))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					type is (complex(kind=8))
						do i=1,length
							Res(i)=A(i)/B(i)
						end do
					class default
						call Array_Divide_Array(Res,A,B,length)
				end select

			class default
				call writemess('ERROR class type')
				call error_stop
		end select
	end subroutine

