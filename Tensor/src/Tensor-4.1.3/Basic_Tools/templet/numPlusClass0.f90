	subroutine numPlusClassFUNCNAME(Res,A,B,length)
		DATATYPE,intent(inout)::Res(:)
		class(*),intent(in)::A,B(:)
		integer,intent(in)::length
		integer::i
		if(size(Res).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call error_stop
		end if
		if(size(B).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call error_stop
		end if
		select type(A)
			type is (integer)
				select type(B)
					type is (integer)
						Res(1:length)=A+B(1:length)
					type is (real(kind=4))
						Res(1:length)=A+B(1:length)
					type is (real(kind=8))
						Res(1:length)=A+B(1:length)
					type is (complex(kind=4))
						Res(1:length)=A+B(1:length)
					type is (complex(kind=8))
						Res(1:length)=A+B(1:length)
					type is (character(len=*))
						do i=1,length
							Res(i)=A+B(i)
						end do
					class default
						call num_Plus_Class(Res,A,B,length)
				end select

			type is (real(kind=4))
				select type(B)
					type is (integer)
						Res(1:length)=A+B(1:length)
					type is (real(kind=4))
						Res(1:length)=A+B(1:length)
					type is (real(kind=8))
						Res(1:length)=A+B(1:length)
					type is (complex(kind=4))
						Res(1:length)=A+B(1:length)
					type is (complex(kind=8))
						Res(1:length)=A+B(1:length)
					type is (character(len=*))
						do i=1,length
							Res(i)=A+B(i)
						end do
					class default
						call num_Plus_Class(Res,A,B,length)
				end select

			type is (real(kind=8))
				select type(B)
					type is (integer)
						Res(1:length)=A+B(1:length)
					type is (real(kind=4))
						Res(1:length)=A+B(1:length)
					type is (real(kind=8))
						Res(1:length)=A+B(1:length)
					type is (complex(kind=4))
						Res(1:length)=A+B(1:length)
					type is (complex(kind=8))
						Res(1:length)=A+B(1:length)
					type is (character(len=*))
						do i=1,length
							Res(i)=A+B(i)
						end do
					class default
						call writemess('ERROR class type')
						call error_stop
				end select

			type is (complex(kind=4))
				select type(B)
					type is (integer)
						Res(1:length)=A+B(1:length)
					type is (real(kind=4))
						Res(1:length)=A+B(1:length)
					type is (real(kind=8))
						Res(1:length)=A+B(1:length)
					type is (complex(kind=4))
						Res(1:length)=A+B(1:length)
					type is (complex(kind=8))
						Res(1:length)=A+B(1:length)
					type is (character(len=*))
						do i=1,length
							Res(i)=A+B(i)
						end do
					class default
						call num_Plus_Class(Res,A,B,length)
				end select

			type is (complex(kind=8))
				select type(B)
					type is (integer)
						Res(1:length)=A+B(1:length)
					type is (real(kind=4))
						Res(1:length)=A+B(1:length)
					type is (real(kind=8))
						Res(1:length)=A+B(1:length)
					type is (complex(kind=4))
						Res(1:length)=A+B(1:length)
					type is (complex(kind=8))
						Res(1:length)=A+B(1:length)
					type is (character(len=*))
						do i=1,length
							Res(i)=A+B(i)
						end do
					class default
						call num_Plus_Class(Res,A,B,length)
				end select
			type is (logical)
				select type(B)
					type is (character(len=*))
						do i=1,length
							Res(i)=A+B(i)
						end do
					class default
						call num_Plus_Class(Res,A,B,length)
				end select
			type is (character(len=*))
				select type(B)
					type is (integer)
						do i=1,length
							Res(i)=A+B(i)
						end do
					type is (real(kind=4))
						do i=1,length
							Res(i)=A+B(i)
						end do
					type is (real(kind=8))
						do i=1,length
							Res(i)=A+B(i)
						end do
					type is (complex(kind=4))
						do i=1,length
							Res(i)=A+B(i)
						end do
					type is (complex(kind=8))
						do i=1,length
							Res(i)=A+B(i)
						end do
					type is (logical)
						do i=1,length
							Res(i)=A+B(i)
						end do
					type is (character(len=*))
						do i=1,length
							Res(i)=A+B(i)
						end do
					class default
						call num_Plus_Class(Res,A,B,length)
				end select

			class default
				call num_Plus_Class(Res,A,B,length)
		end select
	end subroutine