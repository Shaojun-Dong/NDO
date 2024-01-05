	subroutine ClassPlusnumFUNCNAME(Res,A,B,length)
		DATATYPE,target,intent(inout)::Res(:)
		class(*),target,intent(in)::A(:),B
		integer,intent(in)::length
		integer::i
		if(size(Res).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call writemess('size(Res)='+size(Res),-1)
			call writemess('length='+length,-1)
			call error_stop
		end if
		if(size(A).lt.length)then
			call writemess('ERROR in (+), length',-1)
			call writemess('size(A)='+size(A),-1)
			call writemess('length='+length,-1)
			call error_stop
		end if
		select type(A)
			type is (integer)
				select type(B)
					type is (integer)
						Res(1:length)=A(1:length)+B
					type is (real(kind=4))
						Res(1:length)=A(1:length)+B
					type is (real(kind=8))
						Res(1:length)=A(1:length)+B
					type is (complex(kind=4))
						Res(1:length)=A(1:length)+B
					type is (complex(kind=8))
						Res(1:length)=A(1:length)+B
					type is (character(len=*))
						do i=1,length
							Res(i)=A(i)+B
						end do
					class default
						call Class_Plus_num(Res,A,B,length)
				end select

			type is (real(kind=4))
				select type(B)
					type is (integer)
						Res(1:length)=A(1:length)+B
					type is (real(kind=4))
						Res(1:length)=A(1:length)+B
					type is (real(kind=8))
						Res(1:length)=A(1:length)+B
					type is (complex(kind=4))
						Res(1:length)=A(1:length)+B
					type is (complex(kind=8))
						Res(1:length)=A(1:length)+B
					type is (character(len=*))
						do i=1,length
							Res(i)=A(i)+B
						end do
					class default
						call Class_Plus_num(Res,A,B,length)
				end select

			type is (real(kind=8))
				select type(B)
					type is (integer)
						Res(1:length)=A(1:length)+B
					type is (real(kind=4))
						Res(1:length)=A(1:length)+B
					type is (real(kind=8))
						Res(1:length)=A(1:length)+B
					type is (complex(kind=4))
						Res(1:length)=A(1:length)+B
					type is (complex(kind=8))
						Res(1:length)=A(1:length)+B
					type is (character(len=*))
						do i=1,length
							Res(i)=A(i)+B
						end do
					class default
						call Class_Plus_num(Res,A,B,length)
				end select

			type is (complex(kind=4))
				select type(B)
					type is (integer)
						Res(1:length)=A(1:length)+B
					type is (real(kind=4))
						Res(1:length)=A(1:length)+B
					type is (real(kind=8))
						Res(1:length)=A(1:length)+B
					type is (complex(kind=4))
						Res(1:length)=A(1:length)+B
					type is (complex(kind=8))
						Res(1:length)=A(1:length)+B
					type is (character(len=*))
						do i=1,length
							Res(i)=A(i)+B
						end do
					class default
						call Class_Plus_num(Res,A,B,length)
				end select

			type is (complex(kind=8))
				select type(B)
					type is (integer)
						Res(1:length)=A(1:length)+B
					type is (real(kind=4))
						Res(1:length)=A(1:length)+B
					type is (real(kind=8))
						Res(1:length)=A(1:length)+B
					type is (complex(kind=4))
						Res(1:length)=A(1:length)+B
					type is (complex(kind=8))
						Res(1:length)=A(1:length)+B
					type is (character(len=*))
						do i=1,length
							Res(i)=A(i)+B
						end do
					class default
						call Class_Plus_num(Res,A,B,length)
				end select

			type is (logical)
				select type(B)
					type is (character(len=*))
						do i=1,length
							Res(i)=A(i)+B
						end do
				class default
					call Class_Plus_num(Res,A,B,length)
				end select
			type is (character(len=*))
				select type(B)
					type is (integer)
						do i=1,length
							Res(i)=A(i)+B
						end do
					type is (real(kind=4))
						do i=1,length
							Res(i)=A(i)+B
						end do
					type is (real(kind=8))
						do i=1,length
							Res(i)=A(i)+B
						end do
					type is (complex(kind=4))
						do i=1,length
							Res(i)=A(i)+B
						end do
					type is (complex(kind=8))
						do i=1,length
							Res(i)=A(i)+B
						end do
					type is (logical)
						do i=1,length
							Res(i)=A(i)+B
						end do
					type is (character(len=*))
						do i=1,length
							Res(i)=A(i)+B
						end do
					class default
						call Class_Plus_num(Res,A,B,length)
				end select

			class default
				call Class_Plus_num(Res,A,B,length)
		end select
	end subroutine

